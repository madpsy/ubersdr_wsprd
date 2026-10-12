#!/usr/bin/env bash
# Build WSJT-X's jt9 and regenerate the expected outputs in tests/data/fst4w/
# (see the README there): the oracle for our FST4W decoder. jt9 is built by WSJT-X's own CMake in Docker
# (docker/wsjtx-oracle.Dockerfile: Ubuntu 24.04, gfortran 13.3, FFTW 3.3.10,
# as this machine), twice: with WSJT-X's Release flags and with -O2
# -ffp-contract=off; the script fails if the two disagree anywhere. Every
# case is run a second time with the jt9_wisdom.dat its first run left in the
# data dir (as a data dir kept from period to period has it), which must
# change nothing.
#
#   scripts/fst4w-oracle.sh [--build-only | --sim-only] [--no-sim]
#
#   WSJTX=~/repos/wsjtx            the WSJT-X tree (at REF_COMMIT)
#   FST4W_ORACLE_WORK=DIR          where to build and run (default: a temp
#                                  dir; a kept one makes rebuilds incremental)
#
# --no-sim keeps recordings/fst4w-sim/ as it is instead of regenerating it.
# --build-only builds jt9 and stops; --sim-only builds it and makes the
# inputs scripts/jt9-compare.sh runs (recordings/fst4w-sim/, the same file for
# file with the same libc, and recordings/fst4w-real/), nothing else: no
# spread set, no cases, no fixtures.
# Inputs are generated WAVs only (scripts/fst4w-sim.f90) plus recorded
# periods under recordings/; nothing is decoded live. recordings/ is not in
# git; cases on a missing input are skipped by the tests.
set -euo pipefail

REF_COMMIT=567ad29ce
IMAGE=wsjtx-oracle

# ================================================================ inside
# The container runs this script again with --inside: /src is the WSJT-X
# tree (read-only), /repo this repository, /work the work dir.
if [[ ${1:-} == --inside ]]; then
shift
build_only=0; sim=1; sim_only=0
while (( $# )); do
    case $1 in --build-only) build_only=1 ;; --no-sim) sim=0 ;; --sim-only) sim_only=1 ;; esac
    shift
done
repo=/repo; work=/work; src=/src
data=$repo/tests/data/fst4w
simdir=$repo/recordings/fst4w-sim

# ---------------------------------------------------------------- build
# WSJT-X's CMake, Release; the docs, man pages, tests and QMAP left out (they
# do not touch jt9). The second build overrides only the optimisation flags
# CMake's Release type sets: CMake/CompilerFlags.cmake still adds its own.
configure() {  # dir [cmake args]
    local dir=$1; shift
    mkdir -p "$dir"
    cmake -S "$src" -B "$dir" -DCMAKE_BUILD_TYPE=Release -DWSJT_GENERATE_DOCS=OFF \
        -DWSJT_SKIP_MANPAGES=ON -DWSJT_ENABLE_TESTS=OFF -DWSJT_SKIP_QMAP=ON "$@" > "$dir/configure.log" 2>&1
}
# libquadmath: gfortran's on x86-64, none on aarch64.
quadmath=; [[ $(gcc -print-file-name=libquadmath.so) == /* ]] && quadmath=-lquadmath
o2="-O2 -ffp-contract=off"
configure "$work/release"
configure "$work/o2" -DCMAKE_Fortran_FLAGS_RELEASE="$o2" -DCMAKE_C_FLAGS_RELEASE="$o2" \
    -DCMAKE_CXX_FLAGS_RELEASE="$o2"
cmake --build "$work/release" --target jt9 fst4sim -j"$(nproc)" > "$work/release/build.log" 2>&1
cmake --build "$work/o2" --target jt9 fst4sim -j"$(nproc)" > "$work/o2/build.log" 2>&1

# fst4w-sim against the Release build's libraries, as WSJT-X links fst4sim.
(
    cd "$work/release"
    mkdir -p sim
    read -r -a fflags <<< "$(sed -n 's/^Fortran_FLAGS = //p' CMakeFiles/fst4sim.dir/flags.make)"
    gfortran "${fflags[@]}" -Wno-intrinsic-shadow -I fortran_modules -I "$src/lib" -J sim \
        -c "$repo/scripts/fst4w-sim.f90" -o sim/fst4w-sim.o
    c++ -fopenmp -pthread sim/fst4w-sim.o -o fst4w-sim libwsjt_fort.a libwsjt_cxx.a -lm \
        /usr/lib/$(gcc -dumpmachine)/libboost_{log_setup,log,chrono,filesystem,regex,thread,atomic}.so.1.83.0 \
        -lfftw3f_threads -lfftw3f -lgfortran $quadmath
)
# The unpack and LDPC harnesses (scripts/fst4w-unpack-oracle.f90,
# fst4w-ldpc-oracle.f90), the same way, against each build's library (the
# harness itself compiled with the Release flags).
for b in release o2; do (
    cd "$work/$b"
    mkdir -p sim
    read -r -a fflags <<< "$(sed -n 's/^Fortran_FLAGS = //p' "$work/release/CMakeFiles/fst4sim.dir/flags.make")"
    for h in unpack ldpc frontend; do
        gfortran "${fflags[@]}" -I fortran_modules -I "$src/lib" -J sim -c "$repo/scripts/fst4w-$h-oracle.f90" \
            -o "sim/fst4w-$h-oracle.o"
        c++ -fopenmp -pthread "sim/fst4w-$h-oracle.o" -o "fst4w-$h-oracle" libwsjt_fort.a libwsjt_cxx.a -lm \
            /usr/lib/$(gcc -dumpmachine)/libboost_{log_setup,log,chrono,filesystem,regex,thread,atomic}.so.1.83.0 \
            -lfftw3f_threads -lfftw3f -lgfortran $quadmath
    done
) done
echo "fst4w-oracle: built in $work/release and $work/o2" >&2
[[ $build_only == 1 ]] && exit 0

# ---------------------------------------------------------------- inputs
# Names are YYMMDD_HHMM.wav: jt9 takes nutc from them. Periods start on
# their own boundaries (F5 on 5 minutes, F15 on 15, F30 on 30). Signals:
# scripts/fst4w-sim.f90 ("F:MESSAGE@f,dt,fspread,snr[,delay]" FST4W,
# "W:MESSAGE@f,t0,drift,snr" WSPR-2); jt9 looks at 1400-1600 Hz (-f 1500
# -F 100), which on a real band is WSPR's window too.
gen() { "$work/release/fst4w-sim" "$@" > /dev/null; }
slots=()
wspr_slots() {  # period message f drift snr: a WSPR signal every 2 minutes, added to slots
    local p=$1 m=$2 f=$3 d=$4 s=$5 t
    for (( t=1; t+111 <= p; t+=120 )); do slots+=("W:$m@$f,$t,$d,$s"); done
}
if [[ $sim == 1 ]]; then
    rm -rf "$simdir"; mkdir -p "$simdir"
    # FST4W-120 ------------------------------------------------------------
    # Types 1 and 2, Lorentzian fading, two WSPR signals in the window.
    gen -p 120 -S 1 -o "$simdir/f2-mix-261012_1200.wav" \
        "F:K1ABC FN42 37@1450,0.0,0,-24" "F:G4XYZ IO91 23@1520,0.5,-0.3,-26" \
        "F:PJ4/K1ABC 37@1550,-0.3,0,-22" "W:DL1ABC JO62 30@1480,1.0,0.5,-20" \
        "W:W9XYZ EN52 33@1530,1.0,0,-15"
    # SNR through the threshold (about -32.8 dB for FST4W-120).
    args=(); calls=(K1AA W2BB G3CC DL4DD JA5EE VK6FF ZL7GG PY8HH OH9II SM0JJ)
    snrs=(-27 -28 -29 -30 -31 -32 -33 -34 -35 -36)
    for i in "${!calls[@]}"; do args+=("F:${calls[$i]} FN42 37@$(( 1405 + i*19 )),$(( i % 3 - 1 )).0,0,${snrs[$i]}"); done
    gen -p 120 -S 2 -o "$simdir/f2-sweep-261012_1210.wav" "${args[@]}"
    # A busy WSPR band: ten WSPR signals across the window and two FST4W.
    args=(); wcalls=(K1WA W2WB G3WC DL4WD JA5WE VK6WF ZL7WG PY8WH OH9WI SM0WJ)
    for i in "${!wcalls[@]}"; do args+=("W:${wcalls[$i]} FN31 30@$(( 1408 + i*19 )),$(( i % 3 )).$(( i % 7 )),$(( i % 3 - 1 )),$(( -8 - i*2 ))"); done
    args+=("F:K1ABC FN42 37@1462,0.0,0,-28" "F:G4XYZ IO91 23@1547,0.3,0,-31")
    gen -p 120 -S 3 -o "$simdir/f2-busy-261012_1220.wav" "${args[@]}"
    # Doppler spread: Watterson 0.2-2 Hz (one with a 1 ms delay), Lorentzian.
    gen -p 120 -S 4 -o "$simdir/f2-fade-261012_1230.wav" \
        "F:K1AA FN42 37@1410,0,0.2,-24" "F:W2BB FN42 37@1435,0,0.5,-24" \
        "F:G3CC FN42 37@1460,0,1.0,-23" "F:DL4DD FN42 37@1485,0,2.0,-22" \
        "F:JA5EE FN42 37@1510,0,0.5,-24,1.0" "F:VK6FF FN42 37@1535,0,-0.1,-25" \
        "F:ZL7GG FN42 37@1560,0,-0.5,-24" "F:PY8HH FN42 37@1585,0,-1.5,-22"
    # The window's edges and start times off nominal.
    gen -p 120 -S 5 -o "$simdir/f2-edge-261012_1240.wav" \
        "F:K1AA FN42 37@1396,0,0,-20" "F:W2BB FN42 37@1402,0,0,-20" \
        "F:G3CC FN42 37@1597,0,0,-20" "F:DL4DD FN42 37@1603,0,0,-20" \
        "F:JA5EE FN42 37@1440,-1.5,0,-22" "F:VK6FF FN42 37@1480,-0.8,0,-22" \
        "F:ZL7GG FN42 37@1520,1.5,0,-22" "F:PY8HH FN42 37@1560,2.5,0,-22"
    gen -p 120 -S 6 -o "$simdir/f2-noise-261012_1250.wav"
    # Digital silence (all zeros): what a receiver gap gives the decoder.
    gen -p 120 -S 6 -N -o "$simdir/f2-silence-261012_1252.wav"
    # Hashed calls over three periods: type 2 full calls, then their type 3
    # "<CALL> GRID6" (and one no one sent the full call of), then more.
    gen -p 120 -S 7 -o "$simdir/f2-hash-261012_1300.wav" \
        "F:PJ4/K1ABC 37@1430,0,0,-22" "F:K9AN/12 33@1490,0,0,-23" "F:W1XY IO91 23@1550,0,0,-21"
    gen -p 120 -S 8 -o "$simdir/f2-hash-261012_1302.wav" \
        "F:<PJ4/K1ABC> FK52UD@1430,0,0,-22" "F:<G4UNK/P> IO91WM@1490,0,0,-23" "F:DL2ZZ JO62 0@1550,0,0,-21"
    gen -p 120 -S 9 -o "$simdir/f2-hash-261012_1304.wav" \
        "F:<K9AN/12> EN50WC@1430,0,0,-22" "F:K1ABC FN42 37@1490,0,0,-23" "F:<PJ4/K1ABC> FK52UD@1550,0,0,-24"
    # Type 3 of jt9's own default calls (-c K1ABC, -x W9XYZ, when neither
    # is given): fst4_decode's AP setup puts their hashes in the
    # table, so they resolve with no type 2 ever heard.
    gen -p 120 -S 24 -o "$simdir/f2-hash-261012_1306.wav" \
        "F:<K1ABC> FN42MI@1440,0,0,-21" "F:<W9XYZ> EN52AB@1520,0,0,-21"
    # FST4W-300 type 3s of the calls the FST4W-120 periods above sent in
    # full: a hash table carried from those periods resolves them; jt9 run
    # alone on the file does not.
    gen -p 300 -S 25 -o "$simdir/f5-hash-261012_1325.wav" \
        "F:<PJ4/K1ABC> FK52UD@1440,0,0,-26" "F:<K9AN/12> EN50WC@1520,0,0,-26"
    # FST4W on FST4's short periods: fst4_decode takes FST4W (iwspr=1) on
    # 15, 30 and 60 s too, and jt9's -p is 60 unless given (ours, src/jt9_main.cpp, takes it).
    gen -p 15 -S 26 -o "$simdir/p15-mix-261012_1410.wav" \
        "F:K1ABC FN42 37@1440,0,0,-14" "F:PJ4/K1ABC 37@1500,0.3,-0.5,-16" "F:G4XYZ IO91 23@1560,-0.2,0,-20"
    gen -p 30 -S 27 -o "$simdir/p30-mix-261012_1411.wav" \
        "F:K1ABC FN42 37@1430,0,0,-18" "F:K9AN/12 33@1500,0,1.0,-17" "F:G4XYZ IO91 23@1570,0.5,0,-23"
    gen -p 60 -S 28 -o "$simdir/p60-mix-261012_1412.wav" \
        "F:K1ABC FN42 37@1420,0,0,-20" "F:PJ4/K1ABC 37@1480,0,-0.3,-21" "F:W1XY IO91 23@1540,0,0,-25" \
        "F:G3CC FN42 37@1590,0.5,0,-28"
    # Below Keff=66's reach, for deep decoding's Keff=50 (calls already in
    # fst4w_calls.txt only).
    args=(); snrs=(-32.5 -33 -33.5 -34 -34.5 -35 -35.5 -36 -36.5 -37)
    for i in "${!calls[@]}"; do args+=("F:${calls[$i]} FN42 37@$(( 1405 + i*19 )),0,0,${snrs[$i]}"); done
    gen -p 120 -S 10 -o "$simdir/f2-weak-261012_1310.wav" "${args[@]}"
    # Near the threshold with starts between the 17 ms samples the sync
    # works in: some decode only one sample off isbest (jitter +1 or -1).
    args=()
    for i in "${!calls[@]}"; do
        dt=$(awk -v i="$i" 'BEGIN { printf "%.4f", (i*7 + 23) % 17 * 0.0017 }')
        snr=$(awk -v i="$i" 'BEGIN { printf "%.1f", -31.5 - (i % 4) * 0.5 }')
        args+=("F:${calls[$i]} FN42 37@$(( 1405 + i*19 )),$dt,0,$snr")
    done
    gen -p 120 -S 23 -o "$simdir/f2-jitter-261012_1312.wav" "${args[@]}"
    # FST4W-300 -----------------------------------------------------------
    slots=(); wspr_slots 300 "W2WB FN31 30" 1490 0 -18; wspr_slots 300 "OH9WI KP20 23" 1575 0.5 -12
    gen -p 300 -S 11 -o "$simdir/f5-mix-261012_1315.wav" \
        "F:K1ABC FN42 37@1430,0,0,-30" "F:PJ4/K1ABC 37@1470,0.4,0,-32" \
        "F:G4XYZ IO91 23@1510,-0.2,-0.3,-29" "F:DL4DD FN42 37@1550,0,0.5,-31" "${slots[@]}"
    args=(); snrs=(-36 -36.5 -37 -37.5 -38 -38.5 -39 -39.5 -40 -40.5)
    for i in "${!calls[@]}"; do args+=("F:${calls[$i]} FN42 37@$(( 1405 + i*19 )),0,0,${snrs[$i]}"); done
    gen -p 300 -S 12 -o "$simdir/f5-weak-261012_1320.wav" "${args[@]}"
    # FST4W-900 and FST4W-1800 --------------------------------------------
    slots=(); wspr_slots 900 "W2WB FN31 30" 1470 0 -16; wspr_slots 900 "OH9WI KP20 23" 1530 0 -20
    gen -p 900 -S 13 -o "$simdir/f15-mix-261012_1330.wav" \
        "F:K1ABC FN42 37@1440,0,0,-36" "F:G4XYZ IO91 23@1500,0,-0.1,-38" "F:JA5EE PM95 30@1560,0.6,0,-41" "${slots[@]}"
    slots=(); wspr_slots 1800 "W2WB FN31 30" 1470 0 -16
    gen -p 1800 -S 14 -o "$simdir/f30-mix-261012_1400.wav" \
        "F:K1ABC FN42 37@1440,0,0,-39" "F:G4XYZ IO91 23@1500,0,-0.05,-41" "F:JA5EE PM95 30@1560,0.6,0,-44" "${slots[@]}"
    (cd "$simdir" && sha256sum -- *.wav > SHA256SUMS)
fi
if [[ $sim == 1 && $sim_only == 0 ]]; then
    # The spread set (tests/test_fst4w_spread.cpp): 62 periods, eight
    # signals each, no spread, Watterson or Lorentzian spreads of 0.02-3 Hz,
    # SNRs from -16 dB (F2) down to -43 dB (F30), dt jitter.
    rm -rf "$repo/recordings/fst4w-spread"
    out=$repo/recordings/fst4w-spread
    mkdir -p $out
    gen=$work/release/fst4w-sim
    calls=(K1AA W2BB G3CC DL4DD JA5EE VK6FF ZL7GG PY8HH OH9II SM0JJ EA1KK I2LL F3MM ON4NN OZ5OO SP6PP)
    grids=(FN42 FN31 IO91 JO62 PM95 QF56 RF80 GG66 KP20 JO89 IN83 JN45 JN18 JO21 JO55 JO91)
    seed=100
    mk() {  # period count snrlo snrhi hh0
      local p=$1 n=$2 lo=$3 hi=$4 t0=$5
      for f in $(seq 0 $((n-1))); do
        seed=$((seed+1)); RANDOM=$seed
        args=()
        for s in $(seq 0 7); do
          c=${calls[$(( (f*8+s) % 16 ))]}; g=${grids[$(( (f*3+s) % 16 ))]}
          fr=$((1410 + s*22 + RANDOM % 6))
          dt=$(awk -v r=$RANDOM 'BEGIN{printf "%.2f", (r%200)/100.0-0.5}')
          kind=$((RANDOM % 3))
          sp=$(awk -v r=$RANDOM 'BEGIN{v=0.02*exp((r%1000)/1000.0*log(150)); printf "%.3f", v}')
          if [ $kind = 0 ]; then sp=0; elif [ $kind = 2 ]; then sp=-$sp; fi
          snr=$(awk -v r=$RANDOM -v lo=$lo -v hi=$hi 'BEGIN{printf "%.1f", lo+(r%1000)/1000.0*(hi-lo)}')
          args+=("F:$c $g $(( (RANDOM % 4) * 10 ))@$fr,$dt,$sp,$snr")
        done
        hh=$(( t0 + f*p/60 )); name=$(printf "s%d-%02d-261013_%02d%02d" $p $f $(( (hh/60)%24 )) $((hh%60)))
        $gen -p $p -S $seed -o $out/$name.wav "${args[@]}" > /dev/null
      done
    }
    mk 120 40 -30 -16 0
    mk 300 12 -36 -22 600
    mk 900 6 -40 -26 840
    mk 1800 4 -43 -29 960
    (cd "$out" && sha256sum -- *.wav > SHA256SUMS)
fi

# Real FST4W from WSJT-X's own samples: samples/FST4+FST4W/201230_0300.wav,
# an on-air FST4W-1800 period (DL0HOT JO60 30 at -44 dB), copied in as
# recordings/fst4w-real/ (the other sample there, 210115_0058.wav, is FST4-60,
# not FST4W).
mkdir -p "$repo/recordings/fst4w-real"
cp "$src/samples/FST4+FST4W/201230_0300.wav" "$repo/recordings/fst4w-real/201230_0300.wav"
if [[ $sim_only == 1 ]]; then
    echo "fst4w-oracle: inputs made in recordings/fst4w-sim and recordings/fst4w-real" >&2
    exit 0
fi

# ---------------------------------------------------------------- cases
# A case: name|calls seed|jt9's arguments after "-a ."|flags. The seed is
# fst4w_calls.txt put in the data dir first: a file under tests/, a case
# whose final fst4w_calls.txt is copied, or -. Flags: "bare" leaves out
# plotspec and decdata (without them jt9 measures no spread and writes no
# fst4_decodes.dat; the usual run makes both). Arguments starting
# recordings/ are this repository's files.
S=recordings/fst4w-sim
R40=recordings/40m-wspr-20261010-1914
k="-p 120 --fst4w -f 1500 -F 100"   # the usual options, plus -d
cases=()
for f in mix sweep busy fade edge noise; do
    inp=$(ls "$repo/$S" | grep "^f2-$f-")
    for d in 1 2 3; do cases+=("f2-$f/d$d|-|$k -d $d $S/$inp|"); done
done
cases+=("f2-mix/d3-bare|-|$k -d 3 $S/f2-mix-261012_1200.wav|bare")
# Silence: jt9's decoder skips a period whose first 15 s have an RMS under
# 0.5 (no decoded.txt written); alone, and between two periods that decode
# (decoded.txt is left as the first wrote it).
cases+=("f2-silence/d3|-|$k -d 3 $S/f2-silence-261012_1252.wav|")
cases+=("f2-silence/between-d3|-|$k -d 3 $S/f2-mix-261012_1200.wav $S/f2-silence-261012_1252.wav|")
h1=$S/f2-hash-261012_1300.wav h2=$S/f2-hash-261012_1302.wav h3=$S/f2-hash-261012_1304.wav
cases+=("f2-hash/1300-d3|-|$k -d 3 $h1|" "f2-hash/1302-d3|-|$k -d 3 $h2|" "f2-hash/1304-d3|-|$k -d 3 $h3|")
cases+=("f2-hash/1306-d3|-|$k -d 3 $S/f2-hash-261012_1306.wav|")
cases+=("f2-hash/sequence-d3|-|$k -d 3 $h1 $h2 $h3|")
cases+=("f2-hash/1302-d3-hash22|-|-p 120 -Y -f 1500 -F 100 -d 3 $h2|")
for d in 2 3; do cases+=("f2-jitter/d$d|-|$k -d $d $S/f2-jitter-261012_1312.wav|"); done
cases+=("f2-jitter/d3-seeded|tests/data/fst4w/calls-seed.txt|$k -d 3 $S/f2-jitter-261012_1312.wav|")
cases+=("f2-weak/d3|-|$k -d 3 $S/f2-weak-261012_1310.wav|")
cases+=("f2-weak/d3-seeded|tests/data/fst4w/calls-seed.txt|$k -d 3 $S/f2-weak-261012_1310.wav|")
cases+=("f2-weak/d2-seeded|tests/data/fst4w/calls-seed.txt|$k -d 2 $S/f2-weak-261012_1310.wav|")
k5="-p 300 --fst4w -f 1500 -F 100"
for d in 1 3; do cases+=("f5-mix/d$d|-|$k5 -d $d $S/f5-mix-261012_1315.wav|"); done
cases+=("f5-weak/d3-seeded|tests/data/fst4w/calls-seed.txt|$k5 -d 3 $S/f5-weak-261012_1320.wav|")
for d in 1 3; do cases+=("f15-mix/d$d|-|-p 900 --fst4w -f 1500 -F 100 -d $d $S/f15-mix-261012_1330.wav|"); done
for d in 1 3; do cases+=("f30-mix/d$d|-|-p 1800 --fst4w -f 1500 -F 100 -d $d $S/f30-mix-261012_1400.wav|"); done
for d in 1 3; do cases+=("f30-real-0300/d$d|-|-p 1800 --fst4w -f 1500 -F 100 -d $d recordings/fst4w-real/201230_0300.wav|"); done
# FST4W on 15, 30 and 60 s periods, and jt9's default period (60 s: no -p).
for p in 15 30 60; do
    inp=$(ls "$repo/$S" | grep "^p$p-mix-")
    for d in 1 3; do cases+=("p$p-mix/d$d|-|-p $p --fst4w -f 1500 -F 100 -d $d $S/$inp|"); done
done
cases+=("p60-mix/no-p-d3|-|--fst4w -f 1500 -F 100 -d 3 $S/p60-mix-261012_1412.wav|")
# Recorded 40 m WSPR periods (dial 7.0386): no FST4W, only the WSPR signals
# FST4W-120 tries and fails on, the busy-band CPU benchmark. Each alone and
# all five in one run.
seq=()
for t in 1914 1916 1918 1920 1922; do
    for d in 1 3; do cases+=("40m-$t/d$d|-|$k -d $d $R40/261010_$t.wav|"); done
    seq+=("$R40/261010_$t.wav")
done
cases+=("40m-sequence/d3|-|$k -d 3 ${seq[*]}|")

# jt9's outputs kept: stdout, stderr if any, the exit status,
# fst4_decodes.dat, decoded.txt, fst4w_calls.txt, and dopspread's spectra
# (fort.52, gzipped); timer.out (times) and jt9_wisdom.dat (FFTW's) are not.
run_case() {  # build case-spec rundir seeds-dir [wisdom file]
    local build=$1 spec=$2 dir=$3 seeds=$4 wis=${5:-}
    local name seed args flags
    IFS='|' read -r name seed args flags <<< "$spec"
    rm -rf "$dir"; mkdir -p "$dir"
    if [[ $seed == tests/* ]]; then cp "$repo/$seed" "$dir/fst4w_calls.txt"
    elif [[ $seed != - ]]; then cp "$seeds/$seed/fst4w_calls.txt" "$dir/fst4w_calls.txt"; fi
    [[ $seed != - ]] && cp "$dir/fst4w_calls.txt" "$dir/fst4w_calls.in"
    [[ $flags == *bare* ]] || touch "$dir/plotspec" "$dir/decdata"
    [[ -n $wis ]] && cp "$wis" "$dir/jt9_wisdom.dat"
    echo "$args" > "$dir/args.txt"
    [[ -n $flags ]] && echo "$flags" > "$dir/flags.txt"
    local argv=() tok
    for tok in $args; do
        case $tok in recordings/*|tests/*) argv+=("$repo/$tok") ;; *) argv+=("$tok") ;; esac
    done
    set +e
    (cd "$dir" && "$work/$build/jt9" -a . "${argv[@]}" > stdout.txt 2> stderr.txt)
    echo $? > "$dir/status.txt"
    set -e
    [[ -s $dir/stderr.txt ]] || rm -f "$dir/stderr.txt"
    [[ -f $dir/fort.52 ]] && gzip -n "$dir/fort.52"
    rm -f "$dir/timer.out" "$dir/plotspec" "$dir/decdata"
}

missing=0
for spec in "${cases[@]}"; do
    IFS='|' read -r name seed args flags <<< "$spec"
    for tok in $args; do
        if [[ $tok == recordings/* && ! -f $repo/$tok ]]; then
            echo "fst4w-oracle: $name needs $tok, missing" >&2; missing=1
        fi
    done
done
(( missing )) && exit 1

run_all() {  # build out-dir [wisdom-from-dir]
    local build=$1 out=$2 wfrom=${3:-} running=0 pass spec seed wis
    for pass in 1 2; do  # seeded cases after the cases they are seeded from
        for spec in "${cases[@]}"; do
            seed=${spec#*|}; seed=${seed%%|*}
            if [[ ($pass == 1 && ($seed == - || $seed == tests/*)) || ($pass == 2 && $seed != - && $seed != tests/*) ]]; then
                wis=; [[ -n $wfrom ]] && wis=$wfrom/${spec%%|*}/jt9_wisdom.dat
                if (( running >= $(nproc) )); then wait -n; running=$((running-1)); fi
                run_case "$build" "$spec" "$out/${spec%%|*}" "$out" "$wis" &
                running=$((running+1))
            fi
        done
        wait; running=0
    done
}
rm -rf "$work/run"
run_all release "$work/run/release"
run_all o2 "$work/run/o2"
run_all release "$work/run/wisdom" "$work/run/release"

# The builds must agree on every file but dopspread's spectra: Release's -O3
# has gen_fst4wave's cos/sin table (and fst4_baseline's log10/pow) computed
# by glibc's vector maths (libmvec), whose last bits differ from the scalar
# functions -O2 calls; fort.52 prints the spectrum to 1e-6 and shows that,
# nothing else has. Report how much it differs.
ok=1
if ! diff -r -x jt9_wisdom.dat -x fort.52.gz "$work/run/release" "$work/run/o2" > "$work/release-vs-o2.diff"; then
    echo "fst4w-oracle: the Release and $o2 builds differ: $work/release-vs-o2.diff" >&2
    ok=0
else
    echo "fst4w-oracle: Release and $o2 agree on all ${#cases[@]} cases (fort.52 aside)" >&2
fi
(
    cd "$work/run/release"
    for f in */*/fort.52.gz; do
        paste <(zcat "$f") <(zcat "../o2/$f") | awk '{ d = $2 - $4; if (d < 0) d = -d; if ($1 != $3) bad++; if (d > 0) n++; if (d > m) m = d }
            END { printf "%d %d %g\n", n + 0, bad + 0, m + 0 }'
    done | awk '{ n += $1; bad += $2; if ($3 > m) m = $3; files++ }
        END { printf "fst4w-oracle: fort.52, %d files: %d values differ (largest %g), %d frequencies\n", files, n, m, bad }' >&2
)
if ! diff -r -x jt9_wisdom.dat "$work/run/release" "$work/run/wisdom" > "$work/wisdom.diff"; then
    echo "fst4w-oracle: a second run with the first's jt9_wisdom.dat differs: $work/wisdom.diff" >&2
    ok=0
else
    echo "fst4w-oracle: a second run with the first's jt9_wisdom.dat changes nothing" >&2
fi
(( ok )) || exit 1

# ---------------------------------------------------------------- unpack
# The unpack vectors: tests/data/fst4w/unpack/*.ops (written outside by
# scripts/fst4w-unpack-vectors.py), through packjt77 itself.
for ops in "$data"/unpack/*.ops; do
    "$work/release/fst4w-unpack-oracle" < "$ops" > "${ops%.ops}.expected"
done
echo "fst4w-oracle: unpack vectors: $(cat "$data"/unpack/*.ops | grep -c '^U') payloads," \
    "$(cat "$data"/unpack/*.expected | grep -c '^1|') unpacked" >&2

# ---------------------------------------------------------------- ldpc
# The LDPC vectors: tests/data/fst4w/ldpc/vectors.ops (written outside by
# scripts/fst4w-ldpc-vectors.py), through decode240_74 and fastosd240_74
# themselves; and the -O2 build's, which must agree (its float code is not
# vectorised differently, but check).
"$work/release/fst4w-ldpc-oracle" < "$data/ldpc/vectors.ops" > "$data/ldpc/vectors.expected"
"$work/o2/fst4w-ldpc-oracle" < "$data/ldpc/vectors.ops" > "$work/ldpc-o2.expected"
if ! cmp -s "$data/ldpc/vectors.expected" "$work/ldpc-o2.expected"; then
    echo "fst4w-oracle: the -O2 build's LDPC decoding differs: $work/ldpc-o2.expected" >&2
    exit 1
fi
echo "fst4w-oracle: ldpc vectors: $(grep -c '^R' "$data/ldpc/vectors.expected") decodes," \
    "$(grep -c '^R [1-9]' "$data/ldpc/vectors.expected") decoded" \
    "($(grep -c '^R 1 ' "$data/ldpc/vectors.expected") by belief propagation)," \
    "$(grep -c '^S' "$data/ldpc/vectors.expected") OSD alone" >&2

# ---------------------------------------------------------------- frontend
# The front end's traces: each period through fst4_decode's own front-end
# routines (scripts/fst4w-frontend-oracle.f90), tests/data/fst4w/frontend/
# <name>.trace, listed in cases.txt (name, input, period, depth). The -O2
# build's traces are compared too, and differ only where -O3's libmvec
# baseline does (reported, not an error).
fe=$data/frontend
mkdir -p "$fe"
: > "$fe/cases.txt"
for f in "$simdir"/*.wav; do
    b=$(basename "$f" .wav)
    case $b in f2-*) p=120 ;; f5-*) p=300 ;; f15-*) p=900 ;; f30-*) p=1800 ;; esac
    echo "$b recordings/fst4w-sim/$b.wav $p 3" >> "$fe/cases.txt"
done
for t in 1914 1916 1918 1920 1922; do
    echo "40m-$t recordings/40m-wspr-20261010-1914/261010_$t.wav 120 3" >> "$fe/cases.txt"
done
echo "f30-real-0300 recordings/fst4w-real/201230_0300.wav 1800 3" >> "$fe/cases.txt"
sort -u -o "$fe/cases.txt" "$fe/cases.txt"
rm -f "$fe"/*.trace
running=0
while read -r name input p d; do
    if (( running >= $(nproc) )); then wait -n; running=$((running-1)); fi
    ( "$work/release/fst4w-frontend-oracle" "$repo/$input" "$p" "$d" > "$fe/$name.trace"
      "$work/o2/fst4w-frontend-oracle" "$repo/$input" "$p" "$d" > "$work/$name-o2.trace" ) &
    running=$((running+1))
done < "$fe/cases.txt"
wait
ndiff=0
while read -r name input p d; do
    cmp -s "$fe/$name.trace" "$work/$name-o2.trace" || ndiff=$((ndiff+1))
done < "$fe/cases.txt"
echo "fst4w-oracle: frontend traces: $(wc -l < "$fe/cases.txt") periods, $(cat "$fe"/*.trace | grep -c '^J') bit-metric frames;" \
    "the -O2 build's differ on $ndiff (libmvec)" >&2

# ---------------------------------------------------------------- fixtures
rm -rf "$data/cases"
for spec in "${cases[@]}"; do
    name=${spec%%|*}
    mkdir -p "$data/cases/$name"
    (cd "$work/run/release/$name" && cp -- $(ls | grep -v '^jt9_wisdom.dat$') "$data/cases/$name/")
done
{
    echo "WSJT-X $(git -C "$src" rev-parse --short=9 HEAD 2>/dev/null || echo "$REF_COMMIT")"
    gfortran --version | head -1
    echo "FFTW $(dpkg-query -W -f='${Version}' libfftw3-single3) ($(uname -m))"
    echo "Release Fortran: $(sed -n 's/^Fortran_FLAGS = //p' "$work/release/CMakeFiles/wsjt_fort_omp.dir/flags.make")"
    echo "Release C:       $(sed -n 's/^C_FLAGS = //p' "$work/release/CMakeFiles/wsjt_cxx.dir/flags.make")"
    echo "O2 Fortran:      $(sed -n 's/^Fortran_FLAGS = //p' "$work/o2/CMakeFiles/wsjt_fort_omp.dir/flags.make")"
} > "$data/BUILD.txt"
echo "fst4w-oracle: ${#cases[@]} cases in $data/cases" >&2

# ---------------------------------------------------------------- coverage
# fst4_decodes.dat's fields (format 3021): nutc icand itry nsyncoh iaptype
# ijitter npct ntype Keff nsync_qual nharderrors dmin nhp hd sync xsnr xdt
# fsig w50 message.
cat "$data"/cases/*/*/fst4_decodes.dat | awk '
    { keff[$9]++; jit[$6]++; try[$3]++; n++ }
    END {
        printf "decodes %d; Keff:", n; for (k in keff) printf " %s=%d", k, keff[k]
        printf "; jitter:"; for (k in jit) printf " %s=%d", k, jit[k]
        printf "; try:"; for (k in try) printf " %s=%d", k, try[k]; printf "\n"
    }' >&2
for d in 1 2 3; do
    printf 'depth %s: %s decodes\n' $d "$(cat "$data"/cases/*/d$d*/stdout.txt 2>/dev/null | grep -c '^[0-9]' || true)" >&2
done
resolved=$(grep -h '<[A-Z0-9/]*>' "$data"/cases/*/*/stdout.txt | grep -vc '<\.\.\.>' || true)
echo "hashed calls resolved: $resolved" >&2
gaps=()
for kf in 66 50; do awk -v k=$kf '$9==k{f=1} END{exit !f}' "$data"/cases/*/*/fst4_decodes.dat || gaps+=("Keff $kf"); done
for j in 0 1 2; do awk -v j=$j '$6==j{f=1} END{exit !f}' "$data"/cases/*/*/fst4_decodes.dat || gaps+=("jitter $j"); done
(( resolved > 0 )) || gaps+=("a resolved hash")
if (( ${#gaps[@]} )); then
    echo "fst4w-oracle: not covered: ${gaps[*]}" >&2
    exit 1
fi
echo "fst4w-oracle: done" >&2
exit 0
fi

# ================================================================ outside
repo=$(cd "$(dirname "$0")/.." && pwd)
wsjtx=${WSJTX:-$HOME/repos/wsjtx}
work=${FST4W_ORACLE_WORK:-$(mktemp -d "${TMPDIR:-/tmp}/fst4w-oracle.XXXXXX")}
mkdir -p "$work"
work=$(cd "$work" && pwd)
for a in "$@"; do
    case $a in --build-only|--no-sim|--sim-only) ;; *) echo "usage: $0 [--build-only | --sim-only] [--no-sim]" >&2; exit 2 ;; esac
done
full=1
for a in "$@"; do [[ $a == --build-only || $a == --sim-only ]] && full=0; done
head=$(git -C "$wsjtx" rev-parse --short=9 HEAD)
if [[ $head != "$REF_COMMIT" ]]; then
    echo "fst4w-oracle: $wsjtx is at $head, not $REF_COMMIT" >&2
    exit 1
fi
if [[ -n $(git -C "$wsjtx" status --porcelain -- lib CMake CMakeLists.txt) ]]; then
    echo "fst4w-oracle: $wsjtx has local changes in lib/ or its CMake files" >&2
    exit 1
fi
if (( full )); then  # the vectors the full run's harnesses read
    python3 "$repo/scripts/fst4w-unpack-vectors.py" "$repo/tests/data/fst4w/unpack"
    mkdir -p "$repo/tests/data/fst4w/ldpc"
    python3 "$repo/scripts/fst4w-ldpc-vectors.py" "$repo/tests/data/fst4w/ldpc/vectors.ops" "$wsjtx"
fi
docker build -q -t "$IMAGE" -f "$repo/docker/wsjtx-oracle.Dockerfile" "$repo/docker" > /dev/null
echo "fst4w-oracle: working in $work" >&2
exec docker run --rm -u "$(id -u):$(id -g)" -e HOME=/work \
    -v "$wsjtx:/src:ro" -v "$repo:/repo" -v "$work:/work" \
    "$IMAGE" bash /repo/scripts/fst4w-oracle.sh --inside "$@"
