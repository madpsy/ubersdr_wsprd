#!/usr/bin/env bash
# Rebuild WSJT-X's wsprd by hand and regenerate the expected outputs in
# tests/data/wsprd/ (see the README there): the oracle for the in-process
# WSPR decoder. Builds the reference twice, with WSJT-X's own Release flags
# and with -O2 -ffp-contract=off, and fails if the two disagree anywhere.
#
#   scripts/wsprd-oracle.sh [--build-only] [--no-sim] [--wisdom FILE]
#
#   WSJTX=~/repos/wsjtx            the WSJT-X tree (lib/wsprd at REF_COMMIT)
#   WSPRD_ORACLE_WORK=/tmp/...     where to build and run (default: a temp dir)
#
# --no-sim keeps recordings/wsprd-sim/ as it is instead of regenerating it.
# --wisdom FILE: run every case once more with FILE as the data dir's
# wspr_wisdom.dat (which wsprd imports before planning) and report the cases
# whose outputs change; the fixtures are not touched.
# Inputs under recordings/ are not in git; cases on them are skipped by the
# tests when the file is missing.
set -euo pipefail

REF_COMMIT=567ad29ce
repo=$(cd "$(dirname "$0")/.." && pwd)
wsjtx=${WSJTX:-$HOME/repos/wsjtx}
work=${WSPRD_ORACLE_WORK:-$(mktemp -d "${TMPDIR:-/tmp}/wsprd-oracle.XXXXXX")}
build_only=0
sim=1
wisdom=
while (( $# )); do
    case $1 in
        --build-only) build_only=1 ;;
        --no-sim) sim=0 ;;
        --wisdom) wisdom=$(realpath "$2"); shift ;;
        *) echo "usage: $0 [--build-only] [--no-sim] [--wisdom FILE]" >&2; exit 2 ;;
    esac
    shift
done

head=$(git -C "$wsjtx" rev-parse --short=9 HEAD)
if [[ $head != "$REF_COMMIT" ]]; then
    echo "wsprd-oracle: $wsjtx is at $head, not $REF_COMMIT" >&2
    exit 1
fi
if [[ -n $(git -C "$wsjtx" status --porcelain -- lib/wsprd lib/indexx.f90 lib/lookup3.c lib/lookup3.h lib/init_random_seed.c lib/init_random_seed.h) ]]; then
    echo "wsprd-oracle: $wsjtx has local changes in the wsprd sources" >&2
    exit 1
fi

# ---------------------------------------------------------------- build
# CMake/Sources.cmake wsprd_CSRCS + lib/indexx.f90 + lib/wsprd/osdwspr.f90;
# wsprsim_CSRCS for wsprsim.
src=$work/src
rm -rf "$src"
mkdir -p "$src/lib/wsprd" "$work/bin"
cp "$wsjtx"/lib/wsprd/*.[ch] "$wsjtx"/lib/wsprd/*.f90 "$src/lib/wsprd/"
rm -f "$src/lib/wsprd/fftw3.h"  # a stale copy; CMake builds against the system's
cp "$wsjtx"/lib/{indexx.f90,lookup3.c,lookup3.h,init_random_seed.c,init_random_seed.h} "$src/lib/"

wsprd_c=(lib/lookup3.c lib/wsprd/wsprd.c lib/wsprd/wsprsim_utils.c lib/wsprd/wsprd_utils.c
         lib/wsprd/fano.c lib/wsprd/jelinek.c lib/wsprd/tab.c lib/wsprd/nhash.c
         lib/init_random_seed.c lib/wsprd/wsprd_stream.c)
wsprd_f=(lib/indexx.f90 lib/wsprd/osdwspr.f90)
wsprsim_c=(lib/lookup3.c lib/wsprd/wsprsim.c lib/wsprd/wsprsim_utils.c lib/wsprd/wsprd_utils.c
           lib/wsprd/fano.c lib/wsprd/tab.c lib/wsprd/nhash.c)

# CMAKE_BUILD_TYPE=Release (the default) with gcc/gfortran, from
# CMake/CompilerFlags.cmake and CMake's GNU defaults (-O3 -DNDEBUG); C gets
# -fopenmp as CMake finds OpenMP; the link -Wl,--gc-sections.
release_c=(-Wall -Wextra -fopenmp -pthread -O3 -DNDEBUG -fdata-sections -ffunction-sections)
release_f=(-Werror=trampolines -Wa,--noexecstack -O3 -DNDEBUG -funroll-loops -fno-f2c
           -ffpe-summary=invalid,zero,overflow,underflow -Wall -Wno-conversion
           -Wno-c-binding-type -fno-second-underscore)
release_l=(-fopenmp -pthread -Wl,--gc-sections)
o2_c=(-O2 -ffp-contract=off -w)
o2_f=(-O2 -ffp-contract=off -fno-second-underscore -w)
o2_l=()

build() {  # name cflags-array fflags-array lflags-array [extra -D...]
    local name=$1 cf=$2[@] ff=$3[@] lf=$4[@]
    local cflags=("${!cf}") fflags=("${!ff}") lflags=("${!lf}")
    shift 4
    local obj=$work/obj-$name
    rm -rf "$obj"; mkdir -p "$obj"
    local objs=() f
    for f in "${wsprd_c[@]}"; do
        gcc "${cflags[@]}" "$@" -c "$src/$f" -o "$obj/$(basename "$f").o" 2>>"$work/build-$name.log"
        objs+=("$obj/$(basename "$f").o")
    done
    for f in "${wsprd_f[@]}"; do
        gfortran "${fflags[@]}" -J "$obj" -c "$src/$f" -o "$obj/$(basename "$f").o" 2>>"$work/build-$name.log"
        objs+=("$obj/$(basename "$f").o")
    done
    gfortran "${lflags[@]}" "${objs[@]}" -o "$work/bin/wsprd-$name" -lfftw3f -lm 2>>"$work/build-$name.log"
}
: > "$work/build-release.log"; : > "$work/build-o2.log"
build release release_c release_f release_l
build o2 o2_c o2_f o2_l

simobj=$work/obj-sim; rm -rf "$simobj"; mkdir -p "$simobj"
objs=()
for f in "${wsprsim_c[@]}"; do
    gcc "${release_c[@]}" -w -c "$src/$f" -o "$simobj/$(basename "$f").o"
    objs+=("$simobj/$(basename "$f").o")
done
gcc "${release_l[@]}" "${objs[@]}" -o "$work/bin/wsprsim" -lm
objs=()
for f in "${wsprsim_c[@]}"; do
    [[ $f == lib/wsprd/wsprsim.c ]] && continue
    objs+=("$simobj/$(basename "$f").o")
done
gcc -O2 -w -I"$src/lib/wsprd" "$repo/scripts/wsprd-sim.c" "${objs[@]}" -o "$work/bin/wsprd-sim" -lm
echo "wsprd-oracle: built in $work/bin" >&2
[[ $build_only == 1 ]] && exit 0


# ---------------------------------------------------------------- inputs
# recordings/wsprd-sim/: WSJT-X's own sample plus synthetic files (wsprd-sim,
# seeded, so the same bytes each time with the same libm). Names end in
# YYMMDD_HHMM.c2/.wav: wsprd takes the date and time from them. Signals are
# "MESSAGE@f,t0,drift,snr" (scripts/wsprd-sim.c).
simdir=$repo/recordings/wsprd-sim
data=$repo/tests/data/wsprd
mkdir -p "$simdir"
gen() { "$work/bin/wsprd-sim" "$@"; }
sweep_calls=(K1ABC W2DEF G3GHI DL1JKL JA2MNO VK3PQR ZL4STU PY5VWX OH6YZA SM7BCD LU8EFG EA9HIJ)
if [[ $sim == 1 ]]; then
    cp "$wsjtx/samples/WSPR/150426_0918.wav" "$simdir/150426_0918.wav"
    # Every call of the other files, strong: its hashtable.txt seeds the
    # cases that need calls already known (OSD only accepts known calls).
    args=()
    for i in "${!sweep_calls[@]}"; do args+=("${sweep_calls[$i]} FN42 37@$(( -99 + i*9 )),1.0,0,-10"); done
    args+=("PJ4/K1ABC 37@9,1.0,0,-10" "G4XYZ/P 23@18,1.0,0,-10" "K9AN/12 33@27,1.0,0,-10"
           "W1XY IO91 23@36,1.0,0,-10" "DL2ZZ JO62 0@45,1.0,0,-10" "JA1XYZ PM95 60@54,1.0,0,-10"
           "N0XYZ EN34 30@63,1.0,0,-10" "VE3ABC FN03 20@72,1.0,0,-10" "OK1AB JN79 27@81,1.0,0,-10")
    gen -S 1 -F 14.0956 -o "$simdir/sim-strong-261010_1158.c2" "${args[@]}"
    # Type 1, 2 and 3 messages: prefix and suffix calls, hashed calls with
    # 6-character grids, a hash no one has sent (<...>), a 4-character grid
    # in a type 3.
    gen -S 2 -F 14.0956 -o "$simdir/sim-types-261010_1200.c2" \
        "K1ABC FN42 37@-96,1.0,0,-18" "PJ4/K1ABC 37@-78,0.5,0,-16" \
        "<PJ4/K1ABC> FK52UD 37@-60,1.5,0,-19" "G4XYZ/P 23@-42,1.2,0,-20" \
        "K9AN/12 33@-24,1.0,0,-21" "<G4XYZ/P> IO91WM 23@-6,1.0,0,-22" \
        "<W1NEW> FN31PR 30@12,1.0,0,-15" "<PA0SKT/2> JO33 40@30,1.0,0,-17" \
        "DL2ZZ JO62 0@48,0.0,0,-20" "JA1XYZ PM95 60@66,2.0,0,-23" \
        "<K1ABC> FN42MI 37@84,1.0,0,-24" "VE3/W1AW 30@102,1.0,0,-19"
    # SNR sweep through the Fano threshold into OSD's reach and past it.
    args=(); snrs=(-22 -24 -26 -27 -28 -29 -30 -31 -32 -33 -34 -35)
    for i in "${!sweep_calls[@]}"; do args+=("${sweep_calls[$i]} FN42 37@$(( -88 + i*16 )),1.0,0,${snrs[$i]}"); done
    gen -S 7 -F 14.0956 -o "$simdir/sim-sweep-261010_1202.c2" "${args[@]}"
    args=(); snrs=(-26.5 -27 -27.5 -28 -28.5 -29 -29.5 -30 -30.5 -31 -31.5 -32)
    for i in "${!sweep_calls[@]}"; do args+=("${sweep_calls[$i]} FN42 37@$(( -88 + i*16 )),$(( i % 3 )).2,0,${snrs[$i]}"); done
    gen -S 8 -F 10.1387 -o "$simdir/sim-sweep2-261010_1204.c2" "${args[@]}"
    # Drift and start time.
    gen -S 3 -F 7.0386 -o "$simdir/sim-drift-261010_1206.c2" \
        "K1ABC FN42 37@-95,-1.0,-4,-20" "W2DEF FN42 37@-75,-0.3,-3,-21" \
        "G3GHI FN42 37@-55,0.4,-2.5,-22" "DL1JKL FN42 37@-35,1.0,-1,-23" \
        "JA2MNO FN42 37@-15,1.7,0.7,-24" "VK3PQR FN42 37@5,2.3,1.5,-25" \
        "ZL4STU FN42 37@25,2.9,2,-26" "PY5VWX FN42 37@45,3.5,3,-22" \
        "OH6YZA FN42 37@65,0.8,4,-21" "SM7BCD FN42 37@85,1.3,3.3,-27"
    # Crowded: 6-8 Hz apart, a call twice 3 Hz apart (the dupe check), two
    # calls 2 Hz apart, strong and weak side by side.
    args=(); calls=("${sweep_calls[@]}" W1XY DL2ZZ JA1XYZ N0XYZ VE3ABC OK1AB G0AAA F5BBB I2CCC SP3DDD HA4EEE YO5FFF)
    for i in $(seq 0 23); do
        f=$(( -84 + i*7 )); s=$(( -12 - (i*7) % 17 ))
        args+=("${calls[$i]} JO01 33@$f,$(( i % 4 )).0,$(( i % 3 - 1 )),$s")
    done
    args+=("K1ABC JO01 33@-81,1.0,0,-18" "XE1GGG EK09 10@-60.5,1.0,0,-20" "CE2HHH FF46 13@-58.5,1.0,0,-20")
    gen -S 4 -F 14.0956 -o "$simdir/sim-crowd-261010_1210.c2" "${args[@]}"
    # Outside +/-110 Hz (only -w reaches them) and just inside.
    gen -S 5 -F 3.5686 -o "$simdir/sim-wide-261010_1212.c2" \
        "K1ABC FN42 37@-148,1.0,0,-18" "W2DEF FN42 37@-131,1.0,0,-19" "G3GHI FN42 37@-116,1.0,0,-20" \
        "DL1JKL FN42 37@-108,1.0,0,-21" "JA2MNO FN42 37@0,1.0,0,-20" "VK3PQR FN42 37@107,1.0,0,-22" \
        "ZL4STU FN42 37@118,1.0,0,-20" "PY5VWX FN42 37@137,1.0,0,-19" "OH6YZA FN42 37@149,1.0,0,-18"
    # Starts off either end of the file (partial decodes), very strong.
    gen -S 6 -F 1.8366 -o "$simdir/sim-edge-261010_1214.c2" \
        "K1ABC FN42 37@-90,-3.0,0,-15" "W2DEF FN42 37@-60,-1.5,0,-18" "G3GHI FN42 37@-30,4.5,0,-18" \
        "DL1JKL FN42 37@0,6.5,0,-15" "JA2MNO FN42 37@30,1.0,0,5" "VK3PQR FN42 37@60,8.0,0,-12" \
        "ZL4STU FN42 37@90,1.0,0,-3"
    gen -S 9 -F 28.1246 -o "$simdir/sim-noise-261010_1216.c2"
    # WSPR-15 as a .c2 (type 15 in the header): f in audio Hz, as wsprd shows.
    gen -S 10 -m -F 0.4742 -o "$simdir/sim-wspr15-261010_1215.c2" \
        "K1ABC FN42 37@-10,1.0,0,-30" "G3GHI IO91 23@-4,3.0,0.5,-33" "DL1JKL JO62 30@3,1.0,0,-36" \
        "<PJ4/K1ABC> FK52UD 37@9,0.0,0,-31"
    # 12 kHz audio, through readwavfile's FFTs. WSPR-15 is 1/8 the width:
    # wsprd looks within +/-13.75 Hz of 1612.5 Hz.
    gen -S 11 -F 7.0386 -W "$simdir/sim-wav2-261010_1218.wav" \
        "K1ABC FN42 37@-90,1.0,0,-20" "PJ4/K1ABC 37@-50,0.5,1,-22" "<PJ4/K1ABC> FK52UD 37@-10,1.5,-1,-24" \
        "G3GHI IO91 23@20,2.0,0,-26" "DL1JKL JO62 30@50,1.0,2,-28" "JA2MNO PM95 10@80,0.2,0,-30"
    gen -S 12 -m -F 0.4742 -W "$simdir/sim-wav15-261010_1230.wav" \
        "K1ABC FN42 37@-10,1.0,0,-30" "G3GHI IO91 23@-2.5,3.0,0.2,-34" "DL1JKL JO62 30@6,1.0,0,-37"
    (cd "$simdir" && sha256sum -- *.c2 *.wav > SHA256SUMS)
fi
# The few small enough for git, so the tests run without recordings/.
mkdir -p "$data/inputs"
rm -f "$data"/inputs/*.c2 "$data"/inputs/*.wav
for f in sim-types-261010_1200.c2 sim-sweep-261010_1202.c2 sim-crowd-261010_1210.c2 sim-drift-261010_1206.c2; do
    cp "$simdir/$f" "$data/inputs/$f"
done
# .wav files cut short: 10000 samples and a half (fread keeps the half's
# low byte, 0x9a, in the buffer), and 5 samples (fewer than the header's 22
# shorts, which stay in the buffer after them).
head -c $(( 44 + 20001 )) "$simdir/sim-wav2-261010_1218.wav" > "$data/inputs/cut-odd-261010_1218.wav"
head -c 54 "$simdir/sim-wav2-261010_1218.wav" > "$data/inputs/cut-tiny-261010_1218.wav"
(cd "$data/inputs" && sha256sum -- *.c2 *.wav > SHA256SUMS)

# ---------------------------------------------------------------- cases
# A case: name|hashtable seed (a case whose final hashtable.txt is copied in
# first, a repository file under tests/, or -)|runs, ';'-separated, each wsprd's argv without -a[|files].
# Every case runs in a fresh data dir that is also the cwd (-c writes
# there). files, ';'-separated, are put in the data dir first: "SRC NAME"
# copies SRC (repository-relative; - for an empty file) to NAME, "SRC NAME
# ro" makes it read-only too; the case's files.txt lists them, one a line.
dials() {  # input -> its dial in MHz (a .c2 carries its own; -f is given anyway)
    case $1 in
        *150426_0918*|*sweep-*|*types*|*crowd*|*strong*) echo 14.0956 ;;
        *sweep2*) echo 10.1387 ;;
        *drift*|*wav2*) echo 7.0386 ;;
        *wide*) echo 3.5686 ;;
        *edge*) echo 1.8366 ;;
        *noise*) echo 28.1246 ;;
        *wspr15*|*wav15*) echo 0.4742 ;;
    esac
}
T=tests/data/wsprd/inputs
Rd=recordings/wsprd-sim
optsets=(
    "skim|-c -C 500 -o 4 -d"
    "skim-nodeep|-c -C 500 -o 4"
    "default|"
    "wide|-w"
    "quick|-q"
    "single|-s"
    "osd0|-o 0"
    "nohash|-H -o 4"
    "noblock|-B"
    "stack|-J"
    "bias|-z 0.6"
    "dialerr|-e 2.5 -d"
    "verbose|-v -d"
)
cases=()
for inp in $T/sim-types-261010_1200.c2 $T/sim-sweep-261010_1202.c2 $T/sim-crowd-261010_1210.c2 \
          $T/sim-drift-261010_1206.c2 $Rd/150426_0918.wav $Rd/sim-sweep2-261010_1204.c2 \
          $Rd/sim-wide-261010_1212.c2 $Rd/sim-edge-261010_1214.c2 $Rd/sim-noise-261010_1216.c2 \
          $Rd/sim-wspr15-261010_1215.c2 $Rd/sim-strong-261010_1158.c2 $Rd/sim-wav2-261010_1218.wav \
          $Rd/sim-wav15-261010_1230.wav; do
    stem=$(basename "$inp"); stem=${stem%.*}; stem=${stem%-??????_????}
    [[ $stem == 150426_0918 ]] && stem=sample
    m=; [[ $inp == *wav15* ]] && m="-m "
    for o in "${optsets[@]}"; do
        cases+=("$stem/${o%%|*}|-|${m}${o#*|} -f $(dials "$inp") $inp")
    done
done
# OSD accepts only calls already in the hashtable: the sweeps again, with
# the strong file's hashtable.txt as the data dir's.
cases+=("sim-sweep/skim-seeded|sim-strong/skim|-c -C 500 -o 4 -d -f 14.0956 $T/sim-sweep-261010_1202.c2")
cases+=("sim-sweep/osd0-seeded|sim-strong/skim|-o 0 -f 14.0956 $T/sim-sweep-261010_1202.c2")
cases+=("sim-sweep/nohash-seeded|sim-strong/skim|-H -o 4 -f 14.0956 $T/sim-sweep-261010_1202.c2")
cases+=("sim-sweep2/skim-seeded|sim-strong/skim|-c -C 500 -o 4 -d -f 10.1387 $Rd/sim-sweep2-261010_1204.c2")
cases+=("sim-sweep2/stack-seeded|sim-strong/skim|-J -o 4 -d -f 10.1387 $Rd/sim-sweep2-261010_1204.c2")
cases+=("sim-types/skim-seeded|sim-strong/skim|-c -C 500 -o 4 -d -f 14.0956 $T/sim-types-261010_1200.c2")
# One data dir across periods, as a station decoding period after period
# keeps it: hashtable.txt, ALL_WSPR.TXT and
# wspr_wisdom.dat carry over from run to run.
k="-c -C 500 -o 4 -d -f"
cases+=("sequence/skim-tsc|-|$k 14.0956 $T/sim-types-261010_1200.c2;$k 14.0956 $T/sim-sweep-261010_1202.c2;$k 14.0956 $T/sim-crowd-261010_1210.c2")
cases+=("sequence/skim-strong-first|-|$k 14.0956 $Rd/sim-strong-261010_1158.c2;$k 14.0956 $T/sim-sweep-261010_1202.c2;$k 10.1387 $Rd/sim-sweep2-261010_1204.c2;$k 0.4742 $Rd/sim-wspr15-261010_1215.c2")
cases+=("sequence/skim-wav|-|$k 7.0386 $Rd/sim-wav2-261010_1218.wav;$k 14.0956 $T/sim-types-261010_1200.c2;$k 14.0956 $Rd/150426_0918.wav")
# Off the air: five consecutive 40 m periods (recordings/40m-wspr-20261010-1914,
# cut from an UberSDR IQ recording, dial 7.0386): each
# alone, and all five in one data dir, where type 3 messages resolve from
# calls heard in earlier periods. OSD at depth 3 (deep decoding runs 4;
# deeper is too slow to sweep).
R40=recordings/40m-wspr-20261010-1914
k3="-c -C 500 -o 3 -d -f 7.0386"
seq=()
for t in 1914 1916 1918 1920 1922; do
    cases+=("40m-$t/skim-o3|-|$k3 $R40/261010_$t.wav")
    seq+=("$k3 $R40/261010_$t.wav")
done
join() { local IFS=';'; echo "$*"; }
cases+=("sequence/40m-skim-o3|-|$(join "${seq[@]}")")
# Deep decoding's options (no -c, OSD depth 4).
seq=(); for t in 1914 1916 1918 1920 1922; do seq+=("-C 500 -o 4 -d -f 7.0386 $R40/261010_$t.wav"); done
cases+=("sequence/40m-skim-o4|-|$(join "${seq[@]}")")
seq=(); for t in 1914 1916 1918 1920 1922; do seq+=("-f 7.0386 $R40/261010_$t.wav"); done
cases+=("sequence/40m-default|-|$(join "${seq[@]}")")
seq=(); for t in 1914 1916 1918 1920 1922; do seq+=("-J -o 2 -d -w -f 7.0386 $R40/261010_$t.wav"); done
cases+=("sequence/40m-stack-o2|-|$(join "${seq[@]}")")
cases+=("40m-1916/single|-|-s -f 7.0386 $R40/261010_1916.wav")
cases+=("40m-1920/nohash-o2|-|-H -o 2 -q -f 7.0386 $R40/261010_1920.wav")
# .wav files cut short (see inputs above).
cases+=("cut-odd/skim-o3|-|-c -C 500 -o 3 -d -f 7.0386 $T/cut-odd-261010_1218.wav")
cases+=("cut-tiny/osd0|-|-c -o 0 -f 7.0386 $T/cut-tiny-261010_1218.wav")
# Names too short for the date (11 characters before the suffix) or the
# time (4): wsprd reads on into the argument before (argv's strings lie end
# to end, in the order given whatever getopt moves). Each reaches back only
# as far as -a's ".", which run() too is given.
st=$T/sim-types-261010_1200.c2
cases+=("names/short|-|-e 1.5 x.c2;-C 500 a.c2;-f 14.0956 1200.c2;-C 500 x1200.c2|$st x.c2;$st a.c2;$st 1200.c2;$st x1200.c2")
# A file wsprd cannot open for writing: it crashes (fclose or fprintf of
# NULL), losing what stdio had not yet written.
cases+=("readonly/hashtable|-|-c x.c2|$st x.c2;- hashtable.txt ro")
cases+=("readonly/all-wspr|-|-c x.c2|$st x.c2;- ALL_WSPR.TXT ro")
cases+=("readonly/spots|-|-H x.c2|$st x.c2;- wspr_spots.txt ro")
cases+=("readonly/timer|-|-c x.c2|$st x.c2;- wspr_timer.out ro")
# With no decodes the first fprintf of NULL never comes: wsprd crashes at
# fclose(NULL) after the timings, having written wspr_wisdom.dat only.
ct=$T/cut-tiny-261010_1218.wav
cases+=("readonly/all-wspr-none|-|-c -o 0 -f 7.0386 $ct|- ALL_WSPR.TXT ro")
cases+=("readonly/spots-none|-|-o 0 -f 7.0386 $ct|- wspr_spots.txt ro")
# ... and so never rewrites hashtable.txt, odd lines and all.
cases+=("readonly/all-wspr-none-hash|tests/data/wsprd/hashtable-odd.txt|-c -o 0 -f 7.0386 $ct|- ALL_WSPR.TXT ro")
# The stack decoder with a large bias: its metric falls below -6400 and the
# bucket index wraps, a read some 16 GB past its array, a crash for
# certain: in the first pass (-z 20), and after the first pass has written
# the -c file (-z 18). (A large negative bias takes the index past the
# array's end instead, into wsprd's own heap: what that does depends on how
# the heap lies, so no case.)
cases+=("sim-sweep/stack-z20|-|-c -J -z 20 -C 500 -f 14.0956 $T/sim-sweep-261010_1202.c2")
cases+=("sim-sweep/stack-z18|-|-c -J -z 18 -C 500 -f 14.0956 $T/sim-sweep-261010_1202.c2")


run_case() {  # build case-spec rundir [rundir of the seeds] [wisdom]
    local build=$1 spec=$2 dir=$3 seeds=${4:-$work/run/$1} wis=${5:-}
    local name=${spec%%|*} rest=${spec#*|}
    local seed=${rest%%|*} runs=${rest#*|} files=
    if [[ $runs == *"|"* ]]; then files=${runs#*|}; runs=${runs%%|*}; fi
    [[ -d $dir ]] && chmod -R u+w "$dir"
    rm -rf "$dir"; mkdir -p "$dir"
    if [[ $seed == tests/* ]]; then cp "$repo/$seed" "$dir/hashtable.txt"; cp "$dir/hashtable.txt" "$dir/hashtable.in"
    elif [[ $seed != - ]]; then cp "$seeds/$seed/hashtable.txt" "$dir/hashtable.txt"; cp "$dir/hashtable.txt" "$dir/hashtable.in"; fi
    local copied=() f fsrc fname fro
    if [[ -n $files ]]; then
        IFS=';' read -r -a fl <<< "$files"
        : > "$dir/files.txt"
        for f in "${fl[@]}"; do
            echo "$f" >> "$dir/files.txt"
            read -r fsrc fname fro <<< "$f"
            if [[ $fsrc == - ]]; then : > "$dir/$fname"; else cp "$repo/$fsrc" "$dir/$fname"; copied+=("$fname"); fi
            [[ $fro == ro ]] && chmod a-w "$dir/$fname"
        done
    fi
    local i=0 line tok argv
    IFS=';' read -r -a lines <<< "$runs"
    [[ -n $wis ]] && cp "$wis" "$dir/wspr_wisdom.dat"
    : > "$dir/args.txt"; : > "$dir/status.txt"
    for line in "${lines[@]}"; do
        i=$((i+1))
        echo "$line" >> "$dir/args.txt"
        argv=()
        for tok in $line; do
            case $tok in recordings/*|tests/*) argv+=("$repo/$tok") ;; *) argv+=("$tok") ;; esac
        done
        set +e
        (cd "$dir" && "$work/bin/wsprd-$build" -a . "${argv[@]}" > "stdout.$i.txt" 2> "stderr.$i.txt")
        echo $? >> "$dir/status.txt"
        set -e
        [[ -s $dir/stderr.$i.txt ]] || rm -f "$dir/stderr.$i.txt"
        [[ -f $dir/wspr_spots.txt ]] && cp "$dir/wspr_spots.txt" "$dir/wspr_spots.$i.txt"
        if [[ -f $dir/000000_0001.c2 ]]; then
            (cd "$dir" && sha256sum 000000_0001.c2 > "c2.$i.sha256" && rm 000000_0001.c2)
        fi
    done
    rm -f "$dir/wspr_spots.txt" "$dir/wspr_timer.out"
    for f in "${copied[@]}"; do rm -f "$dir/$f"; done
}

jobs_max=$(nproc)
running=0
if [[ -n $wisdom ]]; then
    out=$work/run/wisdom-$(basename "$wisdom")
    for pass in 1 2; do
        for spec in "${cases[@]}"; do
            seed=${spec#*|}; seed=${seed%%|*}
            if [[ ($pass == 1 && $seed == -) || ($pass == 2 && $seed != -) ]]; then
                if (( running >= jobs_max )); then wait -n; running=$((running-1)); fi
                run_case release "$spec" "$out/${spec%%|*}" "$out" "$wisdom" &
                running=$((running+1))
            fi
        done
        wait; running=0
    done
    if diff -r -x wspr_wisdom.dat "$data/cases" "$out" > "$out.diff"; then
        echo "wsprd-oracle: $wisdom changes nothing" >&2
    else
        echo "wsprd-oracle: with $wisdom these differ ($out.diff):" >&2
        diff -rq -x wspr_wisdom.dat "$data/cases" "$out" >&2 || true
    fi
    exit 0
fi
for build in release o2; do
    # seeded cases need their seed's case first
    for pass in 1 2; do
        for spec in "${cases[@]}"; do
            seed=${spec#*|}; seed=${seed%%|*}
            if [[ ($pass == 1 && $seed == -) || ($pass == 2 && $seed != -) ]]; then
                if (( running >= jobs_max )); then wait -n; running=$((running-1)); fi
                run_case "$build" "$spec" "$work/run/$build/${spec%%|*}" &
                running=$((running+1))
            fi
        done
        wait; running=0
    done
done

# The two builds must agree to the byte (wspr_wisdom.dat aside: FFTW's).
if ! diff -r -x wspr_wisdom.dat "$work/run/release" "$work/run/o2" > "$work/release-vs-o2.diff"; then
    echo "wsprd-oracle: the Release and -O2 -ffp-contract=off builds differ: $work/release-vs-o2.diff" >&2
    exit 1
fi
echo "wsprd-oracle: Release and -O2 -ffp-contract=off agree on all ${#cases[@]} cases" >&2

[[ -d $data/cases ]] && chmod -R u+w "$data/cases"
rm -rf "$data/cases"
for spec in "${cases[@]}"; do
    name=${spec%%|*}
    mkdir -p "$data/cases/$name"
    (cd "$work/run/release/$name" && cp -- $(ls | grep -v '^wspr_wisdom.dat$') "$data/cases/$name/")
    chmod -R u+w "$data/cases/$name"
done
echo "wsprd-oracle: ${#cases[@]} cases in $data/cases" >&2
