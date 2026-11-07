#!/usr/bin/env bash
# The reference for MSK144 in our jt9: WSJT-X's own jt9, msk144sim and
# libwsjt_fort.a at 567ad29ce (3.2.0-rc1), built by WSJT-X's CMake in the
# wsjtx-oracle image (docker/wsjtx-oracle.Dockerfile: Ubuntu 24.04, gfortran
# 13.3, FFTW 3.3.10, as scripts/fst4w-oracle.sh uses), and the corpus
# tests/data/msk144 that every piece of the port is checked against.
#
#   scripts/msk144-oracle.sh build         the three builds, the two-build check
#   scripts/msk144-oracle.sh corpus        inputs, the cases run, fixtures stored
#   scripts/msk144-oracle.sh focus         Task 6's cases (repeat-across-files, lastmode-*,
#                                          ntol-wide, freq-edge, off-grid) and Task 7's swl-state and swl-pair, one
#                                          reference process a file
#   scripts/msk144-oracle.sh extra         the final review's cases: odd-trunc, odd-stream (a WAV cut
#                                          mid-data; the streaming size 0xFFFFFFFF) and real-uacnj-* (20
#                                          real slots, meteor-scatter's command line), one reference
#                                          process a file
#   scripts/msk144-oracle.sh all           build, corpus, focus and extra
#   scripts/msk144-oracle.sh --print-env   the paths a later step sources
#   scripts/msk144-oracle.sh fortran FILE.f90 [NAME [omp]]
#                                          compile a Fortran oracle against the
#                                          reference's libwsjt_fort.a; print the
#                                          binary's path
#   scripts/msk144-oracle.sh exec CMD ...  run CMD in the image, the work dir
#                                          and the repository mounted at their
#                                          own paths (the binaries link the
#                                          image's libraries, not the host's)
#
#   WSJTX=~/repos/wsjtx            the WSJT-X tree, at REF_COMMIT
#   MSK144_ORACLE_WORK=DIR         where to build and run (default
#                                  ~/.cache/msk144-oracle; kept, so rebuilds
#                                  are incremental)
#   DOCKER=podman                  the container runner (default: docker if
#                                  there is one, else podman). Under podman
#                                  the runs take --cgroup-manager=cgroupfs and
#                                  --userns=keep-id: a rootless podman with no
#                                  user session cannot use systemd's.
#
# Three builds, all WSJT-X's Release unless said:
#   release   WSJT-X's own flags: the reference ($MSK144_REF)
#   o2        -O2 -ffp-contract=off; must agree with release on the corpus
#   patched   release with decode_msk144.f90's bshmsg and bswl (the MSK40
#             short messages and the SWL mode, which jt9 never turns on) read
#             from the environment: MSK144_SHMSG=1, MSK144_SWL=1. Both set is
#             the two forced true; each alone is what our --msk-shmsg and
#             --msk-swl are to do. With neither set it is release's decoder.
#
# The corpus's variants (tests/data/msk144/README.md): d1, d3, defaults and
# ms (the line meteor-scatter runs) from release; shmsg and swl from patched.
set -euo pipefail

REF_COMMIT=567ad29ce
IMAGE=wsjtx-oracle
repo=$(cd "$(dirname "$0")/.." && pwd)

# ================================================================ inside
# The container runs this script again with --inside: /src is the WSJT-X
# tree (read-only), /repo this repository, /work the work dir.
if [[ ${1:-} == --inside ]]; then
shift
cmd=$1; shift
work=/work; src=/src; repo=/repo
b=$work/release

# libquadmath: gfortran's on x86-64, none on aarch64.
quadmath=; [[ $(gcc -print-file-name=libquadmath.so) == /* ]] && quadmath=-lquadmath

# fcompile SRC.f90 OUT: a Fortran main against the release build's library,
# with the flags CMake gives msk144sim, as WSJT-X links its own simulators.
fcompile() {
    local f=$1 out=$2 variant=${3:-} o mismatch=()
    local lib=libwsjt_fort.a mods=fortran_modules omp=()
    # variant omp: against the OpenMP variant of the library, as the release jt9 links it
    if [[ $variant == omp ]]; then lib=libwsjt_fort_omp.a mods=fortran_modules_omp omp=(-fopenmp); fi
    o=$(mktemp -d "$work/fc.XXXXXX")
    read -r -a fflags <<< "$(sed -n 's/^Fortran_FLAGS = //p' "$b/CMakeFiles/msk144sim.dir/flags.make")"
    # -fallow-argument-mismatch for the msk40 oracle only (the others compile strictly): lib/msk40spd.f90
    # leaves nhasharray undeclared (a scalar by implicit typing) and passes the array mskrtd gave it on to msk40decodeframe, which declares it 10 by 10: the reference's files
    # are separate and compile; an oracle that holds copies of both in one file is told of the mismatch
    [[ $(basename "$f") == msk40-* ]] && mismatch=(-fallow-argument-mismatch)
    (cd "$b" && gfortran "${fflags[@]}" "${omp[@]}" "${mismatch[@]}" -Wno-intrinsic-shadow -I $mods -I "$src/lib" -J "$o" \
            -c "$f" -o "$o/main.o" &&
        c++ -fopenmp -pthread "$o/main.o" -o "$out" $lib libwsjt_cxx.a -lm \
            /usr/lib/$(gcc -dumpmachine)/libboost_{log_setup,log,chrono,filesystem,regex,thread,atomic}.so.1.83.0 \
            -lfftw3f_threads -lfftw3f -lgfortran $quadmath)
    rm -rf "$o"
}

# ---- helpers shared by the focus and extra steps
data=$repo/tests/data/msk144
odd=$work/inputs

# write_sums: SHA256SUMS, every input by its path from the repository's root
write_sums() {
    (cd "$repo" && find recordings/msk144-samples recordings/msk144-real recordings/msk144-uacnj recordings/msk144-sim recordings/msk144-float \
        -name '*.wav' | sort | xargs sha256sum > tests/data/msk144/inputs/SHA256SUMS)
}

# make_odd: the odd WAVs (stereo, 11025 Hz, truncated, ...) are the ones scripts/jt9-compare.sh makes: its
# generator, taken out of the script, so the two cannot differ
make_odd() {
    rm -rf "$odd"; mkdir -p "$odd"
    sed -n '/^JT9_COMPARE_REPO=\$repo python3 -I - /,/^EOF$/p' "$repo/scripts/jt9-compare.sh" | sed '1d;$d' |
        JT9_COMPARE_REPO=/repo python3 -I - "$odd" > /dev/null
}

# one_run DIR ARGS...: the reference, once, in DIR; its stdout with the
# repository's directory taken out of the paths. --msk-shmsg and --msk-swl are our options, not the
# reference's: the patched build runs, with MSK144_SHMSG=1 and MSK144_SWL=1 (as the corpus step does)
one_run() {
    local dir=$1; shift
    rm -rf "$dir"; mkdir -p "$dir"
    local argv=() tok envs=() build=release
    for tok in "$@"; do
        case $tok in
            recordings/*) argv+=("$repo/$tok") ;;
            @I/*) argv+=("$odd/${tok#@I/}") ;;
            --msk-shmsg) envs+=(MSK144_SHMSG=1); build=patched ;;
            --msk-swl) envs+=(MSK144_SWL=1); build=patched ;;
            *) argv+=("$tok") ;;
        esac
    done
    set +e
    (cd "$dir" && env "${envs[@]}" "$work/$build/jt9" "${argv[@]}" > stdout.txt 2> stderr.txt < /dev/null)
    echo $? > "$dir/status.txt"
    set -e
    sed -i "s|$repo/||g; s|$odd/|@I/|g" "$dir/stdout.txt" "$dir/stderr.txt"
    rm -f "$dir/timer.out" "$dir/jt9_wisdom.dat" "$dir/plotspec" "$dir/decdata"
}
# store DIR ARGS-LINE [STDIN-FILE]: the case directory from the run just made in $dir
store() {
    local out=$1 run=$2 args=$3
    rm -rf "$out"; mkdir -p "$out"
    echo "jt9 $args" > "$out/args.txt"
    cp "$run/stdout.txt" "$run/status.txt" "$out/"
    [[ -s $run/stderr.txt ]] && cp "$run/stderr.txt" "$out/"
    [[ -f $run/decoded.txt ]] && cp "$run/decoded.txt" "$out/"
    return 0
}
# single CASE VARIANT FILE ARGS...: one file, one process
single() {
    local case=$1 v=$2 file=$3; shift 3
    one_run "$work/focus-run/$case-$v" "$@" "$file"
    store "$data/cases/$case/jt9-$v" "$work/focus-run/$case-$v" "$* $file"
}

case $cmd in
# ---------------------------------------------------------------- build
build|all)
    # WSJT-X's CMake, Release; the docs, man pages, tests and QMAP left out
    # (they do not touch jt9). The o2 build overrides only the optimisation
    # flags CMake's Release type sets: CMake/CompilerFlags.cmake adds its own.
    configure() {  # srcdir builddir [cmake args]
        local s=$1 dir=$2; shift 2
        mkdir -p "$dir"
        cmake -S "$s" -B "$dir" -DCMAKE_BUILD_TYPE=Release -DWSJT_GENERATE_DOCS=OFF \
            -DWSJT_SKIP_MANPAGES=ON -DWSJT_ENABLE_TESTS=OFF -DWSJT_SKIP_QMAP=ON "$@" > "$dir/configure.log" 2>&1
    }
    o2="-O2 -ffp-contract=off"
    configure "$src" "$work/release"
    configure "$src" "$work/o2" -DCMAKE_Fortran_FLAGS_RELEASE="$o2" -DCMAKE_C_FLAGS_RELEASE="$o2" \
        -DCMAKE_CXX_FLAGS_RELEASE="$o2"
    # The patched copy: the tree without .git and samples, and the one file
    # decode_msk144.f90 changed by sed. The two switches it declares false
    # (bshmsg, bswl) are read from MSK144_SHMSG and MSK144_SWL instead.
    rm -rf "$work/src-patched"; mkdir -p "$work/src-patched"
    (cd "$src" && tar --exclude=.git --exclude=./samples -cf - .) | tar -C "$work/src-patched" -xf -
    f=$work/src-patched/lib/decode_msk144.f90
    sed -i \
        -e 's/^  logical :: bswl = \.false\.$/&\n  character(len=1) :: msk144_env/' \
        -e 's/^  message_count = 0$/&\n  call get_environment_variable("MSK144_SHMSG", msk144_env)\n  bshmsg = msk144_env == "1"\n  call get_environment_variable("MSK144_SWL", msk144_env)\n  bswl = msk144_env == "1"/' \
        "$f"
    [[ $(diff <(cat "$src/lib/decode_msk144.f90") "$f" | grep -c '^>') == 5 ]] ||
        { echo "msk144-oracle: the sed patch did not apply to decode_msk144.f90" >&2; exit 1; }
    diff "$src/lib/decode_msk144.f90" "$f" > "$work/patch.diff" || true
    configure "$work/src-patched" "$work/patched"
    cmake --build "$work/release" --target jt9 msk144sim -j"$(nproc)" > "$work/release/build.log" 2>&1
    cmake --build "$work/o2" --target jt9 msk144sim -j"$(nproc)" > "$work/o2/build.log" 2>&1
    cmake --build "$work/patched" --target jt9 -j"$(nproc)" > "$work/patched/build.log" 2>&1
    # msk144-sim (scripts/msk144-sim.f90) against the release library, and
    # a Fortran oracle's directory
    mkdir -p "$work/oracles"
    fcompile "$repo/scripts/msk144-sim.f90" "$b/msk144-sim"
    echo "msk144-oracle: built in $work/release, $work/o2 and $work/patched" >&2

    # The two-build check on WSJT-X's own samples, and the patched build
    # (neither switch set) against release on them.
    chk=$work/check; rm -rf "$chk"; mkdir -p "$chk"
    bad=0
    for w in 181211_120500 181211_120800; do
        for build in release o2 patched; do
            mkdir -p "$chk/$build"
            (cd "$chk/$build" && "$work/$build/jt9" -a . -k -p 15 -f 1500 -F 50 -d 3 "$src/samples/MSK144/$w.wav" \
                > "$w.out" 2> "$w.err"; echo $? > "$w.status")
        done
        for build in o2 patched; do
            for e in out err status; do
                cmp -s "$chk/release/$w.$e" "$chk/$build/$w.$e" || { echo "msk144-oracle: $build differs from release on $w ($e)" >&2; bad=1; }
            done
        done
        echo "msk144-oracle: $w: $(grep -c '&' "$chk/release/$w.out") decode(s): $(grep '&' "$chk/release/$w.out" | tr -s ' ' | tr '\n' ';')" >&2
    done
    (( bad == 0 )) || exit 1
    echo "msk144-oracle: release, $o2 and patched agree on the two samples" >&2
    [[ $cmd == all ]] || exit 0
    ;&

# ---------------------------------------------------------------- corpus
corpus)
    data=$repo/tests/data/msk144
    samples=$repo/recordings/msk144-samples
    mkdir -p "$samples" "$repo/recordings/msk144-real" "$data/inputs"
    cp "$src/samples/MSK144/181211_120500.wav" "$src/samples/MSK144/181211_120800.wav" "$samples/"
    for f in 160915_113100 160915_113230; do
        [[ -f $samples/$f.wav ]] || echo "msk144-oracle: $samples/$f.wav is not here (README.md says where from); its cases are left out" >&2
    done
    [[ -f $repo/recordings/msk144-real/msk144_261023_105030.wav ]] ||
        echo "msk144-oracle: recordings/msk144-real/msk144_261023_105030.wav is not here; its cases are left out" >&2
    bash "$repo/scripts/msk144-sim-grid.sh" "$b/msk144-sim" "$b/msk144sim" "$repo/recordings/msk144-sim" "$samples"

    make_odd
    write_sums

    # ---- the cases: name|T/R|input|extra
    # input is a path from the repository's root, or @I/ an odd WAV
    cases=()
    while read -r stem tr; do
        [[ -n $stem ]] || continue
        cases+=("$(sed 's/-261106_[0-9]*$//' <<< "$stem")|$tr|recordings/msk144-sim/$stem.wav")
    done < "$repo/recordings/msk144-sim/MANIFEST"
    for f in "$samples"/*.wav; do
        [[ -f $f ]] || continue
        cases+=("sample-$(basename "$f" .wav)|15|recordings/msk144-samples/$(basename "$f")")
    done
    if [[ -f $repo/recordings/msk144-real/msk144_261023_105030.wav ]]; then
        cases+=("real-261023-p30|30|recordings/msk144-real/msk144_261023_105030.wav"
                "real-261023-p15|15|recordings/msk144-real/msk144_261023_105030.wav")
    fi
    cases+=("odd-stereo|15|@I/stereo-261012_1200.wav" "odd-rate11025|15|@I/rate11025-261012_1200.wav"
            "odd-bits8|15|@I/bits8-261012_1200.wav" "odd-header|15|@I/header-261012_1200.wav"
            "odd-tiny|15|@I/tiny-261012_1200.wav")
    # the plain variants for every case; shmsg and swl for the cases below
    flagged='^(sample-|real-|snr-|msg-|pings|msk40-|stock-|odd-silence|odd-noise|odd-clipped)'

    # variant templates: FILE, TR; -a . first as the other modes' args.txt have it
    argv_for() {  # variant tr file
        case $1 in
            d1)       echo "-a . -k -p $2 -f 1500 -F 200 -d 1 $3" ;;
            d3)       echo "-a . -k -p $2 -f 1500 -F 200 -d 3 $3" ;;
            defaults) echo "-a . -k -p $2 $3" ;;
            ms)       echo "-Y --msk144 -p $2 -f 1500 -F 200 -d 3 -a . $3" ;;   # meteor-scatter's line
            shmsg)    echo "-a . -k -p $2 -f 1500 -F 200 -d 3 --msk-shmsg $3" ;;
            swl)      echo "-a . -k -p $2 -f 1500 -F 200 -d 3 --msk-swl $3" ;;
        esac
    }
    # the work list: case|variant|tr|file
    jobs=()
    for spec in "${cases[@]}"; do
        IFS='|' read -r name tr file <<< "$spec"
        vars=(d1 d3 defaults ms)
        [[ $name =~ $flagged ]] && vars+=(shmsg swl)
        for v in "${vars[@]}"; do jobs+=("$name|$v|$tr|$file"); done
    done
    # jt9 with no options but -k (period 60 s, as jt9 runs it): the samples
    for spec in "${cases[@]}"; do
        IFS='|' read -r name tr file <<< "$spec"
        [[ $name == sample-* ]] && jobs+=("$name|bare|$tr|$file")
    done

    run_job() {  # build job rundir
        local build=$1 job=$2 dir=$3 name v tr file
        IFS='|' read -r name v tr file <<< "$job"
        local args envs=()
        if [[ $v == bare ]]; then args="-a . -k $file"; else args=$(argv_for "$v" "$tr" "$file"); fi
        rm -rf "$dir"; mkdir -p "$dir"
        echo "jt9 $args" > "$dir/args.txt"
        local argv=() tok
        for tok in $args; do
            case $tok in
                recordings/*) argv+=("$repo/$tok") ;;
                @I/*) argv+=("$odd/${tok#@I/}") ;;
                --msk-shmsg) envs+=(MSK144_SHMSG=1) ;;
                --msk-swl) envs+=(MSK144_SWL=1) ;;
                *) argv+=("$tok") ;;
            esac
        done
        set +e
        (cd "$dir" && env "${envs[@]}" "$work/$build/jt9" "${argv[@]}" > stdout.txt 2> stderr.txt < /dev/null)
        echo $? > "$dir/status.txt"
        set -e
        sed -i "s|$repo/||g; s|$odd/|@I/|g" "$dir/stdout.txt" "$dir/stderr.txt"
        [[ -s $dir/stderr.txt ]] || rm -f "$dir/stderr.txt"
        rm -f "$dir/timer.out" "$dir/jt9_wisdom.dat" "$dir/plotspec" "$dir/decdata"
    }
    run_all() {  # build out-dir [jobs matching REGEX only]
        local build=$1 out=$2 only=${3:-} running=0 job name v
        for job in "${jobs[@]}"; do
            IFS='|' read -r name v _ <<< "$job"
            [[ -n $only && ! $v =~ $only ]] && continue
            if (( running >= $(nproc) )); then wait -n; running=$((running - 1)); fi
            run_job "$build" "$job" "$out/$name/jt9-$v" &
            running=$((running + 1))
        done
        wait
    }
    rm -rf "$work/run"
    run_all release "$work/run/release" '^(d1|d3|defaults|ms|bare)$'
    run_all o2 "$work/run/o2" '^(d1|d3|defaults|ms|bare)$'
    run_all patched "$work/run/patched" '^(shmsg|swl)$'
    # release and o2 agree on every run but the flagged variants (patched only)
    ok=1; nrun=0; : > "$work/release-vs-o2.diff"
    for d in "$work"/run/o2/*/jt9-*; do
        n=${d#"$work/run/o2/"}; nrun=$((nrun + 1))
        diff -r "$work/run/release/$n" "$d" >> "$work/release-vs-o2.diff" 2>&1 || { echo "differs: $n" >> "$work/release-vs-o2.diff"; ok=0; }
    done
    if (( ok )); then echo "msk144-oracle: release and -O2 agree on all $nrun runs" >&2
    else echo "msk144-oracle: release and -O2 differ: $work/release-vs-o2.diff" >&2; exit 1; fi
    # The patched build's flagged cases come from it; the rest from release.
    for d in "$work"/run/patched/*/jt9-*; do
        n=${d#"$work/run/patched/"}
        rm -rf "$work/run/release/$n"; mkdir -p "$(dirname "$work/run/release/$n")"; cp -a "$d" "$work/run/release/$n"
    done

    # the data dirs keep args, stdout, status, stderr where there is some, and
    # the empty decoded.txt jt9 opens: MSK144 writes nothing to it
    extra=$(cd "$work/run/release" && find . -type f ! -name args.txt ! -name stdout.txt ! -name status.txt \
            ! -name stderr.txt ! \( -name decoded.txt -size 0 \) || true)
    [[ -z $extra ]] || echo "msk144-oracle: files written beside stdout: $(head -5 <<< "$extra")" >&2

    rm -rf "$data/cases"
    mkdir -p "$data/cases"
    cp -a "$work/run/release/." "$data/cases/"
    {
        echo "WSJT-X $(git -c safe.directory='*' -C "$src" rev-parse --short=9 HEAD 2>/dev/null || echo "$REF_COMMIT")"
        gfortran --version | head -1
        echo "FFTW $(dpkg-query -W -f='${Version}' libfftw3-single3) ($(uname -m))"
        echo "Release Fortran: $(sed -n 's/^Fortran_FLAGS = //p' "$work/release/CMakeFiles/wsjt_fort_omp.dir/flags.make")"
        echo "O2 Fortran:      $(sed -n 's/^Fortran_FLAGS = //p' "$work/o2/CMakeFiles/wsjt_fort_omp.dir/flags.make")"
        echo "patch (patched build, decode_msk144.f90):"; sed 's/^/  /' "$work/patch.diff"
    } > "$data/BUILD.txt"
    # what the corpus holds
    echo "msk144-oracle: ${#cases[@]} cases, ${#jobs[@]} runs, in $data/cases" >&2
    for v in d1 d3 defaults ms shmsg swl bare; do
        printf 'msk144-oracle: %-8s %3d runs, %3d with a decode\n' "$v" \
            "$(ls -d "$data"/cases/*/jt9-$v 2>/dev/null | wc -l)" \
            "$(grep -l '^[0-9]' "$data"/cases/*/jt9-$v/stdout.txt 2>/dev/null | wc -l)" >&2
    done
       [[ $cmd == all ]] || exit 0
    ;&


# ---------------------------------------------------------------- focus
# Task 6's cases, stored in the corpus's layout beside the others (the corpus
# step is not run again): the command line's behaviours that each look wrong
# and must stay as they stand. One reference process a file, with OMP_DYNAMIC
# and OMP_THREAD_LIMIT unset: after its first file the reference's decoder.f90
# calls omp_set_dynamic(.true.), which makes msk144sync's team (and so a
# second file's result) depend on the load. A multi-file case's outputs are the
# single-file runs' put end to end: the files' stdout in turn; for --feed, with
# "<FileDone> PATH N" after each (N the decode lines of that file).
focus)
    data=$repo/tests/data/msk144
    unset OMP_DYNAMIC OMP_THREAD_LIMIT
    sim=recordings/msk144-sim
    sample=recordings/msk144-samples/181211_120500.wav
    hashab=$(ls "$repo/$sim"/hash-ab-*.wav | sed "s|^$repo/||")
    hashb=$(ls "$repo/$sim"/hash-b-*.wav | sed "s|^$repo/||")
    m200=$(ls "$repo/$sim"/freq-m200-sm02-*.wav | sed "s|^$repo/||")
    edge=$(ls "$repo/$sim"/freq-edge-f0800-*.wav | sed "s|^$repo/||")
    og6=$(ls "$repo/$sim"/offgrid-f1537-sm06-*.wav | sed "s|^$repo/||")
    og4=$(ls "$repo/$sim"/offgrid-f1537-sm04-*.wav | sed "s|^$repo/||")
    swla=$(ls "$repo/$sim"/swl-state-a-*.wav | sed "s|^$repo/||")
    swlb=$(ls "$repo/$sim"/swl-state-b-*.wav | sed "s|^$repo/||")
    swlp=$(ls "$repo/$sim"/swl-pair-*.wav | sed "s|^$repo/||")
    early=$(ls "$repo/$sim"/early-ping-*.wav | sed "s|^$repo/||")
    B="-a . -k -p 15 -f 1500 -F 200 -d 3"

    # multi CASE VARIANT FILES...: one process a file, the outputs put end to end, as `jt9 B FILES` and
    # as `jt9 --feed B` print them
    multi() {
        local case=$1 v=$2 B=${MULTI_B:-$B}; shift 2
        local files=("$@") i=0 f d n
        local cat=$work/focus-run/$case-$v.files cfeed=$work/focus-run/$case-$v.feed
        rm -rf "$work/focus-run/$case-$v" "$cat" "$cfeed"; mkdir -p "$work/focus-run/$case-$v"
        : > "$work/focus-run/$case-$v/stdout.txt"; : > "$cfeed.stdout"
        : > "$work/focus-run/$case-$v/stdin.txt"
        for f in "${files[@]}"; do
            i=$((i + 1)); d=$work/focus-run/$case-$v/one$i
            one_run "$d" $B "$f"
            [[ $(cat "$d/status.txt") == 0 && ! -s $d/stderr.txt ]] || { echo "msk144-oracle: focus: $f: status/stderr" >&2; exit 1; }
            cat "$d/stdout.txt" >> "$work/focus-run/$case-$v/stdout.txt"
            cat "$d/stdout.txt" >> "$cfeed.stdout"
            n=$(grep -c -v -e '^ EOF' -e '^<DecodeFinished>' "$d/stdout.txt" || true)
            echo "<FileDone> $f $n" >> "$cfeed.stdout"
            echo "$f" >> "$work/focus-run/$case-$v/stdin.txt"
            [[ -f $d/decoded.txt ]] && : > "$work/focus-run/$case-$v/decoded.txt"
        done
        echo 0 > "$work/focus-run/$case-$v/status.txt"
        # the file list on the command line, and --feed's
        store "$data/cases/$case/jt9-files" "$work/focus-run/$case-$v" "$B ${files[*]}"
        mkdir -p "$data/cases/$case/jt9-feed"
        echo "jt9 --feed $B" > "$data/cases/$case/jt9-feed/args.txt"
        cp "$cfeed.stdout" "$data/cases/$case/jt9-feed/stdout.txt"
        cp "$work/focus-run/$case-$v/status.txt" "$work/focus-run/$case-$v/stdin.txt" "$data/cases/$case/jt9-feed/"
        [[ -f $work/focus-run/$case-$v/decoded.txt ]] && cp "$work/focus-run/$case-$v/decoded.txt" "$data/cases/$case/jt9-feed/"
        return 0
    }
    rm -rf "$work/focus-run"; mkdir -p "$work/focus-run"
    rm -rf "$data"/cases/{repeat-across-files,lastmode-Y,lastmode-8,ntol-wide,freq-edge,offgrid-sm06,offgrid-sm04,swl-state,swl-pair}

    # a decode in the first block twice (its SNR, with no noise level yet, shows state kept from the file before),
    # the same decode further in twice, then a hashed call that the first file's earlier message resolves and a
    # file that holds only the hashed message (<...> there), twice over
    multi repeat-across-files all "$early" "$early" "$sample" "$sample" "$hashab" "$hashb" "$hashab" "$hashb"
    # the last mode wins, in both orders
    single lastmode-Y Ymsk144 "$sample" -Y --msk144 -p 15 -f 1500 -F 200 -d 3 -a .
    single lastmode-Y msk144Y "$sample" --msk144 -Y -p 15 -f 1500 -F 200 -d 3 -a .
    single lastmode-8 msk1448 "$sample" --msk144 -8 -p 15 -f 1500 -F 200 -d 3 -a .
    single lastmode-8 8msk144 "$sample" -8 --msk144 -p 15 -f 1500 -F 200 -d 3 -a .
    # -F beyond what jt9.f90 allows (min(ntol, 1000))
    for F in 300 1000 1001 5000; do single ntol-wide "F$F" "$m200" -a . -k -p 15 -f 1500 -F $F -d 3; done
    # -f at the band's edges
    # (f2000-F5000: 1200 Hz from the carrier, past what jt9.f90's limit of 1000 reaches; f1801 and f1802 sit
    # 1001 and 1002 Hz away, where a limit of 1001 would answer differently)
    for a in "200 1000" "200 700" "200 200" "0 1000" "800 1000" "6000 5000" "9000 1000" "2000 5000" "2000 1000" \
             "1801 5000" "1802 5000" "1803 5000"; do
        set -- $a; single freq-edge "f$1-F$2" "$edge" -a . -k -p 15 -f $1 -F $2 -d 3
    done
    # SWL mode across files (Task 7): the MSK40 pair hash of the first decode must not outlive its file. The
    # first and last files hold the MSK40 pings alone; the middle one the standard message that gives the
    # pair. One process a file, the reference's table is always fresh; a port that kept it would print the
    # last file's line another way
    MULTI_B="$B --msk-swl" multi swl-state all "$swlb" "$swla" "$swlb"
    # and within one file: the standard message at 5 s, then the MSK40 pings at 7 s and 11 s that its calls explain
    single swl-pair swl "$swlp" $B --msk-swl
    # an off-grid carrier: -d 1 (short pings only) and -d 3 (the averages too)
    for d in 1 3; do
        single offgrid-sm06 d$d "$og6" -a . -k -p 15 -f 1500 -F 200 -d $d
        single offgrid-sm04 d$d "$og4" -a . -k -p 15 -f 1500 -F 200 -d $d
    done
    echo "msk144-oracle: focus: $(ls -d "$data"/cases/{repeat-across-files,lastmode-Y,lastmode-8,ntol-wide,freq-edge,offgrid-sm06,offgrid-sm04,swl-state,swl-pair}/jt9-* | wc -l) runs stored" >&2
    [[ $cmd == all ]] || exit 0
    ;&

# ---------------------------------------------------------------- extra
# The final review's cases, stored beside the others, one reference process a file with OMP_DYNAMIC and
# OMP_THREAD_LIMIT unset (as focus does; every case here is one file):
#   odd-trunc, odd-stream  a WAV cut mid-data, and one with the streaming size 0xFFFFFFFF (@I/, made by
#                          scripts/jt9-compare.sh's generator): the reference refuses both, 3.0.2 decoded
#                          the samples present (tests/data/msk144/AB-3.0.2.md); variants d1, d3, defaults, ms
#   real-uacnj-*           20 slots from UACNJ's meteor-scatter archive (recordings/msk144-uacnj, 12 kHz
#                          s16 mono, 30 s, each with the JSON sidecar meteor-scatter wrote): variant ms,
#                          the line meteor-scatter runs
extra)
    unset OMP_DYNAMIC OMP_THREAD_LIMIT
    make_odd
    write_sums
    rm -rf "$data"/cases/{odd-trunc,odd-stream} "$data"/cases/real-uacnj-*
    for f in trunc stream; do
        # the corpus step's variants d1, d3, defaults and ms, at T/R 15 s
        single "odd-$f" d1 "@I/$f-261012_1200.wav" -a . -k -p 15 -f 1500 -F 200 -d 1
        single "odd-$f" d3 "@I/$f-261012_1200.wav" -a . -k -p 15 -f 1500 -F 200 -d 3
        single "odd-$f" defaults "@I/$f-261012_1200.wav" -a . -k -p 15
        single "odd-$f" ms "@I/$f-261012_1200.wav" -Y --msk144 -p 15 -f 1500 -F 200 -d 3 -a .
    done
    for w in "$repo"/recordings/msk144-uacnj/*.wav; do
        n=$(basename "$w" .wav)
        single "real-uacnj-$n" ms "recordings/msk144-uacnj/$n.wav" -Y --msk144 -p 30 -f 1500 -F 200 -d 3 -a .
    done
    echo "msk144-oracle: extra: $(ls -d "$data"/cases/{odd-trunc,odd-stream,real-uacnj-*}/jt9-* | wc -l) runs stored" >&2
    ;;

# ---------------------------------------------------------------- fortran
fortran)
    f=$1; name=${2:-$(basename "$f" .f90)}
    [[ $f == /* ]] || f=$PWD/$f
    mkdir -p "$work/oracles"
    fcompile "$f" "$work/oracles/$name" "${3:-}"
    echo "$work/oracles/$name"
    ;;
exec)
    exec "$@"
    ;;
*) echo "msk144-oracle: unknown step $cmd" >&2; exit 2 ;;
esac
exit 0
fi

# ================================================================ outside
wsjtx=${WSJTX:-$HOME/repos/wsjtx}
work=${MSK144_ORACLE_WORK:-${XDG_CACHE_HOME:-$HOME/.cache}/msk144-oracle}
mkdir -p "$work"
work=$(cd "$work" && pwd)

print_env() {
    cat <<EOF
MSK144_ORACLE_WORK=$work
MSK144_REF=$work/release/jt9
MSK144_REF_O2=$work/o2/jt9
MSK144_REF_PATCHED=$work/patched/jt9
MSK144_SIM=$work/release/msk144sim
MSK144_SIM_SEEDED=$work/release/msk144-sim
MSK144_LIB=$work/release/libwsjt_fort.a
MSK144_MODULES=$work/release/fortran_modules
MSK144_IMAGE=$IMAGE
EOF
}
if [[ ${1:-} == --print-env ]]; then print_env; exit 0; fi
step=${1:-}
case $step in build|corpus|focus|extra|all|fortran|exec) shift ;;
    *) sed -n '2,25p' "$0" | sed 's/^# \?//' >&2; exit 2 ;;
esac

# the container runner: docker if there is one, else podman
runner=${DOCKER:-}
if [[ -z $runner ]]; then
    if command -v docker > /dev/null; then runner=docker; else runner=podman; fi
fi
ctr_pre=() ctr_build=() ctr_run=()
if [[ $(basename "$runner") == podman ]]; then
    ctr_pre=(--cgroup-manager=cgroupfs); ctr_build=(--isolation=chroot); ctr_run=(--userns=keep-id)
fi

head=$(git -c safe.directory='*' -C "$wsjtx" rev-parse --short=9 HEAD)
if [[ $head != "$REF_COMMIT" ]]; then
    echo "msk144-oracle: $wsjtx is at $head, not $REF_COMMIT" >&2
    exit 1
fi
if [[ -n $(git -c safe.directory='*' -C "$wsjtx" status --porcelain -- lib CMake CMakeLists.txt) ]]; then
    echo "msk144-oracle: $wsjtx has local changes in lib/ or its CMake files" >&2
    exit 1
fi
"$runner" "${ctr_pre[@]}" build "${ctr_build[@]}" -q -t "$IMAGE" -f "$repo/docker/wsjtx-oracle.Dockerfile" "$repo/docker" > /dev/null
echo "msk144-oracle: working in $work" >&2
if [[ $step == exec ]]; then
    exec "$runner" "${ctr_pre[@]}" run --rm "${ctr_run[@]}" -u "$(id -u):$(id -g)" -e HOME=/work \
        -v "$wsjtx:/src:ro" -v "$repo:$repo" -v "$work:$work" "$IMAGE" "$@"
fi
if [[ $step == fortran ]]; then
    f=$(realpath "$1"); shift
    out=$("$runner" "${ctr_pre[@]}" run --rm "${ctr_run[@]}" -u "$(id -u):$(id -g)" -e HOME=/work \
        -v "$wsjtx:/src:ro" -v "$repo:/repo" -v "$work:/work" -v "$f:/in/$(basename "$f"):ro" "$IMAGE" \
        bash /repo/scripts/msk144-oracle.sh --inside fortran "/in/$(basename "$f")" "$@")
    echo "${out/#\/work/$work}"
    exit 0
fi
exec "$runner" "${ctr_pre[@]}" run --rm "${ctr_run[@]}" -u "$(id -u):$(id -g)" -e HOME=/work \
    -v "$wsjtx:/src:ro" -v "$repo:/repo" -v "$work:/work" "$IMAGE" \
    bash /repo/scripts/msk144-oracle.sh --inside "$step"
