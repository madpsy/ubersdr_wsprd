#!/usr/bin/env bash
# jt9.exe under wine against WSJT-X's own jt9: every case in
# tests/data/fst4w/cases whose inputs are here (recordings/fst4w-sim and
# -real come from scripts/fst4w-oracle.sh --sim-only), run as the case was,
# its stdout, fst4_decodes.dat, decoded.txt and fst4w_calls.txt compared
# byte for byte (line ends aside: Windows' CRLF), fort.52 (the spread's
# debug spectrum) within 1e-5 a value (its last printed digit: FFTW's and
# the maths library's rounding differ on Windows, as on arm64). Then FT8's
# (jt9 -8): the FT8 oracle's jt9 cases (tests/data/ft8/cases/*/jt9-*) and
# the file-mode cases (tests/data/ft8/filemode/cases) whose inputs are here,
# stdout, the exit status, stderr where the case has one and every file
# written byte for byte (line ends aside, and the repository's directory
# in the paths printed). Then FT4's (jt9 -5) the same way: the FT4
# oracle's jt9 cases (tests/data/ft4/cases/*/jt9-*) and its file-mode
# cases (tests/data/ft4/filemode/cases). FT8's and FT4's run several at
# once (JT9_WINE_JOBS, default nproc), each in its own directory.
# JT9_WINE_ONLY=REGEX keeps only the cases whose name (the path under
# tests/data, e.g. ft4/filemode/cases/bad-second/d2 or fst4w/cases/...)
# matches.
#
#   scripts/jt9-wine-check.sh JT9_EXE [WORK]
#
# Prints a line a case that differs, then the counts; exits 1 on any.
set -uo pipefail
(( $# >= 1 )) || { echo "usage: $0 JT9_EXE [WORK]" >&2; exit 2; }
exe=$(realpath "$1")
repo=$(cd "$(dirname "$0")/.." && pwd)
work=$(realpath -m "${2:-$(mktemp -d "${TMPDIR:-/tmp}/jt9-wine.XXXXXX")}")
export WINEPREFIX=$work/prefix WINEDEBUG=-all
mkdir -p "$work"
# Make the prefix first: wine says so on stderr when it creates one, which
# would land in the first case's stderr (FST4W's cases, run first, used to
# hide it; with JT9_WINE_ONLY they may not run).
wineboot --init > /dev/null 2>&1 || true
pass=0 fail=0 skip=0
while IFS= read -r a; do
    d=$(dirname "$a"); name=${d#"$repo/tests/data/fst4w/cases/"}
    [[ -n ${JT9_WINE_ONLY:-} && ! fst4w/cases/$name =~ $JT9_WINE_ONLY ]] && continue
    argv=(-a . --spread exact) missing=
    for t in $(cat "$a"); do
        case $t in
            recordings/*|tests/*) [[ -f $repo/$t ]] || missing=$t; argv+=("$repo/$t") ;;
            *) argv+=("$t") ;;
        esac
    done
    if [[ -n $missing ]]; then skip=$((skip + 1)); continue; fi
    run=$work/run/$name
    rm -rf "$run"; mkdir -p "$run"
    [[ -f $d/fst4w_calls.in ]] && cp "$d/fst4w_calls.in" "$run/fst4w_calls.txt"
    grep -qs bare "$d/flags.txt" || touch "$run/plotspec" "$run/decdata"
    (cd "$run" && wine "$exe" "${argv[@]}" 2> /dev/null | tr -d '\r' > stdout.txt)
    bad=
    cmp -s "$run/stdout.txt" "$d/stdout.txt" || bad+=" stdout"
    for f in fst4_decodes.dat decoded.txt fst4w_calls.txt; do
        if [[ -f $d/$f || -f $run/$f ]]; then
            cmp -s <(tr -d '\r' < "$run/$f" 2> /dev/null) "$d/$f" 2> /dev/null || bad+=" $f"
        fi
    done
    if [[ -f $d/fort.52.gz ]]; then
        python3 -I - "$run/fort.52" "$d/fort.52.gz" <<'PY' || bad+=" fort.52"
import gzip, sys
try:
    a = [l.split() for l in open(sys.argv[1]).read().replace("\r", "").splitlines()]
except OSError:
    sys.exit(1)
b = [l.split() for l in gzip.open(sys.argv[2], "rt").read().splitlines()]
same = len(a) == len(b) and all(
    len(x) == len(y) and all(abs(float(p) - float(q)) <= 1e-5 for p, q in zip(x, y)) for x, y in zip(a, b))
sys.exit(0 if same else 1)
PY
    fi
    if [[ -z $bad ]]; then pass=$((pass + 1)); else fail=$((fail + 1)); echo "differs: $name:$bad"; fi
done < <(find "$repo/tests/data/fst4w/cases" -name args.txt | sort)

# FT8 and FT4: run every case whose inputs are here, then compare.
ft8=() nft4=0
while IFS= read -r a; do
    missing=
    for t in $(cut -d' ' -f2- "$a"); do [[ $t == recordings/* && ! -f $repo/$t ]] && missing=$t; done
    d=$(dirname "$a")
    [[ -n ${JT9_WINE_ONLY:-} && ! ${d#"$repo/tests/data/"} =~ $JT9_WINE_ONLY ]] && continue
    if [[ -n $missing ]]; then skip=$((skip + 1)); else ft8+=("$d"); [[ $d == */ft4/* ]] && nft4=$((nft4 + 1)); fi
done < <(for m in ft8 ft4; do
             find "$repo/tests/data/$m/cases" -path '*/jt9-*' -name args.txt 2> /dev/null | sort
             find "$repo/tests/data/$m/filemode/cases" -name args.txt 2> /dev/null | sort
         done)
run_ft8() {  # case dir
    local d=$1 name run t step steps argv=()
    name=${d#"$repo/tests/data/"}; run=$work/run/$name
    rm -rf "$run"; mkdir -p "$run/cwd"
    if [[ -f $d/prep.txt && $(cat "$d/prep.txt") != - ]]; then
        IFS=';' read -r -a steps < "$d/prep.txt"
        for step in "${steps[@]}"; do
            case $step in
                mkdir:*) mkdir -p "$run/cwd/${step#mkdir:}" ;;
                seed:*) step=${step#seed:}; cp "$repo/${step%%>*}" "$run/cwd/${step#*>}" ;;
            esac
        done
    fi
    for t in $(cut -d' ' -f2- "$d/args.txt"); do
        case $t in recordings/*|tests/*) argv+=("$repo/$t") ;; *) argv+=("$t") ;; esac
    done
    (cd "$run/cwd" && wine "$exe" "${argv[@]}" > ../stdout 2> ../stderr < /dev/null; echo $? > ../status)
    rm -f "$run/cwd/timer.out" "$run/cwd/jt9_wisdom.dat"
}
jobs_max=${JT9_WINE_JOBS:-$(nproc)} running=0
for d in "${ft8[@]}"; do
    if (( running >= jobs_max )); then wait -n; running=$((running - 1)); fi
    run_ft8 "$d" &
    running=$((running + 1))
done
wait
# (the fixtures' stdout and stderr name inputs by their paths in the
# repository, as scripts/ft8-filemode-oracle.sh keeps them)
nocr() { tr -d '\r' < "$1" 2> /dev/null | sed "s|$repo/||g"; }
for d in "${ft8[@]}"; do
    name=${d#"$repo/tests/data/"}; run=$work/run/$name bad=
    exp_out=$d/stdout.txt
    if [[ -f $d/stdout.txt.gz ]]; then gzip -dc "$d/stdout.txt.gz" > "$run/stdout.expected"; exp_out=$run/stdout.expected; fi
    cmp -s <(nocr "$run/stdout") "$exp_out" || bad+=" stdout"
    cmp -s "$run/status" "$d/status.txt" || bad+=" status($(cat "$run/status"))"
    [[ -f $d/stderr.txt ]] && { cmp -s <(nocr "$run/stderr") "$d/stderr.txt" || bad+=" stderr"; }
    # the files written: filemode's under files/, the oracle's beside args.txt
    if [[ -d $d/files || $d == */filemode/* ]]; then
        for f in $( (cd "$d/files" 2> /dev/null && find . -type f; cd "$run/cwd" && find . -type f) | sort -u); do
            cmp -s <(nocr "$run/cwd/$f") "$d/files/$f" 2> /dev/null || bad+=" ${f#./}"
        done
    else
        cmp -s <(nocr "$run/cwd/decoded.txt") "$d/decoded.txt" || bad+=" decoded.txt"
    fi
    if [[ -z $bad ]]; then pass=$((pass + 1)); else fail=$((fail + 1)); echo "differs: $name:$bad"; fi
done

echo "jt9.exe under wine: $pass cases as WSJT-X's jt9 ($((${#ft8[@]} - nft4)) of them FT8, $nft4 FT4), $fail differing, $skip without their inputs"
(( fail == 0 && pass > 0 ))
