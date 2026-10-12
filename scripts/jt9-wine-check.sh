#!/usr/bin/env bash
# jt9.exe under wine against WSJT-X's own jt9: every case in
# tests/data/fst4w/cases whose inputs are here (recordings/fst4w-sim and
# -real come from scripts/fst4w-oracle.sh --sim-only), run as the case was,
# its stdout, fst4_decodes.dat, decoded.txt and fst4w_calls.txt compared
# byte for byte (line ends aside: Windows' CRLF), fort.52 (the spread's
# debug spectrum) within 1e-5 a value (its last printed digit: FFTW's and
# the maths library's rounding differ on Windows, as on arm64).
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
pass=0 fail=0 skip=0
while IFS= read -r a; do
    d=$(dirname "$a"); name=${d#"$repo/tests/data/fst4w/cases/"}
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
echo "jt9.exe under wine: $pass cases as WSJT-X's jt9, $fail differing, $skip without their inputs"
(( fail == 0 && pass > 0 ))
