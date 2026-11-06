#!/usr/bin/env bash
# jt9.exe under wine (or a native jt9: a path not ending .exe, run as it is:
# the quick --check of build.sh and export-wsprd.sh) against WSJT-X's own
# jt9's stored outputs, no WSJT-X built: every case in
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
# cases (tests/data/ft4/filemode/cases). Then FT2's (jt9 --ft2): the FT2
# oracle's jt9 cases (tests/data/ft2/cases/*/jt9-*, made by the hybrid
# ft2-jt9). Then MSK144's (jt9 -k): the cases in tests/data/msk144/cases/
# (the reference's outputs, scripts/msk144-oracle.sh), left out while the jt9
# run does not list --msk144 in its usage, the shmsg and swl variants
# (--msk-shmsg, --msk-swl) while it lists neither (a jt9 listing --spread
# runs them all, and fails when an option is missing), and the cases on an odd
# input (@I/..., made by scripts/jt9-compare.sh) always, counted and
# printed. FT8's, FT4's, FT2's and
# MSK144's run several at once (JT9_WINE_JOBS,
# default nproc), each in its own directory.
# A native jt9 off x86-64 (aarch64: the Pi) may differ where
# tests/data/jt9-arm64-accepted.txt says it does: FFTW's NEON codelets round
# differently from the AVX ones the stored outputs were made with (user
# 2026-10-13, docs/ft8-decoder.md phase 7; tests/data/ft8/arm64-fftw.txt).
# Each line there is a case, a file and the sha256 of the arm64 output
# accepted for it: those are counted "accepted", not differing; any other
# difference, or another output for a listed one, still fails.
# JT9_CHECK_PRINT_ACCEPT=1 prints such a line for each difference (to add
# to the list once a difference is understood and accepted).
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
native=false; [[ $exe == *.exe ]] || native=true
# (the fixtures' stdout and stderr name inputs by their paths in the
# repository, as scripts/ft8-filemode-oracle.sh keeps them)
nocr() { tr -d '\r' < "$1" 2> /dev/null | sed "s|$repo/||g"; }
accepted_list=$repo/tests/data/jt9-arm64-accepted.txt
may_accept=false; $native && [[ $(uname -m) != x86_64 && -f $accepted_list ]] && may_accept=true
# accepted CASE FILE GOT: GOT's output (as compared) listed for CASE's FILE.
accepted() {
    local sha; sha=$(nocr "$3" | sha256sum | cut -d' ' -f1)
    [[ ${JT9_CHECK_PRINT_ACCEPT:-} == 1 ]] && echo "accept-line: $1 $2 $sha" >&2
    $may_accept && grep -q "^$1 $2 $sha\b" "$accepted_list"
}
# The jt9 run: under wine, or as it is.
jt9() { if $native; then "$exe" "$@"; else wine "$exe" "$@"; fi; }
# Make the prefix first: wine says so on stderr when it creates one, which
# would land in the first case's stderr (FST4W's cases, run first, used to
# hide it; with JT9_WINE_ONLY they may not run).
$native || wineboot --init > /dev/null 2>&1 || true
pass=0 fail=0 skip=0 accept=0
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
    (cd "$run" && jt9 "${argv[@]}" 2> /dev/null | tr -d '\r' > stdout.txt)
    bad= acc= cname=fst4w/cases/$name
    differs() {  # FILE GOT
        if accepted "$cname" "$1" "$2"; then acc+=" $1"; else bad+=" $1"; fi
    }
    cmp -s "$run/stdout.txt" "$d/stdout.txt" || differs stdout "$run/stdout.txt"
    for f in fst4_decodes.dat decoded.txt fst4w_calls.txt; do
        if [[ -f $d/$f || -f $run/$f ]]; then
            cmp -s <(tr -d '\r' < "$run/$f" 2> /dev/null) "$d/$f" 2> /dev/null || differs "$f" "$run/$f"
        fi
    done
    if [[ -f $d/fort.52.gz ]]; then
        python3 -I - "$run/fort.52" "$d/fort.52.gz" <<'PY'

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
        [[ $? == 0 ]] || differs fort.52 "$run/fort.52"
    fi
    if [[ -n $bad ]]; then fail=$((fail + 1)); echo "differs: $cname:$bad"
    elif [[ -n $acc ]]; then accept=$((accept + 1)); echo "accepted (arm64, FFTW's NEON rounding): $cname:$acc"
    else pass=$((pass + 1)); fi
done < <(find "$repo/tests/data/fst4w/cases" -name args.txt | sort)

# FT8, FT4 and FT2: run every case whose inputs are here, then compare.
ft8=() nft4=0 nft2=0 nmsk=0 ninput=0 ninput_msk=0
usage_text=$(jt9 -h 2> /dev/null | tr -d '\r' || true)
# Our port lists --spread (the WSJT-X reference does not). With it, every MSK144
# case runs, the shmsg and swl variants (--msk-shmsg, --msk-swl) too, and fails
# when an option is missing: a binary that lists --spread must run all of
# them. Without --spread (the reference as the jt9 under test) the usage
# decides.
msk_ok=false msk_flags=false msk_left=0
[[ $usage_text == *--msk144* ]] && msk_ok=true
[[ $usage_text == *--msk-shmsg* && $usage_text == *--msk-swl* ]] && msk_flags=true
if [[ $usage_text == *--spread* ]]; then msk_ok=true msk_flags=true; fi
while IFS= read -r a; do
    missing=
    d=$(dirname "$a")
    if [[ $d == */msk144/* ]]; then
        if ! $msk_ok || { [[ $d == */jt9-shmsg || $d == */jt9-swl ]] && ! $msk_flags; }; then msk_left=$((msk_left + 1)); continue; fi
    fi
    for t in $(cut -d' ' -f2- "$a"); do
        [[ $t == recordings/* && ! -f $repo/$t ]] && missing=$t
        # @I/ inputs are made by scripts/jt9-compare.sh (odd WAVs), not kept
        # in the repository: such cases cannot run here, in any mode. They are
        # counted (ninput, and ninput_msk for MSK144's) and printed below
        if [[ $t == @I/* ]]; then
            ninput=$((ninput + 1)); [[ $d == */msk144/* ]] && ninput_msk=$((ninput_msk + 1))
            continue 2
        fi
    done
    [[ -n ${JT9_WINE_ONLY:-} && ! ${d#"$repo/tests/data/"} =~ $JT9_WINE_ONLY ]] && continue
    if [[ -n $missing ]]; then skip=$((skip + 1)); else ft8+=("$d"); [[ $d == */ft4/* ]] && nft4=$((nft4 + 1)); [[ $d == */ft2/* ]] && nft2=$((nft2 + 1)); [[ $d == */msk144/* ]] && nmsk=$((nmsk + 1)); fi
done < <(for m in ft8 ft4 ft2 msk144; do
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
    # --feed's stdin (MSK144's cases): stdin.txt, its paths from the repository's root
    local in=/dev/null
    if [[ -f $d/stdin.txt ]]; then sed "s|^|$repo/|" "$d/stdin.txt" > "$run/stdin"; in=$run/stdin; fi
    (cd "$run/cwd" && jt9 "${argv[@]}" > ../stdout 2> ../stderr < "$in"; echo $? > ../status)
    find "$run/cwd" \( -name timer.out -o -name jt9_wisdom.dat \) -delete
}
jobs_max=${JT9_WINE_JOBS:-$(nproc)} running=0
for d in "${ft8[@]}"; do
    if (( running >= jobs_max )); then wait -n; running=$((running - 1)); fi
    run_ft8 "$d" &
    running=$((running + 1))
done
wait
for d in "${ft8[@]}"; do
    name=${d#"$repo/tests/data/"}; run=$work/run/$name bad= acc=
    differs() {  # FILE GOT
        if accepted "$name" "$1" "$2"; then acc+=" $1"; else bad+=" $1"; fi
    }
    exp_out=$d/stdout.txt
    if [[ -f $d/stdout.txt.gz ]]; then gzip -dc "$d/stdout.txt.gz" > "$run/stdout.expected"; exp_out=$run/stdout.expected; fi
    cmp -s <(nocr "$run/stdout") "$exp_out" || differs stdout "$run/stdout"
    cmp -s "$run/status" "$d/status.txt" || bad+=" status($(cat "$run/status"))"
    [[ -f $d/stderr.txt ]] && { cmp -s <(nocr "$run/stderr") "$d/stderr.txt" || bad+=" stderr"; }
    # the files written: filemode's under files/, the oracle's beside args.txt
    if [[ -d $d/files || $d == */filemode/* ]]; then
        for f in $( (cd "$d/files" 2> /dev/null && find . -type f; cd "$run/cwd" && find . -type f) | sort -u); do
            cmp -s <(nocr "$run/cwd/$f") "$d/files/$f" 2> /dev/null || differs "${f#./}" "$run/cwd/$f"
        done
    elif [[ $d == */msk144/* ]]; then  # nothing is written but what is stored (usually nothing)
        for f in $( (cd "$d" && ls | grep -v -x -e args.txt -e stdout.txt -e status.txt -e stderr.txt -e stdin.txt; cd "$run/cwd" && find . -type f | sed 's|^\./||') | sort -u); do
            cmp -s <(nocr "$run/cwd/$f") "$d/$f" 2> /dev/null || differs "$f" "$run/cwd/$f"
        done
    else
        cmp -s <(nocr "$run/cwd/decoded.txt") "$d/decoded.txt" || differs decoded.txt "$run/cwd/decoded.txt"
    fi
    if [[ -n $bad ]]; then fail=$((fail + 1)); echo "differs: $name:$bad"
    elif [[ -n $acc ]]; then accept=$((accept + 1)); echo "accepted (arm64, FFTW's NEON rounding): $name:$acc"
    else pass=$((pass + 1)); fi
done

echo "jt9-wine-check: MSK144: $nmsk cases ran, $msk_left left out for lack of an option (--msk144 in the usage: $msk_ok, the flags: $msk_flags)" >&2
echo "jt9-wine-check: $ninput cases left out for their @I/ inputs, $ninput_msk of them MSK144's (scripts/jt9-compare.sh makes those WAVs; --deep-check runs them)" >&2
if [[ $usage_text == *--spread* ]] && [[ -z ${JT9_WINE_ONLY:-} ]]; then
    # The floor: 347 MSK144 cases that run here, plus 20 for the real UACNJ slots (recordings/msk144-uacnj). The 347 is
    # the 367 runs of the corpus as it stood at Task 8, less the 20 @I/ runs (odd-stereo, -rate11025, -bits8,
    # -header and -tiny, 4 variants each). The final fix wave added odd-trunc and odd-stream (8 more @I/ runs,
    # left out too: ninput_msk counts them) and the 20 real-uacnj runs. Raise it with the corpus. Nothing may
    # be left out for lack of an option, and the @I/ runs left out may not number fewer than the 28 there are.
    (( nmsk < 367 )) && { echo "jt9-wine-check: FAIL: only $nmsk of the 367 MSK144 cases ran" >&2; fail=$((fail + 1)); }
    (( msk_left > 0 )) && { echo "jt9-wine-check: FAIL: $msk_left MSK144 cases left out" >&2; fail=$((fail + 1)); }
    (( ninput_msk != 28 )) && { echo "jt9-wine-check: FAIL: $ninput_msk MSK144 cases left out for their @I/ inputs, 28 expected" >&2; fail=$((fail + 1)); }
fi
label="jt9.exe under wine"; $native && label="$(basename "$exe")"
acc_note=; (( accept )) && acc_note=", $accept as accepted for arm64 (tests/data/jt9-arm64-accepted.txt)"
echo "$label: $pass cases as WSJT-X's jt9 ($((${#ft8[@]} - nft4 - nft2 - nmsk)) of them FT8, $nft4 FT4, $nft2 FT2, $nmsk MSK144)$acc_note, $fail differing, $skip without their inputs"
(( fail == 0 && pass > 0 ))
