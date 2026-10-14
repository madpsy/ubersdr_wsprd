#!/usr/bin/env bash
# A wsprd against WSJT-X's own wsprd's stored outputs, no WSJT-X built (the
# quick --check of build.sh and export-wsprd.sh; scripts/wsprd-compare.sh,
# the deep one, runs WSJT-X's wsprd beside it): every case in
# tests/data/wsprd/cases whose inputs are here, run as
# scripts/wsprd-oracle.sh ran it (a fresh data dir that is also the cwd,
# "-a ." first, seeded with hashtable.in and files.txt's files, each args.txt
# line a run), as tests/test_wsprd.cpp runs them in process: each run's
# stdout, wspr_spots.txt and 000000_0001.c2 (its sha256: a .wav's only on
# x86-64, where FFTW's first FFTs are the reference's), then the statuses,
# ALL_WSPR.TXT and hashtable.txt, byte for byte (line ends aside: Windows'
# CRLF). A wsprd ending .exe runs under wine.
#
#   scripts/wsprd-stored-check.sh WSPRD [WORK]
#
# Prints a line a case that differs, then the counts; exits 1 on any.
# WSPRD_STORED_JOBS=n cases at once (default nproc).
set -uo pipefail
(( $# >= 1 )) || { echo "usage: $0 WSPRD [WORK]" >&2; exit 2; }
exe=$(realpath "$1")
repo=$(cd "$(dirname "$0")/.." && pwd)
work=$(realpath -m "${2:-$(mktemp -d "${TMPDIR:-/tmp}/wsprd-stored.XXXXXX")}")
mkdir -p "$work"
export WINEPREFIX=$work/prefix WINEDEBUG=-all
native=false; [[ $exe == *.exe ]] || native=true
wsprd() { if $native; then "$exe" "$@"; else wine "$exe" "$@"; fi; }
$native || wineboot --init > /dev/null 2>&1 || true
wav_c2=false; [[ $(uname -m) == x86_64 ]] && $native && wav_c2=true
nocr() { tr -d '\r' < "$1"; }
cases=$repo/tests/data/wsprd/cases
results=$work/results; rm -rf "$results"; mkdir -p "$results"
# One case: its verdict in results/ ("pass", "skip" or what differs).
one() {
    local cdir=$1 name t missing dir rostill from to mode bad n status line argv out err c2 f
    name=${cdir#"$cases"/}
    missing=
    for t in $(cat "$cdir/args.txt"); do
        [[ $t == recordings/* || $t == tests/* ]] && [[ ! -f $repo/$t ]] && missing=$t
    done
    local verdict=$results/${name//\//__}
    if [[ -n $missing ]]; then echo skip > "$verdict"; return; fi
    dir=$work/run/$name
    chmod -R u+w "$dir" 2> /dev/null; rm -rf "$dir"; mkdir -p "$dir"
    [[ -f $cdir/hashtable.in ]] && cp "$cdir/hashtable.in" "$dir/hashtable.txt"
    rostill=
    if [[ -f $cdir/files.txt ]]; then
        while read -r from to mode; do
            [[ -z ${to:-} ]] && continue
            if [[ $from == - ]]; then : > "$dir/$to"; else cp "$repo/$from" "$dir/$to"; fi
            if [[ ${mode:-} == ro ]]; then
                chmod a-w "$dir/$to"
                [[ -w $dir/$to ]] && rostill=1  # root: read-only changes nothing
            fi
        done < "$cdir/files.txt"
    fi
    if [[ -n $rostill ]]; then echo skip > "$verdict"; return; fi
    bad= n=0 status=
    while IFS= read -r line; do
        n=$((n + 1))
        argv=(-a .)
        for t in $line; do
            case $t in recordings/*|tests/*) argv+=("$repo/$t") ;; *) argv+=("$t") ;; esac
        done
        (cd "$dir" && wsprd "${argv[@]}" > "../$(basename "$dir").out.$n" 2> "../$(basename "$dir").err.$n" < /dev/null)
        status+="$?"$'\n'
        out=$work/run/$(dirname "$name")/$(basename "$dir").out.$n
        err=$work/run/$(dirname "$name")/$(basename "$dir").err.$n
        cmp -s <(nocr "$out") "$cdir/stdout.$n.txt" || bad+=" stdout.$n"
        [[ -s $err ]] && bad+=" stderr.$n"
        if [[ -f $cdir/wspr_spots.$n.txt || -f $dir/wspr_spots.txt ]]; then
            cmp -s <(nocr "$dir/wspr_spots.txt" 2> /dev/null) "$cdir/wspr_spots.$n.txt" 2> /dev/null || bad+=" wspr_spots.$n"
        fi
        c2=$dir/000000_0001.c2
        if [[ -f $cdir/c2.$n.sha256 || -f $c2 ]]; then
            if [[ ! -f $cdir/c2.$n.sha256 || ! -f $c2 ]]; then
                bad+=" c2.$n(there in one only)"
            elif $wav_c2 || [[ $line != *.wav* ]]; then
                [[ "$(cd "$dir" && sha256sum 000000_0001.c2)" == "$(cat "$cdir/c2.$n.sha256")" ]] || bad+=" c2.$n"
            fi
            rm -f "$c2"
        fi
    done < "$cdir/args.txt"
    [[ "$status" == "$(cat "$cdir/status.txt")"$'\n' ]] || bad+=" status"
    for f in ALL_WSPR.TXT hashtable.txt; do
        if [[ -f $cdir/$f || -f $dir/$f ]]; then
            cmp -s <(nocr "$dir/$f" 2> /dev/null) "$cdir/$f" 2> /dev/null || bad+=" $f"
        fi
    done
    chmod -R u+w "$dir" 2> /dev/null
    if [[ -z $bad ]]; then echo pass > "$verdict"; else echo "differs: $name:$bad" > "$verdict"; fi
}
jobs_max=${WSPRD_STORED_JOBS:-$(nproc)} running=0
for cdir in "$cases"/*/*/; do
    if (( running >= jobs_max )); then wait -n; running=$((running - 1)); fi
    one "${cdir%/}" &
    running=$((running + 1))
done
wait
pass=0 fail=0 skip=0
for v in "$results"/*; do
    case $(cat "$v") in
        pass) pass=$((pass + 1)) ;;
        skip) skip=$((skip + 1)) ;;
        *) fail=$((fail + 1)); cat "$v" ;;
    esac
done
label="$(basename "$exe")"; $native || label="$label under wine"
echo "$label: $pass cases as WSJT-X's wsprd's stored outputs, $fail differing, $skip without their inputs"
(( fail == 0 && pass > 0 ))
