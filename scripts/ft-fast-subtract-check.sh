#!/usr/bin/env bash
# jt9 -8, -5 and --ft2 --fast-subtract, and with -w 2 (measured FFT plans), against jt9 as WSJT-X decodes: neither is
# WSJT-X's to the bit, so they are held to its decodes instead. Each FT8, FT4 and FT2 recording set (recordings/ft8-real,
# -sim, -filemode, the same for FT4, ft2-sim and ft2-live's 40 m set) is decoded as one run (one decoder for its files,
# in name order) at depths 1, 2 and 3, as jt9 decodes it and the two other ways, and each way's decodes (period, message; '?' ones and the a1..a7 marks left out) are compared with
# jt9's: every one must keep at least MIN_SAME percent (default 99) of jt9's decodes and find no more than the rest
# new. On x86-64 every run so far has had exactly jt9's decodes.
#
#   scripts/ft-fast-subtract-check.sh JT9 [WORK]
#
# Prints a line a run, then the counts; exits 1 if any run falls short.
set -uo pipefail
(( $# >= 1 )) || { echo "usage: $0 JT9 [WORK]" >&2; exit 2; }
jt9=$(realpath "$1")
repo=$(cd "$(dirname "$0")/.." && pwd)
work=$(realpath -m "${2:-$(mktemp -d "${TMPDIR:-/tmp}/ft-fast-subtract.XXXXXX")}")
min_same=${MIN_SAME:-99}
rm -rf "$work"; mkdir -p "$work"

# The decodes of a run's stdout (the mode's separator after the frequency: FT8 ~, FT4 +, FT2 *), one
# 'PERIOD# MESSAGE' a line
decodes() {
    awk -v sep="$2" '/^<DecodeFinished>/ { n++; next }
         index($0, " " sep " ") { m = substr($0, index($0, " " sep " ") + 3); sub(/^ +/, "", m); sub(/ +$/, "", m)
                 if (m ~ /[?]$/) next
                 sub(/ +a[0-9]$/, "", m); gsub(/  +/, " ", m); print n, m }' "$1" | sort -u
}

pass=0 fail=0
check_set() {
    local mode=$1 set=$2 sep=$3 f depth way run other nb no same new lost pct line
    local -a files=() opts=()
    for f in "$repo/recordings/$set"/*.wav; do [[ -f $f ]] && files+=("${f##*/}"); done
    (( ${#files[@]} )) || return 0
    for depth in 1 2 3; do
        for way in exact fast fast-w2; do
            run=$work/${set//\//_}-d$depth-$way; mkdir -p "$run"
            cp "$repo/recordings/$set"/*.wav "$run/"
            opts=()
            [[ $way == fast ]] && opts=(--fast-subtract)
            [[ $way == fast-w2 ]] && opts=(-w 2 --fast-subtract)
            (cd "$run" && "$jt9" "$mode" -a . -d "$depth" "${opts[@]}" "${files[@]}" > stdout 2> stderr)
            rm -f "$run"/*.wav
            decodes "$run/stdout" "$sep" > "$run/decodes"
        done
        base=$work/${set//\//_}-d$depth-exact/decodes
        nb=$(wc -l < "$base")
        for way in fast fast-w2; do
            other=$work/${set//\//_}-d$depth-$way/decodes
            no=$(wc -l < "$other")
            same=$(comm -12 "$base" "$other" | wc -l)
            new=$((no - same)) lost=$((nb - same))
            pct=$(awk -v s="$same" -v b="$nb" 'BEGIN { printf "%.1f", b ? 100 * s / b : 100 }')
            line="jt9 $mode recordings/$set -d $depth $([[ $way == fast-w2 ]] && echo '-w 2 ')--fast-subtract: $no decodes against $nb: $same the same ($pct%), $lost lost, $new new"
            if awk -v s="$same" -v b="$nb" -v n="$new" -v m="$min_same" \
                    'BEGIN { exit !(b == 0 || (100 * s >= m * b && 100 * n <= (100 - m) * b)) }'; then
                echo "$line"; pass=$((pass + 1))
            else
                echo "SHORT: $line"; fail=$((fail + 1))
            fi
        done
    done
}

for set in ft8-real ft8-sim ft8-filemode; do check_set -8 "$set" '~'; done
for set in ft4-real ft4-sim ft4-filemode; do check_set -5 "$set" '+'; done
for set in ft2-sim ft2-live/40m-20261014-0226; do check_set --ft2 "$set" '*'; done
echo "FT8/FT4/FT2 --fast-subtract (and -w 2): $((pass + fail)) runs held to ${min_same}% of jt9's decodes, $fail short"
(( fail == 0 ))
