#!/usr/bin/env bash
# jt9 -8 --fast-subtract, and with -w 2 (measured FFT plans), against jt9 -8 as WSJT-X decodes: neither is WSJT-X's to
# the bit, so they are held to its decodes instead. Each FT8 recording set (recordings/ft8-real, -sim, -filemode) is
# decoded as one run (one decoder for its files, in name order) at depths 1, 2 and 3, as jt9 decodes it and the two
# other ways, and each way's decodes (period, message; '?' ones and the a1..a7 marks left out) are compared with
# jt9's: every one must keep at least MIN_SAME percent (default 99) of jt9's decodes and find no more than the rest
# new. On x86-64 every run so far has had exactly jt9's decodes.
#
#   scripts/ft8-fast-subtract-check.sh JT9 [WORK]
#
# Prints a line a run, then the counts; exits 1 if any run falls short.
set -uo pipefail
(( $# >= 1 )) || { echo "usage: $0 JT9 [WORK]" >&2; exit 2; }
jt9=$(realpath "$1")
repo=$(cd "$(dirname "$0")/.." && pwd)
work=$(realpath -m "${2:-$(mktemp -d "${TMPDIR:-/tmp}/ft8-fast-subtract.XXXXXX")}")
min_same=${MIN_SAME:-99}
rm -rf "$work"; mkdir -p "$work"

# The decodes of a run's stdout, one 'PERIOD# MESSAGE' a line
decodes() {
    awk '/^<DecodeFinished>/ { n++; next }
         / ~ / { m = $0; sub(/^.* ~  */, "", m); sub(/ +$/, "", m); if (m ~ /[?]$/) next
                 sub(/ +a[0-9]$/, "", m); gsub(/  +/, " ", m); print n, m }' "$1" | sort -u
}

pass=0 fail=0
for set in ft8-real ft8-sim ft8-filemode; do
    files=()
    for f in "$repo/recordings/$set"/*.wav; do [[ -f $f ]] && files+=("${f##*/}"); done
    (( ${#files[@]} )) || continue
    for depth in 1 2 3; do
        for way in exact fast fast-w2; do
            run=$work/$set-d$depth-$way; mkdir -p "$run"
            cp "$repo/recordings/$set"/*.wav "$run/"
            opts=()
            [[ $way == fast ]] && opts=(--fast-subtract)
            [[ $way == fast-w2 ]] && opts=(-w 2 --fast-subtract)
            (cd "$run" && "$jt9" -8 -a . -d "$depth" "${opts[@]}" "${files[@]}" > stdout 2> stderr)
            rm -f "$run"/*.wav
            decodes "$run/stdout" > "$run/decodes"
        done
        base=$work/$set-d$depth-exact/decodes
        nb=$(wc -l < "$base")
        for way in fast fast-w2; do
            other=$work/$set-d$depth-$way/decodes
            no=$(wc -l < "$other")
            same=$(comm -12 "$base" "$other" | wc -l)
            new=$((no - same)) lost=$((nb - same))
            pct=$(awk -v s="$same" -v b="$nb" 'BEGIN { printf "%.1f", b ? 100 * s / b : 100 }')
            line="recordings/$set -d $depth $([[ $way == fast-w2 ]] && echo '-w 2 ')--fast-subtract: $no decodes against $nb: $same the same ($pct%), $lost lost, $new new"
            if awk -v s="$same" -v b="$nb" -v n="$new" -v m="$min_same" \
                    'BEGIN { exit !(b == 0 || (100 * s >= m * b && 100 * n <= (100 - m) * b)) }'; then
                echo "$line"; pass=$((pass + 1))
            else
                echo "SHORT: $line"; fail=$((fail + 1))
            fi
        done
    done
done
echo "FT8 --fast-subtract (and -w 2): $((pass + fail)) runs held to ${min_same}% of jt9's decodes, $fail short"
(( fail == 0 ))
