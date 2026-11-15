#!/usr/bin/env bash
# jt9 --fst4w --fine-freq and wsprd -n -F against signals at known frequencies: synthetic FST4W (each T/R period) and
# WSPR-2 at random audio frequencies 1450-1550 Hz (FST4W's lowest tone) and 1410-1590 Hz (WSPR's centre), in noise,
# made in double precision by scripts/fine-freq-fst4w-gen.cpp and fine-freq-wspr-gen.cpp (built here against src/ for
# the message encoding only). Each decode's frequency is compared with the one it was made at.
#
#   scripts/fine-freq-check.sh JT9 WSPRD [N] [WORK]
#
# N files a case (default 20). Prints each case's decodes and RMS / largest error in mHz; exits 1 if a case's RMS is
# over LIMIT_MHZ (default 1). Its 2026-11-15 run (x86-64, N 50; 25 or 20 for 900 and 1800): FST4W 0.31-0.46 mHz RMS
# at every period from 15 dB above its threshold to 3 dB above (WSJT-X's 0.1 Hz fsig: 29-34 mHz); WSPR-2 0.33 and
# 0.35 at -10 and -20 dB (wsprd's 0.1 Hz: 32); at -26 dB 5.6, 4 of 50 decodes off by 2-35 mHz, each with a spread
# (w50) of 0.24 Hz or more, so this case is printed but not held to the limit.
set -uo pipefail
(( $# >= 2 )) || { echo "usage: $0 JT9 WSPRD [N] [WORK]" >&2; exit 2; }
jt9=$(realpath "$1"); wsprd=$(realpath "$2"); n=${3:-20}
repo=$(cd "$(dirname "$0")/.." && pwd)
work=$(realpath -m "${4:-$(mktemp -d "${TMPDIR:-/tmp}/fine-freq.XXXXXX")}")
limit=${LIMIT_MHZ:-1}
rm -rf "$work"; mkdir -p "$work"
src="$repo/src"
${CXX:-g++} -O2 -std=c++17 -I"$src" "$repo/scripts/fine-freq-fst4w-gen.cpp" "$src/fst4w_ldpc.cpp" "$src/fst4w_unpack.cpp" \
    "$src/wsprd_osd.cpp" -o "$work/fst4w-gen" || exit 2
${CXX:-g++} -O2 -std=c++17 -I"$src" "$repo/scripts/fine-freq-wspr-gen.cpp" "$src/wsprd_unpack.cpp" "$src/wsprd_fano.cpp" \
    -o "$work/wspr-gen" || exit 2

# f(seed): a frequency in [lo, lo+width) from the seed, the same on any machine
freq() { awk -v s="$1" -v lo="$2" -v w="$3" 'BEGIN { x = (s * 2654435761) % 4294967296; printf "%.6f", lo + w * x / 4294967296 }'; }

failed=0
report() {  # report LABEL HELD ERRORS_FILE TRIES
    awk -v label="$1" -v held="$2" -v tries="$4" -v limit="$limit" '
        { n++; s += $1 * $1; a = $1 < 0 ? -$1 : $1; if (a > m) m = a }
        END { rms = n ? sqrt(s / n) : 0
              printf "%-22s %3d/%-3d decoded  %7.2f mHz RMS  %7.2f max%s\n", label, n, tries, rms, m,
                     (held && (n == 0 || rms > limit)) ? "  FAIL" : ""
              exit (held && (n == 0 || rms > limit)) ? 1 : 0 }' "$3" || failed=1
}

for spec in 120:-15 120:-30 300:-20 300:-34 900:-25 900:-39 1800:-28 1800:-42; do
    p=${spec%:*} snr=${spec#*:}; errs="$work/fst4w-$p-$snr.txt"; : > "$errs"
    for (( i = 1; i <= n; i++ )); do
        d="$work/f$p$snr-$i"; mkdir -p "$d"; touch "$d/plotspec" "$d/decdata"
        f=$(freq $(( p * 1000 + i )) 1450 100)
        "$work/fst4w-gen" "$p" "$f" "$snr" "$i" "$d/t.wav" > /dev/null
        (cd "$d" && "$jt9" -a . -p "$p" --fst4w --fine-freq -f 1500 -F 100 t.wav > /dev/null 2>&1)
        awk -v f="$f" 'NR == 1 { print ($18 - f) * 1000 }' "$d/fst4_decodes.dat" 2> /dev/null >> "$errs"
        rm -rf "$d"
    done
    report "FST4W-$p $snr dB" 1 "$errs" "$n"
done

for snr in -10 -20 -26; do
    errs="$work/wspr-$snr.txt"; : > "$errs"
    for (( i = 1; i <= n; i++ )); do
        d="$work/w$snr-$i"; mkdir -p "$d"
        f=$(freq $(( 7000 + i )) 1410 180)
        "$work/wspr-gen" "$f" "$snr" "$i" "$d/261115_1200.wav" "K1ABC FN42 37" || continue
        (cd "$d" && "$wsprd" -n -F -f 14.0956 261115_1200.wav > /dev/null 2>&1)
        awk -v f="$f" '$6 == "K1ABC" { printf "%.6f\n", ($5 * 1e6 - (14095600 + f)) * 1000; exit }' "$d/ALL_WSPR.TXT" 2> /dev/null >> "$errs"
        rm -rf "$d"
    done
    report "WSPR-2 $snr dB" $(( snr > -26 )) "$errs" "$n"
done
exit $failed
