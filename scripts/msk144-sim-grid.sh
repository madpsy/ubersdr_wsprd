#!/usr/bin/env bash
# The simulated MSK144 inputs of tests/data/msk144 (made by
# scripts/msk144-oracle.sh, which builds the two generators and runs this in
# the oracle's container).
#
#   scripts/msk144-sim-grid.sh SIM STOCK_SIM OUTDIR [SAMPLES_DIR]
#
#   SIM        scripts/msk144-sim.f90's binary (seeded noise, chosen pings)
#   STOCK_SIM  WSJT-X's own msk144sim (a ping every second; unseeded noise,
#              which is rand() from its default seed: the same every run)
#   OUTDIR     where the WAVs and MANIFEST go (recordings/msk144-sim)
#   SAMPLES_DIR  with 181211_120800.wav in it: the source of the clipped case
#
# Every file is named CASE-261106_HHMMSS.wav: jt9 takes nutc from the last six
# digits. MANIFEST lists "stem T/R-seconds", what the oracle runs each with.
# Seeds are fixed: file number n has seed 1000+n, so a case keeps its noise
# when cases are added at the end.
#
# The grid, each axis swept with the others at their base (the message
# "K1ABC W9XYZ EN37", 1500 Hz, one ping 0.2 s wide at 5 s, T/R 15 s):
#   snr      -8 to +10 dB in 2 dB steps
#   freq     offsets 0, +-100, +-200 Hz from 1500 Hz, at +4 and -2 dB
#   width    0.05, 0.2, 1, 5 s and the whole slot, at +4 and 0 dB
#   pings    1 and 3 a slot, at 0 and +6 dB
#   tr       T/R 30 s (two pings), at -2, +4 and +10 dB; and a 2.5 s width
#   msg      six more message types at +6 dB
#   msk40    four MSK40 frames (genmsk40, through the short-message
#            syntax "<CALL CALL> report") at +6 and 0 dB
#   stock    msk144sim's own two usage examples
#   odd      silence, noise alone, a 30 s file (run with -p 15), a 7 s file,
#            and the 181211_120800 sample clipped
set -euo pipefail
(( $# >= 3 )) || { sed -n '2,12p' "$0" >&2; exit 2; }
sim=$(realpath "$1") stock=$(realpath "$2")
out=$3; samples=${4:-}
rm -rf "$out"; mkdir -p "$out"
out=$(realpath "$out")
n=0
manifest=$out/MANIFEST; : > "$manifest"

# gen CASE TR SNR WIDTH PINGS FREQ MESSAGE [more msk144-sim options]
gen() {
    local case=$1 tr=$2 snr=$3 width=$4 pings=$5 freq=$6 msg=$7; shift 7
    n=$((n + 1))
    local stem; stem=$(printf '%s-261106_12%02d%02d' "$case" $((n / 60)) $((n % 60)))
    "$sim" -p "$tr" -S $((1000 + n)) -f "$freq" -s "$snr" -w "$width" -t "$pings" "$@" \
        -o "$out/$stem.wav" "$msg" > /dev/null
    echo "$stem $tr" >> "$manifest"
}
sn() { printf '%s%02d' "$([[ $1 == -* ]] && echo m || echo p)" "${1#-}"; }
M="K1ABC W9XYZ EN37"

for s in -8 -6 -4 -2 0 2 4 6 8 10; do gen "snr-$(sn $s)" 15 $s 0.2 5 1500 "$M"; done
for o in -200 -100 0 100 200; do
    for s in 4 -2; do
        gen "freq-$(sn $o)-s$(sn $s)" 15 $s 0.2 5 $((1500 + o)) "$M"
    done
done
for w in 0.05 0.2 1 5 slot; do
    for s in 4 0; do gen "width-${w}-s$(sn $s)" 15 $s $w 5 1500 "$M"; done
done
for s in 0 6; do
    gen "pings1-s$(sn $s)" 15 $s 0.2 8 1500 "$M"
    gen "pings3-s$(sn $s)" 15 $s 0.2 2,7,11 1500 "$M"
done
for s in -2 4 10; do gen "tr30-s$(sn $s)" 30 $s 0.2 4,19 1500 "$M"; done
gen "tr30-wide-s04" 30 4 2.5 4,19 1500 "$M"
gen "msg-cq" 15 6 0.2 5 1500 "CQ K1ABC FN42"
gen "msg-report" 15 6 0.2 5 1500 "K1ABC W9XYZ R-03"
gen "msg-rrr" 15 6 0.2 5 1500 "K1ABC W9XYZ RRR"
gen "msg-prefix" 15 6 0.2 5 1500 "PJ4/K1ABC W9XYZ"
gen "msg-suffix" 15 6 0.2 5 1500 "K1ABC/P W9XYZ EN37"
gen "msg-freetext" 15 6 0.2 5 1500 "TNX 73 GL"
# MSK40: the hash of the calls in <> is the 12-bit key jt9's defaults (-c K1ABC
# -x W9XYZ) decode against.
gen "msk40-rpt03" 15 6 0.2 5 1500 "<K1ABC W9XYZ> +03"
gen "msk40-r03" 15 6 0.2 5 1500 "<W9XYZ K1ABC> R-03"
gen "msk40-rrr" 15 0 0.2 5 1500 "<K1ABC W9XYZ> RRR"
gen "msk40-73" 15 0 0.2 5 1500 "<K1ABC W9XYZ> 73"

# msk144sim's own two usage examples, in a scratch directory (it writes
# 000000_000001.wav to the current one).
tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
for spec in "15 1500 0.12 2 stock-tr15-s02" "30 1500 2.5 15 stock-tr30-s15"; do
    read -r tr f w s name <<< "$spec"
    n=$((n + 1))
    stem=$(printf '%s-261106_12%02d%02d' "$name" $((n / 60)) $((n % 60)))
    (cd "$tmp" && "$stock" "K1ABC W9XYZ EN37" "$tr" "$f" "$w" "$s" 1 > /dev/null)
    mv "$tmp/000000_000001.wav" "$out/$stem.wav"
    echo "$stem $tr" >> "$manifest"
done

# odd inputs
gen "odd-silence" 15 0 0.2 5 1500 "" -Z
gen "odd-noise" 15 0 0.2 5 1500 "" -n
gen "odd-tr30-run-p15" 30 4 0.2 4,19 1500 "$M"
sed -i '$ s/ 30$/ 15/' "$manifest"        # a 30 s file, run with -p 15
gen "odd-short7" 7 4 0.2 3 1500 "$M"
sed -i '$ s/ 7$/ 15/' "$manifest"         # a 7 s file, run with -p 15
if [[ -n $samples && -f $samples/181211_120800.wav ]]; then
    n=$((n + 1))
    stem=$(printf 'odd-clipped-261106_12%02d%02d' $((n / 60)) $((n % 60)))
    python3 -I - "$samples/181211_120800.wav" "$out/$stem.wav" <<'EOF'
import struct, sys
b = open(sys.argv[1], "rb").read()
i = b.index(b"data")
n = struct.unpack("<I", b[i + 4:i + 8])[0] // 2
s = struct.unpack(f"<{n}h", b[i + 8:i + 8 + 2 * n])
c = [max(-32768, min(32767, 40 * v)) for v in s]       # 40x gain, clipped
open(sys.argv[2], "wb").write(b[:i + 8] + struct.pack(f"<{n}h", *c))
EOF
    echo "$stem 15" >> "$manifest"
fi

# Task 6's inputs (the Review Focus cases of scripts/msk144-oracle.sh focus),
# listed in MANIFEST.focus, not MANIFEST: the corpus runs no variants on them.
#   hash-a, hash-b, hash-ab: a ping of "K1ABC W9XYZ EN37" at 3 s; a ping of
#       "<W9XYZ> PJ4/K1ABC RRR" at 9 s (a hashed call: it prints <...> unless the
#       decoder has seen W9XYZ in full); and the two files' samples added (the
#       hash resolves in the second ping)
#   offgrid-f1537-sm06, offgrid-f1537-sm04: the carrier at 1537 Hz (on none of the
#       averaging searches' frequency grids, as every other carrier is), a 1 s
#       ping at -6 and -4 dB: the -d 1 short-ping decoder misses them, and
#       the -d 3 averages find them, so msk144sync's step (10/navg Hz) shows
#   early-ping: a ping at 0.15 s, +10 dB: decoded in the first block, where mskrtd
#       has no noise level yet (SNR 0), so a state carried over from a file
#       before it would show in the SNR
#   freq-edge-f0800: the carrier at 800 Hz, +8 dB: -f and -F reach it from
#       the band's lower edge (-f 200 -F 700, -f 0, a window cut at 0 Hz)
#   swl-state-a, swl-state-b (Task 7): a ping of "N0ABC W1AW EN34" at +10 dB, and three
#       pings of the MSK40 frame "<N0ABC W1AW> RRR" at +8 dB (jt9's calls are K1ABC
#       and W9XYZ, so the hash is none the decoder is asked for). Run with
#       --msk-swl, files b, a, b: the pair's hash is in the SWL table when the
#       second b runs only if the table survives from a to the b after it
#   swl-pair: the two samples added, as hash-ab is: the standard message at 5 s
#       gives the decoder both calls, and the MSK40 pings at 7 s and 11 s are then
#       "<N0ABC W1AW> RRR" by the pair's hash (the ping at 2 s came first: "<hhhh>")
focus=$out/MANIFEST.focus; : > "$focus"
genf() { gen "$@"; tail -n 1 "$manifest" >> "$focus"; sed -i '$d' "$manifest"; }
genf "hash-a" 15 10 0.2 3 1500 "K1ABC W9XYZ EN37"
genf "hash-b" 15 10 0.2 9 1500 "<W9XYZ> PJ4/K1ABC RRR"
n=$((n + 1))
stem=$(printf 'hash-ab-261106_12%02d%02d' $((n / 60)) $((n % 60)))
python3 -I - "$out" "$(tail -n 2 "$focus" | head -n 1 | cut -d' ' -f1)" "$(tail -n 1 "$focus" | cut -d' ' -f1)" "$stem" <<'PYEOF'
import struct, sys
out, a, b, stem = sys.argv[1:5]
def rd(name):
    d = open(f"{out}/{name}.wav", "rb").read()
    i = d.index(b"data")
    n = struct.unpack("<I", d[i + 4:i + 8])[0] // 2
    return d[:i + 8], struct.unpack(f"<{n}h", d[i + 8:i + 8 + 2 * n])
h, x = rd(a)
_, y = rd(b)
z = [max(-32768, min(32767, p + q)) for p, q in zip(x, y)]
open(f"{out}/{stem}.wav", "wb").write(h + struct.pack(f"<{len(z)}h", *z))
PYEOF
echo "$stem 15" >> "$focus"
genf "freq-edge-f0800" 15 8 0.2 5 800 "$M"
genf "early-ping" 15 10 0.2 0.15 1500 "$M"
genf "offgrid-f1537-sm06" 15 -6 1 5 1537 "$M"
genf "offgrid-f1537-sm04" 15 -4 1 5 1537 "$M"
genf "swl-state-a" 15 10 0.2 5 1500 "N0ABC W1AW EN34"
genf "swl-state-b" 15 8 0.2 2,7,11 1500 "<N0ABC W1AW> RRR"
n=$((n + 1))
stem=$(printf 'swl-pair-261106_12%02d%02d' $((n / 60)) $((n % 60)))
python3 -I - "$out" "$(tail -n 2 "$focus" | head -n 1 | cut -d' ' -f1)" "$(tail -n 1 "$focus" | cut -d' ' -f1)" "$stem" <<'PYEOF'
import struct, sys
out, a, b, stem = sys.argv[1:5]
def rd(name):
    d = open(f"{out}/{name}.wav", "rb").read()
    i = d.index(b"data")
    n = struct.unpack("<I", d[i + 4:i + 8])[0] // 2
    return d[:i + 8], struct.unpack(f"<{n}h", d[i + 8:i + 8 + 2 * n])
h, x = rd(a)
_, y = rd(b)
z = [max(-32768, min(32767, p + q)) for p, q in zip(x, y)]
open(f"{out}/{stem}.wav", "wb").write(h + struct.pack(f"<{len(z)}h", *z))
PYEOF
echo "$stem 15" >> "$focus"
echo "msk144-sim-grid: $(grep -c . "$manifest") files in $out, $(grep -c . "$focus") more in MANIFEST.focus" >&2
