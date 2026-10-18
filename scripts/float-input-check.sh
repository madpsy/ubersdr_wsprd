#!/usr/bin/env bash
# wsprd and jt9 on 32-bit float .wav input, against the same samples as
# 16-bit: each 16-bit mono .wav of recordings/40m-wspr-20261010-1914 and
# recordings/wsprd-sim is written again as WAVE_FORMAT_IEEE_FLOAT, every
# sample s as s/32768 (exact in a float), and each run on the float file must
# leave the files the run on the 16-bit file leaves, byte for byte: wsprd's
# stdout, ALL_WSPR.TXT, wspr_spots.txt, hashtable.txt and 000000_0001.c2 under
# each of its option sets below, and jt9 --fst4w's stdout, decoded.txt,
# fst4_decodes.dat and fst4w_calls.txt: every file but wspr_timer.out, whose
# CPU times differ from run to run. Each run starts from a fresh data dir that
# is also the cwd.
#
#   scripts/float-input-check.sh WSPRD JT9 [WORK]
#
# Prints a line a run that differs, then the counts; exits 1 on any.
set -uo pipefail
(( $# >= 2 )) || { echo "usage: $0 WSPRD JT9 [WORK]" >&2; exit 2; }
wsprd=$(realpath "$1")
jt9=$(realpath "$2")
repo=$(cd "$(dirname "$0")/.." && pwd)
work=$(realpath -m "${3:-$(mktemp -d "${TMPDIR:-/tmp}/float-input.XXXXXX")}")
rm -rf "$work"; mkdir -p "$work/float"

# The float copy: the 16-bit file's samples, each s/32768, under a 44-byte
# WAVE_FORMAT_IEEE_FLOAT header (as sox and radiod write float .wav files,
# less their optional chunks).
to_float() {
    python3 -I - "$1" "$2" <<'EOF'
import struct, sys, wave
with wave.open(sys.argv[1], 'rb') as w:
    assert w.getnchannels() == 1 and w.getsampwidth() == 2
    rate, raw = w.getframerate(), w.readframes(w.getnframes())
n = len(raw) // 2
data = struct.pack(f'<{n}f', *(s / 32768.0 for s in struct.unpack(f'<{n}h', raw)))
with open(sys.argv[2], 'wb') as f:
    f.write(b'RIFF' + struct.pack('<I', 36 + len(data)) + b'WAVE' + b'fmt ' +
            struct.pack('<IHHIIHH', 16, 3, 1, rate, rate * 4, 4, 32) + b'data' + struct.pack('<I', len(data)) + data)
EOF
}

# One run of PROGRAM ARGS... on the 16-bit file and on its float copy; NAME
# labels it. 0 when they leave the same files.
same() {
    local name=$1 wav=$2 prog=$3; shift 3
    local base=${wav##*/} d f bad=
    for d in int float; do
        rm -rf "$work/run/$d"; mkdir -p "$work/run/$d"
        if [[ $d == int ]]; then cp "$wav" "$work/run/$d/$base"; else cp "$work/float/$base" "$work/run/$d/$base"; fi
        (cd "$work/run/$d" && "$prog" -a . "$@" "$base" > stdout 2> stderr; echo $? > status)
        rm -f "$work/run/$d/$base"
    done
    for f in $(cd "$work/run/int" && ls) $(cd "$work/run/float" && ls); do
        [[ $f == wspr_timer.out ]] && continue
        cmp -s "$work/run/int/$f" "$work/run/float/$f" || bad="$bad $f"
    done
    if [[ -n $bad ]]; then echo "differs: $name:$(tr ' ' '\n' <<< "$bad" | sort -u | tr '\n' ' ')"; return 1; fi
    return 0
}

pass=0 fail=0
for wav in "$repo"/recordings/40m-wspr-20261010-1914/*.wav "$repo"/recordings/wsprd-sim/*.wav; do
    [[ -f $wav ]] || continue
    to_float "$wav" "$work/float/${wav##*/}"
    for opts in "" "-c -C 500 -o 4 -d" "-c -n" "-w -s -B"; do
        # shellcheck disable=SC2086
        if same "wsprd $opts ${wav#"$repo"/}" "$wav" "$wsprd" $opts; then pass=$((pass + 1)); else fail=$((fail + 1)); fi
    done
    if same "jt9 --fst4w ${wav#"$repo"/}" "$wav" "$jt9" -p 120 --fst4w -f 1500 -F 100 -d 3; then pass=$((pass + 1)); else fail=$((fail + 1)); fi
done
echo "float input: $((pass + fail)) runs as on 16-bit input, $fail differing"
(( fail == 0 ))
