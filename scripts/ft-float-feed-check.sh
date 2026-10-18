#!/usr/bin/env bash
# jt9's FT8 (-8) and FT4 (-5) on 32-bit float .wav input and through --feed,
# against jt9's file mode on the 16-bit files. Each set of recordings below
# is decoded as one run (one decoder for all its files, in name order, so the
# hashes and a7 tables carried from period to period are checked too), under
# each option set, four ways:
#   int:        jt9 -8 OPTS FILE...            the 16-bit files (jt9's own way)
#   float:      jt9 -8 OPTS FILE...            their float copies, every sample
#                                              s as s/32768 (exact in a float)
#   feed:       FILE... on stdin | jt9 -8 OPTS --feed      the 16-bit files
#   feed-float: the same, the float copies
# float must leave the files int leaves, byte for byte (stdout, decoded.txt,
# the exit status); feed and feed-float the same, less their "<FileDone> PATH
# N" lines, one a file in order. Each run starts from a fresh data dir that
# is also the cwd.
#
#   scripts/ft-float-feed-check.sh JT9 [WORK]
#
# Prints a line a run that differs, then the counts; exits 1 on any.
set -uo pipefail
(( $# >= 1 )) || { echo "usage: $0 JT9 [WORK]" >&2; exit 2; }
jt9=$(realpath "$1")
repo=$(cd "$(dirname "$0")/.." && pwd)
work=$(realpath -m "${2:-$(mktemp -d "${TMPDIR:-/tmp}/ft-float-feed.XXXXXX")}")
rm -rf "$work"; mkdir -p "$work"

# The float copy, as float-input-check.sh writes it.
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

pass=0 fail=0
# MODE SET: each recordings/SET's .wav files as one run.
check_set() {
    local mode=$1 set=$2 opts way f bad
    local -a files=()
    for f in "$repo/recordings/$set"/*.wav; do [[ -f $f ]] && files+=("${f##*/}"); done
    (( ${#files[@]} )) || return 0
    mkdir -p "$work/$set/int" "$work/$set/float"
    for f in "${files[@]}"; do
        cp "$repo/recordings/$set/$f" "$work/$set/int/$f"
        to_float "$work/$set/int/$f" "$work/$set/float/$f"
    done
    for opts in "-d 3" "-d 2 -L 100 -H 4000" "-d 1"; do
        for way in int float feed feed-float; do
            local run=$work/run/$way src=$work/$set/${way#feed-}
            [[ $way == feed ]] && src=$work/$set/int
            rm -rf "$run"; mkdir -p "$run"
            cp "$src"/*.wav "$run/"
            # shellcheck disable=SC2086
            if [[ $way == feed* ]]; then
                (cd "$run" && printf '%s\n' "${files[@]}" | "$jt9" "$mode" -a . $opts --feed > stdout.raw 2> stderr
                 echo $? > status
                 grep -v '^<FileDone> ' stdout.raw > stdout
                 grep '^<FileDone> ' stdout.raw | cut -d' ' -f2 > done.txt)
                if ! printf '%s\n' "${files[@]}" | cmp -s - "$run/done.txt"; then
                    echo "differs: jt9 $mode $opts recordings/$set ($way): its <FileDone> lines"; fail=$((fail + 1))
                fi
                rm -f "$run/stdout.raw" "$run/done.txt"
            else
                (cd "$run" && "$jt9" "$mode" -a . $opts "${files[@]}" > stdout 2> stderr; echo $? > status)
            fi
            rm -f "$run"/*.wav "$run/jt9_wisdom.dat"
        done
        for way in float feed feed-float; do
            bad=
            for f in stdout decoded.txt status; do
                [[ -e $work/run/int/$f || -e $work/run/$way/$f ]] || continue
                cmp -s "$work/run/int/$f" "$work/run/$way/$f" || bad="$bad $f"
            done
            if [[ -n $bad ]]; then
                echo "differs: jt9 $mode $opts recordings/$set ($way):$bad"; fail=$((fail + 1))
            else
                pass=$((pass + 1))
            fi
        done
    done
}

for set in ft8-real ft8-sim ft8-filemode; do check_set -8 "$set"; done
for set in ft4-real ft4-sim ft4-filemode; do check_set -5 "$set"; done
echo "FT8/FT4 float input and --feed: $((pass + fail)) runs as jt9's on the 16-bit files, $fail differing"
(( fail == 0 ))
