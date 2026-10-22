#!/usr/bin/env bash
# The WebAssembly build (wasm/dist) against a native jt9 of the same source,
# byte for byte, on the FT8, FT4 and FT2 recordings here:
#
#   - jt9.cjs (the command line under Node) and the native jt9 on each
#     recording set as one run, at depths 1 to 3: stdout and decoded.txt
#   - ft.mjs (the JavaScript API, wasm/check-api.mjs) given what jt9's file
#     mode reads and its parameters: every period's lines as the native
#     jt9's
#
#   wasm/check.sh NATIVE_JT9 [WORK]
set -euo pipefail

here=$(cd "$(dirname "$0")" && pwd)
root=$(cd "$here/.." && pwd)
native=$(readlink -f "$1")
work=${2:-$here/build/check}
[ -x "$native" ] || { echo "no native jt9 at $1" >&2; exit 2; }
[ -f "$here/dist/jt9.cjs" ] || { echo "no wasm/dist: run wasm/build.sh first" >&2; exit 2; }
command -v node >/dev/null || { echo "node is needed" >&2; exit 2; }

failed=0
run_set() {  # name, mode flag, depth, files...
    local name=$1 flag=$2 d=$3; shift 3
    local dir=$work/$name-d$d
    rm -rf "$dir"
    mkdir -p "$dir/native" "$dir/wasm"
    (cd "$dir/native" && "$native" "$flag" -d "$d" "$@" > stdout)
    (cd "$dir/wasm" && node "$here/dist/jt9.cjs" "$flag" -d "$d" "$@" > stdout)
    local ok=true
    for f in stdout decoded.txt; do
        if [ -f "$dir/native/$f" ] || [ -f "$dir/wasm/$f" ]; then
            cmp -s "$dir/native/$f" "$dir/wasm/$f" || ok=false
        fi
    done
    if $ok && node "$here/check-api.mjs" "$flag" "$d" "$dir/native/stdout" "$@"; then
        echo "  $name -d $d: identical ($(grep -vc '^<DecodeFinished>' "$dir/native/stdout") lines)"
    else
        echo "  $name -d $d: DIFFERS ($dir)"
        failed=1
    fi
}

sets=()
for s in ft8-sim:-8 ft8-real:-8 ft4-sim:-5 ft4-real:-5 ft2-sim:--ft2; do
    [ -n "$(ls "$root/recordings/${s%%:*}"/*.wav 2>/dev/null)" ] && sets+=("$s")
done
[ ${#sets[@]} -gt 0 ] || { echo "no FT8, FT4 or FT2 recordings under recordings/" >&2; exit 2; }
echo "the wasm build against $native"
for s in "${sets[@]}"; do
    mapfile -t wavs < <(ls "$root/recordings/${s%%:*}"/*.wav)
    for d in 1 2 3; do run_set "${s%%:*}" "${s#*:}" "$d" "${wavs[@]}"; done
done
[ "$failed" = 0 ] && echo "the wasm build decodes as the native jt9" || { echo "the wasm build differs" >&2; exit 1; }
