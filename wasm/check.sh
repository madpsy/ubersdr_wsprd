#!/usr/bin/env bash
# The WebAssembly builds (wasm/dist) against native jt9 and wsprd of the
# same source, byte for byte, on the recordings here:
#
#   - jt9.cjs (the command line under Node) and the native jt9 on each FT8,
#     FT4 and FT2 recording set as one run, at depths 1 to 3: stdout and
#     decoded.txt
#   - ft.mjs (the JavaScript API, wasm/check-api.mjs) given what jt9's file
#     mode reads and its parameters: every period's lines as the native
#     jt9's
#   - wsprd.cjs and the native wsprd on each WSPR recording set, its files
#     one after another in one directory, under several option sets: stdout,
#     ALL_WSPR.TXT, wspr_spots.txt and hashtable.txt; and wspr.mjs
#     (wasm/check-wspr-api.mjs) given each period's samples and the same
#     options: its lines as the native wsprd's stdout
#
#   wasm/check.sh NATIVE_JT9 NATIVE_WSPRD [WORK]
set -euo pipefail

here=$(cd "$(dirname "$0")" && pwd)
root=$(cd "$here/.." && pwd)
[ $# -ge 2 ] || { echo "usage: $0 NATIVE_JT9 NATIVE_WSPRD [WORK]" >&2; exit 2; }
native=$(readlink -f "$1")
native_wsprd=$(readlink -f "$2")
work=${3:-$here/build/check}
[ -x "$native" ] || { echo "no native jt9 at $1" >&2; exit 2; }
[ -x "$native_wsprd" ] || { echo "no native wsprd at $2" >&2; exit 2; }
[ -f "$here/dist/jt9/jt9.cjs" ] && [ -f "$here/dist/wsprd/wsprd.cjs" ] ||
    { echo "no wasm/dist: run wasm/build.sh first" >&2; exit 2; }
command -v node >/dev/null || { echo "node is needed" >&2; exit 2; }

failed=0
run_set() {  # name, mode flag, depth, files...
    local name=$1 flag=$2 d=$3; shift 3
    local dir=$work/$name-d$d
    rm -rf "$dir"
    mkdir -p "$dir/native" "$dir/wasm"
    (cd "$dir/native" && "$native" "$flag" -d "$d" "$@" > stdout)
    (cd "$dir/wasm" && node "$here/dist/jt9/jt9.cjs" "$flag" -d "$d" "$@" > stdout)
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

# wsprd: name, "options", files... (each file a run, one directory)
run_wsprd() {
    local name=$1 opts=$2; shift 2
    local dir=$work/wsprd-$name
    rm -rf "$dir"
    mkdir -p "$dir/native" "$dir/wasm"
    local f wavs=()
    for f in "$@"; do
        # shellcheck disable=SC2086
        (cd "$dir/native" && "$native_wsprd" $opts "$f" >> stdout)
        # shellcheck disable=SC2086
        (cd "$dir/wasm" && node "$here/dist/wsprd/wsprd.cjs" $opts "$f" >> stdout)
        case "$f" in *.wav) wavs+=("$f") ;; esac
    done
    local ok=true
    for f in stdout ALL_WSPR.TXT wspr_spots.txt hashtable.txt; do
        if [ -f "$dir/native/$f" ] || [ -f "$dir/wasm/$f" ]; then
            cmp -s "$dir/native/$f" "$dir/wasm/$f" || ok=false
        fi
    done
    if $ok && [ ${#wavs[@]} -gt 0 ] && [ ${#wavs[@]} = $# ]; then
        node "$here/check-wspr-api.mjs" "$dir/native/stdout" "$opts" "${wavs[@]}" || ok=false
    fi
    if $ok; then
        echo "  wsprd $name: identical ($(grep -vc '^<DecodeFinished>' "$dir/native/stdout") decodes)"
    else
        echo "  wsprd $name: DIFFERS ($dir)"
        failed=1
    fi
}

sets=()
for s in ft8-sim:-8 ft8-real:-8 ft4-sim:-5 ft4-real:-5 ft2-sim:--ft2; do
    [ -n "$(ls "$root/recordings/${s%%:*}"/*.wav 2>/dev/null)" ] && sets+=("$s")
done
[ ${#sets[@]} -gt 0 ] || { echo "no FT8, FT4 or FT2 recordings under recordings/" >&2; exit 2; }
echo "the wasm builds against $native and $native_wsprd"
for s in "${sets[@]}"; do
    mapfile -t wavs < <(ls "$root/recordings/${s%%:*}"/*.wav)
    for d in 1 2 3; do run_set "${s%%:*}" "${s#*:}" "$d" "${wavs[@]}"; done
done

# WSPR: the 40 m periods (7.0386 MHz), WSJT-X's sample (14.0956), the
# generated .c2 periods (their own dial; the command line only).
mapfile -t forty < <(ls "$root"/recordings/40m-wspr-*/*.wav 2>/dev/null)
sample=$root/recordings/wspr2-150426_0918-14.0956.wav
mapfile -t c2s < <(ls "$root"/recordings/wsprd-sim/*.c2 2>/dev/null)
for opts in "-C 500 -o 4" "-C 500 -o 4 -d" "" "-q" "-s" "-B" "-n -C 500 -o 4" "-w -H -C 500"; do
    tag=$(tr -d ' ' <<<"${opts:-defaults}")
    [ ${#forty[@]} -gt 0 ] && run_wsprd "40m$tag" "$opts -f 7.0386" "${forty[@]}"
    [ -f "$sample" ] && run_wsprd "sample$tag" "$opts -f 14.0956" "$sample"
    [ ${#c2s[@]} -gt 0 ] && run_wsprd "c2$tag" "$opts" "${c2s[@]}"
done

[ "$failed" = 0 ] && echo "the wasm builds decode as the native jt9 and wsprd" || { echo "the wasm builds differ" >&2; exit 1; }
