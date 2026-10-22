#!/usr/bin/env bash
# The FT8/FT4/FT2 decoders (and the jt9 command line) built to WebAssembly
# with Emscripten, as an npm package in wasm/dist:
#
#   ft.mjs + ft.wasm   the decoders for JavaScript (src/ft_wasm.cpp: new
#                      Decoder("ft8"|"ft4"|"ft2").decode(samples, utc, opts)),
#                      browsers, workers and Node
#   jt9.cjs + jt9.wasm our jt9 (jt9_cli.h) for Node, on the real file system:
#                      node jt9.cjs -8 -d 2 FILE.wav
#
# FFTW 3.3.10 (single precision, no SIMD: none for wasm) is built once into
# wasm/build/fftw. Needs emsdk (EMSDK, or ~/emsdk).
#
#   wasm/build.sh                build
#   wasm/build.sh --pack         build, then npm pack in wasm/dist
#   wasm/build.sh --version TAG  built as release TAG (jt9 -v, package.json)
set -euo pipefail

pack=false
version=""
while [ $# -gt 0 ]; do
    case "$1" in
        --pack) pack=true ;;
        --version) version=$2; shift ;;
        *) echo "unknown argument: $1" >&2; exit 2 ;;
    esac
    shift
done

here=$(cd "$(dirname "$0")" && pwd)
root=$(cd "$here/.." && pwd)
build=$here/build
dist=$here/dist

emsdk=${EMSDK:-$HOME/emsdk}
if ! command -v emcc >/dev/null; then
    # shellcheck disable=SC1091
    source "$emsdk/emsdk_env.sh" >/dev/null 2>&1
fi
command -v emcc >/dev/null || { echo "emcc not found: install emsdk (EMSDK=$emsdk)" >&2; exit 1; }

jobs=$(nproc 2>/dev/null || echo 4)
mkdir -p "$build" "$dist"

# FFTW, as docker/windows.Dockerfile fetches it.
fftw=$build/fftw
if [ ! -f "$fftw/lib/libfftw3f.a" ]; then
    sha=56c932549852cddcfafdab3820b0200c7742675be92179e59e6215b340e26467
    tgz=$build/fftw-3.3.10.tar.gz
    [ -f "$tgz" ] || curl -sSLf -o "$tgz" https://www.fftw.org/fftw-3.3.10.tar.gz
    echo "$sha  $tgz" | sha256sum -c - >/dev/null
    rm -rf "$build/fftw-src" && mkdir "$build/fftw-src"
    tar xzf "$tgz" -C "$build/fftw-src" --strip-components=1
    (cd "$build/fftw-src" &&
        emconfigure ./configure -q --prefix="$fftw" --enable-float --disable-shared --enable-static \
            --disable-fortran --disable-doc --disable-threads CFLAGS="-O3" >/dev/null &&
        emmake make -j"$jobs" >/dev/null && emmake make install >/dev/null)
    rm -rf "$build/fftw-src"
fi

# The decoders' sources (what jt9_cli.cpp reaches), jt9_shm.cpp replaced.
decoders=(
    fft.cpp fst4w.cpp fst4w_frontend.cpp fst4w_jt9.cpp fst4w_ldpc.cpp fst4w_unpack.cpp
    ft2.cpp ft2_jt9.cpp ft4.cpp ft4_frontend.cpp ft4_jt9.cpp
    ft8.cpp ft8_frontend.cpp ft8_jt9.cpp ft8_ldpc.cpp ft8_maths.cpp
    ft_fast_subtract.cpp packjt77.cpp wsprd_osd.cpp
)
for f in ft2_sync.cpp ft4_candidates.cpp ft4_metrics.cpp ft4_sync.cpp ft4_wave.cpp; do
    [ -f "$root/src/$f" ] && decoders+=("$f")
done

# -ffp-contract=off everywhere, as CMakeLists.txt builds the decoders: wasm
# has no FMA, so the floats are scalar x86-64's.
cflags=(-std=c++17 -O3 -ffp-contract=off -fwasm-exceptions -I"$root/src" -I"$fftw/include" -DNDEBUG)
[ -n "$version" ] || version="dev-$(git -C "$root" rev-parse --short HEAD 2>/dev/null || echo unknown)"

objs=$build/obj
mkdir -p "$objs"
compile() {  # source, object
    if [ ! -f "$2" ] || [ "$1" -nt "$2" ] || [ "$0" -nt "$2" ]; then
        em++ "${cflags[@]}" -c "$1" -o "$2"
    fi
}
pids=()
lib_objs=()
for f in "${decoders[@]}"; do
    o=$objs/${f%.cpp}.o
    lib_objs+=("$o")
    compile "$root/src/$f" "$o" &
    pids+=($!)
    if [ ${#pids[@]} -ge "$jobs" ]; then wait "${pids[0]}"; pids=("${pids[@]:1}"); fi
done
compile "$here/src/ft_wasm.cpp" "$objs/ft_wasm.o" &
compile "$here/src/jt9_shm_none.cpp" "$objs/jt9_shm_none.o" &
compile "$root/src/jt9_cli.cpp" "$objs/jt9_cli.o" &
em++ "${cflags[@]}" -DJT9_VERSION="\"$version\"" -c "$root/src/jt9_main.cpp" -o "$objs/jt9_main.o" &
wait

ldflags=(-O3 -fwasm-exceptions -sALLOW_MEMORY_GROWTH=1 -sSTACK_SIZE=4MB "$fftw/lib/libfftw3f.a")

em++ "$objs/ft_wasm.o" "${lib_objs[@]}" "${ldflags[@]}" -lembind \
    -sMODULARIZE=1 -sEXPORT_ES6=1 -sEXPORT_NAME=createFt -sENVIRONMENT=web,worker,node \
    -o "$dist/ft.mjs"

em++ "$objs/jt9_main.o" "$objs/jt9_cli.o" "$objs/jt9_shm_none.o" "${lib_objs[@]}" "${ldflags[@]}" \
    -sNODERAWFS=1 -sENVIRONMENT=node -sEXIT_RUNTIME=1 \
    -o "$dist/jt9.cjs"
# run as a program too (npx jt9-wasm, package.json's bin)
sed -i '1i #!/usr/bin/env node' "$dist/jt9.cjs"
chmod +x "$dist/jt9.cjs"

cp "$here/package/"* "$dist/"
cp "$root/LICENSE" "$dist/LICENSE"
# package.json's: vYYYY.MM.DD[-N] as YYYY.MM.DD[-N], anything else 0.0.0-dev
semver=$(sed -nE 's/^v([0-9]+\.[0-9]+\.[0-9]+(-[0-9]+)?)$/\1/p' <<<"$version")
sed -i "s/\"version\": \"[^\"]*\"/\"version\": \"${semver:-0.0.0-dev}\"/" "$dist/package.json"

ls -l "$dist"
if $pack; then (cd "$dist" && rm -f ./*.tgz && npm pack --silent); fi
