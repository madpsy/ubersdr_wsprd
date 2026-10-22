#!/usr/bin/env bash
# The FT8/FT4/FT2 decoders and the WSPR decoder (and jt9's and wsprd's
# command lines) built to WebAssembly with Emscripten, as two npm packages:
#
#   wasm/dist/jt9 (ubersdr-jt9-wasm)
#     ft.mjs + ft.wasm     the FT decoders for JavaScript (src/ft_wasm.cpp:
#                          new Decoder("ft8"|"ft4"|"ft2").decode(samples,
#                          utc, opts)), browsers, workers and Node
#     jt9.cjs + jt9.wasm   our jt9 (jt9_cli.h) for Node, on the real file
#                          system: node jt9.cjs -8 -d 2 FILE.wav
#   wasm/dist/wsprd (ubersdr-wsprd-wasm)
#     wspr.mjs + wspr.wasm the WSPR decoder for JavaScript (src/wspr_wasm.cpp:
#                          new WsprDecoder().decode(samples, opts))
#     wsprd.cjs + wsprd.wasm  our wsprd (wsprd.h's run) for Node: node
#                          wsprd.cjs -f 14.0956 FILE.wav
#
# FFTW 3.3.10 (single precision, no SIMD: none for wasm) is built once into
# wasm/build/fftw. Needs emsdk (EMSDK, or ~/emsdk).
#
#   wasm/build.sh                build
#   wasm/build.sh --pack         build, then npm pack each package (wasm/dist/*/*.tgz)
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

# The FT decoders' sources (what jt9_cli.cpp reaches), jt9_shm.cpp replaced.
decoders=(
    fft.cpp fst4w.cpp fst4w_frontend.cpp fst4w_jt9.cpp fst4w_ldpc.cpp fst4w_unpack.cpp
    ft2.cpp ft2_jt9.cpp ft4.cpp ft4_frontend.cpp ft4_jt9.cpp
    ft8.cpp ft8_frontend.cpp ft8_jt9.cpp ft8_ldpc.cpp ft8_maths.cpp
    ft_fast_subtract.cpp packjt77.cpp wsprd_osd.cpp
)
for f in ft2_sync.cpp ft4_candidates.cpp ft4_metrics.cpp ft4_sync.cpp ft4_wave.cpp; do
    [ -f "$root/src/$f" ] && decoders+=("$f")
done
# wsprd's (fft_wisdom.cpp: no wisdom file of its own; wsprd_osd.cpp above).
wsprd_sources=(wsprd.cpp wsprd_fano.cpp wsprd_unpack.cpp wsprd_files.cpp)
wisdom=$root/src/fft_wisdom.cpp
[ -f "$root/wsprd/src/fft_wisdom.cpp" ] && wisdom=$root/wsprd/src/fft_wisdom.cpp

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
wsprd_objs=()
for f in "${wsprd_sources[@]}"; do
    o=$objs/${f%.cpp}.o
    wsprd_objs+=("$o")
    compile "$root/src/$f" "$o" &
done
compile "$wisdom" "$objs/fft_wisdom.o" &
compile "$here/src/ft_wasm.cpp" "$objs/ft_wasm.o" &
compile "$here/src/wspr_wasm.cpp" "$objs/wspr_wasm.o" &
compile "$here/src/jt9_shm_none.cpp" "$objs/jt9_shm_none.o" &
compile "$root/src/jt9_cli.cpp" "$objs/jt9_cli.o" &
em++ "${cflags[@]}" -DJT9_VERSION="\"$version\"" -c "$root/src/jt9_main.cpp" -o "$objs/jt9_main.o" &
em++ "${cflags[@]}" -DWSPRD_VERSION="\"$version\"" -c "$root/src/wsprd_main.cpp" -o "$objs/wsprd_main.o" &
wait

ldflags=(-O3 -fwasm-exceptions -sALLOW_MEMORY_GROWTH=1 -sSTACK_SIZE=4MB "$fftw/lib/libfftw3f.a")

es6=(-lembind -sMODULARIZE=1 -sEXPORT_ES6=1 -sENVIRONMENT=web,worker,node)
cli=(-sNODERAWFS=1 -sENVIRONMENT=node -sEXIT_RUNTIME=1)
rm -rf "${dist:?}"/*
mkdir -p "$dist/jt9" "$dist/wsprd"

em++ "$objs/ft_wasm.o" "${lib_objs[@]}" "${ldflags[@]}" "${es6[@]}" -sEXPORT_NAME=createFt -o "$dist/jt9/ft.mjs"
em++ "$objs/jt9_main.o" "$objs/jt9_cli.o" "$objs/jt9_shm_none.o" "${lib_objs[@]}" "${ldflags[@]}" "${cli[@]}" \
    -o "$dist/jt9/jt9.cjs"
em++ "$objs/wspr_wasm.o" "${wsprd_objs[@]}" "$objs/wsprd_osd.o" "$objs/fft.o" "$objs/fft_wisdom.o" \
    "${ldflags[@]}" "${es6[@]}" -sEXPORT_NAME=createWspr -o "$dist/wsprd/wspr.mjs"
em++ "$objs/wsprd_main.o" "${wsprd_objs[@]}" "$objs/wsprd_osd.o" "$objs/fft.o" "$objs/fft_wisdom.o" \
    "${ldflags[@]}" "${cli[@]}" -o "$dist/wsprd/wsprd.cjs"

# package.json's: vYYYY.MM.DD[-N] as YYYY.MM.DD[-N], anything else 0.0.0-dev
semver=$(sed -nE 's/^v([0-9]+\.[0-9]+\.[0-9]+(-[0-9]+)?)$/\1/p' <<<"$version")
for p in jt9 wsprd; do
    # run as a program too (npx jt9-wasm / wsprd-wasm, package.json's bin)
    sed -i '1i #!/usr/bin/env node' "$dist/$p/$p.cjs"
    chmod +x "$dist/$p/$p.cjs"
    cp "$here/package/$p/"* "$dist/$p/"
    cp "$root/LICENSE" "$dist/$p/LICENSE"
    sed -i "s/\"version\": \"[^\"]*\"/\"version\": \"${semver:-0.0.0-dev}\"/" "$dist/$p/package.json"
done

ls -l "$dist"/*
if $pack; then for p in jt9 wsprd; do (cd "$dist/$p" && npm pack --silent); done; fi
