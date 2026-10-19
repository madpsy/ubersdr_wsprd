#!/usr/bin/env bash
# Build wsprd and jt9.
#
#   ./build.sh            build wsprd and jt9 for this machine, Linux on
#                         x86-64 or arm64 (a Raspberry Pi too): build/wsprd,
#                         build/jt9, and stripped, build/dist/wsprd-linux-<arch>
#                         and build/dist/jt9-linux-<arch>
#   ./build.sh --windows  the same, and wsprd.exe and jt9.exe for 64-bit Windows,
#                         cross-compiled here on Linux in Docker (mingw-w64):
#                         build/dist/wsprd-windows-x86_64.exe, checked under
#                         wine when it is installed
#   ./build.sh --check    the same as ./build.sh, then this wsprd and jt9
#                         against WSJT-X's own's stored outputs, byte for byte
#                         (scripts/wsprd-stored-check.sh, scripts/jt9-wine-check.sh
#                         run natively; about four minutes, nothing of WSJT-X
#                         built), and their 32-bit float input against the same
#                         samples as 16-bit (scripts/float-input-check.sh). On arm64 jt9's known differences (FFTW's
#                         NEON rounding: tests/data/jt9-arm64-accepted.txt)
#                         print as "accepted" and pass; only "differs:"
#                         fails
#   ./build.sh --deep-check
#                         the same as ./build.sh, then this wsprd and jt9
#                         run beside WSJT-X's own, built here, on every option
#                         and recording (scripts/wsprd-compare.sh,
#                         scripts/jt9-compare.sh; longer)
#   ./build.sh --version TAG
#                         built as release TAG (wsprd --version, jt9 -v); otherwise
#                         dev-<commit>
#
# Releases (the maintainer's):
#   ./build.sh --release [TAG]
#                         --check and --windows (--deep-check too when given),
#                         then release both as TAG on
#                         GitHub (default vYYYY.MM.DD, with -2, -3... if that
#                         is taken), and as the release "latest" (its tag
#                         moved to TAG's commit, its files replaced): the
#                         newest's files always at
#                         .../releases/download/latest/<file>. Earlier
#                         releases are kept.
#   ./build.sh --add-built TAG DIR
#                         no build: add the wsprd and jt9 another machine built (a
#                         copy of its build/dist, made with --version TAG) to
#                         release TAG, and to "latest" while that is TAG: an
#                         arm64 build from a Pi, for one
#
# The deep check needs WSJT-X's source at the decoder's commit (WSJTX, default
# ~/repos/wsjtx, at 567ad29ce; scripts/wsprd-oracle.sh --build-only builds its
# wsprd with WSJT-X's Release flags; scripts/fst4w-oracle.sh builds its jt9
# with WSJT-X's own CMake in Docker, and makes jt9's inputs) and takes some
# minutes (jt9's the longer, building WSJT-X); WSPRD_SKIP_COMPARE=1 leaves
# wsprd's out, JT9_SKIP_COMPARE=1 jt9's. Publishing needs the GitHub CLI logged
# in (gh auth login), a clean working tree and HEAD pushed, so what is
# released is a commit on GitHub; Windows needs Docker.
#
# Debian or Ubuntu: cmake, g++, pkg-config and libfftw3-dev (and gfortran and
# Docker for the deep check).
set -euo pipefail

cd "$(dirname "$0")"
root=$(pwd)
build=$root/build

release=false
windows=false
check=false
deep=false
added=""
tag=""
version=""
while [ $# -gt 0 ]; do
    case "$1" in
        --release)
            release=true
            if [ $# -gt 1 ] && [ "${2#-}" = "$2" ]; then tag=$2; shift; fi
            ;;
        --windows) windows=true ;;
        --check) check=true ;;
        --deep-check) deep=true ;;
        --version)
            [ $# -gt 1 ] || { echo "--version needs the release's tag" >&2; exit 2; }
            version=$2
            shift
            ;;
        --add-built)
            [ $# -gt 2 ] || { echo "--add-built needs the release's tag and a directory" >&2; exit 2; }
            tag=$2
            added=$3
            shift 2
            ;;
        -h|--help)
            sed -n '2,55p' "$0" | sed 's/^# \{0,1\}//'
            exit 0
            ;;
        *)
            echo "unknown argument: $1 (see --help)" >&2
            exit 2
            ;;
    esac
    shift
done
if $release; then windows=true; check=true; fi
if $release && [ -n "$version" ]; then
    echo "--version is for a plain build; a release is built as its tag" >&2
    exit 2
fi

say() { printf '\n== %s\n' "$*"; }

# A dist directory's files for one platform (linux-x86_64, linux-aarch64,
# windows-x86_64).
dist_assets() {
    local dir=$1 platform=$2
    case "$platform" in
        windows-*) echo "$dir/wsprd-$platform.exe" "$dir/jt9-$platform.exe" "$dir/SHA256SUMS-$platform" ;;
        *) echo "$dir/wsprd-$platform" "$dir/jt9-$platform" "$dir/SHA256SUMS-$platform" ;;
    esac
}

# The notes' section for one platform: where it was built, what it needs,
# and the checksum.
section_for() {
    local dir=$1 platform=$2 needs built_on fence='```'
    case "$platform" in
        linux-x86_64) needs="an x86-64 CPU with AVX2 (Intel since Haswell, AMD since Zen); nothing else to install" ;;
        linux-aarch64) needs="64-bit ARM (a Raspberry Pi 4 or 5 on a 64-bit OS, for one); nothing else to install" ;;
        windows-x86_64) needs="64-bit Windows 10 or 11 on an x86-64 CPU with AVX2; nothing else to install" ;;
        *) needs="an ${platform#*-} CPU" ;;
    esac
    built_on=$dir/built-on
    [ -f "$dir/built-on-$platform" ] && built_on=$dir/built-on-$platform
    printf '### %s\n\nBuilt on %s. Needs %s.\n\n%s\n%s\n%s' "$platform" "$(cat "$built_on")" "$needs" \
        "$fence" "$(cat "$dir/SHA256SUMS-$platform")" "$fence"
}

# "latest": the newest release's files under one tag that does not change.
# Its title names the release it holds.
latest_holds() {
    gh release view latest --json name -q .name 2>/dev/null | sed -n 's/^latest (\(.*\))$/\1/p'
}
make_latest() {
    local tag=$1 commit=$2 notes=$3; shift 3
    git tag -f latest "$commit" >/dev/null
    git push -q -f origin refs/tags/latest
    # Made again each time, not edited: a release's date is when it was
    # first published, and an edited one would go on showing the first.
    if gh release view latest >/dev/null 2>&1; then gh release delete latest -y; fi
    gh release create latest "$@" --verify-tag --title "latest ($tag)" --notes "$notes" --latest
    echo "  $(gh release view latest --json url -q .url)"
}

if [ -n "$added" ]; then
    gh auth status >/dev/null 2>&1 || { echo "gh is not logged in: run gh auth login" >&2; exit 1; }
    gh release view "$tag" >/dev/null 2>&1 || { echo "there is no release $tag" >&2; exit 1; }
    arch=$(cd "$added" && ls wsprd-linux-* 2>/dev/null | head -1 | sed 's/^wsprd-linux-//')
    [ -n "$arch" ] || { echo "$added holds no wsprd-linux-<arch>" >&2; exit 1; }
    [ -f "$added/jt9-linux-$arch" ] || { echo "$added holds no jt9-linux-$arch" >&2; exit 1; }
    for f in $(dist_assets "$added" "linux-$arch") "$added/built-on"; do
        [ -f "$f" ] || { echo "$added has no $(basename "$f")" >&2; exit 1; }
    done
    (cd "$added" && sha256sum --quiet -c "SHA256SUMS-linux-$arch") || { echo "$added: checksums do not match" >&2; exit 1; }
    built=$(cat "$added/version" 2>/dev/null || true)
    if [ "$built" != "$tag" ]; then
        echo "$added was built as ${built:-no version}, not $tag: build it with ./build.sh --version $tag" >&2
        exit 1
    fi
    say "adding linux-$arch from $added to $tag"
    # shellcheck disable=SC2046
    gh release upload "$tag" $(dist_assets "$added" "linux-$arch") --clobber
    body=$(gh release view "$tag" --json body -q .body)
    if ! grep -qx "### linux-$arch" <<<"$body"; then
        gh release edit "$tag" --notes "$body

$(section_for "$added" "linux-$arch")"
    fi
    echo "  $(gh release view "$tag" --json url -q .url)"
    if [ "$(latest_holds)" = "$tag" ]; then
        # shellcheck disable=SC2046
        gh release upload latest $(dist_assets "$added" "linux-$arch") --clobber
        gh release edit latest --notes "$(gh release view "$tag" --json body -q .body)"
        echo "  and latest"
    fi
    exit 0
fi

# --- release preconditions, before spending time on the build ----------------
if $release; then
    say "checking the release can be made"
    gh auth status >/dev/null 2>&1 || { echo "gh is not logged in: run gh auth login" >&2; exit 1; }
    docker info >/dev/null 2>&1 || { echo "a release carries wsprd.exe, built in Docker: Docker is not usable here" >&2; exit 1; }
    if [ -n "$(git status --porcelain)" ]; then
        echo "the working tree has changes; commit or stash them first:" >&2
        git status --short >&2
        exit 1
    fi
    git fetch --quiet origin
    branch=$(git rev-parse --abbrev-ref HEAD)
    if ! git merge-base --is-ancestor HEAD "origin/$branch" 2>/dev/null; then
        echo "HEAD ($(git rev-parse --short HEAD)) is not pushed to origin/$branch; push it first" >&2
        exit 1
    fi
    if [ -z "$tag" ]; then
        # Today's date, numbered past any release or tag of today already
        # made, so the newest always sorts last.
        base="v$(date -u +%Y.%m.%d)"
        last=0
        for t in $(gh release list --limit 200 --json tagName -q '.[].tagName' 2>/dev/null) $(git tag -l "$base*"); do
            case "$t" in
                "$base") [ "$last" -lt 1 ] && last=1 ;;
                "$base"-*) n=${t#"$base"-}; [ "$n" -eq "$n" ] 2>/dev/null && [ "$n" -gt "$last" ] && last=$n ;;
            esac
        done
        if [ "$last" -eq 0 ]; then tag=$base; else tag="$base-$((last + 1))"; fi
    elif gh release view "$tag" >/dev/null 2>&1; then
        echo "release $tag already exists" >&2
        exit 1
    fi
    [ "$tag" != latest ] || { echo "\"latest\" is the newest release's own tag; name the release by its date" >&2; exit 1; }
    echo "will release $(git rev-parse --short HEAD) as $tag"
    version=$tag
fi
[ -n "$version" ] || version="dev-$(git rev-parse --short HEAD 2>/dev/null || echo unknown)"

# --- build (Linux, this machine) ---------------------------------------------
say "building wsprd and jt9 $version (Release)"
rm -rf "$build"
cmake -S "$root" -B "$build" -DCMAKE_BUILD_TYPE=Release -DWSPRD_VERSION="$version"
cmake --build "$build" -j "$(nproc)"
"$build/wsprd" --version | grep -qF "wsprd $version " || { echo "wsprd was not built as $version" >&2; exit 1; }
"$build/jt9" -v | grep -qF "jt9 $version " || { echo "jt9 was not built as $version" >&2; exit 1; }
# FFTW is linked in (CMakeLists.txt, libfftw3-dev's libfftw3f.a): a binary
# that would need it installed is not one to hand out.
for p in wsprd jt9; do
    if ldd "$build/$p" | grep -q libfftw3; then
        echo "$p needs FFTW's shared library: FFTW's static one (libfftw3f.a, in libfftw3-dev) was not found" >&2
        exit 1
    fi
done

arch=$(uname -m)
dist=$build/dist
mkdir -p "$dist"
for p in wsprd jt9; do
    cp "$build/$p" "$dist/$p-linux-$arch"
    strip "$dist/$p-linux-$arch"
done
(cd "$dist" && sha256sum -- "wsprd-linux-$arch" "jt9-linux-$arch" > "SHA256SUMS-linux-$arch")
(. /etc/os-release && echo "$PRETTY_NAME") > "$dist/built-on"
echo "$version" > "$dist/version"
echo "built: $dist/wsprd-linux-$arch, $dist/jt9-linux-$arch ($version)"

# --- the check: WSJT-X's own's stored outputs ---------------------------------
if $check && [ "${WSPRD_SKIP_COMPARE:-0}" != 1 ]; then
    say "checking wsprd against WSJT-X's wsprd's stored outputs"
    "$root/scripts/wsprd-stored-check.sh" "$build/wsprd" "$build/stored-wsprd"
fi
if $check && [ "${JT9_SKIP_COMPARE:-0}" != 1 ]; then
    say "checking jt9 against WSJT-X's jt9's stored outputs"
    "$root/scripts/jt9-wine-check.sh" "$build/jt9" "$build/stored-jt9"
fi
if $check && [ "${WSPRD_SKIP_COMPARE:-0}" != 1 ] && [ "${JT9_SKIP_COMPARE:-0}" != 1 ]; then
    say "checking wsprd and jt9 on 32-bit float input against the same samples as 16-bit"
    "$root/scripts/float-input-check.sh" "$build/wsprd" "$build/jt9" "$build/float-input"
    say "checking jt9's FT8, FT4 and FT2 on 32-bit float input and through --feed against its 16-bit file mode"
    "$root/scripts/ft-float-feed-check.sh" "$build/jt9" "$build/ft-float-feed"
fi

# --- the deep check: WSJT-X's own wsprd, every option and recording -----------
if $deep && [ "${WSPRD_SKIP_COMPARE:-0}" != 1 ]; then
    say "comparing with WSJT-X's wsprd"
    WSPRD_ORACLE_WORK=$build/stock "$root/scripts/wsprd-oracle.sh" --build-only
    WSPRD_COMPARE_WORK=$build/compare "$root/scripts/wsprd-compare.sh" "$build/wsprd" "$build/stock/bin/wsprd-release"
fi
if $deep && [ "${JT9_SKIP_COMPARE:-0}" != 1 ]; then
    say "comparing with WSJT-X's jt9"
    # WSJT-X's jt9 (built in Docker, image wsjtx-oracle) and the inputs.
    FST4W_ORACLE_WORK=$build/stock-jt9 "$root/scripts/fst4w-oracle.sh" --sim-only
    # On arm64 fort.52 (the spread's debug spectrum) differs in its last
    # printed digit: FFTW's NEON code rounds differently from the reference's.
    ignore=; [ "$(uname -m)" = x86_64 ] || ignore=fort.52
    JT9_COMPARE_IGNORE=$ignore JT9_COMPARE_WORK=$build/compare-jt9 \
        "$root/scripts/jt9-compare.sh" "$build/jt9" "$build/stock-jt9/release/jt9" --docker wsjtx-oracle
fi

$windows || exit 0

# --- Windows (Docker, mingw-w64) -----------------------------------------------
command -v docker >/dev/null && docker info >/dev/null 2>&1 ||
    { echo "--windows builds wsprd.exe in Docker: Docker is not usable here" >&2; exit 1; }
say "building wsprd.exe and jt9.exe (mingw-w64)"
image=ubersdr-wsprd-win
wbuild=$build/windows
docker build -q -t "$image" -f docker/windows.Dockerfile docker >/dev/null
mkdir -p "$wbuild"
docker run --rm -u "$(id -u):$(id -g)" -v "$root:$root:ro" -v "$wbuild:/b" "$image" sh -c "
    cmake -S '$root' -B /b -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_TOOLCHAIN_FILE=/win/toolchain.cmake \
          -DWSPRD_VERSION='$version' &&
    cmake --build /b"
for p in wsprd jt9; do
    file "$wbuild/$p.exe" | grep -q 'PE32+ executable.*x86-64' || { echo "$p.exe was not built" >&2; exit 1; }
    # One file: any DLL but Windows' own would have to be shipped beside it.
    other=$(docker run --rm -v "$wbuild:/b" "$image" x86_64-w64-mingw32-objdump -p "/b/$p.exe" |
        awk '/DLL Name:/ {print $3}' |
        grep -viE '^(kernel32|advapi32|bcrypt|crypt32|ws2_32|msvcrt|user32|shell32|ole32)\.dll$' || true)
    [ -z "$other" ] || { echo "$p.exe needs DLLs Windows does not have: $other" >&2; exit 1; }
    docker run --rm -u "$(id -u):$(id -g)" -v "$wbuild:/b" "$image" \
        x86_64-w64-mingw32-strip -o "/b/$p-windows-x86_64.exe" "/b/$p.exe"
done
# Under wine, when it is there: the same decodes as Linux's of a .c2 (whose
# FFT is the same on both; a .wav's first FFTs' last bits are FFTW's build's).
if command -v wine >/dev/null; then
    say "checking wsprd.exe under wine"
    c2=$root/tests/data/wsprd/inputs/sim-types-261010_1200.c2
    for side in linux windows; do
        rm -rf "$build/wine-$side"; mkdir -p "$build/wine-$side"; cp "$c2" "$build/wine-$side/"
    done
    (cd "$build/wine-linux" && "$build/wsprd" -C 500 -o 4 -d "$(basename "$c2")" > stdout)
    (cd "$build/wine-windows" && WINEPREFIX="$build/wineprefix" WINEDEBUG=-all \
        wine "$wbuild/wsprd-windows-x86_64.exe" -C 500 -o 4 -d "$(basename "$c2")" | tr -d '\r' > stdout)
    for f in stdout ALL_WSPR.TXT wspr_spots.txt; do
        cmp -s <(tr -d '\r' < "$build/wine-linux/$f") <(tr -d '\r' < "$build/wine-windows/$f") ||
            { echo "wsprd.exe's $f differs from Linux's" >&2; exit 1; }
    done
    echo "wsprd.exe decodes as Linux's does"
    # jt9.exe: every case of WSJT-X's jt9 whose inputs are here (FST4W's
    # made by --deep-check, or kept from an earlier one) printed and written
    # as WSJT-X's.
    say "checking jt9.exe under wine"
    "$root/scripts/jt9-wine-check.sh" "$wbuild/jt9-windows-x86_64.exe" "$build/wine-jt9"
else
    echo "wine is not installed: wsprd.exe and jt9.exe are built but not checked"
fi
cp "$wbuild/wsprd-windows-x86_64.exe" "$wbuild/jt9-windows-x86_64.exe" "$dist/"
(cd "$dist" && sha256sum -- wsprd-windows-x86_64.exe jt9-windows-x86_64.exe > SHA256SUMS-windows-x86_64)
compiler=$(docker run --rm "$image" x86_64-w64-mingw32-g++ -dumpfullversion)
echo "$(. /etc/os-release && echo "$PRETTY_NAME"), cross-compiled with mingw-w64 GCC $compiler" > "$dist/built-on-windows-x86_64"
echo
echo "built: $dist/wsprd-windows-x86_64.exe, $dist/jt9-windows-x86_64.exe"
$release || exit 0

# --- publish -------------------------------------------------------------------
say "releasing $tag"
commit=$(git rev-parse HEAD)
notes="Built from $commit.

WSJT-X's \`wsprd\` (567ad29ce's decoder), a drop-in for it: the same options, output and files, byte for byte, in a tenth of its CPU on x86-64 and a sixth on a Raspberry Pi 5. One option more, \`-n\`: no search for a drifting signal, and each ALL_WSPR.TXT line ending in the signal's spread in Hz. \`wsprd --version\` gives this release.

And WSJT-X's \`jt9\` for FST4W (FST4W-120 to -1800, and 15, 30 and 60 s), FT8 (\`jt9 -8\`) and FT4 (\`jt9 -5\`), run on WAV files: the same command line, output and files, byte for byte, in about a seventh of its CPU for FST4W, 2.1 to 2.9 times less for FT8 and 2.9 to 4.3 times less for FT4 (x86-64, depths 1 to 3). And FT2 (\`jt9 --ft2\`) as WSJT-X improved 3.2.0's \`jt9\` runs it: improved's FT2 code on the FT4 decoder above, in 4.2 to 6.0 times less CPU than improved's \`jt9\` (x86-64, depths 1 to 3). One option more, \`--spread fast|exact|off\` (FST4W). FST4W, FT8, FT4 and FT2 only so far: keep WSJT-X's own \`jt9\` for WSJT-X and its other modes. \`jt9 -v\` gives this release.

- \`wsprd-linux-x86_64\`, \`wsprd-linux-aarch64\`, \`jt9-linux-x86_64\`, \`jt9-linux-aarch64\`: Linux.
- \`wsprd-windows-x86_64.exe\`, \`jt9-windows-x86_64.exe\`: Windows."
assets=()
for platform in "linux-$arch" windows-x86_64; do
    notes="$notes

$(section_for "$dist" "$platform")"
    # shellcheck disable=SC2207
    assets+=($(dist_assets "$dist" "$platform"))
done
gh release create "$tag" "${assets[@]}" --target "$commit" --title "$tag" --notes "$notes" --latest=false
echo "  $(gh release view "$tag" --json url -q .url)"
make_latest "$tag" "$commit" "$notes" "${assets[@]}"
