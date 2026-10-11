#!/usr/bin/env bash
# Build wsprd.
#
#   ./build.sh            build wsprd for this machine, Linux on x86-64 or
#                         arm64 (a Raspberry Pi too): build/wsprd, and
#                         stripped, build/dist/wsprd-linux-<arch>
#   ./build.sh --windows  the same, and wsprd.exe for 64-bit Windows,
#                         cross-compiled here on Linux in Docker (mingw-w64):
#                         build/dist/wsprd-windows-x86_64.exe, checked under
#                         wine when it is installed
#   ./build.sh --check    the same as ./build.sh, then this wsprd compared
#                         with WSJT-X's own, byte for byte, on every option
#                         and recording (scripts/wsprd-compare.sh)
#   ./build.sh --version TAG
#                         built as release TAG (wsprd --version); otherwise
#                         dev-<commit>
#
# Releases (the maintainer's):
#   ./build.sh --release [TAG]
#                         --check and --windows, then release both as TAG on
#                         GitHub (default vYYYY.MM.DD, with -2, -3... if that
#                         is taken), and as the release "latest" (its tag
#                         moved to TAG's commit, its files replaced): the
#                         newest's files always at
#                         .../releases/download/latest/<file>. Earlier
#                         releases are kept.
#   ./build.sh --add-built TAG DIR
#                         no build: add the wsprd another machine built (a
#                         copy of its build/dist, made with --version TAG) to
#                         release TAG, and to "latest" while that is TAG: an
#                         arm64 build from a Pi, for one
#
# The check needs WSJT-X's source at the decoder's commit (WSJTX, default
# ~/repos/wsjtx, at 567ad29ce; scripts/wsprd-oracle.sh --build-only builds its
# wsprd with WSJT-X's Release flags) and takes some minutes;
# WSPRD_SKIP_COMPARE=1 leaves it out. Publishing needs the GitHub CLI logged
# in (gh auth login), a clean working tree and HEAD pushed, so what is
# released is a commit on GitHub; Windows needs Docker.
#
# Debian or Ubuntu: cmake, g++, pkg-config and libfftw3-dev (and gfortran for
# the check).
set -euo pipefail

cd "$(dirname "$0")"
root=$(pwd)
build=$root/build

release=false
windows=false
check=false
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
            sed -n '2,41p' "$0" | sed 's/^# \{0,1\}//'
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
    local dir=$1 platform=$2 exe=""
    case "$platform" in windows-*) exe=.exe ;; esac
    echo "$dir/wsprd-$platform$exe" "$dir/SHA256SUMS-$platform"
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
    if gh release view latest >/dev/null 2>&1; then
        for a in $(gh release view latest --json assets -q '.assets[].name'); do
            gh release delete-asset latest "$a" -y
        done
        gh release edit latest --title "latest ($tag)" --notes "$notes" --latest
        gh release upload latest "$@" --clobber
    else
        gh release create latest "$@" --title "latest ($tag)" --notes "$notes" --latest
    fi
    echo "  $(gh release view latest --json url -q .url)"
}

if [ -n "$added" ]; then
    gh auth status >/dev/null 2>&1 || { echo "gh is not logged in: run gh auth login" >&2; exit 1; }
    gh release view "$tag" >/dev/null 2>&1 || { echo "there is no release $tag" >&2; exit 1; }
    arch=$(cd "$added" && ls wsprd-linux-* 2>/dev/null | head -1 | sed 's/^wsprd-linux-//')
    [ -n "$arch" ] || { echo "$added holds no wsprd-linux-<arch>" >&2; exit 1; }
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
say "building wsprd $version (Release)"
rm -rf "$build"
cmake -S "$root" -B "$build" -DCMAKE_BUILD_TYPE=Release -DWSPRD_VERSION="$version"
cmake --build "$build" -j "$(nproc)"
"$build/wsprd" --version | grep -qF "wsprd $version " || { echo "wsprd was not built as $version" >&2; exit 1; }
# FFTW is linked in (CMakeLists.txt, libfftw3-dev's libfftw3f.a): a binary
# that would need it installed is not one to hand out.
if ldd "$build/wsprd" | grep -q libfftw3; then
    echo "wsprd needs FFTW's shared library: FFTW's static one (libfftw3f.a, in libfftw3-dev) was not found" >&2
    exit 1
fi

arch=$(uname -m)
dist=$build/dist
mkdir -p "$dist"
cp "$build/wsprd" "$dist/wsprd-linux-$arch"
strip "$dist/wsprd-linux-$arch"
(cd "$dist" && sha256sum -- "wsprd-linux-$arch" > "SHA256SUMS-linux-$arch")
(. /etc/os-release && echo "$PRETTY_NAME") > "$dist/built-on"
echo "$version" > "$dist/version"
echo "built: $dist/wsprd-linux-$arch ($version)"

# --- the check: WSJT-X's own wsprd, every option and recording ---------------
if $check && [ "${WSPRD_SKIP_COMPARE:-0}" != 1 ]; then
    say "comparing with WSJT-X's wsprd"
    WSPRD_ORACLE_WORK=$build/stock "$root/scripts/wsprd-oracle.sh" --build-only
    WSPRD_COMPARE_WORK=$build/compare "$root/scripts/wsprd-compare.sh" "$build/wsprd" "$build/stock/bin/wsprd-release"
fi

$windows || exit 0

# --- Windows (Docker, mingw-w64) -----------------------------------------------
command -v docker >/dev/null && docker info >/dev/null 2>&1 ||
    { echo "--windows builds wsprd.exe in Docker: Docker is not usable here" >&2; exit 1; }
say "building wsprd.exe (mingw-w64)"
image=ubersdr-wsprd-win
wbuild=$build/windows
docker build -q -t "$image" -f docker/windows.Dockerfile docker >/dev/null
mkdir -p "$wbuild"
docker run --rm -u "$(id -u):$(id -g)" -v "$root:$root:ro" -v "$wbuild:/b" "$image" sh -c "
    cmake -S '$root' -B /b -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_TOOLCHAIN_FILE=/win/toolchain.cmake \
          -DWSPRD_VERSION='$version' &&
    cmake --build /b"
file "$wbuild/wsprd.exe" | grep -q 'PE32+ executable.*x86-64' || { echo "wsprd.exe was not built" >&2; exit 1; }
# One file: any DLL but Windows' own would have to be shipped beside it.
other=$(docker run --rm -v "$wbuild:/b" "$image" x86_64-w64-mingw32-objdump -p /b/wsprd.exe |
    awk '/DLL Name:/ {print $3}' |
    grep -viE '^(kernel32|advapi32|bcrypt|crypt32|ws2_32|msvcrt|user32|shell32|ole32)\.dll$' || true)
[ -z "$other" ] || { echo "wsprd.exe needs DLLs Windows does not have: $other" >&2; exit 1; }
docker run --rm -u "$(id -u):$(id -g)" -v "$wbuild:/b" "$image" \
    x86_64-w64-mingw32-strip -o /b/wsprd-windows-x86_64.exe /b/wsprd.exe
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
else
    echo "wine is not installed: wsprd.exe is built but not checked"
fi
cp "$wbuild/wsprd-windows-x86_64.exe" "$dist/"
(cd "$dist" && sha256sum -- wsprd-windows-x86_64.exe > SHA256SUMS-windows-x86_64)
compiler=$(docker run --rm "$image" x86_64-w64-mingw32-g++ -dumpfullversion)
echo "$(. /etc/os-release && echo "$PRETTY_NAME"), cross-compiled with mingw-w64 GCC $compiler" > "$dist/built-on-windows-x86_64"
echo
echo "built: $dist/wsprd-windows-x86_64.exe"
$release || exit 0

# --- publish -------------------------------------------------------------------
say "releasing $tag"
commit=$(git rev-parse HEAD)
notes="Built from $commit.

WSJT-X's \`wsprd\` (567ad29ce's decoder), a drop-in for it: the same options, output and files, byte for byte, in a tenth of its CPU on x86-64 and a sixth on a Raspberry Pi 5. One option more, \`-n\`: no search for a drifting signal, and each ALL_WSPR.TXT line ending in the signal's spread in Hz. \`wsprd --version\` gives this release. WSPR-2 and WSPR-15 only, as WSJT-X's \`wsprd\`: FST4W is decoded by WSJT-X's \`jt9\`, which this does not replace.

- \`wsprd-linux-x86_64\`, \`wsprd-linux-aarch64\`: Linux.
- \`wsprd-windows-x86_64.exe\`: Windows."
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
