cmake_minimum_required(VERSION 3.16)
project(ubersdr_wsprd LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# x86-64-v3 (AVX2, FMA, BMI2: any x86 since Haswell or Zen): the decoder's
# searches run eight lanes wide with AVX2. "native" tunes for the build
# machine alone, an empty value builds for the compiler's baseline (the plain
# code, the same results). Not applied on ARM, where NEON is always there.
set(WSPRD_X86_ARCH "x86-64-v3" CACHE STRING "-march for x86-64 builds (x86-64-v3, native, or empty)")
if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64" AND WSPRD_X86_ARCH)
  add_compile_options(-march=${WSPRD_X86_ARCH})
endif()
# GCC notes, on ARM, that GCC 10 changed how some structs are passed; nothing
# here is called across that line.
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  add_compile_options(-Wno-psabi)
endif()
# Windows, built with mingw-w64 (docker/windows.Dockerfile): strict C++17 hides
# M_PI there, and windows.h's min and max macros would hide std::min and
# std::max.
if(WIN32)
  add_compile_definitions(_USE_MATH_DEFINES _POSIX_THREAD_SAFE_FUNCTIONS=200112L NOMINMAX WIN32_LEAN_AND_MEAN)
endif()

find_package(Threads REQUIRED)
find_package(PkgConfig REQUIRED)
pkg_check_modules(FFTW3F REQUIRED IMPORTED_TARGET fftw3f)
# FFTW linked into the binary (libfftw3-dev's libfftw3f.a), so it runs with
# nothing to install; OFF, or no libfftw3f.a, links the shared library.
option(WSPRD_STATIC_FFTW "Link FFTW's static library into wsprd" ON)
if(WSPRD_STATIC_FFTW)
  find_library(FFTW3F_STATIC_LIBRARY NAMES libfftw3f.a HINTS ${FFTW3F_LIBRARY_DIRS})
endif()
if(WSPRD_STATIC_FFTW AND FFTW3F_STATIC_LIBRARY)
  message(STATUS "FFTW linked in: ${FFTW3F_STATIC_LIBRARY}")
  set(WSPRD_FFTW ${FFTW3F_STATIC_LIBRARY} m)
  include_directories(${FFTW3F_INCLUDE_DIRS})
else()
  set(WSPRD_FFTW PkgConfig::FFTW3F)
endif()

# The decoder must decode exactly as WSJT-X's wsprd does, floats to the bit:
# as wsprd.c built with gcc -O3 and no contraction. On x86-64 that is WSJT-X's
# own CMake Release build (no -march, and plain x86-64 has no FMA to fuse
# with), so here no contraction under x86-64-v3; its vectorising changes no
# result (no reassociation without -ffast-math). Elsewhere (aarch64) WSJT-X's
# Release build fuses multiply-adds; not fused here on any platform or
# compiler, so every multiply and add is rounded as written and the NEON code
# (vmulq/vaddq, never vfmaq/vmlaq) stays the plain code's. Check an aarch64
# build's objects with objdump -d: no fmla/fmls/fmadd/fmsub/fnmadd/fnmsub.
set(WSPRD_SOURCES src/wsprd.cpp src/wsprd_osd.cpp src/wsprd_fano.cpp src/wsprd_unpack.cpp src/wsprd_files.cpp)
set_source_files_properties(${WSPRD_SOURCES} PROPERTIES COMPILE_OPTIONS "-ffp-contract=off")

# The release's tag (build.sh), else dev-<commit>: wsprd --version.
set(WSPRD_VERSION "" CACHE STRING "The version built in (a release's tag; empty: dev-<commit>)")
set(WSPRD_VERSION_BUILT "${WSPRD_VERSION}")
if(NOT WSPRD_VERSION_BUILT)
  set(WSPRD_VERSION_BUILT "dev")
  find_package(Git QUIET)
  if(GIT_FOUND)
    execute_process(COMMAND "${GIT_EXECUTABLE}" rev-parse --short HEAD WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
                    OUTPUT_VARIABLE wsprd_commit OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
    if(wsprd_commit)
      set(WSPRD_VERSION_BUILT "dev-${wsprd_commit}")
    endif()
  endif()
endif()
message(STATUS "wsprd version: ${WSPRD_VERSION_BUILT}")

add_executable(wsprd src/wsprd_main.cpp ${WSPRD_SOURCES} src/fft.cpp src/fft_wisdom.cpp)
target_include_directories(wsprd PRIVATE src)
target_link_libraries(wsprd PRIVATE ${WSPRD_FFTW} Threads::Threads)
target_compile_definitions(wsprd PRIVATE WSPRD_VERSION="${WSPRD_VERSION_BUILT}")
target_compile_options(wsprd PRIVATE -Wall -Wextra)
