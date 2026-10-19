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

# jt9: WSJT-X's jt9 run on WAV files, a mode at a time (FST4W, FT8, FT4 and FT2), its
# command line, output and files as WSJT-X's (src/jt9_cli.h). Both must
# decode exactly as WSJT-X's jt9 does, floats to the bit: no contraction, as
# above. FST4W calls glibc's vector maths (libmvec, through libm) where
# WSJT-X's -O3 build does (on Windows, without glibc, the plain functions:
# the same printed results, scripts/jt9-wine-check.sh). FT8, FT4 and FT2 call no
# maths library function but sqrtf, fabsf, fmodf and lroundf: glibc's own
# (and x86-64 libmvec's) are carried in src/ft8_maths.cpp, so they decode
# alike everywhere. Linked as wsprd is (FFTW built in).
set(FST4W_SOURCES src/fst4w_unpack.cpp src/fst4w_ldpc.cpp src/fst4w_frontend.cpp src/fst4w.cpp src/fst4w_jt9.cpp)
set_source_files_properties(${FST4W_SOURCES} PROPERTIES COMPILE_OPTIONS "-ffp-contract=off")
# GCC before 13 (mingw-w64's 12.2 in docker/windows.Dockerfile, for one)
# fuses complex multiplies into vfmaddsub under -march=x86-64-v3 despite
# -ffp-contract=off: there FT8's, FT4's and FT2's files get no FMA instructions at all
# (-mno-fma). ft8_maths.cpp keeps FMA for its explicit std::fma, the only
# fused operations in it.
set(FT8_FP_OPTIONS "-ffp-contract=off")
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 13
   AND CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64" AND WSPRD_X86_ARCH)
  set(FT8_FP_OPTIONS "-ffp-contract=off;-mno-fma")
endif()
set(FT8_SOURCES src/packjt77.cpp src/ft8_ldpc.cpp src/ft8_frontend.cpp src/ft8.cpp)
# --fast-subtract's smoothing (not WSJT-X's to the bit, so none of its floating-point options)
set(FT_FAST_SOURCES src/ft_fast_subtract.cpp)
set_source_files_properties(src/ft8_maths.cpp PROPERTIES COMPILE_OPTIONS "-ffp-contract=off")
set_source_files_properties(${FT8_SOURCES} PROPERTIES COMPILE_OPTIONS "${FT8_FP_OPTIONS}")
set(FT4_SOURCES src/ft4_frontend.cpp src/ft4_candidates.cpp src/ft4_sync.cpp src/ft4_metrics.cpp src/ft4_wave.cpp
  src/ft4.cpp src/ft4_jt9.cpp)
set_source_files_properties(${FT4_SOURCES} PROPERTIES COMPILE_OPTIONS "${FT8_FP_OPTIONS}")
# FT2: FT4's decoder at FT2's sizes (its own copies of what they change), as
# WSJT-X improved 3.2.0 runs FT2, built as FT4's files are.
set(FT2_SOURCES src/ft2_sync.cpp src/ft2.cpp src/ft2_jt9.cpp)
set_source_files_properties(${FT2_SOURCES} PROPERTIES COMPILE_OPTIONS "${FT8_FP_OPTIONS}")
add_executable(jt9 src/jt9_main.cpp src/jt9_cli.cpp src/jt9_shm.cpp ${FST4W_SOURCES} ${FT8_SOURCES} ${FT_FAST_SOURCES} ${FT4_SOURCES} ${FT2_SOURCES}
  src/ft8_maths.cpp src/ft8_jt9.cpp
  src/wsprd_osd.cpp src/fft.cpp src/fft_wisdom.cpp)
target_include_directories(jt9 PRIVATE src)
target_link_libraries(jt9 PRIVATE ${WSPRD_FFTW} m Threads::Threads)
target_compile_definitions(jt9 PRIVATE JT9_VERSION="${WSPRD_VERSION_BUILT}")
target_compile_options(jt9 PRIVATE -Wall -Wextra)
