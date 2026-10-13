# Cross-compiles wsprd for 64-bit Windows from Linux: mingw-w64, its posix
# variant (the default, win32 threads, has no std::thread or std::mutex).
# FFTW is built here, static, into /win: wsprd.exe then needs no DLL that
# Windows does not have.
FROM debian:bookworm

RUN apt-get update && apt-get install -y --no-install-recommends \
        g++-mingw-w64-x86-64-posix gcc-mingw-w64-x86-64-posix binutils-mingw-w64-x86-64 \
        cmake ninja-build make pkg-config curl ca-certificates file \
    && rm -rf /var/lib/apt/lists/* \
    && for t in g++ gcc; do \
         update-alternatives --set x86_64-w64-mingw32-$t /usr/bin/x86_64-w64-mingw32-$t-posix; \
       done

ENV HOST=x86_64-w64-mingw32 PREFIX=/win
ENV PKG_CONFIG_LIBDIR=/win/lib/pkgconfig

# FFTW, single precision, with the SIMD Ubuntu's package has (SSE2 and AVX,
# no AVX2: FFTW's AVX2 codelets round differently, and jt9's FT8 must decode
# as on Linux): the kernels are chosen at run time. Its own aligned malloc:
# mingw has no memalign for configure to find.
ARG FFTW_SHA256=56c932549852cddcfafdab3820b0200c7742675be92179e59e6215b340e26467
RUN curl -sSLf -o /tmp/fftw.tgz https://www.fftw.org/fftw-3.3.10.tar.gz \
    && echo "${FFTW_SHA256}  /tmp/fftw.tgz" | sha256sum -c - \
    && mkdir /tmp/fftw && tar xzf /tmp/fftw.tgz -C /tmp/fftw --strip-components=1 \
    && cd /tmp/fftw && ./configure -q --host=$HOST --prefix=$PREFIX --enable-float \
         --enable-sse2 --enable-avx --disable-shared --enable-static \
         --disable-fortran --disable-doc --with-our-malloc \
    && make -j"$(nproc)" >/dev/null && make install >/dev/null \
    && rm -rf /tmp/fftw*

COPY windows-toolchain.cmake /win/toolchain.cmake
