# Build environment for the reference jt9 (scripts/fst4w-oracle.sh): WSJT-X's
# own Linux CI packages (.github/images/linux-ci/install-packages.sh, normal),
# minus what only packaging and the GUI tests use, plus the distribution's
# Hamlib in place of the one CI builds (jt9 and fst4sim never link it). Ubuntu
# 24.04 like the build host, so gfortran 13.3 and FFTW 3.3.10 are the ones
# our decoder is built and checked with. The WSJT-X tree is mounted, not
# copied: the image holds no source.
FROM ubuntu:24.04
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates build-essential cmake gfortran pkg-config \
        libfftw3-dev libboost-log-dev \
        qtbase5-dev qttools5-dev qtmultimedia5-dev libqt5serialport5-dev \
        libqt5sql5-sqlite libqt5websockets5-dev \
        libusb-1.0-0-dev libudev-dev portaudio19-dev libhamlib-dev \
    && rm -rf /var/lib/apt/lists/*
