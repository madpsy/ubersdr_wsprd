#!/usr/bin/env bash
# One piece of the MSK144 port against the reference, byte for byte:
#
#   scripts/msk144-piece-check.sh PIECE
#
# builds scripts/msk144-PIECE-oracle.f90 against WSJT-X's own library
# (scripts/msk144-oracle.sh fortran, in the wsjtx-oracle image) and
# scripts/msk144-PIECE-check.cpp against the port's sources with jt9's
# floating-point flags (CMakeLists.txt's FT8_FP_OPTIONS), runs the oracle,
# gives its output to the driver on stdin, and diffs the two outputs. Exit 0
# when they are identical. The oracle makes its own fixed inputs and prints
# them with its results; the driver reads the inputs back out of that output
# and prints the port's results in the same format.
#
# A later piece adds its oracle, its driver, and one line to the table below
# (the port's sources that driver links).
#
# ldpc also checks intermediate values: the oracle gets, appended, a copy of
# the reference's lib/bpdecode128_90.f90 (sed: renamed bpdecode128_90_trace,
# plus a declaration and two calls: at the top of each iteration,
# msk144_trace_print, which prints zn and tov, every value's bits, and
# checksums of them; after the tanh, msk144_trace_msg, which prints every toc
# and tanhtoc, in each iteration whose messages are read again) and
# scripts/msk144-bp-trace.f90 (those printers), and the port's optional hooks
# (bpdecode128_90_trace and bpdecode128_90_msg_trace in src/msk144_ldpc.h,
# through scripts/msk144-bp-trace.h) print the same. Every vector is traced.
# The driver also runs the line-for-line bpdecode128_90_reference beside the
# port's decoder and prints a MISMATCH line (which the oracle never prints) if
# their answers differ. The sed anchors are checked: the script stops if the
# reference moved.
#
# frontend also checks intermediate values: the oracle gets, appended, a copy
# of the reference's lib/analytic.f90 (sed: renamed analytic_trace, plus two
# calls, one before the filter's `if( beq )`, one before `c(1)=0.5*c(1)` and
# one after it, that print the filter response h, the equalizer corr, and c
# after the forward FFT, after the filter and after the half of the DC term), and the port's hook (analytic_trace in
# src/msk_frontend.h) prints the same. It also gets a copy of the reference's
# lib/msk144_freq_search.f90 (sed: renamed msk144_freq_search_trace, plus one
# call after `xcc=abs(cc)` that prints each trial's averaged frames c, its
# correlation cc and its magnitudes xcc), checked against the library's
# routine, and the port's hook (freq_search_trace in src/msk_frontend.h) prints
# the same: every trial of every search, not only the one that wins. It also
# runs msk144_freq_search once per trial frequency, for each trial's
# correlation peak. The driver links FFTW (-lfftw3f).
#
# sync also checks intermediate values: the oracle gets, appended, two copies
# of the reference's routines (sed: lib/msk144sync.f90 renamed
# msk144sync_trace, with calls that print the constants pp and cb at its first
# call and c and xcc before its peak search; lib/msk144spd.f90 renamed
# msk144spd_trace, with calls that print its constants, its spectrum windows,
# each step's peaks and ratios, the detection arrays and median, the
# candidates, each call of msk144sync, and where it calls msk144decodeframe a
# hook that runs the reference's frame decoder, on the frame itself (in place,
# as the library's msk144spd does), and prints a hash of what it was given and
# what it answered). The port's hooks (sync_trace and spd_trace in
# src/msk144_sync.h) print the same, its frame decoder is the real one
# (src/msk144_decode.cpp), and the oracle's ct, softbits and msgreceived after
# msk144spd (success or failure) are compared with the port's. Every case also
# runs the library's routines and the oracle stops if a copy's results differ
# from them (ct included).
#
# The oracle is compiled afresh on every run (scripts/msk144-oracle.sh fortran
# does not cache by name), so an edit to a oracle's source always shows.
#
# Flags: the driver is built with what CMakeLists.txt gives jt9's MSK144
# files: FT8_FP_OPTIONS (-ffp-contract=off, plus -mno-fma on GCC before 13
# for x86-64), on top of CMake's Release -O3 and the -march=x86-64-v3
# default. To verify they still match: run
#   cmake --build build -t jt9 -- VERBOSE=1 (after touching src/msk144_ldpc.cpp)
# and compare the compile line of msk144_ldpc.cpp with `flags` below.
#
#   MSK144_ORACLE_WORK=DIR   as scripts/msk144-oracle.sh (default
#                            ~/.cache/msk144-oracle); needs its `build`
#
# frame also checks intermediate values: the oracle gets, appended, a copy of
# the reference's lib/msk144decodeframe.f90 (sed: renamed msk144decodeframe_trace,
# with calls that print its constants pp and cb, cca/ccb/phase0, cfac, c after the
# phase correction, the raw softbits, the sync errors, the normalisation, llr,
# bpdecode128_90's answer and the message type bits), and its call of bpdecode128_90
# goes to the ldpc piece's traced copy (every iteration's zn, tov, toc and tanhtoc);
# every frame runs the copy and
# the library's routine and the oracle stops if their results differ. The driver's
# hooks (decode_trace in src/msk144_decode.h) print the same. The module packjt77's
# tables are shared state: the driver makes the same pack77 and unpack77 calls in the
# same order (the lines PK, RUN and QS). It also checks encode_128_90, genmsk_128_90,
# fmtmsg and atan2 (the port's atan2f_glibc against the container's glibc 2.39).
#
# msk40 (MSK40, WSJT-X's bshmsg: scripts/msk40-oracle.f90 and scripts/msk40-check.cpp, not msk144-msk40-*)
# checks bpdecode40, encode_msk40, hash, update_msk40_hasharray, msk40decodeframe, msk40sync and msk40spd. The
# oracle gets, appended, instrumented copies of the reference's lib/bpdecode40.f90, msk40decodeframe.f90,
# msk40sync.f90 and msk40spd.f90 (sed: renamed *_trace, with calls that print the iterations' zn and tov, the
# constants pp and cb, cca/phase0, cfac, c after the phase correction, the raw softbits, the sync errors, the
# normalisation, llr, bpdecode40's answer, the hash comparison, msk40sync's constants and its c and xcc before the
# peak search, msk40spd's constants, windows, steps, detection arrays and median, candidates and sync calls); every
# case runs the copy and the library's routine and the oracle stops if their results differ. The driver's hooks
# (msk40_trace in src/msk40.h) print the same. Like sync, the same oracle linked against the OpenMP library
# (msk40sync then runs min(8, OMP_GET_MAX_THREADS()) threads) must print the same at 8 and at 4 threads.
#
# Pieces: ldpc frontend sync frame msk40
set -euo pipefail

repo=$(cd "$(dirname "$0")/.." && pwd)
piece=${1:-}
case $piece in
    ldpc) sources="src/msk144_ldpc.cpp src/ft8_ldpc.cpp src/wsprd_osd.cpp" ;;
    frontend) sources="src/msk_frontend.cpp src/ft8_maths.cpp src/fft.cpp src/fft_wisdom.cpp" ;;
    sync) sources="src/msk144_sync.cpp src/msk144_decode.cpp src/msk144_ldpc.cpp src/ft8_ldpc.cpp src/wsprd_osd.cpp src/packjt77.cpp src/msk_frontend.cpp src/ft8_maths.cpp src/fft.cpp src/fft_wisdom.cpp" ;;
    frame) sources="src/msk144_decode.cpp src/msk144_ldpc.cpp src/ft8_ldpc.cpp src/wsprd_osd.cpp src/packjt77.cpp src/ft8_maths.cpp" ;;
    msk40) sources="src/msk40.cpp src/msk144_sync.cpp src/msk144_decode.cpp src/msk144_ldpc.cpp src/ft8_ldpc.cpp src/wsprd_osd.cpp src/packjt77.cpp src/msk_frontend.cpp src/ft8_maths.cpp src/fft.cpp src/fft_wisdom.cpp" ;;
    *) sed -n '/^# Pieces:/p' "$0" | sed 's/^# //' >&2; exit 2 ;;
esac
# The least number of lines an oracle's output may have (the lines the two outputs are compared on). The
# oracle decides its own case list, and the driver reads that list back out of the oracle's output, so a
# shorter list in the oracle shortens both outputs and still matches: the END line's count is no guard
# either. A floor does guard it. Each is today's line count, 16561, 119319, 230961, 82153 and 54759,
# rounded down to two significant figures, i.e. about 1 to 2 % below: an edit that drops a case group
# (the smallest is some 3 % of its piece) fails, and an added case never does. Raise the floor with the
# count when a piece gains cases.
case $piece in
    ldpc) min_lines=16000 ;;
    frontend) min_lines=119000 ;;
    sync) min_lines=228000 ;;
    frame) min_lines=82000 ;;
    msk40) min_lines=54000 ;;
esac
oracle_src=$repo/scripts/msk144-$piece-oracle.f90
driver_src=$repo/scripts/msk144-$piece-check.cpp
[[ $piece == msk40 ]] && { oracle_src=$repo/scripts/msk40-oracle.f90; driver_src=$repo/scripts/msk40-check.cpp; }
[[ -f $oracle_src && -f $driver_src ]] || { echo "msk144-piece-check: $piece has no oracle or driver" >&2; exit 2; }

work=${MSK144_ORACLE_WORK:-${XDG_CACHE_HOME:-$HOME/.cache}/msk144-oracle}/piece-$piece
mkdir -p "$work"

# jt9's flags: CMake's Release (-O3) and FT8_FP_OPTIONS: -ffp-contract=off,
# and -mno-fma on GCC before 13 for x86-64; -march=x86-64-v3 as the default
# build, where this CPU has AVX2.
cxx=${CXX:-c++}
flags=(-std=c++17 -O3 -DNDEBUG -ffp-contract=off -Wall -Wextra -I "$repo/src")
if [[ $(uname -m) == x86_64 ]]; then
    gccmajor=$("$cxx" -dumpversion | cut -d. -f1)
    [[ $("$cxx" --version | head -1) != *clang* && $gccmajor -lt 13 ]] && flags+=(-mno-fma)
    grep -q avx2 /proc/cpuinfo && flags+=(-march=x86-64-v3)
fi
srcs=(); for s in $sources; do srcs+=("$repo/$s"); done

# bpdecode128_90_trace for the ldpc and frame oracles: the reference's routine, renamed, printing (when
# the common /msktrace/ ltrace is true) its iterations' zn and tov at the top and its toc and tanhtoc after
# the tanh (iter < maxiterations: the last iteration's are never read), and the printers
bp_trace_routine() {  # OUTFILE
    local ref=${WSJTX:-$HOME/repos/wsjtx}/lib/bpdecode128_90.f90
    sed -e 's/^subroutine bpdecode128_90(/subroutine bpdecode128_90_trace(/' \
        -e 's/^end subroutine bpdecode128_90$/end subroutine bpdecode128_90_trace/' \
        -e 's/^  real Tmn$/&\n  logical ltrace\n  common \/msktrace\/ ltrace/' \
        -e 's/^! Check to see if we have a codeword (check before we do any iteration)\.$/    if( ltrace ) call msk144_trace_print(iter,zn,tov)\n&/' \
        -e 's/^    do j=1,N$/    if( ltrace .and. iter.lt.maxiterations ) call msk144_trace_msg(iter,toc,tanhtoc,nrw)\n&/' \
        "$ref" > "$1.tmp"
    [[ $(diff "$ref" "$1.tmp" | grep -c '^>') == 6 ]] ||
        { echo "msk144-piece-check: the sed instrumentation did not apply to $ref" >&2; exit 1; }
    cat "$repo/scripts/msk144-bp-trace.f90" "$1.tmp" > "$1"
}

# the oracle's source: for ldpc, with the instrumented copy of the reference's routine appended
oracle_full=$oracle_src
if [[ $piece == ldpc ]]; then
    oracle_full=$work/msk144-ldpc-oracle-traced.f90
    bp_trace_routine "$work/trace-routine.f90"
    cat "$oracle_src" "$work/trace-routine.f90" > "$oracle_full"
fi

if [[ $piece == frontend ]]; then
    ref=${WSJTX:-$HOME/repos/wsjtx}/lib/analytic.f90
    oracle_full=$work/msk144-frontend-oracle-traced.f90
    sed -e 's/^subroutine analytic(/subroutine analytic_trace(/' \
        -e 's/^end subroutine analytic$/end subroutine analytic_trace/' \
        -e 's/^  if( beq ) then$/  call msk144_dump_hcf(h,corr,c,nfft)\n&/' \
        -e 's/^  c(1)=0\.5\*c(1) .*$/  call msk144_dump_g(c,nfft)\n&\n  call msk144_dump_w(c)/' \
        "$ref" > "$work/trace-routine.f90"
    [[ $(diff "$ref" "$work/trace-routine.f90" | grep -c '^>') == 5 ]] ||
        { echo "msk144-piece-check: the sed instrumentation did not apply to $ref" >&2; exit 1; }
    ref2=${WSJTX:-$HOME/repos/wsjtx}/lib/msk144_freq_search.f90
    sed -e 's/^subroutine msk144_freq_search(/subroutine msk144_freq_search_trace(/' \
        -e 's/^end subroutine msk144_freq_search$/end subroutine msk144_freq_search_trace/' \
        -e 's/^     xcc=abs(cc)$/&\n     call msk144_dump_trial(ifr,c,cc,xcc)/' \
        "$ref2" > "$work/trace-search.f90"
    [[ $(diff "$ref2" "$work/trace-search.f90" | grep -c '^>') == 3 ]] ||
        { echo "msk144-piece-check: the sed instrumentation did not apply to $ref2" >&2; exit 1; }
    cat "$oracle_src" "$work/trace-routine.f90" "$work/trace-search.f90" > "$oracle_full"
fi

if [[ $piece == sync ]]; then
    lib=${WSJTX:-$HOME/repos/wsjtx}/lib
    oracle_full=$work/msk144-sync-oracle-traced.f90
    sed -e 's/^subroutine msk144sync(/subroutine msk144sync_trace(/' \
        -e 's/^end subroutine msk144sync$/end subroutine msk144sync_trace/' \
        -e 's/^     first=\.false\.$/&\n     call msk144_sync_consts(pp,cb)/' \
        -e 's/^! Find npeaks largest peaks$/  call msk144_sync_peak(xmax,fest,c,xcc,nfreqs)\n&/' \
        "$lib/msk144sync.f90" > "$work/trace-sync.f90"
    [[ $(diff "$lib/msk144sync.f90" "$work/trace-sync.f90" | grep -c '^>') == 4 ]] ||
        { echo "msk144-piece-check: the sed instrumentation did not apply to $lib/msk144sync.f90" >&2; exit 1; }
    sed -e 's/^subroutine msk144spd(/subroutine msk144spd_trace(/' \
        -e 's/^end subroutine msk144spd$/end subroutine msk144spd_trace/' \
        -e 's/^     first=\.false\.$/&\n     call msk144_spd_consts(rcw)/' \
        -e 's/^  i4000=nint(nfhi\/df) + 1$/&\n  call msk144_spd_win(ihlo,ihhi,illo,ilhi,i2000,i4000)/' \
        -e 's/^    detfer(istp)=ferr$/&\n    call msk144_spd_step(istp,ihpk,ilpk,deltah,deltal,ah,al,trath,tratl,ferr)/' \
        -e 's/^  xmed=detmet(indices(nstep\/4))$/&\n  call msk144_spd_med(nstep,indices,xmed,detmet,detmet2,detfer)/' \
        -e 's/^  npeaks=2$/  call msk144_spd_cands(ndet,nstart,ferrs,snrs)\n&/' \
        -e 's/call msk144sync(cdat,3,/call msk144sync_trace(cdat,3,/' \
        -e 's/^           nsyncsuccess,xmax,c)$/&\n      call msk144_spd_sync(icand,iav,nsyncsuccess,fest,xmax,npkloc)/' \
        -e 's/call msk144decodeframe(ct,softbits,msgreceived,ndecodesuccess)/call msk144_frame_hook(ct,softbits,msgreceived,ndecodesuccess)/' \
        "$lib/msk144spd.f90" > "$work/trace-spd.f90"
    [[ $(diff "$lib/msk144spd.f90" "$work/trace-spd.f90" | grep -c '^>') == 10 ]] ||
        { echo "msk144-piece-check: the sed instrumentation did not apply to $lib/msk144spd.f90" >&2; exit 1; }
    cat "$oracle_src" "$work/trace-sync.f90" "$work/trace-spd.f90" > "$oracle_full"
fi

# An oracle that fails a self-check dies with `error stop` (non-zero exit), and every oracle ends its
# output with one line, "END <piece> <count>", which the driver prints too (its own count): an output
# cut short, by a crash, a `stop` that exits 0 or a full disk, lacks it and cannot match. Both are
# required here before the two outputs are compared.
run_oracle() {  # BINARY OUTFILE [env args...]
    local bin=$1 out=$2; shift 2
    if ! "$repo/scripts/msk144-oracle.sh" exec "$@" "$bin" "$repo" > "$out"; then
        echo "msk144-piece-check: $piece: the oracle exited non-zero (output $out, $(wc -l < "$out") lines)" >&2
        tail -3 "$out" >&2 || true
        exit 1
    fi
    grep -Eq "^END $piece [0-9]+\$" <(tail -n 1 "$out") ||
        { echo "msk144-piece-check: $piece: the oracle's output ($out) does not end with 'END $piece N': truncated" >&2; exit 1; }
    (( $(wc -l < "$out") >= min_lines )) ||
        { echo "msk144-piece-check: $piece: the oracle's output ($out) has $(wc -l < "$out") lines, under the floor of $min_lines: cases are missing" >&2; exit 1; }
}

if [[ $piece == frame ]]; then
    lib=${WSJTX:-$HOME/repos/wsjtx}/lib
    oracle_full=$work/msk144-frame-oracle-traced.f90
    sed -e 's/^subroutine msk144decodeframe(/subroutine msk144decodeframe_trace(/' \
        -e 's/^end subroutine msk144decodeframe$/end subroutine msk144decodeframe_trace/' \
        -e 's/^     first=\.false\.$/&\n     call msk144_df_consts(pp,cb)/' \
        -e 's/^  phase0=atan2(imag(cca+ccb),real(cca+ccb))$/&\n  call msk144_df_phase(cca,ccb,phase0)/' \
        -e 's/^  c=c\*conjg(cfac)$/&\n  call msk144_df_corr(cfac,c)/' \
        -e 's/^! sync word hard error weight is used as a discriminator for $/  call msk144_df_soft(softbits)\n&/' \
        -e 's/^  nbadsync=nbadsync1+nbadsync2$/&\n  call msk144_df_sync(nbadsync1,nbadsync2)/' \
        -e 's/^  ssig=sqrt(s2av-sav\*sav)$/&\n  call msk144_df_norm(sav,s2av,ssig)/' \
        -e 's/^  llr=2\.0\*llr\/(sigma\*sigma)$/&\n  call msk144_df_llr(llr)/' \
        -e 's/^  call bpdecode128_90(llr,apmask,max_iterations,decoded77,cw,nharderror,niterations)$/  call bpdecode128_90_trace(llr,apmask,max_iterations,decoded77,cw,nharderror,niterations)\n  call msk144_df_bp(decoded77,cw,nharderror,niterations)/' \
        -e "s/^    read(c77(72:77),'(2b3)') n3,i3\$/&\n    call msk144_df_tag(c77,n3,i3)/" \
        "$lib/msk144decodeframe.f90" > "$work/trace-frame.f90"
    [[ $(diff "$lib/msk144decodeframe.f90" "$work/trace-frame.f90" | grep -c '^>') == 12 ]] ||
        { echo "msk144-piece-check: the sed instrumentation did not apply to $lib/msk144decodeframe.f90" >&2; exit 1; }
    grep -q '^  call bpdecode128_90_trace(' "$work/trace-frame.f90" ||
        { echo "msk144-piece-check: the frame copy does not call bpdecode128_90_trace" >&2; exit 1; }
    bp_trace_routine "$work/trace-bp.f90"
    cat "$oracle_src" "$work/trace-frame.f90" "$work/trace-bp.f90" > "$oracle_full"
fi

if [[ $piece == msk40 ]]; then
    lib=${WSJTX:-$HOME/repos/wsjtx}/lib
    oracle_full=$work/msk40-oracle-traced.f90
    # each sed's anchors are lines of the reference; the count of lines added is checked
    sed -e 's/^subroutine bpdecode40(/subroutine bpdecode40_trace(/' \
        -e 's/^end subroutine bpdecode40$/end subroutine bpdecode40_trace/' \
        -e 's/^! Check to see if we have a codeword (check before we do any iteration)\.$/  call msk40_bp_iter(iter,zn,tov)\n&/' \
        "$lib/bpdecode40.f90" > "$work/trace-bp.f90"
    [[ $(diff "$lib/bpdecode40.f90" "$work/trace-bp.f90" | grep -c '^>') == 3 ]] ||
        { echo "msk144-piece-check: the sed instrumentation did not apply to $lib/bpdecode40.f90" >&2; exit 1; }
    sed -e 's/^subroutine msk40decodeframe(/subroutine msk40decodeframe_trace(/' \
        -e 's/^end subroutine msk40decodeframe$/end subroutine msk40decodeframe_trace/' \
        -e 's/^     first=\.false\.$/&\n     call msk40_df_consts(pp,cb)/' \
        -e 's/^  phase0=atan2(imag(cca),real(cca))$/&\n  call msk40_df_phase(cca,phase0)/' \
        -e 's/^  c=c\*conjg(cfac)$/&\n  call msk40_df_corr(cfac,c)/' \
        -e 's/^! Sync word hard error weight is used to reject frames that $/  call msk40_df_soft(softbits)\n&/' \
        -e 's/^  nbadsync=nbadsync1$/&\n  call msk40_df_sync(nbadsync1)/' \
        -e 's/^  ssig=sqrt(s2av-sav\*sav)$/&\n  call msk40_df_norm(sav,s2av,ssig)/' \
        -e 's/^  llr=2\.0\*llr\/(sigma\*sigma)$/&\n  call msk40_df_llr(sigma,llr)/' \
        -e 's/^  call bpdecode40(llr,max_iterations,decoded,niterations)$/&\n  call msk40_df_bp(decoded,niterations)/' \
        -e 's/^    nrxhash=(imsg-nrxrpt)\/16$/&\n    call msk40_df_msg(nhammd,cord,imsg,nrxrpt,nrxhash,ihash)/' \
        "$lib/msk40decodeframe.f90" > "$work/trace-frame.f90"
    [[ $(diff "$lib/msk40decodeframe.f90" "$work/trace-frame.f90" | grep -c '^>') == 11 ]] ||
        { echo "msk144-piece-check: the sed instrumentation did not apply to $lib/msk40decodeframe.f90" >&2; exit 1; }
    sed -e 's/^subroutine msk40sync(/subroutine msk40sync_trace(/' \
        -e 's/^end subroutine msk40sync$/end subroutine msk40sync_trace/' \
        -e 's/^     first=\.false\.$/&\n     call msk40_sync_consts(pp,cb)/' \
        -e 's/^! Find npeaks largest peaks$/  call msk40_sync_peak(xmax,fest,c,xcc,nfreqs)\n&/' \
        "$lib/msk40sync.f90" > "$work/trace-sync.f90"
    [[ $(diff "$lib/msk40sync.f90" "$work/trace-sync.f90" | grep -c '^>') == 4 ]] ||
        { echo "msk144-piece-check: the sed instrumentation did not apply to $lib/msk40sync.f90" >&2; exit 1; }
    sed -e 's/^subroutine msk40spd(/subroutine msk40spd_trace(/' \
        -e 's/^end subroutine msk40spd$/end subroutine msk40spd_trace/' \
        -e 's/^     first=\.false\.$/&\n     call msk40_spd_consts(rcw)/' \
        -e 's/^  i4000=nint(nfhi\/df)+1$/&\n  call msk40_spd_win(ihlo,ihhi,illo,ilhi,i2000,i4000)/' \
        -e 's/^    detfer(istp)=ferr$/&\n    call msk40_spd_step(istp,ihpk,ilpk,deltah,deltal,ah,al,trath,tratl,ferr)/' \
        -e 's/^  xmed=detmet(indices(nstep\/4))$/&\n  call msk40_spd_med(nstep,indices,xmed,detmet,detmet2,detfer)/' \
        -e 's/^  npeaks=2$/  call msk40_spd_cands(ndet,nstart,ferrs,snrs)\n&/' \
        -e 's/call msk40sync(cdat,3,/call msk40sync_trace(cdat,3,/' \
        -e 's/^                     nsyncsuccess,c)$/&\n      call msk40_spd_sync(icand,iav,nsyncsuccess,fest,npkloc)/' \
        "$lib/msk40spd.f90" > "$work/trace-spd.f90"
    [[ $(diff "$lib/msk40spd.f90" "$work/trace-spd.f90" | grep -c '^>') == 9 ]] ||
        { echo "msk144-piece-check: the sed instrumentation did not apply to $lib/msk40spd.f90" >&2; exit 1; }
    cat "$oracle_src" "$work/trace-bp.f90" "$work/trace-frame.f90" "$work/trace-sync.f90" "$work/trace-spd.f90" > "$oracle_full"
fi

# the oracle
oracle=$("$repo/scripts/msk144-oracle.sh" fortran "$oracle_full" "msk144-$piece-oracle")
run_oracle "$oracle" "$work/oracle.out"

# sync: the release jt9 links the OpenMP variant of the library (msk144sync then runs
# min(4, OMP_GET_MAX_THREADS()) threads, each with a slice of the trial frequencies, merged by a strict >).
# The port is serial. The same oracle, linked against libwsjt_fort_omp.a with -fopenmp and run with four
# threads (a full team: OMP_DYNAMIC=false, OMP_THREAD_LIMIT unset), must print the same.
if [[ $piece == sync ]]; then
    oracle_omp=$("$repo/scripts/msk144-oracle.sh" fortran "$oracle_full" "msk144-$piece-oracle-omp" omp)
    run_oracle "$oracle_omp" "$work/oracle-omp.out" env -u OMP_THREAD_LIMIT OMP_NUM_THREADS=4 OMP_DYNAMIC=false
    cmp -s "$work/oracle.out" "$work/oracle-omp.out" ||
        { echo "msk144-piece-check: sync: the OpenMP library at 4 threads DIFFERS: diff $work/oracle.out $work/oracle-omp.out" >&2
          diff "$work/oracle.out" "$work/oracle-omp.out" | head -10 >&2 || true; exit 1; }
fi

# msk40: the same, linked against the OpenMP library (msk40sync's team of min(8, max threads)) at 8 and at 4 threads
if [[ $piece == msk40 ]]; then
    oracle_omp=$("$repo/scripts/msk144-oracle.sh" fortran "$oracle_full" "msk144-$piece-oracle-omp" omp)
    for nt in 8 4; do
        run_oracle "$oracle_omp" "$work/oracle-omp$nt.out" env -u OMP_THREAD_LIMIT OMP_NUM_THREADS=$nt OMP_DYNAMIC=false
        cmp -s "$work/oracle.out" "$work/oracle-omp$nt.out" ||
            { echo "msk144-piece-check: msk40: the OpenMP library at $nt threads DIFFERS: diff $work/oracle.out $work/oracle-omp$nt.out" >&2
              diff "$work/oracle.out" "$work/oracle-omp$nt.out" | head -10 >&2 || true; exit 1; }
    done
fi

# the driver (no output file if the port is missing: the link fails)
rm -f "$work/driver"
"$cxx" "${flags[@]}" "$driver_src" "${srcs[@]}" -lfftw3f -lm -lpthread -o "$work/driver"
"$work/driver" < "$work/oracle.out" > "$work/driver.out" 2> "$work/driver.err" ||
    { echo "msk144-piece-check: $piece: the driver exited non-zero" >&2; cat "$work/driver.err" >&2; tail -3 "$work/driver.out" >&2 || true; exit 1; }
cat "$work/driver.err" >&2
# ldpc: the driver must say the fast tables passed (else bpdecode128_90 ran its reference fallback
# and the identical output proves nothing about the fast path)
if [[ $piece == ldpc ]] && ! grep -q 'fast tables passed' "$work/driver.err"; then
    echo "msk144-piece-check: ldpc: the driver did not report 'fast tables passed'" >&2; exit 1
fi

if cmp -s "$work/oracle.out" "$work/driver.out"; then
    if [[ $piece == frontend ]]; then
        echo "msk144-piece-check: $piece: identical ($(wc -l < "$work/oracle.out") lines: $(grep -c '^A' "$work/oracle.out") analytic, $(grep -c '^T' "$work/oracle.out") tweak1, $(grep -c '^S' "$work/oracle.out") freq_search, $(grep -c '^R' "$work/oracle.out") traced trials, $(grep -c '^Q' "$work/oracle.out") trial peaks)"
        exit 0
    fi
    if [[ $piece == frame ]]; then
        echo "msk144-piece-check: frame: identical ($(wc -l < "$work/oracle.out") lines: $(grep -c '^RUN' "$work/oracle.out") frames, $(grep -c '^R .* [1-9][0-9]* |' "$work/oracle.out") decoded, $(grep -c '^QS' "$work/oracle.out") signal-quality runs, $(grep -c '^KD' "$work/oracle.out") LDPC runs ($(grep -c '^I ' "$work/oracle.out") traced iterations), $(grep -c '^EN' "$work/oracle.out") encodes, $(grep -c '^GM' "$work/oracle.out") genmsk, $(grep -c '^FM' "$work/oracle.out") fmtmsg, $(grep -c '^AV' "$work/oracle.out") atan2 values and $(grep -c '^AT' "$work/oracle.out") atan2 sweeps)"
        exit 0
    fi
    if [[ $piece == msk40 ]]; then
        echo "msk144-piece-check: msk40: identical, and the OpenMP library at 8 and at 4 threads the same ($(wc -l < "$work/oracle.out") lines: $(grep -c '^BI' "$work/oracle.out") bpdecode40 runs ($(grep -c '^BD [0-9]* [0-9]' "$work/oracle.out") decoded), $(grep -c '^RUN' "$work/oracle.out") msk40decodeframe runs ($(grep -c '^R [0-9]* [1-3] ' "$work/oracle.out") decoded), $(grep -c '^YP' "$work/oracle.out") msk40sync runs ($(grep -c '^YA [0-9]* 1 ' "$work/oracle.out") syncs), $(grep -c '^ZP' "$work/oracle.out") msk40spd runs ($(grep -c '^X [0-9]* [1-3] ' "$work/oracle.out") decoded), $(grep -c '^HA' "$work/oracle.out") hashes, $(grep -c '^NH' "$work/oracle.out") hash arrays)"
        exit 0
    fi
    if [[ $piece == sync ]]; then
        echo "msk144-piece-check: sync: identical, and the OpenMP library at 4 threads the same ($(wc -l < "$work/oracle.out") lines: $(grep -c '^SYN' "$work/oracle.out") msk144sync and $(grep -c '^SPD' "$work/oracle.out") msk144spd runs; $(awk '$1 == "X" && $3 == 1 && $4 !~ /^\|/' "$work/oracle.out" | wc -l) msk144sync successes; $(grep -c '^F ' "$work/oracle.out") frame decodes by the port's own decoder ($(grep -c '^F .* [1-9][0-9]* |' "$work/oracle.out") decoded); $(grep -c '^T' "$work/oracle.out") spd steps, $(grep -c '^[LJKW] ' "$work/oracle.out") piece cases)"
        exit 0
    fi
    echo "msk144-piece-check: $piece: identical ($(wc -l < "$work/oracle.out") lines, $(grep -c '^R' "$work/oracle.out") decodes run, $(grep -c '^I' "$work/oracle.out") traced iterations ($(grep -c '^BO' "$work/oracle.out") with their messages), $(grep -c '^P' "$work/oracle.out") platanh and $(grep -c '^T' "$work/oracle.out") tanh values, $(awk '$1 == "R" && $4 >= 0' "$work/oracle.out" | wc -l) decoded, $(grep -c '^C' "$work/oracle.out") CRC checks)"
else
    echo "msk144-piece-check: $piece: DIFFERS: diff $work/oracle.out $work/driver.out" >&2
    diff "$work/oracle.out" "$work/driver.out" | head -10 >&2 || true
    exit 1
fi
