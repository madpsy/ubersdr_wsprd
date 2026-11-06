! WSJT-X's MSK144 front end for scripts/msk144-piece-check.sh frontend:
! generates its own fixed inputs, prints them and what the reference does with
! them, for scripts/msk144-frontend-check.cpp (src/msk_frontend.cpp) to match
! line for line. Built against WSJT-X's own library (libwsjt_fort.a): analytic,
! tweak1, msk144_freq_search, four2a (npatience 1, as jt9's default: the plans
! are FFTW_ESTIMATE_PATIENT). The piece check appends analytic_trace to this
! file: lib/analytic.f90 with two calls added (see scripts/msk144-piece-check.sh)
! that print the filter response and the FFT outputs. The oracle runs both
! routines on every case and stops if their outputs differ.
!
! Every float is its bits, 8 hex digits (a real*8, 16); arrays 16 values a line,
! "TAG k first-index values", complex values as re, im pairs.
!
!   A k npts nfft beq pc(1..5)   an analytic case (pc: real*8 bits)
!   D k i v...                   its input d(1:npts), as mskrtd makes it (the
!                                int16 block over its rms)
!   H k i v...                   h(1:nfft/2+1), the filter response   [analytic's]
!   K k i v...                   corr(1:nfft/2+1), the equalizer      [analytic's]
!   F k i v...                   c(1:nfft) after the forward FFT      [analytic's]
!   G k i v...                   c(1:nfft/2+1) after the filter       [analytic's]
!   W k i v...                   c(1:2) after the half of the DC term [analytic's]
!   O k i v...                   c(1:nfft), the analytic signal
!   T k src jz f0                tweak1 of case src's O, jz points, f0 (float bits)
!   U k i v...                   cb(1:jz)
!   B i v...                     the sync word cb(1:42), as msk144sync makes it
!   S k src nframes navmask(8) fc if1 if2 delf
!                                msk144_freq_search on src's O(1:864*nframes)
!   X k xmax bestf               its xmax and bestf
!   Y k i v...                   cs(1:864)
!   Z k i v...                   xccs(0:863)
!   V k i v...                   cdat2(1:864*nframes), the last trial's tweak
!   Q k ifr xmax bestf           the same call with if1 = if2 = ifr: each
!                                trial's correlation peak, for ifr = if1..if2
!
! Inputs: WSJT-X's two samples and two of our sim files (tests/data/msk144's
! snr-p10 and freq-p200-sp04): in each, the 7168-sample block (steps of 864)
! of most energy and the block of least.
!
! License: GNU GPL v3 (as WSJT-X)
program msk144_frontend_oracle
   implicit none
   integer, parameter :: NZ = 7168, NSPM = 864, NF = 8192
   character(len=512) :: repo, path(4)
   integer :: nf4, nfft, k, i, j, ifile, iblk, npts, nsrc(64), nsamp, ios, istart, ibest, iworst, ntol, nfreqs, if1, if2
   integer*2, allocatable :: wav(:)
   real :: d(NF), rms, e, ebest, eworst, fc, delf
   real*8 :: pc(5)
   logical*1 :: beq
   complex :: c(NF), c2(NF), cb(42), cs(NSPM), cdat2(NSPM*8), cs1(NSPM), cd1(NSPM*8)
   complex, allocatable :: osave(:,:)
   real :: foff(3), xmax, bestf, f0
   integer :: navmask(8), nframes, ifr, nint_
   integer :: npatience, nthreads
   common /patience/ npatience, nthreads
   common /mskcase/ kcase
   integer :: kcase
   real :: pp(12), pi, cbi(42), cbq(42)
   integer :: s8(8)
   data s8/0,1,1,1,0,0,1,0/
   data foff/-1500.0, -1507.3, 1488.9/

   npatience = 1
   nthreads = 1
   call get_command_argument(1, repo)
   path(1) = trim(repo)//'/recordings/msk144-samples/181211_120500.wav'
   path(2) = trim(repo)//'/recordings/msk144-samples/181211_120800.wav'
   path(3) = trim(repo)//'/recordings/msk144-sim/snr-p10-261106_120010.wav'
   path(4) = trim(repo)//'/recordings/msk144-sim/freq-p200-sp04-261106_120019.wav'
   allocate(osave(NF, 64))
   osave = 0.0
   nsrc = 0

   ! the sync word, as msk144sync defines it
   pi = 4.0*atan(1.0)
   do i = 1, 12
      pp(i) = sin((i-1)*pi/12.0)
   enddo
   s8 = 2*s8 - 1
   cbq(1:6) = pp(7:12)*s8(1)
   cbq(7:18) = pp*s8(3)
   cbq(19:30) = pp*s8(5)
   cbq(31:42) = pp*s8(7)
   cbi(1:12) = pp*s8(2)
   cbi(13:24) = pp*s8(4)
   cbi(25:36) = pp*s8(6)
   cbi(37:42) = pp(1:6)*s8(8)
   cb = cmplx(cbi, cbq)
   call prc('B', 0, cb, 42)

   ! ---- analytic: cases 1-8 are (file, hottest block), (file, quietest block) for the four files
   k = 0
   do ifile = 1, 4
      open(10, file=trim(path(ifile)), access='stream', form='unformatted', status='old', iostat=ios)
      if (ios /= 0) error stop 'oracle: cannot open a WAV'
      inquire(10, size=nsamp)
      nsamp = (nsamp - 44)/2
      allocate(wav(nsamp))
      read(10, pos=45) wav
      close(10)
      ebest = -1.0; eworst = huge(1.0)
      ibest = 0; iworst = 0
      istart = 0
      do while (istart + NZ <= nsamp)
         e = 0.0
         do i = 1, NZ
            e = e + real(wav(istart+i))**2
         enddo
         if (e > ebest) then
            ebest = e; ibest = istart
         endif
         if (e < eworst) then
            eworst = e; iworst = istart
         endif
         istart = istart + NSPM*2
      enddo
      do iblk = 1, 2
         istart = ibest
         if (iblk == 2) istart = iworst
         d = 0.0
         d(1:NZ) = wav(istart+1:istart+NZ)
         rms = sqrt(sum(d(1:NZ)*d(1:NZ))/NZ)
         if (rms < 1.0) rms = 1.0
         d(1:NZ) = (1.0/rms)*d(1:NZ)
         d(NZ+1:NF) = 0.0
         k = k + 1
         pc = 0.0d0
         beq = .false.
         call run_analytic(k, d, NZ, NF, pc, beq)
      enddo
      deallocate(wav)
   enddo
   ! 9: the model-waveform shape msk144signalquality uses: 864 points, nfft 1024 (a new nfft: h made again)
   k = 9
   pc = 0.0d0; beq = .false.
   d(1:1024) = 0.0
   d(1:NSPM) = real(osave(1:NSPM, 7))      ! any real block: the real part of case 7's analytic signal
   call run_analytic(k, d, NSPM, 1024, pc, beq)
   ! 10: back to 8192, the equalizer on, phase coefficients set (corr made, used)
   k = 10
   pc = (/0.1d0, -0.2d0, 0.3d0, -0.15d0, 0.07d0/)
   beq = .true.
   d(1:NF) = 0.0
   d(1:NZ) = real(osave(1:NZ, 5))
   call run_analytic(k, d, NZ, NF, pc, beq)
   ! 11: the same pc, the equalizer off (corr kept, not used)
   k = 11
   beq = .false.
   call run_analytic(k, d, NZ, NF, pc, beq)
   ! 12: 4096 points of a 3000-point block
   k = 12
   pc = 0.0d0
   d(1:NF) = 0.0
   d(1:3000) = real(osave(1:3000, 3))
   call run_analytic(k, d, 3000, 4096, pc, beq)

   ! ---- tweak1: three offsets on each hottest block (cases 1, 3, 5, 7)
   k = 12
   do i = 1, 7, 2
      do j = 1, 3
         k = k + 1
         f0 = foff(j)
         call run_tweak(k, i, NZ, f0)
      enddo
   enddo

   ! ---- msk144_freq_search: windows as msk144sync makes them (nfreqs, nstep), and the odd ones
   ! (k src nframes navmask fc if1 if2 delf)
   fc = 1500.0
   ! the four averaging patterns of mskrtd; deltaf = 10/navg; ntol 15, 25, 40, 20 (thread slices after)
   navmask = (/1,1,1,1,0,0,0,0/); call win(1, 8, navmask, 15.0, 1)
   navmask = (/1,1,1,1,0,0,0,0/); call win(1, 8, navmask, 15.0, 3)
   navmask = (/0,0,1,1,1,1,0,0/); call win(3, 8, navmask, 25.0, 1)
   navmask = (/1,1,1,1,1,0,0,0/); call win(5, 8, navmask, 40.0, 1)
   navmask = (/1,1,1,1,1,1,1,0/); call win(7, 8, navmask, 15.0, 1)
   navmask = (/1,1,1,1,1,1,1,0/); call win(7, 8, navmask, 15.0, 4)
   navmask = (/1,1,1,0,0,0,0,0/); call win(5, 3, navmask, 20.0, 1)   ! 3 frames
   navmask = (/1,1,1,1,0,0,0,0/); call win(2, 8, navmask, 15.0, 1)   ! the quietest blocks
   navmask = (/1,1,1,1,1,1,1,0/); call win(6, 8, navmask, 15.0, 1)
   navmask = (/0,0,0,0,0,0,0,0/); call one(1, 8, navmask, 0, 4, 2.5)   ! navg 0: nothing done
   navmask = (/1,1,1,1,0,0,0,0/); call one(1, 8, navmask, 5, -5, 2.5)  ! if1 > if2: nothing done
   navmask = (/1,1,1,1,0,0,0,0/); call one(3, 8, navmask, 0, 0, 2.5)   ! one trial
   navmask = (/1,1,1,1,0,0,0,0/); call one(3, 8, navmask, -3, 3, 0.0)  ! delf 0: the same trial seven times
   navmask = (/1,1,1,1,0,0,0,0/); call one(3, 8, navmask, -3, 3, 100.0) ! far from the signal

   ! ---- analytic again, after the searches (k goes on from the last search)
   ! the stale-corr path: a pc change seen at nfft 1024 with the equalizer off makes corr(1:513) and
   ! sets pclast; at nfft 8192 with the equalizer on and the same pc, h is made again (nfft changed)
   ! but corr is not: its entries from 514 on are the earlier calls' (msk144signalquality and mskrtd
   ! call analytic this way)
   k = k + 1
   pc = (/0.05d0, -0.1d0, 0.8d0, -0.3d0, 0.2d0/)
   beq = .false.
   d(1:NF) = 0.0
   d(1:NSPM) = real(osave(1:NSPM, 7))
   call run_analytic(k, d, NSPM, 1024, pc, beq)
   k = k + 1
   beq = .true.
   d(1:NF) = 0.0
   d(1:NZ) = real(osave(1:NZ, 5))
   call run_analytic(k, d, NZ, NF, pc, beq)
   ! the DC term: with the equalizer, h(1) = 0 times corr(1) = (cos pd, sin pd) with pd = 2.475
   ! (negative real part, positive imaginary part) times a DC term of negative real part gives
   ! c(1) = (+0, -0) before the half: where a real 0.5*c(1) and the complex product with (0.5, 0)
   ! differ in the sign of a zero
   k = k + 1
   pc = (/0.0d0, 0.0d0, 1.1d0, 0.0d0, 0.0d0/)
   d(1:NF) = 0.0
   d(1:NZ) = real(osave(1:NZ, 5)) - 2.0
   call run_analytic(k, d, NZ, NF, pc, beq)
   ! the end marker: a truncated output (a failed self-check, a crash) cannot match without it
   write(*, '(a,1x,i0)') 'END frontend', k

contains

   ! k's window under msk144sync's split: ntol/delf; nthreads slices, slice id (1..nthreads)
   subroutine win(src, nfr, mask, ntol_, nthr)
      integer :: src, nfr, mask(8), nthr, ifrq1, ifrq2, nstep, nfreqs_, navg
      real :: ntol_, delf_
      integer :: idt
      navg = sum(mask(1:8))
      delf_ = 10.0/real(navg)
      nfreqs_ = 2*nint(ntol_/delf_) + 1
      nstep = nfreqs_/nthr
      idt = min(2, nthr)               ! the second slice when there are several
      ifrq1 = -nint(ntol_/delf_) + (idt-1)*nstep
      ifrq2 = ifrq1 + nstep - 1
      if (idt == nthr) ifrq2 = nint(ntol_/delf_)
      call run_search_args(src, nfr, mask, ifrq1, ifrq2, delf_)
   end subroutine win

   subroutine one(src, nfr, mask, ifrq1, ifrq2, delf_)
      integer :: src, nfr, mask(8), ifrq1, ifrq2
      real :: delf_
      call run_search_args(src, nfr, mask, ifrq1, ifrq2, delf_)
   end subroutine one

   subroutine run_search_args(src, nfr, mask, ifrq1, ifrq2, delf_)
      integer :: src, nfr, mask(8), ifrq1, ifrq2, n
      real :: delf_
      real :: ccs(0:NSPM-1), x1, b1
      complex :: cdat(NSPM*8), cdat2l(NSPM*8), csl(NSPM), cs1l(NSPM), cd1l(NSPM*8)
      integer :: kk, ifrq
      k = k + 1
      n = NSPM*nfr
      cdat(1:n) = osave(1:n, src)
      write(*, '(a,1x,i0,1x,i0,1x,i0,8(1x,i0),1x,z8.8,2(1x,i0),1x,z8.8)') 'S', k, src, nfr, mask(1:8), &
         transfer(fc, 0), ifrq1, ifrq2, transfer(delf_, 0)
      cdat2l = 0.0
      call msk144_freq_search(cdat, fc, ifrq1, ifrq2, delf_, nfr, mask, cb, cdat2l, xmax, bestf, csl, ccs)
      write(*, '(a,1x,i0,2(1x,z8.8))') 'X', k, transfer(xmax, 0), transfer(bestf, 0)
      call prc('Y', k, csl, NSPM)
      call prr('Z', k, ccs, NSPM)
      call prc('V', k, cdat2l, n)
      do ifrq = ifrq1, ifrq2
         cs1l = 0.0
         call msk144_freq_search(cdat, fc, ifrq, ifrq, delf_, nfr, mask, cb, cd1l, x1, b1, cs1l, ccs)
         write(*, '(a,1x,i0,1x,i0,2(1x,z8.8))') 'Q', k, ifrq, transfer(x1, 0), transfer(b1, 0)
      enddo
   end subroutine run_search_args

   subroutine run_analytic(kk, dd, np, nfft_, pcc, beqq)
      integer :: kk, np, nfft_
      real :: dd(*)
      real*8 :: pcc(5), pctr(5)
      logical*1 :: beqq
      complex, allocatable :: ca(:), cb2(:)
      allocate(ca(NF), cb2(NF))
      kcase = kk
      write(*, '(a,1x,i0,1x,i0,1x,i0,1x,i0,5(1x,z16.16))') 'A', kk, np, nfft_, merge(1, 0, beqq), &
         (transfer(pcc(i), 0_8), i=1,5)
      call prr('D', kk, dd, np)
      pctr = pcc
      call analytic_trace(dd, np, nfft_, ca, pctr, beqq)    ! prints H, K, F, G
      call analytic(dd, np, nfft_, cb2, pcc, beqq)
      if (any(transfer(ca(1:nfft_), (/0/)) /= transfer(cb2(1:nfft_), (/0/)))) &
         error stop 'oracle: analytic_trace differs from analytic'
      call prc('O', kk, ca, nfft_)
      osave(1:nfft_, kk) = ca(1:nfft_)
      osave(nfft_+1:NF, kk) = 0.0
   end subroutine run_analytic

   subroutine run_tweak(kk, src, jz, ff0)
      integer :: kk, src, jz
      real :: ff0
      complex :: ca(NZ), cbb(NZ)
      ca = osave(1:NZ, src)
      write(*, '(a,1x,i0,1x,i0,1x,i0,1x,z8.8)') 'T', kk, src, jz, transfer(ff0, 0)
      call tweak1(ca, jz, ff0, cbb)
      call prc('U', kk, cbb, jz)
   end subroutine run_tweak

end program msk144_frontend_oracle

! ---- printing: real and complex arrays, 16 values a line
subroutine prr(tag, k, x, n)
   implicit none
   character(len=*) :: tag
   integer :: k, n, i, i2, j
   real :: x(*)
   do i = 1, n, 16
      i2 = min(i + 15, n)
      write(*, '(a,1x,i0,1x,i0,16(1x,z8.8))') tag, k, i, (transfer(x(j), 0), j = i, i2)
   enddo
end subroutine prr

subroutine prc(tag, k, x, n)
   implicit none
   character(len=*) :: tag
   integer :: k, n, i, i2, j
   complex :: x(*)
   real :: v(2*n)
   do j = 1, n
      v(2*j-1) = real(x(j))
      v(2*j) = aimag(x(j))
   enddo
   do i = 1, 2*n, 16
      i2 = min(i + 15, 2*n)
      write(*, '(a,1x,i0,1x,i0,16(1x,z8.8))') tag, k, i, (transfer(v(j), 0), j = i, i2)
   enddo
end subroutine prc

! the trace calls analytic_trace makes
subroutine msk144_dump_hcf(h, corr, c, nfft)
   implicit none
   integer :: nfft, nh, kcase
   real :: h(*)
   complex :: corr(*), c(*)
   common /mskcase/ kcase
   nh = nfft/2
   call prr('H', kcase, h, nh+1)
   call prc('K', kcase, corr, nh+1)
   call prc('F', kcase, c, nfft)
end subroutine msk144_dump_hcf

subroutine msk144_dump_w(c)
   implicit none
   integer :: kcase
   complex :: c(*)
   common /mskcase/ kcase
   call prc('W', kcase, c, 2)
end subroutine msk144_dump_w

subroutine msk144_dump_g(c, nfft)
   implicit none
   integer :: nfft, kcase
   complex :: c(*)
   common /mskcase/ kcase
   call prc('G', kcase, c, nfft/2+1)
end subroutine msk144_dump_g
