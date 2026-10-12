! WSJT-X's FST4W front end on a period, for scripts/fst4w-oracle.sh: runs
! fst4_decode's own routines (blanker, four2a, get_candidates_fst4,
! fst4_downsample, fst4_sync_search, get_fst4_bitmetrics) in fst4_decode's
! order (iwspr=1) on a 12 kHz WAV and prints each stage, for
! src/fst4w_frontend.cpp to match (tests/test_fst4w_frontend.cpp). Built
! against WSJT-X's own library (libwsjt_fort.a, module fst4_decode).
!
!   fst4w-frontend-oracle FILE.wav NTRPERIOD NDEPTH [NFQSO NTOL]
!
! Prints (floats as their bits in hex, arrays as hash32: h = (h xor w) *
! 16777619 mod 2^32 over their 32-bit words, from h = 2166136261):
!   P ntr ndepth nfft1 nfft2 nss
!   B hash(c_bigfft)                     after blanker and the r2c FFT
!   N ncand | N -                        get_candidates_fst4 (- : returned early)
!   C i f1 f2 f3 f4 f5                   candidates0(i,1:5) after the search
!   S icand hash(c2) sbest fcbest isbest  each candidate's downsample and sync
!   D ncand, then C rows                 after the duplicates are removed
!   W icand hash(c2)                     the second downsample
!   J icand ijitter is0 badsync nsync_qual hash(bitmetrics) hash(s4)
!                                        (nsync_qual - when not reached)
!
! License: GNU GPL v3 (as WSJT-X)
program fst4w_frontend_oracle
   use fst4_decode
   use, intrinsic :: iso_c_binding, only: c_int
   implicit none
   interface
      integer(c_int) function fftwf_init_threads() bind(C, name='fftwf_init_threads')
         import c_int
      end function fftwf_init_threads
      subroutine fftwf_plan_with_nthreads(n) bind(C, name='fftwf_plan_with_nthreads')
         import c_int
         integer(c_int), value :: n
      end subroutine fftwf_plan_with_nthreads
   end interface
   integer :: npatience, nthreads
   common/patience/npatience, nthreads
   integer, parameter :: NN = 160   ! fst4_params.f90's
   character(len=512) :: fname, arg
   integer*2, allocatable :: iwave(:)
   complex, allocatable :: c_bigfft(:), c2(:), cframe(:)
   real :: candidates0(200,5), bitmetrics(320,4), s4(0:3,NN)
   integer :: ntrperiod, ndepth, nfqso, ntol, nsps, nmax, ndown, nfft1, nfft2, nss, nh1
   integer :: hmod, nfa, nfb, ncand, icand, ic2, ic, isbest, isbest2, nspsec, ijitter, ioffset, is0, iend
   integer :: jittermax, nsync_qual, nsyncoh, i, nb, npct, ndropmax, inb1, inb2, inb, nexp_decode
   real :: fs, fs2, dt, dt2, tt, baud, sigbw, fa, fb, minsync, fc0, fc, fc2, sbest, fcbest, emedelay, fc_synced
   logical :: badsync
   integer :: nhdr(11), ios, nsamp

   ! jt9's FFTW setup: four2a's patience 1 (FFTW_ESTIMATE_PATIENT), one
   ! thread, FFTW's threads initialised.
   npatience = 1
   nthreads = 1
   ios = fftwf_init_threads()
   call fftwf_plan_with_nthreads(1)

   call get_command_argument(1, fname)
   call get_command_argument(2, arg); read(arg, *) ntrperiod
   call get_command_argument(3, arg); read(arg, *) ndepth
   nfqso = 1500
   ntol = 100
   if (command_argument_count() .ge. 5) then
      call get_command_argument(4, arg); read(arg, *) nfqso
      call get_command_argument(5, arg); read(arg, *) ntol
   end if
   nexp_decode = 3*256
   emedelay = 0.0
   hmod = 1
   fcbest = 0.0
   isbest = 0

   ! The period's samples (16-bit mono after a 44-byte header), zero padded
   ! to jt9's buffer, as the decoder sees them.
   allocate(iwave(30*60*12000))
   iwave = 0
   open(10, file=trim(fname), access='stream', status='old')
   read(10) nhdr
   nsamp = nhdr(11)/2
   read(10, iostat=ios) iwave(1:min(nsamp, size(iwave)))
   close(10)

   if (ntrperiod .eq. 120) then
      nsps = 8200; nmax = 120*12000; ndown = 205; nfft1 = 7200*200
   else if (ntrperiod .eq. 300) then
      nsps = 21504; nmax = 300*12000; ndown = 512; nfft1 = 7020*512
   else if (ntrperiod .eq. 900) then
      nsps = 66560; nmax = 900*12000; ndown = 1664; nfft1 = 6480*1664
   else
      nsps = 134400; nmax = 1800*12000; ndown = 3360; nfft1 = 6426*3360
   end if
   nss = nsps/ndown
   fs = 12000.0
   fs2 = fs/ndown
   nspsec = nint(fs2)
   dt = 1.0/fs
   dt2 = 1.0/fs2
   tt = nsps*dt
   baud = 1.0/tt
   sigbw = 4.0*baud
   nfft2 = nfft1/ndown
   nfft1 = nfft2*ndown
   nh1 = nfft1/2
   allocate(c_bigfft(0:nfft1/2), c2(0:nfft2-1), cframe(0:160*nss-1))
   write(*, '(a,i0,1x,i0,1x,i0,1x,i0,1x,i0)') 'P ', ntrperiod, ndepth, nfft1, nfft2, nss

   jittermax = 2
   if (ndepth .eq. 1) jittermax = 0
   ndropmax = 1
   npct = 0
   nb = nexp_decode/256 - 3
   if (nb .ge. 0) npct = nb
   inb1 = 20
   inb2 = 5
   if (nb .eq. -1) then
      inb2 = 5
   else if (nb .eq. -2) then
      inb2 = 2
   else if (nb .eq. -3) then
      inb2 = 1
   else
      inb1 = 0
   end if
   nfa = max(100, nfqso-ntol-100)
   nfb = min(4800, nfqso+ntol+100)
   fa = max(100, nint(nfqso+1.5*baud-ntol))
   fb = min(4800, nint(nfqso+1.5*baud+ntol))

   do inb = 0, inb1, inb2
      if (nb .lt. 0) npct = inb
      call blanker(iwave, nfft1, ndropmax, npct, c_bigfft)
      call four2a(c_bigfft, nfft1, 1, -1, 0)
      write(*, '(a,z8.8)') 'B ', h32c(c_bigfft, nfft1/2+1)
      nsyncoh = 8
      minsync = 1.20
      if (ntrperiod .eq. 15) minsync = 1.15
      ncand = -1
      call get_candidates_fst4(c_bigfft, nfft1, nsps, hmod, fs, fa, fb, nfa, nfb, minsync, ncand, candidates0)
      if (ncand .lt. 0) then
         write(*, '(a)') 'N -'
         cycle
      end if
      write(*, '(a,i0)') 'N ', ncand
      do i = 1, ncand
         call put_c(i)
      end do
      do icand = 1, ncand
         fc0 = candidates0(icand,1)
         call fst4_downsample(c_bigfft, nfft1, ndown, fc0, sigbw, c2)
         call fst4_sync_search(c2, nfft2, hmod, fs2, nss, ntrperiod, nsyncoh, emedelay, sbest, fcbest, isbest)
         write(*, '(a,i0,1x,z8.8,1x,z8.8,1x,z8.8,1x,i0)') 'S ', icand, h32c(c2, nfft2), &
            transfer(sbest, 0), transfer(fcbest, 0), isbest
         fc_synced = fc0 + fcbest
         candidates0(icand,3) = fc_synced
         candidates0(icand,4) = isbest
      end do
      do icand = 1, ncand
         fc = candidates0(icand,3)
         isbest = nint(candidates0(icand,4))
         do ic2 = icand+1, ncand
            fc2 = candidates0(ic2,3)
            isbest2 = nint(candidates0(ic2,4))
            if (fc2 .gt. 0.0) then
               if (abs(fc2-fc) .lt. 0.10*baud) then
                  if (abs(isbest2-isbest) .le. 2) then
                     candidates0(ic2,3) = -1
                  end if
               end if
            end if
         end do
      end do
      ic = 0
      do icand = 1, ncand
         if (candidates0(icand,3) .gt. 0) then
            ic = ic+1
            candidates0(ic,:) = candidates0(icand,:)
         end if
      end do
      ncand = ic
      write(*, '(a,i0)') 'D ', ncand
      do i = 1, ncand
         call put_c(i)
      end do
      do icand = 1, ncand
         fc_synced = candidates0(icand,3)
         isbest = nint(candidates0(icand,4))
         call fst4_downsample(c_bigfft, nfft1, ndown, fc_synced, sigbw, c2)
         write(*, '(a,i0,1x,z8.8)') 'W ', icand, h32c(c2, nfft2)
         do ijitter = 0, jittermax
            if (ijitter .eq. 0) ioffset = 0
            if (ijitter .eq. 1) ioffset = 1
            if (ijitter .eq. 2) ioffset = -1
            is0 = isbest+ioffset
            iend = is0+160*nss-1
            if (is0 .lt. 0 .or. iend .gt. (nfft2-1)) cycle
            cframe = c2(is0:iend)
            bitmetrics = 0
            nsync_qual = -999
            call get_fst4_bitmetrics(cframe, nss, bitmetrics, s4, nsync_qual, badsync)
            if (nsync_qual .eq. -999) then
               write(*, '(a,i0,1x,i0,1x,i0,1x,l1,a,z8.8,1x,z8.8)') 'J ', icand, ijitter, is0, badsync, ' - ', &
                  h32r(bitmetrics, 320*4), h32r(s4, 4*NN)
            else
               write(*, '(a,i0,1x,i0,1x,i0,1x,l1,1x,i0,1x,z8.8,1x,z8.8)') 'J ', icand, ijitter, is0, badsync, &
                  nsync_qual, h32r(bitmetrics, 320*4), h32r(s4, 4*NN)
            end if
         end do
      end do
   end do

contains

   subroutine put_c(i)
      integer, intent(in) :: i
      write(*, '(a,i0,5(1x,z8.8))') 'C ', i, (transfer(candidates0(i,ic), 0), ic=1,5)
   end subroutine put_c

   integer function h32w(h, w)
      integer(8), intent(in) :: h
      integer, intent(in) :: w
      integer(8) :: x
      x = ieor(h, iand(int(w, 8), 4294967295_8))
      h32w = int(iand(x*16777619_8, 4294967295_8) - merge(4294967296_8, 0_8, iand(x*16777619_8, 4294967295_8) .gt. 2147483647_8))
   end function h32w

   integer function h32r(a, n)
      integer, intent(in) :: n
      real, intent(in) :: a(n)
      integer(8) :: h
      integer :: k
      h = 2166136261_8
      do k = 1, n
         h = iand(int(h32w(h, transfer(a(k), 0)), 8), 4294967295_8)
      end do
      h32r = int(h - merge(4294967296_8, 0_8, h .gt. 2147483647_8))
   end function h32r

   integer function h32c(a, n)
      integer, intent(in) :: n
      complex, intent(in) :: a(n)
      integer(8) :: h
      integer :: k
      h = 2166136261_8
      do k = 1, n
         h = iand(int(h32w(h, transfer(real(a(k)), 0)), 8), 4294967295_8)
         h = iand(int(h32w(h, transfer(aimag(a(k)), 0)), 8), 4294967295_8)
      end do
      h32c = int(h - merge(4294967296_8, 0_8, h .gt. 2147483647_8))
   end function h32c

end program fst4w_frontend_oracle
