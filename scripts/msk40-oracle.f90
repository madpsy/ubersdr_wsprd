! WSJT-X's MSK40 routines for scripts/msk144-piece-check.sh msk40: makes its
! own fixed inputs, prints them and what the reference does with them, for
! scripts/msk40-check.cpp (src/msk40.cpp) to match line for line. Built
! against WSJT-X's own library (libwsjt_fort.a): msk40decodeframe, msk40sync,
! msk40spd, bpdecode40, encode_msk40, update_msk40_hasharray, hash, fmtmsg,
! packjt77 (module). The piece check appends instrumented copies of
! lib/bpdecode40.f90, msk40decodeframe.f90, msk40sync.f90 and msk40spd.f90
! (sed, see scripts/msk144-piece-check.sh): bpdecode40_trace,
! msk40decodeframe_trace, msk40sync_trace and msk40spd_trace, with calls
! added that print intermediate values. Every case runs the copy and the
! library's routine and the oracle stops (error stop) if their results
! differ. The copy of msk40spd calls the copy of msk40sync and the library's
! msk40decodeframe, whose frames this oracle does not print (the direct
! frame cases check that routine); the driver turns its frame traces off
! there.
!
! Every float is its bits, 8 hex digits; arrays 16 values a line, "TAG k
! first-index values", complex values as re, im pairs.
!
! Lines, in the order a case prints them
!   EN k |m16| |cw32|            encode_msk40 (16 explicit messages)
!   EA n h1 h2                   encode_msk40 on all 65536 messages: hashes of the codewords
!   HA k len |text| hash         hash(text, len): the 15 bits
!   PK k |text| i3 n3 c77        pack77_legacy_truncating_fallback of a text; its unpack77 follows
!   UN k |msg|                   unpack77 of that c77 (saves its calls: packjt77's recent_calls)
!   RC k |call|x10               recent_calls after it
!   NH k v...                    update_msk40_hasharray: nhasharray(1:10,1:10), i fastest
!   BI k llr...                  bpdecode40's input llr(1:32)
!   BZ k*100+iter i v...         the trace: zn(1:32) then tov(1:3,1:32)
!   BD k niter |decoded16|       its answer ('-' for decoded when niter < 0)
!   FP k |mycall| |hiscall| bswl xsnr   msk40decodeframe's other arguments (calls character*12)
!   FC k i v...                  the frame c(1:240) it is given
!   RUN k                        msk40decodeframe runs (the copy, then the library)
!   KP/KB k i v...               first call only: pp(1:12), cb(1:42)
!   KH k ccar ccai phase0        cca, phase0
!   KF k cfacr cfaci / KC k ..   cfac; c(1:240) after c=c*conjg(cfac)
!   KS k i v...                  softbits(1:40) after the matched filter
!   KY k nbadsync1               the sync word errors
!   KN k sav s2av ssig           the normalisation
!   KL k sigma i v...            llr(1:32) after 2*llr/sigma**2 (sigma on the line before: KG)
!   KD k niter |decoded16|       bpdecode40's answer, as the frame decoder sees it
!   KM k nhammd cord imsg nrxrpt nrxhash ihash   what the decoder compares
!   R k nsuccess |msgreceived|   the answer; O k i v...: c as the routine left it
!   SY k |...| parameters        msk40sync's cases: SC k i v... its input cdat(1:720) (3 frames)
!   YP k fc navmask(3)           the arguments; YR k: it runs; YA k nsuccess fest npklocs(1:2) / YC k i v...: c(1:240)
!   YK k first-index values      first call only: pp, cb; YX k xmax fest nfreqs, YW k i v...: c, xcc before the peaks
!   ZB k i v...                  msk40spd's input block cbig(1:6912)
!   ZP k |mycall| |hiscall| bswl ntol fc   the arguments; ZR k: it runs
!   N/T/D/H/O/M/E/J              msk40spd's trace (windows, steps, detmet, detmet2, detfer, median, candidates, syncs)
!   X k nsuccess |msgreceived| fret tret navg    its answer
! The end: END msk40 n (the cases run)
!
! four2a runs at npatience 1 (FFTW_ESTIMATE_PATIENT), as jt9's default: the plans of the
! spd's squared-signal spectrum (240 points) are the reference's.
!
! Inputs: a 32-bit LCG (x = x*1664525 + 1013904223 mod 2**32, seed 20261106).
! Frames are the transmitted waveform (the half-sine pulses: c = (xi, xq) of
! the 8 sync bits and 32 codeword bits), rotated by a phase and a frequency
! error, plus Gaussian noise (a sum of 12 uniforms less 6, times sigma). The
! spd blocks are that waveform repeated, as an analytic signal at a carrier,
! in a window of the block, plus noise.
!
! Blind spot: NaNs compare equal whatever their sign and payload, and print
! as the default NaN (cbits below): the instrumented copies and the library
! are different objects, and x86 returns the first operand's NaN. A port that
! gets a NaN's sign or payload wrong passes here.
!
! License: GNU GPL v3 (as WSJT-X)
program msk40_oracle
   use packjt77
   implicit none
   integer, parameter :: NSPM = 240, NB = 6912
   integer*8 :: state
   integer :: kc
   common /m40case/ kc
   integer :: nhasharr(MAXRECENT, MAXRECENT)
   integer :: npatience, nthreads
   common /patience/ npatience, nthreads
   real :: pp(12)
   integer :: i, j, k, im
   integer*1 :: m16(16)
   character(len=12) :: calls(6)
   character(len=37) :: texts(14)

   state = 20261106_8
   kc = 0
   npatience = 1        ! four2a's FFTW_ESTIMATE_PATIENT, as jt9's default
   nthreads = 1
   nhasharr = 0
   do i = 1, 12
      pp(i) = real(sin(dble(i - 1)*3.14159265358979d0/12.0d0))
   enddo

   call direct_encode()
   call direct_hash()
   ! ---- frames before any call is known: nhasharray zero, recent_calls blank
   calls(1) = 'K1ABC'; calls(2) = 'W9XYZ'
   calls(3) = 'PJ4/K1ABC'; calls(4) = 'k1abc'
   calls(5) = ' '; calls(6) = 'N0ABC'
   call frames_first()
   ! ---- the calls the decoder has heard
   texts = ' '
   texts(1) = 'CQ N0ABC FN31'
   texts(2) = 'W1AW N0ABC R-05'
   texts(3) = 'CQ DX VK3ABC QF22'
   texts(4) = 'JA1XYZ K5ABC RR73'
   texts(5) = 'K1ABC W9XYZ EN37'
   texts(6) = 'CQ PJ4/K1ABC'
   texts(7) = '<PJ4/K1ABC> W9XYZ RR73'
   texts(8) = 'TNX BOB 73 GL'
   texts(9) = 'N0ABC W1AW 73'
   texts(10) = 'DL1ABC OH2XYZ -12'
   texts(11) = 'CQ G4ABC IO91'
   texts(12) = 'EA8ABC LU1XYZ GF05'
   texts(13) = 'F5ABC ON4XYZ JO20'
   texts(14) = 'CQ ZL1ABC RF73'
   call hear_calls(texts, 14)
   call bp_cases()
   call frames_swl()
   call sync_cases()
   call spd_cases()
   write(*, '(a,1x,i0)') 'END msk40', kc

contains

   real(8) function nextu()
      state = mod(state*1664525_8 + 1013904223_8, 4294967296_8)
      nextu = real(state, 8)/4294967296.0d0
   end function nextu

   real function gauss()
      real(8) :: t
      integer :: kk
      t = -6.0d0
      do kk = 1, 12
         t = t + nextu()
      enddo
      gauss = real(t)
   end function gauss

   subroutine differ(what)
      character(len=*) :: what
      write(*, '(a)') 'DIFFER ' // what
      error stop 'oracle: an instrumented copy differs from the library'
   end subroutine differ

   ! the 12-bit hash of two calls, as msk40decodeframe makes it
   integer function ihash_of(a, b)
      character(len=*), intent(in) :: a, b
      character*37 :: hashmsg
      integer :: iz, ih
      hashmsg = trim(a) // ' ' // trim(b)
      call fmtmsg(hashmsg, iz)
      call hash(hashmsg, 37, ih)
      ihash_of = iand(ih, 4095)
   end function ihash_of

   ! the 16 message bits of hash*16+report, most significant last (message(i) = bit i-1)
   subroutine mk_msg(ihash, irpt, m)
      integer, intent(in) :: ihash, irpt
      integer*1, intent(out) :: m(16)
      integer :: ig, ii
      ig = 16*ihash + irpt
      do ii = 1, 16
         m(ii) = int(iand(1, ishft(ig, 1 - ii)), 1)
      enddo
   end subroutine mk_msg

   subroutine direct_encode()
      integer*1 :: mm(16), cw(32), cw2(32)
      integer :: it, ii, jj, nn
      integer*8 :: h1, h2, w
      character(len=16) :: sm
      character(len=32) :: sc
      do it = 1, 16
         do ii = 1, 16
            select case (it)
            case (1)
               mm(ii) = 0
            case (2)
               mm(ii) = 1
            case (3)
               mm(ii) = int(mod(ii, 2), 1)
            case (4)
               mm(ii) = int(1 - mod(ii, 2), 1)
            case default
               mm(ii) = int(nextu()*2.0d0, 1)
            end select
         enddo
         call encode_msk40(mm, cw)
         write(sm, '(16i1)') mm
         write(sc, '(32i1)') cw
         write(*, '(a,1x,i0,1x,a,a,a,1x,a,a,a)') 'EN', it, '|', sm, '|', '|', sc, '|'
      enddo
      h1 = 1_8; h2 = 1_8
      do nn = 0, 65535
         do ii = 1, 16
            mm(ii) = int(iand(1, ishft(nn, 1 - ii)), 1)
         enddo
         call encode_msk40(mm, cw)
         do jj = 1, 32
            w = cw(jj)
            h1 = mod(h1*1000003_8 + w + 1_8, 2147483647_8)
            h2 = mod(h2*999983_8 + w*(jj + 1), 2147483629_8)
         enddo
      enddo
      write(*, '(a,1x,i0,1x,i0,1x,i0)') 'EA', 65536, h1, h2
      cw2 = cw
   end subroutine direct_encode

   subroutine direct_hash()
      character(len=60) :: txt(14)
      integer :: lens(14), it, ih, n
      character(len=60) :: buf
      txt(1) = 'K1ABC W9XYZ'
      txt(2) = ' '
      txt(3) = 'A'
      txt(4) = 'PJ4/K1ABC W9XYZ'
      txt(5) = 'K1ABC W9XYZ RRR'
      txt(6) = 'ABCDEFGHIJKL'
      txt(7) = 'ABCDEFGHIJKLM'
      txt(8) = 'ABCDEFGHIJKLMNOPQRSTUVWXY'
      txt(9) = 'N0ABC W1AW'
      txt(10) = 'W1AW N0ABC'
      txt(11) = 'abcdefghijklmnopqrstuvwxyz0123456789'
      txt(12) = '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'
      txt(13) = 'x'
      txt(14) = 'JA1XYZ K5ABC'
      lens = (/ 37, 37, 37, 37, 37, 12, 13, 25, 37, 37, 37, 37, 1, 37 /)
      do it = 1, 14
         n = lens(it)
         buf = txt(it)
         if (n > len_trim(buf)) buf(len_trim(buf) + 1:60) = ' '
         call hash(buf(1:n), n, ih)
         write(*, '(a,1x,i0,1x,i0,1x,a,a,a,1x,i0)') 'HA', it, n, '|', buf(1:n), '|', ih
      enddo
   end subroutine direct_hash

   ! pack, unpack (saving the calls), the recent calls, and nhasharray from them
   subroutine hear_calls(tx, nt)
      character(len=37), intent(in) :: tx(*)
      integer, intent(in) :: nt
      integer :: it, i3, n3, ii
      character*77 :: c77
      character*37 :: msg
      logical :: ok
      type(pack77_options) :: opt
      do it = 1, nt
         i3 = -1; n3 = -1
         c77 = ' '
         opt%record_tx_hashes = .true.
         call pack77_legacy_truncating_fallback(tx(it), i3, n3, c77, opt)
         write(*, '(a,1x,i0,1x,a,a,a,1x,i0,1x,i0,1x,a)') 'PK', it, '|', tx(it), '|', i3, n3, c77
         if (i3 < 0 .or. n3 < 0) cycle
         msg = ' '
         call unpack77(c77, 1, msg, ok)
         write(*, '(a,1x,i0,1x,a,a,a)') 'UN', it, '|', msg, '|'
         write(*, '(a,1x,i0,10(a,a13))') 'RC', it, ('|', recent_calls(ii), ii = 1, MAXRECENT)
      enddo
      call update_msk40_hasharray(nhasharr)
      write(*, '(a,1x,i0,100(1x,i0))') 'NH', 1, nhasharr
   end subroutine hear_calls

   ! ---- bpdecode40
   subroutine bp_cases()
      integer*1 :: mm(16), cw(32), dtr(16), dlib(16)
      real :: llr(32), amp, sig
      integer :: it, ii, ia, is, ntr, nlib, nf, ifl
      real, parameter :: amps(4) = (/0.5, 1.0, 2.0, 4.0/)
      real, parameter :: sigs(5) = (/0.0, 0.5, 1.0, 1.6, 2.5/)
      real :: nanv
      nanv = transfer(int(z'7FC00000'), 0.0)
      it = 0
      do ia = 1, 4
         do is = 1, 5
            do ii = 1, 12
               call rand_msg(mm)
               call encode_msk40(mm, cw)
               amp = amps(ia); sig = sigs(is)
               do nf = 1, 32
                  llr(nf) = real(2*cw(nf) - 1)*amp
                  if (sig > 0.0) llr(nf) = llr(nf) + sig*gauss()
               enddo
               ! a few flipped bits in a third of the cases
               if (mod(ii, 3) == 0) then
                  do nf = 1, 3
                     ifl = 1 + int(nextu()*32.0d0)
                     llr(ifl) = -llr(ifl)
                  enddo
               endif
               it = it + 1
               call bp_case(it, llr)
            enddo
         enddo
      enddo
      ! special inputs
      llr = 0.0; it = it + 1; call bp_case(it, llr)
      llr = 5.0; it = it + 1; call bp_case(it, llr)
      llr = -5.0; it = it + 1; call bp_case(it, llr)
      llr = 1.0e30; it = it + 1; call bp_case(it, llr)
      llr = -1.0e-30; it = it + 1; call bp_case(it, llr)
      call rand_msg(mm); call encode_msk40(mm, cw)
      do nf = 1, 32
         llr(nf) = real(2*cw(nf) - 1)*2.0
      enddo
      llr(7) = nanv; it = it + 1; call bp_case(it, llr)
      llr(7) = 2.0*real(2*cw(7) - 1); llr(20) = transfer(int(z'7F800000'), 0.0); it = it + 1; call bp_case(it, llr)
      llr(20) = 2.0*real(2*cw(20) - 1); llr(3) = -0.0; it = it + 1; call bp_case(it, llr)
      ntr = 0; nlib = 0
   end subroutine bp_cases

   subroutine rand_msg(mm)
      integer*1, intent(out) :: mm(16)
      integer :: ii
      do ii = 1, 16
         mm(ii) = int(nextu()*2.0d0, 1)
      enddo
   end subroutine rand_msg

   subroutine bp_case(it, llr)
      integer, intent(in) :: it
      real, intent(in) :: llr(32)
      integer*1 :: dtr(16), dlib(16)
      integer :: ntr, nlib, ii
      character(len=16) :: sd
      kc = it
      call prr('BI', it, llr, 32)
      dtr = -1; dlib = -1
      ntr = -99; nlib = -99
      call bpdecode40_trace(llr, 5, dtr, ntr)
      call bpdecode40(llr, 5, dlib, nlib)
      if (ntr /= nlib) call differ('bpdecode40 niterations')
      if (nlib >= 0) then
         if (any(dtr /= dlib)) call differ('bpdecode40 decoded')
         write(sd, '(16i1)') dlib
         write(*, '(a,1x,i0,1x,i0,1x,a,a,a)') 'BD', it, nlib, '|', sd, '|'
      else
         write(*, '(a,1x,i0,1x,i0,1x,a)') 'BD', it, nlib, '|-|'
      endif
   end subroutine bp_case

   ! ---- frames
   subroutine frames_first()
      integer*1 :: m(16)
      integer :: ih
      ! nothing known: nhasharray zero, recent_calls blank; the hash of the calls matches
      ih = ihash_of(calls(1), calls(2))
      call mk_msg(ih, 8, m)
      call do_frame(m, 0, 0, 0.0, 0.4, 0.0, 1.0, calls(1), calls(2), .false., 3.0, 0)
      call do_frame(m, 0, 0, 0.0, 0.4, 0.0, 1.0, calls(1), calls(2), .true., 3.0, 0)
      call mk_msg(ih, 14, m)
      call do_frame(m, 0, 0, 0.0, -2.0, 0.0, 1.0, calls(1), calls(2), .false., -3.0, 0)
      call mk_msg(ih, 15, m)
      call do_frame(m, 0, 0, 0.0, 1.0, 0.0, 1.0, calls(1), calls(2), .false., 0.0, 0)
      ! a report below 7 does not decode (the types 0..6 are for the long form)
      call mk_msg(ih, 3, m)
      call do_frame(m, 0, 0, 0.0, 1.0, 0.0, 1.0, calls(1), calls(2), .true., 0.0, 0)
      ! the wrong hash: nothing with bswl false; with bswl and nhasharray zero, hash 0 matches the zeros
      call mk_msg(ih + 1, 9, m)
      call do_frame(m, 0, 0, 0.0, 1.0, 0.0, 1.0, calls(1), calls(2), .false., 0.0, 0)
      call do_frame(m, 0, 0, 0.0, 1.0, 0.0, 1.0, calls(1), calls(2), .true., 0.0, 0)
      call mk_msg(0, 10, m)
      call do_frame(m, 0, 0, 0.0, 0.3, 0.0, 1.0, calls(1), calls(2), .true., 0.0, 0)
      ! other calls: a prefix call, lower case, blank hiscall, blank both (ihash 9999)
      ih = ihash_of(calls(3), calls(2)); call mk_msg(ih, 11, m)
      call do_frame(m, 0, 0, 0.0, 0.3, 0.0, 1.0, calls(3), calls(2), .false., 1.0, 0)
      ih = ihash_of(calls(4), calls(2)); call mk_msg(ih, 12, m)
      call do_frame(m, 0, 0, 0.0, 0.3, 0.0, 1.0, calls(4), calls(2), .false., 1.0, 0)
      call do_frame(m, 0, 0, 0.0, 0.3, 0.0, 1.0, calls(1), calls(2), .false., 1.0, 0)
      call mk_msg(ih, 12, m)
      call do_frame(m, 0, 0, 0.0, 0.3, 0.0, 1.0, calls(1), calls(5), .false., 1.0, 0)
      call mk_msg(4095, 12, m)
      call do_frame(m, 0, 0, 0.0, 0.3, 0.0, 1.0, calls(5), calls(5), .true., 1.0, 0)
      ! the same calls again (the saved ihash), and the calls changed back
      ih = ihash_of(calls(1), calls(2)); call mk_msg(ih, 13, m)
      call do_frame(m, 0, 0, 0.0, 0.7, 0.0, 1.0, calls(1), calls(2), .false., 1.0, 0)
      call do_frame(m, 0, 0, 0.0, 0.7, 0.0, 1.0, calls(1), calls(2), .false., 1.0, 0)
      ! phases round the circle, with a frequency error
      call mk_msg(ih, 8, m)
      do im = 0, 11
         call do_frame(m, 0, 0, 0.0, -3.14159265 + 0.5236*real(im), 0.0005*real(im - 5), 1.0, calls(1), calls(2), .false., 2.0, 0)
      enddo
      ! the phase at an axis and the sync polarity
      call do_frame(m, 0, 0, 0.0, 0.0, 0.0, 1.0, calls(1), calls(2), .false., 2.0, 0)
      call do_frame(m, 0, 0, 0.0, 1.5707963, 0.0, 1.0, calls(1), calls(2), .false., 2.0, 0)
      call do_frame(m, 0, 0, 0.0, 3.1415927, 0.0, 1.0, calls(1), calls(2), .false., 2.0, 0)
      call do_frame(m, 0, 0, 0.0, -1.5707963, 0.0, 1.0, calls(1), calls(2), .false., 2.0, 0)
   end subroutine frames_first

   subroutine frames_swl()
      integer*1 :: m(16)
      integer :: ih, ii, ij, i, j, ir
      real :: ph, sg
      real, parameter :: sgs(6) = (/0.1, 0.2, 0.3, 0.4, 0.5, 0.7/)
      ! the hash of two recent calls (both orders), seen by a receiver with bswl: success 2
      do ii = 1, 4
         ij = ii + 3
         ih = nhasharr(ii, ij)
         if (ih < 0) cycle
         call mk_msg(ih, 7 + ii, m)
         call do_frame(m, 0, 0, 0.0, 0.2*real(ii), 0.0, 1.0, calls(1), calls(2), .true., 1.0, 0)
         call do_frame(m, 0, 0, 0.0, 0.2*real(ii), 0.0, 1.0, calls(1), calls(2), .false., 1.0, 0)
         ih = nhasharr(ij, ii)
         call mk_msg(ih, 14, m)
         call do_frame(m, 0, 0, 0.0, 0.2*real(ii), 0.0, 1.0, calls(1), calls(2), .true., 1.0, 0)
      enddo
      ! a pair that matches more than one entry is possible only by chance; a hash that matches none: success 3
      do ii = 1, 6
         call mk_msg(1000 + 301*ii, 7 + mod(ii*3, 9), m)
         call do_frame(m, 0, 0, 0.0, 0.1*real(ii), 0.0, 1.0, calls(1), calls(2), .true., 2.0, 0)
      enddo
      ! the hash of the calls themselves with bswl: success 1 first
      ih = ihash_of(calls(6), calls(2)); call mk_msg(ih, 9, m)
      call do_frame(m, 0, 0, 0.0, 0.5, 0.0, 1.0, calls(6), calls(2), .true., 2.0, 0)
      ! noise sweep, bit errors, sync word errors, degenerate input
      ih = ihash_of(calls(1), calls(2))
      do ir = 7, 15, 2
         call mk_msg(ih, ir, m)
         do i = 1, 6
            do j = 1, 2
               ph = real(nextu()*6.2831853d0 - 3.1415927d0)
               sg = sgs(i)
               call do_frame(m, 0, 0, sg, ph, real(0.002d0*(nextu() - 0.5d0)), 1.0, calls(1), calls(2), mod(i + j, 2) == 0, &
                             real(i) - 3.0, 0)
            enddo
         enddo
      enddo
      call mk_msg(ih, 9, m)
      do i = 1, 8
         call do_frame(m, i, 0, 0.0, 0.4, 0.0, 1.0, calls(1), calls(2), .false., 1.0, 0)
      enddo
      do i = 1, 5
         call do_frame(m, 0, i, 0.0, 0.9, 0.0, 1.0, calls(1), calls(2), .false., 1.0, 0)
      enddo
      do i = 0, 5
         call do_frame(m, 0, 0, 0.0, 0.0, 0.0, 10.0**(real(i)*10.0 - 30.0), calls(1), calls(2), .true., 1.0, 0)
      enddo
      call do_frame(m, 0, 0, 0.0, 0.0, 0.0, 0.0, calls(1), calls(2), .true., 1.0, 0)
      call do_frame(m, 0, 0, 0.0, 0.0, 0.0, -1.0, calls(1), calls(2), .true., 1.0, 0)
      do i = 1, 5
         call do_frame(m, 0, 0, 0.0, 0.0, 0.0, 1.0, calls(1), calls(2), .true., 1.0, i)
      enddo
      do i = 1, 6
         call do_frame(m, 0, 0, 3.0, 0.0, 0.0, 0.0, calls(1), calls(2), .true., 1.0, 0)
      enddo
      ! many random messages through frames without noise: the hash rarely matches
      do i = 1, 40
         call rand_msg(m)
         call do_frame(m, 0, 0, 0.0, real(nextu()*6.0d0 - 3.0d0), 0.0, 1.0, calls(1), calls(2), .true., 1.0, 0)
      enddo
   end subroutine frames_swl

   ! the waveform: bits(1:40) as +-1: c = xi + j xq
   subroutine waveform(bs, xi, xq)
      integer, intent(in) :: bs(40)
      real, intent(out) :: xi(NSPM), xq(NSPM)
      integer :: ii, is
      xq = 0.0; xi = 0.0
      xq(1:6) = real(bs(1))*pp(7:12)
      do ii = 1, 19
         is = (ii - 1)*12 + 7
         xq(is:is + 11) = real(bs(2*ii + 1))*pp
      enddo
      xq(NSPM - 5:NSPM) = real(bs(1))*pp(1:6)
      do ii = 1, 20
         is = (ii - 1)*12 + 1
         xi(is:is + 11) = real(bs(2*ii))*pp
      enddo
   end subroutine waveform

   ! the 40 bits (+-1) of a message: the sync word then the codeword, nflip codeword bits and nsyncflip sync bits wrong
   subroutine make_bits(m, nflip, nsyncflip, bs)
      integer*1, intent(in) :: m(16)
      integer, intent(in) :: nflip, nsyncflip
      integer, intent(out) :: bs(40)
      integer*1 :: cw(32)
      integer :: ii, ifl
      integer, parameter :: s8r(8) = (/1, 0, 1, 1, 0, 0, 0, 1/)
      integer, parameter :: synpos(5) = (/1, 4, 2, 7, 5/)
      call encode_msk40(m, cw)
      bs(1:8) = s8r
      bs(9:40) = cw
      do ii = 1, nflip
         ifl = 8 + 1 + int(nextu()*32.0d0)
         bs(ifl) = 1 - bs(ifl)
      enddo
      do ii = 1, nsyncflip
         bs(synpos(ii)) = 1 - bs(synpos(ii))
      enddo
      bs = 2*bs - 1
   end subroutine make_bits

   ! one frame through the decoder. nflip, nsyncflip: bits wrong; sigma noise; phase, dwr (rad/sample),
   ! scale; odd: a degenerate point (1 NaN, 2 Inf, 3 NaN in the sync, 4 NaN in the data, 5 -Inf)
   subroutine do_frame(m, nflip, nsyncflip, sigma, phase, dwr, scale, mycall, hiscall, bswl, xsnr, odd)
      integer*1, intent(in) :: m(16)
      integer, intent(in) :: nflip, nsyncflip, odd
      real, intent(in) :: sigma, phase, dwr, scale, xsnr
      character(len=*), intent(in) :: mycall, hiscall
      logical, intent(in) :: bswl
      integer :: bs(40), ii
      real :: xi(NSPM), xq(NSPM), ang, cr, ci, nanv, infv
      complex :: cf(NSPM), ctr(NSPM), clib(NSPM), rot
      character(len=37) :: mtr, mlib
      character(len=12) :: mc, hc
      logical*1 :: lsw
      integer :: ntr, nlib, cbits
      nanv = transfer(int(z'7FC00000'), 0.0)
      infv = transfer(int(z'7F800000'), 0.0)
      kc = kc + 1
      call make_bits(m, nflip, nsyncflip, bs)
      call waveform(bs, xi, xq)
      do ii = 1, NSPM
         ang = phase + dwr*real(ii)
         rot = cmplx(cos(ang), sin(ang))
         cf(ii) = cmplx(xi(ii), xq(ii))*rot
         cr = real(cf(ii)); ci = aimag(cf(ii))
         if (sigma > 0.0) then
            cr = cr + sigma*gauss()
            ci = ci + sigma*gauss()
         endif
         cf(ii) = cmplx(scale*cr, scale*ci)
      enddo
      if (odd == 1) cf(30) = cmplx(nanv, aimag(cf(30)))
      if (odd == 2) cf(60) = cmplx(infv, aimag(cf(60)))
      if (odd == 3) cf(5) = cmplx(real(cf(5)), nanv)
      if (odd == 4) cf(150) = cmplx(nanv, nanv)
      if (odd == 5) cf(100) = cmplx(-infv, real(cf(100)))
      mc = mycall; hc = hiscall
      write(*, '(a,1x,i0,1x,a,a,a,a,a,1x,l1,1x,z8.8)') 'FP', kc, '|', mc, '|', hc, '|', bswl, transfer(xsnr, 0)
      call prc('FC', kc, cf, NSPM)
      write(*, '(a,1x,i0)') 'RUN', kc
      ctr = cf; clib = cf
      mtr = 'unset'; mlib = 'unset'
      ntr = -1; nlib = -1
      lsw = bswl
      call msk40decodeframe_trace(ctr, mc, hc, xsnr, lsw, nhasharr, mtr, ntr)
      call msk40decodeframe(clib, mc, hc, xsnr, lsw, nhasharr, mlib, nlib)
      if (ntr /= nlib .or. mtr /= mlib) call differ('msk40decodeframe result')
      do ii = 1, NSPM
         if (cbits(real(ctr(ii))) /= cbits(real(clib(ii))) .or. cbits(aimag(ctr(ii))) /= cbits(aimag(clib(ii)))) &
            call differ('msk40decodeframe c')
      enddo
      write(*, '(a,1x,i0,1x,i0,1x,a,a,a)') 'R', kc, nlib, '|', mlib, '|'
      call prc('O', kc, clib, NSPM)
   end subroutine do_frame

   ! ---- msk40sync
   subroutine sync_cases()
      integer*1 :: m(16)
      integer :: bs(40), ih, ii, in, ip, ncase, nsucc1, nsucc2, nfr
      integer, parameter :: navp(3, 6) = reshape((/0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1/), (/3, 6/))
      real :: xi(NSPM), xq(NSPM), cdat(2, 3*NSPM), fo, ph, sg
      complex :: cd(3*NSPM), c1(NSPM), c2(NSPM)
      real :: fest1, fest2, ang, cr, ci, dw, nanv
      integer :: nk1(2), nk2(2), jj, cbits
      integer :: navmask(3)
      real, parameter :: foffs(6) = (/0.0, 11.0, -23.0, 40.0, -3.5, 90.0/)
      real, parameter :: sgs(6) = (/0.0, 0.3, 0.6, 1.0, 1.6, 0.0/)
      ih = ihash_of(calls(1), calls(2))
      nanv = transfer(int(z'7FC00000'), 0.0)
      ncase = 0
      do ii = 1, 6
         call mk_msg(ih, 7 + ii, m)
         call make_bits(m, 0, 0, bs)
         call waveform(bs, xi, xq)
         ph = real(nextu()*6.2831853d0 - 3.1415927d0)
         dw = real(2.0d0*3.14159265358979d0*dble(1500.0 + foffs(ii))/12000.0d0)
         sg = sgs(ii)
         do jj = 1, 3*NSPM
            in = mod(jj - 1, NSPM) + 1
            ang = ph + dw*real(jj)
            cd(jj) = cmplx(xi(in), xq(in))*cmplx(cos(ang), sin(ang))
            if (sg > 0.0) then
               cr = real(cd(jj)) + sg*gauss()
               ci = aimag(cd(jj)) + sg*gauss()
               cd(jj) = cmplx(cr, ci)
            endif
         enddo
         if (ii == 6) cd(300) = cmplx(nanv, 0.0)
         ! the carrier passed in is 1500 Hz plus the carrier's error, 0 to 15 Hz off
         do ip = 1, 6
            navmask = navp(:, ip)
            ncase = ncase + 1
            kc = kc + 1
            fo = 1500.0 + foffs(ii) - 7.0*real(mod(ip, 3))
            write(*, '(a,1x,i0,1x,z8.8,3(1x,i0))') 'YP', kc, transfer(fo, 0), navmask
            call prc('SC', kc, cd, 3*NSPM)
            call run_sync(cd, navmask, fo)
         enddo
      enddo
      nsucc1 = 0; nsucc2 = 0; nfr = 0; fest1 = 0.0; fest2 = 0.0
   end subroutine sync_cases

   subroutine run_sync(cd, navmask, fo)
      complex, intent(in) :: cd(3*NSPM)
      integer, intent(in) :: navmask(3)
      real, intent(in) :: fo
      complex :: c1(NSPM), c2(NSPM)
      real :: fest1, fest2
      integer :: nk1(2), nk2(2), ns1, ns2, jj, cbits
      write(*, '(a,1x,i0)') 'YR', kc
      nk1 = -1; nk2 = -1; ns1 = -9; ns2 = -9; fest1 = -1.0; fest2 = -1.0
      c1 = 0.0; c2 = 0.0
      call msk40sync_trace(cd, 3, 29, 7.2, navmask, 2, fo, fest1, nk1, ns1, c1)
      call msk40sync(cd, 3, 29, 7.2, navmask, 2, fo, fest2, nk2, ns2, c2)
      if (ns1 /= ns2 .or. cbits(fest1) /= cbits(fest2) .or. any(nk1 /= nk2)) call differ('msk40sync result')
      do jj = 1, NSPM
         if (cbits(real(c1(jj))) /= cbits(real(c2(jj))) .or. cbits(aimag(c1(jj))) /= cbits(aimag(c2(jj)))) &
            call differ('msk40sync c')
      enddo
      write(*, '(a,1x,i0,1x,i0,1x,z8.8,2(1x,i0))') 'YA', kc, ns2, cbits(fest2), nk2
      call prc('YC', kc, c2, NSPM)
   end subroutine run_sync

   ! ---- msk40spd: a ping of the repeated frame in a block of noise
   subroutine spd_cases()
      integer*1 :: m(16)
      integer :: bs(40), ih, icase, jj, in, n0, n1, ncalls
      real :: xi(NSPM), xq(NSPM), fcar, ph, sg, amp, fc, ang, cr, ci
      complex :: z(NB)
      real, parameter :: amps(14) = (/3.0, 1.5, 1.0, 0.7, 0.5, 0.35, 3.0, 1.0, 0.5, 2.0, 1.0, 0.7, 4.0, 0.3/)
      real, parameter :: sgs(14) = (/0.05, 0.2, 0.3, 0.3, 0.3, 0.3, 0.1, 0.3, 0.3, 0.2, 0.3, 0.3, 0.05, 0.6/)
      real, parameter :: dfs(14) = (/0.0, 8.0, -15.0, 30.0, -45.0, 3.0, 60.0, -80.0, 17.0, 0.0, 5.0, -5.0, 100.0, 0.0/)
      integer, parameter :: ns0(14) = (/200, 700, 1500, 2500, 100, 3000, 4000, 1200, 800, 2000, 600, 3500, 300, 900/)
      integer, parameter :: lens(14) = (/2400, 1800, 1000, 900, 4000, 700, 2500, 3000, 4200, 1500, 2000, 1000, 2900, 5000/)
      integer, parameter :: ntols(14) = (/20, 100, 200, 50, 100, 20, 100, 100, 300, 100, 100, 100, 100, 100/)
      character(len=12) :: mc, hc
      logical :: bsw
      do icase = 1, 14
         if (mod(icase, 2) == 1) then
            mc = calls(1); hc = calls(2)
         else
            mc = calls(3); hc = calls(2)
         endif
         bsw = icase >= 6
         ih = ihash_of(mc, hc)
         if (icase == 9 .or. icase == 12) then
            ih = nhasharr(1, 2)           ! the hash of two recent calls: success 2 with bswl
         endif
         call mk_msg(ih, 7 + mod(icase, 9), m)
         call make_bits(m, 0, 0, bs)
         call waveform(bs, xi, xq)
         fc = 1500.0
         fcar = fc + dfs(icase)
         ph = real(nextu()*6.2831853d0 - 3.1415927d0)
         amp = amps(icase); sg = sgs(icase)
         n0 = ns0(icase); n1 = min(NB, n0 + lens(icase) - 1)
         do jj = 1, NB
            z(jj) = cmplx(0.0, 0.0)
            if (jj >= n0 .and. jj <= n1) then
               in = mod(jj - 1, NSPM) + 1
               ang = real(2.0d0*3.14159265358979d0*dble(fcar)*dble(jj)/12000.0d0 + dble(ph))
               z(jj) = amp*cmplx(xi(in), xq(in))*cmplx(cos(ang), sin(ang))
            endif
            cr = real(z(jj)) + sg*gauss()
            ci = aimag(z(jj)) + sg*gauss()
            z(jj) = cmplx(cr, ci)
         enddo
         kc = kc + 1
         write(*, '(a,1x,i0,1x,a,a,a,a,a,1x,l1,1x,i0,1x,z8.8)') 'ZP', kc, '|', mc, '|', hc, '|', bsw, ntols(icase), transfer(fc, 0)
         call prc('ZB', kc, z, NB)
         call run_spd(z, mc, hc, bsw, ntols(icase), fc)
      enddo
      ! a block of noise alone, a block of zeros, a block with an Inf carrier, and a block shorter than 3 frames
      do jj = 1, NB
         z(jj) = cmplx(0.3*gauss(), 0.3*gauss())
      enddo
      kc = kc + 1
      mc = calls(1); hc = calls(2)
      write(*, '(a,1x,i0,1x,a,a,a,a,a,1x,l1,1x,i0,1x,z8.8)') 'ZP', kc, '|', mc, '|', hc, '|', .true., 100, transfer(1500.0, 0)
      call prc('ZB', kc, z, NB)
      call run_spd(z, mc, hc, .true., 100, 1500.0)
      z = cmplx(0.0, 0.0)
      kc = kc + 1
      write(*, '(a,1x,i0,1x,a,a,a,a,a,1x,l1,1x,i0,1x,z8.8)') 'ZP', kc, '|', mc, '|', hc, '|', .true., 100, transfer(1500.0, 0)
      call prc('ZB', kc, z, NB)
      call run_spd(z, mc, hc, .true., 100, 1500.0)
      ncalls = 0
   end subroutine spd_cases

   subroutine run_spd(z, mc, hc, bsw, ntol, fc)
      complex, intent(in) :: z(NB)
      character(len=12), intent(in) :: mc, hc
      logical, intent(in) :: bsw
      integer, intent(in) :: ntol
      real, intent(in) :: fc
      character(len=37) :: m1, m2
      real :: fr1, tr1, fr2, tr2
      integer :: ns1, ns2, nv1, nv2, cbits
      logical*1 :: lsw
      write(*, '(a,1x,i0)') 'ZR', kc
      lsw = bsw
      m1 = 'unset'; m2 = 'unset'
      fr1 = -5.0; tr1 = -5.0; fr2 = -5.0; tr2 = -5.0
      nv1 = -5; nv2 = -5; ns1 = -9; ns2 = -9
      call msk40spd_trace(z, NB, ntol, mc, hc, lsw, nhasharr, ns1, m1, fc, fr1, tr1, nv1)
      call msk40spd(z, NB, ntol, mc, hc, lsw, nhasharr, ns2, m2, fc, fr2, tr2, nv2)
      if (ns1 /= ns2 .or. m1 /= m2 .or. cbits(fr1) /= cbits(fr2) .or. cbits(tr1) /= cbits(tr2) .or. nv1 /= nv2) &
         call differ('msk40spd result')
      write(*, '(a,1x,i0,1x,i0,1x,a,a,a,2(1x,z8.8),1x,i0)') 'X', kc, ns2, '|', m2, '|', cbits(fr2), cbits(tr2), nv2
   end subroutine run_spd

end program msk40_oracle

! ---- printing: real and complex arrays, 16 values a line
subroutine prr(tag, k, x, n)
   implicit none
   character(len=*) :: tag
   integer :: k, n, i, i2, j
   real :: x(*)
   integer :: cbits
   do i = 1, n, 16
      i2 = min(i + 15, n)
      write(*, '(a,1x,i0,1x,i0,16(1x,z8.8))') tag, k, i, (cbits(x(j)), j = i, i2)
   enddo
end subroutine prr

subroutine prc(tag, k, x, n)
   implicit none
   character(len=*) :: tag
   integer :: k, n, i, i2, j
   complex :: x(*)
   real :: v(2*n)
   integer :: cbits
   do j = 1, n
      v(2*j-1) = real(x(j))
      v(2*j) = aimag(x(j))
   enddo
   do i = 1, 2*n, 16
      i2 = min(i + 15, 2*n)
      write(*, '(a,1x,i0,1x,i0,16(1x,z8.8))') tag, k, i, (cbits(v(j)), j = i, i2)
   enddo
end subroutine prc

! ---- what the instrumented copies call (see scripts/msk144-piece-check.sh)
subroutine msk40_bp_iter(iter, zn, tov)
   implicit none
   integer :: iter, kc
   real :: zn(32), tov(3, 32)
   common /m40case/ kc
   call prr('BZ', kc*100 + iter, zn, 32)
   call prr('BV', kc*100 + iter, tov, 96)
end subroutine msk40_bp_iter

subroutine msk40_df_consts(pp, cb)
   implicit none
   real :: pp(12)
   complex :: cb(42)
   integer :: kc
   common /m40case/ kc
   call prr('KP', kc, pp, 12)
   call prc('KB', kc, cb, 42)
end subroutine msk40_df_consts

subroutine msk40_df_phase(cca, phase0)
   implicit none
   complex :: cca
   real :: phase0
   integer :: kc, cbits
   common /m40case/ kc
   write(*, '(a,1x,i0,3(1x,z8.8))') 'KH', kc, cbits(real(cca)), cbits(aimag(cca)), cbits(phase0)
end subroutine msk40_df_phase

subroutine msk40_df_corr(cfac, c)
   implicit none
   complex :: cfac, c(240)
   integer :: kc, cbits
   common /m40case/ kc
   write(*, '(a,1x,i0,2(1x,z8.8))') 'KF', kc, cbits(real(cfac)), cbits(aimag(cfac))
   call prc('KC', kc, c, 240)
end subroutine msk40_df_corr

subroutine msk40_df_soft(softbits)
   implicit none
   real :: softbits(40)
   integer :: kc
   common /m40case/ kc
   call prr('KS', kc, softbits, 40)
end subroutine msk40_df_soft

subroutine msk40_df_sync(n1)
   implicit none
   integer :: n1, kc
   common /m40case/ kc
   write(*, '(a,1x,i0,1x,i0)') 'KY', kc, n1
end subroutine msk40_df_sync

subroutine msk40_df_norm(sav, s2av, ssig)
   implicit none
   real :: sav, s2av, ssig
   integer :: kc, cbits
   common /m40case/ kc
   write(*, '(a,1x,i0,3(1x,z8.8))') 'KN', kc, cbits(sav), cbits(s2av), cbits(ssig)
end subroutine msk40_df_norm

subroutine msk40_df_llr(sigma, llr)
   implicit none
   real :: sigma, llr(32)
   integer :: kc, cbits
   common /m40case/ kc
   write(*, '(a,1x,i0,1x,z8.8)') 'KG', kc, cbits(sigma)
   call prr('KL', kc, llr, 32)
end subroutine msk40_df_llr

subroutine msk40_df_bp(decoded, niterations)
   implicit none
   integer*1 :: decoded(16)
   integer :: niterations, kc
   character(len=16) :: sd
   common /m40case/ kc
   ! decoded is written only when a codeword is found: undefined otherwise, so not printed
   if (niterations >= 0) then
      write(sd, '(16i1)') decoded
      write(*, '(a,1x,i0,1x,i0,1x,a,a,a)') 'KD', kc, niterations, '|', sd, '|'
   else
      write(*, '(a,1x,i0,1x,i0,1x,a)') 'KD', kc, niterations, '|-|'
   endif
end subroutine msk40_df_bp

subroutine msk40_df_msg(nhammd, cord, imsg, nrxrpt, nrxhash, ihash)
   implicit none
   integer :: nhammd, imsg, nrxrpt, nrxhash, ihash, kc, cbits
   real :: cord
   common /m40case/ kc
   write(*, '(a,1x,i0,1x,i0,1x,z8.8,4(1x,i0))') 'KM', kc, nhammd, cbits(cord), imsg, nrxrpt, nrxhash
   write(*, '(a,1x,i0,1x,i0)') 'KI', kc, ihash
end subroutine msk40_df_msg

subroutine msk40_sync_consts(pp, cb)
   implicit none
   real :: pp(12)
   complex :: cb(42)
   integer :: kc
   common /m40case/ kc
   call prr('YKP', kc, pp, 12)
   call prc('YKB', kc, cb, 42)
end subroutine msk40_sync_consts

subroutine msk40_sync_peak(xmax, fest, c, xcc, nfreqs)
   implicit none
   real :: xmax, fest, xcc(0:239)
   complex :: c(240)
   integer :: nfreqs, kc, cbits
   common /m40case/ kc
   write(*, '(a,1x,i0,2(1x,z8.8),1x,i0)') 'YX', kc, cbits(xmax), cbits(fest), nfreqs
   call prc('YW', kc, c, 240)
   call prr('YV', kc, xcc, 240)
end subroutine msk40_sync_peak

subroutine msk40_spd_consts(rcw)
   implicit none
   real :: rcw(12)
   integer :: kc
   common /m40case/ kc
   call prr('R', kc, rcw, 12)
end subroutine msk40_spd_consts

subroutine msk40_spd_win(ihlo, ihhi, illo, ilhi, i2000, i4000)
   implicit none
   integer :: ihlo, ihhi, illo, ilhi, i2000, i4000, kc
   common /m40case/ kc
   write(*, '(a,1x,i0,6(1x,i0))') 'N', kc, ihlo, ihhi, illo, ilhi, i2000, i4000
end subroutine msk40_spd_win

subroutine msk40_spd_step(istp, ihpk, ilpk, deltah, deltal, ah, al, trath, tratl, ferr)
   implicit none
   integer :: istp, ihpk, ilpk, kc, cbits
   real :: deltah, deltal, ah, al, trath, tratl, ferr
   common /m40case/ kc
   write(*, '(a,1x,i0,3(1x,i0),7(1x,z8.8))') 'T', kc, istp, ihpk, ilpk, cbits(deltah), cbits(deltal), &
      cbits(ah), cbits(al), cbits(trath), cbits(tratl), cbits(ferr)
end subroutine msk40_spd_step

subroutine msk40_spd_med(nstep, indices, xmed, detmet, detmet2, detfer)
   implicit none
   integer :: nstep, indices(150), kc, j, cbits
   real :: xmed, detmet(-2:153), detmet2(-2:153), detfer(150)
   common /m40case/ kc
   call prr('D', kc, detmet(1), nstep)
   call prr('H', kc, detmet2(1), nstep)
   call prr('O', kc, detfer, nstep)
   write(*, '(a,1x,i0,1x,i0,1x,z8.8,150(1x,i0))') 'M', kc, nstep, cbits(xmed), (indices(j), j = 1, nstep)
end subroutine msk40_spd_med

subroutine msk40_spd_cands(ndet, nstart, ferrs, snrs)
   implicit none
   integer :: ndet, nstart(5), kc, j, cbits
   real :: ferrs(5), snrs(5)
   common /m40case/ kc
   write(*, '(a,1x,i0,1x,i0,1x,i0)') 'E', kc, 0, ndet
   do j = 1, ndet
      write(*, '(a,1x,i0,1x,i0,1x,i0,2(1x,z8.8))') 'E', kc, j, nstart(j), cbits(ferrs(j)), cbits(snrs(j))
   enddo
end subroutine msk40_spd_cands

subroutine msk40_spd_sync(icand, iav, nsyncsuccess, fest, npkloc)
   implicit none
   integer :: icand, iav, nsyncsuccess, npkloc(10), kc, cbits
   real :: fest
   common /m40case/ kc
   write(*, '(a,1x,i0,3(1x,i0),1x,z8.8,2(1x,i0))') 'J', kc, icand, iav, nsyncsuccess, cbits(fest), npkloc(1), npkloc(2)
end subroutine msk40_spd_sync

! the bits of a real, a NaN as the default one (the sign and payload of a NaN are not compared)
integer function cbits(x)
   implicit none
   real :: x
   if (x /= x) then
      cbits = int(z'7FC00000')
   else
      cbits = transfer(x, 0)
   endif
end function cbits
