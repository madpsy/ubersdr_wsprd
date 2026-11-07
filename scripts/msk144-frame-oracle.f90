! WSJT-X's MSK144 frame decoder, signal quality, encode_128_90, genmsk_128_90,
! fmtmsg and atan2 for scripts/msk144-piece-check.sh frame: makes its own fixed
! inputs, prints them and what the reference does with them, for
! scripts/msk144-frame-check.cpp (src/msk144_decode.cpp) to match line for
! line. Built against WSJT-X's own library (libwsjt_fort.a): msk144decodeframe,
! msk144signalquality, genmsk_128_90, encode_128_90, fmtmsg, packjt77 (module).
! The piece check appends msk144decodeframe_trace to this file: a copy of
! lib/msk144decodeframe.f90 with calls added (sed, see
! scripts/msk144-piece-check.sh) that print its intermediate values. Every frame
! runs the copy and the library's routine and the oracle stops (error stop) if
! their results differ.
!
! Every float is its bits, 8 hex digits (a real*8, 16); arrays 16 values a
! line, "TAG k first-index values", complex values as re, im pairs. The module
! packjt77's tables are shared state, so the driver makes the same calls in
! the same order (PK, RUN, QS lines); the frame decoder runs twice here (the
! copy, then the library), the port's once: what it saves is the same calls
! again.
!
! Lines, in the order a case prints them
!   PK k |text| i3 n3 c77      pack77_legacy_truncating_fallback of a 37-char text
!   PN k |text| i3 n3 c77      the same, the hashes of its calls not saved (record_tx_hashes false)
!   MB k m77                   the 77 message bits the frame carries
!   C k i v...                 the frame cframe(1:864) the decoder is given
!   RUN k                      msk144decodeframe runs (the copy, then the library)
!   KP k i v... / KB k i v...  first call only: pp(1:12), cb(1:42)
!   KH k ccar ccai ccbr ccbi phase0     cca, ccb, phase0
!   KF k cfacr cfaci           cfac
!   KC k i v...                c(1:864) after c=c*conjg(cfac)
!   KS k i v...                softbits(1:144) after the matched filter
!   KY k nbadsync1 nbadsync2   the sync word errors
!   KN k sav s2av ssig         the normalisation
!   KL k i v...                llr(1:128) after 2*llr/sigma**2
!   I, J, BZ, BV, BO, BH       bpdecode128_90's iterations: zn, tov, toc and tanhtoc
!                              (scripts/msk144-bp-trace.f90)
!   KD k nharderror niter |decoded77| |cw|    bpdecode128_90's answer
!   KT k |c77| n3 i3           the message bits and its type
!   R k nsuccess |msgreceived| the answer
!   S k i v... / O k i v...    softbits(1:144) and c(1:864) as the routine left them
!   QS k |msg| snr freq t0     msk144signalquality on the frame's softbits (the
!                              line is its input); Q k nbiterrors eyeopening
! Direct checks
!   EN k m77 crcmask cw128     encode_128_90 (crcmask 0: the library's own;
!                              the line is checked equal to the library's)
!   GM k |msg| ichk itype |msgsent| tones144    genmsk_128_90
!   FM k len |in| iz |out|     fmtmsg
!   AT k n h1 h2               atan2 of n pairs: a hash of the result bits
!   AV k yhex xhex res         atan2(y, x)
! The end: END frame n (the cases run)
!
! Inputs: a 32-bit LCG (x = x*1664525 + 1013904223 mod 2**32, seed 20261106).
! Frames are the transmitted waveform (the half-sine pulses of genmsk_128_90's
! xi and xq: c = (xi, xq)), rotated by a phase and a frequency error, plus
! Gaussian noise (a sum of 12 uniforms less 6, times sigma).
!
! Blind spot: NaNs compare equal whatever their sign and payload, and print as
! the default NaN (cbits below): the instrumented copy and the library are
! different objects, and x86 returns the first operand's NaN. A port that
! gets a NaN's sign or payload wrong passes here (the atan2 sweeps hash raw bits).
!
! License: GNU GPL v3 (as WSJT-X)
program msk144_frame_oracle
   use packjt77
   implicit none
   integer, parameter :: NSPM = 864
   interface
      subroutine msk144signalquality(cframe, snr, freq, t0, softbits, msg, dxcall, btrain, datadir, nbiterrors, &
                                     eyeopening, pcoeffs)
         complex :: cframe(864)
         real :: snr, freq, t0, softbits(144), eyeopening
         character*37 msg
         character*12 dxcall
         logical*1 btrain
         character*(*) datadir
         integer :: nbiterrors
         real*8 :: pcoeffs(5)
      end subroutine msk144signalquality
   end interface
   integer*8 :: state, state2
   integer :: k, i, j, is, nmsg, imsg, nsynf
   integer*1 :: m77(77), cw(128)
   integer :: kc
   common /frcase/ kc
   logical :: ltrace
   common /msktrace/ ltrace
   character(len=37) :: texts(40), lastmsg
   real :: sig, ph, dw
   integer, parameter :: nfl(6) = (/1, 3, 8, 14, 20, 28/)
   complex :: c(NSPM)
   real, parameter :: sigs(8) = (/0.15, 0.3, 0.45, 0.6, 0.8, 1.0, 1.4, 2.0/)
   logical :: haveframe

   state = 20261106_8
   k = 0
   ltrace = .true.   ! bpdecode128_90_trace prints every iteration (scripts/msk144-bp-trace.f90)
   lastmsg = ' '

   texts = ' '
   texts(1) = 'CQ K1ABC FN42'
   texts(2) = 'K1ABC W9XYZ -05'
   texts(3) = 'W9XYZ K1ABC R-12'
   texts(4) = 'K1ABC W9XYZ RR73'
   texts(5) = 'CQ DX K1ABC FN42'
   texts(6) = '<PJ4/K1ABC> W9XYZ RR73'       ! the hash not known yet: <...>
   texts(7) = 'CQ PJ4/K1ABC'                 ! saves PJ4/K1ABC
   texts(8) = '<PJ4/K1ABC> W9XYZ RR73'       ! now it resolves
   texts(9) = 'W9XYZ <PJ4/K1ABC> RRR'
   texts(10) = 'TNX BOB 73 GL'
   texts(11) = '123456789ABCDEF012'
   texts(12) = 'K1ABC/R W9XYZ/R R FN42'
   texts(13) = 'CQ POTA K1ABC FN42'
   texts(14) = 'CQ K1ABC'
   texts(15) = 'W9XYZ K1ABC 73'
   texts(16) = 'K1ABC W9XYZ +10'
   texts(17) = 'CQ NA K1ABC FN42'
   texts(18) = 'K1ABC RR73; W9XYZ <KH1/KH7Z> -08'  ! DXpedition (i3 1 n3: refused by the decoder)
   texts(19) = 'TU; K1ABC W9XYZ 579 NJ'             ! ARRL RTTY (i3 3: refused)
   texts(20) = '<...> K1ABC RRR'
   texts(21) = 'HELLO WORLD'
   texts(22) = 'K1ABC W9XYZ'
   texts(23) = 'W9XYZ K1ABC EN34'
   texts(24) = 'K1ABC <PJ4/K1ABC> -03'
   nmsg = 24

   call direct_encode()
   call direct_fmtmsg()
   call direct_atan2()
   ! ---- clean frames of each text, phase and noise as drawn
   do imsg = 1, nmsg
      call text_bits(imsg, m77, haveframe, imsg == 6)   ! text 6: its hash not saved by the packing
      if (.not. haveframe) cycle
      call do_frame(m77, 0, 0, 0.0, 0.3 + 0.2*real(imsg), 0.0, 1.0, .true., texts(imsg))
   enddo

   ! ---- phases round the circle, one message, no noise, with a frequency error
   do i = 0, 11
      ph = -3.14159265 + 0.5236*real(i)
      call text_bits(2, m77, haveframe)
      call do_frame(m77, 0, 0, 0.0, ph, 0.0005*real(i - 5), 1.0, .true., texts(2))
   enddo
   call text_bits(4, m77, haveframe)
   call do_frame(m77, 0, 0, 0.0, 0.7853982, 0.0, 1.0, .true., texts(4))      ! the phase at an axis: cca+ccb on a diagonal
   call do_frame(m77, 0, 0, 0.0, 1.5707963, 0.0, 1.0, .true., texts(4))
   call do_frame(m77, 0, 0, 0.0, 3.1415927, 0.0, 1.0, .true., texts(4))
   call do_frame(m77, 0, 0, 0.0, 0.0, 0.0, 1.0, .true., texts(4))
   call do_frame(m77, 0, 0, 0.0, -1.5707963, 0.0, 1.0, .true., texts(4))
   call do_frame(m77, 0, 0, 0.0, -3.1415927, 0.0, 1.0, .true., texts(4))

   ! ---- noise sweeps: 8 levels, 5 messages, two draws: some decode, some do not
   do imsg = 1, 5
      call text_bits(imsg, m77, haveframe)
      do i = 1, 8
         do j = 1, 2
            ph = real(nextu()*6.2831853d0 - 3.1415927d0)
            call do_frame(m77, 0, 0, sigs(i), ph, real(0.002d0*(nextu() - 0.5d0)), 1.0, .true., texts(imsg))
         enddo
      enddo
   enddo

   ! ---- bit errors in the codeword: 1, 3, 8, 14, 20, 28 flips, noise-free: the LDPC corrects some
   call text_bits(1, m77, haveframe)
   do i = 1, 6
      call do_frame(m77, nfl(i), 0, 0.0, 0.4, 0.0, 1.0, .true., texts(1))
   enddo
   ! the CRC wrong with the parity right (the codeword of message+CRC xor the mask): never decodes
   call text_bits(2, m77, haveframe)
   call do_frame(m77, 0, 0, 0.0, -0.4, 0.0, 1.0, .true., texts(2), 1)
   call do_frame(m77, 0, 0, 0.0, -0.4, 0.0, 1.0, .true., texts(2), 4096)
   call do_frame(m77, 0, 0, 0.3, 2.0, 0.0, 1.0, .true., texts(2), 37)
   call do_frame(m77, 0, 0, 0.0, 2.0, 0.0, 1.0, .true., texts(2), 8191)

   ! ---- sync word errors: 0..8 of the 16 sync bits wrong (nbadsync 4 still goes on, 5 returns)
   call text_bits(4, m77, haveframe)
   do nsynf = 1, 8
      call do_frame(m77, 0, nsynf, 0.0, 0.9, 0.0, 1.0, .true., texts(4))
   enddo

   ! ---- every (n3, i3): 64 messages of random bits with the type bits set, clean frames
   do is = 0, 63
      do i = 1, 77
         m77(i) = int(nextu()*2.0d0, 1)
      enddo
      call set_type(m77, is/8, mod(is, 8))
      call do_frame(m77, 0, 0, 0.0, 0.1*real(is) - 3.0, 0.0, 1.0, .false., 'random 77 bits')
   enddo

   ! ---- frames that are not frames: zeros, NaN, Inf, scaled
   call text_bits(1, m77, haveframe)
   call do_frame(m77, 0, 0, 0.0, 0.0, 0.0, 0.0, .true., texts(1))             ! scale 0: all zeros
   call do_frame(m77, 0, 0, 0.0, 0.0, 0.0, 1.0e-30, .true., texts(1))
   call do_frame(m77, 0, 0, 0.0, 0.0, 0.0, 1.0e20, .true., texts(1))
   call do_frame(m77, 0, 0, 0.0, 0.0, 0.0, 1.0e30, .true., texts(1))
   call do_frame(m77, 0, 0, 0.0, 0.0, 0.0, 1.0e-3, .true., texts(1))
   call do_frame(m77, 0, 0, 0.0, 0.0, 0.0, -1.0, .true., texts(1))
   call do_frame(m77, 0, 0, 0.0, 0.0, 0.0, 1.0, .true., texts(1), 0, 1)     ! one NaN point
   call do_frame(m77, 0, 0, 0.0, 0.0, 0.0, 1.0, .true., texts(1), 0, 2)     ! one Inf point
   call do_frame(m77, 0, 0, 0.0, 0.0, 0.0, 1.0, .true., texts(1), 0, 3)     ! a NaN in the sync
   call do_frame(m77, 0, 0, 0.0, 0.0, 0.0, 1.0, .true., texts(1), 0, 4)     ! a NaN in a data point
   call do_frame(m77, 0, 0, 0.0, 0.0, 0.0, 1.0, .true., texts(1), 0, 5)     ! -Inf point in the second sync
   ! noise only
   do i = 1, 6
      call do_frame(m77, 0, 0, 3.0, 0.0, 0.0, 0.0, .true., texts(1))
   enddo

   ! genmsk_128_90 packs 'CQ PJ4/K1ABC' and the '<PJ4/K1ABC>' message: it runs last, so that it does not save the
   ! hash before the frames meet it unknown and then known
   call direct_genmsk()

   write(*, '(a,1x,i0)') 'END frame', kc

contains

   ! the text of case imsg as 77 bits (pack77_legacy_truncating_fallback as genmsk_128_90 calls it); printed
   ! nosave: the encoded calls' hashes are not saved (pack77_options record_tx_hashes false), so that the
   ! decoder meets the hash unknown; the PN line instead of PK
   subroutine text_bits(it, m, ok, nosave)
      integer, intent(in) :: it
      integer*1, intent(out) :: m(77)
      logical, intent(out) :: ok
      logical, intent(in), optional :: nosave
      integer :: i3, n3, jj
      character*77 :: c77
      type(pack77_options) :: opt
      logical :: ns
      i3 = -1; n3 = -1
      c77 = ' '
      ns = .false.
      if (present(nosave)) ns = nosave
      opt%record_tx_hashes = .not. ns
      call pack77_legacy_truncating_fallback(texts(it), i3, n3, c77, opt)
      write(*, '(a,1x,i0,1x,a,a,a,1x,i0,1x,i0,1x,a)') merge('PN', 'PK', ns), it, '|', texts(it), '|', i3, n3, c77
      m = 0
      ok = i3 >= 0 .and. n3 >= 0
      if (ok) then
         do jj = 1, 77
            if (c77(jj:jj) == '1') m(jj) = 1
         enddo
      endif
   end subroutine text_bits

   subroutine set_type(m, n3, i3)
      integer*1, intent(inout) :: m(77)
      integer, intent(in) :: n3, i3
      integer :: jj
      do jj = 1, 3
         m(71 + jj) = int(iand(ishft(n3, jj - 3), 1), 1)
         m(74 + jj) = int(iand(ishft(i3, jj - 3), 1), 1)
      enddo
   end subroutine set_type

   ! encode_128_90 with the CRC xor-ed with crcmask before the parity (a copy of the library's routine)
   subroutine encode_mask(message77, crcmask, codeword)
      use, intrinsic :: iso_c_binding
      use crc
      integer, intent(in) :: crcmask
      integer*1, intent(in) :: message77(77)
      integer*1, intent(out) :: codeword(128)
      integer, parameter :: NN = 128, KK = 90, MM = NN - KK
      character*90 tmpchar
      integer*1 gen(MM, KK)
      integer*1 message(KK)
      integer*1 pchecks(MM)
      integer*1, target :: i1MsgBytes(12)
      integer :: ii, jj, jjj, icol, ibmax, istr, nsum, ncrc13
      include "ldpc_128_90_generator.f90"
      gen = 0
      do ii = 1, MM
         do jj = 1, 23
            read(g(ii)(jj:jj), "(Z1)") istr
            ibmax = 4
            if (jj .eq. 23) ibmax = 2
            do jjj = 1, ibmax
               icol = (jj - 1)*4 + jjj
               if (btest(istr, 4 - jjj)) gen(ii, icol) = 1
            enddo
         enddo
      enddo
      write(tmpchar, '(77i1)') message77
      tmpchar(78:80) = '000'
      i1MsgBytes = 0
      read(tmpchar, '(10b8)') i1MsgBytes(1:10)
      ncrc13 = crc13(c_loc(i1MsgBytes), 12)
      ncrc13 = ieor(ncrc13, crcmask)
      write(tmpchar(78:90), '(b13)') ncrc13
      read(tmpchar, '(90i1)') message
      do ii = 1, MM
         nsum = 0
         do jj = 1, KK
            nsum = nsum + message(jj)*gen(ii, jj)
         enddo
         pchecks(ii) = mod(nsum, 2)
      enddo
      codeword(1:KK) = message
      codeword(KK + 1:NN) = pchecks
   end subroutine encode_mask

   ! one frame through the decoder and (when it decoded, or when asked) the signal quality.
   ! nflip codeword bits wrong; nsyncflip sync bits wrong; sigma noise; phase, frequency error (rad/sample),
   ! scale; crcmask: the CRC xor-ed; odd: a degenerate point (1 NaN, 2 Inf, 3 NaN in the sync, 4 NaN in
   ! the data, 5 -Inf in the second sync word)
   subroutine do_frame(m, nflip_, nsyncflip, sigma, phase, dwr, scale, dosq, text, crcmask, odd)
      integer*1, intent(in) :: m(77)
      integer, intent(in) :: nflip_, nsyncflip
      real, intent(in) :: sigma, phase, dwr, scale
      logical, intent(in) :: dosq
      character(len=*), intent(in) :: text
      integer, intent(in), optional :: crcmask, odd
      integer*1 :: cwk(128), bitseq(144), s8(8)
      real :: pp(12), xi(864), xq(864), ang, cr, ci
      complex :: cf(864), ctr(864), clib(864), rot
      real :: sbtr(144), sblib(144), eye, snr, freq, t0
      real*8 :: pcoef(5)
      logical*1 :: btrain
      character(len=37) :: mtr, mlib
      integer :: ii, ntr, nlib, nbit, nodd, mask, ifl, pos, nsq, is, cbits
      integer, parameter :: synpos(8) = (/1, 57, 3, 60, 5, 62, 7, 64/)
      real :: nanv, infv
      data s8 /0, 1, 1, 1, 0, 0, 1, 0/
      nanv = transfer(int(z'7FC00000'), 0.0)
      infv = transfer(int(z'7F800000'), 0.0)
      kc = kc + 1
      nodd = 0
      if (present(odd)) nodd = odd
      mask = 0
      if (present(crcmask)) mask = crcmask
      write(*, '(a,1x,i0,1x,77i1)') 'MB', kc, m
      call encode_mask(m, mask, cwk)
      if (mask == 0) call check_enc(m, cwk)
      ! codeword bits flipped
      do ii = 1, nflip_
         ifl = 1 + int(nextu()*128.0d0)
         cwk(ifl) = 1 - cwk(ifl)
      enddo
      bitseq(1:8) = s8
      bitseq(9:56) = cwk(1:48)
      bitseq(57:64) = s8
      bitseq(65:144) = cwk(49:128)
      ! sync word bits wrong: positions spread over the two words
      do ii = 1, nsyncflip
         pos = synpos(ii)
         bitseq(pos) = 1 - bitseq(pos)
      enddo
      bitseq = 2*bitseq - 1
      do ii = 1, 12
         pp(ii) = real(sin(dble(ii - 1)*3.14159265358979d0/12.0d0))
      enddo
      xq = 0.0; xi = 0.0
      xq(1:6) = real(bitseq(1))*pp(7:12)
      do ii = 1, 71
         is = (ii - 1)*12 + 7
         xq(is:is + 11) = real(bitseq(2*ii + 1))*pp
      enddo
      xq(864 - 5:864) = real(bitseq(1))*pp(1:6)
      do ii = 1, 72
         is = (ii - 1)*12 + 1
         xi(is:is + 11) = real(bitseq(2*ii))*pp
      enddo
      do ii = 1, 864
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
      if (nodd == 1) cf(100) = cmplx(nanv, aimag(cf(100)))
      if (nodd == 2) cf(200) = cmplx(infv, aimag(cf(200)))
      if (nodd == 3) cf(10) = cmplx(real(cf(10)), nanv)
      if (nodd == 4) cf(500) = cmplx(nanv, nanv)
      if (nodd == 5) cf(350) = cmplx(-infv, real(cf(350)))
      call prc('C', kc, cf, 864)
      write(*, '(a,1x,i0)') 'RUN', kc
      ctr = cf; clib = cf
      sbtr = -3.0; sblib = -3.0
      mtr = 'unset'; mlib = 'unset'
      ntr = -1; nlib = -1
      call msk144decodeframe_trace(ctr, sbtr, mtr, ntr)
      call msk144decodeframe(clib, sblib, mlib, nlib)
      ! (NaNs compare equal whatever their sign and payload: the copy and the library are different objects
      ! and x86 returns the first operand's NaN, so which one's sign survives follows the operand order)
      if (ntr /= nlib .or. mtr /= mlib) call differ()
      do ii = 1, 144
         if (cbits(sbtr(ii)) /= cbits(sblib(ii))) call differ()
      enddo
      do ii = 1, 864
         if (cbits(real(ctr(ii))) /= cbits(real(clib(ii))) .or. cbits(aimag(ctr(ii))) /= cbits(aimag(clib(ii)))) call differ()
      enddo
      write(*, '(a,1x,i0,1x,i0,1x,a,a,a)') 'R', kc, nlib, '|', mlib, '|'
      call prr('S', kc, sblib, 144)
      call prc('O', kc, clib, 864)
      ! signal quality: on a decode; on a failure when asked (the last decoded message, these softbits)
      nsq = 0
      if (nlib > 0) then
         nsq = 1
         lastmsg = mlib
      else if (dosq .and. lastmsg /= ' ' .and. mod(kc, 3) == 0) then
         nsq = 2
      endif
      if (nsq > 0) then
         snr = 5.0 + real(kc)*0.1
         freq = 1500.0 + real(kc)
         t0 = 0.5*real(kc)
         btrain = .false.
         pcoef = 0.0d0
         write(*, '(a,1x,i0,1x,a,a,a,3(1x,z8.8))') 'QS', kc, '|', lastmsg, '|', transfer(snr, 0), transfer(freq, 0), &
            transfer(t0, 0)
         nbit = -99
         eye = -7.0
         call msk144signalquality(clib, snr, freq, t0, sblib, lastmsg, 'K1ABC       ', btrain, '/tmp', nbit, eye, pcoef)
         write(*, '(a,1x,i0,1x,i0,1x,z8.8)') 'Q', kc, nbit, cbits(eye)
      endif
   end subroutine do_frame

   subroutine differ()
      error stop 'oracle: msk144decodeframe_trace differs from msk144decodeframe'
   end subroutine differ

   subroutine check_enc(m, cwk)
      integer*1, intent(in) :: m(77), cwk(128)
      integer*1 :: cwl(128)
      call encode_128_90(m, cwl)
      if (any(cwl /= cwk)) error stop 'oracle: encode_mask differs from encode_128_90'
   end subroutine check_enc

   subroutine direct_encode()
      integer*1 :: mm(77), cwl(128)
      integer :: it, ii
      do it = 1, 14
         do ii = 1, 77
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
         call encode_128_90(mm, cwl)
         write(*, '(a,1x,i0,1x,77i1,1x,i0,1x,128i1)') 'EN', it, mm, 0, cwl
      enddo
   end subroutine direct_encode

   subroutine direct_fmtmsg()
      character(len=100) :: buf
      integer :: it, n, iz
      character(len=100) :: ins(16)
      integer :: lens(16)
      ins(1) = 'cq k1abc  fn42';                 lens(1) = 37
      ins(2) = '  leading blanks';               lens(2) = 37
      ins(3) = 'ALREADY CLEAN';                  lens(3) = 13
      ins(4) = '                    ';           lens(4) = 20
      ins(5) = 'a   b    c     d';               lens(5) = 22
      ins(6) = 'x';                              lens(6) = 1
      ins(7) = ' ';                              lens(7) = 1
      ins(8) = 'Mixed Case With  Two  Spaces  ';  lens(8) = 37
      ins(9) = 'abc';                            lens(9) = 0
      ins(10) = 'a' // repeat(' ', 70) // 'b';   lens(10) = 80
      ins(11) = 'tnx bob  73  gl  ';             lens(11) = 18
      ins(12) = '{}~[\]^_`@:;';                  lens(12) = 11
      ins(13) = 'az AZ09 +-./?';                 lens(13) = 13
      ins(14) = 'a  ';                           lens(14) = 3
      ins(15) = '  ';                            lens(15) = 2
      ins(16) = 'k1abc w9xyz  rr73' // repeat(' ', 40); lens(16) = 60
      do it = 1, 16
         n = lens(it)
         buf = ins(it)
         if (n > 0) then
            call fmtmsg(buf(1:n), iz)
            write(*, '(a,1x,i0,1x,i0,1x,a,a,a,1x,i0,1x,a,a,a)') 'FM', it, n, '|', ins(it)(1:n), '|', iz, '|', buf(1:n), '|'
         else
            call fmtmsg(buf(1:0), iz)
            write(*, '(a,1x,i0,1x,i0,1x,a,1x,i0,1x,a)') 'FM', it, n, '||', iz, '||'
         endif
      enddo
   end subroutine direct_fmtmsg

   subroutine direct_atan2()
      integer :: n, idx
      integer*8 :: h1, h2, w
      real :: y, x, z
      integer :: ib
      real :: sp(14)
      integer :: ii, jj
      integer :: kk
      real*8 :: sc(0:13), u1, u2
      real, volatile :: one      ! not folded: gfortran turns atan2(y, 1.0) into atanf(y)
      data sc /1.0d-7, 1.0d-6, 1.0d-5, 1.0d-4, 1.0d-3, 1.0d-2, 1.0d-1, 1.0d0, 1.0d1, 1.0d2, 1.0d3, 1.0d4, 1.0d5, 1.0d6/
      state2 = 777_8
      one = 1.0
      ! sweep 1: random values over many decades (draws in the order y value, y scale, x value, x scale)
      h1 = 1_8; h2 = 1_8
      do idx = 1, 1500000
         u1 = nextu2(); u2 = nextu2()
         y = real((u1 - 0.5d0)*2.0d0*sc(int(u2*14.0d0)))
         u1 = nextu2(); u2 = nextu2()
         x = real((u1 - 0.5d0)*2.0d0*sc(int(u2*14.0d0)))
         z = atan2(y, x)
         w = iand(int(transfer(z, 0), 8), 4294967295_8)
         h1 = mod(h1*1000003_8 + w, 2147483647_8)
         h2 = mod(h2*999983_8 + w, 2147483629_8)
      enddo
      write(*, '(a,1x,i0,1x,i0,1x,i0,1x,i0)') 'AT', 1, 1500000, h1, h2
      ! sweep 2: x = 1 (atanf of y; x kept out of the compiler's reach) for every 251st float bit pattern
      h1 = 1_8; h2 = 1_8
      n = 0
      do ib = -2147483647 - 1, 2147483647 - 251, 251
         y = transfer(ib, 0.0)
         z = atan2(y, one)
         w = iand(int(transfer(z, 0), 8), 4294967295_8)
         h1 = mod(h1*1000003_8 + w, 2147483647_8)
         h2 = mod(h2*999983_8 + w, 2147483629_8)
         n = n + 1
      enddo
      write(*, '(a,1x,i0,1x,i0,1x,i0,1x,i0)') 'AT', 2, n, h1, h2
      ! sweep 3: |y/x| in each atanf argument-reduction range, both signs of both, y = ratio, x = +-1.5 (not 1)
      h1 = 1_8; h2 = 1_8
      n = 0
      do idx = 1, 1000000
         z = real(nextu2()*12.0d0)       ! a ratio in 0..12
         if (idx > 500000) z = real(nextu2()*0.75d0)
         y = z*1.5
         x = 1.5
         do ii = 1, 4
            if (ii == 2 .or. ii == 4) y = -y
            if (ii >= 3) x = -1.5
            z = atan2(y, x)
            w = iand(int(transfer(z, 0), 8), 4294967295_8)
            h1 = mod(h1*1000003_8 + w, 2147483647_8)
            h2 = mod(h2*999983_8 + w, 2147483629_8)
            n = n + 1
         enddo
      enddo
      write(*, '(a,1x,i0,1x,i0,1x,i0,1x,i0)') 'AT', 3, n, h1, h2
      ! special values: zeros, infinities, NaN, huge, tiny, denormal, 1.0
      sp = (/0.0, -0.0, 1.0, -1.0, 2.0, 1.0e-30, -1.0e30, transfer(int(z'7F800000'), 0.0), transfer(int(z'FF800000'), 0.0), &
             transfer(int(z'7FC00000'), 0.0), transfer(1, 0.0), transfer(int(z'80000001'), 0.0), 1.0e38, 3.0e-39/)
      kk = 0
      do ii = 1, 14
         do jj = 1, 14
            kk = kk + 1
            y = sp(ii); x = sp(jj)
            z = atan2(y, x)
            write(*, '(a,1x,i0,3(1x,z8.8))') 'AV', kk, transfer(y, 0), transfer(x, 0), transfer(z, 0)
         enddo
      enddo
   end subroutine direct_atan2

   subroutine direct_genmsk()
      character*37 :: msgs(12), msgsent
      integer*4 :: i4tone(144)
      integer :: it, itype, ichk
      msgs = ' '
      msgs(1) = 'CQ K1ABC FN42'
      msgs(2) = 'K1ABC W9XYZ -05'
      msgs(3) = '  K1ABC W9XYZ RR73'
      msgs(4) = 'CQ PJ4/K1ABC'
      msgs(5) = '<PJ4/K1ABC> W9XYZ RR73'
      msgs(6) = 'TNX BOB 73 GL'
      msgs(7) = '<...> K1ABC RRR'
      msgs(8) = '!!!!!'
      msgs(9) = 'k1abc w9xyz -05'
      msgs(10) = ' '
      msgs(11) = 'K1ABC W9XYZ -05' // char(0) // 'JUNK'
      msgs(12) = 'THIS MESSAGE IS FAR TOO LONG TO PACK'
      do it = 1, 12
         ichk = 0
         i4tone = -77
         msgsent = 'unset'
         itype = -99
         call genmsk_128_90(msgs(it), ichk, msgsent, i4tone, itype)
         write(*, '(a,1x,i0,1x,a,a,a,1x,i0,1x,i0,1x,a,a,a,1x,144i1)') 'GM', it, '|', msgs(it), '|', ichk, itype, '|', &
            msgsent, '|', max(0, min(9, i4tone))
      enddo
   end subroutine direct_genmsk

   real(8) function nextu2()
      state2 = mod(state2*1664525_8 + 1013904223_8, 4294967296_8)
      nextu2 = real(state2, 8)/4294967296.0d0
   end function nextu2

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

end program msk144_frame_oracle

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

! ---- what the instrumented copy of msk144decodeframe calls (see scripts/msk144-piece-check.sh)
subroutine msk144_df_consts(pp, cb)
   implicit none
   real :: pp(12)
   complex :: cb(42)
   integer :: kc
   common /frcase/ kc
   call prr('KP', kc, pp, 12)
   call prc('KB', kc, cb, 42)
end subroutine msk144_df_consts

subroutine msk144_df_phase(cca, ccb, phase0)
   implicit none
   complex :: cca, ccb
   real :: phase0
   integer :: kc
   common /frcase/ kc
   integer :: cbits
   write(*, '(a,1x,i0,5(1x,z8.8))') 'KH', kc, cbits(real(cca)), cbits(aimag(cca)), cbits(real(ccb)), &
      cbits(aimag(ccb)), cbits(phase0)
end subroutine msk144_df_phase

subroutine msk144_df_corr(cfac, c)
   implicit none
   complex :: cfac, c(864)
   integer :: kc
   common /frcase/ kc
   integer :: cbits
   write(*, '(a,1x,i0,2(1x,z8.8))') 'KF', kc, cbits(real(cfac)), cbits(aimag(cfac))
   call prc('KC', kc, c, 864)
end subroutine msk144_df_corr

subroutine msk144_df_soft(softbits)
   implicit none
   real :: softbits(144)
   integer :: kc
   common /frcase/ kc
   call prr('KS', kc, softbits, 144)
end subroutine msk144_df_soft

subroutine msk144_df_sync(n1, n2)
   implicit none
   integer :: n1, n2, kc
   common /frcase/ kc
   write(*, '(a,1x,i0,2(1x,i0))') 'KY', kc, n1, n2
end subroutine msk144_df_sync

subroutine msk144_df_norm(sav, s2av, ssig)
   implicit none
   real :: sav, s2av, ssig
   integer :: kc
   common /frcase/ kc
   integer :: cbits
   write(*, '(a,1x,i0,3(1x,z8.8))') 'KN', kc, cbits(sav), cbits(s2av), cbits(ssig)
end subroutine msk144_df_norm

subroutine msk144_df_llr(llr)
   implicit none
   real :: llr(128)
   integer :: kc
   common /frcase/ kc
   call prr('KL', kc, llr, 128)
end subroutine msk144_df_llr

subroutine msk144_df_bp(decoded77, cw, nharderror, niterations)
   implicit none
   integer*1 :: decoded77(77), cw(128)
   integer :: nharderror, niterations, kc
   common /frcase/ kc
   ! decoded77 is written only on success: undefined otherwise, so not printed
   if (nharderror >= 0) then
      write(*, '(a,1x,i0,1x,i0,1x,i0,1x,77i1,1x,128i1)') 'KD', kc, nharderror, niterations, decoded77, cw
   else
      write(*, '(a,1x,i0,1x,i0,1x,i0,1x,a,1x,128i1)') 'KD', kc, nharderror, niterations, '-', cw
   endif
end subroutine msk144_df_bp

subroutine msk144_df_tag(c77, n3, i3)
   implicit none
   character*77 :: c77
   integer :: n3, i3, kc
   common /frcase/ kc
   write(*, '(a,1x,i0,1x,a,1x,i0,1x,i0)') 'KT', kc, c77, n3, i3
end subroutine msk144_df_tag

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
