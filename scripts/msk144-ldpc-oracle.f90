! WSJT-X's (128,90) LDPC decoder, CRC-13 check, platanh and tanh, for
! scripts/msk144-piece-check.sh ldpc: generates its own fixed inputs and
! prints them and what the reference does with them, for
! scripts/msk144-ldpc-check.cpp (src/msk144_ldpc.cpp) to match line for line.
! Built against WSJT-X's own library (libwsjt_fort.a): encode_128_90,
! bpdecode128_90, chkcrc13a, platanh. The piece check appends
! bpdecode128_90_trace to this file: lib/bpdecode128_90.f90 with one call
! added at the top of each iteration (see scripts/msk144-piece-check.sh).
!
!   L MAXITER TRACE HEX*128   the soft bits (float bits, 8 hex digits each)
!   I iter zn-sum tov-sum zn(1) zn(65) zn(128) tov(1,1) tov(2,64) tov(3,128)
!                             (hex floats) the trace, at the top of each
!                             iteration, for the vectors with TRACE 1; sums in
!                             a fixed order (zn i=1..128; tov j=1..128, kk=1..3)
!   R maxiter iter nharderror message77 cw128 nbadcrc
!                             bpdecode128_90 with apmask 0 on the L line above;
!                             message77 carries over from vector to vector (it
!                             is written only on success); nbadcrc is chkcrc13a
!                             on cw(1:90)
!   C BITS90 nbadcrc          chkcrc13a on a 90-bit word
!   P x y                     platanh(x, y), hex floats: 2000 x in [-1.1, 1.1]
!                             and each branch boundary and zero, with +-1 and
!                             +-2 ulp
!   T x y                     tanh(x) of a default real as bpdecode128_90 has
!                             it: x from -30..30, dense near 0, at the edges
!                             of tanhf's branches (2**-28, 1, 22, the
!                             saturation near 9.01), with +-1 and +-2 ulp
!
! Inputs: a 32-bit LCG (x = x*1664525 + 1013904223 mod 2**32, seed 20261106);
! 200 random 77-bit messages encoded by encode_128_90, soft bits
! (2/s**2)*(+-1 + s*g) with g a sum of 12 uniforms less 6, at s = 0.55, 0.75
! and 0.95 in turn (maxiter 10, as msk144decodeframe) and every fifth at
! maxiter 30; then 40 of pure noise; traces for vectors 1-6 and 201-203;
! then 90-bit words for chkcrc13a: the first 90 bits of 40 codewords, each
! also with one bit flipped, and 40 random; then the P and T sweeps.
!
! License: GNU GPL v3 (as WSJT-X)
program msk144_ldpc_oracle
   implicit none
   integer*8 :: state
   integer*1 :: apmask(128), msg(77), cw(128), cwout(128), m77(77), word(90), enc(40,128)
   real :: llr(128), sigma(0:2), s
   integer :: ibits(128), iv, i, j, k, maxiter, iter, nharderror, nbadcrc, nenc, ntrace, nrun
   logical :: ltrace
   common /msktrace/ ltrace
   external platanh
   data sigma/0.55, 0.75, 0.95/

   apmask = 0
   state = 20261106_8
   m77 = 0
   nenc = 0
   nrun = 0
   do iv = 1, 240
      if (iv .le. 200) then
         do i = 1, 77
            msg(i) = int(nextu() * 2.0d0, 1)
         enddo
         call encode_128_90(msg, cw)
         if (nenc .lt. 40) then
            nenc = nenc + 1
            enc(nenc, :) = cw
         endif
         s = sigma(mod(iv, 3))
         do i = 1, 128
            llr(i) = (2.0 / (s * s)) * (real(2 * cw(i) - 1) + s * gauss())
         enddo
      else
         do i = 1, 128
            llr(i) = 4.0 * gauss()
         enddo
      endif
      maxiter = 10
      if (mod(iv, 5) .eq. 0) maxiter = 30
      ltrace = iv .le. 6 .or. (iv .ge. 201 .and. iv .le. 203)
      ntrace = 0
      if (ltrace) ntrace = 1
      ibits = transfer(llr, ibits)
      write(*, '(a,i0,1x,i0,1x,128z8.8)') 'L ', maxiter, ntrace, ibits
      nharderror = -99
      iter = -99
      cwout = 0
      call bpdecode128_90_trace(llr, apmask, maxiter, m77, cwout, nharderror, iter)
      word = cwout(1:90)
      call chkcrc13a(word, nbadcrc)
      write(*, '(a,i0,1x,i0,1x,i0,1x,77i1,1x,128i1,1x,i0)') 'R ', maxiter, iter, nharderror, m77, cwout, nbadcrc
      nrun = nrun + 1
   enddo

   do iv = 1, 40
      word = enc(iv, 1:90)
      call chkcrc13a(word, nbadcrc)
      write(*, '(a,90i1,1x,i0)') 'C ', word, nbadcrc
      i = 1 + int(nextu() * 90.0d0)
      word(i) = 1 - word(i)
      call chkcrc13a(word, nbadcrc)
      write(*, '(a,90i1,1x,i0)') 'C ', word, nbadcrc
      do j = 1, 90
         word(j) = int(nextu() * 2.0d0, 1)
      enddo
      call chkcrc13a(word, nbadcrc)
      write(*, '(a,90i1,1x,i0)') 'C ', word, nbadcrc
   enddo

   ! platanh: the random sweep, then the boundaries and neighbours
   do k = 1, 2000
      call pl(real(2.2d0 * nextu() - 1.1d0))
   enddo
   call pl(0.0)
   call pl(-0.0)
   do k = 1, 7
      call pnear(pb(k))
      call pnear(-pb(k))
   enddo

   ! tanh: random over -30..30, dense near 0, and the edges
   do k = 1, 600
      call th(real(60.0d0 * nextu() - 30.0d0))
   enddo
   do k = 1, 600
      call th(real(2.0d0 * nextu() - 1.0d0))
   enddo
   do k = 1, 200
      call th(real(1.0d-3 * (nextu() - 0.5d0)))
   enddo
   call th(0.0)
   call th(-0.0)
   do k = 1, 8
      call tnear(tb(k))
      call tnear(-tb(k))
   enddo
   do k = 1, 200
      call th(real(8.5d0 + 1.0d0 * nextu()))
      call th(real(-8.5d0 - 1.0d0 * nextu()))
   enddo
   ! the end marker: a truncated output (a failed self-check, a crash) cannot match without it
   write(*, '(a,1x,i0)') 'END ldpc', nrun
contains
   real function pb(k)
      integer, intent(in) :: k
      real :: b(7)
      b = (/ 0.664, 0.9217, 0.9951, 0.9998, 1.0, 0.0, 0.83 /)
      pb = b(k)
   end function pb
   real function tb(k)
      integer, intent(in) :: k
      real :: b(8)
      b = (/ 2.0**(-28), 2.0**(-12), 0.5, 1.0, 9.0107, 9.0108, 22.0, 30.0 /)
      tb = b(k)
   end function tb
   subroutine pl(x)
      real, intent(in) :: x
      real :: y
      call platanh(x, y)
      write(*, '(a,z8.8,1x,z8.8)') 'P ', transfer(x, 0), transfer(y, 0)
   end subroutine pl
   subroutine pnear(x)
      real, intent(in) :: x
      real :: v
      integer :: m
      v = x
      do m = 1, 2
         v = nearest(v, -1.0)
      enddo
      do m = 1, 5
         call pl(v)
         v = nearest(v, 1.0)
      enddo
   end subroutine pnear
   subroutine th(x)
      real, intent(in) :: x
      real :: y
      y = tanh(x)
      write(*, '(a,z8.8,1x,z8.8)') 'T ', transfer(x, 0), transfer(y, 0)
   end subroutine th
   subroutine tnear(x)
      real, intent(in) :: x
      real :: v
      integer :: m
      v = x
      do m = 1, 2
         v = nearest(v, -1.0)
      enddo
      do m = 1, 5
         call th(v)
         v = nearest(v, 1.0)
      enddo
   end subroutine tnear
   real(8) function nextu()
      state = mod(state * 1664525_8 + 1013904223_8, 4294967296_8)
      nextu = real(state, 8) / 4294967296.0d0
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
end program msk144_ldpc_oracle

subroutine msk144_trace_print(iter, zn, tov)
   implicit none
   integer, intent(in) :: iter
   real, intent(in) :: zn(128), tov(4,128)
   real :: s1, s2
   integer :: i, j, kk
   s1 = 0.0
   do i = 1, 128
      s1 = s1 + zn(i)
   enddo
   s2 = 0.0
   do j = 1, 128
      do kk = 1, 3
         s2 = s2 + tov(kk, j)
      enddo
   enddo
   write(*, '(a,i0,7(1x,z8.8))') 'I ', iter, transfer(s1, 0), transfer(s2, 0), transfer(zn(1), 0), &
        transfer(zn(65), 0), transfer(zn(128), 0), transfer(tov(1,1), 0), transfer(tov(2,64), 0)
   write(*, '(a,z8.8)') 'J ', transfer(tov(3,128), 0)
end subroutine msk144_trace_print
