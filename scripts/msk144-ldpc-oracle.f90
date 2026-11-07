! WSJT-X's (128,90) LDPC decoder, CRC-13 check, platanh and tanh, for
! scripts/msk144-piece-check.sh ldpc: generates its own fixed inputs and
! prints them and what the reference does with them, for
! scripts/msk144-ldpc-check.cpp (src/msk144_ldpc.cpp) to match line for line.
! Built against WSJT-X's own library (libwsjt_fort.a): encode_128_90,
! bpdecode128_90, chkcrc13a, platanh. The piece check appends
! bpdecode128_90_trace to this file: lib/bpdecode128_90.f90 with two calls
! added in each iteration, and their printers, scripts/msk144-bp-trace.f90
! (the I, J, BZ, BV, BO and BH lines; see scripts/msk144-piece-check.sh).
!
!   L MAXITER TRACE HEX*128   the soft bits (float bits, 8 hex digits each)
!   I, J, BZ, BV, BO, BH      the trace (scripts/msk144-bp-trace.f90), for
!                             the vectors with TRACE 1 (all of them)
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
! maxiter 30; then 40 of pure noise; then 20 at the edges of the float range
! (vectors 241-260: codewords' soft bits scaled down to subnormals and up to
! 1e37, all zeros, signed zeros, +inf and -inf among noise, small integers,
! noise at 0.5); every vector traced;
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
   do iv = 1, 260
      if (iv .gt. 240) then
         call edge_llr(iv - 240)
      else if (iv .le. 200) then
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
      ltrace = .true.
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
   ! the edge cases: k = 1..20
   subroutine edge_llr(k)
      integer, intent(in) :: k
      real, parameter :: sc(9) = (/ 1.0e-40, 1.0e-38, 1.0e-30, 1.0e-6, 1.0e3, 1.0e6, 1.0e20, 1.0e30, 1.0e37 /)
      real :: pinf
      integer :: ii
      pinf = transfer(int(z'7F800000'), 1.0)
      if (k .le. 9) then
         do ii = 1, 128
            llr(ii) = sc(k) * (real(2 * enc(k, ii) - 1) + 0.8 * gauss())
         enddo
      else if (k .eq. 10) then
         llr = 0.0
      else if (k .eq. 11) then
         do ii = 1, 128
            llr(ii) = 0.0
            if (mod(ii, 2) .eq. 0) llr(ii) = -0.0
         enddo
      else if (k .le. 13) then
         do ii = 1, 128
            llr(ii) = 4.0 * (real(2 * enc(k, ii) - 1) + 0.6 * gauss())
         enddo
         llr(3) = pinf
         llr(70) = pinf
         if (k .eq. 13) then
            llr(3) = -pinf
            llr(70) = -pinf
            llr(128) = -pinf
         endif
      else if (k .le. 16) then
         do ii = 1, 128
            llr(ii) = real(nint(3.0 * (real(2 * enc(k, ii) - 1) + 0.9 * gauss())))
         enddo
      else
         do ii = 1, 128
            llr(ii) = 0.5 * gauss()
         enddo
      endif
   end subroutine edge_llr
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
