! The trace printers of bpdecode128_90_trace (the reference's
! lib/bpdecode128_90.f90 with two calls added by scripts/msk144-piece-check.sh),
! shared by the ldpc and frame oracles; scripts/msk144-bp-trace.h prints the same
! from the port's hooks (src/msk144_ldpc.h).
!
!   I iter zn-sum tov-sum zn(1) zn(65) zn(128) tov(1,1) tov(2,64)
!   J tov(3,128)            checksums, at the top of each iteration
!   BZ iter zn(1:128)       every bit's zn, at the top of each iteration
!   BV iter tov(1:3,1:128)  every message from a check to a bit (kk fastest),
!                           there too
!   BO iter toc(1:nrw(j),j), j = 1..38
!                           every message from a bit to a check, and
!   BH iter tanhtoc(1:nrw(j),j), j = 1..38
!                           its tanh(-toc/2), after the tanh, in each iteration
!                           whose tov is read again (iter < maxiterations: the
!                           last iteration's messages are never read)
! All hex floats (their bits); rows past nrw(j) are never read and not printed.
!
! License: GNU GPL v3 (as WSJT-X)
subroutine msk144_trace_print(iter, zn, tov)
   implicit none
   integer, intent(in) :: iter
   real, intent(in) :: zn(128), tov(4,128)
   real :: s1, s2
   integer :: i, j, kk, ib(384)
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
   write(*, '(a,i0,1x,128z8.8)') 'BZ ', iter, transfer(zn, ib(1:128))
   do j = 1, 128
      do kk = 1, 3
         ib(3*(j-1)+kk) = transfer(tov(kk, j), 0)
      enddo
   enddo
   write(*, '(a,i0,1x,384z8.8)') 'BV ', iter, ib
end subroutine msk144_trace_print

subroutine msk144_trace_msg(iter, toc, tanhtoc, nrw)
   implicit none
   integer, intent(in) :: iter, nrw(38)
   real, intent(in) :: toc(11,38), tanhtoc(11,38)
   integer :: i, j, n, ib(418), ih(418)
   n = 0
   do j = 1, 38
      do i = 1, nrw(j)
         n = n + 1
         ib(n) = transfer(toc(i, j), 0)
         ih(n) = transfer(tanhtoc(i, j), 0)
      enddo
   enddo
   write(*, '(a,i0,1x,i0,1x,*(z8.8))') 'BO ', iter, n, ib(1:n)
   write(*, '(a,i0,1x,i0,1x,*(z8.8))') 'BH ', iter, n, ih(1:n)
end subroutine msk144_trace_msg
