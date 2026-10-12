! WSJT-X's (240,74) LDPC decoding on soft-bit vectors, for
! scripts/fst4w-oracle.sh: reads ops from stdin and prints what the
! reference does, for src/fst4w_ldpc.cpp to match
! (tests/test_fst4w_ldpc.cpp). Built against WSJT-X's own library
! (libwsjt_fort.a): decode240_74, fastosd240_74, encode240_74, get_crc24.
!
!   L HEX*240          the soft bits (float bits, 8 hex digits each)
!   V KEFF MAXOSD NORD decode240_74 with apmask 0 and dmin 0.0 first, as
!                      fst4_decode calls it; message74 and cw carry over
!                      from op to op as fst4_decode's own do. Prints
!                      "R ntype nharderror dmin(hex) message74 cw"
!   O K NDEEP          fastosd240_74 on the soft bits (message74, cw carried
!                      over): "S nhardmin dmin(hex) message74 cw"
!   C BITS74           get_crc24: "C ncrc"
!   E BITS74           encode240_74: "E codeword"
!
! License: GNU GPL v3 (as WSJT-X)
program fst4w_ldpc_oracle
   implicit none
   character(len=2048) :: line
   real :: llr(240), dmin
   integer :: ibits(240), ios, keff, maxosd, norder, k, ndeep, ntype, nharderror, nhardmin, ncrc
   integer*1 :: apmask(240), message74(74), cw(240), m74(74)

   apmask = 0
   message74 = 0
   cw = 0
   do
      read(*, '(a)', iostat=ios) line
      if (ios .ne. 0) exit
      select case (line(1:1))
      case ('L')
         read(line(3:), '(240z8)') ibits
         llr = transfer(ibits, llr)
      case ('V')
         read(line(3:), *) keff, maxosd, norder
         dmin = 0.0
         nharderror = -1
         call decode240_74(llr, keff, maxosd, norder, apmask, message74, cw, ntype, nharderror, dmin)
         write(*, '(a,i0,1x,i0,1x,z8.8,1x,74i1,1x,240i1)') 'R ', ntype, nharderror, transfer(dmin, 0), message74, cw
      case ('O')
         read(line(3:), *) k, ndeep
         dmin = 0.0
         call fastosd240_74(llr, k, apmask, ndeep, message74, cw, nhardmin, dmin)
         write(*, '(a,i0,1x,z8.8,1x,74i1,1x,240i1)') 'S ', nhardmin, transfer(dmin, 0), message74, cw
      case ('C')
         read(line(3:76), '(74i1)') m74
         call get_crc24(m74, 74, ncrc)
         write(*, '(a,i0)') 'C ', ncrc
      case ('E')
         read(line(3:76), '(74i1)') m74
         call encode240_74(m74, cw)
         write(*, '(a,240i1)') 'E ', cw
      end select
   end do
end program fst4w_ldpc_oracle
