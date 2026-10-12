! WSJT-X's packjt77 on FST4W payloads, for scripts/fst4w-oracle.sh: reads
! ops from stdin and prints what the reference does, for src/fst4w_unpack.cpp
! to match (tests/test_fst4w_unpack.cpp). Built against WSJT-X's own library
! (libwsjt_fort.a), so the module state is packjt77's.
!
!   I MYCALL HISCALL   fst4_decode's setup with those calls (- for blank):
!                      mycall13/dxcall13 set and, when either changed, the AP
!                      setup's pack77/unpack77 of "MYCALL HISCALL RR73"
!   U BITS             unpack77 of the 50 bits with fst4_decode's tail;
!                      prints "1|msg|" (msg the full character*37) or "0"
!   D                  the 22-bit table: "D nzhash", then "i ihash22 |call|"
!
! After each op, when nzhash or any ihash22/calls22 entry changed, prints
! "H nzhash ihash22(1) |calls22(1)|".
!
! License: GNU GPL v3 (as WSJT-X)
program fst4w_unpack_oracle
   use packjt77
   implicit none
   character(len=200) :: line
   character(len=12) :: mycall, hiscall, mycall0, hiscall0, hiscall0b
   character(len=37) :: msg, msgsent
   character(len=77) :: c77
   character(len=13) :: prev_calls(MAXHASH)
   integer :: prev_ihash(MAXHASH), prev_nz
   integer :: i, i3, n3, l1, ios, k1, k2
   logical :: ok, unpk77_success

   mycall0 = ''
   hiscall0 = ''
   prev_calls = calls22
   prev_ihash = ihash22
   prev_nz = nzhash
   do
      read(*, '(a)', iostat=ios) line
      if (ios .ne. 0) exit
      if (line(1:1) .eq. 'I') then
         k1 = 3
         k2 = index(line(k1:), ' ') + k1 - 1
         mycall = line(k1:k2-1)
         hiscall = adjustl(line(k2+1:))
         if (mycall .eq. '-') mycall = ''
         if (hiscall .eq. '-') hiscall = ''
         ! fst4_decode, from "dxcall13=hiscall" to label 10.
         dxcall13 = hiscall
         mycall13 = mycall
         l1 = index(mycall, char(0))
         if (l1 .ne. 0) mycall(l1:) = " "
         l1 = index(hiscall, char(0))
         if (l1 .ne. 0) hiscall(l1:) = " "
         if (mycall .ne. mycall0 .or. hiscall .ne. hiscall0) then
            if (len(trim(mycall)) .ge. 3) then
               hiscall0b = hiscall
               if (len(trim(hiscall0b)) .lt. 3) hiscall0b = mycall
               msg = trim(mycall)//' '//trim(hiscall0b)//' RR73'
               i3 = -1
               n3 = -1
               call pack77(msg, i3, n3, c77)
               call unpack77(c77, 1, msgsent, unpk77_success)
            end if
            mycall0 = mycall
            hiscall0 = hiscall
         end if
      else if (line(1:1) .eq. 'U') then
         c77 = line(3:52)
         c77(51:77) = '000000000000000000000110000'
         call unpack77(c77, 1, msg, ok)
         if (ok) then
            write(*, '(a)') '1|'//msg//'|'
         else
            write(*, '(a)') '0'
         end if
      else if (line(1:1) .eq. 'D') then
         write(*, '(a,i0)') 'D ', nzhash
         do i = 1, nzhash
            write(*, '(i0,1x,i0,1x,a)') i, ihash22(i), '|'//calls22(i)//'|'
         end do
      end if
      if (nzhash .ne. prev_nz .or. any(ihash22 .ne. prev_ihash) .or. any(calls22 .ne. prev_calls)) then
         write(*, '(a,i0,1x,i0,1x,a)') 'H ', nzhash, ihash22(1), '|'//calls22(1)//'|'
         prev_calls = calls22
         prev_ihash = ihash22
         prev_nz = nzhash
      end if
   end do
end program fst4w_unpack_oracle
