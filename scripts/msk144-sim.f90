! Synthetic MSK144 and MSK40 inputs for scripts/msk144-oracle.sh: WSJT-X's
! msk144sim with the choices it fixes made arguments, so one file can hold
! one ping, three pings, a ping of any width, or a carrier for the whole slot,
! with seeded noise. Built against WSJT-X's own library (genmsk_128_90, and
! through it genmsk40; wavhdr; gran), the way msk144sim is linked. The signal,
! the ping shape and the noise scaling are msk144sim's and makepings':
!   - the ping amplitude sqrt(2)*10^(snr/20) at its peak, a peak at the ping's
!     start plus one width: 2.718*t*exp(-t) with t = (time - start)/width,
!     zero before the start and beyond ten widths;
!   - noise of variance 6000/2500 a sample, the sum times 30 into 16 bits.
! msk144sim sends a ping every second from t = 1 s; here the list is given.
!
! Copyright 2016 Steven Franke, K9AN and Joe Taylor, K1JT (msk144sim)
! License: GNU GPL v3
!
!   msk144-sim -p TR -S SEED [-f HZ] [-s SNR] [-w WIDTH|slot] [-t T1,T2,..]
!              [-N | -Z] -o FILE.wav "MESSAGE"
!
!   -p TR     the file's length in seconds (default 15)
!   -S SEED   seeds the C library's rand(), which gran() draws on; the same
!             arguments give the same file with the same libc (default 1)
!   -f HZ     the centre of the two tones, as msk144sim's freq (default 1500)
!   -s SNR    dB in 2500 Hz (default 0)
!   -w WIDTH  the ping's width in seconds (default 0.2); "slot": no ping, the
!             signal at its peak amplitude for the whole file
!   -t LIST   the pings' start times in seconds, comma separated (default 1)
!   -N        the signal alone, no noise
!   -Z        digital silence: all zeros, whatever else is given
!   MESSAGE   as msk144sim takes it; "<K1ABC W9XYZ> +03" makes an MSK40 frame
!             (genmsk40: the short message, 40 symbols repeated). Without a
!             MESSAGE (or with -n), only noise.
program msk144_sim
   use wavhdr
   use, intrinsic :: iso_c_binding, only: c_int
   interface
      subroutine c_srand(seed) bind(C, name='srand')
         import c_int
         integer(c_int), value :: seed
      end subroutine c_srand
   end interface
   integer, parameter :: NMAX = 1800*12000
   type(hdr) h
   character(len=512) :: arg, outfile, tlist
   character*37 msg, msgsent
   integer*4 itone(144)
   integer*2, allocatable :: iwave(:)
   real, allocatable :: pings(:), waveform(:)
   real*8 twopi, freq, phi, dphi0, dphi1, dphi, t, tstart(64), baud
   real :: snr, width, sig, fac, xx, wave
   integer :: nseed, ntr, npts, nping, i, j, k, nsym, nsps, itype, nrep, jrep, ip, iarg, ios, p
   logical :: nonoise, silence, noiseonly, wholeslot

   nseed = 1; ntr = 15; freq = 1500.d0; snr = 0.0; width = 0.2
   tlist = '1'; outfile = ''; msg = ' '
   nonoise = .false.; silence = .false.; noiseonly = .true.; wholeslot = .false.
   iarg = 1
   do while (iarg <= command_argument_count())
      call get_command_argument(iarg, arg)
      select case (trim(arg))
      case ('-N'); nonoise = .true.
      case ('-Z'); silence = .true.
      case ('-n'); noiseonly = .true.
      case ('-p', '-S', '-f', '-s', '-w', '-t', '-o')
         iarg = iarg + 1
         if (iarg > command_argument_count()) call usage()
         select case (trim(arg))
         case ('-p'); call get_command_argument(iarg, arg); read(arg, *, iostat=ios) ntr
         case ('-S'); call get_command_argument(iarg, arg); read(arg, *, iostat=ios) nseed
         case ('-f'); call get_command_argument(iarg, arg); read(arg, *, iostat=ios) freq
         case ('-s'); call get_command_argument(iarg, arg); read(arg, *, iostat=ios) snr
         case ('-w')
            call get_command_argument(iarg, arg)
            if (trim(arg) == 'slot') then
               wholeslot = .true.
            else
               read(arg, *, iostat=ios) width
            end if
         case ('-t'); call get_command_argument(iarg, tlist)
         case ('-o'); call get_command_argument(iarg, outfile)
         end select
         if (ios /= 0) call usage()
      case default
         msg = trim(arg); noiseonly = .false.
      end select
      iarg = iarg + 1
   end do
   if (len_trim(outfile) == 0 .or. ntr < 1 .or. ntr*12000 > NMAX) call usage()
   if (msg == ' ') noiseonly = .true.

   nping = 0
   p = 1
   do while (p <= len_trim(tlist))
      i = index(tlist(p:), ',')
      if (i == 0) i = len_trim(tlist) - p + 2
      nping = nping + 1
      read(tlist(p:p+i-2), *) tstart(nping)
      p = p + i
   end do

   npts = ntr*12000
   allocate(iwave(0:npts-1), pings(0:npts-1), waveform(0:npts-1))
   iwave = 0; pings = 0.0; waveform = 0.0
   sig = sqrt(2.0)*10.0**(0.05*snr)

   if (.not. noiseonly .and. .not. silence) then
      call genmsk_128_90(msg, 0, msgsent, itone, itype)
      if (itype < 1 .or. itype > 7) then
         write(*, *) 'Illegal message'
         stop 1
      end if
      write(*, '(a,a)') 'Message sent: ', trim(msgsent)
      twopi = 8.d0*atan(1.d0)
      nsym = 144; nsps = 6; baud = 2000.d0
      if (itone(41) < 0) nsym = 40
      dphi0 = twopi*(freq - 0.25d0*baud)/12000.d0
      dphi1 = twopi*(freq + 0.25d0*baud)/12000.d0
      phi = 0.0; k = 0; dphi = dphi0
      nrep = npts/(nsym*nsps)
      do jrep = 1, nrep
         do i = 1, nsym
            dphi = dphi0
            if (itone(i) /= 0) dphi = dphi1
            do j = 1, nsps
               waveform(k) = cos(phi)
               k = k + 1
               phi = mod(phi + dphi, twopi)
            end do
         end do
      end do
      ! the tail the whole repeats do not reach continues the last phase
      do while (k < npts)
         waveform(k) = cos(phi)
         k = k + 1
         phi = mod(phi + dphi, twopi)
      end do
      if (wholeslot) then
         pings = sig
      else
         do i = 0, npts - 1
            do ip = 1, nping
               t = ((i + 1)/12000.d0 - tstart(ip))/width
               if (t >= 0.d0 .and. t <= 10.d0) pings(i) = max(pings(i), real(2.718*t*dexp(-t)*sig))
            end do
         end do
      end if
   end if

   if (.not. silence) then
      call c_srand(int(nseed, c_int))
      fac = sqrt(6000.0/2500.0)
      do i = 0, npts - 1
         xx = 0.0
         if (.not. nonoise) xx = gran()
         wave = pings(i)*waveform(i) + fac*xx
         iwave(i) = int(30.0*wave, kind=2)
      end do
   end if

   h = default_header(12000, npts)
   open(10, file=trim(outfile), status='replace', access='stream')
   write(10) h, iwave
   close(10)

contains
   subroutine usage()
      write(*, '(a)') 'Usage: msk144-sim -p TR -S SEED [-f HZ] [-s SNR] [-w WIDTH|slot] [-t T1,T2,..] [-N|-Z] -o FILE.wav "MESSAGE"'
      stop 2
   end subroutine usage
end program msk144_sim
