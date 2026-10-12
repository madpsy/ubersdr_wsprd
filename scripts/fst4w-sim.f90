! Synthetic FST4W inputs for scripts/fst4w-oracle.sh: several signals per
! file, FST4W and WSPR-2 mixed (on a real band FST4W's window is WSPR's), in
! seeded noise, as a 12 kHz 16-bit WAV of the period's length. Built against
! WSJT-X's own library (genfst4, gen_fst4wave, watterson, lorentzian_fading,
! genwspr, gran) the way its fst4sim is, so signals, fading and noise are
! fst4sim's: unit-variance noise at 12 kHz, a signal of S dB in 2500 Hz an
! amplitude of sqrt(2*2500/6000)*10^(S/20), all times 100 into 16 bits.
!
! Copyright 2020 Steven Franke, K9AN and Joe Taylor, K1JT (fst4sim)
! License: GNU GPL v3
!
!   fst4w-sim -p PERIOD -S SEED [-N] -o FILE.wav [SIGNAL ...]
!
! SIGNAL is one of
!   F:MESSAGE@f,dt,fspread,snr[,delay]  FST4W-PERIOD as fst4sim sends it: f
!       the lowest tone's audio Hz (what jt9 reports), dt from the nominal
!       start (1 s; 0.5 s for 15 s periods), fspread > 0 Watterson Doppler
!       spread (Hz, with delay ms), < 0 Lorentzian fading of -fspread Hz, 0
!       none; snr dB in 2500 Hz
!   W:MESSAGE@f,t0,drift,snr  WSPR-2: f the middle of its four tones (Hz),
!       t0 its start in seconds into the file, drift Hz over the
!       transmission; snr dB in 2500 Hz
! -N: no noise. The noise and the fading come from the C library's rand(),
! seeded with SEED before anything is drawn, in the order the signals are
! given: the same arguments give the same file with the same libc.
program fst4w_sim
   use wavhdr
   use, intrinsic :: iso_c_binding, only: c_int
   include 'fst4/fst4_params.f90'
   interface
      subroutine c_srand(seed) bind(C, name='srand')
         import c_int
         integer(c_int), value :: seed
      end subroutine c_srand
   end interface
   type(hdr) h
   character(len=512) :: arg, outfile, spec
   character*37 msg37, msgsent37
   character*22 msg22, msgsent22
   complex, allocatable :: c0(:)
   real, allocatable :: wave(:), total(:)
   integer*2, allocatable :: iwave(:)
   integer*1 msgbits(101)
   integer itone(NN), itw(162)
   integer hmod, nsec, nseed, nargs, iarg, i, k, kat, nsps, iwspr, nmax, nwave, nz, ios
   logical noise
   real f00, xdt, fspread, delay, snrdb, sig, f0, fs, dt, baud, t0, drift
   real*8 phi, dphi, twopi8, fsym

   nsec = 0
   nseed = 1
   noise = .true.
   outfile = ''
   nargs = command_argument_count()
   iarg = 1
   do while (iarg .le. nargs)
      call get_command_argument(iarg, arg)
      if (arg .eq. '-p') then
         iarg = iarg + 1; call get_command_argument(iarg, arg); read(arg,*) nsec
      else if (arg .eq. '-S') then
         iarg = iarg + 1; call get_command_argument(iarg, arg); read(arg,*) nseed
      else if (arg .eq. '-o') then
         iarg = iarg + 1; call get_command_argument(iarg, outfile)
      else if (arg .eq. '-N') then
         noise = .false.
      else
         exit
      end if
      iarg = iarg + 1
   end do
   nsps = 0
   if (nsec .eq. 15) nsps = 720
   if (nsec .eq. 30) nsps = 1680
   if (nsec .eq. 60) nsps = 3888
   if (nsec .eq. 120) nsps = 8200
   if (nsec .eq. 300) nsps = 21504
   if (nsec .eq. 900) nsps = 66560
   if (nsec .eq. 1800) nsps = 134400
   if (nsps .eq. 0 .or. len_trim(outfile) .eq. 0) then
      write(0,*) 'usage: fst4w-sim -p 15|30|60|120|300|900|1800 -S SEED [-N] -o FILE.wav ', &
         '"F:MESSAGE@f,dt,fspread,snr[,delay]" | "W:MESSAGE@f,t0,drift,snr" ...'
      stop 2
   end if

   call c_srand(int(nseed, c_int))
   fs = 12000.0
   dt = 1.0 / fs
   twopi8 = 8.d0 * atan(1.d0)
   nmax = nsec * 12000
   nz = nsps * NN
   nwave = max(nmax, (NN + 2) * nsps)
   baud = fs / nsps
   hmod = 1
   allocate(c0(0:nwave-1), wave(nwave), total(nmax), iwave(nmax))
   total = 0.0

   do while (iarg .le. nargs)
      call get_command_argument(iarg, spec)
      iarg = iarg + 1
      kat = index(spec, '@', back=.true.)
      if (kat .eq. 0 .or. spec(2:2) .ne. ':') call bad(spec)
      if (spec(1:1) .eq. 'F') then
         delay = 0.0
         read(spec(kat+1:), *, iostat=ios) f00, xdt, fspread, snrdb, delay
         if (ios .ne. 0) then
            delay = 0.0
            read(spec(kat+1:), *, iostat=ios) f00, xdt, fspread, snrdb
            if (ios .ne. 0) call bad(spec)
         end if
         msg37 = spec(3:kat-1)
         if (msg37(1:1) .eq. '<') then
            ! A type 3 ("<CALL> GRID6") packs only once its call's hash is
            ! known: let the encoder see the call first.
            iwspr = 1
            call genfst4(msg37(2:index(msg37, '>')-1)//' 30', 1, msgsent37, msgbits, itone, iwspr)
         end if
         iwspr = 1   ! genfst4 writes it back: (240,74) or not
         call genfst4(msg37, 0, msgsent37, msgbits, itone, iwspr)
         if (iwspr .ne. 1) call bad(spec)
         if (trim(msgsent37) .eq. '*** bad message ***') call bad(spec)
         f0 = f00 + 1.5 * hmod * baud
         call gen_fst4wave(itone, NN, nsps, nwave, fs, hmod, f0, 1, c0, wave)
         k = nint((xdt + 1.0) / dt)
         if (nsec .eq. 15) k = nint((xdt + 0.5) / dt)
         c0 = cshift(c0, -k)
         if (k .gt. 0) c0(0:k-1) = 0.0
         if (k .lt. 0) c0(nmax+k:nmax-1) = 0.0
         if (fspread .gt. 0.0 .or. delay .ne. 0.0) call watterson(c0, nwave, nz, fs, delay, fspread)
         if (fspread .lt. 0.0) call lorentzian_fading(c0, nwave, fs, -fspread)
         sig = sqrt(2 * 2500.0 / 6000.0) * 10.0**(0.05 * snrdb)
         total = total + sig * aimag(c0(0:nmax-1))
         write(*,'(a,a,a,a)') 'F ', trim(msgsent37), '  ', trim(spec(kat+1:))
      else if (spec(1:1) .eq. 'W') then
         read(spec(kat+1:), *, iostat=ios) f00, t0, drift, snrdb
         if (ios .ne. 0) call bad(spec)
         msg22 = spec(3:kat-1)
         call genwspr(msg22, msgsent22, itw)
         sig = sqrt(2 * 2500.0 / 6000.0) * 10.0**(0.05 * snrdb)
         phi = 0.d0
         do i = 0, 162 * 8192 - 1
            k = nint(t0 * fs) + i
            fsym = f00 + drift / 2.0 * (i / 8192 - 81) / 81.0 + (itw(i/8192 + 1) - 1.5d0) * fs / 8192.d0
            dphi = twopi8 * fsym / fs
            if (k .ge. 0 .and. k .lt. nmax) total(k+1) = total(k+1) + sig * sin(phi)
            phi = modulo(phi + dphi, twopi8)
         end do
         write(*,'(a,a,a,a)') 'W ', trim(msgsent22), '  ', trim(spec(kat+1:))
      else
         call bad(spec)
      end if
   end do

   if (noise) then
      do i = 1, nmax
         total(i) = total(i) + gran()
      end do
   end if
   total = 100.0 * total
   if (any(abs(total) .gt. 32767.0)) then
      write(0,*) 'fst4w-sim: signals too strong, the file would clip'
      stop 1
   end if
   iwave = nint(total)
   h = default_header(12000, nmax)
   open(10, file=trim(outfile), status='replace', access='stream')
   write(10) h, iwave
   close(10)

contains

   subroutine bad(s)
      character(len=*), intent(in) :: s
      write(0,*) 'fst4w-sim: bad signal: ', trim(s)
      stop 2
   end subroutine bad

end program fst4w_sim
