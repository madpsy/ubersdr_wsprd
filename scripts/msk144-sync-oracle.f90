! WSJT-X's MSK144 short-ping detector and sync for scripts/msk144-piece-check.sh
! sync: makes its own fixed inputs, prints them and what the reference does
! with them, for scripts/msk144-sync-check.cpp (src/msk144_sync.cpp) to match
! line for line. Built against WSJT-X's own library (libwsjt_fort.a): analytic,
! msk144spd, msk144sync, msk144decodeframe, msk_spectrum (msk_clip_spectrum_window,
! msk_peak_offset), indexx. The piece check appends two instrumented copies of
! the reference's routines to this file (lib/msk144spd.f90 as msk144spd_trace,
! lib/msk144sync.f90 as msk144sync_trace: sed adds calls, see
! scripts/msk144-piece-check.sh) that print their intermediate values; every
! case runs the copy and the library's routine and stops if their results
! differ. The copy of msk144spd calls msk144sync_trace and, for each frame it
! decodes, msk144_frame_hook, which runs the library's msk144decodeframe on the
! frame itself (it corrects the frame in place, as it does in the library's
! msk144spd, so the copy's ct, softbits and msgreceived are the library's) and
! prints a hash of what it was given and what it answered: the port's driver
! runs its own msk144decodeframe (src/msk144_decode.cpp) on the frames its own
! msk144spd cuts, and prints the same.
!
! Every float is its bits, 8 hex digits (a real*8, 16); arrays 16 values a
! line, "TAG k first-index values", complex values as re, im pairs.
!
! Inputs
!   I s i v...                  slice s: 7168 points of the analytic signal of a
!                               block of a recording (mskrtd's: the int16 block
!                               over its rms, analytic, nfft 8192, equalizer on,
!                               pc 0); slices 1..: per file, up to three blocks
!                               that msk144spd decodes, the hottest, the quietest.
!                               Slice 99: 28000 points, blocks end to end.
! Direct checks of the pieces
!   L k kind n lo hi res v...   maxloc: kind 0 maxloc(a(1:n)), kind 1
!                               maxloc(a, mask) with the mask true on lo..hi;
!                               res the index (0: none)
!   J k n v... | idx...         indexx(v, n, idx)
!   K k nb peak res re im ...   msk_peak_offset(spectrum(1:nb), peak) = res
!   W k nfft lo hi res(lo hi u) msk_clip_spectrum_window(nfft, lo, hi, ...) with
!                               the requested bounds real*8 bits
! Sync and short-ping cases
!   SYN k s off nf ntol delf mask(8) npeaks fc nanidx infidx scale
!                               msk144sync on slice s from point off*864+1,
!                               nf frames, as msk144sync's arguments (nanidx
!                               and infidx: a point made NaN or +Inf, 1-based,
!                               0 none (infidx < 0: the real part of point
!                               -infidx made 1e20); scale: every point times it, first;
!                               s 0: zeros)
!   SPD k s n ntol fc nanidx infidx scale
!                               msk144spd on the first n points of slice s
!   P pp(1:12), B cb(1:42)      msk144sync's constants, at its first call
!   R rcw(1:12)                 msk144spd's, at its first call
!   N ihlo ihhi illo ilhi i2000 i4000   spd's windows
!   T istp ihpk ilpk deltah deltal ah al trath tratl ferr   one per step
!   D k i v...                  spd's detmet(1:nstep) as the steps made it
!   H k i v..., O k i v...      detmet2(1:nstep), detfer(1:nstep)
!   M k nstep xmed idx...       the median and indexx's indices
!   E k 0 ndet / E k j nstart ferr snr   the candidates
!   Q k nfreqs xmax fest        sync, before its peak search
!   Y k i v... / Z k i v...     sync's c(1:864), xcc(0:863) there
!   G k icand iav ns fest xmax npkloc(1:10)   spd's call of sync
!   F k j h1 h2 ns |msg|        a call of the frame decoder (ct as two
!                               hashes of its bits)
!   X k ns fest xmax npklocs    msk144sync's results
!   X k ns |msg| fret tret navg msk144spd's; with U k i v... (softbits) and
!                               V k i v... (ct) as it left them
!
! License: GNU GPL v3 (as WSJT-X)
program msk144_sync_oracle
   use msk_spectrum
   use packjt77, only: calls10, calls12, calls22, recent_calls, ihash22, nzhash
   implicit none
   integer, parameter :: NZ = 7168, NSPM = 864, NF = 8192, NSL = 40, NLONG = 28000
   character(len=512) :: repo, path(4)
   integer :: k, i, j, ifile, iblk, nsamp, ios, istart, ibest, iworst, nslice, ndec, nz6912
   integer :: ns1, nsl_file(4), nsdec, sfirst, sdec1, shot1
   integer*2, allocatable :: wav(:)
   real :: d(NF), rms, e, ebest, eworst
   real*8 :: pc(5)
   logical*1 :: beq
   complex :: cd(NF), ct(NSPM)
   complex, allocatable :: sl(:,:), cl(:)
   real :: sb(144), fret, tret
   integer :: navg, nsuc
   character(len=37) :: msg
   integer :: kc, nframe
   common /skcase/ kc, nframe
   integer :: npatience, nthreads
   common /patience/ npatience, nthreads
   real :: nanv, infv
   integer :: ihot(4), iquiet(4), idec(3,4), ndecf(4), s, mask8(8)
   real :: sc, fcx
   integer :: s1x
   integer :: patt(8,4)
   integer*8 :: sst
   character(len=77) :: synbits(3)
   real, parameter :: synsig(3) = (/0.15, 0.3, 0.6/)
   integer, parameter :: synstart(3) = (/201, 437, 650/)
   data patt/1,1,1,1,0,0,0,0, 0,0,1,1,1,1,0,0, 1,1,1,1,1,0,0,0, 1,1,1,1,1,1,1,0/

   npatience = 1
   nthreads = 1
   sst = 4242_8
   synbits(1) = '11100000001100101110001110101111110010011011000110000111000100100010001001100'
   synbits(2) = '11001000001000001110101001101101111001100010101001100011001000110111011101001'
   synbits(3) = '00001001101111011110001101010000011000010100100111011100000111111001110101001'
   nanv = transfer(int(z'7FC00000'), 0.0)
   infv = transfer(int(z'7F800000'), 0.0)
   call get_command_argument(1, repo)
   path(1) = trim(repo)//'/recordings/msk144-samples/181211_120500.wav'
   path(2) = trim(repo)//'/recordings/msk144-samples/181211_120800.wav'
   path(3) = trim(repo)//'/recordings/msk144-sim/snr-p10-261106_120010.wav'
   path(4) = trim(repo)//'/recordings/msk144-sim/freq-p200-sp04-261106_120019.wav'
   allocate(sl(NZ, NSL), cl(NLONG))
   sl = 0.0
   cl = 0.0
   kc = 0
   nframe = 0
   k = 0

   ! ---- the direct checks
   call direct_checks(k)

   ! ---- slices
   nslice = 0
   do ifile = 1, 4
      open(10, file=trim(path(ifile)), access='stream', form='unformatted', status='old', iostat=ios)
      if (ios /= 0) error stop 'oracle: cannot open a WAV'
      inquire(10, size=nsamp)
      nsamp = (nsamp - 44)/2
      allocate(wav(nsamp))
      read(10, pos=45) wav
      close(10)
      ebest = -1.0; eworst = huge(1.0)
      ibest = 0; iworst = 0
      ndec = 0
      istart = 0
      do while (istart + NZ <= nsamp)
         call block_analytic(wav, istart, cd)
         ! the energy of the block, as the frontend's check ranks them
         e = 0.0
         do i = 1, NZ
            e = e + real(wav(istart+i))**2
         enddo
         if (e > ebest) then
            ebest = e; ibest = istart
         endif
         if (e < eworst) then
            eworst = e; iworst = istart
         endif
         if (ndec < 3) then
            kc = 0
            call quiet_spd(cd, 6912, 200, 1500.0, nsuc)
            if (nsuc > 0) then
               ndec = ndec + 1
               nslice = nslice + 1
               call add_slice(nslice, cd)
            endif
         endif
         istart = istart + NSPM*2
      enddo
      nslice = nslice + 1
      call block_analytic(wav, ibest, cd)
      call add_slice(nslice, cd)
      nslice = nslice + 1
      call block_analytic(wav, iworst, cd)
      call add_slice(nslice, cd)
      nsl_file(ifile) = nslice
      ndecf(ifile) = ndec
      if (ifile == 1) then
         ! slice 99: 30000 points, hottest blocks of the files end to end, filled as the files come
         cl(1:NZ) = sl(1:NZ, nslice-1)
      else if (ifile == 2) then
         cl(NZ+1:2*NZ) = sl(1:NZ, nslice-1)
      else if (ifile == 3) then
         cl(2*NZ+1:3*NZ) = sl(1:NZ, nslice-1)
      else
         cl(3*NZ+1:NLONG) = sl(1:NLONG-3*NZ, nslice-1)
      endif
      deallocate(wav)
   enddo
   call print_slice(99, cl, NLONG)
   call check_synbits()
   ! check_synbits runs the library's decoder and leaves the shared packjt77 tables (calls, recent_calls, hashes)
   ! filled: the port starts empty, so start the cases from the same empty state
   calls10 = ''; calls12 = ''; calls22 = ''; recent_calls = ''; ihash22 = -1; nzhash = 0
   ! ---- synthetic blocks: eight repeats of one frame at 1500 Hz, so that the real frame decoder meets frames
   ! it refuses. The 77 bits: one that unpacks to
   ! 'CQ <...>' (refused, with that text left in msgreceived), one a message the decoder returns no success for
   ! (text left too), and 'K1ABC W9XYZ RR73' (decoded).
   ! (message bits A, B, C; noise sigma; the first sample of the first frame)
   do i = 1, 3
      do j = 1, 3
         call synth_slice(synbits(i), synsig(j), synstart(1 + mod(i + j, 3)))
      enddo
   enddo
   call synth_slice(synbits(1), 1.2, 201)

   shot1 = nsl_file(1) - 1
   sdec1 = 1
   sfirst = 1
   sdec1 = 1
   nsdec = ndecf(1)
   ! the first slice that decodes, or the first file's hottest
   if (ndecf(1) == 0) sdec1 = shot1

   ! ---- msk144sync directly: spd's three-frame patterns on each slice's start, mskrtd's four eight-frame ones,
   ! then the odd cases (k src off nf ntol delf mask npeaks fc nan inf scale)
   do s = 1, nslice
      do i = 1, 6
         call navpat(i, mask8)
         call run_sync(k, s, 0, 3, 8, 2.0, mask8, 2, 1500.0, 0, 0, 1.0)
      enddo
      do i = 1, 4
         mask8 = patt(:, i)
         call run_sync(k, s, 0, 8, 15, 10.0/real(sum(mask8)), mask8, 2, 1500.0, 0, 0, 1.0)
      enddo
      ! mskrtd's four patterns with a small user ntol: ntol/deltaf in [0.5, 1.5) gives nint = 1, so nfreqs = 3
      ! (2*1+1) trial frequencies, fewer than the OpenMP library's four threads; ntol 2 and 3 over deltaf 2.5, 2.0
      ! and 1.43 (ratios 0.8, 1.2; 1.0, 1.5; 1.4, 2.1), ntol 1 over them (0.4, 0.5, 0.7)
      do i = 1, 4
         mask8 = patt(:, i)
         do j = 1, 3
            call run_sync(k, s, 0, 8, j, 10.0/real(sum(mask8)), mask8, 2, 1500.0, 0, 0, 1.0)
         enddo
      enddo
   enddo
   s = sdec1
   mask8 = patt(:, 1)
   call run_sync(k, s, 1, 3, 8, 2.0, (/0,1,0,0,0,0,0,0/), 2, 1500.0, 0, 0, 1.0)
   call run_sync(k, s, 2, 3, 8, 2.0, (/1,1,1,0,0,0,0,0/), 1, 1500.0, 0, 0, 1.0)
   call run_sync(k, s, 0, 3, 8, 2.0, (/1,1,1,0,0,0,0,0/), 3, 1503.7, 0, 0, 1.0)
   call run_sync(k, s, 0, 3, 8, 2.0, (/1,1,1,0,0,0,0,0/), 5, 1500.0, 0, 0, 1.0)
   call run_sync(k, s, 0, 3, 0, 2.0, (/1,1,1,0,0,0,0,0/), 2, 1500.0, 0, 0, 1.0)     ! ntol 0: one trial
   call run_sync(k, s, 0, 3, 15, 4.0, (/1,1,1,0,0,0,0,0/), 2, 1500.0, 0, 0, 1.0)    ! nint(3.75) = 4
   call run_sync(k, s, 0, 3, 3, 10.0, (/1,1,1,0,0,0,0,0/), 2, 1500.0, 0, 0, 1.0)    ! nint(0.3) = 0
   call run_sync(k, s, 0, 3, 7, 2.0, (/1,1,1,0,0,0,0,0/), 2, 1500.0, 0, 0, 1.0)     ! nint(3.5) = 4
   call run_sync(k, s, 0, 3, 8, 2.0, (/0,0,0,0,0,0,0,0/), 2, 1500.0, 0, 0, 1.0)     ! navg 0
   call run_sync(k, s, 0, 1, 8, 2.0, (/1,0,0,0,0,0,0,0/), 2, 1500.0, 0, 0, 1.0)     ! one frame
   call run_sync(k, s, 0, 2, 8, 2.0, (/1,1,0,0,0,0,0,0/), 2, 1500.0, 0, 0, 1.0)
   call run_sync(k, s, 0, 6, 8, 2.0, (/1,1,1,1,1,1,0,0/), 2, 1500.0, 0, 0, 1.0)
   call run_sync(k, s, 0, 0, 8, 2.0, (/0,0,0,0,0,0,0,0/), 2, 1500.0, 0, 0, 1.0)     ! no frames
   call run_sync(k, s, 0, 3, 8, 2.0, (/1,1,1,0,0,0,0,0/), 2, 1500.0, 1500, 0, 1.0)  ! a NaN point
   call run_sync(k, s, 0, 3, 8, 2.0, (/1,1,1,0,0,0,0,0/), 2, 1500.0, 0, 7, 1.0)     ! an Inf point
   ! xmax exactly 1.3 (0x3FA66666), where nsuccess turns on (xmax >= 1.3): the scale found by stepping the
   ! float scale (the correlation is linear in it) with the library's msk144sync, then the scales either side
   call find_scale(s1x, fcx, sc)
   call run_sync(k, s1x, 0, 3, 8, 2.0, (/1,1,1,0,0,0,0,0/), 2, fcx, 0, 0, nearest(sc, -1.0))
   call run_sync(k, s1x, 0, 3, 8, 2.0, (/1,1,1,0,0,0,0,0/), 2, fcx, 0, 0, sc)
   call run_sync(k, s1x, 0, 3, 8, 2.0, (/1,1,1,0,0,0,0,0/), 2, fcx, 0, 0, nearest(sc, 1.0))
   call run_sync(k, s, 0, 3, 8, 2.0, (/1,1,1,0,0,0,0,0/), 2, 1500.0, 0, -7, 1.0)    ! a point of 1e20
   call run_sync(k, s, 0, 3, 8, 2.0, (/1,1,1,0,0,0,0,0/), 2, 1500.0, 0, 0, 1.0e20)  ! overflow
   call run_sync(k, s, 0, 3, 8, 2.0, (/1,1,1,0,0,0,0,0/), 2, 1500.0, 0, 0, 0.0)     ! scaled to zero
   call run_sync(k, 0, 0, 3, 8, 2.0, (/1,1,1,0,0,0,0,0/), 2, 1500.0, 0, 0, 1.0)     ! all zeros
   call run_sync(k, 0, 0, 3, 8, 2.0, (/1,0,1,0,0,0,0,0/), 5, 1500.0, 0, 0, 1.0)     ! all zeros, 5 peaks
   call run_sync(k, s, 0, 3, 8, 2.0, (/1,1,1,0,0,0,0,0/), 2, 1.0e30, 0, 0, 1.0)     ! a huge carrier
   call run_sync(k, s, 0, 3, 8, 2.0, (/1,1,1,0,0,0,0,0/), 2, -1500.0, 0, 0, 1.0)    ! a negative one
   call run_sync(k, s, 0, 5, 8, 2.5, (/0,1,1,0,1,0,0,0/), 4, 1500.0, 0, 0, 1.0)     ! a gap in the mask
   call run_sync(k, s, 0, 8, 40, 2.0, (/1,1,1,1,1,0,0,0/), 2, 1500.0, 0, 0, 1.0)    ! wide
   call run_sync(k, s, 0, 8, 15, 0.5, (/1,1,1,1,0,0,0,0/), 2, 1500.0, 0, 0, 1.0)    ! fine steps

   ! ---- msk144spd: every slice as mskrtd runs it (np = 6912, ntol 200, fc 1500), then the odd cases
   ! (k s n ntol fc nan inf scale)
   do s = 1, nslice
      call run_spd(k, s, 6912, 200, 1500.0, 0, 0, 1.0)
   enddo
   s = sdec1
   call run_spd(k, s, 7168, 200, 1500.0, 0, 0, 1.0)
   call run_spd(k, s, 6912, 15, 1500.0, 0, 0, 1.0)
   call run_spd(k, s, 6912, 40, 1500.0, 0, 0, 1.0)
   call run_spd(k, s, 6912, 0, 1500.0, 0, 0, 1.0)         ! ntol 0
   call run_spd(k, s, 6912, -1, 1500.0, 0, 0, 1.0)        ! ntol < 0: nothing done
   call run_spd(k, s, 6912, 100000, 1500.0, 0, 0, 1.0)    ! windows clipped to the spectrum
   call run_spd(k, s, 6912, 1000, 1500.0, 0, 0, 1.0)
   call run_spd(k, s, 2591, 200, 1500.0, 0, 0, 1.0)       ! n < 3*NSPM: nothing done
   call run_spd(k, s, 2592, 200, 1500.0, 0, 0, 1.0)       ! the least n: 8 steps
   call run_spd(k, s, 2593, 200, 1500.0, 0, 0, 1.0)
   call run_spd(k, s, 3000, 200, 1500.0, 0, 0, 1.0)
   call run_spd(k, s, 4321, 200, 1500.0, 0, 0, 1.0)
   call run_spd(k, s, 6912, 200, 1400.0, 0, 0, 1.0)
   call run_spd(k, s, 6912, 200, 1503.3, 0, 0, 1.0)
   call run_spd(k, s, 6912, 200, 5300.0, 0, 0, 1.0)       ! the high window clipped
   call run_spd(k, s, 6912, 200, 5500.0, 0, 0, 1.0)       ! the high window gone
   call run_spd(k, s, 6912, 200, 6000.0, 0, 0, 1.0)
   call run_spd(k, s, 6912, 200, -600.0, 0, 0, 1.0)       ! the low window clipped
   call run_spd(k, s, 6912, 200, -1000.0, 0, 0, 1.0)      ! the low window gone
   call run_spd(k, s, 6912, 200, 500.0, 0, 0, 1.0)
   call run_spd(k, s, 6912, 200, nanv, 0, 0, 1.0)         ! fc NaN
   call run_spd(k, s, 6912, 200, infv, 0, 0, 1.0)
   call run_spd(k, s, 6912, 200, -infv, 0, 0, 1.0)
   call run_spd(k, s, 6912, 200, 1.0e30, 0, 0, 1.0)
   call run_spd(k, s, 6912, 200, 1500.0, 3000, 0, 1.0)    ! a NaN point
   call run_spd(k, s, 6912, 200, 1500.0, 1, 0, 1.0)
   call run_spd(k, s, 6912, 200, 1500.0, 6912, 0, 1.0)
   call run_spd(k, s, 6912, 200, 1500.0, 0, 100, 1.0)     ! an Inf point
   call run_spd(k, s, 6912, 200, 1500.0, 2000, 4000, 1.0)
   call run_spd(k, s, 6912, 200, 1500.0, 0, -3000, 1.0)   ! one point of 1e20: its square is Inf, spread by the FFT
   call run_spd(k, s, 6912, 200, 1500.0, 0, -100, 1.0)
   call run_spd(k, s, 6912, 200, 1500.0, 0, -6000, 1.0)
   call run_spd(k, s, 6912, 200, 1500.0, 0, -2800, 1.0)   ! in the tail edge (the last 12 points) of step 10's window
   call run_spd(k, s, 6912, 200, 1500.0, 0, -860, 1.0)    ! and of step 1's
   call run_spd(k, s, 6912, 200, 1500.0, 0, -2161, 1.0)
   call run_spd(k, s, 6912, 200, 1500.0, 0, -2162, 1.0)
   call run_spd(k, s, 6912, 200, 1500.0, 0, -1297, 1.0)
   call run_spd(k, s, 6912, 200, 1500.0, 0, -4000, 1.0)
   call run_spd(k, s, 6912, 200, 1500.0, 0, 0, 1.0e20)    ! overflow in the squares
   call run_spd(k, s, 6912, 200, 1500.0, 0, 0, 1.0e-20)
   call run_spd(k, s, 6912, 200, 1500.0, 0, 0, 1.0e-30)
   call run_spd(k, s, 6912, 200, 1500.0, 0, 0, 0.0)       ! zeros made by scaling
   call run_spd(k, 0, 6912, 200, 1500.0, 0, 0, 1.0)       ! all zeros
   call run_spd(k, 0, 2592, 200, 1500.0, 0, 0, 1.0)
   call run_spd(k, shot1, 6912, 200, 1500.0, 2500, 0, 1.0)
   call run_spd(k, nsl_file(2) - 1, 6912, 200, 1500.0, 0, 5000, 1.0)
   call run_spd(k, nsl_file(3) - 1, 6912, 200, 1500.0, 0, 0, 1.0e20)
   call run_spd(k, 99, 28000, 200, 1500.0, 0, 0, 1.0)     ! 125 steps asked for, 100 allowed
   call run_spd(k, 99, 22464, 200, 1500.0, 0, 0, 1.0)     ! 99
   call run_spd(k, 99, 22680, 200, 1500.0, 0, 0, 1.0)     ! 101 asked for
   call run_spd(k, 99, 28000, 200, 1500.0, 15000, 0, 1.0)
   ! the end marker: a truncated output (a failed self-check, a crash) cannot match without it
   write(*, '(a,1x,i0)') 'END sync', k

contains

   ! analytic signal of the int16 block, as mskrtd makes it
   subroutine block_analytic(w, i0, c)
      integer*2 :: w(*)
      integer :: i0, ii
      complex :: c(NF)
      real :: dd(NF), r
      real*8 :: pcc(5)
      logical*1 :: bq
      dd = 0.0
      dd(1:NZ) = w(i0+1:i0+NZ)
      r = sqrt(sum(dd(1:NZ)*dd(1:NZ))/NZ)
      if (r < 1.0) r = 1.0
      dd(1:NZ) = (1.0/r)*dd(1:NZ)
      dd(NZ+1:NF) = 0.0
      pcc = 0.0d0
      bq = .true.
      ii = 0
      call analytic(dd, NZ, NF, c, pcc, bq)
   end subroutine block_analytic

   ! the clean MSK frame (864 points, c = xi + j xq, genmsk_128_90's half-sine pulses) of 77 bits '0'/'1'
   subroutine msk_frame(mbits, cf)
      character(len=*) :: mbits
      complex :: cf(NSPM)
      integer*1 :: m(77), cwk(128), bitseq(144), s8(8)
      real :: pp(12), xi(NSPM), xq(NSPM)
      integer :: ii, is
      data s8 /0, 1, 1, 1, 0, 0, 1, 0/
      read(mbits, '(77i1)') m
      call encode_128_90(m, cwk)
      bitseq(1:8) = s8
      bitseq(9:56) = cwk(1:48)
      bitseq(57:64) = s8
      bitseq(65:144) = cwk(49:128)
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
      cf = cmplx(xi, xq)
   end subroutine msk_frame

   ! synbits are the 77 bits of messages with known fates: the library's own decoder, on the clean frame of
   ! each, must give these (the strings were found with scripts/msk144-frame-oracle.f90); the oracle stops if they no longer do
   subroutine check_synbits()
      complex :: cf(NSPM), ctt(NSPM)
      real :: sb(144)
      character(len=37) :: mg
      character(len=37), parameter :: want(3) = (/ 'CQ <...>                             ', &
         'RQ3NXC/R QU1RK R BR49                ', 'K1ABC W9XYZ RR73                     ' /)
      integer, parameter :: wantns(3) = (/ 0, 0, 1 /)
      integer :: ii, nsx
      do ii = 1, 3
         call msk_frame(synbits(ii), cf)
         ctt = cf
         call msk144decodeframe(ctt, sb, mg, nsx)
         if (mg /= want(ii) .or. nsx /= wantns(ii)) error stop 'oracle: a synbits message no longer decodes as recorded'
      enddo
   end subroutine check_synbits

   ! a slice of 7168 points: the analytic signal of the MSK frame of these 77 bits repeated from point `start` on at
   ! 1500 Hz (the frame's samples are the half-sine pulses of genmsk_128_90's xi and xq, c = xi + j xq, times
   ! sqrt 2 for unit rms of the real signal), with Gaussian noise of this sigma on each part; a new slice
   subroutine synth_slice(mbits, sigma, start)
      character(len=*) :: mbits
      real :: sigma
      integer :: start
      real :: ang, zr, zi
      complex :: cf(NSPM), z(NZ)
      integer :: ii, jj
      call msk_frame(mbits, cf)
      do ii = 1, NZ
         jj = mod(ii - start + 864*10, 864) + 1
         ang = real(mod(6.283185307179586d0*1500.0d0*dble(ii - 1)/12000.0d0, 6.283185307179586d0))
         z(ii) = 1.4142135*cf(jj)*cmplx(cos(ang), sin(ang))
         zr = real(z(ii)); zi = aimag(z(ii))
         if (ii < start) then
            zr = 0.0; zi = 0.0
         endif
         zr = zr + sigma*sgauss()
         zi = zi + sigma*sgauss()
         z(ii) = cmplx(zr, zi)
      enddo
      nslice = nslice + 1
      sl(1:NZ, nslice) = z
      call print_slice(nslice, sl(1,nslice), NZ)
   end subroutine synth_slice

   ! the synthetic slices' noise: a sum of 12 uniforms less 6, from a 32-bit LCG (seed 4242)
   real function sgauss()
      real*8 :: t
      integer :: kk
      t = -6.0d0
      do kk = 1, 12
         sst = mod(sst*1664525_8 + 1013904223_8, 4294967296_8)
         t = t + real(sst, 8)/4294967296.0d0
      enddo
      sgauss = real(t)
   end function sgauss

   subroutine add_slice(ks, c)
      integer :: ks
      complex :: c(NF)
      sl(1:NZ, ks) = c(1:NZ)
      call print_slice(ks, sl(1,ks), NZ)
   end subroutine add_slice

   subroutine print_slice(ks, c, n)
      integer :: ks, n
      complex :: c(*)
      call prc('I', ks, c, n)
   end subroutine print_slice

   ! the library's msk144spd on a slice, nothing printed
   subroutine quiet_spd(c, n, ntol, fc, ns)
      complex :: c(*)
      integer :: n, ntol, ns, navgq
      real :: fc, fq, tq, sbq(144)
      complex :: ctq(NSPM)
      character(len=37) :: mq
      ctq = 0.0
      sbq = 0.0
      call msk144spd(c, n, ntol, ns, mq, fc, fq, tq, navgq, ctq, sbq)
   end subroutine quiet_spd

   ! a slice, a carrier and a float scale at which msk144sync's xmax (3 frames from the slice's start, mask 1 1 1,
   ! ntol 8, delf 2) is exactly 1.3: xmax steps by more than one bit per step of the scale, so some of the
   ! variants (slice, carrier) never land on it; the first that does is returned
   subroutine find_scale(ss, fcv, scl)
      integer :: ss, i, j, jj, nk(10), nsu, mm(8), s1, ifc
      real :: scl, xm, fe, g, fcv
      complex :: cd(NSPM*8), c0(NSPM)
      mm = (/1,1,1,0,0,0,0,0/)
      do s1 = 1, nslice
         do ifc = 0, 3
            fcv = 1500.0 + 0.37*real(ifc)
            call make_input(s1, 0, NSPM*3, 0, 0, 1.0, cd)
            call msk144sync(cd, 3, 8, 2.0, mm, 2, fcv, fe, nk, nsu, xm, c0)
            if (xm <= 0.0) cycle
            scl = 1.3/xm
            do i = 0, 60
               do j = -1, 1, 2
                  g = scl
                  do jj = 1, i
                     g = nearest(g, real(j))
                  enddo
                  call make_input(s1, 0, NSPM*3, 0, 0, g, cd)
                  call msk144sync(cd, 3, 8, 2.0, mm, 2, fcv, fe, nk, nsu, xm, c0)
                  if (transfer(xm, 0) == int(z'3FA66666')) then
                     scl = g
                     ss = s1
                     return
                  endif
               enddo
            enddo
         enddo
      enddo
      error stop 'oracle: no scale gives xmax 1.3'
   end subroutine find_scale

   subroutine navpat(ip, m)
      integer :: ip, m(8)
      integer :: np(3, 6)
      data np/0,1,0, 1,0,0, 0,0,1, 1,1,0, 0,1,1, 1,1,1/
      m = 0
      m(1:3) = np(1:3, ip)
   end subroutine navpat

   subroutine make_input(s, off, n, nanidx, infidx, scale, c)
      integer :: s, off, n, nanidx, infidx, i
      real :: scale
      complex :: c(*)
      if (s == 0) then
         c(1:n) = (0.0, 0.0)
      else if (s == 99) then
         c(1:n) = cl(off+1:off+n)
      else
         c(1:n) = sl(off+1:off+n, s)
      endif
      if (s /= 0) then
         do i = 1, n
            c(i) = cmplx(scale*real(c(i)), scale*aimag(c(i)))
         enddo
      endif
      if (nanidx > 0) c(nanidx) = cmplx(nanv, aimag(c(nanidx)))
      if (infidx > 0) c(infidx) = cmplx(infv, aimag(c(infidx)))
      if (infidx < 0) c(-infidx) = cmplx(1.0e20, aimag(c(-infidx)))
   end subroutine make_input

   subroutine run_sync(kk, s, off, nf, ntol, delf, mask, npeaks, fc, nanidx, infidx, scale)
      integer :: kk, s, off, nf, ntol, mask(8), npeaks, nanidx, infidx
      real :: delf, fc, scale
      complex :: cdat(NSPM*8), c1(NSPM), c2(NSPM)
      integer :: nk1(10), nk2(10), nsuc1, nsuc2, m1(8)
      real :: fest1, fest2, xmax1, xmax2
      kk = kk + 1
      kc = kk
      cdat = 0.0
      if (nf > 0) call make_input(s, off*NSPM, NSPM*nf, nanidx, infidx, scale, cdat)
      write(*, '(a,1x,i0,1x,i0,1x,i0,1x,i0,1x,i0,1x,z8.8,8(1x,i0),1x,i0,1x,z8.8,2(1x,i0),1x,z8.8)') 'SYN', kk, s, off, &
         nf, ntol, transfer(delf, 0), mask, npeaks, transfer(fc, 0), nanidx, infidx, transfer(scale, 0)
      nk1 = -1; nk2 = -1
      c1 = (7.0, -7.0); c2 = (7.0, -7.0)
      fest1 = -1.0; fest2 = -1.0; xmax1 = -1.0; xmax2 = -1.0
      nsuc1 = -1; nsuc2 = -1
      m1 = mask
      call msk144sync_trace(cdat, nf, ntol, delf, m1, npeaks, fc, fest1, nk1, nsuc1, xmax1, c1)
      call msk144sync(cdat, nf, ntol, delf, m1, npeaks, fc, fest2, nk2, nsuc2, xmax2, c2)
      if (nsuc1 /= nsuc2 .or. transfer(fest1, 0) /= transfer(fest2, 0) .or. transfer(xmax1, 0) /= transfer(xmax2, 0) &
          .or. any(nk1(1:npeaks) /= nk2(1:npeaks)) .or. &
          any(transfer(c1, (/0/)) /= transfer(c2, (/0/)))) error stop 'oracle: msk144sync_trace differs from msk144sync'
      write(*, '(a,1x,i0,1x,i0,2(1x,z8.8),10(1x,i0))') 'X', kk, nsuc1, transfer(fest1, 0), transfer(xmax1, 0), &
         nk1(1:npeaks)
   end subroutine run_sync

   subroutine run_spd(kk, s, n, ntol, fc, nanidx, infidx, scale)
      integer :: kk, s, n, ntol, nanidx, infidx
      real :: fc, scale
      complex, allocatable :: cbig(:)
      complex :: ct1(NSPM), ct2(NSPM)
      real :: sb1(144), sb2(144), fr1, fr2, tr1, tr2
      integer :: na1, na2, ns1_, ns2_
      character(len=37) :: m1, m2
      kk = kk + 1
      kc = kk
      nframe = 0
      allocate(cbig(max(n, 1)))
      call make_input(s, 0, n, nanidx, infidx, scale, cbig)
      write(*, '(a,1x,i0,1x,i0,1x,i0,1x,i0,1x,z8.8,2(1x,i0),1x,z8.8)') 'SPD', kk, s, n, ntol, transfer(fc, 0), &
         nanidx, infidx, transfer(scale, 0)
      ct1 = (7.0, -7.0); ct2 = (7.0, -7.0)
      sb1 = -3.0; sb2 = -3.0
      fr1 = -1.0; fr2 = -1.0; tr1 = -1.0; tr2 = -1.0
      na1 = -1; na2 = -1
      m1 = 'unset'; m2 = 'unset'
      call msk144spd_trace(cbig, n, ntol, ns1_, m1, fc, fr1, tr1, na1, ct1, sb1)
      call msk144spd(cbig, n, ntol, ns2_, m2, fc, fr2, tr2, na2, ct2, sb2)
      if (ns1_ /= ns2_ .or. m1 /= m2 .or. transfer(fr1, 0) /= transfer(fr2, 0) .or. transfer(tr1, 0) /= transfer(tr2, 0) &
          .or. na1 /= na2 .or. any(transfer(sb1, (/0/)) /= transfer(sb2, (/0/)))) &
         error stop 'oracle: msk144spd_trace differs from msk144spd'
      ! ct: both are corrected in place by the frame decoder (msk144_frame_hook gives it the frame itself), so they
      ! agree on every case, success or failure, whether or not a frame was cut
      if (any(transfer(ct1, (/0/)) /= transfer(ct2, (/0/)))) &
         error stop 'oracle: msk144spd_trace differs from msk144spd (ct)'
      write(*, '(a,1x,i0,1x,i0,a,a,a,2(1x,z8.8),1x,i0)') 'X', kk, ns1_, ' |', m1, '|', transfer(fr1, 0), transfer(tr1, 0), na1
      call prr('U', kk, sb1, 144)
      call prc('V', kk, ct1, NSPM)
      deallocate(cbig)
   end subroutine run_spd

   ! ---- the pieces on their own
   subroutine direct_checks(kk)
      integer :: kk
      call maxloc_cases(kk)
      call indexx_cases(kk)
      call peak_cases(kk)
      call window_cases(kk)
   end subroutine direct_checks

   subroutine maxloc_cases(kk)
      integer :: kk, it, n, il, lo, hi, kind
      real :: a(100)
      integer*8 :: seed
      ! hand-made arrays
      call ml(kk, 0, 8, 1, 8, (/1.0, 5.0, 3.0, 5.0, 2.0, 0.0, 7.0, 7.0/))
      call ml(kk, 0, 8, 1, 8, (/nanv, 5.0, 3.0, nanv, 2.0, 0.0, 7.0, nanv/))
      call ml(kk, 0, 8, 1, 8, (/nanv, nanv, nanv, nanv, nanv, nanv, nanv, nanv/))
      call ml(kk, 0, 8, 1, 8, (/nanv, nanv, -infv, nanv, -infv, -infv, -infv, nanv/))
      call ml(kk, 0, 8, 1, 8, (/-1.0, -2.0, -3.0, -4.0, -5.0, -6.0, -7.0, -8.0/))
      call ml(kk, 0, 8, 1, 8, (/nanv, -1.0, -2.0, -3.0, -4.0, -5.0, -6.0, -8.0/))
      call ml(kk, 0, 8, 1, 8, (/infv, 1.0, infv, nanv, 1.0, 1.0, 1.0, 1.0/))
      call ml(kk, 0, 8, 1, 8, (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/))
      call ml(kk, 0, 8, 1, 8, (/-0.0, 0.0, -0.0, 0.0, -0.0, 0.0, -0.0, 0.0/))
      call ml(kk, 0, 8, 1, 8, (/0.0, -0.0, 0.0, -0.0, 0.0, -0.0, 0.0, -0.0/))
      call ml(kk, 0, 8, 1, 8, (/nanv, 2.0, nanv, 3.0, nanv, 3.0, nanv, 1.0/))
      call ml(kk, 0, 1, 1, 1, (/nanv, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/))
      call ml(kk, 0, 1, 1, 1, (/2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/))
      call ml(kk, 0, 4, 1, 4, (/1.0, 2.0, 3.0, 9.0, 0.0, 0.0, 0.0, 0.0/))
      call ml(kk, 1, 8, 3, 6, (/9.0, 9.0, 1.0, 5.0, 5.0, 2.0, 9.0, 9.0/))
      call ml(kk, 1, 8, 3, 6, (/9.0, 9.0, nanv, nanv, nanv, nanv, 9.0, 9.0/))
      call ml(kk, 1, 8, 3, 6, (/9.0, 9.0, nanv, 4.0, nanv, nanv, 9.0, 9.0/))
      call ml(kk, 1, 8, 3, 6, (/9.0, 9.0, 0.0, 0.0, 0.0, 0.0, 9.0, 9.0/))
      call ml(kk, 1, 8, 3, 6, (/9.0, 9.0, -1.0, -1.0, -3.0, -1.0, 9.0, 9.0/))
      call ml(kk, 1, 8, 3, 6, (/9.0, 9.0, -infv, nanv, -infv, nanv, 9.0, 9.0/))
      call ml(kk, 1, 8, 1, 1, (/nanv, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0/))
      call ml(kk, 1, 8, 8, 8, (/9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, infv/))
      call ml(kk, 1, 8, 4, 3, (/9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 1.0/))
      ! pseudo-random arrays with ties and NaNs
      seed = 12345_8
      do it = 1, 60
         n = 8 + mod(it*7, 90)
         kind = mod(it, 2)
         do j = 1, n
            seed = mod(seed*1103515245_8 + 12345_8, 2147483648_8)
            a(j) = real(mod(seed/65536_8, 7_8)) - 3.0
            if (mod(seed/1024_8, 9_8) == 0) a(j) = nanv
            if (mod(seed/1024_8, 11_8) == 1) a(j) = -infv
         enddo
         lo = 1 + mod(it, 5)
         hi = n - mod(it, 4)
         call ml(kk, kind, n, lo, hi, a(1:n))
      enddo
      il = 0
   end subroutine maxloc_cases

   subroutine ml(kk, kind, n, lo, hi, a)
      integer :: kk, kind, n, lo, hi, il, i
      real :: a(n)
      real :: det(-2:103), ts(864)
      logical :: m(864)
      kk = kk + 1
      if (kind == 0) then
         det = 0.0
         det(1:n) = a
         call mlu(det, n, il)
      else
         ts = 0.0
         ts(1:n) = a
         m = .false.
         if (lo <= hi) m(lo:hi) = .true.
         call mlm(ts, m, il)
      endif
      write(*, '(a,1x,i0,1x,i0,1x,i0,1x,i0,1x,i0,1x,i0,100(1x,z8.8))') 'L', kk, kind, n, lo, hi, il, &
         (transfer(a(i), 0), i = 1, n)
   end subroutine ml

   subroutine indexx_cases(kk)
      integer :: kk, it, n, j
      real :: a(100)
      integer*8 :: seed
      call jx(kk, 8, (/3.0, 1.0, 2.0, 5.0, 4.0, 9.0, 0.0, 7.0/))
      call jx(kk, 8, (/1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0/))
      call jx(kk, 8, (/8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0/))
      call jx(kk, 8, (/1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0/))
      call jx(kk, 8, (/nanv, 2.0, 1.0, nanv, 3.0, 0.0, 4.0, nanv/))
      call jx(kk, 8, (/nanv, nanv, nanv, nanv, nanv, nanv, nanv, nanv/))
      call jx(kk, 8, (/infv, 2.0, -infv, nanv, 3.0, 0.0, infv, -infv/))
      call jx(kk, 1, (/5.0/))
      call jx(kk, 2, (/5.0, 4.0/))
      call jx(kk, 3, (/5.0, 4.0, 9.0/))
      seed = 777_8
      do it = 1, 40
         n = 7 + mod(it*13, 94)
         do j = 1, n
            seed = mod(seed*1103515245_8 + 12345_8, 2147483648_8)
            a(j) = real(mod(seed/65536_8, 12_8)) - 3.0
            if (it > 30 .and. mod(seed/1024_8, 8_8) == 0) a(j) = nanv
         enddo
         call jx(kk, n, a(1:n))
      enddo
   end subroutine indexx_cases

   subroutine jx(kk, n, a)
      integer :: kk, n, idx(100), j
      real :: a(n)
      kk = kk + 1
      call indexx(a, n, idx)
      write(*, '(a,1x,i0,1x,i0)', advance='no') 'J', kk, n
      do j = 1, n
         write(*, '(1x,z8.8)', advance='no') transfer(a(j), 0)
      enddo
      write(*, '(a)', advance='no') ' |'
      do j = 1, n
         write(*, '(1x,i0)', advance='no') idx(j)
      enddo
      write(*, '(a)') ''
   end subroutine jx

   subroutine peak_cases(kk)
      integer :: kk, it, nb, j, pb
      complex :: sp(40)
      integer*8 :: seed
      real :: eps, nzero
      nzero = transfer(int(z'80000000'), 0.0)
      ! a smooth peak, and the awkward ones
      sp = (0.0, 0.0)
      sp(1:5) = (/(1.0, 0.5), (3.0, 1.0), (8.0, -2.0), (2.0, 0.5), (1.0, 0.0)/)
      call pk(kk, 5, 3, sp)
      call pk(kk, 5, 2, sp)
      call pk(kk, 5, 1, sp)      ! the first bin: nothing done
      call pk(kk, 5, 5, sp)      ! the last
      call pk(kk, 5, 0, sp)
      call pk(kk, 5, 6, sp)
      sp = (0.0, 0.0)
      call pk(kk, 5, 3, sp)      ! scale 0
      sp = (1.0, 0.0)
      call pk(kk, 5, 3, sp)      ! flat: no curvature
      sp = (0.0, 0.0)
      sp(1:5) = (/(1.0, 0.0), (3.0, 0.0), (1.0, 0.0), (0.0, 0.0), (0.0, 0.0)/)
      call pk(kk, 5, 2, sp)      ! curvature 4, offset 0
      sp(1:5) = (/(0.0, 1.0), (0.0, 3.0), (0.0, 1.0), (0.0, 0.0), (0.0, 0.0)/)
      call pk(kk, 5, 2, sp)      ! |d| > |c|: the other branch
      sp(1:5) = (/(1.0, 4.0), (1.0, 9.0), (1.5, 3.0), (0.0, 0.0), (0.0, 0.0)/)
      call pk(kk, 5, 2, sp)
      sp(1:5) = (/(4.0, 1.0), (9.0, 1.0), (3.0, 1.5), (0.0, 0.0), (0.0, 0.0)/)
      call pk(kk, 5, 2, sp)
      sp(1:5) = (/(1.0, 1.0), (1.0, 1.0), (1.0, 1.0), (0.0, 0.0), (0.0, 0.0)/)
      call pk(kk, 5, 2, sp)
      eps = epsilon(1.0)
      sp(1:5) = (/(1.0, 0.0), cmplx(1.0 + 4.0*eps, 0.0), (1.0, 0.0), (0.0, 0.0), (0.0, 0.0)/)
      call pk(kk, 5, 2, sp)      ! curvature near epsilon*scale
      sp(1:5) = (/(1.0, 0.0), cmplx(1.0 + 2.0*eps, 0.0), (1.0, 0.0), (0.0, 0.0), (0.0, 0.0)/)
      call pk(kk, 5, 2, sp)
      sp(1:5) = (/(1.0, 0.0), cmplx(1.0 + 1.0*eps, 0.0), (1.0, 0.0), (0.0, 0.0), (0.0, 0.0)/)
      call pk(kk, 5, 2, sp)
      sp(1:5) = (/cmplx(nanv, 0.0), (3.0, 0.0), (1.0, 0.0), (0.0, 0.0), (0.0, 0.0)/)
      call pk(kk, 5, 2, sp)
      sp(1:5) = (/(1.0, 0.0), cmplx(nanv, 0.0), (1.0, 0.0), (0.0, 0.0), (0.0, 0.0)/)
      call pk(kk, 5, 2, sp)
      sp(1:5) = (/(1.0, 0.0), (3.0, 0.0), cmplx(nanv, 0.0), (0.0, 0.0), (0.0, 0.0)/)
      call pk(kk, 5, 2, sp)
      sp(1:5) = (/cmplx(infv, 0.0), (3.0, 0.0), (1.0, 0.0), (0.0, 0.0), (0.0, 0.0)/)
      call pk(kk, 5, 2, sp)
      sp(1:5) = (/(1.0, 0.0), cmplx(infv, 0.0), (1.0, 0.0), (0.0, 0.0), (0.0, 0.0)/)
      call pk(kk, 5, 2, sp)
      sp(1:5) = (/(1.0e30, 1.0e30), (3.0e30, 0.0), (1.0e30, 5.0e30), (0.0, 0.0), (0.0, 0.0)/)
      call pk(kk, 5, 2, sp)
      sp(1:5) = (/(1.0e-30, 1.0e-30), (3.0e-30, 0.0), (1.0e-30, 5.0e-30), (0.0, 0.0), (0.0, 0.0)/)
      call pk(kk, 5, 2, sp)
      sp(1:5) = (/(1.0, 0.0), (1.0, 0.0), (3.0, 0.0), cmplx(2.0e-38, transfer(77777, 0.0)), cmplx(transfer(1, 0.0), 0.0)/)
      call pk(kk, 5, 3, sp)
      sp(1:5) = (/(-0.0, 0.0), (-0.0, -0.0), (0.0, -0.0), (0.0, 0.0), (0.0, 0.0)/)
      call pk(kk, 5, 2, sp)
      sp(1:5) = (/(1.0, 0.0), (-2.0, 0.0), (1.0, 0.0), (0.0, 0.0), (0.0, 0.0)/)
      call pk(kk, 5, 2, sp)      ! a trough
      ! |c| = |d| (a tie in the division's branch), c = -d and c = d
      sp(1:5) = (/cmplx(0.75, 0.25), cmplx(1.0, -0.5), cmplx(0.25, -0.25), (0.0, 0.0), (0.0, 0.0)/)
      call pk(kk, 5, 2, sp)
      sp(1:5) = (/cmplx(0.5, 1.0), cmplx(1.0, 1.0), cmplx(0.5, 0.0), (0.0, 0.0), (0.0, 0.0)/)
      call pk(kk, 5, 2, sp)
      ! signed zeros in the curvature: s0 = (1, -0): 2.0*s0 as a complex product has imaginary part +0
      sp(1:5) = (/cmplx(nzero, 0.0), cmplx(1.0, nzero), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0)/)
      call pk(kk, 5, 2, sp)
      sp(1:5) = (/cmplx(nzero, nzero), cmplx(1.0, nzero), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0)/)
      call pk(kk, 5, 2, sp)
      sp(1:5) = (/cmplx(0.0, nzero), cmplx(2.0, nzero), cmplx(nzero, 0.0), (0.0, 0.0), (0.0, 0.0)/)
      call pk(kk, 5, 2, sp)
      seed = 4242_8
      do it = 1, 60
         nb = 5
         do j = 1, nb
            seed = mod(seed*1103515245_8 + 12345_8, 2147483648_8)
            sp(j) = cmplx(real(mod(seed/65536_8, 2001_8) - 1000_8)/100.0, real(mod(seed/131072_8, 2001_8) - 1000_8)/100.0)
         enddo
         pb = 2 + mod(it, 3)
         call pk(kk, nb, pb, sp)
      enddo
   end subroutine peak_cases

   subroutine pk(kk, nb, pb, sp)
      integer :: kk, nb, pb, j
      complex :: sp(40)
      real :: r
      kk = kk + 1
      r = msk_peak_offset(sp(1:nb), pb)
      write(*, '(a,1x,i0,1x,i0,1x,i0,1x,z8.8,80(1x,z8.8))') 'K', kk, nb, pb, transfer(r, 0), &
         (transfer(real(sp(j)), 0), transfer(aimag(sp(j)), 0), j = 1, nb)
   end subroutine pk

   subroutine window_cases(kk)
      integer :: kk
      real*8 :: nan8, inf8
      nan8 = transfer(int(z'7FF8000000000000', 8), 0.0d0)
      inf8 = transfer(int(z'7FF0000000000000', 8), 0.0d0)
      call wn(kk, 864, 260.0d0, 316.0d0)
      call wn(kk, 864, 1.0d0, 864.0d0)
      call wn(kk, 864, 0.0d0, 100.0d0)
      call wn(kk, 864, -5.0d0, 100.0d0)
      call wn(kk, 864, 800.0d0, 900.0d0)
      call wn(kk, 864, 864.4d0, 900.0d0)
      call wn(kk, 864, 864.5d0, 900.0d0)
      call wn(kk, 864, 865.0d0, 900.0d0)
      call wn(kk, 864, -10.0d0, 0.4d0)
      call wn(kk, 864, -10.0d0, 0.5d0)
      call wn(kk, 864, -10.0d0, 0.6d0)
      call wn(kk, 864, 100.5d0, 100.5d0)
      call wn(kk, 864, 100.4d0, 100.6d0)
      call wn(kk, 864, 100.6d0, 100.4d0)
      call wn(kk, 864, 2.5d0, 3.5d0)
      call wn(kk, 864, 1.5d0, 1.5d0)
      call wn(kk, 864, nan8, 100.0d0)
      call wn(kk, 864, 100.0d0, nan8)
      call wn(kk, 864, -inf8, inf8)
      call wn(kk, 864, 100.0d0, inf8)
      call wn(kk, 864, -inf8, 100.0d0)
      call wn(kk, 864, 1.0d300, 1.0d301)
      call wn(kk, 864, -1.0d300, 1.0d300)
      call wn(kk, 0, 1.0d0, 2.0d0)
      call wn(kk, 1, 1.0d0, 2.0d0)
      call wn(kk, 1, 0.6d0, 1.4d0)
      call wn(kk, -3, 1.0d0, 2.0d0)
      call wn(kk, 2147483647, 1.0d0, 3.0d9)
      call wn(kk, 2147483647, 2.0d9, 3.0d9)
   end subroutine window_cases

   subroutine wn(kk, nfft, rlo, rhi)
      integer :: kk, nfft, lo, hi
      real*8 :: rlo, rhi
      logical :: u
      kk = kk + 1
      call msk_clip_spectrum_window(nfft, rlo, rhi, lo, hi, u)
      write(*, '(a,1x,i0,1x,i0,2(1x,z16.16),3(1x,i0))') 'W', kk, nfft, transfer(rlo, 0_8), transfer(rhi, 0_8), lo, hi, &
         merge(1, 0, u)
   end subroutine wn

end program msk144_sync_oracle

! maxloc as msk144spd writes it
subroutine mlu(detmet, nstep, il)
   implicit none
   integer :: nstep, il
   integer, dimension(1) :: iloc
   real :: detmet(-2:103)
   iloc = maxloc(detmet(1:nstep))
   il = iloc(1)
end subroutine mlu

subroutine mlm(tonespec, ismask, il)
   implicit none
   integer :: il
   integer, dimension(1) :: iloc
   real :: tonespec(864)
   logical :: ismask(864)
   iloc = maxloc(tonespec, ismask)
   il = iloc(1)
end subroutine mlm

! ---- printing: real and complex arrays, 16 values a line
subroutine prr(tag, k, x, n)
   implicit none
   character(len=*) :: tag
   integer :: k, n, i, i2, j
   real :: x(*)
   do i = 1, n, 16
      i2 = min(i + 15, n)
      write(*, '(a,1x,i0,1x,i0,16(1x,z8.8))') tag, k, i, (transfer(x(j), 0), j = i, i2)
   enddo
end subroutine prr

subroutine prc(tag, k, x, n)
   implicit none
   character(len=*) :: tag
   integer :: k, n, i, i2, j
   complex :: x(*)
   real :: v(2*n)
   do j = 1, n
      v(2*j-1) = real(x(j))
      v(2*j) = aimag(x(j))
   enddo
   do i = 1, 2*n, 16
      i2 = min(i + 15, 2*n)
      write(*, '(a,1x,i0,1x,i0,16(1x,z8.8))') tag, k, i, (transfer(v(j), 0), j = i, i2)
   enddo
end subroutine prc

! two hashes of the bits of n reals: h = h*m + w, modulo a prime, w the 32 bits unsigned
subroutine hashr(x, n, h1, h2)
   implicit none
   integer :: n, i
   real :: x(n)
   integer*8 :: h1, h2, w
   h1 = 1_8
   h2 = 1_8
   do i = 1, n
      w = iand(int(transfer(x(i), 0), 8), 4294967295_8)
      h1 = mod(h1*1000003_8 + w, 2147483647_8)
      h2 = mod(h2*999983_8 + w, 2147483629_8)
   enddo
end subroutine hashr

! what the instrumented copy of msk144spd calls where the reference calls msk144decodeframe
subroutine msk144_frame_hook(ct, softbits, msgreceived, nsuccess)
   implicit none
   complex :: ct(864)
   real :: softbits(144), v(1728)
   character(len=37) :: msgreceived
   integer :: nsuccess, kc, nframe, j
   integer*8 :: h1, h2
   common /skcase/ kc, nframe
   ! msk144decodeframe changes the frame it is given (the phase correction, in place), as it does in the library's
   ! msk144spd: hash the frame as it came, then give the decoder the frame itself: ct is left as the decoder left it
   do j = 1, 864
      v(2*j-1) = real(ct(j))
      v(2*j) = aimag(ct(j))
   enddo
   call hashr(v, 1728, h1, h2)
   call msk144decodeframe(ct, softbits, msgreceived, nsuccess)
   nframe = nframe + 1
   write(*, '(a,1x,i0,1x,i0,2(1x,i0),1x,i0,a,a,a)') 'F', kc, nframe, h1, h2, nsuccess, ' |', msgreceived, '|'
end subroutine msk144_frame_hook

! ---- the traces of msk144spd_trace and msk144sync_trace
subroutine msk144_spd_consts(rcw)
   implicit none
   real :: rcw(12)
   integer :: kc, nframe
   common /skcase/ kc, nframe
   call prr('R', kc, rcw, 12)
end subroutine msk144_spd_consts

subroutine msk144_spd_win(ihlo, ihhi, illo, ilhi, i2000, i4000)
   implicit none
   integer :: ihlo, ihhi, illo, ilhi, i2000, i4000, kc, nframe
   common /skcase/ kc, nframe
   write(*, '(a,1x,i0,6(1x,i0))') 'N', kc, ihlo, ihhi, illo, ilhi, i2000, i4000
end subroutine msk144_spd_win

subroutine msk144_spd_step(istp, ihpk, ilpk, deltah, deltal, ah, al, trath, tratl, ferr)
   implicit none
   integer :: istp, ihpk, ilpk, kc, nframe
   real :: deltah, deltal, ah, al, trath, tratl, ferr
   common /skcase/ kc, nframe
   write(*, '(a,1x,i0,3(1x,i0),7(1x,z8.8))') 'T', kc, istp, ihpk, ilpk, transfer(deltah, 0), transfer(deltal, 0), &
      transfer(ah, 0), transfer(al, 0), transfer(trath, 0), transfer(tratl, 0), transfer(ferr, 0)
end subroutine msk144_spd_step

subroutine msk144_spd_med(nstep, indices, xmed, detmet, detmet2, detfer)
   implicit none
   integer :: nstep, indices(100), kc, nframe, j
   real :: xmed, detmet(-2:103), detmet2(-2:103), detfer(100)
   common /skcase/ kc, nframe
   call prr('D', kc, detmet(1), nstep)
   call prr('H', kc, detmet2(1), nstep)
   call prr('O', kc, detfer, nstep)
   write(*, '(a,1x,i0,1x,i0,1x,z8.8,100(1x,i0))') 'M', kc, nstep, transfer(xmed, 0), (indices(j), j = 1, nstep)
end subroutine msk144_spd_med

subroutine msk144_spd_cands(ndet, nstart, ferrs, snrs)
   implicit none
   integer :: ndet, nstart(5), kc, nframe, j
   real :: ferrs(5), snrs(5)
   common /skcase/ kc, nframe
   write(*, '(a,1x,i0,1x,i0,1x,i0)') 'E', kc, 0, ndet
   do j = 1, ndet
      write(*, '(a,1x,i0,1x,i0,1x,i0,2(1x,z8.8))') 'E', kc, j, nstart(j), transfer(ferrs(j), 0), transfer(snrs(j), 0)
   enddo
end subroutine msk144_spd_cands

subroutine msk144_spd_sync(icand, iav, nsyncsuccess, fest, xmax, npkloc)
   implicit none
   integer :: icand, iav, nsyncsuccess, npkloc(10), kc, nframe
   real :: fest, xmax
   common /skcase/ kc, nframe
   write(*, '(a,1x,i0,1x,i0,1x,i0,1x,i0,2(1x,z8.8),2(1x,i0))') 'G', kc, icand, iav, nsyncsuccess, transfer(fest, 0), &
      transfer(xmax, 0), npkloc(1), npkloc(2)
end subroutine msk144_spd_sync

subroutine msk144_sync_consts(pp, cb)
   implicit none
   real :: pp(12)
   complex :: cb(42)
   integer :: kc, nframe
   common /skcase/ kc, nframe
   call prr('P', kc, pp, 12)
   call prc('B', kc, cb, 42)
end subroutine msk144_sync_consts

subroutine msk144_sync_peak(xmax, fest, c, xcc, nfreqs)
   implicit none
   real :: xmax, fest, xcc(0:863)
   complex :: c(864)
   integer :: nfreqs, kc, nframe
   common /skcase/ kc, nframe
   write(*, '(a,1x,i0,1x,i0,2(1x,z8.8))') 'Q', kc, nfreqs, transfer(xmax, 0), transfer(fest, 0)
   call prc('Y', kc, c, 864)
   call prr('Z', kc, xcc, 864)
end subroutine msk144_sync_peak
