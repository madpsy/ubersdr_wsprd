// FT4's front end: what jt9's ft4_decode does from a period's samples to
// each candidate's soft bits and back to the samples with a decode taken
// out (WSJT-X lib/ft4 at 567ad29ce:
// getcandidates4 with ft4_baseline; ft4_downsample; sync4d;
// get_ft4_bitmetrics; subtractft4 with gen_ft4wave; genft4.f90's entry
// get_ft4_tones_from_77bits), giving the reference's numbers to the bit
// on x86-64. ft4_decode itself (the subtraction passes, the candidate
// loop's DT segments and sync search, the decoding passes, a priori,
// decode174_91, the output) is the driver's (ft4.h); it calls these as
// ft4_decode calls them, with one Frontend as one jt9 process's SAVEd
// state.
//
// The reference is 567ad29ce, not WSJT-X improved:
// ft4_params' NMAX is 21*3456 (72,576), so the big FFT is 72,576 points and
// the downsampled period NDMAX 4,032. Its decoder never calls genft4's
// packing half (no a7 in 567ad29ce), only the tones entry.
//
// What it rests on, as FT8's front end (ft8_frontend.h): no multiply-adds
// fused (-ffp-contract=off, FT8_FP_OPTIONS in CMakeLists.txt, for every
// ft4_*.cpp); float sums in the Fortran's order; real*8 where the Fortran
// has it (ft4_baseline is implicit real*8); glibc 2.39's maths where
// gfortran calls it, carried (ft8_maths.h: sincosf, cosf, sinf, cabsf as
// hypotf, powf, log10f, lroundf as nint), and libmvec on exactly the
// elements the reference's -O3 objects send through it (their disassembly,
// site by site); FFTs as four2a
// plans them (fst4w_frontend.h's cached plans) on buffers aligned as the
// reference's (16 bytes or more at every FFT site: the vectors' '@' lines).
// The helpers FT8's decoder uses (nint, maxloc, pctile, polyfit, cis,
// cabs_, cmul, twkfreq1, the trace hashes, ...) are ft8_fortran.h's;
// nuttal_window, normalizebmet and gfsk_pulse are ft8_frontend.h's,
// called as they are (the same source at 567ad29ce).
//
// The work is split by routine group, one source file each, so that each
// can be written and tested on its own (tests/test_ft4_frontend.cpp,
// tags [ft4-fe-candidates], [ft4-fe-sync], [ft4-fe-metrics],
// [ft4-fe-wave], all of them [ft4-frontend]):
//   src/ft4_candidates.cpp  getcandidates4, ft4_baseline
//   src/ft4_sync.cpp        ft4_downsample, sync4d
//   src/ft4_metrics.cpp     get_ft4_bitmetrics
//   src/ft4_wave.cpp        gen_ft4wave, subtractft4,
//                           get_ft4_tones_from_77bits
//   src/ft4_frontend.cpp    the constructor, destructor, heap_bytes
// Each group keeps its SAVEd state in its own struct (Cand, Sync,
// Metrics, Wave), defined in its file only, with that file's make_*,
// del_* and *_bytes: no group's state is declared here, so no group edits
// this header.
//
// getcandidates4.f90, ft4_baseline.f90, ft4_downsample.f90, sync4d.f90,
// get_ft4_bitmetrics.f90, gen_ft4wave.f90, subtractft4.f90, genft4.f90,
// ft4_params.f90: Joe Taylor, K1JT, Steven Franke, K9AN, and the WSJT-X
// developers (GNU GPL v3), and so is this.
#pragma once

#include <complex>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <string>

namespace packjt77 {
class Packjt77;
}

namespace ft4 {

using cfloat = std::complex<float>;

// ft4_params.f90 (567ad29ce)
constexpr int KK = 91;                         // information bits (77 + CRC14)
constexpr int ND = 87;                         // data symbols
constexpr int NS = 16;                         // sync symbols
constexpr int NN = NS + ND;                    // 103, sync and data symbols
constexpr int NN2 = NS + ND + 2;               // 105, with the ramp symbols
constexpr int NSPS = 576;                      // samples per symbol at 12000 S/s
constexpr int NZ = NSPS * NN;                  // 59,328
constexpr int NZ2 = NSPS * NN2;                // 60,480, the shaped waveform
constexpr int NMAX = 21 * 3456;                // 72,576, samples in iwave
constexpr int NFFT1 = 2304, NH1 = NFFT1 / 2;   // symbol spectra
constexpr int NSTEP = NSPS;                    // coarse time-sync step
constexpr int NHSYM = (NMAX - NFFT1) / NSTEP;  // 122 symbol spectra
constexpr int NDOWN = 18;                      // downsample factor
// ft4_decode's, sync4d's, ft4_downsample's
constexpr int NSS = NSPS / NDOWN;              // 32 samples a symbol downsampled
constexpr int NDMAX = NMAX / NDOWN;            // 4,032, the downsampled period (NP, NFFT2)
constexpr int MAXCAND = 200;                   // ft4_decode's
// subtractft4's (its own parameters: NMAX hard-coded 21*3456 there too)
constexpr int NFILT = 1400;
constexpr int NFRAME = NN2 * NSPS;             // 60,480

using Trace = std::function<void(const std::string&)>;

// ft4_baseline(s, nfa, nfb, sbase): s and sbase are NH1 (Fortran's
// s(NH1), element i at [i-1]); s is changed (to dB over ia..ib) as there,
// sbase written over ia..ib only (the rest left as the caller had it).
// nfa and nfb are bins. Its trace lines (l, p) go to trace when given.
// (src/ft4_candidates.cpp)
void ft4_baseline(float* s, int nfa, int nfb, float* sbase, const Trace* trace = nullptr);

// genft4.f90's entry get_ft4_tones_from_77bits(msgbits, i4tone): msgbits
// is scrambled with rvec in place (the caller's array changes, as there:
// ft4_decode's message77 after its subtraction), i4tone NN tones 0..3.
// (src/ft4_wave.cpp)
void get_ft4_tones_from_77bits(int8_t msgbits[77], int i4tone[NN]);

// genft4(msg0, ichk, msgsent, msgbits, i4tone) with pk as the packjt77
// module's state (its hash tables; pack77_legacy_truncating_fallback and
// unpack77 as genft8's): msg0 cut at its first NUL and its leading blanks
// stripped as there (k blanks strip 1+2+... characters). Encoded and
// unpacked: msgsent the unpacked message, msgbits the 77 bits scrambled
// with rvec (as the entry leaves them), i4tone the tones. Not encoded,
// or not unpacked: msgbits and i4tone 0, msgsent "*** bad message ***".
// ichk 1: msgsent only (when encoded; msgbits and i4tone untouched).
// 567ad29ce's decoder never calls it (no a7). (src/ft4_wave.cpp)
void genft4(packjt77::Packjt77& pk, const std::string& msg0, int ichk, std::string& msgsent, int8_t msgbits[77],
            int i4tone[NN]);

// ft4_decode's reading of get_ft4_bitmetrics' result (bitmetrics 2*NN x 3
// as Frontend::get_ft4_bitmetrics writes it): nsync_qual, the hard-sync
// count of bitmetrics(:,1) >= 0 against the four Costas arrays' bits
// (0..32; the decoder goes on at 20 or more), and llra, llrb, llrc (2*ND
// each) from columns 1, 2, 3: the data rows, times scalefac 2.83.
// (src/ft4_metrics.cpp)
int ft4_nsync_qual(const float* bitmetrics);
void ft4_llrs(const float* bitmetrics, float* llra, float* llrb, float* llrc);

class Frontend {
public:
    Frontend();
    ~Frontend();
    Frontend(const Frontend&) = delete;
    Frontend& operator=(const Frontend&) = delete;

    // ---------------------------------------------- src/ft4_candidates.cpp
    // getcandidates4(dd, fa, fb, syncmin, nfqso, maxcand, savg, candidate,
    // ncand, sbase): dd NMAX; savg and sbase NH1 ([i-1]); candidate is
    // 2 x maxcand (Fortran's (i,k) at [(k-1)*2 + i-1]). On its early
    // return (any sbase(nfa:nfb) <= 0) ncand and candidate are left as the
    // caller set them (ft4_decode: 0 and zeros).
    void getcandidates4(const float* dd, float fa, float fb, float syncmin, int nfqso, int maxcand, float* savg,
                        float* candidate, int& ncand, float* sbase);

    // ---------------------------------------------- src/ft4_sync.cpp
    // ft4_downsample(dd, newdata, f0, c): c is NDMAX; the big r2c FFT
    // (x, SAVEd) recomputed when newdata, which is not cleared here (the
    // caller clears its dobigfft).
    void ft4_downsample(const float* dd, bool newdata, float f0, cfloat* c);

    // sync4d(cd0, i0, ctwk, itwk, sync): cd0 NDMAX, ctwk 2*NSS. ft4_decode
    // always passes itwk 1; with 0 csync2 (a local, uninitialised in the
    // reference) is whatever it last held here.
    float sync4d(const cfloat* cd0, int i0, const cfloat* ctwk, int itwk);

    // The sync search's calls, batched: sync[k] = sync4d(cd0,
    // istart0 + k*step, ctwk, 1), k = 0..n-1, bit for bit, with csync2
    // left as the last call leaves it. sync4d_prepare(cd0) first, after
    // cd0's last change (without it, or with a step other than 1 or 4,
    // each start is still exact, by sync4d's scalar sums).
    // tests/test_ft4_frontend.cpp checks it against sync4d, start by start.
    void sync4d_prepare(const cfloat* cd0);
    void sync4d_many(const cfloat* cd0, int istart0, int n, int step, const cfloat* ctwk, float* sync);

    // ---------------------------------------------- src/ft4_metrics.cpp
    // get_ft4_bitmetrics(cd, bitmetrics, badsync): cd NN*NSS; bitmetrics
    // is 2*NN x 3 (Fortran's (i,j) at [(j-1)*2*NN + i-1]), written only
    // when badsync is false (with badsync the reference returns before).
    void get_ft4_bitmetrics(const cfloat* cd, float* bitmetrics, bool& badsync);

    // get_ft4_bitmetrics in its plain form, line for line with the
    // reference: the same results, kept for tests/test_ft4_frontend.cpp's
    // side-by-side check ([ft4-fe-reference]).
    void get_ft4_bitmetrics_reference(const cfloat* cd, float* bitmetrics, bool& badsync);

    // ---------------------------------------------- src/ft4_wave.cpp
    // gen_ft4wave(itone, nsym, nsps, fsample, f0, cwave, wave, icmplx,
    // nwave): cwave or wave (icmplx 1 or 0) nwave long. Its pulse is
    // made on the first call only (SAVEd, from that call's nsps).
    void gen_ft4wave(const int* itone, int nsym, int nsps, float fsample, float f0, cfloat* cwave, float* wave,
                     int icmplx, int nwave);

    // gen_ft4wave in its plain form (without sincosf four at a time, mod
    // without fprem, the pulse products), line for line with the
    // reference: the same results, kept for the side-by-side check.
    void gen_ft4wave_reference(const int* itone, int nsym, int nsps, float fsample, float f0, cfloat* cwave,
                               float* wave, int icmplx, int nwave);

    // subtractft4(dd, itone, f0, dt): dd NMAX, changed in place; its
    // /heap4/ arrays and filter kept here.
    void subtractft4(float* dd, const int itone[NN], float f0, float dt);

    // ---------------------------------------------- src/ft4_frontend.cpp
    // Its buffers on the heap now, bytes (FFTW's plans
    // not counted). Not while it runs.
    size_t heap_bytes() const;

    // The insides as scripts/ft8-frontend-patch.py's instrumented FT4
    // reference writes them (tests/data/ft4/frontend/README.md, the lines
    // each routine writes; tests/test_ft4_frontend.cpp): one call a line,
    // none when unset.
    void set_trace(Trace f) { trace_ = std::move(f); }

    // Each routine group's SAVEd state, defined in its own file.
    struct Cand;
    struct Sync;
    struct Metrics;
    struct Wave;

private:
    // Made, freed and measured by the group's file (make_cand, del_cand,
    // cand_bytes in ft4_candidates.cpp, and so on), so that this header
    // never needs their definitions.
    static Cand* make_cand();
    static void del_cand(Cand*);
    static size_t cand_bytes(const Cand*);
    void sync4d_init();
    static Sync* make_sync();
    static void del_sync(Sync*);
    static size_t sync_bytes(const Sync*);
    static Metrics* make_metrics();
    static void del_metrics(Metrics*);
    static size_t metrics_bytes(const Metrics*);
    static Wave* make_wave();
    static void del_wave(Wave*);
    static size_t wave_bytes(const Wave*);

    Cand* cand_ = nullptr;
    Sync* sync_ = nullptr;
    Metrics* metrics_ = nullptr;
    Wave* wave_ = nullptr;

    Trace trace_;
};

}  // namespace ft4
