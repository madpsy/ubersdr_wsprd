// FT8's front end: what jt9's ft8_decode does from a period's samples to
// each candidate's soft bits and back to the samples with a decode taken
// out (WSJT-X lib/ft8 at 567ad29ce: sync8 with get_spectrum_baseline,
// baseline, nuttal_window, pctile, shell and indexx; ft8_downsample;
// sync8d; ft8b, its fine sync, symbol spectra, bit metrics, normalizebmet,
// a-priori set-up, decode174_91 calls, checks and SNR; subtractft8 with
// gen_ft8wave, gfsk_pulse and peakup; genft8's get_ft8_tones_from_77bits;
// ft8_a7.f90's ft8_a7d), giving the reference's numbers to the bit on
// x86-64. ft8_decode itself (the passes, a7's tables, the 41/47/50
// sequence) is the driver's (ft8.h); it calls these
// as ft8_decode calls them, with one Frontend as one jt9 process's SAVEd
// state.
//
// What it rests on: no multiply-adds fused (-ffp-contract=off,
// CMakeLists.txt); float sums in the Fortran's order; real*8 where the
// Fortran has it (baseline and polyfit); glibc's maths where gfortran calls
// it (sincosf for cos/sin pairs, hypotf as cabsf for a complex abs, lroundf
// for nint, log10f, powf, erff), and glibc's vector maths (libmvec) on the
// elements the reference's -O3 build sends through it (from its
// disassembly): baseline's log10f four at a time from ia, the last
// (ib-ia+1) mod 4 scalar; nuttal_window's three cosf; ft8_downsample's taper
// (1..100; taper(0) folded to 1); subtractft8's window (j = -1999..2000; its
// j = -2000 folded by gcc, bits 0x2709aded); gen_ft8wave's table (entries
// 1..65532 by libmvec cosf and sinf, 0 and 65533..65535 by sincosf) and its
// ramps. erff, log10f and those libmvec functions as glibc 2.39 has them on
// x86-64, carried (ft8_maths.h), so that every platform gets their bits
// (the Pi's glibc and AArch64's libmvec differ in the last bits; its
// reference build sends the same elements through _ZGVnN4v_*).
// FFTs as four2a plans them (FFTW_ESTIMATE_PATIENT, in place;
// fst4w_frontend.h's cached plans), on buffers aligned as the reference's:
// subtractft8's refine-dt r2c runs on an array gfortran places 4 bytes past
// a 16-byte boundary, which FFTW plans without SIMD (other bits), so ours
// is placed the same.
//
// sync8.f90, get_spectrum_baseline.f90, baseline.f90, ft8_downsample.f90,
// sync8d.f90, ft8b.f90, subtractft8.f90, gen_ft8wave.f90, genft8.f90,
// ft8_a7.f90, nuttal_window.f90, gfsk_pulse.f90, peakup.f90, pctile.f90,
// shell.f90, polyfit.f90, determ.f90, indexx.f90: Joe Taylor, K1JT, Steven
// Franke, K9AN, and the WSJT-X developers (GNU GPL v3), and so is this.
#pragma once

#include <array>
#include <complex>
#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "fst4w_frontend.h"  // FftwVec
#include "ft_fast_subtract.h"  // --fast-subtract's smoothing

namespace ft8 {

using cfloat = std::complex<float>;

// ft8_params.f90
constexpr int KK = 91, ND = 58, NS = 21, NN = 79;
constexpr int NSPS = 1920, NZ = NSPS * NN, NMAX = 15 * 12000;
constexpr int NFFT1 = 2 * NSPS, NH1 = NFFT1 / 2, NSTEP = NSPS / 4;
constexpr int NHSYM = NMAX / NSTEP - 3, NDOWN = 60;
constexpr int NP2 = 2812;  // ft8b's and ft8_a7d's

// What ft8b and ft8_a7d ask of packjt77 (one Packjt77 a decoder, the
// driver's): ft8b's "dxcall13=hiscall12; mycall13=mycall12", unpack77,
// stdcall and genft8's packing half.
class Codec {
public:
    virtual ~Codec() = default;
    virtual void set_calls(const std::string& mycall12, const std::string& hiscall12) = 0;
    // unpack77(c77, nrx, msg, success): c77 is 77 '0'/'1' characters;
    // msg37 its character*37.
    virtual bool unpack77(const char* c77, int nrx, std::string& msg37) = 0;
    virtual bool stdcall(const std::string& call12) = 0;
    // genft8(msg, i3, n3, msgsent, msgbits, itone) up to its "go to 2"
    // (ok) or its bad-message branch (msgbits 0).
    struct Packed {
        bool ok = false;
        std::string msgsent;  // 37 characters
        std::array<int8_t, 77> msgbits{};
    };
    virtual Packed genft8(const std::string& msg37) = 0;
};

// get_ft8_tones_from_77bits(msgbits, itone) (genft8.f90's entry).
void get_ft8_tones_from_77bits(const int8_t msgbits[77], int itone[NN]);

// normalizebmet(bmet, n).
void normalizebmet(float* bmet, int n);

// peakup(ym, y0, yp, dx).
float peakup(float ym, float y0, float yp);

// gfsk_pulse(b, t).
float gfsk_pulse(float b, float t);

// nuttal_window(win, n).
void nuttal_window(float* win, int n);

// baseline(s, nfa, nfb, sbase): s and sbase are 1920 (Fortran's s(1920),
// element i at [i-1]); s is changed (to dB from ia to ib) as there.
void baseline(float* s, int nfa, int nfb, float* sbase);

// ft8_decode's xbase for a candidate (a7: max(1, ...) as its line has it).
float xbase_of(const float* sbase, float f1, bool a7);

class Frontend {
public:
    Frontend();
    ~Frontend();
    Frontend(const Frontend&) = delete;
    Frontend& operator=(const Frontend&) = delete;

    // sync8(dd, npts, nfa, nfb, syncmin, nfqso, maxcand, candidate, ncand,
    // sbase): candidate is 3 x maxcand (Fortran's (i,k) at [(k-1)*3 + i-1]),
    // sbase NH1 ([i-1]). nfa and nfb are get_spectrum_baseline's to clamp.
    void sync8(const float* dd, int npts, int& nfa, int& nfb, float syncmin, int nfqso, int maxcand,
               float* candidate, int& ncand, float* sbase);

    // get_spectrum_baseline(dd, nfa, nfb, sbase).
    void get_spectrum_baseline(const float* dd, int& nfa, int& nfb, float* sbase);

    // ft8_downsample(dd, newdat, f0, c1): c1 is 3200; the big FFT kept here
    // (recomputed, and newdat cleared, when newdat).
    void ft8_downsample(const float* dd, bool& newdat, float f0, cfloat* c1);

    // sync8d(cd0, i0, ctwk, itwk, sync).
    float sync8d(const cfloat* cd0, int i0, const cfloat* ctwk, int itwk);

    // ft8b(dd0, newdat, nQSOProgress, nfqso, nftx, ndepth, nzhsym, lapon,
    // lapcqonly, napwid, lsubtract, nagain, ncontest, imetric, iaptype,
    // mycall12, hiscall12, f1, xdt, xbase, apsym, aph10, nharderrors, dmin,
    // nbadcrc, ipass, msg37, xsnr, itone), the arguments as there (ipass is
    // ft8b's pass, iappass in ft8_decode); packjt77 through codec.
    void ft8b(float* dd0, bool& newdat, int nQSOProgress, int nfqso, int nftx, int ndepth, int nzhsym, bool lapon,
              bool lapcqonly, int napwid, bool lsubtract, bool nagain, int ncontest, int imetric, int& iaptype,
              const std::string& mycall12, const std::string& hiscall12, float& f1, float& xdt, float xbase,
              const int apsym[58], const int aph10[10], int& nharderrors, float& dmin, int& nbadcrc, int& ipass,
              std::string& msg37, float& xsnr, int itone[NN], Codec& codec);

    // subtractft8(dd0, itone, f0, dt, lrefinedt).
    void subtractft8(float* dd0, const int itone[NN], float f0, float dt, bool lrefinedt);
    // Ours (jt9 --fast-subtract): sqf's 4001-tap smoothing of dd*conj(cref) computed at 300 S/s and interpolated back
    // (sqf_fast), in place of its two 180,000-point FFTs. Not WSJT-X's to the bit; the same decodes.
    void set_fast_subtract(bool v) { fast_subtract_ = v; }

    // gen_ft8wave(itone, nsym, nsps, bt, fsample, f0, cwave, wave, icmplx,
    // nwave): cwave or wave (icmplx) nwave long; its tables kept here.
    void gen_ft8wave(const int* itone, int nsym, int nsps, float bt, float fsample, float f0, cfloat* cwave,
                     float* wave, int icmplx, int nwave);

    // The line-for-line translations subtractft8 (with its sqf) and
    // gen_ft8wave were optimised from, bit for bit the same, kept for the
    // side-by-side test.
    void subtractft8_reference(float* dd0, const int itone[NN], float f0, float dt, bool lrefinedt);
    void gen_ft8wave_reference(const int* itone, int nsym, int nsps, float bt, float fsample, float f0,
                               cfloat* cwave, float* wave, int icmplx, int nwave);

    // ft8_a7d(dd0, newdat, call_1, call_2, grid4, xdt, f1, xbase,
    // nharderrors, dmin, msg37, xsnr): call_1 and call_2 character*12,
    // grid4 character*4.
    void ft8_a7d(const float* dd0, bool& newdat, const std::string& call_1, const std::string& call_2,
                 const std::string& grid4, float& xdt, float& f1, float xbase, int& nharderrors, float& dmin,
                 std::string& msg37, float& xsnr, Codec& codec);

    // Its buffers on the heap now, bytes (FFTW's plans
    // not counted). Not while it runs.
    size_t heap_bytes() const;

    // The insides as scripts/ft8-frontend-patch.py's instrumented reference
    // writes them (tests/test_ft8_frontend.cpp): one call a line, none
    // when unset.
    void set_trace(std::function<void(const std::string&)> f) { trace_ = std::move(f); }
    // The trace's pending chain of candidates that failed the sync check,
    // written out (its "g" line), as the instrumented reference does at
    // ft8_decode's start and end.
    void trace_flush();

private:
    // ft8b's and ft8_a7d's shared steps: the fine sync from the coarse
    // start (ft8_downsample twice, sync8d) and the symbol spectra.
    struct Sync {
        int ibest0 = 0, ibest = 0, nsync = 0;
        float delfbest = 0, smax = 0;
        float ss[9] = {};
        uint32_t h_cd0 = 0;  // the first downsample's, for the trace
    };
    void fine_sync(const float* dd0, bool& newdat, float& f1, float& xdt, Sync& sy, bool a7);
    // Bit metrics nsym 1..3 into bmeta..bmetd (and bmete when given).
    void bit_metrics(int imetric, float* bmeta, float* bmetb, float* bmetc, float* bmetd, float* bmete);
    // The trace's G line, or the candidate's words into the chain.
    void trace_sync(const Sync& sy, float xdt, bool chain, float f1_in, float xdt_in, float xbase, bool newdat_in);
    // The trace's M (sampled), before normalizebmet, and L/l lines after.
    void trace_metrics(int nll, float* const* v, bool after);
    // The trace's r line, at ft8b's returns past the sync check (xsnr null:
    // every pass failed, not set).
    void trace_ret(int nbadcrc, int ipass, int iaptype, int nharderrors, float f1, float xdt, float dmin,
                   const float* xsnr, const std::string& msg37);

    float sqf(int idt, float* dd0, float f0, float dt, bool ldt);
    float sqf_reference(int idt, float* dd0, float f0, float dt, bool ldt);
    void st_make_filter();       // the filter and end correction, first call
    void st_plan_r2c(float* x);  // the trial's r2c on x, made once
    void gw_tables(int nsym, int nsps, float bt, float fsample);  // pulse, ctab

    // get_spectrum_baseline's
    bool gsb_first_ = true;
    std::vector<float> gsb_window_;
    // sync8's spectra s(NH1, NHSYM)
    std::vector<float> s_;
    fst4w::FftwVec<float> x3840_;
    // ft8_downsample's
    bool ds_first_ = true;
    float taper_[101] = {};
    fst4w::FftwVec<float> ds_x_;  // x(NFFT1+2), cx equivalenced
    fst4w::FftwVec<cfloat> c3200_;
    // sync8d's
    bool s8d_first_ = true;
    float s8d_twopi_ = 0;
    cfloat csync_[7][32];
    // ft8b's and ft8_a7d's (cd0, symbol spectra)
    fst4w::FftwVec<cfloat> cd0_, csymb_;
    cfloat cs_[NN][8];
    float s8_[NN][8];
    // ft8b's a-priori tables (DATA, SAVEd), mapped on its first call and
    // on each change of ncontest
    bool ft8b_first_ = true;
    int ncontest0_ = 0;
    int mcq_[29] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0};
    int mcqru_[29] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0};
    int mcqfd_[29] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0};
    int mcqtest_[29] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0};
    int mcqww_[29] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0};
    int mrrr_[19] = {0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1};
    int m73_[19] = {0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1};
    int mrr73_[19] = {0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1};
    // subtractft8's
    bool st_first_ = true;
    float endcorrection_[2001] = {};
    bool fast_subtract_ = false;
    std::unique_ptr<FastSubtract> fs_;  // fast subtraction's (sqf_fast), made on its first call
    std::vector<float> fs_gain_;         // and the end correction a sample
    float sqf_fast(int idt, float* dd0, float f0, float dt, bool ldt);
    fst4w::FftwVec<cfloat> cref_, camp_, cfilt_, cw_;
    std::vector<float> st_dd_;
    fst4w::FftwVec<float> st_x_;
    void* st_plan_ = nullptr;  // its r2c (fftwf_plan), for st_x_ as sqf places it
    // gen_ft8wave's
    float fchk0_ = 0.0f, gw_twopi_ = 0, gw_dt_ = 0, gw_hmod_ = 0;
    std::vector<float> pulse_ = std::vector<float>(23040);
    std::vector<cfloat> ctab_ = std::vector<cfloat>(65536);
    std::vector<float> dphi_;
    bool ctab_made_ = false;      // ctab_ depends on nothing but twopi
    int gw_npulse_ = 0;           // the pulse_ (re)computations, for pdphi_
    int pdphi_npulse_ = -1, pdphi_nsps_ = 0;
    std::vector<float> pdphi_;    // dphi_peak*pulse(m), m < 3*nsps
    std::vector<float> gw_phi_;   // phi before each output sample
    std::vector<int> gw_ix_;      // its table index
    std::vector<float> gw_c_;     // the envelope ramp

    std::function<void(const std::string&)> trace_;
    int trace_nl_ = 0;     // the trace's M count (every 499th sampled)
    bool trace_smp_ = false;
    int trace_nch_ = 0;    // the pending chain
    uint32_t trace_ch_ = 2166136261u;
};

}  // namespace ft8
