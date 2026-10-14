// FT2's front end: what the FT2 hybrid's
// ft4_decode (tests/data/ft2/README.md) does from a period's stretched
// samples to each candidate's soft bits and back. FT2 is FT4's decoder
// (567ad29ce's, ft4_frontend.h) run on audio stretched 2x, at the sizes
// WSJT-X improved 3.2.0's ft4_params.f90 gives it for FT2: NMAX 23*3456 =
// 79,488 ("sized for FT2 stretched"), NHSYM pinned at 122. So:
//
//  - getcandidates4 (with ft4_baseline), get_ft4_bitmetrics, subtractft4
//    (with gen_ft4wave) and get_ft4_tones_from_77bits are FT4's own,
//    called unchanged through an ft4::Frontend: getcandidates4 reads only
//    dd(1:72,000) (NHSYM 122), the bit metrics know nothing of NMAX, and
//    subtractft4 has its own NMAX 21*3456 in both trees (FT2's stretched
//    samples past 72,576 are never subtracted, as there). Only
//    getcandidates4's trace line S hashes dd over NMAX: here over 79,488.
//  - ft4_downsample and sync4d have FT4's text but run at FT2's sizes (the
//    big r2c FFT 79,488 points, df 12000/79,488, iwt/iwf/iws 69/552/138,
//    NFFT2 = NP = 4,416), so they are FT2's own copies here
//    (src/ft2_sync.cpp), not edits of src/ft4_sync.cpp: FT4's validated
//    code is not touched. sync4d_many and
//    sync4d_prepare are FT4's batched sync search (AVX2, NEON),
//    copied with NP 4,416; sync4d itself is kept, line for line, as the
//    reference they are checked against ([ft2-fe-reference]).
//
// Everything else as ft4_frontend.h says: no fused multiply-adds
// (FT8_FP_OPTIONS), the reference's float order, glibc 2.39's maths
// carried, FFTs planned as four2a plans them.
//
// getcandidates4.f90, ft4_baseline.f90, ft4_downsample.f90, sync4d.f90,
// get_ft4_bitmetrics.f90, gen_ft4wave.f90, subtractft4.f90, ft4_params.f90:
// Joe Taylor, K1JT, Steven Franke, K9AN, and the WSJT-X developers (GNU GPL
// v3); FT2's sizes from WSJT-X improved 3.2.0 (Uwe Risse, DG2YCB, and the
// same); and so is this.
#pragma once

#include <complex>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <string>

#include "ft4_frontend.h"

namespace ft2 {

using cfloat = std::complex<float>;
using Trace = ft4::Trace;

// FT4's (ft4_params.f90, the same in both trees)
using ft4::KK;
using ft4::MAXCAND;
using ft4::ND;
using ft4::NDOWN;
using ft4::NFFT1;
using ft4::NH1;
using ft4::NN;
using ft4::NS;
using ft4::NSPS;
using ft4::NSS;

// improved 3.2.0's ft4_params.f90 (FT2's sizes)
constexpr int NMAX = 23 * 3456;          // 79,488: dd, the big FFT
constexpr int NHSYM = 122;               // pinned (FT4's (NMAX-NFFT1)/NSTEP at 72,576)
constexpr int NDMAX = NMAX / NDOWN;      // 4,416: the downsampled period (NP, NFFT2)
static_assert(NHSYM == (ft4::NMAX - NFFT1) / ft4::NSTEP, "getcandidates4 is FT4's");

// improved's decoder.f90 FT2 branch
constexpr int FT2_NMAX = 23 * 1728;  // 39,744: the period's samples stretched
constexpr int FT2_STRETCH = 2;
static_assert(FT2_STRETCH * FT2_NMAX == NMAX, "the stretch fills dd");

class Frontend {
public:
    Frontend();
    ~Frontend();
    Frontend(const Frontend&) = delete;
    Frontend& operator=(const Frontend&) = delete;

    // FT4's getcandidates4 (ft4_frontend.h) on dd (NMAX here; it reads
    // the first 72,000), its S trace line's hash taken over NMAX.
    void getcandidates4(const float* dd, float fa, float fb, float syncmin, int nfqso, int maxcand, float* savg,
                        float* candidate, int& ncand, float* sbase);

    // ---------------------------------------------- src/ft2_sync.cpp
    // ft4_downsample at FT2's sizes: dd NMAX, c NDMAX; the big r2c FFT
    // (SAVEd) recomputed when newdata, which is not cleared here.
    void ft4_downsample(const float* dd, bool newdata, float f0, cfloat* c);

    // sync4d at NP = NDMAX: cd0 NDMAX, ctwk 2*NSS (ft4_frontend.h's).
    float sync4d(const cfloat* cd0, int i0, const cfloat* ctwk, int itwk);

    // The sync search's calls batched, bit for bit (ft4_frontend.h's
    // contract: sync4d_prepare(cd0) after cd0's last change).
    void sync4d_prepare(const cfloat* cd0);
    void sync4d_many(const cfloat* cd0, int istart0, int n, int step, const cfloat* ctwk, float* sync);

    // ---------------------------------------------- FT4's, unchanged
    void get_ft4_bitmetrics(const cfloat* cd, float* bitmetrics, bool& badsync) {
        fe4_.get_ft4_bitmetrics(cd, bitmetrics, badsync);
    }
    // subtractft4 on dd's first 72,576 samples (its own NMAX).
    void subtractft4(float* dd, const int itone[NN], float f0, float dt) { fe4_.subtractft4(dd, itone, f0, dt); }

    // Its buffers on the heap, bytes (FFTW's plans not counted).
    size_t heap_bytes() const;

    // The instrumented hybrid's front-end lines (tests/data/ft2/frontend,
    // tests/data/ft4/frontend/README.md's list), FT4's routines' included.
    void set_trace(Trace f);

    struct Sync;

private:
    ft4::Frontend fe4_;
    Sync* sync_ = nullptr;
    const float* dd_ = nullptr;  // getcandidates4's dd, for its S line
    Trace trace_;
    void sync4d_init();
};

}  // namespace ft2
