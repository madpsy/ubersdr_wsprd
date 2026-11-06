// MSK144's front end: WSJT-X's analytic (lib/analytic.f90), tweak1
// (lib/tweak1.f90) and msk144_freq_search (lib/msk144_freq_search.f90) at
// 567ad29ce, with the FFTs as four2a plans them, giving the reference's
// numbers to the bit on x86-64.
//
// What it rests on: no multiply-adds fused (-ffp-contract=off,
// CMakeLists.txt); complex products written out as gfortran does (Fortran
// rules, no NaN fix-up); float sums in the Fortran's order; glibc's maths
// where gfortran calls it (ft8_maths.h: cosf and sinf, sincosf for a cos/sin
// pair, hypotf for a complex abs); FFTW in place with FFTW_ESTIMATE_PATIENT
// (jt9's npatience 1), as four2a plans, forward and backward.
//
// analytic.f90, tweak1.f90, msk144_freq_search.f90, four2a.f90: Joe Taylor,
// K1JT, Steven Franke, K9AN, and the WSJT-X team (GNU GPL v3), and so is
// this.
#pragma once

#include <array>
#include <complex>

namespace msk144 {

using cfloat = std::complex<float>;

// four2a(a, nfft, 1, isign, 1): a complex transform of nfft points, in place,
// isign -1 forward and +1 backward, both unnormalised; planned as four2a
// plans (one plan per size, direction and buffer alignment, kept).
void fft_c2c_inplace(cfloat* a, int nfft, int isign);

// analytic(d, npts, nfft, c, pc, beq): the analytic signal of the real
// passband d(1:npts), zero-padded to nfft, in c(1:nfft). pc holds the five
// static phase coefficients; beq turns the static equalizer on. Like the
// Fortran's, it keeps its filter response (and the equalizer's table) from
// call to call, remade when nfft (pc) changes: not thread-safe.
void analytic(const float* d, int npts, int nfft, cfloat* c, std::array<double, 5>& pc, bool beq);

// tweak1(ca, jz, f0, cb): ca(1:jz) shifted by f0 Hz at 12000 Hz sampling
// into cb(1:jz). Not in place.
void tweak1(const cfloat* ca, int jz, float f0, cfloat* cb);

// msk144_freq_search(cdat, fc, if1, if2, delf, nframes, navmask, cb, cdat2,
// xmax, bestf, cs, xccs), with NSPM = 864: cdat and cdat2 are 864*nframes
// points, navmask has nframes entries (1: the frame is averaged), cb the 42
// points of the sync word, cs 864 points, xccs 864 values (the Fortran's
// 0:863). Outputs are as the Fortran leaves them, early returns included.
void msk144_freq_search(const cfloat* cdat, float fc, int if1, int if2, float delf, int nframes, const int* navmask,
                        const cfloat* cb, cfloat* cdat2, float& xmax, float& bestf, cfloat* cs, float* xccs);

// A debug-only trace of analytic, for the piece check: null in production
// (two untaken branches a call). When set, it is called with stage 0 just
// before the filter is applied (h and corr for 1..nfft/2+1, c(1:nfft) after
// the forward FFT), stage 1 just before the half of the DC term (c after
// the filter) and stage 2 just after it (c(1:2)), as scripts/msk144-frontend-oracle.f90's copy of analytic
// prints them.
using AnalyticTrace = void (*)(int stage, const float* h, const cfloat* corr, const cfloat* c, int nfft);
extern AnalyticTrace analytic_trace;

}  // namespace msk144
