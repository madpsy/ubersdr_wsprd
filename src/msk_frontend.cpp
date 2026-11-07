// MSK144's front end (msk_frontend.h): WSJT-X's analytic.f90, tweak1.f90 and
// msk144_freq_search.f90 at 567ad29ce, line for line, with four2a's plans.
// Arrays keep the Fortran's order; the 1-based loops of the Fortran are
// shifted to 0-based here with the Fortran's index in a comment where it helps.
//
// Joe Taylor, K1JT, Steven Franke, K9AN, and the WSJT-X team (GNU GPL v3),
// and so is this.
#include "msk_frontend.h"

#include <fftw3.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <mutex>
#include <new>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#include "fft.h"  // fft_planner_lock, the decoders' wisdom and plan flags
#include "ft8_maths.h"

namespace msk144 {

AnalyticTrace analytic_trace = nullptr;
FreqSearchTrace freq_search_trace = nullptr;

namespace {

constexpr int NSMALL = 16385;  // four2a: the largest "small" FFT, whose data a plan keeps

// four2a's planning: in place, FFTW_ESTIMATE_PATIENT (jt9's npatience 1),
// with the decoders' wisdom (fft.h), the process's set aside and put back
// after. Planning at this patience reads no data, but four2a keeps up to
// NSMALL values across it anyway: so do we.
template <typename Plan>
fftwf_plan make_plan(void* data, size_t keep_bytes, Plan make) {
    std::lock_guard<std::mutex> lock(skimmer::fft_planner_lock());
    char* const saved = fftwf_export_wisdom_to_string();
    fftwf_forget_wisdom();
    skimmer::decoder_wisdom_before_plan();
    std::vector<char> aa;
    if (keep_bytes) aa.assign(static_cast<char*>(data), static_cast<char*>(data) + keep_bytes);
    fftwf_plan p = make();
    if (keep_bytes) std::memcpy(data, aa.data(), keep_bytes);
    skimmer::decoder_wisdom_after_plan();
    fftwf_forget_wisdom();
    if (saved) {
        fftwf_import_wisdom_from_string(saved);
        std::free(saved);
    }
    if (!p) throw std::runtime_error("msk144: FFTW could not plan");
    return p;
}

// One plan per (n, direction) for the process, made on an fftwf_malloc
// buffer of its own (freed once planned) and keeping that buffer's alignment.
struct CachedPlan {
    fftwf_plan plan = nullptr;
    int alignment = 0;
};

CachedPlan cached_plan(int n, int isign) {
    static std::mutex mu;
    static std::vector<std::pair<std::pair<int, int>, CachedPlan>> cache;
    std::lock_guard<std::mutex> lock(mu);
    for (const auto& e : cache)
        if (e.first.first == n && e.first.second == isign) return e.second;
    const size_t bytes = sizeof(cfloat) * static_cast<size_t>(n);
    void* buf = fftwf_malloc(bytes);
    if (!buf) throw std::bad_alloc();
    std::memset(buf, 0, bytes);
    auto* b = static_cast<fftwf_complex*>(buf);
    CachedPlan c;
    // the kept plans' patience: four2a's FFTW_ESTIMATE_PATIENT unless jt9 -w asked for another
    const unsigned flags = skimmer::decoder_plan_flags();
    c.plan = make_plan(buf, 0, [&] { return fftwf_plan_dft_1d(n, b, b, isign < 0 ? FFTW_FORWARD : FFTW_BACKWARD, flags); });
    c.alignment = fftwf_alignment_of(static_cast<float*>(buf));
    fftwf_free(buf);
    cache.push_back({{n, isign}, c});
    return c;
}

}  // namespace

void fft_c2c_inplace(cfloat* a, int nfft, int isign) {
    auto* d = reinterpret_cast<fftwf_complex*>(a);
    const CachedPlan c = cached_plan(nfft, isign);
    if (fftwf_alignment_of(reinterpret_cast<float*>(a)) == c.alignment) {
        fftwf_execute_dft(c.plan, d, d);
        return;
    }
    // a buffer aligned otherwise than the kept plan's: a plan made on it, used once
    const size_t keep = nfft <= NSMALL ? sizeof(cfloat) * static_cast<size_t>(nfft) : 0;
    fftwf_plan p = make_plan(a, keep, [&] {
        return fftwf_plan_dft_1d(nfft, d, d, isign < 0 ? FFTW_FORWARD : FFTW_BACKWARD, FFTW_ESTIMATE_PATIENT);
    });
    fftwf_execute(p);
    std::lock_guard<std::mutex> lock(skimmer::fft_planner_lock());
    fftwf_destroy_plan(p);
}

namespace {

// complex products as gfortran writes them (-fcx-fortran-rules): the plain
// formula, no NaN recovery.
inline cfloat cmul(cfloat a, cfloat b) {
    const float ar = a.real(), ai = a.imag(), br = b.real(), bi = b.imag();
    return cfloat(ar * br - ai * bi, ar * bi + ai * br);
}

constexpr int NFFTMAX = 1024 * 1024;

}  // namespace

void analytic(const float* d, int npts, int nfft, cfloat* c, std::array<double, 5>& pc, bool beq) {
    // save corr, nfft0, h, ac, aclast, pclast, pi, t, beta
    static std::vector<float> h(NFFTMAX / 2 + 1);          // real h(NFFTMAX/2), 1-based: h[i-1]
    static std::vector<cfloat> corr(NFFTMAX / 2 + 1);      // complex corr(NFFTMAX/2)
    static int nfft0 = 0;
    static double aclast[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
    static double pclast[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
    static const double ac[5] = {1.0, 0.0, 0.0, 0.0, 0.0};  // data ac/1.0,0.0,0.0,0.0,0.0/
    static float pi, t, beta;

    const float df = 12000.0f / static_cast<float>(nfft);
    const int nh = nfft / 2;
    if (nfft != nfft0) {
        // pi=4.0*atan(1.0): gfortran folds atan(1.0) to the float nearest pi/4
        pi = 4.0f * 0.785398185253143310546875f;
        t = 1.0f / 2000.0f;
        beta = 0.1f;
        for (int i = 1; i <= nh + 1; ++i) {
            const float ff = static_cast<float>(i - 1) * df;
            const float f = ff - 1500.0f;
            h[i - 1] = 1.0f;
            if (std::fabs(f) > (1.0f - beta) / (2.0f * t) && std::fabs(f) <= (1.0f + beta) / (2.0f * t)) {
                h[i - 1] = h[i - 1] * 0.5f *
                           (1.0f + ft8::cosf_glibc((pi * t / beta) * (std::fabs(f) - (1.0f - beta) / (2.0f * t))));
            } else if (std::fabs(f) > (1.0f + beta) / (2.0f * t)) {
                h[i - 1] = 0.0f;
            }
        }
        nfft0 = nfft;
    }

    bool changed = false;
    for (int i = 0; i < 5; ++i)
        if (aclast[i] != ac[i] || pclast[i] != pc[i]) changed = true;
    if (changed) {
        for (int i = 0; i < 5; ++i) {
            aclast[i] = ac[i];
            pclast[i] = pc[i];
        }
        // The reference's -O3 vectorised this loop: four elements at a time from the first, glibc's
        // _ZGVbN4v_cosf and _ZGVbN4v_sinf for cos(pd) and sin(pd), and the last 1-3 elements by sincosf.
        auto make_corr = [&](int i, float& pd) {  // i is the Fortran's index; returns pd, sets corr's real part
            const float ff = static_cast<float>(i - 1) * df;
            const float f = ff - 1500.0f;
            const double fp = f / 1000.0f;  // f/1000.0 in single, then widened
            corr[i - 1] = cfloat(static_cast<float>(ac[0] + fp * (ac[1] + fp * (ac[2] + fp * (ac[3] + fp * ac[4])))), 0.0f);
            pd = static_cast<float>(fp * fp * (pc[2] + fp * (pc[3] + fp * pc[4])));  // ignore 1st two terms
        };
        const int ncorr = nh + 1;
        const int nvec = ncorr - ncorr % 4;
        for (int i0 = 1; i0 <= nvec; i0 += 4) {
            float pd[4], cs4[4], sn4[4];
            for (int l = 0; l < 4; ++l) make_corr(i0 + l, pd[l]);
            ft8::cosf4_glibc(pd, cs4);
            ft8::sinf4_glibc(pd, sn4);
            for (int l = 0; l < 4; ++l) corr[i0 + l - 1] = cmul(corr[i0 + l - 1], cfloat(cs4[l], sn4[l]));
        }
        for (int i = nvec + 1; i <= ncorr; ++i) {
            float pd, sn, cs;
            make_corr(i, pd);
            ft8::sincosf_glibc(pd, &sn, &cs);
            corr[i - 1] = cmul(corr[i - 1], cfloat(cs, sn));
        }
    }

    const float fac = 2.0f / static_cast<float>(nfft);
    for (int i = 0; i < npts; ++i) c[i] = cfloat(fac * d[i], 0.0f);
    for (int i = npts; i < nfft; ++i) c[i] = cfloat(0.0f, 0.0f);
    fft_c2c_inplace(c, nfft, -1);  // Forward c2c FFT

    if (analytic_trace) analytic_trace(0, h.data(), corr.data(), c, nfft);
    if (beq) {
        for (int i = 0; i <= nh; ++i) {
            // h*corr*c = (h*corr)*c: the real h made complex (h, 0) for each product, as gfortran
            // does (signed zeros show it where h is 0)
            c[i] = cmul(cmul(cfloat(h[i], 0.0f), corr[i]), c[i]);
        }
    } else {
        for (int i = 0; i <= nh; ++i) c[i] = cmul(cfloat(h[i], 0.0f), c[i]);
    }

    if (analytic_trace) analytic_trace(1, h.data(), corr.data(), c, nfft);
    // c(1)=0.5*c(1): the real 0.5 made complex (0.5, 0) for a full complex product, as gfortran does
    // (signed zeros show it)
    c[0] = cmul(cfloat(0.5f, 0.0f), c[0]);  // Half of DC term
    if (analytic_trace) analytic_trace(2, h.data(), corr.data(), c, nfft);
    for (int i = nh + 1; i < nfft; ++i) c[i] = cfloat(0.0f, 0.0f);  // Zero the negative frequencies
    fft_c2c_inplace(c, nfft, 1);                                    // Inverse c2c FFT
}

namespace {

// tweak1's step: wstep=cmplx(cos(dphi),sin(dphi)), a complex*16, with dphi (a default real) =
// twopi*f0/12000.d0. tweak1 and msk144_freq_search share it.
void tweak1_wstep(float f0, double& sr, double& si) {
    constexpr double twopi = 6.283185307179586476925286766559;
    const float dphi = static_cast<float>(twopi * static_cast<double>(f0) / 12000.0);
    float sn, cs;
    ft8::sincosf_glibc(dphi, &sn, &cs);
    sr = cs;
    si = sn;
}

// Four lanes of doubles and of floats (GCC's and Clang's vector extensions: each operation
// acts lane by lane, with the scalar operation's IEEE rounding).
typedef double v4d __attribute__((vector_size(32)));
typedef float v4f __attribute__((vector_size(16)));

}  // namespace

void tweak1(const cfloat* ca, int jz, float f0, cfloat* cb) {
    double wr = 1.0, wi = 0.0;  // w=1.d0
    double sr, si;
    tweak1_wstep(f0, sr, si);
    for (int i = 0; i < jz; ++i) {
        const double nr = wr * sr - wi * si;
        const double ni = wr * si + wi * sr;
        wr = nr;
        wi = ni;
        const cfloat w4(static_cast<float>(wr), static_cast<float>(wi));
        cb[i] = cmul(w4, ca[i]);
    }
}

// msk144_freq_search, in two places in another order than the Fortran's, with the same bits.
//
// 1. The frequency shift. The Fortran calls tweak1 once a trial frequency, then sums the
// averaged frames. tweak1's phasor w steps by a complex*16 product each sample, a chain of
// dependent multiplies and adds that keeps the CPU waiting on each result. Here NL trials
// run side by side: one pass over the samples steps NL phasors, each with its own wstep and
// the same double expressions as tweak1, so each lane's w takes the same values in the same
// order, and multiplies cdat by each lane's w4 with the same complex product. The NL chains
// do not depend on each other, so the CPU (and the compiler, across vector lanes) overlaps
// them. The frames go into each trial's c as the Fortran adds them: from 0, frame by frame in
// ascending order. Only the last trial's shifted data reaches cdat2, as tweak1's last call
// leaves it in the Fortran. Lanes past if2 in the last group repeat if2 and are dropped.
//
// 2. The correlation. The Fortran runs lag by lag: for each ish, a dot product over the 42
// sync points kk. Here kk runs outside and ish inside, over the arrays sr and si, so the
// compiler keeps eight lags in each AVX2 register (four in NEON's) and steps them together. Every lag stays its
// own sum, from 0, over kk = 0..41 in ascending order, with the same expression for each term.
// dr and di hold ct2(1+j)+ct2(337+j) for each j: the same single add the Fortran makes inside
// the dot product, made once a trial instead of 42 times.
//
// Why both are exact: the build forbids contraction (FT8_FP_OPTIONS: -ffp-contract=off, and
// -mno-fma on GCC before 13), so each + and * rounds on its own, as the Fortran's does; IEEE
// rounding depends only on the operands, so a value computed from the same operands in the
// same order has the same bits whichever lane, register or loop order computes it. What moved
// is the order across independent trials and lags, and no trial or lag reads another's
// values. The trials' results then go through the peak search in ascending ifr, as in the
// Fortran. The magnitudes stay glibc's cabsf: hypotf_n_glibc gives hypotf_glibc's bits, four at a
// time, with no squared-magnitude shortcut, so the peak search sees the same values.
void msk144_freq_search(const cfloat* cdat, float fc, int if1, int if2, float delf, int nframes, const int* navmask,
                        const cfloat* cb, cfloat* cdat2, float& xmax, float& bestf, cfloat* cs, float* xccs) {
    constexpr int NSPM = 864;
    constexpr int NSYNC = 42;
    constexpr int ND = NSPM + NSYNC - 1;  // ct2(1+ish+kk), ish 0..863, kk 0..41: 905 sums
    constexpr int NL = 8;                 // trials side by side
    constexpr int NG = NL / 4;            // groups of four trials, one vector of each
    std::vector<cfloat> c(NSPM);
    std::vector<v4f> clr(static_cast<size_t>(NSPM) * NG), cli(clr.size());  // [j*NG+l/4][l%4]: trial l's c(j)
    cfloat ct2[2 * NSPM];
    alignas(32) float dr[ND], di[ND];
    alignas(32) float sr[NSPM], si[NSPM];
    float xcc[NSPM];

    xmax = 0.0f;
    bestf = 0.0f;
    for (int i = 0; i < NSPM; ++i) cs[i] = cfloat(0.0f, 0.0f);
    for (int i = 0; i < NSPM; ++i) xccs[i] = 0.0f;

    int navg = 0;
    for (int i = 0; i < nframes; ++i) navg += navmask[i];
    if (navg <= 0 || if1 > if2) return;
    const float fac = 1.0f / (48.0f * std::sqrt(static_cast<float>(navg)));

    for (int ifr0 = if1; ifr0 <= if2; ifr0 += NL) {  // Find freq that maximizes sync
        const int nl = std::min(NL, if2 - ifr0 + 1);    // trials in this group
        const bool last = ifr0 + nl - 1 == if2;         // the group that holds if2: it fills cdat2

        // tweak1(cdat, n, -(fc+ferr), cdat2) for each trial of the group, and the frames' sums;
        // trial l is lane l%4 of group l/4
        v4d wsr[NG], wsi[NG], wr[NG], wi[NG];
        for (int l = 0; l < NL; ++l) {
            const float ferr = static_cast<float>(ifr0 + std::min(l, nl - 1)) * delf;
            double sr_, si_;
            tweak1_wstep(-(fc + ferr), sr_, si_);
            wsr[l / 4][l % 4] = sr_;
            wsi[l / 4][l % 4] = si_;
            wr[l / 4][l % 4] = 1.0;  // w=1.d0
            wi[l / 4][l % 4] = 0.0;
        }
        for (v4f& x : clr) x = v4f{0.0f, 0.0f, 0.0f, 0.0f};
        for (v4f& x : cli) x = v4f{0.0f, 0.0f, 0.0f, 0.0f};
        // The first ngu groups of four, all a group of nl trials uses: with four trials or
        // fewer (a short last group, or a short call), group 1 holds only repeats of if2,
        // whose sums nothing reads (the peak search reads l < nl, cdat2 comes from lane
        // nl-1, both in group 0), so it is not run; its sums stay 0, unread.
        auto frames = [&](auto ngu_c) {
            constexpr int ngu = decltype(ngu_c)::value;
            for (int fr = 0; fr < nframes; ++fr) {
                const int ib = fr * NSPM;  // ib=(i-1)*NSPM+1, 0-based here
                const bool avg = navmask[fr] == 1;
                for (int j = 0; j < NSPM; ++j) {
                    const float ar = cdat[ib + j].real(), ai = cdat[ib + j].imag();
                    const v4f arv = {ar, ar, ar, ar}, aiv = {ai, ai, ai, ai};
                    v4f xr[ngu], xi[ngu];
                    for (int g = 0; g < ngu; ++g) {
                        const v4d nr = wr[g] * wsr[g] - wi[g] * wsi[g];
                        const v4d ni = wr[g] * wsi[g] + wi[g] * wsr[g];
                        wr[g] = nr;
                        wi[g] = ni;
                        const v4f w4r = __builtin_convertvector(nr, v4f), w4i = __builtin_convertvector(ni, v4f);
                        xr[g] = w4r * arv - w4i * aiv;  // cmul(w4, cdat(i))
                        xi[g] = w4r * aiv + w4i * arv;
                    }
                    if (avg) {
                        v4f* const cr = &clr[static_cast<size_t>(j) * NG];
                        v4f* const ci = &cli[static_cast<size_t>(j) * NG];
                        for (int g = 0; g < ngu; ++g) {
                            cr[g] = cr[g] + xr[g];
                            ci[g] = ci[g] + xi[g];
                        }
                    }
                    if (last) cdat2[ib + j] = cfloat(xr[(nl - 1) / 4][(nl - 1) % 4], xi[(nl - 1) / 4][(nl - 1) % 4]);
                }
            }
        };
        static_assert(NG == 2, "frames() is run for one group or both");
        if (nl > 4)
            frames(std::integral_constant<int, 2>{});
        else
            frames(std::integral_constant<int, 1>{});

        for (int l = 0; l < nl; ++l) {
            const int ifr = ifr0 + l;
            const float ferr = static_cast<float>(ifr) * delf;
            for (int i = 0; i < NSPM; ++i) c[i] = cfloat(clr[static_cast<size_t>(i) * NG + l / 4][l % 4], cli[static_cast<size_t>(i) * NG + l / 4][l % 4]);

            for (int i = 0; i < NSPM; ++i) ct2[i] = c[i];
            for (int i = 0; i < NSPM; ++i) ct2[NSPM + i] = c[i];

            // cc(ish) = dot_product(ct2(1+ish:42+ish)+ct2(337+ish:378+ish), cb(1:42)): the sum of
            // the conjugate of the first times the second; lags inside, sync points outside (see
            // above)
            for (int j = 0; j < ND; ++j) {
                dr[j] = ct2[j].real() + ct2[336 + j].real();
                di[j] = ct2[j].imag() + ct2[336 + j].imag();
            }
            for (int ish = 0; ish < NSPM; ++ish) {
                sr[ish] = 0.0f;
                si[ish] = 0.0f;
            }
            for (int kk = 0; kk < NSYNC; ++kk) {
                const float br = cb[kk].real(), bi = cb[kk].imag();
                const float* const tr = dr + kk;
                const float* const ti = di + kk;
                for (int ish = 0; ish < NSPM; ++ish) {
                    sr[ish] = sr[ish] + (tr[ish] * br - (-ti[ish]) * bi);
                    si[ish] = si[ish] + (tr[ish] * bi + (-ti[ish]) * br);
                }
            }

            // abs of a complex: glibc's cabsf, hypotf_glibc's bits four at a time
            ft8::hypotf_n_glibc(sr, si, xcc, NSPM);
            if (freq_search_trace) freq_search_trace(ifr, c.data(), sr, si, xcc);
            bool finite = true;
            for (int i = 0; i < NSPM; ++i)
                if (!std::isfinite(xcc[i])) finite = false;
            if (!finite) continue;
            float mx = xcc[0];
            for (int i = 1; i < NSPM; ++i)
                if (xcc[i] > mx) mx = xcc[i];
            const float xb = mx * fac;
            if (!std::isfinite(xb)) continue;
            if (ifr == if1 || xb > xmax) {
                xmax = xb;
                bestf = ferr;
                for (int i = 0; i < NSPM; ++i) cs[i] = c[i];
                for (int i = 0; i < NSPM; ++i) xccs[i] = xcc[i];
            }
        }
    }
}

}  // namespace msk144
