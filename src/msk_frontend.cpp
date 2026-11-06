// MSK144's front end (msk_frontend.h): WSJT-X's analytic.f90, tweak1.f90 and
// msk144_freq_search.f90 at 567ad29ce, line for line, with four2a's plans.
// Arrays keep the Fortran's order; the 1-based loops of the Fortran are
// shifted to 0-based here with the Fortran's index in a comment where it helps.
//
// Joe Taylor, K1JT, Steven Franke, K9AN, and the WSJT-X team (GNU GPL v3),
// and so is this.
#include "msk_frontend.h"

#include <fftw3.h>

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <mutex>
#include <new>
#include <stdexcept>
#include <utility>
#include <vector>

#include "fft.h"  // fft_planner_lock, the decoders' wisdom and plan flags
#include "ft8_maths.h"

namespace msk144 {

AnalyticTrace analytic_trace = nullptr;

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

void tweak1(const cfloat* ca, int jz, float f0, cfloat* cb) {
    constexpr double twopi = 6.283185307179586476925286766559;
    // w=1.d0; dphi (a default real) = twopi*f0/12000.d0
    double wr = 1.0, wi = 0.0;
    const float dphi = static_cast<float>(twopi * static_cast<double>(f0) / 12000.0);
    float sn, cs;
    ft8::sincosf_glibc(dphi, &sn, &cs);
    const double sr = cs, si = sn;  // wstep=cmplx(cos(dphi),sin(dphi)), a complex*16
    for (int i = 0; i < jz; ++i) {
        const double nr = wr * sr - wi * si;
        const double ni = wr * si + wi * sr;
        wr = nr;
        wi = ni;
        const cfloat w4(static_cast<float>(wr), static_cast<float>(wi));
        cb[i] = cmul(w4, ca[i]);
    }
}

void msk144_freq_search(const cfloat* cdat, float fc, int if1, int if2, float delf, int nframes, const int* navmask,
                        const cfloat* cb, cfloat* cdat2, float& xmax, float& bestf, cfloat* cs, float* xccs) {
    constexpr int NSPM = 864;
    std::vector<cfloat> c(NSPM);
    cfloat ct2[2 * NSPM];
    cfloat cc[NSPM];
    float xcc[NSPM];

    xmax = 0.0f;
    bestf = 0.0f;
    for (int i = 0; i < NSPM; ++i) cs[i] = cfloat(0.0f, 0.0f);
    for (int i = 0; i < NSPM; ++i) xccs[i] = 0.0f;

    int navg = 0;
    for (int i = 0; i < nframes; ++i) navg += navmask[i];
    if (navg <= 0 || if1 > if2) return;
    const int n = nframes * NSPM;
    const float fac = 1.0f / (48.0f * std::sqrt(static_cast<float>(navg)));

    for (int ifr = if1; ifr <= if2; ++ifr) {  // Find freq that maximizes sync
        const float ferr = static_cast<float>(ifr) * delf;
        tweak1(cdat, n, -(fc + ferr), cdat2);
        for (int i = 0; i < NSPM; ++i) c[i] = cfloat(0.0f, 0.0f);
        for (int i = 1; i <= nframes; ++i) {
            const int ib = (i - 1) * NSPM;  // ib=(i-1)*NSPM+1, 0-based here
            if (navmask[i - 1] == 1)
                for (int j = 0; j < NSPM; ++j) c[j] = cfloat(c[j].real() + cdat2[ib + j].real(), c[j].imag() + cdat2[ib + j].imag());
        }

        for (int i = 0; i < NSPM; ++i) cc[i] = cfloat(0.0f, 0.0f);
        for (int i = 0; i < NSPM; ++i) ct2[i] = c[i];
        for (int i = 0; i < NSPM; ++i) ct2[NSPM + i] = c[i];

        for (int ish = 0; ish <= NSPM - 1; ++ish) {
            // dot_product(ct2(1+ish:42+ish)+ct2(337+ish:378+ish), cb(1:42)): the sum of the
            // conjugate of the first times the second
            float sr = 0.0f, si = 0.0f;
            for (int kk = 0; kk < 42; ++kk) {
                const float tr = ct2[ish + kk].real() + ct2[336 + ish + kk].real();
                const float ti = ct2[ish + kk].imag() + ct2[336 + ish + kk].imag();
                const float br = cb[kk].real(), bi = cb[kk].imag();
                sr = sr + (tr * br - (-ti) * bi);
                si = si + (tr * bi + (-ti) * br);
            }
            cc[ish] = cfloat(sr, si);
        }

        bool finite = true;
        for (int i = 0; i < NSPM; ++i) {
            xcc[i] = ft8::hypotf_glibc(cc[i].real(), cc[i].imag());  // abs of a complex: glibc's cabsf
            if (!std::isfinite(xcc[i])) finite = false;
        }
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

}  // namespace msk144
