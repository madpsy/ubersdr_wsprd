// FST4W's front end (fst4w_frontend.h): WSJT-X lib/fst4_decode.f90's
// get_candidates_fst4, fst4_downsample, fst4_sync_search and sync_fst4,
// lib/blanker.f90, lib/four2a.f90's planning, lib/fst4/fst4_baseline.f90,
// lib/pctile.f90, lib/shell.f90, lib/polyfit.f90, lib/determ.f90,
// lib/fst4/get_fst4_bitmetrics.f90 and lib/ft8/ft8b.f90's normalizebmet, at
// 567ad29ce, line for line. Arrays keep the Fortran's indexing where that
// keeps the two readable side by side.
//
// Steven Franke, K9AN, and Joe Taylor, K1JT, WSJT-X (GNU GPL v3), and so is
// this.
#include "fst4w_frontend.h"

#include <fftw3.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <mutex>
#include <stdexcept>

#include "fft.h"  // fft_planner_lock
#include "fst4w_simd.h"


namespace fst4w {

namespace {

constexpr int NN = 160;  // fst4_params.f90: sync and data symbols

// gfortran's maxloc of a real section: the first maximum; NaNs never win
// unless all are NaN (then the first element).
template <typename F>
int maxloc(int lo, int hi, F at) {
    int best = -1;
    float m = 0.0f;
    for (int i = lo; i <= hi; ++i) {
        const float v = at(i);
        if (best < 0) {
            if (!std::isnan(v)) {
                best = i;
                m = v;
            }
        } else if (v > m) {
            best = i;
            m = v;
        }
    }
    return best < 0 ? lo : best;
}

constexpr int NSMALL = 16385;
constexpr unsigned PATIENCE = FFTW_ESTIMATE_PATIENT;

// four2a's planning: in place, FFTW_ESTIMATE_PATIENT, with the decoders'
// wisdom (fft.h: none in the skimmer, jt9_wisdom.dat's in jt9), the
// process's set aside and put back after. Planning at this patience reads
// no data, but four2a keeps up to NSMALL values across it anyway: so do we.
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
    if (!p) throw std::runtime_error("fst4w: FFTW could not plan");
    return p;
}

// The plans, one per (n, kind) for the process: kind 0 r2c, -1 forward,
// +1 backward. Each is made on an fftwf_malloc buffer of its own (freed
// once planned) and keeps that buffer's alignment.
struct CachedPlan {
    fftwf_plan plan = nullptr;
    int alignment = 0;
};

CachedPlan cached_plan(int n, int kind) {
    static std::mutex mu;
    static std::vector<std::pair<std::pair<int, int>, CachedPlan>> cache;
    std::lock_guard<std::mutex> lock(mu);
    for (const auto& e : cache)
        if (e.first.first == n && e.first.second == kind) return e.second;
    const size_t bytes = kind == 0 ? sizeof(float) * 2 * static_cast<size_t>(n / 2 + 1) : sizeof(cfloat) * static_cast<size_t>(n);
    void* buf = fftwf_malloc(bytes);
    if (!buf) throw std::bad_alloc();
    std::memset(buf, 0, bytes);
    CachedPlan c;
    c.plan = make_plan(buf, 0, [&] {
        if (kind == 0) return fftwf_plan_dft_r2c_1d(n, static_cast<float*>(buf), static_cast<fftwf_complex*>(buf), PATIENCE);
        return fftwf_plan_dft_1d(n, static_cast<fftwf_complex*>(buf), static_cast<fftwf_complex*>(buf),
                                 kind < 0 ? FFTW_FORWARD : FFTW_BACKWARD, PATIENCE);
    });
    c.alignment = fftwf_alignment_of(static_cast<float*>(buf));
    fftwf_free(buf);
    cache.push_back({{n, kind}, c});
    return c;
}

// A plan for this very buffer (one not aligned as the cached plans are),
// used once.
template <typename Plan>
void plan_and_run(void* data, size_t keep_bytes, Plan make) {
    fftwf_plan p = make_plan(data, keep_bytes, make);
    fftwf_execute(p);
    std::lock_guard<std::mutex> lock(skimmer::fft_planner_lock());
    fftwf_destroy_plan(p);
}

// cabsf(x + iy) exactly as glibc's (hypotf: the double sum of the exact
// double squares, its square root rounded to float), for finite values
// (eight at a time: simd::abs_exact).
inline float abs_exact(float x, float y) {
    return static_cast<float>(std::sqrt(static_cast<double>(x) * static_cast<double>(x) +
                                        static_cast<double>(y) * static_cast<double>(y)));
}

// ---------------------------------------------------------------- baseline

// shell(n, a): Shell's sort, ascending, a 1-based.
void shell(int n, float* a) {
    int inc = 1;
    do inc = 3 * inc + 1;
    while (inc <= n);
    do {
        inc = inc / 3;
        for (int i = inc + 1; i <= n; ++i) {
            const float v = a[i];
            int j = i;
            while (a[j - inc] > v) {
                a[j] = a[j - inc];
                j = j - inc;
                if (j <= inc) break;
            }
            a[j] = v;
        }
    } while (inc > 1);
}

// pctile(x, npts, npct, xpct), x 1-based.
float pctile(const float* x, int npts, int npct) {
    if (npts < 0 || npct < 0 || npct > 100) return 1.0f;
    std::vector<float> tmp(static_cast<size_t>(npts) + 1);
    for (int i = 1; i <= npts; ++i) tmp[i] = x[i];
    shell(npts, tmp.data());
    int j = static_cast<int>(std::lround(static_cast<float>(npts) * 0.01f * static_cast<float>(npct)));
    if (j < 1) j = 1;
    if (j > npts) j = npts;
    return tmp[j];
}

// determ(array, norder): array(10,10) column-major, modified as there.
double determ(double* array, int norder) {
    auto A = [&](int i, int j) -> double& { return array[(j - 1) * 10 + (i - 1)]; };
    double det = 1.0;
    for (int k = 1; k <= norder; ++k) {
        if (A(k, k) == 0) {
            int j;
            for (j = k; j <= norder; ++j)
                if (A(k, j) != 0) break;
            if (j > norder) return 0.0;
            for (int i = k; i <= norder; ++i) {
                const double s8 = A(i, j);
                A(i, j) = A(i, k);
                A(i, k) = s8;
            }
            det = -1. * det;
        }
        det = det * A(k, k);
        if (k < norder) {
            const int k1 = k + 1;
            for (int i = k1; i <= norder; ++i)
                for (int j = k1; j <= norder; ++j) A(i, j) = A(i, j) - A(i, k) * A(k, j) / A(k, k);
        }
    }
    return det;
}

// polyfit(x, y, sigmay, npts, nterms, mode=0, a, chisqr): x, y, a 1-based;
// chisqr not computed (nothing reads it).
void polyfit(const double* x, const double* y, int npts, int nterms, double* a) {
    const int nmax = 2 * nterms - 1;
    double sumx[11] = {}, sumy[11] = {}, array[100];
    for (int i = 1; i <= npts; ++i) {
        const double xi = x[i], yi = y[i];
        const double weight = 1;
        double xterm = weight;
        for (int n = 1; n <= nmax; ++n) {
            sumx[n] = sumx[n] + xterm;
            xterm = xterm * xi;
        }
        double yterm = weight * yi;
        for (int n = 1; n <= nterms; ++n) {
            sumy[n] = sumy[n] + yterm;
            yterm = yterm * xi;
        }
    }
    auto A = [&](int i, int j) -> double& { return array[(j - 1) * 10 + (i - 1)]; };
    for (int j = 1; j <= nterms; ++j)
        for (int k = 1; k <= nterms; ++k) A(j, k) = sumx[j + k - 1];
    const double delta = determ(array, nterms);
    if (delta == 0) {
        for (int l = 1; l <= nterms; ++l) a[l] = 0.;
        return;
    }
    for (int l = 1; l <= nterms; ++l) {
        for (int j = 1; j <= nterms; ++j) {
            for (int k = 1; k <= nterms; ++k) A(j, k) = sumx[j + k - 1];
            A(j, l) = sumy[j];
        }
        a[l] = determ(array, nterms) / delta;
    }
}

// fst4_baseline(s, np, ia, ib, npct, sbase): s, sbase 1-based.
void fst4_baseline(const float* s, int np, int ia, int ib, int npct, float* sbase) {
    std::vector<float> sw(static_cast<size_t>(np) + 1, 0.0f);
    // sw(i) = 10.0*log10(s(i)), i = ia..ib: libmvec four at a time, the
    // rest scalar (fewer than 4: all scalar).
    {
        const int count = ib - ia + 1;
        int i = ia;
        if (count > 3) {
            for (; i + 3 <= ib; i += 4) {
                float l[4];
                simd::log10f4(&s[i], l);
                for (int q = 0; q < 4; ++q) sw[i + q] = l[q] * 10.0f;
            }
        }
        for (; i <= ib; ++i) sw[i] = std::log10(s[i]) * 10.0f;
    }

    const int nseg = 8, nterms = 3;
    const int nlen = (ib - ia + 1) / nseg;  // Length of test segment
    const int i0 = (ib - ia + 1) / 2;      // Midpoint
    static thread_local double x[1001], y[1001];
    int k = 0;
    for (int n = 1; n <= nseg; ++n) {
        const int ja = ia + (n - 1) * nlen;
        const int jb = ja + nlen - 1;
        const float base = pctile(&sw[ja] - 1, nlen, npct);  // lowest npct of points
        for (int i = ja; i <= jb; ++i) {
            if (sw[i] <= base) {
                if (k < 1000) ++k;  // Save all "lower envelope" points
                x[k] = static_cast<double>(i - i0);
                y[k] = static_cast<double>(sw[i]);
            }
        }
    }
    const int kz = k;
    double a[6] = {};
    polyfit(x, y, kz, nterms, a);  // Fit a low-order polynomial
    for (int i = 1; i <= np; ++i) sbase[i] = 0.0f;
    for (int i = ia; i <= ib; ++i) {
        const double t = static_cast<double>(i - i0);
        sbase[i] = static_cast<float>(a[1] + t * (a[2] + t * (a[3])) + static_cast<double>(0.2f));
    }
    // sbase = 10**(sbase/10.0) over 1..np: libmvec four at a time from 1.
    {
        int i = 1;
        if (np > 3) {
            for (; i + 3 <= np; i += 4) {
                float e[4];
                for (int q = 0; q < 4; ++q) e[q] = sbase[i + q] / 10.0f;
                simd::powf4(10.0f, e, &sbase[i]);
            }
        }
        for (; i <= np; ++i) sbase[i] = std::pow(10.0f, sbase[i] / 10.0f);
    }
}

}  // namespace

Period Period::of(int ntrperiod) {
    Period p;
    p.ntrperiod = ntrperiod;
    switch (ntrperiod) {
    case 15:
        p.nsps = 720;
        p.nmax = 15 * 12000;
        p.ndown = 18;
        p.nfft1 = (p.nmax / p.ndown) * p.ndown;
        break;
    case 30:
        p.nsps = 1680;
        p.nmax = 30 * 12000;
        p.ndown = 42;
        p.nfft1 = 359856;
        break;
    case 60:
        p.nsps = 3888;
        p.nmax = 60 * 12000;
        p.ndown = 108;
        p.nfft1 = 7500 * 96;
        break;
    case 120:
        p.nsps = 8200;
        p.nmax = 120 * 12000;
        p.ndown = 205;
        p.nfft1 = 7200 * 200;
        break;
    case 300:
        p.nsps = 21504;
        p.nmax = 300 * 12000;
        p.ndown = 512;
        p.nfft1 = 7020 * 512;
        break;
    case 900:
        p.nsps = 66560;
        p.nmax = 900 * 12000;
        p.ndown = 1664;
        p.nfft1 = 6480 * 1664;
        break;
    case 1800:
        p.nsps = 134400;
        p.nmax = 1800 * 12000;
        p.ndown = 3360;
        p.nfft1 = 6426 * 3360;
        break;
    default:
        throw std::invalid_argument("fst4w: FST4W periods are 15, 30, 60, 120, 300, 900 and 1800 s");
    }
    p.nss = p.nsps / p.ndown;
    p.fs = 12000.0f;
    p.fs2 = p.fs / static_cast<float>(p.ndown);
    p.nspsec = static_cast<int>(std::lround(p.fs2));
    p.dt = 1.0f / p.fs;
    p.dt2 = 1.0f / p.fs2;
    p.tt = static_cast<float>(p.nsps) * p.dt;
    p.baud = 1.0f / p.tt;
    p.sigbw = 4.0f * p.baud;
    p.nfft2 = p.nfft1 / p.ndown;
    p.nfft1 = p.nfft2 * p.ndown;
    return p;
}

void blanker(int16_t* iwave, int nz, int ndropmax, int npct, float* c_bigfft) {
    static thread_local std::vector<int> hist;
    hist.assign(32769, 0);
    const float fblank = 0.01f * static_cast<float>(npct);
    for (int i = 0; i < nz; ++i) {
        if (iwave[i] == -32768) iwave[i] = -32767;
        const int n = std::abs(static_cast<int>(iwave[i]));
        hist[n] = hist[n] + 1;
    }
    int n = 0, i;
    const int target =
        static_cast<int>(std::lround(static_cast<float>(nz) * fblank / static_cast<float>(ndropmax)));
    for (i = 32768; i >= 0; --i) {
        n = n + hist[i];
        if (n >= target) break;
    }
    const int nthresh = i;
    int ndrop = 0;
    float xx = 0.0f;
    for (int ii = 1; ii <= nz; ++ii) {
        int i0 = iwave[ii - 1];
        if (ndrop > 0) {
            i0 = 0;
            ndrop = ndrop - 1;
        }
        // Start to apply blanking
        if (std::abs(i0) > nthresh) {
            i0 = 0;
            ndrop = ndropmax;
        }
        // Now copy the data into c_bigfft
        if ((ii & 1) == 1) {
            xx = static_cast<float>(i0);
        } else {
            const int j = ii / 2 - 1;
            c_bigfft[2 * j] = xx;
            c_bigfft[2 * j + 1] = static_cast<float>(i0);
        }
    }
}

void blanker(float* iwave, int nz, int ndropmax, int npct, float* c_bigfft) {
    static thread_local std::vector<int> hist;
    hist.assign(32769, 0);
    const float fblank = 0.01f * static_cast<float>(npct);
    for (int i = 0; i < nz; ++i) {
        const float a = std::fabs(iwave[i]);
        const int n = a >= 32768.0f ? 32768 : static_cast<int>(a);
        hist[n] = hist[n] + 1;
    }
    int n = 0, i;
    const int target =
        static_cast<int>(std::lround(static_cast<float>(nz) * fblank / static_cast<float>(ndropmax)));
    for (i = 32768; i >= 0; --i) {
        n = n + hist[i];
        if (n >= target) break;
    }
    // An int sample's |i0| > nthresh. The top bin holds |x| >= 32768 here,
    // which no int sample reaches: nthresh 32768 (npct 0) blanks nothing.
    const float nthresh = i >= 32768 ? std::numeric_limits<float>::infinity() : static_cast<float>(i);
    int ndrop = 0;
    float xx = 0.0f;
    for (int ii = 1; ii <= nz; ++ii) {
        float i0 = iwave[ii - 1];
        if (ndrop > 0) {
            i0 = 0.0f;
            ndrop = ndrop - 1;
        }
        // Start to apply blanking (an int sample's |i0| > nthresh is |i0| >= nthresh+1)
        if (std::fabs(i0) >= nthresh + 1.0f) {
            i0 = 0.0f;
            ndrop = ndropmax;
        }
        if ((ii & 1) == 1) {
            xx = i0;
        } else {
            const int j = ii / 2 - 1;
            c_bigfft[2 * j] = xx;
            c_bigfft[2 * j + 1] = i0;
        }
    }
}

void* fftw_alloc_bytes(size_t n) { return fftwf_malloc(n); }
void fftw_free_bytes(void* p) { fftwf_free(p); }

void fft_r2c_inplace(float* data, int n) {
    const CachedPlan c = cached_plan(n, 0);
    if (fftwf_alignment_of(data) == c.alignment) {
        fftwf_execute_dft_r2c(c.plan, data, reinterpret_cast<fftwf_complex*>(data));
        return;
    }
    const size_t keep = n <= NSMALL ? sizeof(float) * 2 * static_cast<size_t>(n / 2 + 1) : 0;
    plan_and_run(data, keep, [&] {
        return fftwf_plan_dft_r2c_1d(n, data, reinterpret_cast<fftwf_complex*>(data), PATIENCE);
    });
}

void fft_c2c_inplace(cfloat* data, int n, int isign) {
    auto* d = reinterpret_cast<fftwf_complex*>(data);
    const CachedPlan c = cached_plan(n, isign == -1 ? -1 : 1);
    if (fftwf_alignment_of(reinterpret_cast<float*>(data)) == c.alignment) {
        fftwf_execute_dft(c.plan, d, d);
        return;
    }
    const size_t keep = n <= NSMALL ? sizeof(cfloat) * static_cast<size_t>(n) : 0;
    plan_and_run(data, keep, [&] {
        return fftwf_plan_dft_1d(n, d, d, isign == -1 ? FFTW_FORWARD : FFTW_BACKWARD, PATIENCE);
    });
}

bool get_candidates_fst4(const float* c_bigfft, int nfft1, int nsps, int hmod, float fs, float fa, float fb, int nfa,
                         int nfb, float minsync, int& ncand, float* candidates) {
    static const float xdb[7] = {0.25f, 0.50f, 0.75f, 1.0f, 0.75f, 0.50f, 0.25f};  // xdb(-3:3)
    auto re = [&](int j) { return c_bigfft[2 * j]; };
    auto im = [&](int j) { return c_bigfft[2 * j + 1]; };

    const float df1 = fs / static_cast<float>(nfft1);
    const float baud = fs / static_cast<float>(nsps);  // Keying rate
    const float df2 = baud / 2.0f;
    const int nd = static_cast<int>(df2 / df1);  // s() sums this many bins of big FFT
    const int ndh = nd / 2;
    int ia = static_cast<int>(std::lround(std::max(100.0f, fa) / df2));  // Low frequency search limit
    int ib = static_cast<int>(std::lround(std::min(4800.0f, fb) / df2));  // High frequency limit
    int ina = static_cast<int>(std::lround(std::max(100.0f, static_cast<float>(nfa)) / df2));
    int inb = static_cast<int>(std::lround(std::min(4800.0f, static_cast<float>(nfb)) / df2));
    if (ia < ina) ia = ina;
    if (ib > inb) ib = inb;

    const int nnw = static_cast<int>(std::lround(48000.0f * static_cast<float>(nsps) * 2.0f / fs));
    std::vector<float> s(static_cast<size_t>(nnw) + 1, 0.0f);  // Low resolution power spectrum, 1-based
    for (int i = ina; i <= inb; ++i) {
        const int j0 = static_cast<int>(std::lround(static_cast<float>(i) * df2 / df1));
        for (int j = j0 - ndh; j <= j0 + ndh; ++j) s[i] = s[i] + re(j) * re(j) + im(j) * im(j);
    }

    ina = std::max(ina, 1 + 3 * hmod);  // Don't run off the ends
    inb = std::min(inb, nnw - 3 * hmod);
    std::vector<float> s2(static_cast<size_t>(nnw) + 1, 0.0f), sbase(static_cast<size_t>(nnw) + 1, 0.0f);
    for (int i = ina; i <= inb; ++i) s2[i] = s[i - hmod * 3] + s[i - hmod] + s[i + hmod] + s[i + hmod * 3];
    const int npctile = 30;
    fst4_baseline(s2.data(), nnw, ina + hmod * 3, inb - hmod * 3, npctile, sbase.data());
    for (int i = ina; i <= inb; ++i)
        if (sbase[i] <= 0.0f) return false;
    for (int i = ina; i <= inb; ++i) s2[i] = s2[i] / sbase[i];  // Normalize wrt noise level

    ncand = 0;
    std::fill(candidates, candidates + 200 * 5, 0.0f);
    auto C = [&](int i, int j) -> float& { return candidates[(j - 1) * 200 + (i - 1)]; };
    if (ia < 3) ia = 3;
    if (ib > nnw - 2) ib = nnw - 2;

    // Find candidates, using the CLEAN algorithm to remove a model of each
    // one from s2() after it has been found.
    while (ncand < 200) {
        const int iploc = maxloc(ia, ib, [&](int i) { return s2[i]; });  // Index of CCF peak
        const float pval = s2[iploc];                                    // Peak value
        if (pval < minsync) break;
        for (int i = -3; i <= 3; ++i) {  // Remove 0.9 of a model CCF at this frequency from s2()
            const int k = iploc + 2 * hmod * i;
            if (k >= ia && k <= ib) {
                // max(0., v) as gfortran's MAX: a NaN argument is passed over
                // (on silence the baseline is NaN, and this clears the peak).
                const float v = s2[k] - 0.9f * pval * xdb[i + 3];
                float m = 0.0f;
                if (v > m || std::isnan(m)) m = v;
                s2[k] = m;
            }
        }
        ncand = ncand + 1;
        C(ncand, 1) = df2 * static_cast<float>(iploc);  // Candidate frequency
        C(ncand, 2) = pval;                              // Rough estimate of SNR
        C(ncand, 5) = sbase[iploc];
    }
    return true;
}

void fst4_downsample(const float* c_bigfft, int nfft1, int ndown, float f0, float sigbw, cfloat* c1) {
    auto cb = [&](int j) { return cfloat(c_bigfft[2 * j], c_bigfft[2 * j + 1]); };
    const float df = 12000.0f / static_cast<float>(nfft1);
    const int i0 = static_cast<int>(std::lround(f0 / df));
    const int ih = static_cast<int>(std::lround((f0 + 1.3f * sigbw / 2.0f) / df));
    const int nbw = ih - i0 + 1;
    const int nfft2 = nfft1 / ndown;
    for (int i = 0; i < nfft2; ++i) c1[i] = 0.0f;
    c1[0] = cb(i0);
    for (int i = 1; i <= nbw; ++i) {
        if (i0 + i <= nfft1 / 2) c1[i] = cb(i0 + i);
        if (i0 - i >= 0) c1[nfft2 - i] = cb(i0 - i);
    }
    const float r = static_cast<float>(nfft2);
    for (int i = 0; i < nfft2; ++i) c1[i] = cfloat(c1[i].real() / r, c1[i].imag() / r);
    fft_c2c_inplace(c1, nfft2, 1);  // c2c FFT back to time domain
}

void Frontend::sync_tables(float f0, int hmod, int nss, int ntr, float fs) {
    static const int isyncword1[8] = {0, 1, 3, 2, 1, 0, 2, 3};
    static const int isyncword2[8] = {2, 3, 1, 0, 3, 2, 0, 1};
    const int nz = 8 * nss;
    if (nss != nss0_ || ntr != ntr0_) {
        twopi_ = 8.0f * std::atan(1.0f);
        sdt_ = 1 / fs;
        int k = 1;
        float phi1 = 0.0f, phi2 = 0.0f;
        for (int i = 0; i <= 7; ++i) {
            const float dphi1 = twopi_ * static_cast<float>(hmod) * (static_cast<float>(isyncword1[i]) - 1.5f) /
                                static_cast<float>(nss);
            const float dphi2 = twopi_ * static_cast<float>(hmod) * (static_cast<float>(isyncword2[i]) - 1.5f) /
                                static_cast<float>(nss);
            for (int j = 1; j <= nss; ++j) {
                float s, c;
                sincosf(phi1, &s, &c);
                csync1_[k - 1] = cfloat(c, s);
                sincosf(phi2, &s, &c);
                csync2_[k - 1] = cfloat(c, s);
                phi1 = std::fmod(phi1 + dphi1, twopi_);
                phi2 = std::fmod(phi2 + dphi2, twopi_);
                ++k;
            }
        }
        fac_ = 1.0f / (8.0f * static_cast<float>(nss));
        nss0_ = nss;
        ntr0_ = ntr;
        f0save_ = -1.e30f;
    }

    if (f0 != f0save_) {
        const float dphi = twopi_ * f0 * sdt_;
        float phi = 0.0f;
        for (int i = 1; i <= nz; ++i) {
            float s, c;
            sincosf(phi, &s, &c);
            const cfloat ctwk(c, s);
            csynct1_[i - 1] = ctwk * csync1_[i - 1];
            csynct2_[i - 1] = ctwk * csync2_[i - 1];
            phi = std::fmod(phi + dphi, twopi_);
        }
        f0save_ = f0;
        nss0_ = nss;
    }
}

float Frontend::sync_fst4(const cfloat* cd0, int i0, float f0, int hmod, int ncoh, int np, int nss, int ntr, float fs) {
    const int nz = 8 * nss;
    sync_tables(f0, hmod, nss, ntr, fs);

    const int i1 = i0, i2 = i0 + 38 * nss, i3 = i0 + 76 * nss, i4 = i0 + 114 * nss, i5 = i0 + 152 * nss;
    float s1 = 0.0f, s2 = 0.0f, s3 = 0.0f, s4 = 0.0f, s5 = 0.0f;
    // sum(cd0(a:a+n-1)*conjg(cs(b+1:b+n))), in order.
    auto corr = [&](int a, const std::vector<cfloat>& cs, int b, int n) {
        float zr = 0.0f, zi = 0.0f;
        for (int t = 0; t < n; ++t) {
            const cfloat x = cd0[a + t];
            const cfloat y = cs[b + t];
            // x*conjg(y), as gfortran multiplies: (xr*yr - xi*(-yi)) + i(xr*(-yi) + xi*yr)
            const float yi = -y.imag();
            zr = zr + (x.real() * y.real() - x.imag() * yi);
            zi = zi + (x.real() * yi + x.imag() * y.real());
        }
        return std::abs(cfloat(zr, zi));
    };
    if (ncoh > 0) {
        const int nsec = 8 / ncoh;
        for (int i = 1; i <= nsec; ++i) {
            const int is = (i - 1) * ncoh * nss;
            const int n = ncoh * nss;
            const float az1 = (i1 + is >= 1) ? corr(i1 + is, csynct1_, is, n) : 0.0f;
            const float az2 = corr(i2 + is, csynct2_, is, n);
            const float az3 = corr(i3 + is, csynct1_, is, n);
            const float az4 = corr(i4 + is, csynct2_, is, n);
            const float az5 = (i5 + is + ncoh * nss - 1 <= np) ? corr(i5 + is, csynct1_, is, n) : 0.0f;
            const float fnz = static_cast<float>(nz);
            s1 = s1 + az1 / fnz;
            s2 = s2 + az2 / fnz;
            s3 = s3 + az3 / fnz;
            s4 = s4 + az4 / fnz;
            s5 = s5 + az5 / fnz;
        }
    } else {
        const int nsub = -ncoh;
        const int nps = nss / nsub;
        const float d = static_cast<float>(8 * nss);
        for (int i = 1; i <= 8; ++i) {
            for (int isub = 1; isub <= nsub; ++isub) {
                const int is = (i - 1) * nss + (isub - 1) * nps;
                const float az1 = (i1 + is >= 1) ? corr(i1 + is, csynct1_, is, nps) : 0.0f;
                const float az2 = corr(i2 + is, csynct2_, is, nps);
                const float az3 = corr(i3 + is, csynct1_, is, nps);
                const float az4 = corr(i4 + is, csynct2_, is, nps);
                const float az5 = (i5 + is + ncoh * nss - 1 <= np) ? corr(i5 + is, csynct1_, is, nps) : 0.0f;
                s1 = s1 + az1 / d;
                s2 = s2 + az2 / d;
                s3 = s3 + az3 / d;
                s4 = s4 + az4 / d;
                s5 = s5 + az5 / d;
            }
        }
    }
    return s1 + s2 + s3 + s4 + s5;
}

void Frontend::fst4_sync_search(const cfloat* c2, int nfft2, int hmod, float fs2, int nss, int ntrperiod, int nsyncoh,
                                float emedelay, float& sbest, float& fcbest, int& isbest) {
    const int nspsec = static_cast<int>(fs2);
    const float baud = fs2 / static_cast<float>(nss);
    float fc1 = 0.0f;
    int is0, ishw;
    if (emedelay < 0.1f) {  // search offsets from 0 s to 2 s
        is0 = static_cast<int>(1.5f * static_cast<float>(nspsec));
        ishw = static_cast<int>(1.5f * static_cast<float>(nspsec));
    } else {  // search plus or minus 1.5 s centered on emedelay
        is0 = static_cast<int>(std::lround((emedelay + 1.0f) * static_cast<float>(nspsec)));
        ishw = static_cast<int>(1.5f * static_cast<float>(nspsec));
    }

    const bool lanes = nsyncoh == 8;  // FST4W's: one coherent section
    std::vector<float> syncs;
    auto run = [&](int first, int last, int step, float fc) {
        const int n = first <= last ? (last - first) / step + 1 : 0;
        syncs.assign(static_cast<size_t>(n), 0.0f);
        if (lanes) {
            sync_lanes(c2, first, step, n, fc, hmod, nfft2, nss, ntrperiod, fs2, syncs.data());
        } else {
            for (int q = 0; q < n; ++q) syncs[q] = sync_fst4(c2, first + q * step, fc, hmod, nsyncoh, nfft2, nss, ntrperiod, fs2);
        }
        return n;
    };

    if (lanes) split(c2, nfft2, 4 * hmod);
    sbest = -1.e30f;
    for (int ifr = -12; ifr <= 12; ++ifr) {
        const float fc = fc1 + 0.1f * baud * static_cast<float>(ifr);
        const int first = std::max(1, is0 - ishw), step = 4 * hmod;
        const int n = run(first, is0 + ishw, step, fc);
        for (int q = 0; q < n; ++q) {
            if (syncs[q] > sbest) {
                fcbest = fc;
                isbest = first + q * step;
                sbest = syncs[q];
            }
        }
    }

    fc1 = fcbest;
    is0 = isbest;
    ishw = 4 * hmod;
    const int isst = 1 * hmod;

    if (lanes) split(c2, nfft2, isst);
    sbest = 0.0f;
    for (int ifr = -7; ifr <= 7; ++ifr) {
        const float fc = fc1 + 0.02f * baud * static_cast<float>(ifr);
        const int first = std::max(1, is0 - ishw);
        const int n = run(first, is0 + ishw, isst, fc);
        for (int q = 0; q < n; ++q) {
            if (syncs[q] > sbest) {
                fcbest = fc;
                isbest = first + q * isst;
                sbest = syncs[q];
            }
        }
    }
}

void Frontend::split(const cfloat* cd0, int np, int step) {
    // cd0 split by residue mod step: D_r[q] = cd0[q*step + r], zero padded,
    // so the eight lanes' samples (starts step apart) lie side by side.
    // Padding: the fifth section of lanes past the data reads zeros (those
    // lanes' z5 is not used).
    const int len = (np + 8 * 160 * 40 + 8 * step) / step + 2;
    split_len_ = len;
    dre_.assign(static_cast<size_t>(step) * len, 0.0f);
    dim_.assign(static_cast<size_t>(step) * len, 0.0f);
    for (int i = 0; i < np; ++i) {
        const size_t at = static_cast<size_t>(i % step) * len + i / step;
        dre_[at] = cd0[i].real();
        dim_[at] = cd0[i].imag();
    }
}

void Frontend::sync_lanes(const cfloat* cd0, int first, int step, int n, float f0, int hmod, int np, int nss, int ntr,
                          float fs, float* out) {
    sync_tables(f0, hmod, nss, ntr, fs);
    const int nz = 8 * nss;
    (void)cd0;
    const int len = split_len_;
    const int off[5] = {0, 38 * nss, 76 * nss, 114 * nss, 152 * nss};
    const float fnz = static_cast<float>(nz);
    const simd::V8 vnz = simd::set1(fnz);
    const cfloat* y1 = csynct1_.data();
    const cfloat* y2 = csynct2_.data();
    for (int q0 = 0; q0 < n; q0 += 8) {
        // The five sections' sums side by side (ten independent chains),
        // each lane's in sync_fst4's order: t = 0, 1, ... nz-1.
        simd::V8 zr[5], zi[5];
        int r[5], q[5];
        for (int k = 0; k < 5; ++k) {
            zr[k] = zi[k] = simd::zero();
            const int idx = first + q0 * step + off[k];  // lane 0's first sample
            r[k] = idx % step;
            q[k] = idx / step;
        }
        const float* re = dre_.data();
        const float* im = dim_.data();
        for (int t = 0; t < nz; ++t) {
            const simd::V8 y1r = simd::set1(y1[t].real()), y1i = simd::set1(-y1[t].imag());
            const simd::V8 y2r = simd::set1(y2[t].real()), y2i = simd::set1(-y2[t].imag());
            for (int k = 0; k < 5; ++k) {
                const size_t at = static_cast<size_t>(r[k]) * len + q[k];
                const simd::V8 xr = simd::load(re + at), xi = simd::load(im + at);
                const simd::V8 yr = (k & 1) ? y2r : y1r, yi = (k & 1) ? y2i : y1i;
                // x*conjg(y) as sync_fst4: (xr*yr - xi*(-yi)) + i(xr*(-yi) + xi*yr)
                zr[k] = simd::add(zr[k], simd::sub(simd::mul(xr, yr), simd::mul(xi, yi)));
                zi[k] = simd::add(zi[k], simd::add(simd::mul(xr, yi), simd::mul(xi, yr)));
                if (++r[k] == step) {
                    r[k] = 0;
                    ++q[k];
                }
            }
        }
        float a[5][8];
        for (int k = 0; k < 5; ++k) simd::store(a[k], simd::div(simd::abs_exact(zr[k], zi[k]), vnz));  // 0 + abs(z_k)/nz
        for (int l = 0; l < 8 && q0 + l < n; ++l) {
            const int i0 = first + (q0 + l) * step;
            // z1 needs i1 >= 1 (always: starts are at least 1); z5 only when
            // its samples lie within the data.
            const float s1 = a[0][l], s2 = a[1][l], s3 = a[2][l], s4 = a[3][l];
            const float s5 = (i0 + 152 * nss + nz - 1 <= np) ? a[4][l] : 0.0f;
            out[q0 + l] = s1 + s2 + s3 + s4 + s5;
        }
    }
}

namespace {

// ft8b.f90's normalizebmet(bmet, n).
void normalizebmet(float* bmet, int n) {
    float sum1 = 0.0f, sum2 = 0.0f;
    for (int i = 0; i < n; ++i) sum1 = sum1 + bmet[i];
    for (int i = 0; i < n; ++i) sum2 = sum2 + bmet[i] * bmet[i];
    const float bmetav = sum1 / static_cast<float>(n);
    const float bmet2av = sum2 / static_cast<float>(n);
    const float var = bmet2av - bmetav * bmetav;
    float bmetsig;
    if (var > 0.0f) bmetsig = std::sqrt(var);
    else bmetsig = std::sqrt(bmet2av);
    for (int i = 0; i < n; ++i) bmet[i] = bmet[i] / bmetsig;
}

// maxval(s2(0:n-1), mask): gfortran's, the masked elements' largest
// (-huge when none; NaNs skipped unless all are NaN).
template <typename M>
float maxval_masked(const float* s, int n, M mask) {
    float m = -std::numeric_limits<float>::max();
    bool any = false, nonnan = false;
    for (int i = 0; i < n; ++i) {
        if (!mask(i)) continue;
        any = true;
        const float v = s[i];
        if (!nonnan) {
            if (!std::isnan(v)) {
                nonnan = true;
                m = v;
            }
        } else if (v > m) {
            m = v;
        }
    }
    if (any && !nonnan) return std::numeric_limits<float>::quiet_NaN();
    return m;
}

}  // namespace

void Frontend::get_fst4_bitmetrics_reference(const cfloat* cd, int nss, float* bitmetrics, float* s4, int& nsync_qual,
                                             bool& badsync) {
    static const int isyncword1[8] = {0, 1, 3, 2, 1, 0, 2, 3};
    static const int isyncword2[8] = {2, 3, 1, 0, 3, 2, 0, 1};
    static const int graymap[4] = {0, 1, 3, 2};
    auto BM = [&](int i, int j) -> float& { return bitmetrics[(j - 1) * 2 * NN + (i - 1)]; };
    auto S4 = [&](int itone, int k) -> float& { return s4[(k - 1) * 4 + itone]; };
    auto one = [](int i, int j) { return (i & (1 << j)) != 0; };

    if (bm_first_ || nss != bm_nss0_) {
        ci_.assign(static_cast<size_t>(nss) * 4, cfloat(0.0f, 0.0f));
        const float twopi = 8.0f * std::atan(1.0f);
        const float dphi = twopi / static_cast<float>(nss);
        for (int itone = 0; itone <= 3; ++itone) {
            const float dp = (static_cast<float>(itone) - 1.5f) * dphi;
            float phi = 0.0f;
            for (int j = 1; j <= nss; ++j) {
                float s, c;
                sincosf(phi, &s, &c);
                ci_[static_cast<size_t>(itone) * nss + j - 1] = cfloat(c, s);
                phi = std::fmod(phi + dp, twopi);
            }
        }
        bm_first_ = false;
        bm_nss0_ = nss;
    }

    // cs(0:3, NN)
    cfloat cs[4][NN + 1];
    for (int k = 1; k <= NN; ++k) {
        const int i1 = (k - 1) * nss;
        for (int itone = 0; itone <= 3; ++itone) {
            float zr = 0.0f, zi = 0.0f;
            const cfloat* c = &ci_[static_cast<size_t>(itone) * nss];
            for (int j = 0; j < nss; ++j) {
                const cfloat x = cd[i1 + j];
                const float yr = c[j].real(), yi = -c[j].imag();
                zr = zr + (x.real() * yr - x.imag() * yi);
                zi = zi + (x.real() * yi + x.imag() * yr);
            }
            cs[itone][k] = cfloat(zr, zi);
        }
        for (int itone = 0; itone <= 3; ++itone) {
            const float a = std::abs(cs[itone][k]);
            S4(itone, k) = a * a;
        }
    }

    // Sync quality check
    int is1 = 0, is2 = 0, is3 = 0, is4 = 0, is5 = 0;
    badsync = false;
    auto ploc = [&](int k) { return maxloc(0, 3, [&](int t) { return S4(t, k); }); };
    for (int k = 1; k <= 8; ++k) {
        if (isyncword1[k - 1] == ploc(k)) ++is1;
        if (isyncword2[k - 1] == ploc(k + 38)) ++is2;
        if (isyncword1[k - 1] == ploc(k + 76)) ++is3;
        if (isyncword2[k - 1] == ploc(k + 114)) ++is4;
        if (isyncword1[k - 1] == ploc(k + 152)) ++is5;
    }
    const int nsync = is1 + is2 + is3 + is4 + is5;  // correct hard sync symbols, 0-40
    badsync = false;
    if (nsync < 16) {
        badsync = true;
        return;
    }

    for (int i = 0; i < 2 * NN * 4; ++i) bitmetrics[i] = 0.0f;

    cfloat c1[5][9], c2[17][5], c4[257][3];  // c1(4,8), c2(16,4), c4(256,2), 1-based
    float* s2 = s2_.data();
    for (int k = 1; k <= NN; k += 8) {
        for (int m = 1; m <= 8; ++m) {  // 4 1-symbol correlations for each of 8 symbols
            for (int n = 1; n <= 4; ++n) {
                c1[n][m] = cs[graymap[n - 1]][k + m - 1];
                s2[n - 1] = std::abs(c1[n][m]);
            }
            const int ipt = (k - 1) * 2 + 2 * (m - 1) + 1;
            for (int ib = 0; ib <= 1; ++ib) {
                const float bm = maxval_masked(s2, 4, [&](int i) { return one(i, 1 - ib); }) -
                                 maxval_masked(s2, 4, [&](int i) { return !one(i, 1 - ib); });
                if (ipt + ib > 2 * NN) continue;
                BM(ipt + ib, 1) = bm;
            }
        }

        for (int m = 1; m <= 4; ++m) {  // 16 2-symbol correlations for each of 4 2-symbol groups
            for (int i = 1; i <= 4; ++i)
                for (int j = 1; j <= 4; ++j) {
                    const int is = (i - 1) * 4 + j;
                    c2[is][m] = c1[i][2 * m - 1] - c1[j][2 * m];
                    s2[is - 1] = std::abs(c2[is][m]);
                }
            const int ipt = (k - 1) * 2 + 4 * (m - 1) + 1;
            for (int ib = 0; ib <= 3; ++ib) {
                const float bm = maxval_masked(s2, 16, [&](int i) { return one(i, 3 - ib); }) -
                                 maxval_masked(s2, 16, [&](int i) { return !one(i, 3 - ib); });
                if (ipt + ib > 2 * NN) continue;
                BM(ipt + ib, 2) = bm;
            }
        }

        for (int m = 1; m <= 2; ++m) {  // 256 4-symbol correlations for each of 2 4-symbol groups
            for (int i = 1; i <= 16; ++i)
                for (int j = 1; j <= 16; ++j) {
                    const int is = (i - 1) * 16 + j;
                    c4[is][m] = c2[i][2 * m - 1] + c2[j][2 * m];
                    s2[is - 1] = std::abs(c4[is][m]);
                }
            const int ipt = (k - 1) * 2 + 8 * (m - 1) + 1;
            for (int ib = 0; ib <= 7; ++ib) {
                const float bm = maxval_masked(s2, 256, [&](int i) { return one(i, 7 - ib); }) -
                                 maxval_masked(s2, 256, [&](int i) { return !one(i, 7 - ib); });
                if (ipt + ib > 2 * NN) continue;
                BM(ipt + ib, 3) = bm;
            }
        }

        // 65536 8-symbol correlations for the entire group (the reference's
        // s2=0 before each level is not copied: only what was just written is read)
        for (int i = 1; i <= 256; ++i)
            for (int j = 1; j <= 256; ++j) {
                const int is = (i - 1) * 256 + j;
                s2[is - 1] = std::abs(c4[i][1] + c4[j][2]);
            }
        const int ipt = (k - 1) * 2 + 1;
        for (int ib = 0; ib <= 15; ++ib) {
            const float bm = maxval_masked(s2, 65536, [&](int i) { return one(i, 15 - ib); }) -
                             maxval_masked(s2, 65536, [&](int i) { return !one(i, 15 - ib); });
            if (ipt + ib > 2 * NN) continue;
            BM(ipt + ib, 4) = bm;
        }
    }

    static const int sw1[16] = {0, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0};
    static const int sw2[16] = {1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1};
    auto hbit = [&](int i) { return BM(i, 1) >= 0 ? 1 : 0; };
    auto count = [&](int start, const int* w) {
        int c = 0;
        for (int t = 0; t < 16; ++t)
            if (hbit(start + t) == w[t]) ++c;
        return c;
    };
    nsync_qual = count(1, sw1) + count(77, sw2) + count(153, sw1) + count(229, sw2) + count(305, sw1);
    if (nsync_qual < 46) {
        badsync = true;
        return;
    }

    for (int j = 1; j <= 4; ++j) normalizebmet(&BM(1, j), 2 * NN);
    const float scalefac = 2.83f;
    for (int i = 0; i < 2 * NN * 4; ++i) bitmetrics[i] = scalefac * bitmetrics[i];
}

namespace {

// For values v[(i-1)*n + (j-1)], i, j = 1..n (n a power of two, 2^h): the
// largest of those whose index has bit b set and of those with it clear,
// from the rows' and columns' largest (maxima of exact values: any order).
struct HalfMax {
    float set, clear;
};
inline HalfMax half_max(const float* rowmax, const float* colmax, int n, int h, int b) {
    const float* m = b < h ? colmax : rowmax;
    const int bit = b < h ? b : b - h;
    float s = -std::numeric_limits<float>::max(), c = -std::numeric_limits<float>::max();
    for (int t = 0; t < n; ++t) {
        if ((t >> bit) & 1) s = std::max(s, m[t]);
        else c = std::max(c, m[t]);
    }
    return {s, c};
}

}  // namespace

void Frontend::get_fst4_bitmetrics(const cfloat* cd, int nss, float* bitmetrics, float* s4, int& nsync_qual,
                                   bool& badsync) {
    for (int i = 0; i < NN * nss; ++i)
        if (!std::isfinite(cd[i].real()) || !std::isfinite(cd[i].imag())) {
            get_fst4_bitmetrics_reference(cd, nss, bitmetrics, s4, nsync_qual, badsync);
            return;
        }
    static const int isyncword1[8] = {0, 1, 3, 2, 1, 0, 2, 3};
    static const int isyncword2[8] = {2, 3, 1, 0, 3, 2, 0, 1};
    static const int graymap[4] = {0, 1, 3, 2};
    auto BM = [&](int i, int j) -> float& { return bitmetrics[(j - 1) * 2 * NN + (i - 1)]; };
    auto S4 = [&](int itone, int k) -> float& { return s4[(k - 1) * 4 + itone]; };

    if (bm_first_ || nss != bm_nss0_) {
        ci_.assign(static_cast<size_t>(nss) * 4, cfloat(0.0f, 0.0f));
        const float twopi = 8.0f * std::atan(1.0f);
        const float dphi = twopi / static_cast<float>(nss);
        for (int itone = 0; itone <= 3; ++itone) {
            const float dp = (static_cast<float>(itone) - 1.5f) * dphi;
            float phi = 0.0f;
            for (int j = 1; j <= nss; ++j) {
                float sn, c;
                sincosf(phi, &sn, &c);
                ci_[static_cast<size_t>(itone) * nss + j - 1] = cfloat(c, sn);
                phi = std::fmod(phi + dp, twopi);
            }
        }
        bm_first_ = false;
        bm_nss0_ = nss;
    }

    cfloat cs[4][NN + 1];
    for (int k = 1; k <= NN; ++k) {
        const int i1 = (k - 1) * nss;
        for (int itone = 0; itone <= 3; ++itone) {
            float zr = 0.0f, zi = 0.0f;
            const cfloat* c = &ci_[static_cast<size_t>(itone) * nss];
            for (int j = 0; j < nss; ++j) {
                const cfloat x = cd[i1 + j];
                const float yr = c[j].real(), yi = -c[j].imag();
                zr = zr + (x.real() * yr - x.imag() * yi);
                zi = zi + (x.real() * yi + x.imag() * yr);
            }
            cs[itone][k] = cfloat(zr, zi);
        }
        for (int itone = 0; itone <= 3; ++itone) {
            const float a = abs_exact(cs[itone][k].real(), cs[itone][k].imag());
            S4(itone, k) = a * a;
        }
    }
    // A correlation overflowing to infinity: the line-for-line code.
    for (int k = 1; k <= NN; ++k)
        for (int t = 0; t < 4; ++t)
            if (!std::isfinite(cs[t][k].real()) || !std::isfinite(cs[t][k].imag()) || !std::isfinite(S4(t, k))) {
                get_fst4_bitmetrics_reference(cd, nss, bitmetrics, s4, nsync_qual, badsync);
                return;
            }

    int is1 = 0, is2 = 0, is3 = 0, is4 = 0, is5 = 0;
    auto ploc = [&](int k) { return maxloc(0, 3, [&](int t) { return S4(t, k); }); };
    for (int k = 1; k <= 8; ++k) {
        if (isyncword1[k - 1] == ploc(k)) ++is1;
        if (isyncword2[k - 1] == ploc(k + 38)) ++is2;
        if (isyncword1[k - 1] == ploc(k + 76)) ++is3;
        if (isyncword2[k - 1] == ploc(k + 114)) ++is4;
        if (isyncword1[k - 1] == ploc(k + 152)) ++is5;
    }
    const int nsync = is1 + is2 + is3 + is4 + is5;
    badsync = false;
    if (nsync < 16) {
        badsync = true;
        return;
    }
    for (int i = 0; i < 2 * NN * 4; ++i) bitmetrics[i] = 0.0f;

    // c1(4,8), c2(16,4), c4(256,2) as split real/imaginary arrays, 1-based.
    float c1r[5][9], c1i[5][9], c2r[17][5], c2i[17][5];
    alignas(32) float c4r[3][257 + 7], c4i[3][257 + 7];
    alignas(32) float rowmax[256], colmax[256];
    for (int k = 1; k <= NN; k += 8) {
        for (int m = 1; m <= 8; ++m) {  // 4 1-symbol correlations for each of 8 symbols
            float s2v[4];
            for (int n = 1; n <= 4; ++n) {
                c1r[n][m] = cs[graymap[n - 1]][k + m - 1].real();
                c1i[n][m] = cs[graymap[n - 1]][k + m - 1].imag();
                s2v[n - 1] = abs_exact(c1r[n][m], c1i[n][m]);
            }
            const int ipt = (k - 1) * 2 + 2 * (m - 1) + 1;
            for (int ib = 0; ib <= 1; ++ib) {
                const int b = 1 - ib;
                float s = -std::numeric_limits<float>::max(), c = s;
                for (int t = 0; t < 4; ++t) {
                    if ((t >> b) & 1) s = std::max(s, s2v[t]);
                    else c = std::max(c, s2v[t]);
                }
                if (ipt + ib > 2 * NN) continue;
                BM(ipt + ib, 1) = s - c;
            }
        }
        for (int m = 1; m <= 4; ++m) {  // 16 2-symbol correlations for each of 4 2-symbol groups
            float rm[4], cm[4];
            for (int t = 0; t < 4; ++t) rm[t] = cm[t] = -std::numeric_limits<float>::max();
            for (int i = 1; i <= 4; ++i)
                for (int j = 1; j <= 4; ++j) {
                    const int is = (i - 1) * 4 + j;
                    c2r[is][m] = c1r[i][2 * m - 1] - c1r[j][2 * m];
                    c2i[is][m] = c1i[i][2 * m - 1] - c1i[j][2 * m];
                    const float a = abs_exact(c2r[is][m], c2i[is][m]);
                    rm[i - 1] = std::max(rm[i - 1], a);
                    cm[j - 1] = std::max(cm[j - 1], a);
                }
            const int ipt = (k - 1) * 2 + 4 * (m - 1) + 1;
            for (int ib = 0; ib <= 3; ++ib) {
                const HalfMax h = half_max(rm, cm, 4, 2, 3 - ib);
                if (ipt + ib > 2 * NN) continue;
                BM(ipt + ib, 2) = h.set - h.clear;
            }
        }
        for (int m = 1; m <= 2; ++m) {  // 256 4-symbol correlations for each of 2 4-symbol groups
            float rm[16], cm[16];
            for (int t = 0; t < 16; ++t) rm[t] = cm[t] = -std::numeric_limits<float>::max();
            for (int i = 1; i <= 16; ++i)
                for (int j = 1; j <= 16; ++j) {
                    const int is = (i - 1) * 16 + j;
                    c4r[m][is] = c2r[i][2 * m - 1] + c2r[j][2 * m];
                    c4i[m][is] = c2i[i][2 * m - 1] + c2i[j][2 * m];
                    const float a = abs_exact(c4r[m][is], c4i[m][is]);
                    rm[i - 1] = std::max(rm[i - 1], a);
                    cm[j - 1] = std::max(cm[j - 1], a);
                }
            const int ipt = (k - 1) * 2 + 8 * (m - 1) + 1;
            for (int ib = 0; ib <= 7; ++ib) {
                const HalfMax h = half_max(rm, cm, 16, 4, 7 - ib);
                if (ipt + ib > 2 * NN) continue;
                BM(ipt + ib, 3) = h.set - h.clear;
            }
        }
        // 65536 8-symbol correlations: |c4(i,1) + c4(j,2)|, eight j at a time.
        const simd::V8 lowest = simd::set1(-std::numeric_limits<float>::max());
        for (int j = 0; j < 256; j += 8) simd::store(&colmax[j], lowest);
        for (int i = 1; i <= 256; ++i) {
            const simd::V8 ar = simd::set1(c4r[1][i]), ai = simd::set1(c4i[1][i]);
            simd::V8 rmax = lowest;
            for (int j = 0; j < 256; j += 8) {
                const simd::V8 re = simd::add(ar, simd::load(&c4r[2][j + 1]));
                const simd::V8 im = simd::add(ai, simd::load(&c4i[2][j + 1]));
                const simd::V8 a = simd::abs_exact(re, im);
                rmax = simd::max(rmax, a);
                simd::store(&colmax[j], simd::max(simd::load(&colmax[j]), a));
            }
            float t[8];
            simd::store(t, rmax);
            float r = t[0];
            for (int l = 1; l < 8; ++l) r = std::max(r, t[l]);
            rowmax[i - 1] = r;
        }
        const int ipt = (k - 1) * 2 + 1;
        for (int ib = 0; ib <= 15; ++ib) {
            const HalfMax h = half_max(rowmax, colmax, 256, 8, 15 - ib);
            if (ipt + ib > 2 * NN) continue;
            BM(ipt + ib, 4) = h.set - h.clear;
        }
    }

    static const int sw1[16] = {0, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0};
    static const int sw2[16] = {1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1};
    auto hbit = [&](int i) { return BM(i, 1) >= 0 ? 1 : 0; };
    auto count = [&](int start, const int* w) {
        int c = 0;
        for (int t = 0; t < 16; ++t)
            if (hbit(start + t) == w[t]) ++c;
        return c;
    };
    nsync_qual = count(1, sw1) + count(77, sw2) + count(153, sw1) + count(229, sw2) + count(305, sw1);
    if (nsync_qual < 46) {
        badsync = true;
        return;
    }
    for (int j = 1; j <= 4; ++j) normalizebmet(&BM(1, j), 2 * NN);
    const float scalefac = 2.83f;
    for (int i = 0; i < 2 * NN * 4; ++i) bitmetrics[i] = scalefac * bitmetrics[i];
}

}  // namespace fst4w
