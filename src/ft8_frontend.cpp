// FT8's front end (ft8_frontend.h): WSJT-X lib/ft8/sync8.f90,
// get_spectrum_baseline.f90, baseline.f90, ft8_downsample.f90, sync8d.f90,
// ft8b.f90 (with normalizebmet), subtractft8.f90, gen_ft8wave.f90,
// genft8.f90's get_ft8_tones_from_77bits, ft8_a7.f90's ft8_a7d, and
// lib/nuttal_window.f90, gfsk_pulse.f90, peakup.f90, pctile.f90, shell.f90,
// polyfit.f90, determ.f90 at 567ad29ce, line for line. Arrays keep the
// Fortran's indexing where that keeps the two readable side by side.
//
// Joe Taylor, K1JT, Steven Franke, K9AN, and the WSJT-X developers (GNU GPL
// v3), and so is this.
#include "ft8_frontend.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <cstdlib>

#include <fftw3.h>

#include <mutex>
#include <stdexcept>

#include "fft.h"         // fft_planner_lock
#include "ft8_maths.h"   // glibc 2.39's maths and x86-64 libmvec, carried
#include "ft8_ldpc.h"
#include "ft8_fortran.h"  // the Fortran helpers (shared with FT4)
#include "wsprd_osd.h"   // indexx (lib/indexx.f90)

namespace ft8 {

using namespace fortran;

namespace {

const int icos7[7] = {3, 1, 4, 0, 6, 5, 2};
const int graymap[8] = {0, 1, 3, 2, 5, 6, 4, 7};

}  // namespace

// ---------------------------------------------------------------- free routines

void get_ft8_tones_from_77bits(const int8_t msgbits[77], int itone[NN]) {
    int8_t codeword[174];
    encode174_91(msgbits, codeword);
    for (int i = 0; i < 7; ++i) {
        itone[i] = icos7[i];
        itone[36 + i] = icos7[i];
        itone[NN - 7 + i] = icos7[i];
    }
    int k = 7;
    for (int j = 1; j <= ND; ++j) {
        const int i = 3 * j - 2;
        ++k;
        if (j == 30) k += 7;
        const int indx = codeword[i - 1] * 4 + codeword[i] * 2 + codeword[i + 1];
        itone[k - 1] = graymap[indx];
    }
}

void normalizebmet(float* bmet, int n) {
    float s = 0.0f;
    for (int i = 0; i < n; ++i) s += bmet[i];
    const float bmetav = s / static_cast<float>(n);
    float s2 = 0.0f;
    for (int i = 0; i < n; ++i) s2 += bmet[i] * bmet[i];
    const float bmet2av = s2 / static_cast<float>(n);
    const float var = bmet2av - bmetav * bmetav;
    float bmetsig;
    if (var > 0.0f)
        bmetsig = std::sqrt(var);
    else
        bmetsig = std::sqrt(bmet2av);
    for (int i = 0; i < n; ++i) bmet[i] = bmet[i] / bmetsig;
}

float peakup(float ym, float y0, float yp) {
    const float b = yp - ym;
    const float c = (yp + ym) - 2.0f * y0;
    return -b / (2.0f * c);
}

float gfsk_pulse(float b, float t) {
    const float c = bitsf(0x40AAC42Bu);  // pi*sqrt(2.0/log(2.0)), folded
    const float cb = c * b;
    return 0.5f * (erff_glibc(cb * (t + 0.5f)) - erff_glibc(cb * (t - 0.5f)));
}

void nuttal_window(float* win, int n) {
    const float a0 = 0.3635819f, a1 = -0.4891775f, a2 = 0.1365995f, a3 = -0.0106411f;
    const float c2 = 2.0f * PI, c4 = 4.0f * PI, c6 = bitsf(0x4196CBE4u);  // 6*pi folded
    const float fn = static_cast<float>(n);
    // libmvec four at a time from i = 1, the rest scalar.
    int i = 1;
    for (; i + 3 <= n; i += 4) {
        float x1[4], x2[4], x3[4], y1[4], y2[4], y3[4];
        for (int q = 0; q < 4; ++q) {
            const float im1 = static_cast<float>(i + q - 1);
            x1[q] = (c2 * im1) / fn;
            x2[q] = (c4 * im1) / fn;
            x3[q] = (c6 * im1) / fn;
        }
        cosf4_glibc(x1, y1);
        cosf4_glibc(x2, y2);
        cosf4_glibc(x3, y3);
        for (int q = 0; q < 4; ++q) win[i + q - 1] = ((a0 + a1 * y1[q]) + a2 * y2[q]) + a3 * y3[q];
    }
    for (; i <= n; ++i) {
        const float im1 = static_cast<float>(i - 1);
        win[i - 1] = ((a0 + a1 * cosf_glibc((c2 * im1) / fn)) + a2 * cosf_glibc((c4 * im1) / fn)) +
                     a3 * cosf_glibc((c6 * im1) / fn);
    }
}

void baseline(float* s0, int nfa, int nfb, float* sbase0) {
    float* s = s0 - 1;  // 1-based
    float* sbase = sbase0 - 1;
    const int nseg = 10, npct = 10;
    const double df = 3.125;  // 12000.0/3840.0, a default real
    const int ia = std::max(1, static_cast<int>(std::lround(static_cast<double>(nfa) / df)));
    const int ib = static_cast<int>(std::lround(static_cast<double>(nfb) / df));
    // s(i) = 10.0*log10(s(i)): libmvec four at a time from ia, the rest
    // scalar.
    {
        int i = ia;
        for (; i + 3 <= ib; i += 4) {
            float l[4];
            log10f4_glibc(&s[i], l);
            for (int q = 0; q < 4; ++q) s[i + q] = 10.0f * l[q];
        }
        for (; i <= ib; ++i) s[i] = 10.0f * log10f_glibc(s[i]);
    }
    const int nterms = 5;
    const int nlen = (ib - ia + 1) / nseg;  // Length of test segment
    const int i0 = (ib - ia + 1) / 2;      // Midpoint
    static thread_local double x[1001], y[1001];
    int k = 0;
    for (int n = 1; n <= nseg; ++n) {  // Loop over all segments
        const int ja = ia + (n - 1) * nlen;
        const int jb = ja + nlen - 1;
        const float base = pctile(&s[ja] - 1, nlen, npct);  // Find lowest npct of points
        for (int i = ja; i <= jb; ++i) {
            if (s[i] <= base) {
                if (k < 1000) ++k;  // Save all "lower envelope" points
                x[k] = static_cast<double>(i - i0);
                y[k] = static_cast<double>(s[i]);
            }
        }
    }
    const int kz = k;
    double a[6] = {};
    polyfit(x, y, kz, nterms, a);  // Fit a low-order polynomial
    for (int i = 1; i <= 1920; ++i) sbase[i] = 0.0f;
    for (int i = ia; i <= ib; ++i) {
        const double t = static_cast<double>(i - i0);
        sbase[i] = static_cast<float>(a[1] + t * (a[2] + t * (a[3] + t * (a[4] + t * (a[5])))) +
                                      static_cast<double>(0.65f));
    }
}

float xbase_of(const float* sbase, float f1, bool a7) {
    int i = nint(f1 / 3.125f);
    if (a7) i = std::max(1, i);
    return powf_glibc(10.0f, 0.1f * (sbase[i - 1] - 40.0f));
}

// ---------------------------------------------------------------- Frontend

Frontend::Frontend()
    : s_(static_cast<size_t>(NH1) * NHSYM),
      x3840_(NFFT1 + 2),
      ds_x_(192000 + 2),
      c3200_(3200),
      cd0_(3200),
      csymb_(32),
      cref_(NZ),
      camp_(NMAX),
      cfilt_(NMAX),
      cw_(NMAX),
      st_dd_(NMAX),
      st_x_(NMAX + 2 + 4) {}

size_t Frontend::heap_bytes() const {
    auto b = [](const auto& v) { return v.capacity() * sizeof(v[0]); };
    return b(gsb_window_) + b(s_) + b(x3840_) + b(ds_x_) + b(c3200_) + b(cd0_) + b(csymb_) + b(cref_) + b(camp_) +
           b(cfilt_) + b(cw_) + b(st_dd_) + b(st_x_) + b(pulse_) + b(ctab_) + b(dphi_) + b(pdphi_) + b(gw_phi_) +
           b(gw_ix_) + b(gw_c_);
}

Frontend::~Frontend() {
    if (st_plan_) {
        std::lock_guard<std::mutex> lock(skimmer::fft_planner_lock());
        fftwf_destroy_plan(static_cast<fftwf_plan>(st_plan_));
    }
}

void Frontend::trace_flush() {
    if (trace_ && trace_nch_ > 0) trace_("g " + I(trace_nch_) + " " + X(trace_ch_));
    trace_nch_ = 0;
    trace_ch_ = 2166136261u;
}

void Frontend::get_spectrum_baseline(const float* dd, int& nfa, int& nfb, float* sbase) {
    constexpr int NST = NFFT1 / 2, NF = 93;
    if (gsb_first_) {
        gsb_first_ = false;
        gsb_window_.assign(NFFT1, 0.0f);
        nuttal_window(gsb_window_.data(), NFFT1);
        float sumw = 0.0f;
        for (int i = 0; i < NFFT1; ++i) sumw += gsb_window_[i];
        for (int i = 0; i < NFFT1; ++i)
            gsb_window_[i] = ((gsb_window_[i] / sumw) * static_cast<float>(NSPS)) * 2.0f / 300.0f;
    }
    float savg[NH1];
    std::fill(savg, savg + NH1, 0.0f);
    float* x = x3840_.data();
    for (int j = 1; j <= NF; ++j) {
        const int ia = (j - 1) * NST + 1;
        const int ib = ia + NFFT1 - 1;
        if (ib > NMAX) break;
        for (int i = 0; i < NFFT1; ++i) x[i] = dd[ia - 1 + i] * gsb_window_[i];
        fst4w::fft_r2c_inplace(x, NFFT1);  // r2c FFT
        for (int i = 1; i <= NH1; ++i) {
            const float a = hypotf_glibc(x[2 * i], x[2 * i + 1]);
            savg[i - 1] = savg[i - 1] + a * a;  // Average spectrum
        }
    }
    const int nwin = nfb - nfa;
    if (nfa < 100) {
        nfa = 100;
        if (nwin < 100) nfb = nfa + nwin;  // nagain
    }
    if (nfb > 4910) {
        nfb = 4910;
        if (nwin < 100) nfa = nfb - nwin;
    }
    if (trace_) trace_("b " + X(h_r(savg, NH1)) + " " + I(nfa) + " " + I(nfb));
    baseline(savg, nfa, nfb, sbase);
    if (trace_) trace_("B " + X(h_r(savg, NH1)) + " " + X(h_r(sbase, NH1)));
}

void Frontend::sync8(const float* dd, int npts, int& nfa, int& nfb, float syncmin, int nfqso, int maxcand,
                     float* candidate, int& ncand, float* sbase) {
    constexpr int MAXPRECAND = 1000, JZ = 62;
    if (trace_) {
        trace_flush();
        trace_("S " + X(h_r(dd, static_cast<size_t>(npts))) + " " + I(npts) + " " + I(nfa) + " " + I(nfb) + " " +
               XF(syncmin) + " " + I(nfqso) + " " + I(maxcand));
    }
    auto S = [&](int i, int j) -> float& { return s_[static_cast<size_t>(j - 1) * NH1 + (i - 1)]; };
    // Compute symbol spectra, stepping by NSTEP steps.
    const float tstep = static_cast<float>(NSTEP) / 12000.0f;
    const float df = 12000.0f / static_cast<float>(NFFT1);  // 3.125 Hz
    const float fac = 1.0f / 300.0f;
    float* x = x3840_.data();
    for (int j = 1; j <= NHSYM; ++j) {
        const int ia = (j - 1) * NSTEP + 1;
        for (int i = 0; i < NSPS; ++i) x[i] = fac * dd[ia - 1 + i];
        std::fill(x + NSPS, x + NFFT1 + 2, 0.0f);
        fst4w::fft_r2c_inplace(x, NFFT1);  // r2c FFT
        for (int i = 1; i <= NH1; ++i) S(i, j) = x[2 * i] * x[2 * i] + x[2 * i + 1] * x[2 * i + 1];
    }
    if (trace_) trace_("s " + X(h_r(s_.data(), s_.size())));
    get_spectrum_baseline(dd, nfa, nfb, sbase);

    const int ia = std::max(1, nint(static_cast<float>(nfa) / df));
    const int ib = nint(static_cast<float>(nfb) / df);
    const int nssy = NSPS / NSTEP;  // # steps per symbol
    const int nfos = NFFT1 / NSPS;  // # frequency bin oversampling factor
    const int jstrt = static_cast<int>(0.5f / tstep);
    static thread_local float candidate0[3 * MAXPRECAND];
    auto C0 = [&](int i, int k) -> float& { return candidate0[(k - 1) * 3 + (i - 1)]; };
    std::fill(candidate0, candidate0 + 3 * MAXPRECAND, 0.0f);
    int k = 0;

    std::vector<float> sync2d(static_cast<size_t>(NH1) * (2 * JZ + 1));
    auto SY = [&](int i, int j) -> float& { return sync2d[static_cast<size_t>(j + JZ) * NH1 + (i - 1)]; };
    for (int i = ia; i <= ib; ++i) {
        for (int j = -JZ; j <= JZ; ++j) {
            float ta = 0, tb = 0, tc = 0, t0a = 0, t0b = 0, t0c = 0;
            for (int n = 0; n <= 6; ++n) {
                const int m = j + jstrt + nssy * n;
                if (m >= 1 && m <= NHSYM) {
                    ta = ta + S(i + nfos * icos7[n], m);
                    float sm = 0.0f;
                    for (int q = 0; q <= 6; ++q) sm += S(i + nfos * q, m);
                    t0a = t0a + sm;
                }
                tb = tb + S(i + nfos * icos7[n], m + nssy * 36);
                {
                    float sm = 0.0f;
                    for (int q = 0; q <= 6; ++q) sm += S(i + nfos * q, m + nssy * 36);
                    t0b = t0b + sm;
                }
                if (m + nssy * 72 <= NHSYM) {
                    tc = tc + S(i + nfos * icos7[n], m + nssy * 72);
                    float sm = 0.0f;
                    for (int q = 0; q <= 6; ++q) sm += S(i + nfos * q, m + nssy * 72);
                    t0c = t0c + sm;
                }
            }
            float t = (ta + tb) + tc;
            float t0 = (t0a + t0b) + t0c;
            t0 = (t0 - t) / 6.0f;
            const float sync_abc = t / t0;
            t = tb + tc;
            t0 = t0b + t0c;
            t0 = (t0 - t) / 6.0f;
            const float sync_bc = t / t0;
            // max(sync_abc, sync_bc): maxss, sync_bc when either is NaN.
            SY(i, j) = sync_abc > sync_bc ? sync_abc : sync_bc;
        }
    }

    std::vector<float> red(NH1 + 1, 0.0f), red2(NH1 + 1, 0.0f);  // 1-based
    std::vector<int> jpeak(NH1 + 1, 0), jpeak2(NH1 + 1, 0);
    const int mlag = 13, mlag2 = JZ;
    for (int i = ia; i <= ib; ++i) {
        jpeak[i] = maxloc(-mlag, mlag, [&](int j) { return SY(i, j); });
        red[i] = SY(i, jpeak[i]);
        jpeak2[i] = maxloc(-mlag2, mlag2, [&](int j) { return SY(i, j); });
        red2[i] = SY(i, jpeak2[i]);
    }
    const int iz = ib - ia + 1;
    std::vector<int> indx(NH1 + 1), indx2(NH1 + 1);  // 1-based
    wsprd::indexx(&red[ia], iz, &indx[1]);
    const int npctile = nint(0.40f * static_cast<float>(iz));
    if (npctile < 1) {  // something is wrong; bail out
        ncand = 0;
        return;
    }
    int ibase = indx[npctile] - 1 + ia;
    if (ibase < 1) ibase = 1;
    if (ibase > NH1) ibase = NH1;
    const float base = red[ibase];
    for (int i = 1; i <= NH1; ++i) red[i] = red[i] / base;
    wsprd::indexx(&red2[ia], iz, &indx2[1]);
    int ibase2 = indx2[npctile] - 1 + ia;
    if (ibase2 < 1) ibase2 = 1;
    if (ibase2 > NH1) ibase2 = NH1;
    const float base2 = red2[ibase2];
    for (int i = 1; i <= NH1; ++i) red2[i] = red2[i] / base2;
    if (trace_) {
        uint32_t h = 2166136261u;
        for (int j = -JZ; j <= JZ; ++j)
            for (int i = ia; i <= ib; ++i) h = hw(h, fbits(SY(i, j)));
        trace_("y " + X(h) + " " + X(h_r(&red[1], NH1)) + " " + X(h_r(&red2[1], NH1)) + " " +
               X(h_i(&jpeak[ia], static_cast<size_t>(iz))) + " " + X(h_i(&jpeak2[ia], static_cast<size_t>(iz))) +
               " " + I(ibase) + " " + I(ibase2) + " " + XF(base) + " " + XF(base2));
    }
    for (int i = 1; i <= std::min(MAXPRECAND, iz); ++i) {
        const int n = ia + indx[iz + 1 - i] - 1;
        if (k >= MAXPRECAND) break;
        if (red[n] >= syncmin && !std::isnan(red[n])) {
            ++k;
            C0(1, k) = static_cast<float>(n) * df;
            C0(2, k) = (static_cast<float>(jpeak[n]) - 0.5f) * tstep;
            C0(3, k) = red[n];
        }
        if (std::abs(jpeak2[n] - jpeak[n]) == 0) continue;
        if (k >= MAXPRECAND) break;
        if (red2[n] >= syncmin && !std::isnan(red2[n])) {
            ++k;
            C0(1, k) = static_cast<float>(n) * df;
            C0(2, k) = (static_cast<float>(jpeak2[n]) - 0.5f) * tstep;
            C0(3, k) = red2[n];
        }
    }
    ncand = k;

    // Save only the best of near-dupe freqs.
    for (int i = 1; i <= ncand; ++i) {
        if (i >= 2) {
            for (int j = 1; j <= i - 1; ++j) {
                const float fdiff = std::fabs(C0(1, i)) - std::fabs(C0(1, j));
                const float tdiff = std::fabs(C0(2, i) - C0(2, j));
                if (std::fabs(fdiff) < 4.0f && tdiff < 0.04f) {
                    if (C0(3, i) >= C0(3, j)) C0(3, j) = 0.0f;
                    if (C0(3, i) < C0(3, j)) C0(3, i) = 0.0f;
                }
            }
        }
    }
    if (trace_) trace_("K " + I(ncand) + " " + X(h_r(candidate0, 3 * static_cast<size_t>(ncand))));
    // (fac=20.0/maxval(s); s=fac*s: s is not used again.)

    // Sort by sync
    std::vector<float> candidate_sync(static_cast<size_t>(ncand) + 1);
    for (int i = 1; i <= ncand; ++i) candidate_sync[i] = C0(3, i);
    wsprd::indexx(&candidate_sync[1], ncand, &indx[1]);
    auto CA = [&](int i, int kk) -> float& { return candidate[(kk - 1) * 3 + (i - 1)]; };
    // Place candidates within 10 Hz of nfqso at the top of the list
    k = 1;
    for (int i = 1; i <= ncand; ++i) {
        if (std::fabs(C0(1, i) - static_cast<float>(nfqso)) <= 10.0f && C0(3, i) >= syncmin) {
            for (int q = 1; q <= 3; ++q) CA(q, k) = C0(q, i);
            C0(3, i) = 0.0f;
            ++k;
        }
    }
    for (int i = ncand; i >= 1; --i) {
        const int j = indx[i];
        if (C0(3, j) >= syncmin) {
            CA(2, k) = C0(2, j);
            CA(3, k) = C0(3, j);
            CA(1, k) = std::fabs(C0(1, j));
            ++k;
            if (k > maxcand) break;
        }
    }
    ncand = k - 1;
    if (trace_) trace_("k " + I(ncand) + " " + X(h_r(candidate, 3 * static_cast<size_t>(ncand))));
}

void Frontend::ft8_downsample(const float* dd, bool& newdat, float f0, cfloat* c1) {
    constexpr int NFFT1D = 192000, NFFT2 = 3200;
    if (ds_first_) {
        // taper(0) = 1 (folded); 1..100 by libmvec, four at a time.
        taper_[0] = 1.0f;
        for (int i = 1; i <= 100; i += 4) {
            float a[4], c[4];
            for (int q = 0; q < 4; ++q) a[q] = (static_cast<float>(i + q) * PI) / 100.0f;
            cosf4_glibc(a, c);
            for (int q = 0; q < 4; ++q) taper_[i + q] = (1.0f + c[q]) * 0.5f;
        }
        ds_first_ = false;
    }
    float* x = ds_x_.data();
    if (newdat) {
        // Data in dd have changed, recompute the long FFT
        std::memcpy(x, dd, sizeof(float) * NMAX);
        std::fill(x + NMAX, x + NFFT1D + 2, 0.0f);  // Zero-pad the x array
        fst4w::fft_r2c_inplace(x, NFFT1D);         // r2c FFT to freq domain
        if (trace_) {
            trace_flush();
            trace_("X " + X(h_r(x, NFFT1D + 2)));
        }
        newdat = false;
    }
    const float df = 12000.0f / static_cast<float>(NFFT1D);
    const float baud = 12000.0f / static_cast<float>(NSPS);
    const int i0 = nint(f0 / df);
    const float ft = f0 + 8.5f * baud;
    const int it = std::min(nint(ft / df), NFFT1D / 2);
    const float fb = f0 - 1.5f * baud;
    const int ib = std::max(1, nint(fb / df));
    int k = 0;
    cfloat* c = c3200_.data();
    std::fill(c, c + NFFT2, cfloat(0.0f, 0.0f));
    for (int i = ib; i <= it; ++i) {
        c[k] = cfloat(x[2 * i], x[2 * i + 1]);
        ++k;
    }
    for (int m = 0; m <= 100; ++m) c[m] = cscale(c[m], taper_[100 - m]);
    for (int m = 0; m <= 100; ++m) c[k - 1 - 100 + m] = cscale(c[k - 1 - 100 + m], taper_[m]);
    // c1 = cshift(c1, i0-ib)
    int sh = (i0 - ib) % NFFT2;
    if (sh < 0) sh += NFFT2;
    for (int m = 0; m < NFFT2; ++m) c1[m] = c[(m + sh) % NFFT2];
    fst4w::fft_c2c_inplace(c1, NFFT2, 1);  // c2c FFT back to time domain
    const float fac = bitsf(0x38293696u);   // 1.0/sqrt(float(NFFT1)*NFFT2), folded
    for (int m = 0; m < NFFT2; ++m) c1[m] = cscale(c1[m], fac);
}

float Frontend::sync8d(const cfloat* cd0, int i0, const cfloat* ctwk, int itwk) {
    // Set some constants and compute the csync array.
    if (s8d_first_) {
        s8d_twopi_ = TWOPI;
        for (int i = 0; i <= 6; ++i) {
            float phi = 0.0f;
            const float dphi = s8d_twopi_ * static_cast<float>(icos7[i]) / 32.0f;
            for (int j = 1; j <= 32; ++j) {
                csync_[i][j - 1] = cis(phi);  // Waveform for 7x7 Costas array
                phi = std::fmod(phi + dphi, s8d_twopi_);
            }
        }
        s8d_first_ = false;
    }
    float sync = 0;
    cfloat csync2[32];
    for (int i = 0; i <= 6; ++i) {  // Sum over 7 Costas frequencies and
        const int i1 = i0 + i * 32;   // three Costas arrays
        const int i2 = i1 + 36 * 32;
        const int i3 = i1 + 72 * 32;
        for (int j = 0; j < 32; ++j) csync2[j] = csync_[i][j];
        if (itwk == 1)
            for (int j = 0; j < 32; ++j) csync2[j] = cmul(ctwk[j], csync2[j]);  // Tweak the frequency
        auto zsum = [&](int a) {
            float re = 0.0f, im = 0.0f;
            for (int j = 0; j < 32; ++j) {
                const cfloat p = cmul(cd0[a + j], cfloat(csync2[j].real(), -csync2[j].imag()));
                re = re + p.real();
                im = im + p.imag();
            }
            return cfloat(re, im);
        };
        cfloat z1(0.0f, 0.0f), z2(0.0f, 0.0f), z3(0.0f, 0.0f);
        if (i1 >= 0 && i1 + 31 <= NP2 - 1) z1 = zsum(i1);
        if (i2 >= 0 && i2 + 31 <= NP2 - 1) z2 = zsum(i2);
        if (i3 >= 0 && i3 + 31 <= NP2 - 1) z3 = zsum(i3);
        const float p1 = z1.real() * z1.real() + z1.imag() * z1.imag();
        const float p2 = z2.real() * z2.real() + z2.imag() * z2.imag();
        const float p3 = z3.real() * z3.real() + z3.imag() * z3.imag();
        sync = ((sync + p1) + p2) + p3;
    }
    return sync;
}


void Frontend::fine_sync(const float* dd0, bool& newdat, float& f1, float& xdt, Sync& sy, bool a7) {
    const float fs2 = 12000.0f / static_cast<float>(NDOWN);
    const float dt2 = 1.0f / fs2;
    const float twopi = TWOPI;
    float delfbest = 0.0f;
    int ibest = 0;
    cfloat* cd0 = cd0_.data();
    cfloat ctwk[32];

    ft8_downsample(dd0, newdat, f1, cd0);  // Mix f1 to baseband and downsample
    if (trace_) sy.h_cd0 = h_c(cd0, 3200);

    const int i0 = nint((xdt + 0.5f) * fs2);  // Initial guess for start of signal
    float smax = 0.0f;
    for (int idt = i0 - 10; idt <= i0 + 10; ++idt) {  // Search over +/- one quarter symbol
        const float sync = sync8d(cd0, idt, ctwk, 0);   // (ctwk not used here)
        if (sync > smax) {
            smax = sync;
            ibest = idt;
        }
    }
    sy.ibest0 = ibest;

    // Now peak up in frequency
    smax = 0.0f;
    for (int ifr = -5; ifr <= 5; ++ifr) {  // Search over +/- 2.5 Hz
        const float delf = static_cast<float>(ifr) * 0.5f;
        const float dphi = (twopi * delf) * dt2;
        float phi = 0.0f;
        for (int i = 0; i < 32; ++i) {
            ctwk[i] = cis(phi);
            phi = std::fmod(phi + dphi, twopi);
        }
        const float sync = sync8d(cd0, ibest, ctwk, 1);
        if (sync > smax) {
            smax = sync;
            delfbest = delf;
        }
    }
    // (twkfreq1(cd0, NP2, fs2, a, cd0): its output is overwritten by the
    // downsample below.)
    f1 = f1 + delfbest;  // Improved estimate of DF

    bool no = false;
    ft8_downsample(dd0, no, f1, cd0);  // Mix f1 to baseband and downsample

    for (int idt = -4; idt <= 4; ++idt) sy.ss[idt + 4] = sync8d(cd0, ibest + idt, ctwk, 0);
    smax = maxval(sy.ss, 9);
    const int iloc = maxloc(1, 9, [&](int q) { return sy.ss[q - 1]; });
    ibest = iloc - 5 + ibest;
    xdt = static_cast<float>(ibest - 1) * dt2;
    if (a7) xdt = xdt - 0.5f;

    cfloat* csymb = csymb_.data();
    for (int k = 1; k <= NN; ++k) {
        const int i1 = ibest + (k - 1) * 32;
        if (i1 >= 0 && i1 + 31 <= NP2 - 1)
            for (int q = 0; q < 32; ++q) csymb[q] = cd0[i1 + q];
        else
            for (int q = 0; q < 32; ++q) csymb[q] = cfloat(0.0f, 0.0f);
        fst4w::fft_c2c_inplace(csymb, 32, -1);
        for (int q = 0; q < 8; ++q) {
            cs_[k - 1][q] = cfloat(csymb[q].real() / 1e3f, csymb[q].imag() / 1e3f);
            s8_[k - 1][q] = cabs_(csymb[q]);
        }
    }

    // sync quality check
    int is1 = 0, is2 = 0, is3 = 0;
    for (int k = 1; k <= 7; ++k) {
        int ip = maxloc(1, 8, [&](int q) { return s8_[k - 1][q - 1]; });
        if (icos7[k - 1] == ip - 1) ++is1;
        ip = maxloc(1, 8, [&](int q) { return s8_[k + 36 - 1][q - 1]; });
        if (icos7[k - 1] == ip - 1) ++is2;
        ip = maxloc(1, 8, [&](int q) { return s8_[k + 72 - 1][q - 1]; });
        if (icos7[k - 1] == ip - 1) ++is3;
    }
    sy.nsync = is1 + is2 + is3;  // hard sync sum - max is 21
    sy.ibest = ibest;
    sy.delfbest = delfbest;
    sy.smax = smax;
}

void Frontend::trace_sync(const Sync& sy, float xdt, bool chain, float f1_in, float xdt_in, float xbase,
                          bool newdat_in) {
    const uint32_t h2 = h_c(cd0_.data(), 3200), hss = h_r(sy.ss, 9), hcs = h_c(&cs_[0][0], 8 * NN),
                   hs8 = h_r(&s8_[0][0], 8 * NN);
    if (chain) {
        if (trace_nch_ == 0) trace_ch_ = 2166136261u;
        ++trace_nch_;
        const uint32_t w[] = {fbits(f1_in), fbits(xdt_in),   fbits(xbase), newdat_in ? 1u : 0u,
                              static_cast<uint32_t>(sy.ibest0), fbits(sy.delfbest), static_cast<uint32_t>(sy.ibest),
                              fbits(xdt),   fbits(sy.smax),   sy.h_cd0,     h2,
                              hss,          hcs,              hs8,          static_cast<uint32_t>(sy.nsync)};
        for (uint32_t x : w) trace_ch_ = hw(trace_ch_, x);
        return;
    }
    trace_("G " + I(sy.ibest0) + " " + XF(sy.delfbest) + " " + I(sy.ibest) + " " + XF(xdt) + " " + XF(sy.smax) +
           " " + X(sy.h_cd0) + " " + X(h2) + " " + X(hss) + " " + X(hcs) + " " + X(hs8) + " " + I(sy.nsync));
}

void Frontend::trace_metrics(int nll, float* const* v, bool after) {
    std::string l = after ? "L" : "M";
    if (!after) {
        ++trace_nl_;
        trace_smp_ = trace_nl_ % 499 == 1;
        if (!trace_smp_) return;
    }
    for (int i = 0; i < nll; ++i) l += " " + X(h_r(v[i], 174));
    trace_(l);
    if (after && trace_smp_) {
        std::string f = "l";
        for (int i = 0; i < nll; ++i)
            for (int k = 0; k < 174; ++k) f += " " + XF(v[i][k]);
        trace_(f);
    }
}

void Frontend::bit_metrics(int imetric, float* bmeta, float* bmetb, float* bmetc, float* bmetd, float* bmete) {
    float s2[512];
    for (int nsym = 1; nsym <= 3; ++nsym) {
        const int nt = 1 << (3 * nsym);
        for (int ihalf = 1; ihalf <= 2; ++ihalf) {
            for (int k = 1; k <= 29; k += nsym) {
                const int ks = ihalf == 1 ? k + 7 : k + 43;
                for (int i = 0; i <= nt - 1; ++i) {
                    const int i1 = i / 64;
                    const int i2 = (i & 63) / 8;
                    const int i3 = i & 7;
                    if (nsym == 1) {
                        s2[i] = cabs_(cs_[ks - 1][graymap[i3]]);
                    } else if (nsym == 2) {
                        s2[i] = cabs_(cadd(cs_[ks - 1][graymap[i2]], cs_[ks][graymap[i3]]));
                    } else {
                        s2[i] = cabs_(cadd(cadd(cs_[ks - 1][graymap[i1]], cs_[ks][graymap[i2]]),
                                           cs_[ks + 1][graymap[i3]]));
                    }
                }
                if (imetric == 2)
                    for (int i = 0; i < nt; ++i) s2[i] = s2[i] * s2[i];
                const int i32 = 1 + (k - 1) * 3 + (ihalf - 1) * 87;
                const int ibmax = nsym == 1 ? 2 : nsym == 2 ? 5 : 8;
                for (int ib = 0; ib <= ibmax; ++ib) {
                    const float m1 = maxval_bit(s2, nt, ibmax - ib, true);
                    const float m0 = maxval_bit(s2, nt, ibmax - ib, false);
                    const float bm = m1 - m0;
                    if (i32 + ib > 174) continue;
                    if (nsym == 1) {
                        bmeta[i32 + ib - 1] = bm;
                        const float den = m1 > m0 ? m1 : m0;
                        float cm;
                        if (den > 0.0f)
                            cm = bm / den;
                        else  // erase it
                            cm = 0.0f;
                        bmetd[i32 + ib - 1] = cm;
                    } else if (nsym == 2) {
                        bmetb[i32 + ib - 1] = bm;
                    } else {
                        bmetc[i32 + ib - 1] = bm;
                    }
                }
            }
        }
    }
    if (bmete) {
        for (int i = 0; i < 174; ++i) {
            const float temp[3] = {bmeta[i], bmetb[i], bmetc[i]};
            const int ip = maxloc(1, 3, [&](int q) { return std::fabs(temp[q - 1]); });
            bmete[i] = temp[ip - 1];
        }
    }
}

void Frontend::ft8b(float* dd0, bool& newdat, int nQSOProgress, int nfqso, int nftx, int ndepth, int nzhsym,
                    bool lapon, bool lapcqonly, int napwid, bool lsubtract, bool nagain, int ncontest, int imetric,
                    int& iaptype, const std::string& mycall12, const std::string& hiscall12, float& f1, float& xdt,
                    float xbase, const int apsym[58], const int aph10[10], int& nharderrors, float& dmin,
                    int& nbadcrc, int& ipass, std::string& msg37, float& xsnr, int itone[NN], Codec& codec) {
    // The a-priori tables, mapped from 0/1 to -1/+1 "if(first.or.(ncontest.ne.ncontest0))": SAVEd as
    // they are DATA, so each later change of ncontest maps them again (-1 to -3, ...), as there.
    static const int nappasses[6] = {2, 2, 2, 4, 4, 3};
    static const int naptypes[6][4] = {{1, 2, 0, 0}, {2, 3, 0, 0}, {2, 3, 0, 0}, {3, 4, 5, 6}, {3, 4, 5, 6}, {3, 1, 2, 0}};
    if (ft8b_first_ || ncontest != ncontest0_) {
        for (int* m : {mcq_, mcqfd_, mcqru_, mcqtest_, mcqww_})
            for (int i = 0; i < 29; ++i) m[i] = 2 * m[i] - 1;
        for (int* m : {mrrr_, m73_, mrr73_})
            for (int i = 0; i < 19; ++i) m[i] = 2 * m[i] - 1;
        ft8b_first_ = false;
        ncontest0_ = ncontest;
    }
    const int *mcq = mcq_, *mcqru = mcqru_, *mcqfd = mcqfd_, *mcqtest = mcqtest_, *mcqww = mcqww_;
    const int *mrrr = mrrr_, *m73 = m73_, *mrr73 = mrr73_;
    auto pm = [](int v) { return static_cast<float>(v); };

    codec.set_calls(mycall12, hiscall12);  // dxcall13=hiscall12; mycall13=mycall12

    nharderrors = -1;
    const float f1_in = f1, xdt_in = xdt;
    const bool newdat_in = newdat;
    Sync sy;
    fine_sync(dd0, newdat, f1, xdt, sy, false);
    const int nsync = sy.nsync;
    int syncmin = 6;
    if (imetric == 2) syncmin = 7;
    if (ndepth <= 2) syncmin = 8;
    if (trace_) {
        if (nsync <= syncmin) {
            trace_sync(sy, xdt, true, f1_in, xdt_in, xbase, newdat_in);
        } else {
            trace_flush();
            trace_("F " + XF(f1_in) + " " + XF(xdt_in) + " " + XF(xbase) + " " + L(newdat_in));
            trace_sync(sy, xdt, false, 0, 0, 0, false);
        }
    }
    if (nsync <= syncmin) {  // bail out
        nbadcrc = 1;
        if (trace_) trace_ch_ = hw(trace_ch_, static_cast<uint32_t>(nbadcrc));
        return;
    }

    float bmeta[174], bmetb[174], bmetc[174], bmetd[174], bmete[174];
    bit_metrics(imetric, bmeta, bmetb, bmetc, bmetd, bmete);
    float* const bm[5] = {bmeta, bmetb, bmetc, bmetd, bmete};
    if (trace_) trace_metrics(5, bm, false);
    normalizebmet(bmeta, 174);
    normalizebmet(bmetb, 174);
    normalizebmet(bmetc, 174);
    normalizebmet(bmetd, 174);
    normalizebmet(bmete, 174);

    const float scalefac = 2.83f;
    float llra[174], llrb[174], llrc[174], llrd[174], llre[174];
    for (int i = 0; i < 174; ++i) {
        llra[i] = scalefac * bmeta[i];
        llrb[i] = scalefac * bmetb[i];
        llrc[i] = scalefac * bmetc[i];
        llrd[i] = scalefac * bmetd[i];
        llre[i] = scalefac * bmete[i];
    }
    float* const ll[5] = {llra, llrb, llrc, llrd, llre};
    if (trace_) trace_metrics(5, ll, true);

    int npasses;
    if (lapon || ncontest == 7) {  // Hounds always use AP
        if (!lapcqonly)
            npasses = 5 + 2 * nappasses[nQSOProgress];
        else
            npasses = 7;
    } else {
        npasses = 5;
    }
    if (nzhsym < 50) npasses = 5;

    float llrz[174];
    int8_t apmask[174] = {};
    int8_t message77[77] = {}, message91[91] = {}, cw[174];
    float apmag = 0.0f;
    for (ipass = 1; ipass <= npasses; ++ipass) {
        std::memcpy(llrz, llra, sizeof llrz);
        if (ipass == 2) std::memcpy(llrz, llrb, sizeof llrz);
        if (ipass == 3) std::memcpy(llrz, llrc, sizeof llrz);
        if (ipass == 4) std::memcpy(llrz, llrd, sizeof llrz);
        if (ipass == 5) std::memcpy(llrz, llre, sizeof llrz);
        if (ipass <= 5) {
            std::memset(apmask, 0, sizeof apmask);
            iaptype = 0;
        }
        if (ipass > 5) {
            std::memcpy(llrz, llra, sizeof llrz);
            if ((ipass - 5) % 2 == 1) std::memcpy(llrz, llra, sizeof llrz);
            if ((ipass - 5) % 2 == 0) std::memcpy(llrz, llrc, sizeof llrz);
            float a[174];
            for (int i = 0; i < 174; ++i) a[i] = std::fabs(llrz[i]);
            apmag = maxval(a, 174) * 1.1f;
            if (!lapcqonly)
                iaptype = naptypes[nQSOProgress][(ipass - 4) / 2 - 1];
            else
                iaptype = 1;

            // Conditions that cause us to bail out of AP decoding
            if (ncontest <= 5 && iaptype >= 3 &&
                (std::fabs(f1 - static_cast<float>(nfqso)) > static_cast<float>(napwid) &&
                 std::fabs(f1 - static_cast<float>(nftx)) > static_cast<float>(napwid)))
                continue;
            if (ncontest == 6) continue;                   // No AP for Foxes
            if (ncontest == 7 && f1 > 950.0f) continue;     // Hounds use AP only for signals below 950 Hz
            if (iaptype >= 2 && apsym[0] > 1) continue;     // No, or nonstandard, mycall
            if (ncontest == 7 && iaptype >= 2 && aph10[0] > 1) continue;
            if (iaptype >= 3 && apsym[29] > 1) continue;    // No, or nonstandard, dxcall

            auto setm = [&](int a1, int a2) {  // apmask(a1:a2)=1
                for (int i = a1; i <= a2; ++i) apmask[i - 1] = 1;
            };
            auto ap = [&](int i) { return apmag * static_cast<float>(apsym[i - 1]); };
            if (iaptype == 1) {  // CQ or CQ RU or CQ TEST or CQ FD
                std::memset(apmask, 0, sizeof apmask);
                setm(1, 29);
                const int* m = nullptr;
                if (ncontest == 0) m = mcq;
                if (ncontest == 1) m = mcqtest;
                if (ncontest == 2) m = mcqtest;
                if (ncontest == 3) m = mcqfd;
                if (ncontest == 4) m = mcqru;
                if (ncontest == 5) m = mcqww;
                if (ncontest == 7) m = mcq;
                if (ncontest == 8) m = mcqtest;
                if (m)
                    for (int i = 1; i <= 29; ++i) llrz[i - 1] = apmag * pm(m[i - 1]);
                setm(75, 77);
                llrz[74] = apmag * (-1.0f);
                llrz[75] = apmag * (-1.0f);
                llrz[76] = apmag * (+1.0f);
            }

            if (iaptype == 2) {  // MyCall,???,???
                std::memset(apmask, 0, sizeof apmask);
                if (ncontest == 0 || ncontest == 1 || ncontest == 5 || ncontest == 8) {
                    setm(1, 29);
                    for (int i = 1; i <= 29; ++i) llrz[i - 1] = ap(i);
                    setm(75, 77);
                    llrz[74] = apmag * (-1.0f);
                    llrz[75] = apmag * (-1.0f);
                    llrz[76] = apmag * (+1.0f);
                } else if (ncontest == 2) {
                    setm(1, 28);
                    for (int i = 1; i <= 28; ++i) llrz[i - 1] = ap(i);
                    setm(72, 74);
                    llrz[71] = apmag * (-1.0f);
                    llrz[72] = apmag * (+1.0f);
                    llrz[73] = apmag * (-1.0f);
                    setm(75, 77);
                    for (int i = 75; i <= 77; ++i) llrz[i - 1] = apmag * (-1.0f);
                } else if (ncontest == 3) {
                    setm(1, 28);
                    for (int i = 1; i <= 28; ++i) llrz[i - 1] = ap(i);
                    setm(75, 77);
                    for (int i = 75; i <= 77; ++i) llrz[i - 1] = apmag * (-1.0f);
                } else if (ncontest == 4) {
                    setm(2, 29);
                    for (int i = 2; i <= 29; ++i) llrz[i - 1] = ap(i - 1);
                    setm(75, 77);
                    llrz[74] = apmag * (-1.0f);
                    llrz[75] = apmag * (+1.0f);
                    llrz[76] = apmag * (+1.0f);
                } else if (ncontest == 7) {  // ??? RR73; MyCall <Fox Call hash10> ???
                    setm(29, 56);
                    for (int i = 29; i <= 56; ++i) llrz[i - 1] = ap(i - 28);
                    setm(57, 66);
                    for (int i = 57; i <= 66; ++i) llrz[i - 1] = apmag * static_cast<float>(aph10[i - 57]);
                    setm(72, 77);
                    llrz[71] = apmag * (-1.0f);
                    llrz[72] = apmag * (-1.0f);
                    llrz[73] = apmag * (+1.0f);
                    for (int i = 75; i <= 77; ++i) llrz[i - 1] = apmag * (-1.0f);
                }
            }

            if (iaptype == 3) {  // MyCall,DxCall,???
                std::memset(apmask, 0, sizeof apmask);
                if (ncontest == 0 || ncontest == 1 || ncontest == 2 || ncontest == 5 || ncontest == 7 ||
                    ncontest == 8) {
                    setm(1, 58);
                    for (int i = 1; i <= 58; ++i) llrz[i - 1] = ap(i);
                    setm(75, 77);
                    llrz[74] = apmag * (-1.0f);
                    llrz[75] = apmag * (-1.0f);
                    llrz[76] = apmag * (+1.0f);
                } else if (ncontest == 3) {  // Field Day
                    setm(1, 56);
                    for (int i = 1; i <= 28; ++i) llrz[i - 1] = ap(i);
                    for (int i = 29; i <= 56; ++i) llrz[i - 1] = ap(i + 1);
                    setm(72, 74);
                    setm(75, 77);
                    for (int i = 75; i <= 77; ++i) llrz[i - 1] = apmag * (-1.0f);
                } else if (ncontest == 4) {
                    setm(2, 57);
                    for (int i = 2; i <= 29; ++i) llrz[i - 1] = ap(i - 1);
                    for (int i = 30; i <= 57; ++i) llrz[i - 1] = ap(i);
                    setm(75, 77);
                    llrz[74] = apmag * (-1.0f);
                    llrz[75] = apmag * (+1.0f);
                    llrz[76] = apmag * (+1.0f);
                }
            }

            if (iaptype == 5 && ncontest == 7) continue;  // Hound
            if (iaptype == 4 || iaptype == 5 || iaptype == 6) {
                std::memset(apmask, 0, sizeof apmask);
                if (ncontest <= 5 || (ncontest == 7 && iaptype == 6) || ncontest == 8) {
                    setm(1, 77);  // mycall, hiscall, RRR|73|RR73
                    for (int i = 1; i <= 58; ++i) llrz[i - 1] = ap(i);
                    const int* m = iaptype == 4 ? mrrr : iaptype == 5 ? m73 : mrr73;
                    for (int i = 59; i <= 77; ++i) llrz[i - 1] = apmag * pm(m[i - 59]);
                } else if (ncontest == 7 && iaptype == 4) {  // Hound listens for MyCall RR73;...
                    setm(1, 28);
                    for (int i = 1; i <= 28; ++i) llrz[i - 1] = ap(i);
                    setm(57, 66);
                    for (int i = 57; i <= 66; ++i) llrz[i - 1] = apmag * static_cast<float>(aph10[i - 57]);
                    setm(72, 77);
                    llrz[71] = apmag * (-1.0f);
                    llrz[72] = apmag * (-1.0f);
                    llrz[73] = apmag * (1.0f);
                    for (int i = 75; i <= 77; ++i) llrz[i - 1] = apmag * (-1.0f);
                }
            }
        }

        std::memset(cw, 0, sizeof cw);
        dmin = 0.0f;
        const int norder = 2;
        int maxosd = 2;
        if (ndepth == 1) maxosd = -1;  // BP only
        if (ndepth == 3 && (std::fabs(static_cast<float>(nfqso) - f1) <= static_cast<float>(napwid) ||
                            std::fabs(static_cast<float>(nftx) - f1) <= static_cast<float>(napwid) || ncontest == 7))
            maxosd = 2;
        const int Keff = 91;
        int ntype = 0;
        decode174_91(llrz, Keff, maxosd, norder, apmask, message91, cw, ntype, nharderrors, dmin);
        if (nharderrors >= 0) std::memcpy(message77, message91, 77);
        if (trace_ && (ipass > 5 || nharderrors >= 0)) {
            std::string t = "T " + I(ipass) + " " + I(ntype) + " " + I(nharderrors) + " " + XF(dmin);
            if (ipass > 5) t += " " + I(iaptype) + " " + XF(apmag) + " " + X(h_r(llrz, 174));
            if (nharderrors >= 0) t += " " + X(h_i1(cw, 174));
            trace_(t + " " + I(maxosd) + " " + X(h_i1(apmask, 174)));
        }

        msg37.assign(37, ' ');
        nbadcrc = 1;
        if (nharderrors < 0 || nharderrors > 36) continue;
        if (std::count(cw, cw + 174, 0) == 174) continue;  // Reject the all-zero codeword
        char c77[78];
        for (int i = 0; i < 77; ++i) c77[i] = static_cast<char>('0' + message77[i]);
        c77[77] = 0;
        const int n3 = 4 * message77[71] + 2 * message77[72] + message77[73];
        const int i3 = 4 * message77[74] + 2 * message77[75] + message77[76];
        if (i3 > 5 || (i3 == 0 && n3 > 6)) continue;
        if (i3 == 0 && n3 == 2) continue;
        const bool unpk77_success = codec.unpack77(c77, 1, msg37);
        if (trace_)
            trace_("Q " + I(ipass) + " " + L(unpk77_success) + " " + I(i3) + " " + I(n3) + " |" + msg37 + "|");
        if (!unpk77_success || msg37.find("/R") != std::string::npos || msg37.compare(0, 4, "TU; ") == 0) {
            if (i3 >= 1 && i3 <= 3 && ncontest == 0) continue;
        }
        if (!unpk77_success) continue;
        // If we get this far: valid codeword, valid (i3,n3), nonquirky message.
        nbadcrc = 0;
        get_ft8_tones_from_77bits(message77, itone);
        if (trace_) {
            std::string o = "O " + I(ipass) + " ";
            for (int i = 0; i < 77; ++i) o += static_cast<char>('0' + message77[i]);
            o += ' ';
            for (int i = 0; i < NN; ++i) o += static_cast<char>('0' + itone[i]);
            trace_(o);
        }
        if (lsubtract) subtractft8(dd0, itone, f1, xdt, false);
        float xsig = 0.0f, xnoi = 0.0f;
        for (int i = 1; i <= 79; ++i) {
            xsig = xsig + s8_[i - 1][itone[i - 1]] * s8_[i - 1][itone[i - 1]];
            const int ios = (itone[i - 1] + 4) % 7;
            xnoi = xnoi + s8_[i - 1][ios] * s8_[i - 1][ios];
        }
        xsnr = 0.001f;
        float xsnr2 = 0.001f;
        float arg = xsig / xnoi - 1.0f;
        if (arg > 0.1f) xsnr = arg;
        arg = xsig / xbase / 3.0e6f - 1.0f;
        if (arg > 0.1f) xsnr2 = arg;
        xsnr = 10.0f * log10f_glibc(xsnr) - 27.0f;
        xsnr2 = 10.0f * log10f_glibc(xsnr2) - 27.0f;
        if (trace_) trace_("Z " + XF(xsig) + " " + XF(xnoi) + " " + XF(xsnr) + " " + XF(xsnr2));
        if (!nagain) xsnr = xsnr2;
        if (nsync <= 10 && xsnr < -25.0f) {  // bail out, likely false decode
            nbadcrc = 1;
            if (trace_) trace_ret(nbadcrc, ipass, iaptype, nharderrors, f1, xdt, dmin, &xsnr, msg37);
            return;
        }
        if (xsnr < -25.0f) xsnr = -25.0f;
        if (trace_) trace_ret(nbadcrc, ipass, iaptype, nharderrors, f1, xdt, dmin, &xsnr, msg37);
        return;
    }
    if (trace_) trace_ret(nbadcrc, ipass, iaptype, nharderrors, f1, xdt, dmin, nullptr, msg37);
}

void Frontend::trace_ret(int nbadcrc, int ipass, int iaptype, int nharderrors, float f1, float xdt, float dmin,
                         const float* xsnr, const std::string& msg37) {
    trace_("r " + I(nbadcrc) + " " + I(ipass) + " " + I(iaptype) + " " + I(nharderrors) + " " + XF(f1) + " " +
           XF(xdt) + " " + XF(dmin) + " " + (xsnr ? XF(*xsnr) : std::string("-")) + " |" + msg37 + "|");
}

// gen_ft8wave's SAVEd tables, (re)made when fchk changes as there. ctab
// depends on twopi alone (always TWOPI), so it is made once: the same
// entries the reference makes again each time.
void Frontend::gw_tables(int nsym, int nsps, float bt, float fsample) {
    constexpr int NTAB = 65536;
    const float fchk = (static_cast<float>(nsym + nsps) + bt) + fsample;
    if (fchk != fchk0_) {
        gw_twopi_ = TWOPI;
        gw_dt_ = 1.0f / fsample;
        gw_hmod_ = 1.0f;
        // Compute the frequency-smoothing pulse
        ++gw_npulse_;
        for (int i = 1; i <= 3 * nsps; ++i) {
            const float tt = (static_cast<float>(i) - 1.5f * static_cast<float>(nsps)) / static_cast<float>(nsps);
            pulse_[i - 1] = gfsk_pulse(bt, tt);
        }
        if (!ctab_made_) {
            // ctab: entries 1..65532 by libmvec cosf and sinf, four at a time;
            // 0 and 65533..65535 by sincosf, as the reference's -O3 build.
            const float scale = 1.0f / 65536.0f;
            ctab_[0] = cis((0.0f * gw_twopi_) * scale);
            for (int i = 1; i + 3 <= 65532; i += 4) {
                float a[4], c[4], s[4];
                for (int q = 0; q < 4; ++q) a[q] = (static_cast<float>(i + q) * gw_twopi_) * scale;
                cosf4_glibc(a, c);
                sinf4_glibc(a, s);
                for (int q = 0; q < 4; ++q) ctab_[i + q] = cfloat(c[q], s[q]);
            }
            for (int i = 65533; i < NTAB; ++i) ctab_[i] = cis((static_cast<float>(i) * gw_twopi_) * scale);
            ctab_made_ = true;
        }
        fchk0_ = fchk;
    }
}

// subtractft8's filter and end correction, made on the first call.
void Frontend::st_make_filter() {
    constexpr int NFFT = NMAX, NFILT = 4000;
    if (st_first_) {  // Create and normalize the filter
        const float fac = 1.0f / static_cast<float>(NFFT);
        std::vector<float> window(NFILT + 1);  // window(-NFILT/2:NFILT/2) at [j+2000]
        // window(-2000) folded by gcc; -1999..2000 by libmvec, four at a time.
        window[0] = bitsf(0x2709ADEDu);
        float sumw = 0.0f;
        sumw = sumw + window[0];
        for (int j = -NFILT / 2 + 1; j <= NFILT / 2; j += 4) {
            float a[4], c[4];
            for (int q = 0; q < 4; ++q) a[q] = (PI * static_cast<float>(j + q)) / static_cast<float>(NFILT);
            cosf4_glibc(a, c);
            for (int q = 0; q < 4; ++q) {
                window[j + q + NFILT / 2] = c[q] * c[q];
                sumw = sumw + window[j + q + NFILT / 2];
            }
        }
        cfloat* cw = cw_.data();
        std::vector<cfloat> tmp(NMAX, cfloat(0.0f, 0.0f));
        for (int i = 1; i <= NFILT + 1; ++i) tmp[i - 1] = cfloat(window[i - 1] / sumw, 0.0f);
        for (int m = 0; m < NMAX; ++m) cw[m] = tmp[(m + NFILT / 2 + 1) % NMAX];  // cshift(cw, NFILT/2+1)
        fst4w::fft_c2c_inplace(cw, NFFT, -1);
        for (int m = 0; m < NMAX; ++m) cw[m] = cscale(cw[m], fac);
        st_first_ = false;
        for (int j = 1; j <= NFILT / 2 + 1; ++j) {
            float s = 0.0f;
            for (int q = j - 1; q <= NFILT / 2; ++q) s += window[q + NFILT / 2];
            endcorrection_[j - 1] = 1.0f / (1.0f - s / sumw);
        }
    }
}

// four2a's plan for sqf's x (as fst4w's: FFTW_ESTIMATE_PATIENT, in place,
// the decoders' wisdom, fft.h), made once on the first trial.
void Frontend::st_plan_r2c(float* x) {
    constexpr int NFFT = NMAX;
    if (!st_plan_) {
        std::lock_guard<std::mutex> lock(skimmer::fft_planner_lock());
        char* const saved = fftwf_export_wisdom_to_string();
        fftwf_forget_wisdom();
        skimmer::decoder_wisdom_before_plan();
        std::vector<float> keep(x, x + NFFT + 2);
        st_plan_ = fftwf_plan_dft_r2c_1d(NFFT, x, reinterpret_cast<fftwf_complex*>(x), FFTW_ESTIMATE_PATIENT);
        std::copy(keep.begin(), keep.end(), x);
        skimmer::decoder_wisdom_after_plan();
        fftwf_forget_wisdom();
        if (saved) {
            fftwf_import_wisdom_from_string(saved);
            std::free(saved);
        }
        if (!st_plan_) throw std::runtime_error("ft8: FFTW could not plan");
    }
}

// gen_ft8wave as gen_ft8wave_reference, faster with every value the same.
// dphi_peak*pulse(m) is kept between calls (the same product each time),
// and a symbol whose tone is 0 adds nothing to dphi (its terms are +0, and
// dphi is never -0). The phase recurrence's fmod becomes a subtraction:
// when every dphi it adds is in [0, twopi) phi stays in [0, twopi), each
// sum phi + dphi is below 2*twopi (2*pred(twopi) is a float, so rounding
// cannot reach it), and for twopi <= s < 2*twopi fmod(s, twopi) is s - twopi,
// which Sterbenz makes exact; otherwise fmod as the reference. The table
// index is then computed from the stored phases in a separate pass.
void Frontend::gen_ft8wave(const int* itone, int nsym, int nsps, float bt, float fsample, float f0, cfloat* cwave,
                           float* wave, int icmplx, int nwave) {
    constexpr int NTAB = 65536;
    if (icmplx == 0) {
        gen_ft8wave_reference(itone, nsym, nsps, bt, fsample, f0, cwave, wave, icmplx, nwave);
        return;
    }
    gw_tables(nsym, nsps, bt, fsample);

    const float dphi_peak = (gw_twopi_ * gw_hmod_) / static_cast<float>(nsps);
    if (pdphi_npulse_ != gw_npulse_ || pdphi_nsps_ != nsps) {
        pdphi_.resize(static_cast<size_t>(3 * nsps));
        for (int m = 0; m < 3 * nsps; ++m) pdphi_[m] = dphi_peak * pulse_[m];
        pdphi_npulse_ = gw_npulse_;
        pdphi_nsps_ = nsps;
    }
    const int ndphi = (nsym + 2) * nsps;
    dphi_.assign(static_cast<size_t>(ndphi), 0.0f);
    float* dphi = dphi_.data();
    const float* pd = pdphi_.data();
    for (int j = 1; j <= nsym; ++j) {
        if (itone[j - 1] == 0) continue;
        const int ib = (j - 1) * nsps;
        const float t = static_cast<float>(itone[j - 1]);
        for (int m = 0; m < 3 * nsps; ++m) dphi[ib + m] = dphi[ib + m] + pd[m] * t;
    }
    {
        const float* pulse = pulse_.data();
        const float a = dphi_peak * static_cast<float>(itone[0]);
        for (int m = 0; m < 2 * nsps; ++m) dphi[m] = dphi[m] + a * pulse[nsps + m];
        const float b = dphi_peak * static_cast<float>(itone[nsym - 1]);
        for (int m = 0; m < 2 * nsps; ++m) dphi[nsym * nsps + m] = dphi[nsym * nsps + m] + b * pulse[m];
    }

    const float twopi = gw_twopi_;
    const float shift = (twopi * f0) * gw_dt_;
    for (int m = 0; m < ndphi; ++m) dphi[m] = dphi[m] + shift;
    unsigned out = 0;
    for (int m = 0; m < ndphi; ++m) out |= static_cast<unsigned>(!(dphi[m] >= 0.0f && dphi[m] < twopi));  // NaN too
    const bool in_range = out == 0;
    if (gw_phi_.size() < static_cast<size_t>(nwave)) gw_phi_.resize(static_cast<size_t>(nwave));
    float* ph = gw_phi_.data();
    float phi = 0.0f;
    const float* dp = dphi + nsps;  // Don't include dummy symbols
    if (in_range) {
        for (int k = 0; k < nwave; ++k) {
            ph[k] = phi;
            const float s = phi + dp[k];
            phi = (s >= twopi) ? s - twopi : s;
        }
    } else {
        for (int k = 0; k < nwave; ++k) {
            ph[k] = phi;
            phi = std::fmod(phi + dp[k], twopi);
        }
    }
    // The index (as the reference's: IEEE divide, truncated), then the
    // table's entries copied whole.
    if (gw_ix_.size() < static_cast<size_t>(nwave)) gw_ix_.resize(static_cast<size_t>(nwave));
    int* ix = gw_ix_.data();
    for (int k = 0; k < nwave; ++k) ix[k] = static_cast<int>((ph[k] * static_cast<float>(NTAB)) / twopi);
    const cfloat* ctab = ctab_.data();
    for (int k = 0; k < nwave; ++k) std::memcpy(&cwave[k], &ctab[ix[k]], sizeof(cfloat));

    // Apply envelope shaping to the first and last symbols, as the reference.
    const int nramp = nint(static_cast<float>(nsps) / 8.0f);
    if (gw_c_.size() < static_cast<size_t>(nramp)) gw_c_.resize(static_cast<size_t>(nramp));
    float* c = gw_c_.data();
    {
        const float den = 2.0f * static_cast<float>(nramp);
        int i = 0;
        for (; i + 3 < nramp; i += 4) {
            float a[4];
            for (int q = 0; q < 4; ++q) a[q] = (twopi * static_cast<float>(i + q)) / den;
            cosf4_glibc(a, &c[i]);
        }
        for (; i < nramp; ++i) c[i] = cosf_glibc((twopi * static_cast<float>(i)) / den);
    }
    const int k1 = nsym * nsps - nramp + 1;
    for (int i = 0; i < nramp; ++i) {
        const float up = (1.0f - c[i]) / 2.0f, down = (1.0f + c[i]) / 2.0f;
        cwave[i] = cscale(cwave[i], up);
        cwave[k1 - 1 + i] = cscale(cwave[k1 - 1 + i], down);
    }
}

// sqf as sqf_reference with the same values, less memory traffic: cfilt is
// written directly (camp's zeros only outside the frame's valid j), dd0 is
// read in place of its copy; a trial (ldt) writes only x, the subtraction
// (not ldt) writes dd0 in place and no x (only a trial's r2c reads x; it
// returns 0 as the reference does). x stays where the reference puts it,
// on the unaligned plan's buffer.
float Frontend::sqf(int idt, float* dd0, float f0, float dt, bool ldt) {
    if (fast_subtract_) return sqf_fast(idt, dd0, f0, dt, ldt);
    constexpr int NFRAME = 1920 * 79, NFFT = NMAX, NFILT = 4000;
    const int nstart = static_cast<int>((dt * 12000.0f + 1.0f) + static_cast<float>(idt));
    cfloat* cfilt = cfilt_.data();
    const cfloat* cref = cref_.data();
    const cfloat* cw = cw_.data();
    // i in [ilo, ihi] has j = nstart - 1 + i in [1, NMAX]
    const int ilo = std::max(1, 2 - nstart), ihi = std::min(NFRAME, NMAX - nstart + 1);
    const int n0 = std::min(std::max(ilo - 1, 0), NFRAME), n1 = std::max(ihi, n0);
    std::fill(cfilt, cfilt + n0, cfloat(0.0f, 0.0f));
    {
        const float* d = dd0 + (nstart - 1);  // d[i-1] = dd0[j-1]
        for (int i = n0; i < n1; ++i) cfilt[i] = cfloat(d[i] * cref[i].real(), d[i] * -cref[i].imag());
    }
    std::fill(cfilt + n1, cfilt + NMAX, cfloat(0.0f, 0.0f));
    fst4w::fft_c2c_inplace(cfilt, NFFT, -1);
    for (int i = 0; i < NFFT; ++i) cfilt[i] = cmul(cfilt[i], cw[i]);
    fst4w::fft_c2c_inplace(cfilt, NFFT, 1);
    for (int m = 1; m <= NFILT / 2 + 1; ++m) cfilt[m - 1] = cscale(cfilt[m - 1], endcorrection_[m - 1]);
    for (int m = 1; m <= NFILT / 2 + 1; ++m)
        cfilt[NFRAME - m] = cscale(cfilt[NFRAME - m], endcorrection_[m - 1]);
    if (!ldt) {
        float* d = dd0 + (nstart - 1);
        for (int i = n0; i < n1; ++i) {
            const float zr = cfilt[i].real() * cref[i].real() - cfilt[i].imag() * cref[i].imag();
            d[i] = d[i] - 2.0f * zr;  // Subtract the reconstructed signal
        }
        return 0.0f;
    }
    float* x = st_x_.data() + 1;
    {
        const float* d = dd0 + (nstart - 1);
        std::fill(x, x + n0, 0.0f);
        for (int i = n0; i < n1; ++i) {
            const float zr = cfilt[i].real() * cref[i].real() - cfilt[i].imag() * cref[i].imag();
            x[i] = d[i] - 2.0f * zr;
        }
        std::fill(x + n1, x + NFFT + 2, 0.0f);
    }
    st_plan_r2c(x);
    fftwf_execute(static_cast<fftwf_plan>(st_plan_));  // Forward FFT, r2c
    float sqq = 0.0f;
    const float df = 12000.0f / static_cast<float>(NFFT);
    const int ia = static_cast<int>((f0 - 1.5f * 6.25f) / df);
    const int ib = static_cast<int>((f0 + 8.5f * 6.25f) / df);
    for (int i = ia; i <= ib; ++i) sqq = (sqq + x[2 * i] * x[2 * i]) + x[2 * i + 1] * x[2 * i + 1];
    return sqq;
}

// Fast subtraction (fast_subtract_), sqf's job without its 180,000-point FFTs. sqf smooths c = dd0*conj(cref) over the
// frame with the 4001-tap normalized cos^2 window w: y[i] = sum_j w[j] c[i+1-j], j = -2000..2000 (the +1 is cshift's;
// c zero outside the frame and outside dd0), end-corrected over its first and last 2001 samples, then subtracts
// 2*Re(y*cref). y's band is a few Hz wide, so here c is taken down to 300 S/s through a 2D-sample triangle (two
// cascaded 40-sample block means: their aliases ~80 dB down), convolved with w sampled every 40 samples, and taken back
// to each sample linearly between the block centres; the subtraction is done in the same loop. A trial (ldt) writes
// the residual to x and measures its band power with sqf's own r2c FFT: measuring it more cheaply (block sums of the
// residual mixed down to the band) moved the refined DT enough to change ~5% of the decodes on crowded bands.
float Frontend::sqf_fast(int idt, float* dd0, float f0, float dt, bool ldt) {
    constexpr int NFRAME = 1920 * 79, NFILT = 4000, D = 40, M = NFRAME / D, Q = NFILT / 2 / D + 1;
    if (fs_kernel_.empty()) {
        std::vector<double> w(NFILT + 1);
        double sumw = 0;
        for (int j = -NFILT / 2; j <= NFILT / 2; ++j) {
            const double a = std::cos(3.14159265358979323846 * j / NFILT);
            w[j + NFILT / 2] = a * a;
            sumw += a * a;
        }
        std::vector<double> k(2 * Q + 1, 0.0);  // K[q]: the window at lag q*D, D of them per decimated sample
        for (int q = -Q; q <= Q; ++q) {
            const int j = q * D;
            if (j >= -NFILT / 2 && j <= NFILT / 2) k[q + Q] = D * w[j + NFILT / 2] / sumw;
        }
        fs_kernel_.assign(k.begin(), k.end());
        fs_b_.resize(M);
        fs_y_.resize(M);
    }
    const int nstart = static_cast<int>((dt * 12000.0f + 1.0f) + static_cast<float>(idt));
    const int ilo = std::max(1, 2 - nstart), ihi = std::min(NFRAME, NMAX - nstart + 1);
    const int n0 = std::min(std::max(ilo - 1, 0), NFRAME), n1 = std::max(ihi, n0);
    float* d = dd0 + (nstart - 1);
    const cfloat* cref = cref_.data();
    // c decimated to 300 S/s through a 2D-sample triangle (two cascaded block means: its aliases are down ~80 dB where
    // a single block mean's were ~40 dB, which mattered on crowded bands), centred on the block centres t_m
    for (int m = 0; m < M; ++m) {
        const float tc = m * D + (D - 1) / 2.0f;  // the triangle's centre
        const int i0 = std::max(m * D - D / 2, n0), i1 = std::min(m * D + D + D / 2, n1);
        float sr = 0, si = 0;
        for (int i = i0; i < i1; ++i) {
            const float g = 1.0f - std::fabs(static_cast<float>(i) - tc) / D;
            if (g <= 0) continue;
            const float v = g * d[i];
            sr += v * cref[i].real();
            si -= v * cref[i].imag();
        }
        fs_b_[m] = cfloat(sr / D, si / D);
    }
    // y at the block centres t_m = m*D + (D-1)/2 - 1
    const float* kq = fs_kernel_.data();
    for (int m = 0; m < M; ++m) {
        const int qlo = std::max(-Q, m - (M - 1)), qhi = std::min(Q, m);
        float yr = 0, yi = 0;
        for (int q = qlo; q <= qhi; ++q) {
            yr += kq[q + Q] * fs_b_[m - q].real();
            yi += kq[q + Q] * fs_b_[m - q].imag();
        }
        fs_y_[m] = cfloat(yr, yi);
    }
    // Each sample: y linear between centres, end-corrected, and the residual d - 2*Re(y*cref): subtracted in dd0, or,
    // for a trial, written to sqf's x for its r2c FFT and band power, as sqf measures it
    float* x = ldt ? st_x_.data() + 1 : nullptr;
    if (ldt) std::fill(x, x + n0, 0.0f);
    for (int i = n0; i < n1; ++i) {
        const int u = 2 * i - (D - 3);  // 2*(i - t_0), t_0 = (D-1)/2 - 1
        int m = u >= 0 ? u / (2 * D) : -((-u + 2 * D - 1) / (2 * D));
        float f = static_cast<float>(u - 2 * D * m) / (2 * D);
        if (m < 0) m = 0, f = 0;
        if (m >= M - 1) m = M - 2, f = 1;
        float yr = fs_y_[m].real() + f * (fs_y_[m + 1].real() - fs_y_[m].real());
        float yi = fs_y_[m].imag() + f * (fs_y_[m + 1].imag() - fs_y_[m].imag());
        float ec = 1.0f;
        if (i <= NFILT / 2) ec = endcorrection_[i];
        else if (i >= NFRAME - NFILT / 2 - 1) ec = endcorrection_[NFRAME - 1 - i];
        const float r = d[i] - 2.0f * ((yr * ec) * cref[i].real() - (yi * ec) * cref[i].imag());
        if (ldt) x[i] = r;
        else d[i] = r;
    }
    if (!ldt) return 0.0f;
    std::fill(x + n1, x + NMAX + 2, 0.0f);
    st_plan_r2c(x);
    fftwf_execute(static_cast<fftwf_plan>(st_plan_));  // Forward FFT, r2c
    float sqq = 0.0f;
    const float df = 12000.0f / static_cast<float>(NMAX);
    const int ia = static_cast<int>((f0 - 1.5f * 6.25f) / df);
    const int ib = static_cast<int>((f0 + 8.5f * 6.25f) / df);
    for (int i = ia; i <= ib; ++i) sqq = (sqq + x[2 * i] * x[2 * i]) + x[2 * i + 1] * x[2 * i + 1];
    return sqq;
}

// subtractft8 as subtractft8_reference: sqf subtracts in dd0 itself, so no
// copy back.
void Frontend::subtractft8(float* dd0, const int itone[NN], float f0, float dt, bool lrefinedt) {
    constexpr int NFRAME = 1920 * 79;
    if (trace_) {
        trace_flush();
        std::string t = "U " + XF(f0) + " " + XF(dt) + " " + L(lrefinedt) + " ";
        for (int i = 0; i < NN; ++i) t += static_cast<char>('0' + itone[i]);
        trace_(t + " " + X(h_r(dd0, NMAX)));
    }
    st_make_filter();

    // Generate complex reference waveform cref
    gen_ft8wave(itone, 79, 1920, 2.0f, 12000.0f, f0, cref_.data(), nullptr, 1, NFRAME);
    if (trace_) trace_("V " + X(h_c(cref_.data(), NFRAME)));

    if (lrefinedt) {  // Are we refining DT ?
        const float sqa = sqf(-90, dd0, f0, dt, true);
        const float sqb = sqf(+90, dd0, f0, dt, true);
        const float sq0 = sqf(0, dd0, f0, dt, true);
        const float dx = peakup(sqa, sq0, sqb);
        if (trace_) trace_("R " + XF(sqa) + " " + XF(sqb) + " " + XF(sq0) + " " + XF(dx));
        if (std::fabs(dx) > 1.0f) return;  // No acceptable minimum: do not subtract
        const int i2 = nint(90.0f * dx);   // Best estimate of idt
        sqf(i2, dd0, f0, dt, false);       // Do the subtraction with idt=i2
    } else {
        sqf(0, dd0, f0, dt, false);  // Do the subtraction with idt=0
    }
    if (trace_) trace_("u " + X(h_r(dd0, NMAX)));
}

void Frontend::gen_ft8wave_reference(const int* itone, int nsym, int nsps, float bt, float fsample, float f0,
                                     cfloat* cwave, float* wave, int icmplx, int nwave) {
    constexpr int NTAB = 65536;
    gw_tables(nsym, nsps, bt, fsample);

    // Compute the smoothed frequency waveform.
    // Length = (nsym+2)*nsps samples, first and last symbols extended
    const float dphi_peak = (gw_twopi_ * gw_hmod_) / static_cast<float>(nsps);
    const int ndphi = (nsym + 2) * nsps;
    dphi_.assign(static_cast<size_t>(ndphi), 0.0f);
    float* dphi = dphi_.data();
    const float* pulse = pulse_.data();
    for (int j = 1; j <= nsym; ++j) {
        const int ib = (j - 1) * nsps;
        const float t = static_cast<float>(itone[j - 1]);
        for (int m = 0; m < 3 * nsps; ++m) dphi[ib + m] = dphi[ib + m] + (dphi_peak * pulse[m]) * t;
    }
    // Add dummy symbols at beginning and end with tone values equal to 1st
    // and last symbol, respectively
    {
        const float a = dphi_peak * static_cast<float>(itone[0]);
        for (int m = 0; m < 2 * nsps; ++m) dphi[m] = dphi[m] + a * pulse[nsps + m];
        const float b = dphi_peak * static_cast<float>(itone[nsym - 1]);
        for (int m = 0; m < 2 * nsps; ++m) dphi[nsym * nsps + m] = dphi[nsym * nsps + m] + b * pulse[m];
    }

    // Calculate and insert the audio waveform
    float phi = 0.0f;
    const float shift = (gw_twopi_ * f0) * gw_dt_;  // Shift frequency up by f0
    for (int m = 0; m < ndphi; ++m) dphi[m] = dphi[m] + shift;
    if (icmplx == 0) std::fill(wave, wave + nwave, 0.0f);
    if (icmplx != 0) std::fill(cwave, cwave + nwave, cfloat(0.0f, 0.0f));

    int k = 0;
    for (int j = nsps; j <= nsps + nwave - 1; ++j) {  // Don't include dummy symbols
        ++k;
        if (icmplx == 0) {
            wave[k - 1] = sinf_glibc(phi);
        } else {
            const int i = static_cast<int>((phi * static_cast<float>(NTAB)) / gw_twopi_);
            cwave[k - 1] = ctab_[i];
        }
        phi = std::fmod(phi + dphi[j], gw_twopi_);
    }

    // Apply envelope shaping to the first and last symbols: cos by libmvec
    // four at a time from i = 0, the rest scalar.
    const int nramp = nint(static_cast<float>(nsps) / 8.0f);
    std::vector<float> c(static_cast<size_t>(nramp));
    {
        const float den = 2.0f * static_cast<float>(nramp);
        int i = 0;
        for (; i + 3 < nramp; i += 4) {
            float a[4];
            for (int q = 0; q < 4; ++q) a[q] = (gw_twopi_ * static_cast<float>(i + q)) / den;
            cosf4_glibc(a, &c[i]);
        }
        for (; i < nramp; ++i) c[i] = cosf_glibc((gw_twopi_ * static_cast<float>(i)) / den);
    }
    const int k1 = nsym * nsps - nramp + 1;
    for (int i = 0; i < nramp; ++i) {
        const float up = (1.0f - c[i]) / 2.0f, down = (1.0f + c[i]) / 2.0f;
        if (icmplx == 0) {
            wave[i] = wave[i] * up;
            wave[k1 - 1 + i] = wave[k1 - 1 + i] * down;
        } else {
            cwave[i] = cscale(cwave[i], up);
            cwave[k1 - 1 + i] = cscale(cwave[k1 - 1 + i], down);
        }
    }
}

float Frontend::sqf_reference(int idt, float* dd0, float f0, float dt, bool ldt) {
    constexpr int NFRAME = 1920 * 79, NFFT = NMAX, NFILT = 4000;
    const int nstart = static_cast<int>((dt * 12000.0f + 1.0f) + static_cast<float>(idt));
    float* dd = st_dd_.data();
    cfloat* camp = camp_.data();
    cfloat* cfilt = cfilt_.data();
    const cfloat* cref = cref_.data();
    const cfloat* cw = cw_.data();
    std::fill(camp, camp + NMAX, cfloat(0.0f, 0.0f));
    std::memcpy(dd, dd0, sizeof(float) * NMAX);
    for (int i = 1; i <= NFRAME; ++i) {
        const int j = nstart - 1 + i;
        if (j >= 1 && j <= NMAX) camp[i - 1] = cfloat(dd[j - 1] * cref[i - 1].real(), dd[j - 1] * -cref[i - 1].imag());
    }
    std::memcpy(cfilt, camp, sizeof(cfloat) * NFRAME);
    std::fill(cfilt + NFRAME, cfilt + NMAX, cfloat(0.0f, 0.0f));
    fst4w::fft_c2c_inplace(cfilt, NFFT, -1);
    for (int i = 0; i < NFFT; ++i) cfilt[i] = cmul(cfilt[i], cw[i]);
    fst4w::fft_c2c_inplace(cfilt, NFFT, 1);
    for (int m = 1; m <= NFILT / 2 + 1; ++m) cfilt[m - 1] = cscale(cfilt[m - 1], endcorrection_[m - 1]);
    for (int m = 1; m <= NFILT / 2 + 1; ++m)
        cfilt[NFRAME - m] = cscale(cfilt[NFRAME - m], endcorrection_[m - 1]);
    // x as sqf has it: in the internal function's frame, 4 bytes past a
    // 16-byte boundary (its r2c's plan is FFTW's for unaligned data, whose
    // bits differ from the aligned plan's).
    float* x = st_x_.data() + 1;
    std::fill(x, x + NFFT + 2, 0.0f);
    for (int i = 1; i <= NFRAME; ++i) {
        const int j = nstart + i - 1;
        if (j >= 1 && j <= NMAX) {
            const float zr = cfilt[i - 1].real() * cref[i - 1].real() - cfilt[i - 1].imag() * cref[i - 1].imag();
            dd[j - 1] = dd[j - 1] - 2.0f * zr;  // Subtract the reconstructed signal
            x[i - 1] = dd[j - 1];
        }
    }
    float sqq = 0.0f;
    if (ldt) {
        st_plan_r2c(x);
        fftwf_execute(static_cast<fftwf_plan>(st_plan_));  // Forward FFT, r2c
        const float df = 12000.0f / static_cast<float>(NFFT);
        const int ia = static_cast<int>((f0 - 1.5f * 6.25f) / df);
        const int ib = static_cast<int>((f0 + 8.5f * 6.25f) / df);
        for (int i = ia; i <= ib; ++i) sqq = (sqq + x[2 * i] * x[2 * i]) + x[2 * i + 1] * x[2 * i + 1];
    }
    return sqq;
}

void Frontend::subtractft8_reference(float* dd0, const int itone[NN], float f0, float dt, bool lrefinedt) {
    constexpr int NFRAME = 1920 * 79;
    if (trace_) {
        trace_flush();
        std::string t = "U " + XF(f0) + " " + XF(dt) + " " + L(lrefinedt) + " ";
        for (int i = 0; i < NN; ++i) t += static_cast<char>('0' + itone[i]);
        trace_(t + " " + X(h_r(dd0, NMAX)));
    }
    st_make_filter();

    // Generate complex reference waveform cref
    gen_ft8wave_reference(itone, 79, 1920, 2.0f, 12000.0f, f0, cref_.data(), nullptr, 1, NFRAME);
    if (trace_) trace_("V " + X(h_c(cref_.data(), NFRAME)));

    bool ldt = lrefinedt;
    float sq0;
    if (ldt) {  // Are we refining DT ?
        const float sqa = sqf_reference(-90, dd0, f0, dt, ldt);
        const float sqb = sqf_reference(+90, dd0, f0, dt, ldt);
        sq0 = sqf_reference(0, dd0, f0, dt, ldt);
        const float dx = peakup(sqa, sq0, sqb);
        if (trace_) trace_("R " + XF(sqa) + " " + XF(sqb) + " " + XF(sq0) + " " + XF(dx));
        if (std::fabs(dx) > 1.0f) return;  // No acceptable minimum: do not subtract
        const int i2 = nint(90.0f * dx);   // Best estimate of idt
        ldt = false;
        sq0 = sqf_reference(i2, dd0, f0, dt, ldt);  // Do the subtraction with idt=i2
    } else {
        sq0 = sqf_reference(0, dd0, f0, dt, ldt);  // Do the subtraction with idt=0
    }
    (void)sq0;
    std::memcpy(dd0, st_dd_.data(), sizeof(float) * NMAX);  // Return dd0 with this signal subtracted
    if (trace_) trace_("u " + X(h_r(dd0, NMAX)));
}

void Frontend::ft8_a7d(const float* dd0, bool& newdat, const std::string& call_1in, const std::string& call_2in,
                       const std::string& grid4in, float& xdt, float& f1, float xbase, int& nharderrors, float& dmin,
                       std::string& msg37, float& xsnr, Codec& codec) {
    auto fix = [](std::string s, size_t n) {
        s.resize(n, ' ');
        return s;
    };
    auto trim = [](const std::string& s) {
        size_t n = s.size();
        while (n > 0 && s[n - 1] == ' ') --n;
        return s.substr(0, n);
    };
    const std::string call_1 = fix(call_1in, 12), call_2 = fix(call_2in, 12), grid4 = fix(grid4in, 4);

    bool std_1 = codec.stdcall(call_1);
    if (call_1.compare(0, 3, "CQ ") == 0) std_1 = true;
    const bool std_2 = codec.stdcall(call_2);
    if (trace_) {
        trace_flush();
        trace_("a " + XF(f1) + " " + XF(xdt) + " " + XF(xbase) + " " + L(newdat) + " " + L(std_1) + " " + L(std_2) +
               " |" + call_1 + "|" + call_2 + "|" + grid4 + "|");
    }

    Sync sy;
    fine_sync(dd0, newdat, f1, xdt, sy, true);
    if (trace_) trace_sync(sy, xdt, false, 0, 0, 0, false);
    // (no bail-out on nsync here)

    float bmeta[174], bmetb[174], bmetc[174], bmetd[174];
    bit_metrics(1, bmeta, bmetb, bmetc, bmetd, nullptr);
    float* const bm[4] = {bmeta, bmetb, bmetc, bmetd};
    if (trace_) trace_metrics(4, bm, false);
    normalizebmet(bmeta, 174);
    normalizebmet(bmetb, 174);
    normalizebmet(bmetc, 174);
    normalizebmet(bmetd, 174);
    const float scalefac = 2.83f;
    float llra[174], llrb[174], llrc[174], llrd[174];
    for (int i = 0; i < 174; ++i) {
        llra[i] = scalefac * bmeta[i];
        llrb[i] = scalefac * bmetb[i];
        llrc[i] = scalefac * bmetc[i];
        llrd[i] = scalefac * bmetd[i];
    }
    float* const ll[4] = {llra, llrb, llrc, llrd};
    if (trace_) trace_metrics(4, ll, true);

    const int MAXMSG = 206;
    float dmm[206];
    float pbest = 0.0f;
    dmin = 1.e30f;
    nharderrors = -1;
    std::string msgbest(37, ' ');
    // The message's character*37, with Fortran's assignment and substring
    // rules.
    auto set = [&](std::string& m, const std::string& v) { m = fix(v, 37); };
    auto put = [&](std::string& m, int j0, const std::string& v) {  // m(j0:j0+len-1) = v
        for (size_t q = 0; q < v.size() && j0 - 1 + q < 37; ++q) m[j0 - 1 + q] = v[q];
    };
    auto i32 = [](int v) {  // (i3.2)
        char b[8];
        std::snprintf(b, sizeof b, "%s%02d", v < 0 ? "-" : "", v < 0 ? -v : v);
        std::string s = b;
        if (s.size() < 3) s = std::string(3 - s.size(), ' ') + s;
        return s.size() > 3 ? std::string("***") : s;
    };
    for (int imsg = 1; imsg <= MAXMSG; ++imsg) {
        std::string msg;
        set(msg, trim(call_1) + " " + trim(call_2));
        const int i = imsg;
        if (call_1.compare(0, 3, "CQ ") == 0 && i != 5) set(msg, "QU1RK " + trim(call_2));
        if (!std_1) {
            if (i == 1 || i >= 6) set(msg, "<" + trim(call_1) + "> " + trim(call_2));
            if (i >= 2 && i <= 4) set(msg, trim(call_1) + " <" + trim(call_2) + ">");
        } else if (!std_2) {
            if (i <= 4 || i == 6) set(msg, "<" + trim(call_1) + "> " + trim(call_2));
            if (i >= 7) set(msg, trim(call_1) + " <" + trim(call_2) + ">");
        }
        const int j0 = static_cast<int>(trim(msg).size()) + 2;
        if (i == 2) put(msg, j0, "RRR");
        if (i == 3) put(msg, j0, "RR73");
        if (i == 4) put(msg, j0, "73");
        if (i == 5) {
            if (std_2) {
                set(msg, "CQ " + trim(call_2));
                if (call_1[2] == '_') set(msg, trim(call_1) + " " + trim(call_2));
                if (grid4 != "RR73") set(msg, trim(msg) + " " + grid4);
            }
            if (!std_2) set(msg, "CQ " + trim(call_2));
        }
        if (i == 6 && std_2) put(msg, j0, grid4);
        if (i >= 7) {
            const int isnr = -50 + (i - 7) / 2;
            if ((i & 1) == 1) {
                put(msg, j0, i32(isnr));
                if (msg[j0 - 1] == ' ') msg[j0 - 1] = '+';
            } else {
                put(msg, j0, "R" + i32(isnr));
                if (msg[j0] == ' ') msg[j0] = '+';
            }
        }

        const Codec::Packed p = codec.genft8(msg);  // Source-encode this message
        int itone[NN];
        int8_t msgbits[77];
        for (int q = 0; q < 77; ++q) msgbits[q] = p.msgbits[q];
        if (p.ok)
            get_ft8_tones_from_77bits(msgbits, itone);
        else
            std::fill(itone, itone + NN, 0);
        int8_t cw[174];
        encode174_91(msgbits, cw);  // Get codeword for this message
        float pow = 0.0f;
        for (int q = 1; q <= 79; ++q) pow = pow + s8_[q - 1][itone[q - 1]] * s8_[q - 1][itone[q - 1]];

        auto dist = [&](const float* llr) {
            float d = 0.0f;
            for (int q = 0; q < 174; ++q) {
                const int hdec = llr[q] >= 0.0f ? 1 : 0;
                const int nxor = hdec ^ cw[q];
                d = d + static_cast<float>(nxor) * std::fabs(llr[q]);
            }
            return d;
        };
        const float da = dist(llra), dbb = dist(llrb), dc = dist(llrc), dd = dist(llrd);
        float dm = da;
        if (dbb < dm) dm = dbb;
        if (dc < dm) dm = dc;
        if (dd < dm) dm = dd;
        dmm[imsg - 1] = dm;
        if (dm < dmin) {
            dmin = dm;
            msgbest = p.msgsent;
            pbest = pow;
            auto count = [&](const float* llr) {
                int n = 0;
                for (int q = 0; q < 174; ++q)
                    if (static_cast<float>(2 * cw[q] - 1) * llr[q] < 0.0f) ++n;
                return n;
            };
            if (dm == da)
                nharderrors = count(llra);
            else if (dm == dbb)
                nharderrors = count(llrb);
            else if (dm == dc)
                nharderrors = count(llrc);
            else if (dm == dd)
                nharderrors = count(llrd);
        }
    }
    if (trace_)
        trace_("D " + X(h_r(dmm, 206)) + " " + XF(pbest) + " " + XF(dmin) + " " + I(nharderrors) + " |" + msgbest +
               "|");

    auto minloc = [&]() {
        int best = -1;
        float m = 0.0f;
        for (int q = 0; q < 206; ++q) {
            const float v = dmm[q];
            if (best < 0) {
                if (!std::isnan(v)) {
                    best = q;
                    m = v;
                }
            } else if (v < m) {
                best = q;
                m = v;
            }
        }
        return best < 0 ? 0 : best;
    };
    dmm[minloc()] = 1.e30f;
    const float dmin2 = dmm[minloc()];
    xsnr = -25.0f;
    const float arg = pbest / xbase / 3.0e6f - 1.0f;
    if (arg > 0.0f) {
        float db = -99.0f;
        if (arg > 1.259e-10f) db = 10.0f * log10f_glibc(arg);
        const float v = db - 27.0f;
        xsnr = v > -25.0f ? v : -25.0f;
    }
    if (dmin > 100.0f || dmin2 / (dmin + 1e-6f) < 1.3f) nharderrors = -1;
    msg37 = msgbest;
    if (msg37.compare(0, 3, "CQ ") == 0 && std_2 && grid4 == "    ") nharderrors = -1;
    if (msg37.compare(0, 6, "QU1RK ") == 0) nharderrors = -1;
    if (trace_)
        trace_("H " + I(nharderrors) + " " + XF(dmin) + " " + XF(dmin2) + " " + XF(xsnr) + " |" + msg37 + "|");
}

}  // namespace ft8
