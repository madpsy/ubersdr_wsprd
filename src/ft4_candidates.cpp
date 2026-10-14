// FT4's front end (ft4_frontend.h): lib/ft4/getcandidates4.f90 and
// ft4_baseline.f90 at 567ad29ce, line for line (no receiver filter-edge
// trim: that block is improved 3.2.0's only).
//
// As the Release object computes them (from its disassembly): getcandidates4's x = (fac*dd)*window, fac = 1/300 folded;
// abs(cx)**2 as cabsf squared; savg/NHSYM a division by 122.0; savsm each
// 15-term sum in index order from 0, over 15.0; nfa = fa/df truncated,
// df = 12000/2304 single; the peak's parabola grouped as written.
// ft4_baseline (implicit real*8): df the single 12000/2304 widened;
// 10*log10 by libmvec log10f four at a time from ia (scalar log10f for the
// last (ib-ia+1) mod 4); the polynomial in double, + 0.65 (the single 0.65
// widened), rounded to single; then 10**(sbase/10.0) by libmvec's powf
// four at a time from ia (x = 10 in every lane), the rest scalar powf.
//
// Joe Taylor, K1JT, Steven Franke, K9AN, and the WSJT-X developers (GNU GPL
// v3), and so is this.
#include "ft4_frontend.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <vector>

#include "ft8_fortran.h"   // the Fortran helpers, FT8's
#include "ft8_frontend.h"  // nuttal_window, normalizebmet, gfsk_pulse
#include "ft8_maths.h"     // glibc 2.39's maths and x86-64 libmvec, carried

namespace ft4 {

using namespace ft8::fortran;
using ft8::log10f4_glibc;
using ft8::log10f_glibc;
using ft8::powf4_glibc;
using ft8::powf_glibc;

namespace {

// A real*8's bits as the trace's z16.16.
std::string X64(double d) {
    uint64_t u;
    std::memcpy(&u, &d, sizeof u);
    char b[24];
    std::snprintf(b, sizeof b, "%016llX", static_cast<unsigned long long>(u));
    return b;
}

}  // namespace

// This group's SAVEd state: getcandidates4's first and window; its x/cx
// (a stack local there, 16-byte aligned as FFTW's buffers are) and the
// '@ cand-x' line's last alignment.
struct Frontend::Cand {
    bool first = true;
    std::vector<float> window;
    fst4w::FftwVec<float> x = fst4w::FftwVec<float>(NFFT1 + 2);  // x(NFFT1) equivalenced with cx(0:NH1)
    int al = -1;
};

Frontend::Cand* Frontend::make_cand() { return new Cand(); }
void Frontend::del_cand(Cand* p) { delete p; }
size_t Frontend::cand_bytes(const Cand* p) {
    return p->window.capacity() * sizeof(float) + p->x.capacity() * sizeof(float);
}

void ft4_baseline(float* s0, int nfa, int nfb, float* sbase0, const Trace* trace) {
    float* s = s0 - 1;  // 1-based
    float* sbase = sbase0 - 1;
    const int nseg = 10, npct = 10;
    const double df = static_cast<double>(12000.0f / static_cast<float>(NFFT1));  // 5.21 Hz
    const int ia = std::max(static_cast<int>(std::lround(static_cast<double>(200.0f) / df)), nfa);
    const int ib = std::min(NH1, nfb);
    // s(i) = 10.0*log10(s(i)) (Convert to dB scale): libmvec four at a
    // time from ia, the rest scalar.
    {
        int i = ia;
        for (; i + 3 <= ib; i += 4) {
            float l[4];
            log10f4_glibc(&s[i], l);
            for (int q = 0; q < 4; ++q) s[i + q] = 10.0f * l[q];
        }
        for (; i <= ib; ++i) s[i] = 10.0f * log10f_glibc(s[i]);
    }
    if (trace && *trace) (*trace)("l " + I(ia) + " " + I(ib) + " " + X(h_r(&s[ia], static_cast<size_t>(std::max(0, ib - ia + 1)))));

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
    if (trace && *trace)
        (*trace)("p " + I(kz) + " " + X64(a[1]) + " " + X64(a[2]) + " " + X64(a[3]) + " " + X64(a[4]) + " " +
                 X64(a[5]));
    for (int i = ia; i <= ib; ++i) {
        const double t = static_cast<double>(i - i0);
        sbase[i] = static_cast<float>(a[1] + t * (a[2] + t * (a[3] + t * (a[4] + t * (a[5])))) +
                                      static_cast<double>(0.65f));
    }
    // sbase(i) = 10**(sbase(i)/10.0): libmvec's powf four at a time from
    // ia, the rest scalar.
    {
        const float ten[4] = {10.0f, 10.0f, 10.0f, 10.0f};
        int i = ia;
        for (; i + 3 <= ib; i += 4) {
            float e[4], p[4];
            for (int q = 0; q < 4; ++q) e[q] = sbase[i + q] / 10.0f;
            powf4_glibc(ten, e, p);
            for (int q = 0; q < 4; ++q) sbase[i + q] = p[q];
        }
        for (; i <= ib; ++i) sbase[i] = powf_glibc(10.0f, sbase[i] / 10.0f);
    }
}

void Frontend::getcandidates4(const float* dd, float fa, float fb, float syncmin, int nfqso, int maxcand,
                              float* savg0, float* candidate, int& ncand, float* sbase0) {
    Cand& st = *cand_;
    float* savg = savg0 - 1;  // 1-based
    float* sbase = sbase0 - 1;
    const bool tr = static_cast<bool>(trace_);
    if (tr)
        trace_("S " + X(h_r(dd, NMAX)) + " " + XF(fa) + " " + XF(fb) + " " + XF(syncmin) + " " + I(nfqso) + " " +
               I(maxcand));
    if (st.first) {
        st.first = false;
        st.window.assign(NFFT1, 0.0f);
        ft8::nuttal_window(st.window.data(), NFFT1);
    }

    // Compute symbol spectra, stepping by NSTEP steps.
    for (int i = 1; i <= NH1; ++i) savg[i] = 0.0f;
    const float df = 12000.0f / static_cast<float>(NFFT1);
    const float fac = 1.0f / 300.0f;
    float* x = st.x.data();
    uint32_t hs = 2166136261u;  // h(s(1:NH1,1:NHSYM)), s itself never read again
    for (int j = 1; j <= NHSYM; ++j) {
        const int ia = (j - 1) * NSTEP + 1;
        const int ib = ia + NFFT1 - 1;
        if (ib > NMAX) break;
        for (int i = 0; i < NFFT1; ++i) x[i] = (fac * dd[ia - 1 + i]) * st.window[i];
        if (tr) {
            const int al = static_cast<int>(reinterpret_cast<uintptr_t>(x) % 16);
            if (al != st.al) {
                st.al = al;
                trace_("@ cand-x " + I(al));
            }
        }
        fst4w::fft_r2c_inplace(x, NFFT1);  // r2c FFT
        float mag[NH1 + 1];  // abs(cx(1:NH1)), cabsf's bits four at a time
        ft8::hypotf_c_glibc(x + 2, mag + 1, NH1);
        for (int i = 1; i <= NH1; ++i) {
            const float c = mag[i];
            const float s = c * c;
            if (tr) hs = hw(hs, fbits(s));
            savg[i] = savg[i] + s;  // Average spectrum
        }
    }
    for (int i = 1; i <= NH1; ++i) savg[i] = savg[i] / static_cast<float>(NHSYM);
    if (tr) trace_("s " + X(hs) + " " + X(h_r(savg0, NH1)));
    float savsm[NH1 + 1];  // savsm(1:NH1) at [1..NH1]
    for (int i = 1; i <= NH1; ++i) savsm[i] = 0.0f;
    for (int i = 8; i <= NH1 - 7; ++i) {
        float sum = 0.0f;
        for (int q = i - 7; q <= i + 7; ++q) sum = sum + savg[q];
        savsm[i] = sum / 15.0f;
    }

    int nfa = static_cast<int>(fa / df);
    if (nfa < nint(200.0f / df)) nfa = nint(200.0f / df);
    int nfb = static_cast<int>(fb / df);
    if (nfb > nint(4910.0f / df)) nfb = nint(4910.0f / df);
    if (tr) trace_("a " + I(nfa) + " " + I(nfb) + " " + X(h_r(&savsm[1], NH1)));
    ft4_baseline(savg0, nfa, nfb, sbase0, tr ? &trace_ : nullptr);
    if (tr)
        trace_("A " + X(h_r(savg0, NH1)) + " " +
               X(h_r(&sbase[nfa], static_cast<size_t>(std::max(0, nfb - nfa + 1)))));
    for (int i = nfa; i <= nfb; ++i)
        if (sbase[i] <= 0.0f) {
            if (tr) trace_("k E");
            return;
        }
    for (int i = nfa; i <= nfb; ++i) savsm[i] = savsm[i] / sbase[i];
    if (tr) trace_("n " + X(h_r(&savsm[1], NH1)));
    const float f_offset = -1.5f * 12000.0f / static_cast<float>(NSPS);
    ncand = 0;
    std::vector<float> candt(static_cast<size_t>(2) * static_cast<size_t>(std::max(0, maxcand)), 0.0f);
    for (int i = nfa + 1; i <= nfb - 1; ++i) {
        if (savsm[i] >= savsm[i - 1] && savsm[i] >= savsm[i + 1] && savsm[i] >= syncmin) {
            const float den = (savsm[i - 1] - 2.0f * savsm[i]) + savsm[i + 1];
            float del = 0.0f;
            if (den != 0.0f) del = (0.5f * (savsm[i - 1] - savsm[i + 1])) / den;
            const float fpeak = (static_cast<float>(i) + del) * df + f_offset;
            if (fpeak < 200.0f || fpeak > 4910.0f) continue;
            const float speak = savsm[i] - (0.25f * (savsm[i - 1] - savsm[i + 1])) * del;
            ncand = ncand + 1;
            candt[static_cast<size_t>(ncand - 1) * 2] = fpeak;
            candt[static_cast<size_t>(ncand - 1) * 2 + 1] = speak;
            if (ncand == maxcand) break;
        }
    }
    if (tr) trace_("K " + I(ncand) + " " + X(h_r(candt.data(), static_cast<size_t>(2 * ncand))));
    for (int i = 0; i < 2 * maxcand; ++i) candidate[i] = 0.0f;
    const float fq = static_cast<float>(nfqso);
    int nq = 0;
    for (int i = 0; i < ncand; ++i)
        if (std::fabs(candt[static_cast<size_t>(i) * 2] - fq) <= 20.0f) ++nq;
    int n1 = 1;
    int n2 = nq + 1;
    for (int i = 1; i <= ncand; ++i) {
        const float* c = &candt[static_cast<size_t>(i - 1) * 2];
        int& n = std::fabs(c[0] - fq) <= 20.0f ? n1 : n2;
        candidate[(n - 1) * 2] = c[0];
        candidate[(n - 1) * 2 + 1] = c[1];
        n = n + 1;
    }
    if (tr) {
        std::string l = "k " + I(ncand);
        for (int i = 0; i < 2 * ncand; ++i) l += " " + XF(candidate[i]);
        trace_(l);
    }
}

}  // namespace ft4
