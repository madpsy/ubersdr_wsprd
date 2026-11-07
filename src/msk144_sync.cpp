// MSK144's short-ping detector and sync (msk144_sync.h): WSJT-X's
// msk144spd.f90, msk144sync.f90, msk_spectrum.f90 and indexx.f90 at 567ad29ce,
// line for line. Arrays keep the Fortran's order; the 1-based loops of the
// Fortran are shifted to 0-based here with the Fortran's index in a comment
// where it helps.
//
// Joe Taylor, K1JT, Steven Franke, K9AN, and the WSJT-X team (GNU GPL v3),
// and so is this.
#include "msk144_sync.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <vector>

#include "ft8_maths.h"

namespace msk144 {

SyncTrace sync_trace;
SpdTrace spd_trace;

namespace {

constexpr int NSPM = 864;

// complex products as gfortran writes them (-fcx-fortran-rules): the plain
// formula, no NaN recovery.
inline cfloat cmul(cfloat a, cfloat b) {
    const float ar = a.real(), ai = a.imag(), br = b.real(), bi = b.imag();
    return cfloat(ar * br - ai * bi, ar * bi + ai * br);
}

// nint of a float (the reference calls lroundf) and of a double (lround), as the default integer
// (llround: a 64-bit result everywhere, so NaN and infinity truncate to 0 on Windows, where long is 32 bits, as on Linux)
inline int nint_f(float x) { return static_cast<int>(llroundf(x)); }
inline int nint_d(double x) { return static_cast<int>(llround(x)); }

// MAX(a, b) as the reference's msk144spd leaves it (objdump: maxss with b, the later argument, in the
// destination): b if it is greater than a, else a: a when either is NaN
inline float max_gf(float a, float b) { return b > a ? b : a; }

// ---- msk144sync's constants: pp (the half-sine pulse) and cb (the sync word), made at the first call.
// (msk144sync.f90:30-50.) The reference's compiler vectorised the loop for pp(i)=sin((i-1)*pi/12.0),
// pi=4.0*atan(1.0) a default real: libmvec's _ZGVbN4v_sinf for i = 2..9, the first (sin 0) and the last
// three folded to the correctly rounded constants. Provenance of the folded ones: gfortran 13.3 in the
// wsjtx-oracle image, read from its printout (the oracle's P and B lines); re-derive them as the
// correctly rounded sinf of the float angle of i = 10, 11, 12. The oracle's P/B/R lines check them.
struct SyncConsts {
    float pp[12];
    cfloat cb[42];
    SyncConsts() {
        const float pi = 4.0f * 0.785398185253143310546875f;  // 4.0*atan(1.0): atan(1.0) folded to the float nearest pi/4
        float ang[12];
        for (int i = 1; i <= 12; ++i) ang[i - 1] = static_cast<float>(i - 1) * pi / 12.0f;
        pp[0] = 0.0f;
        for (int g = 0; g < 2; ++g) ft8::sinf4_glibc(ang + 1 + 4 * g, pp + 1 + 4 * g);
        const uint32_t folded[3] = {0x3F3504F3u, 0x3EFFFFFAu, 0x3E8483EAu};  // sin of the angles of i = 10, 11, 12
        for (int i = 0; i < 3; ++i) std::memcpy(&pp[9 + i], &folded[i], 4);

        int s8[8] = {0, 1, 1, 1, 0, 0, 1, 0};
        for (int& v : s8) v = 2 * v - 1;
        float cbi[42], cbq[42];
        for (int i = 0; i < 6; ++i) cbq[i] = pp[6 + i] * static_cast<float>(s8[0]);
        for (int i = 0; i < 12; ++i) cbq[6 + i] = pp[i] * static_cast<float>(s8[2]);
        for (int i = 0; i < 12; ++i) cbq[18 + i] = pp[i] * static_cast<float>(s8[4]);
        for (int i = 0; i < 12; ++i) cbq[30 + i] = pp[i] * static_cast<float>(s8[6]);
        for (int i = 0; i < 12; ++i) cbi[i] = pp[i] * static_cast<float>(s8[1]);
        for (int i = 0; i < 12; ++i) cbi[12 + i] = pp[i] * static_cast<float>(s8[3]);
        for (int i = 0; i < 12; ++i) cbi[24 + i] = pp[i] * static_cast<float>(s8[5]);
        for (int i = 0; i < 6; ++i) cbi[36 + i] = pp[i] * static_cast<float>(s8[7]);
        for (int i = 0; i < 42; ++i) cb[i] = cfloat(cbi[i], cbq[i]);
        if (sync_trace.consts) sync_trace.consts(pp, cb);
    }
};

// ---- msk144spd's constants: rcw, the raised-cosine edge window (msk144spd.f90:62-80).
// angle=(i-1)*pi/12.0 with pi a real*8, rcw(i)=(1-cos(angle))/2: the same vectorisation (libmvec's
// _ZGVbN4v_cosf for i = 2..9; i = 1 is cos 0, i = 10..12 folded: as for pp, gfortran 13.3's constants,
// re-derived as float (1 - c)/2 with c the correctly rounded cosf of the float angle; the oracle's R line
// checks them).
struct SpdConsts {
    float rcw[12];
    SpdConsts() {
        const double pi = 4.0 * std::atan(1.0);
        float ang[12], cs[12];
        for (int i = 1; i <= 12; ++i) ang[i - 1] = static_cast<float>((static_cast<double>(i - 1) * pi) / 12.0);
        cs[0] = 1.0f;
        for (int g = 0; g < 2; ++g) ft8::cosf4_glibc(ang + 1 + 4 * g, cs + 1 + 4 * g);
        const uint32_t folded[3] = {0x3F5A827Au, 0x3F6ED9ECu, 0x3F7BA375u};  // (1-cos)/2 for i = 10, 11, 12
        for (int i = 0; i < 9; ++i) rcw[i] = (1.0f - cs[i]) / 2.0f;
        for (int i = 0; i < 3; ++i) std::memcpy(&rcw[9 + i], &folded[i], 4);
        if (spd_trace.consts) spd_trace.consts(rcw);
    }
};

}  // namespace

// ---- indexx.f90:1-91 (Numerical Recipes' quicksort of indices), the gotos kept
void indexx(const float* arr, int n, int* indx) {
    constexpr int M = 7, NSTACK = 50;
    int istack[NSTACK + 1];
    auto A = [&](int i) -> float { return arr[i - 1]; };
    auto X = [&](int i) -> int& { return indx[i - 1]; };
    int i, indxt, ir, itemp, j, jstack, k, l;
    float a;

    for (j = 1; j <= n; ++j) X(j) = j;

    jstack = 0;
    l = 1;
    ir = n;
L1:
    if (ir - l < M) {
        for (j = l + 1; j <= ir; ++j) {
            indxt = X(j);
            a = A(indxt);
            for (i = j - 1; i >= 1; --i) {
                if (A(X(i)) <= a) goto L2;
                X(i + 1) = X(i);
            }
            i = 0;
        L2:
            X(i + 1) = indxt;
        }
        if (jstack == 0) return;

        ir = istack[jstack];
        l = istack[jstack - 1];
        jstack = jstack - 2;
    } else {
        k = (l + ir) / 2;
        itemp = X(k);
        X(k) = X(l + 1);
        X(l + 1) = itemp;

        if (A(X(l + 1)) > A(X(ir))) {
            itemp = X(l + 1);
            X(l + 1) = X(ir);
            X(ir) = itemp;
        }

        if (A(X(l)) > A(X(ir))) {
            itemp = X(l);
            X(l) = X(ir);
            X(ir) = itemp;
        }

        if (A(X(l + 1)) > A(X(l))) {
            itemp = X(l + 1);
            X(l + 1) = X(l);
            X(l) = itemp;
        }

        i = l + 1;
        j = ir;
        indxt = X(l);
        a = A(indxt);
    L3:
        i = i + 1;
        if (A(X(i)) < a) goto L3;

    L4:
        j = j - 1;
        if (A(X(j)) > a) goto L4;
        if (j < i) goto L5;
        itemp = X(i);
        X(i) = X(j);
        X(j) = itemp;
        goto L3;

    L5:
        X(l) = X(j);
        X(j) = indxt;
        jstack = jstack + 2;
        if (jstack > NSTACK) throw std::runtime_error("NSTACK too small in indexx");
        if (ir - i + 1 >= j - l) {
            istack[jstack] = ir;
            istack[jstack - 1] = i;
            ir = j - 1;
        } else {
            istack[jstack] = j - 1;
            istack[jstack - 1] = l;
            l = i;
        }
    }
    goto L1;
}

// ---- MAXLOC as gfortran -O3 writes it for a real array: elements that are NaN are passed over, the first
// other element starts the search, later ones replace it if greater; with nothing but NaN, the first
// element (of the mask)
int maxloc_f(const float* a, int n) {
    int i = 0, pos = 0;
    float limit = 0.0f;
    for (; i < n; ++i) {
        if (!std::isnan(a[i])) {
            pos = i + 1;
            limit = a[i];
            ++i;
            break;
        }
    }
    if (pos == 0) return n > 0 ? 1 : 0;
    for (; i < n; ++i) {
        if (a[i] > limit) {
            limit = a[i];
            pos = i + 1;
        }
    }
    return pos;
}

int maxloc_masked_f(const float* a, int lo, int hi) {
    int i = lo, pos = 0;
    float limit = 0.0f;
    for (; i <= hi; ++i) {
        if (!std::isnan(a[i - 1])) {
            pos = i;
            limit = a[i - 1];
            ++i;
            break;
        }
    }
    if (pos == 0) return lo <= hi ? lo : 0;
    for (; i <= hi; ++i) {
        if (a[i - 1] > limit) {
            limit = a[i - 1];
            pos = i;
        }
    }
    return pos;
}

// ---- msk_spectrum.f90: msk_clip_spectrum_window lines 12-36, msk_peak_offset lines 38-57
void msk_clip_spectrum_window(int nfft, double requested_lo, double requested_hi, int& lo, int& hi, bool& usable) {
    lo = 1;
    hi = 0;
    usable = false;
    if (nfft < 1) return;
    if (!std::isfinite(requested_lo) || !std::isfinite(requested_hi)) return;
    if (requested_lo > requested_hi) return;
    const double dn = static_cast<double>(nfft);
    if (requested_hi < 0.5 || requested_lo >= dn + 0.5) return;

    lo = nint_d(std::max(1.0, std::min(dn, requested_lo)));
    hi = nint_d(std::max(1.0, std::min(dn, requested_hi)));
    usable = lo <= hi;
    if (!usable) {
        lo = 1;
        hi = 0;
    }
}

float msk_peak_offset(const cfloat* spectrum, int n, int peak_bin) {
    float offset = 0.0f;
    // lbound 1, ubound n
    if (peak_bin <= 1 || peak_bin >= n) return offset;

    const cfloat sm = spectrum[peak_bin - 2], s0 = spectrum[peak_bin - 1], sp = spectrum[peak_bin];
    // curvature=2.0*spectrum(peak_bin)-spectrum(peak_bin-1)-spectrum(peak_bin+1): the real 2.0 made complex
    cfloat curvature = cmul(cfloat(2.0f, 0.0f), s0);
    curvature = cfloat(curvature.real() - sm.real(), curvature.imag() - sm.imag());
    curvature = cfloat(curvature.real() - sp.real(), curvature.imag() - sp.imag());
    // scale=max(abs(a),abs(b),abs(c)): as the reference's code leaves it (maxss: the later operand if either is NaN)
    const float ab = ft8::hypotf_glibc(sm.real(), sm.imag());
    const float a0 = ft8::hypotf_glibc(s0.real(), s0.imag());
    const float ap = ft8::hypotf_glibc(sp.real(), sp.imag());
    float scale = ab > a0 ? ab : a0;
    scale = scale > ap ? scale : ap;
    if (scale <= 0.0f) return offset;
    constexpr float eps = 1.1920928955078125e-07f;  // epsilon(scale)
    if (ft8::hypotf_glibc(curvature.real(), curvature.imag()) <= eps * scale) return offset;

    // candidate=-real((spectrum(peak_bin-1)-spectrum(peak_bin+1))/curvature): gfortran's complex division
    // with Fortran rules (range reduction by the larger of |c|, |d|), the real part only
    const float a = sm.real() - sp.real(), b = sm.imag() - sp.imag();
    const float c = curvature.real(), d = curvature.imag();
    float tr;
    if (std::fabs(d) > std::fabs(c)) {
        const float ratio = c / d;
        const float div = ratio * c + d;
        tr = (a * ratio + b) / div;
    } else {
        const float ratio = d / c;
        const float div = d * ratio + c;
        tr = (a + b * ratio) / div;
    }
    const float candidate = -tr;
    if (std::isfinite(candidate)) offset = candidate;
    return offset;
}

// ---- msk144sync.f90:1-101
void msk144sync(const cfloat* cdat, int nframes, int ntol, float delf, const int* navmask, int npeaks, float fc, float& fest,
                int* npklocs, int& nsuccess, float& xmax, cfloat* c) {
    static const SyncConsts consts;
    std::vector<cfloat> cdat2(static_cast<size_t>(NSPM) * static_cast<size_t>(std::max(nframes, 0)) + 1);
    cfloat cs[NSPM];
    float xcc[NSPM];  // xcc(0:NSPM-1)
    float xccs[NSPM];
    float xm = 0.0f, bf = 0.0f;

    const int nfreqs = 2 * nint_f(static_cast<float>(ntol) / delf) + 1;
    // nthreads = 1: id = 1
    const int nstep = nfreqs / 1;
    const int if1 = -nint_f(static_cast<float>(ntol) / delf) + (1 - 1) * nstep;
    // if2=if1+nstep-1, then, as id == nthreads, if2=nint(ntol/delf)
    const int if2 = nint_f(static_cast<float>(ntol) / delf);
    msk144_freq_search(cdat, fc, if1, if2, delf, nframes, navmask, consts.cb, cdat2.data(), xm, bf, cs, xccs);

    xmax = xm;
    fest = fc + bf;
    for (int i = 0; i < NSPM; ++i) c[i] = cs[i];
    for (int i = 0; i < NSPM; ++i) xcc[i] = xccs[i];

    if (sync_trace.peak) sync_trace.peak(xmax, fest, c, xcc, nfreqs);

    // Find npeaks largest peaks
    for (int ipk = 1; ipk <= npeaks; ++ipk) {
        const int ic2 = maxloc_f(xcc, NSPM);
        npklocs[ipk - 1] = ic2;
        // pkamps(ipk)=xcc(ic2-1) is never used
        for (int j = std::max(0, ic2 - 7); j <= std::min(NSPM - 1, ic2 + 7); ++j) xcc[j] = 0.0f;
    }

    nsuccess = 0;
    if (xmax >= 1.3f) nsuccess = 1;
}

// ---- msk144spd.f90:1-208 (the calls to msk144decodeframe at lines 194 go to `decode`)
void msk144spd(const cfloat* cbig, int n, int ntol, int& nsuccess, char* msgreceived, float fc, float& fret, float& tret,
               int& navg, cfloat* ct, float* softbits, const DecodeSink& decode) {
    constexpr int MAXSTEPS = 100, NFFT = NSPM, MAXCAND = 5, NPATTERNS = 6;
    static const int navpatterns[NPATTERNS][3] = {{0, 1, 0}, {1, 0, 0}, {0, 0, 1}, {1, 1, 0}, {0, 1, 1}, {1, 1, 1}};
    // data tpat/1.5,0.5,2.5,1.0,2.0,1.5/ is never used
    constexpr double fs = 12000.0;
    constexpr double df = fs / NFFT;

    nsuccess = 0;
    std::memset(msgreceived, ' ', 37);
    if (n < 3 * NSPM || ntol < 0 || !std::isfinite(fc)) return;

    static const SpdConsts consts;
    const float* rcw = consts.rcw;

    cfloat cdat[3 * NSPM];
    cfloat c[NSPM];
    alignas(32) cfloat ctmp[NFFT];
    int indices[MAXSTEPS];
    int npkloc[10] = {0};
    int navmask[3];
    int nstart[MAXCAND];
    float detmet[MAXSTEPS + 6];   // detmet(-2:MAXSTEPS+3): Fortran index k is [k+2]
    float detmet2[MAXSTEPS + 6];
    float detfer[MAXSTEPS];
    float ferrs[MAXCAND];
    float snrs[MAXCAND];
    float tonespec[NFFT];
    auto DM = [&](int k) -> float& { return detmet[k + 2]; };
    auto DM2 = [&](int k) -> float& { return detmet2[k + 2]; };

    // fill the detmet, detferr arrays
    const int nstep = std::min(MAXSTEPS, (n - NSPM) / 216);  // 72ms/4=18ms steps
    for (float& v : detmet) v = 0.0f;
    for (float& v : detmet2) v = 0.0f;
    for (float& v : detfer) v = -999.99f;
    const double nfhi = 2.0 * (static_cast<double>(fc) + 500.0);
    const double nflo = 2.0 * (static_cast<double>(fc) - 500.0);
    const double ihlo_requested = (nfhi - 2.0 * static_cast<double>(ntol)) / df + 1.0;
    const double ihhi_requested = (nfhi + 2.0 * static_cast<double>(ntol)) / df + 1.0;
    const double illo_requested = (nflo - 2.0 * static_cast<double>(ntol)) / df + 1.0;
    const double ilhi_requested = (nflo + 2.0 * static_cast<double>(ntol)) / df + 1.0;
    int ihlo, ihhi, illo, ilhi;
    bool high_window_usable, low_window_usable;
    msk_clip_spectrum_window(NFFT, ihlo_requested, ihhi_requested, ihlo, ihhi, high_window_usable);
    msk_clip_spectrum_window(NFFT, illo_requested, ilhi_requested, illo, ilhi, low_window_usable);
    if (!high_window_usable || !low_window_usable) return;
    const int i2000 = nint_d(nflo / df) + 1;
    const int i4000 = nint_d(nfhi / df) + 1;
    if (spd_trace.windows) spd_trace.windows(ihlo, ihhi, illo, ilhi, i2000, i4000);
    for (int istp = 1; istp <= nstep; ++istp) {
        const int ns = 1 + 216 * (istp - 1);
        const int ne = ns + NSPM - 1;
        if (ne > n) break;
        for (int i = 0; i < NFFT; ++i) ctmp[i] = cfloat(0.0f, 0.0f);
        for (int i = 0; i < NSPM; ++i) ctmp[i] = cbig[ns - 1 + i];

        // Coarse carrier frequency sync - seek tones at 2000 Hz and 4000 Hz in
        // squared signal spectrum.
        for (int i = 0; i < NFFT; ++i) ctmp[i] = cmul(ctmp[i], ctmp[i]);
        for (int i = 0; i < 12; ++i) ctmp[i] = cmul(ctmp[i], cfloat(rcw[i], 0.0f));
        for (int i = 0; i < 12; ++i) ctmp[NSPM - 12 + i] = cmul(ctmp[NSPM - 12 + i], cfloat(rcw[11 - i], 0.0f));
        fft_c2c_inplace(ctmp, NFFT, -1);
        for (int i = 0; i < NFFT; ++i) {
            const float ab = ft8::hypotf_glibc(ctmp[i].real(), ctmp[i].imag());
            tonespec[i] = ab * ab;
        }

        // high tone search window ihlo:ihhi
        const int ihpk = maxloc_masked_f(tonespec, ihlo, ihhi);
        const float deltah = msk_peak_offset(ctmp, NFFT, ihpk);
        const float ah = tonespec[ihpk - 1];
        float sumh = 0.0f;
        for (int i = ihlo; i <= ihhi; ++i) sumh = sumh + tonespec[i - 1];
        const float ahavp = (sumh - ah) / static_cast<float>(ihhi - ihlo + 1);
        const float trath = ah / (ahavp + 0.01f);
        // window for low tone illo:ilhi
        const int ilpk = maxloc_masked_f(tonespec, illo, ilhi);
        const float deltal = msk_peak_offset(ctmp, NFFT, ilpk);
        const float al = tonespec[ilpk - 1];
        float suml = 0.0f;
        for (int i = illo; i <= ilhi; ++i) suml = suml + tonespec[i - 1];
        const float alavp = (suml - al) / static_cast<float>(ilhi - illo + 1);
        const float tratl = al / (alavp + 0.01f);
        // fdiff=(ihpk+deltah-ilpk-deltal)*df is never used
        // (ihpk+deltah-i4000)*df/2.0: a single sum, then the double df, the 2.0 a double
        const float ferrh = static_cast<float>(static_cast<double>((static_cast<float>(ihpk) + deltah) - static_cast<float>(i4000)) * df / 2.0);
        const float ferrl = static_cast<float>(static_cast<double>((static_cast<float>(ilpk) + deltal) - static_cast<float>(i2000)) * df / 2.0);
        float ferr;
        if (ah >= al) {
            ferr = ferrh;
        } else {
            ferr = ferrl;
        }
        DM(istp) = max_gf(ah, al);
        DM2(istp) = max_gf(trath, tratl);
        detfer[istp - 1] = ferr;
        if (spd_trace.step) spd_trace.step(istp, ihpk, ilpk, deltah, deltal, ah, al, trath, tratl, ferr);
    }  // end of detection-metric and frequency error estimation loop

    indexx(&DM(1), nstep, indices);  // find median of detection metric vector
    const float xmed = DM(indices[nstep / 4 - 1]);
    if (spd_trace.median) spd_trace.median(nstep, indices, xmed, &DM(1), &DM2(1), detfer);
    if (xmed <= 0.0f || !std::isfinite(xmed)) return;
    for (float& v : detmet) v = v / xmed;  // noise floor of detection metric is 1.0
    int ndet = 0;

    for (int ip = 1; ip <= MAXCAND; ++ip) {  // Find candidates
        const int il = maxloc_f(&DM(1), nstep);
        if (DM(il) < 3.0f) break;
        if (std::fabs(detfer[il - 1]) <= static_cast<float>(ntol)) {
            ndet = ndet + 1;
            nstart[ndet - 1] = 1 + (il - 1) * 216 + 1;
            ferrs[ndet - 1] = detfer[il - 1];
            snrs[ndet - 1] = 12.0f * ft8::log10f_glibc(DM(il)) / 2.0f - 9.0f;
        }
        DM(il) = 0.0f;
    }

    if (ndet < 3) {
        const int nloop = MAXCAND - ndet;  // a Fortran DO's trip count is made once, at the loop's start
        for (int ip = 1; ip <= nloop; ++ip) {  // Find candidates
            const int il = maxloc_f(&DM2(1), nstep);
            if (DM2(il) < 12.0f) break;
            if (std::fabs(detfer[il - 1]) <= static_cast<float>(ntol)) {
                ndet = ndet + 1;
                nstart[ndet - 1] = 1 + (il - 1) * 216 + 1;
                ferrs[ndet - 1] = detfer[il - 1];
                snrs[ndet - 1] = 12.0f * ft8::log10f_glibc(DM2(il)) / 2.0f - 9.0f;
            }
            DM2(il) = 0.0f;
        }
    }

    if (spd_trace.candidates) spd_trace.candidates(ndet, nstart, ferrs, snrs);
    constexpr int npeaks = 2;
    constexpr int ntol0 = 8;
    constexpr float deltaf = 2.0f;
    for (int icand = 1; icand <= ndet; ++icand) {  // Try to sync/demod/decode each candidate.
        int ib = std::max(1, nstart[icand - 1] - NSPM);
        int ie = ib - 1 + 3 * NSPM;
        if (ie > n) {
            ie = n;
            ib = ie - 3 * NSPM + 1;
        }
        for (int i = 0; i < 3 * NSPM; ++i) cdat[i] = cbig[ib - 1 + i];
        const float fo = fc + ferrs[icand - 1];
        for (int iav = 1; iav <= NPATTERNS; ++iav) {
            for (int i = 0; i < 3; ++i) navmask[i] = navpatterns[iav - 1][i];
            float fest, xmax;
            int nsyncsuccess;
            msk144sync(cdat, 3, ntol0, deltaf, navmask, npeaks, fo, fest, npkloc, nsyncsuccess, xmax, c);
            if (spd_trace.sync) spd_trace.sync(icand, iav, nsyncsuccess, fest, xmax, npkloc);

            if (nsyncsuccess == 0) continue;

            for (int ipk = 1; ipk <= npeaks; ++ipk) {
                for (int is = 1; is <= 3; ++is) {
                    int ic0 = npkloc[ipk - 1];
                    if (is == 2) ic0 = std::max(1, ic0 - 1);
                    if (is == 3) ic0 = std::min(NSPM, ic0 + 1);
                    // ct=cshift(c,ic0-1)
                    for (int i = 0; i < NSPM; ++i) ct[i] = c[(i + ic0 - 1) % NSPM];
                    int ndecodesuccess = 0;
                    decode(ct, softbits, msgreceived, ndecodesuccess);
                    if (ndecodesuccess > 0) {
                        tret = static_cast<float>(static_cast<double>(nstart[icand - 1] + NSPM / 2) / fs);
                        fret = fest;
                        navg = navmask[0] + navmask[1] + navmask[2];
                        nsuccess = 1;
                        return;
                    }
                }
            }
        }
    }  // candidate loop
}

}  // namespace msk144
