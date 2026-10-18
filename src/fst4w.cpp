// FST4W decoding (fst4w.h): WSJT-X lib/fst4_decode.f90's decode with
// iwspr=1 and its dopspread, lib/fst4/gen_fst4wave.f90, lib/gfsk_pulse.f90
// and genfst4.f90's get_fst4_tones_from_bits, at 567ad29ce, line for line.
// The FST4 (iwspr=0) branches and AP decoding are not here: FST4W never
// reaches them.
//
// Steven Franke, K9AN, and Joe Taylor, K1JT, WSJT-X (GNU GPL v3), and so is
// this.
#include "fst4w.h"

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <limits>

#include "fst4w_frontend.h"
#include "fst4w_ldpc.h"
#include "fst4w_simd.h"
#include "fst4w_unpack.h"


namespace fst4w {

namespace {

constexpr int NN = 160, ND = 120;
constexpr int MAXWCALLS = 100;
constexpr int NTAB = 65536;

// Fortran's 37-character message, blank padded.
std::string pad37(const std::string& s) { return (s + std::string(37, ' ')).substr(0, 37); }
std::string rtrim(std::string s) {
    while (!s.empty() && s.back() == ' ') s.pop_back();
    return s;
}
// index(s, t) (1-based, 0 when absent).
int findex(const std::string& s, const std::string& t) {
    const size_t p = s.find(t);
    return p == std::string::npos ? 0 : static_cast<int>(p) + 1;
}

// gfortran's MAX(a, b) and MIN(a, b) of reals: a NaN argument is passed over.
float fmax_g(float a, float b) {
    float r = a;
    if (b > r || std::isnan(r)) r = b;
    return r;
}
float fmin_g(float a, float b) {
    float r = a;
    if (b < r || std::isnan(r)) r = b;
    return r;
}

// gfsk_pulse(b, t), c = pi*sqrt(2/log(2)) as gfortran folds it.
float gfsk_pulse(float b, float t) {
    static const float c = [] {
        const uint32_t u = 0x40aac42bu;
        float f;
        std::memcpy(&f, &u, 4);
        return f;
    }();
    const float cb = c * b;
    return 0.5f * (erff(cb * (t + 0.5f)) - erff(cb * (t - 0.5f)));
}

// get_fst4_tones_from_bits(msgbits, i4tone, iwspr=1).
void tones_from_bits(const int8_t message74[74], int itone[NN]) {
    static const int isyncword1[8] = {0, 1, 3, 2, 1, 0, 2, 3};
    static const int isyncword2[8] = {2, 3, 1, 0, 3, 2, 0, 1};
    int8_t codeword[240];
    encode240_74(message74, codeword);
    int itmp[ND];
    for (int i = 1; i <= ND; ++i) {  // Grayscale mapping: 00 0, 01 1, 11 2, 10 3
        const int is = codeword[2 * i - 1] + 2 * codeword[2 * i - 2];
        if (is <= 1) itmp[i - 1] = is;
        if (is == 2) itmp[i - 1] = 3;
        if (is == 3) itmp[i - 1] = 2;
    }
    for (int i = 0; i < 8; ++i) {
        itone[i] = isyncword1[i];
        itone[38 + i] = isyncword2[i];
        itone[76 + i] = isyncword1[i];
        itone[114 + i] = isyncword2[i];
        itone[152 + i] = isyncword1[i];
    }
    for (int i = 0; i < 30; ++i) {
        itone[8 + i] = itmp[i];
        itone[46 + i] = itmp[30 + i];
        itone[84 + i] = itmp[60 + i];
        itone[122 + i] = itmp[90 + i];
    }
}

// The Fortran's x = x*(1.0 -/+ cos(twopi*i/real(nsps/2)))/2.0 over cwave's
// ramps: libmvec's cosf four at a time from the first element, the rest
// scalar, as the reference's vectorised loops.
void ramp(cfloat* cw, int n, float twopi, float den, bool up) {
    int i = 0;
    if (n > 3) {
        for (; i + 3 < n; i += 4) {
            float x[4], c[4];
            for (int l = 0; l < 4; ++l) x[l] = twopi * static_cast<float>(i + l) / den;
            simd::cosf4(x, c);
            for (int l = 0; l < 4; ++l) {
                const float f = (up ? (1.0f - c[l]) : (1.0f + c[l])) / 2.0f;
                cw[i + l] = cfloat(cw[i + l].real() * f, cw[i + l].imag() * f);
            }
        }
    }
    for (; i < n; ++i) {
        const float c = std::cos(twopi * static_cast<float>(i) / den);
        const float f = (up ? (1.0f - c) : (1.0f + c)) / 2.0f;
        cw[i] = cfloat(cw[i].real() * f, cw[i].imag() * f);
    }
}

// The fast spread's low-pass: a Kaiser-windowed sinc at 12 kHz, 1201 taps,
// cut-off 50 Hz, beta 10 (about 100 dB), unit DC gain.
const std::vector<float>& spread_filter() {
    static const std::vector<float> h = [] {
        const int L = 1201, half = L / 2;
        const double fcut = 50.0 / 12000.0, beta = 10.0;
        auto i0 = [](double x) {
            double sum = 1, term = 1;
            for (int k = 1; k < 60; ++k) {
                term *= (x / (2 * k)) * (x / (2 * k));
                sum += term;
            }
            return sum;
        };
        std::vector<double> d(L);
        double total = 0;
        for (int n = 0; n < L; ++n) {
            const double t = n - half;
            const double sinc = t == 0 ? 2 * fcut : std::sin(2 * M_PI * fcut * t) / (M_PI * t);
            const double r = t / half;
            d[n] = sinc * i0(beta * std::sqrt(std::max(0.0, 1 - r * r))) / i0(beta);
            total += d[n];
        }
        std::vector<float> f(L);
        for (int n = 0; n < L; ++n) f[n] = static_cast<float>(d[n] / total);
        return f;
    }();
    return h;
}

}  // namespace

struct Decoder::Impl {
    Unpacker unpacker;
    Ldpc ldpc;
    Frontend fe;
    bool first = true;
    std::vector<std::string> wcalls;  // character*20 entries, trimmed
    std::vector<std::string> pending_calls;
    bool have_pending = false;
    float fcbest = 0.0f;  // data fcbest/0.0/ (saved)
    int isbest_saved = 0;
    float w50_saved = 0.0f;  // w50 when no dopspread ran (never set: 0)
    int ncall = 0;           // dopspread's call count (saved)
    const std::atomic<bool>* stop = nullptr;  // set_stop()
    bool stopped = false;                     // the last decode was
    bool halt() { return stopped || (stop && stop->load(std::memory_order_relaxed) && (stopped = true)); }
    // gen_fst4wave's saved state.
    bool gw_first = true, lshape = true;
    int gw_nsps0 = -99;
    float gw_twopi = 0, gw_dt = 0, gw_tsym = 0;
    std::vector<cfloat> ctab;
    std::vector<float> pulse;

    void gen_fst4wave(const int* itone, int nsym, int nsps, int nwave, float fsample, int hmod, float f0, cfloat* cwave);
    template <class Sample>
    void dopspread(const int* itone, const Sample* iwave, int nsps, int nmax, int ndown, int hmod, int i0, float fc,
                   float& fmid, float& w50, std::vector<std::pair<float, float>>* spectrum, bool fast = false,
                   bool count = true);
    template <class Sample>
    void decode(Sample* iwave, int nutc, const Options& o, const std::function<void(const Result&)>& on_result,
                bool* calls_changed);
};

void Decoder::Impl::gen_fst4wave(const int* itone, int nsym, int nsps, int nwave, float fsample, int hmod, float f0,
                                 cfloat* cwave) {
    if (gw_first) {
        gw_twopi = 8.0f * std::atan(1.0f);
        ctab.resize(NTAB);
        for (int i = 0; i < NTAB; i += 4) {  // cos and sin by libmvec, four at a time
            float phi[4], c[4], s[4];
            for (int l = 0; l < 4; ++l) phi[l] = static_cast<float>(i + l) * gw_twopi / static_cast<float>(NTAB);
            simd::cosf4(phi, c);
            simd::sinf4(phi, s);
            for (int l = 0; l < 4; ++l) ctab[i + l] = cfloat(c[l], s[l]);
        }
        // gfortran's -O3 folded the loop's first and last three iterations
        // at compile time (correctly rounded, not libmvec's): its constants.
        auto fb = [](uint32_t u) {
            float f;
            std::memcpy(&f, &u, 4);
            return f;
        };
        ctab[0] = cfloat(1.0f, 0.0f);
        ctab[65533] = cfloat(fb(0x3f7fffffu), fb(0xb996a888u));
        ctab[65534] = cfloat(fb(0x3f800000u), fb(0xb948d111u));
        ctab[65535] = cfloat(fb(0x3f800000u), fb(0xb8c8a221u));
        const char* v = std::getenv("FST4_NOSHAPING");
        if (v && std::strlen(v) == 1 && v[0] == '1') lshape = false;
    }
    if (gw_first || nsps != gw_nsps0) {
        pulse.assign(3 * static_cast<size_t>(nsps), 0.0f);
        gw_dt = 1.0f / fsample;
        gw_tsym = static_cast<float>(nsps) / fsample;
        // Compute the smoothed frequency-deviation pulse
        for (int i = 1; i <= 3 * nsps; ++i) {
            const float tt = (static_cast<float>(i) - 1.5f * static_cast<float>(nsps)) / static_cast<float>(nsps);
            pulse[i - 1] = gfsk_pulse(2.0f, tt);
        }
        gw_first = false;
        gw_nsps0 = nsps;
    }
    const float twopi = gw_twopi;
    // Generate one symbol at a time, retaining phase across symbol boundaries.
    std::vector<float> dphi(static_cast<size_t>(nsps));
    const float dphi_peak = twopi * static_cast<float>(hmod) / static_cast<float>(nsps);
    std::vector<float> scaled_pulse(3 * static_cast<size_t>(nsps));
    for (size_t i = 0; i < scaled_pulse.size(); ++i) scaled_pulse[i] = dphi_peak * pulse[i];
    const float carrier = twopi * (f0 - 1.5f * static_cast<float>(hmod) / gw_tsym) * gw_dt;
    for (int i = nsym * nsps; i < nwave; ++i) cwave[i] = 0.0f;
    float phi = 0.0f;
    int k = 0;
    for (int iblock = 1; iblock <= nsym; ++iblock) {
        std::fill(dphi.begin(), dphi.end(), 0.0f);
        // Accumulate overlapping pulses in tone order to preserve phase rounding.
        for (int j = std::max(1, iblock - 1); j <= std::min(nsym, iblock + 1); ++j) {
            const int ip = (iblock - j + 1) * nsps + 1;
            const float t = static_cast<float>(itone[j - 1]);
            for (int m = 0; m < nsps; ++m) dphi[m] = dphi[m] + scaled_pulse[ip - 1 + m] * t;
        }
        for (int m = 0; m < nsps; ++m) dphi[m] = dphi[m] + carrier;
        for (int j = 0; j < nsps; ++j) {
            ++k;
            int i = static_cast<int>(phi * static_cast<float>(NTAB) / twopi);
            i = i & (NTAB - 1);
            cwave[k - 1] = ctab[i];
            phi = phi + dphi[j];
            if (phi > twopi) phi = phi - twopi;
        }
    }
    // Compute the ramp-up and ramp-down symbols
    if (lshape) {
        const float den = static_cast<float>(nsps / 2);
        ramp(cwave, nsps / 4, twopi, den, true);
        const int k1 = (nsym - 1) * nsps + 3 * nsps / 4 + 1;
        ramp(cwave + (k1 - 1), nsps / 4 + 1, twopi, den, false);
    }
}

template <class Sample>
void Decoder::Impl::dopspread(const int* itone, const Sample* iwave, int nsps, int nmax, int ndown, int hmod, int i0,
                              float fc, float& fmid, float& w50, std::vector<std::pair<float, float>>* spectrum,
                              bool fast, bool count) {
    if (count) ++ncall;
    const int nfft = 2 * nmax;
    const int nwave = std::max(nmax, (NN + 2) * nsps);
    std::vector<cfloat> cwave(static_cast<size_t>(nwave));
    FftwVec<cfloat> g(fast ? 0 : static_cast<size_t>(nfft));
    FftwVec<cfloat> gd;
    int nwrap = nfft;
    const float fsample = 12000.0f;
    gen_fst4wave(itone, NN, nsps, nwave, fsample, hmod, fc, cwave.data());
    // cwave = cshift(cwave, -i0*ndown): element k of the shifted wave is
    // cwave[(k + sh) mod nwave]; read in place.
    long sh = -static_cast<long>(i0) * ndown;
    sh %= nwave;
    if (sh < 0) sh += nwave;
    auto cw = [&](int k) -> const cfloat& {
        long j = k + sh;
        if (j >= nwave) j -= nwave;
        return cwave[static_cast<size_t>(j)];
    };
    const float fac = 1.0f / 32768;
    if (!fast) {
        for (int k = 0; k < nmax; ++k) {
            const float s = fac * static_cast<float>(iwave[k]);
            g[k] = cfloat(s * cw(k).real(), s * -cw(k).imag());
        }
        for (int k = nmax; k < nfft; ++k) g[k] = 0.0f;
        fft_c2c_inplace(g.data(), nfft, -1);  // Forward c2c FFT
    } else {
        // Not jt9's: the same product low-passed (flat to 1e-5 within +/-10.2
        // Hz, 100 dB down past 89.8 Hz) and decimated by 120 to 100 Hz, its
        // FFT 120 times smaller with the same bins (nfft/120 points at
        // 100 Hz). Only the spectrum's shape is used (it is normalised by
        // its peak), so the filter's gain does not matter.
        const int D = 120, n2 = nfft / D;
        const std::vector<float>& h = spread_filter();
        const int half = static_cast<int>(h.size()) / 2;
        // The product, zero padded by the filter's half length each side.
        std::vector<float> gr(static_cast<size_t>(nmax + 2 * half + 16), 0.0f), gi(gr.size(), 0.0f);
        float* pr = gr.data() + half;
        float* pi = gi.data() + half;
        for (int k = 0; k < nmax; ++k) {
            const float s = fac * static_cast<float>(iwave[k]);
            pr[k] = s * cw(k).real();
            pi[k] = s * -cw(k).imag();
        }
        gd.assign(static_cast<size_t>(n2), cfloat(0.0f, 0.0f));
        const int L = static_cast<int>(h.size());
        const int m0 = -(half / D), m1 = (nmax - 1 + half) / D;
        for (int m = m0; m <= m1; ++m) {
            const int c = m * D;  // the output's centre sample
            const int lo = std::max(-half, c - half);
            const int hi = std::min(nmax - 1 + half, c + half);
            const float* hr = h.data() + (lo - (c - half));
            simd::V8 ar = simd::zero(), ai = simd::zero();
            int t = lo;
            for (; t + 7 <= hi; t += 8) {
                const simd::V8 w = simd::load(hr + (t - lo));
                ar = simd::fma_loose(w, simd::load(pr + t), ar);
                ai = simd::fma_loose(w, simd::load(pi + t), ai);
            }
            float sr = simd::hsum(ar), si = simd::hsum(ai);
            for (; t <= hi; ++t) {
                sr += hr[t - lo] * pr[t];
                si += hr[t - lo] * pi[t];
            }
            (void)L;
            gd[static_cast<size_t>((m % n2 + n2) % n2)] = cfloat(sr, si);
        }
        fft_c2c_inplace(gd.data(), n2, -1);
        nwrap = n2;
    }

    const float df = 12000.0f / static_cast<float>(nfft);
    const cfloat* spec = fast ? gd.data() : g.data();
    auto pw = [&](int i) {
        const int j = i < 0 ? i + nwrap : i;
        return spec[j].real() * spec[j].real() + spec[j].imag() * spec[j].imag();
    };
    int ia = static_cast<int>(1.0f / df);
    float smax = 0.0f;
    for (int i = -ia; i <= ia; ++i) smax = fmax_g(pw(i), smax);  // smax in +/- 1 Hz around 0

    ia = static_cast<int>(10.1f / df);
    std::vector<float> ssv(2 * static_cast<size_t>(ia) + 1);
    float* ss = ssv.data() + ia;  // ss(-ia:ia)
    float sum1 = 0.0f, sum2 = 0.0f;
    int nns = 0;
    for (int i = -ia; i <= ia; ++i) {
        ss[i] = pw(i) / smax;
        const float f = static_cast<float>(i) * df;
        if (f >= -4.0f && f <= -2.0f) {
            sum1 = sum1 + ss[i];  // Power between -2 and -4 Hz
            ++nns;
        } else if (f >= 2.0f && f <= 4.0f) {
            sum2 = sum2 + ss[i];  // Power between +2 and +4 Hz
        }
    }
    const float avg = fmin_g(sum1 / static_cast<float>(nns), sum2 / static_cast<float>(nns));

    sum1 = 0.0f;
    for (int i = -ia; i <= ia; ++i) {
        const float f = static_cast<float>(i) * df;
        if (std::fabs(f) <= 1.0f) sum1 = sum1 + ss[i] - avg;  // Power in abs(f) < 1 Hz
    }

    ia = static_cast<int>(std::lround(1.0f / df)) + 1;
    sum2 = 0.0f;
    float xi1 = -999, xi2 = -999, xi3 = -999, sum2z = 0.0f;
    for (int i = -ia; i <= ia; ++i) {  // Find freq range that has 50% of signal power
        sum2 = sum2 + ss[i] - avg;
        if (sum2 >= 0.25f * sum1 && xi1 == -999.0f) xi1 = static_cast<float>(i - 1) + (0.25f * sum1 - sum2) / (sum2 - sum2z);
        if (sum2 >= 0.50f * sum1 && xi2 == -999.0f) xi2 = static_cast<float>(i - 1) + (sum2 - 0.50f * sum1) / (sum2 - sum2z);
        if (sum2 >= 0.75f * sum1) {
            xi3 = static_cast<float>(i - 1) + (sum2 - 0.75f * sum1) / (sum2 - sum2z);
            break;
        }
        sum2z = sum2;
    }
    const float xdiff = std::sqrt(1.0f + (xi3 - xi1) * (xi3 - xi1));  // Keep small values from fluctuating too widely
    w50 = xdiff * df;  // Doppler spread
    fmid = xi2 * df;   // Frequency midpoint of signal power

    if (spectrum) {  // the spectrum for plotting (fort.52)
        spectrum->clear();
        const int nx = static_cast<int>(std::lround(xi2));
        for (int i = -ia; i <= ia; ++i) {
            float y = static_cast<float>(ncall - 1);
            const int j = i + nx;
            if (std::fabs(static_cast<float>(j) * df) < 10.0f) y = 0.99f * ss[i + nx] + static_cast<float>(ncall) - 1;
            spectrum->emplace_back(static_cast<float>(i) * df, y);
        }
    }
}

template <class Sample>
void Decoder::Impl::decode(Sample* iwave, int nutc, const Options& o,
                           const std::function<void(const Result&)>& on_result, bool* calls_changed) {
    if (calls_changed) *calls_changed = false;
    stopped = false;
    if (first) {
        if (have_pending) {  // fst4w_calls.txt: up to MAXWCALLS lines, to the first empty one
            for (const std::string& l : pending_calls) {
                if (static_cast<int>(wcalls.size()) >= MAXWCALLS) break;
                std::string w = l.substr(0, std::min<size_t>(20, l.size()));
                size_t k = 0;
                while (k < w.size() && w[k] == ' ') ++k;
                w = rtrim(w.substr(k));
                if (w.empty()) break;
                wcalls.push_back(w);
            }
        }
        first = false;
    }
    unpacker.set_calls(o.mycall, o.hiscall);

    const Period P = Period::of(o.ntrperiod);
    const int hmod = 1;
    const int nsps = P.nsps, nmax = P.nmax, ndown = P.ndown, nfft1 = P.nfft1, nfft2 = P.nfft2, nss = P.nss;
    const float fs2 = P.fs2, baud = P.baud, sigbw = P.sigbw;
    const int nspsec = P.nspsec;

    FftwVec<float> c_bigfft(2 * static_cast<size_t>(nfft1 / 2 + 1), 0.0f);
    FftwVec<cfloat> c2(static_cast<size_t>(nfft2));
    std::vector<cfloat> cframe(160 * static_cast<size_t>(nss));

    const int ndepth = o.ndepth & 3;
    int jittermax = 2;
    bool do_k50_decode = false;
    if (ndepth == 3) {
        jittermax = 2;
        do_k50_decode = true;
    } else if (ndepth == 2) {
        jittermax = 2;
    } else if (ndepth == 1) {
        jittermax = 0;
    }

    // Noise blanker setup
    const int ndropmax = 1;
    int npct = 0;
    const int nb = o.nexp_decode / 256 - 3;
    if (nb >= 0) npct = nb;
    int inb1 = 20, inb2 = 5;
    if (nb == -1) inb2 = 5;
    else if (nb == -2) inb2 = 2;
    else if (nb == -3) inb2 = 1;
    else inb1 = 0;

    // FST4W: noise baseline window and signal search window.
    const int nfa = std::max(100, o.nfqso - o.ntol - 100);
    const int nfb = std::min(4800, o.nfqso + o.ntol + 100);
    const float fa = static_cast<float>(std::max(
        100, static_cast<int>(std::lround(static_cast<float>(o.nfqso) + 1.5f * baud - static_cast<float>(o.ntol)))));
    const float fb = static_cast<float>(std::min(
        4800, static_cast<int>(std::lround(static_cast<float>(o.nfqso) + 1.5f * baud + static_cast<float>(o.ntol)))));

    std::vector<std::string> decodes;  // decodes(100)
    bool new_callsign = false;
    std::vector<float> candidates0(200 * 5), bitmetrics(320 * 4), s4(4 * NN), llrs(240 * 4), llr(240);
    auto C = [&](int i, int j) -> float& { return candidates0[(j - 1) * 200 + (i - 1)]; };
    auto BM = [&](int i, int j) { return bitmetrics[(j - 1) * 320 + (i - 1)]; };
    auto LL = [&](int i, int j) -> float& { return llrs[(j - 1) * 240 + (i - 1)]; };
    auto S4 = [&](int itone, int k) { return s4[(k - 1) * 4 + itone]; };
    int8_t apmask[240] = {}, cw[240] = {}, message74[74] = {};
    int ncand = 0;
    float w50 = w50_saved;

    for (int inb = 0; inb <= inb1 && !halt(); inb += inb2) {
        if (nb < 0) npct = inb;  // looping over blanker settings
        blanker(iwave, nfft1, ndropmax, npct, c_bigfft.data());
        // The big fft: the smoothed spectrum, and downconverting each candidate.
        fft_r2c_inplace(c_bigfft.data(), nfft1);
        const int nsyncoh = 8;
        float minsync = 1.20f;
        if (o.ntrperiod == 15) minsync = 1.15f;

        // First approximation of candidate frequencies
        if (!get_candidates_fst4(c_bigfft.data(), nfft1, nsps, hmod, P.fs, fa, fb, nfa, nfb, minsync, ncand,
                                 candidates0.data()))
            ncand = 0;  // the reference's ncand is then whatever it held: baseline <= 0 cannot follow its RMS check
        int isbest = 0;
        for (int icand = 1; icand <= ncand && !halt(); ++icand) {
            const float fc0 = C(icand, 1);
            fst4_downsample(c_bigfft.data(), nfft1, ndown, fc0, sigbw, c2.data());
            float sbest = 0.0f;
            fe.fst4_sync_search(c2.data(), nfft2, hmod, fs2, nss, o.ntrperiod, nsyncoh, o.emedelay, sbest, fcbest,
                                isbest);
            C(icand, 3) = fc0 + fcbest;
            C(icand, 4) = static_cast<float>(isbest);
        }

        // remove duplicate candidates
        for (int icand = 1; icand <= ncand; ++icand) {
            const float fc = C(icand, 3);
            isbest = static_cast<int>(std::lround(C(icand, 4)));
            for (int ic2 = icand + 1; ic2 <= ncand; ++ic2) {
                const float fc2 = C(ic2, 3);
                const int isbest2 = static_cast<int>(std::lround(C(ic2, 4)));
                if (fc2 > 0.0f && std::fabs(fc2 - fc) < 0.10f * baud && std::abs(isbest2 - isbest) <= 2)
                    C(ic2, 3) = -1;
            }
        }
        int ic = 0;
        for (int icand = 1; icand <= ncand; ++icand) {
            if (C(icand, 3) > 0) {
                ++ic;
                for (int j = 1; j <= 5; ++j) C(ic, j) = C(icand, j);
            }
        }
        ncand = ic;
        const std::vector<float>& candidates = candidates0;  // FST4W: candidates = candidates0
        auto CC = [&](int i, int j) { return candidates[(j - 1) * 200 + (i - 1)]; };

        for (int icand = 1; icand <= ncand && !halt(); ++icand) {
            const float csync = CC(icand, 2);
            const float fc_synced = CC(icand, 3);
            isbest = static_cast<int>(std::lround(CC(icand, 4)));
            float xdt = static_cast<float>(isbest - nspsec) / fs2;
            if (o.ntrperiod == 15) xdt = (static_cast<float>(isbest) - static_cast<float>(nspsec) / 2.0f) / fs2;
            fst4_downsample(c_bigfft.data(), nfft1, ndown, fc_synced, sigbw, c2.data());

            bool next_candidate = false;
            for (int ijitter = 0; ijitter <= jittermax && !next_candidate; ++ijitter) {
                const int ioffset = ijitter == 0 ? 0 : ijitter == 1 ? 1 : -1;
                const int is0 = isbest + ioffset;
                const int iend = is0 + 160 * nss - 1;
                if (is0 < 0 || iend > nfft2 - 1) continue;
                std::copy(c2.begin() + is0, c2.begin() + iend + 1, cframe.begin());
                std::fill(bitmetrics.begin(), bitmetrics.end(), 0.0f);
                int nsync_qual = 0;
                bool badsync = false;
                fe.get_fst4_bitmetrics(cframe.data(), nss, bitmetrics.data(), s4.data(), nsync_qual, badsync);
                if (badsync) continue;

                for (int il = 1; il <= 4; ++il) {
                    for (int i = 0; i < 60; ++i) {
                        LL(1 + i, il) = BM(17 + i, il);
                        LL(61 + i, il) = BM(93 + i, il);
                        LL(121 + i, il) = BM(169 + i, il);
                        LL(181 + i, il) = BM(245 + i, il);
                    }
                }
                const int nblock = 4, ntmax = 4;  // 50-bit messages, no ap decoding
                for (int itry = 1; itry <= ntmax; ++itry) {
                    for (int i = 1; i <= 240; ++i) llr[i - 1] = LL(i, itry);
                    int iaptype = 0;
                    (void)nblock;
                    float dmin = 0.0f;
                    int nharderrors = -1, ntype = 0, keff = 66;
                    bool unpk77_success = false;
                    std::string msg = pad37("");
                    bool skip66 = false;
                    // Try decoding with Keff=66
                    ldpc.decode240_74(llr.data(), 66, 2, 3, apmask, message74, cw, ntype, nharderrors, dmin);
                    if (nharderrors < 0) {
                        skip66 = true;
                    } else {
                        int ones = 0;
                        for (int i = 0; i < 240; ++i) ones += cw[i] == 1;
                        if (ones == 0) {
                            nharderrors = -nharderrors;
                            continue;
                        }
                        const auto r = unpacker.unpack(reinterpret_cast<const uint8_t*>(message74));
                        msg = r.msg;
                        unpk77_success = r.ok;
                        if (o.lprinthash22 && unpk77_success && findex(msg, "<...>") > 0) {
                            int n22tmp = 0;
                            for (int i = 0; i < 22; ++i) n22tmp = 2 * n22tmp + message74[i];
                            const int i1 = findex(msg, " ");
                            const std::string wpart = rtrim(msg.substr(static_cast<size_t>(i1)).substr(0, 20));
                            char b[16];
                            std::snprintf(b, sizeof b, "<%07d>", n22tmp);
                            msg = pad37(std::string(b) + " " + wpart);
                        }
                        if (unpk77_success && do_k50_decode) {
                            // Keff=66 decode: save call/grid in fst4w_calls.txt if not there already.
                            const int i1 = findex(msg, " ");
                            const int i2 = i1 + findex(msg.substr(static_cast<size_t>(i1)), " ");
                            const std::string wpart = rtrim(msg.substr(0, static_cast<size_t>(i2)).substr(0, 20));
                            // Only save callsigns/grids from type 1 messages
                            if (wpart.find('/') == std::string::npos && wpart.find('<') == std::string::npos) {
                                bool ifound = false;
                                for (const std::string& w : wcalls)
                                    if (w == wpart) ifound = true;  // index(wcalls(i), wpart) of character*20s
                                if (!ifound) {  // This is a new callsign
                                    new_callsign = true;
                                    if (static_cast<int>(wcalls.size()) < MAXWCALLS) {
                                        wcalls.push_back(wpart);
                                    } else {
                                        wcalls.erase(wcalls.begin());
                                        wcalls.push_back(wpart);
                                    }
                                }
                            }
                        }
                    }
                    (void)skip66;
                    // If no decode then try Keff=50
                    iaptype = 0;
                    if (!unpk77_success && do_k50_decode) {
                        keff = 50;
                        ldpc.decode240_74(llr.data(), 50, 1, 4, apmask, message74, cw, ntype, nharderrors, dmin);
                        int ones = 0;
                        for (int i = 0; i < 240; ++i) ones += cw[i] == 1;
                        if (ones == 0) {
                            nharderrors = -nharderrors;
                            continue;
                        }
                        const auto r = unpacker.unpack(reinterpret_cast<const uint8_t*>(message74));
                        msg = r.msg;
                        unpk77_success = r.ok;
                        // No CRC in this mode: only call/grids seen before
                        if (unpk77_success) {
                            unpk77_success = false;
                            for (const std::string& w : wcalls)
                                if (findex(msg, w) > 0) unpk77_success = true;
                        }
                    }

                    if (nharderrors >= 0 && unpk77_success) {
                        bool dupe = false;
                        for (const std::string& d : decodes)
                            if (d == msg) dupe = true;
                        if (dupe) {
                            next_candidate = true;
                            break;
                        }
                        decodes.push_back(msg);

                        int itone[NN];
                        tones_from_bits(message74, itone);
                        float fmid = -999.0f;
                        Result res;
                        if (o.plotspec && o.spread != Options::Spread::Off) {
                            const bool fast = o.spread == Options::Spread::Fast;
                            dopspread(itone, iwave, nsps, nmax, ndown, hmod, isbest, fc_synced, fmid, w50,
                                      &res.spectrum, fast);
                            if (o.spread == Options::Spread::Both)
                                dopspread(itone, iwave, nsps, nmax, ndown, hmod, isbest, fc_synced, res.fmid_fast,
                                          res.w50_fast, nullptr, true, false);
                        }
                        float xsig = 0.0f;
                        for (int i = 1; i <= NN; ++i) xsig = xsig + S4(itone[i - 1], i);
                        const float base = CC(icand, 5);
                        float snr_calfac;
                        switch (o.ntrperiod) {
                        case 15: snr_calfac = 800.0f; break;
                        case 30: snr_calfac = 600.0f; break;
                        case 60: snr_calfac = 430.0f; break;
                        case 120: snr_calfac = 390.0f; break;
                        case 300: snr_calfac = 340.0f; break;
                        case 900: snr_calfac = 320.0f; break;
                        case 1800: snr_calfac = 320.0f; break;
                        default: snr_calfac = 430.0f;
                        }
                        const float arg = snr_calfac * xsig / base - 1.0f;
                        float xsnr;
                        if (arg > 0.0f) {
                            // 10*log10(1.46/2500) is folded at compile time; computed
                            // here the same way, from the float quotient.
                            static const float c1 = 10.0f * static_cast<float>(std::log10(static_cast<double>(1.46f / 2500.0f)));
                            xsnr = 10.0f * std::log10(arg) + c1 + 10.0f * std::log10(8200.0f / static_cast<float>(nsps));
                        } else {
                            xsnr = -99.9f;
                        }
                        res.nutc = nutc;
                        res.sync = 0.0f;
                        res.nsnr = static_cast<int>(std::lround(xsnr));
                        res.xdt = xdt;
                        res.fsig = fc_synced - 1.5f * baud;
                        res.msg = msg;
                        res.nap = iaptype;
                        res.qual = 0.0f;
                        res.ntrperiod = o.ntrperiod;
                        res.fmid = fmid;
                        res.w50 = w50;
                        if (o.decdata) {
                            int nhp = 0;
                            float hd = 0.0f;
                            for (int i = 1; i <= 240; ++i) {
                                const int hdec = LL(i, 1) >= 0.0f ? 1 : 0;
                                const int x = hdec ^ cw[i - 1];
                                nhp += x != 0;  // hard errors wrt N=1 soft symbols
                                hd = hd + static_cast<float>(x) * std::fabs(LL(i, 1));
                            }
                            res.icand = icand;
                            res.itry = itry;
                            res.nsyncoh = nsyncoh;
                            res.iaptype = iaptype;
                            res.ijitter = ijitter;
                            res.npct = npct;
                            res.ntype = ntype;
                            res.keff = keff;
                            res.nsync_qual = nsync_qual;
                            res.nharderrors = nharderrors;
                            res.dmin = dmin;
                            res.nhp = nhp;
                            res.hd = hd;
                            res.csync = csync;
                            res.xsnr = xsnr;
                        }
                        on_result(res);
                        next_candidate = true;
                        break;
                    }
                }  // metrics
            }  // istart jitter
        }  // candidate list
    }  // noise blanker loop
    w50_saved = w50;
    if (new_callsign && do_k50_decode && calls_changed) *calls_changed = true;
}

Decoder::Decoder() : d_(new Impl) {}
Decoder::~Decoder() = default;

void Decoder::load_calls(const std::vector<std::string>& lines) {
    d_->pending_calls = lines;
    d_->have_pending = true;
}

std::vector<std::string> Decoder::calls() const { return d_->wcalls; }

std::string Decoder::save_hashes() const { return d_->unpacker.save(); }
bool Decoder::load_hashes(const std::string& text) { return d_->unpacker.load(text); }
void Decoder::set_stop(const std::atomic<bool>* stop) { d_->stop = stop; }
bool Decoder::stopped() const { return d_->stopped; }

size_t Decoder::memory_needed(const Options& o) {
    // decode()'s own buffers at their most: the big FFT's, the downsampled
    // period's and a frame's, and dopspread's (made per decode, beside
    // them): its waveform and the exact spread's FFT of twice the period,
    // or the fast one's filtered product. The caller's samples apart.
    const Period P = Period::of(o.ntrperiod);
    const size_t c = sizeof(cfloat), f = sizeof(float);
    size_t n = f * 2 * (static_cast<size_t>(P.nfft1) / 2 + 1) + c * static_cast<size_t>(P.nfft2) + c * 160 * static_cast<size_t>(P.nss);
    if (o.plotspec && o.spread != Options::Spread::Off) {
        n += c * std::max<size_t>(P.nmax, (NN + 2) * static_cast<size_t>(P.nsps));
        const size_t exact = c * 2 * static_cast<size_t>(P.nmax);
        const size_t fast = 2 * f * (static_cast<size_t>(P.nmax) + 1216) + c * 2 * static_cast<size_t>(P.nmax) / 120;
        n += o.spread == Options::Spread::Exact ? exact : o.spread == Options::Spread::Fast ? fast : exact + fast;
    }
    return n;
}

void Decoder::decode(int16_t* iwave, int nutc, const Options& o, const std::function<void(const Result&)>& on_result,
                     bool* calls_changed) {
    d_->decode(iwave, nutc, o, on_result, calls_changed);
}

void Decoder::decode(float* iwave, int nutc, const Options& o, const std::function<void(const Result&)>& on_result,
                     bool* calls_changed) {
    d_->decode(iwave, nutc, o, on_result, calls_changed);
}

}  // namespace fst4w
