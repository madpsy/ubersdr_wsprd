// FT4's front end (ft4_frontend.h): lib/ft4/gen_ft4wave.f90,
// subtractft4.f90 and genft4.f90 (with its entry
// get_ft4_tones_from_77bits) at 567ad29ce, line for line.
//
// As the Release object computes them (from its disassembly):
// - gen_ft4wave: the pulse made on the first call only (SAVEd, from that
//   call's nsps), tt = (real(i) - 1.5*real(nsps))/real(nsps);
//   dphi_peak = (twopi*hmod)/real(nsps); each symbol adds
//   (dphi_peak*pulse(k))*real(itone(j)); then dphi + (twopi*f0)*dt. Each
//   sample by scalar cosf and sinf (not sincosf; sin(0) of the first
//   folded, the same +0), phi = mod(phi + dphi(j), twopi) by x87 fprem on
//   the single sum (exactly fmodf). The ramps' cos by libmvec four at a
//   time from i = 0, the nsps mod 4 left by scalar cosf, argument
//   (twopi*real(i))/(2.0*nsps); cwave*(1-+cos)/2.0 as complex*real (the
//   full product with a zero imaginary part) then complex/real by Smith's
//   rule (ratio 0, /2.0 folded to *0.5).
// - subtractft4: nstart = int((dt*12000.0 + 1.0) - 576.0) in single;
//   camp = dd*conjg(cref) as complex*real; the window cos(pi*j/1400)**2
//   by libmvec for j = -699..700 four at a time, window(-700) folded
//   (0x2709ADED, FT8's bits), sum in index order; cw = window/sum,
//   cshift'ed, its FFT times fac (folded 0x37672AE4) as complex*real; the
//   last subtraction dd - 2*real(cfilt*cref) with the doubling as x + x.
// - genft4: the NUL cut and the blank strip (strips 1+2+... characters,
//   as there), pack77_legacy_truncating_fallback, the i3/n3 check,
//   unpack77; the entry scrambles the caller's msgbits with rvec in place.
//
// Joe Taylor, K1JT, Steven Franke, K9AN, and the WSJT-X developers (GNU GPL
// v3), and so is this.
#include "ft4_frontend.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "fst4w_frontend.h"  // FftwVec, fft_c2c_inplace (four2a's plans)
#include "ft8_fortran.h"     // the Fortran helpers, FT8's
#include "ft8_frontend.h"    // gfsk_pulse
#include "ft8_ldpc.h"        // encode174_91
#include "ft8_maths.h"       // glibc 2.39's maths and x86-64 libmvec, carried
#include "ft_fast_subtract.h"  // --fast-subtract's smoothing
#include "packjt77.h"

namespace ft4 {

using namespace ft8::fortran;
using ft8::cosf4_glibc;
using ft8::cosf_glibc;
using ft8::sinf_glibc;

// This group's SAVEd state: gen_ft4wave's SAVE and subtractft4's
// /heap4/ common and first.
struct Frontend::Wave {
    // gen_ft4wave
    bool gw_first = true;
    float twopi = 0.0f, dt = 0.0f, hmod = 0.0f;
    std::vector<float> pulse = std::vector<float>(6912, 0.0f);  // pulse(6912), static: zeros
    std::vector<float> pp = std::vector<float>(6912, 0.0f);     // dphi_peak*pulse
    std::vector<float> dphi;                                    // allocatable, (nsym+2)*nsps
    // subtractft4: /heap4/ cref(NFRAME), cfilt(NMAX), cw(NMAX),
    // xjunk(NFRAME), made on the first call (camp(NMAX) made in cfilt)
    bool st_first = true;
    fst4w::FftwVec<cfloat> cref, cfilt, cw;
    std::vector<float> xjunk;
    int al_cw = -1, al_cfilt = -1;  // the trace's '@' lines (loc() mod 16)
    std::unique_ptr<ft8::FastSubtract> fs;  // fast subtraction's (subtractft4_fast), made on its first call
};

Frontend::Wave* Frontend::make_wave() { return new Wave(); }
void Frontend::del_wave(Wave* p) { delete p; }
size_t Frontend::wave_bytes(const Wave* w) {
    if (!w) return 0;
    return (w->pulse.capacity() + w->pp.capacity()) * sizeof(float) + w->dphi.capacity() * sizeof(float) +
           (w->cref.capacity() + w->cfilt.capacity() + w->cw.capacity()) * sizeof(cfloat) +
           w->xjunk.capacity() * sizeof(float) + (w->fs ? w->fs->bytes() : 0);
}

namespace {

// genft4.f90's DATA
const int8_t rvec[77] = {0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0,
                         1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1,
                         1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1};
const int icos4a[4] = {0, 1, 3, 2};
const int icos4b[4] = {1, 0, 2, 3};
const int icos4c[4] = {2, 3, 1, 0};
const int icos4d[4] = {3, 2, 0, 1};

// cwave(k) = cwave(k)*w/2.0 as the object has it: complex*real, the full
// product with the real's zero imaginary part, then complex/real by
// Smith's rule with ratio 0 (/2.0 as *0.5).
inline cfloat ramp_c(cfloat a, float w) {
    const float pr = w * a.real() - a.imag() * 0.0f;
    const float pi = w * a.imag() + a.real() * 0.0f;
    return cfloat((pr + pi * 0.0f) * 0.5f, (pi - pr * 0.0f) * 0.5f);
}

// The ramps' cos(twopi*i/(2.0*nsps)), i = 0..nsps-1: libmvec four at a
// time from i = 0, the nsps mod 4 left by scalar cosf.
void ramp_cos(float twopi, int nsps, std::vector<float>& c) {
    c.resize(static_cast<size_t>(std::max(nsps, 0)));
    const float fnsps = static_cast<float>(nsps);
    const float den = fnsps + fnsps;  // 2.0*nsps
    int i = 0;
    if (nsps > 3) {
        for (; i + 3 < nsps; i += 4) {
            float a[4];
            for (int q = 0; q < 4; ++q) a[q] = (static_cast<float>(i + q) * twopi) / den;
            cosf4_glibc(a, &c[i]);
        }
    }
    for (; i < nsps; ++i) c[i] = cosf_glibc((static_cast<float>(i) * twopi) / den);
}

}  // namespace

// genft4.f90's entry get_ft4_tones_from_77bits(msgbits, i4tone), from
// label 2.
void get_ft4_tones_from_77bits(int8_t msgbits[77], int i4tone[NN]) {
    for (int i = 0; i < 77; ++i) msgbits[i] = static_cast<int8_t>((msgbits[i] + rvec[i]) % 2);
    int8_t codeword[2 * ND];
    ft8::encode174_91(msgbits, codeword);
    // Grayscale mapping: 00 0, 01 1, 11 2, 10 3
    int itmp[ND];
    for (int i = 1; i <= ND; ++i) {
        const int is = codeword[2 * i - 1] + 2 * codeword[2 * i - 2];
        if (is <= 1) itmp[i - 1] = is;
        if (is == 2) itmp[i - 1] = 3;
        if (is == 3) itmp[i - 1] = 2;
    }
    for (int i = 0; i < 4; ++i) {
        i4tone[i] = icos4a[i];
        i4tone[33 + i] = icos4b[i];
        i4tone[66 + i] = icos4c[i];
        i4tone[99 + i] = icos4d[i];
    }
    for (int i = 0; i < 29; ++i) {
        i4tone[4 + i] = itmp[i];
        i4tone[37 + i] = itmp[29 + i];
        i4tone[70 + i] = itmp[58 + i];
    }
}

void genft4(packjt77::Packjt77& pk, const std::string& msg0, int ichk, std::string& msgsent, int8_t msgbits[77],
            int i4tone[NN]) {
    std::string message = fix(msg0, 37);
    for (int i = 1; i <= 37; ++i) {
        if (message[i - 1] == '\0') {
            for (int k = i - 1; k < 37; ++k) message[k] = ' ';
            break;
        }
    }
    for (int i = 1; i <= 37; ++i) {  // Strip leading blanks
        if (message[0] != ' ') break;
        message = fix(message.substr(static_cast<size_t>(i)), 37);
    }
    bool bad = true;
    const packjt77::Packjt77::Packed packed = pk.pack77_legacy_truncating_fallback(message);
    if (packed.i3 >= 0 && packed.n3 >= 0) {
        const packjt77::Packjt77::Unpacked u = pk.unpack77(packed.c77.data(), 0);
        msgsent = fix(u.msg, 37);
        if (ichk == 1) return;
        // read(c77,'(77i1)',err=1) msgbits: blanks read as 0, any other
        // non-digit an error.
        bool err = false;
        int8_t b[77];
        for (int i = 0; i < 77; ++i) {
            const char c = packed.c77[static_cast<size_t>(i)];
            if (c == ' ') {
                b[i] = 0;
            } else if (c >= '0' && c <= '9') {
                b[i] = static_cast<int8_t>(c - '0');
            } else {
                err = true;
                break;
            }
        }
        if (!err) {
            std::memcpy(msgbits, b, 77);
            if (u.ok) bad = false;
        }
    }
    if (bad) {
        std::memset(msgbits, 0, 77);
        std::fill(i4tone, i4tone + NN, 0);
        msgsent = "*** bad message ***                  ";
        return;
    }
    get_ft4_tones_from_77bits(msgbits, i4tone);
}

void Frontend::gen_ft4wave(const int* itone, int nsym, int nsps, float fsample, float f0, cfloat* cwave, float* wave,
                           int icmplx, int nwave) {
    Wave& w = *wave_;
    const int n = (nsym + 2) * nsps;
    w.dphi.assign(static_cast<size_t>(std::max(n, 0)), 0.0f);
    float* dphi = w.dphi.data();

    if (w.gw_first) {
        w.twopi = TWOPI;  // 8.0*atan(1.0), folded
        w.dt = 1.0f / fsample;
        w.hmod = 1.0f;
        // Compute the smoothed frequency-deviation pulse
        for (int i = 1; i <= 3 * nsps; ++i) {
            const float tt = (static_cast<float>(i) - 1.5f * static_cast<float>(nsps)) / static_cast<float>(nsps);
            w.pulse[i - 1] = ft8::gfsk_pulse(1.0f, tt);
        }
        if (trace_) trace_("Q " + X(h_r(w.pulse.data(), static_cast<size_t>(3 * nsps))));
        w.gw_first = false;
    }

    // Compute the smoothed frequency waveform.
    // Length = (nsym+2)*nsps samples, zero-padded
    const float dphi_peak = (w.twopi * w.hmod) / static_cast<float>(nsps);
    // dphi_peak*pulse(m) once a call (the same product each symbol
    // makes), and, when every product is +0 or more and finite, tone 0's
    // symbols skipped: their terms are then +0, and dphi, from +0 up by
    // such terms, is never -0, so adding +0 changes nothing.
    float* pp = w.pp.data();
    bool skip0 = true;
    for (int m = 0; m < 3 * nsps; ++m) {
        pp[m] = dphi_peak * w.pulse[m];
        skip0 = skip0 && !std::signbit(pp[m]) && std::isfinite(pp[m]);
    }
    for (int j = 1; j <= nsym; ++j) {
        const int ib = (j - 1) * nsps;
        if (skip0 && itone[j - 1] == 0) continue;
        const float t = static_cast<float>(itone[j - 1]);
        for (int m = 0; m < 3 * nsps; ++m) dphi[ib + m] = dphi[ib + m] + pp[m] * t;
    }

    // Calculate and insert the audio waveform
    uint32_t h0 = 0;
    if (trace_) h0 = h_r(dphi, static_cast<size_t>(n));
    float phi = 0.0f;
    const float shift = (w.twopi * f0) * w.dt;
    for (int j = 0; j < n; ++j) dphi[j] = dphi[j] + shift;
    if (trace_) trace_("q " + X(h0) + " " + X(h_r(dphi, static_cast<size_t>(n))));
    if (wave) std::fill(wave, wave + nwave, 0.0f);
    if (icmplx == 1) std::fill(cwave, cwave + nwave, cfloat(0.0f, 0.0f));
    // The complex wave (subtractft4's) faster, the same bits: cos
    // and sin of each phi by sincosf (glibc's sincosf, sinf and cosf share
    // the reduction and polynomials: equal on every float, [ft8-maths]), and mod(phi + dphi, twopi) without fprem: with phi and
    // dphi in [0, twopi) the sum s is below 2*twopi, so the remainder is s
    // or s - twopi, exact (Sterbenz) as fmod's; s rounded to 2*twopi
    // itself goes to fmod. Any dphi outside [0, twopi) (or NaN): the
    // reference's loop below.
    bool fast = icmplx == 1;
    for (int j = 0; fast && j < n; ++j) fast = dphi[j] >= 0.0f && dphi[j] < w.twopi;
    const float twopi2 = w.twopi + w.twopi;
    // The phases a block at a time, then their sines and cosines four at
    // a time (sincosf_n_glibc, the same bits).
    constexpr int BLK = 256;
    for (int j0 = 0; fast && j0 < n; j0 += BLK) {
        const int nb = std::min(BLK, n - j0);
        float ph[BLK], sn[BLK], cs[BLK];
        for (int j = 0; j < nb; ++j) {
            ph[j] = phi;
            const float t = phi + dphi[j0 + j];
            phi = t < w.twopi ? t : (t < twopi2 ? t - w.twopi : std::fmod(t, w.twopi));
        }
        ft8::sincosf_n_glibc(ph, sn, cs, nb);
        for (int j = 0; j < nb; ++j) cwave[j0 + j] = cfloat(cs[j], sn[j]);
    }
    int k = 0;
    for (int j = 0; !fast && j < n; ++j) {
        ++k;
        if (icmplx == 0) {
            wave[k - 1] = sinf_glibc(phi);
        } else {
            const float c = cosf_glibc(phi);
            const float s = sinf_glibc(phi);
            cwave[k - 1] = cfloat(c, s);
        }
        phi = std::fmod(phi + dphi[j], w.twopi);
    }

    // Compute the ramp-up and ramp-down symbols
    std::vector<float> c;
    ramp_cos(w.twopi, nsps, c);
    const int k1 = (nsym + 1) * nsps + 1;
    if (icmplx == 0) {
        for (int i = 0; i < nsps; ++i) wave[i] = (wave[i] * (1.0f - c[i])) * 0.5f;
        for (int i = 0; i < nsps; ++i) wave[k1 - 1 + i] = (wave[k1 - 1 + i] * (1.0f + c[i])) * 0.5f;
    } else {
        for (int i = 0; i < nsps; ++i) cwave[i] = ramp_c(cwave[i], 1.0f - c[i]);
        for (int i = 0; i < nsps; ++i) cwave[k1 - 1 + i] = ramp_c(cwave[k1 - 1 + i], 1.0f + c[i]);
    }
}

// gen_ft4wave in its plain form, line for line with the reference
// (each sample's cosf and sinf, fmod), kept to check gen_ft4wave against
// (tests/test_ft4_frontend.cpp, [ft4-fe-reference]).
void Frontend::gen_ft4wave_reference(const int* itone, int nsym, int nsps, float fsample, float f0, cfloat* cwave, float* wave,
                           int icmplx, int nwave) {
    Wave& w = *wave_;
    const int n = (nsym + 2) * nsps;
    w.dphi.assign(static_cast<size_t>(std::max(n, 0)), 0.0f);
    float* dphi = w.dphi.data();

    if (w.gw_first) {
        w.twopi = TWOPI;  // 8.0*atan(1.0), folded
        w.dt = 1.0f / fsample;
        w.hmod = 1.0f;
        // Compute the smoothed frequency-deviation pulse
        for (int i = 1; i <= 3 * nsps; ++i) {
            const float tt = (static_cast<float>(i) - 1.5f * static_cast<float>(nsps)) / static_cast<float>(nsps);
            w.pulse[i - 1] = ft8::gfsk_pulse(1.0f, tt);
        }
        if (trace_) trace_("Q " + X(h_r(w.pulse.data(), static_cast<size_t>(3 * nsps))));
        w.gw_first = false;
    }

    // Compute the smoothed frequency waveform.
    // Length = (nsym+2)*nsps samples, zero-padded
    const float dphi_peak = (w.twopi * w.hmod) / static_cast<float>(nsps);
    for (int j = 1; j <= nsym; ++j) {
        const int ib = (j - 1) * nsps;
        const float t = static_cast<float>(itone[j - 1]);
        for (int m = 0; m < 3 * nsps; ++m) dphi[ib + m] = dphi[ib + m] + (dphi_peak * w.pulse[m]) * t;
    }

    // Calculate and insert the audio waveform
    uint32_t h0 = 0;
    if (trace_) h0 = h_r(dphi, static_cast<size_t>(n));
    float phi = 0.0f;
    const float shift = (w.twopi * f0) * w.dt;
    for (int j = 0; j < n; ++j) dphi[j] = dphi[j] + shift;
    if (trace_) trace_("q " + X(h0) + " " + X(h_r(dphi, static_cast<size_t>(n))));
    if (wave) std::fill(wave, wave + nwave, 0.0f);
    if (icmplx == 1) std::fill(cwave, cwave + nwave, cfloat(0.0f, 0.0f));
    int k = 0;
    for (int j = 0; j < n; ++j) {
        ++k;
        if (icmplx == 0) {
            wave[k - 1] = sinf_glibc(phi);
        } else {
            const float c = cosf_glibc(phi);
            const float s = sinf_glibc(phi);
            cwave[k - 1] = cfloat(c, s);
        }
        phi = std::fmod(phi + dphi[j], w.twopi);
    }

    // Compute the ramp-up and ramp-down symbols
    std::vector<float> c;
    ramp_cos(w.twopi, nsps, c);
    const int k1 = (nsym + 1) * nsps + 1;
    if (icmplx == 0) {
        for (int i = 0; i < nsps; ++i) wave[i] = (wave[i] * (1.0f - c[i])) * 0.5f;
        for (int i = 0; i < nsps; ++i) wave[k1 - 1 + i] = (wave[k1 - 1 + i] * (1.0f + c[i])) * 0.5f;
    } else {
        for (int i = 0; i < nsps; ++i) cwave[i] = ramp_c(cwave[i], 1.0f - c[i]);
        for (int i = 0; i < nsps; ++i) cwave[k1 - 1 + i] = ramp_c(cwave[k1 - 1 + i], 1.0f + c[i]);
    }
}

void Frontend::subtractft4(float* dd, const int itone[NN], float f0, float dt) {
    constexpr int NFFT = NMAX;
    Wave& w = *wave_;
    if (w.cref.empty()) {
        w.cref.assign(NFRAME, cfloat(0.0f, 0.0f));
        w.cfilt.assign(NMAX, cfloat(0.0f, 0.0f));
        w.cw.assign(NMAX, cfloat(0.0f, 0.0f));
        w.xjunk.assign(NFRAME, 0.0f);
    }
    cfloat* cref = w.cref.data();
    cfloat* cfilt = w.cfilt.data();
    cfloat* cw = w.cw.data();

    const int nstart = static_cast<int>((dt * 12000.0f + 1.0f) - 576.0f);
    if (trace_) {
        std::string t;
        for (int i = 0; i < NN; ++i) t += static_cast<char>('0' + itone[i]);
        trace_("U " + XF(f0) + " " + XF(dt) + " " + I(nstart) + " " + t + " " + X(h_r(dd, NMAX)));
    }
    gen_ft4wave(itone, 103, NSPS, 12000.0f, f0, cref, w.xjunk.data(), 1, NFRAME);
    if (trace_) trace_("V " + X(h_c(cref, NFRAME)));
    if (fast_subtract_) {
        subtractft4_fast(dd, nstart);
        return;
    }
    // camp = dd*conjg(cref), zero elsewhere, made where cfilt = camp puts
    // it (cfilt(1:NFRAME) = camp(1:NFRAME), the rest zero, as
    // camp(NFRAME+1:NMAX) is; camp itself is read nowhere else).
    std::fill(cfilt, cfilt + NMAX, cfloat(0.0f, 0.0f));
    for (int i = 1; i <= NFRAME; ++i) {
        const int id = nstart - 1 + i;
        if (id >= 1 && id <= NMAX) {
            const cfloat cr = cref[i - 1];
            cfilt[i - 1] = cmul(cfloat(dd[id - 1], 0.0f), cfloat(cr.real(), -cr.imag()));
        }
    }
    if (trace_) trace_("v " + X(h_c(cfilt, NMAX)));

    if (w.st_first) {
        // Create and normalize the filter
        const float fac = bitsf(0x37672AE4u);  // 1.0/float(nfft), folded
        float window[NFILT + 1];                // window(-NFILT/2:NFILT/2) at [j+700]
        // window(-700) folded by gcc; -699..700 by libmvec, four at a time.
        window[0] = bitsf(0x2709ADEDu);
        float sum = 0.0f;
        sum = sum + window[0];
        for (int j = -NFILT / 2 + 1; j <= NFILT / 2; j += 4) {
            float a[4], c[4];
            for (int q = 0; q < 4; ++q) a[q] = (PI * static_cast<float>(j + q)) / static_cast<float>(NFILT);
            cosf4_glibc(a, c);
            for (int q = 0; q < 4; ++q) {
                window[j + q + NFILT / 2] = c[q] * c[q];
                sum = sum + window[j + q + NFILT / 2];
            }
        }
        std::vector<cfloat> tmp(NMAX, cfloat(0.0f, 0.0f));
        for (int i = 1; i <= NFILT + 1; ++i) tmp[i - 1] = cfloat(window[i - 1] / sum, 0.0f);
        for (int m = 0; m < NMAX; ++m) cw[m] = tmp[(m + NFILT / 2 + 1) % NMAX];  // cshift(cw, NFILT/2+1)
        if (trace_) {
            const int a = static_cast<int>(reinterpret_cast<uintptr_t>(cw) % 16);
            if (a != w.al_cw) {
                w.al_cw = a;
                trace_("@ sub-cw " + I(a));
            }
        }
        fst4w::fft_c2c_inplace(cw, NFFT, -1);
        for (int m = 0; m < NMAX; ++m) cw[m] = cmul(cw[m], cfloat(fac, 0.0f));
        if (trace_) trace_("F " + X(h_r(window, NFILT + 1)) + " " + XF(sum) + " " + X(h_c(cw, NMAX)));
        w.st_first = false;
    }

    if (trace_) {
        const int a = static_cast<int>(reinterpret_cast<uintptr_t>(cfilt) % 16);
        if (a != w.al_cfilt) {
            w.al_cfilt = a;
            trace_("@ sub-cfilt " + I(a));
        }
    }
    fst4w::fft_c2c_inplace(cfilt, NFFT, -1);
    for (int m = 0; m < NFFT; ++m) cfilt[m] = cmul(cfilt[m], cw[m]);
    fst4w::fft_c2c_inplace(cfilt, NFFT, 1);
    if (trace_) trace_("f " + X(h_c(cfilt, NMAX)));

    // Subtract the reconstructed signal
    for (int i = 1; i <= NFRAME; ++i) {
        const int j = nstart + i - 1;
        if (j >= 1 && j <= NMAX) {
            const float z = cfilt[i - 1].real() * cref[i - 1].real() - cfilt[i - 1].imag() * cref[i - 1].imag();
            dd[j - 1] = dd[j - 1] - (z + z);
        }
    }
    if (trace_) trace_("u " + X(h_r(dd, NMAX)));
}

// Fast subtraction (fast_subtract_): subtractft4's smoothing (the 1401-tap window over the 60,480-sample frame, no
// end correction) by ft_fast_subtract.h's, at 857 S/s (14-sample blocks: the window is 100 blocks a side, as FT8's
// 4000 taps are at 40), and subtracted in dd.
void Frontend::subtractft4_fast(float* dd, int nstart) {
    Wave& w = *wave_;
    if (!w.fs) w.fs = std::make_unique<ft8::FastSubtract>(NFRAME, NFILT, 14);
    // i in [ilo, ihi] (1-based) has j = nstart - 1 + i in [1, NMAX]
    const int ilo = std::max(1, 2 - nstart), ihi = std::min(NFRAME, NMAX - nstart + 1);
    const int n0 = std::min(std::max(ilo - 1, 0), NFRAME), n1 = std::max(ihi, n0);
    float* d = dd + (nstart - 1);
    w.fs->subtract(d, w.cref.data(), n0, n1, nullptr, d);
    if (trace_) trace_("u " + X(h_r(dd, NMAX)));
}

}  // namespace ft4
