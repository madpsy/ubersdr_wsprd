// MSK40, the short-message form of MSK144 (msk40.h): WSJT-X's msk40spd.f90,
// msk40sync.f90, msk40_freq_search.f90, msk40decodeframe.f90, bpdecode40.f90,
// encode_msk40.f90, update_msk40_hasharray.f90 and hash.f90 at 567ad29ce, line
// for line. Arrays keep the Fortran's order; the 1-based loops of the Fortran
// are shifted to 0-based here with the Fortran's index in a comment where it
// helps (bpdecode40 keeps the Fortran's 1-based indices, element 0 unused).
//
// Joe Taylor, K1JT, Steven Franke, K9AN, and the WSJT-X team (GNU GPL v3),
// and so is this. lookup3.c: Bob Jenkins, 2006, public domain.
#include "msk40.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <vector>

#include "ft8_maths.h"
#include "msk144_decode.h"  // fmtmsg
#include "msk144_ldpc.h"    // platanh
#include "msk144_sync.h"    // indexx, maxloc_f, msk_clip_spectrum_window, msk_peak_offset
#include "ft8_ldpc.h"       // ft8::tanhf_glibc

namespace msk144 {

Msk40Trace msk40_trace;

namespace {

constexpr int NSPM = 240;  // samples per MSK40 frame

// complex products as gfortran writes them (-fcx-fortran-rules): the plain
// formula, no NaN recovery.
inline cfloat cmul(cfloat a, cfloat b) {
    const float ar = a.real(), ai = a.imag(), br = b.real(), bi = b.imag();
    return cfloat(ar * br - ai * bi, ar * bi + ai * br);
}
inline cfloat conjg(cfloat a) { return cfloat(a.real(), -a.imag()); }

// nint of a float (the reference calls lroundf) and of a double (lround), as the default integer
inline int nint_f(float x) { return static_cast<int>(llroundf(x)); }
inline int nint_d(double x) { return static_cast<int>(llround(x)); }

// MAX(a, b) as gfortran -O3 writes it for the reference's msk40spd: b if it is greater than a, else a
inline float max_gf(float a, float b) { return b > a ? b : a; }

// ---- lib/lookup3.c's wsjt_hashlittle, lines 160-305, its byte-at-a-time path (the aligned word paths give the same value)
#define rot(x, k) (((x) << (k)) | ((x) >> (32 - (k))))
#define mix(a, b, c)                                                                                     \
    {                                                                                                    \
        a -= c; a ^= rot(c, 4);  c += b;                                                                 \
        b -= a; b ^= rot(a, 6);  a += c;                                                                 \
        c -= b; c ^= rot(b, 8);  b += a;                                                                 \
        a -= c; a ^= rot(c, 16); c += b;                                                                 \
        b -= a; b ^= rot(a, 19); a += c;                                                                 \
        c -= b; c ^= rot(b, 4);  b += a;                                                                 \
    }
#define final_(a, b, c)                                                                                  \
    {                                                                                                    \
        c ^= b; c -= rot(b, 14);                                                                         \
        a ^= c; a -= rot(c, 11);                                                                         \
        b ^= a; b -= rot(a, 25);                                                                         \
        c ^= b; c -= rot(b, 16);                                                                         \
        a ^= c; a -= rot(c, 4);                                                                          \
        b ^= a; b -= rot(a, 14);                                                                         \
        c ^= b; c -= rot(b, 24);                                                                         \
    }

uint32_t hashlittle(const void* key, size_t length, uint32_t initval) {
    uint32_t a, b, c;
    a = b = c = 0xdeadbeef + static_cast<uint32_t>(length) + initval;
    const uint8_t* k = static_cast<const uint8_t*>(key);
    while (length > 12) {  // all but the last block
        a += k[0];
        a += static_cast<uint32_t>(k[1]) << 8;
        a += static_cast<uint32_t>(k[2]) << 16;
        a += static_cast<uint32_t>(k[3]) << 24;
        b += k[4];
        b += static_cast<uint32_t>(k[5]) << 8;
        b += static_cast<uint32_t>(k[6]) << 16;
        b += static_cast<uint32_t>(k[7]) << 24;
        c += k[8];
        c += static_cast<uint32_t>(k[9]) << 8;
        c += static_cast<uint32_t>(k[10]) << 16;
        c += static_cast<uint32_t>(k[11]) << 24;
        mix(a, b, c);
        length -= 12;
        k += 12;
    }
    switch (length) {  // the last block: all the cases fall through
        case 12: c += static_cast<uint32_t>(k[11]) << 24; [[fallthrough]];
        case 11: c += static_cast<uint32_t>(k[10]) << 16; [[fallthrough]];
        case 10: c += static_cast<uint32_t>(k[9]) << 8; [[fallthrough]];
        case 9: c += k[8]; [[fallthrough]];
        case 8: b += static_cast<uint32_t>(k[7]) << 24; [[fallthrough]];
        case 7: b += static_cast<uint32_t>(k[6]) << 16; [[fallthrough]];
        case 6: b += static_cast<uint32_t>(k[5]) << 8; [[fallthrough]];
        case 5: b += k[4]; [[fallthrough]];
        case 4: a += static_cast<uint32_t>(k[3]) << 24; [[fallthrough]];
        case 3: a += static_cast<uint32_t>(k[2]) << 16; [[fallthrough]];
        case 2: a += static_cast<uint32_t>(k[1]) << 8; [[fallthrough]];
        case 1: a += k[0]; break;
        case 0: return c;
    }
    final_(a, b, c);
    return c;
}
#undef rot
#undef mix
#undef final_

#include "msk40_tables.inc"

// ---- msk40sync's and msk40decodeframe's constants, made at the first call: pp, the half-sine pulse, and
// cb, the sync word waveforms (s8r = 1,0,1,1,0,0,0,1, then 2*s8r-1: the reverse of MSK144's).
// The reference's compiler vectorised the pp loop (libmvec's _ZGVbN4v_sinf for i = 2..9, i = 1 is sin 0,
// i = 10..12 folded to constants), and folded different constants in the two routines, because
// msk40sync's pi is a default real (4.0*atan(1.0)) and msk40decodeframe's a real*8 (4d0*datan(1d0)).
// The constants are the correctly rounded sinf of each routine's float angles of i = 10, 11, 12, as in
// src/msk144_sync.cpp and src/msk144_decode.cpp (their comments say how they were read); the oracle's
// pp and cb lines check them.
struct Consts {
    float pp[12];
    cfloat cb[42];
    Consts(bool single_pi) {
        float ang[12];
        if (single_pi) {
            const float pi = 4.0f * 0.785398185253143310546875f;  // 4.0*atan(1.0): atan(1.0) folded to the float nearest pi/4
            for (int i = 1; i <= 12; ++i) ang[i - 1] = static_cast<float>(i - 1) * pi / 12.0f;
        } else {
            const double pi = 4.0 * 0.78539816339744830961566084581988;
            for (int i = 1; i <= 12; ++i) ang[i - 1] = static_cast<float>((static_cast<double>(i - 1) * pi) / 12.0);
        }
        pp[0] = 0.0f;
        for (int g = 0; g < 2; ++g) ft8::sinf4_glibc(ang + 1 + 4 * g, pp + 1 + 4 * g);
        const uint32_t folded_single[3] = {0x3F3504F3u, 0x3EFFFFFAu, 0x3E8483EAu};
        const uint32_t folded_double[3] = {0x3F3504F3u, 0x3F000001u, 0x3E8483F1u};
        const uint32_t* folded = single_pi ? folded_single : folded_double;
        for (int i = 0; i < 3; ++i) std::memcpy(&pp[9 + i], &folded[i], 4);

        int s8r[8] = {1, 0, 1, 1, 0, 0, 0, 1};
        for (int& v : s8r) v = 2 * v - 1;
        float cbi[42], cbq[42];
        for (int i = 0; i < 6; ++i) cbq[i] = pp[6 + i] * static_cast<float>(s8r[0]);       // cbq(1:6)=pp(7:12)*s8r(1)
        for (int i = 0; i < 12; ++i) cbq[6 + i] = pp[i] * static_cast<float>(s8r[2]);      // cbq(7:18)=pp*s8r(3)
        for (int i = 0; i < 12; ++i) cbq[18 + i] = pp[i] * static_cast<float>(s8r[4]);     // cbq(19:30)=pp*s8r(5)
        for (int i = 0; i < 12; ++i) cbq[30 + i] = pp[i] * static_cast<float>(s8r[6]);     // cbq(31:42)=pp*s8r(7)
        for (int i = 0; i < 12; ++i) cbi[i] = pp[i] * static_cast<float>(s8r[1]);          // cbi(1:12)=pp*s8r(2)
        for (int i = 0; i < 12; ++i) cbi[12 + i] = pp[i] * static_cast<float>(s8r[3]);     // cbi(13:24)=pp*s8r(4)
        for (int i = 0; i < 12; ++i) cbi[24 + i] = pp[i] * static_cast<float>(s8r[5]);     // cbi(25:36)=pp*s8r(6)
        for (int i = 0; i < 6; ++i) cbi[36 + i] = pp[i] * static_cast<float>(s8r[7]);      // cbi(37:42)=pp(1:6)*s8r(8)
        for (int i = 0; i < 42; ++i) cb[i] = cfloat(cbi[i], cbq[i]);
    }
};

// ---- msk40spd's constants: rcw, the raised-cosine edge window (msk40spd.f90:62-80, as msk144spd's)
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
        if (msk40_trace.spd_consts) msk40_trace.spd_consts(rcw);
    }
};

std::string trim_right(const std::string& s) {
    size_t n = s.size();
    while (n > 0 && s[n - 1] == ' ') --n;
    return s.substr(0, n);
}

// write(msgreceived, ...) into a character*37: the text, blank padded
void put37(char* msgreceived, const std::string& text) {
    const size_t n = std::min<size_t>(text.size(), 37);
    std::memcpy(msgreceived, text.data(), n);
    if (n < 37) std::memset(msgreceived + n, ' ', 37 - n);
}

}  // namespace

// ---- hash.f90:1-10: ihash=nhash(c_loc(string),len,146), ihash=iand(i,MASK15)
int hash_f90(const char* string, int len) {
    // nhash.c's own mask (0x7fff, for a key that is not empty) is the same as MASK15
    uint32_t h = hashlittle(string, static_cast<size_t>(len), 146u);
    if (len != 0) h &= 0x7fffu;
    return static_cast<int>(h & 32767u);
}

// ---- encode_msk40.f90:1-46
void encode_msk40(const signed char message[16], signed char codeword[32]) {
    // data g/Z'4428',Z'5a6b',Z'1b04',Z'2c12',Z'60c4',Z'1071',Z'be6a',Z'36dd', Z'c580',Z'ad9a',Z'eca2',Z'7843',
    //        Z'332e',Z'a685',Z'5906',Z'1efe'/
    static const int g[16] = {0x4428, 0x5a6b, 0x1b04, 0x2c12, 0x60c4, 0x1071, 0xbe6a, 0x36dd,
                              0xc580, 0xad9a, 0xeca2, 0x7843, 0x332e, 0xa685, 0x5906, 0x1efe};
    // fill the generator matrix: gen40(i,j)=1 if btest(g(i),16-j)
    struct Gen {
        signed char gen40[16][16];  // gen40[i-1][j-1]
        Gen() {
            std::memset(gen40, 0, sizeof gen40);
            for (int i = 1; i <= 16; ++i)
                for (int j = 1; j <= 16; ++j)
                    if ((g[i - 1] >> (16 - j)) & 1) gen40[i - 1][j - 1] = 1;
        }
    };
    static const Gen G;
    signed char pchecks[16];
    signed char itmp[32];
    for (int i = 1; i <= 16; ++i) {
        int nsum = 0;
        for (int j = 1; j <= 16; ++j) nsum = nsum + message[j - 1] * G.gen40[i - 1][j - 1];
        pchecks[i - 1] = static_cast<signed char>(nsum % 2);
    }
    for (int i = 0; i < 16; ++i) itmp[i] = pchecks[i];       // itmp(1:16)=pchecks
    for (int i = 0; i < 16; ++i) itmp[16 + i] = message[i];  // itmp(17:32)=message(1:16)
    for (int i = 0; i < 32; ++i) codeword[kColorder[i]] = itmp[i];  // codeword(colorder+1)=itmp(1:32)
}

// ---- bpdecode40.f90:1-148. The Fortran's Nm(i,j), Mn(kk,j), nrw(j) and the arrays tov(3,N), toc(7,M),
// tanhtoc(7,M), zn(N), llr(N) keep their 1-based indices (element 0 unused).
void bpdecode40(const float llr_in[32], int maxiterations, signed char decoded[16], int& niterations) {
    constexpr int N = 32, K = 16, M = N - K;
    float tov[4][N + 1];
    float toc[8][M + 1];
    float tanhtoc[8][M + 1];
    float zn[N + 1];
    float llr[N + 1];
    int cw[N + 1];
    float Tmn;
    const int ncw = 3;
    auto Nm = [](int i, int j) { return static_cast<int>(kNm[j - 1][i - 1]); };
    auto Mn = [](int kk, int j) { return static_cast<int>(kMn[j - 1][kk - 1]); };
    auto nrw = [](int j) { return static_cast<int>(kNrw[j - 1]); };
    for (int i = 1; i <= N; ++i) llr[i] = llr_in[i - 1];

    std::memset(toc, 0, sizeof toc);
    std::memset(tov, 0, sizeof tov);
    std::memset(tanhtoc, 0, sizeof tanhtoc);

    // initialize messages to checks
    for (int j = 1; j <= M; ++j)
        for (int i = 1; i <= nrw(j); ++i) toc[i][j] = llr[Nm(i, j)];

    for (int iter = 0; iter <= maxiterations; ++iter) {
        // Update bit log likelihood ratios (tov=0 in iteration 0).
        for (int i = 1; i <= N; ++i) {
            float s = 0.0f;  // sum(tov(1:ncw,i)), in order
            for (int kk = 1; kk <= ncw; ++kk) s += tov[kk][i];
            zn[i] = llr[i] + s;
        }

        if (msk40_trace.bp) msk40_trace.bp(iter, zn, tov);  // the piece check's hook, not in the Fortran

        // Check to see if we have a codeword (check before we do any iteration).
        for (int i = 1; i <= N; ++i) cw[i] = zn[i] > 0.f ? 1 : 0;
        int ncheck = 0;
        for (int i = 1; i <= M; ++i) {
            int synd = 0;
            for (int k = 1; k <= nrw(i); ++k) synd += cw[Nm(k, i)];
            if (synd % 2 != 0) ++ncheck;
        }

        if (ncheck == 0) {  // we have a codeword - reorder the columns and return it
            niterations = iter;
            // codeword=cw(colorder+1); decoded=codeword(M+1:N)
            for (int i = 1; i <= K; ++i) decoded[i - 1] = static_cast<signed char>(cw[kColorder[M + i - 1] + 1]);
            return;
        }

        // Send messages from bits to check nodes
        for (int j = 1; j <= M; ++j)
            for (int i = 1; i <= nrw(j); ++i) {
                const int ibj = Nm(i, j);
                toc[i][j] = zn[ibj];
                for (int kk = 1; kk <= ncw; ++kk)  // subtract off what the bit had received from the check
                    if (Mn(kk, ibj) == j) toc[i][j] = toc[i][j] - tov[kk][ibj];
            }

        // send messages from check nodes to variable nodes
        for (int i = 1; i <= M; ++i)
            for (int r = 1; r <= 7; ++r) tanhtoc[r][i] = ft8::tanhf_glibc(-toc[r][i] / 2);

        for (int j = 1; j <= N; ++j)
            for (int i = 1; i <= ncw; ++i) {
                const int ichk = Mn(i, j);  // Mn(:,j) are the checks that include bit j
                // product(tanhtoc(1:nrw(ichk),ichk), mask=Nm(1:nrw(ichk),ichk).ne.j)
                Tmn = 1.0f;
                for (int r = 1; r <= nrw(ichk); ++r)
                    if (Nm(r, ichk) != j) Tmn = Tmn * tanhtoc[r][ichk];
                const float y = platanh(-Tmn);
                tov[i][j] = 2 * y;
            }
    }
    niterations = -1;
}

// ---- update_msk40_hasharray.f90:1-25
void update_msk40_hasharray(Msk40State& ms, const Packjt77& pk) {
    constexpr int MAXRECENT = Packjt77::MAXRECENT;
    for (int i = 0; i < MAXRECENT; ++i)
        for (int j = 0; j < MAXRECENT; ++j) ms.nhasharray[i][j] = -1;  // nhasharray=-1
    for (int i = 1; i <= MAXRECENT; ++i)
        for (int j = i + 1; j <= MAXRECENT; ++j) {
            const std::string ci = pk.recent_call(i), cj = pk.recent_call(j);
            if (ci[0] != ' ' && cj[0] != ' ') {
                // hashmsg=trim(recent_calls(i))//' '//trim(recent_calls(j)), a character*37
                std::string hashmsg = trim_right(ci) + ' ' + trim_right(cj);
                hashmsg.resize(37, ' ');
                int iz;
                fmtmsg(hashmsg, iz);
                int ihash = hash_f90(hashmsg.data(), 37);
                ihash = ihash & 4095;
                ms.nhasharray[i - 1][j - 1] = ihash;
                hashmsg = trim_right(cj) + ' ' + trim_right(ci);
                hashmsg.resize(37, ' ');
                fmtmsg(hashmsg, iz);
                ihash = hash_f90(hashmsg.data(), 37);
                ihash = ihash & 4095;
                ms.nhasharray[j - 1][i - 1] = ihash;
            }
        }
}

// ---- msk40decodeframe.f90:1-163
void msk40decodeframe(Msk40State& ms, const Packjt77& pk, cfloat* c, const std::string& mycall_in,
                      const std::string& hiscall_in, float xsnr, bool bswl, char* msgreceived, int& nsuccess) {
    static const char* const rpt[16] = {"-03 ", "+00 ", "+03 ", "+06 ", "+10 ", "+13 ", "+16 ", "R-03",
                                        "R+00", "R+03", "R+06", "R+10", "R+13", "R+16", "RRR ", "73  "};
    static const Consts K(false);  // first=.true. branch, once
    const float* pp = K.pp;
    const cfloat* cb = K.cb;
    static const bool traced_consts = (msk40_trace.df_consts ? (msk40_trace.df_consts(K.pp, K.cb), true) : true);
    (void)traced_consts;
    static const int s8r[8] = {1, -1, 1, 1, -1, -1, -1, 1};  // s8r=2*s8r-1

    std::string mycall = mycall_in, hiscall = hiscall_in;  // character*12
    mycall.resize(12, ' ');
    hiscall.resize(12, ' ');

    if (mycall != ms.mycall0 || hiscall != ms.hiscall0) {
        std::string hashmsg = trim_right(mycall) + ' ' + trim_right(hiscall);  // a character*37
        hashmsg.resize(37, ' ');
        if (hashmsg != std::string(37, ' ') && hiscall != std::string(12, ' ')) {  // protect against blank mycall/hiscall
            int iz;
            fmtmsg(hashmsg, iz);
            int ihash = hash_f90(hashmsg.data(), 37);
            ms.ihash = ihash & 4095;
        } else {
            ms.ihash = 9999;  // so that it can never match a received hash
        }
        ms.mycall0 = mycall;
        ms.hiscall0 = hiscall;
    }

    nsuccess = 0;
    std::memset(msgreceived, ' ', 37);

    // Estimate carrier phase.
    // cca=sum(c(1:1+41)*conjg(cb)): a running sum from zero, in order
    cfloat cca(0.0f, 0.0f);
    for (int i = 0; i < 42; ++i) {
        const cfloat p = cmul(c[i], conjg(cb[i]));
        cca = cfloat(cca.real() + p.real(), cca.imag() + p.imag());
    }
    const float phase0 = ft8::atan2f_glibc(cca.imag(), cca.real());
    if (msk40_trace.df_phase) msk40_trace.df_phase(cca, phase0);

    // Remove phase error - want constellation rotated so that sample points lie on I/Q axes
    float sphase, cphase;
    ft8::sincosf_glibc(phase0, &sphase, &cphase);  // cfac=cmplx(cos(phase0),sin(phase0))
    const cfloat cfac(cphase, sphase);
    const cfloat cfacc = conjg(cfac);
    for (int i = 0; i < NSPM; ++i) c[i] = cmul(c[i], cfacc);  // c=c*conjg(cfac)
    if (msk40_trace.df_corrected) msk40_trace.df_corrected(cfac, c);

    // Matched filter.
    float softbits[40];
    {
        float s1 = 0.0f, s2 = 0.0f;
        for (int j = 0; j < 6; ++j) s1 = s1 + c[j].imag() * pp[6 + j];                 // imag(c(1:6))*pp(7:12)
        for (int j = 0; j < 6; ++j) s2 = s2 + c[NSPM - 6 + j].imag() * pp[j];          // imag(c(NSPM-5:NSPM))*pp(1:6)
        softbits[0] = s1 + s2;
        float s3 = 0.0f;
        for (int j = 0; j < 12; ++j) s3 = s3 + c[j].real() * pp[j];                    // real(c(1:12))*pp
        softbits[1] = s3;
    }
    for (int i = 2; i <= 20; ++i) {
        float sa = 0.0f, sb = 0.0f;
        // imag(c(1+(i-1)*12-6:1+(i-1)*12+5))*pp: 0-based start (i-1)*12-6
        for (int j = 0; j < 12; ++j) sa = sa + c[(i - 1) * 12 - 6 + j].imag() * pp[j];
        // real(c(7+(i-1)*12-6:7+(i-1)*12+5))*pp: 0-based start (i-1)*12
        for (int j = 0; j < 12; ++j) sb = sb + c[(i - 1) * 12 + j].real() * pp[j];
        softbits[2 * i - 2] = sa;  // softbits(2*i-1)
        softbits[2 * i - 1] = sb;  // softbits(2*i)
    }
    if (msk40_trace.df_raw) msk40_trace.df_raw(softbits);

    // Sync word hard error weight is used to reject frames that are unlikely to decode.
    int hardbits[40];
    for (int i = 0; i < 40; ++i) hardbits[i] = softbits[i] >= 0.0f ? 1 : 0;
    int sum1 = 0;
    for (int i = 0; i < 8; ++i) sum1 += (2 * hardbits[i] - 1) * s8r[i];
    const int nbadsync1 = (8 - sum1) / 2;
    const int nbadsync = nbadsync1;
    if (msk40_trace.df_sync) msk40_trace.df_sync(nbadsync1);
    if (nbadsync > 3) return;

    // Normalize the softsymbols before submitting to decoder.
    float ssum = 0.0f, s2sum = 0.0f;
    for (int i = 0; i < 40; ++i) ssum = ssum + softbits[i];
    for (int i = 0; i < 40; ++i) s2sum = s2sum + softbits[i] * softbits[i];
    const float sav = ssum / 40.0f;
    const float s2av = s2sum / 40.0f;
    const float ssig = std::sqrt(s2av - sav * sav);
    if (msk40_trace.df_norm) msk40_trace.df_norm(sav, s2av, ssig);
    for (int i = 0; i < 40; ++i) softbits[i] = softbits[i] / ssig;

    float sigma = 0.75f;
    if (xsnr < 0.0f) sigma = 0.75f - 0.11f * xsnr;
    float llr[32];
    for (int i = 0; i < 32; ++i) llr[i] = softbits[8 + i];  // llr(1:32)=softbits(9:40)
    for (int i = 0; i < 32; ++i) llr[i] = 2.0f * llr[i] / (sigma * sigma);
    if (msk40_trace.df_llr) msk40_trace.df_llr(sigma, llr);

    const int max_iterations = 5;
    signed char decoded[16] = {};  // integer*1 decoded(16): not set unless a codeword is found
    int niterations = 0;
    bpdecode40(llr, max_iterations, decoded, niterations);
    if (msk40_trace.df_bp) msk40_trace.df_bp(decoded, niterations);
    if (niterations >= 0) {
        signed char cw[32];
        encode_msk40(decoded, cw);
        int nhammd = 0;
        float cord = 0.0f;
        for (int i = 1; i <= 32; ++i) {
            if (cw[i - 1] != hardbits[i + 8 - 1]) {
                nhammd = nhammd + 1;
                cord = cord + std::fabs(softbits[i + 8 - 1]);
            }
        }

        int imsg = 0;
        for (int i = 1; i <= 16; ++i) imsg = imsg * 2 + (decoded[17 - i - 1] & 1);  // imsg=ishft(imsg,1)+iand(1_1,decoded(17-i))
        const int nrxrpt = imsg & 15;
        const int nrxhash = (imsg - nrxrpt) / 16;
        if (msk40_trace.df_msg) msk40_trace.df_msg(nhammd, cord, imsg, nrxrpt, nrxhash, ms.ihash);
        if (nhammd <= 4 && cord < 0.65f && nrxhash == ms.ihash && nrxrpt >= 7) {
            nsuccess = 1;
            // write(msgreceived,'(a1,a,1x,a,a1,1x,a4)') "<",trim(mycall),hiscall,">",rpt(nrxrpt)
            put37(msgreceived, "<" + trim_right(mycall) + " " + trim_right(hiscall) + "> " + rpt[nrxrpt]);
            return;
        } else if (bswl && nhammd <= 4 && cord < 0.65f && nrxrpt >= 7) {
            constexpr int MAXRECENT = Packjt77::MAXRECENT;
            for (int i = 1; i <= MAXRECENT; ++i) {
                for (int j = i + 1; j <= MAXRECENT; ++j) {
                    if (nrxhash == ms.nhasharray[i - 1][j - 1]) {
                        nsuccess = 2;
                        put37(msgreceived, "<" + trim_right(pk.recent_call(i)) + " " + trim_right(pk.recent_call(j)) +
                                               "> " + rpt[nrxrpt]);
                    } else if (nrxhash == ms.nhasharray[j - 1][i - 1]) {
                        nsuccess = 2;
                        put37(msgreceived, "<" + trim_right(pk.recent_call(j)) + " " + trim_right(pk.recent_call(i)) +
                                               "> " + rpt[nrxrpt]);
                    }
                }
            }
            if (nsuccess == 0) {
                nsuccess = 3;
                // write(msgreceived,'(a1,i4.4,a1,1x,a4)') "<",nrxhash,">",rpt(nrxrpt)
                char h[8];
                std::snprintf(h, sizeof h, "%04d", nrxhash);
                put37(msgreceived, std::string("<") + h + "> " + rpt[nrxrpt]);
            }
        }
    }
}

// ---- msk40_freq_search.f90:1-63
void msk40_freq_search(const cfloat* cdat, float fc, int if1, int if2, float delf, int nframes, const int* navmask,
                       const cfloat* cb, cfloat* cdat2, float& xmax, float& bestf, cfloat* cs, float* xccs) {
    cfloat c[NSPM];
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
    const float fac = 1.0f / (24.0f * std::sqrt(static_cast<float>(navg)));

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
            // dot_product(ct2(1+ish:42+ish), cb(1:42)): the sum of the conjugate of the first times the second
            float sr = 0.0f, si = 0.0f;
            for (int kk = 0; kk < 42; ++kk) {
                const float tr = ct2[ish + kk].real(), ti = ct2[ish + kk].imag();
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

// ---- msk40sync.f90:1-113
void msk40sync(const cfloat* cdat, int nframes, int ntol, float delf, const int* navmask, int npeaks, float fc,
               float& fest, int* npklocs, int& nsuccess, cfloat* c) {
    static const Consts K(true);
    static const bool traced_consts = (msk40_trace.sync_consts ? (msk40_trace.sync_consts(K.pp, K.cb), true) : true);
    (void)traced_consts;
    std::vector<cfloat> cdat2(static_cast<size_t>(NSPM) * static_cast<size_t>(std::max(nframes, 0)) + 1);
    cfloat cs[NSPM];
    float xcc[NSPM];  // xcc(0:NSPM-1)
    float xccs[NSPM];
    float xm = 0.0f, bf = 0.0f;

    const int nfreqs = 2 * nint_f(static_cast<float>(ntol) / delf) + 1;
    // nthreads = 1: id = 1, nstep = nfreqs, if1 = -nint(ntol/delf), and as id == nthreads, if2 = nint(ntol/delf)
    const int if1 = -nint_f(static_cast<float>(ntol) / delf);
    const int if2 = nint_f(static_cast<float>(ntol) / delf);
    msk40_freq_search(cdat, fc, if1, if2, delf, nframes, navmask, K.cb, cdat2.data(), xm, bf, cs, xccs);

    const float xmax = xm;
    fest = fc + bf;
    for (int i = 0; i < NSPM; ++i) c[i] = cs[i];
    for (int i = 0; i < NSPM; ++i) xcc[i] = xccs[i];

    if (msk40_trace.sync_peak) msk40_trace.sync_peak(xmax, fest, c, xcc, nfreqs);

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

// ---- msk40spd.f90:1-210 (msk40decodeframe is called at lines 194-195)
void msk40spd(Msk40State& ms, const Packjt77& pk, const cfloat* cbig, int n, int ntol, const std::string& mycall12,
              const std::string& hiscall12, bool bswl, int& nsuccess, char* msgreceived, float fc, float& fret,
              float& tret, int& navg) {
    constexpr int MAXSTEPS = 150, NFFT = NSPM, MAXCAND = 5, NPATTERNS = 6;
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
    cfloat ct[NSPM];
    alignas(32) cfloat ctmp[NFFT];
    int indices[MAXSTEPS];
    int npkloc[10] = {0};
    int navmask[3];
    int nstart[MAXCAND];
    float detmet[MAXSTEPS + 6];  // detmet(-2:MAXSTEPS+3): Fortran index k is [k+2]
    float detmet2[MAXSTEPS + 6];
    float detfer[MAXSTEPS];
    float ferrs[MAXCAND];
    float snrs[MAXCAND];
    float tonespec[NFFT];
    auto DM = [&](int k) -> float& { return detmet[k + 2]; };
    auto DM2 = [&](int k) -> float& { return detmet2[k + 2]; };

    // fill the detmet, detferr arrays
    const int nstep = std::min(MAXSTEPS, (n - NSPM) / 60);  // 20ms/4=5ms steps
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
    if (msk40_trace.spd_windows) msk40_trace.spd_windows(ihlo, ihhi, illo, ilhi, i2000, i4000);
    for (int istp = 1; istp <= nstep; ++istp) {
        const int ns = 1 + 60 * (istp - 1);
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
        if (msk40_trace.spd_step) msk40_trace.spd_step(istp, ihpk, ilpk, deltah, deltal, ah, al, trath, tratl, ferr);
    }  // end of detection-metric and frequency error estimation loop

    indexx(&DM(1), nstep, indices);  // find median of detection metric vector
    const float xmed = DM(indices[nstep / 4 - 1]);
    if (msk40_trace.spd_median) msk40_trace.spd_median(nstep, indices, xmed, &DM(1), &DM2(1), detfer);
    if (xmed <= 0.0f || !std::isfinite(xmed)) return;
    for (float& v : detmet) v = v / xmed;  // noise floor of detection metric is 1.0
    int ndet = 0;

    for (int ip = 1; ip <= MAXCAND; ++ip) {  // Find candidates
        const int il = maxloc_f(&DM(1), nstep);
        if (DM(il) < 3.5f) break;
        if (std::fabs(detfer[il - 1]) <= static_cast<float>(ntol)) {
            ndet = ndet + 1;
            nstart[ndet - 1] = 1 + (il - 1) * 60 + 1;
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
                nstart[ndet - 1] = 1 + (il - 1) * 60 + 1;
                ferrs[ndet - 1] = detfer[il - 1];
                snrs[ndet - 1] = 12.0f * ft8::log10f_glibc(DM2(il)) / 2.0f - 9.0f;
            }
            DM2(il) = 0.0f;
        }
    }

    if (msk40_trace.spd_candidates) msk40_trace.spd_candidates(ndet, nstart, ferrs, snrs);
    constexpr int npeaks = 2;
    constexpr int ntol0 = 29;
    constexpr float deltaf = 7.2f;
    for (int icand = 1; icand <= ndet; ++icand) {  // Try to sync/demod/decode each candidate.
        int ib = std::max(1, nstart[icand - 1] - NSPM);
        int ie = ib - 1 + 3 * NSPM;
        if (ie > n) {
            ie = n;
            ib = ie - 3 * NSPM + 1;
        }
        for (int i = 0; i < 3 * NSPM; ++i) cdat[i] = cbig[ib - 1 + i];
        const float fo = fc + ferrs[icand - 1];
        const float xsnr = snrs[icand - 1];
        for (int iav = 1; iav <= NPATTERNS; ++iav) {
            for (int i = 0; i < 3; ++i) navmask[i] = navpatterns[iav - 1][i];
            float fest;
            int nsyncsuccess;
            msk40sync(cdat, 3, ntol0, deltaf, navmask, npeaks, fo, fest, npkloc, nsyncsuccess, c);
            if (msk40_trace.spd_sync) msk40_trace.spd_sync(icand, iav, nsyncsuccess, fest, npkloc);
            if (nsyncsuccess == 0) continue;

            for (int ipk = 1; ipk <= npeaks; ++ipk) {
                for (int is = 1; is <= 3; ++is) {
                    int ic0 = npkloc[ipk - 1];
                    if (is == 2) ic0 = std::max(1, ic0 - 1);
                    if (is == 3) ic0 = std::min(NSPM, ic0 + 1);
                    // ct=cshift(c,ic0-1)
                    for (int i = 0; i < NSPM; ++i) ct[i] = c[(i + ic0 - 1) % NSPM];
                    int ndecodesuccess = 0;
                    msk40decodeframe(ms, pk, ct, mycall12, hiscall12, xsnr, bswl, msgreceived, ndecodesuccess);
                    if (ndecodesuccess > 0) {
                        tret = static_cast<float>(static_cast<double>(nstart[icand - 1] + NSPM / 2) / fs);
                        fret = fest;
                        navg = navmask[0] + navmask[1] + navmask[2];
                        nsuccess = ndecodesuccess;
                        return;
                    }
                }
            }
        }
    }  // candidate loop
}

}  // namespace msk144
