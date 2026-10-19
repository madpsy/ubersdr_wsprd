// The built-in FT2 decoder (ft2.h): improved 3.2.0's FT2-own code
// (decoder.f90's FT2 branch, ft2_decoded, jt9.f90's --ft2 file mode) on
// 567ad29ce's FT4 (decoder.f90's RMS gate and decoded.txt, ft4_decode.f90
// whole, decode_completion.f90), line for line, as the FT2 hybrid of
// tests/data/ft2/README.md has them; and the oracle harness's FT2 modes
// (scripts/ft8-oracle.f90 -2, --file -2) as run_periods.
//
// ft4_decode here is src/ft4.cpp's, copied (FT4's validated code is not
// touched) with FT2's differences: NMAX 79,488 (dd, iwave, the trace's
// hashes), NDMAX 4,416 (cd2, cb, sum2's divisor 4416.0, cd's cut at
// min(NDMAX-1, ...)), and syncmin 1.0 (improved's ft4_decode.f90:246,
// "if(tperiod.eq.3.75) syncmin=1").
//
// Joe Taylor, K1JT, Steven Franke, K9AN, and the WSJT-X developers (GNU GPL
// v3); FT2's code from WSJT-X improved 3.2.0 (Uwe Risse, DG2YCB, and the
// same); and so is this.
#include "ft2.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <type_traits>

#include "fst4w_jt9.h"    // fortran_f, fortran_i
#include "ft8_fortran.h"  // fix, trim, findex, nint, maxval, twkfreq1, fbits, the trace hashes
#include "ft8_ldpc.h"     // decode174_91, encode174_91
#include "ft8_maths.h"    // log10f_glibc

namespace ft2 {

using ft8::fortran::fbits;
using ft8::fortran::findex;
using ft8::fortran::fix;
using ft8::fortran::h_c;
using ft8::fortran::h_i;
using ft8::fortran::h_r;
using ft8::fortran::hw;
using ft8::fortran::L;
using ft8::fortran::maxval;
using ft8::fortran::nint;
using ft8::fortran::trim;
using ft8::fortran::twkfreq1;
using ft8::fortran::X;
using ft8::fortran::XF;

namespace {

constexpr int NPTS = 180000;  // id2's length (the RMS gate's)

// ft4_decode's frame from decodes(100) to the end of cd: the hybrid's
// Release object has 567ad29ce's offsets (decodes 0x2980(%rsp), savg
// 0x3800, sbase 0x4a00, cd 0x5c00), ft4.h's decodes(100) note.
constexpr size_t DLEN = 37;
constexpr size_t SAVG_AT = 0x3800 - 0x2980;
constexpr size_t SBASE_AT = 0x4a00 - 0x2980;
constexpr size_t CD_AT = 0x5c00 - 0x2980;
constexpr size_t FRAME_BYTES = CD_AT + NN * NSS * sizeof(cfloat);
static_assert(SAVG_AT == 100 * DLEN + 12 && SBASE_AT == SAVG_AT + NH1 * sizeof(float) &&
              CD_AT == SBASE_AT + NH1 * sizeof(float) && CD_AT + NN * NSS * sizeof(cfloat) == 0xc300 - 0x2980);
static_assert(3 * 200 * DLEN <= FRAME_BYTES, "the most decodes one call can make stay inside cd");

// The front-end trace's raw dumps: every 499th get_ft4_bitmetrics call
// from the first (scripts/ft8-frontend-patch.py's FT4_NSAMP).
constexpr long FE_NSAMP = 499;

// ft4_decode's DATA tables (ft4.cpp's).
const int8_t RVEC[77] = {0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0,
                         1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1,
                         1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1};
const int MCQ[29] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0};
const int MCQRU[29] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0};
const int MCQFD[29] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0};
const int MCQTEST[29] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0};
const int MCQWW[29] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0};

struct Ap29 {
    int v[29];
    explicit Ap29(const int* m) {
        for (int i = 0; i < 29; ++i) v[i] = 2 * ((m[i] + RVEC[i]) % 2) - 1;
    }
};
const Ap29& mcq() {
    static const Ap29 a(MCQ);
    return a;
}
const Ap29& mcqru() {
    static const Ap29 a(MCQRU);
    return a;
}
const Ap29& mcqfd() {
    static const Ap29 a(MCQFD);
    return a;
}
const Ap29& mcqtest() {
    static const Ap29 a(MCQTEST);
    return a;
}
const Ap29& mcqww() {
    static const Ap29 a(MCQWW);
    return a;
}

const int NAPPASSES[6] = {2, 2, 2, 2, 2, 3};
const int NAPTYPES[6][4] = {{1, 2, 0, 0}, {2, 3, 0, 0}, {2, 3, 0, 0}, {3, 6, 0, 0}, {3, 6, 0, 0}, {3, 1, 2, 0}};

std::string hex8(float v) {
    char b[12];
    std::snprintf(b, sizeof b, "%08X", fbits(v));
    return b;
}

std::string I(long v) { return std::to_string(v); }

// The hash of integer*2 samples, sign-extended (the trace's feh_i2).
uint32_t h_i2(const int16_t* a, size_t n) {
    uint32_t h = 2166136261u;
    for (size_t i = 0; i < n; ++i) h = hw(h, static_cast<uint32_t>(static_cast<int>(a[i])));
    return h;
}
// Float samples' (ours, no reference trace): their bits.
uint32_t h_i2(const float* a, size_t n) {
    uint32_t h = 2166136261u;
    for (size_t i = 0; i < n; ++i) h = hw(h, fbits(a[i]));
    return h;
}

// The float samples' RMS gate (call()): only all zeros is silence.
bool all_zero(const float* x, int n) {
    for (int i = 0; i < n; ++i)
        if (x[i] != 0.0f) return false;
    return true;
}
bool all_zero(const int16_t*, int) { return false; }  // not asked: the 16-bit gate is jt9's

// sum2 = sum(cd2*conjg(cd2))/(real(NMAX)/real(NDOWN)) as ft4.cpp's
// sum2_cd2, over NDMAX 4,416, divided by 4416.0.
float sum2_cd2(const cfloat* c) {
    float sr = 0.0f, si = 0.0f;
    for (int k = 0; k < NDMAX; ++k) {
        const float re = c[k].real(), im = c[k].imag(), nim = -im;
        const float p = re * nim;
        sr = sr + (re * re - im * nim);
        si = si + (p - p);
    }
    return (si * 0.0f + sr) / 4416.0f;
}
static_assert(static_cast<float>(NMAX) / static_cast<float>(NDOWN) == 4416.0f, "sum2's divisor");

// sum(abs(cb)**2)/(real(NSS)*NN) over NDMAX, as ft4.cpp's sum2_cb.
float sum2_cb(const cfloat* c) {
    float a[NDMAX];
    ft8::hypotf_c_glibc(reinterpret_cast<const float*>(c), a, NDMAX);
    float s = 0.0f;
    for (int k = 0; k < NDMAX; ++k) s = s + a[k] * a[k];
    return s / 3296.0f;
}

// c = c/sqrt(sum2) over NDMAX, as ft4.cpp's normalise.
void normalise(cfloat* c, float sum2) {
    const float s = std::sqrt(sum2);
    const float r = 0.0f / s;
    const float d = s + 0.0f * r;
    for (int k = 0; k < NDMAX; ++k) {
        const float a = c[k].real(), b = c[k].imag();
        c[k] = cfloat((a + b * r) / d, (b - a * r) / d);
    }
}

// jt9_input_validation.f90's is_real_token (as ft4.cpp's).
bool real_token(const std::string& t) {
    const size_t n = t.size();
    if (n == 0) return false;
    size_t i = 0;
    if (t[i] == '+' || t[i] == '-') ++i;
    if (i >= n) return false;
    int digits = 0;
    while (i < n && t[i] >= '0' && t[i] <= '9') ++digits, ++i;
    if (i < n && t[i] == '.') {
        ++i;
        while (i < n && t[i] >= '0' && t[i] <= '9') ++digits, ++i;
    }
    if (digits == 0) return false;
    if (i < n) {
        if (t[i] != 'e' && t[i] != 'E') return false;
        ++i;
        if (i < n && (t[i] == '+' || t[i] == '-')) ++i;
        if (i >= n) return false;
        digits = 0;
        while (i < n) {
            if (t[i] < '0' || t[i] > '9') return false;
            ++digits, ++i;
        }
        if (digits == 0) return false;
    }
    return i >= n;
}

std::string digits(const int* v, int n) {
    std::string s;
    for (int i = 0; i < n; ++i) s += static_cast<char>('0' + v[i]);
    return s;
}
std::string digits(const int8_t* v, int n) {
    std::string s;
    for (int i = 0; i < n; ++i) s += static_cast<char>('0' + v[i]);
    return s;
}

}  // namespace

// ---------------------------------------------------------------- improved's FT2 code

Halved halve(int nfa, int nfb, int nfqso) {
    Halved h;
    h.nfa = std::max(0, nfa / FT2_STRETCH);
    h.nfb = std::max(h.nfa + 1, nfb / FT2_STRETCH);
    h.nfqso = std::max(0, nfqso / FT2_STRETCH);
    return h;
}

namespace {
template <class S>
void stretch_t(const S* id2, S* id2_ft4) {
    for (int i = 1; i <= FT2_NMAX; ++i) {
        id2_ft4[FT2_STRETCH * i - FT2_STRETCH] = id2[i - 1];
        id2_ft4[FT2_STRETCH * i - 1] = id2[i - 1];
    }
}
}  // namespace

void stretch(const int16_t* id2, int16_t* id2_ft4) { stretch_t(id2, id2_ft4); }
void stretch(const float* id2, float* id2_ft4) { stretch_t(id2, id2_ft4); }

Decode decoded_lines(int nutc, float sync, int snr, float dt, float freq, const std::string& decoded, int nap,
                     float qual) {
    std::string decoded0 = fix(decoded, 37);
    const float dt_ft2 = 0.5f * dt - 0.25f;
    const float freq_ft2 = 2.0f * freq;
    std::string annot = "  ";
    if (nap != 0) {
        // write(annot,'(a1,i1)') 'a',nap: i1 overflows to '*'
        annot = nap >= 0 && nap <= 9 ? std::string("a") + static_cast<char>('0' + nap) : std::string("a*");
        if (qual < 0.17f) decoded0[36] = '?';
    }
    Decode d;
    d.sync = sync;
    d.snr = snr;
    d.dt = dt_ft2;
    d.freq = freq_ft2;
    d.msg37 = decoded0;
    d.nap = nap;
    d.qual = qual;
    // 1001 format(i6.6,i4,f5.1,i5,' * ',1x,a37,1x,a2)
    d.line = fst4w::fortran_i(nutc, 6, 6) + fst4w::fortran_i(snr, 4) + fst4w::fortran_f(dt_ft2, 5, 1) +
             fst4w::fortran_i(nint(freq_ft2), 5) + " * " + " " + decoded0 + " " + annot;
    // 1002 format(i6.6,i4,i5,f6.1,f8.0,i4,3x,a37,' FT2')
    d.decoded_line = fst4w::fortran_i(nutc, 6, 6) + fst4w::fortran_i(nint(sync), 4) + fst4w::fortran_i(snr, 5) +
                     fst4w::fortran_f(dt_ft2, 6, 1) + fst4w::fortran_f(freq_ft2, 8, 0) + fst4w::fortran_i(0, 4) +
                     "   " + decoded0 + " FT2";
    return d;
}

// ---------------------------------------------------------------- the Decoder

Decoder::Decoder()
    : mycall0_(12, ' '),
      hiscall0_(12, ' '),
      dd_(NMAX, 0.0f),
      id2_ft4_(NMAX, 0),
      candidate_(2 * MAXCAND, 0.0f),
      cd2_(NDMAX),
      cb_(NDMAX),
      bitmetrics_(2 * NN * 3, 0.0f),
      frame_(FRAME_BYTES / sizeof(float), 0.0f) {}

Decoder::~Decoder() = default;

size_t Decoder::memory() const {
    auto b = [](const auto& v) { return v.capacity() * sizeof(v[0]); };
    return sizeof(*this) + fe_.heap_bytes() + b(dd_) + b(id2_ft4_) + b(id2_ft4f_) + b(candidate_) + b(cd2_) + b(cb_) + b(frame_) +
           b(bitmetrics_) + b(syncs_) + decoded_text_.capacity() + decoded_pending_.capacity() +
           static_cast<size_t>(NPTS) * sizeof(int16_t);  // period()'s id2
}

void Decoder::set_frontend_trace(std::function<void(const std::string&)> f) {
    fetrace_ = std::move(f);
    fe_.set_trace(fetrace_);
}

// decoder.f90: the RMS gate, decoded.txt, then improved's FT2 branch.
template <class S>
Decoder::Call Decoder::call_t(const S* id2, int nutc, const Params& p) {
    Call c;
    // rms=sqrt(dot_product(float(id2(1:180000)),float(id2(1:180000)))/180000.0)
    float sum = 0.0f;
    for (int i = 0; i < NPTS; ++i) {
        const float x = static_cast<float>(id2[i]);
        sum = sum + x * x;
    }
    const float rms = std::sqrt(sum / 180000.0f);
    if (std::is_same<S, float>::value ? !all_zero(id2, NPTS) : !(rms < 0.5f)) {
        c.gated = true;
        decoded_exists_ = true;
        decoded_pending_.clear();
        wrote_ = 0;
        const int ncontest = p.nexp_decode & 7;
        // if(params%nmode.eq.52): nfa/nfb/nfqso halved, id2(1:23*1728)
        // stretched into id2_ft4, FT4's decode with them.
        const Halved h = halve(p.nfa, p.nfb, p.nfqso);
        S* id2_ft4;
        if constexpr (std::is_same<S, float>::value) {
            if (id2_ft4f_.empty()) id2_ft4f_.assign(NMAX, 0.0f);
            id2_ft4 = id2_ft4f_.data();
        } else {
            id2_ft4 = id2_ft4_.data();
        }
        stretch(id2, id2_ft4);
        if (trace_ || fetrace_) {
            const std::string s = "S " + I(p.nfa) + " " + I(p.nfb) + " " + I(p.nfqso) + " " + I(h.nfa) + " " +
                                  I(h.nfb) + " " + I(h.nfqso) + " " + X(h_i2(id2_ft4, NMAX));
            if (trace_) trace_(s);
            if (fetrace_) fetrace_(s);
        }
        Params q = p;
        q.nfa = h.nfa;
        q.nfb = h.nfb;
        q.nfqso = h.nfqso;
        decode(id2_ft4, q, nutc, ncontest, c);
        if (wrote_ > 0) decoded_text_ = p.nagain ? decoded_text_ + decoded_pending_ : decoded_pending_;
    }
    c.ndecoded = static_cast<int>(c.decodes.size());
    return c;
}

Decoder::Call Decoder::call(const int16_t* id2, int nutc, const Params& p) { return call_t(id2, nutc, p); }
Decoder::Call Decoder::call(const float* id2, int nutc, const Params& p) { return call_t(id2, nutc, p); }

template <class S>
Decoder::Period Decoder::period_t(const S* samples, size_t n, int nutc, const Params& p, size_t take) {
    Period r;
    std::vector<S> id2(NPTS, 0);
    std::copy(samples, samples + std::min({n, take, static_cast<size_t>(NPTS)}), id2.begin());
    r.c = call(id2.data(), nutc, p);
    for (const Decode& d : r.c.decodes) r.out += d.line + "\n";
    r.ndecoded = r.c.ndecoded;
    // decode_completion.f90: ('<DecodeFinished>',2i4,i9)
    r.out += "<DecodeFinished>" + fst4w::fortran_i(0, 4) + fst4w::fortran_i(r.ndecoded, 4) + fst4w::fortran_i(0, 9) +
             "\n";
    return r;
}

// ---------------------------------------------------------------- ft4_decode.f90 (FT2's sizes)

Decoder::Period Decoder::period(const int16_t* samples, size_t n, int nutc, const Params& p, size_t take) {
    return period_t(samples, n, nutc, p, take);
}
Decoder::Period Decoder::period(const float* samples, size_t n, int nutc, const Params& p, size_t take) {
    return period_t(samples, n, nutc, p, take);
}

template <class S>
void Decoder::decode(const S* iwave, const Params& p, int nutc, int ncontest, Call& c) {
    std::string mycall = fix(p.mycall, 12), hiscall = fix(p.hiscall, 12);
    const int nfqso = p.nfqso, ndepth = p.ndepth, nqso = p.nqso_progress;
    const bool lapcqonly = p.lapcqonly;
    pk_.set_calls(mycall, hiscall);  // dxcall13=hiscall; mycall13=mycall
    fe_.set_fast_subtract(p.fast_subtract);
    const bool fe = static_cast<bool>(fetrace_);

    if (first_) {
        const float fs = 12000.0f / static_cast<float>(NDOWN);
        cfloat ctwk[2 * NSS];
        for (int idf = -16; idf <= 16; ++idf) {
            const float a[5] = {static_cast<float>(idf), 0.0f, 0.0f, 0.0f, 0.0f};
            for (cfloat& z : ctwk) z = cfloat(1.0f, 0.0f);
            twkfreq1(ctwk, 2 * NSS, fs / 2.0f, a, ctwk2_[idf + 16]);
        }
        mycall0_.assign(12, ' ');
        hiscall0_.assign(12, ' ');
        first_ = false;
    }

    int l1 = findex(mycall, std::string(1, '\0'));
    if (l1 != 0) std::fill(mycall.begin() + (l1 - 1), mycall.end(), ' ');
    l1 = findex(hiscall, std::string(1, '\0'));
    if (l1 != 0) std::fill(hiscall.begin() + (l1 - 1), hiscall.end(), ' ');
    if (mycall != mycall0_ || hiscall != hiscall0_) {
        std::fill(apbits_, apbits_ + 2 * ND, 0);
        apbits_[0] = 99;
        apbits_[29] = 99;
        std::fill(apmy_ru_, apmy_ru_ + 28, 0);
        std::fill(aphis_fd_, aphis_fd_ + 28, 0);
        if (trim(mycall).size() >= 3) {
            bool nohiscall = false;
            std::string h0 = hiscall;
            if (trim(h0).size() < 3) {
                h0 = mycall;  // use mycall for dummy hiscall - mycall won't be hashed.
                nohiscall = true;
            }
            const std::string message = fix(trim(mycall) + " " + trim(h0) + " RR73", 37);
            const packjt77::Packjt77::Packed pk = pk_.pack77(message);
            const packjt77::Packjt77::Unpacked u = pk_.unpack77(pk.c77.c_str(), 1);
            if (!(pk.i3 != 1 || message != fix(u.msg, 37) || !u.ok)) {
                int8_t message77[77], cw[174];
                for (int i = 0; i < 77; ++i) message77[i] = static_cast<int8_t>(pk.c77[i] - '0');
                for (int i = 0; i < 28; ++i) {
                    apmy_ru_[i] = 2 * ((message77[i] + RVEC[i + 1]) % 2) - 1;
                    aphis_fd_[i] = 2 * ((message77[i + 29] + RVEC[i + 28]) % 2) - 1;
                }
                for (int i = 0; i < 77; ++i) message77[i] = static_cast<int8_t>((message77[i] + RVEC[i]) % 2);
                ft8::encode174_91(message77, cw);
                for (int i = 0; i < 2 * ND; ++i) apbits_[i] = 2 * cw[i] - 1;
                if (nohiscall) apbits_[29] = 99;
            }
        }
        mycall0_ = mycall;
        hiscall0_ = hiscall;
    }
    unsigned char* const frame = reinterpret_cast<unsigned char*>(frame_.data());
    float* const savg = frame_.data() + SAVG_AT / sizeof(float);
    float* const sbase = frame_.data() + SBASE_AT / sizeof(float);
    cfloat* const cd = reinterpret_cast<cfloat*>(frame_.data() + CD_AT / sizeof(float));
    std::memset(frame, ' ', 100 * DLEN);  // decodes=' '
    int ndecodes = 0;
    const float fa = static_cast<float>(p.nfa), fb = static_cast<float>(p.nfb);
    for (int i = 0; i < NMAX; ++i) dd_[i] = static_cast<float>(iwave[i]);
    if (fe)
        fetrace_("C " + I(ndepth) + " " + I(nfqso) + " " + I(p.nfa) + " " + I(p.nfb) + " " + I(ncontest) + " " +
                 I(nqso) + " " + L(lapcqonly) + " " + X(h_i2(iwave, NMAX)) + " " + X(h_i(apbits_, 2 * ND)) + " |" +
                 mycall + "|" + hiscall + "|");

    // ndepth=3: 3 passes, bp+osd; ndepth=2: 3 passes, bp only; ndepth=1: 1
    // pass, no subtraction. syncmin: improved's FT2 threshold (1.18 for FT4).
    const float syncmin = 1.0f;
    bool dosubtract = true, doosd = true;
    int nsp = 3;
    if (ndepth == 2) doosd = false;
    if (ndepth == 1) {
        nsp = 1;
        dosubtract = false;
        doosd = false;
    }
    if (trace_)
        trace_("C " + I(ndepth) + " " + I(nsp) + " " + (dosubtract ? "T" : "F") + " " + (doosd ? "T" : "F"));

    float smax1 = 0.0f;
    int nd1 = 0;
    int nharderror = -1, ntype = 0, iaptype = 0;
    float dmin = 0.0f;
    float llr[2 * ND] = {}, llra[2 * ND], llrb[2 * ND], llrc[2 * ND], llrd[2 * ND] = {};
    int8_t apmask[2 * ND] = {}, message91[91] = {}, cw[2 * ND] = {}, message77[77] = {};
    int i4tone[NN];
    std::string message(37, ' ');

    for (int isp = 1; isp <= nsp; ++isp) {
        if (stopped()) break;  // the caller's watchdog (set_stop): no more passes
        if (isp == 2) {
            if (ndecodes == 0) break;
            nd1 = ndecodes;
        } else if (isp == 3) {
            const int nd2 = ndecodes - nd1;
            if (nd2 == 0) break;
        }

        if (fe) fetrace_("P " + I(isp) + " " + X(h_r(dd_.data(), NMAX)));
        std::fill(candidate_.begin(), candidate_.end(), 0.0f);
        int ncand = 0;
        fe_.getcandidates4(dd_.data(), fa, fb, syncmin, nfqso, MAXCAND, savg, candidate_.data(), ncand, sbase);
        if (trace_) trace_("P " + I(isp) + " " + I(ncand));
        bool dobigfft = true;
        for (int icand = 1; icand <= ncand; ++icand) {
            if (stopped()) break;
            const float f0 = candidate_[(icand - 1) * 2];
            const float snr = candidate_[(icand - 1) * 2 + 1] - 1.0f;
            if (fe) fetrace_("E " + I(icand) + " " + XF(f0) + " " + XF(snr) + " " + L(dobigfft));
            fe_.ft4_downsample(dd_.data(), dobigfft, f0, cd2_.data());  // Downsample to 32 Sam/Sym
            if (dobigfft) dobigfft = false;
            {
                const float sum2 = sum2_cd2(cd2_.data());
                if (sum2 > 0.0f) normalise(cd2_.data(), sum2);
                if (fe) fetrace_("e " + XF(sum2) + " " + X(h_c(cd2_.data(), NDMAX)));
            }
            fe_.sync4d_prepare(cd2_.data());  // for sync4d_many
            int ibest = -1, idfbest = 0;
            float smax = -99.0f;
            for (int iseg = 1; iseg <= 3; ++iseg) {  // DT search is done over 3 segments
                for (int isync = 1; isync <= 2; ++isync) {
                    int idfmin, idfmax, idfstp, ibmin, ibmax, ibstp;
                    if (isync == 1) {
                        idfmin = -12;
                        idfmax = 12;
                        idfstp = 3;
                        ibmin = -344;
                        ibmax = 1012;
                        if (iseg == 1) {
                            ibmin = 108;
                            ibmax = 560;
                        } else if (iseg == 2) {
                            smax1 = smax;
                            ibmin = 560;
                            ibmax = 1012;
                        } else if (iseg == 3) {
                            ibmin = -344;
                            ibmax = 108;
                        }
                        ibstp = 4;
                    } else {
                        idfmin = idfbest - 4;
                        idfmax = idfbest + 4;
                        idfstp = 1;
                        ibmin = ibest - 5;
                        ibmax = ibest + 5;
                        ibstp = 1;
                    }
                    ibest = -1;
                    idfbest = 0;
                    smax = -99.0f;
                    uint32_t hy = 2166136261u;
                    const int nstart = std::max(0, (ibmax - ibmin + ibstp) / ibstp);
                    for (int idf = idfmin; idf <= idfmax; idf += idfstp) {
                        fe_.sync4d_many(cd2_.data(), ibmin, nstart, ibstp, ctwk2_[idf + 16], syncs_.data());
                        for (int k = 0; k < nstart; ++k) {
                            const int istart = ibmin + k * ibstp;
                            const float sync = syncs_[k];  // sync power
                            if (fe) hy = hw(hy, fbits(sync));
                            if (sync > smax) {
                                smax = sync;
                                ibest = istart;
                                idfbest = idf;
                            }
                        }
                    }
                    if (fe)
                        fetrace_("Y " + I(iseg) + " " + I(isync) + " " + X(hy) + " " + XF(smax) + " " + I(ibest) +
                                 " " + I(idfbest));
                }
                if (iseg == 1) smax1 = smax;
                if (smax < 1.2f) continue;
                if (iseg > 1 && smax < smax1) continue;
                const float f1 = f0 + static_cast<float>(idfbest);
                if (f1 <= 10.0f || f1 >= 4990.0f) continue;
                if (fe) fetrace_("G " + I(iseg) + " " + XF(f1) + " " + L(dobigfft));
                fe_.ft4_downsample(dd_.data(), dobigfft, f1, cb_.data());  // Final downsample, corrected f0
                float sum2b;
                {
                    sum2b = sum2_cb(cb_.data());
                    if (sum2b > 0.0f) normalise(cb_.data(), sum2b);
                }
                std::fill(cd, cd + NN * NSS, cfloat(0.0f, 0.0f));
                if (ibest >= 0) {
                    const int it = std::min(NDMAX - 1, ibest + NN * NSS - 1);
                    const int np = it - ibest + 1;
                    for (int k = 0; k < np; ++k) cd[k] = cb_[ibest + k];
                } else {
                    for (int k = 0; k < NN * NSS + 2 * ibest; ++k) cd[-ibest + k] = cb_[k];
                }
                bool smp = false;
                if (fe) {
                    fetrace_("g " + XF(sum2b) + " " + X(h_c(cb_.data(), NDMAX)) + " " + X(h_c(cd, NN * NSS)));
                    ++nbm_;
                    smp = nbm_ % FE_NSAMP == 1;
                    if (smp) {
                        std::string l = "I";
                        l.reserve(1 + 9 * 2 * NN * NSS);
                        for (int k = 0; k < NN * NSS; ++k) l += " " + hex8(cd[k].real()) + " " + hex8(cd[k].imag());
                        fetrace_(l);
                    }
                }
                bool badsync = false;
                fe_.get_ft4_bitmetrics(cd, bitmetrics_.data(), badsync);
                if (fe && smp && !badsync) {
                    std::string l = "J";
                    l.reserve(1 + 9 * 3 * 2 * NN);
                    for (int k = 0; k < 3 * 2 * NN; ++k) l += " " + hex8(bitmetrics_[k]);
                    fetrace_(l);
                }
                if (badsync) continue;
                const int nsync_qual = ft4::ft4_nsync_qual(bitmetrics_.data());
                if (fe) fetrace_("H " + I(nsync_qual));
                if (nsync_qual < 20) continue;

                ft4::ft4_llrs(bitmetrics_.data(), llra, llrb, llrc);  // scalefac 2.83

                float absa[2 * ND];
                for (int i = 0; i < 2 * ND; ++i) absa[i] = std::fabs(llra[i]);
                const float apmag = maxval(absa, 2 * ND) * 1.1f;
                int npasses = 3 + NAPPASSES[nqso];
                if (lapcqonly) npasses = 4;
                if (ndepth == 1) npasses = 3;
                if (ncontest >= 6) npasses = 3;  // Don't support Fox and Hound
                for (int ipass = 1; ipass <= npasses; ++ipass) {
                    if (ipass == 1) std::memcpy(llr, llra, sizeof llr);
                    if (ipass == 2) std::memcpy(llr, llrb, sizeof llr);
                    if (ipass == 3) std::memcpy(llr, llrc, sizeof llr);
                    if (ipass <= 3) {
                        std::memset(apmask, 0, sizeof apmask);
                        iaptype = 0;
                    }

                    if (ipass > 3) {
                        std::memcpy(llrd, llrc, sizeof llrd);
                        iaptype = NAPTYPES[nqso][ipass - 3 - 1];
                        if (lapcqonly) iaptype = 1;

                        // Conditions that cause us to bail out of AP decoding
                        // (napwid=50: folded in the reference's object)
                        if (ncontest <= 5 && iaptype >= 3 && std::fabs(f1 - static_cast<float>(nfqso)) > 50.0f)
                            continue;
                        if (iaptype >= 2 && apbits_[0] > 1) continue;   // No, or nonstandard, mycall
                        if (iaptype >= 3 && apbits_[29] > 1) continue;  // No, or nonstandard, dxcall

                        auto apm = [&](int v) { return apmag * static_cast<float>(v); };
                        if (iaptype == 1) {  // CQ or CQ TEST or CQ FD or CQ RU or CQ WW
                            std::memset(apmask, 0, sizeof apmask);
                            for (int i = 0; i < 29; ++i) apmask[i] = 1;
                            const Ap29* t = nullptr;
                            if (ncontest == 0) t = &mcq();
                            if (ncontest == 1) t = &mcqtest();
                            if (ncontest == 2) t = &mcqtest();
                            if (ncontest == 3) t = &mcqfd();
                            if (ncontest == 4) t = &mcqru();
                            if (ncontest == 5) t = &mcqww();
                            if (t)
                                for (int i = 0; i < 29; ++i) llrd[i] = apm(t->v[i]);
                        }

                        if (iaptype == 2) {  // MyCall,???,???
                            std::memset(apmask, 0, sizeof apmask);
                            if (ncontest == 0 || ncontest == 1 || ncontest == 5) {
                                for (int i = 0; i < 29; ++i) apmask[i] = 1;
                                for (int i = 0; i < 29; ++i) llrd[i] = apm(apbits_[i]);
                            } else if (ncontest == 2 || ncontest == 3) {
                                for (int i = 0; i < 28; ++i) apmask[i] = 1;
                                for (int i = 0; i < 28; ++i) llrd[i] = apm(apbits_[i]);
                            } else if (ncontest == 4) {
                                for (int i = 1; i < 29; ++i) apmask[i] = 1;
                                for (int i = 1; i < 29; ++i) llrd[i] = apm(apmy_ru_[i - 1]);
                            }
                        }

                        if (iaptype == 3) {  // MyCall,DxCall,???
                            std::memset(apmask, 0, sizeof apmask);
                            if (ncontest == 0 || ncontest == 1 || ncontest == 2 || ncontest == 5) {
                                for (int i = 0; i < 58; ++i) apmask[i] = 1;
                                for (int i = 0; i < 58; ++i) llrd[i] = apm(apbits_[i]);
                            } else if (ncontest == 3) {  // Field Day
                                for (int i = 0; i < 56; ++i) apmask[i] = 1;
                                for (int i = 0; i < 28; ++i) llrd[i] = apm(apbits_[i]);
                                for (int i = 28; i < 56; ++i) llrd[i] = apm(aphis_fd_[i - 28]);
                            } else if (ncontest == 4) {
                                for (int i = 1; i < 57; ++i) apmask[i] = 1;
                                for (int i = 1; i < 29; ++i) llrd[i] = apm(apmy_ru_[i - 1]);
                                for (int i = 29; i < 57; ++i) llrd[i] = apm(apbits_[i]);
                            }
                        }

                        if (iaptype == 4 || iaptype == 5 || iaptype == 6) {
                            std::memset(apmask, 0, sizeof apmask);
                            if (ncontest <= 5) {
                                for (int i = 0; i < 77; ++i) apmask[i] = 1;  // mycall, hiscall, RRR|73|RR73
                                if (iaptype == 6)
                                    for (int i = 0; i < 77; ++i) llrd[i] = apm(apbits_[i]);
                            }
                        }

                        std::memcpy(llr, llrd, sizeof llr);
                    }
                    std::memset(message77, 0, sizeof message77);
                    dmin = 0.0f;

                    const int ndeep = 2;
                    int maxosd = 2;
                    if (std::fabs(static_cast<float>(nfqso) - f1) <= 50.0f) maxosd = 3;
                    if (!doosd) maxosd = -1;
                    ft8::decode174_91(llr, 91, maxosd, ndeep, apmask, message91, cw, ntype, nharderror, dmin);
                    std::memcpy(message77, message91, 77);

                    int sum77 = 0;
                    for (int i = 0; i < 77; ++i) sum77 += message77[i];
                    if (sum77 == 0) {
                        if (trace_ && nharderror >= 0)
                            trace_("Z " + I(isp) + " " + I(icand) + " " + I(iseg) + " " + I(ipass) + " " +
                                   I(nharderror));
                        continue;
                    }
                    if (nharderror >= 0) {
                        for (int i = 0; i < 77; ++i)
                            message77[i] = static_cast<int8_t>((message77[i] + RVEC[i]) % 2);  // remove rvec scrambling
                        char c77[78];
                        for (int i = 0; i < 77; ++i) c77[i] = static_cast<char>('0' + message77[i]);
                        c77[77] = '\0';
                        const packjt77::Packjt77::Unpacked u = pk_.unpack77(c77, 1);
                        message = fix(u.msg, 37);
                        if (!u.ok) {
                            if (trace_)
                                trace_("F " + I(isp) + " " + I(icand) + " " + I(iseg) + " " + I(ipass) + " " +
                                       I(iaptype) + " |" + std::string(c77, 77) + "|");
                            break;
                        }
                        if (dosubtract) {
                            int8_t m77in[77];
                            if (fe) std::memcpy(m77in, message77, sizeof m77in);
                            ft4::get_ft4_tones_from_77bits(message77, i4tone);
                            if (fe) fetrace_("O " + digits(m77in, 77) + " " + digits(i4tone, NN) + " " + digits(message77, 77));
                            const float dt = static_cast<float>(ibest) / 666.67f;
                            if (trace_) trace_("U " + hex8(f1) + " " + hex8(dt));
                            fe_.subtractft4(dd_.data(), i4tone, f1, dt);
                        }
                        auto dline = [&](const std::string& qual, const char* what) {
                            const int n3 = 4 * (c77[71] - '0') + 2 * (c77[72] - '0') + (c77[73] - '0');
                            const int i3 = 4 * (c77[74] - '0') + 2 * (c77[75] - '0') + (c77[76] - '0');
                            trace_("D " + I(isp) + " " + I(icand) + " " + I(iseg) + " " + I(ipass) + " " +
                                   I(iaptype) + " " + I(nharderror) + " " + I(ntype) + " " + hex8(dmin) + " " +
                                   I(ibest) + " " + I(idfbest) + " " + hex8(f1) + " " + hex8(smax) + " " + qual +
                                   " " + I(i3) + " " + I(n3) + " " + what + " |" + message + "|");
                        };
                        bool idupe = false;
                        for (int i = 0; i < ndecodes; ++i)
                            if (std::memcmp(frame + static_cast<size_t>(i) * DLEN, message.data(), DLEN) == 0)
                                idupe = true;
                        if (idupe) {
                            if (trace_) dline("-", "dupe");
                            break;
                        }
                        ++ndecodes;
                        std::memcpy(frame + static_cast<size_t>(ndecodes - 1) * DLEN, message.data(), DLEN);
                        float xsnr;
                        if (snr > 0.0f)
                            xsnr = 10.0f * ft8::log10f_glibc(snr) - 14.8f;
                        else
                            xsnr = -21.0f;
                        const int nsnr = nint(std::max(-21.0f, xsnr));
                        const float xdt = static_cast<float>(ibest) / 666.67f - 0.5f;
                        const float qual = 1.0f - (static_cast<float>(nharderror) + dmin) / 60.0f;
                        if (trace_) dline(hex8(qual), "new");
                        // decoder_callbacks.f90's ft2_decoded (improved's)
                        Decode d = decoded_lines(nutc, smax, nsnr, xdt, f1, message, iaptype, qual);
                        decoded_pending_ += d.decoded_line + "\n";
                        ++wrote_;
                        c.decodes.push_back(std::move(d));
                        break;
                    }
                }  // Sequence estimation
                if (fe) fetrace_("R " + I(iseg) + " " + I(nharderror));
                if (nharderror >= 0) break;
            }  // 3 DT segments
        }  // Candidate list
    }  // Subtraction loop
    if (trace_) trace_("N " + I(ndecodes));
}

// ---------------------------------------------------------------- the harness

int run_periods(const std::vector<std::string>& args, const std::string& cwd, RunOutput& r,
                const std::function<void(const std::string&)>& fetrace) {
    auto at = [&](const std::string& p) { return !p.empty() && p[0] == '/' ? p : cwd + "/" + p; };
    auto usage = [&]() {
        r.err += "usage: ft2-oracle [-2] -d 1|2|3 WAV...\n"
                 "       ft2-oracle --file [-2] [-d D -L HZ -H HZ -f HZ -F HZ -c CALL -x CALL -g GRID -G GRID -Q N "
                 "-X N -p S] WAV...\n";
        return 2;
    };
    // The harness's options (and jt9's, for "jt9 -a DIR --ft2 ...").
    auto ends = [](const std::string& a, const std::string& t) {
        return a.size() >= t.size() && a.compare(a.size() - t.size(), t.size(), t) == 0;
    };
    bool filemode = !args.empty() && (args[0] == "jt9" || ends(args[0], "/jt9") || args[0] == "ft2-jt9" ||
                                      ends(args[0], "/ft2-jt9"));
    bool ft2 = false;
    int ndepth = -1, flow = 200, fhigh = 4000, nrxfreq = 1500, nqso = 0, nexp_decode = 0;
    std::string mycall = "K1ABC", hiscall = "W9XYZ";
    size_t iarg = 1;
    auto word = [](const std::string& a) {
        const size_t e = a.find_first_of(" ,/");
        return a.substr(0, e);
    };
    auto num = [&](int& v) {
        if (++iarg >= args.size()) return false;
        try {
            v = std::stoi(args[iarg]);
        } catch (...) {
            return false;
        }
        return true;
    };
    for (; iarg < args.size(); ++iarg) {
        const std::string& a = args[iarg];
        if (a.empty() || a[0] != '-') break;
        if (a == "--file") {
            filemode = true;
        } else if (a == "-2" || a == "--ft2") {
            ft2 = true;
        } else if (a == "-a") {
            if (++iarg >= args.size()) return usage();
        } else if (a == "-d") {
            if (!num(ndepth)) return usage();
        } else if (a == "-L") {
            if (!num(flow)) return usage();
        } else if (a == "-H") {
            if (!num(fhigh)) return usage();
        } else if (a == "-f") {
            if (!num(nrxfreq)) return usage();
        } else if (a == "-F") {
            int ntol;
            if (!num(ntol)) return usage();  // (FT2 does not read ntol)
        } else if (a == "-Q") {
            if (!num(nqso)) return usage();
        } else if (a == "-X") {
            if (!num(nexp_decode)) return usage();
        } else if (a == "-p") {
            if (++iarg >= args.size()) return usage();
            const std::string& t = args[iarg];
            if (!real_token(t) || !std::isfinite(std::strtod(t.c_str(), nullptr))) {
                r.err += "jt9: invalid value for -p: " + trim(t) + " (must be a finite number)\n";
                return 2;
            }
            const double trperiod = std::strtod(t.c_str(), nullptr);
            if (trperiod <= 0.0 || trperiod > 1800.0) {
                r.err += "jt9: invalid value for -p: " + trim(t) + " (must be between 0 and 1800 seconds)\n";
                return 2;
            }
        } else if (a == "-c" || a == "-x" || a == "-G" || a == "-g") {
            if (++iarg >= args.size()) return usage();
            const std::string v = word(args[iarg]);
            if (a == "-c") mycall = v;
            if (a == "-x") hiscall = v;
        } else {
            return usage();
        }
    }
    if (!ft2) return usage();  // (this is FT2's harness)
    if (!filemode && (ndepth < 1 || ndepth > 3)) return usage();
    if (filemode && ndepth < 0) ndepth = 1;
    if (iarg >= args.size()) return usage();

    Params p;
    p.ndepth = ndepth;
    if (filemode) {
        // jt9_params_init's for mode 52 (mode 5's)
        if (mycall == "b") mycall.clear();
        if (hiscall == "b") hiscall.clear();
        p.nfa = flow;
        p.nfb = fhigh;
        p.nfqso = nrxfreq;
        p.nqso_progress = nqso;
        p.nexp_decode = nexp_decode;
        p.mycall = mycall;
        p.hiscall = hiscall;
    }
    if (p.nqso_progress < 0 || p.nqso_progress > 5) return usage();  // nappasses(0:5)
    Decoder dec;
    if (r.want_trace) dec.set_trace([&r](const std::string& l) { r.trace += l + "\n"; });
    if (fetrace) dec.set_frontend_trace(fetrace);
    auto stop2 = [&]() {
        r.decoded_exists = dec.decoded_exists();
        r.decoded = dec.decoded_text();
        return 2;
    };
    for (; iarg < args.size(); ++iarg) {
        const std::string& infile = args[iarg];
        std::vector<int16_t> s;
        int rate = 0;
        std::string why;
        if (!ft8::read_wav(at(infile), s, rate, why)) {
            r.err += (filemode ? "jt9" : "ft2-oracle") + std::string(": cannot read WAV file ") + infile + ": " +
                     why + "\n";
            return stop2();
        }
        if (rate != 12000) {
            if (filemode)
                r.err += "jt9: unsupported " + std::to_string(rate) + " Hz WAV sample rate for the selected mode\n";
            else
                r.err += "ft2-oracle: " + infile + " is not 12000 Hz\n";
            return stop2();
        }
        const int nutc = ft8::nutc_of(infile);
        const std::string w = "W " + fst4w::fortran_i(nutc, 6, 6) + " " + I(ndepth) + " " + trim(infile);
        if (r.want_trace) r.trace += w + "\n";
        if (fetrace) fetrace(w);
        if (filemode) {
            // jt9 --ft2 reads npts = 21*1728 in blocks of 3456: 10 blocks,
            // 34,560 samples, and says so when the file ends first.
            const size_t want = (21 * 1728 / 3456) * 3456;
            if (s.size() < want) r.out += " EOF on input file " + trim(infile) + "\n";
            const Decoder::Period per = dec.period(s.data(), s.size(), nutc, p, want);
            r.out += per.out;
        } else {
            // the GUI's request: the cycle's 45,000 samples (3.75 s).
            const Decoder::Period per = dec.period(s.data(), s.size(), nutc, p, 45000);
            r.out += per.out;
            r.decoded_calls += "== " + fst4w::fortran_i(nutc, 6, 6) + " 23\n";
            if (dec.decoded_exists()) {
                size_t a = 0;
                const std::string& t = dec.decoded_text();
                while (a < t.size()) {
                    size_t e = t.find('\n', a);
                    if (e == std::string::npos) e = t.size();
                    r.decoded_calls += trim(t.substr(a, e - a)) + "\n";
                    a = e + 1;
                }
            }
        }
    }
    r.decoded_exists = dec.decoded_exists();
    r.decoded = dec.decoded_text();
    return 0;
}

}  // namespace ft2
