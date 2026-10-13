// The built-in FT8 decoder (ft8.h): WSJT-X lib/ft8_decode.f90 whole,
// decoder.f90's FT8 branch, decoder_callbacks.f90's ft8_decoded,
// decode_completion.f90, jt9a.f90's 41/47/50 sequence, ft8/ft8_a7.f90's
// ft8_a7_save and getmsg, ft8/ft8_a8d.f90, ft8/twkfreq1.f90, smo121.f90 and
// db.f90 at 567ad29ce, line for line; and scripts/ft8-oracle.f90 (the
// harness that runs them as jt9a and jt9 -8 do) as run_periods.
//
// Joe Taylor, K1JT, Steven Franke, K9AN, and the WSJT-X developers (GNU GPL
// v3), and so is this.
#include "ft8.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>

#include "fst4w_jt9.h"  // fortran_f, fortran_i
#include "ft8_maths.h"  // log10f_glibc, sincosf_glibc, logf_glibc

namespace ft8 {

namespace {

constexpr int MAXCAND = 1000, MAX_EARLY = 200, NPTS = 15 * 12000;

// Fortran's character*n assignment: cut or blank padded.
std::string fix(std::string s, size_t n) {
    s.resize(n, ' ');
    return s;
}

// trim(s): trailing blanks only (NULs stay).
std::string trim(const std::string& s) {
    size_t n = s.size();
    while (n > 0 && s[n - 1] == ' ') --n;
    return s.substr(0, n);
}

// index(s, sub): 1-based, 0 when absent, 1 for an empty sub.
int findex(const std::string& s, const std::string& sub) {
    const size_t p = s.find(sub);
    return p == std::string::npos ? 0 : static_cast<int>(p) + 1;
}

// s(a:b), 1-based, empty when b < a (b cut to the string's length).
std::string substr(const std::string& s, int a, int b) {
    b = std::min(b, static_cast<int>(s.size()));
    if (b < a) return std::string();
    return s.substr(static_cast<size_t>(a - 1), static_cast<size_t>(b - a + 1));
}

// Fortran's nint of a default real.
inline int nint(float x) { return static_cast<int>(std::lround(x)); }

// jt9_input_validation.f90's is_real_token.
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

// A float's bits as the trace's z8.8.
std::string hex8(float v) {
    uint32_t u;
    std::memcpy(&u, &v, 4);
    char b[12];
    std::snprintf(b, sizeof b, "%08X", u);
    return b;
}

bool blank(const std::string& s) { return s.find_first_not_of(' ') == std::string::npos; }

// Fortran's complex product (no NaN recovery), as the front end's.
inline cfloat cmul(cfloat a, cfloat b) {
    return cfloat(a.real() * b.real() - a.imag() * b.imag(), a.real() * b.imag() + a.imag() * b.real());
}

// db.f90
float db(float x) {
    float d = -99.0f;
    if (x > 1.259e-10f) d = 10.0f * log10f_glibc(x);
    return d;
}

// smo121.f90 on x(1:nz).
void smo121(float* x, int nz) {
    float x0 = x[0];
    for (int i = 2; i <= nz - 1; ++i) {
        const float x1 = x[i - 1];
        x[i - 1] = 0.5f * x[i - 1] + 0.25f * (x0 + x[i]);
        x0 = x1;
    }
}

// twkfreq1.f90 (ca and cb may be the same array, as ft8_a8d passes them).
void twkfreq1(const cfloat* ca, int npts, float fsample, const float a[5], cfloat* cb) {
    const float twopi = 6.283185307f;
    cfloat w(1.0f, 0.0f);
    const float x0 = 0.5f * static_cast<float>(npts + 1);
    const float s = 2.0f / static_cast<float>(npts);
    for (int i = 1; i <= npts; ++i) {
        const float x = s * (static_cast<float>(i) - x0);
        const float p2 = 1.5f * x * x - 0.5f;
        const float p3 = 2.5f * (x * x * x) - 1.5f * x;
        const float p4 = 4.375f * ((x * x) * (x * x)) - 3.75f * (x * x) + 0.375f;
        const float dphi = (a[0] + x * a[1] + p2 * a[2] + p3 * a[3] + p4 * a[4]) * (twopi / fsample);
        float sn, cs;
        sincosf_glibc(dphi, &sn, &cs);
        w = cmul(w, cfloat(cs, sn));
        cb[i - 1] = cmul(w, ca[i - 1]);
    }
}

// (i3.2)
std::string i32(int v) {
    char b[8];
    std::snprintf(b, sizeof b, "%s%02d", v < 0 ? "-" : "", v < 0 ? -v : v);
    std::string s = b;
    if (s.size() < 3) s = std::string(3 - s.size(), ' ') + s;
    return s.size() > 3 ? std::string("***") : s;
}

// msg(j0:j0+len-1) = v on a character*37.
void put(std::string& m, int j0, const std::string& v) {
    for (size_t q = 0; q < v.size() && j0 - 1 + q < m.size(); ++q) m[j0 - 1 + q] = v[q];
}

// ft8_a7.f90's getmsg(i, mycall, hiscall, hisgrid, msg).
std::string getmsg(int i, const std::string& mycall, const std::string& hiscall, const std::string& hisgrid) {
    const bool my_std = packjt77::stdcall(mycall);
    const bool his_std = packjt77::stdcall(hiscall);
    int isnr = 0;
    std::string msg = fix(trim(mycall) + " " + trim(hiscall), 37);
    if (!my_std) {
        if (i == 1 || i >= 6) msg = fix("<" + trim(mycall) + "> " + trim(hiscall), 37);
        if (i >= 2 && i <= 4) msg = fix(trim(mycall) + " <" + trim(hiscall) + ">", 37);
    } else if (!his_std) {
        if (i <= 4 || i == 6) msg = fix("<" + trim(mycall) + "> " + trim(hiscall), 37);
        if (i >= 7) msg = fix(trim(mycall) + " <" + trim(hiscall) + ">", 37);
    }
    const int j0 = static_cast<int>(trim(msg).size()) + 2;
    if (i == 2) put(msg, j0, "RRR");
    if (i == 3) put(msg, j0, "RR73");
    if (i == 4) put(msg, j0, "73");
    if (i == 5) {
        if (his_std) msg = fix("CQ " + trim(hiscall) + " " + hisgrid.substr(0, 4), 37);
        if (!his_std) msg = fix("CQ " + trim(hiscall), 37);
    }
    if (i == 6 && his_std) put(msg, j0, hisgrid.substr(0, 4));
    if (i >= 7 && i <= 206) {
        isnr = -50 + (i - 7) / 2;
        if ((i & 1) == 1) {
            put(msg, j0, i32(isnr));
            if (msg[j0 - 1] == ' ') msg[j0 - 1] = '+';
        } else {
            put(msg, j0, "R" + i32(isnr));
            if (msg[j0] == ' ') msg[j0] = '+';
        }
    }
    if (std::abs(isnr) > 30) msg = std::string(37, ' ');
    return msg;
}

// isgrid4 (ft8_a7_save's statement function) on g4, its character*4.
bool isgrid4(const std::string& w) {
    const std::string g4 = fix(w, 4);
    return trim(g4).size() == 4 && g4[0] >= 'A' && g4[0] <= 'R' && g4[1] >= 'A' && g4[1] <= 'R' && g4[2] >= '0' &&
           g4[2] <= '9' && g4[3] >= '0' && g4[3] <= '9';
}

}  // namespace

// packjt77 for the front end: the decoder's own Packjt77.
class Decoder::PackCodec : public Codec {
public:
    explicit PackCodec(packjt77::Packjt77& p) : p_(p) {}
    void set_calls(const std::string& my, const std::string& his) override { p_.set_calls(my, his); }
    bool unpack77(const char* c77, int nrx, std::string& msg37) override {
        const packjt77::Packjt77::Unpacked u = p_.unpack77(c77, nrx);
        msg37 = fix(u.msg, 37);
        return u.ok;
    }
    bool stdcall(const std::string& call12) override { return packjt77::stdcall(call12); }
    Packed genft8(const std::string& msg37) override {
        const packjt77::Packjt77::Genft8 g = p_.genft8(msg37);
        Packed r;
        r.ok = g.ok;
        r.msgsent = g.msgsent;
        for (int i = 0; i < 77; ++i) r.msgbits[i] = static_cast<int8_t>(g.msgbits[i]);
        return r;
    }

private:
    packjt77::Packjt77& p_;
};

Decoder::Decoder()
    : codec_(std::make_unique<PackCodec>(pk_)),
      dd_(NPTS, 0.0f),
      dd1_(NPTS, 0.0f),
      itone_save_(static_cast<size_t>(NN) * MAX_EARLY, 0),
      allmessages_(MAX_EARLY, std::string(37, ' ')),
      sbase_(NH1, 0.0f),
      candidate_(3 * MAXCAND, 0.0f),
      a7mem_(A7_END, 0),
      a8_cd0_(3200),
      a8_csymb_(32) {}

Decoder::~Decoder() = default;

size_t Decoder::memory() const {
    auto b = [](const auto& v) { return v.capacity() * sizeof(v[0]); };
    size_t n = sizeof(*this) + fe_.heap_bytes() + b(dd_) + b(dd1_) + b(itone_save_) + b(sbase_) + b(candidate_) +
               b(a7mem_) + b(a8_cd0_) + b(a8_csymb_) + decoded_text_.capacity() + decoded_pending_.capacity();
    for (const std::string& m : allmessages_) n += sizeof m + m.capacity();
    return n + 2 * static_cast<size_t>(NPTS) * sizeof(int16_t);  // period()'s id2 and id2a
}

// ---------------------------------------------------------------- ft8_a7's storage

float Decoder::a7_real(size_t base, int e) const {
    const size_t o = base + 4 * static_cast<size_t>(e);
    float v = 0.0f;
    if (o + 4 <= A7_END) std::memcpy(&v, &a7mem_[o], 4);
    return v;
}

void Decoder::set_a7_real(size_t base, int e, float v) {
    const size_t o = base + 4 * static_cast<size_t>(e);
    if (o + 4 <= A7_END) std::memcpy(&a7mem_[o], &v, 4);
}

std::string Decoder::msg0(int e) const {
    std::string m(37, '\0');
    const size_t o = A7_MSG0 + 37 * static_cast<size_t>(e);
    for (size_t q = 0; q < 37 && o + q < A7_END; ++q) m[q] = static_cast<char>(a7mem_[o + q]);
    return m;
}

void Decoder::set_msg0(int e, const std::string& m) {
    const size_t o = A7_MSG0 + 37 * static_cast<size_t>(e);
    for (size_t q = 0; q < 37 && o + q < A7_END; ++q) a7mem_[o + q] = static_cast<unsigned char>(m[q]);
}

int Decoder::jseq() const {
    int32_t v;
    std::memcpy(&v, &a7mem_[A7_JSEQ], 4);
    return v;
}

void Decoder::set_jseq(int v) {
    const int32_t w = v;
    std::memcpy(&a7mem_[A7_JSEQ], &w, 4);
}

// ---------------------------------------------------------------- decoder.f90

Decoder::Call Decoder::call(const int16_t* id2, int nzhsym, int nutc, const Params& p) {
    Call c;
    // rms=sqrt(dot_product(float(id2(1:180000)),float(id2(1:180000)))/180000.0)
    float sum = 0.0f;
    for (int i = 0; i < 180000; ++i) {
        const float x = static_cast<float>(id2[i]);
        sum = sum + x * x;
    }
    const float rms = std::sqrt(sum / 180000.0f);
    if (!(rms < 0.5f)) {
        c.gated = true;
        // FT8 defaults at start to enable a8 decoding
        if (nzhsym == 41) ltry_a8_ = true;
        decoded_exists_ = true;  // open(13, ..., status='unknown')
        decoded_pending_.clear();
        wrote_ = 0;
        const int ncontest = p.nexp_decode & 7;
        if (!(p.nfa > p.nfb)) decode(id2, p, nutc, nzhsym, ncontest, c);
        if (wrote_ > 0) decoded_text_ = decoded_pending_;  // close(13): truncated after the writes
    }
    // 800: ndecoded
    int ndecoded = static_cast<int>(c.decodes.size());
    if (nzhsym == 41)
        ndec41_ = ndecoded;
    else if (nzhsym == 47)
        ndec47_ = ndecoded;
    else if (nzhsym == 50)
        ndecoded = ndec41_ + ndec47_ + ndecoded;
    c.ndecoded = ndecoded;
    return c;
}

Decoder::Period Decoder::period(const int16_t* samples, size_t n, int nutc, const Params& p, size_t fifty,
                                const std::function<void(int)>& after_call) {
    Period r;
    std::vector<int16_t> id2(NPTS, 0), id2a(NPTS, 0);
    std::copy(samples, samples + std::min<size_t>(n, NPTS), id2.begin());
    auto lines = [&](const Call& c) {
        for (const Decode& d : c.decodes) r.out += d.line + "\n";
    };
    // lib/jt9a.f90: the two early calls on id2a, the third on id2.
    for (int nearly : {41, 47}) {
        const size_t k = static_cast<size_t>(nearly) * 3456;
        std::fill(id2a.begin(), id2a.end(), 0);
        std::copy(id2.begin(), id2.begin() + static_cast<long>(k), id2a.begin());
        Call c = call(id2a.data(), nearly, nutc, p);
        lines(c);
        if (after_call) after_call(nearly);
        (nearly == 41 ? r.c41 : r.c47) = std::move(c);
    }
    if (fifty < static_cast<size_t>(NPTS)) {
        std::fill(id2a.begin(), id2a.end(), 0);
        std::copy(id2.begin(), id2.begin() + static_cast<long>(fifty), id2a.begin());
        r.c50 = call(id2a.data(), 50, nutc, p);
    } else {
        r.c50 = call(id2.data(), 50, nutc, p);
    }
    lines(r.c50);
    if (after_call) after_call(50);
    r.ndecoded = r.c50.ndecoded;
    // decode_completion.f90: ('<DecodeFinished>',2i4,i9)
    r.out += "<DecodeFinished>" + fst4w::fortran_i(0, 4) + fst4w::fortran_i(r.ndecoded, 4) + fst4w::fortran_i(0, 9) +
             "\n";
    return r;
}

// decoder_callbacks.f90's ft8_decoded (ncontest 6's Fox branch not reached).
void Decoder::callback(float sync, int snr, float dt, float freq, const std::string& decoded, int nap, float qual,
                       int nutc, int nfqso, Call& c) {
    std::string decoded0 = fix(decoded, 37);
    const float fdiff = freq - static_cast<float>(nfqso);
    if (std::fabs(fdiff) < 3.0f) ltry_a8_ = false;
    std::string annot = "  ";
    if (nap != 0) {
        annot = nap >= 0 && nap <= 9 ? std::string("a") + static_cast<char>('0' + nap) : std::string("a*");
        if (qual < 0.17f) decoded0[36] = '?';
    }
    Decode d;
    d.sync = sync;
    d.snr = snr;
    d.dt = dt;
    d.freq = freq;
    d.msg37 = decoded0;
    d.nap = nap;
    d.qual = qual;
    // 1001 format(i6.6,i4,f5.1,i5,' ~ ',1x,a37,1x,a2)
    d.line = fst4w::fortran_i(nutc, 6, 6) + fst4w::fortran_i(snr, 4) + fst4w::fortran_f(dt, 5, 1) +
             fst4w::fortran_i(nint(freq), 5) + " ~ " + " " + decoded0 + " " + annot;
    // 1002 format(i6.6,i4,i5,f6.1,f8.0,i4,3x,a37,' FT8')
    d.decoded_line = fst4w::fortran_i(nutc, 6, 6) + fst4w::fortran_i(nint(sync), 4) + fst4w::fortran_i(snr, 5) +
                     fst4w::fortran_f(dt, 6, 1) + fst4w::fortran_f(freq, 8, 0) + fst4w::fortran_i(0, 4) + "   " +
                     decoded0 + " FT8";
    decoded_pending_ += d.decoded_line + "\n";
    ++wrote_;
    c.decodes.push_back(std::move(d));
}

// ---------------------------------------------------------------- ft8_decode.f90

void Decoder::decode(const int16_t* iwave16, const Params& p, int nutc, int nzhsym, int ncontest, Call& c) {
    const std::string mycall12 = fix(p.mycall, 12), hiscall12 = fix(p.hiscall, 12), hisgrid = fix(p.hisgrid, 6);
    const int ndepth = p.ndepth, nfqso = p.nfqso, nftx = p.nftx;
    const bool lft8apon = p.lft8apon, nagain = p.nagain;
    auto done = [&]() {
        if (after_decode_) after_decode_(nutc, nzhsym, pk_);
    };
    auto load = [&](std::vector<float>& d) {  // dd=iwave
        for (int i = 0; i < NPTS; ++i) d[i] = static_cast<float>(iwave16[i]);
    };

    bool la8 = true;
    if (nzhsym == 41) set_jseq((nutc / 5) % 2);
    if (nutc0_ == -1) {
        for (int e = 0; e < 800; ++e) {
            set_msg0(e, std::string(37, ' '));
            set_dt0(e, 0.0f);
            set_f0(e, 0.0f);
        }
    }
    // jseq as the statements below read it. Only a rotation of more than 600
    // entries (600 saves in one period) writes over it (msg0's copy past the
    // array's end); then the reference indexes ndec with what was copied
    // there, which is not emulated (ft8.h, a7mem_).
    const int jseq_ = jseq();

    // Added 41==nzhsym to force a reset if the same wav file is processed
    // twice or more in a row, in which case nutc.eq.nutc0 and ndec(jseq,1)
    // doesn't get reset
    if (nzhsym == 41 || nutc != nutc0_) {
        // New UTC.  Move previously saved 'a7' data from k=1 to k=0 (past
        // 200 entries the copy runs on into the next columns and, for jseq
        // 1, past the arrays' end, as the reference's).
        const int iz = ndec_[jseq_][1];
        for (int i = 1; i <= iz; ++i) set_dt0(a7(i, jseq_, 0), dt0(a7(i, jseq_, 1)));
        for (int i = 1; i <= iz; ++i) set_f0(a7(i, jseq_, 0), f0(a7(i, jseq_, 1)));
        for (int i = 1; i <= iz; ++i) set_msg0(a7(i, jseq_, 0), msg0(a7(i, jseq_, 1)));
        ndec_[jseq_][0] = iz;
        ndec_[jseq_][1] = 0;
        nutc0_ = nutc;
        for (int i = 1; i <= 200; ++i) set_dt0(a7(i, jseq_, 1), 0.0f);
        for (int i = 1; i <= 200; ++i) set_f0(a7(i, jseq_, 1), 0.0f);
    }
    if (trace_)
        trace_("C " + std::to_string(nzhsym) + " early " + std::to_string(ndec_early_) + " jseq " +
               std::to_string(jseq_));

    if (ndepth == 1 && nzhsym < 50) {
        ndec_early_ = 0;
        done();
        return;
    }
    if (ndepth == 1 && nzhsym == 50) load(dd_);
    std::array<int, 58> apsym2;
    std::array<int, 10> aph10;
    packjt77::ft8apset(pk_, mycall12, hiscall12, ncontest, apsym2, aph10);

    if (nzhsym <= 47) {
        load(dd_);
        dd1_ = dd_;
    }

    int ndecodes;
    if (nzhsym == 41) {
        ndecodes = 0;
        std::fill(allmessages_.begin(), allmessages_.end(), std::string(37, ' '));
    } else {
        ndecodes = ndec_early_;
    }

    std::string msg37(37, ' ');
    float xdt = 0, f1 = 0, xsnr = 0, dmin = 0;
    int nharderrors = 0, iaptype = 0, nsnr = 0;
    bool newdat = true;

    bool to900 = false;
    if (nzhsym == 47 && ndec_early_ == 0) {
        dd1_ = dd_;
        goto l800;
    }

    if (nzhsym == 47 && ndec_early_ >= 1) {
        std::fill(lsubtracted_, lsubtracted_ + MAX_EARLY, false);
        bool lrefinedt = true;
        if (ndepth <= 2) lrefinedt = false;
        for (int i = 1; i <= ndec_early_; ++i) {
            if (xdt_save_[i - 1] - 0.5f < 0.396f) {
                fe_.subtractft8(dd_.data(), &itone_save_[static_cast<size_t>(i - 1) * NN], f1_save_[i - 1],
                                xdt_save_[i - 1], lrefinedt);
                lsubtracted_[i - 1] = true;
            }
        }
        dd1_ = dd_;
        to900 = true;
        goto l900;
    }

    if (nzhsym == 50 && ndec_early_ >= 1 && !nagain) {
        const int n = 47 * 3456;
        std::copy(dd1_.begin(), dd1_.begin() + n, dd_.begin());
        for (int i = n; i < NPTS; ++i) dd_[i] = static_cast<float>(iwave16[i]);
        for (int i = 1; i <= ndec_early_; ++i) {
            if (lsubtracted_[i - 1]) continue;
            fe_.subtractft8(dd_.data(), &itone_save_[static_cast<size_t>(i - 1) * NN], f1_save_[i - 1],
                            xdt_save_[i - 1], true);
        }
    }

    {
        int ifa = p.nfa, ifb = p.nfb;
        if (nzhsym == 50 && nagain) {
            load(dd_);
            ifa = nfqso - 20;
            ifb = nfqso + 20;
        }

        // ndepth=1: 1 pass, bp; ndepth=2: subtraction, 3 passes, bp+osd (no
        // subtract refinement); ndepth=3: subtraction, 3 passes, bp+osd
        int npass = 3;
        int imetric = 1;
        bool lsubtract = true;
        if (ndepth == 1) npass = 2;
        for (int ipass = 1; ipass <= npass; ++ipass) {
            if (stopped()) break;
            newdat = true;
            float syncmin = 1.3f;
            if (ndepth <= 2) syncmin = 2.1f;
            if (ipass == 1) {
                lsubtract = true;
                imetric = 1;
            } else if (ipass == 2) {
                imetric = 2;
                lsubtract = true;
            } else if (ipass == 3) {
                imetric = 2;
                if (ndecodes == 0) continue;
                lsubtract = true;
            }
            int ncand = 0;
            fe_.sync8(dd_.data(), NPTS, ifa, ifb, syncmin, nfqso, MAXCAND, candidate_.data(), ncand, sbase_.data());
            if (trace_)
                trace_("P " + std::to_string(ipass) + " " + std::to_string(imetric) + " " + std::to_string(ncand) +
                       (ncand >= 1 ? " " + hex8(candidate_[(ncand - 1) * 3 + 2]) : std::string(" -")));
            for (int icand = 1; icand <= ncand; ++icand) {
                if (stopped()) break;
                sync_ = candidate_[(icand - 1) * 3 + 2];
                f1 = candidate_[(icand - 1) * 3];
                xdt = candidate_[(icand - 1) * 3 + 1];
                const float xbase = xbase_of(sbase_.data(), f1, false);
                msg37.assign(37, ' ');
                int nbadcrc = 0, iappass = 0;
                int itone[NN];
                fe_.ft8b(dd_.data(), newdat, p.nqso_progress, nfqso, nftx, ndepth, nzhsym, lft8apon, p.lapcqonly,
                         p.napwid, lsubtract, nagain, ncontest, imetric, iaptype, mycall12, hiscall12, f1, xdt, xbase,
                         apsym2.data(), aph10.data(), nharderrors, dmin, nbadcrc, iappass, msg37, xsnr, itone,
                         *codec_);
                nsnr = nint(xsnr);
                xdt = xdt - 0.5f;
                if (nbadcrc == 0) {
                    bool ldupe = false;
                    for (int id = 1; id <= ndecodes; ++id)
                        if (msg37 == allmessages_[id - 1]) ldupe = true;
                    if (!ldupe) {
                        if (ndecodes >= MAX_EARLY) continue;
                        ++ndecodes;
                        allmessages_[ndecodes - 1] = msg37;
                        f1_save_[ndecodes - 1] = f1;
                        xdt_save_[ndecodes - 1] = xdt + 0.5f;
                        std::copy(itone, itone + NN, &itone_save_[static_cast<size_t>(ndecodes - 1) * NN]);
                    }
                    if (!ldupe) {
                        const float qual = 1.0f - (static_cast<float>(nharderrors) + dmin) / 60.0f;
                        callback(sync_, nsnr, xdt, f1, msg37, iaptype, qual, nutc, nfqso, c);
                        ft8_a7_save(jseq_, xdt, f1, msg37);  // Enter decode in table
                    }
                }
            }
        }
    }

l800:
    ndec_early_ = 0;
    if (nzhsym < 50) ndec_early_ = ndecodes;

l900:
    (void)to900;
    if (lft8apon && ncontest != 6 && ncontest != 7 && nzhsym == 50 && ndec_[jseq_][0] >= 1) {
        newdat = true;
        const int nz = ndec_[jseq_][0];
        for (int i = 1; i <= nz; ++i) {
            if (stopped()) break;
            const int e = a7(i, jseq_, 0);
            if (f0(e) == -99.0f) break;
            if (f0(e) == -98.0f) continue;
            if (findex(msg0(e), "<") >= 1) continue;  // ### Temporary ###
            msg37 = msg0(e);
            const int i1 = findex(msg37, " ");
            const int i2 = findex(msg37.substr(static_cast<size_t>(i1)), " ") + i1;
            const std::string call_1 = fix(substr(msg37, 1, i1 - 1), 12);
            const std::string call_2 = fix(substr(msg37, i1 + 1, i2 - 1), 12);
            std::string grid4 = fix(substr(msg37, i2 + 1, i2 + 4), 4);
            if (grid4 == "RR73" || findex(grid4, "+") > 0 || findex(grid4, "-") > 0) grid4 = "    ";
            xdt = dt0(e);
            f1 = f0(e);
            const float xbase = xbase_of(sbase_.data(), f1, true);
            msg37.assign(37, ' ');
            if (a7d_first_) {  // ft8_a7d's one(0:511,0:8) on its first call, in its storage
                for (int b = 0; b < 512; ++b)
                    for (int j = 0; j < 9; ++j) {
                        const int32_t one = (b & (1 << j)) != 0 ? 1 : 0;
                        std::memcpy(&a7mem_[A7_ONE + 4 * static_cast<size_t>(b + 512 * j)], &one, 4);
                    }
                a7d_first_ = false;
            }
            fe_.ft8_a7d(dd_.data(), newdat, call_1, call_2, grid4, xdt, f1, xbase, nharderrors, dmin, msg37, xsnr,
                        *codec_);
            if (trace_)
                trace_("A " + std::to_string(i) + " " + std::to_string(nharderrors) + " " + hex8(dmin) + " " +
                       hex8(xsnr) + " " + hex8(sync_) + " |" + msg37 + "|");
            if (nharderrors >= 0) {
                nsnr = static_cast<int>(xsnr);  // nsnr=xsnr
                iaptype = 7;
                const float qual = 1.0f;
                if (findex(msg37, trim(hiscall12)) > 0) la8 = false;
                callback(sync_, nsnr, xdt, f1, msg37, iaptype, qual, nutc, nfqso, c);
                ft8_a7_save(jseq_, xdt, f1, msg37);  // Enter decode in table
            }
        }
    }

    if (lft8apon && ncontest != 6 && ncontest != 7 && nzhsym == 50 && la8 && trim(hiscall12).size() >= 3 &&
        trim(hisgrid).size() >= 4 && ltry_a8_ && !stopped()) {
        // Try for an a8 decode at nfqso
        f1 = static_cast<float>(nfqso);
        float fbest = 0, plog = 0;
        ft8_a8d(dd_.data(), mycall12, hiscall12, hisgrid, f1, xdt, fbest, xsnr, plog, msg37);
        if (msg37[0] != ' ') {
            sync_ = 10.0f;  // ### ???
            nsnr = nint(xsnr);
            iaptype = 8;
            float qual = 1.0f;
            if (plog < -147.0f) qual = 0.16f;
            callback(sync_, nsnr, xdt, fbest, msg37, iaptype, qual, nutc, nfqso, c);
            ft8_a7_save(jseq_, xdt, f1, msg37);  // Enter decode in the a7 table
        }
    }
    done();
}

// ---------------------------------------------------------------- ft8_a7.f90

void Decoder::ft8_a7_save(int jseq, float dt, float f, std::string& msg) {
    if (findex(msg, "/") >= 1 || findex(msg, "<") >= 1) return;
    packjt77::Split77 s = packjt77::split77(msg);  // Parse msg into words
    msg = s.msg;                                    // (split77 rewrites its argument)
    if (s.nwords < 1) return;
    if (s.w[0].compare(0, 3, "CQ_") == 0) return;
    const int j = jseq;

    // Add this decode to current table for this sequence
    ndec_[j][1] = ndec_[j][1] + 1;  // Number of decodes in this sequence
    const int i = ndec_[j][1];      // i is index of a new table entry
    if (i > 200) return;            // Prevent table overflow (indexes start at 1)

    const int e = a7(i, j, 1);
    set_dt0(e, dt);  // Save dt in table
    set_f0(e, f);    // Save f in table
    std::string m = fix(trim(s.w[0]) + " " + trim(s.w[1]), 37);  // Save "call_1 call_2"
    if (s.w[0].compare(0, 3, "CQ ") == 0 && s.nw[1] <= 2)
        m = fix("CQ " + trim(s.w[1]) + " " + trim(s.w[2]), 37);  // Save "CQ DX Call_2"
    // Include grid as part of message
    if (isgrid4(s.w[s.nwords - 1])) m = fix(trim(m) + " " + trim(s.w[s.nwords - 1]), 37);

    // If a transmission at this frequency with message fragment "call_1
    // call_2" was decoded in the previous sequence, flag it as "DO NOT USE"
    // because we have already decoded and subtracted that station's next
    // transmission.
    const packjt77::Split77 s2 = packjt77::split77(m);  // Parse msg into words
    set_msg0(e, s2.msg);
    const std::string key = " " + trim(s2.w[1]);
    for (int ii = 1; ii <= ndec_[j][0]; ++ii) {
        const int q = a7(ii, j, 0);
        if (f0(q) <= -98.0f) continue;
        const int i2 = findex(msg0(q), key);
        if (std::fabs(f - f0(q)) <= 3.0f && i2 >= 3) set_f0(q, -98.0f);  // Flag as "do not use"
    }
}

// ---------------------------------------------------------------- ft8_a8d.f90

void Decoder::ft8_a8d(const float* dd, const std::string& mycall, const std::string& dxcall,
                      const std::string& dxgrid, float f1a, float& xdt, float& fbest, float& xsnr, float& plog,
                      std::string& msgbest) {
    constexpr int NWAVE = NN * NSPS / NDOWN;  // Length of generated cwave = 2528
    constexpr int NZZ = 3200;                 // Length of downsampled arrays
    constexpr int NFFT = NZZ, NH = NFFT / 2;  // FFT length
    constexpr int NMSGS = 206;                // Number of messages to be tried
    std::vector<float> s(2 * NH + 1, 0.0f), s0(2 * NH + 1, 0.0f), s1(2 * NH + 1, 0.0f);  // (-NH:NH) at [j+NH]
    std::vector<cfloat> cwave(NZZ), cd(NZZ);
    std::vector<float> xjunk(NZZ);
    cfloat* cd0 = a8_cd0_.data();
    cfloat* csymb = a8_csymb_.data();
    int itone[NN], itone_best[NN];

    const float f1 = f1a;
    bool newdat = true;
    // Mix from f1 to baseband and downsample from dd into cd (complex, 200
    // Hz sampling)
    fe_.ft8_downsample(dd, newdat, f1, cd.data());
    const float fac = 1.e-6f;

    const float fsd = 200.0f;  // Downsampled rate (Hz)
    const float bt = 2.0f;     // Gaussian smoothing parameter
    const float dt = 1.0f / fsd;
    const float df = 200.0f / static_cast<float>(NFFT);

    // Find the message that fits the received waveform best
    float sbest = 0.0f, tbest = 0.0f;
    fbest = 0.0f;
    msgbest.assign(37, ' ');
    std::fill(itone_best, itone_best + NN, 0);
    for (int imsg = 1; imsg <= NMSGS; ++imsg) {
        const std::string msg = getmsg(imsg, mycall, dxcall, dxgrid);

        // Source-encode the message, get itone(), and generate complex FT8
        // waveform.
        pk_.pack77(msg);
        const packjt77::Packjt77::Genft8 g = pk_.genft8(msg);
        if (g.ok) {
            int8_t msgbits[77];
            for (int q = 0; q < 77; ++q) msgbits[q] = static_cast<int8_t>(g.msgbits[q]);
            get_ft8_tones_from_77bits(msgbits, itone);
        } else {
            std::fill(itone, itone + NN, 0);
        }
        fe_.gen_ft8wave(itone, NN, 32, bt, fsd, 0.0f, cwave.data(), xjunk.data(), 1, NWAVE);
        std::fill(cwave.begin() + NWAVE, cwave.end(), cfloat(0.0f, 0.0f));

        // For this message, find the best lag
        float spk = 0.0f, fpk = 0.0f, tpk = 0.0f;
        int lagpk = 0, lag1 = -200, lag2 = 200, lagstep = 4;
        for (int iter = 1; iter <= 2; ++iter) {
            if (iter == 2) {
                lag1 = lagpk - 8;
                lag2 = lagpk + 8;
                lagstep = 1;
            }
            for (int lag = lag1; lag <= lag2; lag += lagstep) {
                for (int i = 0; i <= NWAVE - 1; ++i) {
                    const int j = i + lag + 100;
                    cd0[i] = cfloat(0.0f, 0.0f);
                    if (j >= 0 && j <= NWAVE - 1) cd0[i] = cmul(cd[j], std::conj(cwave[i]));
                }
                std::fill(cd0 + NWAVE, cd0 + NFFT, cfloat(0.0f, 0.0f));
                fst4w::fft_c2c_inplace(cd0, NFFT, -1);  // Forward c2c FFT

                for (int i = 0; i <= NFFT - 1; ++i) {
                    int j = i;
                    if (i > NH) j = i - NFFT;
                    s[j + NH] = fac * (cd0[i].real() * cd0[i].real() + cd0[i].imag() * cd0[i].imag());
                }

                smo121(s.data(), NFFT + 1);

                for (int j = -NH; j <= NH; ++j) {
                    if (s[j + NH] > spk) {
                        spk = s[j + NH];
                        fpk = static_cast<float>(j) * df + f1;
                        lagpk = lag;
                        tpk = static_cast<float>(lag) * dt;
                        s0 = s;  // Save best-fit spectrum for this msg
                    }
                }
            }
        }

        if (spk > sbest) {
            sbest = spk;
            fbest = fpk;
            tbest = tpk;
            msgbest = msg;
            std::copy(itone, itone + NN, itone_best);
            s1 = s0;  // Save overall best-fit spectrum
        }
    }

    float a[5] = {0, 0, 0, 0, 0};
    a[0] = f1 - fbest;
    if (std::fabs(a[0]) > 5.0f) {  // Effective Ftol = 5 Hz
        msgbest.assign(37, ' ');
        return;
    }
    twkfreq1(cd.data(), NZZ, 200.0f, a, cd.data());  // Re-center the spectrum

    xdt = tbest;
    float sa = 0.0f, sb = 0.0f;
    for (int j = -200; j <= -100; ++j) sa = sa + s1[j + NH];
    for (int j = 100; j <= 200; ++j) sb = sb + s1[j + NH];
    const float ave = (sa + sb) / 202.0f;
    for (float& v : s1) v = v / ave - 1.0f;
    float s1pk = s1[-32 + NH];
    for (int j = -31; j <= 32; ++j)
        if (s1[j + NH] > s1pk || std::isnan(s1pk)) s1pk = s1[j + NH];
    float sig = 0.0f;
    int nsig = 0;
    for (int i = -32; i <= 32; ++i) {
        if (s1[i + NH] < 0.5f * s1pk) continue;
        sig = sig + s1[i + NH];
        ++nsig;
    }
    sig = sig / static_cast<float>(nsig);
    xsnr = db(sig) - 35.0f;
    if (xsnr < -30.0f) xsnr = -30.0f;  // Set min SNR for 'a8' decodes

    if (!blank(msgbest)) {
        // Compute probability of successful decode.
        plog = 0.0f;
        int nhard = 0, nsum = 0;
        float sum_sync = 0.0f, sum_sig = 0.0f, sum_big = 0.0f;
        for (int k = 1; k <= NN; ++k) {
            const int i0 = 32 * (k - 1) + nint((tbest + 0.5f) / 0.005f);
            std::fill(csymb, csymb + 32, cfloat(0.0f, 0.0f));
            for (int i = 0; i <= 31; ++i)
                if (i0 + i >= 0 && i0 + i <= NZZ - 1) csymb[i] = cd[i0 + i];
            fst4w::fft_c2c_inplace(csymb, 32, -1);  // c2c
            float s8[8];
            for (int i = 0; i <= 7; ++i) s8[i] = csymb[i].real() * csymb[i].real() + csymb[i].imag() * csymb[i].imag();
            float s8sum = 0.0f;
            for (int i = 0; i <= 7; ++i) s8sum = s8sum + s8[i];
            if (s8sum > 0.0f) {
                const float pr = s8[itone_best[k - 1]] / s8sum;
                plog = plog + logf_glibc(pr);
                ++nsum;
            }
            int ipk = 0;
            for (int i = 1; i <= 7; ++i)
                if (s8[i] > s8[ipk] || (std::isnan(s8[ipk]) && !std::isnan(s8[i]))) ipk = i;
            if (ipk != itone_best[k - 1]) ++nhard;
            if (k <= 7 || (k >= 37 && k <= 43) || k >= 73)
                sum_sync = sum_sync + s8[itone_best[k - 1]];
            else
                sum_sig = sum_sig + s8[itone_best[k - 1]];
            sum_big = sum_big + s8[ipk];
        }
        if (nsum < NN) plog = plog + static_cast<float>(NN - nsum) * logf_glibc(0.125f);
        const float sigobig = (sum_sync + sum_sig) / sum_big;
        if (nhard > 54 || plog < -159.0f || sigobig < 0.71f) msgbest.assign(37, ' ');
    }
}

// ---------------------------------------------------------------- the harness

bool read_wav(const std::string& path, std::vector<int16_t>& samples, int& rate, std::string& why) {
    std::ifstream f(path, std::ios::binary);
    if (!f) {
        why = "cannot open";
        return false;
    }
    char hdr[12];
    if (!f.read(hdr, 12) || std::memcmp(hdr, "RIFF", 4) != 0 || std::memcmp(hdr + 8, "WAVE", 4) != 0) {
        why = "not a WAV file";
        return false;
    }
    rate = 0;
    bool fmt = false;
    for (;;) {
        char id[4];
        uint32_t len;
        if (!f.read(id, 4) || !f.read(reinterpret_cast<char*>(&len), 4)) break;
        if (std::memcmp(id, "fmt ", 4) == 0) {
            std::vector<char> b(len);
            if (!f.read(b.data(), static_cast<std::streamsize>(len)) || len < 16) break;
            uint16_t format, channels, bits;
            uint32_t sr;
            std::memcpy(&format, &b[0], 2);
            std::memcpy(&channels, &b[2], 2);
            std::memcpy(&sr, &b[4], 4);
            std::memcpy(&bits, &b[14], 2);
            if (format != 1 || channels != 1 || bits != 16) {
                why = "not 16-bit mono PCM";
                return false;
            }
            rate = static_cast<int>(sr);
            fmt = true;
            if (len & 1) f.seekg(1, std::ios::cur);
        } else if (std::memcmp(id, "data", 4) == 0) {
            if (!fmt) break;
            samples.assign(len / 2, 0);
            f.read(reinterpret_cast<char*>(samples.data()), static_cast<std::streamsize>(samples.size() * 2));
            samples.resize(static_cast<size_t>(f.gcount()) / 2);
            return true;
        } else {
            f.seekg(len + (len & 1), std::ios::cur);
        }
    }
    why = "no data";
    return false;
}

int nutc_of(const std::string& filename) {
    std::string name = trim(filename);
    if (name.empty()) return 0;
    const size_t dot = name.rfind('.');
    if (dot == std::string::npos) return 0;
    const std::string ext = name.substr(dot);
    if (ext != ".wav" && ext != ".WAV") return 0;
    const std::string stem = name.substr(0, dot);
    auto digits = [](const std::string& t) {
        if (t.empty()) return false;
        for (char ch : t)
            if (ch < '0' || ch > '9') return false;
        return true;
    };
    if (stem.size() >= 6 && digits(stem.substr(stem.size() - 6))) return std::stoi(stem.substr(stem.size() - 6));
    if (stem.size() >= 5 && stem[stem.size() - 5] == '_' && digits(stem.substr(stem.size() - 4)))
        return std::stoi(stem.substr(stem.size() - 4));
    return 0;
}

int run_periods(const std::vector<std::string>& args, const std::string& cwd, RunOutput& r) {
    auto at = [&](const std::string& p) { return !p.empty() && p[0] == '/' ? p : cwd + "/" + p; };
    auto usage = [&]() {
        r.err += "usage: ft8-oracle -d 1|2|3 WAV...\n"
                 "       ft8-oracle --file [-d D -L HZ -H HZ -f HZ -F HZ -c CALL -x CALL -g GRID -G GRID -Q N -X N "
                 "-p S] WAV...\n";
        return 2;
    };
    // The harness's options (and jt9's, for "jt9 -a DIR -8 ...").
    bool filemode = !args.empty() && (args[0] == "jt9" || (args[0].size() >= 4 &&
                                                               args[0].compare(args[0].size() - 4, 4, "/jt9") == 0));
    double trperiod = 60.0;  // jt9's TRperiod (-p)
    int ndepth = -1, flow = 200, fhigh = 4000, nrxfreq = 1500, nqso = 0, nexp_decode = 0;
    std::string mycall = "K1ABC", hiscall = "W9XYZ", hisgrid = "EN37";
    size_t iarg = 1;
    // list-directed read(arg,*) of a character variable: the first value
    // (up to a blank, comma or slash).
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
        } else if (a == "-8") {
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
            if (!num(ntol)) return usage();  // (FT8 does not read ntol)
        } else if (a == "-Q") {
            if (!num(nqso)) return usage();
        } else if (a == "-X") {
            if (!num(nexp_decode)) return usage();
        } else if (a == "-p") {
            // require_real and jt9's range check (TRperiod, real*8)
            if (++iarg >= args.size()) return usage();
            const std::string& t = args[iarg];
            if (!real_token(t) || !std::isfinite(std::strtod(t.c_str(), nullptr))) {
                r.err += "jt9: invalid value for -p: " + trim(t) + " (must be a finite number)\n";
                return 2;
            }
            trperiod = std::strtod(t.c_str(), nullptr);
            if (trperiod <= 0.0 || trperiod > 1800.0) {
                r.err += "jt9: invalid value for -p: " + trim(t) + " (must be between 0 and 1800 seconds)\n";
                return 2;
            }
        } else if (a == "-c" || a == "-x" || a == "-G" || a == "-g") {
            if (++iarg >= args.size()) return usage();
            const std::string v = word(args[iarg]);
            if (a == "-c") mycall = v;
            if (a == "-x") hiscall = v;
            if (a == "-g") hisgrid = v;
        } else {
            return usage();
        }
    }
    if (!filemode && (ndepth < 1 || ndepth > 3)) return usage();
    if (filemode && ndepth < 0) ndepth = 1;
    if (iarg >= args.size()) return usage();

    Params p;
    p.ndepth = ndepth;
    if (filemode) {
        // jt9_params_init's: nftx not set (zero), napwid 75.
        if (mycall == "b") mycall.clear();
        if (hiscall == "b") {
            hiscall.clear();
            hisgrid.clear();
        }
        p.nfa = flow;
        p.nfb = fhigh;
        p.nfqso = nrxfreq;
        p.nftx = 0;
        p.napwid = 75;
        p.nqso_progress = nqso;
        p.nexp_decode = nexp_decode;
        p.mycall = mycall;
        p.hiscall = hiscall;
        p.hisgrid = hisgrid;
    }
    Decoder dec;
    if (r.want_trace) dec.set_trace([&r](const std::string& l) { r.trace += l + "\n"; });
    for (; iarg < args.size(); ++iarg) {
        const std::string& infile = args[iarg];
        std::vector<int16_t> s;
        int rate = 0;
        std::string why;
        if (!read_wav(at(infile), s, rate, why)) {
            r.err += (filemode ? "jt9" : "ft8-oracle") + std::string(": cannot read WAV file ") + infile + ": " +
                     why + "\n";
            return 2;
        }
        if (rate != 12000) {
            if (filemode)
                r.err += "jt9: unsupported " + std::to_string(rate) + " Hz WAV sample rate for the selected mode\n";
            else
                r.err += "ft8-oracle: " + infile + " is not 12000 Hz\n";
            return 2;
        }
        const int nutc = nutc_of(infile);
        if (filemode) {
            // jt9 -8 reads whole blocks of 3456 samples, npts/3456 of them
            // (npts=TRperiod*12000), up to 179,712, and says so when the file
            // ends first.
            const int npts = static_cast<int>(trperiod * 12000.0);
            const size_t want = std::min<size_t>(static_cast<size_t>(npts / 3456) * 3456, 179712);
            if (s.size() < want) r.out += " EOF on input file " + trim(infile) + "\n";
            s.resize(std::min(s.size(), want));
            const Decoder::Period per = dec.period(s.data(), s.size(), nutc, p, 50 * 3456);
            r.out += per.out;
        } else {
            auto snapshot = [&](int nz) {
                r.decoded_calls += "== " + fst4w::fortran_i(nutc, 6, 6) + " " + std::to_string(nz) + "\n";
                if (dec.decoded_exists()) {
                    // read(79,'(a)') line; write(80,'(a)') trim(line)
                    size_t a = 0;
                    const std::string& t = dec.decoded_text();
                    while (a < t.size()) {
                        size_t e = t.find('\n', a);
                        if (e == std::string::npos) e = t.size();
                        r.decoded_calls += trim(t.substr(a, e - a)) + "\n";
                        a = e + 1;
                    }
                }
            };
            const Decoder::Period per = dec.period(s.data(), s.size(), nutc, p, NPTS, snapshot);
            r.out += per.out;
        }
    }
    r.decoded_exists = dec.decoded_exists();
    r.decoded = dec.decoded_text();
    return 0;
}

}  // namespace ft8
