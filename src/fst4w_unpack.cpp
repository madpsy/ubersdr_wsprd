// FST4W's message unpacking and callsign hash table (fst4w_unpack.h):
// WSJT-X's lib/77bit/packjt77.f90, packjt77_schema.f90 and
// packjt77_grammar.f90 at 567ad29ce, the parts an FST4W decode reaches:
// unpack77 for i3=0 n3=6 (WSPR types 1, 2 and 3) with unpack28, hash22,
// to_grid4, to_grid, pack77_wspr_affix_text, callok, save_hash_call,
// ihashcall and the effects unpack77 applies, and fst4_decode's AP setup.
//
// Copyright 2018-2025 Joe Taylor, K1JT, Steven Franke, K9AN, and the WSJT-X
// developers; GNU GPL v3, and so is this.
//
// Strings are Fortran's fixed-length ones (Call13 and friends, blank
// padded); assignments truncate and pad as Fortran's do, comparisons ignore
// trailing blanks as Fortran's do.
#include "fst4w_unpack.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <stdexcept>

namespace fst4w {

namespace {

constexpr int PACK77_NTOKENS = 2063592;
constexpr int PACK77_MAX22 = 4194304;
constexpr int PACK77_WSPR_NZZZ = 46656;
constexpr const char* PACK77_BASE36 = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
constexpr const char* PACK77_BASE38 = " 0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ/";

// fst4_decode's fixed tail: c77(51:77).
constexpr const char* FST4W_TAIL = "000000000000000000000110000";

// Fortran character(len=N) helpers.
template <size_t N>
void assign(std::array<char, N>& dst, const std::string& src) {
    for (size_t i = 0; i < N; ++i) dst[i] = i < src.size() ? src[i] : ' ';
}

template <size_t N>
std::array<char, N> blank() {
    std::array<char, N> a;
    a.fill(' ');
    return a;
}

template <size_t N>
int len_trim(const std::array<char, N>& a) {
    int n = static_cast<int>(N);
    while (n > 0 && a[n - 1] == ' ') --n;
    return n;
}

template <size_t N>
std::string trim(const std::array<char, N>& a) {
    return std::string(a.data(), static_cast<size_t>(len_trim(a)));
}

template <size_t N>
std::string full(const std::array<char, N>& a) {
    return std::string(a.data(), N);
}

// Fortran's a .eq. b for strings: the shorter blank padded.
bool feq(const std::string& a, const std::string& b) {
    const size_t n = std::max(a.size(), b.size());
    for (size_t i = 0; i < n; ++i) {
        const char x = i < a.size() ? a[i] : ' ';
        const char y = i < b.size() ? b[i] : ' ';
        if (x != y) return false;
    }
    return true;
}

// adjustl: leading blanks moved to the end.
template <size_t N>
std::array<char, N> adjustl(const std::array<char, N>& a) {
    size_t k = 0;
    while (k < N && a[k] == ' ') ++k;
    std::array<char, N> r = blank<N>();
    for (size_t i = k; i < N; ++i) r[i - k] = a[i];
    return r;
}

// index(s, c): 1-based, 0 when absent.
template <size_t N>
int index_of(const std::array<char, N>& a, char c) {
    for (size_t i = 0; i < N; ++i)
        if (a[i] == c) return static_cast<int>(i) + 1;
    return 0;
}

int bits_value(const char* c77, int start_bit, int width) {
    int v = 0;
    for (int i = 0; i < width; ++i) v = 2 * v + (c77[start_bit - 1 + i] == '1');
    return v;
}

bool isdig(char c) { return c >= '0' && c <= '9'; }
bool islet(char c) { return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z'); }

// packjt77.f90's callok(trim(c13)).
bool callok(const std::string& w) {
    const int n = static_cast<int>(w.size());  // len(trim(w)) of an already trimmed w
    if (n < 3) return false;
    if (w[0] == 'Q') return false;
    int i;
    for (i = n; i >= 1; --i)
        if (isdig(w[i - 1])) break;
    const int i0 = i;
    if (i0 != 2 && i0 != 3) return false;
    // pfx is character*2: w(1:i0-1), blank padded; len_trim of it.
    std::string pfx = w.substr(0, static_cast<size_t>(i0 - 1));
    while (!pfx.empty() && pfx.back() == ' ') pfx.pop_back();
    int nlp = 0, ndp = 0;
    const int np = static_cast<int>(pfx.size());
    for (int k = 0; k < np; ++k) {
        if (isdig(pfx[k])) ++ndp;
        if (islet(pfx[k])) ++nlp;
    }
    if (nlp + ndp != np) return false;
    if (nlp == 0) return false;
    int nls = 0;
    const int ns = n - i0;
    if (ns < 1 || ns > 3) return false;
    for (int k = i0 + 1; k <= n; ++k)
        if (islet(w[k - 1])) ++nls;
    if (nls < ns) return false;
    return true;
}

// packjt77_grammar.f90's pack77_c28_standard_shape.
bool c28_standard_shape(const std::string& token) {
    std::string t = token;
    while (!t.empty() && t.back() == ' ') t.pop_back();
    const int n = static_cast<int>(t.size());
    if (n < 3 || n > 6) return false;
    for (char c : t)
        if (!((c >= 'A' && c <= 'Z') || isdig(c))) return false;
    int iarea = 0;
    for (int i = n; i >= 2; --i)
        if (isdig(t[i - 1])) {
            iarea = i;
            break;
        }
    if (iarea != 2 && iarea != 3) return false;
    int npdig = 0, nplet = 0;
    for (int i = 1; i <= iarea - 1; ++i) {
        const char c = t[i - 1];
        if (isdig(c)) ++npdig;
        else if (c >= 'A' && c <= 'Z') ++nplet;
        else return false;
    }
    if (nplet == 0 || npdig >= iarea - 1) return false;
    int nslet = 0;
    for (int i = iarea + 1; i <= n; ++i) {
        const char c = t[i - 1];
        if (c >= 'A' && c <= 'Z') ++nslet;
        else return false;
    }
    return nslet <= 3;
}

// packjt77_grammar.f90's pack77_valid_hash_call_token.
bool valid_hash_call_token(const Call13& token) {
    const int n = len_trim(token);
    if (n < 5 || n > 13) return false;
    if (token[0] != '<' || token[n - 1] != '>') return false;
    for (int i = 1; i < n - 1; ++i)
        if (token[i] == '<' || token[i] == '>') return false;
    return true;
}

// packjt77.f90's normalize_hash_call: the call a "<CALL>" token or a plain
// call stages, or none.
bool normalize_hash_call(const Call13& c13, Call13& cw) {
    cw = c13;
    if (cw[0] == ' ' || std::memcmp(cw.data(), "<...>", 5) == 0) return false;
    if (cw[0] == '<') {
        if (!valid_hash_call_token(cw)) return false;
        Call13 s = blank<13>();
        for (int i = 1; i < 13; ++i) s[i - 1] = cw[i];
        cw = s;
        const int i = index_of(cw, '>');
        if (i > 0)
            for (int k = i - 1; k < 13; ++k) cw[k] = ' ';
    }
    return len_trim(cw) >= 3;
}

// packjt77_grammar.f90's to_grid4.
bool to_grid4(int n, std::array<char, 4>& grid4) {
    const int j1 = n / (18 * 10 * 10);
    if (j1 < 0 || j1 > 17) return false;
    n = n - j1 * 18 * 10 * 10;
    const int j2 = n / (10 * 10);
    if (j2 < 0 || j2 > 17) return false;
    n = n - j2 * 10 * 10;
    const int j3 = n / 10;
    if (j3 < 0 || j3 > 9) return false;
    const int j4 = n - j3 * 10;
    if (j4 < 0 || j4 > 9) return false;
    grid4 = {static_cast<char>(j1 + 'A'), static_cast<char>(j2 + 'A'), static_cast<char>(j3 + '0'),
             static_cast<char>(j4 + '0')};
    return true;
}

// packjt77_grammar.f90's to_grid (4- or 6-character).
bool to_grid(int n, std::array<char, 6>& grid6) {
    const int j1 = n / (18 * 10 * 10 * 25 * 25);
    if (j1 < 0 || j1 > 17) return false;
    n = n - j1 * 18 * 10 * 10 * 25 * 25;
    const int j2 = n / (10 * 10 * 25 * 25);
    if (j2 < 0 || j2 > 17) return false;
    n = n - j2 * 10 * 10 * 25 * 25;
    const int j3 = n / (10 * 25 * 25);
    if (j3 < 0 || j3 > 9) return false;
    n = n - j3 * 10 * 25 * 25;
    const int j4 = n / (25 * 25);
    if (j4 < 0 || j4 > 9) return false;
    n = n - j4 * 25 * 25;
    const int j5 = n / 25;
    if (j5 < 0 || j5 > 24) return false;
    const int j6 = n - j5 * 25;
    if (j6 < 0 || j6 > 24) return false;
    grid6 = blank<6>();
    grid6[0] = static_cast<char>(j1 + 'A');
    grid6[1] = static_cast<char>(j2 + 'A');
    grid6[2] = static_cast<char>(j3 + '0');
    grid6[3] = static_cast<char>(j4 + '0');
    if (j5 != 24 || j6 != 24) {
        grid6[4] = static_cast<char>(j5 + 'A');
        grid6[5] = static_cast<char>(j6 + 'A');
    }
    return true;
}

// packjt77_grammar.f90's pack77_wspr_affix_text.
bool wspr_affix_text(int npfx, std::array<char, 3>& cpfx, bool& is_prefix) {
    int n = npfx;
    cpfx = blank<3>();
    is_prefix = n < PACK77_WSPR_NZZZ;
    if (is_prefix) {
        for (int i = 3; i >= 1; --i) {
            const int j = n % 36 + 1;
            cpfx[i - 1] = PACK77_BASE36[j - 1];
            n = n / 36;
            if (n == 0) break;
        }
    } else {
        n = n - PACK77_WSPR_NZZZ;
        if (n <= 35) {
            cpfx[0] = PACK77_BASE36[n];
        } else if (n > 35 && n <= 1295) {
            cpfx[0] = PACK77_BASE36[n / 36];
            cpfx[1] = PACK77_BASE36[n % 36];
        } else if (n > 1295 && n <= 12959) {
            cpfx[0] = PACK77_BASE36[n / 360];
            cpfx[1] = PACK77_BASE36[(n / 10) % 36];
            cpfx[2] = PACK77_BASE36[n % 10];
        } else {
            return false;
        }
    }
    return true;
}

// nint(idbm*10.0/3.0) in default real, then write(crpt,'(i3)'): trimmed.
std::string dbm_text(int idbm) {
    const float x = static_cast<float>(idbm) * 10.0f / 3.0f;
    const long v = std::lround(x);
    char b[16];
    std::snprintf(b, sizeof b, "%3ld", v);
    std::string s(b);
    if (s.size() > 3) s = "***";  // i3 overflow (cannot happen: 0..60)
    size_t k = 0;
    while (k < s.size() && s[k] == ' ') ++k;
    return s.substr(k);
}

}  // namespace

Unpacker::Unpacker() {
    for (auto& c : calls22_) c = blank<13>();
    ihash22_.fill(-1);
    mycall13_ = blank<13>();
    dxcall13_ = blank<13>();
    mycall0_ = blank<12>();
    hiscall0_ = blank<12>();
    mycall13_0_ = blank<13>();
}

std::string Unpacker::calls22(int i) const { return trim(calls22_[i - 1]); }

int Unpacker::ihashcall(const Call13& c0, int m) {
    int64_t n8 = 0;
    for (int i = 0; i < 11; ++i) {
        const char* p = std::strchr(PACK77_BASE38, c0[i]);
        // Fortran's index() - 1: -1 when absent (and for NUL, never in a call).
        const int64_t j = (c0[i] != '\0' && p) ? (p - PACK77_BASE38) : -1;
        n8 = 38 * n8 + j;
    }
    // ihashcall_from_n8: the top m bits of the low 64 bits of 47055833459*n8,
    // in 16-bit limbs with Fortran's (truncating) mod and division.
    static const int64_t hash_factor[4] = {2419, 62655, 10, 0};
    int64_t b[4], p[4];
    int64_t x = n8;
    for (int i = 0; i < 4; ++i) {
        b[i] = x % 65536;
        x = x / 65536;
    }
    int64_t carry = 0;
    for (int i = 0; i < 4; ++i) {
        int64_t term = carry;
        for (int j = 0; j <= i; ++j) term = term + hash_factor[j] * b[i - j];
        p[i] = term % 65536;
        carry = term / 65536;
    }
    if (m <= 16) return static_cast<int>(p[3] / (int64_t(1) << (16 - m)));
    return static_cast<int>(p[3] * (int64_t(1) << (m - 16)) + p[2] / (int64_t(1) << (32 - m)));
}

// packjt77.f90's save_hash_call, its 22-bit table (the 10- and 12-bit ones
// are not kept: FST4W never reads them).
void Unpacker::save_hash_call(const Call13& c13) {
    Call13 cw = c13;
    if (cw[0] == ' ' || std::memcmp(cw.data(), "<...>", 5) == 0) return;
    if (cw[0] == '<') {
        Call13 s = blank<13>();
        for (int i = 1; i < 13; ++i) s[i - 1] = cw[i];
        cw = s;
    }
    const int i = index_of(cw, '>');
    if (i > 0)
        for (int k = i - 1; k < 13; ++k) cw[k] = ' ';
    if (len_trim(cw) < 3) return;
    const int n22 = ihashcall(cw, 22);
    // any(ihash22.eq.n22) over the whole array, unused (-1) slots too.
    bool exists = false;
    for (int k = 0; k < MAXHASH; ++k)
        if (ihash22_[k] == n22) exists = true;
    if (exists) {
        for (int k = 0; k < MAXHASH; ++k)
            if (ihash22_[k] == n22) calls22_[k] = cw;
        return;
    }
    for (int k = MAXHASH - 1; k >= 1; --k) {
        ihash22_[k] = ihash22_[k - 1];
        calls22_[k] = calls22_[k - 1];
    }
    ihash22_[0] = n22;
    calls22_[0] = cw;
    if (nzhash_ < MAXHASH) ++nzhash_;
}

Call13 Unpacker::hash22(int n22) const {
    Call13 c13;
    assign(c13, "<...>");
    for (int i = 0; i < nzhash_; ++i) {
        if (ihash22_[i] == n22) {
            assign(c13, "<" + trim(calls22_[i]) + ">");
            break;
        }
    }
    return c13;
}

// packjt77.f90's unpack28 (unpack28_core, two-digit prefixes allowed).
bool Unpacker::unpack28(int n28_0, Call13& c13) const {
    static const char* c1 = " 0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    static const char* c2 = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    static const char* c3 = "0123456789";
    static const char* c4 = " ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    bool success = true;
    int n28 = n28_0;
    bool done = false;
    if (n28 < PACK77_NTOKENS) {
        if (n28 == 0) assign(c13, "DE");
        if (n28 == 1) assign(c13, "QRZ");
        if (n28 == 2) assign(c13, "CQ");
        if (n28 <= 2) {
            done = true;
        } else if (n28 <= 1002) {
            char b[16];
            std::snprintf(b, sizeof b, "CQ_%03d", n28 - 3);
            assign(c13, b);
            done = true;
        } else if (n28 <= 532443) {
            int n = n28 - 1003;
            const int i1 = n / (27 * 27 * 27);
            n = n - 27 * 27 * 27 * i1;
            const int i2 = n / (27 * 27);
            n = n - 27 * 27 * i2;
            const int i3 = n / 27;
            const int i4 = n - 27 * i3;
            Call13 t;
            assign(t, std::string{c4[i1], c4[i2], c4[i3], c4[i4]});
            t = adjustl(t);
            assign(c13, "CQ_" + full(t).substr(0, 10));
            done = true;
        }
    }
    if (!done) {
        n28 = n28 - PACK77_NTOKENS;
        if (n28 < PACK77_MAX22) {
            c13 = hash22(n28);
        } else {
            int n = n28 - PACK77_MAX22;
            const int i1 = n / (36 * 10 * 27 * 27 * 27);
            n = n - 36 * 10 * 27 * 27 * 27 * i1;
            const int i2 = n / (10 * 27 * 27 * 27);
            n = n - 10 * 27 * 27 * 27 * i2;
            const int i3 = n / (27 * 27 * 27);
            n = n - 27 * 27 * 27 * i3;
            const int i4 = n / (27 * 27);
            n = n - 27 * 27 * i4;
            const int i5 = n / 27;
            const int i6 = n - 27 * i5;
            assign(c13, std::string{c1[i1], c2[i2], c3[i3], c4[i4], c4[i5], c4[i6]});
            c13 = adjustl(c13);
            if (!callok(trim(c13))) {
                assign(c13, "QU1RK");
                success = false;
            }
        }
    }
    const int i0 = index_of(c13, ' ');
    if (i0 != 0 && i0 < len_trim(c13)) {
        assign(c13, "QU1RK");
        success = false;
    }
    return success;
}

void Unpacker::set_calls(const std::string& mycall_in, const std::string& hiscall_in) {
    std::array<char, 12> mycall, hiscall;
    assign(mycall, mycall_in);
    assign(hiscall, hiscall_in);
    // dxcall13=hiscall; mycall13=mycall (character*12 into *13).
    assign(dxcall13_, full(hiscall));
    assign(mycall13_, full(mycall));
    // Cut at a NUL, as fst4_decode does.
    for (auto* c : {&mycall, &hiscall}) {
        const int l1 = index_of(*c, '\0');
        if (l1 != 0)
            for (int k = l1 - 1; k < 12; ++k) (*c)[k] = ' ';
    }
    if (full(mycall) == full(mycall0_) && full(hiscall) == full(hiscall0_)) return;
    if (len_trim(mycall) >= 3) {
        std::array<char, 12> his0 = hiscall;
        if (len_trim(his0) < 3) his0 = mycall;  // mycall as a dummy hiscall
        // pack77("MYCALL HISCALL RR73"): a type 1 message whose two calls
        // pack28 stages, committed in that order; then unpack77 of it, whose
        // entry saves mycall13 (when it changed) and whose type 1 effects
        // save the same two calls again: no change.
        const std::string a = trim(mycall), b = trim(his0);
        if (!c28_standard_shape(a) || !callok(a) || !c28_standard_shape(b) || !callok(b))
            throw std::invalid_argument("fst4w: only standard calls are supported for mycall and hiscall");
        Call13 ca, cb;
        assign(ca, a);
        assign(cb, b);
        save_hash_call(ca);
        save_hash_call(cb);
        if (!feq(full(mycall13_), full(mycall13_0_)) && len_trim(mycall13_) > 2) {
            mycall13_0_ = mycall13_;
            save_hash_call(mycall13_);
        }
        save_hash_call(ca);
        save_hash_call(cb);
    }
    mycall0_ = mycall;
    hiscall0_ = hiscall;
}

Unpacker::Result Unpacker::unpack(const uint8_t* payload) {
    // unpack77's entry: mycall13 saved when it changed.
    if (!feq(full(mycall13_), full(mycall13_0_)) && len_trim(mycall13_) > 2) {
        mycall13_0_ = mycall13_;
        save_hash_call(mycall13_);
    }

    char c77[78];
    for (int i = 0; i < 50; ++i) c77[i] = payload[i] ? '1' : '0';
    std::memcpy(c77 + 50, FST4W_TAIL, 27);
    c77[77] = '\0';

    std::array<char, 37> msg = blank<37>();
    bool success = true;
    Call13 staged = blank<13>();
    bool have_staged = false;
    auto set_msg = [&](const std::string& s) { assign(msg, s); };

    // decode_pack77_tag: i3 = 0, n3 = 6 always with this tail. Then the
    // payload type from bits 49 and 50.
    int wspr_type = 0;
    if (c77[49] == '1') wspr_type = 2;
    else if (c77[48] == '0') wspr_type = 1;
    else if (c77[48] == '1') wspr_type = 3;

    if (wspr_type == 2) {
        // n28 1-28, npfx 29-44, idbm 45-49 (0..18), j50 = 1.
        const int n28 = bits_value(c77, 1, 28);
        const int npfx = bits_value(c77, 29, 16);
        const int idbm_f = bits_value(c77, 45, 5);
        if (idbm_f > 18) return {false, full(msg)};
        Call13 call_1;
        if (!unpack28(n28, call_1)) success = false;
        const std::string crpt = dbm_text(idbm_f);
        std::array<char, 3> cpfx;
        bool is_prefix;
        if (!wspr_affix_text(npfx, cpfx, is_prefix)) return {false, full(msg)};
        Call13 call_1a;
        if (is_prefix) assign(call_1a, trim(adjustl(cpfx)) + "/" + trim(call_1));
        else assign(call_1a, trim(call_1) + "/" + trim(adjustl(cpfx)));
        set_msg(trim(call_1a) + " " + crpt);
        have_staged = normalize_hash_call(call_1a, staged);
    } else if (wspr_type == 1) {
        // n28 1-28, igrid4 29-43, idbm 44-48 (0..18), j49 = j50 = 0.
        const int n28 = bits_value(c77, 1, 28);
        const int igrid4 = bits_value(c77, 29, 15);
        const int idbm_f = bits_value(c77, 44, 5);
        if (idbm_f > 18) return {false, full(msg)};
        Call13 call_1;
        if (!unpack28(n28, call_1)) success = false;
        std::array<char, 4> grid4 = blank<4>();  // the reference's is undefined when to_grid4 fails
        if (!to_grid4(igrid4, grid4)) success = false;
        set_msg(trim(call_1) + " " + full(grid4) + " " + dbm_text(idbm_f));
        have_staged = normalize_hash_call(call_1, staged);
    } else if (wspr_type == 3) {
        // n22 1-22, igrid6 23-47, selector 48-50 = 2.
        const int n22 = bits_value(c77, 1, 22);
        const int igrid6 = bits_value(c77, 23, 25);
        if (bits_value(c77, 48, 3) != 2) return {false, full(msg)};
        Call13 call_1;
        if (!unpack28(n22 + PACK77_NTOKENS, call_1)) success = false;
        std::array<char, 6> grid6 = blank<6>();  // undefined in the reference when to_grid fails
        if (!to_grid(igrid6, grid6)) success = false;
        set_msg(trim(call_1) + " " + full(grid6));
    } else {
        return {false, full(msg)};
    }
    if (std::memcmp(msg.data(), "CQ <", 4) == 0) success = false;
    if (success && have_staged) save_hash_call(staged);
    return {success, full(msg)};
}

std::string Unpacker::save() const {
    std::string t = "# FST4W callsign hashes: packjt77's 22-bit table, newest first\n";
    auto field = [](const std::string& v) { return v.empty() ? std::string("-") : v; };
    t += "mycall0 " + field(trim(mycall0_)) + "\n";
    t += "hiscall0 " + field(trim(hiscall0_)) + "\n";
    t += "mycall13_0 " + field(trim(mycall13_0_)) + "\n";
    for (int k = 0; k < nzhash_; ++k) t += "hash " + std::to_string(ihash22_[k]) + " " + trim(calls22_[k]) + "\n";
    return t;
}

bool Unpacker::load(const std::string& text) {
    Unpacker u;
    std::istringstream in(text);
    // A call field: no blanks, at most its length ("-" empty).
    auto take = [](const std::string& v, size_t n, auto& dst) {
        if (v == "-") return true;
        if (v.empty() || v.size() > n) return false;
        assign(dst, v);
        return true;
    };
    for (std::string l; std::getline(in, l);) {
        if (!l.empty() && l.back() == '\r') l.pop_back();
        if (l.empty() || l[0] == '#') continue;
        std::istringstream w(l);
        std::string key, a, b, more;
        w >> key >> a;
        if (key == "hash") {
            w >> b;
            if (w >> more || b.empty() || u.nzhash_ >= MAXHASH) return false;
            Call13 c = blank<13>();
            if (!take(b, 13, c)) return false;
            char* end = nullptr;
            const long h = std::strtol(a.c_str(), &end, 10);
            if (a.empty() || *end || h != ihashcall(c, 22)) return false;
            u.ihash22_[u.nzhash_] = static_cast<int>(h);
            u.calls22_[u.nzhash_] = c;
            ++u.nzhash_;
            continue;
        }
        if (w >> more) return false;
        if (key == "mycall0") {
            if (!take(a, 12, u.mycall0_)) return false;
        } else if (key == "hiscall0") {
            if (!take(a, 12, u.hiscall0_)) return false;
        } else if (key == "mycall13_0") {
            if (!take(a, 13, u.mycall13_0_)) return false;
        } else {
            return false;
        }
    }
    *this = u;
    return true;
}

}  // namespace fst4w
