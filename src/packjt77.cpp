// WSJT-X's 77-bit message packing (packjt77.h): lib/77bit/packjt77.f90,
// packjt77_schema.f90 and packjt77_grammar.f90 at 567ad29ce, translated
// routine by routine in the reference's order of effects; and genft8's
// packing half (lib/ft8/genft8.f90), ft8apset (lib/ft8/ft8apset.f90),
// split77, chkcall (lib/chkcall.f90) and stdcall (lib/qra/q65/
// q65_set_list.f90).
//
// Copyright 2018-2025 Joe Taylor, K1JT, Steven Franke, K9AN, and the WSJT-X
// developers; GNU GPL v3, and so is this.
//
// Strings are Fortran's CHARACTER values (FStr): a variable's length is
// fixed, assignment (set) cuts or blank pads to it, comparison ignores
// trailing blanks, substrings and index are 1-based. The routine names are
// the reference's; where one is left out it is said why.
#include "packjt77.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

namespace packjt77 {

namespace {

// ---------------------------------------------------------------- strings

// A Fortran CHARACTER value of length 0..CAP. Copy construction makes a
// value; Fortran assignment into a variable is set(), which keeps the
// variable's length (so plain = is not allowed).
class FStr {
public:
    static constexpr int CAP = 80;

    FStr() = default;
    FStr(const char* s) : FStr(s, static_cast<int>(std::strlen(s))) {}
    FStr(const char* s, int n) {
        if (n < 0 || n > CAP) std::abort();
        n_ = n;
        std::memcpy(c_, s, static_cast<size_t>(n));
    }
    explicit FStr(const std::string& s) : FStr(s.data(), static_cast<int>(s.size())) {}
    explicit FStr(const Call13& c) : FStr(c.data(), 13) {}
    FStr(const FStr&) = default;
    FStr& operator=(const FStr&) = delete;

    // A variable's declaration: character(len=n), blank.
    void declare(int n) {
        if (n < 0 || n > CAP) std::abort();
        n_ = n;
        std::memset(c_, ' ', static_cast<size_t>(n));
    }

    static FStr blank(int n) {
        FStr s;
        if (n < 0 || n > CAP) std::abort();
        s.n_ = n;
        std::memset(s.c_, ' ', static_cast<size_t>(n));
        return s;
    }

    int len() const { return n_; }
    char operator()(int i) const { return c_[i - 1]; }  // s(i:i)
    char& at(int i) { return c_[i - 1]; }
    const char* data() const { return c_; }
    std::string str() const { return std::string(c_, static_cast<size_t>(n_)); }

    // s(i:j); empty when j < i.
    FStr sub(int i, int j) const {
        if (j < i) return FStr();
        return FStr(c_ + i - 1, j - i + 1);
    }
    FStr from(int i) const { return sub(i, n_); }  // s(i:)

    // s = v
    void set(const FStr& v) {
        const int k = v.n_ < n_ ? v.n_ : n_;
        std::memmove(c_, v.c_, static_cast<size_t>(k));
        if (k < n_) std::memset(c_ + k, ' ', static_cast<size_t>(n_ - k));
    }
    // s(i:j) = v
    void put(int i, int j, const FStr& v) {
        if (j < i) return;
        const int n = j - i + 1;
        const int k = v.n_ < n ? v.n_ : n;
        std::memmove(c_ + i - 1, v.c_, static_cast<size_t>(k));
        if (k < n) std::memset(c_ + i - 1 + k, ' ', static_cast<size_t>(n - k));
    }

    friend FStr operator+(const FStr& a, const FStr& b) {
        FStr r;
        if (a.n_ + b.n_ > CAP) std::abort();
        r.n_ = a.n_ + b.n_;
        std::memcpy(r.c_, a.c_, static_cast<size_t>(a.n_));
        std::memcpy(r.c_ + a.n_, b.c_, static_cast<size_t>(b.n_));
        return r;
    }
    friend FStr operator+(const FStr& a, const char* b) { return a + FStr(b); }
    friend FStr operator+(const char* a, const FStr& b) { return FStr(a) + b; }

    // Fortran's a .eq. b: the shorter blank padded.
    friend bool operator==(const FStr& a, const FStr& b) {
        const int n = a.n_ > b.n_ ? a.n_ : b.n_;
        for (int i = 0; i < n; ++i) {
            const char x = i < a.n_ ? a.c_[i] : ' ';
            const char y = i < b.n_ ? b.c_[i] : ' ';
            if (x != y) return false;
        }
        return true;
    }
    friend bool operator!=(const FStr& a, const FStr& b) { return !(a == b); }
    friend bool operator==(const FStr& a, const char* b) { return a == FStr(b); }
    friend bool operator!=(const FStr& a, const char* b) { return !(a == FStr(b)); }

private:
    int n_ = 0;
    char c_[CAP];
};

int len_trim(const FStr& s) {
    int n = s.len();
    while (n > 0 && s(n) == ' ') --n;
    return n;
}

FStr trim(const FStr& s) { return s.sub(1, len_trim(s)); }

FStr adjustl(const FStr& s) {
    int k = 1;
    while (k <= s.len() && s(k) == ' ') ++k;
    FStr r = FStr::blank(s.len());
    r.set(s.from(k));
    return r;
}

FStr adjustr(const FStr& s) {
    const int n = len_trim(s);
    FStr r = FStr::blank(s.len());
    r.put(s.len() - n + 1, s.len(), s.sub(1, n));
    return r;
}

// index(s, t): 1-based, 0 when absent; index(s, '') is 1.
int index(const FStr& s, const FStr& t) {
    for (int i = 1; i + t.len() - 1 <= s.len(); ++i)
        if (std::memcmp(s.data() + i - 1, t.data(), static_cast<size_t>(t.len())) == 0) return i;
    return 0;
}
int index(const FStr& s, const char* t) { return index(s, FStr(t)); }
int index(const FStr& s, char c) { return index(s, FStr(&c, 1)); }

FStr chr(char c) { return FStr(&c, 1); }

Call13 to13(const FStr& s) {
    FStr t = FStr::blank(13);
    t.set(s);
    Call13 c;
    std::memcpy(c.data(), t.data(), 13);
    return c;
}

// write(c, '(iw.m)') v: right justified, at least m digits, asterisks when
// it does not fit.
FStr fmt_i(int v, int w, int m) {
    char b[32];
    std::snprintf(b, sizeof b, "%s%0*d", v < 0 ? "-" : "", m, v < 0 ? -v : v);
    const int n = static_cast<int>(std::strlen(b));
    if (n > w) return FStr(std::string(static_cast<size_t>(w), '*'));
    FStr r = FStr::blank(w);
    r.put(w - n + 1, w, FStr(b));
    return r;
}

// Fortran's nint on a default real: lroundf.
int nint(float x) { return static_cast<int>(std::lround(x)); }

// ---------------------------------------------------------------- grammar
// packjt77_grammar.f90

constexpr int PACK77_NSEC = 86;
constexpr int PACK77_NUSCAN = 171;
constexpr int PACK77_WSPR_NZZZ = 46656;
constexpr int PACK77_NTOKENS = 2063592;
constexpr int PACK77_MAX22 = 4194304;
constexpr int PACK77_MAXGRID4 = 32400;
const FStr PACK77_BASE36("0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ");
const FStr PACK77_BASE38(" 0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ/");

const char* const PACK77_ARRL_SECTIONS[PACK77_NSEC] = {
    "AB ", "AK ", "AL ", "AR ", "AZ ", "BC ", "CO ", "CT ", "DE ", "EB ",  //
    "EMA", "ENY", "EPA", "EWA", "GA ", "GH ", "IA ", "ID ", "IL ", "IN ",  //
    "KS ", "KY ", "LA ", "LAX", "NS ", "MB ", "MDC", "ME ", "MI ", "MN ",  //
    "MO ", "MS ", "MT ", "NC ", "ND ", "NE ", "NFL", "NH ", "NL ", "NLI",  //
    "NM ", "NNJ", "NNY", "TER", "NTX", "NV ", "OH ", "OK ", "ONE", "ONN",  //
    "ONS", "OR ", "ORG", "PAC", "PR ", "QC ", "RI ", "SB ", "SC ", "SCV",  //
    "SD ", "SDG", "SF ", "SFL", "SJV", "SK ", "SNJ", "STX", "SV ", "TN ",  //
    "UT ", "VA ", "VI ", "VT ", "WCF", "WI ", "WMA", "WNY", "WPA", "WTX",  //
    "WV ", "WWA", "WY ", "DX ", "PE ", "NB "};

const char* const PACK77_RTTY_MULTIPLIERS_STATES[72] = {
    "AL ", "AK ", "AZ ", "AR ", "CA ", "CO ", "CT ", "DE ", "FL ", "GA ",  //
    "HI ", "ID ", "IL ", "IN ", "IA ", "KS ", "KY ", "LA ", "ME ", "MD ",  //
    "MA ", "MI ", "MN ", "MS ", "MO ", "MT ", "NE ", "NV ", "NH ", "NJ ",  //
    "NM ", "NY ", "NC ", "ND ", "OH ", "OK ", "OR ", "PA ", "RI ", "SC ",  //
    "SD ", "TN ", "TX ", "UT ", "VT ", "VA ", "WA ", "WV ", "WI ", "WY ",  //
    "NB ", "NS ", "QC ", "ON ", "MB ", "SK ", "AB ", "BC ", "NWT", "NF ",  //
    "LB ", "NU ", "YT ", "PEI", "DC ", "DR ", "FR ", "GD ", "GR ", "OV ",  //
    "ZH ", "ZL "};

// PACK77_RTTY_MULTIPLIERS(i): the 72 above, then X01..X99.
FStr rtty_multiplier(int i) {
    if (i <= 72) return FStr(PACK77_RTTY_MULTIPLIERS_STATES[i - 1]);
    char b[16];
    std::snprintf(b, sizeof b, "X%02d", i - 72);
    return FStr(b);
}

int pack77_arrl_section_index(const FStr& section) {
    if (section.len() < 3) return -1;
    for (int i = 1; i <= PACK77_NSEC; ++i)
        if (FStr(PACK77_ARRL_SECTIONS[i - 1]) == section.sub(1, 3)) return i;
    return -1;
}

FStr pack77_arrl_section_name(int isec) {
    if (isec >= 1 && isec <= PACK77_NSEC) return FStr(PACK77_ARRL_SECTIONS[isec - 1]);
    return FStr("   ");
}

int pack77_rtty_multiplier_index(const FStr& multiplier) {
    if (multiplier.len() < 3) return -1;
    for (int i = 1; i <= PACK77_NUSCAN; ++i)
        if (rtty_multiplier(i) == multiplier.sub(1, 3)) return i;
    return -1;
}

FStr pack77_rtty_multiplier_name(int imult) {
    if (imult >= 1 && imult <= PACK77_NUSCAN) return rtty_multiplier(imult);
    return FStr("   ");
}

int pack77_report_index_from_snr(int isnr) {
    if (isnr < -50 || isnr > 50) return -1;
    // irpt 1..85 is -30..+50 dB; -50..-31 dB fold to irpt 86..105.
    int irpt = isnr;
    if (irpt >= -50 && irpt <= -31) irpt = irpt + 101;
    return irpt + 35;
}

int pack77_snr_from_report_index(int irpt) {
    int isnr = irpt - 35;
    if (isnr > 50) isnr = isnr - 101;
    return isnr;
}

// write(crpt,'(i3.2)') isnr, a leading blank made '+'.
FStr pack77_format_snr_report(int isnr) {
    FStr crpt = fmt_i(isnr, 3, 2);
    if (crpt(1) == ' ') crpt.at(1) = '+';
    return crpt;
}

int pack77_rtty_report_index(const FStr& report) {
    if (len_trim(report) != 3) return -1;
    if (report(1) != '5') return -1;
    if (report(2) < '2' || report(2) > '9') return -1;
    if (report(3) != '9') return -1;
    return report(2) - '2';
}

// format('5',i1,'9')
FStr pack77_format_rtty_report(int irpt) { return "5" + fmt_i(irpt + 2, 1, 1) + "9"; }

// format(i2,i4.4)
FStr pack77_format_vhf_exchange(int irpt, int iserial) { return fmt_i(52 + irpt, 2, 1) + fmt_i(iserial, 4, 4); }

FStr pack77_format_exchange_message(const FStr& call_1, const FStr& call_2, int ir, const FStr& exchange,
                                    bool reply_has_separate_token) {
    FStr msg = FStr::blank(37);
    msg.set(trim(call_1) + " " + trim(call_2) + " " + exchange);
    if (ir == 1) {
        if (reply_has_separate_token) msg.set(trim(call_1) + " " + trim(call_2) + " R " + exchange);
        else msg.set(trim(call_1) + " " + trim(call_2) + " R" + exchange);
    }
    return msg;
}

FStr pack77_format_type3_message(int itu, int ir, const FStr& call_1, const FStr& call_2, const FStr& crpt,
                                 const FStr& exchange) {
    FStr body = FStr::blank(37);
    body.set(pack77_format_exchange_message(call_1, call_2, ir, crpt + " " + exchange, true));
    FStr msg = FStr::blank(37);
    msg.set(body);
    if (itu == 1) msg.set("TU; " + trim(body));
    return msg;
}

int pack77_type3_message_len(int itu, int ir, int call_1_len, int call_2_len, int exchange_len) {
    int nmsg = call_1_len + 1 + call_2_len + 1 + 3 + 1 + exchange_len;
    if (ir == 1) nmsg = call_1_len + 1 + call_2_len + 3 + 3 + 1 + exchange_len;
    if (itu == 1) nmsg = nmsg + 4;
    return nmsg;
}

int pack77_type3_display_call_len(int n28, const FStr& call_text) {
    int ncall = len_trim(call_text);
    if (n28 >= PACK77_NTOKENS && n28 < PACK77_NTOKENS + PACK77_MAX22) {
        if (ncall > 11) ncall = -1;
        else ncall = ncall + 2;
    }
    return ncall;
}

bool pack77_type3_display_fits(int itu, int ir, int n28a, int n28b, const FStr& call_1, const FStr& call_2,
                               const FStr& exchange) {
    const int ncall_1 = pack77_type3_display_call_len(n28a, call_1);
    if (ncall_1 < 0) return false;
    const int ncall_2 = pack77_type3_display_call_len(n28b, call_2);
    if (ncall_2 < 0) return false;
    return pack77_type3_message_len(itu, ir, ncall_1, ncall_2, len_trim(exchange)) <= 37;
}

bool pack77_type3_render_fits(int itu, int ir, const FStr& call_1, const FStr& call_2, const FStr& exchange) {
    const int ncall_1 = len_trim(call_1);
    const int ncall_2 = len_trim(call_2);
    if (call_1(1) == '<' && index(call_1, '>') < 1) return false;
    if (call_2(1) == '<' && index(call_2, '>') < 1) return false;
    return pack77_type3_message_len(itu, ir, ncall_1, ncall_2, len_trim(exchange)) <= 37;
}

FStr pack77_format_qso_tail(const FStr& call_1, const FStr& call_2, int itail) {
    FStr msg = FStr::blank(37);
    msg.set(trim(call_1) + " " + trim(call_2));
    if (itail == 1) msg.set(trim(call_1) + " " + trim(call_2) + " RRR");
    if (itail == 2) msg.set(trim(call_1) + " " + trim(call_2) + " RR73");
    if (itail == 3) msg.set(trim(call_1) + " " + trim(call_2) + " 73");
    return msg;
}

int pack77_qso_tail_index(const FStr& token) {
    int itail = 0;
    if (trim(token) == "RRR") itail = 1;
    if (trim(token) == "RR73") itail = 2;
    if (trim(token) == "73") itail = 3;
    return itail;
}

int pack77_type12_irpt_from_tail(int itail) { return itail >= 1 && itail <= 3 ? itail + 1 : 0; }
int pack77_tail_from_type12_irpt(int irpt) { return irpt >= 2 && irpt <= 4 ? irpt - 1 : 0; }
int pack77_type4_nrpt_from_tail(int itail) { return itail >= 1 && itail <= 3 ? itail : 0; }
int pack77_tail_from_type4_nrpt(int nrpt) { return nrpt >= 1 && nrpt <= 3 ? nrpt : 0; }

bool pack77_is_digit_char(char c) { return c >= '0' && c <= '9'; }

bool pack77_all_digits(const FStr& token) {
    const int n = len_trim(token);
    if (n < 1) return false;
    for (int i = 1; i <= n; ++i)
        if (!pack77_is_digit_char(token(i))) return false;
    return true;
}

// read(token,*,err=10) value on a token of digits only (pack77_all_digits):
// gfortran's list-directed read fails past huge(0).
void pack77_parse_uint(const FStr& token, int& value, bool& ok) {
    value = 0;
    ok = pack77_all_digits(token);
    if (!ok) return;
    long long v = 0;
    for (int i = 1; i <= len_trim(token); ++i) {
        v = 10 * v + (token(i) - '0');
        if (v > 2147483647LL) {
            ok = false;
            return;
        }
    }
    value = static_cast<int>(v);
}

void pack77_parse_signed_report(const FStr& token, int& value, bool& ok) {
    value = 0;
    ok = false;
    if (len_trim(token) != 3) return;
    if (token(1) != '+' && token(1) != '-') return;
    if (!pack77_is_digit_char(token(2))) return;
    if (!pack77_is_digit_char(token(3))) return;
    value = 10 * (token(2) - '0') + (token(3) - '0');
    if (token(1) == '-') value = -value;
    ok = true;
}

void pack77_parse_dxpedition_report(const FStr& token, int& n5, bool& ok) {
    n5 = 0;
    int isnr;
    pack77_parse_signed_report(token, isnr, ok);
    if (!ok) return;
    if (isnr < -30 || isnr > 32) {
        ok = false;
        return;
    }
    n5 = (isnr + 30) / 2;
    ok = (2 * n5 - 30) == isnr;
}

void pack77_parse_type3_serial(const FStr& token, int& nserial, bool& ok) {
    nserial = 0;
    ok = false;
    if (len_trim(token) != 4) return;
    pack77_parse_uint(token, nserial, ok);
    if (!ok) return;
    ok = nserial >= 1 && nserial <= 7999;
}

void pack77_parse_wspr_dbm(const FStr& token, int& idbm, bool& ok) {
    idbm = 0;
    int dbm;
    pack77_parse_uint(token, dbm, ok);
    if (!ok) return;
    if (len_trim(token) > 1 && token(1) == '0') {
        ok = false;
        return;
    }
    if (dbm < 0 || dbm > 60) {
        ok = false;
        return;
    }
    idbm = nint(0.3f * static_cast<float>(dbm));
    const int decoded_dbm = nint(static_cast<float>(idbm) * 10.0f / 3.0f);
    ok = decoded_dbm == dbm;
}

void pack77_parse_field_day_exchange(const FStr& token, int& ntx, int& nclass, bool& ok) {
    ntx = -1;
    nclass = -1;
    ok = false;
    const int n = len_trim(token);
    if (n < 2) return;
    if (n > 2 && token(1) == '0') return;
    pack77_parse_uint(token.sub(1, n - 1), ntx, ok);
    if (!ok) return;
    if (ntx < 1 || ntx > 32) {
        ok = false;
        return;
    }
    nclass = static_cast<unsigned char>(token(n)) - 'A';
    ok = nclass >= 0 && nclass <= 5;
}

bool pack77_valid_hash_call_token(const FStr& token) {
    const int n = len_trim(token);
    if (n < 5 || n > 13) return false;
    if (token(1) != '<' || token(n) != '>') return false;
    if (index(token.sub(2, n - 1), '<') > 0) return false;
    if (index(token.sub(2, n - 1), '>') > 0) return false;
    return true;
}

bool pack77_valid_cq_modifier(const FStr& token) {
    const int n = len_trim(token);
    if (n < 1 || n > 4) return false;
    int nletters = 0, ndigits = 0;
    for (int i = 1; i <= n; ++i) {
        if (token(i) >= 'A' && token(i) <= 'Z') ++nletters;
        else if (pack77_is_digit_char(token(i))) ++ndigits;
        else return false;
    }
    return (ndigits == 3 && nletters == 0) || (nletters >= 1 && nletters <= 4 && ndigits == 0);
}

// pack77_split_source_tokens(msg, max_words, nwords, w, ok): w has size_w
// elements, each of its own (fixed) length.
void pack77_split_source_tokens(const FStr& msg, int max_words, int& nwords, FStr* w, int size_w, bool& ok) {
    nwords = 0;
    int n = 0;
    char c0 = ' ';
    for (int k = 0; k < size_w; ++k) w[k].set(FStr::blank(37));
    ok = true;
    const int iz = len_trim(msg);
    for (int i = 1; i <= iz; ++i) {
        char c = msg(i);
        if (c == '\0') c = ' ';
        if (c >= 'a' && c <= 'z') c = static_cast<char>(c - 32);
        if (c == ' ' && c0 == ' ') continue;
        if (c != ' ' && c0 == ' ') {
            ++nwords;
            n = 0;
            if (nwords > max_words || nwords > size_w) {
                ok = false;
                return;
            }
        }
        if (c != ' ') {
            ++n;
            if (n > w[0].len()) {
                ok = false;
                return;
            }
            w[nwords - 1].at(n) = c;
        }
        c0 = c;
    }
}

bool pack77_gate_staged_call_matches(const FStr& input_body, const Call13* staged_calls, int nstaged) {
    const int n = std::min(std::max(nstaged, 0), 32);
    for (int i = 1; i <= n; ++i)
        if (trim(FStr(staged_calls[i - 1])) == trim(input_body)) return true;
    return false;
}

void pack77_gate_message_checks(const FStr& input_msg, const FStr& decoded_msg, const Call13* staged_calls,
                                int nstaged, bool& messages_match, bool& resolved_render_fits) {
    FStr input_tokens[19], decoded_tokens[19];
    for (int k = 0; k < 19; ++k) {
        input_tokens[k].declare(37);
        decoded_tokens[k].declare(37);
    }
    FStr input_body = FStr::blank(37);
    int input_nwords, decoded_nwords;
    bool input_ok, decoded_ok;
    messages_match = false;
    resolved_render_fits = false;
    pack77_split_source_tokens(input_msg, 19, input_nwords, input_tokens, 19, input_ok);
    pack77_split_source_tokens(decoded_msg, 19, decoded_nwords, decoded_tokens, 19, decoded_ok);
    if (!input_ok || !decoded_ok) return;
    if (input_nwords != decoded_nwords) return;
    messages_match = true;
    resolved_render_fits = true;
    // nrender: the display length after hash resolution (a decoded <...>
    // renders as the staged call and its brackets).
    int nrender = 0;
    for (int i = 1; i <= input_nwords; ++i) {
        const FStr& dt = decoded_tokens[i - 1];
        const FStr& it = input_tokens[i - 1];
        if (i > 1) nrender = nrender + 1;
        bool token_match = false;
        if (trim(dt) == trim(it)) {
            token_match = true;
            nrender = nrender + len_trim(dt);
        } else {
            input_body.set(FStr::blank(37));
            input_body.set("<" + trim(it) + ">");
            if (trim(dt) == trim(input_body)) {
                token_match = true;
                nrender = nrender + len_trim(dt);
            } else if (trim(dt) == "<...>") {
                input_body.set(it);
                const int nbody = len_trim(input_body);
                if (nbody >= 2 && input_body(1) == '<' && input_body(nbody) == '>') {
                    input_body.set(FStr::blank(37));
                    if (nbody > 2) input_body.set(it.sub(2, nbody - 1));
                }
                token_match = pack77_gate_staged_call_matches(input_body, staged_calls, nstaged);
                if (len_trim(input_body) > 11) resolved_render_fits = false;
                nrender = nrender + len_trim(input_body) + 2;
            } else {
                nrender = nrender + len_trim(dt);
            }
        }
        if (!token_match) messages_match = false;
        if (trim(dt) == "<...>" && !token_match) resolved_render_fits = false;
        if (nrender > 37) resolved_render_fits = false;
    }
}

bool pack77_c28_standard_shape(const FStr& token) {
    const int n = len_trim(token);
    if (n < 3 || n > 6) return false;
    for (int i = 1; i <= n; ++i) {
        const char c = token(i);
        if (!((c >= 'A' && c <= 'Z') || pack77_is_digit_char(c))) return false;
    }
    int iarea = 0;
    for (int i = n; i >= 2; --i)
        if (pack77_is_digit_char(token(i))) {
            iarea = i;
            break;
        }
    if (iarea != 2 && iarea != 3) return false;
    int npdig = 0, nplet = 0;
    for (int i = 1; i <= iarea - 1; ++i) {
        const char c = token(i);
        if (pack77_is_digit_char(c)) ++npdig;
        else if (c >= 'A' && c <= 'Z') ++nplet;
        else return false;
    }
    if (nplet == 0 || npdig >= iarea - 1) return false;
    int nslet = 0;
    for (int i = iarea + 1; i <= n; ++i) {
        const char c = token(i);
        if (c >= 'A' && c <= 'Z') ++nslet;
        else return false;
    }
    return nslet <= 3;
}

bool pack77_type12_standard_call(const FStr& token) {
    if (!pack77_c28_standard_shape(token)) return false;
    // Q-prefix rejection is a source-grammar policy, not a c28 wire shape fact.
    return token(1) != 'Q';
}

enum { PACK77_TYPE12_SUFFIX_NONE = 0, PACK77_TYPE12_SUFFIX_R = 1, PACK77_TYPE12_SUFFIX_P = 2 };

struct Type12CallSource {
    FStr c28_token = FStr::blank(13);
    int suffix = PACK77_TYPE12_SUFFIX_NONE;
    bool is_hash = false;
    bool is_cq = false;
};

void pack77_parse_type12_call(const FStr& token, bool allow_special, Type12CallSource& call_source, bool& ok) {
    call_source.c28_token.set(FStr::blank(13));
    call_source.suffix = PACK77_TYPE12_SUFFIX_NONE;
    call_source.is_hash = false;
    call_source.is_cq = false;
    ok = false;
    const int n = len_trim(token);
    if (n < 1 || n > 13) return;

    if (allow_special) {
        if (trim(token) == "CQ") {
            call_source.c28_token.set("CQ");
            call_source.is_cq = true;
            ok = true;
            return;
        }
        if (trim(token) == "DE" || trim(token) == "QRZ") {
            call_source.c28_token.set(trim(token));
            ok = true;
            return;
        }
        if (n >= 4 && token.sub(1, 3) == "CQ_") {
            FStr modifier = FStr::blank(37);
            modifier.set(token.sub(4, n));
            if (!pack77_valid_cq_modifier(modifier)) return;
            call_source.c28_token.set(trim(token));
            call_source.is_cq = true;
            ok = true;
            return;
        }
    }

    if (token(1) == '<' || index(token.sub(1, n), '>') > 0) {
        if (!pack77_valid_hash_call_token(token)) return;
        call_source.c28_token.set(trim(token));
        call_source.is_hash = true;
        ok = true;
        return;
    }

    const int slash = index(token.sub(1, n), '/');
    FStr base = FStr::blank(37);
    if (slash > 0) {
        if (slash < 2 || slash != n - 1) return;
        if (index(token.sub(slash + 1, n), '/') > 0) return;
        if (token.sub(slash, n) == "/R") call_source.suffix = PACK77_TYPE12_SUFFIX_R;
        else if (token.sub(slash, n) == "/P") call_source.suffix = PACK77_TYPE12_SUFFIX_P;
        else return;
        base.set(token.sub(1, slash - 1));
    } else {
        base.set(token.sub(1, n));
    }

    if (!pack77_type12_standard_call(base)) return;
    call_source.c28_token.set(trim(base));
    ok = true;
}

bool pack77_type12_call_has_suffix(const Type12CallSource& s) { return s.suffix != PACK77_TYPE12_SUFFIX_NONE; }

bool pack77_is_grid4(const FStr& grid4) {
    if (grid4.len() < 4) return false;
    if (len_trim(grid4) != 4) return false;
    return grid4(1) >= 'A' && grid4(1) <= 'R' && grid4(2) >= 'A' && grid4(2) <= 'R' && grid4(3) >= '0' &&
           grid4(3) <= '9' && grid4(4) >= '0' && grid4(4) <= '9';
}

bool pack77_is_grid6(const FStr& grid6, bool allow_grid4) {
    if (grid6.len() < 4) return false;
    const int ngrid = len_trim(grid6);
    if (ngrid == 4) return allow_grid4 && pack77_is_grid4(grid6.sub(1, 4));
    if (grid6.len() < 6 || ngrid != 6) return false;
    return pack77_is_grid4(grid6.sub(1, 4)) && grid6(5) >= 'A' && grid6(5) <= 'X' && grid6(6) >= 'A' &&
           grid6(6) <= 'X';
}

int pack77_grid4_index(const FStr& grid4) {
    return (grid4(1) - 'A') * 18 * 10 * 10 + (grid4(2) - 'A') * 10 * 10 + (grid4(3) - '0') * 10 + (grid4(4) - '0');
}

int pack77_grid6_index(const FStr& grid6) {
    return pack77_grid4_index(grid6.sub(1, 4)) * 24 * 24 + (grid6(5) - 'A') * 24 + (grid6(6) - 'A');
}

struct Type12Source {
    Type12CallSource call_1, call_2;
    int ir = 0;
    int igrid4 = 0;
    int i3 = 1;
};

void pack77_parse_type12_source(const FStr& msg, Type12Source& source, bool& ok) {
    ok = false;
    FStr w[5] = {FStr::blank(37), FStr::blank(37), FStr::blank(37), FStr::blank(37), FStr::blank(37)};
    int nwords;
    bool split_ok;
    pack77_split_source_tokens(msg, 5, nwords, w, 5, split_ok);
    if (!split_ok) return;
    if (nwords < 2) return;

    FStr call_1_token = FStr::blank(37);
    call_1_token.set(w[0]);
    int call_2_index = 2, tail_start = 3;
    if (trim(w[0]) == "CQ" && nwords >= 3 && pack77_valid_cq_modifier(w[1])) {
        call_1_token.set(FStr::blank(37));
        call_1_token.set("CQ_" + trim(w[1]));
        call_2_index = 3;
        tail_start = 4;
    }

    bool call_ok;
    pack77_parse_type12_call(call_1_token, true, source.call_1, call_ok);
    if (!call_ok) return;
    pack77_parse_type12_call(w[call_2_index - 1], false, source.call_2, call_ok);
    if (!call_ok) return;
    if (source.call_1.is_cq && source.call_2.is_hash) return;
    if (source.call_1.is_hash &&
        (pack77_type12_call_has_suffix(source.call_2) || index(source.call_2.c28_token, '/') > 0))
        return;
    if (source.call_2.is_hash &&
        (pack77_type12_call_has_suffix(source.call_1) || index(source.call_1.c28_token, '/') > 0))
        return;

    if (source.call_1.suffix == PACK77_TYPE12_SUFFIX_R && source.call_2.suffix == PACK77_TYPE12_SUFFIX_P) return;
    if (source.call_1.suffix == PACK77_TYPE12_SUFFIX_P && source.call_2.suffix == PACK77_TYPE12_SUFFIX_R) return;

    source.i3 = 1;
    if (source.call_1.suffix == PACK77_TYPE12_SUFFIX_P || source.call_2.suffix == PACK77_TYPE12_SUFFIX_P)
        source.i3 = 2;

    const int ntail = nwords - tail_start + 1;
    if (ntail < 0 || ntail > 2) return;
    if (ntail == 0) {
        source.ir = 0;
        source.igrid4 = PACK77_MAXGRID4 + 1;
        ok = true;
        return;
    }

    if (ntail == 2) {
        if (trim(w[tail_start - 1]) != "R") return;
        if (source.call_1.is_cq) return;
        if (len_trim(w[tail_start]) != 4) return;
        if (!pack77_is_grid4(w[tail_start])) return;
        source.ir = 1;
        source.igrid4 = pack77_grid4_index(w[tail_start].sub(1, 4));
        ok = true;
        return;
    }

    FStr tail = FStr::blank(37);
    tail.set(w[tail_start - 1]);
    bool tail_is_grid = false;
    if (len_trim(tail) == 4) tail_is_grid = pack77_is_grid4(tail);
    if (tail_is_grid) {
        source.ir = 0;
        source.igrid4 = pack77_grid4_index(tail.sub(1, 4));
        ok = true;
        return;
    }
    if (source.call_1.is_cq) return;

    int isnr, irpt;
    if (tail(1) == '+' || tail(1) == '-') {
        pack77_parse_signed_report(tail, isnr, ok);
        if (!ok) return;
        irpt = pack77_report_index_from_snr(isnr);
        if (irpt < 1) {
            ok = false;
            return;
        }
        source.ir = 0;
        source.igrid4 = PACK77_MAXGRID4 + irpt;
        ok = true;
        return;
    }
    if (tail.sub(1, 2) == "R+" || tail.sub(1, 2) == "R-") {
        pack77_parse_signed_report(tail.from(2), isnr, ok);
        if (!ok) return;
        irpt = pack77_report_index_from_snr(isnr);
        if (irpt < 1) {
            ok = false;
            return;
        }
        source.ir = 1;
        source.igrid4 = PACK77_MAXGRID4 + irpt;
        ok = true;
        return;
    }

    const int itail = pack77_qso_tail_index(tail);
    if (itail < 1) return;
    source.ir = 0;
    source.igrid4 = PACK77_MAXGRID4 + pack77_type12_irpt_from_tail(itail);
    ok = true;
}

bool pack77_valid_c28_call_slot(const FStr& token) {
    const int n = len_trim(token);
    if (n < 3 || n > 13) return false;
    if (pack77_valid_hash_call_token(token)) return false;
    for (int i = 1; i <= n; ++i) {
        const char c = token(i);
        if (!((c >= 'A' && c <= 'Z') || pack77_is_digit_char(c) || c == '/')) return false;
    }
    if (pack77_type12_standard_call(token)) return true;
    const int slash = index(token.sub(1, n), '/');
    if (slash < 2 || slash >= n) return false;
    if (index(token.sub(slash + 1, n), '/') > 0) return false;
    return true;
}

struct DxpeditionSource {
    FStr call_1 = FStr::blank(13), call_2 = FStr::blank(13), hash_token = FStr::blank(13);
    int n5 = 0;
};

void pack77_parse_dxpedition_source(const FStr& msg, DxpeditionSource& source, bool& ok) {
    ok = false;
    FStr w[5] = {FStr::blank(37), FStr::blank(37), FStr::blank(37), FStr::blank(37), FStr::blank(37)};
    int nwords;
    bool split_ok;
    pack77_split_source_tokens(msg, 5, nwords, w, 5, split_ok);
    if (!split_ok) return;
    if (nwords != 5) return;
    if (trim(w[1]) != "RR73;") return;
    if (!pack77_valid_c28_call_slot(w[0])) return;
    if (!pack77_valid_c28_call_slot(w[2])) return;
    if (!pack77_valid_hash_call_token(w[3])) return;
    pack77_parse_dxpedition_report(w[4], source.n5, ok);
    if (!ok) return;
    source.call_1.set(trim(w[0]));
    source.call_2.set(trim(w[2]));
    source.hash_token.set(trim(w[3]));
}

struct FieldDaySource {
    FStr call_1 = FStr::blank(13), call_2 = FStr::blank(13);
    int ir = 0, ntx = 0, nclass = 0, isec = -1;
};

void pack77_parse_field_day_source(const FStr& msg, FieldDaySource& source, bool& ok) {
    ok = false;
    FStr w[5] = {FStr::blank(37), FStr::blank(37), FStr::blank(37), FStr::blank(37), FStr::blank(37)};
    int nwords;
    bool split_ok;
    pack77_split_source_tokens(msg, 5, nwords, w, 5, split_ok);
    if (!split_ok) return;
    if (nwords < 4 || nwords > 5) return;
    const int section_index = nwords;
    source.isec = pack77_arrl_section_index(w[section_index - 1]);
    if (source.isec == -1) return;
    const FStr section = pack77_arrl_section_name(source.isec);
    if (len_trim(w[section_index - 1]) != len_trim(section)) return;
    if (!pack77_valid_c28_call_slot(w[0])) return;
    if (!pack77_valid_c28_call_slot(w[1])) return;
    if (nwords == 5) {
        if (trim(w[2]) != "R") return;
        source.ir = 1;
    }
    const int exchange_index = nwords - 1;
    pack77_parse_field_day_exchange(w[exchange_index - 1], source.ntx, source.nclass, ok);
    if (!ok) return;
    source.call_1.set(trim(w[0]));
    source.call_2.set(trim(w[1]));
    ok = true;
}

bool pack77_wspr_affix_text(int npfx, FStr& cpfx, bool& is_prefix);

int pack77_wspr_affix_index(const FStr& token, int slash_index, int token_len, bool& ok) {
    ok = true;
    int npfx = 0;
    if (slash_index <= 4 && token_len - slash_index < 3) {
        ok = false;
        return npfx;
    }
    if (slash_index > 4 && token_len - slash_index > 3) {
        ok = false;
        return npfx;
    }
    for (int i = 1; i <= slash_index - 1; ++i)
        if (index(PACK77_BASE36, token(i)) < 1) {
            ok = false;
            return npfx;
        }
    for (int i = slash_index + 1; i <= token_len; ++i)
        if (index(PACK77_BASE36, token(i)) < 1) {
            ok = false;
            return npfx;
        }
    auto b36 = [&](int i) { return index(PACK77_BASE36, token(i)) - 1; };
    if (slash_index <= 4) {
        npfx = b36(1);
        if (slash_index >= 3) npfx = 36 * npfx + b36(2);
        if (slash_index == 4) npfx = 36 * npfx + b36(3);
    } else {
        if (token_len - slash_index == 1) npfx = b36(slash_index + 1);
        if (token_len - slash_index == 2) npfx = 36 * b36(slash_index + 1) + b36(slash_index + 2);
        if (token_len - slash_index == 3) {
            if (!pack77_is_digit_char(token(slash_index + 3))) {
                ok = false;
                return npfx;
            }
            npfx = 36 * 10 * b36(slash_index + 1) + 10 * b36(slash_index + 2) + b36(slash_index + 3);
        }
        npfx = npfx + PACK77_WSPR_NZZZ;
    }
    return npfx;
}

bool pack77_wspr_affix_text(int npfx, FStr& cpfx, bool& is_prefix) {
    int n = npfx;
    cpfx.set("   ");
    is_prefix = n < PACK77_WSPR_NZZZ;
    if (is_prefix) {
        for (int i = 3; i >= 1; --i) {
            const int j = n % 36 + 1;
            cpfx.at(i) = PACK77_BASE36(j);
            n = n / 36;
            if (n == 0) break;
        }
    } else {
        n = n - PACK77_WSPR_NZZZ;
        if (n <= 35) {
            cpfx.at(1) = PACK77_BASE36(n + 1);
        } else if (n > 35 && n <= 1295) {
            cpfx.at(1) = PACK77_BASE36(n / 36 + 1);
            cpfx.at(2) = PACK77_BASE36(n % 36 + 1);
        } else if (n > 1295 && n <= 12959) {
            cpfx.at(1) = PACK77_BASE36(n / 360 + 1);
            cpfx.at(2) = PACK77_BASE36((n / 10) % 36 + 1);
            cpfx.at(3) = PACK77_BASE36(n % 10 + 1);
        } else {
            return false;
        }
    }
    return true;
}

struct WsprSource {
    FStr call_token = FStr::blank(13), hash_token = FStr::blank(13);
    FStr base_call = FStr::blank(6), grid6 = FStr::blank(6);
    FStr grid4 = FStr::blank(4);
    int subtype = 0, idbm = 0, npfx = 0;
};

void pack77_parse_wspr_source(const FStr& msg, WsprSource& source, bool& ok) {
    ok = false;
    FStr w[3] = {FStr::blank(37), FStr::blank(37), FStr::blank(37)};
    int nwords;
    bool split_ok;
    pack77_split_source_tokens(msg, 3, nwords, w, 3, split_ok);
    if (!split_ok) return;
    if (nwords < 2 || nwords > 3) return;
    const int m1 = len_trim(w[0]);
    const int m2 = len_trim(w[1]);

    bool grid4_ok = false;
    if (nwords == 3) grid4_ok = m2 == 4 && pack77_is_grid4(w[1]);

    if (nwords == 3) {
        if (!pack77_type12_standard_call(w[0])) return;
        if (!grid4_ok) return;
        pack77_parse_wspr_dbm(w[2], source.idbm, ok);
        if (!ok) return;
        source.subtype = 1;
        source.call_token.set(trim(w[0]));
        source.grid4.set(w[1].sub(1, 4));
        ok = true;
        return;
    }

    if (m1 >= 5 && m1 <= 10 && m2 <= 2) {
        const int slash = index(w[0].sub(1, m1), '/');
        if (slash >= 2 && slash < m1 && index(w[0].sub(slash + 1, m1), '/') == 0) {
            FStr base = FStr::blank(37);
            if (slash <= 4) base.set(w[0].sub(slash + 1, m1));
            else base.set(w[0].sub(1, slash - 1));
            if (pack77_type12_standard_call(base)) {
                pack77_parse_wspr_dbm(w[1], source.idbm, ok);
                if (!ok) return;
                source.npfx = pack77_wspr_affix_index(w[0], slash, m1, ok);
                if (!ok) return;
                FStr cpfx = FStr::blank(3);
                bool is_prefix;
                ok = pack77_wspr_affix_text(source.npfx, cpfx, is_prefix);
                if (!ok) return;
                ok = false;
                cpfx.set(adjustl(cpfx));
                if (slash <= 4) {
                    if (!is_prefix) return;
                    if (trim(cpfx) != w[0].sub(1, slash - 1)) return;
                } else {
                    if (is_prefix) return;
                    if (trim(cpfx) != w[0].sub(slash + 1, m1)) return;
                }
                source.subtype = 2;
                source.call_token.set(trim(w[0]));
                source.base_call.set(trim(base));
                ok = true;
                return;
            }
        }
    }

    if (m1 >= 5 && m1 <= 12 && m2 <= 6) {
        if (!pack77_valid_hash_call_token(w[0])) return;
        if (!pack77_is_grid6(w[1], true)) return;
        source.subtype = 3;
        source.hash_token.set(trim(w[0]));
        source.grid6.set(w[1].sub(1, 6));
        ok = true;
    }
}

struct Type3Source {
    FStr call_1 = FStr::blank(13), call_2 = FStr::blank(13);
    int itu = 0, ir = 0, irpt = -1, nexch = 0;
};

void pack77_parse_type3_source(const FStr& msg, Type3Source& source, bool& ok) {
    ok = false;
    FStr w[6] = {FStr::blank(37), FStr::blank(37), FStr::blank(37),
                 FStr::blank(37), FStr::blank(37), FStr::blank(37)};
    int nwords;
    bool split_ok;
    pack77_split_source_tokens(msg, 6, nwords, w, 6, split_ok);
    if (!split_ok) return;
    if (nwords < 4 || nwords > 6) return;

    if (trim(w[0]) == "TU;") source.itu = 1;
    if (trim(w[3 + source.itu - 1]) == "R") source.ir = 1;
    const int report_index = 3 + source.itu + source.ir;
    if (!pack77_valid_c28_call_slot(w[1 + source.itu - 1])) return;
    if (!pack77_valid_c28_call_slot(w[2 + source.itu - 1])) return;
    source.call_1.set(trim(w[1 + source.itu - 1]));
    source.call_2.set(trim(w[2 + source.itu - 1]));

    const int exchange_index = report_index + 1;
    if (exchange_index != nwords) return;

    source.irpt = pack77_rtty_report_index(w[report_index - 1]);
    if (source.irpt < 0) return;

    int nserial;
    bool serial_ok;
    pack77_parse_type3_serial(w[exchange_index - 1], nserial, serial_ok);
    if (serial_ok) {
        source.nexch = nserial;
        ok = true;
        return;
    }

    const int imult = pack77_rtty_multiplier_index(w[exchange_index - 1]);
    if (imult <= 0) return;
    const FStr exchange = pack77_rtty_multiplier_name(imult);
    if (len_trim(w[exchange_index - 1]) != len_trim(exchange)) return;
    source.nexch = 8000 + imult;
    ok = true;
}

bool pack77_valid_type4_c11(const FStr& token) {
    const int n = len_trim(token);
    if (n < 1 || n > 11) return false;
    for (int i = 1; i <= n; ++i)
        if (index(PACK77_BASE38, token(i)) < 1) return false;
    return true;
}

struct Type4Source {
    FStr hash_token = FStr::blank(13);
    FStr c11 = FStr::blank(11);
    int icq = 0, iflip = 0, itail = 0;
};

void pack77_parse_type4_source(const FStr& msg, Type4Source& source, bool& ok) {
    ok = false;
    FStr w[3] = {FStr::blank(37), FStr::blank(37), FStr::blank(37)};
    int nwords;
    bool split_ok;
    pack77_split_source_tokens(msg, 3, nwords, w, 3, split_ok);
    if (!split_ok) return;
    if (nwords < 2 || nwords > 3) return;

    const bool hash_1 = pack77_valid_hash_call_token(w[0]);
    const bool hash_2 = pack77_valid_hash_call_token(w[1]);

    if (trim(w[0]) == "CQ") {
        if (nwords != 2) return;
        if (hash_2) return;
        const FStr& c11_token = w[1];
        if (len_trim(c11_token) < 5) return;
        if (!pack77_valid_type4_c11(c11_token)) return;
        if (pack77_type12_standard_call(c11_token)) return;
        source.c11.set(c11_token.sub(1, 11));
        source.icq = 1;
        ok = true;
        return;
    }

    if (nwords == 3) {
        source.itail = pack77_qso_tail_index(w[2]);
        if (source.itail == 0) return;
    }

    if (hash_1 && !hash_2) {
        if (!pack77_valid_type4_c11(w[1])) return;
        if (pack77_type12_standard_call(w[1])) return;
        source.hash_token.set(trim(w[0]));
        source.c11.set(w[1].sub(1, 11));
        source.iflip = 0;
        ok = true;
        return;
    }

    if (hash_2 && !hash_1) {
        if (!pack77_valid_type4_c11(w[0])) return;
        if (pack77_type12_standard_call(w[0])) return;
        source.hash_token.set(trim(w[1]));
        source.c11.set(w[0].sub(1, 11));
        source.iflip = 1;
        ok = true;
        return;
    }
}

struct Type5Source {
    FStr hash_token_1 = FStr::blank(13), hash_token_2 = FStr::blank(13);
    FStr grid6 = FStr::blank(6);
    int ir = 0, irpt = 0, iserial = 0;
};

void pack77_valid_type5_exchange(const FStr& token, int& nx, bool& ok) {
    nx = -1;
    ok = false;
    if (len_trim(token) != 6) return;
    pack77_parse_uint(token, nx, ok);
}

void pack77_parse_type5_source(const FStr& msg, Type5Source& source, bool& ok) {
    ok = false;
    FStr w[5] = {FStr::blank(37), FStr::blank(37), FStr::blank(37), FStr::blank(37), FStr::blank(37)};
    int nwords;
    bool split_ok;
    pack77_split_source_tokens(msg, 5, nwords, w, 5, split_ok);
    if (!split_ok) return;
    if (nwords < 4 || nwords > 5) return;
    if (!pack77_valid_hash_call_token(w[0])) return;
    if (!pack77_valid_hash_call_token(w[1])) return;
    if (nwords == 5) {
        if (trim(w[2]) != "R") return;
        source.ir = 1;
    }
    const int exchange_index = nwords - 1;
    const int grid_index = nwords;
    int nx;
    pack77_valid_type5_exchange(w[exchange_index - 1], nx, ok);
    if (!ok) return;
    if (nx < 520001 || nx > 592047) {
        ok = false;
        return;
    }
    source.iserial = nx % 10000;
    if (source.iserial < 1 || source.iserial > 2047) {
        ok = false;
        return;
    }
    if (len_trim(w[grid_index - 1]) != 6) {
        ok = false;
        return;
    }
    if (!pack77_is_grid6(w[grid_index - 1], false)) {
        ok = false;
        return;
    }
    source.hash_token_1.set(trim(w[0]));
    source.hash_token_2.set(trim(w[1]));
    source.irpt = nx / 10000 - 52;
    source.grid6.set(w[grid_index - 1].sub(1, 6));
    ok = true;
}

FStr pack77_type4_call_text_from_n58(int64_t n58) {
    FStr c11 = FStr::blank(11);
    int64_t n = n58;
    for (int i = 11; i >= 1; --i) {
        const int j = static_cast<int>(n % 38) + 1;
        c11.at(i) = PACK77_BASE38(j);
        n = n / 38;
    }
    return c11;
}

int64_t pack77_type4_n58_from_call_text(const FStr& c11) {
    int64_t n58 = 0;
    for (int i = 1; i <= 11; ++i) n58 = n58 * 38 + index(PACK77_BASE38, c11(i)) - 1;
    return n58;
}

bool to_grid4(int n, FStr& grid4) {
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
    grid4.at(1) = static_cast<char>(j1 + 'A');
    grid4.at(2) = static_cast<char>(j2 + 'A');
    grid4.at(3) = static_cast<char>(j3 + '0');
    grid4.at(4) = static_cast<char>(j4 + '0');
    return true;
}

bool to_grid6(int n, FStr& grid6) {
    const int j1 = n / (18 * 10 * 10 * 24 * 24);
    if (j1 < 0 || j1 > 17) return false;
    n = n - j1 * 18 * 10 * 10 * 24 * 24;
    const int j2 = n / (10 * 10 * 24 * 24);
    if (j2 < 0 || j2 > 17) return false;
    n = n - j2 * 10 * 10 * 24 * 24;
    const int j3 = n / (10 * 24 * 24);
    if (j3 < 0 || j3 > 9) return false;
    n = n - j3 * 10 * 24 * 24;
    const int j4 = n / (24 * 24);
    if (j4 < 0 || j4 > 9) return false;
    n = n - j4 * 24 * 24;
    const int j5 = n / 24;
    if (j5 < 0 || j5 > 23) return false;
    const int j6 = n - j5 * 24;
    if (j6 < 0 || j6 > 23) return false;
    grid6.at(1) = static_cast<char>(j1 + 'A');
    grid6.at(2) = static_cast<char>(j2 + 'A');
    grid6.at(3) = static_cast<char>(j3 + '0');
    grid6.at(4) = static_cast<char>(j4 + '0');
    grid6.at(5) = static_cast<char>(j5 + 'A');
    grid6.at(6) = static_cast<char>(j6 + 'A');
    return true;
}

bool to_grid(int n, FStr& grid6) {
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
    grid6.set(FStr::blank(6));
    grid6.at(1) = static_cast<char>(j1 + 'A');
    grid6.at(2) = static_cast<char>(j2 + 'A');
    grid6.at(3) = static_cast<char>(j3 + '0');
    grid6.at(4) = static_cast<char>(j4 + '0');
    if (j5 != 24 || j6 != 24) {
        grid6.at(5) = static_cast<char>(j5 + 'A');
        grid6.at(6) = static_cast<char>(j6 + 'A');
    }
    return true;
}

// ---------------------------------------------------------------- schema
// packjt77_schema.f90: the 77-bit layouts, MSB first, 1-based. get_int's
// range checks are the decode gates; put_int's the encode ones.

using C77 = FStr;  // character*77 of '0' and '1'

int64_t field(const C77& c77, int start, int width) {
    int64_t v = 0;
    for (int i = 0; i < width; ++i) v = 2 * v + (c77(start + i) == '1');
    return v;
}

// get_int(ok, c77, field, value) with the field's min/max.
void get_int(bool& ok, const C77& c77, int start, int width, int64_t vmin, int64_t vmax, int& value) {
    value = 0;
    if (!ok) return;
    const int64_t v = field(c77, start, width);
    ok = v >= vmin && v <= vmax;
    if (ok) value = static_cast<int>(v);
}

void require_int(bool& ok, const C77& c77, int start, int width, int expected) {
    if (!ok) return;
    ok = field(c77, start, width) == expected;
}

void put_int(bool& ok, C77& c77, int start, int width, int64_t vmin, int64_t vmax, int64_t value) {
    if (!ok) return;
    ok = false;
    if (value < vmin || value > vmax) return;
    for (int i = width; i >= 1; --i) {
        const int64_t place = int64_t(1) << (width - i);
        c77.at(start + i - 1) = ((value / place) % 2 == 1) ? '1' : '0';
    }
    ok = true;
}

bool bits_are_binary(const C77& c77) {
    for (int i = 1; i <= c77.len(); ++i)
        if (c77(i) != '0' && c77(i) != '1') return false;
    return true;
}

C77 begin_message() { return C77(std::string(77, '0')); }

constexpr int64_t M28 = 268435455;

bool encode_pack77_free_text(const FStr& text_bits, C77& c77) {
    c77.set(begin_message());
    // put_bits: width 71, all '0'/'1' (packtext77's output always is).
    c77.put(1, 71, text_bits);
    return true;
}

bool encode_pack77_telemetry(int ntel1, int ntel2, int ntel3, C77& c77) {
    c77.set(begin_message());
    bool ok = true;
    put_int(ok, c77, 1, 23, 0, 8388607, ntel1);
    put_int(ok, c77, 24, 24, 0, 16777215, ntel2);
    put_int(ok, c77, 48, 24, 0, 16777215, ntel3);
    put_int(ok, c77, 72, 3, 5, 5, 5);
    put_int(ok, c77, 75, 3, 0, 0, 0);
    return ok;
}

bool encode_pack77_dxpedition(int n28a, int n28b, int n10, int n5, C77& c77) {
    c77.set(begin_message());
    bool ok = true;
    put_int(ok, c77, 1, 28, 0, M28, n28a);
    put_int(ok, c77, 29, 28, 0, M28, n28b);
    put_int(ok, c77, 57, 10, 0, 1023, n10);
    put_int(ok, c77, 67, 5, 0, 31, n5);
    put_int(ok, c77, 72, 3, 1, 1, 1);
    put_int(ok, c77, 75, 3, 0, 0, 0);
    return ok;
}

bool encode_pack77_field_day(int n3, int n28a, int n28b, int ir, int intx, int nclass, int isec, C77& c77) {
    c77.set(begin_message());
    bool ok = true;
    put_int(ok, c77, 1, 28, 0, M28, n28a);
    put_int(ok, c77, 29, 28, 0, M28, n28b);
    put_int(ok, c77, 57, 1, 0, 1, ir);
    put_int(ok, c77, 58, 4, 0, 15, intx);
    put_int(ok, c77, 62, 3, 0, 7, nclass);
    put_int(ok, c77, 65, 7, 0, 127, isec);
    put_int(ok, c77, 72, 3, n3, n3, n3);
    put_int(ok, c77, 75, 3, 0, 0, 0);
    return ok;
}

bool encode_pack77_wspr_type1(int n28, int igrid4, int idbm, C77& c77) {
    c77.set(begin_message());
    bool ok = true;
    put_int(ok, c77, 1, 28, 0, M28, n28);
    put_int(ok, c77, 29, 15, 0, 32767, igrid4);
    put_int(ok, c77, 44, 5, 0, 18, idbm);
    put_int(ok, c77, 49, 1, 0, 0, 0);
    put_int(ok, c77, 50, 1, 0, 0, 0);
    put_int(ok, c77, 51, 21, 0, 0, 0);
    put_int(ok, c77, 72, 3, 6, 6, 6);
    put_int(ok, c77, 75, 3, 0, 0, 0);
    return ok;
}

bool encode_pack77_wspr_type2(int n28, int npfx, int idbm, C77& c77) {
    c77.set(begin_message());
    bool ok = true;
    put_int(ok, c77, 1, 28, 0, M28, n28);
    put_int(ok, c77, 29, 16, 0, 65535, npfx);
    put_int(ok, c77, 45, 5, 0, 18, idbm);
    put_int(ok, c77, 50, 1, 1, 1, 1);
    put_int(ok, c77, 51, 21, 0, 0, 0);
    put_int(ok, c77, 72, 3, 6, 6, 6);
    put_int(ok, c77, 75, 3, 0, 0, 0);
    return ok;
}

bool encode_pack77_type12(int i3, int n28a, int ipa, int n28b, int ipb, int ir, int igrid4, C77& c77) {
    c77.set(begin_message());
    bool ok = true;
    put_int(ok, c77, 1, 28, 0, M28, n28a);
    put_int(ok, c77, 29, 1, 0, 1, ipa);
    put_int(ok, c77, 30, 28, 0, M28, n28b);
    put_int(ok, c77, 58, 1, 0, 1, ipb);
    put_int(ok, c77, 59, 1, 0, 1, ir);
    put_int(ok, c77, 60, 15, 0, 32767, igrid4);
    put_int(ok, c77, 75, 3, i3, i3, i3);
    return ok;
}

bool encode_pack77_type3(int itu, int n28a, int n28b, int ir, int irpt, int nexch, C77& c77) {
    c77.set(begin_message());
    bool ok = true;
    put_int(ok, c77, 1, 1, 0, 1, itu);
    put_int(ok, c77, 2, 28, 0, M28, n28a);
    put_int(ok, c77, 30, 28, 0, M28, n28b);
    put_int(ok, c77, 58, 1, 0, 1, ir);
    put_int(ok, c77, 59, 3, 0, 7, irpt);
    put_int(ok, c77, 62, 13, 0, 8191, nexch);
    put_int(ok, c77, 75, 3, 3, 3, 3);
    return ok;
}

bool encode_pack77_type4(int n12, int64_t n58, int iflip, int nrpt, int icq, C77& c77) {
    c77.set(begin_message());
    bool ok = true;
    put_int(ok, c77, 1, 12, 0, 4095, n12);
    put_int(ok, c77, 13, 58, 0, 288230376151711743LL, n58);
    put_int(ok, c77, 71, 1, 0, 1, iflip);
    put_int(ok, c77, 72, 2, 0, 3, nrpt);
    put_int(ok, c77, 74, 1, 0, 1, icq);
    put_int(ok, c77, 75, 3, 4, 4, 4);
    return ok;
}

bool encode_pack77_type5(int n12, int n22, int ir, int irpt, int iserial, int igrid6, C77& c77) {
    c77.set(begin_message());
    bool ok = true;
    put_int(ok, c77, 1, 12, 0, 4095, n12);
    put_int(ok, c77, 13, 22, 0, 4194303, n22);
    put_int(ok, c77, 35, 1, 0, 1, ir);
    put_int(ok, c77, 36, 3, 0, 7, irpt);
    put_int(ok, c77, 39, 11, 1, 2047, iserial);
    put_int(ok, c77, 50, 25, 0, 18662399, igrid6);
    put_int(ok, c77, 75, 3, 5, 5, 5);
    return ok;
}

// decode_pack77_tag: i3 from bits 75-77, n3 (i3 = 0) from 72-74; fails
// unless every character is '0' or '1'.
bool decode_pack77_tag(const C77& c77, int& i3, int& n3) {
    i3 = -1;
    n3 = -1;
    if (!bits_are_binary(c77)) return false;
    i3 = static_cast<int>(field(c77, 75, 3));
    n3 = i3 == 0 ? static_cast<int>(field(c77, 72, 3)) : 0;
    return true;
}

// pack77_wspr_payload_type
int pack77_wspr_payload_type(const C77& c77) {
    if (c77(50) == '1') return 2;
    if (c77(49) == '0') return 1;
    if (c77(49) == '1') return 3;  // invalid selector 110 too: the schema rejects it
    return 0;
}

// ---------------------------------------------------------------- text
// packtext77 / unpacktext77: 13 characters of ' 0-9A-Z+-./?' as a 71-bit
// base-42 number (the reference's mp_short_* byte arithmetic is exactly
// this; 42^13 < 2^71).

const FStr TEXT_ALPHABET(" 0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ+-./?");

FStr packtext77(const FStr& c13) {
    const FStr w = adjustr(c13);
    unsigned __int128 v = 0;
    for (int i = 1; i <= 13; ++i) {
        int j = index(TEXT_ALPHABET, w(i)) - 1;
        if (j < 0) j = 0;
        v = 42 * v + static_cast<unsigned>(j);
    }
    FStr c71 = FStr::blank(71);
    for (int i = 71; i >= 1; --i) {
        c71.at(i) = (v & 1) ? '1' : '0';
        v >>= 1;
    }
    return c71;
}

FStr unpacktext77(const FStr& c71) {
    unsigned __int128 v = 0;
    for (int i = 1; i <= 71; ++i) v = 2 * v + (c71(i) == '1');
    FStr c13 = FStr::blank(13);
    for (int i = 13; i >= 1; --i) {
        const int ir = static_cast<int>(v % 42);
        c13.at(i) = TEXT_ALPHABET(ir + 1);
        v = v / 42;
    }
    return c13;
}

// pack77_normalized_message: NUL to blank, upper case, blanks collapsed
// (as pack77_split_source_tokens canonicalises).
FStr pack77_normalized_message(const FStr& msg) {
    FStr normalized = FStr::blank(37);
    char c0 = ' ';
    int j = 0;
    const int iz = len_trim(msg);
    for (int i = 1; i <= std::min(iz, 37); ++i) {
        char c = msg(i);
        if (c == '\0') c = ' ';
        if (c == ' ' && c0 == ' ') continue;
        if (c >= 'a' && c <= 'z') c = static_cast<char>(c - 32);
        ++j;
        normalized.at(j) = c;
        c0 = c;
    }
    return normalized;
}

bool pack77_free_text_alphabet_ok(const FStr& msg) {
    for (int i = 1; i <= len_trim(msg); ++i)
        if (index(TEXT_ALPHABET, msg(i)) == 0) return false;
    return true;
}

bool pack77_exact_free_text_ok(const FStr& msg) {
    if (len_trim(msg) > 13) return false;
    return pack77_free_text_alphabet_ok(msg.sub(1, len_trim(msg)));
}

int pack77_free_text_reject_status(const FStr& msg) {
    if (!pack77_free_text_alphabet_ok(msg)) return PACK77_STATUS_FREE_TEXT_INVALID;
    if (len_trim(msg) > 13) return PACK77_STATUS_FREE_TEXT_TOO_LONG;
    return PACK77_STATUS_FREE_TEXT_INVALID;
}

// ---------------------------------------------------------------- results

struct PackResult {
    bool encoded = false;
    int status = PACK77_STATUS_NOT_ENCODED;
    int i3 = -1, n3 = -1;
    C77 c77 = FStr::blank(77);
};

PackResult pack77_no_match() { return PackResult(); }

PackResult pack77_reject(int status) {
    PackResult r;
    r.status = status;
    return r;
}

PackResult pack77_accept(int i3, int n3, const C77& c77) {
    PackResult r;
    r.encoded = true;
    r.status = PACK77_STATUS_ENCODED;
    r.i3 = i3;
    r.n3 = n3;
    r.c77.set(c77);
    return r;
}

void copy_result(PackResult& dst, const PackResult& src) {
    dst.encoded = src.encoded;
    dst.status = src.status;
    dst.i3 = src.i3;
    dst.n3 = src.n3;
    dst.c77.set(src.c77);
}

// The calls an accepted encode makes resolvable (pack77_hash_facts).
struct HashFacts {
    int ncalls = 0;
    Call13 calls[32];
};

struct Candidate {
    PackResult result;
    HashFacts hash_facts;
};

struct EncodeCore {
    PackResult result;
    HashFacts hash_facts;
};

// Hash and recent-call side effects a decode collects, applied only when
// it succeeds (unpack77_effects).
struct Effects {
    int nhash = 0, nrecent = 0;
    Call13 hash_calls[32];
    Call13 recent_calls[32];
};

struct Context {
    int nrx = 0;
    bool mycall_set = false, dxcall_set = false;
    FStr mycall = FStr::blank(13), dxcall = FStr::blank(13);
    int hashmy10 = -1, hashmy12 = -1, hashmy22 = -1, hashdx10 = -1;
};

struct CoreResult {
    bool success = false;
    int i3 = -1, n3 = -1;
    FStr msg = FStr::blank(37);
    Effects effects;
};

// normalize_hash_call: the call a "<CALL>" token or a plain call stages.
bool normalize_hash_call(const FStr& c13, FStr& cw) {
    cw.set(c13);
    if (cw(1) == ' ' || cw.sub(1, 5) == "<...>") return false;
    if (cw(1) == '<') {
        if (!pack77_valid_hash_call_token(cw)) return false;
        cw.set(FStr(cw).from(2));
        const int i = index(cw, '>');
        if (i > 0) cw.put(i, 13, FStr("         "));
    }
    return len_trim(cw) >= 3;
}

void stage_hash_call(HashFacts& hash_facts, const FStr& c13) {
    FStr cw = FStr::blank(13);
    if (!normalize_hash_call(c13, cw)) return;
    if (hash_facts.ncalls >= 32) return;
    hash_facts.calls[hash_facts.ncalls++] = to13(cw);
}

void stage_unpack_call(Call13* calls, int& ncalls, const FStr& c13) {
    FStr cw = FStr::blank(13);
    if (!normalize_hash_call(c13, cw)) return;
    if (ncalls >= 32) return;
    calls[ncalls++] = to13(cw);
}
void stage_unpack_hash_call(Effects& e, const FStr& c13) { stage_unpack_call(e.hash_calls, e.nhash, c13); }
void stage_unpack_recent_call(Effects& e, const FStr& c13) { stage_unpack_call(e.recent_calls, e.nrecent, c13); }

int ihashcall_f(const FStr& c0, int m) {
    int64_t n8 = 0;
    for (int i = 1; i <= 11; ++i) {
        const int64_t j = index(PACK77_BASE38, c0(i)) - 1;
        n8 = 38 * n8 + j;
    }
    // ihashcall_from_n8: the top m bits of the low 64 bits of
    // 47055833459*n8, in 16-bit limbs with Fortran's truncating mod and
    // division.
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

void pack77_hash_index(HashFacts& hash_facts, const FStr& token, int nbits, int& n, bool& ok) {
    FStr hash_call = FStr::blank(13);
    ok = normalize_hash_call(token, hash_call);
    if (!ok) return;
    n = ihashcall_f(hash_call, nbits);
    stage_hash_call(hash_facts, hash_call);
}

// pack28_core: a special token, a 22-bit hash or a standard call as c28.
// c13 is changed as the reference changes its argument (a "<CALL>" token
// left as CALL). The '3DA0' and '3X' workarounds set a variable every path
// below overwrites or never reads: left out.
void pack28_core(FStr& c13, int& n28, HashFacts& hash_facts) {
    static const FStr a1(" 0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ");
    static const FStr a2("0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ");
    static const FStr a3("0123456789");
    static const FStr a4(" ABCDEFGHIJKLMNOPQRSTUVWXYZ");
    n28 = -1;
    if (c13.sub(1, 3) == "DE ") {
        n28 = 0;
    } else if (c13.sub(1, 4) == "QRZ ") {
        n28 = 1;
    } else if (c13.sub(1, 3) == "CQ ") {
        n28 = 2;
    } else {
        bool done = false;
        if (c13.sub(1, 3) == "CQ_") {
            const int n = len_trim(c13);
            if (n >= 4 && pack77_valid_cq_modifier(c13.sub(4, n))) {
                if (pack77_all_digits(c13.sub(4, n))) {
                    const int nqsy = 100 * (c13(4) - '0') + 10 * (c13(5) - '0') + (c13(6) - '0');
                    n28 = 3 + nqsy;
                } else {
                    FStr c4 = FStr::blank(4);
                    c4.set(c13.sub(4, n));
                    c4.set(adjustr(c4));
                    int m = 0;
                    for (int i = 1; i <= 4; ++i) {
                        int j = 0;
                        const char c = c4(i);
                        if (c >= 'A' && c <= 'Z') j = c - 'A' + 1;
                        m = 27 * m + j;
                    }
                    n28 = 3 + 1000 + m;
                }
                done = true;
            }
        }
        if (!done) {
            if (c13(1) == '<') {
                if (!pack77_valid_hash_call_token(c13)) {
                    n28 = PACK77_NTOKENS;
                } else {
                    stage_hash_call(hash_facts, c13);
                    const int i2 = index(c13, '>');
                    c13.set(FStr(c13).sub(2, i2 - 1));
                    n28 = PACK77_NTOKENS + ihashcall_f(c13, 22);
                }
            } else if (!pack77_c28_standard_shape(c13)) {
                // A nonstandard call: its 22-bit hash.
                stage_hash_call(hash_facts, c13);
                n28 = PACK77_NTOKENS + ihashcall_f(c13, 22);
            } else {
                const int n = len_trim(c13);
                int i;
                for (i = n; i >= 2; --i)
                    if (pack77_is_digit_char(c13(i))) break;
                const int iarea = i;
                stage_hash_call(hash_facts, c13);
                FStr callsign = FStr::blank(6);
                if (iarea == 2) callsign.set(" " + c13.sub(1, 5));
                if (iarea == 3) callsign.set(c13.sub(1, 6));
                const int i1 = index(a1, callsign(1)) - 1;
                const int i2 = index(a2, callsign(2)) - 1;
                const int i3 = index(a3, callsign(3)) - 1;
                const int i4 = index(a4, callsign(4)) - 1;
                const int i5 = index(a4, callsign(5)) - 1;
                const int i6 = index(a4, callsign(6)) - 1;
                if (i1 < 0 || i2 < 0 || i3 < 0 || i4 < 0 || i5 < 0 || i6 < 0) {
                    n28 = PACK77_NTOKENS + ihashcall_f(c13, 22);
                } else {
                    n28 = 36 * 10 * 27 * 27 * 27 * i1 + 10 * 27 * 27 * 27 * i2 + 27 * 27 * 27 * i3 + 27 * 27 * i4 +
                          27 * i5 + i6;
                    n28 = n28 + PACK77_NTOKENS + PACK77_MAX22;
                }
            }
        }
    }
    n28 = n28 & ((1 << 28) - 1);
}

// lower_c28_slot: special tokens, explicit <CALL> hashes and standard calls;
// other calls hashed only where the field allows it.
void lower_c28_slot(const FStr& token, bool allow_auto_hash, int& n28, HashFacts& hash_facts, bool& ok) {
    ok = true;
    n28 = -1;
    FStr cw = FStr::blank(13);
    cw.set(token);
    if (cw.sub(1, 3) == "DE " || cw.sub(1, 4) == "QRZ " || cw.sub(1, 3) == "CQ " || cw.sub(1, 3) == "CQ_") {
        HashFacts trial_facts = hash_facts;
        pack28_core(cw, n28, trial_facts);
        if (n28 < PACK77_NTOKENS) return;
    }
    cw.set(token);
    if (pack77_valid_hash_call_token(cw) || pack77_c28_standard_shape(cw)) {
        pack28_core(cw, n28, hash_facts);
        return;
    }
    if (!allow_auto_hash) {
        ok = false;
        return;
    }
    pack28_core(cw, n28, hash_facts);
}

int pack77_c28_render_len(const FStr& token) {
    const int m = len_trim(token);
    int n = m;
    if (pack77_valid_hash_call_token(token)) return n;
    if (!pack77_c28_standard_shape(token)) n = m > 11 ? -1 : m + 2;
    return n;
}

bool pack77_dxpedition_render_fits(const DxpeditionSource& source) {
    int n = pack77_c28_render_len(source.call_1);
    if (n < 0) return false;
    const int ncall = pack77_c28_render_len(source.call_2);
    if (ncall < 0) return false;
    n = n + 1 + 5 + 1 + ncall + 1 + len_trim(source.hash_token) + 1 + 3;
    return n <= 37;
}

bool pack77_field_day_render_fits(const FieldDaySource& source) {
    const FStr section = pack77_arrl_section_name(source.isec);
    int nexch = 2;
    if (source.ntx >= 10) nexch = 3;
    int n = pack77_c28_render_len(source.call_1);
    if (n < 0) return false;
    const int ncall = pack77_c28_render_len(source.call_2);
    if (ncall < 0) return false;
    n = n + 1 + ncall + 1 + nexch + 1 + len_trim(section);
    if (source.ir == 1) n = n + 2;
    return n <= 37;
}

void pack77_accept_candidate_encoded(Candidate& candidate, int i3, int n3, const C77& c77, bool ok) {
    copy_result(candidate.result, pack77_no_match());
    if (ok) copy_result(candidate.result, pack77_accept(i3, n3, c77));
}

PackResult pack77_05(const FStr& msg) {
    FStr text = FStr::blank(37);
    text.set(adjustl(msg));
    const int n = len_trim(text);
    if (n < 1 || n > 18) return pack77_no_match();
    for (int i = 1; i <= n; ++i) {
        char ch = text(i);
        if (ch >= 'a' && ch <= 'z') ch = static_cast<char>(ch - 32);
        if (index(FStr("0123456789ABCDEF"), ch) == 0) return pack77_no_match();
        text.at(i) = ch;
    }
    FStr c18 = FStr::blank(18);
    c18.set(text.sub(1, n));
    c18.set(adjustr(c18));
    // read(c18,'(3z6)') ntel: blanks in a field are ignored (zero when all
    // blank); only hex digits are left, so it cannot fail.
    int ntel[3];
    for (int k = 0; k < 3; ++k) {
        int v = 0;
        for (int i = 1; i <= 6; ++i) {
            const char ch = c18(6 * k + i);
            if (ch == ' ') continue;
            v = 16 * v + (ch <= '9' ? ch - '0' : ch - 'A' + 10);
        }
        ntel[k] = v;
    }
    C77 c77 = FStr::blank(77);
    if (ntel[0] >= 0 && ntel[1] >= 0 && ntel[2] >= 0) {
        if (encode_pack77_telemetry(ntel[0], ntel[1], ntel[2], c77)) return pack77_accept(0, 5, c77);
    }
    return pack77_no_match();
}

PackResult pack77_free_text_accept(const FStr& msg, int status) {
    FStr msg0 = FStr::blank(37);
    msg0.set(msg);
    msg0.put(14, 37, FStr("                        "));
    const FStr text_bits = packtext77(msg0.sub(1, 13));
    C77 c77 = FStr::blank(77);
    if (!encode_pack77_free_text(text_bits, c77)) return pack77_no_match();
    PackResult r = pack77_accept(0, 0, c77);
    r.status = status;
    return r;
}

PackResult pack77_free_text_exact(const FStr& msg) {
    const FStr msg0 = pack77_normalized_message(msg);
    if (!pack77_exact_free_text_ok(msg0)) return pack77_reject(pack77_free_text_reject_status(msg0));
    return pack77_free_text_accept(msg0, PACK77_STATUS_ENCODED);
}

}  // namespace

// ---------------------------------------------------------------- the module

struct Packjt77::Impl {
    static void hash10(const Packjt77& p, int n10, FStr& c13) {
        c13.set("<...>");
        if (n10 < 0 || n10 > 1023) return;
        if (len_trim(FStr(p.calls10_[n10])) > 0) {
            c13.set(FStr(p.calls10_[n10]));
            c13.set("<" + trim(c13) + ">");
        }
    }

    static void hash12(const Packjt77& p, int n12, FStr& c13) {
        c13.set("<...>");
        if (n12 < 0 || n12 > 4095) return;
        if (len_trim(FStr(p.calls12_[n12])) > 0) {
            c13.set(FStr(p.calls12_[n12]));
            c13.set("<" + trim(c13) + ">");
        }
    }

    static void hash22(const Packjt77& p, int n22, FStr& c13) {
        c13.set("<...>");
        for (int i = 1; i <= p.nzhash_; ++i) {
            if (p.ihash22_[i - 1] == n22) {
                c13.set(FStr(p.calls22_[i - 1]));
                c13.set("<" + trim(c13) + ">");
                return;
            }
        }
    }

    static void save_hash_call(Packjt77& p, const FStr& c13, int& n10, int& n12, int& n22) {
        FStr cw = FStr::blank(13);
        cw.set(c13);
        if (cw(1) == ' ' || cw.sub(1, 5) == "<...>") return;
        if (cw(1) == '<') cw.set(FStr(cw).from(2));
        const int i = index(cw, '>');
        if (i > 0) cw.put(i, 13, FStr("         "));
        if (len_trim(cw) < 3) return;

        const FStr mycall13(p.mycall13_);
        n10 = ihashcall_f(cw, 10);
        if (n10 >= 0 && n10 <= 1023 && cw != mycall13) p.calls10_[n10] = to13(cw);
        n12 = ihashcall_f(cw, 12);
        if (n12 >= 0 && n12 <= 4095 && cw != mycall13) p.calls12_[n12] = to13(cw);

        n22 = ihashcall_f(cw, 22);
        const Call13 c = to13(cw);
        // any(ihash22.eq.n22) over all MAXHASH entries (unused ones are -1):
        // the entry gets the latest spelling, in place.
        bool exists = false;
        for (int k = 0; k < MAXHASH; ++k)
            if (p.ihash22_[k] == n22) exists = true;
        if (exists) {
            for (int k = 0; k < MAXHASH; ++k)
                if (p.ihash22_[k] == n22) p.calls22_[k] = c;
            return;
        }
        // A new entry at the top, the table moved down.
        for (int k = MAXHASH - 1; k >= 1; --k) {
            p.ihash22_[k] = p.ihash22_[k - 1];
            p.calls22_[k] = p.calls22_[k - 1];
        }
        p.ihash22_[0] = n22;
        p.calls22_[0] = c;
        if (p.nzhash_ < MAXHASH) ++p.nzhash_;
    }

    static void add_call_to_recent_calls(Packjt77& p, const Call13& callsign) {
        // Added unless already on the list (a call last on it is added again).
        bool ladd = true;
        for (int i = 1; i <= MAXRECENT - 1; ++i)
            if (p.recent_calls_[i - 1] == callsign) ladd = false;
        if (ladd) {
            for (int i = MAXRECENT; i >= 2; --i) p.recent_calls_[i - 1] = p.recent_calls_[i - 2];
            p.recent_calls_[0] = callsign;
        }
    }

    static void commit_hash_facts(Packjt77& p, const HashFacts& hash_facts) {
        int n10, n12, n22;
        for (int i = 0; i < hash_facts.ncalls; ++i) save_hash_call(p, FStr(hash_facts.calls[i]), n10, n12, n22);
    }

    static void apply_unpack77_effects(Packjt77& p, const Effects& effects) {
        int n10, n12, n22;
        for (int i = 0; i < effects.nhash; ++i) save_hash_call(p, FStr(effects.hash_calls[i]), n10, n12, n22);
        for (int i = 0; i < effects.nrecent; ++i) add_call_to_recent_calls(p, effects.recent_calls[i]);
    }

    // unpack28 (unpack28_core, two-digit prefixes allowed).
    static bool unpack28(const Packjt77& p, int n28_0, FStr& c13) {
        static const FStr c1(" 0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ");
        static const FStr c2("0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ");
        static const FStr c3("0123456789");
        static const FStr c4(" ABCDEFGHIJKLMNOPQRSTUVWXYZ");
        bool success = true;
        int n28 = n28_0;
        bool done = false;
        if (n28 < PACK77_NTOKENS) {
            // Special tokens DE, QRZ, CQ, CQ_nnn, CQ_aaaa.
            if (n28 == 0) c13.set("DE           ");
            if (n28 == 1) c13.set("QRZ          ");
            if (n28 == 2) c13.set("CQ           ");
            if (n28 <= 2) {
                done = true;
            } else if (n28 <= 1002) {
                c13.set("CQ_" + fmt_i(n28 - 3, 3, 3));
                done = true;
            } else if (n28 <= 532443) {
                int n = n28 - 1003;
                const int i1 = n / (27 * 27 * 27);
                n = n - 27 * 27 * 27 * i1;
                const int i2 = n / (27 * 27);
                n = n - 27 * 27 * i2;
                const int i3 = n / 27;
                const int i4 = n - 27 * i3;
                c13.set(chr(c4(i1 + 1)) + chr(c4(i2 + 1)) + chr(c4(i3 + 1)) + chr(c4(i4 + 1)));
                c13.set(adjustl(c13));
                c13.set("CQ_" + FStr(c13).sub(1, 10));
                done = true;
            }
        }
        if (!done) {
            n28 = n28 - PACK77_NTOKENS;
            if (n28 < PACK77_MAX22) {
                // A 22-bit hash (negative when n28_0 is 532444..NTOKENS-1).
                hash22(p, n28, c13);
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
                c13.set(chr(c1(i1 + 1)) + chr(c2(i2 + 1)) + chr(c3(i3 + 1)) + chr(c4(i4 + 1)) + chr(c4(i5 + 1)) +
                        chr(c4(i6 + 1)));
                c13.set(adjustl(c13));
                if (!callok(trim(c13))) {
                    c13.set("QU1RK");
                    success = false;
                }
            }
        }
        const int i0 = index(c13, ' ');
        if (i0 != 0 && i0 < len_trim(c13)) {
            c13.set("QU1RK");
            success = false;
        }
        return success;
    }

    // packjt77.f90's callok (outside the module).
    static bool callok(const FStr& w) {
        auto islet = [](char c) {
            const int k = static_cast<unsigned char>(c);
            return (k >= 65 && k <= 90) || (k >= 97 && k <= 122);
        };
        auto isdig = [](char c) {
            const int k = static_cast<unsigned char>(c);
            return k >= 48 && k <= 57;
        };
        const int n = len_trim(w);
        if (n < 3) return false;
        if (w(1) == 'Q') return false;
        int i;
        for (i = n; i >= 1; --i)
            if (isdig(w(i))) break;
        const int i0 = i;
        if (i0 != 2 && i0 != 3) return false;
        FStr pfx = FStr::blank(2);
        pfx.set(w.sub(1, i0 - 1));
        int nlp = 0, ndp = 0;
        const int np = len_trim(pfx);
        for (int k = 1; k <= np; ++k) {
            if (isdig(pfx(k))) ++ndp;
            if (islet(pfx(k))) ++nlp;
        }
        if (nlp + ndp != np) return false;
        if (nlp == 0) return false;
        int nls = 0;
        const int ns = n - i0;
        if (ns < 1 || ns > 3) return false;
        for (int k = i0 + 1; k <= n; ++k)
            if (islet(w(k))) ++nls;
        if (nls < ns) return false;
        return true;
    }

    // ------------------------------------------------------------ unpack77

    static void decode_i3_0(const Packjt77& p, const C77& c77, const Context& context, CoreResult& decoded) {
        FStr call_1 = FStr::blank(13), call_2 = FStr::blank(13), call_3 = FStr::blank(13);
        FStr call_1a = FStr::blank(13);
        FStr crpt = FStr::blank(3), cntx = FStr::blank(3), cpfx = FStr::blank(3);
        FStr grid6 = FStr::blank(6), grid4 = FStr::blank(4);
        bool ok = true;
        const int n3 = decoded.n3;

        if (n3 == 0) {
            // 0.0 Free text. decode_pack77_free_text: n3 and i3 are 0.
            decoded.msg.put(1, 13, unpacktext77(c77.sub(1, 71)));
            decoded.msg.put(14, 37, FStr("                        "));
            decoded.msg.set(adjustl(decoded.msg));
            if (decoded.msg(1) == ' ') decoded.success = false;

        } else if (n3 == 1) {
            // 0.1 K1ABC RR73; W9XYZ <KH1/KH7Z> -11     DXpedition Mode
            const int n28a = static_cast<int>(field(c77, 1, 28));
            const int n28b = static_cast<int>(field(c77, 29, 28));
            const int n10 = static_cast<int>(field(c77, 57, 10));
            const int n5 = static_cast<int>(field(c77, 67, 5));
            const int irpt = 2 * n5 - 30;
            crpt.set(pack77_format_snr_report(irpt));
            bool unpk28_success = unpack28(p, n28a, call_1);
            if (!unpk28_success || n28a <= 2) decoded.success = false;
            unpk28_success = unpack28(p, n28b, call_2);
            if (!unpk28_success || n28b <= 2) decoded.success = false;
            hash10(p, n10, call_3);
            if (context.nrx == 1) {
                if (context.dxcall_set && len_trim(context.dxcall) >= 3 && context.hashdx10 == n10)
                    call_3.set("<" + trim(context.dxcall) + ">");
            }
            if (context.nrx == 0 && context.mycall_set && n10 == context.hashmy10)
                call_3.set("<" + trim(context.mycall) + ">");
            decoded.msg.set(trim(call_1) + " RR73; " + trim(call_2) + " " + trim(call_3) + " " + crpt);

        } else if (n3 == 2) {
            decoded.success = false;

        } else if (n3 == 3 || n3 == 4) {
            // 0.3/0.4 WA9XYZ KA1ABC R 16A EMA          ARRL Field Day
            const int n28a = static_cast<int>(field(c77, 1, 28));
            const int n28b = static_cast<int>(field(c77, 29, 28));
            const int ir = static_cast<int>(field(c77, 57, 1));
            const int intx = static_cast<int>(field(c77, 58, 4));
            const int nclass = static_cast<int>(field(c77, 62, 3));
            int isec = static_cast<int>(field(c77, 65, 7));
            if (isec > PACK77_NSEC || isec < 1) {
                decoded.success = false;
                isec = 1;
            }
            bool unpk28_success = unpack28(p, n28a, call_1);
            if (!unpk28_success || n28a <= 2) decoded.success = false;
            unpk28_success = unpack28(p, n28b, call_2);
            if (!unpk28_success || n28b <= 2) decoded.success = false;
            int ntx = intx + 1;
            if (n3 == 4) ntx = ntx + 16;
            cntx.put(1, 2, fmt_i(ntx, 2, 1));
            cntx.at(3) = static_cast<char>('A' + nclass);
            decoded.msg.set(pack77_format_exchange_message(
                call_1, call_2, ir, trim(adjustl(cntx)) + " " + pack77_arrl_section_name(isec), true));

        } else if (n3 == 5) {
            // 0.5 Telemetry (18 hex digits): write(msg,'(3z6.6)') ntel.
            char b[32];
            std::snprintf(b, sizeof b, "%06X%06X%06X", static_cast<unsigned>(field(c77, 1, 23)),
                          static_cast<unsigned>(field(c77, 24, 24)), static_cast<unsigned>(field(c77, 48, 24)));
            decoded.msg.set(FStr(b));

        } else if (n3 == 6) {
            switch (pack77_wspr_payload_type(c77)) {
                case 2: {
                    // WSPR Type 2. decode_pack77_wspr_type2's gates: idbm
                    // <= 18, j50 = 1, pad 51-71 = 0.
                    int n28 = 0, npfx = 0, idbm = 0;
                    get_int(ok, c77, 1, 28, 0, M28, n28);
                    get_int(ok, c77, 29, 16, 0, 65535, npfx);
                    get_int(ok, c77, 45, 5, 0, 18, idbm);
                    require_int(ok, c77, 50, 1, 1);
                    require_int(ok, c77, 51, 21, 0);
                    if (!ok) {
                        decoded.success = false;
                        return;
                    }
                    idbm = nint(static_cast<float>(idbm) * 10.0f / 3.0f);
                    if (!unpack28(p, n28, call_1)) decoded.success = false;
                    crpt.set(fmt_i(idbm, 3, 1));
                    bool is_prefix;
                    if (!pack77_wspr_affix_text(npfx, cpfx, is_prefix)) {
                        decoded.success = false;
                        return;
                    }
                    if (is_prefix) call_1a.set(trim(adjustl(cpfx)) + "/" + trim(call_1));
                    else call_1a.set(trim(call_1) + "/" + trim(adjustl(cpfx)));
                    decoded.msg.set(trim(call_1a) + " " + trim(adjustl(crpt)));
                    stage_unpack_hash_call(decoded.effects, call_1a);
                    break;
                }
                case 1: {
                    // WSPR Type 1: idbm <= 18, j49 = j50 = 0, pad = 0.
                    int n28 = 0, igrid4 = 0, idbm = 0;
                    get_int(ok, c77, 1, 28, 0, M28, n28);
                    get_int(ok, c77, 29, 15, 0, 32767, igrid4);
                    get_int(ok, c77, 44, 5, 0, 18, idbm);
                    require_int(ok, c77, 49, 1, 0);
                    require_int(ok, c77, 50, 1, 0);
                    require_int(ok, c77, 51, 21, 0);
                    if (!ok) {
                        decoded.success = false;
                        return;
                    }
                    idbm = nint(static_cast<float>(idbm) * 10.0f / 3.0f);
                    if (!unpack28(p, n28, call_1)) decoded.success = false;
                    // grid4 is undefined in the reference when to_grid4
                    // fails (and the message is then not used).
                    if (!to_grid4(igrid4, grid4)) decoded.success = false;
                    crpt.set(fmt_i(idbm, 3, 1));
                    decoded.msg.set(trim(call_1) + " " + grid4 + " " + trim(adjustl(crpt)));
                    stage_unpack_hash_call(decoded.effects, call_1);
                    break;
                }
                case 3: {
                    // WSPR Type 3: selector 48-50 = 2, pad = 0.
                    int n22 = 0, igrid6 = 0;
                    get_int(ok, c77, 1, 22, 0, 4194303, n22);
                    get_int(ok, c77, 23, 25, 0, 33554431, igrid6);
                    require_int(ok, c77, 48, 3, 2);
                    require_int(ok, c77, 51, 21, 0);
                    if (!ok) {
                        decoded.success = false;
                        return;
                    }
                    const int n28 = n22 + PACK77_NTOKENS;
                    if (!unpack28(p, n28, call_1)) decoded.success = false;
                    if (!to_grid(igrid6, grid6)) decoded.success = false;
                    decoded.msg.set(trim(call_1) + " " + grid6);
                    break;
                }
                default:
                    decoded.success = false;
                    return;
            }

        } else if (n3 > 6) {
            decoded.success = false;
        }
    }

    static void decode_type12(const Packjt77& p, const C77& c77, const Context& context, CoreResult& decoded) {
        FStr call_1 = FStr::blank(13), call_2 = FStr::blank(13);
        FStr crpt = FStr::blank(3), grid4 = FStr::blank(4);
        // Type 1 (standard message) or Type 2 ("/P"): every field is in
        // range, the i3 field is the tag already read.
        const int n28a = static_cast<int>(field(c77, 1, 28));
        const int ipa = static_cast<int>(field(c77, 29, 1));
        const int n28b = static_cast<int>(field(c77, 30, 28));
        const int ipb = static_cast<int>(field(c77, 58, 1));
        const int ir = static_cast<int>(field(c77, 59, 1));
        const int igrid4 = static_cast<int>(field(c77, 60, 15));
        bool unpk28_success = unpack28(p, n28a, call_1);
        if (context.nrx == 1 && context.mycall_set && context.hashmy22 == (n28a - PACK77_NTOKENS)) {
            call_1.set("<" + trim(context.mycall) + ">");
            unpk28_success = true;
        }
        if (!unpk28_success) decoded.success = false;
        unpk28_success = unpack28(p, n28b, call_2);
        if (!unpk28_success) decoded.success = false;
        if (call_1.sub(1, 3) == "CQ_") call_1.at(3) = ' ';
        if (index(call_1, '<') <= 0) {
            const int i = index(call_1, ' ');
            if (i >= 4 && ipa == 1 && decoded.i3 == 1) call_1.put(i, i + 1, FStr("/R"));
            if (i >= 4 && ipa == 1 && decoded.i3 == 2) call_1.put(i, i + 1, FStr("/P"));
            // call_1 is hashed only by configured (MTft8) decodes.
            if (i >= 4) stage_unpack_recent_call(decoded.effects, call_1);
        }
        if (index(call_2, '<') <= 0) {
            const int i = index(call_2, ' ');
            if (i >= 4 && ipb == 1 && decoded.i3 == 1) call_2.put(i, i + 1, FStr("/R"));
            if (i >= 4 && ipb == 1 && decoded.i3 == 2) call_2.put(i, i + 1, FStr("/P"));
            if (i >= 4) {
                stage_unpack_recent_call(decoded.effects, call_2);
                stage_unpack_hash_call(decoded.effects, call_2);
            }
        }
        if (igrid4 <= PACK77_MAXGRID4) {
            if (!to_grid4(igrid4, grid4)) decoded.success = false;
            decoded.msg.set(pack77_format_exchange_message(call_1, call_2, ir, grid4, true));
            if (decoded.msg.sub(1, 3) == "CQ " && ir == 1) decoded.success = false;
        } else {
            const int irpt = igrid4 - PACK77_MAXGRID4;
            if (irpt >= 1 && irpt <= 4)
                decoded.msg.set(pack77_format_qso_tail(call_1, call_2, pack77_tail_from_type12_irpt(irpt)));
            if (irpt >= 5) {
                crpt.set(pack77_format_snr_report(pack77_snr_from_report_index(irpt)));
                decoded.msg.set(pack77_format_exchange_message(call_1, call_2, ir, crpt, false));
            }
            if (decoded.msg.sub(1, 3) == "CQ " && irpt >= 2) decoded.success = false;
        }
    }

    static void decode_type3(const Packjt77& p, const C77& c77, CoreResult& decoded) {
        FStr call_1 = FStr::blank(13), call_2 = FStr::blank(13);
        FStr exchange = FStr::blank(4), crpt = FStr::blank(3);
        // Type 3: ARRL RTTY Contest.
        const int itu = static_cast<int>(field(c77, 1, 1));
        const int n28a = static_cast<int>(field(c77, 2, 28));
        const int n28b = static_cast<int>(field(c77, 30, 28));
        const int ir = static_cast<int>(field(c77, 58, 1));
        const int irpt = static_cast<int>(field(c77, 59, 3));
        const int nexch = static_cast<int>(field(c77, 62, 13));
        crpt.set(pack77_format_rtty_report(irpt));
        if (!unpack28(p, n28a, call_1)) decoded.success = false;
        if (!unpack28(p, n28b, call_2)) decoded.success = false;
        if (nexch > 8000 && nexch - 8000 <= PACK77_NUSCAN) {
            exchange.set(pack77_rtty_multiplier_name(nexch - 8000));
            if (!pack77_type3_render_fits(itu, ir, call_1, call_2, exchange)) {
                decoded.success = false;
                return;
            }
            decoded.msg.set(pack77_format_type3_message(itu, ir, call_1, call_2, crpt, exchange));
        } else if (nexch >= 1 && nexch <= 7999) {
            exchange.set(fmt_i(nexch, 4, 4));
            if (!pack77_type3_render_fits(itu, ir, call_1, call_2, exchange)) {
                decoded.success = false;
                return;
            }
            decoded.msg.set(pack77_format_type3_message(itu, ir, call_1, call_2, crpt, exchange));
        } else {
            decoded.success = false;
        }
    }

    static void decode_type4(const Packjt77& p, const C77& c77, const Context& context, CoreResult& decoded) {
        FStr call_1 = FStr::blank(13), call_2 = FStr::blank(13), call_3 = FStr::blank(13);
        const int n12 = static_cast<int>(field(c77, 1, 12));
        const FStr c11 = pack77_type4_call_text_from_n58(field(c77, 13, 58));
        const int iflip = static_cast<int>(field(c77, 71, 1));
        const int nrpt = static_cast<int>(field(c77, 72, 2));
        const int icq = static_cast<int>(field(c77, 74, 1));
        hash12(p, n12, call_3);
        if (iflip == 0) {
            // 12-bit hash for the TO call.
            call_1.set(call_3);
            call_2.set(adjustl(c11) + "  ");
            stage_unpack_recent_call(decoded.effects, call_2);
            stage_unpack_hash_call(decoded.effects, call_2);
            if (context.nrx == 1 && context.dxcall_set && context.mycall_set && call_2 == context.dxcall &&
                n12 == context.hashmy12)
                call_1.set("<" + trim(context.mycall) + ">");
            if (context.nrx == 1 && context.mycall_set && index(call_1, "<...>") > 0 && n12 == context.hashmy12)
                call_1.set("<" + trim(context.mycall) + ">");
        } else {
            // 12-bit hash for the DE call (call_1 hashed only by configured decodes).
            call_1.set(adjustl(c11));
            call_2.set(call_3);
            stage_unpack_recent_call(decoded.effects, call_1);
            if (context.nrx == 0 && context.mycall_set && n12 == context.hashmy12)
                call_2.set("<" + trim(context.mycall) + ">");
        }
        if (icq == 0) decoded.msg.set(pack77_format_qso_tail(call_1, call_2, pack77_tail_from_type4_nrpt(nrpt)));
        else decoded.msg.set("CQ " + trim(call_2));
    }

    static void decode_type5(const Packjt77& p, const C77& c77, const Context& context, CoreResult& decoded) {
        FStr call_1 = FStr::blank(13), call_2 = FStr::blank(13);
        FStr cexch = FStr::blank(6), grid6 = FStr::blank(6);
        // Type 5 <PA3XYZ> <G4ABC/P> R 590003 IO91NP    EU VHF contest.
        // decode_pack77_type5's gates: iserial >= 1, igrid6 <= 18662399.
        bool ok = true;
        int n12 = 0, n22 = 0, ir = 0, irpt = 0, iserial = 0, igrid6 = 0;
        get_int(ok, c77, 1, 12, 0, 4095, n12);
        get_int(ok, c77, 13, 22, 0, 4194303, n22);
        get_int(ok, c77, 35, 1, 0, 1, ir);
        get_int(ok, c77, 36, 3, 0, 7, irpt);
        get_int(ok, c77, 39, 11, 1, 2047, iserial);
        get_int(ok, c77, 50, 25, 0, 18662399, igrid6);
        if (!ok) {
            decoded.success = false;
            return;
        }
        hash12(p, n12, call_1);
        // No mycall_set guard here, so before mycall13 is ever set the
        // hash is unpack77's zero static.
        if (n12 == context.hashmy12) call_1.set("<" + trim(context.mycall) + ">");
        hash22(p, n22, call_2);
        cexch.set(pack77_format_vhf_exchange(irpt, iserial));
        decoded.success = to_grid6(igrid6, grid6);
        decoded.msg.set(pack77_format_exchange_message(call_1, call_2, ir, cexch + " " + grid6, true));
    }

    static CoreResult unpack77_core(const Packjt77& p, const C77& c77, const Context& context) {
        CoreResult decoded;
        if (!decode_pack77_tag(c77, decoded.i3, decoded.n3)) {
            decoded.msg.set("failed unpack");
            return decoded;
        }
        decoded.success = true;
        if (decoded.i3 == 0) decode_i3_0(p, c77, context, decoded);
        else if (decoded.i3 == 1 || decoded.i3 == 2) decode_type12(p, c77, context, decoded);
        else if (decoded.i3 == 3) decode_type3(p, c77, decoded);
        else if (decoded.i3 == 4) decode_type4(p, c77, context, decoded);
        else if (decoded.i3 == 5) decode_type5(p, c77, context, decoded);
        else if (decoded.i3 >= 6) decoded.success = false;
        // CQ addressed to a hashed call is no message in any type.
        if (decoded.msg.sub(1, 4) == "CQ <") decoded.success = false;
        return decoded;
    }

    // ------------------------------------------------------------ pack77

    // take_candidate: an encoded candidate wins only if the neutral decoder
    // (no calls set, the live tables) gives it back.
    static bool take_candidate(const Packjt77& p, const Candidate& candidate, const FStr& msg, EncodeCore& encoded) {
        bool won = candidate.result.encoded;
        if (won) {
            const CoreResult decoded = unpack77_core(p, candidate.result.c77, Context());
            bool messages_match, resolved_render_fits;
            pack77_gate_message_checks(msg, decoded.msg, candidate.hash_facts.calls, candidate.hash_facts.ncalls,
                                       messages_match, resolved_render_fits);
            if (!(decoded.success && decoded.i3 == candidate.result.i3 && decoded.n3 == candidate.result.n3 &&
                  messages_match && resolved_render_fits)) {
                won = false;
                copy_result(encoded.result, pack77_reject(PACK77_STATUS_INTERNAL_ROUNDTRIP_REJECTED));
            }
        }
        if (!won) return false;
        copy_result(encoded.result, candidate.result);
        encoded.hash_facts = candidate.hash_facts;
        return true;
    }

    static void pack77_01_candidate(const FStr& msg, Candidate& candidate) {
        // Type 0.1 K1ABC RR73; W9XYZ <KH1/KH7Z> -11   28 28 10 5   DXpedition
        DxpeditionSource source;
        bool ok;
        pack77_parse_dxpedition_source(msg, source, ok);
        if (!ok) return;
        if (!pack77_dxpedition_render_fits(source)) return;
        int n28a, n28b, n10 = 0;
        lower_c28_slot(source.call_1, true, n28a, candidate.hash_facts, ok);
        if (!ok) return;
        lower_c28_slot(source.call_2, true, n28b, candidate.hash_facts, ok);
        if (!ok) return;
        pack77_hash_index(candidate.hash_facts, source.hash_token, 10, n10, ok);
        if (!ok) return;
        C77 c77 = FStr::blank(77);
        ok = encode_pack77_dxpedition(n28a, n28b, n10, source.n5, c77);
        pack77_accept_candidate_encoded(candidate, 0, 1, c77, ok);
    }

    static void pack77_03_candidate(const FStr& msg, Candidate& candidate) {
        // 0.3/0.4 WA9XYZ KA1ABC R 16A EMA   28 28 1 4 3 7   ARRL Field Day
        FieldDaySource source;
        bool ok;
        pack77_parse_field_day_source(msg, source, ok);
        if (!ok) return;
        if (!pack77_field_day_render_fits(source)) return;
        int n3 = 3;
        int intx = source.ntx - 1;
        if (intx >= 16) {
            n3 = 4;
            intx = source.ntx - 17;
        }
        int n28a, n28b;
        lower_c28_slot(source.call_1, true, n28a, candidate.hash_facts, ok);
        if (!ok) return;
        lower_c28_slot(source.call_2, true, n28b, candidate.hash_facts, ok);
        if (!ok) return;
        C77 c77 = FStr::blank(77);
        ok = encode_pack77_field_day(n3, n28a, n28b, source.ir, intx, source.nclass, source.isec, c77);
        pack77_accept_candidate_encoded(candidate, 0, n3, c77, ok);
    }

    static void pack77_06_candidate(const FStr& msg, Candidate& candidate) {
        WsprSource source;
        bool ok;
        pack77_parse_wspr_source(msg, source, ok);
        if (!ok) return;
        C77 c77 = FStr::blank(77);
        int n28;
        if (source.subtype == 1) {
            // WSPR Type 1
            lower_c28_slot(source.call_token, false, n28, candidate.hash_facts, ok);
            if (!ok) return;
            const int igrid4 = pack77_grid4_index(source.grid4);
            ok = encode_pack77_wspr_type1(n28, igrid4, source.idbm, c77);
            pack77_accept_candidate_encoded(candidate, 0, 6, c77, ok);
            return;
        }
        if (source.subtype == 2) {
            // WSPR Type 2
            lower_c28_slot(source.base_call + "       ", false, n28, candidate.hash_facts, ok);
            if (!ok) return;
            ok = encode_pack77_wspr_type2(n28, source.npfx, source.idbm, c77);
            pack77_accept_candidate_encoded(candidate, 0, 6, c77, ok);
            return;
        }
        // WSPR Type 3 only with pack77_options' prefer_wspr_50bit (FST4W's
        // transmit side): never for FT8.
    }

    static void pack77_1_candidate(const FStr& msg, Candidate& candidate) {
        Type12Source source;
        bool ok;
        pack77_parse_type12_source(msg, source, ok);
        if (!ok) return;
        int n28a, n28b;
        lower_c28_slot(source.call_1.c28_token, false, n28a, candidate.hash_facts, ok);
        if (!ok) return;
        lower_c28_slot(source.call_2.c28_token, false, n28b, candidate.hash_facts, ok);
        if (!ok) return;
        const int ipa = pack77_type12_call_has_suffix(source.call_1) ? 1 : 0;
        const int ipb = pack77_type12_call_has_suffix(source.call_2) ? 1 : 0;
        C77 c77 = FStr::blank(77);
        ok = encode_pack77_type12(source.i3, n28a, ipa, n28b, ipb, source.ir, source.igrid4, c77);
        pack77_accept_candidate_encoded(candidate, source.i3, 0, c77, ok);
    }

    static void pack77_3_candidate(const FStr& msg, Candidate& candidate) {
        // 3 TU; W9XYZ K1ABC R 579 MA    1 28 28 1 3 13   ARRL RTTY contest
        Type3Source source;
        bool ok;
        pack77_parse_type3_source(msg, source, ok);
        if (!ok) return;
        int n28a, n28b;
        lower_c28_slot(source.call_1, true, n28a, candidate.hash_facts, ok);
        if (!ok) return;
        lower_c28_slot(source.call_2, true, n28b, candidate.hash_facts, ok);
        if (!ok) return;
        FStr exchange = FStr::blank(4);
        exchange.set(" ");
        if (source.nexch > 0 && source.nexch < 8000) exchange.set(fmt_i(source.nexch, 4, 4));
        if (source.nexch > 8000) exchange.set(pack77_rtty_multiplier_name(source.nexch - 8000));
        if (!pack77_type3_display_fits(source.itu, source.ir, n28a, n28b, source.call_1, source.call_2, exchange))
            return;
        C77 c77 = FStr::blank(77);
        ok = encode_pack77_type3(source.itu, n28a, n28b, source.ir, source.irpt, source.nexch, c77);
        pack77_accept_candidate_encoded(candidate, 3, 0, c77, ok);
    }

    static void pack77_4_candidate(const FStr& msg, Candidate& candidate) {
        // <WA9XYZ> PJ4/KA1ABC RR73   12 58 1 2 1   one nonstandard call, one hashed
        Type4Source source;
        bool ok;
        pack77_parse_type4_source(msg, source, ok);
        if (!ok) return;
        FStr hash_token = FStr::blank(13);
        hash_token.set(source.hash_token);
        if (source.icq == 1 && len_trim(hash_token) == 0) hash_token.set(source.c11);
        int n12 = 0;
        pack77_hash_index(candidate.hash_facts, hash_token, 12, n12, ok);
        if (!ok) return;
        const FStr c11 = adjustr(source.c11);
        const int64_t n58 = pack77_type4_n58_from_call_text(c11);
        int nrpt = pack77_type4_nrpt_from_tail(source.itail);
        if (source.icq == 1) nrpt = 0;
        C77 c77 = FStr::blank(77);
        ok = encode_pack77_type4(n12, n58, source.iflip, nrpt, source.icq, c77);
        pack77_accept_candidate_encoded(candidate, 4, 0, c77, ok);
    }

    static void pack77_5_candidate(const FStr& msg, Candidate& candidate) {
        // <PA3XYZ> <G4ABC/P> R 590003 IO91NP   h12 h22 r1 s3 s11 g25   EU VHF
        Type5Source source;
        bool ok;
        pack77_parse_type5_source(msg, source, ok);
        if (!ok) return;
        int n12 = 0, n22 = 0;
        pack77_hash_index(candidate.hash_facts, source.hash_token_1, 12, n12, ok);
        if (!ok) return;
        pack77_hash_index(candidate.hash_facts, source.hash_token_2, 22, n22, ok);
        if (!ok) return;
        const int igrid6 = pack77_grid6_index(source.grid6);
        C77 c77 = FStr::blank(77);
        ok = encode_pack77_type5(n12, n22, source.ir, source.irpt, source.iserial, igrid6, c77);
        pack77_accept_candidate_encoded(candidate, 5, 0, c77, ok);
    }

    // pack77_encode_core (no options: no DX macro, no WSPR 50-bit preference).
    static void pack77_encode_core(const Packjt77& p, const FStr& msg0, EncodeCore& encoded) {
        FStr msg = FStr::blank(37);
        msg.set(msg0);
        // Upper case, blanks collapsed, words (19 at most, 13 characters each).
        FStr source_tokens[19] = {
            FStr::blank(13), FStr::blank(13), FStr::blank(13), FStr::blank(13), FStr::blank(13),
            FStr::blank(13), FStr::blank(13), FStr::blank(13), FStr::blank(13), FStr::blank(13),
            FStr::blank(13), FStr::blank(13), FStr::blank(13), FStr::blank(13), FStr::blank(13),
            FStr::blank(13), FStr::blank(13), FStr::blank(13), FStr::blank(13)};
        int nwords;
        bool split_ok;
        pack77_split_source_tokens(msg, 19, nwords, source_tokens, 19, split_ok);
        const bool word_overflow = !split_ok;
        const bool starts_with_special_token =
            msg.sub(1, 3) == "CQ " || msg.sub(1, 3) == "DE " || msg.sub(1, 4) == "QRZ ";
        if (word_overflow) {
            if (!starts_with_special_token) {
                copy_result(encoded.result, pack77_05(msg));
                if (encoded.result.encoded) return;
            }
            copy_result(encoded.result, pack77_free_text_exact(msg));
            return;
        }
        if (!starts_with_special_token) {
            {
                Candidate candidate;  // 0.1 DXpedition
                pack77_01_candidate(msg, candidate);
                if (take_candidate(p, candidate, msg, encoded)) return;
            }
            {
                Candidate candidate;  // 0.3 and 0.4 ARRL Field Day
                pack77_03_candidate(msg, candidate);
                if (take_candidate(p, candidate, msg, encoded)) return;
            }
            if (nwords < 2) {
                // 0.5 telemetry
                copy_result(encoded.result, pack77_05(msg));
                if (encoded.result.encoded) return;
            }
        }
        {
            Candidate candidate;  // 0.6 WSPR
            pack77_06_candidate(msg, candidate);
            if (take_candidate(p, candidate, msg, encoded)) return;
        }
        {
            Candidate candidate;  // Type 1 or 2, with optional /R or /P
            pack77_1_candidate(msg, candidate);
            if (take_candidate(p, candidate, msg, encoded)) return;
        }
        {
            Candidate candidate;  // Type 3 ARRL RTTY contest
            pack77_3_candidate(msg, candidate);
            if (take_candidate(p, candidate, msg, encoded)) return;
        }
        {
            Candidate candidate;  // Type 4 one nonstandard call, one hashed
            pack77_4_candidate(msg, candidate);
            if (take_candidate(p, candidate, msg, encoded)) return;
        }
        {
            Candidate candidate;  // Type 5 EU VHF contest
            pack77_5_candidate(msg, candidate);
            if (take_candidate(p, candidate, msg, encoded)) return;
        }
        // Exactly representable free text, or nothing.
        copy_result(encoded.result, pack77_free_text_exact(msg));
    }

    static PackResult pack77_encode_result(Packjt77& p, const FStr& msg0, bool record_tx_hashes) {
        EncodeCore core;
        pack77_encode_core(p, msg0, core);
        PackResult encoded;
        copy_result(encoded, core.result);
        // (pack77_failed_preferred_wspr_shape needs prefer_wspr_50bit.)
        if (record_tx_hashes) commit_hash_facts(p, core.hash_facts);
        return encoded;
    }

    static PackResult pack77_legacy_truncating_result(Packjt77& p, const FStr& msg0) {
        PackResult encoded = pack77_encode_result(p, msg0, true);
        if (encoded.encoded) return encoded;
        if (encoded.status == PACK77_STATUS_PREFERRED_FAMILY_REJECTED) return encoded;
        if (encoded.status == PACK77_STATUS_INTERNAL_ROUNDTRIP_REJECTED) return encoded;
        FStr msg = FStr::blank(37);
        msg.set(msg0);
        msg.set(pack77_normalized_message(msg));
        if (len_trim(msg) < 1) return encoded;
        FStr payload = FStr::blank(13);
        payload.set(msg.sub(1, 13));
        if (!pack77_free_text_alphabet_ok(payload)) return encoded;
        const PackResult fallback = pack77_free_text_accept(payload, PACK77_STATUS_ENCODED);
        if (fallback.encoded) copy_result(encoded, fallback);
        return encoded;
    }

    static Packed packed(const PackResult& r) {
        Packed out;
        out.i3 = r.i3;
        out.n3 = r.n3;
        out.status = r.status;
        out.c77 = r.c77.str();
        return out;
    }
};

// ---------------------------------------------------------------- public

Packjt77::Packjt77() {
    Call13 b;
    b.fill(' ');
    calls10_.fill(b);
    calls12_.fill(b);
    calls22_.fill(b);
    ihash22_.fill(-1);
    recent_calls_.fill(b);
    mycall13_ = b;
    dxcall13_ = b;
    mycall13_0_ = b;
    dxcall13_0_ = b;
}

void Packjt77::set_calls(const std::string& mycall12, const std::string& hiscall12) {
    FStr my = FStr::blank(12), his = FStr::blank(12);
    my.set(FStr(mycall12.substr(0, 12)));
    his.set(FStr(hiscall12.substr(0, 12)));
    dxcall13_ = to13(his);
    mycall13_ = to13(my);
}

Packjt77::Unpacked Packjt77::unpack77(const char* c77in, int nrx) {
    // The entry's bookkeeping: mycall13's hashes saved (and the call
    // itself saved) when it changed; dxcall13's hashes computed.
    const FStr mycall13(mycall13_), dxcall13(dxcall13_);
    if (mycall13 != FStr(mycall13_0_)) {
        if (len_trim(mycall13) > 2) {
            mycall13_set_ = true;
            mycall13_0_ = mycall13_;
            Impl::save_hash_call(*this, mycall13, hashmy10_, hashmy12_, hashmy22_);
        } else {
            mycall13_set_ = false;
        }
    }
    if (dxcall13 != FStr(dxcall13_0_)) {
        if (len_trim(dxcall13) > 2) {
            dxcall13_set_ = true;
            dxcall13_0_ = dxcall13_;
            hashdx10_ = ihashcall_f(dxcall13, 10);
            hashdx12_ = ihashcall_f(dxcall13, 12);
            hashdx22_ = ihashcall_f(dxcall13, 22);
        } else {
            dxcall13_set_ = false;
        }
    }

    Context context;
    context.nrx = nrx;
    context.mycall_set = mycall13_set_;
    context.dxcall_set = dxcall13_set_;
    context.mycall.set(FStr(mycall13_));
    context.dxcall.set(FStr(dxcall13_));
    context.hashmy10 = hashmy10_;
    context.hashmy12 = hashmy12_;
    context.hashmy22 = hashmy22_;
    context.hashdx10 = hashdx10_;

    // run_unpack77: the effects applied only on success, hashes first.
    const CoreResult decoded = Impl::unpack77_core(*this, FStr(c77in, 77), context);
    if (decoded.success) Impl::apply_unpack77_effects(*this, decoded.effects);
    return {decoded.success, decoded.msg.str()};
}

Packjt77::Packed Packjt77::pack77(const std::string& msg37, bool record_tx_hashes) {
    FStr msg = FStr::blank(37);
    msg.set(FStr(msg37.substr(0, 37)));
    return Impl::packed(Impl::pack77_encode_result(*this, msg, record_tx_hashes));
}

Packjt77::Packed Packjt77::pack77_legacy_truncating_fallback(const std::string& msg37) {
    FStr msg = FStr::blank(37);
    msg.set(FStr(msg37.substr(0, 37)));
    return Impl::packed(Impl::pack77_legacy_truncating_result(*this, msg));
}

Packjt77::Genft8 Packjt77::genft8(const std::string& msg37) {
    Genft8 g;
    const Packed packed = pack77_legacy_truncating_fallback(msg37);
    g.i3 = packed.i3;
    g.n3 = packed.n3;
    bool unpk77_success = false;
    if (g.i3 >= 0) {
        const Unpacked u = unpack77(packed.c77.data(), 0);
        unpk77_success = u.ok;
        g.msgsent = u.msg;
    }
    // read(c77,'(77i1)') msgbits: blanks read as 0.
    if (g.i3 >= 0 && unpk77_success) {
        for (int i = 0; i < 77; ++i) g.msgbits[i] = packed.c77[i] == '1';
        g.ok = true;
    } else {
        g.msgbits.fill(0);
        g.msgsent = "*** bad message ***                  ";
    }
    return g;
}

void Packjt77::save_hash_call(const std::string& c13, int& n10, int& n12, int& n22) {
    FStr c = FStr::blank(13);
    c.set(FStr(c13.substr(0, 13)));
    Impl::save_hash_call(*this, c, n10, n12, n22);
}

int Packjt77::ihashcall(const std::string& c13, int m) {
    FStr c = FStr::blank(13);
    c.set(FStr(c13.substr(0, 13)));
    return ihashcall_f(c, m);
}

Split77 split77(const std::string& msg37) {
    FStr msg = FStr::blank(37);
    msg.set(FStr(msg37.substr(0, 37)));
    FStr w[19] = {FStr::blank(13), FStr::blank(13), FStr::blank(13), FStr::blank(13), FStr::blank(13),
                  FStr::blank(13), FStr::blank(13), FStr::blank(13), FStr::blank(13), FStr::blank(13),
                  FStr::blank(13), FStr::blank(13), FStr::blank(13), FStr::blank(13), FStr::blank(13),
                  FStr::blank(13), FStr::blank(13), FStr::blank(13), FStr::blank(13)};
    int nw[19] = {};
    int iz = len_trim(msg);
    int j = 0, k = 0, n = 0;
    char c0 = ' ';
    for (int i = 1; i <= iz; ++i) {
        if (msg(i) == '\0') msg.at(i) = ' ';
        char c = msg(i);
        if (c == ' ' && c0 == ' ') continue;  // leading and repeated blanks
        if (c != ' ' && c0 == ' ') {
            ++k;  // a new word
            n = 0;
        }
        ++j;
        ++n;
        if (c >= 'a' && c <= 'z') c = static_cast<char>(c - 32);
        msg.at(j) = c;
        if (k <= 19) {
            nw[k - 1] = n;
            if (n <= 13) w[k - 1].at(n) = c;
        }
        c0 = c;
    }
    iz = j;
    int nwords = k;
    if (nwords > 0) {
        for (int i = 1; i <= std::min(nwords, 19); ++i)
            if (nw[i - 1] <= 13) nw[i - 1] = len_trim(w[i - 1]);
        msg.put(iz + 1, 37, FStr("                                     "));
        if (nwords >= 3) {
            std::string bcall_1;
            const bool ok1 = chkcall(w[2].str(), bcall_1);
            if (ok1 && w[0].sub(1, 3) == "CQ " && pack77_valid_cq_modifier(w[1])) {
                w[0].set("CQ_" + w[1].sub(1, 10));  // "CQ " made "CQ_"
                for (int i = 2; i <= 12; ++i) w[i - 1].set(FStr(w[i]));
                nwords = nwords - 1;
                for (int i = 1; i <= nwords; ++i) nw[i - 1] = len_trim(w[i - 1]);
            }
        }
    }
    Split77 s;
    s.msg = msg.str();
    s.nwords = nwords;
    for (int i = 0; i < 19; ++i) {
        s.nw[i] = nw[i];
        s.w[i] = w[i].str();
    }
    return s;
}

bool chkcall(const std::string& w13, std::string& bc6) {
    FStr w = FStr::blank(13);
    w.set(FStr(w13.substr(0, 13)));
    FStr bc = FStr::blank(6);
    auto isdigit = [](char c) { return c >= '0' && c <= '9'; };
    auto isletter = [](char c) { return c >= 'A' && c <= 'Z'; };
    bool cok = false;
    bc.set(w.sub(1, 6));
    const int n1 = len_trim(w);
    [&] {
        if (n1 > 11) return;
        if (index(w, '.') >= 1) return;
        if (index(w, '+') >= 1) return;
        if (index(w, '-') >= 1) return;
        if (index(w, '?') >= 1) return;
        if (n1 > 6 && index(w, '/') <= 0) return;
        const int i0 = index(w, '/');
        if (std::max(i0 - 1, n1 - i0) > 6) return;  // the base call under 7 characters
        if (i0 >= 2 && i0 <= n1 - 1) {              // the base call of a compound call
            if (i0 - 1 <= n1 - i0) bc.set(w.sub(i0 + 1, n1) + "   ");
            if (i0 - 1 > n1 - i0) bc.set(w.sub(1, i0 - 1) + "   ");
        }
        const int nbc = len_trim(bc);
        if (nbc > 6) return;
        // One of the first two characters a letter; Q only for QU1RK.
        if (!isletter(bc(1)) && !isletter(bc(2))) return;
        if (bc(1) == 'Q' && bc.sub(1, 5) != "QU1RK") return;
        // A digit second or third.
        int i1 = 0;
        if (isdigit(bc(2))) i1 = 2;
        if (isdigit(bc(3))) i1 = 3;
        if (i1 == 0) return;
        // A suffix of 1-3 letters.
        if (i1 == nbc) return;
        int n = 0;
        for (int i = i1 + 1; i <= nbc; ++i) {
            const char c = bc(i);
            if (c < 'A' || c > 'Z') return;
            ++n;
        }
        if (n >= 1 && n <= 3) cok = true;
    }();
    bc6 = bc.str();
    return cok;
}

bool stdcall(const std::string& call12) {
    FStr callsign = FStr::blank(12);
    callsign.set(FStr(call12.substr(0, 12)));
    auto is_digit = [](char c) { return c >= '0' && c <= '9'; };
    auto is_letter = [](char c) { return c >= 'A' && c <= 'Z'; };
    const int n = len_trim(callsign);
    // The DO variable after the loop: the exit's i, else 1 (or n when the
    // loop runs zero times).
    int i;
    if (n < 2) {
        i = n;
    } else {
        for (i = n; i >= 2; --i)
            if (is_digit(callsign(i))) break;
    }
    const int iarea = i;  // right-most digit (call area)
    int npdig = 0, nplet = 0;
    for (int k = 1; k <= iarea - 1; ++k) {
        if (is_digit(callsign(k))) ++npdig;
        if (is_letter(callsign(k))) ++nplet;
    }
    int nslet = 0;
    for (int k = iarea + 1; k <= n; ++k)
        if (is_letter(callsign(k))) ++nslet;
    return !(iarea < 2 || iarea > 3 || nplet == 0 || npdig >= iarea - 1 || nslet > 3);
}

void ft8apset(Packjt77& p, const std::string& mycall12_in, const std::string& hiscall12_in, int ncontest,
              std::array<int, 58>& apsym, std::array<int, 10>& aph10) {
    FStr mycall12 = FStr::blank(12), hiscall12 = FStr::blank(12);
    mycall12.set(FStr(mycall12_in.substr(0, 12)));
    hiscall12.set(FStr(hiscall12_in.substr(0, 12)));
    apsym.fill(0);
    apsym[0] = 99;
    apsym[29] = 99;
    aph10.fill(0);
    aph10[0] = 99;
    if (len_trim(mycall12) < 3) return;

    bool nohiscall = false;
    FStr hiscall = FStr::blank(12);
    hiscall.set(hiscall12);
    if (len_trim(hiscall) < 3) {
        hiscall.set("KA1ABC");  // a dummy hiscall
        nohiscall = true;
    } else {
        int n10 = 0, n12 = 0, n22 = 0;
        p.save_hash_call(hiscall.str(), n10, n12, n22);
        // write(c10,'(b10.10)') iand(n10,Z'3FF'); read(c10,'(10i1.1)') aph10
        const int v = n10 & 0x3FF;
        for (int i = 0; i < 10; ++i) aph10[i] = 2 * ((v >> (9 - i)) & 1) - 1;
    }

    // A dummy standard message: i3=1, 28 1 28 1 1 15.
    FStr msg = FStr::blank(37);
    msg.set(trim(mycall12) + " " + trim(hiscall) + " RRR");
    if (!stdcall(mycall12.str())) msg.set("<" + trim(mycall12) + "> " + trim(hiscall) + " RRR");
    const Packjt77::Packed packed = p.pack77(msg.str());
    const Packjt77::Unpacked msgchk = p.unpack77(packed.c77.data(), 1);
    if (ncontest == 7 && (packed.i3 != 1 || !msgchk.ok)) return;
    if (ncontest <= 5 && (packed.i3 != 1 || msg != FStr(msgchk.msg) || !msgchk.ok)) return;

    // read(c77,'(58i1)') apsym(1:58): blanks read as 0.
    for (int i = 0; i < 58; ++i) apsym[i] = 2 * (packed.c77[i] == '1') - 1;
    if (nohiscall) {
        apsym[29] = 99;
        aph10[0] = 99;
    }
}

}  // namespace packjt77
