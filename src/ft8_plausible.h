// Ours (jt9 -8 --nn only): message plausibility rules. A decode that breaks them is not a real
// transmission: no confirmed decode on the validation recordings breaks them, false decodes often do.
//
// plausible(m77, deep), on the 77 message bits (m77[0] first):
//   always: a standard message (i3 = 1, 2) with a locator past RR99 ("SA00"), a report code no encoder makes (> 135)
//   or a CQ with the R flag; an EU VHF message (i3 = 5) with a 6-character locator past RR99XX.
//   deep (for OSD and a-priori decodes, which can turn noise into a CRC-valid word): also every full callsign must be
//   well-formed and start with an amateur call block in use (the call fields of standard, DXpedition (i3 0, n3 1) and
//   RTTY Roundup (i3 3) messages, and the non-standard call of i3 = 4 messages).
#pragma once

#include <cstdint>
#include <string>

namespace ft8::plaus {

// Amateur call-sign blocks in use: first characters
// that are a block on their own, then two-character blocks.
inline const char* kBlocks1 = "FGIKMNRUW";
inline const char* kBlocks2 =
    "1A2A2D2E2I2J2M2O2Q2R2U2V2W3A3B3C3D3E3F3G3H3I3J3K3L3M3N3O3P3Q3R3S3T3U3V3W3X3Y3Z4A4B4C4D4E4F4G4H4I4J4K"
    "4L4M4O4P4Q4R4S4T4U4V4W4X4Y4Z5A5B5C5D5E5F5G5H5I5J5K5L5M5N5O5P5Q5R5S5T5U5V5W5X5Y5Z6A6B6C6D6E6F6G6H6I6J"
    "6K6L6M6N6O6P6Q6R6S6T6U6V6W6X6Y6Z7A7B7C7D7E7F7G7H7I7J7K7L7M7N7O7P7Q7R7S7T7U7V7W7X7Y7Z8A8B8C8D8E8F8G8H"
    "8I8J8K8L8M8N8O8P8Q8R8S8T8U8V8W8X8Y8Z9A9B9C9D9E9F9G9H9I9J9K9L9M9N9O9P9Q9R9S9T9U9V9W9X9Y9ZA2A3A4A5A6A7"
    "A8A9AAABACADAEAFAGAHAIAJAKALAMANAOAPAQARASATAUAVAWAXAYAZB0B1B2B3B4B5B6B7B8B9BABDBGBHBIBJBLBMBNBOBPBQ"
    "BSBTBUBVBWBXBYBZC2C3C4C5C6C7C8C9CACBCCCDCECFCGCHCICJCKCLCMCNCOCPCQCRCSCTCUCVCWCXCYCZD2D3D4D5D6D7D8D9"
    "DADBDCDDDEDFDGDHDIDJDKDLDMDNDODPDQDRDSDTDUDVDWDXDYDZE2E3E4E5E6E7EAEBECEDEEEFEGEHEIEJEKELEMENEOEPEQER"
    "ESETEUEVEWEXEYEZH2H3H4H5H6H7H8H9HAHBHCHDHEHFHGHHHIHJHKHLHNHOHPHQHRHSHTHUHVHWHXHYHZJ2J3J4J5J6J7J8JAJD"
    "JEJFJGJHJIJJJKJLJMJNJOJPJQJRJSJTJUJVJWJXJYL1L2L3L4L5L6L7L8L9LALBLCLDLELFLGLHLILJLKLLLMLNLOLPLQLRLSLT"
    "LULVLWLXLYLZOAOBOCODOEOFOGOHOIOJOKOLOMONOOOPOQOROSOTOUOVOWOXOYOZP2P3P4P5P6P7P8P9PAPBPCPDPEPFPGPHPIPJ"
    "PKPLPMPNPOPPPQPRPSPTPUPVPWPXPYPZS0S2S3S4S5S6S7S8S9SASBSCSDSESFSGSHSISJSKSLSMSNSOSPSQSRSSSTSUSVSWSXSY"
    "SZT2T3T4T5T6T7T8TATBTCTDTETFTGTHTITJTKTLTMTNTOTPTQTRTSTTTUTVTWTXTYTZV2V3V4V5V6V7V8V9VAVBVCVDVEVFVGVI"
    "VJVKVLVOVPVQVRVTVUVVVWVXVYVZXAXBXCXDXEXFXGXHXIXJXKXLXMXNXOXPXQXRXSXTXUXVXWXXXYXZY2Y3Y4Y5Y6Y7Y8Y9YAYB"
    "YCYDYEYFYGYHYIYJYKYLYMYNYOYPYQYRYSYTYUYVYWYXYYZ2Z3Z6Z8ZAZBZCZDZFZGZKZLZMZPZQZRZSZTZUZVZWZXZYZZ";

constexpr int64_t kNtokens = 2063592, kMax22 = 4194304, kMaxGrid4 = 32400;
constexpr int64_t kMaxGrid6 = 18LL * 18 * 10 * 10 * 24 * 24;

inline uint64_t bits(const int8_t* m, int a, int b) {  // m[a..b) as an unsigned integer, m[a] most significant
    uint64_t v = 0;
    for (int i = a; i < b; ++i) v = (v << 1) | static_cast<uint64_t>(m[i] & 1);
    return v;
}

inline bool in_blocks(const std::string& base) {
    if (base.empty()) return false;
    for (const char* p = kBlocks1; *p; ++p)
        if (*p == base[0]) return true;
    if (base.size() < 2) return false;
    for (const char* p = kBlocks2; p[0] && p[1]; p += 2)
        if (p[0] == base[0] && p[1] == base[1]) return true;
    return false;
}

inline bool call_ok(uint64_t n28) {  // a 28-bit call field: tokens and hashes pass; a full call must be in a block
    if (static_cast<int64_t>(n28) < kNtokens + kMax22) return true;
    uint64_t n = n28 - kNtokens - kMax22;
    const uint64_t i1 = n / (36ULL * 10 * 27 * 27 * 27);
    if (i1 > 36) return false;
    static const char c1[] = " 0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ", c2[] = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ",
                      c3[] = "0123456789", c4[] = " ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    n %= 36ULL * 10 * 27 * 27 * 27;
    const uint64_t i2 = n / (10 * 27 * 27 * 27);
    n %= 10 * 27 * 27 * 27;
    const uint64_t i3 = n / (27 * 27 * 27);
    n %= 27 * 27 * 27;
    const uint64_t i4 = n / (27 * 27), i5 = n / 27 % 27, i6 = n % 27;
    std::string s{c1[i1], c2[i2], c3[i3], c4[i4], c4[i5], c4[i6]};
    const size_t a = s.find_first_not_of(' '), b = s.find_last_not_of(' ');
    if (a == std::string::npos) return false;
    s = s.substr(a, b - a + 1);
    return s.find(' ') == std::string::npos && in_blocks(s);
}

inline bool c58_ok(uint64_t n58) {  // i3 = 4 non-standard call: <= 11 chars, <= 2 '/', longest part a real call
    static const char c38[] = " 0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ/";
    uint64_t lim = 1;
    for (int i = 0; i < 11; ++i) lim *= 38;
    if (n58 >= lim) return false;
    std::string s(11, ' ');
    for (int i = 10; i >= 0; --i) {
        s[static_cast<size_t>(i)] = c38[n58 % 38];
        n58 /= 38;
    }
    const size_t a = s.find_first_not_of(' '), b = s.find_last_not_of(' ');
    if (a == std::string::npos) return false;
    s = s.substr(a, b - a + 1);
    if (s.find(' ') != std::string::npos) return false;
    std::string parts[4], base;
    int np = 0;
    for (char c : s) {
        if (c == '/') {
            if (++np > 2) return false;
        } else {
            parts[np] += c;
        }
    }
    for (int i = 0; i <= np; ++i) {
        if (parts[i].empty()) return false;
        if (parts[i].size() > base.size()) base = parts[i];
    }
    bool digit = false, alpha = false;
    for (char c : base) {
        digit = digit || (c >= '0' && c <= '9');
        alpha = alpha || (c >= 'A' && c <= 'Z');
    }
    return digit && alpha && in_blocks(base);
}

inline bool plausible(const int8_t* m, bool deep) {
    const uint64_t i3 = bits(m, 74, 77);
    if (i3 == 5) return static_cast<int64_t>(bits(m, 49, 74)) < kMaxGrid6;
    if (!deep && i3 != 1 && i3 != 2) return true;
    if (i3 == 4) return c58_ok(bits(m, 12, 70));
    if (i3 == 3) return call_ok(bits(m, 1, 29)) && call_ok(bits(m, 29, 57));
    if (i3 == 0) {
        if (bits(m, 71, 74) == 1) return call_ok(bits(m, 0, 28)) && call_ok(bits(m, 28, 56));
        return true;
    }
    if (i3 != 1 && i3 != 2) return true;
    const uint64_t n28a = bits(m, 0, 28), n28b = bits(m, 29, 57), ir = bits(m, 58, 59);
    const int64_t igrid4 = static_cast<int64_t>(bits(m, 59, 74));
    if (igrid4 == kMaxGrid4) return false;
    if (igrid4 - kMaxGrid4 > 135) return false;
    if (n28a >= 2 && static_cast<int64_t>(n28a) < kNtokens && ir) return false;
    return !deep || (call_ok(n28a) && call_ok(n28b));
}

// Contest message formats (ARRL RTTY Roundup i3 3, ARRL Field Day i3 0 n3 3/4, EU VHF i3 5): for OSD and a-priori
// decodes outside contest mode, where nearly all of them are false.
inline bool contest_format(const int8_t* m) {
    const uint64_t i3 = bits(m, 74, 77), n3 = bits(m, 71, 74);
    return i3 == 3 || i3 == 5 || (i3 == 0 && (n3 == 3 || n3 == 4));
}

}  // namespace ft8::plaus
