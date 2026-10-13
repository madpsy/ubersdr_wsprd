// WSJT-X's 77-bit message packing (lib/77bit/packjt77.f90,
// packjt77_schema.f90 and packjt77_grammar.f90 at 567ad29ce) as FT8 reaches
// it: pack77 (and pack77_legacy_truncating_fallback, as genft8 uses it) and
// unpack77 for every i3/n3, the module's callsign tables (calls10, calls12,
// calls22 with ihash22, the recent calls), unpack77's SAVEd mycall/dxcall
// bookkeeping, and the routines FT8's driver calls beside them: split77,
// chkcall, stdcall, genft8's packing half and ft8apset. One Packjt77 is one
// jt9 process's packjt77 module: every pack and unpack reads and writes its
// tables in the order the reference does, so hashes saved by one call
// resolve "<...>" in later ones.
//
// Only the non-configured path (jt9's ordinary decode, not the MTft8
// "configured/var" one): unpack77_configured, the per-thread queues and
// pack77's options other than record_tx_hashes are not reached by FT8.
//
// Strings in and out are Fortran's fixed-length ones as std::string: a
// message is its character*37, a call its character*13, blank padded, and
// any byte may occur (NUL included), as in the reference.
//
// packjt77.f90, packjt77_schema.f90, packjt77_grammar.f90, genft8.f90,
// ft8apset.f90, chkcall.f90, q65_set_list.f90 (stdcall): Copyright 2018-2025
// Joe Taylor, K1JT, Steven Franke, K9AN, and the WSJT-X developers; GNU GPL
// v3, and so is this.
#pragma once

#include <array>
#include <cstdint>
#include <string>

namespace packjt77 {

// pack77's status values (packjt77.f90's PACK77_STATUS_*).
enum {
    PACK77_STATUS_NOT_ENCODED = 0,
    PACK77_STATUS_ENCODED = 1,
    PACK77_STATUS_FREE_TEXT_TOO_LONG = 3,
    PACK77_STATUS_FREE_TEXT_INVALID = 4,
    PACK77_STATUS_PREFERRED_FAMILY_REJECTED = 5,
    PACK77_STATUS_INTERNAL_ROUNDTRIP_REJECTED = 6,
};

// A callsign table entry: character*13, blank padded.
using Call13 = std::array<char, 13>;

class Packjt77 {
public:
    enum { MAXHASH = 1000, MAXRECENT = 10 };

    Packjt77();

    // The module's mycall13 and dxcall13, as ft8b sets them from its
    // character*12 mycall12 and hiscall12 ("dxcall13=hiscall12;
    // mycall13=mycall12"): the values are taken as character*12 (cut or
    // blank padded to 12), then padded to 13.
    void set_calls(const std::string& mycall12, const std::string& hiscall12);

    // unpack77(c77, nrx, msg, success): c77 is 77 characters (normally '0'
    // and '1'; anything else fails as the reference does). msg is the
    // character*37, meaningful when ok (the reference's msg after a failure
    // is not used by FT8 and is partly undefined there).
    struct Unpacked {
        bool ok = false;
        std::string msg;  // 37 characters
    };
    Unpacked unpack77(const char* c77, int nrx);

    // pack77(msg, i3, n3, c77, options, status) on a character*37 message
    // (cut or blank padded to 37). record_tx_hashes is pack77_options'
    // (true: the encoded calls are saved, as every FT8 caller has it).
    // c77 is 77 characters, blanks when nothing was encoded (i3 = n3 = -1).
    struct Packed {
        int i3 = -1, n3 = -1, status = PACK77_STATUS_NOT_ENCODED;
        std::string c77;  // 77 characters
    };
    Packed pack77(const std::string& msg37, bool record_tx_hashes = true);
    Packed pack77_legacy_truncating_fallback(const std::string& msg37);

    // genft8(msg, i3, n3, msgsent, msgbits, itone) up to the tones:
    // pack77_legacy_truncating_fallback, unpack77(c77, 0) when it encoded,
    // the 77 message bits. ok is the reference's "go to 2" (encoded and
    // unpacked); otherwise msgbits are 0 and msgsent "*** bad message ***".
    struct Genft8 {
        int i3 = -1, n3 = -1;
        bool ok = false;
        std::string msgsent;  // 37 characters
        std::array<uint8_t, 77> msgbits{};
    };
    Genft8 genft8(const std::string& msg37);

    // save_hash_call(c13, n10, n12, n22): n10/n12/n22 are left as passed in
    // when the call is not saved (blank, "<...>", under 3 characters).
    void save_hash_call(const std::string& c13, int& n10, int& n12, int& n22);

    // ihashcall(c0, m) on a character*13.
    static int ihashcall(const std::string& c13, int m);

    // The tables, 1-based or 0-based as the reference's: each entry its
    // full character*13.
    std::string calls10(int n10) const { return str(calls10_[n10]); }   // 0..1023
    std::string calls12(int n12) const { return str(calls12_[n12]); }   // 0..4095
    std::string calls22(int i) const { return str(calls22_[i - 1]); }   // 1..MAXHASH
    int ihash22(int i) const { return ihash22_[i - 1]; }                // 1..MAXHASH
    int nzhash() const { return nzhash_; }
    std::string recent_call(int i) const { return str(recent_calls_[i - 1]); }  // 1..MAXRECENT
    std::string mycall13() const { return str(mycall13_); }
    std::string dxcall13() const { return str(dxcall13_); }

    // Internal: the reference's routines on its own string type (packjt77.cpp).
    struct Impl;

private:
    static std::string str(const Call13& c) { return std::string(c.data(), 13); }

    std::array<Call13, 1024> calls10_;
    std::array<Call13, 4096> calls12_;
    std::array<Call13, MAXHASH> calls22_;
    std::array<int, MAXHASH> ihash22_;
    int nzhash_ = 0;
    std::array<Call13, MAXRECENT> recent_calls_;
    Call13 mycall13_, dxcall13_;
    // unpack77's SAVEd locals: the hashes have no initialiser (static, so
    // zero); mycall13_0 is '' and dxcall13_0 ' ' (both 13 blanks).
    int hashmy10_ = 0, hashmy12_ = 0, hashmy22_ = 0;
    int hashdx10_ = 0, hashdx12_ = 0, hashdx22_ = 0;
    bool mycall13_set_ = false, dxcall13_set_ = false;
    Call13 mycall13_0_, dxcall13_0_;
};

// split77(msg, nwords, nw, w): msg (character*37) as split77 leaves it
// (upper case, blanks collapsed), the words (character*13 each) and their
// lengths as it computes them.
struct Split77 {
    std::string msg;  // 37 characters
    int nwords = 0;
    std::array<int, 19> nw{};
    std::array<std::string, 19> w;  // 13 characters each
};
Split77 split77(const std::string& msg37);

// chkcall(w, bc, cok) on a character*13; bc is its character*6.
bool chkcall(const std::string& w13, std::string& bc6);

// stdcall(callsign, std) on a character*12 (lib/qra/q65/q65_set_list.f90).
bool stdcall(const std::string& call12);

// ft8apset(mycall12, hiscall12, ncontest, apsym, aph10) (lib/ft8/ft8apset.f90):
// the a-priori symbols for MyCall passes, with its packjt77 effects.
void ft8apset(Packjt77& p, const std::string& mycall12, const std::string& hiscall12, int ncontest,
              std::array<int, 58>& apsym, std::array<int, 10>& aph10);

}  // namespace packjt77
