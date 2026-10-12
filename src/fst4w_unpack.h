// FST4W's message unpacking and callsign hash table: what WSJT-X's jt9 does
// with a decoded 50-bit FST4W payload (lib/77bit/packjt77.f90 at 567ad29ce:
// unpack77 on the payload with fst4_decode's fixed tail, i3=0 n3=6, WSPR
// types 1, 2 and 3), and the session state that unpacking reads and writes:
// packjt77's 22-bit hash table (calls22, ihash22, nzhash) and the mycall /
// hiscall bookkeeping fst4_decode and unpack77 keep. One Unpacker is one jt9
// process's packjt77: hashes saved by one decode resolve type 3 messages in
// later ones, as they do across the periods of one long-running jt9.
//
// Only what FST4W reaches. packjt77's 10- and 12-bit tables and its recent
// calls are written by the reference too but read only by 77-bit message
// types, so they are not kept.
//
// packjt77.f90, packjt77_schema.f90, packjt77_grammar.f90: Copyright 2018-2025
// Joe Taylor, K1JT, Steven Franke, K9AN, and the WSJT-X developers; GNU GPL v3,
// and so is this.
#pragma once

#include <array>
#include <cstdint>
#include <string>

namespace fst4w {

// A callsign field as Fortran holds it: character*13, blank padded.
using Call13 = std::array<char, 13>;

class Unpacker {
public:
    enum { MAXHASH = 1000 };

    Unpacker();

    // fst4_decode's per-call setup with jt9's -c and -x (jt9 defaults:
    // K1ABC, W9XYZ, when neither is given): mycall13/dxcall13 set,
    // and when either call changed since the last time and mycall has 3+
    // characters, "MYCALL HISCALL RR73" packed and unpacked as fst4_decode
    // does for its AP bits, which saves both calls' hashes. Only standard
    // calls (pack77's c28 shape, passing callok) are supported for that
    // seeding: anything else throws std::invalid_argument.
    void set_calls(const std::string& mycall, const std::string& hiscall);

    // unpack77(c77, 1, msg, success) on payload[0..49] (each 0 or 1) with
    // fst4_decode's tail (c77(51:77) = '000000000000000000000110000').
    // msg is the reference's character*37, blank padded; it is meaningful
    // only when ok. On success the message's call is saved in the table
    // (types 1 and 2), as the reference's effects are.
    struct Result {
        bool ok = false;
        std::string msg;  // 37 characters
    };
    Result unpack(const uint8_t* payload);

    // The table as packjt77 holds it, for tests.
    int nzhash() const { return nzhash_; }
    int ihash22(int i) const { return ihash22_[i - 1]; }      // 1-based, as Fortran
    std::string calls22(int i) const;                        // trailing blanks trimmed

    // The state carried from period to period (the 22-bit table in its
    // order, and fst4_decode's and unpack77's saved calls) as text, one
    // item a line, so a table outlives the process: load(save()) leaves
    // an Unpacker as this one is. load() takes nothing from text it does
    // not wholly accept (a hash not its call's, more than MAXHASH, an
    // unknown line) and says false.
    std::string save() const;
    bool load(const std::string& text);

    // packjt77's ihashcall(c13, m): the call's m-bit hash.
    static int ihashcall(const Call13& c, int m);

private:
    void save_hash_call(const Call13& c13);
    Call13 hash22(int n22) const;
    bool unpack28(int n28, Call13& c13) const;

    std::array<Call13, MAXHASH> calls22_;
    std::array<int, MAXHASH> ihash22_;
    int nzhash_ = 0;
    Call13 mycall13_, dxcall13_;
    // fst4_decode's saved mycall0/hiscall0 (character*12) and unpack77's
    // mycall13_0.
    std::array<char, 12> mycall0_, hiscall0_;
    Call13 mycall13_0_;
};

}  // namespace fst4w
