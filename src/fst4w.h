// FST4W decoding in-process: what WSJT-X's jt9 does with a period of
// FST4W-120/300/900/1800 (lib/fst4_decode.f90's decode with iwspr=1 at
// 567ad29ce, and the code under it: fst4w_frontend.h, fst4w_ldpc.h,
// fst4w_unpack.h), giving the reference's results to the bit on x86-64:
// the same decodes, SNR, dt, frequency, sync, hard errors, dmin, spread
// (w50) and messages, and the same callsign hash table and fst4w_calls.txt
// list carried from period to period as one long-running jt9 carries them.
// fst4w_jt9.cpp runs it as jt9's file mode does (jt9's lines and files).
//
// What it rests on is listed in fst4w_frontend.h and fst4w_ldpc.h; besides
// those, dopspread's waveform (gen_fst4wave) has its tone table and ramps
// from glibc's libmvec cosf/sinf four at a time, as the reference's -O3
// computes them (but the table's first and last three entries, which gcc
// folded at compile time: their constants), and its pulse from glibc's erff.
//
// fst4_decode.f90, gen_fst4wave.f90, genfst4.f90: Steven Franke, K9AN, and
// Joe Taylor, K1JT, WSJT-X (GNU GPL v3), and so is this.
#pragma once

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace fst4w {

// fst4_decode's arguments that matter for FST4W, with jt9's file-mode
// values (jt9 -p P --fst4w -f 1500 -F 100 -d D, no -X, -c, -x).
struct Options {
    int ntrperiod = 120;      // 15, 30, 60, 120, 300, 900 or 1800
    int ndepth = 1;           // 1..3 (jt9's -d; iand(ndepth,3))
    int nfqso = 1500;         // -f
    int ntol = 100;           // -F, at most 100 for FST4W
    int nexp_decode = 3 * 256;
    float emedelay = 0.0f;
    std::string mycall = "K1ABC", hiscall = "W9XYZ";  // jt9's -c, -x defaults
    bool lprinthash22 = false;  // jt9 -Y
    bool plotspec = true;       // a 'plotspec' file: dopspread runs (w50, fort.52)
    // dopspread's spectrum (with plotspec): Fast (the default), jt9's mixing
    // low-passed and decimated by 120 before a 120 times smaller FFT with the
    // same bins: w50 and fmid within about 4e-5 Hz of jt9's (printed w50 the
    // same on 463 of 465 test decodes), F30 several times cheaper; Exact,
    // jt9's own (an FFT of twice the period's length), bit for bit; Off, no
    // spread at all, as jt9 without plotspec; Both, Exact's with Fast's beside
    // it (Result::w50_fast), for comparing.
    enum class Spread { Fast, Exact, Off, Both };
    Spread spread = Spread::Fast;
    bool decdata = true;        // a 'decdata' file: fst4_decodes.dat lines
};

// One decode, as fst4_decode hands it to its callback and writes it to
// fst4_decodes.dat.
struct Result {
    int nutc = 0;
    float sync = 0.0f;  // the callback's sync: the reference's smax1, never set (0)
    int nsnr = 0;
    float xdt = 0.0f, fsig = 0.0f;
    std::string msg;  // character*37
    int nap = 0;      // iaptype
    float qual = 0.0f;
    int ntrperiod = 0;
    float fmid = -999.0f, w50 = 0.0f;
    float fmid_fast = -999.0f, w50_fast = 0.0f;  // Spread::Both: the fast spread's
    // fst4_decodes.dat's fields (written when Options::decdata).
    int icand = 0, itry = 0, nsyncoh = 0, iaptype = 0, ijitter = 0, npct = 0, ntype = 0, keff = 0, nsync_qual = 0,
        nharderrors = 0, nhp = 0;
    float dmin = 0.0f, hd = 0.0f, csync = 0.0f, xsnr = 0.0f;
    // dopspread's spectrum as it writes it to fort.52 (f12.6 pairs), when
    // Options::plotspec.
    std::vector<std::pair<float, float>> spectrum;
};

class Decoder {
public:
    Decoder();
    ~Decoder();
    Decoder(const Decoder&) = delete;
    Decoder& operator=(const Decoder&) = delete;

    // fst4w_calls.txt as fst4_decode reads it on its first call (its lines;
    // call before the first decode, as jt9 has the file then), and the list
    // it holds now, for writing back when calls_changed() says so.
    void load_calls(const std::vector<std::string>& lines);
    std::vector<std::string> calls() const;

    // fst4_decode's decode on iwave (jt9's id2: at least ntrperiod*12000
    // samples; its -32768s become -32767, as the reference's blanker
    // leaves them). on_result is called for each decode in order, as the
    // reference's callback. calls_changed: the period added a call to the
    // list (the reference rewrites fst4w_calls.txt then).
    void decode(int16_t* iwave, int nutc, const Options& o, const std::function<void(const Result&)>& on_result,
                bool* calls_changed);
    // The same on float samples at the 16-bit scale (full scale 32768), as
    // from a 32-bit float .wav: no rounding to 16 bits, nor clipping.
    void decode(float* iwave, int nutc, const Options& o, const std::function<void(const Result&)>& on_result,
                bool* calls_changed);

    // The callsign hash table (fst4w_unpack.h's Unpacker::save/load): kept
    // across runs, it goes on resolving type 3 messages as one jt9 running
    // all along would. load_hashes() false: the text refused, the table as
    // it was.
    std::string save_hashes() const;
    bool load_hashes(const std::string& text);

    // A flag decode() looks at between candidates (and blanker passes):
    // when set, it ends early, its decodes so far given; stopped() then
    // says so. nullptr: none.
    void set_stop(const std::atomic<bool>* stop);
    bool stopped() const;

    // What decode() allocates for o at its most, bytes (iwave apart).
    static size_t memory_needed(const Options& o);

private:
    struct Impl;
    std::unique_ptr<Impl> d_;
};

}  // namespace fst4w
