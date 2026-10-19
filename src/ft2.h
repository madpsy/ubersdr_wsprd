// The built-in FT2 decoder: what one
// long-running jt9 does for one FT2 dial, from a period's 12 kHz samples to
// its stdout lines and decoded.txt, bit for bit on x86-64 Linux with the
// FT2 reference: WSJT-X improved 3.2.0's FT2-own logic on 567ad29ce's FT4
// (the hybrid of tests/data/ft2/README.md: FT2 is built on 567ad29ce's
// FT4; improved's FT4 changes are not taken).
//
// FT2 in improved 3.2.0 is not a decoder of its own: its decoder.f90 takes
// the period's first 23*1728 = 39,744 samples, copies each twice into
// 79,488 (FT2's 41.67 Hz tones and baud become FT4's 20.83), halves nfa,
// nfb and nfqso, and runs FT4's decoder on that (tperiod 3.75: syncmin 1.0
// in ft4_decode, NMAX 23*3456 in ft4_params); its callback ft2_decoded maps
// dt and frequency back (dt/2 - 0.25, 2*freq) and prints " * " and " FT2".
// Here: the RMS gate, decoded.txt and the counts as FT4's decoder.f90
// branch; the branch (ft2::halve, ft2::stretch); ft4_decode whole as
// src/ft4.cpp has it, at FT2's sizes and threshold (its own copy: FT4's
// validated code is not touched); ft2_decoded (ft2::decoded_lines); the
// front end FT4's where its routines serve unchanged, FT2's own copies where
// the sizes differ (ft2_frontend.h).
//
// One Decoder is one jt9 process, as ft4::Decoder (ft4.h): ft4_decode's
// SAVEd state, packjt77's tables, the front end's SAVEd state and plans.
// decodes(100) and the stack frame after it as ft4.h's note: the hybrid's
// Release object lays decodes, savg, sbase and cd out at the same offsets
// as 567ad29ce's (0x2980, 0x3800, 0x4a00, 0x5c00 from %rsp), so an overflow
// past 100 decodes in a call is the reference's (not reached by any FT2
// fixture: at most 50 decodes a call were found possible).
//
// decoder.f90, decoder_callbacks.f90, decode_completion.f90,
// ft4_decode.f90, jt9.f90: Joe Taylor, K1JT, Steven Franke, K9AN, and the
// WSJT-X developers (GNU GPL v3); FT2's branch, ft2_decoded and jt9 --ft2
// from WSJT-X improved 3.2.0 (Uwe Risse, DG2YCB, and the same); and so is
// this.
#pragma once

#include <atomic>
#include <cstdint>
#include <functional>
#include <string>
#include <vector>

#include "ft2_frontend.h"
#include "ft8.h"  // ft8::Params, ft8::Decode, ft8::RunOutput, read_wav, nutc_of
#include "packjt77.h"

namespace ft2 {

// The params block fields FT2 reads, FT4's (ft4.h), unhalved: nfa, nfb and
// nfqso as the GUI or jt9 fill them; the decoder halves them.
using Params = ft8::Params;

// One decode as ft2_decoded gets it: sync, snr, nap and qual as
// ft4_decode passes them; dt and freq FT2's (dt_ft2 = 0.5*dt - 0.25,
// freq_ft2 = 2.0*freq, as printed); its two lines (stdout's format 1001
// with " * ", decoded.txt's 1002 ending " FT2").
using Decode = ft8::Decode;

// improved's decoder.f90 FT2 branch: nfa_ft2 = max(0, nfa/2), nfb_ft2 =
// max(nfa_ft2+1, nfb/2), nfqso_ft2 = max(0, nfqso/2) (integer division,
// toward zero as Fortran's).
struct Halved {
    int nfa = 0, nfb = 0, nfqso = 0;
};
Halved halve(int nfa, int nfb, int nfqso);
// id2_ft4(2i-1) = id2_ft4(2i) = id2(i), i = 1..FT2_NMAX: id2 FT2_NMAX
// samples in, id2_ft4 NMAX out.
void stretch(const int16_t* id2, int16_t* id2_ft4);
void stretch(const float* id2, float* id2_ft4);

// ft2_decoded's two lines from what ft4_decode passes it (FT4's domain:
// dt and freq stretched) and the period's nutc: d's fields filled, dt and
// freq FT2's.
Decode decoded_lines(int nutc, float sync, int snr, float dt, float freq, const std::string& decoded, int nap,
                     float qual);

class Decoder {
public:
    Decoder();
    ~Decoder();
    Decoder(const Decoder&) = delete;
    Decoder& operator=(const Decoder&) = delete;

    // One multimode_decoder_core call (nmode 52): id2 is the call's 180,000
    // samples (the RMS gate reads them all; the branch the first
    // FT2_NMAX). gated: false when the RMS gate skipped the call
    // (decoded.txt not opened). ndecoded: the call's count.
    struct Call {
        bool gated = false;
        int ndecoded = 0;
        std::vector<Decode> decodes;  // in the order printed
    };
    Call call(const int16_t* id2, int nutc, const Params& p);
    // The same on float samples in 16-bit units (a 16-bit sample s is s: a
    // float .wav file's x as fst4w::Wav reads it, x * 32768), so float copies
    // of 16-bit samples decode as they do. One difference, as FT8's and FT4's:
    // the RMS gate passes any call whose samples are not all zero, since a
    // float file's level is its own (an SDR's unscaled output sits far below
    // RMS 0.5).
    Call call(const float* id2, int nutc, const Params& p);

    // A period: one call on samples (n of them, at most `take` used, zeros
    // past them in the 180,000: the GUI's request has the cycle's 45,000,
    // jt9 --ft2 reads 34,560), then <DecodeFinished>.
    struct Period {
        Call c;
        int ndecoded = 0;  // <DecodeFinished>'s
        std::string out;   // stdout: the decode lines and <DecodeFinished>
    };
    Period period(const int16_t* samples, size_t n, int nutc, const Params& p, size_t take = 45000);
    // The same on float samples, as call()'s.
    Period period(const float* samples, size_t n, int nutc, const Params& p, size_t take = 45000);

    // decoded.txt as decoder.f90 leaves it (as ft4::Decoder's).
    bool decoded_exists() const { return decoded_exists_; }
    const std::string& decoded_text() const { return decoded_text_; }

    // Stops a decode in progress (ft4::Decoder's): for a watchdog only.
    void set_stop(bool v) { stop_.store(v, std::memory_order_relaxed); }

    const packjt77::Packjt77& packer() const { return pk_; }

    // The memory it holds now, bytes, FFTW's plans not counted. Between
    // periods only.
    size_t memory() const;

    // The oracle trace (tests/data/ft2/README.md): the branch's S line and
    // ft4_decode's C, P, Z, F, U, D, N lines, without the newline.
    void set_trace(std::function<void(const std::string&)> f) { trace_ = std::move(f); }
    // The front-end trace (tests/data/ft2/frontend): the front end's lines
    // and, from the decoder, the lines the instrumented ft4_decode and the
    // hybrid's branch write (S, C, P, E, e, Y, G, g, I, J, H, O, R:
    // tests/data/ft4/frontend/README.md), all in the reference's order.
    void set_frontend_trace(std::function<void(const std::string&)> f);

private:
    // call() and period() for both sample types (S: int16_t or float).
    template <class S>
    Call call_t(const S* id2, int nutc, const Params& p);
    template <class S>
    Period period_t(const S* samples, size_t n, int nutc, const Params& p, size_t take);
    template <class S>
    void decode(const S* iwave, const Params& p, int nutc, int ncontest, Call& c);
    bool stopped() const { return stop_.load(std::memory_order_relaxed); }

    Frontend fe_;
    packjt77::Packjt77 pk_;

    // ft4_decode's SAVEd state (ft4.h)
    bool first_ = true;
    cfloat ctwk2_[33][2 * NSS];
    int apbits_[2 * ND] = {};
    int apmy_ru_[28] = {}, aphis_fd_[28] = {};
    std::string mycall0_, hiscall0_;
    std::vector<float> dd_;
    std::vector<int16_t> id2_ft4_;
    std::vector<float> id2_ft4f_;  // float calls' stretch: made on the first one
    std::vector<float> candidate_;
    std::vector<cfloat> cd2_, cb_;
    std::vector<float> syncs_ = std::vector<float>(128);
    std::vector<float> bitmetrics_;
    std::vector<float> frame_;  // decodes(100) .. cd, as ft4.h's

    bool decoded_exists_ = false;
    std::string decoded_text_;
    std::string decoded_pending_;
    int wrote_ = 0;

    std::function<void(const std::string&)> trace_;
    std::function<void(const std::string&)> fetrace_;
    long nbm_ = 0;  // get_ft4_bitmetrics calls in the process (the front-end trace's sampling)
    std::atomic<bool> stop_{false};
};

// The oracle's harness in its FT2 modes, as a function, for the tests and
// the standalone jt9 --ft2:
//   ft2-oracle -2 -d D WAV...                 the GUI's decoding (45,000
//                                             samples a period, nzhsym 23),
//                                             decoded-calls
//   ft2-oracle --file -2 [jt9 options] WAV... jt9 --ft2's file mode
//   ft2-jt9 -a DIR --ft2 [options] WAV...     the same (also "jt9 ...")
// Paths relative to cwd. Returns the exit status (2 with a message on err
// for bad usage or input). fetrace, when set, gets the front-end trace
// (the harness's W lines and Decoder::set_frontend_trace's).
using RunOutput = ft8::RunOutput;
int run_periods(const std::vector<std::string>& args, const std::string& cwd, RunOutput& r,
                const std::function<void(const std::string&)>& fetrace = {});

}  // namespace ft2
