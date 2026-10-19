// The built-in FT4 decoder: what one long-running WSJT-X jt9 does for one
// FT4 dial, from a period's 12 kHz samples to its stdout lines and
// decoded.txt (WSJT-X at 567ad29ce: lib/jt9a.f90's one call a request,
// decoder.f90's FT4 branch with its RMS gate, decoded.txt and counts,
// ft4_decode.f90 whole and decoder_callbacks.f90's ft4_decoded), bit for
// bit on x86-64 Linux. The front end (ft4_frontend.h), the LDPC code
// (ft8_ldpc.h) and packjt77 (packjt77.h) are the validated pieces it drives.
//
// One Decoder is one jt9 process: everything jt9 keeps from call to call
// lives in it (ft4_decode's SAVEd first-call state, apbits, mycall0 and
// hiscall0, ctwk2; packjt77's tables; the front end's SAVEd state and FFT
// plans). Unlike FT8, FT4 keeps no decodes from period to period. Not
// thread-safe; one per dial.
//
// ft4_decode's decodes(100), the call's messages for the duplicate check,
// has no bound on its count there (no bounds checks in 567ad29ce's Release
// build): a 101st decode in one call writes past it, over 12 unused bytes,
// then savg, sbase and cd, which getcandidates4 and the next DT segment
// rewrite, so those entries stop matching and a later pass prints the
// message again. jt9 -5 with a wide -H reaches it (the over fixtures:
// a quirk of the reference, kept). Here decodes, savg, sbase and cd share one
// buffer laid out as that stack frame (frame_), so every overflow is the
// reference's; the most one call can decode (3 passes of 200 candidates,
// 22,200 bytes) stays inside cd.
//
// jt9a.f90, decoder.f90, decoder_callbacks.f90, decode_completion.f90,
// ft4_decode.f90, jt9.f90, jt9_params_init.f90: Joe Taylor, K1JT, Steven
// Franke, K9AN, and the WSJT-X developers (GNU GPL v3), and so is this.
#pragma once

#include <atomic>
#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "ft4_frontend.h"
#include "ft8.h"  // ft8::Params, ft8::Decode, ft8::RunOutput, read_wav, nutc_of
#include "packjt77.h"

namespace ft4 {

// The params block fields FT4 reads, as FT8's (ndepth, nfa, nfb, nfqso,
// nqso_progress, nexp_decode, lapcqonly, nagain, mycall, hiscall; the rest
// unused), the GUI's by default: blank calls.
using Params = ft8::Params;

// One decode as ft4_decoded gets it, with its two lines (stdout's format
// 1001 with " + ", decoded.txt's 1002 ending " FT4").
using Decode = ft8::Decode;

class Decoder {
public:
    Decoder();
    ~Decoder();
    Decoder(const Decoder&) = delete;
    Decoder& operator=(const Decoder&) = delete;

    // One multimode_decoder_core call (nmode 5): id2 is the call's 180,000
    // samples (the RMS gate reads them all; ft4_decode the first NMAX).
    // gated: false when the RMS gate skipped the call (decoded.txt not
    // opened). ndecoded: the call's count (decoder.f90's ndecoded).
    struct Call {
        bool gated = false;
        int ndecoded = 0;
        std::vector<Decode> decodes;  // in the order printed
    };
    Call call(const int16_t* id2, int nutc, const Params& p);
    // The same on float samples in 16-bit units (a 16-bit sample s is s: a
    // float .wav file's x as fst4w::Wav reads it, x * 32768), so float copies
    // of 16-bit samples decode as they do. One difference, as FT8's and FST4W's: the
    // RMS gate passes any call whose samples are not all zero, since a float
    // file's level is its own (an SDR's unscaled output sits far below RMS 0.5).
    Call call(const float* id2, int nutc, const Params& p);

    // jt9a's period for FT4: one call on samples (n of them, at most
    // `take` used, zeros past them in the 180,000: jt9a's request has the
    // cycle's 90,000, jt9 -5 reads 72,576), then <DecodeFinished>.
    struct Period {
        Call c;
        int ndecoded = 0;  // <DecodeFinished>'s
        std::string out;   // stdout: the decode lines and <DecodeFinished>
    };
    Period period(const int16_t* samples, size_t n, int nutc, const Params& p, size_t take = 90000);
    // The same on float samples, as call()'s.
    Period period(const float* samples, size_t n, int nutc, const Params& p, size_t take = 90000);

    // decoded.txt as decoder.f90 leaves it (as ft8::Decoder's): absent
    // until a call passes the RMS gate, then the lines of the last call
    // that wrote any.
    bool decoded_exists() const { return decoded_exists_; }
    const std::string& decoded_text() const { return decoded_text_; }

    // Stops a decode in progress: looked at between candidates and passes;
    // the call then ends as if nothing more were found. For the caller's
    // watchdog only (a stopped period's output is not jt9's).
    void set_stop(bool v) { stop_.store(v, std::memory_order_relaxed); }

    const packjt77::Packjt77& packer() const { return pk_; }

    // The memory it holds now, bytes (itself, packjt77's tables, its
    // buffers and the front end's), FFTW's plans not counted. Between
    // periods only.
    size_t memory() const;

    // The oracle trace's lines ft4_decode writes (scripts/ft8-trace-patch.py
    // ft4: C, P, Z, F, U, D, N; tests/data/ft4/README.md), without the
    // newline. For the tests.
    void set_trace(std::function<void(const std::string&)> f) { trace_ = std::move(f); }
    // The front end's trace (ft4::Frontend::set_trace).
    void set_frontend_trace(std::function<void(const std::string&)> f) { fe_.set_trace(std::move(f)); }

private:
    // call() and period() for both sample types (S: int16_t or float).
    template <class S>
    Call call_t(const S* id2, int nutc, const Params& p);
    template <class S>
    Period period_t(const S* samples, size_t n, int nutc, const Params& p, size_t take);
    // ft4_decode's decode(), with the callback.
    template <class S>
    void decode(const S* iwave, const Params& p, int nutc, int ncontest, Call& c);
    // decoder_callbacks.f90's ft4_decoded.
    void callback(float sync, int snr, float dt, float freq, const std::string& decoded, int nap, float qual,
                  int nutc, Call& c);
    bool stopped() const { return stop_.load(std::memory_order_relaxed); }

    Frontend fe_;
    packjt77::Packjt77 pk_;

    // ft4_decode's SAVEd state
    bool first_ = true;
    cfloat ctwk2_[33][2 * NSS];  // ctwk2(:, -16:16)
    int apbits_[2 * ND] = {};    // SAVEd without DATA: zero until set
    // apmy_ru and aphis_fd (ncontest 3 and 4 only): not SAVEd there, so
    // after a call whose calls did not change they hold whatever the stack
    // had; here what the last set-up left (not reached with ncontest 0).
    int apmy_ru_[28] = {}, aphis_fd_[28] = {};
    std::string mycall0_, hiscall0_;
    std::vector<float> dd_;
    // its locals (on the stack there; only what is written before it is
    // read is ever used)
    std::vector<float> candidate_;
    std::vector<cfloat> cd2_, cb_;
    std::vector<float> syncs_ = std::vector<float>(128);  // one idf's sync4d values (at most 115)
    std::vector<float> bitmetrics_;
    // decodes(100), 12 unused bytes, savg(NH1), sbase(NH1) and cd(0:NN*NSS-1)
    // as they lie in 567ad29ce's Release frame (0x2980, 0x37f4, 0x3800,
    // 0x4a00, 0x5c00 from %rsp): see the decodes(100) note above.
    std::vector<float> frame_;

    // decoder.f90's
    bool decoded_exists_ = false;
    std::string decoded_text_;
    std::string decoded_pending_;  // this call's lines to decoded.txt
    int wrote_ = 0;

    std::function<void(const std::string&)> trace_;
    std::atomic<bool> stop_{false};
};

// The oracle's harness (scripts/ft8-oracle.f90) in its FT4 modes, as a
// function, for the tests and the standalone jt9 -5:
//   ft8-oracle -5 -d D WAV...                    the GUI's decoding (jt9a),
//                                                decoded-calls
//   ft8-oracle --file -5 [jt9 options] WAV...    jt9 -5's file mode
// and WSJT-X's own "jt9 -a DIR -5 [options] WAV..." (as --file). Paths are
// relative to cwd. Fills what the harness leaves: stdout, decoded.txt,
// decoded-calls.txt (GUI mode), and with want_trace the oracle trace (its
// W lines and the decoder's). Returns the exit status (2 with a message on
// err for bad usage or input).
using RunOutput = ft8::RunOutput;
int run_periods(const std::vector<std::string>& args, const std::string& cwd, RunOutput& r);

}  // namespace ft4
