// The built-in FT8 decoder: what one long-running
// WSJT-X jt9 does for one FT8 dial, from a period's 180,000 12 kHz samples
// to its stdout lines and decoded.txt (WSJT-X at 567ad29ce: lib/jt9a.f90's
// 41/47/50 sequence, decoder.f90's FT8 branch with its RMS gate, decoded.txt
// and counts, ft8_decode.f90 whole, ft8_a7.f90's ft8_a7_save, ft8_a8d.f90
// with ft8_a7.f90's getmsg, twkfreq1, smo121 and db, and
// decoder_callbacks.f90's ft8_decoded), bit for bit on x86-64 Linux. The
// front end (ft8_frontend.h), the LDPC code (ft8_ldpc.h) and packjt77
// (packjt77.h) are the validated pieces it drives.
//
// One Decoder is one jt9 process: everything jt9 keeps from call to call and
// period to period lives in it (ft8_decode's SAVEd dd, dd1, nutc0,
// ndec_early, the early decodes' tones, frequencies and times, lsubtracted
// and allmessages; ft8_a7's tables and jseq; decoder.f90's ltry_a8, ndec41
// and ndec47; packjt77's tables; the front end's SAVEd state and FFT plans).
// Not thread-safe; one per dial.
//
// jt9a.f90, decoder.f90, decoder_callbacks.f90, decode_completion.f90,
// ft8_decode.f90, ft8_a7.f90, ft8_a8d.f90, twkfreq1.f90, smo121.f90, db.f90,
// jt9.f90, jt9_params_init.f90, jt9_input_validation.f90: Joe Taylor,
// K1JT, Steven Franke, K9AN, and the WSJT-X developers (GNU GPL v3), and so
// is this.
#pragma once

#include <atomic>
#include <cstdint>
#include <functional>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "ft8_frontend.h"
#include "packjt77.h"

namespace ft8 {

// The params block fields FT8 reads (jt9com.f90's params_block), as the
// GUI fills them (blank calls and grids, HHMMSS nutc) by default.
struct Params {
    int ndepth = 2;
    int nfa = 100, nfb = 3200;
    int nfqso = 1500, nftx = 1500;
    int napwid = 50;
    int nqso_progress = 0;
    int nexp_decode = 0;  // ncontest = iand(nexp_decode, 7)
    bool lft8apon = true, lapcqonly = false, nagain = false;
    // Ours: the subtraction's smoothing computed at a low rate (FT8: 300 S/s, Frontend::set_fast_subtract; FT4 and
    // FT2: 857 S/s, ft4::Frontend::set_fast_subtract): the same decodes as WSJT-X's (scripts/ft8-fast-subtract-check.sh),
    // not to the bit, for less CPU.
    bool fast_subtract = false;
    // character*12, *12 and *6 (cut or blank padded; blank: none).
    std::string mycall, hiscall, hisgrid;
};

// One decode as ft8_decoded gets it, with its two lines.
struct Decode {
    float sync = 0;
    int snr = 0;
    float dt = 0, freq = 0;
    std::string msg37;  // as printed ('?' in column 37 when it has one)
    int nap = 0;
    float qual = 0;
    std::string line;          // stdout's format 1001, no newline
    std::string decoded_line;  // decoded.txt's format 1002, no newline
};

class Decoder {
public:
    Decoder();
    ~Decoder();
    Decoder(const Decoder&) = delete;
    Decoder& operator=(const Decoder&) = delete;

    // One multimode_decoder_core call (FT8, lmultift8 false): id2 is the
    // call's 180,000 samples (jt9a's: the first 41*3456 or 47*3456 and
    // zeros for the early calls). Returns the call's ndecoded as decoder.f90
    // counts it (at nzhsym 50 the period's: ndec41 + ndec47 + this call's).
    // gated: false when the RMS gate skipped the call (decoded.txt not
    // opened).
    struct Call {
        bool gated = false;
        int ndecoded = 0;
        std::vector<Decode> decodes;  // in the order printed
    };
    Call call(const int16_t* id2, int nzhsym, int nutc, const Params& p);
    // The same on float samples in 16-bit units (a 16-bit sample s is s: a
    // float .wav file's x as fst4w::Wav reads it, x * 32768), so float copies
    // of 16-bit samples decode as they do. One difference, as FST4W's: the
    // RMS gate passes any call whose samples are not all zero, since a float
    // file's level is its own (an SDR's unscaled output sits far below RMS 0.5).
    Call call(const float* id2, int nzhsym, int nutc, const Params& p);

    // jt9a's period: the calls at nzhsym 41, 47 and 50 on samples (n of them,
    // the first 180,000 used, zeros past n) and <DecodeFinished>'s count.
    // fifty: how many samples the 50 call sees (jt9a's 180,000; jt9 -8's
    // 50*3456). decoded.txt is followed in decoded_text (see below).
    // after_call(nzhsym), when given, runs after each call.
    struct Period {
        Call c41, c47, c50;
        int ndecoded = 0;  // <DecodeFinished>'s
        std::string out;   // stdout: the decode lines and <DecodeFinished>
    };
    Period period(const int16_t* samples, size_t n, int nutc, const Params& p, size_t fifty = 180000,
                  const std::function<void(int)>& after_call = {});
    // The same on float samples, as call()'s.
    Period period(const float* samples, size_t n, int nutc, const Params& p, size_t fifty = 180000,
                  const std::function<void(int)>& after_call = {});

    // decoded.txt as decoder.f90 leaves it in its temp dir: absent until a
    // call passes the RMS gate (opened with status 'unknown'), then the
    // lines of the last call that wrote any (gfortran truncates only after
    // a WRITE; a call that writes none leaves the file as it was).
    bool decoded_exists() const { return decoded_exists_; }
    const std::string& decoded_text() const { return decoded_text_; }

    // Called at the end of each ft8_decode call (after the RMS gate passed
    // and the call ran), with the decoder's packjt77: its state then (the
    // tests' D lines).
    void set_after_decode(std::function<void(int nutc, int nzhsym, const packjt77::Packjt77&)> f) {
        after_decode_ = std::move(f);
    }

    // Stops a decode in progress: looked at between candidates, passes and
    // a7 entries; the call then ends as if nothing more were found. For the
    // caller's watchdog only (a stopped period's output is not jt9's).
    void set_stop(bool v) { stop_.store(v, std::memory_order_relaxed); }

    const packjt77::Packjt77& packer() const { return pk_; }

    // The memory it holds now, bytes: itself (packjt77's tables, the a7
    // tables and the front end's fixed arrays) and its buffers on the heap,
    // plus the two sample copies period() makes while it runs; FFTW's plans
    // not counted. Between periods (not thread-safe).
    size_t memory() const;

    // The oracle trace's lines this decoder can tell (scripts/ft8-trace-patch.py:
    // C at each decode() call after the a7 rotation, P after each sync8, A
    // after each ft8_a7d), without the newline. For the tests.
    void set_trace(std::function<void(const std::string&)> f) { trace_ = std::move(f); }
    // The front end's trace (Frontend::set_trace: sync8, ft8b, subtractft8
    // and ft8_a7d as float bits and hashes) and its pending chain written
    // out, for tests/ft8_bench.cpp's self-comparison of the insides.
    void set_frontend_trace(std::function<void(const std::string&)> f) { fe_.set_trace(std::move(f)); }
    void flush_frontend_trace() { fe_.trace_flush(); }

private:
    // call() and period() for both sample types (S: int16_t or float).
    template <class S>
    Call call_t(const S* id2, int nzhsym, int nutc, const Params& p);
    template <class S>
    Period period_t(const S* samples, size_t n, int nutc, const Params& p, size_t fifty,
                    const std::function<void(int)>& after_call);
    // ft8_decode's decode(), with the callback.
    template <class S>
    void decode(const S* iwave, const Params& p, int nutc, int nzhsym, int ncontest, Call& c);
    // ft8_a7.f90's ft8_a7_save (msg is ft8_decode's msg37, rewritten in place
    // by split77 as the reference's argument is).
    void ft8_a7_save(int jseq, float dt, float f, std::string& msg);
    // ft8_a8d.f90 (with getmsg).
    void ft8_a8d(const float* dd, const std::string& mycall, const std::string& dxcall, const std::string& dxgrid,
                 float f1a, float& xdt, float& fbest, float& xsnr, float& plog, std::string& msgbest);
    // decoder_callbacks.f90's ft8_decoded.
    void callback(float sync, int snr, float dt, float freq, const std::string& decoded, int nap, float qual,
                  int nutc, int nfqso, Call& c);
    bool stopped() const { return stop_.load(std::memory_order_relaxed); }

    Frontend fe_;
    packjt77::Packjt77 pk_;
    class PackCodec;
    std::unique_ptr<PackCodec> codec_;

    // ft8_decode's SAVEd state
    std::vector<float> dd_, dd1_;
    int nutc0_ = -1;
    int ndec_early_ = 0;
    std::vector<int> itone_save_;  // itone_save(NN, 200)
    float f1_save_[200] = {}, xdt_save_[200] = {};
    bool lsubtracted_[200] = {};
    std::vector<std::string> allmessages_;
    // ft8_decode's local sync: on the stack, never initialised, so the value
    // the a7 decodes print is whatever the last candidate (or a8) left there
    // (a quirk of the reference, kept).
    float sync_ = 0.0f;
    std::vector<float> sbase_;  // ft8_decode's local sbase(NH1), the last sync8's
    std::vector<float> candidate_;

    // ft8_a7's module storage as the reference's object lays it out (its .bss,
    // the same in jt9 and the oracle's harnesses): msg0 (200, 0:1, 0:1),
    // character*37, then jseq, itone_a7 (never set: zero), f0 and dt0 (200,
    // 0:1, 0:1), then ft8_a7d's SAVEd logical one(0:511, 0:8). The tables are
    // flat in Fortran's column-major order (i, j, k) at (i-1) + 200*j + 400*k,
    // so indexes past 200 alias the next column (the a7 overflow) and indexes
    // past an array's 800 entries read and write what follows it, as there:
    // jseq 1's overflow (a quirk of the reference, kept). Offsets from msg0's
    // (nm -n -S of ft8_a7.f90.o: msg0 at 0x20, 0x73a0 bytes; jseq 0x73c0;
    // itone_a7 0x73e0; f0 0x7520 and dt0 0x81a0, 0xc80 bytes each; one.2
    // 0x8e20, 0x4800 bytes; the gaps zero). Exact up to 600 saves in a
    // period: past that the rotation's msg0 copy writes over jseq, which the
    // statements after it index ndec with (not modelled; decode() keeps the
    // jseq it read), and past 881 it reads beyond the block (the next
    // object's storage, not modelled: reads give zeros, writes are dropped).
    // Nothing here throws.
    static constexpr size_t A7_MSG0 = 0, A7_JSEQ = 0x73a0, A7_F0 = 0x7500, A7_DT0 = 0x8180, A7_ONE = 0x8e00,
                            A7_END = 0xd600;
    std::vector<unsigned char> a7mem_;
    bool a7d_first_ = true;  // ft8_a7d's first: one is filled by its first call
    int ndec_[2][2] = {};    // ndec(j, k)
    static int a7(int i, int j, int k) { return (i - 1) + 200 * j + 400 * k; }
    float a7_real(size_t base, int e) const;
    void set_a7_real(size_t base, int e, float v);
    float f0(int e) const { return a7_real(A7_F0, e); }
    float dt0(int e) const { return a7_real(A7_DT0, e); }
    void set_f0(int e, float v) { set_a7_real(A7_F0, e, v); }
    void set_dt0(int e, float v) { set_a7_real(A7_DT0, e, v); }
    std::string msg0(int e) const;
    void set_msg0(int e, const std::string& m);
    int jseq() const;
    void set_jseq(int v);

    // decoder.f90's
    bool ltry_a8_ = false;
    int ndec41_ = 0, ndec47_ = 0;
    bool decoded_exists_ = false;
    std::string decoded_text_;
    std::string decoded_pending_;  // this call's lines to decoded.txt
    int wrote_ = 0;

    // ft8_a8d's FFT buffers
    fst4w::FftwVec<cfloat> a8_cd0_, a8_csymb_;

    std::function<void(int, int, const packjt77::Packjt77&)> after_decode_;
    std::function<void(const std::string&)> trace_;
    std::atomic<bool> stop_{false};
};

// A WAV's 16-bit mono samples (the data chunk; 12000 Hz checked by the
// caller). False with why when it cannot be read.
bool read_wav(const std::string& path, std::vector<int16_t>& samples, int& rate, std::string& why);

// jt9's parse_wav_filename_nutc: the last six digits of the name before
// ".wav"/".WAV", else "_HHMM", else 0.
int nutc_of(const std::string& path);

// The oracle's harness (scripts/ft8-oracle.f90) as a function, for the tests
// and the standalone jt9 -8 (ft8_jt9.h):
//   -d D WAV...                       the GUI's decoding (jt9a), decoded-calls
//   --file [jt9 -8 options] WAV...    jt9 -8's file mode
// and WSJT-X's own "jt9 -a DIR -8 [options] WAV..." (as --file). Paths are
// relative to cwd. Fills what the harness leaves: stdout, decoded.txt (as
// the decoder leaves it), decoded-calls.txt (GUI mode), and with want_trace
// the oracle trace's C, P and A lines (Decoder::set_trace). Returns the exit
// status (2 with a message on err for bad usage or input).
struct RunOutput {
    bool want_trace = false;  // in: fill trace (the decoder's C, P and A lines)
    std::string trace;
    std::string out, err;
    bool decoded_exists = false;
    std::string decoded;
    std::string decoded_calls;
};
int run_periods(const std::vector<std::string>& args, const std::string& cwd, RunOutput& r);

}  // namespace ft8
