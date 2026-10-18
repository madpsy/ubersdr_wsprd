// WSPR decoding in-process: WSJT-X's wsprd (lib/wsprd at 567ad29ce, with
// lib/indexx.f90), giving the reference's results to the bit: the same
// decodes, SNR, DT, frequency, drift and sync from the same samples, and,
// through run(), the same stdout, ALL_WSPR.TXT, wspr_spots.txt,
// hashtable.txt and -c file as the wsprd binary, byte for byte
// (tests/data/wsprd holds the reference's outputs; scripts/wsprd-compare.sh
// compares). wsprd (wsprd_main.cpp) is that binary, with -n besides.
//
// What it rests on, outside the source: the wsprd files are compiled with
// no multiply-adds fused (-ffp-contract=off, CMakeLists.txt), as WSJT-X's
// Release build of wsprd.c is on x86_64 (plain x86-64 has no FMA), under
// -march=x86-64-v3 (whose vectorising, with no
// reassociation, changes no result), and on every other platform too: on
// aarch64 that is the reference built with -ffp-contract=off (checked on a
// Pi 5; WSJT-X's own Release build fuses there, which moves sync values in
// the last bits but changed no decode on the fixtures); its maths library calls
// are the reference's (double cos/sin/log10/pow/sqrt of float arguments);
// and its FFTW plans are what FFTW_ESTIMATE makes with only wsprd's own
// wisdom in the planner (Decoder sets the process's aside while it plans).
// The .c2 path's only FFT is the 512-point one; a .wav's first two FFTs'
// bits depend on the FFTW build and CPU (its codelets) as the reference's
// do.
//
// wsprd.c, wsprd_utils.c, wsprsim_utils.c: Copyright 2001-2018, Joe Taylor,
// K1JT, and 2014-2018, Steven Franke, K9AN; fano.c: 1994, Phil Karn, KA9Q;
// WSJT-X is GNU GPL v3, and so is this.
#pragma once

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace wsprd {

// wsprd's command line: every option, with wsprd's defaults. Fields hold
// what its getopt switch stores, in the reference's types.
struct Options {
    std::string data_dir = ".";       // -a <path>: where its files live (run() only)
    int npasses = 3;                  // -B: 2 (no block demodulation pass); -s: 1
    bool subtraction = true;          // -s: false
    bool writec2 = false;             // -c: the IQ after the first pass (run(): 000000_0001.c2)
    unsigned int maxcycles = 10000;   // -C n: Fano/stack cycles per bit
    bool more_candidates = false;     // -d: deeper search
    double dialfreq_error = 0.0;      // -e x: dial error in Hz
    double dialfreq = 0.0;            // -f x: dial in MHz (a .wav's; a .c2 carries its own)
    bool usehashtable = true;         // -H: false (run(): no hashtable.txt; Decoder: a fresh table each decode)
    bool stackdecoder = false;        // -J: the stack decoder instead of Fano
    int wspr_type = 2;                // -m: 15 (WSPR-15 .wav)
    int ndepth = -1;                  // -o n: OSD depth, -1 none
    bool quickmode = false;           // -q
    bool verbose = false;             // -v: dupes too
    float fmin = -110, fmax = 110;    // -w: -150, 150
    float bias = 0.45f;               // -z x: Fano metric bias
    bool stream_mode = false;         // -0: framed PCM on stdin, NDJSON out (run() only)
    std::string infile;               // the first operand: a .wav or .c2 (run() only)
    int infile_arg = -1;              // its index in parse_args' args (run() only)
    // Not wsprd's (no option sets it, nothing it prints reads it): each
    // decode's spread measured (Decode::spread).
    bool spread = true;
    // Not wsprd's: -n, taken only with_n (the wsprd binary): the candidates'
    // coarse drift search held at 0 in every pass (the first passes' +-0.5
    // Hz drift refinement still runs, so a decode's drift may be +-0.5,
    // printed 0); and ALL_WSPR.TXT's lines end in the decode's spread
    // (" %6.3f").
    bool nodrift = false;
};

// usage() as wsprd prints it (to stdout); with_n with an -n line after
// -z's.
std::string usage(bool with_n = false);

// wsprd's getopt over args (its argv after the program name, prog that name
// for getopt's messages), as glibc's getopt runs it: clustered flags (-wd),
// values attached or next (-C500, -C 500), operands anywhere (permuted to
// the end, unless POSIXLY_CORRECT is set) and -- ending the options, the
// values converted as wsprd does (strtoul, strtod, strtol: "-o x" is 0).
// Fills *opt from a default Options. On an unknown option or a missing
// value returns false with getopt's message ("prog: invalid option -- 'x'",
// newline included) in *error: wsprd then prints usage() and exits 1. A
// missing operand is not an error here (infile stays empty): wsprd checks
// that later, after its data directory. with_n: -n taken too (nodrift),
// else it is an invalid option as in wsprd.
bool parse_args(const std::vector<std::string>& args, Options* opt, std::string* error, const std::string& prog = "wsprd",
                bool with_n = false);

// One decode: a line of ALL_WSPR.TXT, all wsprd keeps of it.
struct Decode {
    std::string date;          // YYMMDD and HHMM, from the input's name
    std::string time;
    float sync = 0;
    float snr = 0;
    float dt = 0;
    double freq = 0;           // MHz
    std::string message;       // "K1ABC FN42 37", "PJ4/K1ABC 37", "<K1ABC> FN42MI 37", "<...> FN42MI 37"
    float drift = 0;
    unsigned int cycles = 0;   // ALL_WSPR.TXT prints cycles/81
    int jitter = 0;
    int blocksize = 0;         // 1-3; 4 the bit-by-bit metric
    unsigned int metric = 0;
    int nhardmin = 0;
    int ipass = 0;             // 0-based; the files print ipass+1
    int decodetype = 0;        // 1 for an OSD decode
    // Not wsprd's: the signal's spread, as WSJT-X measures FST4W's (jt9's
    // w50, dopspread in lib/fst4_decode.f90): the width holding the middle
    // half of the power of its channel gain (the samples times the decoded
    // signal's conjugate, before it was subtracted), Hz, at least one bin
    // (1/240 Hz, WSPR-15's 1/1800). -1 none (Options::spread off).
    float spread = -1;
    // The message's fields: callsign ("<...>" for an unknown hash), grid
    // (empty for type 2, 6 characters for type 3), power in dBm.
    std::string callsign;
    std::string grid;
    int power = 0;
};

// The lines wsprd prints for a decode (each with its newline): stdout's,
// ALL_WSPR.TXT's and wspr_spots.txt's, and -0's NDJSON. with_spread: the
// ALL_WSPR.TXT line ends in " %6.3f" of d.spread (-n's).
std::string stdout_line(const Decode& d);
std::string all_wspr_line(const Decode& d, bool with_spread = false);
std::string spots_line(const Decode& d);
std::string stream_line(const Decode& d);

// The call hash table: hashtable.txt, which type 3 messages (<CALL>) and OSD
// decodes are resolved through.
class Hashtable {
public:
    Hashtable();
    void clear();
    // As wsprd reads hashtable.txt: fgets() 80 bytes at a time, each line to
    // wsprd_load_hash_line (a malformed line is skipped).
    void load(const std::string& text);
    bool load_line(const char* line);
    // As wsprd writes it: "%5d %s %s\n" for each slot with a call.
    std::string save() const;
    // The call and grid in slot ihash (0..32767), empty if none.
    std::string call(int ihash) const;
    std::string grid(int ihash) const;
    // nhash(call) & 0x7fff: the slot a call is kept in.
    static int hash(const std::string& call);
    char* hashtab() { return hashtab_.data(); }
    char* loctab() { return loctab_.data(); }

private:
    std::vector<char> hashtab_;  // 32768 x 13
    std::vector<char> loctab_;   // 32768 x 5
};

// wsprd's timers (wspr_timer.out), seconds of clock().
struct Timings {
    float treadwav = 0, tcandidates = 0, tsync0 = 0, tsync1 = 0, tsync2 = 0, tfano = 0, tosd = 0, ttotal = 0;
};

// A .c2 file's contents (readc2file): 45000 complex samples at 375 Hz, q
// already negated as wsprd's reader does, the dial in MHz and the type (2
// or 15) from its header.
struct C2 {
    std::vector<float> idat, qdat;
    double freq = 0;
    int wspr_type = 2;
};

struct Result {
    std::vector<Decode> decodes;  // as wsprd prints them, by frequency
    // With writec2: the samples after the first pass, as writec2file writes
    // them (45000 of each, q not yet negated), and the header's dial.
    bool has_c2 = false;
    std::vector<float> c2_idat, c2_qdat;
    double c2_freq = 0;
    // Where wsprd would stop instead, its exit status: 1 for an input it
    // refuses (nothing decoded); 139 where it crashes mid-decode (the stack
    // decoder's metric out of its bucket array: -J with a large -z), with
    // no decodes, has_c2 if the first pass had ended (its file written),
    // and the Decoder left as the crashed run leaves its data directory:
    // hash table and wisdom as before the decode, timers at 0.
    int status = 0;
    // Ended early by the Decoder's stop flag (set_stop()), as a killed
    // wsprd ends: it writes ALL_WSPR.TXT and hashtable.txt only after its
    // last pass, so no decodes, and the Decoder as a crash leaves it (hash
    // table and wisdom as before the decode, timers at 0).
    bool stopped = false;
};

// What wsprd keeps across its run, and across runs in one data directory:
// the hash table and FFTW's wisdom; and its timers. Each decode*() is one
// wsprd run from its input on. A Decoder decodes one period at a time;
// separate Decoders may decode at once on different threads.
class Decoder {
public:
    explicit Decoder(const Options& opt = Options());
    ~Decoder();
    Decoder(const Decoder&) = delete;
    Decoder& operator=(const Decoder&) = delete;

    Options& options() { return opt_; }
    const Options& options() const { return opt_; }
    // Used and updated by each decode unless !options().usehashtable.
    Hashtable& hashtable() { return hashtable_; }
    // FFTW wisdom (FFTW's text) imported before planning and replaced by
    // what the planner holds after: wspr_wisdom.dat. Empty: none.
    std::string& wisdom() { return wisdom_; }
    // Added to by each decode, as wsprd adds to wspr_timer.out's: carried
    // from one decode to the next as that file carries them ("%7.2f" read
    // back), zeroed by a refused or crashed decode (wsprd has emptied the
    // file by then). The calling thread's CPU time (see cpu_clock()).
    Timings& timings() { return timings_; }
    // A flag (not wsprd's) another thread may set to end a decode early:
    // looked at before planning, before the passes, and at each pass,
    // candidate and try, so within a try's time (some milliseconds; but
    // WSPR-15's planning and first FFT, about a second, run whole); the
    // decode then ends as a killed wsprd does (Result::stopped). Never set,
    // results are wsprd's. nullptr: none.
    void set_stop(const std::atomic<bool>* stop) { stop_ = stop; }
    // Bytes held: the Decoder's buffers (kept from decode to decode), its
    // hash table, and the FFTW buffers and memos its last decode made (gone
    // once it returned): what a decode holds at most, near enough (FFTW's
    // plans' own tables and the stack aside). Not while a decode runs.
    size_t memory() const;
    // What its next decode will hold at most, near enough, by its options
    // (the mode, -J, OSD): memory() as that decode would leave it, for a
    // Decoder that has not decoded yet. Never less than memory().
    size_t memory_needed() const;
    // Frees the buffers a decode makes again (the samples, spectra and
    // scratch: a WSPR period's 5 MB, WSPR-15's 25), keeping the hash table,
    // wisdom and timers: for a Decoder kept idle between periods. Changes
    // no result.
    void release();

    // A period of 12 kHz audio (readwavfile): samples the .wav's data after
    // its 44-byte header, n of them (wsprd reads 114 s, 912 s for WSPR-15,
    // zero-filling what is short; more are ignored); the dial
    // options().dialfreq, WSPR-15 if options().wspr_type is 15. date and
    // time are YYMMDD and HHMM (only their first 6 and 4 characters are
    // kept, as wsprd's fixed fields hold). Refused (status 1) as wsprd
    // refuses it: no samples, or a type other than 2 or 15.
    Result decode_wav(const int16_t* samples, size_t n, const std::string& date, const std::string& time);
    // The same from a 32-bit float .wav's samples, full scale 1.0 (a 16-bit
    // sample s is s/32768 of it): decoded as the 16-bit samples are from
    // the FFT on, without their rounding to 16 bits.
    Result decode_wav(const float* samples, size_t n, const std::string& date, const std::string& time);
    // A .c2 file's samples (readc2file), its own dial and type. Refused
    // (status 1) with fewer than 45000 of either, as wsprd refuses a short
    // file.
    Result decode_c2(const C2& c2, const std::string& date, const std::string& time);

private:
    struct Impl;
    struct Kept;
    template <class Sample>
    Result decode_samples(const Sample* samples, size_t n, const std::string& date, const std::string& time);
    Options opt_;
    Hashtable hashtable_;
    std::string wisdom_;
    Timings timings_;
    bool decoded_ = false;  // a decode has been made: timings_ to carry
    const std::atomic<bool>* stop_ = nullptr;
    void carry_timings();
    bool stopping(Result& res) const;  // the stop flag set: res.stopped
    std::unique_ptr<Impl> impl_;
};

// readc2file: false if the file cannot be read (*error says why when wsprd
// prints something: "Cannot open data file ...") or is short.
bool read_c2(const std::string& path, C2* c2, std::string* error);
// writec2file's bytes: its 14-byte name, type, dial, then the samples.
std::string c2_bytes(const char* name14, int type, double freq, const float* idat, const float* qdat);

// wsprd's main() in full, its files under options' -a, the -c file in cwd,
// and a relative input or -a taken as cwd's (empty: the current directory). out and err get what it writes to
// stdout and stderr, in order (flush: where wsprd flushes stdout); in is
// stdin for -0. Returns wsprd's exit status.
struct Io {
    std::function<void(const std::string&)> out;
    std::function<void(const std::string&)> err;
    std::function<void()> flush;
    std::FILE* in = nullptr;
    std::string cwd;
    // Where wsprd crashes (SIGSEGV): called after all that wsprd wrote to
    // stdout and its files is written, to end the process there without
    // flushing stdio, so what the crash loses is what the C library had
    // not yet written (if it returns, run() returns 139). Unset: run()
    // returns 139, and out gets only what a fully buffered stdout (4096
    // bytes) would have written; so do the files.
    std::function<void()> crash;
    // Not wsprd's: -n taken (parse_args' with_n) and usage() with its line,
    // as the wsprd binary's; -n given, nodrift on and ALL_WSPR.TXT's spread
    // column printed. Off: wsprd's command line exactly.
    bool with_n = false;
};
int run(const std::vector<std::string>& argv, const Io& io);

}  // namespace wsprd
