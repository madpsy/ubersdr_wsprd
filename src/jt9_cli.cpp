// jt9's command line (jt9_cli.h): lib/jt9.f90's option loop and checks,
// lib/options.f90's getopt and print_opt, lib/jt9_input_validation.f90's
// number checks, at 567ad29ce, with WSJT-X improved 3.2.0's --ft2; then
// each mode's file mode (FST4W's: fst4w_jt9.h; FT8's: ft8_jt9.h; FT4's:
// ft4_jt9.h; FT2's: ft2_jt9.h).
#include "jt9_cli.h"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>

#include "fst4w_jt9.h"
#include "ft2_jt9.h"
#include "ft4_jt9.h"
#include "ft8_jt9.h"
#include "fft.h"  // decoder_wisdom_*
#include "jt9_shm.h"

namespace jt9_cli {

namespace {

struct Opt {
    const char* name;  // long, without "--"
    bool has_arg;
    char chr;          // short; 0: a long option alone (ours)
    const char* descr;
    const char* argname;
};

// jt9's own options that the modes so far read, in jt9's order and words
// (improved 3.2.0's --ft2 where its jt9 lists it: before q65, so before
// ft4 here), and ours.
const Opt OPTS[] = {
    {"help", false, 'h', "Display this help message", ""},
    {"version", false, 'v', "Display version and build revision", ""},
    {"shmem", true, 's', "Use shared memory for sample data", "KEY"},
    {"tr-period", true, 'p', "Tx/Rx period, default SECONDS=60", "SECONDS"},
    {"executable-path", true, 'e', "Location of subordinate executables, default PATH=\".\"", "PATH"},
    {"data-path", true, 'a', "Location of writeable data files, default PATH=\".\"", "PATH"},
    {"temp-path", true, 't', "Temporary files path, default PATH=\".\"", "PATH"},
    {"lowest", true, 'L', "Lowest frequency decoded (JT65), default HERTZ=200", "HERTZ"},
    {"highest", true, 'H', "Highest frequency decoded, default HERTZ=4000", "HERTZ"},
    {"rx-frequency", true, 'f', "Receive frequency offset, default HERTZ=1500", "HERTZ"},
    {"freq-tolerance", true, 'F', "Receive frequency tolerance, default HERTZ=20", "HERTZ"},
    {"patience", true, 'w', "FFTW3 planing patience (0-4), default PATIENCE=1", "PATIENCE"},
    {"fft-threads", true, 'm', "Number of threads to process large FFTs, default THREADS=1", "THREADS"},
    {"ft2", false, '2', "FT2 mode", ""},
    {"ft4", false, '5', "FT4 mode", ""},
    {"fst4w", false, 'W', "FST4W mode", ""},
    {"fst4w", false, 'Y', "FST4W mode, print hash22 values", ""},
    {"ft8", false, '8', "FT8 mode", ""},
    {"quiet", false, 'q', "Quiet mode, no <DecodeFinished>", ""},
    {"QSOprog", true, 'Q', "QSO progress (0-5), default PROGRESS=1", "QSOprogress"},
    {"depth", true, 'd', "Decoding depth (1-3), default DEPTH=1", "DEPTH"},
    {"my-call", true, 'c', "my callsign", "CALL"},
    {"my-grid", true, 'G', "my grid locator", "GRID"},
    {"his-call", true, 'x', "his callsign", "CALL"},
    {"his-grid", true, 'g', "his grid locator", "GRID"},
    {"experience-decode", true, 'X', "experience based decoding flags (1..n), default FLAGS=0", "FLAGS"},
    {"fast-subtract", false, 0,
     "Ours (FT8, FT4, FT2): the decoded signals' subtraction smoothed at a low rate, not with two big FFTs: the "
     "same decodes as WSJT-X's, not to the bit, for less CPU",
     ""},
    {"feed", false, 0,
     "Ours (FT8, FT4, FT2): read WAV paths from stdin, one a line (PATH or PATH HHMMSS), and decode each as it arrives "
     "with one decoder for all, as one jt9 -s does period after period; after each file's lines and "
     "<DecodeFinished>, <FileDone> PATH N",
     ""},
    {"spread", true, 0,
     "FST4W's spread (w50): fast (default: within 0.0001 Hz of jt9's, several times cheaper), exact (jt9's, "
     "bit for bit) or off (none)",
     "fast|exact|off"},
};
// getopt's option letters, ':' after those taking a value. No '2':
// improved's getopt string has none, so only --ft2 selects FT2 and "-2" is
// an unrecognised option there (the usage, status 0), as here.
// -s, -e, -r, -m and -w as 567ad29ce's getopt string has them (-r, its
// share path, is not in its list of options).
const char* const SHORT = "hvs:e:a:r:m:p:d:f:F:w:t:L:H:5WY8qQ:c:G:x:g:X:";

std::string rtrim(std::string s) {
    while (!s.empty() && s.back() == ' ') s.pop_back();
    return s;
}

// options.f90's getopt over args (1-based positions, as Fortran's).
class Getopt {
public:
    explicit Getopt(const std::vector<std::string>& a) : args_(a), cnt_(static_cast<int>(a.size())) {}
    struct Got {
        int st = 0;        // 0 an option, 1 none left (or "--"), -1 unknown, -2 its value missing
        char chr = 0;
        const Opt* opt = nullptr;  // the long option it was given as, or none
        std::string optarg;
        int offset = 0, remain = 0;
    };
    Got next(std::ostream& err) {
        Got g;
        bool is_long = false;
        std::string arg, longname;
        char chr = 0;
        auto get = [&](int p) { return args_[static_cast<size_t>(p - 1)]; };
        auto at = [&](size_t i) { return i < arg.size() ? arg[i] : ' '; };  // Fortran's blank padding
        auto done = [&](int st) {
            if (st == -1) err << " error: unrecognised option: " << rtrim(arg) << "\n";
            if (st == -2) {
                if (!is_long) err << " error: option -" << chr << " requires an argument\n";
                else err << " error: option --" << longname << " requires an argument\n";
            }
            g.st = st;
            g.chr = chr;
            g.optarg = arg;
            g.offset = pos_;
            g.remain = cnt_ - pos_;
            ++pos_;
            return g;
        };
        if (pos_ > cnt_) {
            --pos_;
            return done(1);
        }
        arg = get(pos_);
        if (at(0) != '-') {
            --pos_;  // the first file: back to the last option or value
            return done(1);
        }
        chr = at(1);
        if (chr != '-' && rtrim(arg).size() > 2) return done(-1);
        if (chr == '-' && at(2) == ' ') return done(1);  // "--"
        if (chr == '-') {
            is_long = true;
            const std::string name = rtrim(arg.substr(2));
            for (const Opt& o : OPTS)
                if (name == o.name) {
                    g.opt = &o;
                    break;
                }
            if (!g.opt) return done(-1);
            chr = g.opt->chr;
            longname = g.opt->name;
            if (!g.opt->has_arg) return done(0);
            if (pos_ == cnt_) return done(-2);
            ++pos_;
            arg = get(pos_);
            if (at(0) == '-') {
                --pos_;
                return done(-2);
            }
            if (!chr) return done(0);  // ours: its value, once
            // jt9's getopt now checks the letter as for a short option,
            // so a value taken goes and the word after it is taken.
        }
        const char* p = chr ? std::strchr(SHORT, chr) : nullptr;
        if (!p) return done(-1);
        if (p[1] == ':') {
            if (pos_ == cnt_) return done(-2);
            ++pos_;
            arg = get(pos_);
            if (at(0) == '-') {
                --pos_;
                return done(-2);
            }
        }
        return done(0);
    }

private:
    const std::vector<std::string>& args_;
    int cnt_;
    int pos_ = 1;
};

// jt9_input_validation's tokens.
bool integer_token(const std::string& t) {
    size_t i = !t.empty() && (t[0] == '+' || t[0] == '-') ? 1 : 0;
    if (i >= t.size()) return false;
    for (; i < t.size(); ++i)
        if (t[i] < '0' || t[i] > '9') return false;
    return true;
}
bool real_token(const std::string& t) {
    size_t i = 0, n = t.size();
    if (n == 0) return false;
    if (t[i] == '+' || t[i] == '-') ++i;
    if (i >= n) return false;
    int digits = 0;
    while (i < n && t[i] >= '0' && t[i] <= '9') ++digits, ++i;
    if (i < n && t[i] == '.') {
        ++i;
        while (i < n && t[i] >= '0' && t[i] <= '9') ++digits, ++i;
    }
    if (digits == 0) return false;
    if (i < n) {
        if (t[i] != 'e' && t[i] != 'E') return false;
        ++i;
        if (i < n && (t[i] == '+' || t[i] == '-')) ++i;
        if (i >= n) return false;
        digits = 0;
        while (i < n) {
            if (t[i] < '0' || t[i] > '9') return false;
            ++digits, ++i;
        }
        if (digits == 0) return false;
    }
    return i >= n;
}

// Fortran's list-directed read of one character value into a character*len
// (12 for a call, 6 for a grid): blanks skipped, a quoted value to its
// closing quote, else to a blank, comma or slash. Empty or a separator
// first: none read (the value as it was).
bool read_list_chars(const std::string& t, std::string& v, size_t len = 12) {
    size_t i = 0;
    while (i < t.size() && t[i] == ' ') ++i;
    if (i >= t.size() || t[i] == ',' || t[i] == '/') return false;
    std::string r;
    if (t[i] == '\'' || t[i] == '"') {
        const char q = t[i++];
        for (; i < t.size(); ++i) {
            if (t[i] == q) {
                if (i + 1 < t.size() && t[i + 1] == q) {
                    r += q;
                    ++i;
                    continue;
                }
                break;
            }
            r += t[i];
        }
    } else {
        for (; i < t.size() && t[i] != ' ' && t[i] != ',' && t[i] != '/'; ++i) r += t[i];
    }
    v = rtrim(r.substr(0, len));
    return true;
}

// print_opt: "-c ARG", "--name ARG", then the description in lines of at
// most 76 characters after 4 blanks (broken at a blank; each line its full
// width, the blanks after the words too, as Fortran writes a substring).
void print_opt(std::ostream& out, const Opt& o) {
    if (o.chr) out << " -" << o.chr << (o.has_arg ? std::string(" ") + o.argname : "") << "\n";
    out << " --" << o.name << (o.has_arg ? std::string(" ") + o.argname : "") << "\n";
    std::string d(o.descr);
    const size_t l = d.size();
    d.resize(500, ' ');
    size_t c1 = 1;
    while (c1 <= l) {
        size_t c2 = std::min<size_t>(c1 + 76, 500);
        if (c2 != 500)
            while (d[c2 - 1] != ' ') --c2;
        out << "    " << d.substr(c1 - 1, c2 - c1) << "\n";
        c1 = c2 + 1;
    }
}

void usage(std::ostream& out) {
    out << " Usage: jt9 [OPTIONS] file1 [file2 ...]\n";
    out << "        Reads data from *.wav files.\n";
    out << " \n";
    out << "        jt9 -s <key> [-w patience] [-m threads] [-e path] [-a path] [-t path] [-r path]\n";
    out << "        Gets data from shared memory region with key==<key>\n";
    out << " \n";
    out << " OPTIONS:\n";
    out << " \n";
    for (const Opt& o : OPTS) print_opt(out, o);
}

// What jt9.f90 does once its options are read: -s, or each mode's files.
int decode(bool read_files, const std::string& shm_key, int mode, fst4w::FileMode& fm, ft8::FileMode& f8,
           const std::vector<std::string>& args, const Getopt::Got& g, const std::string& cwd, std::ostream& out,
           std::ostream& err, std::ostream* live_err) {
    if (!read_files) {
        jt9_shm::Options so;
        so.key = shm_key;
        so.temp_dir = f8.temp_dir;
        return jt9_shm::run(so, cwd, out, live_err ? *live_err : err);
    }
    std::vector<std::string> files(args.begin() + g.offset, args.begin() + g.offset + g.remain);
    // The modes so far.
    if (mode == 241 || mode == 242) {
        fm.hash22 = mode == 242;
        return fst4w::run_files(fm, files, cwd, out, err);
    }
    if (mode == 8) {
        // Fox (ncontest 6) writes houndcallers.txt, which this jt9 does not.
        if ((f8.nexp_decode & 7) == 6) {
            err << "jt9: -X " << f8.nexp_decode << ": FT8's Fox mode (ncontest 6) is not taken by this jt9\nSTOP 2\n";
            return 2;
        }
        return ft8::run_files(f8, files, cwd, out, err);
    }
    if (mode == 5) return ft4::run_files(f8, files, cwd, out, err);
    if (mode == 52) return ft2::run_files(f8, files, cwd, out, err);
    err << "jt9: no mode given: this jt9 decodes FST4W (--fst4w, -W, or -Y), FT8 (-8), FT4 (-5) and FT2 (--ft2) so "
           "far\nSTOP 2\n";
    return 2;
}

}  // namespace

int run(const std::vector<std::string>& args, const std::string& cwd, std::ostream& out, std::ostream& err,
        const std::string& version, std::ostream* live_err) {
    fst4w::FileMode fm;
    ft8::FileMode f8;  // what FT8, FT4 and FT2 read, beside FST4W's
    int mode = 0;  // jt9's: 241 --fst4w (-W), 242 -Y, 8 -8, 5 -5, 52 --ft2 (improved's)
    bool display_help = false, display_version = false;
    bool read_files = true;  // false: -s
    bool feed = false;       // --feed: the files' names from stdin
    std::string shm_key;
    int nthreads = 1, npatience = 1;  // -m read; -w the decoders' kept FFT plans' patience (fft.h)
    Getopt getopt(args);
    Getopt::Got g;
    // jt9.f90's invalid_option: said, status 2.
    auto invalid = [&](char o, const std::string& text, const char* reason) {
        err << "jt9: invalid value for -" << o << ": " << rtrim(text) << " (" << reason << ")\nSTOP 2\n";
        return 2;
    };
    auto integer = [&](char o, const std::string& text, int& v) {
        if (!integer_token(text)) return false;
        errno = 0;
        char* end = nullptr;
        const long x = std::strtol(text.c_str(), &end, 10);
        if (errno || x < INT32_MIN || x > INT32_MAX) return false;  // the read fails
        v = static_cast<int>(x);
        (void)o;
        return true;
    };
    for (;;) {
        g = getopt.next(err);
        if (g.st != 0) break;
        const std::string& v = g.optarg;
        if (g.opt && !g.opt->chr && std::strcmp(g.opt->name, "fast-subtract") == 0) {
            f8.fast_subtract = true;
            continue;
        }
        if (g.opt && !g.opt->chr && std::strcmp(g.opt->name, "feed") == 0) {
            feed = true;
            continue;
        }
        if (g.opt && !g.opt->chr) {  // --spread
            if (v == "fast") fm.spread = fst4w::Options::Spread::Fast;
            else if (v == "exact") fm.spread = fst4w::Options::Spread::Exact;
            else if (v == "off") fm.spread = fst4w::Options::Spread::Off;
            else if (v == "both") fm.spread = fst4w::Options::Spread::Both;  // for comparing: exact's, fast's to stderr
            else {
                err << "jt9: invalid value for --spread: " << v << " (fast, exact or off)\nSTOP 2\n";
                return 2;
            }
            continue;
        }
        switch (g.chr) {
        case 'h': display_help = true; break;
        case 'v': display_version = true; break;
        case 's':
            read_files = false;
            shm_key = v;
            break;
        case 'e': break;  // exe_dir: no subordinate programs here
        case 'r': break;  // share_dir: no shipped data read here
        case 'm':
            if (!integer('m', v, nthreads)) return invalid('m', v, "must be an integer");
            break;
        case 'w':
            if (!integer('w', v, npatience)) return invalid('w', v, "must be an integer");
            break;
        case 'a': fm.data_dir = f8.data_dir = v; break;
        case 't': fm.temp_dir = f8.temp_dir = v; break;
        case 'p': {
            if (!real_token(v)) return invalid('p', v, "must be a finite number");
            const double p = std::strtod(v.c_str(), nullptr);
            if (!std::isfinite(p)) return invalid('p', v, "must be a finite number");
            if (p <= 0 || p > 1800) return invalid('p', v, "must be between 0 and 1800 seconds");
            fm.trperiod = f8.trperiod = p;
            break;
        }
        case 'd':
            if (!integer('d', v, fm.ndepth)) return invalid('d', v, "must be an integer");
            f8.ndepth = fm.ndepth;
            break;
        case 'f':
            if (!integer('f', v, fm.nrxfreq)) return invalid('f', v, "must be an integer");
            f8.nrxfreq = fm.nrxfreq;
            break;
        case 'L':
            if (!integer('L', v, f8.flow)) return invalid('L', v, "must be an integer");
            break;
        case 'H':
            if (!integer('H', v, f8.fhigh)) return invalid('H', v, "must be an integer");
            break;
        case 'Q':
            if (!integer('Q', v, f8.nqso_progress)) return invalid('Q', v, "must be an integer");
            if (f8.nqso_progress < 0 || f8.nqso_progress > 5) return invalid('Q', v, "must be between 0 and 5");
            break;
        case 'F':
            if (!integer('F', v, fm.ntol)) return invalid('F', v, "must be an integer");
            break;
        case 'X':
            if (!integer('X', v, fm.nexp_decode)) return invalid('X', v, "must be an integer");
            fm.nexp_decode_set = true;
            f8.nexp_decode = fm.nexp_decode;
            break;
        case 'W': mode = 241; break;
        case 'Y': mode = 242; break;
        case '8': mode = 8; break;
        case '5': mode = 5; break;
        case '2': mode = 52; break;  // --ft2 only (getopt refuses -2)
        case 'q': fm.quiet = f8.quiet = true; break;
        case 'c':
        case 'x':
        case 'G':
        case 'g': {
            // read (optarg(:arglen), *) at jt9.f90's line for it: an empty
            // or blank value is the end of the internal file, a runtime
            // error (said as gfortran says it, its backtrace left out).
            if (v.find_first_not_of(' ') == std::string::npos) {
                const int line = g.chr == 'c' ? 233 : g.chr == 'G' ? 235 : g.chr == 'x' ? 237 : 239;
                err << "At line " << line << " of file jt9.f90\nFortran runtime error: End of file\n";
                return 2;
            }
            if (g.chr == 'c') read_list_chars(v, fm.mycall), f8.mycall = fm.mycall;
            if (g.chr == 'x') read_list_chars(v, fm.hiscall), f8.hiscall = fm.hiscall;
            if (g.chr == 'g') read_list_chars(v, f8.hisgrid, 6);
            if (g.chr == 'G') {  // mygrid: read (FT8 decodes without it)
                std::string mygrid;
                read_list_chars(v, mygrid, 6);
            }
            break;
        }
        default: break;
        }
    }
    if (display_version) {
        out << "jt9 " << version
            << " (FST4W, FT8 and FT4: WSJT-X 567ad29ce's decoders; FT2: WSJT-X improved 3.2.0's FT2 on that FT4)\n";
        return 0;
    }
    if (feed && (!read_files || g.remain > 0 || (mode != 8 && mode != 5 && mode != 52))) {
        err << "jt9: --feed takes FT8 (-8), FT4 (-5) or FT2 (--ft2), no -s and no files\nSTOP 2\n";
        return 2;
    }
    if (display_help || g.st < 0 || (!read_files && g.remain > 0) || (read_files && g.remain < 1 && !feed)) {
        usage(out);
        return 0;
    }
    // jt9.f90: FFTW's wisdom imported from <data_dir>/jt9_wisdom.dat
    // before decoding, exported to it at the end; a STOP (status 2) on the
    // way ends jt9 without exporting it. Ours: a -w other than 1 keeps its
    // own, jt9_wisdom-wN.dat: FFTW plans with any wisdom made at the same
    // patience or more, so a measured plan kept in jt9_wisdom.dat would be
    // taken by every later -w 1 run in that directory, which would then no
    // longer decode as WSJT-X's does to the bit.
    const std::string wisfile = ft8::jt9_join(
        ft8::jt9_at(cwd, f8.data_dir),
        npatience == 1 ? std::string("jt9_wisdom.dat") : "jt9_wisdom-w" + std::to_string(npatience) + ".dat");
    skimmer::set_decoder_plan_patience(npatience);
    skimmer::decoder_wisdom_on(wisfile);
    int status;
    if (feed) {
        // Runs until stdin ends: what it says goes out as it happens.
        std::ostream& e = live_err ? *live_err : err;
        status = mode == 8   ? ft8::run_feed(f8, std::cin, cwd, out, e)
                 : mode == 5 ? ft4::run_feed(f8, std::cin, cwd, out, e)
                             : ft2::run_feed(f8, std::cin, cwd, out, e);
    } else {
        status = decode(read_files, shm_key, mode, fm, f8, args, g, cwd, out, err, live_err);
    }
    if (status == 0) skimmer::decoder_wisdom_export(wisfile);
    skimmer::decoder_wisdom_off();
    return status;
}

std::string stderr_as_file(const std::string& err) {
    static const std::string stop = "STOP 2\n";
    if (err.size() < stop.size() || err.compare(err.size() - stop.size(), stop.size(), stop) != 0) return err;
    return stop + err.substr(0, err.size() - stop.size());
}

}  // namespace jt9_cli
