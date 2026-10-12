// jt9's command line (jt9_cli.h): lib/jt9.f90's option loop and checks,
// lib/options.f90's getopt and print_opt, lib/jt9_input_validation.f90's
// number checks, at 567ad29ce; then each mode's file mode (FST4W's:
// fst4w_jt9.h).
#include "jt9_cli.h"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <string>

#include "fst4w_jt9.h"

namespace jt9_cli {

namespace {

struct Opt {
    const char* name;  // long, without "--"
    bool has_arg;
    char chr;          // short; 0: a long option alone (ours)
    const char* descr;
    const char* argname;
};

// jt9's own options that the modes so far read, in jt9's order and words,
// and ours.
const Opt OPTS[] = {
    {"help", false, 'h', "Display this help message", ""},
    {"version", false, 'v', "Display version and build revision", ""},
    {"tr-period", true, 'p', "Tx/Rx period, default SECONDS=60", "SECONDS"},
    {"data-path", true, 'a', "Location of writeable data files, default PATH=\".\"", "PATH"},
    {"temp-path", true, 't', "Temporary files path, default PATH=\".\"", "PATH"},
    {"rx-frequency", true, 'f', "Receive frequency offset, default HERTZ=1500", "HERTZ"},
    {"freq-tolerance", true, 'F', "Receive frequency tolerance, default HERTZ=20", "HERTZ"},
    {"fst4w", false, 'W', "FST4W mode", ""},
    {"fst4w", false, 'Y', "FST4W mode, print hash22 values", ""},
    {"quiet", false, 'q', "Quiet mode, no <DecodeFinished>", ""},
    {"depth", true, 'd', "Decoding depth (1-3), default DEPTH=1", "DEPTH"},
    {"my-call", true, 'c', "my callsign", "CALL"},
    {"his-call", true, 'x', "his callsign", "CALL"},
    {"experience-decode", true, 'X', "experience based decoding flags (1..n), default FLAGS=0", "FLAGS"},
    {"spread", true, 0,
     "FST4W's spread (w50): fast (default: within 0.0001 Hz of jt9's, several times cheaper), exact (jt9's, "
     "bit for bit) or off (none)",
     "fast|exact|off"},
};
// getopt's option letters, ':' after those taking a value.
const char* const SHORT = "hva:p:d:f:F:t:WYqc:x:X:";

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

// Fortran's list-directed read of one character value into len 12: blanks
// skipped, a quoted value to its closing quote, else to a blank, comma or
// slash. Empty or a separator first: none read (the value as it was).
bool read_list_chars(const std::string& t, std::string& v) {
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
    v = rtrim(r.substr(0, 12));
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
    out << " OPTIONS:\n";
    out << " \n";
    for (const Opt& o : OPTS) print_opt(out, o);
}

}  // namespace

int run(const std::vector<std::string>& args, const std::string& cwd, std::ostream& out, std::ostream& err,
        const std::string& version) {
    fst4w::FileMode fm;
    int mode = 0;  // jt9's: 241 --fst4w (-W), 242 -Y
    bool display_help = false, display_version = false;
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
        case 'a': fm.data_dir = v; break;
        case 't': fm.temp_dir = v; break;
        case 'p': {
            if (!real_token(v)) return invalid('p', v, "must be a finite number");
            const double p = std::strtod(v.c_str(), nullptr);
            if (!std::isfinite(p)) return invalid('p', v, "must be a finite number");
            if (p <= 0 || p > 1800) return invalid('p', v, "must be between 0 and 1800 seconds");
            fm.trperiod = p;
            break;
        }
        case 'd':
            if (!integer('d', v, fm.ndepth)) return invalid('d', v, "must be an integer");
            break;
        case 'f':
            if (!integer('f', v, fm.nrxfreq)) return invalid('f', v, "must be an integer");
            break;
        case 'F':
            if (!integer('F', v, fm.ntol)) return invalid('F', v, "must be an integer");
            break;
        case 'X':
            if (!integer('X', v, fm.nexp_decode)) return invalid('X', v, "must be an integer");
            fm.nexp_decode_set = true;
            break;
        case 'W': mode = 241; break;
        case 'Y': mode = 242; break;
        case 'q': fm.quiet = true; break;
        case 'c': read_list_chars(v, fm.mycall); break;
        case 'x': read_list_chars(v, fm.hiscall); break;
        default: break;
        }
    }
    if (display_version) {
        out << "jt9 " << version << " (FST4W: WSJT-X 567ad29ce's decoder)\n";
        return 0;
    }
    if (display_help || g.st < 0 || g.remain < 1) {
        usage(out);
        return 0;
    }
    std::vector<std::string> files(args.begin() + g.offset, args.begin() + g.offset + g.remain);
    // The modes so far.
    if (mode == 241 || mode == 242) {
        fm.hash22 = mode == 242;
        return fst4w::run_files(fm, files, cwd, out, err);
    }
    err << "jt9: no mode given: this jt9 decodes FST4W (--fst4w, -W, or -Y) so far\nSTOP 2\n";
    return 2;
}

}  // namespace jt9_cli
