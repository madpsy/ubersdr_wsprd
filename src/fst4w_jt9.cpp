// jt9's file mode for FST4W (fst4w_jt9.h): lib/jt9.f90's WAV loop,
// decoder.f90's multimode_decoder_core (its RMS check, decoded.txt, the
// FST4W call, <DecodeFinished>), decoder_callbacks.f90's fst4_decoded and
// fst4_decode.f90's file writes, at 567ad29ce.
#include "fst4w_jt9.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <memory>
#include <sstream>
#include <stdexcept>

#include "fst4w.h"

namespace fst4w {

namespace {

constexpr int NTMAX = 30 * 60;  // jt9's id2: NTMAX*12000 samples

std::string join(const std::string& dir, const std::string& name) {
    if (dir.empty() || dir == ".") return name;
    return dir.back() == '/' ? dir + name : dir + "/" + name;
}
std::string at(const std::string& cwd, const std::string& path) {
    if (!path.empty() && path[0] == '/') return path;
    return join(cwd, path);
}
bool exists(const std::string& path) { return static_cast<bool>(std::ifstream(path)); }

// parse_wav_filename_nutc: the last 6 characters of the name before
// .wav/.WAV if all digits, else '_' and the last 4, else 0.
int nutc_of(const std::string& filename) {
    std::string f = filename;
    while (!f.empty() && f.back() == ' ') f.pop_back();
    const size_t dot = f.rfind('.');
    if (dot == std::string::npos) return 0;
    const std::string ext = f.substr(dot);
    if (ext != ".wav" && ext != ".WAV") return 0;
    const std::string stem = f.substr(0, dot);
    auto digits = [](const std::string& s) {
        if (s.empty()) return false;
        for (char c : s)
            if (c < '0' || c > '9') return false;
        return true;
    };
    if (stem.size() >= 6 && digits(stem.substr(stem.size() - 6))) return std::atoi(stem.substr(stem.size() - 6).c_str());
    if (stem.size() >= 5 && stem[stem.size() - 5] == '_' && digits(stem.substr(stem.size() - 4)))
        return std::atoi(stem.substr(stem.size() - 4).c_str());
    return 0;
}

// lib/readwav.f90: a RIFF WAV's header read and checked as it reads it
// (its checks in its order, its messages), then its 16-bit samples from the
// data chunk on demand.
struct Wav {
    std::ifstream f;
    int rate = 0, channels = 0, bits = 0;
    long remaining = 0;  // bytes of data left
    bool open(const std::string& path, std::string& why) {
        f.open(path, std::ios::binary);
        if (!f) {
            why = "cannot open file";
            return false;
        }
        f.seekg(0, std::ios::end);
        const int64_t filesize = static_cast<int64_t>(f.tellg());
        f.seekg(0);
        auto fail = [&](const char* text) {
            why = text;
            f.close();
            return false;
        };
        if (filesize < 12) return fail("file is too short for a RIFF/WAVE header");
        char riff[12];
        if (!f.read(riff, 12)) return fail("cannot read RIFF/WAVE header");
        if (std::memcmp(riff, "RIFF", 4) != 0 || std::memcmp(riff + 8, "WAVE", 4) != 0)
            return fail("not a RIFF/WAVE file");
        int64_t filepos = 13;  // Fortran's, 1-based
        bool format_found = false;
        for (;;) {
            if (filepos + 8 > filesize + 1) return fail("missing data chunk");
            char id[4];
            int32_t size;
            f.seekg(filepos - 1);
            if (!f.read(id, 4) || !f.read(reinterpret_cast<char*>(&size), 4)) return fail("cannot read WAV chunk header");
            filepos += 8;
            if (size < 0) return fail("WAV chunk has an invalid size");
            if (filepos + size > filesize + 1) return fail("WAV chunk extends beyond end of file");
            const int64_t nextpos = filepos + size + (size % 2);
            if (std::memcmp(id, "fmt ", 4) == 0) {
                if (size < 16) return fail("fmt chunk is too short");
                char b[16];
                if (!f.read(b, 16)) return fail("cannot read fmt chunk");
                int16_t fmt_tag, ch, bps;
                int32_t sr;
                std::memcpy(&fmt_tag, &b[0], 2);
                std::memcpy(&ch, &b[2], 2);
                std::memcpy(&sr, &b[4], 4);
                std::memcpy(&bps, &b[14], 2);
                rate = sr;
                channels = ch;
                bits = bps;
                if (fmt_tag != 1 || ch != 1 || bps != 16) return fail("unsupported WAV format; expected mono 16-bit PCM");
                if (sr != 11025 && sr != 12000) return fail("unsupported WAV sample rate");
                format_found = true;
            } else if (std::memcmp(id, "data", 4) == 0) {
                if (!format_found) return fail("data chunk appears before fmt chunk");
                remaining = size;
                f.seekg(filepos - 1);
                return true;
            }
            if (nextpos > filesize + 1) return fail("WAV chunk extends beyond end of file");
            filepos = nextpos;
        }
    }
    // Up to n samples; how many were read.
    long read(int16_t* out, long n) {
        const long want = std::min(n, remaining / 2);
        f.read(reinterpret_cast<char*>(out), want * 2);
        const long got = static_cast<long>(f.gcount()) / 2;
        remaining -= got * 2;
        return got;
    }
};

}  // namespace

std::string fortran_i(long v, int w, int m) {
    char b[64];
    if (m > 0) std::snprintf(b, sizeof b, "%0*ld", m + (v < 0 ? 1 : 0), v);
    else std::snprintf(b, sizeof b, "%ld", v);
    std::string s(b);
    if (static_cast<int>(s.size()) > w) return std::string(static_cast<size_t>(w), '*');
    return std::string(static_cast<size_t>(w) - s.size(), ' ') + s;
}

std::string fortran_f(float x, int w, int d) {
    if (std::isnan(x)) {
        const std::string s = "NaN";
        return static_cast<int>(s.size()) > w ? std::string(static_cast<size_t>(w), '*')
                                               : std::string(static_cast<size_t>(w) - s.size(), ' ') + s;
    }
    if (std::isinf(x)) {
        const std::string s = x < 0 ? (w >= 9 ? "-Infinity" : "-Inf") : (w >= 8 ? "Infinity" : "Inf");
        return static_cast<int>(s.size()) > w ? std::string(static_cast<size_t>(w), '*')
                                               : std::string(static_cast<size_t>(w) - s.size(), ' ') + s;
    }
    char b[128];
    std::snprintf(b, sizeof b, "%.*f", d, static_cast<double>(x));
    std::string s(b);
    if (d == 0) s += '.';
    if (static_cast<int>(s.size()) > w) {
        // The optional leading zero goes first: "0.5" -> ".5", "-0.5" -> "-.5".
        if (s.rfind("0.", 0) == 0) s.erase(0, 1);
        else if (s.rfind("-0.", 0) == 0) s.erase(1, 1);
    }
    if (static_cast<int>(s.size()) > w) return std::string(static_cast<size_t>(w), '*');
    return std::string(static_cast<size_t>(w) - s.size(), ' ') + s;
}

bool jt9_samples(const int16_t* samples, size_t n, int trperiod, std::vector<int16_t>& id2) {
    const size_t npts = static_cast<size_t>(trperiod) * 12000;
    const size_t kstep = 6912 / 2;
    const size_t take = std::min(n, npts / kstep * kstep);
    id2.assign(npts, 0);
    std::copy(samples, samples + take, id2.begin());
    float sum = 0.0f;
    for (size_t i = 0; i < 180000 && i < npts; ++i) {
        const float x = static_cast<float>(id2[i]);
        sum = sum + x * x;
    }
    return !(std::sqrt(sum / 180000.0f) < 0.5f);
}

std::string decodes_dat_line(const Result& r) {
    std::string m = r.msg;
    while (!m.empty() && m.back() == ' ') m.pop_back();
    return fortran_i(r.nutc, 6, 6) + fortran_i(r.icand, 4) + fortran_i(r.itry, 3) + fortran_i(r.nsyncoh, 3) +
           fortran_i(r.iaptype, 3) + fortran_i(r.ijitter, 3) + fortran_i(r.npct, 3) + fortran_i(r.ntype, 3) +
           fortran_i(r.keff, 4) + fortran_i(r.nsync_qual, 4) + fortran_i(r.nharderrors, 4) + fortran_f(r.dmin, 6, 1) +
           fortran_i(r.nhp, 4) + fortran_f(r.hd, 6, 1) + fortran_f(r.csync, 9, 2) + fortran_f(r.xsnr, 6, 1) +
           fortran_f(r.xdt, 6, 2) + fortran_f(r.fsig, 7, 1) + fortran_f(r.w50, 7, 3) + " " + m;
}

std::string calls_file(const std::vector<std::string>& calls) {
    std::string t;
    // write(42,'(a20)') trim(wcalls(i)): A20 right-justifies a shorter string.
    for (const std::string& w : calls) t += (w.size() < 20 ? std::string(20 - w.size(), ' ') + w : w.substr(0, 20)) + "\n";
    return t;
}

int run_files(const FileMode& fm, const std::vector<std::string>& files, const std::string& cwd, std::ostream& out,
              std::ostream& err) {
    std::string data_dir = fm.data_dir, temp_dir = fm.temp_dir;
    const double trperiod = fm.trperiod;
    const int mode = fm.hash22 ? 242 : 241, ndepth = fm.ndepth, nrxfreq = fm.nrxfreq;
    const bool quiet = fm.quiet;
    const Options::Spread spread = fm.spread;
    std::string mycall = fm.mycall, hiscall = fm.hiscall;
    if (mycall == "b") mycall.clear();
    if (hiscall == "b") hiscall.clear();
    const int ntol = std::min(fm.ntol, 100);
    const int nexp_decode = fm.nexp_decode_set ? fm.nexp_decode : 3 * 256;  // single decode off and nb=0

    data_dir = at(cwd, data_dir);
    temp_dir = at(cwd, temp_dir);
    Decoder dec;
    bool first = true;
    std::vector<int16_t> id2(static_cast<size_t>(NTMAX) * 12000);
    for (const std::string& infile : files) {
        Wav wav;
        std::string why;
        if (!wav.open(at(cwd, infile), why)) {
            err << fm.program << ": cannot read WAV file " << infile << ": " << why << "\nSTOP 2\n";
            return 2;
        }
        if (wav.rate != 12000) {
            err << fm.program << ": unsupported " << wav.rate << " Hz WAV sample rate for the selected mode\nSTOP 2\n";
            return 2;
        }
        const int nutc = nutc_of(infile);
        const long npts = static_cast<long>(trperiod * 12000.0);
        const long kstep = 6912 / 2;
        std::fill(id2.begin(), id2.end(), 0);
        const long nblocks = npts / kstep;
        for (long iblk = 1; iblk <= nblocks; ++iblk) {
            const long k0 = (iblk - 1) * kstep;
            const long k = std::min(iblk * kstep, npts);
            const long got = wav.read(&id2[static_cast<size_t>(k0)], k - k0);
            if (got == 0) {
                out << " EOF on input file " << infile << "\n";
                break;
            }
            if (got < k - k0) {
                out << " EOF on input file " << infile << "\n";
                break;
            }
        }

        // multimode_decoder_core: no decode when the first 15 s are quiet.
        int ndecoded = 0;
        float sum = 0.0f;
        for (int i = 0; i < 180000; ++i) {
            const float x = static_cast<float>(id2[static_cast<size_t>(i)]);
            sum = sum + x * x;
        }
        const float rms = std::sqrt(sum / 180000.0f);
        if (!(rms < 0.5f)) {
            // decoded.txt: open(status='unknown') creates it empty, or keeps
            // it as it was until a line is written, which makes it the
            // file's first.
            std::ostringstream decoded;
            if (first) {
                std::ifstream cf(join(data_dir, "fst4w_calls.txt"));
                if (cf) {
                    std::vector<std::string> lines;
                    for (std::string l; std::getline(cf, l);) lines.push_back(l);
                    dec.load_calls(lines);
                }
                first = false;
            }
            Options o;
            o.ntrperiod = static_cast<int>(trperiod);
            o.ndepth = ndepth;
            o.nfqso = nrxfreq;
            o.ntol = ntol;
            o.nexp_decode = nexp_decode;
            o.mycall = mycall;
            o.hiscall = hiscall;
            o.lprinthash22 = mode == 242;
            o.plotspec = exists(join(cwd, "plotspec"));
            o.decdata = exists(join(data_dir, "decdata"));
            o.spread = spread;
            bool calls_changed = false;
            // fst4_decode has its sizes for these periods alone; on any
            // other its decode divides by a size never set, and jt9 dies:
            // of SIGFPE on x86-64, of SIGSEGV on arm64 (where an integer
            // division by zero does not trap, and it goes on to read out of
            // bounds). Said as gfortran says it, the status a shell's for
            // that signal (stdout so far kept: ours is not lost).
            if (o.ntrperiod != 15 && o.ntrperiod != 30 && o.ntrperiod != 60 && o.ntrperiod != 120 &&
                o.ntrperiod != 300 && o.ntrperiod != 900 && o.ntrperiod != 1800) {
                if (!exists(join(temp_dir, "decoded.txt"))) std::ofstream(join(temp_dir, "decoded.txt"));  // opened before the decode
#if defined(__aarch64__)
                err << "\nProgram received signal SIGSEGV: Segmentation fault - invalid memory reference.\n";
                return 128 + 11;
#else
                err << "\nProgram received signal SIGFPE: Floating-point exception - erroneous arithmetic operation.\n";
                return 128 + 8;
#endif
            }
            try {
            dec.decode(id2.data(), nutc, o, [&](const Result& r) {
                if (o.spread == Options::Spread::Both && r.fmid != -999.0f) {
                    std::string m = r.msg;
                    while (!m.empty() && m.back() == ' ') m.pop_back();
                    char b[512];
                    std::snprintf(b, sizeof b, "SPREAD %s %d %.9g %.9g %.9g %.9g %s|%s|%s\n", infile.c_str(), r.nsnr,
                                  static_cast<double>(r.w50), static_cast<double>(r.w50_fast), static_cast<double>(r.fmid),
                                  static_cast<double>(r.fmid_fast), fortran_f(r.w50, 7, 3).c_str(),
                                  fortran_f(r.w50_fast, 7, 3).c_str(), m.c_str());
                    err << b;
                }
                if (r.spectrum.size()) {
                    std::ofstream f52(join(cwd, "fort.52"), std::ios::app);
                    for (const auto& p : r.spectrum) f52 << fortran_f(p.first, 12, 6) << fortran_f(p.second, 12, 6) << "\n";
                }
                if (o.decdata) {
                    std::ofstream dat(join(data_dir, "fst4_decodes.dat"), std::ios::app);
                    dat << decodes_dat_line(r) << "\n";
                }
                // fst4_decoded (ntrperiod >= 60).
                std::string decoded0 = r.msg;
                const std::string annot = "  ";
                // Under 60 s the time is HHMMSS (i6.6), else HHMM (i4.4).
                const int tw = r.ntrperiod < 60 ? 6 : 4;
                std::string line = fortran_i(r.nutc, tw, tw) + fortran_i(r.nsnr, 4) + fortran_f(r.xdt, 5, 1) +
                                   fortran_i(std::lround(r.fsig), 5) + " ` " + " " + decoded0 + " " + annot;
                line.resize(70, ' ');
                if (r.fmid != -999.0f) {
                    const std::string w = r.w50 < 0.95f ? fortran_f(r.w50, 6, 3) : fortran_f(r.w50, 6, 2);
                    line.replace(64, 6, w);
                }
                out << line << "\n";
                if (decoded)
                    decoded << fortran_i(r.nutc, tw, tw) << fortran_i(std::lround(r.sync), 4) << fortran_i(r.nsnr, 5)
                            << fortran_f(r.xdt, 6, 1) << fortran_f(r.fsig, 8, 0) << fortran_i(0, 4) << "   "
                            << decoded0 << " FST4\n";
                ++ndecoded;
            }, &calls_changed);
            } catch (const std::invalid_argument& e) {
                // Not jt9's: -c or -x not a standard call (fst4w_unpack.h).
                err << fm.program << ": " << e.what() << "\n";
                return 2;
            }
            if (!decoded.str().empty() || !exists(join(temp_dir, "decoded.txt"))) {
                std::ofstream f(join(temp_dir, "decoded.txt"), std::ios::trunc);  // (created empty if new)
                f << decoded.str();
            }
            if (calls_changed) {  // re-write the fst4w_calls.txt file
                std::ofstream cf(join(data_dir, "fst4w_calls.txt"), std::ios::trunc);
                cf << calls_file(dec.calls());
            }
        }
        if (!quiet) {
            char b[64];
            std::snprintf(b, sizeof b, "<DecodeFinished>%4d%4d%9d", 0, ndecoded, 0);
            out << b << "\n";
        }
    }
    return 0;
}

}  // namespace fst4w
