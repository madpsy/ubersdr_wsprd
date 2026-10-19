// jt9's file mode for FT8 (ft8_jt9.h): lib/jt9.f90's WAV loop and its
// FT8 branch (the 41/47/50 calls on id2a), jt9_params_init.f90's params,
// decoder.f90's decoded.txt and <DecodeFinished>, at 567ad29ce. The decode
// is ft8::Decoder's, unchanged.
#include "ft8_jt9.h"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <string>

#include "fst4w_jt9.h"  // jt9's WAV reader (fst4w::Wav), jt9_nutc
#include "ft8.h"

namespace ft8 {

std::string jt9_join(const std::string& dir, const std::string& name) {
    if (dir.empty() || dir == ".") return name;
    return dir.back() == '/' ? dir + name : dir + "/" + name;
}
std::string jt9_at(const std::string& cwd, const std::string& path) {
    // Absolute: /..., and on Windows also C:\..., \\server\... and \...
    if (!path.empty() && (path[0] == '/' || path[0] == '\\')) return path;
    if (path.size() > 1 && path[1] == ':') return path;
    return jt9_join(cwd, path);
}
void jt9_write_decoded(bool exists, const std::string& text, const std::string& path) {
    if (!exists) return;
    if (!text.empty())
        std::ofstream(path, std::ios::trunc) << text;
    else if (!std::ifstream(path))
        std::ofstream(path, std::ios::app);
}

namespace {

// init_default_params for mode 8, and the WAV path's overrides. nftx is
// not set there: zero (the allocated shared_data is fresh memory).
Params params_of(const FileMode& fm) {
    // jt9.f90: "b" for none (hiscall's takes hisgrid with it).
    std::string mycall = fm.mycall, hiscall = fm.hiscall, hisgrid = fm.hisgrid;
    if (mycall == "b") mycall.clear();
    if (hiscall == "b") {
        hiscall.clear();
        hisgrid.clear();
    }
    Params p;
    p.ndepth = fm.ndepth;
    p.nfa = fm.flow;
    p.nfb = fm.fhigh;
    p.nfqso = fm.nrxfreq;
    p.nftx = 0;
    p.napwid = 75;
    p.nqso_progress = fm.nqso_progress;
    p.nexp_decode = fm.nexp_decode;
    p.lft8apon = true;
    p.lapcqonly = false;
    p.nagain = false;
    p.mycall = mycall;
    p.hiscall = hiscall;
    p.hisgrid = hisgrid;
    p.fast_subtract = fm.fast_subtract;
    return p;
}

// One file as jt9 -8 decodes it: read (a 32-bit float file too, ours: its
// samples in 16-bit units, ft8.h's float call()), the three calls, the
// lines printed, decoded.txt. false, with why, for a file it cannot read or
// of another rate (nothing printed).
bool decode_file(Decoder& dec, const Params& p, const FileMode& fm, const std::string& infile, int nutc,
                 const std::string& cwd, const std::string& decoded_path, std::ostream& out, std::string& why,
                 int& ndecoded) {
    fst4w::Wav wav;
    if (!wav.open(jt9_at(cwd, infile), why, true)) {
        why = "cannot read WAV file " + infile + ": " + why;
        return false;
    }
    if (wav.rate != 12000) {
        why = "unsupported " + std::to_string(wav.rate) + " Hz WAV sample rate for the selected mode";
        return false;
    }
    constexpr long KSTEP = 6912 / 2, LAST = 52 * KSTEP;  // FT8 reads no block ending past 179,712
    std::vector<int16_t> id2;
    std::vector<float> id2f;
    if (wav.is_float)
        id2f.assign(15 * 12000, 0.0f);
    else
        id2.assign(15 * 12000, 0);
    const long npts = static_cast<long>(fm.trperiod * 12000.0);
    const long nblocks = npts / KSTEP;
    for (long iblk = 1; iblk <= nblocks; ++iblk) {
        const long k0 = (iblk - 1) * KSTEP;
        const long k = std::min(iblk * KSTEP, npts);
        if (k > LAST) break;
        const long got = wav.is_float ? wav.read(&id2f[static_cast<size_t>(k0)], k - k0)
                                      : wav.read(&id2[static_cast<size_t>(k0)], k - k0);
        if (got == 0 || got < k - k0) {
            out << " EOF on input file " << infile << "\n";
            break;
        }
    }

    // The three calls: 41 and 47 on their first blocks, 50 on 50 blocks.
    const Decoder::Period per = wav.is_float ? dec.period(id2f.data(), id2f.size(), nutc, p, 50 * KSTEP)
                                             : dec.period(id2.data(), id2.size(), nutc, p, 50 * KSTEP);
    std::string text = per.out;
    if (fm.quiet) text.erase(text.rfind("<DecodeFinished>"));
    out << text;
    ndecoded = per.ndecoded;

    // decoded.txt: opened (status 'unknown') by every call past the RMS
    // gate, so created then; its lines those of the last call that wrote
    // any (a call writing none leaves it as it was, also as it was
    // before the run).
    jt9_write_decoded(dec.decoded_exists(), dec.decoded_text(), decoded_path);
    return true;
}

}  // namespace

int run_files(const FileMode& fm, const std::vector<std::string>& files, const std::string& cwd, std::ostream& out,
              std::ostream& err) {
    const Params p = params_of(fm);
    const std::string decoded_path = jt9_join(jt9_at(cwd, fm.temp_dir), "decoded.txt");
    Decoder dec;
    for (const std::string& infile : files) {
        std::string why;
        int n = 0;
        if (!decode_file(dec, p, fm, infile, fst4w::jt9_nutc(infile), cwd, decoded_path, out, why, n)) {
            err << fm.program << ": " << why << "\nSTOP 2\n";
            return 2;
        }
    }
    return 0;
}

int pcmrecord_nutc(const std::string& path) {
    const std::string name = path.substr(path.find_last_of('/') + 1);
    // YYYYMMDDTHHMMSSZ_...
    if (name.size() < 17 || name[8] != 'T' || name[15] != 'Z' || name[16] != '_') return -1;
    for (int i : {0, 1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 12, 13, 14})
        if (name[static_cast<size_t>(i)] < '0' || name[static_cast<size_t>(i)] > '9') return -1;
    return std::atoi(name.substr(9, 6).c_str());
}

bool feed_line(const std::string& line, std::string& path, int& nutc) {
    std::string l = line;
    while (!l.empty() && (l.back() == '\r' || l.back() == ' ' || l.back() == '\t')) l.pop_back();
    const size_t b = l.find_first_not_of(" \t");
    if (b == std::string::npos) return false;
    l = l.substr(b);
    // "PATH HHMMSS": the last word, when six digits after a blank.
    const size_t sp = l.find_last_of(" \t");
    if (sp != std::string::npos && l.size() - sp - 1 == 6 &&
        l.find_first_not_of("0123456789", sp + 1) == std::string::npos) {
        path = l.substr(0, l.find_last_not_of(" \t", sp) + 1);
        nutc = std::atoi(l.substr(sp + 1).c_str());
        return true;
    }
    path = l;
    nutc = fst4w::jt9_nutc(path);
    if (nutc == 0) {
        const int t = pcmrecord_nutc(path);
        if (t >= 0) nutc = t;
    }
    return true;
}

int run_feed(const FileMode& fm, std::istream& in, const std::string& cwd, std::ostream& out, std::ostream& err) {
    const Params p = params_of(fm);
    const std::string decoded_path = jt9_join(jt9_at(cwd, fm.temp_dir), "decoded.txt");
    Decoder dec;
    std::string line, path, why;
    int nutc = 0;
    while (std::getline(in, line)) {
        if (!feed_line(line, path, nutc)) continue;
        int n = 0;
        if (!decode_file(dec, p, fm, path, nutc, cwd, decoded_path, out, why, n)) {
            err << fm.program << ": " << why << "\n";
            err.flush();
            n = 0;
        }
        out << "<FileDone> " << path << " " << n << "\n";
        out.flush();
    }
    return 0;
}

}  // namespace ft8
