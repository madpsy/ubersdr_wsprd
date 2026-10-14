// jt9's file mode for FT8 (ft8_jt9.h): lib/jt9.f90's WAV loop and its
// FT8 branch (the 41/47/50 calls on id2a), jt9_params_init.f90's params,
// decoder.f90's decoded.txt and <DecodeFinished>, at 567ad29ce. The decode
// is ft8::Decoder's, unchanged.
#include "ft8_jt9.h"

#include <algorithm>
#include <cstdint>
#include <fstream>

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

int run_files(const FileMode& fm, const std::vector<std::string>& files, const std::string& cwd, std::ostream& out,
              std::ostream& err) {
    // jt9.f90: "b" for none (hiscall's takes hisgrid with it).
    std::string mycall = fm.mycall, hiscall = fm.hiscall, hisgrid = fm.hisgrid;
    if (mycall == "b") mycall.clear();
    if (hiscall == "b") {
        hiscall.clear();
        hisgrid.clear();
    }
    // init_default_params for mode 8, and the WAV path's overrides. nftx is
    // not set there: zero (the allocated shared_data is fresh memory).
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
    const std::string decoded_path = jt9_join(jt9_at(cwd, fm.temp_dir), "decoded.txt");

    Decoder dec;
    constexpr long KSTEP = 6912 / 2, LAST = 52 * KSTEP;  // FT8 reads no block ending past 179,712
    std::vector<int16_t> id2(15 * 12000);
    for (const std::string& infile : files) {
        fst4w::Wav wav;
        std::string why;
        if (!wav.open(jt9_at(cwd, infile), why)) {
            err << fm.program << ": cannot read WAV file " << infile << ": " << why << "\nSTOP 2\n";
            return 2;
        }
        if (wav.rate != 12000) {
            err << fm.program << ": unsupported " << wav.rate << " Hz WAV sample rate for the selected mode\nSTOP 2\n";
            return 2;
        }
        const int nutc = fst4w::jt9_nutc(infile);
        const long npts = static_cast<long>(fm.trperiod * 12000.0);
        std::fill(id2.begin(), id2.end(), 0);
        const long nblocks = npts / KSTEP;
        for (long iblk = 1; iblk <= nblocks; ++iblk) {
            const long k0 = (iblk - 1) * KSTEP;
            const long k = std::min(iblk * KSTEP, npts);
            if (k > LAST) break;
            const long got = wav.read(&id2[static_cast<size_t>(k0)], k - k0);
            if (got == 0 || got < k - k0) {
                out << " EOF on input file " << infile << "\n";
                break;
            }
        }

        // The three calls: 41 and 47 on their first blocks, 50 on 50 blocks.
        const Decoder::Period per = dec.period(id2.data(), id2.size(), nutc, p, 50 * KSTEP);
        std::string text = per.out;
        if (fm.quiet) text.erase(text.rfind("<DecodeFinished>"));
        out << text;

        // decoded.txt: opened (status 'unknown') by every call past the RMS
        // gate, so created then; its lines those of the last call that wrote
        // any (a call writing none leaves it as it was, also as it was
        // before the run).
        jt9_write_decoded(dec.decoded_exists(), dec.decoded_text(), decoded_path);
    }
    return 0;
}

}  // namespace ft8
