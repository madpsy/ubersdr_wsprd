// jt9's file mode for FT4 (ft4_jt9.h): lib/jt9.f90's WAV loop (npts =
// 21*3456 for mode 5) and its one multimode_decoder call,
// jt9_params_init.f90's params, decoder.f90's decoded.txt and
// <DecodeFinished>, at 567ad29ce. The decode is ft4::Decoder's, unchanged;
// the rest is FT8's file mode (ft8_jt9.cpp) with FT4's reading.
#include "ft4_jt9.h"

#include <algorithm>
#include <cstdint>
#include <string>

#include "fst4w_jt9.h"  // jt9's WAV reader (fst4w::Wav), jt9_nutc
#include "ft4.h"

namespace ft4 {

namespace {

// init_default_params for mode 5: what ft4_decode reads.
Params params_of(const ft8::FileMode& fm) {
    // jt9.f90: "b" for none.
    std::string mycall = fm.mycall, hiscall = fm.hiscall;
    if (mycall == "b") mycall.clear();
    if (hiscall == "b") hiscall.clear();
    Params p;
    p.ndepth = fm.ndepth;
    p.nfa = fm.flow;
    p.nfb = fm.fhigh;
    p.nfqso = fm.nrxfreq;
    p.nqso_progress = fm.nqso_progress;
    p.nexp_decode = fm.nexp_decode;
    p.lapcqonly = false;
    p.nagain = false;
    p.mycall = mycall;
    p.hiscall = hiscall;
    return p;
}

// One file as jt9 -5 decodes it (a 32-bit float file too, ours: ft4.h's
// float call()); false, with why, for one it cannot read or of another rate.
bool decode_file(Decoder& dec, const Params& p, const ft8::FileMode& fm, const std::string& infile, int nutc,
                 const std::string& cwd, const std::string& decoded_path, std::ostream& out, std::string& why,
                 int& ndecoded) {
    fst4w::Wav wav;
    if (!wav.open(ft8::jt9_at(cwd, infile), why, true)) {
        why = "cannot read WAV file " + infile + ": " + why;
        return false;
    }
    if (wav.rate != 12000) {
        why = "unsupported " + std::to_string(wav.rate) + " Hz WAV sample rate for the selected mode";
        return false;
    }
    constexpr long KSTEP = 6912 / 2, NPTS = 21 * KSTEP;  // mode 5: npts = 21*3456, nblocks 21
    std::vector<int16_t> id2;
    std::vector<float> id2f;
    if (wav.is_float)
        id2f.assign(static_cast<size_t>(NPTS), 0.0f);
    else
        id2.assign(static_cast<size_t>(NPTS), 0);
    for (long iblk = 1; iblk <= NPTS / KSTEP; ++iblk) {
        const long k0 = (iblk - 1) * KSTEP, k = iblk * KSTEP;
        const long got = wav.is_float ? wav.read(&id2f[static_cast<size_t>(k0)], k - k0)
                                      : wav.read(&id2[static_cast<size_t>(k0)], k - k0);
        if (got == 0 || got < k - k0) {
            out << " EOF on input file " << infile << "\n";
            break;
        }
    }

    const Decoder::Period per =
        wav.is_float ? dec.period(id2f.data(), id2f.size(), nutc, p, static_cast<size_t>(NPTS))
                     : dec.period(id2.data(), id2.size(), nutc, p, static_cast<size_t>(NPTS));
    std::string text = per.out;
    if (fm.quiet) text.erase(text.rfind("<DecodeFinished>"));
    out << text;
    ndecoded = per.ndecoded;
    ft8::jt9_write_decoded(dec.decoded_exists(), dec.decoded_text(), decoded_path);
    return true;
}

}  // namespace

int run_files(const ft8::FileMode& fm, const std::vector<std::string>& files, const std::string& cwd,
              std::ostream& out, std::ostream& err) {
    const Params p = params_of(fm);
    const std::string decoded_path = ft8::jt9_join(ft8::jt9_at(cwd, fm.temp_dir), "decoded.txt");
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

int run_feed(const ft8::FileMode& fm, std::istream& in, const std::string& cwd, std::ostream& out,
             std::ostream& err) {
    const Params p = params_of(fm);
    const std::string decoded_path = ft8::jt9_join(ft8::jt9_at(cwd, fm.temp_dir), "decoded.txt");
    Decoder dec;
    std::string line, path, why;
    int nutc = 0;
    while (std::getline(in, line)) {
        if (!ft8::feed_line(line, path, nutc)) continue;
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

}  // namespace ft4
