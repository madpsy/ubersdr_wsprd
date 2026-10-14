// jt9's file mode for FT2 (ft2_jt9.h): lib/jt9.f90's WAV loop with
// improved 3.2.0's mode 52 (npts = 21*1728 in blocks of 3456: 10 blocks)
// and its one multimode_decoder call, jt9_params_init.f90's params,
// decoder.f90's decoded.txt and <DecodeFinished>. The decode is
// ft2::Decoder's, unchanged; the rest is FT4's file mode (ft4_jt9.cpp)
// with FT2's reading.
#include "ft2_jt9.h"

#include <algorithm>
#include <cstdint>

#include "fst4w_jt9.h"  // jt9's WAV reader (fst4w::Wav), jt9_nutc
#include "ft2.h"

namespace ft2 {

int run_files(const ft8::FileMode& fm, const std::vector<std::string>& files, const std::string& cwd,
              std::ostream& out, std::ostream& err) {
    // jt9.f90: "b" for none.
    std::string mycall = fm.mycall, hiscall = fm.hiscall;
    if (mycall == "b") mycall.clear();
    if (hiscall == "b") hiscall.clear();
    // init_default_params for mode 52 (mode 5's): what the FT2 branch and
    // ft4_decode read; the decoder halves nfa, nfb and nfqso.
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
    const std::string decoded_path = ft8::jt9_join(ft8::jt9_at(cwd, fm.temp_dir), "decoded.txt");

    Decoder dec;
    // mode 52: npts = 21*1728, kstep 3456, nblocks = npts/kstep = 10 (the
    // rest of npts never read).
    constexpr long KSTEP = 6912 / 2, NPTS = 21 * 1728, NBLOCKS = NPTS / KSTEP;
    std::vector<int16_t> id2(static_cast<size_t>(NBLOCKS * KSTEP));
    for (const std::string& infile : files) {
        fst4w::Wav wav;
        std::string why;
        if (!wav.open(ft8::jt9_at(cwd, infile), why)) {
            err << fm.program << ": cannot read WAV file " << infile << ": " << why << "\nSTOP 2\n";
            return 2;
        }
        if (wav.rate != 12000) {
            err << fm.program << ": unsupported " << wav.rate << " Hz WAV sample rate for the selected mode\nSTOP 2\n";
            return 2;
        }
        const int nutc = fst4w::jt9_nutc(infile);
        std::fill(id2.begin(), id2.end(), 0);
        for (long iblk = 1; iblk <= NBLOCKS; ++iblk) {
            const long k0 = (iblk - 1) * KSTEP, k = std::min(iblk * KSTEP, NPTS);
            const long got = wav.read(&id2[static_cast<size_t>(k0)], k - k0);
            if (got == 0 || got < k - k0) {
                out << " EOF on input file " << infile << "\n";
                break;
            }
        }

        const Decoder::Period per = dec.period(id2.data(), id2.size(), nutc, p, id2.size());
        std::string text = per.out;
        if (fm.quiet) text.erase(text.rfind("<DecodeFinished>"));
        out << text;
        ft8::jt9_write_decoded(dec.decoded_exists(), dec.decoded_text(), decoded_path);
    }
    return 0;
}

}  // namespace ft2
