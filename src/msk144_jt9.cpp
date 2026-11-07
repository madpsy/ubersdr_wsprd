// jt9's file mode for MSK144 (msk144_jt9.h): lib/jt9.f90's WAV loop for mode
// 144, jt9_params_init.f90's params, decode_msk144.f90 (src/msk144.h),
// decoder.f90's decoded.txt and <DecodeFinished>, at 567ad29ce.
#include "msk144_jt9.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <string>

#include "fst4w_jt9.h"  // jt9's WAV reader (fst4w::Wav), jt9_nutc
#include "msk144.h"
#include "msk144_decode.h"

namespace msk144 {

namespace {

// What jt9.f90 and init_default_params hand to decode_msk144.
struct Params {
    int ntr, ntol, nfqso, ndepth;
    std::string mycall, hiscall;
    bool shmsg, swl;  // ours: --msk-shmsg, --msk-swl (bshmsg, bswl; false in WSJT-X's jt9)
};

Params params_of(const ft8::FileMode& fm) {
    // jt9.f90: "b" for none (hiscall's takes hisgrid with it).
    std::string mycall = fm.mycall, hiscall = fm.hiscall;
    if (mycall == "b") mycall.clear();
    if (hiscall == "b") hiscall.clear();
    Params p;
    p.ntr = static_cast<int>(fm.trperiod);  // params%ntr = int(TRperiod)
    // jt9.f90:359-366 and init_default_params: ntol = min(ntol, 1000) for the modes without a rule of their own
    p.ntol = std::min(fm.ntol, 1000);
    p.nfqso = fm.nrxfreq;
    p.ndepth = fm.ndepth;
    p.mycall = mycall;
    p.hiscall = hiscall;
    p.shmsg = fm.msk_shmsg;
    p.swl = fm.msk_swl;
    return p;
}

// One file as jt9 -k decodes it; false, with why, for one it cannot read or
// of another rate. n: the decode lines printed.
bool decode_file(RtdState& st, const Params& p, const ft8::FileMode& fm, const std::string& infile, int nutc,
                 const std::string& cwd, const std::string& decoded_path, std::ostream& out, std::string& why,
                 int& n) {
    n = 0;
    fst4w::Wav wav;
    if (!wav.open(ft8::jt9_at(cwd, infile), why, true)) {
        why = "cannot read WAV file " + infile + ": " + why;
        return false;
    }
    if (wav.rate != 12000) {
        why = "unsupported " + std::to_string(wav.rate) + " Hz WAV sample rate for the selected mode";
        return false;
    }
    // jt9.f90:394-447: nsps 6912, kstep 3456, npts = TRperiod*12000, nblocks = npts/kstep
    constexpr long KSTEP = 6912 / 2;
    const long npts = static_cast<long>(fm.trperiod * 12000.0);
    const long nblocks = npts / KSTEP;
    // shared_data%id2=0: what is not read stays zero. A float file (ours) is read into id2f, each
    // sample x as x*32768.0f (Wav::read, exact): the 16-bit counts mskrtd's float(id2) would hold.
    std::vector<int16_t> id2;
    std::vector<float> id2f;
    long nhsym0 = -999;
    for (long iblk = 1; iblk <= nblocks; ++iblk) {
        const long k0 = (iblk - 1) * KSTEP;
        const long k = std::min(iblk * KSTEP, npts);
        long got;
        if (wav.is_float) {
            id2f.resize(static_cast<size_t>(k), 0.0f);
            got = wav.read(&id2f[static_cast<size_t>(k0)], k - k0);
        } else {
            id2.resize(static_cast<size_t>(k), 0);
            got = wav.read(&id2[static_cast<size_t>(k0)], k - k0);
        }
        if (got == 0) {
            out << " EOF on input file " << infile << "\n";
            break;
        }
        const long nhsym = (k - 2048) / KSTEP;
        if (nhsym >= 1 && nhsym != nhsym0) {
            nhsym0 = nhsym;
            if (nhsym >= 181) break;  // mode 144 is none of 240, 241, 242, 66
        }
        if (got < k - k0) {
            out << " EOF on input file " << infile << "\n";
            break;
        }
    }

    // jt9.f90's mode 144: decode_msk144, from a fresh process's state
    Packjt77& pk = msk144_packjt77();
    st.reset();
    msk144_reset_packjt77();
    if (wav.is_float)
        n = decode_msk144(st, pk, id2f.data(), static_cast<long>(id2f.size()), nutc, p.ntr, p.ntol, p.nfqso, p.ndepth,
                          p.mycall, p.hiscall, p.shmsg, p.swl, out);
    else
        n = decode_msk144(st, pk, id2.data(), static_cast<long>(id2.size()), nutc, p.ntr, p.ntol, p.nfqso, p.ndepth,
                          p.mycall, p.hiscall, p.shmsg, p.swl, out);

    // The "normal decoding pass": multimode_decoder, which for MSK144 opens decoded.txt once the rms of
    // id2(1:180000) reaches 0.5 and says <DecodeFinished> 0 0 0 (decoder.f90:104-107)
    float sum = 0.0f;
    for (long i = 0; i < 180000; ++i) {
        const long have = wav.is_float ? static_cast<long>(id2f.size()) : static_cast<long>(id2.size());
        const float x = i >= have ? 0.0f
                        : wav.is_float ? id2f[static_cast<size_t>(i)] : static_cast<float>(id2[static_cast<size_t>(i)]);
        sum = sum + x * x;
    }
    const float rms = std::sqrt(sum / 180000.0f);
    if (!(rms < 0.5f)) ft8::jt9_write_decoded(true, std::string(), decoded_path);
    if (!fm.quiet) out << "<DecodeFinished>" << fst4w::fortran_i(0, 4) << fst4w::fortran_i(0, 4) << fst4w::fortran_i(0, 9) << "\n";
    return true;
}

}  // namespace

int run_files(const ft8::FileMode& fm, const std::vector<std::string>& files, const std::string& cwd,
              std::ostream& out, std::ostream& err) {
    const Params p = params_of(fm);
    const std::string decoded_path = ft8::jt9_join(ft8::jt9_at(cwd, fm.temp_dir), "decoded.txt");
    RtdState st;
    for (const std::string& infile : files) {
        std::string why;
        int n = 0;
        if (!decode_file(st, p, fm, infile, fst4w::jt9_nutc(infile), cwd, decoded_path, out, why, n)) {
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
    RtdState st;
    std::string line, path, why;
    int nutc = 0;
    while (std::getline(in, line)) {
        if (!ft8::feed_line(line, path, nutc)) continue;
        int n = 0;
        if (!decode_file(st, p, fm, path, nutc, cwd, decoded_path, out, why, n)) {
            err << fm.program << ": " << why << "\n";
            err.flush();
            n = 0;
        }
        out << "<FileDone> " << path << " " << n << "\n";
        out.flush();
    }
    return 0;
}

}  // namespace msk144
