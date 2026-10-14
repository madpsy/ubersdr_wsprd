// jt9 -s (jt9_shm.h): lib/jt9a.f90 at WSJT-X v2.7.0, its loop exactly,
// and decoder.f90's reading of the params for FT8, FT4 and FT2 (at
// 567ad29ce, whose decoders ours are).
#include "jt9_shm.h"

#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <set>
#include <thread>
#include <vector>

#ifndef _WIN32
#include <unistd.h>
#endif

#include "fst4w_jt9.h"  // fortran_i
#include "ft2.h"
#include "ft4.h"
#include "ft8.h"
#include "ft8_jt9.h"  // jt9_join, jt9_at, jt9_write_decoded
#include "jt9_commons.h"
#include "qt_shm.h"

namespace jt9_shm {

namespace {

using skimmer::jt9::dec_data_t;
using Params = decltype(dec_data_t::params);

constexpr int MSDELAY = 30;     // jt9a's msdelay, ms
constexpr size_t NPTS = 180000;  // the id2 decoder.f90 is given (jt9a: shared_data%id2)

// decoder.f90: mycall=transfer(params%mycall,mycall) and the like: the
// bytes as they are (a C caller's NULs too).
std::string chars(const char* p, size_t n) { return std::string(p, n); }

// decoder.f90's reading of the params for FT8, FT4 and FT2 (ft4::Params
// and ft2::Params are ft8::Params).
ft8::Params params_of(const Params& q) {
    ft8::Params p;
    p.ndepth = q.ndepth;
    p.nfa = q.nfa;
    p.nfb = q.nfb;
    p.nfqso = q.nfqso;
    p.nftx = q.nftx;
    p.napwid = q.napwid;
    p.nqso_progress = q.nQSOProgress;
    p.nexp_decode = q.nexp_decode;
    p.lft8apon = q.lft8apon;
    p.lapcqonly = q.lapcqonly;
    p.nagain = q.nagain;
    p.mycall = chars(q.mycall, sizeof q.mycall);
    p.hiscall = chars(q.hiscall, sizeof q.hiscall);
    p.hisgrid = chars(q.hisgrid, sizeof q.hisgrid);
    return p;
}

// decode_completion's ('<DecodeFinished>',2i4,i9): nsynced and navg 0 for
// these modes.
std::string finished(int ndecoded) {
    return "<DecodeFinished>" + fst4w::fortran_i(0, 4) + fst4w::fortran_i(ndecoded, 4) + fst4w::fortran_i(0, 9) + "\n";
}

// What this jt9 does not decode as jt9 would: why, or "".
std::string refused(const Params& q) {
    if (q.nmode != 8 && q.nmode != 5 && q.nmode != 52)
        return "mode " + std::to_string(q.nmode) + " (this jt9 decodes FT8, FT4 and FT2 from shared memory)";
    if (q.nmode == 8 && q.lmultift8) return "FT8's multithreaded decoder (lmultift8)";
    const int ncontest = q.nexp_decode & 7;
    if (q.nmode == 8 && ncontest == 6) return "FT8's Fox mode (ncontest 6)";
    if (q.nmode == 8 && ncontest == 7 && q.b_superfox && q.b_even_seq) return "FT8's SuperFox decoder";
    if (q.nmode == 8 && q.emedelay != 0) return "FT8 with an EME delay (emedelay)";
    return "";
}

}  // namespace

int run(const Options& o, const std::string& cwd, std::ostream& out, std::ostream& err) {
    skimmer::qtcompat::SharedMemory shm(o.key);
    if (!shm.attach()) {
        // jt9a: "if(.not.ok) call abort".
        err << o.program << ": cannot attach to shared memory '" << o.key << "': " << shm.errorString() << "\n";
        err.flush();
        std::abort();
    }
    auto* dd = static_cast<dec_data_t*>(shm.data());
    const std::string decoded_path = ft8::jt9_join(ft8::jt9_at(cwd, o.temp_dir), "decoded.txt");
    std::unique_ptr<ft8::Decoder> ft8dec;
    std::unique_ptr<ft4::Decoder> ft4dec;
    std::unique_ptr<ft2::Decoder> ft2dec;
    std::set<std::string> said;  // each refusal once
    std::vector<int16_t> id2(NPTS);
    const auto nap = [] { std::this_thread::sleep_for(std::chrono::milliseconds(MSDELAY)); };
#ifndef _WIN32
    // Ours: a caller that dies without saying 999 leaves its memory
    // attached to us, and jt9a would wait for it for ever.
    const pid_t parent = getppid();
    const auto orphaned = [parent] { return getppid() != parent; };
#else
    const auto orphaned = [] { return false; };
#endif

    for (;;) {
        // 10: wait for ipc(2) (ipc[1]) to be 1, or 999.
        if (!shm.lock()) std::abort();
        if (dd->ipc[1] == 999) {
            shm.unlock();
            shm.detach();
            return 0;
        }
        if (dd->ipc[1] != 1) {
            if (!shm.unlock()) std::abort();
            if (orphaned()) return 0;
            nap();
            continue;
        }
        dd->ipc[1] = 0;
        if (shm.size() < sizeof(dec_data_t)) {
            shm.unlock();
            shm.detach();
            err << " jt9a: Incompatible shared-memory layout.\n";
            return 2;
        }
        const Params q = dd->params;  // a copy: the caller carries on writing to it
        if (!shm.unlock()) std::abort();
        // jt9a reads shared_data%id2 as it decodes, unlocked: copied here.
        std::memcpy(id2.data(), dd->d2, NPTS * sizeof(int16_t));

        std::string text;
        const std::string why = refused(q);
        if (!why.empty()) {
            if (said.insert(why).second) err << o.program << ": not decoded here: " << why << "\n";
            text = finished(0);
        } else {
            const ft8::Params p = params_of(q);
            if (q.nmode == 8) {
                if (!ft8dec) ft8dec = std::make_unique<ft8::Decoder>();
                if (q.ndiskdat && !q.nagain) {
                    // The early passes at nzhsym 41 and 47, then 50.
                    text = ft8dec->period(id2.data(), NPTS, q.nutc, p).out;
                } else {
                    // One call at the caller's nzhsym (WSJT-X's live
                    // decodes, or a decode again); <DecodeFinished> after
                    // each, live, else only at 50.
                    const ft8::Decoder::Call c = ft8dec->call(id2.data(), q.nzhsym, q.nutc, p);
                    for (const auto& d : c.decodes) text += d.line + "\n";
                    if (q.nzhsym == 50 || !q.ndiskdat) text += finished(c.ndecoded);
                }
                ft8::jt9_write_decoded(ft8dec->decoded_exists(), ft8dec->decoded_text(), decoded_path);
            } else if (q.nmode == 5) {
                if (!ft4dec) ft4dec = std::make_unique<ft4::Decoder>();
                text = ft4dec->period(id2.data(), NPTS, q.nutc, p, NPTS).out;
                ft8::jt9_write_decoded(ft4dec->decoded_exists(), ft4dec->decoded_text(), decoded_path);
            } else {
                if (!ft2dec) ft2dec = std::make_unique<ft2::Decoder>();
                text = ft2dec->period(id2.data(), NPTS, q.nutc, p, NPTS).out;
                ft8::jt9_write_decoded(ft2dec->decoded_exists(), ft2dec->decoded_text(), decoded_path);
            }
        }
        out << text;
        out.flush();

        // 100: wait for ipc(3) (ipc[2]) to be 1.
        for (;;) {
            if (!shm.lock()) std::abort();
            if (dd->ipc[2] == 1) {
                dd->ipc[2] = 0;
                if (!shm.unlock()) std::abort();
                break;
            }
            if (!shm.unlock()) std::abort();
            if (orphaned()) return 0;
            nap();
        }
    }
}

}  // namespace jt9_shm
