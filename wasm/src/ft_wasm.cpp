// The FT8, FT4, FT2 and MSK144 decoders for JavaScript (embind): one Decoder
// is one jt9 process for one dial, as ft8::Decoder, ft4::Decoder and
// ft2::Decoder are (their state carried from period to period: the a7
// tables, packjt77's hashes). A period's 12 kHz samples in, its decodes out,
// as WSJT-X's GUI has jt9 decode them (ft8.h's, ft4.h's, ft2.h's period()).
// MSK144 decodes each period afresh, as `jt9 -k` decodes each file
// (msk144_jt9.cpp: the samples jt9 reads of it, mskrtd's state and
// packjt77's tables new for each).
//
// Samples: a Float32Array (-1..1, scaled by 32768 as a float .wav file is
// read) or an Int16Array (as they are).
#include <emscripten/bind.h>
#include <emscripten/val.h>

#include <algorithm>
#include <cstdint>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "ft2.h"
#include "ft4.h"
#include "ft8.h"
#include "msk144.h"
#include "msk144_decode.h"

using emscripten::val;

namespace {

// ft8.h's Params (the GUI's: depth 2, 100-3200 Hz, rx/tx 1500, blank calls)
// with the options given.
ft8::Params params_from(const val& o) {
    ft8::Params p;
    if (o.isUndefined() || o.isNull()) return p;
    auto num = [&](const char* k, int& v) {
        if (o.hasOwnProperty(k)) v = o[k].as<int>();
    };
    auto str = [&](const char* k, std::string& v) {
        if (o.hasOwnProperty(k)) v = o[k].as<std::string>();
    };
    auto flag = [&](const char* k, bool& v) {
        if (o.hasOwnProperty(k)) v = o[k].as<bool>();
    };
    num("depth", p.ndepth);
    num("freqLow", p.nfa);
    num("freqHigh", p.nfb);
    num("rxFreq", p.nfqso);
    num("txFreq", p.nftx);
    num("apWidth", p.napwid);
    num("qsoProgress", p.nqso_progress);
    num("contest", p.nexp_decode);
    flag("apDecode", p.lft8apon);
    flag("apCqOnly", p.lapcqonly);
    flag("fastSubtract", p.fast_subtract);
    str("myCall", p.mycall);
    str("hisCall", p.hiscall);
    str("hisGrid", p.hisgrid);
    return p;
}

std::string trim(const std::string& s) {
    const size_t b = s.find_first_not_of(' ');
    if (b == std::string::npos) return {};
    return s.substr(b, s.find_last_not_of(' ') - b + 1);
}

void push_decodes(val& out, const std::vector<ft8::Decode>& ds) {
    for (const ft8::Decode& d : ds) {
        val o = val::object();
        o.set("snr", d.snr);
        o.set("dt", d.dt);
        o.set("freq", d.freq);
        std::string m = d.msg37;
        const bool low_confidence = m.size() >= 37 && m[36] == '?';
        if (low_confidence) m.resize(36);
        o.set("message", trim(m));
        o.set("lowConfidence", low_confidence);
        o.set("ap", d.nap);
        o.set("quality", d.qual);
        o.set("sync", d.sync);
        o.set("line", d.line);
        out.call<void>("push", o);
    }
}

// What `jt9 -k` reads (msk144_jt9.cpp): -p, -f, -F (limited to 1000), -d, -c,
// -x, and ours, --msk-shmsg and --msk-swl. Defaults: a 15 s period, 1500 Hz,
// 200 Hz either side, depth 3, blank calls.
struct MskParams {
    int ntr = 15, nfqso = 1500, ntol = 200, ndepth = 3;
    std::string mycall, hiscall;
    bool shmsg = false, swl = false;
};

MskParams msk_params_from(const val& o) {
    MskParams p;
    if (o.isUndefined() || o.isNull()) return p;
    auto num = [&](const char* k, int& v) {
        if (o.hasOwnProperty(k)) v = o[k].as<int>();
    };
    auto str = [&](const char* k, std::string& v) {
        if (o.hasOwnProperty(k)) v = o[k].as<std::string>();
    };
    auto flag = [&](const char* k, bool& v) {
        if (o.hasOwnProperty(k)) v = o[k].as<bool>();
    };
    num("period", p.ntr);
    num("rxFreq", p.nfqso);
    num("tolerance", p.ntol);
    num("depth", p.ndepth);
    str("myCall", p.mycall);
    str("hisCall", p.hiscall);
    flag("shortMessages", p.shmsg);
    flag("swl", p.swl);
    p.ntol = std::min(p.ntol, 1000);
    return p;
}

// One MSK144 decode line (mskrtd's format 1021, i6.6,i4,f5.1,i5,a4,a37,
// written a80) as its fields.
val msk_decode_of(const std::string& line) {
    val o = val::object();
    o.set("snr", std::stoi(line.substr(6, 4)));
    o.set("dt", std::stod(line.substr(10, 5)));
    o.set("freq", std::stoi(line.substr(15, 5)));
    o.set("message", trim(line.size() > 24 ? line.substr(24, 37) : std::string()));
    o.set("line", line);
    return o;
}

class Decoder {
public:
    explicit Decoder(const std::string& mode) {
        if (mode == "ft8") {
            mode_ = 8;
            ft8_ = std::make_unique<ft8::Decoder>();
        } else if (mode == "ft4") {
            mode_ = 4;
            ft4_ = std::make_unique<ft4::Decoder>();
        } else if (mode == "ft2") {
            mode_ = 2;
            ft2_ = std::make_unique<ft2::Decoder>();
        } else if (mode == "msk144") {
            mode_ = 144;
            msk_ = std::make_unique<msk144::RtdState>();
        } else {
            throw std::invalid_argument("mode must be ft8, ft4, ft2 or msk144");
        }
    }

    std::string mode() const {
        return mode_ == 8 ? "ft8" : mode_ == 4 ? "ft4" : mode_ == 2 ? "ft2" : "msk144";
    }

    // samples: the period's 12 kHz audio from its start (FT8 15 s, FT4
    // 7.5 s, FT2 3.75 s; shorter is zero-filled). utc: HHMMSS. options:
    // depth (1-3), freqLow, freqHigh, rxFreq, txFreq, apWidth, myCall,
    // hisCall, hisGrid, qsoProgress, contest, apDecode, apCqOnly, fastSubtract.
    // MSK144: the period's audio (period, 15 or 30 s), options period,
    // rxFreq, tolerance, depth, myCall, hisCall, shortMessages, swl.
    // Returns { decodes: [...], count, text } (text: jt9's stdout lines).
    val decode(const val& samples, int utc, const val& options) {
        if (mode_ == 144) return decode_msk(samples, utc, msk_params_from(options));
        const ft8::Params p = params_from(options);
        // An Int16Array through the decoders' 16-bit path (jt9's 16-bit
        // files: the RMS gate as theirs); anything else as floats -1..1 in
        // 16-bit units (jt9's float files: the gate passes any sound).
        if (samples["constructor"]["name"].as<std::string>() == "Int16Array") {
            const std::vector<int16_t> s = emscripten::convertJSArrayToNumberVector<int16_t>(samples);
            return run(s.data(), s.size(), utc, p);
        }
        std::vector<float> s = emscripten::convertJSArrayToNumberVector<float>(samples);
        for (float& x : s) x *= 32768.0f;
        return run(s.data(), s.size(), utc, p);
    }

private:
    // jt9 -k's file mode on one period's samples (msk144_jt9.cpp's decode_file)
    val decode_msk(const val& samples, int utc, const MskParams& p) {
        // The samples jt9 reads of a file: whole blocks of 3456 within the
        // period, none past the block that reaches nhsym 181.
        constexpr long KSTEP = 3456;
        const long npts = static_cast<long>(p.ntr) * 12000L;
        long take = std::max(0L, npts / KSTEP) * KSTEP;
        take = std::min(take, 182L * KSTEP);
        std::string out;
        msk144::Packjt77& pk = msk144::msk144_packjt77();
        msk_->reset();
        msk144::msk144_reset_packjt77();
        std::ostringstream os;
        if (samples["constructor"]["name"].as<std::string>() == "Int16Array") {
            std::vector<int16_t> s = emscripten::convertJSArrayToNumberVector<int16_t>(samples);
            if (static_cast<long>(s.size()) > take) s.resize(static_cast<size_t>(take));
            msk144::decode_msk144(*msk_, pk, s.data(), static_cast<long>(s.size()), utc, p.ntr, p.ntol, p.nfqso,
                                  p.ndepth, p.mycall, p.hiscall, p.shmsg, p.swl, os);
        } else {
            std::vector<float> s = emscripten::convertJSArrayToNumberVector<float>(samples);
            if (static_cast<long>(s.size()) > take) s.resize(static_cast<size_t>(take));
            for (float& x : s) x *= 32768.0f;
            msk144::decode_msk144(*msk_, pk, s.data(), static_cast<long>(s.size()), utc, p.ntr, p.ntol, p.nfqso,
                                  p.ndepth, p.mycall, p.hiscall, p.shmsg, p.swl, os);
        }
        val r = val::object();
        val ds = val::array();
        std::istringstream in(os.str());
        std::string line;
        int n = 0;
        while (std::getline(in, line)) {
            ds.call<void>("push", msk_decode_of(line));
            ++n;
        }
        // <DecodeFinished> as jt9 -k prints it: MSK144 counts nothing there
        r.set("decodes", ds);
        r.set("count", n);
        r.set("text", os.str() + "<DecodeFinished>   0   0        0\n");
        return r;
    }

    template <class S>
    val run(const S* s, size_t n, int utc, const ft8::Params& p) {
        val r = val::object();
        val ds = val::array();
        if (mode_ == 8) {
            const ft8::Decoder::Period per = ft8_->period(s, n, utc, p);
            push_decodes(ds, per.c41.decodes);
            push_decodes(ds, per.c47.decodes);
            push_decodes(ds, per.c50.decodes);
            r.set("count", per.ndecoded);
            r.set("text", per.out);
        } else if (mode_ == 4) {
            const ft4::Decoder::Period per = ft4_->period(s, n, utc, p);
            push_decodes(ds, per.c.decodes);
            r.set("count", per.ndecoded);
            r.set("text", per.out);
        } else {
            const ft2::Decoder::Period per = ft2_->period(s, n, utc, p);
            push_decodes(ds, per.c.decodes);
            r.set("count", per.ndecoded);
            r.set("text", per.out);
        }
        r.set("decodes", ds);
        return r;
    }

    int mode_ = 8;
    std::unique_ptr<ft8::Decoder> ft8_;
    std::unique_ptr<ft4::Decoder> ft4_;
    std::unique_ptr<ft2::Decoder> ft2_;
    std::unique_ptr<msk144::RtdState> msk_;
};

}  // namespace

EMSCRIPTEN_BINDINGS(ft_wasm) {
    emscripten::class_<Decoder>("Decoder")
        .constructor<std::string>()
        .function("decode", &Decoder::decode)
        .property("mode", &Decoder::mode);
}
