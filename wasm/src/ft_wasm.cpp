// The FT8, FT4 and FT2 decoders for JavaScript (embind): one Decoder is one
// jt9 process for one dial, as ft8::Decoder, ft4::Decoder and ft2::Decoder
// are (their state carried from period to period: the a7 tables, packjt77's
// hashes). A period's 12 kHz samples in, its decodes out, as WSJT-X's GUI
// has jt9 decode them (ft8.h's, ft4.h's, ft2.h's period()).
//
// Samples: a Float32Array (-1..1, scaled by 32768 as a float .wav file is
// read) or an Int16Array (as they are).
#include <emscripten/bind.h>
#include <emscripten/val.h>

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "ft2.h"
#include "ft4.h"
#include "ft8.h"

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
        } else {
            throw std::invalid_argument("mode must be ft8, ft4 or ft2");
        }
    }

    std::string mode() const { return mode_ == 8 ? "ft8" : mode_ == 4 ? "ft4" : "ft2"; }

    // samples: the period's 12 kHz audio from its start (FT8 15 s, FT4
    // 7.5 s, FT2 3.75 s; shorter is zero-filled). utc: HHMMSS. options:
    // depth (1-3), freqLow, freqHigh, rxFreq, txFreq, apWidth, myCall,
    // hisCall, hisGrid, qsoProgress, contest, apDecode, apCqOnly, fastSubtract.
    // Returns { decodes: [...], count, text } (text: jt9's stdout lines).
    val decode(const val& samples, int utc, const val& options) {
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
};

}  // namespace

EMSCRIPTEN_BINDINGS(ft_wasm) {
    emscripten::class_<Decoder>("Decoder")
        .constructor<std::string>()
        .function("decode", &Decoder::decode)
        .property("mode", &Decoder::mode);
}
