// The WSPR decoder for JavaScript (embind): one WsprDecoder is wsprd run on
// one dial's periods in one data directory (wsprd.h's Decoder): its hash
// table (hashtable.txt) carried from period to period. A period's 12 kHz
// samples in, its decodes out, as wsprd prints and writes them.
//
// Samples: an Int16Array (a .wav's) or a Float32Array (-1..1, a float
// .wav's).
#include <emscripten/bind.h>
#include <emscripten/val.h>

#include <cstdint>
#include <string>
#include <vector>

#include "wsprd.h"

using emscripten::val;

namespace {

// wsprd's Options (as WSJT-X runs wsprd: -C 500 -o 4) with the options given.
wsprd::Options options_from(const val& o) {
    wsprd::Options opt;
    opt.maxcycles = 500;
    opt.ndepth = 4;
    if (o.isUndefined() || o.isNull()) return opt;
    auto has = [&](const char* k) { return o.hasOwnProperty(k) && !o[k].isUndefined(); };
    if (has("dialFreq")) opt.dialfreq = o["dialFreq"].as<double>();
    if (has("dialError")) opt.dialfreq_error = o["dialError"].as<double>();
    if (has("mode")) opt.wspr_type = o["mode"].as<int>();
    if (has("cycles")) opt.maxcycles = o["cycles"].as<unsigned>();
    if (has("osdDepth")) opt.ndepth = o["osdDepth"].as<int>();
    if (has("deep")) opt.more_candidates = o["deep"].as<bool>();
    if (has("quick")) opt.quickmode = o["quick"].as<bool>();
    if (has("subtraction") && !o["subtraction"].as<bool>()) {
        opt.subtraction = false;  // -s: one pass
        opt.npasses = 1;
    }
    if (has("blockDemod") && !o["blockDemod"].as<bool>()) opt.npasses = 2;  // -B
    if (has("stackDecoder")) opt.stackdecoder = o["stackDecoder"].as<bool>();
    if (has("wide") && o["wide"].as<bool>()) {
        opt.fmin = -150;
        opt.fmax = 150;
    }
    if (has("bias")) opt.bias = o["bias"].as<float>();
    if (has("noDrift")) opt.nodrift = o["noDrift"].as<bool>();
    if (has("hashtable")) opt.usehashtable = o["hashtable"].as<bool>();
    return opt;
}

val decode_object(const wsprd::Decode& d, bool nodrift) {
    val o = val::object();
    o.set("date", d.date);
    o.set("time", d.time);
    o.set("snr", d.snr);
    o.set("dt", d.dt);
    o.set("freq", d.freq);
    o.set("drift", d.drift);
    o.set("message", d.message);
    o.set("callsign", d.callsign);
    o.set("grid", d.grid);
    o.set("power", d.power);
    o.set("sync", d.sync);
    o.set("spread", d.spread);
    o.set("osd", d.decodetype == 1);
    o.set("pass", d.ipass + 1);
    o.set("line", wsprd::stdout_line(d));
    o.set("allWsprLine", wsprd::all_wspr_line(d, nodrift));
    o.set("spotsLine", wsprd::spots_line(d));
    return o;
}

class WsprDecoder {
public:
    // samples: the period's 12 kHz audio from its start (114 s read,
    // WSPR-15 912 s; shorter is zero-filled). options: date (YYMMDD), time
    // (HHMM), dialFreq (MHz), mode (2 or 15), cycles (-C, default 500),
    // osdDepth (-o, default 4; -1 none), deep (-d), quick (-q),
    // subtraction (false: -s), blockDemod (false: -B), stackDecoder (-J),
    // wide (-w), bias (-z), dialError (-e), noDrift (-n), hashtable
    // (false: -H). Returns { decodes, text } (text: wsprd's stdout), or
    // throws for an input wsprd refuses.
    val decode(const val& samples, const val& options) {
        dec_.options() = options_from(options);
        std::string date = "000000", time = "0000";
        if (!options.isUndefined() && !options.isNull()) {
            if (options.hasOwnProperty("date")) date = options["date"].as<std::string>();
            if (options.hasOwnProperty("time")) time = options["time"].as<std::string>();
        }
        wsprd::Result res;
        if (samples["constructor"]["name"].as<std::string>() == "Int16Array") {
            const std::vector<int16_t> s = emscripten::convertJSArrayToNumberVector<int16_t>(samples);
            res = dec_.decode_wav(s.data(), s.size(), date, time);
        } else {
            const std::vector<float> s = emscripten::convertJSArrayToNumberVector<float>(samples);
            res = dec_.decode_wav(s.data(), s.size(), date, time);
        }
        if (res.status != 0)
            throw std::runtime_error("wsprd refused the period (status " + std::to_string(res.status) + ")");
        val r = val::object();
        val ds = val::array();
        std::string text;
        for (const wsprd::Decode& d : res.decodes) {
            ds.call<void>("push", decode_object(d, dec_.options().nodrift));
            text += wsprd::stdout_line(d);
        }
        r.set("decodes", ds);
        r.set("text", text + "<DecodeFinished>\n");
        return r;
    }

    // hashtable.txt as wsprd writes and reads it.
    std::string hashtable() { return dec_.hashtable().save(); }
    void load_hashtable(const std::string& text) { dec_.hashtable().load(text); }

    // Frees the period's buffers between decodes (keeps the hash table).
    void release() { dec_.release(); }

private:
    wsprd::Decoder dec_;
};

}  // namespace

EMSCRIPTEN_BINDINGS(wspr_wasm) {
    emscripten::class_<WsprDecoder>("WsprDecoder")
        .constructor<>()
        .function("decode", &WsprDecoder::decode)
        .function("hashtable", &WsprDecoder::hashtable)
        .function("loadHashtable", &WsprDecoder::load_hashtable)
        .function("release", &WsprDecoder::release);
}
