// Not WSJT-X's: carrier_iq.h.
#include "carrier_iq.h"

#include <algorithm>
#include <cctype>
#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>

namespace carrier {

std::vector<std::complex<float>> decimate(const std::complex<float>* x, long n, double fs_in, int factor,
                                          double pass_hz, double stop_hz) {
    // Blackman's transition is about 5.5/L of the rate: L taps for pass..stop.
    int half = static_cast<int>(std::ceil(2.75 * fs_in / (stop_hz - pass_hz)));
    const int len = 2 * half + 1;
    const double fc = 0.5 * (pass_hz + stop_hz) / fs_in;  // cutoff, cycles a sample
    std::vector<double> h(static_cast<size_t>(len));
    double sum = 0.0;
    for (int i = 0; i < len; ++i) {
        const int k = i - half;
        const double s = k == 0 ? 2.0 * fc : std::sin(2.0 * M_PI * fc * k) / (M_PI * k);
        const double w = 0.42 - 0.5 * std::cos(2.0 * M_PI * i / (len - 1)) + 0.08 * std::cos(4.0 * M_PI * i / (len - 1));
        h[static_cast<size_t>(i)] = s * w;
        sum += h[static_cast<size_t>(i)];
    }
    for (double& v : h) v /= sum;  // unity gain at DC
    const long nout = n > 0 ? (n - 1) / factor + 1 : 0;
    std::vector<std::complex<float>> y(static_cast<size_t>(nout));
    for (long m = 0; m < nout; ++m) {
        const long c = m * factor;
        double re = 0.0, im = 0.0;
        const long lo = std::max(0L, c - half), hi = std::min(n - 1, c + half);
        for (long k = lo; k <= hi; ++k) {
            const double t = h[static_cast<size_t>(k - c + half)];
            re += t * x[k].real();
            im += t * x[k].imag();
        }
        y[static_cast<size_t>(m)] = {static_cast<float>(re), static_cast<float>(im)};
    }
    return y;
}

namespace {
void put32(std::string& s, uint32_t v) {
    for (int i = 0; i < 4; ++i) s += static_cast<char>((v >> (8 * i)) & 0xff);
}
void put16(std::string& s, uint16_t v) {
    s += static_cast<char>(v & 0xff);
    s += static_cast<char>(v >> 8);
}
}  // namespace

bool write_wav(const std::string& path, const std::vector<std::complex<float>>& iq, unsigned rate,
               const std::string& comment, std::string* err) {
    std::string icmt = comment;
    icmt += '\0';
    if (icmt.size() % 2) icmt += '\0';
    std::string list = "INFO";
    list += "ICMT";
    put32(list, static_cast<uint32_t>(icmt.size()));
    list += icmt;
    const uint32_t data_bytes = static_cast<uint32_t>(iq.size() * 8);
    std::string h = "RIFF";
    put32(h, static_cast<uint32_t>(4 + (8 + 16) + (8 + list.size()) + (8 + data_bytes)));
    h += "WAVEfmt ";
    put32(h, 16);
    put16(h, 3);  // WAVE_FORMAT_IEEE_FLOAT
    put16(h, 2);  // I, Q
    put32(h, rate);
    put32(h, rate * 8);
    put16(h, 8);
    put16(h, 32);
    h += "LIST";
    put32(h, static_cast<uint32_t>(list.size()));
    h += list;
    h += "data";
    put32(h, data_bytes);
    const std::string tmp = path + ".tmp";
    FILE* f = std::fopen(tmp.c_str(), "wb");
    if (!f) {
        if (err) *err = "cannot open " + tmp + ": " + std::strerror(errno);
        return false;
    }
    bool ok = std::fwrite(h.data(), 1, h.size(), f) == h.size();
    for (const auto& c : iq) {  // little-endian float32 pairs (the hosts we build for)
        const float v[2] = {c.real(), c.imag()};
        ok = ok && std::fwrite(v, sizeof v, 1, f) == 1;
    }
    ok = (std::fclose(f) == 0) && ok;
    if (!ok || std::rename(tmp.c_str(), path.c_str()) != 0) {
        if (err) *err = "cannot write " + path + ": " + std::strerror(errno);
        std::remove(tmp.c_str());
        return false;
    }
    return true;
}

std::string file_name(const std::string& message, double ref_hz) {
    std::string call = message.substr(0, message.find(' '));
    if (call.size() >= 2 && call.front() == '<' && call.back() == '>') call = call.substr(1, call.size() - 2);
    if (call.empty() || call == "...") return "";
    for (char& c : call)
        if (!std::isalnum(static_cast<unsigned char>(c))) c = '_';
    char b[64];
    std::snprintf(b, sizeof b, "_%.3f.wav", ref_hz);
    return call + b;
}

std::string json_string(const std::string& s) {
    std::string o = "\"";
    for (char c : s) {
        if (c == '"' || c == '\\') o += '\\', o += c;
        else if (static_cast<unsigned char>(c) < 0x20) {
            char b[8];
            std::snprintf(b, sizeof b, "\\u%04x", c);
            o += b;
        } else o += c;
    }
    return o + "\"";
}

}  // namespace carrier
