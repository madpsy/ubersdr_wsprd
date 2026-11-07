// The port's side of scripts/msk144-piece-check.sh msk40: reads the oracle's
// output (scripts/msk40-oracle.f90) on stdin, takes its EN, HA, PK, NH, BI,
// FP, FC, RUN, YP, SC, YR, ZP, ZB and ZR lines as inputs, and prints every
// line as the reference printed it, the results from src/msk40.cpp instead
// (the trace lines through the port's hooks). The two outputs must be
// identical byte for byte.
//
// The module packjt77's state is shared by the reference's calls (the recent
// calls): the driver makes the same pack77 and unpack77 calls on one Packjt77
// in the order the oracle's lines come. nhasharray is the port's Msk40State,
// updated where the oracle's NH line is.
//
// While msk40spd runs, the frame decoder's traces are off (the oracle's copy
// of msk40spd calls the library's frame decoder, which prints nothing); the
// sync's stay on.
//
// Blind spot: NaNs compare equal whatever their sign and payload (bits_of
// prints the default NaN), as in the oracle's header.
//
// License: GNU GPL v3 (as WSJT-X)
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "msk144_decode.h"
#include "msk40.h"

namespace {

using msk144::cfloat;

uint32_t bits_of(float x) {
    if (x != x) return 0x7FC00000u;
    uint32_t b;
    std::memcpy(&b, &x, 4);
    return b;
}
float float_of(const std::string& hex) {
    const uint32_t b = static_cast<uint32_t>(std::strtoul(hex.c_str(), nullptr, 16));
    float x;
    std::memcpy(&x, &b, 4);
    return x;
}

// prr and prc of the oracle: 16 values a line
void prr(const char* tag, int k, const float* x, int n) {
    for (int i = 1; i <= n; i += 16) {
        std::printf("%s %d %d", tag, k, i);
        for (int j = i; j <= std::min(i + 15, n); ++j) std::printf(" %08X", bits_of(x[j - 1]));
        std::printf("\n");
    }
}
void prc(const char* tag, int k, const cfloat* x, int n) {
    std::vector<float> v(2 * static_cast<size_t>(std::max(n, 0)));
    for (int j = 0; j < n; ++j) {
        v[2 * j] = x[j].real();
        v[2 * j + 1] = x[j].imag();
    }
    prr(tag, k, v.data(), 2 * n);
}

int kc = 0;

// the traces, as the oracle's hooks print
void t_bp(int iter, const float zn[33], const float tov[4][33]) {
    prr("BZ", kc * 100 + iter, zn + 1, 32);
    float v[96];
    for (int i = 1; i <= 32; ++i)
        for (int kk = 1; kk <= 3; ++kk) v[(i - 1) * 3 + (kk - 1)] = tov[kk][i];
    prr("BV", kc * 100 + iter, v, 96);
}
void t_df_consts(const float* pp, const cfloat* cb) {
    prr("KP", kc, pp, 12);
    prc("KB", kc, cb, 42);
}
void t_df_phase(cfloat cca, float phase0) {
    std::printf("KH %d %08X %08X %08X\n", kc, bits_of(cca.real()), bits_of(cca.imag()), bits_of(phase0));
}
void t_df_corrected(cfloat cfac, const cfloat* c) {
    std::printf("KF %d %08X %08X\n", kc, bits_of(cfac.real()), bits_of(cfac.imag()));
    prc("KC", kc, c, 240);
}
void t_df_raw(const float* softbits) { prr("KS", kc, softbits, 40); }
void t_df_sync(int n1) { std::printf("KY %d %d\n", kc, n1); }
void t_df_norm(float sav, float s2av, float ssig) {
    std::printf("KN %d %08X %08X %08X\n", kc, bits_of(sav), bits_of(s2av), bits_of(ssig));
}
void t_df_llr(float sigma, const float* llr) {
    std::printf("KG %d %08X\n", kc, bits_of(sigma));
    prr("KL", kc, llr, 32);
}
void t_df_bp(const signed char* decoded, int niterations) {
    if (niterations >= 0) {
        std::string sd;
        for (int i = 0; i < 16; ++i) sd += static_cast<char>('0' + decoded[i]);
        std::printf("KD %d %d |%s|\n", kc, niterations, sd.c_str());
    } else {
        std::printf("KD %d %d |-|\n", kc, niterations);
    }
}
void t_df_msg(int nhammd, float cord, int imsg, int nrxrpt, int nrxhash, int ihash) {
    std::printf("KM %d %d %08X %d %d %d\n", kc, nhammd, bits_of(cord), imsg, nrxrpt, nrxhash);
    std::printf("KI %d %d\n", kc, ihash);
}
void t_sync_consts(const float* pp, const cfloat* cb) {
    prr("YKP", kc, pp, 12);
    prc("YKB", kc, cb, 42);
}
void t_sync_peak(float xmax, float fest, const cfloat* c, const float* xcc, int nfreqs) {
    std::printf("YX %d %08X %08X %d\n", kc, bits_of(xmax), bits_of(fest), nfreqs);
    prc("YW", kc, c, 240);
    prr("YV", kc, xcc, 240);
}
void t_spd_consts(const float* rcw) { prr("R", kc, rcw, 12); }
void t_spd_windows(int ihlo, int ihhi, int illo, int ilhi, int i2000, int i4000) {
    std::printf("N %d %d %d %d %d %d %d\n", kc, ihlo, ihhi, illo, ilhi, i2000, i4000);
}
void t_spd_step(int istp, int ihpk, int ilpk, float deltah, float deltal, float ah, float al, float trath, float tratl,
                float ferr) {
    std::printf("T %d %d %d %d %08X %08X %08X %08X %08X %08X %08X\n", kc, istp, ihpk, ilpk, bits_of(deltah),
                bits_of(deltal), bits_of(ah), bits_of(al), bits_of(trath), bits_of(tratl), bits_of(ferr));
}
void t_spd_median(int nstep, const int* indices, float xmed, const float* detmet, const float* detmet2,
                  const float* detfer) {
    prr("D", kc, detmet, nstep);
    prr("H", kc, detmet2, nstep);
    prr("O", kc, detfer, nstep);
    std::printf("M %d %d %08X", kc, nstep, bits_of(xmed));
    for (int j = 0; j < nstep; ++j) std::printf(" %d", indices[j]);
    std::printf("\n");
}
void t_spd_candidates(int ndet, const int* nstart, const float* ferrs, const float* snrs) {
    std::printf("E %d %d %d\n", kc, 0, ndet);
    for (int j = 0; j < ndet; ++j)
        std::printf("E %d %d %d %08X %08X\n", kc, j + 1, nstart[j], bits_of(ferrs[j]), bits_of(snrs[j]));
}
void t_spd_sync(int icand, int iav, int nsyncsuccess, float fest, const int* npkloc) {
    std::printf("J %d %d %d %d %08X %d %d\n", kc, icand, iav, nsyncsuccess, bits_of(fest), npkloc[0], npkloc[1]);
}

// one oracle line's values from its fourth field on: TAG k i v...
void store(const std::string& line, std::vector<float>& into) {
    std::istringstream is(line);
    std::string tag, hex;
    int k, i;
    is >> tag >> k >> i;
    size_t at = static_cast<size_t>(i - 1);
    while (is >> hex) {
        if (into.size() <= at) into.resize(at + 1);
        into[at++] = float_of(hex);
    }
}

std::vector<cfloat> as_complex(const std::vector<float>& f, int n) {
    std::vector<cfloat> c(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) c[static_cast<size_t>(i)] = cfloat(f.at(2 * i), f.at(2 * i + 1));
    return c;
}

// the text between the n-th bar and the next, as Fortran's a(len) wrote it
std::string between_bars(const std::string& line, size_t& pos, size_t len) {
    const size_t bar = line.find('|', pos);
    std::string f = line.substr(bar + 1, len);
    pos = bar + 1 + len;  // at the closing '|'
    return f;
}

void set_df_traces(bool on) {
    msk144::Msk40Trace& t = msk144::msk40_trace;
    t.df_consts = on ? t_df_consts : nullptr;
    t.df_phase = on ? t_df_phase : nullptr;
    t.df_corrected = on ? t_df_corrected : nullptr;
    t.df_raw = on ? t_df_raw : nullptr;
    t.df_sync = on ? t_df_sync : nullptr;
    t.df_norm = on ? t_df_norm : nullptr;
    t.df_llr = on ? t_df_llr : nullptr;
    t.df_bp = on ? t_df_bp : nullptr;
    t.df_msg = on ? t_df_msg : nullptr;
}

}  // namespace

int main() {
    msk144::Packjt77 pk;
    msk144::Msk40State ms;
    std::map<int, std::vector<float>> data;  // case -> its input (re/im pairs, or llr)
    struct FrameArgs {
        std::string mc, hc;
        bool bswl = false;
        float xsnr = 0.0f;
    };
    std::map<int, FrameArgs> fargs;
    struct SyncArgs {
        float fo = 0.0f;
        int navmask[3] = {0, 0, 0};
    };
    std::map<int, SyncArgs> sargs;
    struct SpdArgs {
        std::string mc, hc;
        bool bswl = false;
        int ntol = 0;
        float fc = 0.0f;
    };
    std::map<int, SpdArgs> pargs;

    msk144::Msk40Trace& t = msk144::msk40_trace;
    // (t.bp only around the direct bpdecode40 cases: the oracle's copy of the frame decoder calls the library's)
    set_df_traces(true);
    t.sync_consts = t_sync_consts;
    t.sync_peak = t_sync_peak;
    t.spd_consts = t_spd_consts;
    t.spd_windows = t_spd_windows;
    t.spd_step = t_spd_step;
    t.spd_median = t_spd_median;
    t.spd_candidates = t_spd_candidates;
    t.spd_sync = t_spd_sync;

    std::string line;
    while (std::getline(std::cin, line)) {
        std::istringstream is(line);
        std::string tag;
        is >> tag;
        if (tag == "EN") {
            int it;
            is >> it;
            size_t pos = 0;
            const std::string sm = between_bars(line, pos, 16);
            pos += 1;
            between_bars(line, pos, 32);
            signed char m[16], cw[32];
            for (int i = 0; i < 16; ++i) m[i] = static_cast<signed char>(sm[static_cast<size_t>(i)] - '0');
            msk144::encode_msk40(m, cw);
            std::string sc;
            for (int i = 0; i < 32; ++i) sc += static_cast<char>('0' + cw[i]);
            std::printf("EN %d |%s| |%s|\n", it, sm.c_str(), sc.c_str());
        } else if (tag == "EA") {
            uint64_t h1 = 1, h2 = 1;
            for (int nn = 0; nn <= 65535; ++nn) {
                signed char m[16], cw[32];
                for (int ii = 1; ii <= 16; ++ii) m[ii - 1] = static_cast<signed char>((nn >> (ii - 1)) & 1);
                msk144::encode_msk40(m, cw);
                for (int jj = 1; jj <= 32; ++jj) {
                    const uint64_t w = static_cast<uint64_t>(cw[jj - 1]);
                    h1 = (h1 * 1000003ULL + w + 1ULL) % 2147483647ULL;
                    h2 = (h2 * 999983ULL + w * static_cast<uint64_t>(jj + 1)) % 2147483629ULL;
                }
            }
            std::printf("EA %d %llu %llu\n", 65536, static_cast<unsigned long long>(h1), static_cast<unsigned long long>(h2));
        } else if (tag == "HA") {
            int it, n;
            is >> it >> n;
            const size_t bar = line.find('|');
            const std::string text = line.substr(bar + 1, static_cast<size_t>(n));
            std::printf("HA %d %d |%s| %d\n", it, n, text.c_str(), msk144::hash_f90(text.data(), n));
        } else if (tag == "PK") {
            int it;
            is >> it;
            size_t pos = 0;
            const std::string text = between_bars(line, pos, 37);
            const msk144::Packjt77::Packed p = pk.pack77_legacy_truncating_fallback(text);
            std::printf("PK %d |%s| %d %d %s\n", it, text.c_str(), p.i3, p.n3, p.c77.c_str());
            if (p.i3 < 0 || p.n3 < 0) continue;
            const msk144::Packjt77::Unpacked u = pk.unpack77(p.c77.data(), 1);
            std::string msg = u.msg;
            msg.resize(37, ' ');
            std::printf("UN %d |%s|\n", it, msg.c_str());
            std::printf("RC %d", it);
            for (int i = 1; i <= msk144::Packjt77::MAXRECENT; ++i) std::printf("|%s", pk.recent_call(i).c_str());
            std::printf("\n");
        } else if (tag == "NH") {
            msk144::update_msk40_hasharray(ms, pk);
            std::printf("NH %d", 1);
            for (int j = 0; j < msk144::Packjt77::MAXRECENT; ++j)
                for (int i = 0; i < msk144::Packjt77::MAXRECENT; ++i) std::printf(" %d", ms.nhasharray[i][j]);
            std::printf("\n");
        } else if (tag == "BI") {
            int k;
            is >> k;
            kc = k;
            std::vector<float>& v = data[k];
            v.clear();
            store(line, v);
            // the second line of the 32 values
            std::string next;
            std::getline(std::cin, next);
            store(next, v);
            prr("BI", k, v.data(), 32);
            signed char dec[16] = {};
            int niter = -99;
            msk144::msk40_trace.bp = t_bp;
            msk144::bpdecode40(v.data(), 5, dec, niter);
            msk144::msk40_trace.bp = nullptr;
            if (niter >= 0) {
                std::string sd;
                for (int i = 0; i < 16; ++i) sd += static_cast<char>('0' + dec[i]);
                std::printf("BD %d %d |%s|\n", k, niter, sd.c_str());
            } else {
                std::printf("BD %d %d |-|\n", k, niter);
            }
        } else if (tag == "FP") {
            int k;
            is >> k;
            FrameArgs fa;
            // |mycall|hiscall| then " bswl xsnr": the middle bar closes the one and opens the other
            const size_t b0 = line.find('|');
            fa.mc = line.substr(b0 + 1, 12);
            fa.hc = line.substr(b0 + 14, 12);
            std::istringstream is2(line.substr(b0 + 27));
            std::string sw, sx;
            is2 >> sw >> sx;
            fa.bswl = sw == "T";
            fa.xsnr = float_of(sx);
            fargs[k] = fa;
            data[k].clear();
            std::puts(line.c_str());
        } else if (tag == "FC" || tag == "SC" || tag == "ZB") {
            int k;
            is >> k;
            store(line, data[k]);
            // (the printed lines are the oracle's, given back with the run that uses them)
        } else if (tag == "RUN") {
            int k;
            is >> k;
            kc = k;
            const FrameArgs& fa = fargs.at(k);
            std::vector<cfloat> c = as_complex(data.at(k), 240);
            prc("FC", k, c.data(), 240);
            std::puts(line.c_str());
            char msg[37];
            std::memset(msg, ' ', 37);
            std::memcpy(msg, "unset", 5);
            int ns = -1;
            msk144::msk40decodeframe(ms, pk, c.data(), fa.mc, fa.hc, fa.xsnr, fa.bswl, msg, ns);
            std::printf("R %d %d |%.37s|\n", k, ns, msg);
            prc("O", k, c.data(), 240);  // c as the decoder left it
        } else if (tag == "YP") {
            int k;
            is >> k;
            std::string sf;
            SyncArgs sa;
            is >> sf >> sa.navmask[0] >> sa.navmask[1] >> sa.navmask[2];
            sa.fo = float_of(sf);
            sargs[k] = sa;
            data[k].clear();
            std::puts(line.c_str());
        } else if (tag == "YR") {
            int k;
            is >> k;
            kc = k;
            const SyncArgs& sa = sargs.at(k);
            const std::vector<cfloat> cd = as_complex(data.at(k), 720);
            prc("SC", k, cd.data(), 720);
            std::puts(line.c_str());
            cfloat c[240];
            int npk[2] = {-1, -1}, ns = -9;
            float fest = -1.0f;
            msk144::msk40sync(cd.data(), 3, 29, 7.2f, sa.navmask, 2, sa.fo, fest, npk, ns, c);
            std::printf("YA %d %d %08X %d %d\n", k, ns, bits_of(fest), npk[0], npk[1]);
            prc("YC", k, c, 240);
        } else if (tag == "ZP") {
            int k;
            is >> k;
            SpdArgs pa;
            const size_t b0 = line.find('|');
            pa.mc = line.substr(b0 + 1, 12);
            pa.hc = line.substr(b0 + 14, 12);
            std::istringstream is2(line.substr(b0 + 27));
            std::string sw, sf;
            is2 >> sw >> pa.ntol >> sf;
            pa.bswl = sw == "T";
            pa.fc = float_of(sf);
            pargs[k] = pa;
            data[k].clear();
            std::puts(line.c_str());
        } else if (tag == "ZR") {
            int k;
            is >> k;
            kc = k;
            const SpdArgs& pa = pargs.at(k);
            const std::vector<cfloat> z = as_complex(data.at(k), 6912);
            prc("ZB", k, z.data(), 6912);
            std::puts(line.c_str());
            char msg[37];
            std::memset(msg, ' ', 37);
            std::memcpy(msg, "unset", 5);
            float fr = -5.0f, tr = -5.0f;
            int nv = -5, ns = -9;
            set_df_traces(false);
            msk144::msk40spd(ms, pk, z.data(), 6912, pa.ntol, pa.mc, pa.hc, pa.bswl, ns, msg, pa.fc, fr, tr, nv);
            set_df_traces(true);
            std::printf("X %d %d |%.37s| %08X %08X %d\n", k, ns, msg, bits_of(fr), bits_of(tr), nv);
        } else if (tag == "KP" || tag == "KB" || tag == "KH" || tag == "KF" || tag == "KC" || tag == "KS" || tag == "KY" ||
                   tag == "KN" || tag == "KG" || tag == "KL" || tag == "KD" || tag == "KM" || tag == "KI" || tag == "R" ||
                   tag == "O" || tag == "BZ" || tag == "BV" || tag == "BD" || tag == "YKP" || tag == "YKB" || tag == "YX" ||
                   tag == "YW" || tag == "YV" || tag == "YA" || tag == "YC" || tag == "N" || tag == "T" || tag == "D" ||
                   tag == "H" || tag == "M" || tag == "E" || tag == "J" || tag == "X" || tag == "UN" || tag == "RC" ||
                   tag == "END") {
            continue;
        } else {
            std::fprintf(stderr, "msk40-check: unknown line: %s\n", line.c_str());
            return 1;
        }
    }
    std::printf("END msk40 %d\n", kc);
    return 0;
}
