// The port's side of scripts/msk144-piece-check.sh sync: reads the oracle's
// output (scripts/msk144-sync-oracle.f90) on stdin, takes its I, L, J, K, W,
// SYN and SPD lines as inputs, and prints every line as the reference printed
// it, the results from src/msk144_sync.cpp instead (the trace lines through
// the port's hooks). The two outputs must be identical byte for byte.
//
// msk144spd's frame decoder is the port's own (src/msk144_decode.cpp): it
// decodes the frames the port's msk144spd cuts, in place as the reference's
// does, and the F lines (a hash of the frame it was given, and what the
// decoder answered), the X, U (softbits) and V (ct) lines after each
// msk144spd are all the port's: a wrong frame, answer, softbits, ct or
// msgreceived shows. One Packjt77 serves the whole run, as the reference's
// module packjt77 does.
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
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "msk144_decode.h"
#include "msk144_sync.h"

namespace {

using msk144::cfloat;
constexpr int NSPM = 864;

uint32_t bits_of(float x) {
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
double double_of(const std::string& hex) {
    const uint64_t b = std::strtoull(hex.c_str(), nullptr, 16);
    double x;
    std::memcpy(&x, &b, 8);
    return x;
}
uint64_t dbits_of(double x) {
    uint64_t b;
    std::memcpy(&b, &x, 8);
    return b;
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

// hashr of the oracle
void hashr(const float* x, int n, uint64_t& h1, uint64_t& h2) {
    h1 = 1;
    h2 = 1;
    for (int i = 0; i < n; ++i) {
        const uint64_t w = bits_of(x[i]);
        h1 = (h1 * 1000003ULL + w) % 2147483647ULL;
        h2 = (h2 * 999983ULL + w) % 2147483629ULL;
    }
}

int kcase = 0;

// the traces, as the oracle's hooks print
void t_sync_consts(const float* pp, const cfloat* cb) {
    prr("P", kcase, pp, 12);
    prc("B", kcase, cb, 42);
}
void t_sync_peak(float xmax, float fest, const cfloat* c, const float* xcc, int nfreqs) {
    std::printf("Q %d %d %08X %08X\n", kcase, nfreqs, bits_of(xmax), bits_of(fest));
    prc("Y", kcase, c, NSPM);
    prr("Z", kcase, xcc, NSPM);
}
void t_spd_consts(const float* rcw) { prr("R", kcase, rcw, 12); }
void t_spd_win(int ihlo, int ihhi, int illo, int ilhi, int i2000, int i4000) {
    std::printf("N %d %d %d %d %d %d %d\n", kcase, ihlo, ihhi, illo, ilhi, i2000, i4000);
}
void t_spd_step(int istp, int ihpk, int ilpk, float deltah, float deltal, float ah, float al, float trath, float tratl, float ferr) {
    std::printf("T %d %d %d %d %08X %08X %08X %08X %08X %08X %08X\n", kcase, istp, ihpk, ilpk, bits_of(deltah), bits_of(deltal),
                bits_of(ah), bits_of(al), bits_of(trath), bits_of(tratl), bits_of(ferr));
}
void t_spd_med(int nstep, const int* indices, float xmed, const float* detmet, const float* detmet2, const float* detfer) {
    prr("D", kcase, detmet, nstep);
    prr("H", kcase, detmet2, nstep);
    prr("O", kcase, detfer, nstep);
    std::printf("M %d %d %08X", kcase, nstep, bits_of(xmed));
    for (int j = 0; j < nstep; ++j) std::printf(" %d", indices[j]);
    std::printf("\n");
}
void t_spd_cands(int ndet, const int* nstart, const float* ferrs, const float* snrs) {
    std::printf("E %d %d %d\n", kcase, 0, ndet);
    for (int j = 0; j < ndet; ++j) std::printf("E %d %d %d %08X %08X\n", kcase, j + 1, nstart[j], bits_of(ferrs[j]), bits_of(snrs[j]));
}
void t_spd_sync(int icand, int iav, int nsyncsuccess, float fest, float xmax, const int* npkloc) {
    std::printf("G %d %d %d %d %08X %08X %d %d\n", kcase, icand, iav, nsyncsuccess, bits_of(fest), bits_of(xmax), npkloc[0], npkloc[1]);
}

// one oracle line's values from its fourth field on: TAG k i v... (the same for several lines, stored at i-1...)
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

std::map<int, std::vector<cfloat>> slices;  // s -> the slice's points

std::vector<cfloat> pairs_to_complex(const std::vector<float>& v) {
    std::vector<cfloat> c(v.size() / 2);
    for (size_t i = 0; i < c.size(); ++i) c[i] = cfloat(v[2 * i], v[2 * i + 1]);
    return c;
}

// make_input of the oracle: points off+1..off+n of slice s (s 0: zeros), each
// times scale (not for zeros), then a NaN real part at nanidx, +Inf at infidx (1e20 at -infidx if negative)
void make_input(int s, int off, int n, int nanidx, int infidx, float scale, cfloat* c) {
    if (s == 0) {
        for (int i = 0; i < n; ++i) c[i] = cfloat(0.0f, 0.0f);
    } else {
        const std::vector<cfloat>& sl = slices.at(s);
        for (int i = 0; i < n; ++i) c[i] = sl.at(static_cast<size_t>(off + i));
        for (int i = 0; i < n; ++i) c[i] = cfloat(scale * c[i].real(), scale * c[i].imag());
    }
    if (nanidx > 0) c[nanidx - 1] = cfloat(std::numeric_limits<float>::quiet_NaN(), c[nanidx - 1].imag());
    if (infidx > 0) c[infidx - 1] = cfloat(std::numeric_limits<float>::infinity(), c[infidx - 1].imag());
    if (infidx < 0) c[-infidx - 1] = cfloat(1.0e20f, c[-infidx - 1].imag());
}

}  // namespace

int main() {
    // the oracle's whole output (its F, U and V lines, and the like, are what the port's own lines are compared with)
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(std::cin, line)) lines.push_back(line);
    msk144::Packjt77 pk;

    msk144::sync_trace.consts = t_sync_consts;
    msk144::sync_trace.peak = t_sync_peak;
    msk144::spd_trace.consts = t_spd_consts;
    msk144::spd_trace.windows = t_spd_win;
    msk144::spd_trace.step = t_spd_step;
    msk144::spd_trace.median = t_spd_med;
    msk144::spd_trace.candidates = t_spd_cands;
    msk144::spd_trace.sync = t_spd_sync;

    std::map<int, std::vector<float>> islice;
    int kmax = 0;  // the last case number: the oracle's END count
    for (const std::string& ln : lines) {
        std::istringstream is(ln);
        std::string tag;
        is >> tag;
        if (tag == "I") {
            std::puts(ln.c_str());
            int s, i;
            std::istringstream is2(ln);
            std::string t2;
            is2 >> t2 >> s >> i;
            store(ln, islice[s]);
            slices[s] = pairs_to_complex(islice[s]);
        } else if (tag == "L") {
            int k, kind, n, lo, hi, res;
            is >> k >> kind >> n >> lo >> hi >> res;
            kmax = std::max(kmax, k);
            std::vector<float> a(static_cast<size_t>(std::max(n, 1)));
            std::string hex;
            for (int i = 0; i < n; ++i) {
                is >> hex;
                a[i] = float_of(hex);
            }
            int il;
            if (kind == 0) {
                il = msk144::maxloc_f(a.data(), n);
            } else {
                il = lo <= hi ? msk144::maxloc_masked_f(a.data(), lo, hi) : 0;
            }
            std::printf("L %d %d %d %d %d %d", k, kind, n, lo, hi, il);
            for (int i = 0; i < n; ++i) std::printf(" %08X", bits_of(a[i]));
            std::printf("\n");
        } else if (tag == "J") {
            int k, n;
            is >> k >> n;
            kmax = std::max(kmax, k);
            std::vector<float> a(static_cast<size_t>(n));
            std::string hex;
            for (int i = 0; i < n; ++i) {
                is >> hex;
                a[i] = float_of(hex);
            }
            std::vector<int> idx(static_cast<size_t>(n));
            msk144::indexx(a.data(), n, idx.data());
            std::printf("J %d %d", k, n);
            for (int i = 0; i < n; ++i) std::printf(" %08X", bits_of(a[i]));
            std::printf(" |");
            for (int i = 0; i < n; ++i) std::printf(" %d", idx[i]);
            std::printf("\n");
        } else if (tag == "K") {
            int k, nb, pb;
            std::string reshex, hex;
            is >> k >> nb >> pb >> reshex;
            kmax = std::max(kmax, k);
            std::vector<cfloat> sp(static_cast<size_t>(nb));
            for (int i = 0; i < nb; ++i) {
                is >> hex;
                const float re = float_of(hex);
                is >> hex;
                sp[i] = cfloat(re, float_of(hex));
            }
            const float r = msk144::msk_peak_offset(sp.data(), nb, pb);
            std::printf("K %d %d %d %08X", k, nb, pb, bits_of(r));
            for (int i = 0; i < nb; ++i) std::printf(" %08X %08X", bits_of(sp[i].real()), bits_of(sp[i].imag()));
            std::printf("\n");
        } else if (tag == "W") {
            int k, nfft;
            std::string lohex, hihex;
            is >> k >> nfft >> lohex >> hihex;
            kmax = std::max(kmax, k);
            int lo, hi;
            bool u;
            msk144::msk_clip_spectrum_window(nfft, double_of(lohex), double_of(hihex), lo, hi, u);
            std::printf("W %d %d %016llX %016llX %d %d %d\n", k, nfft, static_cast<unsigned long long>(dbits_of(double_of(lohex))),
                        static_cast<unsigned long long>(dbits_of(double_of(hihex))), lo, hi, u ? 1 : 0);
        } else if (tag == "SYN") {
            std::puts(ln.c_str());
            int k, s, off, nf, ntol, mask[8], npeaks, nanidx, infidx;
            std::string delfhex, fchex, scalehex;
            is >> k >> s >> off >> nf >> ntol >> delfhex;
            kmax = std::max(kmax, k);
            for (int i = 0; i < 8; ++i) is >> mask[i];
            is >> npeaks >> fchex >> nanidx >> infidx >> scalehex;
            std::vector<cfloat> cdat(static_cast<size_t>(NSPM) * 8, cfloat(0.0f, 0.0f));
            if (nf > 0) make_input(s, off * NSPM, NSPM * nf, nanidx, infidx, float_of(scalehex), cdat.data());
            std::vector<cfloat> c(NSPM, cfloat(7.0f, -7.0f));
            int nk[10];
            for (int& v : nk) v = -1;
            float fest = -1.0f, xmax = -1.0f;
            int nsuc = -1;
            kcase = k;
            msk144::msk144sync(cdat.data(), nf, ntol, float_of(delfhex), mask, npeaks, float_of(fchex), fest, nk, nsuc, xmax, c.data());
            std::printf("X %d %d %08X %08X", k, nsuc, bits_of(fest), bits_of(xmax));
            for (int i = 0; i < npeaks; ++i) std::printf(" %d", nk[i]);
            std::printf("\n");
        } else if (tag == "SPD") {
            std::puts(ln.c_str());
            int k, s, n, ntol, nanidx, infidx;
            std::string fchex, scalehex;
            is >> k >> s >> n >> ntol >> fchex >> nanidx >> infidx >> scalehex;
            kmax = std::max(kmax, k);
            std::vector<cfloat> cbig(static_cast<size_t>(std::max(n, 1)));
            make_input(s, 0, n, nanidx, infidx, float_of(scalehex), cbig.data());
            std::vector<cfloat> ct(NSPM, cfloat(7.0f, -7.0f));
            std::vector<float> sb(144, -3.0f);
            float fret = -1.0f, tret = -1.0f;
            int navg = -1, nsuccess = -1;
            char msg[37];
            std::memset(msg, ' ', 37);
            std::memcpy(msg, "unset", 5);
            kcase = k;
            int nframe = 0;
            const msk144::DecodeSink sink = [&](cfloat* c, float* softbits, char* msgreceived, int& ns) {
                std::vector<float> v(2 * NSPM);
                for (int j = 0; j < NSPM; ++j) {
                    v[2 * j] = c[j].real();
                    v[2 * j + 1] = c[j].imag();
                }
                uint64_t h1, h2;
                hashr(v.data(), 2 * NSPM, h1, h2);
                ++nframe;
                msk144::msk144decodeframe(pk, c, softbits, msgreceived, ns);
                std::printf("F %d %d %llu %llu %d |%.37s|\n", kcase, nframe, static_cast<unsigned long long>(h1),
                            static_cast<unsigned long long>(h2), ns, msgreceived);
            };
            msk144::msk144spd(cbig.data(), n, ntol, nsuccess, msg, float_of(fchex), fret, tret, navg, ct.data(), sb.data(), sink);
            std::printf("X %d %d |%.37s| %08X %08X %d\n", k, nsuccess, msg, bits_of(fret), bits_of(tret), navg);
            prr("U", k, sb.data(), 144);
            prc("V", k, ct.data(), NSPM);
        } else if (tag == "P" || tag == "B" || tag == "R" || tag == "N" || tag == "T" || tag == "D" || tag == "H" || tag == "O" ||
                   tag == "M" || tag == "E" || tag == "G" || tag == "Q" || tag == "Y" || tag == "Z" || tag == "F" || tag == "X" ||
                   tag == "U" || tag == "V" || tag == "END") {
            continue;
        } else {
            std::fprintf(stderr, "msk144-sync-check: unknown line: %s\n", ln.c_str());
            return 1;
        }
    }
    std::printf("END sync %d\n", kmax);
    return 0;
}
