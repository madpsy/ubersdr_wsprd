// FST4W's (240,74) LDPC decoding (fst4w_ldpc.h): WSJT-X's lib/fst4
// decode240_74.f90, fastosd240_74.f90 (with mrbencode74, partial_syndrome,
// nextpat74), encode240_74.f90, get_crc24.f90 and lib/platanh.f90 at
// 567ad29ce, line for line. Arrays keep the Fortran's 1-based indices
// (element 0 unused) where that keeps the two readable side by side.
//
// Steven Franke, K9AN, and Joe Taylor, K1JT, WSJT-X (GNU GPL v3), and so is
// this.
#include "fst4w_ldpc.h"

#include <algorithm>
#include <cmath>
#include <cstring>

#include "fst4w_simd.h"
#include "wsprd_osd.h"  // indexx

namespace fst4w {

namespace {

#include "fst4w_ldpc_tables.inc"

constexpr int N = 240, K = 74, M = N - K;

// encode240_74's generator matrix, from the hex rows (first call: filled
// once, as the Fortran's SAVEd gen).
struct Generator {
    int8_t gen[M][K];
    Generator() {
        std::memset(gen, 0, sizeof gen);
        for (int i = 0; i < M; ++i) {
            for (int j = 1; j <= 19; ++j) {
                const char c = kGen[i][j - 1];
                const int istr = (c >= '0' && c <= '9') ? c - '0' : (c >= 'a' && c <= 'f') ? c - 'a' + 10 : c - 'A' + 10;
                const int ibmax = j == 19 ? 2 : 4;
                for (int jj = 1; jj <= ibmax; ++jj) {
                    const int icol = (j - 1) * 4 + jj;
                    if ((istr >> (4 - jj)) & 1) gen[i][icol - 1] = 1;
                }
            }
        }
    }
};

const Generator& generator() {
    static const Generator g;
    return g;
}

// platanh(x, y): a piecewise-linear atanh.
float platanh(float x) {
    int isign = +1;
    float z = x;
    if (x < 0) {
        isign = -1;
        z = std::fabs(x);
    }
    if (z <= 0.664f) return x / 0.83f;
    if (z <= 0.9217f) return static_cast<float>(isign) * (z - 0.4064f) / 0.322f;
    if (z <= 0.9951f) return static_cast<float>(isign) * (z - 0.8378f) / 0.0524f;
    if (z <= 0.9998f) return static_cast<float>(isign) * (z - 0.9914f) / 0.0012f;
    return static_cast<float>(isign) * 7.0f;
}

// mrbencode74(me, codeword, g2, N, K): g2 is N x k, here column c of g2 is
// row c of rows (g2(:,c) = genmrb(c,:)).
void mrbencode74(const int8_t* me, int8_t* codeword, const int8_t* rows, int k) {
    std::memset(codeword + 1, 0, N);
    for (int i = 1; i <= k; ++i)
        if (me[i] == 1) {
            const int8_t* g = rows + static_cast<size_t>(i - 1) * N;
            for (int r = 1; r <= N; ++r) codeword[r] ^= g[r - 1];
        }
}

// partial_syndrome(me, sp, np, g2, N, K): sp = g2(K:K+np-1, i) summed.
void partial_syndrome(const int8_t* me, int8_t* sp, int np, const int8_t* rows, int k) {
    std::memset(sp + 1, 0, static_cast<size_t>(np));
    for (int i = 1; i <= k; ++i)
        if (me[i] == 1) {
            const int8_t* g = rows + static_cast<size_t>(i - 1) * N;
            for (int r = 1; r <= np; ++r) sp[r] ^= g[k + r - 2];
        }
}

// nextpat74(mi, k, iorder, iflag): the next test pattern of weight iorder.
void nextpat74(int8_t* mi, int k, int iorder, int& iflag) {
    int ind = -1;
    for (int i = 1; i <= k - 1; ++i)
        if (mi[i] == 0 && mi[i + 1] == 1) ind = i;
    if (ind < 0) {
        iflag = ind;
        return;
    }
    int8_t ms[N + 1];
    std::memset(ms, 0, sizeof ms);
    for (int i = 1; i <= ind - 1; ++i) ms[i] = mi[i];
    ms[ind] = 1;
    ms[ind + 1] = 0;
    if (ind + 1 < k) {
        int s = 0;
        for (int i = 1; i <= k; ++i) s += ms[i];
        const int nz = iorder - s;
        for (int i = k - nz + 1; i <= k; ++i) ms[i] = 1;
    }
    std::memcpy(mi + 1, ms + 1, static_cast<size_t>(k));
    for (int i = 1; i <= k; ++i)
        if (mi[i] == 1) {
            iflag = i;
            break;
        }
}

}  // namespace

int get_crc24(const int8_t* mc, int len) {
    static const int8_t p[25] = {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1};
    int8_t r[25];
    std::memcpy(r, mc, 25);
    for (int i = 0; i <= len - 25; ++i) {
        r[24] = mc[i + 24];
        const int8_t r1 = r[0];
        for (int j = 0; j < 25; ++j) r[j] = static_cast<int8_t>((r[j] + r1 * p[j]) % 2);
        const int8_t first = r[0];  // cshift(r, 1)
        for (int j = 0; j < 24; ++j) r[j] = r[j + 1];
        r[24] = first;
    }
    int ncrc = 0;
    for (int j = 0; j < 24; ++j) ncrc = 2 * ncrc + r[j];
    return ncrc;
}

void encode240_74(const int8_t message[74], int8_t codeword[240]) {
    const auto& g = generator().gen;
    int8_t pchecks[M];
    for (int i = 0; i < M; ++i) {
        int nsum = 0;
        for (int j = 0; j < K; ++j) nsum += message[j] * g[i][j];
        pchecks[i] = static_cast<int8_t>(nsum % 2);
    }
    std::memcpy(codeword, message, K);
    std::memcpy(codeword + K, pchecks, M);
}

void Ldpc::decode240_74(const float llr[240], int Keff, int maxosd, int norder, const int8_t apmask[240],
                        int8_t message74[74], int8_t cw[240], int& ntype, int& nharderror, float& dmin) {
    int maxiterations = 30;
    if (Keff == 50) maxiterations = 1;

    float zsave[3][N];  // the Fortran's local zsave(N,3)
    int nosd = 0;
    if (maxosd > 3) maxosd = 3;
    if (maxosd == 0) {  // osd with channel llrs
        nosd = 1;
        std::memcpy(zsave[0], llr, sizeof(float) * N);
    } else if (maxosd > 0) {
        nosd = maxosd;
    } else {  // just bp
        nosd = 0;
    }

    if (maxosd != 0) {
        float toc[M][5], tov[N][3], tanhtoc[M][5], zn[N], zsum[N];
        std::memset(toc, 0, sizeof toc);
        std::memset(tov, 0, sizeof tov);
        std::memset(tanhtoc, 0, sizeof tanhtoc);
        // initialize messages to checks
        for (int j = 0; j < M; ++j)
            for (int i = 0; i < kNrw[j]; ++i) toc[j][i] = llr[kNm[j][i] - 1];

        int ncnt = 0, nclast = 0;
        for (int i = 0; i < N; ++i) zsum[i] = 0.0f;
        for (int iter = 0; iter <= maxiterations; ++iter) {
            // Update bit log likelihood ratios (tov=0 in iteration 0).
            for (int i = 0; i < N; ++i) {
                if (apmask[i] != 1) {
                    float s = 0.0f;
                    for (int kk = 0; kk < kNcw; ++kk) s += tov[i][kk];
                    zn[i] = llr[i] + s;
                } else {
                    zn[i] = llr[i];
                }
            }
            for (int i = 0; i < N; ++i) zsum[i] = zsum[i] + zn[i];
            if (iter > 0 && iter <= maxosd) std::memcpy(zsave[iter - 1], zsum, sizeof zsum);

            // Check to see if we have a codeword (check before we do any iteration).
            for (int i = 0; i < N; ++i) cw[i] = zn[i] > 0.0f ? 1 : 0;
            int ncheck = 0;
            for (int i = 0; i < M; ++i) {
                int synd = 0;
                for (int r = 0; r < kNrw[i]; ++r) synd += cw[kNm[i][r] - 1];
                if (synd % 2 != 0) ++ncheck;
            }
            if (ncheck == 0) {  // we have a codeword - if crc is good, return it
                int8_t m74[74];
                std::memcpy(m74, cw, 74);
                if (get_crc24(m74, 74) == 0) {
                    std::memcpy(message74, cw, 74);
                    int nh = 0;
                    float d = 0.0f;
                    for (int i = 0; i < N; ++i) {
                        const int hdec = llr[i] >= 0 ? 1 : 0;
                        const int nxor = hdec ^ cw[i];
                        nh += nxor;
                        d += static_cast<float>(nxor) * std::fabs(llr[i]);
                    }
                    nharderror = static_cast<int8_t>(nh);  // sum of an integer*1 array: wraps past 127
                    dmin = d;
                    ntype = 1;
                    return;
                }
            }

            if (iter > 0) {  // an early stopping criterion
                const int nd = ncheck - nclast;
                if (nd < 0) ncnt = 0;  // # of unsatisfied parity checks decreased
                else ++ncnt;
                if (ncnt >= 5 && iter >= 10 && ncheck > 15) {
                    nharderror = -1;
                    break;
                }
            }
            nclast = ncheck;

            // Send messages from bits to check nodes
            for (int j = 0; j < M; ++j) {
                for (int i = 0; i < kNrw[j]; ++i) {
                    const int ibj = kNm[j][i];
                    toc[j][i] = zn[ibj - 1];
                    for (int kk = 0; kk < kNcw; ++kk)  // subtract off what the bit had received from the check
                        if (kMn[ibj - 1][kk] == j + 1) toc[j][i] = toc[j][i] - tov[ibj - 1][kk];
                }
            }

            // send messages from check nodes to variable nodes
            for (int i = 0; i < M; ++i)
                for (int r = 0; r < 5; ++r) tanhtoc[i][r] = std::tanh(-toc[i][r] / 2);

            for (int j = 0; j < N; ++j) {
                for (int i = 0; i < kNcw; ++i) {
                    const int ichk = kMn[j][i];  // Mn(:,j) are the checks that include bit j
                    float Tmn = 1.0f;
                    for (int r = 0; r < kNrw[ichk - 1]; ++r)
                        if (kNm[ichk - 1][r] != j + 1) Tmn = Tmn * tanhtoc[ichk - 1][r];
                    const float y = platanh(-Tmn);
                    tov[j][i] = 2 * y;
                }
            }
        }  // bp iterations
    }

    for (int i = 1; i <= nosd; ++i) {
        float zn[N];
        std::memcpy(zn, zsave[i - 1], sizeof zn);
        float dminosd = 0.0f;
        fastosd240_74(zn, Keff, apmask, norder, message74, cw, nharderror, dminosd);
        if (nharderror > 0) {
            int nh = 0;
            float d = 0.0f;
            for (int r = 0; r < N; ++r) {
                const int hdec = llr[r] >= 0 ? 1 : 0;
                const int nxor = hdec ^ cw[r];
                nh += nxor;  // nharderror based on input llrs
                d += static_cast<float>(nxor) * std::fabs(llr[r]);
            }
            nharderror = static_cast<int8_t>(nh);  // integer*1 sum, as above
            dmin = d;
            ntype = 1 + i;
            return;
        }
    }

    ntype = 0;
    nharderror = -1;
}

void Ldpc::fastosd240_74_reference(const float llr[240], int k, const int8_t apmask[240], int ndeep,
                                   int8_t message74[74], int8_t cw_out[240], int& nhardmin, float& dmin) {
    (void)apmask;  // apmaskr is permuted with the columns but never read
    if (first_ || k != ksave_) {  // fill the generator matrix
        // Create generator matrix for partial CRC cascaded with LDPC code:
        // the last 74-k CRC bits cascaded, the first k-50 for detection.
        gen_.assign(static_cast<size_t>(k) * N, 0);
        for (int i = 1; i <= k; ++i) {
            std::memset(message74, 0, 74);
            message74[i - 1] = 1;
            if (i <= 50) {
                const int ncrc24 = get_crc24(message74, 74);
                for (int b = 0; b < 24; ++b) message74[50 + b] = static_cast<int8_t>((ncrc24 >> (23 - b)) & 1);
                for (int b = 51; b <= k; ++b) message74[b - 1] = 0;
            }
            encode240_74(message74, cw_out);
            std::memcpy(&gen_[static_cast<size_t>(i - 1) * N], cw_out, N);
        }
        first_ = false;
        ksave_ = k;
    }

    // 1-based working arrays, as the Fortran's.
    int8_t hdec[N + 1], cw[N + 1], c0[N + 1], ce[N + 1], temp[N + 1], temprow[N + 1];
    int8_t m0[N + 1], me[N + 1], mi[N + 1], sp[33];
    int indices[N + 1], indices2[N + 1], indx[N + 1], nxor[N + 1];
    float absrx[N + 1], rx[N + 1];
    std::vector<int8_t> genmrb(static_cast<size_t>(k) * N);  // genmrb(r,c) = genmrb[(r-1)*N + c-1]
    std::vector<int8_t> tmp(static_cast<size_t>(k) * N);
    auto G = [&](int r, int c) -> int8_t& { return genmrb[static_cast<size_t>(r - 1) * N + c - 1]; };

    const int nswap = 20;
    for (int ibasis = 1; ibasis <= 2; ++ibasis) {
        for (int i = 1; i <= N; ++i) rx[i] = llr[i - 1];
        // Hard decisions on the received word.
        for (int i = 1; i <= N; ++i) hdec[i] = rx[i] >= 0 ? 1 : 0;
        // Use magnitude of received symbols as a measure of reliability.
        for (int i = 1; i <= N; ++i) absrx[i] = std::fabs(llr[i - 1]);
        wsprd::indexx(absrx + 1, N, indx + 1);

        // Re-order the columns of the generator matrix in order of decreasing reliability.
        for (int i = 1; i <= N; ++i) {
            for (int r = 1; r <= k; ++r) G(r, i) = gen_[static_cast<size_t>(r - 1) * N + indx[N + 1 - i] - 1];
            indices[i] = indx[N + 1 - i];
        }

        if (ibasis == 2) {
            for (int i = k - nswap + 1; i <= k; ++i) {
                for (int r = 1; r <= k; ++r) temp[r] = G(r, i);
                for (int r = 1; r <= k; ++r) G(r, i) = G(r, i + nswap);
                for (int r = 1; r <= k; ++r) G(r, i + nswap) = temp[r];
                const int itmp = indices[i];
                indices[i] = indices[i + nswap];
                indices[i + nswap] = itmp;
            }
        }

        // Gaussian elimination: the most reliable received bits in positions
        // 1:k, in order of decreasing reliability (more or less).
        int icol = 1;
        for (int i = 1; i <= N; ++i) indices2[i] = 0;
        int nskipped = 0;
        for (int id = 1; id <= k; ++id) {
            int iflag = 0;
            while (iflag == 0) {
                if (G(id, icol) != 1) {
                    for (int j = id + 1; j <= k; ++j) {
                        if (G(j, icol) == 1) {
                            for (int c = 1; c <= N; ++c) temprow[c] = G(id, c);
                            for (int c = 1; c <= N; ++c) G(id, c) = G(j, c);
                            for (int c = 1; c <= N; ++c) G(j, c) = temprow[c];
                            iflag = 1;
                        }
                    }
                    if (iflag == 0) {  // skip this column
                        ++nskipped;
                        indices2[k + nskipped] = icol;  // put icol where skipped columns go
                        ++icol;                          // look at the next column
                    }
                } else {
                    iflag = 1;
                }
            }
            indices2[id] = icol;
            for (int j = 1; j <= k; ++j)
                if (id != j && G(j, icol) == 1)
                    for (int c = 1; c <= N; ++c) G(j, c) = static_cast<int8_t>(G(id, c) ^ G(j, c));
            ++icol;
        }
        for (int i = k + nskipped + 1; i <= N; ++i) indices2[i] = i;
        // genmrb(1:k,:) = genmrb(1:k,indices2); indices = indices(indices2)
        tmp = genmrb;
        for (int r = 1; r <= k; ++r)
            for (int c = 1; c <= N; ++c) G(r, c) = tmp[static_cast<size_t>(r - 1) * N + indices2[c] - 1];
        {
            int t[N + 1];
            for (int c = 1; c <= N; ++c) t[c] = indices[indices2[c]];
            std::memcpy(indices + 1, t + 1, sizeof(int) * N);
        }
        // g2 = transpose(genmrb): column c of g2 is row c of genmrb.
        const int8_t* rows = genmrb.data();

        // The hard decisions for the k MRB bits define the order 0 message, m0.
        {
            int8_t t[N + 1];
            for (int i = 1; i <= N; ++i) t[i] = hdec[indices[i]];
            std::memcpy(hdec + 1, t + 1, N);
        }
        for (int i = 1; i <= k; ++i) m0[i] = hdec[i];
        for (int i = 1; i <= N; ++i) absrx[i] = std::fabs(llr[indices[i] - 1]);

        mrbencode74(m0, c0, rows, k);
        nhardmin = 0;
        dmin = 0.0f;
        for (int i = 1; i <= N; ++i) {
            nxor[i] = c0[i] ^ hdec[i];
            nhardmin += nxor[i];
            dmin = dmin + static_cast<float>(nxor[i]) * absrx[i];
        }
        const int np = 32;

        std::memcpy(cw + 1, c0 + 1, N);

        if (ndeep != 0) {
            if (ndeep > 4) ndeep = 4;
            int nord = 0, nsyndmax = 0;
            float xlambda = 0.0f;
            if (ndeep == 1) {
                nord = 1;
                xlambda = 0.0f;
                nsyndmax = np;
            } else if (ndeep == 2) {
                nord = 2;
                xlambda = 0.0f;
                nsyndmax = np;
            } else if (ndeep == 3) {
                nord = 3;
                xlambda = 4.0f;
                nsyndmax = 11;
            } else if (ndeep == 4) {
                nord = 4;
                xlambda = 3.4f;
                nsyndmax = 12;
            }

            float s1 = 0.0f, s2 = 0.0f;
            for (int i = 1; i <= k; ++i) s1 = s1 + absrx[i];
            for (int i = k + 1; i <= N; ++i) s2 = s2 + absrx[i];
            const float rho = s1 / (s1 + xlambda * s2);
            float rhodmin = rho * dmin;
            for (int iorder = 1; iorder <= nord; ++iorder) {
                for (int i = 1; i <= k - iorder; ++i) mi[i] = 0;
                for (int i = k - iorder + 1; i <= k; ++i) mi[i] = 1;
                int iflag = k - iorder + 1;
                while (iflag >= 0) {
                    for (int i = 1; i <= k; ++i) me[i] = static_cast<int8_t>(m0[i] ^ mi[i]);
                    float d1 = 0.0f;
                    for (int i = 1; i <= k; ++i) d1 = d1 + static_cast<float>(mi[i]) * absrx[i];
                    if (d1 > rhodmin) break;
                    partial_syndrome(me, sp, np, rows, k);
                    int nwhsp = 0;
                    for (int r = 1; r <= np; ++r) nwhsp += sp[r] ^ hdec[k + r - 1];
                    if (nwhsp <= nsyndmax) {
                        mrbencode74(me, ce, rows, k);
                        float dd = 0.0f;
                        int nh = 0;
                        for (int i = 1; i <= N; ++i) {
                            nxor[i] = ce[i] ^ hdec[i];
                            nh += nxor[i];
                            dd = dd + static_cast<float>(nxor[i]) * absrx[i];
                        }
                        if (dd < dmin) {
                            dmin = dd;
                            rhodmin = rho * dmin;
                            std::memcpy(cw + 1, ce + 1, N);
                            nhardmin = nh;
                        }
                    }
                    // the next test error pattern; iflag goes negative after
                    // the last pattern of weight iorder
                    nextpat74(mi, k, iorder, iflag);
                }
            }
        }
        // Re-order the codeword to [message bits][parity bits] format.
        {
            int8_t t[N + 1];
            std::memcpy(t + 1, cw + 1, N);
            for (int i = 1; i <= N; ++i) cw[indices[i]] = t[i];
        }
        std::memcpy(cw_out, cw + 1, N);
        std::memcpy(message74, cw + 1, 74);
        if (get_crc24(message74, 74) == 0) break;
        nhardmin = -nhardmin;
    }
}

namespace {

// 240 bits as four words: position p (1-based) is bit p-1.
struct B240 {
    uint64_t w[4] = {0, 0, 0, 0};
    bool get(int p) const { return (w[(p - 1) >> 6] >> ((p - 1) & 63)) & 1; }
    void set(int p) { w[(p - 1) >> 6] |= uint64_t(1) << ((p - 1) & 63); }
    void flip(int p) { w[(p - 1) >> 6] ^= uint64_t(1) << ((p - 1) & 63); }
    void operator^=(const B240& o) {
        w[0] ^= o.w[0];
        w[1] ^= o.w[1];
        w[2] ^= o.w[2];
        w[3] ^= o.w[3];
    }
    int count() const {
        return __builtin_popcountll(w[0]) + __builtin_popcountll(w[1]) + __builtin_popcountll(w[2]) +
               __builtin_popcountll(w[3]);
    }
    // Bits p..p+31 (1-based p) as a word.
    uint32_t bits32(int p) const {
        const int b = p - 1, i = b >> 6, s = b & 63;
        uint64_t v = w[i] >> s;
        if (s && i < 3) v |= w[i + 1] << (64 - s);
        return static_cast<uint32_t>(v);
    }
};

inline B240 bx(const B240& a, const B240& b) {
    B240 r;
    r.w[0] = a.w[0] ^ b.w[0];
    r.w[1] = a.w[1] ^ b.w[1];
    r.w[2] = a.w[2] ^ b.w[2];
    r.w[3] = a.w[3] ^ b.w[3];
    return r;
}

// sum(x(i)*absrx(i)), i = 1..N, x the set bits of v: only the set terms, in
// order (the others add +0.0 to a sum of non-negative finite values: no
// change).
inline float sum_set(const B240& v, const float* absrx) {
    float s = 0.0f;
    for (int i = 0; i < 4; ++i) {
        uint64_t x = v.w[i];
        while (x) {
            const int b = __builtin_ctzll(x);
            s = s + absrx[1 + 64 * i + b];
            x &= x - 1;
        }
    }
    return s;
}

// sum_set(v, absrx) < limit, and the sum when it is: the running sum of
// non-negative terms never falls, so once it reaches limit the answer is
// no (the sum itself is then not needed).
inline bool sum_set_below(const B240& v, const float* absrx, float limit, float& sum) {
    float s = 0.0f;
    for (int i = 0; i < 4; ++i) {
        uint64_t x = v.w[i];
        while (x) {
            const int b = __builtin_ctzll(x);
            s = s + absrx[1 + 64 * i + b];
            if (!(s < limit)) return false;
            x &= x - 1;
        }
    }
    sum = s;
    return true;
}

// Lane masks for a byte of bits: kByteMask[b] lane l all ones when bit l of b is set.
struct ByteMasks {
    alignas(32) uint32_t m[256][8];
    ByteMasks() {
        for (int b = 0; b < 256; ++b)
            for (int l = 0; l < 8; ++l) m[b][l] = ((b >> l) & 1) ? 0xffffffffu : 0u;
    }
};
const ByteMasks& byte_masks() {
    static const ByteMasks t;
    return t;
}

// The sum of a240[p-1] over the set bits p of v, in no particular order
// (bytes before first_byte known to be zero).
inline float masked_sum(const B240& v, const float* a240, int first_byte) {
    const ByteMasks& bm = byte_masks();
    simd::V8 acc = simd::zero();
    const uint8_t* bytes = reinterpret_cast<const uint8_t*>(v.w);  // little-endian: byte i is positions 8i+1..8i+8
    for (int i = first_byte; i < 30; ++i) acc = simd::add_masked_byte(acc, a240 + 8 * i, bytes[i], &bm.m[0][0]);
    return simd::hsum(acc);
}

// What the order-1..nord pattern search reads and updates.
struct Search {
    const B240* rows;  // rows[1..k]: the reduced generator (g2's columns)
    const B240* c0;    // the order-0 codeword
    const B240* hdec;
    const B240* parity;
    const float* absrx;  // 1-based, in the reordered positions
    const float* a240;   // the same, 0-based and 32-byte aligned (a240[p-1] = absrx[p])
    int k, nord, nsyndmax;
    uint32_t hsyn;  // hdec(k:k+31)
    float rho;
    float* rhodmin;
    float* dmin;
    B240* cw;
    int* nhardmin;
};

// fastosd240_74's search over test patterns of weight 1..nord, in
// nextpat74's order, a pattern held as its positions p[0] < ... < p[w-1]:
// nextpat74 finds the last i with mi(i) = 0 and mi(i+1) = 1, sets mi(i),
// clears mi(i+1) and puts the rest of the weight at the end (k-nz+1..k);
// that i is p[j]-1 for the last j whose left neighbour is free, so: p[j]
// moves left by one and p[j+1..] go to the end.
void search_patterns(Search& S) {
    const int k = S.k;
    // The rows' and the order-0 codeword's bits k..k+31 (the partial
    // syndrome's), so a pattern's syndrome needs no 240-bit codeword.
    uint32_t rowsyn[75];
    for (int r = 1; r <= k; ++r) rowsyn[r] = S.rows[r].bits32(k);
    const uint32_t c0syn = S.c0->bits32(k) ^ S.hsyn;
    for (int w = 1; w <= S.nord; ++w) {
        int p[4];
        for (int t = 0; t < w; ++t) p[t] = k - w + 1 + t;  // mi(k-iorder+1:k) = 1
        for (;;) {
            // d1 = sum(mi(1:k)*absrx(1:k)): the set positions, in order
            // (the zero terms add +0.0 to a sum of non-negative values).
            float d1 = 0.0f;
            uint32_t syn = c0syn;
            for (int t = 0; t < w; ++t) {
                d1 = d1 + S.absrx[p[t]];
                syn ^= rowsyn[p[t]];
            }
            if (d1 > *S.rhodmin) break;
            // partial syndrome against hdec(k:k+np-1)
            const int nwhsp = __builtin_popcount(syn);
            if (nwhsp <= S.nsyndmax && d1 < *S.dmin) {
                B240 ce = *S.c0;
                for (int t = 0; t < w; ++t) ce ^= S.rows[p[t]];
                // dd = sum(nxor*absrx): positions 1..k of nxor are mi's (the
                // reduced generator is the identity there and hdec(1:k) is
                // m0), so the sum runs on from d1 over the parity positions,
                // in order; it stops once it can no longer be below dmin.
                const B240 nxe = bx(ce, *S.hdec);
                B240 par = nxe;
                for (int q = 0; q < 4; ++q) par.w[q] &= S.parity->w[q];
                // A quick no, exact: P, the parity terms summed in any order
                // (eight lanes), is at most 1.5e-5 above their exact sum, and
                // the sequential sum from d1 at least 1.5e-5 below it
                // (positive terms, at most 244 roundings): once
                // (d1 + P*(1 - 1e-4))*(1 - 3e-5) reaches dmin, so does the
                // sequential sum.
                if ((d1 + masked_sum(par, S.a240, k >> 3) * (1.0f - 1e-4f)) * (1.0f - 3e-5f) < *S.dmin) {
                    float sum = d1;
                    bool below = true;
                    for (int q = 0; q < 4 && below; ++q) {
                        uint64_t x = par.w[q];
                        while (x) {
                            sum = sum + S.absrx[1 + 64 * q + __builtin_ctzll(x)];
                            if (!(sum < *S.dmin)) {
                                below = false;
                                break;
                            }
                            x &= x - 1;
                        }
                    }
                    if (below) {
                        *S.dmin = sum;
                        *S.rhodmin = S.rho * *S.dmin;
                        *S.cw = ce;
                        *S.nhardmin = nxe.count();
                    }
                }
            }
            // nextpat74: the last j whose left neighbour is free
            int j = w - 1;
            while (j >= 0 && (j == 0 ? p[0] < 2 : p[j - 1] == p[j] - 1)) --j;
            if (j < 0) break;
            const int q = p[j];
            p[j] = q - 1;
            if (q < k)  // ind+1 < k: the rest of the weight at the end
                for (int t = j + 1; t < w; ++t) p[t] = k - (w - 1 - t);
        }
    }
}

}  // namespace

void Ldpc::fastosd240_74(const float llr[240], int k, const int8_t apmask[240], int ndeep, int8_t message74[74],
                         int8_t cw_out[240], int& nhardmin, float& dmin) {
    for (int i = 0; i < N; ++i)
        if (!std::isfinite(llr[i])) {
            fastosd240_74_reference(llr, k, apmask, ndeep, message74, cw_out, nhardmin, dmin);
            return;
        }
    if (k < 1 || k > 74) {
        fastosd240_74_reference(llr, k, apmask, ndeep, message74, cw_out, nhardmin, dmin);
        return;
    }
    // The generator for partial CRC cascaded with the LDPC code, per k (the
    // reference rebuilds it whenever k changes; the same rows).
    std::vector<uint64_t>& gb = genbits_[k];
    if (gb.empty()) {
        gb.assign(static_cast<size_t>(k) * 4, 0);
        int8_t m74[74], c[240];
        for (int i = 1; i <= k; ++i) {
            std::memset(m74, 0, 74);
            m74[i - 1] = 1;
            if (i <= 50) {
                const int ncrc24 = get_crc24(m74, 74);
                for (int b = 0; b < 24; ++b) m74[50 + b] = static_cast<int8_t>((ncrc24 >> (23 - b)) & 1);
                for (int b = 51; b <= k; ++b) m74[b - 1] = 0;
            }
            encode240_74(m74, c);
            for (int p = 1; p <= N; ++p)
                if (c[p - 1]) gb[static_cast<size_t>(i - 1) * 4 + ((p - 1) >> 6)] |= uint64_t(1) << ((p - 1) & 63);
        }
    }
    auto genbit = [&](int r, int p) { return (gb[static_cast<size_t>(r - 1) * 4 + ((p - 1) >> 6)] >> ((p - 1) & 63)) & 1; };

    int indx[N + 1], indices[N + 1], indices2[N + 1];
    float absrx[N + 1];
    B240 rows[75];
    const int nswap = 20;
    for (int ibasis = 1; ibasis <= 2; ++ibasis) {
        for (int i = 1; i <= N; ++i) absrx[i] = std::fabs(llr[i - 1]);
        wsprd::indexx(absrx + 1, N, indx + 1);
        for (int i = 1; i <= N; ++i) indices[i] = indx[N + 1 - i];
        if (ibasis == 2)
            for (int i = k - nswap + 1; i <= k; ++i) std::swap(indices[i], indices[i + nswap]);
        // genmrb(r, i) = gen(r, indices(i))
        for (int r = 1; r <= k; ++r) {
            rows[r] = B240();
            for (int i = 1; i <= N; ++i)
                if (genbit(r, indices[i])) rows[r].set(i);
        }
        // Gaussian elimination, as the reference does it (its swaps and all).
        int icol = 1, nskipped = 0;
        for (int i = 1; i <= N; ++i) indices2[i] = 0;
        for (int id = 1; id <= k; ++id) {
            int iflag = 0;
            while (iflag == 0) {
                if (!rows[id].get(icol)) {
                    for (int j = id + 1; j <= k; ++j)
                        if (rows[j].get(icol)) {
                            std::swap(rows[id], rows[j]);
                            iflag = 1;
                        }
                    if (iflag == 0) {
                        ++nskipped;
                        indices2[k + nskipped] = icol;
                        ++icol;
                    }
                } else {
                    iflag = 1;
                }
            }
            indices2[id] = icol;
            for (int j = 1; j <= k; ++j)
                if (id != j && rows[j].get(icol)) rows[j] ^= rows[id];
            ++icol;
        }
        for (int i = k + nskipped + 1; i <= N; ++i) indices2[i] = i;
        // genmrb(1:k,:) = genmrb(1:k,indices2); indices = indices(indices2)
        for (int r = 1; r <= k; ++r) {
            B240 t;
            for (int c = 1; c <= N; ++c)
                if (rows[r].get(indices2[c])) t.set(c);
            rows[r] = t;
        }
        {
            int t[N + 1];
            for (int c = 1; c <= N; ++c) t[c] = indices[indices2[c]];
            std::memcpy(indices + 1, t + 1, sizeof(int) * N);
        }
        // Hard decisions and reliabilities in the new order.
        B240 hdec;
        for (int i = 1; i <= N; ++i) {
            if (llr[indices[i] - 1] >= 0) hdec.set(i);
            absrx[i] = std::fabs(llr[indices[i] - 1]);
        }
        // The order 0 codeword: the k MRB hard decisions encoded.
        B240 c0;
        for (int r = 1; r <= k; ++r)
            if (hdec.get(r)) c0 ^= rows[r];
        B240 nx = bx(c0, hdec);
        nhardmin = nx.count();
        dmin = sum_set(nx, absrx);
        const int np = 32;
        B240 cw = c0;

        if (ndeep != 0) {
            if (ndeep > 4) ndeep = 4;
            int nord = 0, nsyndmax = 0;
            float xlambda = 0.0f;
            if (ndeep == 1) {
                nord = 1;
                nsyndmax = np;
            } else if (ndeep == 2) {
                nord = 2;
                nsyndmax = np;
            } else if (ndeep == 3) {
                nord = 3;
                xlambda = 4.0f;
                nsyndmax = 11;
            } else if (ndeep == 4) {
                nord = 4;
                xlambda = 3.4f;
                nsyndmax = 12;
            }
            float s1 = 0.0f, s2 = 0.0f;
            for (int i = 1; i <= k; ++i) s1 = s1 + absrx[i];
            for (int i = k + 1; i <= N; ++i) s2 = s2 + absrx[i];
            const float rho = s1 / (s1 + xlambda * s2);
            float rhodmin = rho * dmin;
            const uint32_t hsyn = hdec.bits32(k);  // hdec(k:k+np-1)
            // The parity positions k+1..N of a 240-bit word.
            B240 parity;
            for (int p = k + 1; p <= N; ++p) parity.set(p);
            alignas(32) float a240[240];
            for (int p = 1; p <= N; ++p) a240[p - 1] = absrx[p];
            Search srch{rows, &c0, &hdec, &parity, absrx, a240, k, nord, nsyndmax, hsyn, rho, &rhodmin, &dmin, &cw, &nhardmin};
            search_patterns(srch);
        }
        // Re-order the codeword to [message bits][parity bits] format.
        for (int i = 1; i <= N; ++i) cw_out[indices[i] - 1] = static_cast<int8_t>(cw.get(i));
        std::memcpy(message74, cw_out, 74);
        if (get_crc24(message74, 74) == 0) break;
        nhardmin = -nhardmin;
    }
}

}  // namespace fst4w
