// MSK144's real-time decoder (msk144.h): WSJT-X's mskrtd.f90 and
// decode_msk144.f90 at 567ad29ce, line for line. The Fortran's 1-based loops
// are shifted to 0-based with the Fortran's index in a comment where it
// helps.
//
// Joe Taylor, K1JT, Steven Franke, K9AN, and the WSJT-X team (GNU GPL v3),
// and so is this.
#include "msk144.h"

#include <cmath>
#include <cstring>

#include "fst4w_jt9.h"  // fortran_f, fortran_i
#include "ft8_maths.h"
#include "msk144_decode.h"
#include "msk144_sync.h"

namespace msk144 {

namespace {

constexpr int NZ = 7168;       // block size
constexpr int NSPM = 864;      // samples per message frame
constexpr int NFFT1 = 8192;    // FFT size for making the analytic signal
constexpr int NPATTERNS = 4;   // number of frame averaging patterns to try

inline int nint_f(float x) { return static_cast<int>(llroundf(x)); }

std::string blank37() { return std::string(37, ' '); }

}  // namespace

// ---- mskrtd.f90:243-261 (update_recent_shmsgs)
static void update_recent_shmsgs(const std::string& message, std::array<std::string, RtdState::NSHMEM>& msgs) {
    bool seen = false;
    for (const std::string& m : msgs)
        if (m == message) seen = true;
    if (!seen) {
        for (int i = RtdState::NSHMEM; i >= 2; --i) msgs[static_cast<size_t>(i - 1)] = msgs[static_cast<size_t>(i - 2)];
        msgs[0] = message;
    }
}

// ---- mskrtd.f90:1-241 (btrain false; msk40spd at lines 129-134, the bswl branch at 221-237)
std::string mskrtd(RtdState& st, packjt77::Packjt77& pk, const int16_t* id2, int nutc0, float tsec, int ntol,
                   int nrxfreq, int ndepth, const std::string& mycall12, const std::string& hiscall12, bool bshmsg,
                   bool bswl) {
    // data iavpatterns/ 1,1,1,1,0,0,0,0, 0,0,1,1,1,1,0,0, 1,1,1,1,1,0,0,0, 1,1,1,1,1,1,1,0/, (8,NPATTERNS)
    static const int iavpatterns[NPATTERNS][8] = {{1, 1, 1, 1, 0, 0, 0, 0},
                                                  {0, 0, 1, 1, 1, 1, 0, 0},
                                                  {1, 1, 1, 1, 1, 0, 0, 0},
                                                  {1, 1, 1, 1, 1, 1, 1, 0}};
    // data xmc/2.0,4.5,2.5,3.5/     !Used to set time at center of averaging mask
    static const float xmc[NPATTERNS] = {2.0f, 4.5f, 2.5f, 3.5f};

    // mycall and hiscall are character*12
    std::string mycall = mycall12, hiscall = hiscall12;
    mycall.resize(12, ' ');
    hiscall.resize(12, ' ');

    cfloat c[NSPM];
    cfloat ct[NSPM];
    float softbits[144];
    char msgreceived[37];
    int iavmask[8];
    int npkloc[10] = {};
    float d[NFFT1];
    float pow[8];
    int ndecodesuccess = 0;
    float fest = 0.0f, tdec = 0.0f;
    int navg = 0;
    bool bshdecode = false;
    std::string line;

    if (st.first) {  // lines 59-80
        st.tsec0 = tsec;
        st.nutc00 = nutc0;
        st.pnoise = -1.0f;
        pk.blank_recent_calls();  // recent_calls(i)(1:13)=' '
        st.msglast = blank37();
        for (std::string& m : st.recent_shmsgs) m = blank37();  // recent_shmsgs(i)(1:37)=' '
        st.msglastswl = blank37();
        st.nsnrlast = -99;
        st.nsnrlastswl = -99;
        pk.set_calls(mycall, hiscall);  // mycall13=' '; mycall13(1:12)=mycall(1:12), and dxcall13
        st.first = false;
    }

    const float fc = static_cast<float>(nrxfreq);

    // Reset if mycall or dxcall changes
    if (pk.mycall13().substr(0, 12) != mycall || pk.dxcall13().substr(0, 12) != hiscall) st.first = true;

    // Dupe checking setup
    if (st.nutc00 != nutc0 || tsec < st.tsec0) {  // reset dupe checker
        st.msglast = blank37();
        st.msglastswl = blank37();
        st.nsnrlast = -99;
        st.nsnrlastswl = -99;
        st.nutc00 = nutc0;
    }

    const float tframe = static_cast<float>(NSPM) / 12000.0f;
    std::memset(msgreceived, ' ', 37);
    const int max_iterations = 10;
    (void)max_iterations;
    // d(1:NZ)=id2; rms=sqrt(sum(d(1:NZ)*d(1:NZ))/NZ)
    for (int i = 0; i < NZ; ++i) d[i] = static_cast<float>(id2[i]);
    float sumsq = 0.0f;
    for (int i = 0; i < NZ; ++i) sumsq = sumsq + d[i] * d[i];
    const float rms = std::sqrt(sumsq / static_cast<float>(NZ));
    if (rms < 1.0f) {
        st.tsec0 = tsec;  // go to 999
        return line;
    }
    const float fac = 1.0f / rms;
    for (int i = 0; i < NZ; ++i) d[i] = fac * d[i];
    for (int i = NZ; i < NFFT1; ++i) d[i] = 0.0f;
    const bool bvar = true;  // btrain false
    analytic(d, NZ, NFFT1, st.cdat.data(), st.pcoeffs, bvar);
    cfloat* cdat = st.cdat.data();

    // Calculate average power for each frame and for the entire block.
    // pow(i)=real(dot_product(cdat(ib:ie),cdat(ib:ie)))*rms**2
    float pmax = -99.0f;
    for (int i = 1; i <= 8; ++i) {
        const int ib = (i - 1) * NSPM;
        float s = 0.0f;
        for (int j = 0; j < NSPM; ++j) {
            // conjg(a)*a, real part: ar*ar - (-ai)*ai
            const float ar = cdat[ib + j].real(), ai = cdat[ib + j].imag();
            s = s + (ar * ar - (-ai) * ai);
        }
        pow[i - 1] = s * (rms * rms);
        pmax = pmax > pow[i - 1] || std::isnan(pow[i - 1]) ? pmax : pow[i - 1];  // pmax=max(pmax,pow(i))
    }
    float psum = 0.0f;
    for (int i = 0; i < 8; ++i) psum = psum + pow[i];
    const float pavg = psum / 8.0f;

    // Short ping decoder: msk144spd on the 8 frames, with the frame decoder
    const int np = 8 * NSPM;
    const DecodeSink sink = [&pk](cfloat* cc, float* sb, char* mg, int& ns) { msk144decodeframe(pk, cc, sb, mg, ns); };
    msk144spd(cdat, np, ntol, ndecodesuccess, msgreceived, fc, fest, tdec, navg, ct, softbits, sink);
    bool decoded = false;
    if (ndecodesuccess == 0 && (bshmsg || bswl)) {
        msk40spd(st.m40, pk, cdat, np, ntol, mycall, hiscall, bswl, ndecodesuccess, msgreceived, fc, fest, tdec, navg);
        if (ndecodesuccess >= 1) bshdecode = true;
    }
    if (ndecodesuccess >= 1) {
        tdec = tsec + tdec;
        decoded = true;
    } else {
        // Fast - short ping decoder only. Normal - try 4-frame averages. Deep - try 4-, 5- and 7-frame averages.
        int npat = NPATTERNS;
        if (ndepth == 1) npat = 0;
        if (ndepth == 2) npat = 2;
        for (int iavg = 1; iavg <= npat && !decoded; ++iavg) {
            for (int i = 0; i < 8; ++i) iavmask[i] = iavpatterns[iavg - 1][i];
            navg = 0;
            for (int i = 0; i < 8; ++i) navg += iavmask[i];
            const float deltaf = 10.0f / static_cast<float>(navg);  // search increment for frequency sync
            const int npeaks = 2;
            int nsyncsuccess = 0;
            float xmax = 0.0f;
            msk144sync(cdat, 8, ntol, deltaf, iavmask, npeaks, fc, fest, npkloc, nsyncsuccess, xmax, c);
            if (nsyncsuccess == 0) continue;

            for (int ipk = 1; ipk <= npeaks && !decoded; ++ipk) {
                for (int is = 1; is <= 3 && !decoded; ++is) {  // Slicer dither
                    int ic0 = npkloc[ipk - 1];
                    if (is == 2) ic0 = std::max(1, ic0 - 1);
                    if (is == 3) ic0 = std::min(NSPM, ic0 + 1);
                    // ct=cshift(c,ic0-1)
                    for (int i = 0; i < NSPM; ++i) ct[i] = c[(i + ic0 - 1) % NSPM];
                    msk144decodeframe(pk, ct, softbits, msgreceived, ndecodesuccess);
                    if (ndecodesuccess > 0) {
                        tdec = tsec + xmc[iavg - 1] * tframe;
                        decoded = true;
                    }
                }
            }
        }
    }

    if (!decoded) {
        std::memset(msgreceived, ' ', 37);
        // no decode - update noise level used for calculating displayed snr.
        if (st.pnoise < 0.0f) {  // initialize noise level
            st.pnoise = pavg;
        } else if (pavg > st.pnoise) {  // noise level is slow to rise
            st.pnoise = 0.9f * st.pnoise + 0.1f * pavg;
        } else if (pavg < st.pnoise) {  // and quick to fall
            st.pnoise = pavg;
        }
        st.tsec0 = tsec;  // go to 999
        return line;
    }

    // 900: Successful decode - estimate snr
    float snr0;
    if (st.pnoise > 0.0f) {
        snr0 = 10.0f * ft8::log10f_glibc(pmax / st.pnoise - 1.0f);
    } else {
        snr0 = 0.0f;
    }
    int nsnr = nint_f(snr0);

    if (!bshdecode) {
        int ncorrected;
        float eyeopening;
        msk144signalquality(pk, ct, snr0, fest, tdec, softbits, msgreceived, hiscall.c_str(), false, "", ncorrected,
                            eyeopening, st.pcoeffs.data());
    }

    const char* decsym = " &  ";

    if (nsnr < -8) nsnr = -8;
    if (nsnr > 24) nsnr = 24;

    // Dupe check.
    const std::string msgstr(msgreceived, 37);
    const bool bflag = ndecodesuccess == 1 && (msgstr != st.msglast || nsnr > st.nsnrlast || tsec < st.tsec0);
    if (bflag) {
        st.msglast = msgstr;
        st.nsnrlast = nsnr;
        if (!bshdecode) update_msk40_hasharray(st.m40, pk);
        // 1021 format(i6.6,i4,f5.1,i5,a4,a37,a1)
        line = fst4w::fortran_i(nutc0, 6, 6) + fst4w::fortran_i(nsnr, 4) + fst4w::fortran_f(tdec, 5, 1) +
               fst4w::fortran_i(nint_f(fest), 5) + decsym + msgstr;
    } else if (bswl && ndecodesuccess >= 2) {
        bool seenb4 = false;
        for (const std::string& m : st.recent_shmsgs)
            if (msgstr == m) seenb4 = true;
        update_recent_shmsgs(msgstr, st.recent_shmsgs);
        const bool bflag2 = seenb4 && (msgstr != st.msglastswl || nsnr > st.nsnrlastswl || tsec < st.tsec0) && nsnr > -6;
        if (bflag2) {
            st.msglastswl = msgstr;
            st.nsnrlastswl = nsnr;
            line = fst4w::fortran_i(nutc0, 6, 6) + fst4w::fortran_i(nsnr, 4) + fst4w::fortran_f(tdec, 5, 1) +
                   fst4w::fortran_i(nint_f(fest), 5) + decsym + msgstr;
        }
    }
    st.tsec0 = tsec;  // 999
    return line;
}

// ---- decode_msk144.f90:1-88
int decode_msk144(RtdState& st, packjt77::Packjt77& pk, const int16_t* audio, long nsamples, int nutc, int ntr, int ntol,
                  int nfqso, int ndepth, const std::string& mycall12, const std::string& hiscall12, bool bshmsg, bool bswl,
                  std::ostream& out) {
    constexpr int SAMPLING_RATE = 12000, BLOCK_SIZE = 7168, STEP_SIZE = BLOCK_SIZE / 2;
    int message_count = 0;
    // decode in 0.3s blocks
    const long sample_count = static_cast<long>(ntr) * SAMPLING_RATE;
    std::array<int16_t, BLOCK_SIZE> buffer;
    for (long position = 1; position <= sample_count - BLOCK_SIZE + 1; position += STEP_SIZE) {
        // buffer = audio_samples(position:position+BLOCK_SIZE-1): what was never read is zero
        for (int i = 0; i < BLOCK_SIZE; ++i) {
            const long j = position - 1 + i;
            buffer[static_cast<size_t>(i)] = j < nsamples ? audio[j] : static_cast<int16_t>(0);
        }
        const float tsec = static_cast<float>(position) / static_cast<float>(SAMPLING_RATE);
        std::string line = mskrtd(st, pk, buffer.data(), nutc, tsec, ntol, nfqso, ndepth, mycall12, hiscall12, bshmsg,
                                    bswl);
        // line(1:1).ne.char(0): a line was made. line = line(1:index(line,char(0))-1), write(*,1001) line, a80
        if (!line.empty()) {
            line.push_back('\0');  // the format's a1 char(0), then blank to the 80 of `line`
            line.resize(80, ' ');
            line.resize(line.find('\0'));
            line.resize(80, ' ');
            out << line << "\n";
            ++message_count;
        }
    }
    return message_count;
}

}  // namespace msk144
