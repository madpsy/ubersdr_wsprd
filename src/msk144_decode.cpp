// MSK144's frame decoder and signal quality (msk144_decode.h): WSJT-X's
// msk144decodeframe.f90, msk144signalquality.f90, genmsk_128_90.f90 and
// fmtmsg.f90 at 567ad29ce, line for line. Arrays keep the Fortran's order;
// the 1-based loops of the Fortran are shifted to 0-based here with the
// Fortran's index in a comment where it helps.
//
// Joe Taylor, K1JT, Steven Franke, K9AN, and the WSJT-X team (GNU GPL v3),
// and so is this.
#include "msk144_decode.h"

#include <cmath>
#include <cstdint>
#include <cstring>

#include "ft8_maths.h"
#include "msk144_ldpc.h"

namespace msk144 {

DecodeTrace decode_trace;

namespace {

constexpr int NSPM = 864;

// complex products as gfortran writes them (-fcx-fortran-rules): the plain
// formula, no NaN recovery.
inline cfloat cmul(cfloat a, cfloat b) {
    const float ar = a.real(), ai = a.imag(), br = b.real(), bi = b.imag();
    return cfloat(ar * br - ai * bi, ar * bi + ai * br);
}
inline cfloat conjg(cfloat a) { return cfloat(a.real(), -a.imag()); }

// ---- msk144decodeframe's constants, made at the first call (msk144decodeframe.f90:22-46): pp, the half-sine
// pulse, and cb, the sync word waveforms. angle=(i-1)*pi/12.0 with pi a real*8 (4d0*datan(1d0)): a double
// rounded to the float angle; pp(i)=sin(angle). The reference's compiler vectorised the loop: libmvec's
// _ZGVbN4v_sinf for i = 2..9 (two calls), i = 1 is sin 0, i = 10..12 folded to constants (objdump of the
// release build's libwsjt_fort_omp.a, msk144decodeframe.f90.o: the call sites, the angle vectors, and
// the stores of 0x3f3504f3, 0x3f000001, 0x3e8483f1); the constants are the correctly rounded sinf of the float
// angles of i = 10, 11, 12. The oracle's KP and KB lines check them. s8 becomes 2*s8-1.
struct DecodeConsts {
    float pp[12];
    cfloat cb[42];
    int s8[8];
    DecodeConsts() {
        const double pi = 4.0 * 0.78539816339744830961566084581988;
        float ang[12];
        for (int i = 1; i <= 12; ++i) ang[i - 1] = static_cast<float>((static_cast<double>(i - 1) * pi) / 12.0);
        pp[0] = 0.0f;
        for (int g = 0; g < 2; ++g) ft8::sinf4_glibc(ang + 1 + 4 * g, pp + 1 + 4 * g);
        const uint32_t folded[3] = {0x3F3504F3u, 0x3F000001u, 0x3E8483F1u};  // sin of the angles of i = 10, 11, 12
        for (int i = 0; i < 3; ++i) std::memcpy(&pp[9 + i], &folded[i], 4);

        const int s8d[8] = {0, 1, 1, 1, 0, 0, 1, 0};
        for (int i = 0; i < 8; ++i) s8[i] = 2 * s8d[i] - 1;
        float cbi[42], cbq[42];
        for (int i = 0; i < 6; ++i) cbq[i] = pp[6 + i] * static_cast<float>(s8[0]);       // cbq(1:6)=pp(7:12)*s8(1)
        for (int i = 0; i < 12; ++i) cbq[6 + i] = pp[i] * static_cast<float>(s8[2]);      // cbq(7:18)=pp*s8(3)
        for (int i = 0; i < 12; ++i) cbq[18 + i] = pp[i] * static_cast<float>(s8[4]);     // cbq(19:30)=pp*s8(5)
        for (int i = 0; i < 12; ++i) cbq[30 + i] = pp[i] * static_cast<float>(s8[6]);     // cbq(31:42)=pp*s8(7)
        for (int i = 0; i < 12; ++i) cbi[i] = pp[i] * static_cast<float>(s8[1]);          // cbi(1:12)=pp*s8(2)
        for (int i = 0; i < 12; ++i) cbi[12 + i] = pp[i] * static_cast<float>(s8[3]);     // cbi(13:24)=pp*s8(4)
        for (int i = 0; i < 12; ++i) cbi[24 + i] = pp[i] * static_cast<float>(s8[5]);     // cbi(25:36)=pp*s8(6)
        for (int i = 0; i < 6; ++i) cbi[36 + i] = pp[i] * static_cast<float>(s8[7]);      // cbi(37:42)=pp(1:6)*s8(8)
        for (int i = 0; i < 42; ++i) cb[i] = cfloat(cbi[i], cbq[i]);
        if (decode_trace.consts) decode_trace.consts(pp, cb);
    }
};

}  // namespace

Packjt77& msk144_packjt77() {
    static Packjt77 pk;
    return pk;
}
void msk144_reset_packjt77() { msk144_packjt77() = Packjt77(); }

// ---- msk144decodeframe.f90:1-112
void msk144decodeframe(Packjt77& pk, cfloat* c, float* softbits, char* msgreceived, int& nsuccess) {
    static const DecodeConsts K;  // first=.true. branch, once
    const float* pp = K.pp;
    const cfloat* cb = K.cb;
    const int* s8 = K.s8;

    nsuccess = 0;
    std::memset(msgreceived, ' ', 37);  // msgreceived=' '

    // Estimate carrier phase.
    // cca=sum(c(1:1+41)*conjg(cb)); ccb=sum(c(1+56*6:1+56*6+41)*conjg(cb)): each a running sum from zero, in order
    cfloat cca(0.0f, 0.0f), ccb(0.0f, 0.0f);
    for (int i = 0; i < 42; ++i) {
        const cfloat p = cmul(c[i], conjg(cb[i]));
        cca = cfloat(cca.real() + p.real(), cca.imag() + p.imag());
    }
    for (int i = 0; i < 42; ++i) {
        const cfloat p = cmul(c[56 * 6 + i], conjg(cb[i]));
        ccb = cfloat(ccb.real() + p.real(), ccb.imag() + p.imag());
    }
    // (cfac=ccb*conjg(cca) is overwritten below)
    const cfloat ccasum(cca.real() + ccb.real(), cca.imag() + ccb.imag());  // cca+ccb
    const float phase0 = ft8::atan2f_glibc(ccasum.imag(), ccasum.real());
    if (decode_trace.phase) decode_trace.phase(cca, ccb, phase0);

    // Remove phase error - want constellation rotated so that sample points lie on I/Q axes
    float sphase, cphase;
    ft8::sincosf_glibc(phase0, &sphase, &cphase);  // cfac=cmplx(cos(phase0),sin(phase0))
    const cfloat cfac(cphase, sphase);
    const cfloat cfacc = conjg(cfac);
    for (int i = 0; i < NSPM; ++i) c[i] = cmul(c[i], cfacc);  // c=c*conjg(cfac)
    if (decode_trace.corrected) decode_trace.corrected(cfac, c);

    // matched filter -
    {
        float s1 = 0.0f, s2 = 0.0f;
        for (int j = 0; j < 6; ++j) s1 = s1 + c[j].imag() * pp[6 + j];                 // imag(c(1:6))*pp(7:12)
        for (int j = 0; j < 6; ++j) s2 = s2 + c[864 - 6 + j].imag() * pp[j];           // imag(c(864-5:864))*pp(1:6)
        softbits[0] = s1 + s2;
        float s3 = 0.0f;
        for (int j = 0; j < 12; ++j) s3 = s3 + c[j].real() * pp[j];                    // real(c(1:12))*pp
        softbits[1] = s3;
    }
    for (int i = 2; i <= 72; ++i) {
        float sa = 0.0f, sb = 0.0f;
        // imag(c(1+(i-1)*12-6:1+(i-1)*12+5))*pp: 0-based start (i-1)*12-6
        for (int j = 0; j < 12; ++j) sa = sa + c[(i - 1) * 12 - 6 + j].imag() * pp[j];
        // real(c(7+(i-1)*12-6:7+(i-1)*12+5))*pp: 0-based start (i-1)*12
        for (int j = 0; j < 12; ++j) sb = sb + c[(i - 1) * 12 + j].real() * pp[j];
        softbits[2 * i - 2] = sa;  // softbits(2*i-1)
        softbits[2 * i - 1] = sb;  // softbits(2*i)
    }
    if (decode_trace.raw) decode_trace.raw(softbits);

    // sync word hard error weight is used as a discriminator for
    // frames that have reasonable probability of decoding
    int hardbits[144];
    for (int i = 0; i < 144; ++i) hardbits[i] = softbits[i] >= 0.0f ? 1 : 0;
    int sum1 = 0, sum2 = 0;
    for (int i = 0; i < 8; ++i) sum1 += (2 * hardbits[i] - 1) * s8[i];
    for (int i = 0; i < 8; ++i) sum2 += (2 * hardbits[56 + i] - 1) * s8[i];
    const int nbadsync1 = (8 - sum1) / 2;
    const int nbadsync2 = (8 - sum2) / 2;
    const int nbadsync = nbadsync1 + nbadsync2;
    if (decode_trace.sync) decode_trace.sync(nbadsync1, nbadsync2);
    if (nbadsync > 4) return;

    // normalize the softsymbols before submitting to decoder
    float ssum = 0.0f, s2sum = 0.0f;
    for (int i = 0; i < 144; ++i) ssum = ssum + softbits[i];
    for (int i = 0; i < 144; ++i) s2sum = s2sum + softbits[i] * softbits[i];
    const float sav = ssum / 144.0f;
    const float s2av = s2sum / 144.0f;
    const float ssig = std::sqrt(s2av - sav * sav);
    if (decode_trace.norm) decode_trace.norm(sav, s2av, ssig);
    for (int i = 0; i < 144; ++i) softbits[i] = softbits[i] / ssig;

    const float sigma = 0.60f;
    float llr[128];
    for (int i = 0; i < 48; ++i) llr[i] = softbits[8 + i];        // llr(1:48)=softbits(9:9+47)
    for (int i = 0; i < 80; ++i) llr[48 + i] = softbits[64 + i];  // llr(49:128)=softbits(65:65+80-1)
    for (int i = 0; i < 128; ++i) llr[i] = 2.0f * llr[i] / (sigma * sigma);
    if (decode_trace.llr) decode_trace.llr(llr);

    const int max_iterations = 10;
    signed char decoded77[77] = {};
    signed char cw[128] = {};
    int nharderror = -1, niterations = 0;
    bpdecode128_90(llr, max_iterations, decoded77, cw, nharderror, niterations);  // apmask=0
    if (decode_trace.bp) decode_trace.bp(decoded77, cw, nharderror, niterations);
    if (nharderror >= 0 && nharderror < 18) {
        nsuccess = 1;
        char c77[78];  // write(c77,'(77i1)') decoded77
        for (int i = 0; i < 77; ++i) c77[i] = static_cast<char>('0' + decoded77[i]);
        c77[77] = '\0';
        // read(c77(72:77),'(2b3)') n3,i3
        const int n3 = (c77[71] - '0') * 4 + (c77[72] - '0') * 2 + (c77[73] - '0');
        const int i3 = (c77[74] - '0') * 4 + (c77[75] - '0') * 2 + (c77[76] - '0');
        if (decode_trace.tag) decode_trace.tag(c77, n3, i3);
        if ((i3 == 0 && (n3 == 1 || n3 == 3 || n3 == 4 || n3 > 5)) || i3 == 3 || i3 > 5) {
            nsuccess = 0;
        } else {
            const Packjt77::Unpacked u = pk.unpack77(c77, 1);  // call unpack77(c77,1,msgreceived,unpk77_success)
            const size_t n = u.msg.size() < 37 ? u.msg.size() : 37;
            std::memcpy(msgreceived, u.msg.data(), n);
            if (n < 37) std::memset(msgreceived + n, ' ', 37 - n);
            if (!u.ok) nsuccess = 0;
        }
    }
}

void msk144decodeframe(cfloat* c, float* softbits, char* msgreceived, int& nsuccess) {
    msk144decodeframe(msk144_packjt77(), c, softbits, msgreceived, nsuccess);
}

// ---- genmsk_128_90.f90:1-129 for a decoded message (ichk 0 gives the tones; 1 only msgsent): the '@' and genmsk40
// branches (lines 51-55, 71-80) are not ported; the tones are lines 112-120.
void genmsk_128_90(Packjt77& pk, const char* msg0, int ichk, char* msgsent, int* i4tone, int& itype) {
    static const char kBad[38] = "*** bad message ***                  ";
    const int s8d[8] = {0, 1, 1, 1, 0, 0, 1, 0};
    char message[37];
    auto bad = [&]() {  // 998: msgbits=0, bitseq=0, i4tone=0, msgsent='*** bad message ***', itype=-1
        for (int i = 0; i < 144; ++i) i4tone[i] = 0;
        std::memcpy(msgsent, kBad, 37);
        itype = -1;
    };
    itype = 1;
    if (msg0[0] == '@') {  // the fixed tone: not reachable from a decoded message
        bad();
        return;
    }
    std::memcpy(message, msg0, 37);
    for (int i = 0; i < 37; ++i) {  // a NUL ends the message: blanks from there
        if (message[i] == '\0') {
            std::memset(message + i, ' ', 37 - i);
            break;
        }
    }
    for (int i = 1; i <= 37; ++i) {  // Strip leading blanks: message=message(i+1:), the Fortran's i, not 1
        if (message[0] != ' ') break;
        std::memmove(message, message + i, 37 - i);
        std::memset(message + 37 - i, ' ', i);
    }
    if (message[0] == '<') {
        int i2 = 0;  // index(message,'>')
        for (int i = 0; i < 37; ++i)
            if (message[i] == '>') {
                i2 = i + 1;
                break;
            }
        int i1 = 0;
        if (i2 > 0) {  // index(message(1:i2),' ')
            for (int i = 0; i < i2; ++i)
                if (message[i] == ' ') {
                    i1 = i + 1;
                    break;
                }
        }
        if (i1 > 0) {  // genmsk40: an MSK40 short message, not reachable from a decoded message
            bad();
            return;
        }
    }
    const Packjt77::Packed packed = pk.pack77_legacy_truncating_fallback(std::string(message, 37));
    if (packed.i3 < 0 || packed.n3 < 0) {
        bad();
        return;
    }
    const Packjt77::Unpacked sent = pk.unpack77(packed.c77.data(), 0);  // Unpack to get msgsent
    {
        const size_t n = sent.msg.size() < 37 ? sent.msg.size() : 37;
        std::memcpy(msgsent, sent.msg.data(), n);
        if (n < 37) std::memset(msgsent + n, ' ', 37 - n);
    }
    if (ichk == 1) return;
    signed char msgbits[77];  // read(c77,"(77i1)",err=998) msgbits: blanks read as 0
    for (int i = 0; i < 77; ++i) msgbits[i] = packed.c77[i] == '1' ? 1 : 0;
    signed char codeword[128];
    encode_128_90(msgbits, codeword);

    // Create 144-bit channel vector: 8-bit sync word + 48 bits + 8-bit sync word + 80 bits
    int bitseq[144];
    for (int i = 0; i < 8; ++i) bitseq[i] = s8d[i];
    for (int i = 0; i < 48; ++i) bitseq[8 + i] = codeword[i];
    for (int i = 0; i < 8; ++i) bitseq[56 + i] = s8d[i];
    for (int i = 0; i < 80; ++i) bitseq[64 + i] = codeword[48 + i];
    for (int i = 0; i < 144; ++i) bitseq[i] = 2 * bitseq[i] - 1;
    // Map I and Q to tones (the waveforms xi, xq the Fortran also makes are not returned)
    for (int i = 1; i <= 144; ++i) i4tone[i - 1] = 0;
    for (int i = 1; i <= 72; ++i) {
        // integer division truncates toward zero, as Fortran's
        i4tone[2 * i - 1 - 1] = (bitseq[2 * i - 1] * bitseq[2 * i - 2] + 1) / 2;
        i4tone[2 * i - 1] = -(bitseq[2 * i - 1] * bitseq[(2 * i) % 144] - 1) / 2;
    }
    // Flip polarity
    for (int i = 0; i < 144; ++i) i4tone[i] = -i4tone[i] + 1;
}

// ---- msk144signalquality.f90:1-213: the non-training part, lines 44-55 and 83-129 (the btrain-gated
// lines 56-82 and the training, 131-212, are compiled out).
// The state the training keeps (navg, the averaged cross spectrum, currently_training) is never touched with
// btrain false.
void msk144signalquality(Packjt77& pk, const cfloat* /*cframe*/, float /*snr*/, float /*freq*/, float /*t0*/,
                         const float* softbits, const char* msg, const char* /*dxcall*/, bool /*btrain*/,
                         const char* /*datadir*/, int& nbiterrors, float& eyeopening, const double* /*pcoeffs*/) {
    int hardbits[144], msgbits[144], i4tone[144];
    char msgsent[37];
    int itype;

    // use decoded message to figure out how many bit errors in the frame
    for (int i = 0; i < 144; ++i) {
        hardbits[i] = 0;
        if (softbits[i] > 0) hardbits[i] = 1;
    }

    // generate tones from decoded message
    const int ichk = 0;
    genmsk_128_90(pk, msg, ichk, msgsent, i4tone, itype);

    // reconstruct message bits from tones
    msgbits[0] = 0;
    for (int i = 1; i <= 143; ++i) {
        if (i4tone[i - 1] == 0) {
            if (i % 2 == 1) {
                msgbits[i] = msgbits[i - 1];
            } else {
                msgbits[i] = (msgbits[i - 1] + 1) % 2;
            }
        } else {
            if (i % 2 == 1) {
                msgbits[i] = (msgbits[i - 1] + 1) % 2;
            } else {
                msgbits[i] = msgbits[i - 1];
            }
        }
    }

    nbiterrors = 0;
    for (int i = 0; i < 144; ++i) {
        if (hardbits[i] != msgbits[i]) nbiterrors = nbiterrors + 1;
    }

    float eyetop = 1;
    float eyebot = -1;
    for (int i = 0; i < 144; ++i) {
        if (msgbits[i] == 1) {
            if (softbits[i] < eyetop) eyetop = softbits[i];
        } else {
            if (softbits[i] > eyebot) eyebot = softbits[i];
        }
    }
    eyeopening = eyetop - eyebot;
}

void msk144signalquality(const cfloat* cframe, float snr, float freq, float t0, const float* softbits, const char* msg,
                         const char* dxcall, bool btrain, const char* datadir, int& nbiterrors, float& eyeopening,
                         const double* pcoeffs) {
    msk144signalquality(msk144_packjt77(), cframe, snr, freq, t0, softbits, msg, dxcall, btrain, datadir, nbiterrors,
                        eyeopening, pcoeffs);
}

// ---- fmtmsg.f90:1-21
void fmtmsg(std::string& msg, int& iz) {
    const int len = static_cast<int>(msg.size());
    // Convert all letters to upper case
    iz = len;
    for (int i = 1; i <= len; ++i) {
        char& ch = msg[i - 1];
        if (ch >= 'a' && ch <= 'z') ch = static_cast<char>(ch + 'A' - 'a');
        if (ch != ' ') iz = i;
    }
    for (int iter = 1; iter <= 37; ++iter) {  // Collapse multiple blanks into one
        // ib2=index(msg(1:iz),'  ')
        int ib2 = 0;
        for (int i = 1; i + 1 <= iz; ++i)
            if (msg[i - 1] == ' ' && msg[i] == ' ') {
                ib2 = i;
                break;
            }
        if (ib2 < 1) return;
        // msg=msg(1:ib2)//msg(ib2+2:): one blank fewer, blank padded to len
        msg.erase(static_cast<size_t>(ib2), 1);
        msg.push_back(' ');
        iz = iz - 1;
    }
}

}  // namespace msk144
