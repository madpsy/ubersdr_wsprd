// MSK144's frame decoder and signal quality: WSJT-X's msk144decodeframe.f90,
// msk144signalquality.f90 (without its training paths), genmsk_128_90.f90
// (the 77-bit message to the 144 tones, as signal quality uses it) and
// fmtmsg.f90 at 567ad29ce (3.2.0-rc1), line for line, giving the reference's
// numbers to the bit on x86-64. It rests on what src/msk144_ldpc.h and
// src/msk144_sync.h rest on (no multiply-adds fused, complex products written
// out as gfortran does, glibc's maths where gfortran calls it: here sincosf,
// atan2f (ft8::atan2f_glibc, added to src/ft8_maths.h for it) and libmvec's
// sinf four at a time for the half-sine pulse), and on what the reference's
// compiler did with these routines (the comments in src/msk144_decode.cpp say
// where).
//
// The module packjt77's tables are shared state in the reference: unpack77
// (and pack77, in genmsk_128_90) save the calls they see, and a later
// "<...>" resolves from them. Here that state is a Packjt77 (src/packjt77.h):
// the caller's, or, without one, msk144_packjt77()'s, the module's single
// copy. A file's decoding starts it afresh (msk144_reset_packjt77).
//
// Joe Taylor, K1JT, Steven Franke, K9AN, and the WSJT-X team (GNU GPL v3),
// and so is this.
#pragma once

#include <string>

#include "msk_frontend.h"
#include "packjt77.h"

namespace msk144 {

using packjt77::Packjt77;

// msk144decodeframe(c, softbits, msgreceived, nsuccess): c the NSPM = 864
// points of a frame, which the decoder corrects in place (the phase
// correction, msk144decodeframe.f90:59, as the reference's `c=c*conjg(cfac)`
// does: msk144signalquality and msk144spd's caller read the corrected frame),
// softbits 144 values out (the matched filter's, normalised when the frame got
// past the sync check), msgreceived 37 characters (blank padded, as
// Fortran's character*37; written whatever happens), nsuccess out (> 0:
// decoded). Matches msk144::DecodeSink (src/msk144_sync.h), which msk144spd
// calls it through.
void msk144decodeframe(Packjt77& pk, cfloat* c, float* softbits, char* msgreceived, int& nsuccess);
// The same with the module's one Packjt77, msk144_packjt77().
void msk144decodeframe(cfloat* c, float* softbits, char* msgreceived, int& nsuccess);

// The module packjt77's state in the reference (one copy for the process)
// and its restart.
Packjt77& msk144_packjt77();
void msk144_reset_packjt77();

// msk144signalquality(cframe, snr, freq, t0, softbits, msg, dxcall, btrain,
// datadir, nbiterrors, eyeopening, pcoeffs): the Fortran's argument list.
// The training paths (btrain true, and what they write to the data directory)
// are compiled out, as file mode has btrain false for good: cframe, snr,
// freq, t0, dxcall, btrain, datadir and pcoeffs are read by nothing here, and
// btrain must be false. msg is 37 characters, softbits 144 values; the result
// is nbiterrors (the hard decisions of softbits against the bits the message
// makes) and eyeopening. pk is the genmsk_128_90's Packjt77 (the module's:
// pack77 and unpack77 save calls).
void msk144signalquality(Packjt77& pk, const cfloat* cframe, float snr, float freq, float t0, const float* softbits,
                         const char* msg, const char* dxcall, bool btrain, const char* datadir, int& nbiterrors,
                         float& eyeopening, const double* pcoeffs);
void msk144signalquality(const cfloat* cframe, float snr, float freq, float t0, const float* softbits, const char* msg,
                         const char* dxcall, bool btrain, const char* datadir, int& nbiterrors, float& eyeopening,
                         const double* pcoeffs);

// genmsk_128_90(msg0, ichk, msgsent, i4tone, itype) for a decoded message:
// msg0 and msgsent are 37 characters, i4tone 144 values (0/1). The two
// branches the reference sends elsewhere are not ported, since a message the
// decoder returns cannot reach them: a message starting '@' (the fixed tone;
// the decoder's alphabet has no '@') and one starting '<' with a blank
// before its '>' (an MSK40 short message, genmsk40; a decoded call has no
// blank inside its brackets). Both give itype -1 and the all-zero tones here.
void genmsk_128_90(Packjt77& pk, const char* msg0, int ichk, char* msgsent, int* i4tone, int& itype);

// fmtmsg(msg, iz): msg.size() is Fortran's len(msg): the upper case of msg,
// iz its last non-blank character (len when all are blank) less the blanks
// collapsed (up to 37 times, each taking one from iz).
void fmtmsg(std::string& msg, int& iz);

// Debug-only traces for the piece check: null in production (an untaken
// branch at each place). Each is called where scripts/msk144-piece-check.sh
// frame's instrumented copy of the reference calls its hook, with what that
// hook prints.
struct DecodeTrace {
    void (*consts)(const float* pp /*12*/, const cfloat* cb /*42*/) = nullptr;  // first call
    void (*phase)(cfloat cca, cfloat ccb, float phase0) = nullptr;
    void (*corrected)(cfloat cfac, const cfloat* c /*864*/) = nullptr;  // after c=c*conjg(cfac)
    void (*raw)(const float* softbits /*144*/) = nullptr;                // after the matched filter
    void (*sync)(int nbadsync1, int nbadsync2) = nullptr;
    void (*norm)(float sav, float s2av, float ssig) = nullptr;
    void (*llr)(const float* llr /*128*/) = nullptr;
    void (*bp)(const signed char* decoded77 /*77*/, const signed char* cw /*128*/, int nharderror, int niterations) = nullptr;
    void (*tag)(const char* c77 /*77*/, int n3, int i3) = nullptr;
};
extern DecodeTrace decode_trace;

}  // namespace msk144
