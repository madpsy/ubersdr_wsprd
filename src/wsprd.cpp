// WSJT-X's wsprd.c (lib/wsprd at 567ad29ce), the decoder: see wsprd.h.
//
// wsprd.c: Copyright 2001-2018, Joe Taylor, K1JT; much of it based on work
// by Steven Franke, K9AN, Copyright 2014-2018; GNU GPL v3 (or later), and
// so is this.
//
// The reference's code, kept as it was, names, layout and comments, so the
// two read side by side; Decoder::Impl::decode() is main() from the input
// read to the sorted decodes. What had to change, and nothing else:
//   - main()'s locals the run keeps (hash table, wisdom, timers) are the
//     Decoder's, its buffers Impl's; the file handling, getopt and printing
//     are in wsprd_files.cpp (run());
//   - PLAN1/2/3 are per decode, made in one go before any is run, under the
//     process's FFTW planner lock (fft.h) with its wisdom set aside and
//     wsprd's own (wspr_wisdom.dat) in its place: the same plans in the same
//     order from the same wisdom, so the same plans;
//   - the static fplast/dt/df/pi of sync_and_demodulate and
//     noncoherent_sequence_detection are locals (fplast never carries
//     over: i==0 always recomputes);
//   - every maths call on a float names its double, as C's <math.h> calls
//     it: in C++, cos(float) would be cosf;
//   - ps[512][nffts] (a VLA) is a vector, ps[j][i] at j*nffts+i;
//   - subtract_signal() (unused by wsprd) is left out;
//   - osdwspr is wsprd_osd.h's; it takes ndepth by value (the Fortran's 5
//     for more than 5 written back changes nothing, 5 and up decode alike);
//   - the coarse search takes ps's square roots once per pass, in place,
//     and each (ifr, idrift)'s ifd once rather than for every time offset;
//   - with AVX2 (x86-64) or NEON (aarch64), sync_and_demodulate's
//     searches and noncoherent_sequence_detection's sums run eight lanes at
//     a time (sync_search, noncoherent_sums), as do the coarse search's time
//     offsets (coarse_sync), every float as the reference's, and a candidate's sums
//     already made for the same (f0, drift, lag) are taken again rather
//     than redone (SyncMemo, NcMemo); with a drift the turn tables are made
//     row by row as the sums take them, never stored (live_sums), a
//     table's step cos/sin kept by fp for the decode (TurnSteps), with a
//     drift a try's next four shifts summed at once, one table per symbol
//     for all four (noncoherent_lags), with no drift a lag 256 below
//     another's takes that one's sums of the symbol before (one table for
//     every symbol, the same samples), and subtract_signal2's low-pass
//     filter runs 32 outputs at a time, only those the subtraction reads
//     (with NEON, aarch64, 16 at a time; in the plain code each output's
//     sum is kept in a local, not in cfi/cfq, which gcc must store each
//     term to);
//   - subtract_signal2's arrays are kept from call to call (SubScratch),
//     not allocated again (only ci and cq zeroed);
//   - WSPRD_TRACE=1 prints each stage to stderr, as the oracle's
//     instrumented reference does (tests/data/wsprd/README.md), to diff.
#include "wsprd.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <initializer_list>
#include <memory>
#include <mutex>
#include <new>
#include <vector>
#include <fftw3.h>

#if defined(__AVX2__) && (defined(__x86_64__) || defined(_M_X64))
#include <immintrin.h>
#define WSPRD_AVX2 1
#elif defined(__aarch64__) && defined(__ARM_NEON)
// Products and sums only (vmulq, vaddq, vsubq), never vfmaq/vmlaq: each
// rounded as the plain code's are (no contraction, CMakeLists.txt)
#include <arm_neon.h>
#define WSPRD_NEON 1
#endif
#if defined(WSPRD_AVX2) || defined(WSPRD_NEON)
#define WSPRD_SIMD 1  // the sums of sync_and_demodulate and noncoherent_sequence_detection in V8s
#endif

#include "fft.h"  // fft_planner_lock(): every FFTW plan in the process is made under it
#include "wsprd_internal.h"
#include "wsprd_osd.h"

namespace wsprd {

clock_t detail::cpu_clock() {
#ifdef _WIN32
    return clock();
#else
    struct timespec ts;
    if (clock_gettime(CLOCK_THREAD_CPUTIME_ID, &ts) != 0) return (clock_t)-1;
    return ts.tv_sec * CLOCKS_PER_SEC + ts.tv_nsec / (1000000000 / CLOCKS_PER_SEC);
#endif
}

using namespace detail;

namespace {

// WSPRD_TRACE=1: stage-by-stage state on stderr. Floats as %a (exact);
// arrays as FNV-1a 64 of their bytes.
bool trace() {
    static const bool on = [] {
        const char* e = std::getenv("WSPRD_TRACE");
        return e && *e && *e != '0';
    }();
    return on;
}
unsigned long long fnv64(const void* p, size_t n) {
    const unsigned char* b = static_cast<const unsigned char*>(p);
    unsigned long long h = 0xcbf29ce484222325ULL;
    for (size_t i = 0; i < n; i++) {
        h ^= b[i];
        h *= 0x100000001b3ULL;
    }
    return h;
}
#define TR(...) do { if (trace()) fprintf(stderr, __VA_ARGS__); } while (0)

// Possible PATIENCE options: FFTW_ESTIMATE, FFTW_ESTIMATE_PATIENT,
// FFTW_MEASURE, FFTW_PATIENT, FFTW_EXHAUSTIVE
#define PATIENCE FFTW_ESTIMATE

const unsigned char pr3[162]=
{1,1,0,0,0,0,0,0,1,0,0,0,1,1,1,0,0,0,1,0,
    0,1,0,1,1,1,1,0,0,0,0,0,0,0,1,0,0,1,0,1,
    0,0,0,0,0,0,1,0,1,1,0,0,1,1,0,1,0,0,0,1,
    1,0,1,0,0,0,0,1,1,0,1,0,1,0,1,0,1,0,0,1,
    0,0,1,0,1,1,0,0,0,1,1,0,1,0,1,0,0,0,1,0,
    0,0,0,0,1,0,0,1,0,0,1,1,1,0,1,1,0,0,1,1,
    0,1,0,0,0,1,1,1,0,0,0,0,0,1,0,1,0,0,1,1,
    0,0,0,0,0,0,0,1,1,0,1,0,1,1,0,0,0,1,1,0,
    0,0};

// PLAN1 (readwavfile's r2c), PLAN2 (its c2c back) and PLAN3 (the 512-point
// spectra), with their buffers. make() plans them as wsprd would in a
// process that had imported wisdom: under the planner lock, the process's
// own wisdom exported, forgotten and put back after, wisdom (what the
// Decoder has, wspr_wisdom.dat's) imported first and replaced by what the
// planner holds once they are made (what wsprd exports at its end).
struct Plans {
    int nfft1 = 0;  // 0: no PLAN1/PLAN2 (a .c2)
    float* realin = nullptr;
    fftwf_complex* fftout1 = nullptr;
    fftwf_complex* fftin2 = nullptr;
    fftwf_complex* fftout2 = nullptr;
    fftwf_complex* fftin3 = nullptr;
    fftwf_complex* fftout3 = nullptr;
    fftwf_plan PLAN1 = nullptr, PLAN2 = nullptr, PLAN3 = nullptr;
    // The ticks PLAN1 and PLAN2 took to make, buffers included: readwavfile
    // makes them, so they are its timer's (PLAN3, the wisdom swaps and the
    // wait for the lock are not).
    clock_t plan12 = 0;

    Plans() = default;
    Plans(const Plans&) = delete;
    Plans& operator=(const Plans&) = delete;

    void make(int nfft1_, int nfft2, std::string& wisdom) {
        std::lock_guard<std::mutex> lock(skimmer::fft_planner_lock());
        char* const saved = fftwf_export_wisdom_to_string();
        fftwf_forget_wisdom();
        if (!wisdom.empty()) fftwf_import_wisdom_from_string(wisdom.c_str());
        nfft1 = nfft1_;
        bool ok = true;
        if (nfft1) {
            const clock_t t0 = cpu_clock();
            realin = (float*)fftwf_malloc(sizeof(float) * nfft1);
            fftout1 = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * (nfft1 / 2 + 1));
            fftin2 = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * nfft2);
            fftout2 = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * nfft2);
            ok = realin && fftout1 && fftin2 && fftout2;
            if (ok) PLAN1 = fftwf_plan_dft_r2c_1d(nfft1, realin, fftout1, PATIENCE);
            if (ok) PLAN2 = fftwf_plan_dft_1d(nfft2, fftin2, fftout2, FFTW_BACKWARD, PATIENCE);
            ok = ok && PLAN1 && PLAN2;
            plan12 = cpu_clock() - t0;
        }
        fftin3 = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * 512);
        fftout3 = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * 512);
        ok = ok && fftin3 && fftout3;
        if (ok) PLAN3 = fftwf_plan_dft_1d(512, fftin3, fftout3, FFTW_FORWARD, PATIENCE);
        ok = ok && PLAN3;
        if (ok) {
            char* const now = fftwf_export_wisdom_to_string();
            if (now) {
                wisdom = now;
                std::free(now);
            }
        }
        fftwf_forget_wisdom();
        if (saved) {
            fftwf_import_wisdom_from_string(saved);
            std::free(saved);
        }
        if (!ok) throw std::bad_alloc();
    }

    ~Plans() {
        std::lock_guard<std::mutex> lock(skimmer::fft_planner_lock());
        if (PLAN1) fftwf_destroy_plan(PLAN1);
        if (PLAN2) fftwf_destroy_plan(PLAN2);
        if (PLAN3) fftwf_destroy_plan(PLAN3);
        fftwf_free(realin);
        fftwf_free(fftout1);
        fftwf_free(fftin2);
        fftwf_free(fftout2);
        fftwf_free(fftin3);
        fftwf_free(fftout3);
    }
};

//***************************************************************************
// readwavfile() from its buffer on: buf2 the file's samples, npoints of them
// (zero past what the file had).
// A sample at full scale 1.0: a 16-bit one as wsprd scales it (/32768.0),
// a 32-bit float one as it is.
inline float wav_sample(short int s) { return s/32768.0; }
inline float wav_sample(float s) { return s; }

template <class Sample>
unsigned long readwavfile(const Sample* buf2, size_t npoints, int ntrmin, float* idat, float* qdat, Plans& plans)
{
    size_t i, j;
    int nfft1, nfft2, nh2, i0;
    double df;

    nfft2=46080; //this is the number of downsampled points that will be returned
    nh2=nfft2/2;

    if( ntrmin == 2 ) {
        nfft1=nfft2*32;      //need to downsample by a factor of 32
        df=12000.0/nfft1;
        i0=1500.0/df+0.5;
    } else if ( ntrmin == 15 ) {
        nfft1=nfft2*8*32;
        df=12000.0/nfft1;
        i0=(1500.0+112.5)/df+0.5;
    } else {
        fprintf(stderr,"This should not happen\n");
        return 1;
    }

    float *realin = plans.realin;
    fftwf_complex *fftin, *fftout = plans.fftout1;

    for (i=0; i<npoints; i++) {
        realin[i]=wav_sample(buf2[i]);
    }

    for (i=npoints; i<(size_t)nfft1; i++) {
        realin[i]=0.0;
    }

    fftwf_execute(plans.PLAN1);

    fftin=plans.fftin2;

    for (i=0; i<(size_t)nfft2; i++) {
        j=i0+i;
        if( i>(size_t)nh2 ) j=j-nfft2;
        fftin[i][0]=fftout[j][0];
        fftin[i][1]=fftout[j][1];
    }

    fftout=plans.fftout2;
    fftwf_execute(plans.PLAN2);

    for (i=0; i<(size_t)nfft2; i++) {
        idat[i]=fftout[i][0]/1000.0;
        qdat[i]=fftout[i][1]/1000.0;
    }

    return nfft2;
}

uint32_t float_bits(float x) { uint32_t u; memcpy(&u, &x, sizeof u); return u; }

// The cos and sin of a turn table's four steps by its fp to the bit, kept
// for a whole decode: a candidate's searches and tries make the tables of
// the same fps again and again (each of the 17 shifts a try sums, the
// frequency search's centre, ...). Direct-mapped, a slot overwritten on a
// clash. Only the SIMD sums use it.
struct TurnSteps {
    enum { SLOTS=4096 };
    struct Slot { uint32_t fp; bool ok; float c[4], s[4]; };
    std::vector<Slot> slot;
    Slot* at(float fp) {
        if (slot.empty()) slot.resize(SLOTS);
        return &slot[(float_bits(fp)*2654435761u) >> 20];
    }
};

// The syncs one candidate's refinement has already summed, by (f0, drift,
// lag) to the bit: each is a function of those and the data alone, and the
// next search often starts from the last one's best (a lag search's best is
// the frequency search's centre, and so on). Only the SIMD searches use it.
struct SyncMemo {
    struct Case { uint32_t f0, drift; int lag; float sync; };
    Case c[32];
    int n=0;
    TurnSteps* steps=nullptr;
    static uint32_t bits(float x) { return float_bits(x); }
    bool find(float f0, float drift, int lag, float* sync) const {
        for (int i=0; i<n; i++)
            if (c[i].lag == lag && c[i].f0 == bits(f0) && c[i].drift == bits(drift)) {
                *sync=c[i].sync;
                return true;
            }
        return false;
    }
    void add(float f0, float drift, int lag, float sync) {
        if (n < 32) c[n++]={bits(f0), bits(drift), lag, sync};
    }
};

// noncoherent_sequence_detection's is/qs/cf/sf of one candidate's tries,
// by (f0, drift, lag) to the bit, as SyncMemo: the block sizes after the
// first try the same jittered shifts again. Only the SIMD sums use it.
struct NcMemo {
    enum { MAX=17 };  // the shifts a block size tries
    struct Sums { uint32_t f0, drift; int lag; float is[4][162], qs[4][162], cf[4][162], sf[4][162]; };
    std::vector<Sums> s;
    int n=0;
    TurnSteps* steps=nullptr;
    const Sums* find(float f0, float drift, int lag) const {
        for (int i=0; i<n; i++)
            if (s[i].lag == lag && s[i].f0 == float_bits(f0) && s[i].drift == float_bits(drift)) return &s[i];
        return nullptr;
    }
    Sums* add(float f0, float drift, int lag) {
        if (n == MAX) return nullptr;
        if (s.empty()) s.resize(MAX);
        Sums* e=&s[n++];
        e->f0=float_bits(f0); e->drift=float_bits(drift); e->lag=lag;
        return e;
    }
};

#ifdef WSPRD_SIMD
// sync_and_demodulate's and noncoherent_sequence_detection's sums with AVX2
// (x86-64) or NEON (aarch64); elsewhere the plain code. Every float is the
// one the reference computes, by the same operations in the same order: the
// eight vector lanes are the eight sums (or turns) the reference keeps side
// by side, never one sum split up.
//
// V8 is those eight lanes: an __m256, or two NEON vectors (lanes 0..3 and
// 4..7). Its operations are each lane's own IEEE product or sum, never a
// fused one.
#ifdef WSPRD_AVX2
typedef __m256 V8;
inline V8 v8_zero() { return _mm256_setzero_ps(); }
inline V8 v8_dup(const float* p) { return _mm256_broadcast_ss(p); }
inline V8 v8_load(const float* p) { return _mm256_loadu_ps(p); }
inline void v8_store(float* p, V8 x) { _mm256_storeu_ps(p, x); }
inline V8 v8_add(V8 a, V8 b) { return _mm256_add_ps(a, b); }
inline V8 v8_mul(V8 a, V8 b) { return _mm256_mul_ps(a, b); }
// lanes 4..7 then 0..3
inline V8 v8_swap(V8 x) { return _mm256_permute2f128_ps(x, x, 1); }
// lanes 4..7 negated
inline V8 v8_neg_hi(V8 x) {
    const __m256 neg=_mm256_setr_ps(0.0f, 0.0f, 0.0f, 0.0f, -0.0f, -0.0f, -0.0f, -0.0f);
    return _mm256_xor_ps(x, neg);
}
#else
struct V8 { float32x4_t lo, hi; };
inline V8 v8_zero() { return {vdupq_n_f32(0.0f), vdupq_n_f32(0.0f)}; }
inline V8 v8_dup(const float* p) { const float32x4_t x=vld1q_dup_f32(p); return {x, x}; }
inline V8 v8_load(const float* p) { return {vld1q_f32(p), vld1q_f32(p+4)}; }
inline void v8_store(float* p, V8 x) { vst1q_f32(p, x.lo); vst1q_f32(p+4, x.hi); }
inline V8 v8_add(V8 a, V8 b) { return {vaddq_f32(a.lo, b.lo), vaddq_f32(a.hi, b.hi)}; }
inline V8 v8_mul(V8 a, V8 b) { return {vmulq_f32(a.lo, b.lo), vmulq_f32(a.hi, b.hi)}; }
inline V8 v8_swap(V8 x) { return {x.hi, x.lo}; }
#endif

// One fp's turn tables: row j is the reference's c0[j]..c3[j], s0[j]..s3[j]
// as X=[c,s], kept as the two vectors a sum takes, cs and sc=[s,c]
// (swap(X)). With AVX2 cs=[c,-s], and a row from the one before is
// X*[cd,cd] + swap(X)*[-sd,sd]: c*cd + s*(-sd) is c*cd - s*sd, and
// s*cd + c*sd is c*sd + s*cd, exactly. With NEON, whose V8 is the halves c
// and s anyway, cs is X itself and sc its halves swapped, so a row is
// stored once (8 floats, not 16), a row from the one before is the
// reference's c*cd - s*sd and c*sd + s*cd, and the sums subtract id*s as
// it does: no negation, the same floats.
struct Turns {
#ifdef WSPRD_AVX2
    struct Row { float cs[8], sc[8]; };
#else
    struct Row { float x[8]; };
#endif
    alignas(64) Row x[257];
    float fp;
};

// Table t's row j: its cs and sc, its c[i] and s[i]
#ifdef WSPRD_AVX2
inline V8 row_cs(const Turns* t, int j) { return v8_load(t->x[j].cs); }
inline V8 row_sc(const Turns* t, int j) { return v8_load(t->x[j].sc); }
inline float row_c(const Turns* t, int j, int i) { return t->x[j].cs[i]; }
inline float row_s(const Turns* t, int j, int i) { return t->x[j].sc[i]; }
#else
inline V8 row_cs(const Turns* t, int j) { return v8_load(t->x[j].x); }
inline V8 row_sc(const Turns* t, int j) { return v8_swap(v8_load(t->x[j].x)); }
inline float row_c(const Turns* t, int j, int i) { return t->x[j].x[i]; }
inline float row_s(const Turns* t, int j, int i) { return t->x[j].x[4+i]; }
#endif

// Lanes i0..i3, q0..q3 plus one j's terms: each i + id*c + qd*s and
// q - id*s + qd*c in the reference's order (with AVX2 q + id*(-s), which is
// q - id*s exactly)
inline V8 turn_add(V8 a, V8 vi, V8 vq, V8 cs, V8 sc) {
#ifdef WSPRD_AVX2
    return v8_add(v8_add(a, v8_mul(vi, cs)), v8_mul(vq, sc));
#else
    return {vaddq_f32(vaddq_f32(a.lo, vmulq_f32(vi.lo, cs.lo)), vmulq_f32(vq.lo, sc.lo)),
            vaddq_f32(vsubq_f32(a.hi, vmulq_f32(vi.hi, cs.hi)), vmulq_f32(vq.hi, sc.hi))};
#endif
}

void turn_steps(float fp, float cd[4], float sd[4], TurnSteps* steps) {
    TurnSteps::Slot* e=steps ? steps->at(fp) : nullptr;
    if (e && e->ok && e->fp == float_bits(fp)) {
        memcpy(cd, e->c, sizeof e->c); memcpy(sd, e->s, sizeof e->s);
        return;
    }
    const float dt=1.0/375.0, df=375.0/256.0;
    const float pi=3.14159265358979323846;
    float twopidt=2*pi*dt, df15=df*1.5, df05=df*0.5;
    float dphi0=twopidt*(fp-df15);
    cd[0]=cos((double)dphi0);
    sd[0]=sin((double)dphi0);
    float dphi1=twopidt*(fp-df05);
    cd[1]=cos((double)dphi1);
    sd[1]=sin((double)dphi1);
    float dphi2=twopidt*(fp+df05);
    cd[2]=cos((double)dphi2);
    sd[2]=sin((double)dphi2);
    float dphi3=twopidt*(fp+df15);
    cd[3]=cos((double)dphi3);
    sd[3]=sin((double)dphi3);
    if (e) {
        e->fp=float_bits(fp); e->ok=true;
        memcpy(e->c, cd, sizeof e->c); memcpy(e->s, sd, sizeof e->s);
    }
}

// A table's row 0, X=[1,1,1,1,0,0,0,0], and its step [cd,cd], [-sd,sd]
// (with NEON [sd,sd]).
void turn_start(float fp, V8* x, V8* cd, V8* sd, TurnSteps* steps) {
    float c[4], s[4];
    turn_steps(fp, c, s, steps);
#ifdef WSPRD_AVX2
    const __m128 c4=_mm_loadu_ps(c), s4=_mm_loadu_ps(s);
    *cd=_mm256_set_m128(c4, c4);
    *sd=_mm256_set_m128(s4, _mm_sub_ps(_mm_setzero_ps(), s4));
    *x=_mm256_setr_ps(1, 1, 1, 1, 0, 0, 0, 0);
#else
    const float32x4_t c4=vld1q_f32(c), s4=vld1q_f32(s);
    *cd={c4, c4};
    *sd={s4, s4};
    *x={vdupq_n_f32(1.0f), vdupq_n_f32(0.0f)};
#endif
}

// Row j's cs and sc from its X, which then steps to row j+1's.
inline void turn_row(V8* x, V8 cd, V8 sd, V8* cs, V8* sc) {
    *sc=v8_swap(*x);
#ifdef WSPRD_AVX2
    *cs=v8_neg_hi(*x);
    *x=v8_add(v8_mul(*x, cd), v8_mul(*sc, sd));
#else
    *cs=*x;
    *x={vsubq_f32(vmulq_f32(cs->lo, cd.lo), vmulq_f32(cs->hi, sd.lo)),
        vaddq_f32(vmulq_f32(cs->lo, sd.hi), vmulq_f32(cs->hi, cd.hi))};
#endif
}

// n rows (256, or 257 for noncoherent_sequence_detection's cf/sf) of the N
// tables t[], side by side so no step waits on the one before.
template <int N>
void build_turns(Turns* const* t, int n, TurnSteps* steps) {
    V8 x[N], cd[N], sd[N];
    for (int m=0; m<N; m++) turn_start(t[m]->fp, &x[m], &cd[m], &sd[m], steps);
    for (int j=0; j<n; j++) {
        for (int m=0; m<N; m++) {
            V8 cs, sc;
            turn_row(&x[m], cd[m], sd[m], &cs, &sc);
#ifdef WSPRD_AVX2
            v8_store(t[m]->x[j].cs, cs);
            v8_store(t[m]->x[j].sc, sc);
#else
            v8_store(t[m]->x[j].x, cs);  // sc is its halves swapped
#endif
        }
    }
}

void build_turns(Turns* const* t, int n, int rows, TurnSteps* steps) {
    for (; n >= 4; n -= 4, t += 4) build_turns<4>(t, rows, steps);
    switch (n) {
        case 3: build_turns<3>(t, rows, steps); break;
        case 2: build_turns<2>(t, rows, steps); break;
        case 1: build_turns<1>(t, rows, steps); break;
    }
}

// N symbol sums for j in [jlo,jhi): sum m turns samples k0[m]+j by table
// x[m], by turn_add. share: the N
// sums have one table (SAME_TABLE, a lag search) or one k0 (SAME_LAG, a
// frequency search), its loads made once for all N.
enum { SHARE_NONE, SAME_TABLE, SAME_LAG };
template <int N, int share=SHARE_NONE>
void turn_sums(const float* id, const float* qd, const long* k0, const Turns* const* x,
               int jlo, int jhi, V8* acc) {
    V8 a[N];
    for (int m=0; m<N; m++) a[m]=v8_zero();
    for (int j=jlo; j<jhi; j++) {
        V8 cs0, sc0, vi0, vq0;
        if (share == SAME_TABLE) {
            cs0=row_cs(x[0], j); sc0=row_sc(x[0], j);
        }
        if (share == SAME_LAG) {
            vi0=v8_dup(id+k0[0]+j); vq0=v8_dup(qd+k0[0]+j);
        }
        for (int m=0; m<N; m++) {
            const V8 cs=share == SAME_TABLE ? cs0 : row_cs(x[m], j);
            const V8 sc=share == SAME_TABLE ? sc0 : row_sc(x[m], j);
            const V8 vi=share == SAME_LAG ? vi0 : v8_dup(id+k0[m]+j);
            const V8 vq=share == SAME_LAG ? vq0 : v8_dup(qd+k0[m]+j);
            a[m]=turn_add(a[m], vi, vq, cs, sc);
        }
    }
    for (int m=0; m<N; m++) acc[m]=a[m];
}

// turn_sums<N> over a batch of up to five, sharing what it can.
template <int share>
void turn_batch(const float* id, const float* qd, const long* k0, const Turns* const* x,
                int n, V8* acc) {
    switch (n) {
        case 5: turn_sums<5, share>(id, qd, k0, x, 0, 256, acc); break;
        case 4: turn_sums<4, share>(id, qd, k0, x, 0, 256, acc); break;
        case 3: turn_sums<3, share>(id, qd, k0, x, 0, 256, acc); break;
        case 2: turn_sums<2, share>(id, qd, k0, x, 0, 256, acc); break;
        default: turn_sums<1>(id, qd, k0, x, 0, 256, acc); break;
    }
}

// turn_sums<N, share> over all 256 j with live tables: each row made as
// the sums take it, never stored, the rows build_turns stores to the bit.
// For a drift, where fp changes every symbol and a stored table serves one
// symbol only, building it is a chain of steps each waiting on the last;
// here the sums' work overlaps it. fp[m]: sum m's table (fp[0] for all with
// SAME_TABLE); x256, if given, gets each table's row 256 X (cf, sf).
template <int N, int share=SHARE_NONE>
void live_sums(const float* id, const float* qd, const long* k0, const float* fp,
               TurnSteps* steps, V8* acc, V8* x256=nullptr) {
    enum { NT = share == SAME_TABLE ? 1 : N };
    V8 x[NT], cd[NT], sd[NT], a[N];
    for (int t=0; t<NT; t++) turn_start(fp[t], &x[t], &cd[t], &sd[t], steps);
    for (int m=0; m<N; m++) a[m]=v8_zero();
    for (int j=0; j<256; j++) {
        V8 cs[NT], sc[NT], vi0, vq0;
        for (int t=0; t<NT; t++) turn_row(&x[t], cd[t], sd[t], &cs[t], &sc[t]);
        if (share == SAME_LAG) {
            vi0=v8_dup(id+k0[0]+j); vq0=v8_dup(qd+k0[0]+j);
        }
        for (int m=0; m<N; m++) {
            const int t = share == SAME_TABLE ? 0 : m;
            const V8 vi=share == SAME_LAG ? vi0 : v8_dup(id+k0[m]+j);
            const V8 vq=share == SAME_LAG ? vq0 : v8_dup(qd+k0[m]+j);
            a[m]=turn_add(a[m], vi, vq, cs[t], sc[t]);
        }
    }
    for (int m=0; m<N; m++) acc[m]=a[m];
    if (x256) for (int t=0; t<NT; t++) x256[t]=x[t];
}

template <int share>
void live_batch(const float* id, const float* qd, const long* k0, const float* fp,
                TurnSteps* steps, int n, V8* acc) {
    switch (n) {
        case 5: live_sums<5, share>(id, qd, k0, fp, steps, acc); break;
        case 4: live_sums<4, share>(id, qd, k0, fp, steps, acc); break;
        case 3: live_sums<3, share>(id, qd, k0, fp, steps, acc); break;
        case 2: live_sums<2, share>(id, qd, k0, fp, steps, acc); break;
        default: live_sums<1>(id, qd, k0, fp, steps, acc); break;
    }
}

// The j the reference sums for symbol samples from k0: 0 < k0+j < np.
void valid_j(long k0, long np, int* jlo, int* jhi) {
    long lo=1-k0, hi=np-k0;
    *jlo=(int)(lo < 0 ? 0 : lo > 256 ? 256 : lo);
    *jhi=(int)(hi < 0 ? 0 : hi > 256 ? 256 : hi);
    if (*jhi < *jlo) *jhi=*jlo;
}

// The reference's (float)sqrt((double)(i*i + q*q)) of the four tones.
void tone_powers(V8 acc, float p[4]) {
#ifdef WSPRD_AVX2
    const __m128 i=_mm256_castps256_ps128(acc), q=_mm256_extractf128_ps(acc, 1);
    const __m128 pw=_mm_add_ps(_mm_mul_ps(i, i), _mm_mul_ps(q, q));
    _mm_storeu_ps(p, _mm256_cvtpd_ps(_mm256_sqrt_pd(_mm256_cvtps_pd(pw))));
#else
    const float32x4_t pw=vaddq_f32(vmulq_f32(acc.lo, acc.lo), vmulq_f32(acc.hi, acc.hi));
    const float64x2_t lo=vsqrtq_f64(vcvt_f64_f32(vget_low_f32(pw)));
    const float64x2_t hi=vsqrtq_f64(vcvt_high_f64_f32(pw));
    vst1q_f32(p, vcvt_high_f32_f64(vcvt_f32_f64(lo), hi));
#endif
}

// sync_and_demodulate modes 0 and 1: every (ifreq, lag) of the search summed
// at once, symbol by symbol, each with the table of its own fp (rebuilt, as
// the reference does, only when fp changes); a case in memo is taken from
// there, the others added. false: a search too wide for the fixed buffers
// here (none wsprd makes), left to the plain code.
bool sync_search(float *id, float *qd, long np, float *f1, int ifmin, int ifmax, float fstep,
                 int *shift1, int lagmin, int lagmax, int lagstep, float *drift1, float *sync, int mode,
                 SyncMemo* memo)
{
    enum { MAXF=5, MAXC=25, BATCH=5 };
    if( mode == 0 ) {ifmin=0; ifmax=0; fstep=0.0;}
    if( mode == 1 ) {lagmin=*shift1;lagmax=*shift1;}
    if (lagstep <= 0 && lagmin <= lagmax) return false;
    const int nf = ifmax >= ifmin ? ifmax-ifmin+1 : 0;
    const int nlag = lagmax >= lagmin ? (lagmax-lagmin)/lagstep+1 : 0;
    if (nf > MAXF || nf*nlag > MAXC) return false;
    const int nc = nf*nlag;
    TurnSteps* steps = memo ? memo->steps : nullptr;

    Turns tabs[MAXF];
    float f0s[MAXF];
    bool built[MAXF] = {}, used[MAXF] = {};
    int cf[MAXC], clag[MAXC], todo[MAXC], ntodo=0;
    float ss[MAXC], totp[MAXC], known[MAXC];
    bool memoized[MAXC];
    for (int f=0; f<nf; f++) f0s[f]=*f1+(ifmin+f)*fstep;
    for (int c=0; c<nc; c++) {
        cf[c]=c/nlag; clag[c]=lagmin+(c%nlag)*lagstep;
        ss[c]=0.0; totp[c]=0.0;
        memoized[c]=memo && memo->find(f0s[cf[c]], *drift1, clag[c], &known[c]);
        if (!memoized[c]) { todo[ntodo++]=c; used[cf[c]]=true; }
    }

    // symbol i's sums of case c added to its ss and totp, in symbol order
    auto finish = [&](int c, int i, V8 acc) {
        float p[4];
        tone_powers(acc, p);
        float p0=p[0], p1=p[1], p2=p[2], p3=p[3];
        totp[c]=totp[c]+p0+p1+p2+p3;
        float cmet=(p1+p3)-(p0+p2);
        ss[c] = (pr3[i] == 1) ? ss[c]+cmet : ss[c]-cmet;
    };
    auto symbol_fp = [&](int f, int i) -> float {
        float fp = f0s[f] + (*drift1/2.0)*((float)i-81.0)/81.0;
        return fp;
    };
    // a drift: live tables, a lone case (the drift probes) five symbols at a time
    const bool drifting = *drift1 != 0;
    // No drift: every symbol has its frequency's one table, so a case whose
    // lag is 256 below another's (a lag search's first and last) has the
    // same samples and table at symbol i as that one had at i-1: its sums.
    int ahead[MAXC];
    for (int t=0; t<ntodo; t++) {
        const int c=todo[t];
        ahead[c]=-1;
        for (int u=0; u<ntodo && !drifting; u++)
            if (cf[todo[u]] == cf[c] && clag[todo[u]] == clag[c]+256) ahead[c]=todo[u];
    }
    V8 last[MAXC];
    for (int i=0; i<162 && ntodo > 0; i++) {
        long k0[MAXC];
        int full[MAXC], nfull=0;
        for (int t=0; t<ntodo; t++) {
            const int c=todo[t];
            k0[c]=clag[c]+i*256;
            int jlo, jhi;
            valid_j(k0[c], np, &jlo, &jhi);
            if (jlo == 0 && jhi == 256 && !(i > 0 && ahead[c] >= 0)) full[nfull++]=c;
        }
        if (drifting && nfull == ntodo) {
            for (int f=0; f<nf; f++) built[f]=false;
            const int c=todo[0];
            int jlo, jhi;
            valid_j(clag[c]+(i+BATCH-1)*256L, np, &jlo, &jhi);
            if (ntodo == 1 && i+BATCH <= 162 && jlo == 0 && jhi == 256) {
                long bk[BATCH];
                float bf[BATCH];
                V8 ba[BATCH];
                for (int m=0; m<BATCH; m++) { bk[m]=clag[c]+(i+m)*256; bf[m]=symbol_fp(cf[c], i+m); }
                live_sums<BATCH>(id, qd, bk, bf, steps, ba);
                for (int m=0; m<BATCH; m++) finish(c, i+m, ba[m]);
                i+=BATCH-1;
                continue;
            }
            V8 acc[MAXC];
            for (int b=0; b<nfull; b+=BATCH) {
                const int n = nfull-b < BATCH ? nfull-b : BATCH;
                long bk[BATCH];
                float bf[BATCH];
                V8 ba[BATCH];
                bool same_f=true, same_k=true;
                for (int m=0; m<n; m++) {
                    const int cm=full[b+m];
                    bk[m]=k0[cm]; bf[m]=symbol_fp(cf[cm], i);
                    same_f=same_f && cf[cm] == cf[full[b]]; same_k=same_k && bk[m] == bk[0];
                }
                if (same_f) live_batch<SAME_TABLE>(id, qd, bk, bf, steps, n, ba);
                else if (same_k) live_batch<SAME_LAG>(id, qd, bk, bf, steps, n, ba);
                else live_batch<SHARE_NONE>(id, qd, bk, bf, steps, n, ba);
                for (int m=0; m<n; m++) acc[full[b+m]]=ba[m];
            }
            for (int t=0; t<ntodo; t++) finish(todo[t], i, acc[todo[t]]);
            continue;
        }

        Turns* stale[MAXF];
        int nstale=0;
        for (int f=0; f<nf; f++) {
            if (!used[f]) continue;
            float fp = symbol_fp(f, i);
            if( !built[f] || (fp != tabs[f].fp) ) {
                tabs[f].fp=fp;
                built[f]=true;
                stale[nstale++]=&tabs[f];
            }
        }
        build_turns(stale, nstale, 256, steps);

        V8 acc[MAXC];
        const Turns* tx[MAXC];
        for (int t=0; t<ntodo; t++) {
            const int c=todo[t];
            tx[c]=&tabs[cf[c]];
            if (i > 0 && ahead[c] >= 0) {
                acc[c]=last[ahead[c]];
                continue;
            }
            int jlo, jhi;
            valid_j(k0[c], np, &jlo, &jhi);
            if (jlo != 0 || jhi != 256) turn_sums<1>(id, qd, &k0[c], &tx[c], jlo, jhi, &acc[c]);
        }
        for (int b=0; b<nfull; b+=BATCH) {
            const int n = nfull-b < BATCH ? nfull-b : BATCH;
            long bk[BATCH];
            const Turns* bx[BATCH];
            V8 ba[BATCH];
            bool same_x=true, same_k=true;
            for (int m=0; m<n; m++) {
                bk[m]=k0[full[b+m]]; bx[m]=tx[full[b+m]];
                same_x=same_x && bx[m] == bx[0]; same_k=same_k && bk[m] == bk[0];
            }
            if (same_x) turn_batch<SAME_TABLE>(id, qd, bk, bx, n, ba);
            else if (same_k) turn_batch<SAME_LAG>(id, qd, bk, bx, n, ba);
            else turn_batch<SHARE_NONE>(id, qd, bk, bx, n, ba);
            for (int m=0; m<n; m++) acc[full[b+m]]=ba[m];
        }

        for (int t=0; t<ntodo; t++) {
            finish(todo[t], i, acc[todo[t]]);
            last[todo[t]]=acc[todo[t]];
        }
    }

    float syncmax=-1e30, fbest=0.0;
    int best_shift=0;
    for (int c=0; c<nc; c++) {
        float s=memoized[c] ? known[c] : ss[c]/totp[c];
        if (memo && !memoized[c]) memo->add(f0s[cf[c]], *drift1, clag[c], s);
        if( s > syncmax ) {
            syncmax=s;
            best_shift=clag[c];
            fbest=f0s[cf[c]];
        }
    }
    *sync=syncmax;
    *shift1=best_shift;
    *f1=fbest;
    return true;
}

// noncoherent_sequence_detection's is/qs/cf/sf: five symbols at a time,
// each summed with the table of its fp (one table for a run of symbols
// with the same fp, as the reference keeps it).
void noncoherent_sums(float *id, float *qd, long np, float f0, int lag, float drift1,
                      float is[4][162], float qs[4][162], float cf[4][162], float sf[4][162],
                      TurnSteps* steps)
{
    enum { BATCH=5, SLOTS=BATCH+1 };
    Turns tabs[SLOTS];
    int cur=-1;
    for (int i0=0; i0<162; i0+=BATCH) {
        const int n = 162-i0 < BATCH ? 162-i0 : BATCH;
        V8 acc[BATCH];
        long k0[BATCH];
        bool full=true;
        for (int m=0; m<n; m++) {
            k0[m]=lag+(i0+m)*256;
            int jlo, jhi;
            valid_j(k0[m], np, &jlo, &jhi);
            if (jlo != 0 || jhi != 256) full=false;
        }
        if (full && n == BATCH && drift1 != 0) {  // a drift: live tables
            float fps[BATCH];
            V8 x256[BATCH];
            for (int m=0; m<n; m++) {
                const int i=i0+m;
                float fp = f0 + (drift1/2.0)*((float)i-81.0)/81.0;
                fps[m]=fp;
            }
            live_sums<BATCH>(id, qd, k0, fps, steps, acc, x256);
            cur=-1;
            for (int m=0; m<n; m++) {
                const int i=i0+m;
                float a[8], x[8];
                v8_store(a, acc[m]);
                v8_store(x, x256[m]);
                for (int t=0; t<4; t++) {
                    cf[t][i]=x[t]; sf[t][i]=x[4+t];
                    is[t][i]=a[t]; qs[t][i]=a[4+t];
                }
            }
            continue;
        }
        Turns* stale[BATCH];
        int nstale=0;
        const Turns* tx[BATCH];
        for (int m=0; m<n; m++) {
            const int i=i0+m;
            float fp = f0 + (drift1/2.0)*((float)i-81.0)/81.0;
            if( cur < 0 || (fp != tabs[cur].fp) ) {
                cur=(cur+1)%SLOTS;
                tabs[cur].fp=fp;
                stale[nstale++]=&tabs[cur];
            }
            tx[m]=&tabs[cur];
        }
        build_turns(stale, nstale, 257, steps);

        if (full && n == BATCH) {
            turn_sums<BATCH>(id, qd, k0, tx, 0, 256, acc);
        } else {
            for (int m=0; m<n; m++) {
                int jlo, jhi;
                valid_j(k0[m], np, &jlo, &jhi);
                turn_sums<1>(id, qd, &k0[m], &tx[m], jlo, jhi, &acc[m]);
            }
        }
        for (int m=0; m<n; m++) {
            const int i=i0+m;
            float a[8];
            v8_store(a, acc[m]);
            for (int t=0; t<4; t++) {
                cf[t][i]=row_c(tx[m], 256, t); sf[t][i]=row_s(tx[m], 256, t);
                is[t][i]=a[t]; qs[t][i]=a[4+t];
            }
        }
    }
}

// live_sums for S symbols side by side, each with its own table shared by
// N lags: sum (s, m) turns samples k0[m]+s*256+j by table s.
template <int S, int N>
void live_lags(const float* id, const float* qd, const long* k0, const float* fp,
               TurnSteps* steps, V8* acc, V8* x256) {
    V8 x[S], cd[S], sd[S], a[S*N];
    for (int t=0; t<S; t++) turn_start(fp[t], &x[t], &cd[t], &sd[t], steps);
    for (int m=0; m<S*N; m++) a[m]=v8_zero();
    for (int j=0; j<256; j++) {
        for (int t=0; t<S; t++) {
            V8 cs, sc;
            turn_row(&x[t], cd[t], sd[t], &cs, &sc);
            for (int m=0; m<N; m++) {
                const long k=k0[m]+t*256L+j;
                const V8 vi=v8_dup(id+k), vq=v8_dup(qd+k);
                a[t*N+m]=turn_add(a[t*N+m], vi, vq, cs, sc);
            }
        }
    }
    for (int m=0; m<S*N; m++) acc[m]=a[m];
    for (int t=0; t<S; t++) x256[t]=x[t];
}

// noncoherent_sums of up to four lags at once (the next shifts a try steps
// through) into the memo, one table per symbol for them all: with a drift
// each lag alone makes all 162 tables again. A lag the memo has, or no
// room for, is left out.
void noncoherent_lags(float *id, float *qd, long np, float f0, const int* lag, int nlag,
                      float drift1, NcMemo* memo)
{
    enum { MAXL=4 };
    NcMemo::Sums* out[MAXL];
    long k0[MAXL];
    int n=0;
    for (int l=0; l<nlag && n<MAXL; l++) {
        if (memo->find(f0, drift1, lag[l])) continue;
        bool dup=false;
        for (int m=0; m<n; m++) dup=dup || out[m]->lag == lag[l];
        if (dup) continue;
        NcMemo::Sums* e=memo->add(f0, drift1, lag[l]);
        if (!e) break;
        out[n++]=e;
    }
    if (n == 0) return;
    Turns tab;  // an edge symbol's
    // all of symbol i's lags summed whole
    auto full = [&](int i) {
        for (int m=0; m<n; m++) {
            int jlo, jhi;
            valid_j(out[m]->lag+i*256L, np, &jlo, &jhi);
            if (jlo != 0 || jhi != 256) return false;
        }
        return true;
    };
    auto put = [&](int m, int i, V8 acc, const float x[8]) {
        float a[8];
        v8_store(a, acc);
        for (int t=0; t<4; t++) {
            out[m]->cf[t][i]=x[t]; out[m]->sf[t][i]=x[4+t];
            out[m]->is[t][i]=a[t]; out[m]->qs[t][i]=a[4+t];
        }
    };
    for (int i=0; i<162; i++) {
        for (int m=0; m<n; m++) k0[m]=out[m]->lag+i*256L;
        if (i+1 < 162 && full(i) && full(i+1)) {  // two symbols, a table each
            float fp[2];
            for (int t=0; t<2; t++) fp[t] = f0 + (drift1/2.0)*((float)(i+t)-81.0)/81.0;
            V8 acc[2*MAXL], x256[2];
            switch (n) {
                case 4: live_lags<2, 4>(id, qd, k0, fp, memo->steps, acc, x256); break;
                case 3: live_lags<2, 3>(id, qd, k0, fp, memo->steps, acc, x256); break;
                case 2: live_lags<2, 2>(id, qd, k0, fp, memo->steps, acc, x256); break;
                default: live_lags<2, 1>(id, qd, k0, fp, memo->steps, acc, x256); break;
            }
            for (int t=0; t<2; t++) {
                float x[8];
                v8_store(x, x256[t]);
                for (int m=0; m<n; m++) put(m, i+t, acc[t*n+m], x);
            }
            i++;
            continue;
        }
        float fp = f0 + (drift1/2.0)*((float)i-81.0)/81.0;
        Turns* t=&tab;
        const Turns* tx=&tab;
        tab.fp=fp;
        build_turns(&t, 1, 257, memo->steps);
        float x[8];
        for (int t4=0; t4<4; t4++) { x[t4]=row_c(&tab, 256, t4); x[4+t4]=row_s(&tab, 256, t4); }
        for (int m=0; m<n; m++) {
            int jlo, jhi;
            V8 acc;
            valid_j(k0[m], np, &jlo, &jhi);
            turn_sums<1>(id, qd, &k0[m], &tx, jlo, jhi, &acc);
            put(m, i, acc, x);
        }
    }
}
#endif

#ifdef WSPRD_AVX2
// The coarse search's sync of one (ifr, idrift) for all 32 time offsets
// k0=-10..21 at once, lane l of vector v being k0=-10+8v+l: every lane sums
// its symbols in the reference's order, a symbol whose kindex is off the
// spectrogram leaving its sums as they were (masked, never read). ps holds
// the square roots already; ifd[k] is the reference's ifd of symbol k.
// Symbol k's step of coarse_sync; EDGE: some vector may reach off the
// spectrogram (the symbols in between have all four on it).
template <bool EDGE>
inline void coarse_step(const float* ps, int nffts, const int* ifd, int k, __m256* ss, __m256* pw) {
    const __m256i lane=_mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7);
    const __m256i n=_mm256_set1_epi32(nffts), neg1=_mm256_set1_epi32(-1);
    const float* r0=ps+(long)(ifd[k]-3)*nffts;
    const float* r1=ps+(long)(ifd[k]-1)*nffts;
    const float* r2=ps+(long)(ifd[k]+1)*nffts;
    const float* r3=ps+(long)(ifd[k]+3)*nffts;
    const __m256 sign=_mm256_set1_ps((float)(2*pr3[k]-1));
    for (int v=0; v<4; v++) {
        const int kk=-10+8*v+2*k;  // lane 0's kindex
        if (EDGE && (kk >= nffts || kk+8 <= 0)) continue;
        __m256 p0, p1, p2, p3;
        const bool full=!EDGE || (kk >= 0 && kk+8 <= nffts);
        __m256i m=neg1;
        if (full) {
            p0=_mm256_loadu_ps(r0+kk); p1=_mm256_loadu_ps(r1+kk);
            p2=_mm256_loadu_ps(r2+kk); p3=_mm256_loadu_ps(r3+kk);
        } else {
            const __m256i ki=_mm256_add_epi32(_mm256_set1_epi32(kk), lane);
            m=_mm256_and_si256(_mm256_cmpgt_epi32(ki, neg1), _mm256_cmpgt_epi32(n, ki));
            // the pointers of off-spectrogram lanes are never dereferenced
            p0=_mm256_maskload_ps(r0+kk, m); p1=_mm256_maskload_ps(r1+kk, m);
            p2=_mm256_maskload_ps(r2+kk, m); p3=_mm256_maskload_ps(r3+kk, m);
        }
        const __m256 d=_mm256_sub_ps(_mm256_add_ps(p1, p3), _mm256_add_ps(p0, p2));
        const __m256 s=_mm256_add_ps(ss[v], _mm256_mul_ps(sign, d));
        const __m256 w=_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(pw[v], p0), p1), p2), p3);
        if (full) {
            ss[v]=s; pw[v]=w;
        } else {
            const __m256 mf=_mm256_castsi256_ps(m);
            ss[v]=_mm256_blendv_ps(ss[v], s, mf);
            pw[v]=_mm256_blendv_ps(pw[v], w, mf);
        }
    }
}

void coarse_sync(const float* ps, int nffts, const int* ifd, float sync[32]) {
    __m256 ss[4], pw[4];
    for (int v=0; v<4; v++) ss[v]=pw[v]=_mm256_setzero_ps();
    // all four vectors on the spectrogram, kindex -10+2k..2k+21: k from
    // 5 to (nffts-22)/2
    const int klo=5, khi=nffts >= 32 ? (nffts-22)/2+1 : klo;
    int k=0;
    for (; k<162 && k<klo; k++) coarse_step<true>(ps, nffts, ifd, k, ss, pw);
    for (; k<162 && k<khi; k++) coarse_step<false>(ps, nffts, ifd, k, ss, pw);
    for (; k<162; k++) coarse_step<true>(ps, nffts, ifd, k, ss, pw);
    for (int v=0; v<4; v++) _mm256_storeu_ps(sync+8*v, _mm256_div_ps(ss[v], pw[v]));
}
#elif defined(WSPRD_NEON)
// As the AVX2 coarse_sync, in eight NEON vectors: lane l of vector v is
// k0=-10+4v+l. An edge vector's lanes off the spectrogram load 0 and keep
// their sums as they were (vbslq).
template <bool EDGE>
inline void coarse_step(const float* ps, int nffts, const int* ifd, int k,
                        float32x4_t* ss, float32x4_t* pw) {
    const float* r[4]={ps+(long)(ifd[k]-3)*nffts, ps+(long)(ifd[k]-1)*nffts,
                       ps+(long)(ifd[k]+1)*nffts, ps+(long)(ifd[k]+3)*nffts};
    const float32x4_t sign=vdupq_n_f32((float)(2*pr3[k]-1));
    for (int v=0; v<8; v++) {
        const int kk=-10+4*v+2*k;  // lane 0's kindex
        if (EDGE && (kk >= nffts || kk+4 <= 0)) continue;
        float32x4_t p0, p1, p2, p3;
        const bool full=!EDGE || (kk >= 0 && kk+4 <= nffts);
        uint32x4_t m=vdupq_n_u32(~0u);
        if (full) {
            p0=vld1q_f32(r[0]+kk); p1=vld1q_f32(r[1]+kk);
            p2=vld1q_f32(r[2]+kk); p3=vld1q_f32(r[3]+kk);
        } else {
            float t[4][4];
            uint32_t on[4];
            for (int l=0; l<4; l++) {
                const int ki=kk+l;
                on[l]=ki >= 0 && ki < nffts ? ~0u : 0u;
                for (int q=0; q<4; q++) t[q][l]=on[l] ? r[q][ki] : 0.0f;
            }
            p0=vld1q_f32(t[0]); p1=vld1q_f32(t[1]); p2=vld1q_f32(t[2]); p3=vld1q_f32(t[3]);
            m=vld1q_u32(on);
        }
        const float32x4_t d=vsubq_f32(vaddq_f32(p1, p3), vaddq_f32(p0, p2));
        const float32x4_t s=vaddq_f32(ss[v], vmulq_f32(sign, d));
        const float32x4_t w=vaddq_f32(vaddq_f32(vaddq_f32(vaddq_f32(pw[v], p0), p1), p2), p3);
        if (full) {
            ss[v]=s; pw[v]=w;
        } else {
            ss[v]=vbslq_f32(m, s, ss[v]);
            pw[v]=vbslq_f32(m, w, pw[v]);
        }
    }
}

void coarse_sync(const float* ps, int nffts, const int* ifd, float sync[32]) {
    float32x4_t ss[8], pw[8];
    for (int v=0; v<8; v++) ss[v]=pw[v]=vdupq_n_f32(0.0f);
    const int klo=5, khi=nffts >= 32 ? (nffts-22)/2+1 : klo;  // as the AVX2 one's
    int k=0;
    for (; k<162 && k<klo; k++) coarse_step<true>(ps, nffts, ifd, k, ss, pw);
    for (; k<162 && k<khi; k++) coarse_step<false>(ps, nffts, ifd, k, ss, pw);
    for (; k<162; k++) coarse_step<true>(ps, nffts, ifd, k, ss, pw);
    for (int v=0; v<8; v++) vst1q_f32(sync+4*v, vdivq_f32(ss[v], pw[v]));
}
#endif

//***************************************************************************
void sync_and_demodulate(float *id, float *qd, long np,
                         unsigned char *symbols, float *f1, int ifmin, int ifmax, float fstep,
                         int *shift1, int lagmin, int lagmax, int lagstep,
                         float *drift1, int symfac, float *sync, int mode, SyncMemo* memo)
{
    /***********************************************************************
     * mode = 0: no frequency or drift search. find best time lag.          *
     *        1: no time lag or drift search. find best frequency.          *
     *        2: no frequency or time lag search. calculate soft-decision   *
     *           symbols using passed frequency and shift.                  *
     ************************************************************************/

    float fplast=-10000.0;
    const float dt=1.0/375.0, df=375.0/256.0;
    const float pi=3.14159265358979323846;
    float twopidt, df15=df*1.5, df05=df*0.5;

    int i, j, k, lag;
    float i0[162],q0[162],i1[162],q1[162],i2[162],q2[162],i3[162],q3[162];
    float p0,p1,p2,p3,cmet,totp,syncmax,fac;
    float c0[256],s0[256],c1[256],s1[256],c2[256],s2[256],c3[256],s3[256];
    float dphi0, cdphi0, sdphi0, dphi1, cdphi1, sdphi1, dphi2, cdphi2, sdphi2,
    dphi3, cdphi3, sdphi3;
    float f0=0.0, fp, ss, fbest=0.0, fsum=0.0, f2sum=0.0, fsymb[162];
    int best_shift = 0, ifreq;

#ifdef WSPRD_SIMD
    if( mode <= 1 && sync_search(id, qd, np, f1, ifmin, ifmax, fstep, shift1, lagmin, lagmax,
                                 lagstep, drift1, sync, mode, memo) ) return;
#else
    (void)memo;
#endif
    syncmax=-1e30;
    if( mode == 0 ) {ifmin=0; ifmax=0; fstep=0.0; f0=*f1;}
    if( mode == 1 ) {lagmin=*shift1;lagmax=*shift1;f0=*f1;}
    if( mode == 2 ) {lagmin=*shift1;lagmax=*shift1;ifmin=0;ifmax=0;f0=*f1;}

    twopidt=2*pi*dt;
    for(ifreq=ifmin; ifreq<=ifmax; ifreq++) {
        f0=*f1+ifreq*fstep;
        for(lag=lagmin; lag<=lagmax; lag=lag+lagstep) {
            ss=0.0;
            totp=0.0;
            for (i=0; i<162; i++) {
                fp = f0 + (*drift1/2.0)*((float)i-81.0)/81.0;
                if( i==0 || (fp != fplast) ) {  // only calculate sin/cos if necessary
                    dphi0=twopidt*(fp-df15);
                    cdphi0=cos((double)dphi0);
                    sdphi0=sin((double)dphi0);

                    dphi1=twopidt*(fp-df05);
                    cdphi1=cos((double)dphi1);
                    sdphi1=sin((double)dphi1);

                    dphi2=twopidt*(fp+df05);
                    cdphi2=cos((double)dphi2);
                    sdphi2=sin((double)dphi2);

                    dphi3=twopidt*(fp+df15);
                    cdphi3=cos((double)dphi3);
                    sdphi3=sin((double)dphi3);

                    c0[0]=1; s0[0]=0;
                    c1[0]=1; s1[0]=0;
                    c2[0]=1; s2[0]=0;
                    c3[0]=1; s3[0]=0;

                    for (j=1; j<256; j++) {
                        c0[j]=c0[j-1]*cdphi0 - s0[j-1]*sdphi0;
                        s0[j]=c0[j-1]*sdphi0 + s0[j-1]*cdphi0;
                        c1[j]=c1[j-1]*cdphi1 - s1[j-1]*sdphi1;
                        s1[j]=c1[j-1]*sdphi1 + s1[j-1]*cdphi1;
                        c2[j]=c2[j-1]*cdphi2 - s2[j-1]*sdphi2;
                        s2[j]=c2[j-1]*sdphi2 + s2[j-1]*cdphi2;
                        c3[j]=c3[j-1]*cdphi3 - s3[j-1]*sdphi3;
                        s3[j]=c3[j-1]*sdphi3 + s3[j-1]*cdphi3;
                    }
                    fplast = fp;
                }

                i0[i]=0.0; q0[i]=0.0;
                i1[i]=0.0; q1[i]=0.0;
                i2[i]=0.0; q2[i]=0.0;
                i3[i]=0.0; q3[i]=0.0;

                for (j=0; j<256; j++) {
                    k=lag+i*256+j;
                    if( (k>0) && (k<np) ) {
                        i0[i]=i0[i] + id[k]*c0[j] + qd[k]*s0[j];
                        q0[i]=q0[i] - id[k]*s0[j] + qd[k]*c0[j];
                        i1[i]=i1[i] + id[k]*c1[j] + qd[k]*s1[j];
                        q1[i]=q1[i] - id[k]*s1[j] + qd[k]*c1[j];
                        i2[i]=i2[i] + id[k]*c2[j] + qd[k]*s2[j];
                        q2[i]=q2[i] - id[k]*s2[j] + qd[k]*c2[j];
                        i3[i]=i3[i] + id[k]*c3[j] + qd[k]*s3[j];
                        q3[i]=q3[i] - id[k]*s3[j] + qd[k]*c3[j];
                    }
                }
                p0=i0[i]*i0[i] + q0[i]*q0[i];
                p1=i1[i]*i1[i] + q1[i]*q1[i];
                p2=i2[i]*i2[i] + q2[i]*q2[i];
                p3=i3[i]*i3[i] + q3[i]*q3[i];

                p0=sqrt((double)p0);
                p1=sqrt((double)p1);
                p2=sqrt((double)p2);
                p3=sqrt((double)p3);

                totp=totp+p0+p1+p2+p3;
                cmet=(p1+p3)-(p0+p2);
                ss = (pr3[i] == 1) ? ss+cmet : ss-cmet;
                if( mode == 2) {                 //Compute soft symbols
                    if(pr3[i]==1) {
                        fsymb[i]=p3-p1;
                    } else {
                        fsymb[i]=p2-p0;
                    }
                }
            }
            ss=ss/totp;
            if( ss > syncmax ) {          //Save best parameters
                syncmax=ss;
                best_shift=lag;
                fbest=f0;
            }
        } // lag loop
    } //freq loop

    if( mode <=1 ) {                       //Send best params back to caller
        *sync=syncmax;
        *shift1=best_shift;
        *f1=fbest;
        return;
    }

    if( mode == 2 ) {
        *sync=syncmax;
        for (i=0; i<162; i++) {              //Normalize the soft symbols
            fsum=fsum+fsymb[i]/162.0;
            f2sum=f2sum+fsymb[i]*fsymb[i]/162.0;
        }
        fac=sqrt((double)(f2sum-fsum*fsum));
        for (i=0; i<162; i++) {
            fsymb[i]=symfac*fsymb[i]/fac;
            if( fsymb[i] > 127) fsymb[i]=127.0;
            if( fsymb[i] < -128 ) fsymb[i]=-128.0;
            symbols[i]=fsymb[i] + 128;
        }
        return;
    }
    return;
}

inline void noncoherent_sequence_detection(float *id, float *qd, long np,
                                    unsigned char *symbols, float *f1,  int *shift1,
                                    float *drift1, int symfac, int *nblocksize, int *bitmetric,
                                    NcMemo* memo)
{
    /************************************************************************
     *  Noncoherent sequence detection for wspr.                            *
     *  Allowed block lengths are nblock=1,2,3,6, or 9 symbols.             *
     *  Longer block lengths require longer channel coherence time.         *
     *  The whole block is estimated at once.                               *
     *  nblock=1 corresponds to noncoherent detection of individual symbols *
     *     like the original wsprd symbol demodulator.                      *
     ************************************************************************/
#ifndef WSPRD_SIMD  // the sums' (noncoherent_sums has its own)
    float fplast=-10000.0;
    const float dt=1.0/375.0, df=375.0/256.0;
    const float pi=3.14159265358979323846;
    float twopidt, df15=df*1.5, df05=df*0.5;
    int k;
    float c0[257],s0[257],c1[257],s1[257],c2[257],s2[257],c3[257],s3[257];
    float dphi0, cdphi0, sdphi0, dphi1, cdphi1, sdphi1, dphi2, cdphi2, sdphi2,
    dphi3, cdphi3, sdphi3;
    float fp;

    twopidt=2*pi*dt;
#endif
    int i, j, lag, itone, ib, b, nblock, nseq, imask;
    float xi[512],xq[512];
    float is[4][162],qs[4][162],cf[4][162],sf[4][162],cm,sm,cmp,smp;
    float p[512],fac,xm1,xm0;
    float f0, fsum=0.0, f2sum=0.0, fsymb[162];

    f0=*f1;
    lag=*shift1;
    nblock=*nblocksize;
    nseq=1<<nblock;
    int bitbybit=*bitmetric;

#ifdef WSPRD_SIMD
    if (const NcMemo::Sums* m = memo ? memo->find(f0, *drift1, lag) : nullptr) {
        memcpy(is, m->is, sizeof is); memcpy(qs, m->qs, sizeof qs);
        memcpy(cf, m->cf, sizeof cf); memcpy(sf, m->sf, sizeof sf);
    } else {
        noncoherent_sums(id, qd, np, f0, lag, *drift1, is, qs, cf, sf, memo ? memo->steps : nullptr);
        if (NcMemo::Sums* e = memo ? memo->add(f0, *drift1, lag) : nullptr) {
            memcpy(e->is, is, sizeof is); memcpy(e->qs, qs, sizeof qs);
            memcpy(e->cf, cf, sizeof cf); memcpy(e->sf, sf, sizeof sf);
        }
    }
#else
    (void)memo;
    for (i=0; i<162; i++) {
        fp = f0 + (*drift1/2.0)*((float)i-81.0)/81.0;
        if( i==0 || (fp != fplast) ) {  // only calculate sin/cos if necessary
            dphi0=twopidt*(fp-df15);
            cdphi0=cos((double)dphi0);
            sdphi0=sin((double)dphi0);

            dphi1=twopidt*(fp-df05);
            cdphi1=cos((double)dphi1);
            sdphi1=sin((double)dphi1);

            dphi2=twopidt*(fp+df05);
            cdphi2=cos((double)dphi2);
            sdphi2=sin((double)dphi2);

            dphi3=twopidt*(fp+df15);
            cdphi3=cos((double)dphi3);
            sdphi3=sin((double)dphi3);

            c0[0]=1; s0[0]=0;
            c1[0]=1; s1[0]=0;
            c2[0]=1; s2[0]=0;
            c3[0]=1; s3[0]=0;

            for (j=1; j<257; j++) {
                c0[j]=c0[j-1]*cdphi0 - s0[j-1]*sdphi0;
                s0[j]=c0[j-1]*sdphi0 + s0[j-1]*cdphi0;
                c1[j]=c1[j-1]*cdphi1 - s1[j-1]*sdphi1;
                s1[j]=c1[j-1]*sdphi1 + s1[j-1]*cdphi1;
                c2[j]=c2[j-1]*cdphi2 - s2[j-1]*sdphi2;
                s2[j]=c2[j-1]*sdphi2 + s2[j-1]*cdphi2;
                c3[j]=c3[j-1]*cdphi3 - s3[j-1]*sdphi3;
                s3[j]=c3[j-1]*sdphi3 + s3[j-1]*cdphi3;
            }

            fplast = fp;
        }

        cf[0][i]=c0[256]; sf[0][i]=s0[256];
        cf[1][i]=c1[256]; sf[1][i]=s1[256];
        cf[2][i]=c2[256]; sf[2][i]=s2[256];
        cf[3][i]=c3[256]; sf[3][i]=s3[256];

        is[0][i]=0.0; qs[0][i]=0.0;
        is[1][i]=0.0; qs[1][i]=0.0;
        is[2][i]=0.0; qs[2][i]=0.0;
        is[3][i]=0.0; qs[3][i]=0.0;

        for (j=0; j<256; j++) {
            k=lag+i*256+j;
            if( (k>0) && (k<np) ) {
                is[0][i]=is[0][i] + id[k]*c0[j] + qd[k]*s0[j];
                qs[0][i]=qs[0][i] - id[k]*s0[j] + qd[k]*c0[j];
                is[1][i]=is[1][i] + id[k]*c1[j] + qd[k]*s1[j];
                qs[1][i]=qs[1][i] - id[k]*s1[j] + qd[k]*c1[j];
                is[2][i]=is[2][i] + id[k]*c2[j] + qd[k]*s2[j];
                qs[2][i]=qs[2][i] - id[k]*s2[j] + qd[k]*c2[j];
                is[3][i]=is[3][i] + id[k]*c3[j] + qd[k]*s3[j];
                qs[3][i]=qs[3][i] - id[k]*s3[j] + qd[k]*c3[j];
            }
        }
    }
#endif

    for (i=0; i<162; i=i+nblock) {
        for (j=0;j<nseq;j++) {
            xi[j]=0.0; xq[j]=0.0;
            cm=1; sm=0;
            for (ib=0; ib<nblock; ib++) {
                b=(j&(1<<(nblock-1-ib)))>>(nblock-1-ib);
                itone=pr3[i+ib]+2*b;
                xi[j]=xi[j]+is[itone][i+ib]*cm + qs[itone][i+ib]*sm;
                xq[j]=xq[j]+qs[itone][i+ib]*cm - is[itone][i+ib]*sm;
                cmp=cf[itone][i+ib]*cm - sf[itone][i+ib]*sm;
                smp=sf[itone][i+ib]*cm + cf[itone][i+ib]*sm;
                cm=cmp; sm=smp;
            }
            p[j]=xi[j]*xi[j]+xq[j]*xq[j];
            p[j]=sqrt((double)p[j]);
        }
        for (ib=0; ib<nblock; ib++) {
            imask=1<<(nblock-1-ib);
            xm1=0.0; xm0=0.0;
            for (j=0; j<nseq; j++) {
                if((j & imask)!=0) {
                    if(p[j] > xm1) xm1=p[j];
                }
                if((j & imask)==0) {
                    if(p[j]>xm0) xm0=p[j];
                }
            }
            fsymb[i+ib]=xm1-xm0;
            if( bitbybit == 1 ) {
                fsymb[i+ib]=fsymb[i+ib]/(xm1 > xm0 ? xm1 : xm0);
            }
        }
    }
    for (i=0; i<162; i++) {              //Normalize the soft symbols
        fsum=fsum+fsymb[i]/162.0;
        f2sum=f2sum+fsymb[i]*fsymb[i]/162.0;
    }
    fac=sqrt((double)(f2sum-fsum*fsum));
    for (i=0; i<162; i++) {
        fsymb[i]=symfac*fsymb[i]/fac;
        if( fsymb[i] > 127) fsymb[i]=127.0;
        if( fsymb[i] < -128 ) fsymb[i]=-128.0;
        symbols[i]=fsymb[i] + 128;
    }
    return;
}

/******************************************************************************
  Subtract the coherent component of a signal
 *******************************************************************************/
// subtract_signal2's arrays, kept for the next call (made and zeroed once,
// ci and cq zeroed again each call: the rest is written before it is read).
struct SubScratch {
    std::vector<float> refi, refq, ci, cq, cfi, cfq;
    // subtract_signal2's filter length: s(t)*conjugate(r(t)) for the
    // reference's sample i (the data's shift0+i) is at ci/cq[nfilt+i]
    // (spread_gain() reads it there).
    static constexpr int nfilt = 360;
};

void subtract_signal2(float *id, float *qd, long np,
                      float f0, int shift0, float drift0, unsigned char* channel_symbols, SubScratch& sc)
{
    float dt=1.0/375.0, df=375.0/256.0;
    float pi=4.*atan(1.0), twopidt, phi=0, dphi, cs;
    int i, j, k, ii, nsym=162, nspersym=256;
    enum { nfilt=360 }; //nfilt must be even number.
    static_assert(nfilt == SubScratch::nfilt, "spread_gain's layout");
    int nsig=nsym*nspersym;
    int nc2=45000;

    if (sc.ci.empty()) {
        for (auto* v : {&sc.refi, &sc.refq, &sc.ci, &sc.cq, &sc.cfi, &sc.cfq}) v->assign(nc2, 0.0f);
    } else {
        std::fill(sc.ci.begin(), sc.ci.end(), 0.0f);
        std::fill(sc.cq.begin(), sc.cq.end(), 0.0f);
    }
    float *refi=sc.refi.data(), *refq=sc.refq.data(), *ci=sc.ci.data(), *cq=sc.cq.data(),
          *cfi=sc.cfi.data(), *cfq=sc.cfq.data();

    twopidt=2.0*pi*dt;

    /******************************************************************************
     Measured signal:                    s(t)=a(t)*exp( j*theta(t) )
     Reference is:                       r(t) = exp( j*phi(t) )
     Complex amplitude is estimated as:  c(t)=LPF[s(t)*conjugate(r(t))]
     so c(t) has phase angle theta-phi
     Multiply r(t) by c(t) and subtract from s(t), i.e. s'(t)=s(t)-c(t)r(t)
     *******************************************************************************/

    // create reference wspr signal vector, centered on f0.
    //
    for (i=0; i<nsym; i++) {

        cs=(float)channel_symbols[i];

        dphi=twopidt*
        (
         f0 + (drift0/2.0)*((float)i-(float)nsym/2.0)/((float)nsym/2.0)
         + (cs-1.5)*df
         );

        for ( j=0; j<nspersym; j++ ) {
            ii=nspersym*i+j;
            refi[ii]=cos((double)phi); //cannot precompute sin/cos because dphi is changing
            refq[ii]=sin((double)phi);
            phi=phi+dphi;
        }
    }

    float w[nfilt], norm=0, partialsum[nfilt];
    //lowpass filter and remove startup transient
    for (i=0; i<nfilt; i++) partialsum[i]=0.0;
    for (i=0; i<nfilt; i++) {
        w[i]=sin((double)(pi*(float)i/(float)(nfilt-1)));
        norm=norm+w[i];
    }
    for (i=0; i<nfilt; i++) {
        w[i]=w[i]/norm;
    }
    for (i=1; i<nfilt; i++) {
        partialsum[i]=partialsum[i-1]+w[i];
    }

    // s(t) * conjugate(r(t))
    // beginning of first symbol in reference signal is at i=0
    // beginning of first symbol in received data is at shift0.
    // filter transient lasts nfilt samples
    // leave nfilt zeros as a pad at the beginning of the unfiltered reference signal
    for (i=0; i<nsym*nspersym; i++) {
        k=shift0+i;
        if( (k>0) && (k<np) ) {
            ci[i+nfilt] = id[k]*refi[i] + qd[k]*refq[i];
            cq[i+nfilt] = qd[k]*refi[i] - id[k]*refq[i];
        }
    }

    // LPF
#ifdef WSPRD_AVX2
    // 32 outputs at a time, a lane each: every lane's sum the reference's,
    // its j in order; only the outputs the subtraction below reads (its
    // j=nfilt..nfilt+nsig-1, 1296 times 32), not the 3168 others
    static_assert(162*256%32 == 0, "LPF lanes");
    for (i=nfilt; i<nfilt+nsig; i+=32) {
        __m256 ai[4], aq[4];
        for (int v=0; v<4; v++) ai[v]=aq[v]=_mm256_setzero_ps();
        const float* pi_=&ci[i-nfilt/2];
        const float* pq=&cq[i-nfilt/2];
        for (j=0; j<nfilt; j++) {
            const __m256 wj=_mm256_broadcast_ss(&w[j]);
            for (int v=0; v<4; v++) {
                ai[v]=_mm256_add_ps(ai[v], _mm256_mul_ps(wj, _mm256_loadu_ps(pi_+j+8*v)));
                aq[v]=_mm256_add_ps(aq[v], _mm256_mul_ps(wj, _mm256_loadu_ps(pq+j+8*v)));
            }
        }
        for (int v=0; v<4; v++) {
            _mm256_storeu_ps(&cfi[i+8*v], ai[v]);
            _mm256_storeu_ps(&cfq[i+8*v], aq[v]);
        }
    }
#elif defined(WSPRD_NEON)
    // As the AVX2 filter, 16 outputs at a time, a lane each
    static_assert(162*256%16 == 0, "LPF lanes");
    for (i=nfilt; i<nfilt+nsig; i+=16) {
        float32x4_t ai[4], aq[4];
        for (int v=0; v<4; v++) ai[v]=aq[v]=vdupq_n_f32(0.0f);
        const float* pi_=&ci[i-nfilt/2];
        const float* pq=&cq[i-nfilt/2];
        for (j=0; j<nfilt; j++) {
            const float wj=w[j];
            for (int v=0; v<4; v++) {
                ai[v]=vaddq_f32(ai[v], vmulq_n_f32(vld1q_f32(pi_+j+4*v), wj));
                aq[v]=vaddq_f32(aq[v], vmulq_n_f32(vld1q_f32(pq+j+4*v), wj));
            }
        }
        for (int v=0; v<4; v++) {
            vst1q_f32(&cfi[i+4*v], ai[v]);
            vst1q_f32(&cfq[i+4*v], aq[v]);
        }
    }
#else
    // Each sum in locals: gcc can't tell that cfi and cfq don't overlap ci
    // and cq (as wsprd.c's calloc'd arrays), so would store every term
    for (i=nfilt/2; i<45000-nfilt/2; i++) {
        float ai=0.0, aq=0.0;
        for (j=0; j<nfilt; j++) {
            ai=ai+w[j]*ci[i-nfilt/2+j];
            aq=aq+w[j]*cq[i-nfilt/2+j];
        }
        cfi[i]=ai; cfq[i]=aq;
    }
#endif

    // subtract c(t)*r(t) here
    // (ci+j*cq)(refi+j*refq)=(ci*refi-cq*refq)+j(ci*refq+cq*refi)
    // beginning of first symbol in reference signal is at i=nfilt
    // beginning of first symbol in received data is at shift0.
    for (i=0; i<nsig; i++) {
        if( i<nfilt/2 ) {        // take care of the end effect (LPF step response) here
            norm=partialsum[nfilt/2+i];
        } else if( i>(nsig-1-nfilt/2) ) {
            norm=partialsum[nfilt/2+nsig-1-i];
        } else {
            norm=1.0;
        }
        k=shift0+i;
        j=i+nfilt;
        if( (k>0) && (k<np) ) {
            id[k]=id[k] - (cfi[j]*refi[i]-cfq[j]*refq[i])/norm;
            qd[k]=qd[k] - (cfi[j]*refq[i]+cfq[j]*refi[i])/norm;
        }
    }

    return;
}

/******************************************************************************
  Not wsprd's: a decode's spread, as WSJT-X measures FST4W's (dopspread in
  lib/fst4_decode.f90, its w50): nothing wsprd prints or keeps reads it.
 *******************************************************************************/
// The channel's gain g(t) = s(t)*conjugate(r(t)) (dopspread's g, the QEX
// paper's), r the decoded signal's reference as subtract_signal2 makes it
// (f0, shift0, drift0), at the data's own sample indices (dopspread's
// cwave cshifted to the signal), zero elsewhere, in a buffer of nfft. From
// subtract_signal2's ci/cq when it ran for this signal (the same products,
// before it subtracted), else made here as it makes them.
struct SpreadScratch {
    std::unique_ptr<skimmer::Fft> fft;  // made at first use (quick: fft_load_wisdom()'s plan, else FFTW_ESTIMATE)
};

void spread_gain(const float *id, const float *qd, long np, float f0, int shift0, float drift0,
                 const unsigned char* channel_symbols, const SubScratch* sub, std::complex<float>* g)
{
    const int nsym=162, nspersym=256, nsig=nsym*nspersym;
    if (sub) {
        for (int i=0; i<nsig; i++) {
            const int k=shift0+i;
            if( (k>0) && (k<np) ) g[k]={sub->ci[SubScratch::nfilt+i], sub->cq[SubScratch::nfilt+i]};
        }
        return;
    }
    // subtract_signal2's reference, sample by sample, and its products.
    float dt=1.0/375.0, df=375.0/256.0;
    float pi=4.*atan(1.0), twopidt=2.0*pi*dt, phi=0, dphi, cs;
    for (int i=0; i<nsym; i++) {
        cs=(float)channel_symbols[i];
        dphi=twopidt*
        (
         f0 + (drift0/2.0)*((float)i-(float)nsym/2.0)/((float)nsym/2.0)
         + (cs-1.5)*df
         );
        for (int j=0; j<nspersym; j++) {
            const float refi=cos((double)phi), refq=sin((double)phi);
            const int k=shift0+nspersym*i+j;
            if( (k>0) && (k<np) ) g[k]={id[k]*refi + qd[k]*refq, qd[k]*refi - id[k]*refq};
            phi=phi+dphi;
        }
    }
}

// dopspread from its FFT on, every float as the Fortran's: g in the
// buffer (nfft points, fs the samples' rate), transformed in place.
float dopspread(skimmer::Fft& fft, float fs)
{
    const int nfft=(int)fft.size();
    fft.forward();
    const std::complex<float>* g=fft.data();
    const auto power=[&](int i) {
        const std::complex<float> x=g[i<0 ? i+nfft : i];
        return x.real()*x.real() + x.imag()*x.imag();
    };

    const float df=fs/(float)nfft;
    int ia=(int)(1.0f/df);
    float smax=0.;
    for (int i=-ia; i<=ia; i++) {                  // Find smax in +/- 1 Hz around 0.
        const float s=power(i);
        smax=s > smax ? s : smax;
    }

    ia=(int)(10.1f/df);
    std::vector<float> ssv(2*ia+1);                // Allocate space for +/- 10 Hz
    float* ss=ssv.data()+ia;
    float sum1=0., sum2=0.;
    int nns=0;
    for (int i=-ia; i<=ia; i++) {
        ss[i]=power(i)/smax;
        const float f=(float)i*df;
        if (f >= -4.0f && f <= -2.0f) {
            sum1=sum1 + ss[i];                     // Power between -2 and -4 Hz
            nns=nns+1;
        } else if (f >= 2.0f && f <= 4.0f) {
            sum2=sum2 + ss[i];                     // Power between +2 and +4 Hz
        }
    }
    const float avg=std::min(sum1/(float)nns, sum2/(float)nns);  // Compute avg from smaller sum

    sum1=0.;
    for (int i=-ia; i<=ia; i++) {
        const float f=(float)i*df;
        if (std::fabs(f) <= 1.0f) sum1=sum1 + ss[i]-avg;  // Power in abs(f) < 1 Hz
    }

    ia=(int)std::lround(1.0f/df) + 1;
    sum2=0.0;
    float xi1=-999, xi3=-999, sum2z=0.;
    for (int i=-ia; i<=ia; i++) {                  // Find freq range that has 50% of signal power
        sum2=sum2 + ss[i]-avg;
        if (sum2 >= 0.25f*sum1 && xi1 == -999.0f) xi1=(float)(i - 1) + (0.25f*sum1-sum2)/(sum2-sum2z);
        // (xi2, the midpoint, is only plotted)
        if (sum2 >= 0.75f*sum1) {
            xi3=(float)(i - 1) + (sum2-0.75f*sum1)/(sum2-sum2z);
            break;
        }
        sum2z=sum2;
    }
    const float xdiff=std::sqrt(1.0f+(xi3-xi1)*(xi3-xi1));  // Keep small values from fluctuating too widely
    return xdiff*df;                               // Compute Doppler spread
}

// A decode's spread (w50), Hz: its gain over the data's samples (np at
// 375 Hz, 46.875 for WSPR-15), zero-padded to twice its period as
// dopspread pads nmax: so the bins are FST4W-120's (WSPR-2, 1/240 Hz) and
// FST4W-900's (WSPR-15, 1/1800 Hz). sub: subtract_signal2's arrays, if it
// has just run for this signal.
float decode_spread(SpreadScratch& sc, int wspr_type, const float *id, const float *qd, long np, float f0, int shift0,
                    float drift0, const unsigned char* channel_symbols, const SubScratch* sub)
{
    const float fs = wspr_type == 15 ? 375.0f/8 : 375.0f;
    const size_t nfft = wspr_type == 15 ? 2*900*375/8 : 2*120*375;
    if (!sc.fft || sc.fft->size() != nfft) {
        sc.fft.reset();
        sc.fft = std::make_unique<skimmer::Fft>(nfft, true);
    }
    std::complex<float>* g = sc.fft->data();
    std::fill(g, g + nfft, std::complex<float>());
    spread_gain(id, qd, std::min<long>(np, (long)nfft), f0, shift0, drift0, channel_symbols, sub, g);
    return dopspread(*sc.fft, fs);
}

unsigned int count_hard_errors( unsigned char *symbols, unsigned char *channel_symbols)
{
    int i,is;
    unsigned char cw[162];
    unsigned int nerrors;
    for (i=0; i<162; i++) {
        cw[i] = channel_symbols[i] >=2 ? 1:0;
    }
    deinterleave(cw);
    nerrors=0;
    for (i=0; i<162; i++) {
        is = symbols[i] > 127 ? 1:0;
        nerrors = nerrors + (is == cw[i] ? 0:1);
    }
    return nerrors;
}

struct cand { float freq; float snr; int shift; float drift; float sync; };

struct result { char date[7]; char time[5]; float sync; float snr;
    float dt; double freq; char message[23]; float drift;
    unsigned int cycles; int jitter; int blocksize; unsigned int metric;
    int nhardmin; int ipass; int decodetype;
    float spread;};  // not wsprd's: Decode::spread

// The message's fields for Decode.
void split_message(Decode& d) {
    std::vector<std::string> tok;
    std::string cur;
    for (char c : d.message) {
        if (c == ' ') {
            if (!cur.empty()) tok.push_back(cur);
            cur.clear();
        } else {
            cur += c;
        }
    }
    if (!cur.empty()) tok.push_back(cur);
    if (tok.size() == 3) {
        d.callsign = tok[0];
        d.grid = tok[1];
        d.power = std::atoi(tok[2].c_str());
    } else if (tok.size() == 2) {
        d.callsign = tok[0];
        d.power = std::atoi(tok[1].c_str());
    }
}

}  // namespace

// ------------------------------------------------------------------- Hashtable

Hashtable::Hashtable() : hashtab_(WSPRD_HASH_COUNT * WSPRD_CALLSIGN_SIZE), loctab_(WSPRD_HASH_COUNT * WSPRD_GRID4_SIZE) {}

void Hashtable::clear() {
    std::fill(hashtab_.begin(), hashtab_.end(), 0);
    std::fill(loctab_.begin(), loctab_.end(), 0);
}

bool Hashtable::load_line(const char* line) { return wsprd_load_hash_line(line, hashtab_.data(), loctab_.data()) != 0; }

void Hashtable::load(const std::string& text) {
    // fgets(line, 80, f): up to 79 bytes, ending after a newline; a NUL ends
    // what wsprd_load_hash_line sees.
    size_t at = 0;
    while (at < text.size()) {
        size_t n = 0;
        while (n < 79 && at + n < text.size()) {
            ++n;
            if (text[at + n - 1] == '\n') break;
        }
        std::string line = text.substr(at, n);
        load_line(line.c_str());
        at += n;
    }
}

std::string Hashtable::save() const {
    std::string out;
    char buf[64];
    for (int i=0; i<WSPRD_HASH_COUNT; i++) {
        const char* call = hashtab_.data() + i * WSPRD_CALLSIGN_SIZE;
        if( strncmp(call,"\0",1) != 0 ) {
            snprintf(buf, sizeof buf, "%5d %s %s\n", i, call, loctab_.data() + i * WSPRD_GRID4_SIZE);
            out += buf;
        }
    }
    return out;
}

std::string Hashtable::call(int ihash) const {
    if (ihash < 0 || ihash >= WSPRD_HASH_COUNT) return {};
    return std::string(hashtab_.data() + ihash * WSPRD_CALLSIGN_SIZE);
}

std::string Hashtable::grid(int ihash) const {
    if (ihash < 0 || ihash >= WSPRD_HASH_COUNT) return {};
    return std::string(loctab_.data() + ihash * WSPRD_GRID4_SIZE);
}

int Hashtable::hash(const std::string& call) { return (int)nhash(call.data(), call.size(), (uint32_t)146); }

// --------------------------------------------------------------------- Decoder

struct Decoder::Impl {
    std::vector<float> idat, qdat;  // maxpts each
    std::vector<float> ps;          // ps[512][nffts]
    std::vector<snode> stack;       // -J
    std::vector<int16_t> wav;       // readwavfile's buf2
    std::vector<float> wavf;        // readwavfile's buf2, from a float .wav
    SubScratch sub;                 // subtract_signal2's arrays
    SpreadScratch spread;           // decode_spread's FFT (not wsprd's)
    OsdTables osd_tables;           // osdwspr's boxit/fetchit tables
    size_t plan_bytes = 0;          // the last decode's FFTW buffers (memory())

    static constexpr int maxpts = 65536;

    Impl() : idat(maxpts), qdat(maxpts) {}

    void decode(const Options& opt, char* hashtab, char* loctab, unsigned int npoints, int wspr_type,
                double dialfreq, const char* date, const char* uttime, Plans& plans, Timings& tm, Result& res,
                const std::atomic<bool>* stop);
};

void Decoder::Impl::decode(const Options& opt, char* hashtab, char* loctab, unsigned int npoints, int wspr_type,
                           double dialfreq, const char* date, const char* uttime, Plans& plans, Timings& tm,
                           Result& res, const std::atomic<bool>* stop)
{
    int i,j,k;
    unsigned char symbols[162], decdata[11], channel_symbols[162], apmask[162], cw[162];
    signed char message[]={-9,13,-35,123,57,-39,64,0,0,0,0};
    char callsign[13], grid[5], call_loc_pow[23];
    int verbose=opt.verbose, quickmode=opt.quickmode, more_candidates=opt.more_candidates;
    int ipass, nblocksize;
    int nhardmin,ihash;
    int writec2=opt.writec2,maxdrift;
    int shift1, lagmin, lagmax, lagstep, ifmin, ifmax, not_decoded;
    unsigned int nbits=81, stacksize=200000;
    struct snode * stack=NULL;
    unsigned int cycles, maxnp, metric;
    float df=375.0/256.0/2;
    float fsymbs[162];
    float dt=1.0/375.0, dt_print;
    double freq_print;
    double dialfreq_error=opt.dialfreq_error;
    float fmin=opt.fmin, fmax=opt.fmax;
    float f1, fstep, sync1, drift1;
    float dmin;
    float psavg[512];
    float *idat = this->idat.data(), *qdat = this->qdat.data();
    clock_t t0;

    struct cand candidates[200];

    struct result decodes[50];

    memset(apmask,0,sizeof apmask);
    memset(cw,0,sizeof cw);
    memset(decdata,0,sizeof decdata);
    memset(channel_symbols,0,sizeof channel_symbols);
    float allfreqs[100];
    char allcalls[100][13];
    for (i=0; i<100; i++) allfreqs[i]=0.0;
    memset(allcalls,0,sizeof(char)*100*13);

    int uniques=0, noprint=0, ndecodes_pass=0;

    // Parameters used for performance-tuning:
    unsigned int maxcycles=opt.maxcycles;    //Decoder timeout limit
    float minsync1=0.10;                     //First sync limit
    float minsync2=0.12;                     //Second sync limit
    int iifac=8;                             //Step size in final DT peakup
    int symfac=50;                           //Soft-symbol normalizing factor
    int subtraction=opt.subtraction;
    int npasses=opt.npasses;
    int ndepth=opt.ndepth;                   //Depth for OSD

    float minrms=52.0 * (symfac/64.0);      //Final test for plausible decoding
    int delta=60;                            //Fano threshold step
    float bias=opt.bias;                    //Fano metric bias (used for both Fano and stack algorithms)

    fftwf_complex *fftin, *fftout;

    int mettab[2][256];

    if( opt.stackdecoder ) {
        if (this->stack.empty()) this->stack.resize(stacksize);
        stack=this->stack.data();
    }

    // setup metric table
    for(i=0; i<256; i++) {
        mettab[0][i]=round( (double)(10*(metric_tables[2][i]-bias)) );
        mettab[1][i]=round( (double)(10*(metric_tables[2][255-i]-bias)) );
    }

    TR("IQ npoints=%u type=%d dial=%a i=%016llx q=%016llx\n", npoints, wspr_type, dialfreq,
       fnv64(idat, sizeof(float)*npoints), fnv64(qdat, sizeof(float)*npoints));
    // Do windowed ffts over 2 symbols, stepped by half symbols
    int nffts=4*floor(npoints/512)-1;
    fftin=plans.fftin3;
    fftout=plans.fftout3;

    this->ps.assign((size_t)512*nffts, 0.0f);
    float *ps=this->ps.data();
    float w[512];
    for(i=0; i<512; i++) {
        w[i]=sin(0.006147931*i);
    }

    TurnSteps steps;  // the SIMD sums' (allocated on first use)

    // Decoder::set_stop()'s flag, looked at each pass, candidate and try
    // (not the reference's: unset, or never set, it changes nothing). Set,
    // the decode ends there (its caller drops what it has, as a killed wsprd).
    auto stopping = [&] {
        if (stop && stop->load(std::memory_order_relaxed)) res.stopped = true;
        return res.stopped;
    };

    //*************** main loop starts here *****************
    for (ipass=0; ipass<npasses; ipass++) {
        if (stopping()) goto stopped;
        if(ipass==1 && ndecodes_pass == 0 && npasses>2) ipass=2;
        if(ipass < 2) {
            nblocksize=1;
            maxdrift=4;
            minsync2=0.12;
        }
        else {
            nblocksize=4;  // try 3 blocksizes plus bitbybit normalization
            maxdrift=0;    // no drift for smaller frequency estimator variance
            minsync2=0.10;
        }
        if (opt.nodrift) maxdrift=0;  // -n (not wsprd's): the candidates' drift search only
        ndecodes_pass=0;   // still needed?
        TR("PASS %d nblocksize=%d maxdrift=%d minsync2=%a fmin=%a fmax=%a i=%016llx q=%016llx\n", ipass, nblocksize, maxdrift, minsync2,
           fmin, fmax, fnv64(idat, sizeof(float)*npoints), fnv64(qdat, sizeof(float)*npoints));

        for (i=0; i<nffts; i++) {
            for(j=0; j<512; j++ ) {
                k=i*128+j;
                fftin[j][0]=idat[k] * w[j];
                fftin[j][1]=qdat[k] * w[j];
            }
            fftwf_execute(plans.PLAN3);
            for (j=0; j<512; j++ ) {
                k=j+256;
                if( k>511 )
                    k=k-512;
                ps[j*nffts+i]=fftout[k][0]*fftout[k][0]+fftout[k][1]*fftout[k][1];
            }
        }

        TR("PS %d nffts=%d ps=%016llx\n", ipass, nffts, fnv64(ps, sizeof(float)*512*nffts));
        // Compute average spectrum
        for (i=0; i<512; i++) psavg[i]=0.0;
        for (i=0; i<nffts; i++) {
            for (j=0; j<512; j++) {
                psavg[j]=psavg[j]+ps[j*nffts+i];
            }
        }

        // Smooth with 7-point window and limit spectrum to +/-150 Hz
        int window[7]={1,1,1,1,1,1,1};
        float smspec[411];
        for (i=0; i<411; i++) {
            smspec[i]=0.0;
            for(j=-3; j<=3; j++) {
                k=256-205+i+j;
                smspec[i]=smspec[i]+window[j+3]*psavg[k];
            }
        }

        // Sort spectrum values, then pick off noise level as a percentile
        float tmpsort[411];
        for (j=0; j<411; j++) {
            tmpsort[j]=smspec[j];
        }
        qsort(tmpsort, 411, sizeof(float), floatcomp);

        // Noise level of spectrum is estimated as 123/411= 30'th percentile
        float noise_level = tmpsort[122];
        TR("NOISE %d psavg=%016llx smspec=%016llx noise=%a\n", ipass, fnv64(psavg, sizeof psavg), fnv64(smspec, sizeof smspec), noise_level);

        /* Renormalize spectrum so that (large) peaks represent an estimate of snr.
         * We know from experience that threshold snr is near -7dB in wspr bandwidth,
         * corresponding to -7-26.3=-33.3dB in 2500 Hz bandwidth.
         * The corresponding threshold is -42.3 dB in 2500 Hz bandwidth for WSPR-15. */

        float min_snr, snr_scaling_factor;
        min_snr = pow(10.0,-8.0/10.0); //this is min snr in wspr bw
        if( wspr_type == 2 ) {
            snr_scaling_factor=26.3;
        } else {
            snr_scaling_factor=35.3;
        }
        for (j=0; j<411; j++) {
            smspec[j]=smspec[j]/noise_level - 1.0;
            if( smspec[j] < min_snr) smspec[j]=0.1*min_snr;
            continue;
        }

        // Find all local maxima in smoothed spectrum.
        for (i=0; i<200; i++) {
            candidates[i].freq=0.0;
            candidates[i].snr=0.0;
            candidates[i].drift=0.0;
            candidates[i].shift=0;
            candidates[i].sync=0.0;
        }

        int npk=0;
        unsigned char candidate;
        for(j=1; j<410; j++) {
            candidate = (smspec[j]>smspec[j-1]) &&
            (smspec[j]>smspec[j+1]) &&
            (npk<200);
            if ( candidate ) {
                candidates[npk].freq = (j-205)*df;
                candidates[npk].snr  = 10*log10((double)smspec[j])-snr_scaling_factor;
                npk++;
            }
        }
        if( more_candidates ) {
            for(j=0; j<411; j=j+3) {
                candidate = (smspec[j]>min_snr) && (npk<200);
                if ( candidate ) {
                    candidates[npk].freq = (j-205)*df;
                    candidates[npk].snr  = 10*log10((double)smspec[j])-snr_scaling_factor;
                    npk++;
                }
            }
        }

        for (i=0; i<npk; i++)
            TR("CAND0 %d %d freq=%a snr=%a\n", ipass, i, candidates[i].freq, candidates[i].snr);
        // Compute corrected fmin, fmax, accounting for dial frequency error
        fmin += dialfreq_error;    // dialfreq_error is in units of Hz
        fmax += dialfreq_error;

        // Don't waste time on signals outside of the range [fmin,fmax].
        i=0;
        for( j=0; j<npk; j++) {
            if( candidates[j].freq >= fmin && candidates[j].freq <= fmax ) {
                candidates[i]=candidates[j];
                i++;
            }
        }
        npk=i;

        // bubble sort on snr
        int pass;
        struct cand tmp;
        for (pass = 1; pass <= npk - 1; pass++) {
            for (k = 0; k < npk - pass ; k++) {
                if (candidates[k].snr < candidates[k+1].snr) {
                    tmp = candidates[k];
                    candidates[k]=candidates[k+1];
                    candidates[k+1] = tmp;
                }
            }
        }

        for (i=0; i<npk; i++)
            TR("CAND1 %d %d freq=%a snr=%a\n", ipass, i, candidates[i].freq, candidates[i].snr);
        t0=cpu_clock();

        /* Make coarse estimates of shift (DT), freq, and drift

         * Look for time offsets up to +/- 8 symbols (about +/- 5.4 s) relative
         to nominal start time, which is 2 seconds into the file

         * Calculates shift relative to the beginning of the file

         * Negative shifts mean that signal started before start of file

         * The program prints DT = shift-2 s

         * Shifts that cause sync vector to fall off of either end of the data
         vector are accommodated by "partial decoding", such that missing
         symbols produce a soft-decision symbol value of 128

         * The frequency drift model is linear, deviation of +/- drift/2 over the
         span of 162 symbols, with deviation equal to 0 at the center of the
         signal vector.
         */

        int idrift,ifr,if0,k0;
        float smax;
#ifndef WSPRD_SIMD
        int kindex;
        float ss,pow,p0,p1,p2,p3;
#endif
        // Every float as the reference's, with less work: ps's square roots
        // are taken once, in place (the search took (float)sqrt((double)ps)
        // of each power it read, and nothing reads ps again before the next
        // pass makes it anew), and each (ifr, idrift)'s ifd of every symbol
        // once, by the reference's expression, rather than for every k0;
        // with AVX2 or NEON the 32 k0 are summed side by side (coarse_sync).
        if (npk > 0)
            for (i=0; i<512*nffts; i++) ps[i]=sqrt((double)ps[i]);
        int ifds[9][162];
        for(j=0; j<npk; j++) {                              //For each candidate...
            smax=-1e30;
            if0=candidates[j].freq/df+256;
            for (ifr=if0-2; ifr<=if0+2; ifr++) {                      //Freq search
                for (idrift=-maxdrift; idrift<=maxdrift; idrift++)
                    for (k=0; k<162; k++)
                        ifds[idrift+maxdrift][k]=ifr+((float)k-81.0)/81.0*( (float)idrift )/(2.0*df);
#ifdef WSPRD_SIMD
                float syncs[9][32];
                for (idrift=0; idrift<2*maxdrift+1; idrift++)
                    coarse_sync(ps, nffts, ifds[idrift], syncs[idrift]);
#endif
                for( k0=-10; k0<22; k0++) {                             //Time search
                    for (idrift=-maxdrift; idrift<=maxdrift; idrift++) {  //Drift search
#ifdef WSPRD_SIMD
                        sync1=syncs[idrift+maxdrift][k0+10];
#else
                        const int *ifd=ifds[idrift+maxdrift];
                        ss=0.0;
                        pow=0.0;
                        for (k=0; k<162; k++) {                             //Sum over symbols
                            kindex=k0+2*k;
                            if( kindex >= 0 && kindex < nffts ) {
                                p0=ps[(ifd[k]-3)*nffts+kindex];
                                p1=ps[(ifd[k]-1)*nffts+kindex];
                                p2=ps[(ifd[k]+1)*nffts+kindex];
                                p3=ps[(ifd[k]+3)*nffts+kindex];

                                ss=ss+(2*pr3[k]-1)*((p1+p3)-(p0+p2));
                                pow=pow+p0+p1+p2+p3;
                            }
                        }
                        sync1=ss/pow;
#endif
                        if( sync1 > smax ) {                  //Save coarse parameters
                            smax=sync1;
                            candidates[j].shift=128*(k0+1);
                            candidates[j].drift=idrift;
                            candidates[j].freq=(ifr-256)*df;
                            candidates[j].sync=sync1;
                        }
                    }
                }
            }
        }
        tm.tcandidates += (float)(cpu_clock()-t0)/CLOCKS_PER_SEC;
        for (i=0; i<npk; i++)
            TR("CAND2 %d %d freq=%a shift=%d drift=%a sync=%a\n", ipass, i, candidates[i].freq,
               candidates[i].shift, candidates[i].drift, candidates[i].sync);

        /*
         Refine the estimates of freq, shift using sync as a metric.
         Sync is calculated such that it is a float taking values in the range
         [0.0,1.0].

         Function sync_and_demodulate has three modes of operation
         mode is the last argument:

         0 = no frequency or drift search. find best time lag.
         1 = no time lag or drift search. find best frequency.
         2 = no frequency or time lag search. Calculate soft-decision
         symbols using passed frequency and shift.

         NB: best possibility for OpenMP may be here: several worker threads
         could each work on one candidate at a time.
         */
        for (j=0; j<npk; j++) {
            if (stopping()) goto stopped;

            f1=candidates[j].freq;
            drift1=candidates[j].drift;
            shift1=candidates[j].shift;
            sync1=candidates[j].sync;
            SyncMemo memo;
            memo.steps=&steps;

            // coarse-grid lag and freq search, then if sync>minsync1 continue
            fstep=0.0; ifmin=0; ifmax=0;
            lagmin=shift1-128;
            lagmax=shift1+128;
            lagstep=64;
            t0 = cpu_clock();
            sync_and_demodulate(idat, qdat, npoints, symbols, &f1, ifmin, ifmax, fstep, &shift1,
                                lagmin, lagmax, lagstep, &drift1, symfac, &sync1, 0, &memo);
            tm.tsync0 += (float)(cpu_clock()-t0)/CLOCKS_PER_SEC;
            TR("REF %d %d lag f=%a shift=%d drift=%a sync=%a\n", ipass, j, f1, shift1, drift1, sync1);

            fstep=0.25; ifmin=-2; ifmax=2;
            t0 = cpu_clock();
            sync_and_demodulate(idat, qdat, npoints, symbols, &f1, ifmin, ifmax, fstep, &shift1,
                                lagmin, lagmax, lagstep, &drift1, symfac, &sync1, 1, &memo);
            TR("REF %d %d freq f=%a shift=%d drift=%a sync=%a\n", ipass, j, f1, shift1, drift1, sync1);

            if(ipass < 2) {
                // refine drift estimate
                fstep=0.0; ifmin=0; ifmax=0;
                float driftp,driftm,syncp,syncm;
                driftp=drift1+0.5;
                sync_and_demodulate(idat, qdat, npoints, symbols, &f1, ifmin, ifmax, fstep, &shift1,
                                    lagmin, lagmax, lagstep, &driftp, symfac, &syncp, 1, &memo);

                driftm=drift1-0.5;
                sync_and_demodulate(idat, qdat, npoints, symbols, &f1, ifmin, ifmax, fstep, &shift1,
                                    lagmin, lagmax, lagstep, &driftm, symfac, &syncm, 1, &memo);

                TR("REF %d %d drift p=%a syncp=%a m=%a syncm=%a\n", ipass, j, driftp, syncp, driftm, syncm);
                if(syncp>sync1) {
                    drift1=driftp;
                    sync1=syncp;
                } else if (syncm>sync1) {
                    drift1=driftm;
                    sync1=syncm;
                }
            }
            tm.tsync1 += (float)(cpu_clock()-t0)/CLOCKS_PER_SEC;

            // fine-grid lag and freq search
            if( sync1 > minsync1 ) {

                lagmin=shift1-32; lagmax=shift1+32; lagstep=16;
                t0 = cpu_clock();
                sync_and_demodulate(idat, qdat, npoints, symbols, &f1, ifmin, ifmax, fstep, &shift1,
                                    lagmin, lagmax, lagstep, &drift1, symfac, &sync1, 0, &memo);
                tm.tsync0 += (float)(cpu_clock()-t0)/CLOCKS_PER_SEC;
                TR("REF %d %d finelag f=%a shift=%d drift=%a sync=%a\n", ipass, j, f1, shift1, drift1, sync1);

                // fine search over frequency
                fstep=0.05; ifmin=-2; ifmax=2;
                t0 = cpu_clock();
                sync_and_demodulate(idat, qdat, npoints, symbols, &f1, ifmin, ifmax, fstep, &shift1,
                                    lagmin, lagmax, lagstep, &drift1, symfac, &sync1, 1, &memo);
                tm.tsync1 += (float)(cpu_clock()-t0)/CLOCKS_PER_SEC;
                TR("REF %d %d finefreq f=%a shift=%d drift=%a sync=%a\n", ipass, j, f1, shift1, drift1, sync1);

                candidates[j].freq=f1;
                candidates[j].shift=shift1;
                candidates[j].drift=drift1;
                candidates[j].sync=sync1;
            }
        }

        for (i=0; i<npk; i++)
            TR("CAND3 %d %d freq=%a shift=%d drift=%a sync=%a\n", ipass, i, candidates[i].freq,
               candidates[i].shift, candidates[i].drift, candidates[i].sync);
        int nwat=0;
        int idupe;
        for ( j=0; j<npk; j++) {
            idupe=0;
            for (k=0;k<nwat;k++) {
                if( fabsf(candidates[j].freq - candidates[k].freq) < 0.05  &&
                   abs(candidates[j].shift - candidates[k].shift) < 16 ) {
                    idupe=1;
                    break;
                }
            }
            if( idupe == 1 ) {
                if(candidates[j].sync > candidates[k].sync) candidates[k]=candidates[j];
            } else if ( candidates[j].sync > minsync2 ) {
                candidates[nwat]=candidates[j];
                nwat++;
            }
        }

        for (i=0; i<nwat; i++)
            TR("CAND4 %d %d freq=%a shift=%d drift=%a sync=%a snr=%a\n", ipass, i, candidates[i].freq,
               candidates[i].shift, candidates[i].drift, candidates[i].sync, candidates[i].snr);
        int idt, ii, jittered_shift;
        float y,sq,rms;
        int ib, blocksize, bitmetric;
        int n1,n2,n3,nadd,nu,ntype;
        int osd_decode;
        NcMemo ncmemo;
        ncmemo.steps=&steps;
        for (j=0; j<nwat; j++) {
            if (stopping()) goto stopped;
            memset(symbols,0,sizeof(char)*nbits*2);
            memset(callsign,0,sizeof(char)*13);
            memset(grid,0,sizeof(char)*5);
            memset(call_loc_pow,0,sizeof(char)*23);
            f1=candidates[j].freq;
            shift1=candidates[j].shift;
            drift1=candidates[j].drift;
            not_decoded=1;
            osd_decode=0;
            ncmemo.n=0;
            bool subtracted=false;  // subtract_signal2 ran for this decode: its ci/cq the spread's gain

            ib=1;
            while( ib <= nblocksize && not_decoded ) {
                if (ib < 4) { blocksize=ib; bitmetric=0; }
                if (ib == 4) { blocksize=1; bitmetric=1; }

                idt=0; ii=0;
                while ( not_decoded && idt<=(128/iifac)) {
                    if (stopping()) goto stopped;
                    ii=(idt+1)/2;
                    if( idt%2 == 1 ) ii=-ii;
                    ii=iifac*ii;
                    jittered_shift=shift1+ii;
                    nhardmin=0; dmin=0.0;

                    // Get soft-decision symbols
                    t0 = cpu_clock();
#ifdef WSPRD_SIMD
                    if (ib == 1 && idt%4 == 1 && drift1 != 0) {
                        // this shift and the next three: one table per symbol for all
                        int lags[4], nl=0;
                        for (int d=idt; d<idt+4 && d<=(128/iifac); d++) {
                            int jj=(d+1)/2;
                            if( d%2 == 1 ) jj=-jj;
                            lags[nl++]=shift1+iifac*jj;
                        }
                        noncoherent_lags(idat, qdat, npoints, f1, lags, nl, drift1, &ncmemo);
                    }
#endif
                    noncoherent_sequence_detection(idat, qdat, npoints, symbols, &f1,
                                                   &jittered_shift, &drift1, symfac, &blocksize, &bitmetric, &ncmemo);
                    tm.tsync2 += (float)(cpu_clock()-t0)/CLOCKS_PER_SEC;

                    sq=0.0;
                    for(i=0; i<162; i++) {
                        y=(float)symbols[i] - 128.0;
                        sq += y*y;
                    }
                    rms=sqrt(sq/162.0);
                    TR("TRY %d %d ib=%d idt=%d shift=%d sym=%016llx rms=%a%s\n", ipass, j, ib, idt, jittered_shift,
                       fnv64(symbols, 162), rms, rms > minrms ? "" : " low");

                    if(rms > minrms) {
                        deinterleave(symbols);
                        t0 = cpu_clock();

                        if ( stack ) {
                            not_decoded = jelinek(&metric, &cycles, decdata, symbols, nbits,
                                                  stacksize, stack, mettab,maxcycles);
                        } else {
                            not_decoded = fano(&metric,&cycles,&maxnp,decdata,symbols,nbits,
                                               mettab,delta,maxcycles);
                        }

                        tm.tfano += (float)(cpu_clock()-t0)/CLOCKS_PER_SEC;
                        TR("FANO %d %d %s fail=%d metric=%u cycles=%u maxnp=%u data=%016llx\n", ipass, j, stack ? "stack" : "fano",
                           not_decoded, metric, cycles, stack ? 0u : maxnp, fnv64(decdata, 11));

                        if( (ndepth >= 0) && not_decoded ) {
                            for(i=0; i<162; i++) {
                                fsymbs[i]=symbols[i]-128.0;
                            }
                            t0 = cpu_clock();
                            osdwspr(fsymbs,apmask,ndepth,cw,&nhardmin,&dmin,&osd_tables);
                            tm.tosd += (float)(cpu_clock()-t0)/CLOCKS_PER_SEC;
                            TR("OSD %d %d nhardmin=%d dmin=%a cw=%016llx\n", ipass, j, nhardmin, dmin, fnv64(cw, 162));

                            for(i=0; i<162; i++) {
                                symbols[i]=255*cw[i];
                            }
                            fano(&metric,&cycles,&maxnp,decdata,symbols,nbits,
                                 mettab,delta,maxcycles);
                            for(i=0; i<11; i++) {
                                if( decdata[i]>127 ) {
                                    message[i]=decdata[i]-256;
                                } else {
                                    message[i]=decdata[i];
                                }
                            }
                            unpack50(message,&n1,&n2);
                            TR("OSDFANO %d %d metric=%u cycles=%u n1=%d n2=%d\n", ipass, j, metric, cycles, n1, n2);
                            if( !unpackcall(n1,callsign) ) break;
                            callsign[12]=0;
                            if( !unpackgrid(n2, grid) ) break;
                            grid[4]=0;
                            ntype = (n2&127) - 64;
                            int itype;
                            if( (ntype >= 0) && (ntype <= 62) ) {
                                nu = ntype%10;
                                itype=1;
                                if( !(nu == 0 || nu == 3 || nu == 7) ) {
                                    nadd=nu;
                                    if( nu > 3 ) nadd=nu-3;
                                    if( nu > 7 ) nadd=nu-7;
                                    n3=n2/128+32768*(nadd-1);
                                    if( !unpackpfx(n3,callsign) ) {
                                        break;
                                    }
                                    itype=2;
                                }
                                ihash=nhash(callsign,strlen(callsign),(uint32_t)146);
                                TR("OSDHASH %d %d call=%s grid=%s ihash=%d known=%d\n", ipass, j, callsign, grid, ihash,
                                   strncmp(hashtab+ihash*13,callsign,13)==0);
                                if(strncmp(hashtab+ihash*13,callsign,13)==0) {
                                    if( (itype==1 && strncmp(loctab+ihash*5,grid,5)==0) ||
                                        (itype==2) ) {
                                       not_decoded=0;
                                       osd_decode =1;
                                    }
                                }
                            }
                        }

                    }
                    idt++;
                    if( quickmode ) break;
                }
                ib++;
            }

            if( !not_decoded ) {
                ndecodes_pass++;

                for(i=0; i<11; i++) {

                    if( decdata[i]>127 ) {
                        message[i]=decdata[i]-256;
                    } else {
                        message[i]=decdata[i];
                    }

                }

                // Unpack the decoded message, update the hashtable, apply
                // sanity checks on grid and power, and return
                // call_loc_pow string and also callsign (for de-duping).
                noprint=unpk_(message,hashtab,loctab,call_loc_pow,callsign);
                TR("DEC %d %d f=%a shift=%d drift=%a ib=%d ii=%d osd=%d noprint=%d msg=[%s] call=[%s]\n", ipass, j, f1, shift1, drift1,
                   ib, ii, osd_decode, noprint, call_loc_pow, callsign);
                if( subtraction && !noprint ) {
                    if( get_wspr_channel_symbols(call_loc_pow, hashtab, loctab, channel_symbols) ) {
                        subtract_signal2(idat, qdat, npoints, f1, shift1, drift1, channel_symbols, sub);
                        subtracted=true;
                        if(!osd_decode) nhardmin=count_hard_errors(symbols,channel_symbols);
                        TR("SUB %d %d chan=%016llx nhardmin=%d i=%016llx q=%016llx\n", ipass, j, fnv64(channel_symbols, 162), nhardmin,
                           fnv64(idat, sizeof(float)*npoints), fnv64(qdat, sizeof(float)*npoints));
                    } else {
                        TR("NOCHAN %d %d\n", ipass, j);
                        break;
                    }
                }

                // Remove dupes (same callsign and freq within 4 Hz)
                int dupe=0;
                for (i=0; i<uniques; i++) {
                    if(!strcmp(callsign,allcalls[i]) &&
                       (fabs((double)(f1-allfreqs[i])) <4.0)) dupe=1;
                }
                if( (verbose || !dupe) && !noprint &&
                    uniques < (int)(sizeof decodes / sizeof decodes[0]) ) {
                    strcpy(allcalls[uniques],callsign);
                    allfreqs[uniques]=f1;
                    uniques++;

                    // Add an extra space at the end of each line so that wspr-x doesn't
                    // truncate the power (TNX to DL8FCL!)

                    if( wspr_type == 15 ) {
                        freq_print=dialfreq+(1500+112.5+f1/8.0)/1e6;
                        dt_print=shift1*8*dt-1.0;
                    } else {
                        freq_print=dialfreq+(1500+f1)/1e6;
                        dt_print=shift1*dt-1.0;
                    }

                    strcpy(decodes[uniques-1].date,date);
                    strcpy(decodes[uniques-1].time,uttime);
                    decodes[uniques-1].sync=candidates[j].sync;
                    decodes[uniques-1].snr=candidates[j].snr;
                    decodes[uniques-1].dt=dt_print;
                    decodes[uniques-1].freq=freq_print;
                    strcpy(decodes[uniques-1].message,call_loc_pow);
                    decodes[uniques-1].drift=drift1;
                    decodes[uniques-1].cycles=cycles;
                    decodes[uniques-1].jitter=ii;
                    decodes[uniques-1].blocksize=blocksize+3*bitmetric;
                    decodes[uniques-1].metric=metric;
                    decodes[uniques-1].nhardmin=nhardmin;
                    decodes[uniques-1].ipass=ipass;
                    decodes[uniques-1].decodetype=osd_decode;
                    // Not wsprd's: its spread, from the samples as they were
                    // before it was subtracted (those before it gone).
                    decodes[uniques-1].spread=-1;
                    if (opt.spread) {
                        unsigned char own[162];  // without subtraction: wsprd's channel_symbols left as they were
                        if (subtracted)
                            decodes[uniques-1].spread=decode_spread(this->spread, wspr_type, idat, qdat, npoints, f1, shift1,
                                                                    drift1, channel_symbols, &sub);
                        else if (get_wspr_channel_symbols(call_loc_pow, hashtab, loctab, own))
                            decodes[uniques-1].spread=decode_spread(this->spread, wspr_type, idat, qdat, npoints, f1, shift1,
                                                                    drift1, own, nullptr);
                    }
                }
            }
        }

        TR("ENDPASS %d ndecodes=%d uniques=%d\n", ipass, ndecodes_pass, uniques);
        if( ipass == 0 && writec2 ) {
            double carrierfreq=dialfreq;
            res.has_c2 = true;
            res.c2_idat.assign(idat, idat + 45000);
            res.c2_qdat.assign(qdat, qdat + 45000);
            res.c2_freq = carrierfreq;
        }
    }
stopped:

    // sort the result in order of increasing frequency
    struct result temp;
    for (j = 1; j <= uniques - 1; j++) {
        for (k = 0; k < uniques - j ; k++) {
            if (decodes[k].freq > decodes[k+1].freq) {
                temp = decodes[k];
                decodes[k]=decodes[k+1];;
                decodes[k+1] = temp;
            }
        }
    }

    for (i=0; i<uniques; i++) {
        Decode d;
        d.date = decodes[i].date;
        d.time = decodes[i].time;
        d.sync = decodes[i].sync;
        d.snr = decodes[i].snr;
        d.dt = decodes[i].dt;
        d.freq = decodes[i].freq;
        d.message = decodes[i].message;
        d.drift = decodes[i].drift;
        d.cycles = decodes[i].cycles;
        d.jitter = decodes[i].jitter;
        d.blocksize = decodes[i].blocksize;
        d.metric = decodes[i].metric;
        d.nhardmin = decodes[i].nhardmin;
        d.ipass = decodes[i].ipass;
        d.decodetype = decodes[i].decodetype;
        d.spread = decodes[i].spread;
        split_message(d);
        res.decodes.push_back(std::move(d));
    }
}

Decoder::Decoder(const Options& opt) : opt_(opt), impl_(std::make_unique<Impl>()) {}

// What a crashed wsprd leaves in its data directory: hashtable.txt and
// wspr_wisdom.dat as they were (both are written after the decode), and
// wspr_timer.out empty (opened "w" before it), so the next run's timers
// start from 0.
struct Decoder::Kept {
    Hashtable hashtable;
    std::string wisdom;
    explicit Kept(const Decoder& d) : hashtable(d.hashtable_), wisdom(d.wisdom_) {}
    void restore(Decoder& d) const {
        d.hashtable_ = hashtable;
        d.wisdom_ = wisdom;
        d.timings_ = Timings();
    }
    // A stopped decode as a killed wsprd's: it writes ALL_WSPR.TXT and
    // hashtable.txt only after its last pass, so nothing decoded and
    // nothing learned (has_c2 as the crash's).
    Result stopped(Decoder& d, Result& res) const {
        res.stopped = true;
        res.decodes.clear();
        restore(d);
        return std::move(res);
    }
};

bool Decoder::stopping(Result& res) const {
    if (stop_ && stop_->load(std::memory_order_relaxed)) res.stopped = true;
    return res.stopped;
}

Decoder::~Decoder() = default;

void Decoder::release() {
    Impl& m = *impl_;
    // Each is made anew by the next decode that needs it (assign(), or
    // made when empty), so nothing it decodes changes.
    std::vector<float>().swap(m.ps);
    std::vector<snode>().swap(m.stack);
    std::vector<int16_t>().swap(m.wav);
    std::vector<float>().swap(m.wavf);
    for (auto* v : {&m.sub.refi, &m.sub.refq, &m.sub.ci, &m.sub.cq, &m.sub.cfi, &m.sub.cfq}) std::vector<float>().swap(*v);
    m.osd_tables.reset();
    m.spread.fft.reset();
    m.plan_bytes = 0;
}

size_t Decoder::memory() const {
    const Impl& m = *impl_;
    size_t n = (m.idat.capacity() + m.qdat.capacity() + m.ps.capacity()) * sizeof(float) +
               m.stack.capacity() * sizeof(snode) + m.wav.capacity() * sizeof(int16_t) +
               m.wavf.capacity() * sizeof(float) + m.plan_bytes;
    for (const auto* v : {&m.sub.refi, &m.sub.refq, &m.sub.ci, &m.sub.cq, &m.sub.cfi, &m.sub.cfq})
        n += v->capacity() * sizeof(float);
    if (m.spread.fft) n += m.spread.fft->size() * sizeof(std::complex<float>);
    // The hash table, osdwspr's tables (some 330 KB once used), and a
    // decode's own memos (made each decode).
    n += 32768 * (13 + 5) + (m.osd_tables ? 330u << 10 : 0) + sizeof(TurnSteps::Slot) * TurnSteps::SLOTS +
         sizeof(NcMemo::Sums) * NcMemo::MAX;
    return n;
}

size_t Decoder::memory_needed() const {
    const Impl& m = *impl_;
    // decode_wav's (decode_c2's are smaller): readwavfile's samples and
    // FFTW buffers, the spectra (at most idat's 65536 points' worth), and
    // what -J, subtraction and OSD make when used.
    const size_t nfft2 = 46080, nfft1 = opt_.wspr_type == 15 ? nfft2 * 8 * 32 : nfft2 * 32;
    const size_t npoints = opt_.wspr_type == 15 ? 8 * 114 * 12000 : 114 * 12000;
    size_t n = (m.idat.capacity() + m.qdat.capacity()) * sizeof(float) + npoints * sizeof(int16_t) +
               sizeof(float) * nfft1 + sizeof(fftwf_complex) * (nfft1 / 2 + 1 + 2 * nfft2 + 2 * 512) +
               sizeof(float) * 512 * (4 * (Impl::maxpts / 512) - 1);
    if (opt_.stackdecoder) n += 200000 * sizeof(snode);
    if (opt_.subtraction) n += 6 * 45000 * sizeof(float);
    if (opt_.spread) n += 2 * 120 * 375 * sizeof(std::complex<float>);  // decode_spread's FFT (WSPR-15's is smaller)
    n += 32768 * (13 + 5) + (opt_.ndepth >= 0 ? 330u << 10 : 0) + sizeof(TurnSteps::Slot) * TurnSteps::SLOTS +
         sizeof(NcMemo::Sums) * NcMemo::MAX;
    return std::max(n, memory());
}

// Between runs wsprd writes its timers "%7.2f" to wspr_timer.out and reads
// them back with fscanf's %f.
void Decoder::carry_timings() {
    if (decoded_) {
        for (float* t : {&timings_.treadwav, &timings_.tcandidates, &timings_.tsync0, &timings_.tsync1,
                         &timings_.tsync2, &timings_.tfano, &timings_.tosd, &timings_.ttotal}) {
            char buf[64];
            snprintf(buf, sizeof buf, "%7.2f", *t);
            *t = strtof(buf, nullptr);
        }
    }
    decoded_ = true;
}

Result Decoder::decode_wav(const int16_t* samples, size_t n, const std::string& date, const std::string& time) {
    return decode_samples(samples, n, date, time);
}

Result Decoder::decode_wav(const float* samples, size_t n, const std::string& date, const std::string& time) {
    return decode_samples(samples, n, date, time);
}

namespace {
// readwavfile's buf2 for the samples' type: Impl's wav (16-bit) or wavf.
template <class Sample> struct WavBuffer;
template <> struct WavBuffer<int16_t> {
    static std::vector<int16_t>& of(std::vector<int16_t>& wav, std::vector<float>&) { return wav; }
};
template <> struct WavBuffer<float> {
    static std::vector<float>& of(std::vector<int16_t>&, std::vector<float>& wavf) { return wavf; }
};
}  // namespace

template <class Sample>
Result Decoder::decode_samples(const Sample* samples, size_t n, const std::string& date, const std::string& time) {
    Result res;
    const int wspr_type = opt_.wspr_type;
    carry_timings();
    if (n == 0 || (wspr_type != 2 && wspr_type != 15)) {
        timings_ = Timings();
        res.status = 1;
        return res;
    }
    // readwavfile's sizes: 114 s at 12 kHz (912 s for WSPR-15) into nfft1.
    const int nfft2 = 46080;
    const int nfft1 = wspr_type == 15 ? nfft2 * 8 * 32 : nfft2 * 32;
    const size_t npoints = wspr_type == 15 ? 8 * 114 * 12000 : 114 * 12000;
    const Kept kept(*this);
    // Stopped before the passes too: WSPR-15's planning and its 11.8M-point
    // FFT take a second.
    if (stopping(res)) return kept.stopped(*this, res);
    try {
        Plans plans;
        plans.make(nfft1, nfft2, wisdom_);
        if (stopping(res)) return kept.stopped(*this, res);
        clock_t t0 = cpu_clock();
        Impl& m = *impl_;
        m.plan_bytes = sizeof(float) * nfft1 + sizeof(fftwf_complex) * (nfft1 / 2 + 1 + 2 * nfft2 + 2 * 512);
        std::vector<Sample>& wav = WavBuffer<Sample>::of(m.wav, m.wavf);
        wav.assign(npoints, 0);
        std::memcpy(wav.data(), samples, sizeof(Sample) * (n < npoints ? n : npoints));
        std::fill(m.idat.begin(), m.idat.end(), 0.0f);
        std::fill(m.qdat.begin(), m.qdat.end(), 0.0f);
        const unsigned int np = readwavfile(wav.data(), npoints, wspr_type, m.idat.data(), m.qdat.data(), plans);
        timings_.treadwav += (float)(cpu_clock()-t0+plans.plan12)/CLOCKS_PER_SEC;
        const double dialfreq = opt_.dialfreq - (opt_.dialfreq_error*1.0e-06);
        Hashtable scratch_table;
        Hashtable& table = opt_.usehashtable ? hashtable_ : scratch_table;
        m.decode(opt_, table.hashtab(), table.loctab(), np, wspr_type, dialfreq, date.substr(0, 6).c_str(),
                 time.substr(0, 4).c_str(), plans, timings_, res, stop_);
    } catch (const detail::Crash&) {
        res.decodes.clear();
        res.status = 139;
        kept.restore(*this);
    } catch (...) {
        kept.restore(*this);  // as a wsprd that died of it leaves its directory
        throw;
    }
    if (res.stopped) return kept.stopped(*this, res);
    return res;
}

Result Decoder::decode_c2(const C2& c2, const std::string& date, const std::string& time) {
    Result res;
    carry_timings();
    if (c2.idat.size() < 45000 || c2.qdat.size() < 45000) {
        timings_ = Timings();
        res.status = 1;
        return res;
    }
    const Kept kept(*this);
    if (stopping(res)) return kept.stopped(*this, res);
    try {
        Plans plans;
        plans.make(0, 0, wisdom_);
        Impl& m = *impl_;
        m.plan_bytes = sizeof(fftwf_complex) * 2 * 512;
        std::fill(m.idat.begin(), m.idat.end(), 0.0f);
        std::fill(m.qdat.begin(), m.qdat.end(), 0.0f);
        std::copy(c2.idat.begin(), c2.idat.begin() + 45000, m.idat.begin());
        std::copy(c2.qdat.begin(), c2.qdat.begin() + 45000, m.qdat.begin());
        double dialfreq = c2.freq;
        dialfreq -= (opt_.dialfreq_error*1.0e-06);
        Hashtable scratch_table;
        Hashtable& table = opt_.usehashtable ? hashtable_ : scratch_table;
        m.decode(opt_, table.hashtab(), table.loctab(), 45000, c2.wspr_type, dialfreq, date.substr(0, 6).c_str(),
                 time.substr(0, 4).c_str(), plans, timings_, res, stop_);
    } catch (const detail::Crash&) {
        res.decodes.clear();
        res.status = 139;
        kept.restore(*this);
    } catch (...) {
        kept.restore(*this);
        throw;
    }
    if (res.stopped) return kept.stopped(*this, res);
    return res;
}

}  // namespace wsprd
