/*
 Synthetic WSPR inputs for scripts/wsprd-oracle.sh: several signals per
 file, each with its own frequency, start time, drift and SNR, in WSPR-2 or
 WSPR-15, as a .c2 file (wsprd's 375 Hz complex format) or a 12 kHz .wav.
 Built from WSJT-X's lib/wsprd/wsprsim.c, linked against its
 wsprsim_utils.c (get_wspr_channel_symbols) and what that needs.

 Copyright 2015, Steven Franke, K9AN
 License: GNU GPL v3

 The noise comes from a seeded xorshift64* generator, not rand(), so a seed
 gives the same file on any C library. The signal model is wsprsim's (phase
 continuous, tones 375/256 Hz apart at 375/256 baud) plus wsprd's drift
 model: the frequency of symbol i is f + drift/2*(i-81)/81.

   wsprd-sim [-S seed] [-m] [-F dial_MHz] [-N] (-o file.c2 | -W file.wav)
             "MESSAGE@f,t0,drift,snr" ...

 f is Hz from the centre (1500 Hz audio, 1612.5 Hz for WSPR-15 with -m),
 t0 the start in seconds into the file (1 s is DT 0), drift Hz over the
 transmission (in a WSPR-15 .c2, wsprd's own units: 8x the audio's), snr dB
 in 2500 Hz. -N: no noise. The .c2 header's name is
 the last 14 characters of the file name, as wsprsim writes "YYMMDD_HHMM.c2".
 */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "wsprsim_utils.h"
#include "wsprd_utils.h"
#include "fano.h"

int printdata=0;

static uint64_t rng_state = 1;

static double uniform01(void)
{
    rng_state ^= rng_state >> 12;
    rng_state ^= rng_state << 25;
    rng_state ^= rng_state >> 27;
    uint64_t x = rng_state * UINT64_C(2685821657736338717);
    return (double)(x >> 11) / 9007199254740992.0;  // [0,1)
}

// wsprsim's gaussrand, on uniform01()
static double gaussrand(void)
{
    static double V1, V2, S;
    static int phase = 0;
    double X;

    if(phase == 0) {
        do {
            double U1 = uniform01();
            double U2 = uniform01();

            V1 = 2 * U1 - 1;
            V2 = 2 * U2 - 1;
            S = V1 * V1 + V2 * V2;
        } while(S >= 1 || S == 0);

        X = V1 * sqrt(-2 * log(S) / S);
    } else
        X = V2 * sqrt(-2 * log(S) / S);

    phase = 1 - phase;

    return X;
}

struct sig { char message[23]; double f, t0, drift, snr; unsigned char symbols[162]; };

// Complex baseband at 375 Hz (decim 1) or real audio at 12 kHz: nsps
// samples per symbol, tone spacing fs/nsps; f, drift are already in Hz of
// that signal. Samples outside [0,n) are skipped.
static void add_complex(double *isig, double *qsig, long n, double fs, int nsps,
                        double f0, double t0, double drift, double amp,
                        const unsigned char *symbols)
{
    int i, j;
    long ii, idelay = (long)floor(t0*fs + 0.5);
    double phi=0.0, twopidt=8.0*atan(1.0)/fs, df=fs/nsps;
    for (i=0; i<162; i++) {
        double fp = f0 + (drift/2.0)*((double)i-81.0)/81.0;
        double dphi=twopidt*(fp + ((double)symbols[i]-1.5)*df);
        for (j=0; j<nsps; j++) {
            ii=idelay+(long)nsps*i+j;
            if (ii >= 0 && ii < n) {
                isig[ii]+=amp*cos(phi);
                if (qsig) qsig[ii]+=amp*sin(phi);
            }
            phi=phi+dphi;
        }
    }
}

static int write_c2(const char *path, int trmin, double freq, const double *isig, const double *qsig)
{
    static float buffer[2*45000];
    char name[14];
    size_t len = strlen(path);
    const char *base = len >= 14 ? path + len - 14 : path;
    memset(name, 0, sizeof name);
    memcpy(name, base, len >= 14 ? 14 : len);
    FILE *fp = fopen(path, "wb");
    if (!fp) { fprintf(stderr, "wsprd-sim: cannot write %s\n", path); return 0; }
    fwrite(name, 1, 14, fp);
    fwrite(&trmin, sizeof(int), 1, fp);
    fwrite(&freq, sizeof(double), 1, fp);
    for (int i=0; i<45000; i++) {
        buffer[2*i]=isig[i];
        buffer[2*i+1]=-qsig[i];
    }
    size_t n = fwrite(buffer, sizeof(float), 2*45000, fp);
    fclose(fp);
    return n == 2*45000;
}

static void put16(FILE *fp, unsigned v) { fputc(v & 255, fp); fputc((v >> 8) & 255, fp); }
static void put32(FILE *fp, unsigned v) { put16(fp, v & 0xffff); put16(fp, v >> 16); }

static int write_wav(const char *path, const double *x, long n)
{
    FILE *fp = fopen(path, "wb");
    if (!fp) { fprintf(stderr, "wsprd-sim: cannot write %s\n", path); return 0; }
    fwrite("RIFF", 1, 4, fp); put32(fp, 36 + 2*n);
    fwrite("WAVEfmt ", 1, 8, fp); put32(fp, 16); put16(fp, 1); put16(fp, 1);
    put32(fp, 12000); put32(fp, 24000); put16(fp, 2); put16(fp, 16);
    fwrite("data", 1, 4, fp); put32(fp, 2*n);
    for (long i=0; i<n; i++) {
        double v = floor(x[i] + 0.5);
        if (v > 32767) v = 32767;
        if (v < -32768) v = -32768;
        put16(fp, (unsigned)(int)v & 0xffff);
    }
    fclose(fp);
    return 1;
}

static void usage(void)
{
    fprintf(stderr, "usage: wsprd-sim [-S seed] [-m] [-F dial_MHz] [-N] (-o file.c2 | -W file.wav)\n"
                    "                 \"MESSAGE@f,t0,drift,snr\" ...\n");
}

int main(int argc, char *argv[])
{
    uint64_t seed = 1;
    int wspr15 = 0, noise = 1, i;
    double dial = 10.1387;
    const char *c2 = NULL, *wav = NULL;
    static struct sig sigs[64];
    int nsig = 0;
    static char hashtab[32768*13], loctab[32768*5];

    for (i=1; i<argc; i++) {
        if (!strcmp(argv[i], "-S") && i+1 < argc) seed = strtoull(argv[++i], NULL, 10);
        else if (!strcmp(argv[i], "-m")) wspr15 = 1;
        else if (!strcmp(argv[i], "-N")) noise = 0;
        else if (!strcmp(argv[i], "-F") && i+1 < argc) dial = strtod(argv[++i], NULL);
        else if (!strcmp(argv[i], "-o") && i+1 < argc) c2 = argv[++i];
        else if (!strcmp(argv[i], "-W") && i+1 < argc) wav = argv[++i];
        else if (argv[i][0] == '-') { usage(); return 2; }
        else {
            const char *at = strrchr(argv[i], '@');
            if (!at || nsig == 64 || at - argv[i] > 22) { usage(); return 2; }
            struct sig *s = &sigs[nsig++];
            memset(s->message, 0, sizeof s->message);
            memcpy(s->message, argv[i], at - argv[i]);
            if (sscanf(at+1, "%lf,%lf,%lf,%lf", &s->f, &s->t0, &s->drift, &s->snr) != 4) { usage(); return 2; }
            char msg[23];
            memcpy(msg, s->message, sizeof msg);
            if (!get_wspr_channel_symbols(msg, hashtab, loctab, s->symbols)) {
                fprintf(stderr, "wsprd-sim: cannot encode '%s'\n", s->message);
                return 2;
            }
        }
    }
    if ((!c2) == (!wav)) { usage(); return 2; }
    rng_state = seed ? seed : 1;

    if (c2) {
        // wsprsim's scaling: unit-variance noise in I and Q at 375 Hz; the
        // SNR in 375 Hz is 8.2 dB above that in 2500 Hz. For WSPR-15 the
        // file is the 8x decimated signal, so 1/8 the bandwidth again.
        static double isig[45000], qsig[45000];
        if (noise) for (i=0; i<45000; i++) { isig[i]=gaussrand(); qsig[i]=gaussrand(); }
        for (int k=0; k<nsig; k++) {
            struct sig *s = &sigs[k];
            double snr = s->snr + 8.2 + (wspr15 ? 9.03 : 0.0);
            double amp = noise ? pow(10, snr/20.0)*pow(2, 0.5) : 1.0;
            // A WSPR-15 .c2 is the audio decimated 8x further, so its
            // samples look like WSPR-2's: wsprd shows f/8 and t*8 for it,
            // hence 8f and t0/8 here (drift stays in wsprd's units).
            add_complex(isig, qsig, 45000, 375.0, 256, wspr15 ? 8*s->f : s->f,
                        wspr15 ? s->t0/8 : s->t0, s->drift, amp, s->symbols);
        }
        return write_c2(c2, wspr15 ? 15 : 2, dial, isig, qsig) ? 0 : 1;
    }

    // 12 kHz audio: noise sigma 300 (in int16 units) over 6 kHz, so the
    // noise in 2500 Hz is 300^2*2500/6000; a tone of amplitude A has A^2/2.
    long n = wspr15 ? 8L*114*12000 : 114L*12000;
    int nsps = wspr15 ? 8*8192 : 8192;
    double centre = wspr15 ? 1500.0 + 112.5 : 1500.0;
    double sigma = 300.0;
    double *x = calloc(n, sizeof(double));
    if (!x) return 1;
    if (noise) for (long j=0; j<n; j++) x[j]=sigma*gaussrand();
    for (int k=0; k<nsig; k++) {
        struct sig *s = &sigs[k];
        double pn = sigma*sigma*2500.0/6000.0;
        double amp = noise ? sqrt(2.0*pn*pow(10, s->snr/10.0)) : 8000.0/nsig;
        add_complex(x, NULL, n, 12000.0, nsps, centre + s->f, s->t0, s->drift, amp, s->symbols);
    }
    int ok = write_wav(wav, x, n);
    free(x);
    return ok ? 0 : 1;
}
