// wsprd's command line, files and output (WSJT-X's lib/wsprd/wsprd.c main(),
// readc2file, writec2file, usage, and wsprd_stream.c's -0 protocol, at
// 567ad29ce): see wsprd.h. run() is main() with the decoding left to
// Decoder, and stdout, stderr and stdin made Io's.
//
// wsprd.c: Copyright 2001-2018, Joe Taylor, K1JT, and 2014-2018, Steven
// Franke, K9AN; WSJT-X is GNU GPL v3, and so is this.
//
// Kept as the reference has them: the order files are opened and written
// (ALL_WSPR.TXT appended to and wspr_spots.txt truncated before the input
// is read, hashtable.txt created if missing only once it has been), each
// printf's format, the files' stdio modes (text on Windows, as the
// reference's), and -c's 000000_0001.c2 in the current directory rather
// than the data directory. getopt is parse_args() here; the input's FFTs
// (readwavfile, pcm_to_iq) are Decoder::decode_wav's.
#include "wsprd.h"
#include "wsprd_internal.h"

#include <cerrno>
#include <climits>
#include <cstdarg>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <string>
#include <vector>
#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#else
#include <unistd.h>
#endif

namespace wsprd {

namespace {

std::string format(const char* fmt, ...) {
    char buf[512];
    va_list args;
    va_start(args, fmt);
    int n = vsnprintf(buf, sizeof buf, fmt, args);
    va_end(args);
    if (n < 0) return std::string();
    if ((size_t)n < sizeof buf) return std::string(buf, (size_t)n);
    std::string s((size_t)n + 1, '\0');
    va_start(args, fmt);
    vsnprintf(&s[0], s.size(), fmt, args);
    va_end(args);
    s.resize((size_t)n);
    return s;
}

}  // namespace

//***************************************************************************
std::string usage(bool with_n)
{
    std::string s;
    s += "Usage: wsprd [options...] infile\n";
    s += "       infile must have suffix .wav or .c2\n";
    s += "\n";
    s += "Options:\n";
    s += "       -a <path> path to writeable data files, default=\".\"\n";
    s += "       -B disable block demodulation - use single-symbol noncoherent demod\n";
    s += "       -c write .c2 file at the end of the first pass\n";
    s += "       -C maximum number of decoder cycles per bit, default 10000\n";
    s += "       -d deeper search. Slower, a few more decodes\n";
    s += "       -e x (x is transceiver dial frequency error in Hz)\n";
    s += "       -f x (x is transceiver dial frequency in MHz)\n";
    s += "       -H do not use (or update) the hash table\n";
    s += "       -J use the stack decoder instead of Fano decoder\n";
    s += "       -m decode wspr-15 .wav file\n";
    s += "       -o n (0<=n<=5), decoding depth for OSD, default is disabled\n";
    s += "       -q quick mode - doesn't dig deep for weak signals\n";
    s += "       -s single pass mode, no subtraction (same as original wsprd)\n";
    s += "       -v verbose mode (shows dupes)\n";
    s += "       -w wideband mode - decode signals within +/- 150 Hz of center\n";
    s += "       -z x (x is fano metric table bias, default is 0.45)\n";
    if (with_n) s += "       -n no drift, does not attempt to follow a drifting frequency\n";
    s += "       -0 streaming mode: read PCM audio from stdin per the\n";
    s += "          streaming framing protocol; emit decodes as NDJSON\n";
    s += "          on stdout. Replaces the positional\n";
    s += "          input file argument.\n";
    return s;
}

bool parse_args(const std::vector<std::string>& args, Options* opt, std::string* error, const std::string& prog,
                bool with_n) {
    const char* optstring = with_n ? "a:BcC:de:f:HJmno:qstwvz:0" : "a:BcC:de:f:HJmo:qstwvz:0";
    *opt = Options();
    if (error) error->clear();
    auto apply = [opt](char c, const char* optarg) {
        switch (c) {
            case 'a':
                opt->data_dir = optarg;
                break;
            case 'B':
                opt->npasses=2;
                break;
            case 'c':
                opt->writec2=1;
                break;
            case 'C':
                opt->maxcycles=(unsigned int) strtoul(optarg,NULL,10);
                break;
            case 'd':
                opt->more_candidates=1;
                break;
            case 'e':
                opt->dialfreq_error = strtod(optarg,NULL);   // units of Hz
                // dialfreq_error = dial reading - actual, correct frequency
                break;
            case 'f':
                opt->dialfreq = strtod(optarg,NULL); // units of MHz
                break;
            case 'H':
                opt->usehashtable = 0;
                break;
            case 'J': //Stack (Jelinek) decoder, Fano decoder is the default
                opt->stackdecoder = 1;
                break;
            case 'n':  // not wsprd's: no coarse drift search (with_n only)
                opt->nodrift = 1;
                break;
            case 'm':  //15-minute wspr mode
                opt->wspr_type = 15;
                break;
            case 'o':  //use ordered-statistics-decoder
                opt->ndepth=(int) strtol(optarg,NULL,10);
                break;
            case 'q':  //no shift jittering
                opt->quickmode = 1;
                break;
            case 's':  //single pass mode
                opt->subtraction = 0;
                opt->npasses = 1;
                break;
            case 'v':
                opt->verbose = 1;
                break;
            case 'w':
                opt->fmin=-150.0;
                opt->fmax=150.0;
                break;
            case 'z':
                opt->bias=strtod(optarg,NULL); //fano metric bias (default is 0.45)
                break;
            case '0':  // --stream (long-name parity with jt9 --stream)
                opt->stream_mode = 1;
                break;
        }
    };
    // glibc's getopt: operands are moved after the options (so the first is
    // argv[optind]) unless POSIXLY_CORRECT stops it at the first; "--" ends
    // the options; "-" is an operand.
    const bool posix = std::getenv("POSIXLY_CORRECT") != nullptr;
    std::vector<std::string> operands;
    std::vector<int> operand_args;
    size_t a = 0;
    for (; a < args.size(); ++a) {
        const std::string& s = args[a];
        if (s == "--") {
            ++a;
            break;
        }
        if (s.size() < 2 || s[0] != '-') {
            if (posix) break;
            operands.push_back(s);
            operand_args.push_back((int)a);
            continue;
        }
        for (size_t c = 1; c < s.size(); ++c) {
            const char ch = s[c];
            const char* p = (ch == '\0' || ch == ':' || ch == ';') ? nullptr : std::strchr(optstring, ch);
            if (!p) {
                if (error) *error = prog + ": invalid option -- '" + ch + "'\n";
                return false;
            }
            if (p[1] != ':') {
                apply(ch, "");
                continue;
            }
            if (c + 1 < s.size()) {
                apply(ch, s.c_str() + c + 1);
            } else if (a + 1 < args.size()) {
                ++a;
                apply(ch, args[a].c_str());
            } else {
                if (error) *error = prog + ": option requires an argument -- '" + ch + "'\n";
                return false;
            }
            break;
        }
    }
    for (; a < args.size(); ++a) {
        operands.push_back(args[a]);
        operand_args.push_back((int)a);
    }
    if (!operands.empty()) {
        opt->infile = operands[0];
        opt->infile_arg = operand_args[0];
    }
    return true;
}

std::string stdout_line(const Decode& d) {
    return format("%4s %3.0f %4.1f %10.6f %2d  %-s \n",
                  d.time.c_str(), d.snr, d.dt, d.freq,
                  (int)d.drift, d.message.c_str());
}

std::string all_wspr_line(const Decode& d, bool with_spread) {
    std::string s = format("%6s %4s %3.0f %5.2f %11.7f  %-22s %2d %5.2f %2d %2d %4d %2d %3d %5u %5d",
                  d.date.c_str(), d.time.c_str(), d.snr,
                  d.dt, d.freq, d.message.c_str(),
                  (int)d.drift, d.sync,
                  d.ipass+1,d.blocksize,d.jitter,
                  d.decodetype,d.nhardmin,d.cycles/81,
                  d.metric);
    if (with_spread) s += format(" %6.3f", d.spread);
    return s + "\n";
}

std::string spots_line(const Decode& d) {
    return format("%6s %4s %3d %3.0f %4.1f %10.6f  %-22s %2d %5u %4d\n",
                  d.date.c_str(), d.time.c_str(), (int)(10*d.sync),
                  d.snr, d.dt, d.freq,
                  d.message.c_str(), (int)d.drift, d.cycles/81,
                  d.jitter);
}

// ---------------------------------------------------------------- wsprd_stream.c

namespace {

#define WSJT_MAGIC      "WSJT"
#define WSJT_FMT_PCM    0x00
#define FRAME_AUDIO     0x01
#define FRAME_CONTROL   0x02

// Buffer the producer's audio chunks. Sized for WSPR-15 worst case
// (8 * 114 * 12000 = 10,944,000 samples). WSPR-2 needs only 1,368,000.
#define MAX_PCM_SAMPLES (8 * 114 * 12000)
#define MAX_CONTROL_FRAME 1024U

struct wsprd_stream_config {
    char date[7];        // YYMMDD\0
    char uttime[5];      // HHMM\0
    double dialfreq;     // MHz
    int wspr_type;       // 2 or 15
    char mycall[13];     // up to 12-char callsign + NUL
    char mygrid[7];      // up to 6-char grid + NUL
    int configure_received;
};

const char *skip_space(const char *p) {
    while (*p == ' ' || *p == '\t' || *p == '\r' || *p == '\n') p++;
    return p;
}

const char *skip_json_string(const char *p) {
    if (*p != '"') return NULL;
    for (p++; *p; p++) {
        if (*p == '\\' && p[1]) {
            p++;
        } else if (*p == '"') {
            return p + 1;
        }
    }
    return NULL;
}

const char *json_find_key(const char *json, const char *key) {
    const size_t key_len = strlen(key);
    const char *p = skip_space(json);
    if (*p != '{') return NULL;
    p++;

    for (;;) {
        p = skip_space(p);
        if (*p == ',') {
            p = skip_space(p + 1);
        }
        if (*p == '}') return NULL;
        if (*p != '"') return NULL;

        const char *name = p + 1;
        const char *after_name = skip_json_string(p);
        if (!after_name) return NULL;
        const size_t name_len = (size_t)(after_name - name - 1);
        p = skip_space(after_name);
        if (*p != ':') return NULL;
        const char *value = skip_space(p + 1);
        if (name_len == key_len && memcmp(name, key, key_len) == 0) {
            return *value ? value : NULL;
        }

        if (*value == '"') {
            p = skip_json_string(value);
            if (!p) return NULL;
        } else {
            p = value;
            while (*p && *p != ',' && *p != '}') p++;
        }
        if (*p == '}') return NULL;
    }
}

int json_get_string(const char *json, const char *key,
                           char *out, size_t out_size) {
    const char *p = json_find_key(json, key);
    if (!p || *p != '"') return -1;
    p++;
    size_t i = 0;
    while (*p && *p != '"' && i + 1 < out_size) {
        out[i++] = *p++;
    }
    out[i] = '\0';
    return (*p == '"') ? 0 : -1;
}

int json_get_int(const char *json, const char *key, int *out) {
    const char *p = json_find_key(json, key);
    if (!p) return -1;
    char *end = NULL;
    errno = 0;
    long v = strtol(p, &end, 10);
    if (end == p || errno == ERANGE || v < INT_MIN || v > INT_MAX) return -1;
    end = (char *)skip_space(end);
    if (*end != ',' && *end != '}') return -1;
    *out = (int) v;
    return 0;
}

int json_get_double(const char *json, const char *key, double *out) {
    const char *p = json_find_key(json, key);
    if (!p) return -1;
    char *end = NULL;
    double v = strtod(p, &end);
    if (end == p) return -1;
    *out = v;
    return 0;
}

// JSON-escape a string into out. Best-effort: handles ASCII printable
// and the must-escape control chars per RFC 8259.
void json_escape_into(const char *src, char *out, size_t out_size) {
    size_t o = 0;
    for (const char *p = src; *p && o + 6 < out_size; p++) {
        unsigned char c = (unsigned char) *p;
        if (c == '"' || c == '\\') {
            if (o + 2 < out_size) { out[o++] = '\\'; out[o++] = c; }
        } else if (c < 0x20) {
            o += snprintf(out + o, out_size - o, "\\u%04x", c);
        } else {
            out[o++] = c;
        }
    }
    out[o] = '\0';
}

void wsprd_stream_emit_ready(const Io& io) {
    io.out("{\"v\":1,\"t\":\"ready\",\"modes\":[\"WSPR\"],\"protocol\":1}\n");
    io.flush();
}

std::string stream_decode_line(const char *date, const char *uttime,
                              float snr, float dt, double freq,
                              int drift, const char *message) {
    char esc_msg[80];
    json_escape_into(message ? message : "", esc_msg, sizeof esc_msg);
    return format("{\"v\":1,\"t\":\"decode\",\"mode\":\"WSPR\","
            "\"date\":\"%s\",\"time\":\"%s\","
            "\"snr\":%.0f,\"dt\":%.1f,\"freq\":%.6f,"
            "\"drift\":%d,\"message\":\"%s\"}\n",
            date ? date : "",
            uttime ? uttime : "",
            snr, dt, freq, drift, esc_msg);
}

void wsprd_stream_emit_decode_finished(const Io& io, const char *date, const char *uttime) {
    io.out(format("{\"v\":1,\"t\":\"decode_finished\","
            "\"date\":\"%s\",\"time\":\"%s\"}\n",
            date ? date : "",
            uttime ? uttime : ""));
    io.flush();
}

void wsprd_stream_emit_error_code(const Io& io, const char *code, const char *detail) {
    char esc_code[80];
    char esc_detail[256];
    json_escape_into(code ? code : "", esc_code, sizeof esc_code);
    json_escape_into(detail ? detail : "", esc_detail, sizeof esc_detail);
    io.out(format("{\"v\":1,\"t\":\"error\",\"code\":\"%s\",\"detail\":\"%s\"}\n",
            esc_code, esc_detail));
    io.flush();
}

void set_stdin_binary(FILE* in) {
#ifdef _WIN32
    _setmode(_fileno(in), _O_BINARY);
#else
    (void)in;
#endif
}

int wsjt_header_read(const Io& io, FILE *fp) {
    unsigned char hdr[8];
    if (fread(hdr, 1, 8, fp) != 8) {
        io.err("wsprd --stream: short read on WSJT header\n");
        return -1;
    }
    if (memcmp(hdr, WSJT_MAGIC, 4) != 0) {
        io.err("wsprd --stream: bad magic; expected 'WSJT'\n");
        return -1;
    }
    unsigned char fmt = hdr[4];
    unsigned char ch = hdr[5];
    uint16_t rate_khz = (uint16_t)hdr[6] | ((uint16_t)hdr[7] << 8);
    if (fmt != WSJT_FMT_PCM) {
        io.err(format(
                "wsprd --stream: unsupported fmt=0x%02x (only 0x00 PCM today)\n",
                fmt));
        return -1;
    }
    if (ch != 1) {
        io.err(format(
                "wsprd --stream: unsupported channels=%u (only mono supported)\n",
                ch));
        return -1;
    }
    if (rate_khz != 12) {
        io.err(format(
                "wsprd --stream: unsupported rate=%u kHz (only 12 kHz supported)\n",
                rate_khz));
        return -1;
    }
    io.err(format(
            "wsprd --stream: header ok (fmt=%u ch=%u rate=%u kHz)\n",
            fmt, ch, rate_khz));
    return 0;
}

// Reads one frame envelope. *body is malloc'd on success; caller frees.
// Returns 0 ok, -1 err, +1 EOF (clean producer-disconnect).
int frame_read(const Io& io, FILE *fp, unsigned char *type, void **body, uint32_t *len) {
    unsigned char hdr[5];
    size_t nr = fread(hdr, 1, 5, fp);
    if (nr == 0) return 1;  // clean EOF
    if (nr != 5) {
        io.err("wsprd --stream: short read on frame header\n");
        return -1;
    }
    *type = hdr[0];
    *len = (uint32_t)hdr[1] | ((uint32_t)hdr[2] << 8) |
           ((uint32_t)hdr[3] << 16) | ((uint32_t)hdr[4] << 24);
    const uint32_t max_len = *type == FRAME_CONTROL
        ? MAX_CONTROL_FRAME
        : (uint32_t)(MAX_PCM_SAMPLES * sizeof(int16_t));
    if (*len > max_len) {
        char detail[160];
        snprintf(detail, sizeof detail,
                 "frame length %u exceeds type capacity %u", *len, max_len);
        wsprd_stream_emit_error_code(io, "frame_too_large", detail);
        io.err(format("wsprd --stream: %s\n", detail));
        return -1;
    }
    if (*len == 0) {
        *body = NULL;
        return 0;
    }
    *body = malloc(*len);
    if (!*body) {
        io.err(format("wsprd --stream: alloc %u bytes failed\n", *len));
        wsprd_stream_emit_error_code(io, "allocation_failed",
                                     "unable to allocate frame body");
        return -1;
    }
    if (fread(*body, 1, *len, fp) != *len) {
        io.err(format("wsprd --stream: short read on frame body (%u bytes)\n",
                *len));
        free(*body);
        *body = NULL;
        wsprd_stream_emit_error_code(io, "frame_read_error",
                                     "short read on frame body");
        return -1;
    }
    return 0;
}

// A RIFF/WAVE file of 32-bit IEEE float samples (WAVE_FORMAT_IEEE_FLOAT, or
// WAVE_FORMAT_EXTENSIBLE with that subformat), as radiod and sox write them:
// true, with its data chunk's offset and length in bytes, when fp is one
// (mono). Any other file is false and read as wsprd reads it: 16-bit
// samples after a 44-byte header, whatever its header says.
bool wav_float_data(FILE* fp, long& offset, long& length, int& channels) {
    unsigned char h[12];
    if (fseek(fp, 0, SEEK_SET) != 0 || fread(h, 1, 12, fp) != 12) return false;
    if (memcmp(h, "RIFF", 4) != 0 || memcmp(h + 8, "WAVE", 4) != 0) return false;
    auto u16 = [](const unsigned char* b) { return (unsigned)b[0] | ((unsigned)b[1] << 8); };
    auto u32 = [](const unsigned char* b) {
        return (unsigned long)b[0] | ((unsigned long)b[1] << 8) | ((unsigned long)b[2] << 16) |
               ((unsigned long)b[3] << 24);
    };
    bool is_float = false;
    long pos = 12;
    for (;;) {
        unsigned char c[8];
        if (fseek(fp, pos, SEEK_SET) != 0 || fread(c, 1, 8, fp) != 8) return false;
        const unsigned long size = u32(c + 4);
        if (memcmp(c, "fmt ", 4) == 0) {
            unsigned char f[40] = {};
            if (size < 16 || fread(f, 1, size < 40 ? size : 40, fp) < 16) return false;
            unsigned tag = u16(f);
            channels = (int)u16(f + 2);
            const unsigned bits = u16(f + 14);
            if (tag == 0xFFFE && size >= 40) tag = u16(f + 24);  // the subformat GUID's first two bytes
            is_float = tag == 3 && bits == 32;
            if (!is_float) return false;
        } else if (memcmp(c, "data", 4) == 0) {
            if (!is_float) return false;
            offset = pos + 8;
            // A writer that never came back to fill in the size (0, or
            // 0xFFFFFFFF while it streams) or one larger than the file: its
            // samples run to the end of the file.
            long end = offset;
            if (fseek(fp, 0, SEEK_END) == 0) end = ftell(fp);
            length = (size == 0 || size == 0xFFFFFFFFul || offset + (long)size > end) ? end - offset : (long)size;
            return true;
        }
        pos += 8 + (long)size + (long)(size & 1);
    }
}

// pcm_to_iq()'s checks; its FFTs are Decoder::decode_wav's. Returns the
// samples to decode (expected_npcm), or 0.
size_t pcm_check(const Io& io, size_t npcm, int ntrmin) {
    size_t expected_npcm;

    if (ntrmin == 2) {
        expected_npcm = 114 * 12000;
    } else if (ntrmin == 15) {
        expected_npcm = 8 * 114 * 12000;
    } else {
        io.err(format("wsprd --stream: invalid wspr_type=%d\n", ntrmin));
        wsprd_stream_emit_error_code(io, "invalid_wspr_type",
                                     "wspr_type must be 2 or 15");
        return 0;
    }

    if (npcm < expected_npcm) {
        io.err(format(
                "wsprd --stream: short PCM (got %zu, expected %zu for WSPR-%d)\n",
                npcm, expected_npcm, ntrmin));
        wsprd_stream_emit_error_code(io, "insufficient_audio",
                                     "not enough PCM samples for the configured WSPR type");
        return 0;
    }
    return expected_npcm;
}

// wsprd_stream_read() up to its pcm_to_iq(): the configuration and the
// PCM. Returns the samples to decode, 0 on a protocol error.
size_t wsprd_stream_read(const Io& io, struct wsprd_stream_config *cfg, std::vector<int16_t>& pcm) {
    FILE* in = io.in ? io.in : stdin;
    set_stdin_binary(in);
    if (wsjt_header_read(io, in) != 0) return 0;

    pcm.assign(MAX_PCM_SAMPLES, 0);

    memset(cfg, 0, sizeof *cfg);
    cfg->wspr_type = 2;
    int have_date = 0;
    int have_time = 0;
    int have_dialfreq = 0;
    int halted = 0;
    int wspr_type_value = 0;
    size_t pcm_n = 0;

    while (!halted) {
        unsigned char type;
        void *body = NULL;
        uint32_t len = 0;
        int rc = frame_read(io, in, &type, &body, &len);
        if (rc < 0) {
            return 0;
        }
        if (rc > 0) break;  // EOF without halt — treat as drain-and-decode

        if (type == FRAME_CONTROL) {
            char *json = (char *)malloc((size_t)len + 1);
            if (!json) {
                free(body);
                wsprd_stream_emit_error_code(io, "allocation_failed",
                                             "unable to allocate control JSON buffer");
                return 0;
            }
            if (len != 0) {
                memcpy(json, body, len);
            }
            json[len] = '\0';
            free(body);

            char event_type[16];
            if (json_get_string(json, "t", event_type, sizeof event_type) != 0) {
                free(json);
                wsprd_stream_emit_error_code(io, "control_parse_error",
                                             "control frame is missing a string t field");
                continue;
            }
            if (strcmp(event_type, "halt") == 0) {
                halted = 1;
                free(json);
                continue;
            }
            if (strcmp(event_type, "configure") == 0) {
                if (json_get_string(json, "date", cfg->date,
                                    sizeof cfg->date) == 0) {
                    have_date = 1;
                }
                if (json_get_string(json, "time", cfg->uttime,
                                    sizeof cfg->uttime) == 0) {
                    have_time = 1;
                }
                if (json_get_double(json, "dialfreq", &cfg->dialfreq) == 0) {
                    have_dialfreq = 1;
                }
                if (json_find_key(json, "wspr_type")) {
                    if (json_get_int(json, "wspr_type", &wspr_type_value) != 0 ||
                        (wspr_type_value != 2 && wspr_type_value != 15)) {
                        free(json);
                        wsprd_stream_emit_error_code(io, "invalid_wspr_type",
                                                     "wspr_type must be 2 or 15");
                        return 0;
                    }
                    cfg->wspr_type = wspr_type_value;
                }
                json_get_string(json, "mycall", cfg->mycall, sizeof cfg->mycall);
                json_get_string(json, "mygrid", cfg->mygrid, sizeof cfg->mygrid);
                cfg->configure_received = 1;
            }
            free(json);
            continue;
        }

        if (type == FRAME_AUDIO) {
            // body is raw int16 LE PCM. len must be even (whole samples).
            if (len % sizeof(int16_t) != 0) {
                free(body);
                wsprd_stream_emit_error_code(io, "odd_audio_frame",
                                             "audio frame length must contain whole int16 samples");
                return 0;
            }
            size_t nsamples = len / sizeof(int16_t);
            if (pcm_n + nsamples > MAX_PCM_SAMPLES) {
                nsamples = MAX_PCM_SAMPLES - pcm_n;
            }
            if (nsamples != 0) {
                memcpy(pcm.data() + pcm_n, body, nsamples * sizeof(int16_t));
            }
            pcm_n += nsamples;
            free(body);
            continue;
        }

        // Unknown frame type — ignore (forward-compat).
        free(body);
    }

    if (!cfg->configure_received) {
        wsprd_stream_emit_error_code(io, "missing_configure",
                                     "no configure frame received before halt or EOF");
        return 0;
    }
    if (!have_date || !have_time || !have_dialfreq) {
        wsprd_stream_emit_error_code(io, "incomplete_configure",
                                     "configure frame missing date, time, or dialfreq");
        return 0;
    }

    return pcm_check(io, pcm_n, cfg->wspr_type);
}

}  // namespace

std::string stream_line(const Decode& d) {
    return stream_decode_line(d.date.c_str(), d.time.c_str(), d.snr, d.dt, d.freq, (int)d.drift, d.message.c_str());
}

// ------------------------------------------------------------ readc2file etc.

namespace {

// readc2file, opening path but naming name in its message.
bool read_c2_as(const std::string& path, const std::string& name, C2* c2, std::string* error)
{
    float *buffer;
    double dfreq;
    int i,ntrmin;
    char c2file[15];
    size_t nr;
    FILE* fp;
    if (error) error->clear();

    fp = fopen(path.c_str(),"rb");
    if (fp == NULL) {
        if (error) *error = format("Cannot open data file '%s': %s\n", name.c_str(), strerror(errno));
        return false;
    }
    nr=fread(c2file,sizeof(char),14,fp);
    nr=fread(&ntrmin,sizeof(int),1,fp);
    nr=fread(&dfreq,sizeof(double),1,fp);
    c2->freq=dfreq;

    buffer=(float *)calloc(2*65536,sizeof(float));
    nr=fread(buffer,sizeof(float),2*45000,fp);
    fclose(fp);

    c2->wspr_type=ntrmin;

    c2->idat.assign(45000, 0.0f);
    c2->qdat.assign(45000, 0.0f);
    for(i=0; i<45000; i++) {
        c2->idat[i]=buffer[2*i];
        c2->qdat[i]=-buffer[2*i+1];
    }
    free(buffer);

    return nr == 2*45000;
}

}  // namespace

bool read_c2(const std::string& path, C2* c2, std::string* error)
{
    return read_c2_as(path, path, c2, error);
}

std::string c2_bytes(const char* c2filename, int trmin, double freq, const float* idat, const float* qdat)
{
    std::string out;
    out.append(c2filename, 14);
    out.append(reinterpret_cast<const char*>(&trmin), sizeof(int));
    out.append(reinterpret_cast<const char*>(&freq), sizeof(double));
    std::vector<float> buffer(2*45000);
    for(int i=0; i<45000; i++) {
        buffer[2*i]=idat[i];
        buffer[2*i+1]=-qdat[i];
    }
    out.append(reinterpret_cast<const char*>(buffer.data()), sizeof(float)*buffer.size());
    return out;
}

// ----------------------------------------------------------------------- main

namespace {

int build_data_path(const Io& io, std::string& destination, size_t capacity,
                           const char *data_dir, const char *filename)
{
    const char *base = data_dir ? data_dir : ".";
    std::string written = std::string(base) + "/" + filename;
    if (written.size() >= capacity) {
        io.err(format("Error: data path for '%s' exceeds %zu bytes\n",
                filename, capacity - 1));
        return -1;
    }
    destination = written;
    return 0;
}

int data_dir_accessible(const char* dir) {
#ifdef _WIN32
    return _access(dir, 6);
#else
    return access(dir, R_OK | W_OK);
#endif
}

}  // namespace

int run(const std::vector<std::string>& argv, const Io& io)
{
    clock_t t0,t00;
    t00=detail::cpu_clock();

    Options opt;
    std::string error;
    std::vector<std::string> args;
    if (argv.size() > 1) args.assign(argv.begin() + 1, argv.end());
    if (!parse_args(args, &opt, &error, argv.empty() ? std::string("wsprd") : argv[0], io.with_n)) {
        io.err(error);
        io.out(usage(io.with_n));
        return 1;
    }
    const char* data_dir = opt.data_dir.c_str();
    // A relative path (the input, -a's) is io.cwd's; "" stays "", which no
    // call opens, as in any directory.
    auto local = [&](const std::string& path) {
#ifdef _WIN32
        const bool absolute = path[0] == '/' || path[0] == '\\' || (path.size() > 1 && path[1] == ':');
#else
        const bool absolute = path[0] == '/';
#endif
        return io.cwd.empty() || path.empty() || absolute ? path : io.cwd + "/" + path;
    };

    if( data_dir_accessible(local(data_dir).c_str())) {
      io.err(format("Error: inaccessible data directory: '%s'\n", data_dir));
      io.out(usage(io.with_n));
      return EXIT_FAILURE;
    }

    if( !opt.stream_mode && opt.infile_arg < 0) {
        io.out(usage(io.with_n));
        return 1;
    }
    // argv's strings back to back, each with its NUL, as the kernel lays
    // them out: a name too short for the date or time makes the reference
    // read the end of the argument before it (in the order given, as getopt
    // moves only the pointers). Before argv[0] reads as NUL.
    std::string argv_area;
    size_t infile_at = 0;
    for (size_t i = 0; i < argv.size(); ++i) {
        if ((int)i == opt.infile_arg + 1) infile_at = argv_area.size();
        argv_area += argv[i];
        argv_area += '\0';
    }
    if (opt.infile_arg < 0) {
        infile_at = argv_area.size();
        argv_area += opt.infile;
        argv_area += '\0';
    }
    const char* ptr_to_infile = argv_area.c_str() + infile_at;
    const char* ptr_to_infile_suffix = NULL;
    const std::string infile_path = local(ptr_to_infile);

    Decoder decoder(opt);
    Timings& tm = decoder.timings();

    std::string wisdom_fname,all_fname,spots_fname,timer_fname,hash_fname;
    FILE *fp_fftwf_wisdom_file, *fall_wspr, *fwsprd, *fhash, *ftimer;
    if (build_data_path(io, wisdom_fname, 200, data_dir,
                        "wspr_wisdom.dat") != 0 ||
        build_data_path(io, all_fname, 200, data_dir,
                        "ALL_WSPR.TXT") != 0 ||
        build_data_path(io, spots_fname, 200, data_dir,
                        "wspr_spots.txt") != 0 ||
        build_data_path(io, timer_fname, 200, data_dir,
                        "wspr_timer.out") != 0 ||
        build_data_path(io, hash_fname, 200, data_dir,
                        "hashtable.txt") != 0) {
        return EXIT_FAILURE;
    }
    for (std::string* f : {&wisdom_fname, &all_fname, &spots_fname, &timer_fname, &hash_fname}) *f = local(*f);
    if ((fp_fftwf_wisdom_file = fopen(wisdom_fname.c_str(), "r"))) {  //Open FFTW wisdom
        std::string text;
        char buf[4096];
        size_t n;
        while ((n = fread(buf, 1, sizeof buf, fp_fftwf_wisdom_file)) > 0) text.append(buf, n);
        decoder.wisdom() = text;
        fclose(fp_fftwf_wisdom_file);
    }

    fall_wspr=fopen(all_fname.c_str(),"a");
    fwsprd=fopen(spots_fname.c_str(),"w");

    if((ftimer=fopen(timer_fname.c_str(),"r"))) {
        //Accumulate timing data
        int nr=fscanf(ftimer,"%f %f %f %f %f %f %f %f",
               &tm.treadwav,&tm.tcandidates,&tm.tsync0,&tm.tsync1,&tm.tsync2,&tm.tfano,&tm.tosd,&tm.ttotal);
        fclose(ftimer);
        if(nr == 0) io.err(format("Empty timer file: '%s'\n", timer_fname.c_str()));
    }
    ftimer=fopen(timer_fname.c_str(),"w");

    auto close_files = [&] {
        if (fall_wspr) fclose(fall_wspr);
        if (fwsprd) fclose(fwsprd);
        if (ftimer) fclose(ftimer);
    };

    char uttime[5],date[7];
    // strncpy(date,ptr_to_infile_suffix-11,6) and (uttime,...-4,4), into
    // argv_area (see above).
    auto parse_date = [&] {
        const char* start = argv_area.c_str();
        auto copy = [start](char* dst, const char* src, int n) {
            int k = 0;
            for (; k < n; ++k) {
                const char c = src + k < start ? '\0' : src[k];
                if (c == '\0') break;
                dst[k] = c;
            }
            for (; k < n; ++k) dst[k] = '\0';
        };
        copy(date, ptr_to_infile_suffix-11, 6);
        copy(uttime, ptr_to_infile_suffix-4, 4);
        date[6]='\0';
        uttime[4]='\0';
    };

    enum { NONE, WAV, C2IN } kind = NONE;
    std::vector<int16_t> pcm;
    std::vector<float> pcmf;  // a float .wav's samples (pcm then unused)
    size_t npcm = 0;
    C2 c2;
    if (opt.stream_mode) {
        // Read PCM from stdin per the streaming protocol. Configure frame
        // supplies date/time/dialfreq/wspr_type; audio frames carry int16
        // PCM @ 12 kHz mono.
        wsprd_stream_config stream_cfg = {};
        wsprd_stream_emit_ready(io);
        t0 = detail::cpu_clock();
        npcm = wsprd_stream_read(io, &stream_cfg, pcm);
        tm.treadwav += (float)(detail::cpu_clock()-t0)/CLOCKS_PER_SEC;
        if (npcm == 0) { close_files(); return 1; }
        decoder.options().wspr_type = stream_cfg.wspr_type;
        decoder.options().dialfreq = stream_cfg.dialfreq;
        strncpy(date, stream_cfg.date, 6); date[6] = '\0';
        strncpy(uttime, stream_cfg.uttime, 4); uttime[4] = '\0';
        kind = WAV;
    } else if( strstr(ptr_to_infile,".wav") ) {
        ptr_to_infile_suffix=strstr(ptr_to_infile,".wav");

        // readwavfile's reading
        t0 = detail::cpu_clock();
        const size_t npoints = opt.wspr_type == 15 ? 8*114*12000 : 114*12000;
        FILE* fp = fopen(infile_path.c_str(),"rb");
        if (fp == NULL) {
            io.err(format("Cannot open data file '%s'\n", ptr_to_infile));
            close_files();
            return 1;
        }
        long data_offset = 0, data_length = 0;
        int channels = 0;
        if (wav_float_data(fp, data_offset, data_length, channels)) {
            if (channels != 1) {
                io.err(format("Float .wav file '%s' has %d channels, not 1\n", ptr_to_infile, channels));
                fclose(fp);
                close_files();
                return 1;
            }
            pcmf.assign(npoints, 0.0f);
            size_t want = (size_t)data_length / sizeof(float);
            if (want > npoints) want = npoints;
            size_t nr = 0;
            if (fseek(fp, data_offset, SEEK_SET) == 0) nr = fread(pcmf.data(), sizeof(float), want, fp);
            fclose(fp);
            if( nr == 0 ) {
                io.err(format("No data in file '%s'\n", ptr_to_infile));
                close_files();
                return 1;
            }
            npcm = npoints;
            tm.treadwav += (float)(detail::cpu_clock()-t0)/CLOCKS_PER_SEC;
            parse_date();
            kind = WAV;
        } else {
        fseek(fp, 0, SEEK_SET);
        // As calloc and glibc's fread leave buf2: the header's 22 shorts
        // where the data does not reach, and the low byte of a last, half
        // sample (fread stores it but does not count it). Read as bytes so
        // that holds with any C library.
        pcm.assign(npoints, 0);
        char* bytes = reinterpret_cast<char*>(pcm.data());
        size_t nr=fread(bytes,1,2*22,fp)/2;      //Read and ignore header
        nr=fread(bytes,1,2*npoints,fp)/2; //Read raw data
        fclose(fp);
        if( nr == 0 ) {
            io.err(format("No data in file '%s'\n", ptr_to_infile));
            close_files();
            return 1;
        }
        npcm = npoints;
        tm.treadwav += (float)(detail::cpu_clock()-t0)/CLOCKS_PER_SEC;
        parse_date();
        kind = WAV;
        }
    } else if ( strstr(ptr_to_infile,".c2") !=0 )  {
        ptr_to_infile_suffix=strstr(ptr_to_infile,".c2");
        if (!read_c2_as(infile_path, ptr_to_infile, &c2, &error)) {
            io.err(error);
            close_files();
            return 1;
        }
        parse_date();
        kind = C2IN;
    } else {
        io.out(format("Error: Failed to open %s\n",ptr_to_infile));
        io.out("WSPR file must have suffix .wav or .c2\n");
        close_files();
        return 1;
    }

    if( opt.usehashtable ) {
        char line[80];
        if( (fhash=fopen(hash_fname.c_str(),"r+")) ) {
            while (fgets(line, sizeof(line), fhash) != NULL) {
                decoder.hashtable().load_line(line);
            }
        } else {
            fhash=fopen(hash_fname.c_str(),"w+");
        }
        // fclose(NULL) when hashtable.txt can be neither updated nor made:
        // the reference crashes before printing anything.
        if (!fhash) {
            close_files();
            if (io.crash) io.crash();
            return 139;
        }
        fclose(fhash);
    }

    Result res = kind != WAV    ? decoder.decode_c2(c2, date, uttime)
                 : pcmf.empty() ? decoder.decode_wav(pcm.data(), npcm, date, uttime)
                                : decoder.decode_wav(pcmf.data(), npcm, date, uttime);

    // run() refuses what the Decoder would (as wsprd does) before here.
    if (res.status == 1) {
        close_files();
        return 1;
    }

    // Where the reference crashes (SIGSEGV, 139): in the stack decoder (see
    // Result::status), and where it writes to a file it could not open
    // (fprintf or fclose of NULL): at the first decode without
    // ALL_WSPR.TXT or wspr_spots.txt, at the timings without
    // wspr_timer.out, at closing those two (with no decodes), at the end
    // without hashtable.txt. What stdio still held is then lost: io.crash
    // ends the process there to lose it as wsprd does; without it, stdout
    // (unless -0, which flushes as it goes) and the files are taken as
    // fully buffered, written 4096 bytes at a time. hashtable.txt is opened
    // "w" (emptied) at the very end, as wsprd does: till then whether it
    // crashes there is not known, so stdout is held when it might.
    enum { NO_CRASH, AT_DECODE, AT_PRINT, AT_TIMER, AT_CLOSE, AT_HASH } crash = NO_CRASH;
    if (res.status == 139) crash = AT_DECODE;
    else if (!res.decodes.empty() && (!fall_wspr || !fwsprd)) crash = AT_PRINT;
    else if (!ftimer) crash = AT_TIMER;
    else if (!fall_wspr || !fwsprd) crash = AT_CLOSE;
    const bool hold = (crash != NO_CRASH || opt.usehashtable) && !io.crash;
    std::string held_out, held_all, held_spots, held_timer;
    auto out = [&](const std::string& s) {
        if (hold && !opt.stream_mode) held_out += s;
        else io.out(s);
    };
    auto put = [&](FILE* f, std::string& held, const std::string& s) {
        if (hold && (crash == AT_TIMER || crash == AT_CLOSE)) held += s;
        else fputs(s.c_str(), f);
    };
    auto crashed = [&] {
        if (io.crash) {
            io.crash();
        } else {
            auto flushed = [](const std::string& s) { return s.substr(0, s.size() / 4096 * 4096); };
            if (!opt.stream_mode) io.out(flushed(held_out));
            if (fall_wspr) fputs(flushed(held_all).c_str(), fall_wspr);
            if (fwsprd) fputs(flushed(held_spots).c_str(), fwsprd);
            if (ftimer) fputs(flushed(held_timer).c_str(), ftimer);
        }
        close_files();
        return 139;
    };

    if( res.has_c2 ) {
        char c2filename[15];
        strcpy(c2filename,"000000_0001.c2");
        out(format("Writing %s\n",c2filename));
        const std::string path = io.cwd.empty() ? std::string(c2filename) : io.cwd + "/" + c2filename;
        FILE* fp = fopen(path.c_str(),"wb");
        if( fp == NULL ) {
            io.err(format("Could not open c2 file '%s'\n", c2filename));
        } else {
            const std::string bytes = c2_bytes(c2filename, 2, res.c2_freq, res.c2_idat.data(), res.c2_qdat.data());
            fwrite(bytes.data(), 1, bytes.size(), fp);
            fclose(fp);
        }
    }
    if (crash == AT_DECODE) return crashed();

    for (const Decode& d : res.decodes) {
        if (opt.stream_mode) {
            io.out(stream_decode_line(date, d.time.c_str(), d.snr, d.dt, d.freq, (int)d.drift, d.message.c_str()));
            io.flush();
        } else {
            out(stdout_line(d));
        }
        if (crash == AT_PRINT) return crashed();
        put(fall_wspr, held_all, all_wspr_line(d, opt.nodrift));
        put(fwsprd, held_spots, spots_line(d));
    }
    if (opt.stream_mode) {
        wsprd_stream_emit_decode_finished(io, date, uttime);
    } else {
        out("<DecodeFinished>\n");
    }

    if ((fp_fftwf_wisdom_file = fopen(wisdom_fname.c_str(), "w"))) {
        fputs(decoder.wisdom().c_str(), fp_fftwf_wisdom_file);
        fclose(fp_fftwf_wisdom_file);
    }

    tm.ttotal += (float)(detail::cpu_clock()-t00)/CLOCKS_PER_SEC;

    if (crash == AT_TIMER) return crashed();
    std::string timer = format("%7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f\n\n",
            tm.treadwav,tm.tcandidates,tm.tsync0,tm.tsync1,tm.tsync2,tm.tfano,tm.tosd,tm.ttotal);

    timer += "Code segment        Seconds   Frac\n";
    timer += "-----------------------------------\n";
    timer += format("readwavfile        %7.2f %7.2f\n",tm.treadwav,tm.treadwav/tm.ttotal);
    timer += format("Coarse DT f0 f1    %7.2f %7.2f\n",tm.tcandidates,
            tm.tcandidates/tm.ttotal);
    timer += format("sync_and_demod(0)  %7.2f %7.2f\n",tm.tsync0,tm.tsync0/tm.ttotal);
    timer += format("sync_and_demod(1)  %7.2f %7.2f\n",tm.tsync1,tm.tsync1/tm.ttotal);
    timer += format("sync_and_demod(2)  %7.2f %7.2f\n",tm.tsync2,tm.tsync2/tm.ttotal);
    timer += format("Stack/Fano decoder %7.2f %7.2f\n",tm.tfano,tm.tfano/tm.ttotal);
    timer += format("OSD        decoder %7.2f %7.2f\n",tm.tosd,tm.tosd/tm.ttotal);
    timer += "-----------------------------------\n";
    timer += format("Total              %7.2f %7.2f\n",tm.ttotal,1.0);
    put(ftimer, held_timer, timer);

    if (crash == AT_CLOSE) {
        // fclose(fall_wspr), then fclose(fwsprd): the first NULL crashes
        if (fall_wspr) {
            fclose(fall_wspr);
            fall_wspr = NULL;
        }
        return crashed();
    }
    // the files are closed (and so written in full) before hashtable.txt
    close_files();
    fall_wspr = fwsprd = ftimer = NULL;

    if( opt.usehashtable ) {
        if (!(fhash = fopen(hash_fname.c_str(), "w"))) {
            crash = AT_HASH;
            return crashed();
        }
        fputs(decoder.hashtable().save().c_str(), fhash);
        fclose(fhash);
    }
    if (hold && !opt.stream_mode) io.out(held_out);

    return 0;
}

}  // namespace wsprd
