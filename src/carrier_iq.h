// Not WSJT-X's: a decoded signal's carrier, its modulation taken off.
//
// Once a WSPR or FST4W signal is decoded its tones are known, so the samples
// times the conjugate of a reference made from them (at a constant frequency,
// no drift) leave a single carrier: what remains is what the path did to it,
// its Doppler shift and spread over the whole transmission, not just one
// frequency at its middle. wsprd -I DIR and jt9 --fst4w --carrier-dir DIR
// write that product for each decode, low-passed and decimated to RATE_HZ
// (complex, +-7.5 Hz around the reference), as a 2-channel 32-bit float WAV
// (I, Q) from the start of the transmission to its end, with the decode's
// particulars as JSON in the WAV's LIST/INFO/ICMT chunk.
//
// The carrier's frequency at time t0_s + n/rate_hz is ref_hz plus the
// baseband frequency of sample n (Hz, the same audio scale as the decoder's
// own: add the dial frequency for RF). ref_hz is the decoder's frequency:
// FST4W's lowest tone (jt9's fsig), WSPR's centre (wsprd's).
#pragma once

#include <complex>
#include <string>
#include <vector>

namespace carrier {

constexpr unsigned RATE_HZ = 15;

// x (n samples at fs_in) low-passed (flat to pass_hz, down by stop_hz:
// Blackman-windowed sinc) and kept every factor'th sample from x[0], each
// output centred on its input sample (no delay; zeros taken outside x).
std::vector<std::complex<float>> decimate(const std::complex<float>* x, long n, double fs_in, int factor,
                                          double pass_hz, double stop_hz);

// The IQ as a 2-channel IEEE-float WAV at rate Hz, comment (JSON) in
// LIST/INFO/ICMT. False, why in *err, if it cannot be written.
bool write_wav(const std::string& path, const std::vector<std::complex<float>>& iq, unsigned rate,
               const std::string& comment, std::string* err);

// A file name for a decode: the message's callsign ("<K1ABC>" as K1ABC, '/'
// as '_'; empty for an unknown hash "<...>") and its frequency, "K1ABC_1500.124.wav".
std::string file_name(const std::string& message, double ref_hz);

// s as a JSON string (quoted, escaped).
std::string json_string(const std::string& s);

}  // namespace carrier
