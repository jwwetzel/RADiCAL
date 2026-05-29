// ============================================================================
// WaveformUtils.h — waveform signal extraction for CAEN DT5742 (DRS4) data
// ============================================================================
//
// DT5742 pulses are negative-going (voltage drops when signal arrives).
// The pedestal is estimated from early, pre-signal samples.
// All amplitudes returned are positive (pedestal-subtracted and inverted).
// ============================================================================

#ifndef WAVEFORMUTILS_H
#define WAVEFORMUTILS_H

#include <cmath>
#include <cfloat>

static const float kNoTime = -1e6f;  // sentinel for "not found"

struct Pulse {
    float pedestal;      // baseline [mV], averaged over early samples
    float peak;          // amplitude above pedestal (>0) [mV]
    float peakTime;      // time of peak sample [ns]
    float crossingTime;  // CFD rising-edge crossing time [ns], kNoTime if not found
    float charge;        // integral of pedestal-subtracted signal over the pulse window [mV*sample]
    bool  valid;         // false if peak is below noise floor
};

// ---------------------------------------------------------------------------
// ExtractPulse
//   time[0..1023] and amp[0..1023] are pointers into the flat event arrays
//   at the appropriate offsets (see ChannelConfig.h).
//
//   fraction  — CFD threshold expressed as a fraction of the peak amplitude.
//               0.20 (20%) is a good default for timing; gives earlier crossing
//               than 50% and is less susceptible to time-walk at low amplitude.
//   minPeak   — minimum peak amplitude required to declare a valid pulse [mV].
// ---------------------------------------------------------------------------
inline Pulse ExtractPulse(const float* time, const float* amp,
                           float fraction = 0.20f, float minPeak = 5.0f)
{
    Pulse p;
    p.crossingTime = kNoTime;
    p.peakTime     = kNoTime;
    p.charge       = 0.f;
    p.valid        = false;

    // Pedestal: mean of samples 3..52  (avoids the very first samples which
    // can be noisy in DRS4, and stops well before any typical signal arrival)
    double ped = 0.0;
    for (int i = 3; i < 53; ++i) ped += amp[i];
    p.pedestal = static_cast<float>(ped / 50.0);

    // Find the minimum sample (peak of negative-going pulse) in range 3..1023
    float vmin = amp[3];
    int   imin = 3;
    for (int i = 4; i < 1024; ++i) {
        if (amp[i] < vmin) { vmin = amp[i]; imin = i; }
    }

    p.peak     = p.pedestal - vmin;   // positive value
    p.peakTime = time[imin];
    p.valid    = (p.peak >= minPeak);

    if (!p.valid) return p;

    // CFD: find the first sample (scanning forward from sample 3 to the peak)
    // where the rising signal crosses fraction * peak above the pedestal.
    // Linear interpolation between adjacent samples for sub-sample precision.
    float target = fraction * p.peak;
    for (int i = 3; i < imin; ++i) {
        float sig_i   = p.pedestal - amp[i];
        float sig_ip1 = p.pedestal - amp[i+1];
        if (sig_i < target && sig_ip1 >= target) {
            float slope = (sig_ip1 - sig_i) / (time[i+1] - time[i]);
            if (slope > 0.f)
                p.crossingTime = time[i] + (target - sig_i) / slope;
            break;
        }
    }

    // Charge integral (energy proxy) over a window around the peak; positive
    // (signal) contributions only. Matches ExtractPulseMulti's definition.
    {
        int clo = imin - 15;  if (clo < 3)    clo = 3;
        int chi = imin + 200; if (chi > 1023) chi = 1023;
        double q = 0.;
        for (int s = clo; s <= chi; ++s) {
            double v = p.pedestal - amp[s];
            if (v > 0.) q += v;
        }
        p.charge = static_cast<float>(q);
    }

    return p;
}

// ---------------------------------------------------------------------------
// ExtractPulseMultiCFD
//   Same as ExtractPulse but returns CFD crossing times at four fractions:
//   10%, 20%, 30%, 50% of the peak.  crossingTime[0..3] correspond to those
//   four fractions in order.  kNoTime is stored if a crossing is not found.
//   All four fractions are extracted in a single forward scan for efficiency.
// ---------------------------------------------------------------------------
struct PulseMulti {
    float pedestal;
    float peak;
    float peakTime;
    float cfd03;   // CFD crossing at 3% of peak [ns]  (Ledovskoy sim optimum)
    float cfd05;   // CFD crossing at 5% of peak [ns]
    float cfd10;   // CFD crossing at 10% of peak [ns]
    float cfd20;   // CFD crossing at 20% of peak [ns]
    float cfd30;   // CFD crossing at 30% of peak [ns]
    float cfd50;   // CFD crossing at 50% of peak [ns]
    float ledTime; // Fixed-threshold LED crossing [ns] (threshold = ledThresh)
    float totTime; // Time over threshold [ns]  (same ledThresh used for edges)
    float charge;  // integral of pedestal-subtracted signal over the pulse window [mV*sample]
    bool  valid;
    bool  saturated;  // true if peak >= satThresh (default 950 mV) -- waveform may be clipped
    bool  spike;      // true if any pedestal sample (3..52) deviates > 5xped_RMS from mean
};

inline PulseMulti ExtractPulseMulti(const float* time, const float* amp,
                                     float ledThresh = 20.0f,
                                     float minPeak   =  5.0f,
                                     float satThresh = 950.0f)
{
    PulseMulti p;
    p.pedestal = kNoTime;   // will be overwritten
    p.peak     = 0.f;
    p.peakTime = kNoTime;
    p.cfd03    = kNoTime;
    p.cfd05    = kNoTime;
    p.cfd10    = kNoTime;
    p.cfd20    = kNoTime;
    p.cfd30    = kNoTime;
    p.cfd50    = kNoTime;
    p.ledTime   = kNoTime;
    p.totTime   = kNoTime;
    p.charge    = 0.f;
    p.valid     = false;
    p.saturated = false;
    p.spike     = false;

    // --- Pedestal: mean of samples 3..52 ---
    double ped = 0.0;
    for (int i = 3; i < 53; ++i) ped += amp[i];
    p.pedestal = static_cast<float>(ped / 50.0);

    // --- Find peak (minimum of negative-going pulse) in 3..1023 ---
    float vmin = amp[3];
    int   imin = 3;
    for (int i = 4; i < 1024; ++i) {
        if (amp[i] < vmin) { vmin = amp[i]; imin = i; }
    }
    p.peak     = p.pedestal - vmin;
    p.peakTime = time[imin];
    p.valid    = (p.peak >= minPeak);
    if (!p.valid) return p;

    // Targets for each CFD fraction and for LED
    float tgt03  = 0.03f * p.peak;
    float tgt05  = 0.05f * p.peak;
    float tgt10  = 0.10f * p.peak;
    float tgt20  = 0.20f * p.peak;
    float tgt30  = 0.30f * p.peak;
    float tgt50  = 0.50f * p.peak;
    // LED threshold is absolute (in mV above pedestal); cap at peak to avoid
    // missing pulses that don't reach a high threshold.
    float tgtLED = (ledThresh < p.peak) ? ledThresh : 0.5f * p.peak;

    // Flags — stop searching once each crossing is found
    bool f03 = false, f05 = false;
    bool f10 = false, f20 = false, f30 = false, f50 = false, fLED = false;

    // Forward scan: rising edge (sample 3 → peak)
    float tLead = kNoTime;   // LED leading-edge crossing for TOT
    for (int i = 3; i < imin; ++i) {
        float si  = p.pedestal - amp[i];
        float si1 = p.pedestal - amp[i+1];
        float dt  = time[i+1] - time[i];
        if (dt <= 0.f) continue;

        // Linear interpolation helper (inline lambda replacement for ACLiC)
        // t_cross = time[i] + (target - si) / slope
        float slope = (si1 - si) / dt;

        if (!f03  && si < tgt03  && si1 >= tgt03)  { if (slope>0.f) p.cfd03  = time[i] + (tgt03 -si)/slope; f03  = true; }
        if (!f05  && si < tgt05  && si1 >= tgt05)  { if (slope>0.f) p.cfd05  = time[i] + (tgt05 -si)/slope; f05  = true; }
        if (!f10  && si < tgt10  && si1 >= tgt10)  { if (slope>0.f) p.cfd10  = time[i] + (tgt10 -si)/slope; f10  = true; }
        if (!f20  && si < tgt20  && si1 >= tgt20)  { if (slope>0.f) p.cfd20  = time[i] + (tgt20 -si)/slope; f20  = true; }
        if (!f30  && si < tgt30  && si1 >= tgt30)  { if (slope>0.f) p.cfd30  = time[i] + (tgt30 -si)/slope; f30  = true; }
        if (!f50  && si < tgt50  && si1 >= tgt50)  { if (slope>0.f) p.cfd50  = time[i] + (tgt50 -si)/slope; f50  = true; }
        if (!fLED && si < tgtLED && si1 >= tgtLED) { if (slope>0.f) { p.ledTime = time[i] + (tgtLED-si)/slope; tLead = p.ledTime; } fLED = true; }
    }

    // --- TOT: trailing-edge crossing (falling edge, sample imin → 1023) ---
    // The signal drops back below the LED threshold after the peak.
    // We scan forward from the peak to find where sig drops below tgtLED.
    if (tLead > kNoTime + 1.f) {
        for (int i = imin; i < 1023; ++i) {
            float si  = p.pedestal - amp[i];
            float si1 = p.pedestal - amp[i+1];
            float dt  = time[i+1] - time[i];
            if (dt <= 0.f) continue;
            if (si >= tgtLED && si1 < tgtLED) {
                float slope = (si1 - si) / dt;
                if (slope < 0.f) {
                    float tTrail = time[i] + (tgtLED - si) / slope;
                    p.totTime = tTrail - tLead;
                }
                break;
            }
        }
    }

    p.saturated = (p.peak >= satThresh);

    {
        double ped_sum2 = 0.;
        for (int i = 3; i < 53; ++i) { double d = amp[i]-p.pedestal; ped_sum2 += d*d; }
        float ped_rms = static_cast<float>(std::sqrt(ped_sum2 / 50.0));
        float thresh5 = 5.0f * ped_rms;
        for (int i = 3; i < 53 && !p.spike; ++i)
            if (std::fabs(amp[i] - p.pedestal) > thresh5) p.spike = true;
    }

    // --- Charge: integral of the pedestal-subtracted signal over a window
    // around the peak (energy proxy, less shape-sensitive than peak amplitude).
    // Window [imin-15, imin+200] (~ -3 to +40 ns) captures the rise + most of
    // the LYSO decay; only positive (signal) contributions are summed.  Units
    // are mV*sample (calibratable; relative use only).
    {
        int clo = imin - 15;  if (clo < 3)    clo = 3;
        int chi = imin + 200; if (chi > 1023) chi = 1023;
        double q = 0.;
        for (int s = clo; s <= chi; ++s) {
            double v = p.pedestal - amp[s];
            if (v > 0.) q += v;
        }
        p.charge = static_cast<float>(q);
    }

    return p;
}

// ---------------------------------------------------------------------------
// Convenience: extract just the pedestal-subtracted, inverted waveform at
// sample index s.  Useful for filling waveform profiles.
// ---------------------------------------------------------------------------
inline float WaveformSample(const float* amp, float pedestal, int s) {
    return pedestal - amp[s];
}

#endif // WAVEFORMUTILS_H
