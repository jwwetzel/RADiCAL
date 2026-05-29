// ============================================================================
// DRS4Calibration.h — DT5742 / DRS4 time-base recovery and stop-cell timing
//                     correction for the RADiCAL test-beam analysis
// ============================================================================
//
// WHY THIS EXISTS
// ---------------
// The May-2023 raw data was written with a *nominal* DRS4 time axis: every
// sampling cell is assigned the same width (1/5 GS/s = 0.200 ns).  The true
// DRS4 cell widths vary at the ~5-15% level, and the chip's per-cell width
// calibration was NOT applied (verified: cell-width RMS across the buffer is
// 0.84 ps == float32 precision on 0.2 ns, i.e. literally i*0.2).
//
// The hardware trigger is asynchronous to the free-running domino wave, so the
// DRS4 stop cell rotates uniformly over [0,1023] (verified: stop-cell RMS 298
// ≈ uniform 1024/sqrt(12)).  The stop cell is recoverable from the data: the
// nominal axis carries a single zero-width step at the stop cell.
//
// Because the beam signal lands at a fixed readout position but the stop cell
// rotates, a given signal edge samples different *physical* cells event to
// event.  The accumulated (uncalibrated) cell-width error from the stop cell to
// the crossing cell is therefore a reproducible function of the stop cell — a
// fixed pattern that can be measured from data and subtracted.  That is the
// purpose of StopCellCorrection below.
//
// IMPORTANT CANCELLATIONS (measured, not assumed)
// -----------------------------------------------
//   * MCP1-MCP2: both on DRS0, near-identical stop cell + readout position →
//     cell-width error cancels.  σ(MCP1-MCP2) is genuine reference jitter.
//   * (DW-UP)/2 corner estimator: both channels on DRS0 G0, same stop cell and
//     near-identical crossing cell → cell-width error largely cancels.
//   * t_HG - t_MCP (per-channel, and MCP-referenced combinations): HG crossing
//     sits ~130 cells from the MCP crossing → the error does NOT cancel and is
//     correctable here.
//
// DESIGN
// ------
// All functions/methods are inline so the header is ODR-safe when included by
// multiple ACLiC-compiled translation units (each gets its own copy), matching
// the convention used by PlotUtils.h / RADiCALStyle.h.
// ============================================================================

#ifndef DRS4CALIBRATION_H
#define DRS4CALIBRATION_H

#include <cmath>
#include <vector>

namespace drs4 {

// ---------------------------------------------------------------------------
// Number of physical sampling cells in one DT5742 / DRS4 group.
// ---------------------------------------------------------------------------
static const int kNCells = 1024;

// ---------------------------------------------------------------------------
// FindStopCell — recover the DRS4 stop (trigger) cell from a group time axis.
//
// The nominal axis is monotonically increasing at 0.200 ns/cell EXCEPT for a
// single zero/negative step at the stop cell (the readout wrap marker).  The
// stop cell is the index i minimising t[i+1] - t[i].
//
//   timeAxis : pointer to the 1024-sample time axis of one DRS group [ns]
//   returns  : stop-cell index in [0, 1022], or -1 if the axis is unusable
// ---------------------------------------------------------------------------
inline int FindStopCell(const float* timeAxis)
{
    if (!timeAxis) return -1;

    int    iMin  = -1;
    double dtMin = 1e30;
    for (int i = 0; i < kNCells - 1; ++i) {
        const double dt = static_cast<double>(timeAxis[i + 1]) - timeAxis[i];
        if (dt < dtMin) { dtMin = dt; iMin = i; }
    }
    return iMin;
}

// ---------------------------------------------------------------------------
// CellWidthRMS — RMS of the per-cell width [ps], excluding the stop-cell step.
//
// A genuinely cell-width-calibrated axis shows tens of ps of variation here.
// A nominal (uncalibrated) axis returns ~float32 precision (< 1 ps), which is
// the diagnostic signature that no calibration was applied.
// ---------------------------------------------------------------------------
inline double CellWidthRMS(const float* timeAxis, int stopCell)
{
    if (!timeAxis) return -1.;

    double sum = 0., sum2 = 0.;
    int    n   = 0;
    for (int i = 0; i < kNCells - 1; ++i) {
        if (i == stopCell) continue;             // skip the wrap discontinuity
        const double dt = static_cast<double>(timeAxis[i + 1]) - timeAxis[i];
        sum += dt; sum2 += dt * dt; ++n;
    }
    if (n < 2) return -1.;
    const double mean = sum / n;
    const double var  = sum2 / n - mean * mean;
    return (var > 0. ? std::sqrt(var) : 0.) * 1000.;   // ns → ps
}

// ---------------------------------------------------------------------------
// StopCellCorrection — reproducible, stop-cell-indexed timing offset.
//
// Trains a correction table: the mean timing residual in each stop-cell bin.
// Subtracting Offset(stopCell) from a raw time removes the fixed cell-width
// pattern while preserving the global mean (so absolute timing is unchanged).
//
// One instance per channel (each channel's crossing samples a different
// physical cell for a given stop cell, hence a different pattern).
//
// Usage:
//   StopCellCorrection corr(64);
//   for (events in TRAINING set) corr.Accumulate(stopCell, hg_cfd);
//   corr.Finalize();
//   double t_corrected = hg_cfd - corr.Offset(stopCell);   // for ANY event
// ---------------------------------------------------------------------------
class StopCellCorrection {
public:
    explicit StopCellCorrection(int nBins = 64)
        : nBins_(nBins > 0 ? nBins : 64),
          ready_(false),
          globalMean_(0.),
          sum_(nBins_, 0.),
          cnt_(nBins_, 0),
          offset_(nBins_, 0.) {}

    // Add one training sample.  stopCell in [0,1023]; residual is the raw
    // timing quantity (e.g. hg_cfd) whose stop-cell dependence we remove.
    inline void Accumulate(int stopCell, double residual)
    {
        const int b = Bin(stopCell);
        if (b < 0) return;
        sum_[b] += residual;
        cnt_[b] += 1;
    }

    // Compute the per-bin offsets after all training samples are added.
    // Bins with fewer than minCount entries fall back to zero offset (no
    // correction) so a sparsely-populated bin cannot inject noise.
    inline void Finalize(long minCount = 30)
    {
        double gSum = 0.; long gCnt = 0;
        for (int b = 0; b < nBins_; ++b) { gSum += sum_[b]; gCnt += cnt_[b]; }
        globalMean_ = (gCnt > 0) ? gSum / gCnt : 0.;

        for (int b = 0; b < nBins_; ++b) {
            if (cnt_[b] >= minCount)
                offset_[b] = sum_[b] / cnt_[b] - globalMean_;
            else
                offset_[b] = 0.;
        }
        ready_ = true;
    }

    // The correction to SUBTRACT from a raw time for this stop cell.
    inline double Offset(int stopCell) const
    {
        if (!ready_) return 0.;
        const int b = Bin(stopCell);
        return (b < 0) ? 0. : offset_[b];
    }

    // Apply the correction to a raw time (convenience).
    inline double Apply(int stopCell, double rawTime) const
    {
        return rawTime - Offset(stopCell);
    }

    bool   Ready()      const { return ready_; }
    int    NBins()      const { return nBins_; }
    double GlobalMean() const { return globalMean_; }

    // RMS of the per-bin offsets [same units as residual] — the size of the
    // reproducible cell-width pattern this correction removes.
    inline double OffsetRMS() const
    {
        if (!ready_) return 0.;
        double s = 0., s2 = 0.; int n = 0;
        for (int b = 0; b < nBins_; ++b) {
            if (cnt_[b] <= 0) continue;
            s += offset_[b]; s2 += offset_[b] * offset_[b]; ++n;
        }
        if (n < 2) return 0.;
        const double m = s / n;
        const double v = s2 / n - m * m;
        return v > 0. ? std::sqrt(v) : 0.;
    }

private:
    inline int Bin(int stopCell) const
    {
        if (stopCell < 0 || stopCell >= kNCells) return -1;
        int b = stopCell * nBins_ / kNCells;
        if (b < 0) b = 0;
        if (b >= nBins_) b = nBins_ - 1;
        return b;
    }

    int                 nBins_;
    bool                ready_;
    double              globalMean_;
    std::vector<double> sum_;
    std::vector<long>   cnt_;
    std::vector<double> offset_;
};

} // namespace drs4

#endif // DRS4CALIBRATION_H
