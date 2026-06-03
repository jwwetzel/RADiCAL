# RADiCAL Data Lab — student files

A hands-on, no-prior-experience walkthrough for opening the CERN May-2023 test-beam
data, writing an event loop, and making the calorimeter-shadow plots.

**Start with the web walkthrough:** open `data_lab.html` at the repo root
(or live at https://jwwetzel.github.io/RADiCAL/data_lab.html). It explains every
step — installing ROOT on a Mac, getting the data from CERNBox, and running each
macro below.

## Files

| File | What it does | Run |
|---|---|---|
| `radlab.h` | Shared helpers: slot→array offsets, `findPeak`, `findPeakTime`, and the channel map | (included by the macros) |
| `drawWaveform.C` | Draw one raw waveform (one channel, one event) | `root -l 'student/drawWaveform.C(2, 344)'` |
| `peakSpectrum.C` | Loop all events, histogram one capillary's peak height | `root -l 'student/peakSpectrum.C+(2)'` |
| `beamProfile.C` | Reconstruct beam x,y from the wire chamber → 2-D beam spot | `root -l 'student/beamProfile.C+'` |
| `shadowMap.C` | Mean amplitude vs position → connector rings + calorimeter shadow | `root -l 'student/shadowMap.C+'` |

Run everything from the **top of the `RADiCAL` folder**. Each macro defaults to
`Data/RUN1258_150_GeV.root`; pass a different path as the last argument if needed.

`figs/` holds the expected-output images shown in `data_lab.html`.

## The teaching vs. production code

These macros use simple, readable peak/time finders so the logic is visible. The
real analysis uses the more careful versions in `Analysis/WaveformUtils.h`
(`ExtractPulse`, with a constant-fraction crossing time) and the full channel map
in `Analysis/ChannelConfig.h`. `Analysis/transverseMaps.C` is the production version
of `shadowMap.C`. Read those once these make sense.
