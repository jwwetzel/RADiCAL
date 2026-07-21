# analyze/studies/INDEX.md — one line per study macro

Recorded 2026-07-21 (REORG_PLAN Phase B3) from each macro's own header. Status key:
**LOAD-BEARING** = generates a manuscript figure (only two qualify; everything else
published comes from `papers/scripts/` gates); **DEPRECATED** = in-file banner naming
its gate replacement; **RETIRED** upstream methods are flagged where a macro drives
them. All others are the preserved exploration record (`analyze/README.md`).
`[FIG]` = produces figures (into `figures/` or `output/`, both regenerable).

| Macro | Purpose | Status |
|---|---|---|
| `absShowerTime.C` | absolute shower time from reduced data, done correctly (no double MCP subtraction) | exploratory |
| `absTime.C` | absolute shower time vs MCP via leading-edge + amplitude-walk correction | exploratory |
| `adaptiveTiming.C` | per-channel adaptive timing source: clipped→lgcfd, un-clipped→led [FIG] | exploratory |
| `alignmentAnalysis.C` | data-driven transverse alignment of RADiCAL, MCP, beam in track frame [FIG] | exploratory |
| `analyzeResolution.C` | multi-energy energy/timing/hit-map resolution analysis, key physics plots [FIG] | exploratory |
| `anomaly150.C` | why DSB1 150 GeV reads worse than 125; run-to-run drift check | exploratory |
| `averageWaveforms.C` | per-channel average DRS4 waveforms aligned on CFD-20% crossing [FIG] | exploratory |
| `brightCoreScan.C` | (DW−UP)/2 σ_t vs fiducial radius tightening onto the bright core | exploratory |
| `build_chain_html.py` | assemble the chain-of-evidence HTML page from produced figures | site generator |
| `cfdFractionDSB1.C` | (DW−UP)/2 σ_t vs CFD fraction on saturated DSB1 channels | exploratory |
| `cfdInterp.C` | full-stat test of a fixed absolute timing threshold via CFD interpolation | exploratory |
| `cfdScan.C` | is CFD-5% the optimal timing threshold, or does clipping break it | exploratory |
| `channelCalibration.C` | inter-channel relative gain calibration from EM-shower LG amplitudes [FIG] | exploratory |
| `channelCombinationScan.C` | brute-force all 255 HG channel-subset timing combinations per energy | exploratory |
| `channelIntegrity.C` | Layer 1 hardware-integrity report, five-page inter-channel consistency PDF [FIG] | exploratory |
| `chargeProfiles.C` | HG amplitude TProfile2D maps for all 8 capillaries vs position [FIG] | exploratory |
| `checkCenter.C` | does RadView.beamCenter match the real ScanRunCenters | exploratory |
| `compareConfigsPlot.C` | cross-config timing & energy resolution vs energy, first-look comparison [FIG] | exploratory |
| `compareEnergies.C` | cross-energy collated plots, four multi-page summary PDFs [FIG] | exploratory |
| `configBestBin.C` | best-contained-bin timing on a reduced config dataset, first look | exploratory |
| `configBestBinDSB1.C` | best-bin (DW−UP)/2 on DSB1 processRun ntuples to reproduce headline | exploratory |
| `configBestBinHGLG.C` | best-bin (DW−UP)/2 run-folded OOS for HG vs LG chains | exploratory |
| `configCapDiag.C` | identify E-type vs T-type capillaries per build from the data | exploratory |
| `configResolution.C` | first-look timing + energy resolution for a reduced config dataset | exploratory |
| `configResolutionDSB1.C` | same simple method on DSB1 processRun ntuples for cross-config baseline | exploratory |
| `configResolutionFull.C` | per-config σ_t(E) + σ_E/E, OOS best-bin, all four builds [FIG] | exploratory (follow-on paper) |
| `configSystematics.C` | per-config systematic on OOS best-bin σ_t at 150 GeV | exploratory (follow-on paper) |
| `crossBuild.C` | full-spectrum (DW−UP)/2 timing for all four builds, de-biased [FIG] | exploratory |
| `cutScan.C` | time resolution vs cumulative amplitude threshold, per beam energy [FIG] | exploratory |
| `depthCorr.C` | depth-independent/corrected shower-time estimator to break the timing floor | exploratory |
| `desaturateCFD.C` | LG-referenced de-saturation of CFD threshold on raw DSB1 waveforms | exploratory |
| `discoverChannels.C` | infer the slot→role channel map from a raw run file | exploratory |
| `discoverReduced.C` | derive a config's slot→role map from a reduced ntuple | exploratory |
| `drs4Diagnostics.C` | DRS4/DT5742 hardware-integrity diagnostics, six-page PDF [FIG] | exploratory |
| `drs4TimeBase.C` | DT5742/DRS4 time-base verification & stop-cell correction validation [FIG] | exploratory |
| `dsb1OOSandCorr.C` | DSB1 HG best-bin OOS + LG→HG de-saturation viability check | exploratory |
| `edgeMechanism.C` | the cause behind the CFD-fraction timing trend (edge shape/jitter) [FIG] | exploratory |
| `edgeSlope.C` | what limits timing after the DT5742 HG channel clips [FIG] | exploratory |
| `elbowFractionTrend.C` | publication-quality per-channel σ_t vs CFD fraction, all six energies [FIG] | exploratory |
| `elbowInvestigation.C` | diagnose bimodal elbow in Down-capillary CFD-20% timing distributions | exploratory |
| `elbowShape.C` | dump the signed shoulder shape for SE-D at 50 GeV | exploratory |
| `etypeChar.C` | Paper-2 in-beam characterization of E-type vs T-type capillaries [FIG] | exploratory (Paper 2) |
| `etypeEnergy.C` | energy resolution from the E-type (full-length WLS) energy capillary [FIG] | exploratory |
| `fiducialTimingScan.C` | is the timing fiducial radius optimal, all energies (referee-proof) [FIG] | exploratory |
| `fixCFD.C` | does CFD vs LG-predicted true peak recover energy-timing improvement | exploratory |
| `floorErr.C` | the timing floor via free fits with proper per-point errors [FIG] | exploratory |
| `floorFit.C` | is the (DW−UP)/2 floor b physical or a fit artifact [FIG] | exploratory |
| `floorModel.C` | timing floor is model-dependent: photostat vs slew fit [FIG] | exploratory |
| `fourDdemo.C` | Paper-2 capstone: one module measuring time, energy, position simultaneously [FIG] | exploratory (Paper 2) |
| `gateDSB1.C` | read full-stat DSB1 file via canonical reader, sanity σ_t | exploratory |
| `harvestResults.C` | single source of truth for every report number → results.json | report pipeline |
| `headlineLG.C` | the headline (DW−UP)/2 best-bin ladder, cfd05 vs lgcfd [FIG] | exploratory |
| `hgLgClean.C` | off-line populations in HG_peak vs LG_peak, clean-event selection [FIG] | exploratory |
| `hgLgLinearity.C` | HG-vs-LG dual-readout linearity, basis for LG-referenced de-saturation [FIG] | exploratory |
| `hgLgPlot.C` | 8-channel HG_peak vs LG_peak panels with spike-cut + main-line fit [FIG] | **LOAD-BEARING** (paper Fig. 3) |
| `idealUniform.C` | "what if every capillary followed the trend" intrinsic-potential projection [FIG] | exploratory |
| `inferMixed.C` | infer which MIXED corners are DSB1 vs LuAG from data | exploratory |
| `investigatePbGlass.C` | PbGlass dual-band structure, four-population classification | exploratory |
| `layer1Summary.C` … `layer6Summary.C` | hero figures for report Layers 1–6 (hardware, MCP reference, beam, calibration, physics, systematics) [FIG] | report pipeline |
| `lgCFD.C` | CFD on LG-predicted true peak from raw, best-bin comparison [FIG] | exploratory |
| `lightYieldThesis.C` | Paper-1 thesis figure: brightest-1000 σ_t(E) fits, all builds [FIG] | exploratory |
| `lightYieldTiming.C` | follow-on-paper money plot: single-channel σ_t vs light yield [FIG] | exploratory (follow-on) |
| `luagDiag.C` | why LuAG is non-monotonic; diagnose how each energy point builds [FIG] | exploratory |
| `mcpJitter.C` | inter-group (mezzanine) DRS4 timing-reference jitter, not intrinsic MCP resolution | exploratory |
| `mcpJitterCanon.C` | what the "100 ps MCP jitter" really is, on canonical data | exploratory |
| `methodCompare.C` | per-build motivation for timing-estimator choice on identical events | **DEPRECATED** → `papers/scripts/method_gain_postfix/` |
| `methodDist.C` | (DW−UP)/2 distribution + Gaussian-core fit behind every timing point | exploratory |
| `mixedHeadToHead.C` | clean in-event DSB1 vs LuAG:Ce WLS comparison | **DEPRECATED** → `papers/scripts/mixed_killshot_bootstrap/` |
| `mixedLGcheck.C` | is MIXED low-gain readout dead or just dim (HG-LG linearity) | exploratory |
| `mixedSeparate.C` | measure MIXED's two materials separately via same-material diagonal corners | exploratory |
| `moduleCenter.C` | data-driven shashlik centre from calorimeter edge half-max crossings | exploratory |
| `monotonicityEvidence.C` | two supporting plots (saturation trap, outlier mechanism) for the σ_t(E) fix | exploratory |
| `monotonicityFix.C` | before/after of the σ_t(E) estimator fix, per build | exploratory |
| `narrativeFigs.C` | the two crux timing-narrative figures: clip-wall spectrum and clipped pulse [FIG] | **LOAD-BEARING** (paper Figs. 2 + 4) |
| `narrativeLadder.C` | payoff figures with fixed estimators: cfd05 vs hg_lgcfd ladder [FIG] | exploratory |
| `oosFraction.C` | fair held-out test for deterministic brightest-fraction selection | exploratory |
| `oosValidate.C` | held-out fold-by-run validation of lgcfd timing | exploratory |
| `optAppendix.C` | Paper-1 selection/optimization appendix, all four builds in one figure [FIG] | exploratory (Fig. A.10 ancestry; curve re-verification = Phase B item) |
| `optLadder.C` | remove equal-width-bin artifact via consistent brightest-K slice per energy | exploratory |
| `optScan.C` | expose every value the best-bin estimator selects, at every energy | exploratory |
| `outlierPeek.C` | is the kurtosis tail broken-timing events or physics | exploratory |
| `paperSystematics.C` | Paper-1 systematic-uncertainty budget, all four builds | **DEPRECATED** → `papers/scripts/systematics_postfix/` |
| `pbLeakage.C` | is high-E LG sub-linearity longitudinal leakage into PbGlass | exploratory |
| `peaksHGLG125.C` | reproduce legacy 125 GeV HG-vs-LG peaks plot, measure unsaturated scatter | exploratory |
| `perRun150.C` | per-run MCP/SiPM leading-edge timing for the four 150 GeV runs | exploratory |
| `positionCorrection.C` | charge-sharing position reconstruction + position-binned walk correction | exploratory |
| `pubFig.C` | two figures per build (distribution soundness + resolution/floor) [FIG] | exploratory |
| `qualityPlots.C` | data quality, channel performance, cut-optimisation diagnostics [FIG] | exploratory |
| `reductionQA.C` | 16-panel per-build QA of what the reducer extracted from raw waveforms [FIG] | exploratory |
| `satLinearity.C` | HG clips at readout rail while LG stays linear (readout-window ceiling) | exploratory |
| `showerLocalization.C` | Paper 2: reconstruct transverse shower impact from four corner amplitudes [FIG] | exploratory (Paper 2) |
| `sigmaProbe.C` | exhaustive σ_t(E) diagnostic isolating fit-vs-selection non-monotonicity cause | exploratory |
| `sigmaVsAmp.C` | is DSB1 timing slew-limited (1/amplitude)+floor; energies collapse onto one curve | exploratory |
| `sigmaVsBrightness.C` | why 150 GeV reads worse than 125 (best-bin window past minimum) | exploratory |
| `slewTest.C` | σ_t vs light yield (SumLG), all energies pooled, per config | exploratory |
| `slopeVsThresh.C` | where on the rising edge we time; does slope rise with energy | exploratory |
| `systematicUncertainties.C` | Layer 6: A2-combo timing systematic via cut variations | report pipeline |
| `tenergyClean.C` | recompute TENERGY σ_t excluding the contaminating NW energy capillary | exploratory |
| `testConfig.C` | phase-1 gate: JSON config must reproduce the hardcoded channel map | exploratory |
| `threshScan.C` | does timing higher on the edge improve σ_t and track energy | exploratory |
| `timingAllMethods.C` | brightest-1000 (DW−UP)/2 σ_t(E) for all timing methods per build | exploratory |
| `timingContainmentScan.C` | PbGlass containment-cut optimisation vs timing resolution | exploratory |
| `timingEnergyBins.C` | energy-binned timing analysis (arXiv:2401.01747 Sec 5.3) with CFD/walk extensions | exploratory |
| `timingFloorComparison.C` | the "22 vs 17.5 ps floor" (extrapolated) explained in one figure [FIG] | exploratory |
| `timingMethods.C` | compare eight timing reconstruction methods for RADiCAL HG channels | exploratory |
| `timingRegression.C` | validate the σ-monotonicity fix using production RadTiming.h functions | exploratory |
| `timingResolution.C` | combined timing resolution for five corner-combination strategies | exploratory |
| `timingTailAnalysis.C` | characterise non-Gaussian tails of the (DW−UP)/2 distribution | exploratory |
| `transverseMaps.C` | per-channel transverse (x,y) amplitude maps before/after 1×1-cm trigger cut | exploratory |
| `trSync.C` | TR0 as inter-group sync and absolute timing reference | exploratory |
| `uniformityScan.C` | spatial uniformity of A2-weighted combo timing σ_t vs position | exploratory |
| `walkCorrTest.C` | measure CFD-20% amplitude-walk curve; test data-driven correction | exploratory |
| `waveformProfiles.C` | per-channel HG+LG average DRS4 waveform profiles for any build [FIG] | exploratory |
| `wireChamberResolution.C` | data-driven delay-line wire-chamber spatial resolution without a second tracker | exploratory |

Notes: `sigmaT.C`, `timingLadder.C`, `timingHeadline.C` (top-level `analyze/`) drive the
**RETIRED** best-bin method (`rad::timingBestBin`, retirement declared at its definition
in `lib/physics/RadTiming.h` and in `sigmaT.C`'s header). The production headline path is
`rad::timingBrightestK` + `kLGCFD` via `papers/scripts/timing_fit_summary/`.
