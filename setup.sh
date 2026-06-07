#!/bin/bash
# ============================================================================
# setup.sh — source ONCE per shell to drive RADiCAL analysis from anywhere:
#     source setup.sh
# Sets the ROOT include path to the four lib/ domain dirs and RAD_DATA to the
# repo root, so every macro's bare #include and radReduced()/radRaw() resolve.
# (Reduction runs on Argon — see reduce/hpc/. This is for local analysis.)
# ============================================================================
_RAD_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
export ROOT_INCLUDE_PATH="$_RAD_ROOT/lib/waveform:$_RAD_ROOT/lib/io:$_RAD_ROOT/lib/physics:$_RAD_ROOT/lib/viz${ROOT_INCLUDE_PATH:+:$ROOT_INCLUDE_PATH}"
export RAD_DATA="$_RAD_ROOT"
echo "RADiCAL env ready:"
echo "  ROOT_INCLUDE_PATH = lib/{waveform,io,physics,viz}"
echo "  RAD_DATA          = $RAD_DATA"
echo "Try:  root -l -b -q 'reduce/verify.C+(\"data/2023/reduced\")'"
