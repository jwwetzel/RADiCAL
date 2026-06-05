// ============================================================================
// reduceRun.C — HPC per-run driver for the canonical config-driven reduction.
// Compiles radcore/Reducer.C and reduces ONE raw file to the canonical schema.
// (Thin wrapper so `root 'reduceRun.C+(args)'` dispatches cleanly in the SGE
//  array job; one ACLiC translation unit, prebuilt by compile.sh.)
//
//   ROOT_INCLUDE_PATH=radcore:Analysis root -l -b -q \
//     'radcore/reduceRun.C+("datasets/2023/configs/LUAG.json","<raw>.root",150,"<out>.root")'
// ============================================================================
#include "Reducer.C"

void reduceRun(const char* configPath, const char* rawFile, double energy, const char* outFile) {
    ReduceFile(configPath, rawFile, energy, outFile);
}
