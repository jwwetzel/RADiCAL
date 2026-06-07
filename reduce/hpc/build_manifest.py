#!/usr/bin/env python3
# ============================================================================
# build_manifest.py — parse the CERN May-2023 run logbook into an HPC manifest.
#
# Produces a CSV of physics runs (one row per run) with the columns the Argon
# job scripts need:  run, energy_GeV, label, v_up, v_down, beam, capillary,
# status.  Run it locally (the logbook lives on your laptop), commit the CSV,
# then the Argon-side discover_tasklist.sh joins it against the actual files.
#
# Default selection (the validated DSB1 electron analysis):
#   * capillary == "DSB1"
#   * beam      == "<E> GeV e-"        (electrons only; pions excluded)
#   * status    == "DONE"
#   * bias      == NOMINAL_BIAS (42.25 V)  — drops the 125 GeV 41.25 V bias scan
#                  so every energy shares one SiPM gain.  Pass --all-bias to keep
#                  every bias.
#
# Usage:
#   python3 build_manifest.py LOGBOOK.csv -o manifest_dsb1.csv
#   python3 build_manifest.py LOGBOOK.csv --capillary DSB1 --all-bias
# ============================================================================
import csv, re, argparse, collections, sys

NOMINAL_BIAS = "42.25"

def is_int(s): return bool(re.fullmatch(r'\d+', s.strip()))
def parse_energy(beam):
    m = re.search(r'(\d+)\s*GeV', beam); return int(m.group(1)) if m else None
def is_electron(beam): return bool(re.search(r'\be-?\b|electron', beam, re.I)) and not re.search(r'pi|\+|muon|mu-', beam, re.I)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("logbook")
    ap.add_argument("-o", "--out", default="manifest_dsb1.csv")
    ap.add_argument("--capillary", default="DSB1")
    ap.add_argument("--all-bias", action="store_true")
    ap.add_argument("--status", default="DONE")
    a = ap.parse_args()

    rows = list(csv.reader(open(a.logbook, newline='')))
    out = []
    for r in rows:
        if len(r) < 10: r = r + ['']*(10-len(r))
        if not is_int(r[1]): continue
        run, status, beam, vup, vdn, cap = (r[1].strip(), r[2].strip(), r[5].strip(),
                                            r[6].strip(), r[7].strip(), r[9].strip())
        e = parse_energy(beam)
        if e is None: continue                          # not a physics run
        if a.capillary and cap != a.capillary: continue
        if not is_electron(beam): continue              # electrons only
        if a.status and status.upper() != a.status.upper(): continue
        if not a.all_bias and vup and vup != NOMINAL_BIAS: continue
        out.append(dict(run=int(run), energy_GeV=e, label=f"{e}GeV",
                        v_up=vup, v_down=vdn, beam=beam, capillary=cap, status=status))

    out.sort(key=lambda x:(x['energy_GeV'], x['run']))
    with open(a.out, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=["run","energy_GeV","label","v_up","v_down","beam","capillary","status"])
        w.writeheader(); w.writerows(out)

    # summary to stderr
    by = collections.Counter(x['label'] for x in out)
    print(f"wrote {a.out}: {len(out)} runs  (capillary={a.capillary!r}, "
          f"bias={'ALL' if a.all_bias else NOMINAL_BIAS+'V'}, e- only, status={a.status})", file=sys.stderr)
    for lab in sorted(by, key=lambda l:int(l[:-3])):
        runs = sorted(x['run'] for x in out if x['label']==lab)
        print(f"  {lab:7s} {by[lab]:3d} runs  [{runs[0]}..{runs[-1]}]", file=sys.stderr)

if __name__ == "__main__":
    main()
