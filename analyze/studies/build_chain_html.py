#!/usr/bin/env python3
"""build_chain_html.py — assemble the full chain-of-evidence HTML page.
Walks the figures produced by waveformProfiles / reductionQA / methodDist /
timingAllMethods / methodCompare and lays them out build-by-build, stage-by-stage:
waveforms -> reduction -> fits -> timing. Only existing files are linked.
  python3 analyze/studies/build_chain_html.py
Output: chain_of_evidence.html  (open in a browser)
"""
import os
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
FIG  = "figures/2023/narrative"
def ex(p): return os.path.exists(os.path.join(ROOT, p))

BUILDS = [
    ("DSB1",   "DSB1 organic WLS capillary (bright)",          "25–150 GeV",  "all 6 energies"),
    ("LUAG",   "LuAG:Ce ceramic WLS capillary (dim)",          "50–150 GeV",  "50 & 150 GeV only"),
    ("MIXED",  "mixed module: NE,SW = DSB1 · NW,SE = LuAG",     "50–150 GeV",  "50,75,100,125,150 GeV"),
    ("TENERGY","segmented 'Energy' module",                    "50–150 GeV",  "50 & 150 GeV only"),
]
NICE = {"DSB1":"DSB1","LUAG":"LuAG","MIXED":"MIXED","TENERGY":"TENERGY"}

def img(path, cap):
    if not ex(path): return f'<div class="missing">— not available: {cap}</div>'
    return f'<figure><a href="{path}" target="_blank"><img src="{path}" loading="lazy"></a><figcaption>{cap}</figcaption></figure>'

def stage(title, sub, figs):
    body = "\n".join(img(p, c) for p, c in figs)
    return f'<div class="stage"><h3>{title}</h3><p class="sub">{sub}</p><div class="grid">{body}</div></div>'

def build_section(b, desc, erange, wcov):
    s  = f'<section id="{b}"><div class="bhead"><h2>{NICE[b]}</h2><span class="tag">{desc}</span><span class="tag2">{erange}</span></div>'
    s += stage("1 · Waveforms — raw DRS4, mean ± RMS, energies overlaid",
               f"per-channel average pulses straight off the digitizer ({wcov}). HG shows the ~820 mV clip; LG is the slow energy pulse.",
               [(f"output/{b}/hg_waveforms.png", f"{NICE[b]} — 8 high-gain timing channels"),
                (f"output/{b}/lg_waveforms.png", f"{NICE[b]} — 8 low-gain energy channels")])
    s += stage("2 · Reduction — extracted features",
               "what the reducer pulled from every event: per-channel HG peaks (clip pile-up at 820 mV), sum_lg shower spectrum, MCP reference, beam track, fiducial counts.",
               [(f"{FIG}/reduction_{b}.png", f"{NICE[b]} — 16-panel reduction QA, energies 25(violet)→150(red)")])
    s += stage("3 · Fits — the σ extraction, every method × energy",
               "the brightest-1000 (DW−UP)/2 distribution with its Gaussian-core fit behind every timing point; G = Gaussian-core, rob = robust truncated-RMS.",
               [(f"{FIG}/method_dist_{b}.png", f"{NICE[b]} — (DW−UP)/2 distributions + fits (rows: cfd05 / led / lgcfd; cols: energies)")])
    s += stage("4 · Time resolution — σ_t(E)",
               "brightest-1000 (DW−UP)/2 σ_t(E) via the fixed production estimator; (NM)=non-monotonic. Left: all eight methods. Right: adopted-method comparison.",
               [(f"{FIG}/timing_allmethods_{b}.png", f"{NICE[b]} — all methods, a/√E ⊕ b fits"),
                (f"{FIG}/method_compare_{b}.png",    f"{NICE[b]} — adopted-method comparison")])
    s += "</section>"
    return s

def validation_section():
    s  = '<section id="VALID"><div class="bhead"><h2>Estimator validation</h2><span class="tag">σ_t(E) monotonicity fix</span></div>'
    s += stage("Before / after the estimator fix",
               "grey = old (5σ-RMS fit window, no veto); blue = fixed (robust window + in-event broken-timing veto). Adopted methods become monotonic; cfd05 headline preserved.",
               [(f"{FIG}/monotonicity_fix.png", "DSB1 lgcfd · LuAG led · MIXED lgcfd · DSB1 cfd05 (headline)")])
    s += stage("Why it was a fit/selection bug",
               "(a) ranking by HG-pulse sum selects the most-saturated events → σ rises (saturation trap, LG ranking kept). (b) one broken-timing event drove kurtosis 888→0.6 after the veto.",
               [(f"{FIG}/monotonicity_evidence.png", "saturation trap + outlier mechanism")])
    s += "</section>"
    return s

nav = " · ".join(f'<a href="#{b}">{NICE[b]}</a>' for b,_,_,_ in BUILDS) + ' · <a href="#VALID">Validation</a>'
sections = "\n".join(build_section(*b) for b in BUILDS) + "\n" + validation_section()

html = f"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">
<title>RADiCAL 2023 — Chain of Evidence</title>
<style>
 body{{font-family:-apple-system,Segoe UI,Roboto,sans-serif;margin:0;background:#0f1115;color:#e8eaed}}
 header{{padding:22px 30px;background:#161a22;border-bottom:2px solid #2a3140}}
 header h1{{margin:0;font-size:22px}} header p{{margin:6px 0 0;color:#9aa4b2;font-size:14px}}
 nav{{position:sticky;top:0;z-index:9;background:#11151c;padding:12px 30px;border-bottom:1px solid #2a3140;font-size:15px}}
 nav a{{color:#6cb6ff;text-decoration:none;margin-right:6px}} nav a:hover{{text-decoration:underline}}
 section{{padding:26px 30px;border-bottom:8px solid #0a0c10}}
 .bhead{{display:flex;align-items:baseline;gap:14px;flex-wrap:wrap;margin-bottom:6px}}
 .bhead h2{{margin:0;font-size:30px;color:#fff}}
 .tag{{background:#243049;color:#bcd;padding:3px 10px;border-radius:12px;font-size:13px}}
 .tag2{{color:#8a93a3;font-size:13px}}
 .stage{{margin:18px 0}} .stage h3{{margin:0 0 2px;font-size:17px;color:#dfe6f0}}
 .sub{{margin:0 0 10px;color:#8a93a3;font-size:13px;max-width:1100px}}
 .grid{{display:flex;flex-wrap:wrap;gap:16px}}
 figure{{margin:0;background:#161a22;border:1px solid #2a3140;border-radius:8px;padding:8px;max-width:760px}}
 figure img{{max-width:100%;height:auto;border-radius:4px;display:block}}
 figcaption{{color:#9aa4b2;font-size:12px;margin-top:6px}}
 .missing{{color:#c97;font-size:13px;font-style:italic;padding:10px}}
</style></head><body>
<header><h1>RADiCAL 2023 — Full Chain of Evidence</h1>
<p>W/LYSO:Ce shashlik EM calorimeter · raw waveforms → reduction → σ fits → time resolution, build by build. Click any figure to open full size.</p></header>
<nav>{nav}</nav>
{sections}
</body></html>"""

out = os.path.join(ROOT, "chain_of_evidence.html")
open(out, "w").write(html)
print("wrote", out)
print("builds:", ", ".join(NICE[b] for b,_,_,_ in BUILDS))
