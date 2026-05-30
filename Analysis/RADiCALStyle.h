// ============================================================================
// RADiCALStyle.h — unified ROOT plot style for the RADiCAL test-beam analysis
// ============================================================================
//
// Call ApplyRADiCALStyle() once at the start of each macro (before booking
// any histograms) to install the global gStyle settings and register the
// custom color constants declared below.
//
// ── Available helpers ────────────────────────────────────────────────────────
//
//   Color constants (Int_t, initialised by ApplyRADiCALStyle):
//     kRData              — CERN dark navy; default for single-series plots
//     kRBlue, kRRed, kRGreen, kROrange, kRPurple, kRTeal, kRIndigo
//                         — Apple iOS system-color palette
//     kREnergyCols[6]     — cold→hot ramp: Indigo/Blue/Teal/Green/Orange/Red
//                           maps to 25/50/75/100/125/150 GeV
//     kRChannelCols[8]    — eight visually distinct channel colors
//
//   Pad helpers:
//     StylePad(rightPalette, sidebar)
//         Standard margin / tick / grid setup for a sub-pad.
//         rightPalette=true : wide right margin for a COLZ colour bar (0.18)
//         sidebar=true      : extra-wide right margin for a legend sidebar (0.27)
//                             (ignored when rightPalette is true)
//
//     SetupSidebar(canvas, plotPad, legPad, legendFrac=0.24)
//         Splits a fresh single-panel canvas into a plot TPad and a legend
//         TPad.  Preferred for summary / overlay plots with multi-entry legends.
//
//     MakeSidebarLegend(legPad, nEntries=6)
//         Returns a pre-styled TLegend sized for nEntries inside legPad.
//         The caller adds entries and calls leg->Draw().
//
//     MakeLegend(nEntries=6)
//         Returns a TLegend placed in the right-margin sidebar of the *current*
//         pad.  Requires StylePad(false, true) to have been called first.
//
//   Title helpers:
//     DrawPadTitle(text, textSize=0.060)
//         Centred near the top of the current pad (NDC).  Automatically
//         accounts for left/right margins so the label stays over the frame
//         even with a wide sidebar margin.
//
//     DrawPageTitle(text, textSize=0.022)
//         Left-aligned at the very top of the canvas.
//         Call after canvas->cd(0).
// ============================================================================

#ifndef RADICALSTYLE_H
#define RADICALSTYLE_H

#include "TStyle.h"
#include "TColor.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TROOT.h"
#include "TString.h"

// ---------------------------------------------------------------------------
// Custom color IDs
//
// Declared static so each ACLiC-compiled translation unit gets its own copy,
// eliminating ODR conflicts when multiple macros are compiled in the same
// ROOT session.  Values are ROOT fallbacks until ApplyRADiCALStyle() is called.
// ---------------------------------------------------------------------------

/// Single-series data — CERN dark navy (#003087)
static Int_t kRData   = kBlue+2;
/// Apple iOS Blue   (#007AFF)
static Int_t kRBlue   = kBlue+1;
/// Apple iOS Red    (#FF3B30)
static Int_t kRRed    = kRed+1;
/// Apple iOS Green  (#34C759)
static Int_t kRGreen  = kGreen+2;
/// Apple iOS Orange (#FF9500)
static Int_t kROrange = kOrange+2;
/// Apple iOS Purple (#AF52DE)
static Int_t kRPurple = kViolet+1;
/// Apple iOS Teal   (#5AC8FA)
static Int_t kRTeal   = kCyan+1;
/// Apple iOS Indigo (#5856D6)
static Int_t kRIndigo = kBlue+2;

/// Per-energy color ramp (cold→hot, 25 → 150 GeV)
static Int_t kREnergyCols[6] = {
    kBlue+2, kBlue+1, kCyan+1, kGreen+2, kOrange+2, kRed+1
};

/// Per-channel palette (8 visually distinct iOS-inspired colors)
static Int_t kRChannelCols[8] = {
    kBlue+1, kRed+1, kGreen+2, kOrange+2,
    kViolet+1, kCyan+1, kMagenta+1, kYellow+3
};

// ---------------------------------------------------------------------------
// ApplyRADiCALStyle — install gStyle settings and register custom colors
//
// Call once at the top of each macro's main function, before any histogram
// booking or canvas creation.
// ---------------------------------------------------------------------------
static void ApplyRADiCALStyle()
{
    // ── Register Apple iOS + CERN colors ────────────────────────────────────
    kRData   = TColor::GetColor(0x00, 0x30, 0x87);  // CERN dark navy
    kRBlue   = TColor::GetColor(0x00, 0x7A, 0xFF);  // iOS Blue
    kRRed    = TColor::GetColor(0xFF, 0x3B, 0x30);  // iOS Red
    kRGreen  = TColor::GetColor(0x34, 0xC7, 0x59);  // iOS Green
    kROrange = TColor::GetColor(0xFF, 0x95, 0x00);  // iOS Orange
    kRPurple = TColor::GetColor(0xAF, 0x52, 0xDE);  // iOS Purple
    kRTeal   = TColor::GetColor(0x5A, 0xC8, 0xFA);  // iOS Teal
    kRIndigo = TColor::GetColor(0x58, 0x56, 0xD6);  // iOS Indigo

    kREnergyCols[0] = kRIndigo;  // 25 GeV
    kREnergyCols[1] = kRBlue;    // 50 GeV
    kREnergyCols[2] = kRTeal;    // 75 GeV
    kREnergyCols[3] = kRGreen;   // 100 GeV
    kREnergyCols[4] = kROrange;  // 125 GeV
    kREnergyCols[5] = kRRed;     // 150 GeV

    kRChannelCols[0] = kRBlue;
    kRChannelCols[1] = kRRed;
    kRChannelCols[2] = kRGreen;
    kRChannelCols[3] = kROrange;
    kRChannelCols[4] = kRPurple;
    kRChannelCols[5] = kRTeal;
    kRChannelCols[6] = kRIndigo;
    kRChannelCols[7] = TColor::GetColor(0xFF, 0xCC, 0x00);  // iOS Yellow

    // ── Global TStyle ────────────────────────────────────────────────────────
    TStyle* s = gStyle;

    // Suppress statistics / fit / title boxes — we annotate manually
    s->SetOptStat(0);
    s->SetOptFit(0);
    s->SetOptTitle(0);

    // Grid: OFF globally (Ledovskoy house style — gridlines are clutter; the
    // plot reads faster without them).  Colour/style retained for the rare pad
    // that opts back in via gPad->SetGridx(1).
    const Int_t gridColor = TColor::GetColor(0xCC, 0xCC, 0xCC);
    s->SetGridColor(gridColor);
    s->SetGridStyle(3);   // dotted
    s->SetGridWidth(1);
    s->SetPadGridX(kFALSE);
    s->SetPadGridY(kFALSE);

    // Tick marks on all four sides
    s->SetPadTickX(1);
    s->SetPadTickY(1);

    // Frame
    s->SetFrameLineWidth(1);
    s->SetFrameBorderMode(0);

    // 2D colour palette (kFall: blue → warm yellow/red)
    s->SetPalette(kFall);
    s->SetNumberContours(99);

    // Font: 42 = Helvetica (clean, sans-serif)
    const Int_t    font  = 42;
    const Double_t tsize = 0.048;

    s->SetTextFont(font);
    s->SetTextSize(tsize);

    s->SetLabelFont(font, "X");
    s->SetLabelFont(font, "Y");
    s->SetLabelFont(font, "Z");
    s->SetLabelSize(tsize * 0.90, "X");
    s->SetLabelSize(tsize * 0.90, "Y");
    s->SetLabelSize(tsize * 0.90, "Z");

    s->SetTitleFont(font, "X");
    s->SetTitleFont(font, "Y");
    s->SetTitleFont(font, "Z");
    s->SetTitleFont(font, "T");
    s->SetTitleSize(tsize,        "X");
    s->SetTitleSize(tsize,        "Y");
    s->SetTitleSize(tsize,        "Z");

    s->SetTitleOffset(1.20, "X");
    s->SetTitleOffset(1.45, "Y");
    s->SetTitleOffset(1.20, "Z");

    // Default pad margins — overridden per-pad by StylePad()
    s->SetPadLeftMargin  (0.13);
    s->SetPadBottomMargin(0.12);
    s->SetPadRightMargin (0.04);
    s->SetPadTopMargin   (0.10);

    // Canvas
    s->SetCanvasBorderMode(0);
    s->SetCanvasColor(kWhite);

    // Legend defaults
    s->SetLegendBorderSize(0);
    s->SetLegendFillColor(kWhite);
    s->SetLegendFont(font);
    s->SetLegendTextSize(0.042);

    // Marker defaults
    s->SetMarkerStyle(20);
    s->SetMarkerSize(0.9);

    // Error bar end-caps: none (cleaner look)
    s->SetEndErrorSize(0);

    gROOT->ForceStyle();
}

// ---------------------------------------------------------------------------
// StylePad — set per-pad margins, tick marks and grid
//
//   rightPalette  true  → wide right margin for a COLZ colour bar (0.18)
//   sidebar       true  → extra-wide right margin for a legend sidebar (0.27)
//                         Ignored when rightPalette is true.
// ---------------------------------------------------------------------------
static void StylePad(bool rightPalette = false, bool sidebar = false)
{
    gPad->SetLeftMargin  (0.13);
    gPad->SetBottomMargin(0.12);
    gPad->SetTopMargin   (0.10);

    if      (rightPalette) gPad->SetRightMargin(0.18);
    else if (sidebar)      gPad->SetRightMargin(0.27);
    else                   gPad->SetRightMargin(0.04);

    gPad->SetTickx(1);
    gPad->SetTicky(1);
    // Grid OFF everywhere (house style) — see ApplyRADiCALStyle.
    gPad->SetGridx(0);
    gPad->SetGridy(0);
}

// ---------------------------------------------------------------------------
// SetupSidebar — split a single-panel canvas into a plot pad and a legend pad
//
//   c           fresh TCanvas (should have just been created)
//   plotPad     out: the main plot TPad  — cd() here before drawing
//   legPad      out: the legend TPad     — cd() here before drawing legend
//   legendFrac  fraction of canvas width dedicated to the legend [0.18–0.30]
//
// Typical usage:
//   TCanvas c("c","",1200,600);
//   TPad *pPlot, *pLeg;
//   SetupSidebar(c, pPlot, pLeg);
//   pPlot->cd();  h->Draw("HIST");
//   TLegend* leg = MakeSidebarLegend(pLeg, 6);
//   leg->AddEntry(h, "150 GeV", "l");
//   pLeg->cd();  leg->Draw();
// ---------------------------------------------------------------------------
static void SetupSidebar(TCanvas& c,
                          TPad*& plotPad, TPad*& legPad,
                          float legendFrac = 0.24f)
{
    c.cd();
    const float split = 1.0f - legendFrac;

    plotPad = new TPad("rpad_plot", "", 0.f, 0.f, split, 1.f);
    plotPad->SetLeftMargin  (0.15f);
    plotPad->SetBottomMargin(0.13f);
    plotPad->SetRightMargin (0.03f);
    plotPad->SetTopMargin   (0.10f);
    plotPad->SetTickx(1);
    plotPad->SetTicky(1);
    plotPad->SetGridx(0);
    plotPad->SetGridy(0);
    plotPad->Draw();

    legPad = new TPad("rpad_leg", "", split, 0.f, 1.f, 1.f);
    legPad->SetLeftMargin  (0.05f);
    legPad->SetRightMargin (0.05f);
    legPad->SetTopMargin   (0.06f);
    legPad->SetBottomMargin(0.06f);
    legPad->SetFillStyle(0);
    legPad->SetBorderSize(0);
    legPad->Draw();
}

// ---------------------------------------------------------------------------
// MakeSidebarLegend — pre-styled TLegend for use inside a legPad
//
// legPad must be cd()'d before this call (or pass it and the function cds).
// The returned legend is sized for nEntries entries.  Call AddEntry + Draw.
// ---------------------------------------------------------------------------
static TLegend* MakeSidebarLegend(TPad* legPad, int nEntries = 6)
{
    legPad->cd();

    const float rowH = 0.11f;
    const float top  = 0.90f;
    float bot = top - static_cast<float>(nEntries) * rowH;
    if (bot < 0.04f) bot = 0.04f;

    TLegend* leg = new TLegend(0.03f, bot, 0.96f, top);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.10f);   // large — rendered inside a narrow pad
    return leg;
}

// ---------------------------------------------------------------------------
// MakeLegend — TLegend placed in the right-margin sidebar of the current pad
//
// Requires StylePad(false, true) (sidebar=true, right margin = 0.27) to have
// been called first.  The legend fills the top portion of the margin area.
// ---------------------------------------------------------------------------
static TLegend* MakeLegend(int nEntries = 6)
{
    const float rm   = static_cast<float>(gPad->GetRightMargin());
    const float x1   = 1.0f - rm + 0.015f;
    const float x2   = 0.985f;
    const float rowH = 0.075f;
    const float top  = 0.90f;
    float bot = top - static_cast<float>(nEntries) * rowH;
    if (bot < 0.08f) bot = 0.08f;

    TLegend* leg = new TLegend(x1, bot, x2, top);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.040f);
    return leg;
}

// ---------------------------------------------------------------------------
// DrawPadTitle — centred title just above the frame of the current pad
//
// The x position is the midpoint of the *frame* (accounts for left and right
// margins), so the label stays centred over the plot area even with a wide
// sidebar right margin.
// ---------------------------------------------------------------------------
static void DrawPadTitle(const char* text, float textSize = 0.060f)
{
    const float lm = static_cast<float>(gPad->GetLeftMargin());
    const float rm = static_cast<float>(gPad->GetRightMargin());
    const float xc = lm + (1.0f - lm - rm) * 0.5f;
    const float yc = 1.0f - static_cast<float>(gPad->GetTopMargin()) * 0.42f;

    TLatex ptit;
    ptit.SetNDC();
    ptit.SetTextFont(52);    // italic — Ledovskoy-style descriptive pad caption
    ptit.SetTextSize(textSize);
    ptit.SetTextColor(kGray + 3);
    ptit.SetTextAlign(22);   // centre-centre
    ptit.DrawLatex(xc, yc, text);
}

// ---------------------------------------------------------------------------
// DrawPageTitle — left-aligned title at the very top of the canvas
//
// Call after canvas->cd(0) (the canvas frame, not a sub-pad).
// ---------------------------------------------------------------------------
static void DrawPageTitle(const char* text, float textSize = 0.022f)
{
    TLatex lab;
    lab.SetNDC();
    lab.SetTextFont(42);
    lab.SetTextSize(textSize);
    lab.SetTextAlign(12);   // left-centre
    lab.DrawLatex(0.01f, 0.975f, text);   // y<0.99 so tall glyphs aren't clipped at the canvas top
}

// ===========================================================================
// SQUARE-FRAME CANVASES
//
// House rule: the data frame is SQUARE; the canvas is rectangular, made wider/
// taller than the frame by *pixel* margins that are guaranteed large enough to
// hold the axis titles and labels.  We size the margins (not the title offset)
// to fit the title — never crank SetTitleOffset to shove a title that then
// clips off the canvas edge.
//
// Margins are given in PIXELS and converted to the fractional margins ROOT
// wants, so frame_px = canvas_px − margins_px is exactly square.
// ===========================================================================

// Default pixel margins for a single square panel.  Left holds a rotated y-title
// + numeric labels; bottom holds the x-title + labels; top a little breathing
// room (or a pad caption); right a thin margin.
static const int kSqML = 108, kSqMB = 94, kSqMT = 54, kSqMR = 36;

// ROOT's batch canvas/PNG output comes out a few pixels smaller than the
// requested TCanvas size (≈ −4 px wide, −28 px tall on this build — the latter
// is the suppressed title-bar allowance).  We pad the request by these so the
// *actual* canvas, and hence the frame, comes out exactly square.
static const int kCanvasShrinkW = 4, kCanvasShrinkH = 28;

// ---------------------------------------------------------------------------
// NewSquareCanvas — one square-frame plot on a rectangular canvas.
//   frame = side of the square data frame in pixels.
// Returns the canvas; draw straight onto it (it is the pad).
// ---------------------------------------------------------------------------
static TCanvas* NewSquareCanvas(const char* name, int frame = 660,
                                int mL = kSqML, int mB = kSqMB,
                                int mT = kSqMT, int mR = kSqMR)
{
    const int W = frame + mL + mR + kCanvasShrinkW;
    const int H = frame + mT + mB + kCanvasShrinkH;
    TCanvas* c = new TCanvas(name, name, W, H);
    c->SetFillColor(kWhite);
    // Fractions are relative to the ACTUAL (shrunk) canvas, so frame = exactly
    // (frame×frame) px and the margins hold the full titles/labels.
    const double Wa = W - kCanvasShrinkW, Ha = H - kCanvasShrinkH;
    c->SetLeftMargin  (static_cast<double>(mL) / Wa);
    c->SetRightMargin (static_cast<double>(mR) / Wa);
    c->SetTopMargin   (static_cast<double>(mT) / Ha);
    c->SetBottomMargin(static_cast<double>(mB) / Ha);
    c->SetTickx(1); c->SetTicky(1);
    c->SetGridx(0); c->SetGridy(0);
    // Title offsets matched to the pixel margins (so titles land inside them).
    gStyle->SetTitleOffset(1.15, "X");
    gStyle->SetTitleOffset(1.55, "Y");
    return c;
}

// ---------------------------------------------------------------------------
// NewSquareGrid — nx×ny grid of SQUARE-frame pads, with a top band reserved for
// a page title (DrawPageTitle on c->cd(0) never overprints a pad).
//   frame = side of each square pad frame in pixels.
//   returns the grid TPad — call grid->cd(i+1), i = 0..nx*ny−1.
// Per-pad margins are smaller (titles are smaller in a small pad) but still in
// pixels, so every pad frame is exactly square and titles stay on-canvas.
// ---------------------------------------------------------------------------
static TPad* NewSquareGrid(TCanvas*& c, const char* name, int nx, int ny,
                           int frame = 300, int mL = 74, int mB = 66,
                           int mT = 46, int mR = 24, int band = 48)
{
    const int cellW = frame + mL + mR;
    const int cellH = frame + mT + mB;
    const int W = nx * cellW + kCanvasShrinkW;
    const int H = ny * cellH + band + kCanvasShrinkH;
    c = new TCanvas(name, name, W, H);
    c->SetFillColor(kWhite);
    c->cd();
    // Grid pad spans the actual (shrunk) canvas below the title band; each
    // sub-pad then comes out exactly cellW×cellH ⇒ square frames.
    const double Ha = H - kCanvasShrinkH;            // = ny*cellH + band
    const double top = static_cast<double>(ny * cellH) / Ha;
    TPad* grid = new TPad(Form("%s_grid", name), "", 0., 0., 1., top);
    grid->SetBorderMode(0);
    grid->SetFillStyle(0);
    grid->Draw();
    grid->Divide(nx, ny, 0., 0.);   // no inter-pad gap; per-pad margins separate frames
    for (int i = 1; i <= nx * ny; ++i) {
        grid->cd(i);
        gPad->SetLeftMargin  (static_cast<double>(mL) / cellW);
        gPad->SetRightMargin (static_cast<double>(mR) / cellW);
        gPad->SetTopMargin   (static_cast<double>(mT) / cellH);
        gPad->SetBottomMargin(static_cast<double>(mB) / cellH);
        gPad->SetTickx(1); gPad->SetTicky(1);
        gPad->SetGridx(0); gPad->SetGridy(0);
    }
    gStyle->SetTitleOffset(1.05, "X");
    gStyle->SetTitleOffset(1.25, "Y");
    c->cd();
    return grid;
}

// ---------------------------------------------------------------------------
// MakeCornerLegend — small, borderless, transparent legend in a corner of the
// current pad (Ledovskoy style).  corner: "tr"(default)/"tl"/"br"/"bl"
// (1st char = top/bottom, 2nd char = left/right).  Call AddEntry + Draw.
// ---------------------------------------------------------------------------
static TLegend* MakeCornerLegend(int nEntries, const char* corner = "tr",
                                 double textSize = 0.040)
{
    const double pad = 0.02, w = 0.30;
    double h = nEntries * 0.055 + 0.012;
    if (h > 0.55) h = 0.55;

    const double L = gPad->GetLeftMargin(),       R = 1.0 - gPad->GetRightMargin();
    const double B = gPad->GetBottomMargin(),      T = 1.0 - gPad->GetTopMargin();

    const bool left   = (corner[0] && corner[1] == 'l');
    const bool bottom = (corner[0] == 'b');

    double x1 = left ? (L + pad)         : (R - pad - w);
    double x2 = x1 + w;
    double y1 = bottom ? (B + pad)       : (T - pad - h);
    double y2 = y1 + h;

    TLegend* leg = new TLegend(x1, y1, x2, y2);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(textSize);
    return leg;
}

// ---------------------------------------------------------------------------
// StyleColz — one-call COLZ styling for the current pad: kFall palette, wide
// right margin for the colour bar, no grid, optional log-z.  Call before
// Draw("COLZ").  (Set the z-axis title on the histogram: h->GetZaxis()->SetTitle("events").)
// ---------------------------------------------------------------------------
static void StyleColz(bool logz = false)
{
    gPad->SetLeftMargin  (0.13);
    gPad->SetBottomMargin(0.12);
    gPad->SetTopMargin   (0.10);
    gPad->SetRightMargin (0.16);
    gPad->SetTickx(1);
    gPad->SetTicky(1);
    gPad->SetGridx(0);
    gPad->SetGridy(0);
    if (logz) gPad->SetLogz(1);
    gStyle->SetPalette(kFall);
}

// ---------------------------------------------------------------------------
// PrintClean — print a canvas to PDF/PNG with the page sized to the canvas's
// own aspect ratio.  Fixes ROOT's default behaviour of letterboxing every
// canvas into a fixed A4-proportioned page (the "half-empty canvas" bug):
// without this, a wide canvas renders into the top strip of a portrait page.
// ---------------------------------------------------------------------------
static void PrintClean(TCanvas* c, const char* path)
{
    if (!c) return;
    const double w = static_cast<double>(c->GetWw());
    const double h = static_cast<double>(c->GetWh());
    if (w > 0. && h > 0.) {
        const double longest = (w > h) ? w : h;
        const double scale   = 20.0 / longest;   // longest side → 20 cm
        gStyle->SetPaperSize(w * scale, h * scale);
    }
    c->Print(path);
}

#endif // RADICALSTYLE_H
