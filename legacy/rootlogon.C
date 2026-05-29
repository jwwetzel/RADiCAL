{
    // Set the default style when ROOT starts.
//    gROOT->SetStyle("Plain");
//    gStyle->SetCanvasBorderMode(0);
//    gStyle->SetPadBorderMode(0);
//    gStyle->SetPadColor(kWhite);
//    gStyle->SetCanvasColor(kWhite);
//    gStyle->SetTitleColor(kBlack);
//    gStyle->SetStatColor(kWhite);
//    gStyle->SetTitleX(0.45);
//    gStyle->SetTitleAlign(23);
//    gStyle->SetTitleBorderSize(0); 
//
//    // Set the default title of histograms
//    gStyle->SetTitleFont(62, "XYZ");
//    gStyle->SetTitleSize(0.05, "XYZ");
//    gStyle->SetTitleXOffset(0.9);
//    gStyle->SetTitleYOffset(0.9);
//
//    // Set the default font and size for axis labels
//    gStyle->SetLabelFont(62, "XYZ");
//    gStyle->SetLabelSize(0.04, "XYZ");
//
//    // Set the default font and size for axis titles
//    gStyle->SetTitleFont(62, "XYZ");
//    gStyle->SetTitleSize(0.05, "XYZ");
//
//    // Set the number of divisions on the axis
//    gStyle->SetNdivisions(505, "XYZ");
//
//    // Set the tick mark style
//    gStyle->SetPadTickX(1);
//    gStyle->SetPadTickY(1);
//
//    // Set the default line color, style, and width
//    gStyle->SetLineColor(kBlack);
//    gStyle->SetLineStyle(1);
//    gStyle->SetLineWidth(1);
//
//    // Set the default fill color for histograms
//    gStyle->SetHistFillColor(kWhite);
//
//    // Set the legend border size
//    gStyle->SetLegendBorderSize(0);
//
//    // Set the default marker style
//    gStyle->SetMarkerStyle(20);
//    gStyle->SetMarkerSize(1.0);
//    gStyle->SetMarkerColor(kBlack);
//
//    // Optimize for PostScript output
//    gStyle->SetPaperSize(TStyle::kUSLetter);
//    gStyle->SetHatchesLineWidth(5);
//    gStyle->SetHatchesSpacing(0.05);
//
//    // Turn off title and statistics box
////    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    
    gStyle->SetGridColor(kGray);
    gStyle->SetGridStyle(3);
    gStyle->SetGridWidth(1);
    gStyle->SetPadGridX(kTRUE);
    gStyle->SetPadGridY(kTRUE);
    gStyle->SetPadTickX(kTRUE);
    gStyle->SetPadTickY(kTRUE);
//
//    // Set the canvas margins
//    gStyle->SetPadTopMargin(0.08);
    gStyle->SetPadRightMargin(0.2);
//    gStyle->SetPadBottomMargin(0.1);
//    gStyle->SetPadLeftMargin(0.1);
//
//    // Set the color palette for 2D histograms
//    //gStyle->SetPalette(kBird);
    gStyle->SetPalette(kFall);
}

