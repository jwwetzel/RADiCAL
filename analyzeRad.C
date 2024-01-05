#include "TFile.h"
#include "TH2D.h"
#include "TProfile2D.h"
#include "TCanvas.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

//Variable container
struct AmplitudeTimePedestal {
    double peak;
    double peakTime;
    double pedestal;
};

// This function takes the amplitude (x) and time data (y) (waveform data)
// First determines the 'pedestal' by averaging the amplitudes of the 3rd - 53rd time slices
// Then determines 'when' the threshold of 50% of the max amplitude is reached
AmplitudeTimePedestal FindAmplitudeAndTime(float *time, float *amplitude, double fraction=0.5) 
{
    // Initialize containers
    AmplitudeTimePedestal result;
    result.peak         = 0;
    result.peakTime     = -1e+6;
    result.pedestal     = 0;

    double  ymin        = +1e+6;
    int     imin        = -1;
    double  sx          = 0;
    // Loop over the waveform data, skipping the first 3 noisy time slices
    for (int i = 3; i != 1024; ++i)
    {
        //find the peak (the minimum value) and save the ith is
        //pulses are negative so peak pulse is minimum value
        if (amplitude[i] < ymin)
        {
            ymin = amplitude[i];
            imin = i;
        }
        //Add up the first 50 time slices to get average pedestal
        if (i < 53)
        {
            sx += amplitude[i];
        }
    }
    
    // average the 50 time slices and subtract peak value to get the absolute peak in mV.
    result.pedestal     = sx / 50.;
    result.peak         = result.pedestal - ymin;

    double threshold = result.pedestal - fraction * result.peak;
    for (int i = imin - 2; i > 5; --i)
    {
        if (amplitude[i] >= threshold && amplitude[i + 1] < threshold) 
        {
            result.peakTime = time[i] + (threshold - amplitude[i]) / (amplitude[i + 1] - amplitude[i]) * (time[i + 1] - time[i]);
            break;
        }
    }
    
    return result;
}


void analyzeRad(const char* filename)
{
//    gStyle->SetPalette(kCividis);
    gStyle->SetPalette(kRust);
    TColor::InvertPalette();
    
    TFile *drs_file         = TFile::Open(filename);
    
    // Check if the file is opened successfully
    if (!drs_file || drs_file->IsZombie())
    {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }
    
    TString outputFileName  = TString(filename) + "_OUTPUT.root";
    TFile *outputFile       = new TFile(outputFileName, "RECREATE");
    
    TTreeReader reader("pulse",drs_file);
    TTreeReaderValue<int> run(reader, "run");
    TTreeReaderValue<int> event0(reader, "event0");
    TTreeReaderValue<int> event1(reader, "event1");
    TTreeReaderValue<int> trigger(reader, "trigger");
    TTreeReaderArray<float> timevalue(reader, "timevalue");
    TTreeReaderArray<float> amplitude(reader, "amplitude");
    
    // for a channel drs/group/ch:
    // channel index = (1024*9*2) * drs + (1024*9) * group + 1024 * ch
    // time index = (1024*2) * drs + 1024 * group
    
    
    // 1x1 cm scintillator is connected to both DRS modules
    // drs=0, group=1, ch=5
    // drs=1, group=1, ch=6
    //  int ich1 = (1024 * 9 * 2) * 0 + (1024 * 9) * 1 + 1024 * 5;
    //  int ich2 = (1024 * 9 * 2) * 1 + (1024 * 9) * 1 + 1024 * 6;
    int it00  = (1024 * 2) * 0 + 1024 * 0; // DRS 0, channels 0 to 8
    int it01  = (1024 * 2) * 0 + 1024 * 1; // DRS 0, channels 9 to 17
    int it10  = (1024 * 2) * 1 + 1024 * 0; // DRS 1, channels 0 to 8
    int it11  = (1024 * 2) * 1 + 1024 * 1; // DRS 1, channels 9 to 17
    
    const int totalGroups = 2; // Number of groups
    const int totalChannelsPerGroup = 18; // Number of channels per group
    int ich[totalGroups][totalChannelsPerGroup];
    
    for (int group = 0; group < totalGroups; ++group)
    {
        for (int channel = 0; channel < totalChannelsPerGroup; ++channel)
        {
            ich[group][channel] = (1024 * 9 * 2) * group + (1024 * 9) * 0 + 1024 * channel;
        }
    }
    
    
    // Wire Chambers channels:
    // Left:  drs=1, group=1, ch=2
    // Right: drs=1, group=1, ch=1
    // Up:    drs=1, group=1, ch=5
    // Down:  drs=1, group=1, ch=3
    int ichL = (1024 * 9 * 2) * 1 + (1024 * 9) * 1 + 1024 * 2;
    int ichR = (1024 * 9 * 2) * 1 + (1024 * 9) * 1 + 1024 * 1;
    int ichU = (1024 * 9 * 2) * 1 + (1024 * 9) * 1 + 1024 * 5;
    int ichD = (1024 * 9 * 2) * 1 + (1024 * 9) * 1 + 1024 * 3;
    int itL = (1024 * 2) * 1 + 1024 * 1;
    int itR = (1024 * 2) * 1 + 1024 * 1;
    int itU = (1024 * 2) * 1 + 1024 * 1;
    int itD = (1024 * 2) * 1 + 1024 * 1;
    
    //Book the histograms
    TProfile2D* histograms[totalGroups][totalChannelsPerGroup];
    TProfile*   waveFormProfiles[totalGroups][totalChannelsPerGroup];
    TH1F *hChPeak[totalGroups][totalChannelsPerGroup];
    TH2D *hWaveProfiles[totalGroups][totalChannelsPerGroup];
    
    TH2D *h41       = new TH2D("h41",";S1 amplitude in DRS0 (mV) ;S1 amplitude in DRS1 (mV)",
                               200,0,1000,200,0,1000);
    
    TH2D *hCaloSum  = new TH2D("hSum", "Sum of Channels;Sum of Low Gain Channels;Sum of Pb Glass Channels",
                               100, 0, 5000, 100, 0, 5000);
    
    TH1F *hPbGlass  = new TH1F("hPbGlass", "Sum of the 4 Pb Glass; Total Sum (mV); # of Events",
                               500,0,5000);
    TH1F *hLGRad    = new TH1F("hLGRad", "Sum of the 8 Low Gain RADiCAL; Total Sum (mV); # of Events",
                               500,0,5000);
    
    TH1F *hWCX      = new TH1F("hWCX", "Wire Chamber Track X; x Track (mm); # of Events",
                               100,-50,50);
    TH1F *hWCY      = new TH1F("hWCY", "Wire Chamber Track Y; x Track (mm); # of Events",
                               100,-50,50);
    
    TH1F *hMCP      = new TH1F("hMCP", "MCP Timing; Time (ns); # of Events",
                               1000,0,200);
    
    
    std::vector<TString> detectorNames =
    {
        //DRS 0
        "SW DH1",       //Channel 0
        "NW DH2",       //Channel 1
        "NE DH3",       //Channel 2
        "SE DH4",       //Channel 3
        "NE UH1",       //Channel 4
        "NW UH2",       //Channel 5
        "SE UH3",       //Channel 6
        "MCP",          //Channel 7
        "Trigger?",     //Channel 8
        "SW UH4",       //Channel 9
        "PB2 SW",       //Channel 10
        "PB1 NW",       //Channel 11
        "PB3 NE",       //Channel 12
        "PB4 SE",       //Channel 13
        "1x1 cm",       //Channel 14
        "2x2 cm",       //Channel 15
        "MCP",          //Channel 16
        "Trigger?",     //Channel 17
        //DRS1
        "SW DL5",       //Channel 0
        "NW DL6",       //Channel 1
        "NE DL7",       //Channel 2
        "SE DL8",       //Channel 3
        "NE DL5",       //Channel 4
        "NW DL6",       //Channel 5
        "SE DL7",       //Channel 6
        "SW DL8",       //Channel 7
        "Trigger?",     //Channel 8
        "WC-X-Anode",   //Channel 9
        "WC-X-Right",   //Channel 10
        "WC-X-Left",    //Channel 11
        "WC-Y-Down",    //Channel 12
        "WC-Y-Anode",   //Channel 13
        "WC-Y-Up",      //Channel 14
        "1x1 cm",       //Channel 15
        "Empty",        //Channel 16
        "Trigger?"      //Channel 17
    };
    
    for ( int group = 0; group < totalGroups; ++group )
    {
        for ( int channel = 0; channel < totalChannelsPerGroup; ++channel )
        {
            int     index   = group * totalChannelsPerGroup + channel;
            TString name    = Form("h%d_DRS%d",         channel, group);
            TString name2   = Form("prof%d_DRS%d",      channel, group);
            TString name3   = Form("h%d_DRS%d_peak",    channel, group);
            TString name4   = Form("h%d_DRS%d_WF",      channel, group);
            
            TString title   = Form("Ch %d - DRS %d - %s;x Track (mm) ;y Track (mm); Mean S%d Amplitude (mV)",   channel, group, detectorNames[index].Data(), channel);
            TString title2  = Form("Ch %d - DRS %d - %s;Time (ns) ;Amplitude (mV)",                             channel, group, detectorNames[index].Data());
            TString title3  = Form("Ch %d - DRS %d - %s;Amplitude (mV); # of Events;",                          channel, group, detectorNames[index].Data());
            //            histograms[group][channel] = new TProfile2D(name, title, 100, -4, 16, 100, -6, 14);
            histograms[group][channel]          = new TProfile2D(name, title, 200, -10, 30, 200, -20, 20);
            waveFormProfiles[group][channel]    = new TProfile(name2, title2, 1024,0,204.8);
            hChPeak[group][channel]             = new TH1F(name3,title3,900,-50,850);
            hWaveProfiles[group][channel]       = new TH2D(name4,title2,1024,0,204.8,1200,-300,900);
        }
    }
    
    int event = 0;
    while ( reader.Next() )
    {
        // Channel data
        AmplitudeTimePedestal reco[totalGroups][totalChannelsPerGroup];
        
        // Beam position
        AmplitudeTimePedestal recoL = FindAmplitudeAndTime( &timevalue[itL], &amplitude[ichL] );
        AmplitudeTimePedestal recoR = FindAmplitudeAndTime( &timevalue[itR], &amplitude[ichR] );
        AmplitudeTimePedestal recoU = FindAmplitudeAndTime( &timevalue[itU], &amplitude[ichU] );
        AmplitudeTimePedestal recoD = FindAmplitudeAndTime( &timevalue[itD], &amplitude[ichD] );
        
        double x_trk = -1e+6;
        double y_trk = -1e+6;
        bool isOK = true;
        
        isOK = isOK && (recoL.peak > 20.);
        isOK = isOK && (recoR.peak > 20.);
        isOK = isOK && (recoU.peak > 20.);
        isOK = isOK && (recoD.peak > 20.);
        if(isOK)
        {
            isOK = isOK && (recoL.peakTime > 0.);
            isOK = isOK && (recoR.peakTime > 0.);
            isOK = isOK && (recoU.peakTime > 0.);
            isOK = isOK && (recoD.peakTime > 0.);
        }
        if(isOK)
        {
            x_trk = 7. / 36. * (recoR.peakTime - recoL.peakTime);
            y_trk = 7. / 36. * (recoD.peakTime - recoU.peakTime);
            hWCX->Fill(x_trk);
            hWCY->Fill(y_trk);
        }
        
        if(isOK)// && (FindAmplitudeAndTime( &timevalue[it01], &amplitude[ich[0][15]] ).peak < 60.))
        {
            for (int group = 0; group < totalGroups; ++group)
            {
                for(int channel = 0; channel < totalChannelsPerGroup; ++channel)
                {
                    int time;
                    if (group == 0 && channel <= 8)
                    {
                        time = it00;
                    } else if (group == 0 && channel > 8)
                    {
                        time = it01;
                    } else if (group == 1 && channel <= 8)
                    {
                        time = it10;
                    } else if (group == 1 && channel > 8)
                    {
                        time = it11;
                    }
                    reco[group][channel] = FindAmplitudeAndTime( &timevalue[time], &amplitude[ich[group][channel]] );
                    histograms[group][channel]->Fill( x_trk, y_trk, reco[group][channel].peak);
                    hChPeak[group][channel]->Fill(reco[group][channel].peak);
                    
                    for (int timeSlice = 0; timeSlice < 1024; ++timeSlice)
                    {
                        waveFormProfiles[group][channel]->Fill(timeSlice*0.2, -1*(amplitude[ich[group][channel]+timeSlice] - reco[group][channel].pedestal));
                        hWaveProfiles[group][channel]->Fill(timeSlice*0.2, -1*(amplitude[ich[group][channel]+timeSlice] - reco[group][channel].pedestal));
                    }
                }
            }
            double sumPbGlass = reco[0][10].peak + reco[0][11].peak +
            reco[0][12].peak + reco[0][13].peak;
            
            double sumLowGain = reco[1][0].peak + reco[1][1].peak +
            reco[1][2].peak + reco[1][3].peak +
            reco[1][4].peak + reco[1][5].peak +
            reco[1][6].peak + reco[1][7].peak;
            hCaloSum->Fill(sumLowGain, sumPbGlass);
            hPbGlass->Fill(sumPbGlass);
            hLGRad->Fill(sumLowGain);
            
            hMCP->Fill(reco[0][7].peakTime);
        }
        
        event++;
        if(event % 1000 == 0)
        {
//                        break;
            cout << " Event " << event << endl;
        }
    }
    
    TDirectory *channelAnalysisDir = outputFile->mkdir("CHANNEL_ANALYSIS");
    channelAnalysisDir->cd();
    TCanvas *canvas = new TCanvas("canvas","canvas",1600,1500);
    canvas->Divide(2,2);
    canvas->Print(TString(filename)+"_.pdf[");
    
    for (int group = 0; group < totalGroups; ++group)
    {
        for (int channel = 0; channel < totalChannelsPerGroup; ++channel)
        {
            canvas->cd(1);
            histograms[group][channel]->Draw("colz");
            histograms[group][channel]->GetZaxis()->SetRangeUser(0,1000);
            histograms[group][channel]->GetZaxis()->SetTitleOffset(1.5);
            
            TEllipse *sensitiveCircle = new TEllipse(6.5, 4.5, 2.5, 2.5);
            sensitiveCircle->SetLineColor(kRed); // Set the line color to red
            sensitiveCircle->SetLineWidth(2);    // Set the line width
            sensitiveCircle->SetFillStyle(0);    // Set fill style to 0 (hollow)
            sensitiveCircle->Draw("SAME");
            
            TEllipse *photoCathode = new TEllipse(6.5, 4.5, 5.5, 5.5);
            photoCathode->SetLineColor(kRed); // Set the line color to red
            photoCathode->SetLineWidth(2);    // Set the line width
            photoCathode->SetFillStyle(0);    // Set fill style to 0 (hollow)
            photoCathode->Draw("SAME");
            
            TArc *mcpHousing = new TArc(6.5, 4.5, 45/2,-137,43);
            mcpHousing->SetLineColor(kRed); // Set the line color to red
            mcpHousing->SetLineWidth(2);    // Set the line width
            mcpHousing->SetFillStyle(0);    // Set fill style to 0 (hollow)
            mcpHousing->Draw("SAMEONLY");
            
            canvas->cd(2);
            waveFormProfiles[group][channel]->GetYaxis()->SetTitleOffset(1.5);
            waveFormProfiles[group][channel]->GetYaxis()->SetRangeUser(-200,800);
            waveFormProfiles[group][channel]->Draw();
            
            canvas->cd(3);
            gPad->SetLogy();
            hChPeak[group][channel]->Draw();
            
            canvas->cd(4);
            gPad->SetLogz();
            hWaveProfiles[group][channel]->Draw("COLZ");
            
            canvas->Print(TString(filename)+"_.pdf");
            canvas->Write();
        }
    }
    canvas->Print(TString(filename)+"_.pdf]");

    outputFile->cd();
    
    TCanvas *canvas2 = new TCanvas("canvas2", "Canvas2", 800*2, 600*2);
    canvas2->Divide(2,2);
    
    canvas2->cd(1);
    hCaloSum->Draw("colz");
    
    canvas2->cd(2);
    hPbGlass->Draw();
    gPad->SetLogy();
    
    canvas2->cd(3);
    gPad->SetLogy();
    hLGRad->Draw();
//    hWCY->Draw();
    
    canvas2->cd(4);
//    hWCX->Draw();
    hMCP->Draw();
    
    canvas2->Print(TString(filename)+"_morehistograms.pdf");
    canvas2->Write();
    
    outputFile->Close();
    
}
