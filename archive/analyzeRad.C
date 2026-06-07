#include "TFile.h"
#include "TH2D.h"
#include "TProfile2D.h"
#include "TCanvas.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TStyle.h"           // For gStyle
#include "TColor.h"           // For TColor and palette settings
#include "TF1.h"              // For TF1 (function fitting)
#include "TEllipse.h"         // For TEllipse
#include "TArc.h"             // For TArc
#include <fstream>            // For std::ifstream and std::ofstream
#include <iostream>           // For std::cerr and std::cout

//Variable container
struct AmplitudeTimePedestal 
{
    double peak;
    double peakTime;
    double pedestal;
    double crossingTime;
};

// This function takes the amplitude (x) and time data (y) (waveform data)
// First determines the 'pedestal' by averaging the amplitudes of the 3rd - 53rd time slices
// Then determines 'when' the threshold of the fraction of the max amplitude is reached
AmplitudeTimePedestal FindAmplitudeAndTime(float *time, float *amplitude, double fraction = 0.5, double threshold = -1)
{
    // Initialize containers
    AmplitudeTimePedestal result;
    result.peak         = 0;
    result.peakTime     = -1e+6;
    result.pedestal     = 0;
    result.crossingTime = -1e+6;

    double  ymin    = +1e+6;
    int     imin    = -1;
    double  sx      = 0;
    
    // Calculate pedestal
    for ( int i = 3; i < 53; ++i )
    { // Loop over initial 50 time slices for pedestal
        sx += amplitude[i];
    }
    result.pedestal = sx / 50.0;

    // Find the peak amplitude and its index
    for (int i = 3; i < 1024; ++i) 
    {
        if (amplitude[i] < ymin) 
        {
            ymin = amplitude[i];
            imin = i;
        }
    }

    // Calculate peak value and set peak time
    result.peak     = result.pedestal - ymin; // Absolute peak value
    result.peakTime = time[imin]; // Time of peak amplitude

    // Determine effective threshold for crossing
    double effectiveThreshold = (threshold != -1) ? threshold : result.pedestal - fraction * result.peak;
//    std::cout << effectiveThreshold << std::endl;

    // Find crossing time based on effective threshold
    for (int i = 3; i < imin; ++i) // Search for crossing from start up to peak (imin)
    {
        if (result.pedestal - amplitude[i] >= effectiveThreshold)
        {
            // Linear interpolation for more accurate crossing time
            double slope = ((result.pedestal - amplitude[i + 1]) - (result.pedestal - amplitude[i])) / (time[i + 1] - time[i]);
            result.crossingTime = time[i] + (effectiveThreshold - (result.pedestal - amplitude[i])) / slope;
            break;
        }
    }

    return result;
}




void analyzeRad(const char* filename, const char* beamEnergy)
{
    gStyle->SetPalette(kRust);
    TColor::InvertPalette();
    
    TFile *drs_file         = TFile::Open(filename);
    
    std::string paramFilename = TString("Output/" + TString(beamEnergy)+"_fit_parameters.txt").Data(); // Name of the file to save/read parameters
    
    double existingParams[8]; // Array to hold the parameters
    bool paramsExist = false;
    
    std::ifstream paramFile(paramFilename);
    if (paramFile.is_open()) {
        for (int i = 0; i < 8; ++i) {
            if (paramFile >> existingParams[i]) {
                paramsExist = true;
            } else {
                // Handle the case where there are fewer than 8 numbers in the file
                std::cerr << "Warning: Less than 8 parameters found in the file." << std::endl;
                break;
            }
        }
        paramFile.close();
    } else {
        std::cerr << "Unable to open file: " << paramFilename << std::endl;
    }
    
    
    // Check if the file is opened successfully
    if (!drs_file || drs_file->IsZombie())
    {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }
    
    TString outputFileName  = "Output/" + TString(beamEnergy) + "_Analyzed.root";
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
    TH1F*       hChPeak[totalGroups][totalChannelsPerGroup];
    TH1F*       hPeakTimes[totalGroups][totalChannelsPerGroup];
    TH2D*       hWaveProfiles[totalGroups][totalChannelsPerGroup];
    TH2D*       hHighVsLowGain[8];
    TH1F*       hTime[8];
    TH1F*       hAvgTimes[4];
    
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
    TF1  *pFits     = new TF1("linFit", "[0]*x", 0, 140); // Define linear function
    
    
    //        int radEdw0 = (1024 * 9 * 2) * 1 + (1024 * 9) * 0 + 1024 * 1; // DRS1 GRP0 CHN1 NW D
    //        int radEdw1 = (1024 * 9 * 2) * 1 + (1024 * 9) * 0 + 1024 * 2; // DRS1 GRP0 CHN2 NE D
    //        int radEdw2 = (1024 * 9 * 2) * 1 + (1024 * 9) * 0 + 1024 * 3; // DRS1 GRP0 CHN3 SE D
    //        int radEdw3 = (1024 * 9 * 2) * 1 + (1024 * 9) * 0 + 1024 * 0; // DRS1 GRP0 CHN0 SW D
    //        int radEup0 = (1024 * 9 * 2) * 1 + (1024 * 9) * 0 + 1024 * 5; // DRS1 GRP0 CHN5 NW U
    //        int radEup1 = (1024 * 9 * 2) * 1 + (1024 * 9) * 0 + 1024 * 4; // DRS1 GRP0 CHN4 NE U
    //        int radEup2 = (1024 * 9 * 2) * 1 + (1024 * 9) * 0 + 1024 * 6; // DRS1 GRP0 CHN6 SE U
    //        int radEup3 = (1024 * 9 * 2) * 1 + (1024 * 9) * 0 + 1024 * 7; // DRS1 GRP0 CHN7 SW U
    //        int itradE = (1024 * 2) * 1 + 1024 * 0;
    
    //        int radTdw0 = (1024 * 9 * 2) * 0 + (1024 * 9) * 0 + 1024 * 1; // DRS0 GRP0 CHN1 NW D
    //        int radTdw1 = (1024 * 9 * 2) * 0 + (1024 * 9) * 0 + 1024 * 2; // DRS0 GRP0 CHN2 NE D
    //        int radTdw2 = (1024 * 9 * 2) * 0 + (1024 * 9) * 0 + 1024 * 3; // DRS0 GRP0 CHN3 SE D
    //        int radTdw3 = (1024 * 9 * 2) * 0 + (1024 * 9) * 0 + 1024 * 0; // DRS0 GRP0 CHN0 SW D
    //        int itradT0 = (1024 * 2) * 0 + 1024 * 0;
    //        int radTup0 = (1024 * 9 * 2) * 0 + (1024 * 9) * 0 + 1024 * 5; // DRS0 GRP0 CHN5 NW U
    //        int radTup1 = (1024 * 9 * 2) * 0 + (1024 * 9) * 0 + 1024 * 4; // DRS0 GRP0 CHN4 NE U
    //        int radTup2 = (1024 * 9 * 2) * 0 + (1024 * 9) * 0 + 1024 * 6; // DRS0 GRP0 CHN6 SE U
    //        int radTup3 = (1024 * 9 * 2) * 0 + (1024 * 9) * 1 + 1024 * 0; // DRS0 GRP1 CHN0 SW U
    //        int itradT1 = (1024 * 2) * 0 + 1024 * 1;
    
    
    
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
        "NE UL5",       //Channel 4
        "NW UL6",       //Channel 5
        "SE UL7",       //Channel 6
        "SW UL8",       //Channel 7
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
    
    TString name    = "";
    TString title   = "";
    
    float ixlow     = 0.;
    float ixhigh    = 1024.;
    
    for ( int group = 0; group < totalGroups; ++group )
    {
        for ( int channel = 0; channel < totalChannelsPerGroup; ++channel )
        {
            int     index   = group * totalChannelsPerGroup + channel;
            
            ixhigh = (group == 0) ? 204.8 : 1024;
            
            title   = Form("Ch %d - DRS %d - %s;x Track (mm) ;y Track (mm); Mean S%d Amplitude (mV)",   channel, group, detectorNames[index].Data(), channel);
            name    = Form("h%d_DRS%d",         channel, group);
            histograms[group][channel]          = new TProfile2D(name, title, 200, -10, 30, 200, -20, 20);
            //            histograms[group][channel] = new TProfile2D(name, title, 100, -4, 16, 100, -6, 14);
            
            title  = Form("Ch %d - DRS %d - %s;Time (ns) ;Amplitude (mV)",                              channel, group, detectorNames[index].Data());
            name   = Form("prof%d_DRS%d",      channel, group);
            
            waveFormProfiles[group][channel]    = new TProfile(name, title, 1024,0,ixhigh);
            
            name   = Form("h%d_DRS%d_WF",      channel, group);
            hWaveProfiles[group][channel]       = new TH2D(name,title,1024,0,ixhigh,1200,-300,900);
            
            title   = Form("Ch %d - DRS %d - %s;Time (ns) ; # of Events",                               channel, group, detectorNames[index].Data());
            name   = Form("peakTimes_f%d_DRS%d",      channel, group);
            hPeakTimes[group][channel]          = new TH1F(name, title, 1024,0,ixhigh);
            
            title  = Form("Ch %d - DRS %d - %s;Amplitude (mV); # of Events;",                           channel, group, detectorNames[index].Data());
            name   = Form("h%d_DRS%d_peak",    channel, group);
            hChPeak[group][channel]             = new TH1F(name,title,900,-50,850);
        }
    }
    
    const char* titles[8] =
    {
        "NW DH vs NW DL",
        "NE DH vs NE DL",
        "SE DH vs SE DL",
        "SW DH vs SW DL",
        "NW UH vs NW UL",
        "NE UH vs NE UL",
        "SE UH vs SE UL",
        "SW UH vs SW UL"
    };
    
    for (int i = 0; i < 8; ++i)
    {
        TString histogramName = titles[i];
        TString histogramTitle = histogramName + "; Peak Low Gain; Peak High Gain";
        hHighVsLowGain[i]   = new TH2D(Form("hHighVsLowGain_%d", i), histogramTitle, 800, 0, 800, 1000, 0, 1000);
        
        hTime[i]            = new TH1F(Form("timing_%d",i), Form("Ch %d Timing; T_{%d} - T_{MCP};# Events;",i,i),200,25.0,27.0);
    }
    
    // Assuming 'hAvgTimes' is an array of TH1F*
    const char* titlesTimes[] =
    {
        "NW U - NW D",
        "NE U - NE D",
        "SE U - SE D",
        "SW U - SW D"
    };

    for (int i = 0; i < 4; ++i) 
    {
        hAvgTimes[i] = new TH1F(Form("hAvgTimes_%d", i), titlesTimes[i], /* number of bins */ 1000, /* x-axis min */ -2.0, /* x-axis max */ 2.0);
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
        
        AmplitudeTimePedestal MCP1 = FindAmplitudeAndTime( &timevalue[it00], &amplitude[ich[0][7]]);
        AmplitudeTimePedestal MCP2 = FindAmplitudeAndTime( &timevalue[it01], &amplitude[ich[0][16]]);
//        isOK = isOK && (sqrt(pow(x_trk - 6.5, 2) + pow(y_trk - 4.5, 2)) <= 2.0);
//        isOK = isOK && MCP1.peakTime > 80 && MCP1.peakTime < 100;
//        isOK = isOK && MCP1.peak > 100;
        
        int time = 0;
        
        if(isOK)// && (FindAmplitudeAndTime( &timevalue[it01], &amplitude[ich[0][15]] ).peak < 100.))
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
                    hPeakTimes[group][channel]->Fill( reco[group][channel].peakTime );
                    
                    for (int timeSlice = 0; timeSlice < 1024; ++timeSlice)
                    {
                        waveFormProfiles[group][channel]->Fill(timevalue[time+timeSlice], reco[group][channel].pedestal - amplitude[ich[group][channel]+timeSlice]);
                        hWaveProfiles[group][channel]->Fill(timevalue[time+timeSlice], -1*(amplitude[ich[group][channel]+timeSlice] - reco[group][channel].pedestal));
                    }
                }
            }
            
            int highGainMappings[8] = {1, 2, 3, 0, 5, 4, 6, 9};
            int lowGainMappings[8] = {1, 2, 3, 0, 5, 4, 6, 7};
            double times[8] = {0.0};
            double crossingTimes[8] = {0.0};
            AmplitudeTimePedestal MCP;
            for ( int channel = 0; channel != 8; ++channel)
            {
                time = (highGainMappings[channel] != 9) ? it00 : it01;
                MCP = (highGainMappings[channel] != 9) ? MCP1 : MCP2;
//                if ( reco[1][lowGainMappings[channel]].peak > 20 && reco[0][highGainMappings[channel]].peak > 20 && reco[0][7].peakTime > 80 && reco[0][7].peakTime < 100 && reco[0][7].peak < 750)
                {
                    hHighVsLowGain[channel]->Fill(reco[1][lowGainMappings[channel]].peak, reco[0][highGainMappings[channel]].peak);
                    times[channel] = FindAmplitudeAndTime(&timevalue[time], &amplitude[ich[0][highGainMappings[channel]]], 0.5, 0.2*reco[1][lowGainMappings[channel]].peak*existingParams[channel]).crossingTime;
                    if (highGainMappings[channel] != 9)
                    {
                        crossingTimes[channel] = FindAmplitudeAndTime(&timevalue[time], &amplitude[ich[0][7]],0.5,0.2*MCP.peak).crossingTime - times[channel];
                        hTime[channel]->Fill(FindAmplitudeAndTime(&timevalue[time], &amplitude[ich[0][7]],0.5,0.2*MCP.peak).crossingTime - times[channel]);
                    }
                    else
                    {
                        crossingTimes[channel] = FindAmplitudeAndTime(&timevalue[time], &amplitude[ich[0][16]],0.5,0.2*MCP.peak).crossingTime - times[channel];
                        hTime[channel]->Fill(FindAmplitudeAndTime(&timevalue[time], &amplitude[ich[0][16]],0.5,0.2*MCP.peak).crossingTime - times[channel]);
                    }
                }
            }
            

            if (crossingTimes[0] !=0 && crossingTimes[1] !=0 && crossingTimes[2] !=0 && crossingTimes[3] !=0 && crossingTimes[4] !=0 && crossingTimes[5] !=0 && crossingTimes[6] !=0 && crossingTimes[7] !=0)
            {
                
                hAvgTimes[0]->Fill(((crossingTimes[0]+crossingTimes[1]+crossingTimes[2]+crossingTimes[3])/4-(crossingTimes[4]+crossingTimes[5]+crossingTimes[6]+crossingTimes[7])/4)/2);
                hAvgTimes[1]->Fill(((crossingTimes[0]+crossingTimes[1]+crossingTimes[2]+crossingTimes[3])/4-(crossingTimes[4]+crossingTimes[5]+crossingTimes[6]+crossingTimes[7])/4)/2);
                hAvgTimes[2]->Fill(((crossingTimes[0]+crossingTimes[1]+crossingTimes[2]+crossingTimes[3])/4-(crossingTimes[4]+crossingTimes[5]+crossingTimes[6]+crossingTimes[7])/4)/2);
                hAvgTimes[3]->Fill(((crossingTimes[0]+crossingTimes[1]+crossingTimes[2]+crossingTimes[3])/4-(crossingTimes[4]+crossingTimes[5]+crossingTimes[6]+crossingTimes[7])/4)/2);
            }
            
            
            double sumPbGlass = reco[0][10].peak + reco[0][11].peak + reco[0][12].peak + reco[0][13].peak;
            double sumLowGain = reco[1][0].peak + reco[1][1].peak + reco[1][2].peak + reco[1][3].peak + reco[1][4].peak + reco[1][5].peak + reco[1][6].peak + reco[1][7].peak;
            
            hCaloSum->Fill(sumLowGain, sumPbGlass);
            hPbGlass->Fill(sumPbGlass);
            hLGRad->Fill(sumLowGain);
            
            hMCP->Fill(reco[0][7].peakTime);
        }
        
        event++;
        if(event % 1000 == 0)
        {
            //break;
            cout << " Event " << event << endl;
        }
    }
    
    TDirectory *channelAnalysisDir = outputFile->mkdir("CHANNEL_ANALYSIS");
    channelAnalysisDir->cd();
    TCanvas *canvas = new TCanvas("canvas","canvas",1600,1500);
    canvas->Divide(2,3);
    canvas->Print("Output/" + TString(beamEnergy)+"_.pdf[");
    
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
            waveFormProfiles[group][channel]->Write();
            
            canvas->cd(3);
            gPad->SetLogy();
            hChPeak[group][channel]->Draw();
            
            canvas->cd(4);
            gPad->SetLogz();
            hWaveProfiles[group][channel]->Draw("COLZ");
            
            canvas->cd(5);
            hPeakTimes[group][channel]->Draw();
            
            canvas->Print("Output/" + TString(beamEnergy)+"_.pdf");
            canvas->Write();
        }
    }
    canvas->Print("Output/" + TString(beamEnergy)+"_.pdf]");
    
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
    
    canvas2->Print("Output/" + TString(beamEnergy)+"_morehistograms.pdf");
    canvas2->Write();
    
    TCanvas *peakCanvas = new TCanvas("peakCanvas", "peakCanvas", 1440, 720);
    peakCanvas->Divide(4, 2);
    
    std::ofstream outFile(paramFilename, std::ios::out); // Open the file to save new parameters
    
    for (int i = 1; i <= 8; ++i)
    {
        peakCanvas->cd(i);
        // Example adjustments: Increase the left margin, decrease the right margin
        gPad->SetLeftMargin(0.15);  // Increase from 0.1 to 0.15, for example
        gPad->SetRightMargin(0.05); // Decrease from 0.1 to 0.05
        
        // Keep top and bottom margins the same unless you also need to adjust these
        gPad->SetTopMargin(0.1);
        gPad->SetBottomMargin(0.1);
        
        if (hHighVsLowGain[i-1])
        {
            if (paramsExist)
            {
                pFits->SetParameters(0,existingParams[i]);
            }
            hHighVsLowGain[i-1]->Fit(pFits,"R");
            
            double p0 = pFits->GetParameter(0); // slope
            outFile << p0 << std::endl; // Save parameters for each fit
            
            hHighVsLowGain[i-1]->Draw();
        }
        if (i < 5)
        {
            hAvgTimes[i-1]->Write();
        }
    }
    peakCanvas->Update();
    peakCanvas->Print("Output/" + TString(beamEnergy)+"_peaks.pdf");
    
    for (int i = 1; i <= 8; ++i)
    {
        peakCanvas->cd(i);
        hTime[i-1]->Draw();
    }
    
    
    
    peakCanvas->Update();
    peakCanvas->Print("Output/" + TString(beamEnergy)+"_times.pdf");

    outFile.close(); // Close the file after saving parameters
    outputFile->Close();
    
}
