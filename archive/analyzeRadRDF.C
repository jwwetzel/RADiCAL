#include "ROOT/RDataFrame.hxx"
#include "TProfile.h"
#include "TCanvas.h"
#include "TFile.h"

constexpr size_t samplesPerChannel = 1024;
constexpr size_t totalGroups = 2;
constexpr size_t totalModules = 2;
constexpr size_t totalChannelsPerGroupModule = 8; // 8 channels per group/module

void analyzeRadRDF(const char* filename)
{
    ROOT::RDataFrame df("pulse", filename);

    // Create a TProfile for each channel
    TProfile* profiles[totalGroups][totalModules][totalChannelsPerGroupModule];
    for (size_t group = 0; group < totalGroups; ++group)
    {
        for (size_t module = 0; module < totalModules; ++module)
        {
            for (size_t channel = 0; channel < totalChannelsPerGroupModule; ++channel)
            {
                TString profileName = Form("Adjusted_Profile_Group%d_Module%d_Channel%d", group, module, channel);
                profiles[group][module][channel] = new TProfile(profileName, profileName, samplesPerChannel, 0, samplesPerChannel);
            }
        }
    }

    // Process the DataFrame
    auto process = [&](unsigned int slot, const ROOT::RVec<float>& amplitude)
    {
        for (size_t group = 0; group < totalGroups; ++group)
        {
            for (size_t module = 0; module < totalModules; ++module)
            {
                for (size_t channel = 0; channel < totalChannelsPerGroupModule; ++channel)
                {
                    // Calculate indices for amplitude
                    size_t startIdxAmp = 1024 * 9 * 2 * group + 1024 * 9 * module + 1024 * channel;

                    // Calculate the pedestal for this channel
                    double pedestal = 0;
                    for (size_t i = 3; i <= 53; ++i)
                    { // Using 1-based indexing for clarity
                        pedestal += amplitude[startIdxAmp + i - 1]; // Convert to 0-based indexing
                    }
                    pedestal /= 51.0; // Average of the 3rd through 53rd values

                    // Fill the profile, subtracting the pedestal and inverting the waveform
                    for (size_t i = 0; i < samplesPerChannel; ++i)
                    {
                        double adjustedValue = -(amplitude[startIdxAmp + i] - pedestal);
                        profiles[group][module][channel]->Fill(i, adjustedValue);
                    }
                }
            }
        }
    };

    // Apply the processing function on each event
    df.ForeachSlot(process, {"amplitude"});

    // Optionally, save the profiles to a file
    TFile outFile("adjusted_channel_profiles.root", "RECREATE");
    for (size_t group = 0; group < totalGroups; ++group)
    {
        for (size_t module = 0; module < totalModules; ++module)
        {
            for (size_t channel = 0; channel < totalChannelsPerGroupModule; ++channel)
            {
                profiles[group][module][channel]->Write();
            }
        }
    }
}
