# RADiCAL Analysis for CERN May 2023 Test Beam

## Step One - Background

There were two DT5742 - 16+1 CAEN digitzers that recorded data.

These are referred to as DRS units because they use the DRS4 chip.

In the first DRS unit, there are 16 channels for detectors amd 1 channel for the trigger.

To save the data, the DRS puts it into two groups: channels 0-7 and 8-15.

But each channel has the trigger associated with it.  So in the reconstructed data you will analyze:

There are channels 0 to 8 where channel 8 is the trigger, and channels 9 through 17, where channel 17 is the trigger.

So even though the DRS outputs 17 channels, there are 18 channels in this dataset because the trigger is duplicated.

## Step Two - Reco format

To reconstruct the data, each channel's data was piggy backed on linearly.

A channel can read 1024 time samples, and each event saves 1024 time slices - waveform data - for each channel.

All of this data is stored end to end in a ROOT file in the ROOT format.

So if you want to access the waveform saved for channel 0 an event, just read the first 1024 numbers (0 to 1023) from the amplitude tree. 

If you want to access the waveform for channel 1 on the SAME event, read numbers 1024 to 2047.

    // for a given channel for a given event: drs/group/ch

    // where drs can be 0 or 1, group can be 0 or 1, and ch can be 0 to 8.

    // channel index = (1024*9*2) * drs + (1024*9) * group + 1024 * ch

    // time index = (1024*2) * drs + 1024 * group

This list gives the following indices for channel data in the event tree using that formula:

Channels:

DRS0:

0

1024

2048

3072

4096

5120

6144

7168

8192

9216

10240

11264

12288

13312

14336

15360

16384

17408

18432



DRS1:



19456

20480

21504

22528

23552

24576

25600

26624

27648

28672

29696

30720

31744

32768

33792

34816

35840

## Step Three - reading the data

To read the data, use a TTreeReader:

    TTreeReader reader("pulse",drs_file);

    TTreeReaderValue<int> run(reader, "run");

    TTreeReaderValue<int> event0(reader, "event0");

    TTreeReaderValue<int> event1(reader, "event1");

    TTreeReaderValue<int> trigger(reader, "trigger");

    TTreeReaderArray<float> timevalue(reader, "timevalue");

    TTreeReaderArray<float> amplitude(reader, "amplitude");

What we are mainly interested in is:

timevalue

and

amplitude


