
# Task: init
Implement Labelling2D:
1. read in a set of tagged traces (2D sparse image) from the input IFrame::pointer following this /exp/dune/app/users/yuhw/wct-ref/img/src/MaskSlice.cxx, denote this input as "reco"
2. read in SimChannel with a configurable artROOT label, which is essentially also a 2D sparse image, denote this input as "simchannel"
3. try to match each pixel in "reco" with "simchannel" using this function: https://internal.dunescience.org/doxygen/classsim_1_1SimChannel.html#aadcd58ce655a71ec74af4fd77b48a813, denoted the output as "matched"
4. make an output IFrame using the matched[0].TrackID for the matched pixels and -1 for the non-matched pixels, an example making an IFrame output could be larwirecell/Components/CookedFrameSource.cxx

# Task: Update Labelling2D, PID trace
1. use the trackID and ParticleInventoryService to get the PID for each matched pixel. following: /exp/dune/app/users/yuhw/larreco/larreco/WireCell/CellTree_module.cc
2. assign 0 for the non-matched pixels
3. tag the trackID traces and PID traces separately with different configurable labels.