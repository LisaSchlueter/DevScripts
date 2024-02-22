# understand simple_calibration and do some plots along the way
using LegendDataManagement
using LegendSpecFits
using HDF5, LegendHDF5IO
using TypedTables
using Plots
using Distributions # for fit()
using StatsBase # for Histogram
using RadiationSpectra # for peakfinder()
include("../../PrettyPlotArg.jl")
include("../utils.jl")

# load meta data 
l200 = LegendData(:l200)
ch_geds     = channelinfo(l200,(:p03,:r000,:cal),system = :geds, only_processable = true).channel;
ch_dets     = channelinfo(l200,(:p03,:r000,:cal),system = :geds, only_processable = true).detector;

# processing work-qround, some channels didnt make the cut
path       = l200.tier[:jlhitch,:cal,DataPeriod(3),DataRun(0)] .* "/" # get data path
filenames  = path .* readdir(path);
pattern =  r"(?<=-ch)(.*)(?=-tier_jlhit\.lh5)";  # regular expression, look for channel name in file
ch_all = [match(pattern, filename).match for filename in filenames]; # extract all channel number from filename
ch_geds_str = map(x -> string(x)[3:end],ch_geds) # convert into strings
sel = map(x -> x in ch_all, ch_geds_str)
ch_geds = ch_geds[sel]
ch_dets = ch_dets[sel]
nChannel    =  length(ch_geds)

# select channel and load data 
Idx = 1
data        = load_hitch(l200,DataPeriod(3),DataRun(0); channel=ch_geds[Idx]);
e_uncal = data.e_cusp; 

# plot uncalibrated data 
n_bins = 15000 # hardcoded in simple_calibration.jl 
h_uncal = fit(Histogram, e_uncal, nbins=n_bins) #advantage over histogram(), maybe that Plots.jl is not needed (large package)? 
h_uncal2= histogram(e_uncal,nbins=n_bins)
plot!(h_uncal,alpha=0.3)

# peakfinder. peakpos is the peak position in ADC counts. peak finder settings are hardcoded
_, peakpos = RadiationSpectra.peakfinder(h_uncal, σ=5.0, backgroundRemove=true, threshold=10)

# simple calibration: just a linear function without offset
fep_guess = sort(peakpos)[end] #full energy peak(= peak with highest energy) in ADC counts
c = 2614.5 / fep_guess
e_simple = e_uncal .* c # roughly calibrated spectrum in keV
histogram(e_simple,nbins=n_bins)


bin_window_cut = 2103.5 - 10 .< e_simple .< 2103.5 + 10
# get optimal binning for simple calibration
bin_width  = LegendSpecFits.get_friedman_diaconis_bin_width(e_simple[bin_window_cut])
# create histogram for simple calibration
h_calsimple = fit(Histogram, e_simple, 0:bin_width:3000)
# get histograms around calibration lines and peakstats
peakhists = LegendSpecFits.subhist.(Ref(h_calsimple), [(peak-first(window), peak+last(window)) for (peak, window) in zip(th228_lines, window_sizes)])
# peakhists = LegendSpecFits.subhist.([e_simple[peak-window .< e_simple .< peak+window] for (peak, window) in zip(th228_lines, window_sizes)])
peakstats = StructArray(estimate_single_peak_stats.(peakhists))

