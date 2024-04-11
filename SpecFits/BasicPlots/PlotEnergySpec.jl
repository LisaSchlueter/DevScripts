using LegendDataManagement
using HDF5, LegendHDF5IO
using TypedTables
using Plots
using Printf
include("../../PrettyPlotArg.jl")
include("../utils.jl")

path_plot = "PackageDevScripts/DevelopSpecFits/BasicPlots/plots"
path_rel_plot = relpath(path_plot,pwd())
path_abs_plot = pwd() * "/" * path_rel_plot * "/"

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
data_energy = data.e_cusp; 

# find all energies
symb_all =columnnames(data)
symb_e = filter(x -> string(x)[1]=='e',names)

# plot uncalibrated spectrum with different EnergyUncal_Binning
# on 1500bin spectrum - peakfinding is performed
energies = data.e_cusp 
plotlyjs()
histogram(energies,xlabel="Energy (ADC counts)",linestyle=:solid, ylabel="Counts",nbins=5000,label="5000 bins",fillcolor=:orange,linecolor=:orange,fillalpha=0.2)
histogram!(energies,xlabel="Energy (ADC counts)",ylabel="Counts",nbins=10000,label="10000 bins",fillcolor=:dodgerblue,linecolor=:dodgerblue,fillalpha=0.3)
histogram!(energies,xlabel="Energy (ADC counts)",linestyle=:dot,ylabel="Counts",nbins=15000,label="15000 bins",fillcolor=:indianred,linecolor=:indianred,fillalpha=0.4)
histogram!(energies,xlabel="Energy (ADC counts)",linestyle=:dash,ylabel="Counts",nbins=20000,label="20000 bins",fillcolor=:forestgreen,linecolor=:forestgreen,fillalpha=0.2,xlims=(0,16000); PlotBasic()...,PrettyLeg()...,PlotFontSize()...)
xticks!(0:5e3:1.5e4,["0","5k","10k","15k"]) #xformatter= x -> @sprintf("%.0e",x)
pname = path_abs_plot * "EnergyUncal_Binning.pdf"
savefig(pname)

# plot uncal. energy spectrum 
plotlyjs()
histogram(energies,xlabel="Energy (ADC counts)",linestyle=:solid,ylabel="Counts",nbins=15000,label="p3, r0, $(string(ch_dets[Idx])), 15,000 bins",fillcolor=:indianred,linecolor=:indianred,fillalpha=0.4,xlims=(0,16000); PlotBasic()...,PrettyLeg()...,PlotFontSize()...)
xticks!(0:5e3:1.5e4,["0","5k","10k","15k"]) #xformatter= x -> @sprintf("%.0e",x)
pname = path_abs_plot * "EnergyUncal_$(string(ch_dets[Idx])).pdf"
savefig(pname)

