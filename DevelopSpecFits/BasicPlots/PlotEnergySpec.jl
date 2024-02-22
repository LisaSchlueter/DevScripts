using LegendDataManagement
using HDF5, LegendHDF5IO
using TypedTables
using Plots
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


energies = data.e_cusp 
plotlyjs()
histogram(energies,xlabel="Energy (ADC counts)",linestyle=:solid, ylabel="Counts",nbins=500,label="500 bins")
histogram!(energies,xlabel="Energy (ADC counts)",linestyle=:dot,ylabel="Counts",nbins=1000,label="1000 bins")
histogram!(energies,xlabel="Energy (ADC counts)",linestyle=:solid,ylabel="Counts",nbins=1500,label="1500 bins")
histogram!(energies,xlabel="Energy (ADC counts)",linestyle=:dash,ylabel="Counts",nbins=21000,label="2000 bins",xlims=(0,2000); PlotBasic()...,PrettyLeg()...,PlotFontSize()...)
pname = path_abs_plot * "EnergyUncal_Binning.pdf"
savefig(pname)
