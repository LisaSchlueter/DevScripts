using LegendSpecFits
using LegendDataManagement

include("plot_th228fits.jl")
include("../utils.jl")
include("utils.jl")
# fit some data to test plot function 
l200 = LegendData(:l200)
ch_geds     = channelinfo(l200,(:p03,:r000,:cal),system = :geds, only_processable = true).channel;
ch_dets     = channelinfo(l200,(:p03,:r000,:cal),system = :geds, only_processable = true).detector;
nChannel = length(ch_geds)
chIdx = 89
period = DataPeriod(3)
run = DataRun(0)
ch = ch_geds[chIdx]

data_det        = load_hitch(l200,period,run; channel=ch);
result_simple, report_simple, th228_lines, th228_names = do_simple_calibration([data_det.e_cusp...]); # simple calibration, get histograms around the calibration lines and peakstats
result,report = LegendSpecFits.fit_peaks_th228(result_simple.peakhists,result_simple.peakstats,th228_lines) # do to: put this

# test plot function 
plt_bf = plot_th228fits(result,report,th228_names ;residuals=true,fontsize = 37)
pname = _get_bestfit_pltname(period,run,ch_dets[chIdx],e_type)
savefig(plt_bf,pname)