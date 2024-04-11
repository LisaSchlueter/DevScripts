# peak fits of calibration data 
# test pseudo prior on skew widths
using LegendDataManagement
using LegendDataManagement.LDMUtils
using LegendSpecFits
# using HDF5, LegendHDF5IO
# using TypedTables
# using Plots
# using Measures
# using JLD2 # saving files 
# using Printf

path_plt = "PackageDevScripts/DevelopSpecFits/BasicPlots/plots"

path_rel_plt = relpath(path_plt,pwd())
path_abs_plt = pwd() * "/" * path_rel_plt * "/"

l200 = LegendData(:l200)
ch_geds     = channelinfo(l200,(:p03,:r000,:cal),system = :geds, only_processable = true).channel;
ch_dets     = channelinfo(l200,(:p03,:r000,:cal),system = :geds, only_processable = true).detector;
nChannel = length(ch_geds)

load_hitchfile(l200,DataPeriod(3),DataRun(0); channel=ch_geds[chIdx]))
chIdx = 2
data_det        = load_hitch(l200,DataPeriod(3),DataRun(0); channel=ch_geds[chIdx]);
result_simple, report_simple, th228_lines, th228_names = do_simple_calibration([data_det.e_cusp...]); # simple calibration, get histograms around the calibration lines and peakstats
result,report = LegendSpecFits.fit_peaks_th228(result_simple.peakhists,result_simple.peakstats,th228_lines) # do to: put this
plt_bf = PlotBf(result,report,chIdx)

p = plot(report[th228_lines[1]])
xl = xlims()
xlabel!("")
gof = result[th228_lines[1]].gof
pband = plot(range(xl[1],xl[2],10),zeros(10),ribbon = 3 .*ones(10),fillalpha=0.5,fillcolor=:lightgrey,linecolor=:darkgrey,legend=false)
pband = plot!(range(xl[1],xl[2],10),zeros(10),ribbon = ones(10),fillalpha=1,fillcolor=:lightgrey,linecolor=:darkgrey,legend=false)
pres = plot!(gof.bin_centers,gof.residuals_norm,size=(500,100),marker = :dot,markersize=3,color=:red,markeredgecolor = :transparent,linecolor = :transparent,legend=false)
ylims!(-5,5)
xlabel!("Energy (keV)")
ylabel!("Residuals (σ)")

# Create a custom layout where the first subplot is twice as large as the second
p_all = plot(p,pres,layout=layout,size=(600,400),top_margin=-2mm;PlotBasic()...,PlotFontSize(fontsize=10)...)
p_all2 = plot(p,pres,layout=layout,size=(600,400),top_margin=-3mm;PlotBasic()...,PlotFontSize(fontsize=10)...)

xlims!(xl)
# ---test ende 

failIdx = [] 
#for chIdx=85#nChannel
chIdx=84#84
    @info "start idx $chIdx $(ch_dets[chIdx])"
  #  if !isfile(path_abs_plt * "singlepeakfit_p3r0_det$(string(ch_dets[chIdx])).pdf")
        #try
                data_det        = load_hitch(l200,DataPeriod(3),DataRun(0); channel=ch_geds[chIdx]);
                result_simple, report_simple, th228_lines, th228_names = do_simple_calibration([data_det.e_cusp...]); # simple calibration, get histograms around the calibration lines and peakstats
                result,report = LegendSpecFits.fit_peaks_th228(result_simple.peakhists,result_simple.peakstats,th228_lines) # do to: put this
                plt_bf = PlotBf(result,report,chIdx;residuals=true)
                fname_plt = path_abs_plt * "singlepeakfit_p3r0_det$(string(ch_dets[chIdx])).pdf"
                savefig(plt_bf,fname_plt)
                fname = path_abs * "singlepeakfit_p3r0_det$(string(ch_dets[chIdx])).jld2"
                @save fname result plt_bf th228_lines th228_names result_simple
       # catch err
       #     @info "$err"
      #      push!(failIdx,chIdx)
       # end
  #  end
#end
