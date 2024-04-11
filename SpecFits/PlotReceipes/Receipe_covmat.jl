using LegendDataManagement
using LegendSpecFits
using LegendHDF5IO, HDF5
using TypedTables
using Plots, ColorSchemes
include("../utils.jl")

# path = "PackageDevScripts/DevelopSpecFits/BasicPlots/results"
# path_rel = relpath(path,pwd())
# path_abs = pwd() * "/" * path_rel * "/"

# path_plt = "PackageDevScripts/DevelopSpecFits/BasicPlots/plots"
# path_rel_plt = relpath(path_plt,pwd())
# path_abs_plt = pwd() * "/" * path_rel_plt * "/"

# get some data for test 
l200 = LegendData(:l200)
ch_geds     = channelinfo(l200,(:p03,:r000,:cal),system = :geds, only_processable = true).channel;
ch_dets     = channelinfo(l200,(:p03,:r000,:cal),system = :geds, only_processable = true).detector;

# do fit 
chIdx = 2
data_det        = load_hitch(l200,DataPeriod(3),DataRun(0); channel=ch_geds[chIdx]);
result_simple, report_simple, th228_lines, th228_names = do_simple_calibration([data_det.e_cusp...]); # simple calibration, get histograms around the calibration lines and peakstats
result,report = LegendSpecFits.fit_peaks_th228(result_simple.peakhists,result_simple.peakstats,th228_names) # do to: put this

v_ml = map(par -> result["Bi212a"][par],[keys(result["Bi212a"])[1:7]...])
par_names = [keys(result["Bi212a"])[1:7]...]
par_vals =  map(par -> result["Bi212a"][par],[keys(result["Bi212a"])[1:7]...])
v_ml = NamedTuple()
for i=1:length(par_names)
    v_ml = merge(v_ml, (par_names[i] = result["Bi212a"][par_names[i]]))
end
v_ml = @NamedTuple{µ, σ, n, step_amplitude, skew_fraction, skew_width, background}(par_vals[1:7])
LegendSpecFits.p_value_MC(LegendSpecFits.th228_fit_functions.f_fit, result_simple.peakhists[1], v_ml)


peaks =keys(result)
# plot
covmat = result["Tl208DEP"].gof.covmat
heatmap(covmat)

