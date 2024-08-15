using LegendDataManagement
using LegendDataManagement.LDMUtils
using LegendSpecFits
using LegendHDF5IO
using TypedTables
using Unitful, Measures
using Measurements: value as mvalue, uncertainty as muncert
using StructArrays, PropDicts
using Revise
using Plots
using Distributions 
using ValueShapes
using LaTeXStrings
include("./utils_ecal.jl")
# data selection 
l200 = LegendData(:l200)
runsel = (DataPeriod(3), DataRun(4), :cal)
chinfo = channelinfo(l200, runsel; system=:geds, only_processable=true)
dets_ged  = chinfo.detector
det_types = reverse(unique(detector_type.(Ref(l200), dets_ged)))
peakshape = :f_fit#Slope

det = dets_ged[17] 
result_simple, _, th228_names, _ = simple_cal(l200, runsel, det)
result, report = fit_peaks(result_simple.peakhists, result_simple.peakstats, th228_names; e_unit=result_simple.unit, calib_type=:th228, fit_func = peakshape);

# plots 
p = plot(broadcast(k -> plot(report[k], left_margin=20mm, top_margin=-5mm, bottom_margin=-2mm, title=string(k), ms=2), th228_names)..., layout=(length(report), 1), size=(1000,710*length(report)) , thickness_scaling=1.8, titlefontsize = 10, legendfontsize = 8, yguidefontsize = 9, xguidefontsize=11)
[result[peak].skew_fraction for peak in th228_names]
[result[peak].gof.converged for peak in th228_names]
# plot(report_fit[th228_names[4]], legend = :topleft, thickness_scaling=1.1, title =  "$(runsel[1]), $(runsel[2]), $det: " * string(th228_names[4]), titlefontsize = 12, ms = 3)
# pltpath = "$(@__DIR__)/plots/specfit/"
# if !isdir(pltpath)
#     mkdir(pltpath)
# end
# # figname = pltpath *  "specfit_hypermet_vs_$(replace(string(peakshape),"f_fit_" => ""))_$(runsel[1])_$(runsel[2])_$(det_type)_$(det).png"
# # savefig(figname)
 

runperiods = [search_disk(DataRun,l200.tier[:raw, :cal, DataPeriod(period)]) for period in [3,4]]
#= fits with skew_frac above 1.0 --> [det, run, peak]
run 1 --> period 3, run 0 
CartesianIndex(17, 5, 1)
CartesianIndex(87, 9, 1)
CartesianIndex(83, 3, 2)
CartesianIndex(71, 4, 2)
CartesianIndex(10, 6, 2)
CartesianIndex(26, 10, 3)
CartesianIndex(11, 1, 4)
CartesianIndex(44, 2, 4)
CartesianIndex(75, 3, 4)
CartesianIndex(88, 7, 4)
CartesianIndex(14, 9, 4)
CartesianIndex(18, 9, 4)
CartesianIndex(11, 10, 4)
CartesianIndex(13, 4, 5)
CartesianIndex(11, 1, 6)
CartesianIndex(61, 2, 6)
CartesianIndex(13, 3, 6)
CartesianIndex(54, 11, 6)
CartesianIndex(13, 5, 7)
CartesianIndex(74, 7, 7)
CartesianIndex(59, 9, 7)
=# 

