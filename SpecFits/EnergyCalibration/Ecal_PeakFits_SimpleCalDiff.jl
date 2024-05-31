#= energy-calibration residuals plot 
plot difference between simple calibrated µ and calibrated µ for FEP 
=#
using Measures
using Plots, Printf, LaTeXStrings
using JLD2
using LegendDataManagement
using LegendDataManagement.LDMUtils
using Unitful
using Measurements
using Measurements: value as mvalue, uncertainty as muncert
using TypedTables
using StatsBase
include("$(@__DIR__)/utils_ecal.jl")

path_plot = "$(@__DIR__)/plots/p$partition/Residuals/"
if !ispath(path_plot)
    mkpath("$path_plot")
end 

#select data and dsp output 
partition = 1
e_type = :e_cusp_ctc


FitPars, MetaData = get_peakfit_rpars_partition(DataPartition(partition); reload = false, e_type = e_type, cal = true);
µ_cal = FitPars.µ
fwhm_cal = FitPars.fwhm

FitPars_simple, _ = get_peakfit_rpars_partition(DataPartition(partition); reload = false, e_type = e_type, cal = false);
µ_scal = FitPars_simple.µ
fwhm_scal = FitPars_simple.fwhm

µ_diff_FEP = mvalue.(µ_cal[:,:,7] .- µ_scal[:,:,7])
fwhm_diff_FEP = mvalue.(fwhm_cal[:,:,7] .- fwhm_scal[:,:,7])

nruns = length(MetaData.partinfo)
x_str = ["$(MetaData.partinfo[i].period)-$(MetaData.partinfo[i].run)" for i=1:nruns]

#difference in µ for FEP
plt = hline([0]*u"keV", color = :black, linestyle = :dash, linewidth = 1, label = false)
scatter!(1:nruns,µ_diff_FEP', xticks = (1:nruns, x_str), xrotation = 45, 
         xlabel = "Partition$partition: period-run",
         ylabel = "cal - simple cal",
         title = "Difference between calibrated and simple calibrated µ for FEP",
         size = (700,450), dpi = 250,
         label = false,
         markersize = 5, color = :dodgerblue, alpha = 0.2, markerstrokecolor = :auto,
         bottom_margin = 5mm,
         framestyle = :box,
         ylims = (-2.2,2.2),
         xguidefontsize = 14, yguidefontsize = 14)
scatter!([1], [NaN]* u"keV", color = :dodgerblue, markerstrokecolor = :auto, label = "Peak-fit µ of 208TlFEP",
        legend = :topleft, foreground_color_legend = :silver, background_color_legend = :white, legendfontsize = 10) 
fname = path_plot * "Ecal_PeakFit_SimpleCalDiff_part$(partition)_$(e_type).png"
savefig(plt, fname)
@info "save plot to $fname"

#difference in fwhm for FEP
plt2 = hline([0]*u"keV", color = :black, linestyle = :dash, linewidth = 1, label = false)
scatter!(1:nruns,fwhm_diff_FEP', xticks = (1:nruns, x_str), xrotation = 45, 
         xlabel = "Partition$partition: period-run",
         ylabel = "cal - simple cal",
         title = "Difference between calibrated and simple calibrated FWHM for FEP",
         size = (700,450), dpi = 250,
         label = false,
         markersize = 5, color = :darkorange, alpha = 0.2, markerstrokecolor = :auto,
         bottom_margin = 5mm,
         framestyle = :box,
         ylims = (-2.2,2.2),
         xguidefontsize = 14, yguidefontsize = 14)
scatter!([1], [NaN]* u"keV", color = :darkorange, markerstrokecolor = :auto, label = "Peak-fit FWHM of 208TlFEP",
        legend = :topleft, foreground_color_legend = :silver, background_color_legend = :white, legendfontsize = 10) 
# fname2 = path_plot * "Ecal_PeakFit_SimpleCalDiff_fwhm_part$(partition)_$(e_type).png"
# savefig(plt2, fname2)
# @info "save plot to $fname2"