#= p-value overview 
plot 2: p-value of calibration curve fit (linear/quadratic fit)

=#
# using Distributions, StatsBase, DataFrames
using LegendDataManagement
using LegendDataManagement.LDMUtils
using LegendSpecFits
using Measures
using Plots
using PDFmerger, Printf, LaTeXStrings, ColorSchemes
using TypedTables
using Unitful
using Measurements
using Measurements: value as mvalue, uncertainty as muncert
using StatsBase
include("$(@__DIR__)/utils_ecal.jl")

#settings: select data and dsp (energy) output 
partition = 2
e_type = :e_cusp_ctc

FWHMPars, MetaData = get_fwhmcurvefit_rpars_partition(DataPartition(partition); reload = false, e_type = e_type);

path_plot = "$(@__DIR__)/plots/p$partition/FWHMCurve/"
if !ispath(path_plot)
    mkdir("$path_plot")
end

pval_threshold = 0.01
fname_fwhmfit = path_plot * "Ecal_fwhmcurvefit_pLess$(pval_threshold)_$(partition)_$(e_type).pdf"

BadFits = findall(FWHMPars.pvalue .< pval_threshold)

for i in  1:1:length(BadFits)
    @info "$(MetaData.partinfo[BadFits[i][2]]), detector $(MetaData.dets_ged[BadFits[i][1]]): $(MetaData.dets_type[BadFits[i][1]])/$(BadFits[i][1])"
    fwhm_fit  = FWHMPars.fwhm_fit[BadFits[i],:]
    pp_fit = FWHMPars.peaks_literature_keV[BadFits[i],:]
    
    result_fwhm, report_fwhm = fit_fwhm(pp_fit, fwhm_fit, ; pol_order = 1,  e_expression = "Energy_ADC")
   
    residuals_abs = mvalue.(fwhm_fit .- report_fwhm.f_fit.(ustrip.(pp_fit)) .* unit.(pp_fit))
    residuals_norm = result_fwhm.gof.residuals_norm
    ylres = round(Int,maximum(abs.(residuals_norm))+2)
    if ylres > 5
        ylt = 6
    else
        ylt = 3
    end
    plttitle = "$(MetaData.partinfo[BadFits[i][2]].period)-$(MetaData.partinfo[BadFits[i][2]].run), detector $(MetaData.dets_ged[BadFits[i][1]]) ($(MetaData.dets_type[BadFits[i][1]]), $(BadFits[i][1]))\n" * "max. |residual| =  $(round(typeof(residuals_abs[1]),maximum(abs.(residuals_abs)), digits = 2)) ($(round(maximum(abs.(residuals_abs./pp_fit .* 100)),digits =2))%)"
    
    plt = plot(report_fwhm, size = (600, 370), plot_title = plttitle, plot_titlefontsize = 7, thickness_scaling = 1.2, bottom_margin = 5mm, xtickfontsize = 6, ms = 3, ylims = (-ylres, ylres), yticks = -21:ylt:21)
    savefig(plt, path_plot .* "tmp.pdf")
    append_pdf!(fname_fwhmfit, path_plot .* "tmp.pdf", cleanup = true)
end