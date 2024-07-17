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
partition = 1
e_type = :e_cusp_ctc
cal_pol_order = 1
CalPars, MetaData = get_calcurvefit_rpars_partition(DataPartition(partition); reload = false, e_type = e_type);

path_plot = "$(@__DIR__)/plots/p$partition/CalibrationCurve/"
if !ispath(path_plot)
    mkdir("$path_plot")
end
pval_threshold = 0.01
fname_calfit = path_plot * "Ecal_calcurvefit_pLess$(pval_threshold)_$(partition)_$(e_type).pdf"

BadFits = findall(CalPars.pvalue .< pval_threshold)

for i in  1:1:25
    @info "$(MetaData.partinfo[BadFits[i][2]]), detector $(MetaData.dets_ged[BadFits[i][1]]): $(MetaData.dets_type[BadFits[i][1]])/$(BadFits[i][1])"
    µ_fit  = CalPars.µ_ADC[BadFits[i],:]
    pp_fit = CalPars.peaks_literature_keV[BadFits[i],:]
    result_calib, report_calib = fit_calibration(cal_pol_order, μ_fit, pp_fit; e_expression = "Energy_ADC")
                    
    residuals = mvalue.(pp_fit .- report_calib.f_fit.(µ_fit) .* unit.(pp_fit))
    residuals_norm = result_calib.gof.residuals_norm
    y_err_pred = residuals ./ residuals_norm

    ylres = round(Int,maximum(abs.(residuals_norm))+2)
    if ylres > 5
        ylt = 6
    else
        ylt = 3
    end
    plttitle = "$(MetaData.partinfo[BadFits[i][2]].period)-$(MetaData.partinfo[BadFits[i][2]].run), detector $(MetaData.dets_ged[BadFits[i][1]]) ($(MetaData.dets_type[BadFits[i][1]]), $(BadFits[i][1]))\n" * "max. |residual| =  $(round(typeof(residuals[1]),maximum(abs.(residuals)), digits = 2)) ($(round(maximum(abs.(residuals./pp_fit .* 100)),digits =2))%)"
    local plt = plot(report_calib, size = (600, 370), plot_title = plttitle, plot_titlefontsize = 7, thickness_scaling = 1.2, bottom_margin = 5mm, xtickfontsize = 6, ms = 3, ylims = (-ylres, ylres), yticks = -21:ylt:21)
    scatter!(mvalue.(µ_fit), ustrip.(pp_fit), xerr = 1e3 .* muncert.(µ_fit), ms = 1.5, color = :white, markerstrokecolor = :slategrey, label = false) 
    scatter!(mvalue.(µ_fit), ustrip.(pp_fit), yerr = 1e3 .* ustrip.(y_err_pred), ms = 1.5, color = :white, markerstrokecolor = :darkorange, label = "peak-fit uncertainty x 1e3") 
    scatter!(mvalue.(µ_fit), ustrip.(pp_fit), ms = 3, color = :black, markerstrokecolor = :black, label = false) 
    savefig(plt, path_plot .* "tmp.pdf")
    append_pdf!(fname_calfit, path_plot .* "tmp.pdf", cleanup = true)
end