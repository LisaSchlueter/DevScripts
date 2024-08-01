#= comparison gamma peak fits with and without linear background slope. 
--> for all dets or only ICPC, Coax, Bege,....
--> peakwise or all peaks 
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
include("$(@__DIR__)/utils_plot.jl")


fontsize = 16 

# data selection 
e_type = :e_cusp_ctc
periods = [3, 4, 6, 7, 8, 9]
exclPeaks = []# [:Tl208DEP, :Tl208SEP] # exclude peaks from plot

# load fit results 
l200 = LegendData(:l200)
FitPars, MetaData = get_rpars_peakfit(l200, periods; reload = false, e_type = e_type, cal_type = :ecal);
FitPars_tails, MetaData_tails = get_rpars_peakfit(l200, periods; reload = false, e_type = e_type, cal_type = :ecal_tails);
FitPars_bslope, MetaData_bslope = get_rpars_peakfit(l200, periods; reload = false, e_type = e_type, cal_type = :ecal_bslope);
MetaData.ecal_config.th228_names[4]
FitPars1 = FitPars_tails#bslope
FitPars2 = FitPars
MetaData1 = MetaData_tails
MetaData2 = MetaData
skew_frac_diff = ustrip.(mvalue.(FitPars1.skew_frac)) .- ustrip.(mvalue.(FitPars2.skew_frac))
chi2_diff = ustrip.(mvalue.(FitPars1.chi2)) .- ustrip.(mvalue.(FitPars2.chi2))

_def_plot(; fs = 14)
colors = _get_def_peakcolors()

# plot low energy tail fraction vs chi2 difference. 
pall = scatter(reshape(skew_frac_diff,:), reshape(chi2_diff,:), 
        color = :darkgrey, alpha = 0.1, markerstrokewidth=0, 
        margins = 3mm,
        label = "all peaks, model: slope - default", 
        ylabel = "Δ χ²",
        xlabel = "Δ Low-energy tail fraction")
ylims!(-20, 20)
xlims!(-0.105,0.105)
peakplt =  p = Vector(undef, length(MetaData.ecal_config.th228_names))
for (peak, peakname) in enumerate(MetaData.ecal_config.th228_names)
    if peakname in exclPeaks
        continue
    end
    peakplt[peak] = scatter(reshape(skew_frac_diff[:, :, peak],:), reshape(chi2_diff[:,:, peak],:), 
        color = colors[Symbol(peakname)], alpha = 0.1, markerstrokewidth=0, 
        label = peakname, 
        ylabel = "Δ χ²",
        framestyle = :off)
        ylims!(-20, 20)
        xlims!(-0.105,0.105)
        if peak == 7
            xlabel!("Δ Low-energy tail fraction")
        end
end
ptot = plot(pall, peakplt..., layout = (8, 1), size = (600, 300*8), legend = :topright, leftmargin = 20mm)


path_plot = "$(@__DIR__)/plots/period$(join(string.(periods)))/FitParDiff/"
fname = path_plot * "Ecal_PeakFitDiff_all_chi2skewfrac_period$(join(string.(periods)))_$(e_type)_$(MetaData1.cal_type)minus$(MetaData2.cal_type).png"
savefig(ptot, fname)
@info "save plot to $fname"


# plot only DEP
stephist()