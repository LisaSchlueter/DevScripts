#= p-value overview 
plot: p-values of gamma-peak fits 
- 1 panel: overview all peaks 
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

#settings: select data and dsp (energy) output 
partition = 2
e_type = :e_cusp_ctc

path_plot1 = "$(@__DIR__)/plots/p$partition/FitPar/"
if !ispath(path_plot1)
    mkdir("$path_plot1")
end 

FitPars, MetaData = get_peakfit_rpars_partition(DataPartition(partition); reload = false, e_type = e_type);
keys(FitPars)
pvalues_peakfit = FitPars.pvalue

# prepare data for plot 
Tbl_lbl = ["E = $(round(typeof(1u"keV"), MetaData.th228_literature[i], digits = 0))" for i=1:MetaData.npeaks]

# common plot arguments 
fs = 14
colors = [:dodgerblue, :darkorange, :lightseagreen, :midnightblue, :crimson , :olive ,:violet]
HistArg = Dict(:normalize => :probability,
               :fill => :true,
               :foreground_color_legend => :silver,
               :background_color_legend => :white,
               :grid => :off,
               :xlims => (0,1),
               :xguidefontsize => fs+2, :xtickfontsize => fs-2,
               :yguidefontsize => fs+2, :ytickfontsize => fs-2,
               :legendfontsize => fs-2)

xlbl = L"\textit{E}_\textrm{cal.} - \textit{E}_\textrm{lit.}" * "(keV)"
nbins = 400
data = pvalues_peakfit

################################################################
# plot pvalues of peak fits. 1 per detector per run per peak 
# plot overview for all detectors : histograms below each other 
median_tot = median(filter(isfinite,reshape(data,:)))
pall = stephist(reshape(data,:), bins = 100, 
                color = :darkgrey, fillalpha = 0.5, 
                framestyle = :box, 
                legend = :topleft,
                ylims = (0,:auto),
                label = "All peaks",
                ylabel = "Probability",
                xlabel  = "p-value",
                ; HistArg...)
stat_str = "median = " * @sprintf("%.2f",median_tot)
annotate!(0.95, 0.095, text(stat_str, :right, fs-2, :grey))
hline!([0.01], color = :black, label = "Uniform (expection)", linestyle = :dash, linewidth = 1.5)

medians = [median(filter(isfinite, reshape(data[:,:,i],:))) for i = 1:MetaData.npeaks]
p = Vector(undef,MetaData.npeaks)
for i = 1:MetaData.npeaks
    stat_str = "median = " * @sprintf("%.2f",medians[i])
    lbl = "E = $(round(typeof(1u"keV"), MetaData.th228_literature[i], digits = 0))" 
    if MetaData.th228_literature[i] == 1592.0u"keV"
        f_alpha = 0.3
        lbl = lbl * "\n(DEP, excluded)"
    elseif MetaData.th228_literature[i] == 2103.0u"keV"
       f_alpha = 0.3
       lbl = lbl * "\n(SEP, excluded)"
    else
        f_alpha = 0.8
    end

    p[i] = stephist(reshape(data[:,:,i],:), bins = 100, 
            color = colors[i], 
            fillalpha = f_alpha,  
            legend = :topleft,
            yaxis=false, yticks = [0,0.05],
            ylims = (0,:auto),
            label=lbl;
            HistArg...)  
    annotate!(0.95, 0.085, text(stat_str, :right, fs-2, colors[i]))
   # hline!([0.01], color = :black, label = false, linestyle = :dash, linewidth = 1.5)
    if i==MetaData.npeaks
        plot!(p[i], xlabel  =  "p-value" )
    else
        # myxticks = xticks(p[i])
        # Set the labels to empty string
        plot!(p[i],bottom_margin = -7mm, xformatter=_->"",ylabel = " ")#, xticks = (myxticks, fill(" ", length(myxticks))))
    end
end

l = @layout([a{0.2h}; grid(7, 1)])
ptot = plot(pall,p...,layout =l, size = (600,1200), xlims = (0,1), ylims = (0,0.11), left_margin = 10mm, right_margin = 5mm)

fname = path_plot1 * "Ecal_PeakFit_pvalueDist_part$(partition)_$(e_type).png"
savefig(ptot,fname)
@info "save plot to $fname"
################################################################
# plot pvalues of calibration curve fit 
# plot overview for all detectors : histograms below each other 
# nruns = MetaData.nruns
# pvalues_calfit_50 = [median(filter(isfinite,pvalues_calfit[:,i]))  for i = 1:nruns]  
# pvalues_calfit_95 = [quantile(filter(isfinite,pvalues_calfit[:,i]),0.95)  for i = 1:nruns] 
# pvalues_calfit_5 = [quantile(filter(isfinite,pvalues_calfit[:,i]),0.05)  for i = 1:nruns] 
# pvalues_calfit_sigmah = [quantile(filter(isfinite,pvalues_calfit[:,i]),1-(1-0.6827)/2)  for i = 1:nruns] 
# pvalues_calfit_sigmal = [quantile(filter(isfinite,pvalues_calfit[:,i]),(1-0.6827)/2)  for i = 1:nruns] 

# x = 1:11
# default(framestyle = :box, size = (650,400), xguidefontsize = 14, yguidefontsize = 14, xtickfontsize = 12, ytickfontsize = 12, legendfontsize = 12, grid = :off, dpi = 300)
# plt = plot(x, pvalues_calfit_5, fillrange = pvalues_calfit_95, fillalpha = 1, fillcolor = :lightgray, c = :transparent, label = "90% band", linewidth = 0)
# plot!(x, pvalues_calfit_sigmal, fillrange = pvalues_calfit_sigmah, fillalpha = 1, fillcolor = :deepskyblue, c = :transparent, label = "68.3% band", linewidth = 0)
# plot!(x, pvalues_calfit_50, c = :darkorange, label = "median ", linewidth = 2)
# plot!(plt,ylims = (0,1), 
#     xlims = (1,11), 
#     xlabel = "Calibration runs", xticks = 1:1:11,
#     ylabel = "p-value", 
#     left_margin = 10mm, bottom_margin = 10mm,
#     title = "Calibration curve fit, all dets, partition $(partition), $e_type", titlefontsize = 12)

# fname_calfit = path_plot2 * "Ecal_CalFit_pvalueDist_part$(partition)_$(e_type).png"
# savefig(plt,fname_calfit)