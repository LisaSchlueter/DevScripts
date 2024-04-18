#= comparison µ peak fits vs literature values of calibration peaks 
plot:
- µ with uncertainty and literature as a function of run. at the right hand side a histogram. for each peak. 
- measure of goodness: pvalue
- mean and error of the mean
- open question: for all detectors together AND for each detector type separately. 
=#
using Distributions, StatsBase, DataFrames
using LegendDataManagement
using LegendDataManagement.LDMUtils
using LegendSpecFits
using Measurements
using Measurements: value as mvalue, uncertainty as muncert
using Measures
using Plots, Printf, LaTeXStrings
using PropDicts
using TypedTables
using Unitful
using DataFrames
include("../SanityPlots/utils.jl")

path_plot = "$(@__DIR__)/plots/"

l200 = LegendData(:l200)
#select data and dsp output 
partition = 1
e_type = :e_cusp_ctc

# open data
partinfo = partitioninfo(l200)[DataPartition(partition)]
filekey = start_filekey(l200, (partinfo[1].period, partinfo[1].run, :cal)) 
chinfo = Table(channelinfo(l200, filekey; system=:geds, only_processable=true))
dets_ged = chinfo.detector
# dets_icpc = dets_ged[chinfo.det_type .== :icpc]

# load all ProbDicts (for all period-run combination in selected partition). speed up load probdict for partition 
pd_ecal_p1 = [l200.par.rpars.ecal[entr.period, entr.run] for entr in partinfo] # takes a while 
nruns = length(pd_ecal_p1)
th228_literature =sort([pd_ecal_p1[1][Symbol(dets_ged[1])][e_type].cal.peaks..., 1592u"keV", 2103u"keV"]) # get literature values mit denen gefittet wurde. interpolation st Qbb
Qbb = 2039.0u"keV"
npeaks = length(th228_literature)
peak_keys = Symbol.(collect(keys(pd_ecal_p1[1][dets_ged[1]][e_type].fit)))

# all calibration fit results for a specific detector
pvalues_calfit = NaN .* zeros((length(dets_ged), nruns))
pvalues_peakfit = ones(length(dets_ged), nruns, npeaks) .* NaN

for i = 1:nruns
    for (d, det) in enumerate(dets_ged)
        pvalues_calfit[d,i] = pd_ecal_p1[i][det][e_type].cal.gof.pvalue # no this is the pvalue of the calibration curve fit. 
        pvalues = [pd_ecal_p1[i][det][e_type].fit[peak].gof.pvalue for peak in peak_keys if haskey(pd_ecal_p1[i][det][e_type].fit,peak)]
       
        if length(pvalues) !== npeaks
            µ_tmp = [pd_ecal_p1[i][det][e_type].fit[peak].µ for peak in peak_keys if haskey(pd_ecal_p1[i][det][e_type].fit,peak)]
            # find missing peak and put NaN
            pvalues_tmp = Vector{Float64}(undef,npeaks)
            for p = 1:npeaks
                idx = findfirst(abs.(µ_tmp .- th228_literature[p]) .<= 10u"keV")
                if !isnothing(idx) 
                    pvalues_tmp[p] = pvalues[idx]
                else
                    @info "peak missing: run $i, det $det ($d), peak $(th228_literature[p]), replace with NaN"
                    pvalues_tmp[p] =  NaN
                end
            end
            pvalues = pvalues_tmp
        end

        pvalues_peakfit[d,i,:] =  pvalues
    end
end
   
# prepare data for plot 
Tbl_lbl = ["E = $(round(typeof(1u"keV"), th228_literature[i], digits = 0))" for i=1:npeaks]

# common plot arguments 
fs = 14
colors = [:dodgerblue, :darkorange, :lightseagreen, :midnightblue, :crimson , :olive ,:violet]
HistArg = Dict(:normalize => :probability,
               :fill => :true,
               :foreground_color_legend => :silver,
               :background_color_legend => :white,
               :grid => :off,
               :xlims => (0,1),
               :xguidefontsize => fs, :xtickfontsize => fs-2,
               :yguidefontsize => fs, :ytickfontsize => fs-2,
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

medians = [median(filter(isfinite, reshape(data[:,:,i],:))) for i = 1:npeaks]
p = Vector(undef,npeaks)
for i = 1:npeaks
    stat_str = "median = " * @sprintf("%.2f",medians[i])
    lbl = "E = $(round(typeof(1u"keV"), th228_literature[i], digits = 0))" 
    if th228_literature[i] == 1592.0u"keV"
        f_alpha = 0.3
        lbl = lbl * "\n(DEP, excluded)"
    elseif th228_literature[i] == 2103.0u"keV"
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
    if i==npeaks
        plot!(p[i], xlabel  =  "p-value" )
    else
        # myxticks = xticks(p[i])
        # Set the labels to empty strings
        plot!(p[i],bottom_margin = -7mm, xformatter=_->"",ylabel = " ")#, xticks = (myxticks, fill(" ", length(myxticks))))
    end
end
l = @layout([a{0.2h}; grid(7, 1)])
ptot = plot(pall,p...,layout =l, size = (600,1200), xlims = (0,1), ylims = (0,0.11), left_margin = 10mm, right_margin = 5mm)

fname = path_plot * "Ecal_PeakFit_pvalueDist_part$(partition)_$(e_type).png"
savefig(ptot,fname)

################################################################
# plot pvalues of calibration curve fit 
# plot overview for all detectors : histograms below each other 
pvalues_calfit_50 = [median(filter(isfinite,pvalues_calfit[:,i]))  for i = 1:nruns]  
pvalues_calfit_95 = [quantile(filter(isfinite,pvalues_calfit[:,i]),0.95)  for i = 1:nruns] 
pvalues_calfit_5 = [quantile(filter(isfinite,pvalues_calfit[:,i]),0.05)  for i = 1:nruns] 
pvalues_calfit_sigmah = [quantile(filter(isfinite,pvalues_calfit[:,i]),1-(1-0.6827)/2)  for i = 1:nruns] 
pvalues_calfit_sigmal = [quantile(filter(isfinite,pvalues_calfit[:,i]),(1-0.6827)/2)  for i = 1:nruns] 

x = 1:11
default(framestyle = :box, size = (650,400), xguidefontsize = 14, yguidefontsize = 14, xtickfontsize = 12, ytickfontsize = 12, legendfontsize = 12, grid = :off, dpi = 300)
plt = plot(x, pvalues_calfit_5, fillrange = pvalues_calfit_95, fillalpha = 1, fillcolor = :lightgray, c = :transparent, label = "90% band", linewidth = 0)
plot!(x, pvalues_calfit_sigmal, fillrange = pvalues_calfit_sigmah, fillalpha = 1, fillcolor = :deepskyblue, c = :transparent, label = "68.3% band", linewidth = 0)
plot!(x, pvalues_calfit_50, c = :darkorange, label = "median ", linewidth = 2)
plot!(plt,ylims = (0,1), 
    xlims = (1,11), 
    xlabel = "Calibration runs", xticks = 1:1:11,
    ylabel = "p-value", 
    left_margin = 10mm, bottom_margin = 10mm,
    title = "Calibration curve fit, all dets, partition $(partition), $e_type", titlefontsize = 12)

fname_calfit = path_plot * "Ecal_CalFit_pvalueDist_part$(partition)_$(e_type).png"
savefig(plt,fname_calfit)