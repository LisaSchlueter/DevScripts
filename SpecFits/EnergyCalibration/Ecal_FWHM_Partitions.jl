# 1. plot fwhm at Qbb for each run in each partition for each detector? 
# 2. plot position of FEP - literature value for each run in each partition for each detector?

# load per run, detector : fwhm, µ(FEP)
using DataFrames
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
using Statistics, Distributions
include("../SanityPlots/utils.jl")

#select data and dsp output 
partition = 1
e_type = :e_cusp_ctc

l200 = LegendData(:l200)

# open data
partinfo = partitioninfo(l200)[DataPartition(partition)]
filekey = start_filekey(l200, (partinfo[1].period, partinfo[1].run, :cal)) 
chinfo = Table(channelinfo(l200, filekey; system=:geds, only_processable=true))
dets_ged = chinfo.detector
dets_type = chinfo.det_type

# load all ProbDicts (for all period-run combination in selected partition). speed up load probdict for partition 
pd_ecal_p1 = [l200.par.rpars.ecal[entr.period, entr.run] for entr in partinfo] # takes a while
nruns = length(pd_ecal_p1)

# plot:  are FWHM at Qbb stable within a partition (for a given detector)?
function Plot_Evolution_fwhm(det_idx; saveplot = false)
    det_lbl = "$(dets_ged[det_idx]) ($(uppercase(String(dets_type[det_idx])))) "
    run_lbl = replace.(replace.(["$(entr.period) $(entr.run)" for entr in partinfo], " r00" => ",r"),"p0" => "p")

    good_runs = [(haskey(pd_ecal_p1[i], dets_ged[det_idx]) ? i : NaN ) for i in 1:nruns]
    good_runs = filter(x -> !isnan(x), good_runs)
    fwhm_runs =  [pd_ecal_p1[i][dets_ged[det_idx]][e_type].fwhm.qbb for i in good_runs]
    fwhm_wmean = weightedmean(ustrip.(fwhm_runs))u"keV"

    # evaluate compatibility between median of fwhm at Qbb
    chi2 = sum(mvalue.(fwhm_runs .- fwhm_wmean).^2 ./ muncert.(fwhm_runs).^2)
    pval = round.(ccdf.(Ref(Chisq(nruns -1)), chi2), digits = 2)
    
    plt = scatter(1:nruns, fwhm_runs, linecolor = :black,  markercolor = :white,
        ylabel = "FWHM of Qbb", xlabel  = "Period, run",
        xticks = (1:nruns,run_lbl), xrotation = 0 ,
        label = false ,
        legend = :best,
        grid = false, framestyle = :box,
        xguidefontsize = 14, xtickfontsize = 10,
        yguidefontsize = 14, ytickfontsize = 10,
        legendfontsize = 12,
        dpi = 300, 
        foreground_color_legend = :silver,
        background_color_legend = :white
        )
    scatter!(1:nruns, mvalue.(fwhm_runs), linecolor = :transparent,  markercolor = :black, label = "Run-wise fit results " * det_lbl)
    hline!([mvalue(fwhm_wmean)], linewidth = 2, color = :dodgerblue, linestyle = :dash, label = "Weighted mean; p = $pval")

    if saveplot
        if pval < 0.05
            path_plot = "$(@__DIR__)/plots/p$partition/Ecal_FWHM/LowP/"
        else
            path_plot = "$(@__DIR__)/plots/p$partition/Ecal_FWHM/OtherP/"
        end
        if !ispath(path_plot)
            mkpath("$path_plot")
        end 
        fname = path_plot * "Ecal_partition$(partition)_FWHMrunwise_$(dets_type[det_idx])_$(dets_ged[det_idx]).png"
        savefig(plt, fname)
        @info "save plot to $fname "
    end 
    title!("Partition $partition", fontsize = 12)
    return plt
end

# plot all detectors 
for (i, det) in enumerate(dets_ged)
    @info "det = $det ($i)"
    if all([~haskey(pd_ecal_p1[i], det) for i in 1:nruns])
        continue
    end
    Plot_Evolution_fwhm(i, saveplot = true)
end


# overall stability: pvalue distribution over all detectors
pval = []
for (det_idx, det) in enumerate(dets_ged)
    @info "det = $det ($det_idx)"
    if all([~haskey(pd_ecal_p1[i], det) for i in 1:nruns])
        continue
    end
    good_runs = [(haskey(pd_ecal_p1[i], dets_ged[det_idx]) ? i : NaN ) for i in 1:nruns]
    good_runs = filter(x -> !isnan(x), good_runs)
    fwhm_runs =  [pd_ecal_p1[i][dets_ged[det_idx]][e_type].fwhm.qbb for i in good_runs]
    fwhm_wmean = weightedmean(ustrip.(fwhm_runs))u"keV"
    chi2 = sum(mvalue.(fwhm_runs .- fwhm_wmean).^2 ./ muncert.(fwhm_runs).^2)
    local pval_det = round.(ccdf.(Ref(Chisq(nruns -1)), chi2), digits = 2)
    push!(pval, pval_det)
end
default(guidefontsize=16, legendfontsize = 12, xtickfontsize = 12, ytickfontsize = 12, titlefontsize = 12)
stephist(pval, nbins = 100, xlabel = "p-value", fill = true, label = "FWHM stability for each detector", title = "FWHM at Qbb fluctuation around median fwhm", framestyle = :box, grid = false)
title!("Partition $partition")

