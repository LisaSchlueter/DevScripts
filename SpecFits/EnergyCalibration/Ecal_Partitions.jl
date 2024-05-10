# 1. plot fwhm at Qbb for each run in each partition for each detector? 
# 2. plot position of FEP - literature value for each run in each partition for each detector?

# load per run, detector : fwhm, µ(FEP)
using DataFrames #Distributions, StatsBase,
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
partition = 1#[1,2]
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

# load fit parameter 
µ_FEP         = ones(length(dets_ged), nruns) .* NaN * u"keV"  .± NaN  * u"keV"
fwhm_FEP      = ones(length(dets_ged), nruns) .* NaN * u"keV"  .± NaN  * u"keV"
fwhm_Qbb  = ones(length(dets_ged), nruns) .* NaN * u"keV"  .± NaN  * u"keV"
pname =  :Tl208FEP
for i = 1:nruns
    for (d, det) in enumerate(dets_ged)  
        if haskey(pd_ecal_p1[i][det][e_type].fit,pname)
            µ_ADC  = pd_ecal_p1[i][det][e_type].fit[pname].µ ./ pd_ecal_p1[i][det][e_type].m_cal_simple
            fwhm_ADC = pd_ecal_p1[i][det][e_type].fit[pname].fwhm ./ pd_ecal_p1[i][det][e_type].m_cal_simple
            Tbl = Table(e_cusp = [µ_ADC, fwhm_ADC], qdrift = [0, 0])
            cal_func = ljl_propfunc(pd_ecal_p1[i][det][e_type].cal.func)
            µ_FEP[d,i] = collect(cal_func.(Tbl))[1] # calibrated peak position for period given run and detector 
            fwhm_FEP[d,i] = collect(cal_func.(Tbl))[2] .- pd_ecal_p1[i][det][e_type].cal.par[1] # calibrated (w/o offset)
        else
            µ_FEP[d,i] = NaN* u"keV"  .± NaN  * u"keV"
            fwhm_FEP[d,i] = NaN* u"keV"  .± NaN  * u"keV"
        end
        fwhm_Qbb_ADC = pd_ecal_p1[i][det][e_type].fwhm.qbb ./ pd_ecal_p1[i][det][e_type].m_cal_simple
        Tbl_fwhm = Table(e_cusp = [fwhm_Qbb_ADC, 1000], qdrift = [0, 0])
        cal_func = ljl_propfunc(pd_ecal_p1[i][det][e_type].cal.func)
        fwhm_Qbb[d,i] =  collect(cal_func.(Tbl_fwhm))[1] .- pd_ecal_p1[i][det][e_type].cal.par[1] # calibrated (w/o offset)
    end
end  
##########################################################################################################################################################################################################################################
# plot 1: how "good are partitions".
# calculate median fwhm_qbb for partition (for each det)
# look by how much values disagree

fwhm_median = mvalue.([mean(fwhm_Qbb[i,:]) for i = 1:length(dets_ged)]) # median over all runs 
# evaluate compatibility between median of fwhm at Qbb
chi2 = [sum(mvalue.(fwhm_Qbb[i,:] .- fwhm_median[i]).^2 ./ muncert.(fwhm_Qbb[i,:]).^2) for i = 1:length(dets_ged)]
pval = round.(ccdf.(Ref(Chisq(nruns -1)), chi2),digits = 2 )
stephist(pval, nbins = 100, xlabel = "p-value", fill = true, label = "partition $partition, all dets ", title = "FWHM at Qbb fluctuation around median fwhm")

# example plot
det_idx = 1
run_lbl = replace.(replace.(["$(entr.period) $(entr.run)" for entr in partinfo], " r00" => ",r"),"p0" => "p")

for i in eachindex(dets_ged)
    Plot_Evolution_fwhm(i)
end



function Plot_Evolution_fwhm(det_idx)
    title_str = "Calibration: partition $partition, $(dets_ged[det_idx]) ($(uppercase(String(dets_type[det_idx])))) "
    plt = hline([fwhm_median[det_idx]], linewidth = 2, color = :darkorange, linestyle = :dash, 
        label = "Median, p = $(pval[det_idx])")
    scatter!(1:nruns,fwhm_Qbb[det_idx,:], linecolor = :black,  markercolor = :white,
    ylabel = "FWHM of Qbb", xlabel  = "Partition $partition: period, run", xticks = (1:nruns,run_lbl),
    label = "Run-wise fit result (calibrated)" ,
    title = title_str,
    legend = :best,
    grid = false, framestyle = :box,
    xguidefontsize = 14, xtickfontsize = 10,
    yguidefontsize = 14, ytickfontsize = 10,
    legendfontsize = 12,
    dpi = 300, 
    )
    scatter!(1:nruns,mvalue.(fwhm_Qbb[det_idx,:]), linecolor = :transparent,  markercolor = :black, label = false)
    if pval[det_idx] < 0.05
        path_plot = "$(@__DIR__)/plots/p$partition/Ecal_FWHM/LowP/"
    else
        path_plot = "$(@__DIR__)/plots/p$partition/Ecal_FWHM/OtherP/"
    end
    if !ispath(path_plot)
        mkdir("$path_plot")
    end 
    fname = path_plot * "Ecal_partition$(partition)_FWHMrunwise_$(dets_type[det_idx])_$(dets_ged[det_idx]).png"
    savefig(plt, fname)
    @info "save plot to $fname "
end