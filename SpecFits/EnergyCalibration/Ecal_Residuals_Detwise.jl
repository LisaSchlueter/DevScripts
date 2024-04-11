#= comparison µ peak fits vs literature values of calibration peaks =#
using Distributions, StatsBase, DataFrames
using LegendDataManagement
using LegendDataManagement.LDMUtils
using LegendSpecFits
using Measurements
using Measurements: value as mvalue, uncertainty as muncert
using Measures
using Plots
using Printf, LaTeXStrings
using PropDicts
using TypedTables
using Unitful
using DataFrames
using StatsPlots
include("../SanityPlots/utils.jl")

path_plot = "PackageDevScripts/DevelopSpecFits/EnergyCalibration/plots/Ecal_Residuals_dets"
path_rel_plot = relpath(path_plot,pwd())
path_abs_plot = pwd() * "/" * path_rel_plot * "/"

l200 = LegendData(:l200)

#select data and dsp output 
partition = 1
partinfo = partitioninfo(l200)[DataPartition(partition)]
e_type = :e_cusp_ctc

# open data
filekey = start_filekey(l200, (partinfo[1].period, partinfo[1].run, :cal)) 
chinfo = Table(channelinfo(l200, filekey; system=:geds, only_processable=true))
dets_ged = chinfo.detector
dets_type = chinfo.det_type
# dets_icpc = dets_ged[chinfo.det_type .== :icpc]

# load all ProbDicts (for all period-run combination in selected partition). speed up load probdict for partition 
pd_ecal_p1 = [l200.par.rpars.ecal[entr.period, entr.run] for entr in partinfo] # takes a while 
nruns = length(pd_ecal_p1)
#literature values 
th228_literature = sort(pd_ecal_p1[1][Symbol(dets_ged[1])][e_type].cal.peaks) # get literature values mit denen gefittet wurde. interpolation st Qbb
Qbb = 2039.0u"keV"
npeaks = length(th228_literature)

# all calibration fit results for a specific detector
µ_all = ones(length(dets_ged), nruns, npeaks) .* NaN * 1u"keV" .± NaN * 1u"keV"
pvalues = zeros((length(dets_ged), nruns))
for i = 1:nruns
    for (d, det) in enumerate(dets_ged)
        pvalues[d,i] = pd_ecal_p1[i][det][e_type].cal.gof.pvalue 
        µ_ADC = pd_ecal_p1[i][det][e_type].cal.µ
        Tbl = Table(e_cusp = µ_ADC, qdrift = fill(0.0,length(µ_ADC)))
        cal_func = ljl_propfunc(pd_ecal_p1[i][det][e_type].cal.func)
        peakpos_cal = sort(collect(cal_func.(Tbl))) # calibrated peak position for period given run and detector 
        if length(peakpos_cal) !== npeaks
            # find missing peak and put NaN
            peakpos_cal_tmp = Vector{Quantity{Measurement{Float64}}}(undef,npeaks)
            for p = 1:npeaks
                idx = findfirst(abs.(peakpos_cal .- th228_literature[p]) .<= 10u"keV")
                if !isnothing(idx)
                    @info "peak missing: run $i, det $det ($d), peak $(th228_literature[p]), replace with NaN"
                    peakpos_cal_tmp[p] = peakpos_cal[idx]
                else
                    peakpos_cal_tmp[p] =  NaN * 1u"keV" ±  NaN * 1u"keV"
                end
            end
            peakpos_cal = peakpos_cal_tmp
        end
        µ_all[d,i,:] = peakpos_cal
    end
end
    
# plot overview for a single detector
# det_idx = 1
for det_idx in eachindex(dets_ged)
    det = dets_ged[det_idx]
    pval_det = pvalues[det_idx,:]
    residual = µ_all[det_idx,:,:]' .- th228_literature
    residual_norm = mvalue.(residual) ./ muncert.(residual)

    Tbl_lbl = ["E = $(round(typeof(1u"keV"), th228_literature[i], digits = 0))" for i=1:npeaks]
    colors = [:dodgerblue, :darkorange, :lightseagreen, :crimson, :violet]
    Tbl = Table(residual_norm = reshape(residual_norm,:), peak = reshape(repeat(round.(typeof(1u"keV"), th228_literature, digits=0), outer=(1,nruns)),:))
    xmax1 = ceil(Int,maximum(abs.(filter(isfinite,Tbl.residual_norm))))
    hist = @df Tbl groupedhist(:residual_norm, group = :peak, bar_position = :stack,
                        nbins = 50, 
                        xlims = (-xmax1, xmax1),
                        color = [:dodgerblue :darkorange :lightseagreen :crimson :violet],
                        xlabel = L"\textit{E}_\textrm{cal.} - \textit{E}_\textrm{lit.} (σ)", ylabel = "Occurence",
                        legendtitle = "228Th fits", legend = :topleft,  foreground_color_legend = :silver, background_color_legend = :white, legendtitlefontsize = fs-2,
                        size = (600,400), grid = :off)
    median(Tbl.residual_norm)          
    annotate!(xmax1-0.1, ylims()[2]-0.5, text(uppercase("$det ($(dets_type[det_idx]))"), :right, fs-3, :grey))
    xmin = minimum(mvalue.(reshape(residual,:)) .- muncert.(reshape(residual,:))) - 0.15u"keV"
    xmax = maximum(mvalue.(reshape(residual,:)) .+ muncert.(reshape(residual,:)))
    plt = vline([0u"keV"],label = :none, color = :black, linestyle = :dot, linewidth = 2)
    for i = 1:npeaks 
        scatter!(residual[i,:], range(1,11) .+ (-0.125 + (i-1) *0.05),
                ms = 0, linewidth = 1.5, linecolor = colors[i], markercolor = :white, alpha = 1,
                label = :none, legend = :none,
                yticks = (1:11, map(x->"$(x.period), $(x.run)",partinfo)))
        scatter!(mvalue(residual[i,:]), range(1,11) .+ (-0.125 + (i-1) *0.05),
                ms = 4, markerstrokecolor = colors[i], markercolor = colors[i],
                label = "E = $(round(typeof(1u"keV"), th228_literature[i], digits = 0))",
                foreground_color_legend = :silver,
                background_color_legend = :white,
                left_margin = 10mm, xlabel = L"\textit{E}_\textrm{cal.} - \textit{E}_\textrm{lit.}",
                xlims = (xmin,xmax), grid = :on)     
    end
    fs = 12
    for r = 1:nruns
        annotate!(ustrip(xmin)+0.05, r, text("p = "  * @sprintf("%.2g",pval_det[r]), :left, fs-3, :grey))
    end

    plt_all = plot(hist,plt,layout =@layout([a{0.35h}; b{0.65h}]),size = (500,700),
        framestyle = :box,
        xlabelfontsize = fs, ylabelfontsize = fs, xtickfontsize = fs-3, ytickfontsize = fs-3, 
        legendfontsize = fs-3, legendtitlefontsize = fs-3,
        dpi = 300)

    fname = path_abs_plot * "Ecal_Residuals_part$(partition)_$(e_type)_det$(det).png"
    savefig(plt_all,fname)
 end