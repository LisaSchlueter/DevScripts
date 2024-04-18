#= comparison µ peak fits vs literature values of calibration peaks 
plot :
- residual histogram: µ - literature . for all peaks and breakdown of single peaks (used for claibration)
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
using StatsPlots
include("../SanityPlots/utils.jl")

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
    
# prepare data for plot 
residual = permutedims(permutedims(µ_all,(3,1,2)) .- th228_literature,(2,3,1))
residual_norm = mvalue.(residual) ./ muncert.(residual)
Tbl_lbl = ["E = $(round(typeof(1u"keV"), th228_literature[i], digits = 0))" for i=1:npeaks]
colors = [:dodgerblue, :darkorange, :lightseagreen, :crimson, :violet]
# Tbl = Table(residual_norm = reshape(permutedims(residual_norm,(3,1,2)),:), peak = reshape(repeat(round.(typeof(1u"keV"), th228_literature, digits=0), outer=(1,nruns)),:))
# xmax1 = ceil(Int,maximum(abs.(filter(isfinite,Tbl.residual_norm))))

# common plot arguments 
fs = 14
colors = [:dodgerblue, :darkorange, :lightseagreen, :crimson, :violet]
HistArg = Dict(:normalize => :none,
               :fill => :true,
               :foreground_color_legend => :silver,
               :background_color_legend => :white,
               :grid => :off,
               :xlims => (-10,10),
               :xguidefontsize => fs, :xtickfontsize => fs-2,
               :yguidefontsize => fs, :ytickfontsize => fs-2,
               :legendfontsize => fs-2)

Mode = :keV#:norm#
if Mode == :norm
    data = residual_norm
    stat_x = 14
    xl = (-15,15)
    xlbl = L"\textit{E}_\textrm{cal.} - \textit{E}_\textrm{lit.} (σ)"
    xunit = "σ"
    nbins = 300
elseif Mode == :keV
    data = ustrip.(mvalue.(residual))
    stat_x = 0.9
    xl = (-1,1)
    xlbl = L"\textit{E}_\textrm{cal.} - \textit{E}_\textrm{lit.}" * "(keV)"
    xunit = "keV"
    nbins = 400
end
# plot overview for all detectors : histograms below each other 
bin_center, bin_edges, _, counts = get_histstats(ustrip.(mvalue.(data)); nbins = nbins)
median_tot = median(filter(isfinite,reshape(data,:)))
fwhm_tot = get_fwhm(bin_center,counts)
pall = stephist(reshape(data,:), bins = bin_edges, 
                color = :darkgrey, fillalpha = 0.5, 
                framestyle = :box, 
                legend = :topleft,
                ylims = (0,:auto),
                label = "All peaks",
                ylabel = "Occurrence",
                xlabel  = xlbl 
                ; HistArg...)
stat_str = "median = " * @sprintf("%.2f%s\n",ustrip(mvalue(median_tot)),xunit) *  @sprintf("fwhm = %.2f%s",fwhm_tot,xunit)
annotate!(stat_x, ylims()[2]*0.9, text(stat_str, :right, fs-2, :grey))
vline!([0],label = :none, color = :black, linestyle = :dash, linewidth = 2)
medians = fill(0.0,npeaks)
fwhm = fill(0.0,npeaks)
for i = 1:npeaks
    bin_center, counts = nothing, nothing
    bin_center, _ , _, counts = get_histstats(mvalue.(ustrip.(data[:,:,i])); nbins = 300)
    medians[i] = median(filter(isfinite, mvalue.(ustrip.(reshape(data[:,:,i],:)))))
    fwhm[i] = get_fwhm(bin_center,counts)
    @info i
end
p = Vector(undef,npeaks)
for i = 1:npeaks
    stat_str = "median = " * @sprintf("%.2f%s\n",medians[i],xunit) *  @sprintf("fwhm = %.2f%s",fwhm[i],xunit)
    lbl = "E = $(round(typeof(1u"keV"), th228_literature[i], digits = 0))" 
    p[i] = stephist(reshape(data[:,:,i],:), bins = bin_edges, 
            color = colors[i], 
            fillalpha = 0.8,  
            legend = :left,
            yaxis=false, 
            ylims = (0,:auto),
            label=lbl;
            HistArg...)  
    annotate!(stat_x, ylims()[2]*0.5, text(stat_str, :right, fs-2, colors[i]))
    vline!([0],label = :none, color = :black, linestyle = :dash, linewidth = 2)
    if i==npeaks
        plot!(p[i], xlabel  = xlbl )
    else
        # myxticks = xticks(p[i])
        # Set the labels to empty strings
        plot!(p[i],bottom_margin = -7mm, xformatter=_->"")#, xticks = (myxticks, fill(" ", length(myxticks))))
    end
end
l = @layout([a{0.35h}; grid(5, 1)])
ptot = plot(pall,p...,layout =l, size = (650,1000), xlims = xl, left_margin = 10mm, right_margin = 5mm)

path_plot = "$(@__DIR__)/plots/Ecal_FitPar_Overview/"
fname = path_plot * "Ecal_Residuals_$(Mode)_part$(partition)_$(e_type).png"
savefig(ptot,fname)

function get_histstats(data; nbins::Int = 300, bin_edges::Vector= [])
    if length(size(data)) > 1
        data = reshape(data,:)
    end
    if any(isnan.(data))
        data = filter(isfinite,data)
    end
    if isempty(bin_edges)
        h = fit(Histogram, data, nbins = nbins) 
    else
        h = fit(Histogram, data, edges)
    end 
    counts = h.weights
    edges = h.edges[1]
    width = diff(edges)[1]
    center = edges .+ width/2
    return center, edges, width, counts
end


function get_fwhm(centers::StepRangeLen, counts::Vector)
    max_idx = argmax(counts)
    max_val = counts[max_idx]
    half_max = max_val / 2
    left_idx = findfirst(counts[1:max_idx] .>= half_max)
    right_idx = findfirst(counts[max_idx:end] .<= half_max) + (max_idx-1)
    try
        interp_left = LinearInterpolation(counts[left_idx-1:left_idx+1],centers[left_idx-1:left_idx+1])
        interp_right = LinearInterpolation(counts[right_idx+1:-1:right_idx-1],centers[right_idx+1:-1:right_idx-1])
        fwhm = interp_right(half_max) - interp_left(half_max)
        return fwhm 
    catch
        fwhm = centers[right_idx] - centers[left_idx]
        return fwhm
    end
end