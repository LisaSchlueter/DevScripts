#= 
fit paramater of peaks fit (calibration data)
look peak fits with bad pvalue. look at fit parameter , especially µ and \sigma 
different "modes" available: modes = [:norm, :keV, :sigma, :fwhm, :skew_frac]
=#
using Distributions, StatsBase, DataFrames, Statistics
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
# using StatsPlots
using Colors, ColorSchemes
include("../SanityPlots/utils.jl")

#select data and dsp output
l200 = LegendData(:l200) 
partition = 1
e_type = :e_cusp_ctc

# plotting path 
path_plot = "$(@__DIR__)/plots/p$partition/FitPar/"
if !ispath(path_plot)
    mkdir("$path_plot")
end 

# open data
partinfo = partitioninfo(l200)[DataPartition(partition)]
filekey = start_filekey(l200, (partinfo[1].period, partinfo[1].run, :cal)) 
chinfo = Table(channelinfo(l200, filekey; system=:geds, only_processable=true))
dets_ged = chinfo.detector
energy_config = dataprod_config(l200).energy(filekey)

# load all ProbDicts (for all period-run combination in selected partition). speed up load probdict for partition 
pd_ecal_p1 = [l200.par.rpars.ecal[entr.period, entr.run] for entr in partinfo] # takes a while 
nruns = length(pd_ecal_p1)

# load ALL peaks 
th228_names = Symbol.(keys(pd_ecal_p1[1][Symbol(dets_ged[1])][e_type].fit))
npeaks = length(th228_names)

# load fit parameter 
skew_frac = ones(length(dets_ged), nruns, npeaks) .* NaN  .± NaN 
skew_width = ones(length(dets_ged), nruns, npeaks) .* NaN  .± NaN   
µ  = ones(length(dets_ged), nruns, npeaks) .* NaN * u"keV"  .± NaN  * u"keV"
σ  = ones(length(dets_ged), nruns, npeaks) .* NaN * u"keV"  .± NaN  * u"keV"
fwhm  = ones(length(dets_ged), nruns, npeaks) .* NaN * u"keV"  .± NaN  * u"keV"
background = ones(length(dets_ged), nruns, npeaks) .* NaN  .± NaN 
n = ones(length(dets_ged), nruns, npeaks) .* NaN  .± NaN 
step_amplitude = ones(length(dets_ged), nruns, npeaks) .* NaN  .± NaN 
pvalues_peakfit = ones(length(dets_ged), nruns, npeaks) .* NaN

for i = 1:nruns
    for (d, det) in enumerate(dets_ged)  
        for (p, pname) in enumerate(th228_names)
            if haskey(pd_ecal_p1[i][det][e_type].fit,pname)
                skew_frac[d,i,p] = pd_ecal_p1[i][det][e_type].fit[pname].skew_fraction
                skew_width[d,i,p] = pd_ecal_p1[i][det][e_type].fit[pname].skew_width
                # µ[d,i,p] = pd_ecal_p1[i][det][e_type].fit[pname].µ #they are only simply energy calibrated 
                # σ[d,i,p] = pd_ecal_p1[i][det][e_type].fit[pname].σ
                # fwhm[d,i,p] = pd_ecal_p1[i][det][e_type].fit[pname].fwhm
                µ_ADC  = pd_ecal_p1[i][det][e_type].fit[pname].µ ./ pd_ecal_p1[i][det][e_type].m_cal_simple
                σ_ADC = pd_ecal_p1[i][det][e_type].fit[pname].σ ./ pd_ecal_p1[i][det][e_type].m_cal_simple
                fwhm_ADC = pd_ecal_p1[i][det][e_type].fit[pname].fwhm ./ pd_ecal_p1[i][det][e_type].m_cal_simple
                Tbl = Table(e_cusp = [µ_ADC, σ_ADC, fwhm_ADC], qdrift = [0,0,0])
                cal_func = ljl_propfunc(pd_ecal_p1[i][det][e_type].cal.func)
                µ[d,i,p] = collect(cal_func.(Tbl))[1] # calibrated peak position for period given run and detector 
                σ[d,i,p] = collect(cal_func.(Tbl))[2] .- pd_ecal_p1[i][det][e_type].cal.par[1] # calibrated (w/o offset)
                fwhm[d,i,p] = collect(cal_func.(Tbl))[3] .- pd_ecal_p1[i][det][e_type].cal.par[1] # calibrated (w/o offset)

                background[d,i,p] = pd_ecal_p1[i][det][e_type].fit[pname].background
                n[d,i,p] = pd_ecal_p1[i][det][e_type].fit[pname].n
                step_amplitude[d,i,p] = pd_ecal_p1[i][det][e_type].fit[pname].step_amplitude
                pvalues_peakfit[d,i,p] = pd_ecal_p1[i][det][e_type].fit[pname].gof.pvalue 
            else
                skew_frac[d,i,p] = NaN ± NaN
                skew_width[d,i,p] = NaN ± NaN
                µ[d,i,p] = NaN* u"keV"  .± NaN  * u"keV"
                σ[d,i,p] = NaN* u"keV"  .± NaN  * u"keV"
                fwhm[d,i,p] = NaN* u"keV"  .± NaN  * u"keV"
                background[d,i,p] = NaN ± NaN
                n[d,i,p] = NaN ± NaN
                step_amplitude[d,i,p] = NaN ± NaN
                pvalues_peakfit[d,i,p] = NaN
            end
        end
    end
end  

##########################################################################################################################################################################################################################################
# look at outliers == fit with p < 1 %
th228_literature =sort([pd_ecal_p1[1][Symbol(dets_ged[1])][e_type].cal.peaks..., 1592u"keV", 2103u"keV"]) # get literature values mit denen gefittet wurde. interpolation st Qbb

# sort 
IdxSort = sortperm(µ[:1,1,:])
µ_sort= µ[:, :, IdxSort]
σ_sort= σ[:, :, IdxSort]
fwhm_sort= fwhm[:, :, IdxSort]
skew_frac_sort = skew_frac[:, :, IdxSort]
pval_sort = pvalues_peakfit[:, :, IdxSort]
residual = permutedims(permutedims(µ_sort,(3,1,2)) .- th228_literature,(2,3,1))
residual_norm = mvalue.(residual) ./ muncert.(residual)
pval = reshape(pval_sort,:)

# plot
modes = [:norm, :keV, :sigma, :fwhm, :skew_frac]
for Mode in modes 

if Mode == :norm
    data = reshape(residual_norm,:)
    stat_x = 14
    xl = (-15,15)
    xlbl = L"\textit{E}_\textrm{cal.} - \textit{E}_\textrm{lit.} (σ)"
    xunit = "σ"
    nbins = 300
    legloc = :topleft
elseif Mode == :keV
    data = reshape(ustrip.(mvalue.(residual)),:)
    stat_x = 0.9
    xl = (-1,1)
    xlbl = L"\textit{E}_\textrm{cal.} - \textit{E}_\textrm{lit.}" * "(keV)"
    xunit = "keV"
    nbins = 400
    legloc = :topleft
elseif Mode ==:sigma
    σ_mean =  [median(filter(!isnan,σ_sort[:,:,peak])) for peak in eachindex(th228_literature)]
    σ_res = permutedims(permutedims(σ_sort,(3,1,2)) .- σ_mean,(2,3,1))  
    data = reshape(ustrip.(mvalue.(σ_res)),:)
    stat_x = 1
    xl = (-1,2)
    xlbl = L"\sigma_\textrm{peak} - \langle \sigma_\textrm{peak} \rangle" * " (keV)"
    xunit = "keV"
    nbins = 300
    legloc = :topright
elseif Mode ==:fwhm
    fwhm_mean =  [median(filter(!isnan,fwhm_sort[:,:,peak])) for peak in eachindex(th228_literature)]
    fwhm_res = permutedims(permutedims(fwhm_sort,(3,1,2)) .- fwhm_mean,(2,3,1))  
    data = reshape(ustrip.(mvalue.(fwhm_res)),:)
    stat_x = 1
    xl = (-2,3.5)
    xlbl = L"FWHM  - \langle FWHM_\textrm{peak} \rangle" * " (keV)"
    xunit = "keV"
    nbins = 300
    legloc = :topright
elseif Mode ==:skew_frac
    skew_frac_mean =  [median(filter(!isnan,skew_frac_sort[:,:,peak])) for peak in eachindex(th228_literature)]
    skew_frac_res = permutedims(permutedims(skew_frac_sort,(3,1,2)) .- skew_frac_mean,(2,3,1))  
    data = reshape(ustrip.(mvalue.(skew_frac_res)),:)
    stat_x = 1
    xl = (-0.1,0.25)
    xlbl = L"\textrm{tail}_\textrm{frac} - \langle \textrm{tail}_\textrm{frac} \rangle"
    xunit = ""
    nbins = 300
    legloc = :topright
end
bin_center, bin_edges, _, counts = get_histstats(ustrip.(mvalue.(data)); nbins = nbins)

fs = 14
HistArg = Dict(:normalize => :none,
                :fill => :true,
                :foreground_color_legend => :silver,
                :background_color_legend => :white,
                :grid => :off,
                :xguidefontsize => fs, :xtickfontsize => fs-2,
                :yguidefontsize => fs, :ytickfontsize => fs-2,
                :legendfontsize => fs-2)

plt = stephist(data, bins = bin_edges, 
                color = :darkgrey, 
                fillalpha = 0.5, 
                framestyle = :box, 
                legend = legloc,
                ylims = (0,:auto),
                label = "All peaks",
                ylabel = "Occurrence",
                xlabel  = xlbl,
                xlims = xl, 
                ; HistArg...)
              
stephist!(data[pval .< 0.01], 
                bins = bin_edges, 
                color = :red, 
                fillalpha = 0.3, 
                label = "All peaks with p < 0.01",
                ; HistArg...)

plot!(plt,title = "Calibration peak fits, partition $partition, $e_type, all $(length(dets_ged)) dets", 
        titlefontsize = fs-4, dpi = 300)

plt_name = path_plot * "Ecal_PeakFits_FitParOutliers_part$(partition)_$(e_type)_$Mode.png"
savefig(plt_name)
@info "save plot to $plt_name"
end

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
                