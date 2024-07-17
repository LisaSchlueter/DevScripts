#= 
fit paramater of peaks fit (calibration data)
look peak fits with bad pvalue. look at fit parameter , especially µ and \sigma 
different "modes" available: modes = [:norm, :keV, :sigma, :fwhm, :skew_frac]
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

#select data and dsp output
partition = 1
e_type = :e_cusp_ctc
FitPars, MetaData = get_peakfit_rpars_partition(DataPartition(partition); reload = false, e_type = e_type);

# plotting path 
path_plot = "$(@__DIR__)/plots/p$partition/FitPar/"
if !ispath(path_plot)
    mkpath("$path_plot")
end 

# plot
modes = [:norm, :keV, :sigma, :fwhm, :skew_frac]
for Mode in modes 

    if Mode == :norm
        data = reshape(mvalue.(FitPars.residual) ./ muncert.(FitPars.residual),:)
        stat_x = 14
        xl = (-15,15)
        xlbl = L"\textit{E}_\textrm{cal.} - \textit{E}_\textrm{lit.} (σ)"
        xunit = "σ"
        nbins = 300
        legloc = :topleft
    elseif Mode == :keV
        data = reshape(ustrip.(mvalue.(FitPars.residual)),:)
        stat_x = 0.9
        xl = (-1,1)
        xlbl = L"\textit{E}_\textrm{cal.} - \textit{E}_\textrm{lit.}" * "(keV)"
        xunit = "keV"
        nbins = 400
        legloc = :topleft
    elseif Mode ==:sigma
        σ_mean =  [median(filter(!isnan,FitPars.σ[:,:,peak])) for peak in eachindex(MetaData.th228_literature)]
        σ_res = permutedims(permutedims(FitPars.σ,(3,1,2)) .- σ_mean,(2,3,1))  
        data = reshape(ustrip.(mvalue.(σ_res)),:)
        stat_x = 1
        xl = (-1,2)
        xlbl = L"\sigma_\textrm{peak} - \langle \sigma_\textrm{peak} \rangle" * " (keV)"
        xunit = "keV"
        nbins = 300
        legloc = :topright
    elseif Mode ==:fwhm
        fwhm_mean =  [median(filter(!isnan,FitPars.fwhm[:,:,peak])) for peak in eachindex(MetaData.th228_literature)]
        fwhm_res = permutedims(permutedims(FitPars.fwhm,(3,1,2)) .- fwhm_mean,(2,3,1))  
        data = reshape(ustrip.(mvalue.(fwhm_res)),:)
        stat_x = 1
        xl = (-2,3.5)
        xlbl = L"FWHM  - \langle FWHM_\textrm{peak} \rangle" * " (keV)"
        xunit = "keV"
        nbins = 300
        legloc = :topright
    elseif Mode ==:skew_frac
        skew_frac_mean =  [median(filter(!isnan,FitPars.skew_frac[:,:,peak])) for peak in eachindex(MetaData.th228_literature)]
        skew_frac_res = permutedims(permutedims(FitPars.skew_frac,(3,1,2)) .- skew_frac_mean,(2,3,1))  
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

    plot!(plt,title = "Calibration peak fits, partition $partition, $e_type, all $(length(MetaData.dets_ged)) dets", 
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
                