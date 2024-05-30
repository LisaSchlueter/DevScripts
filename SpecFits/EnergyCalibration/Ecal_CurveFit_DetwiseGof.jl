#= 
script to evaluate the goodness of fit of the calibration curve for each detector (and each run within partition)
plot:
    x-axis: detecotr id 
    y-axis: different options "modes" available for y-axis 
        mode = :pval: p-value of calibration curve fit
        mode = :maxres : maximum residual in keV
        mode = :maxres_prct : maximum residual in percent
=#
# using Distributions, StatsBase, DataFrames
using LegendDataManagement
using LegendDataManagement.LDMUtils
using LegendSpecFits
using Measures
using Plots, Printf, LaTeXStrings, ColorSchemes
using TypedTables
using Unitful
using Measurements
using Measurements: value as mvalue, uncertainty as muncert
using StatsBase
include("$(@__DIR__)/utils_ecal.jl")

#settings: select data and dsp (energy) output 
partition = 4
e_type = :e_cusp_ctc
# Mode = :maxres_prct 
Mode = :maxres
# Mode = :pval 
CalPars, MetaData = get_calcurvefit_rpars_partition(DataPartition(partition); reload = false, e_type = e_type);

path_plot = "$(@__DIR__)/plots/p$partition/CalibrationCurve/"
if !ispath(path_plot)
    mkdir("$path_plot")
end
 
good_det_Idx = findall(MetaData.good_det .== true)
pvalue = CalPars.pvalue[good_det_Idx,:]
residuals = CalPars.peaks_literature_keV[good_det_Idx,:,:] .- mvalue.(CalPars.y_fit[good_det_Idx,:,:])
residuals_max =  reshape(maximum(abs.(residuals), dims = 3 ), size(pvalue))
residuals_max_prct =  100 .* reshape(maximum(abs.(residuals) ./ CalPars.peaks_literature_keV[good_det_Idx,:,:], dims = 3), size(pvalue)) 

if Mode == :pval
    data = pvalue
    ylbl = "p-value"
elseif Mode == :maxres
    data = residuals_max
    ylbl = "max. residual "
elseif Mode == :maxres_prct
    data = residuals_max_prct
    ylbl = "max. residual (%)"
end


dets_ged  = MetaData.dets_ged[good_det_Idx]
dets_type = MetaData.dets_type[good_det_Idx]
det_types = unique(dets_type)

# some info 
badpMatrix = pvalue .< 0.05 
BadCalFits = pvalue .< 0.05 
BadCalFits_Prct = round(100*sum(BadCalFits)/length(badpMatrix), digits =1)
BadCalFits_bege_Prct_det = [round(100*sum(BadCalFits[dets_type .== dt ,:])/length(badpMatrix[dets_type .== dt ,:]), digits =1) for dt in det_types]
@info "All dets: $(sum(BadCalFits)) out of $(length(badpMatrix)) ($BadCalFits_Prct%) calibration curve fits have p < 0.05"
for (i, dt) in enumerate(det_types)
    @info "$dt: $(sum(BadCalFits[dets_type .== dt ,:])) out of $(length(badpMatrix[dets_type .== dt ,:])) ($(BadCalFits_bege_Prct_det[i])%) calibration curve fits have p < 0.05"
end

# common plot arguments 
fs = 14
PltArg = Dict( :framestyle => :box,
               :foreground_color_legend => :silver,
               :background_color_legend => :white,
               :xguidefontsize => fs, :xtickfontsize => fs-4,
               :yguidefontsize => fs, :ytickfontsize => fs-2,
               :legendfontsize => fs-2, 
               :markerstrokewidth => 0.5,
               :markeralpha => 0.8,
               :markerstrokealpha => 1 )

# sort after detector type
detSortIdx = [findall(dets_type .== :bege); findall(dets_type .== :icpc); findall(dets_type .== :ppc);  findall(dets_type .== :coax)]
plt_data = data[detSortIdx,:]
dets_type = dets_type[detSortIdx]
dets_ged = dets_ged[detSortIdx]

pltcolors = ColorSchemes.distinguishable_colors(MetaData.nruns)
pltmarkers = [:circle, :square, :diamond, :cross, :xcross, :utriangle, :pentagon, :dtriangle, :rtriangle, :ltriangle, :pentagon, :hexagon, :heptagon, :octagon, :star4, :star5, :star6, :star7, :star8]

plt = scatter(1:length(dets_ged), plt_data[:,1], xticks = (1:length(dets_ged),string.(dets_ged)), 
        marker = pltmarkers[1],
        xrotation = 90,
        ylabel = ylbl,
        xlabel = "Detector",
        size = (1500,500),
        legend = :outertop,
        legend_columns = 6,
        bottom_margin = 15mm,
        left_margin = 10mm,
        color = pltcolors[1],
        xlims = (0.5,length(dets_ged)+0.5),
        markerstrokecolor =  pltcolors[1],
        label = "$(MetaData.partinfo[1].period)-$(MetaData.partinfo[1].run)",
        ; PltArg...)

for i = 2:MetaData.nruns
    scatter!(1:length(dets_ged), 
            data[:,i],
            marker = pltmarkers[i],
            color = pltcolors[i],  
            markerstrokecolor =  pltcolors[i],
            label = "$(MetaData.partinfo[i].period)-$(MetaData.partinfo[i].run)",
            ; PltArg...)
end

plot!(plt,title = "Partition $partition: calibration curve p-values", fontsize = fs-4, dpi = 300)

fname_calfit = path_plot * "Ecal_calcurve_$(Mode)_overview_part$(partition)_$(e_type).png"
savefig(plt,fname_calfit)