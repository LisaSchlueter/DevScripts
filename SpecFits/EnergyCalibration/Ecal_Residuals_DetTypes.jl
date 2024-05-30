#= energy-calibration residuals plot 
residuals = difference of calibrated peak positions (fit) ans literature values 
comparison of residuals (all peaks) between the different detector types
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
cal = true # calibration data

FitPars, MetaData = get_peakfit_rpars_partition(DataPartition(partition); reload = false, e_type = e_type, cal = cal);
#plot only those used in calibration fit 
CalPeakIdx = findall([(MetaData.th228_names[peak] !==  :Tl208DEP) && (MetaData.th228_names[peak] !==  :Tl208SEP) for peak in eachindex(MetaData.th228_names)])
npeaks = length(CalPeakIdx)
th228_literature = MetaData.th228_literature[CalPeakIdx]

residual = reshape(FitPars.residual[:,:,CalPeakIdx], length(MetaData.dets_ged),:)
residual_norm = reshape(mvalue.(residual) ./ muncert.(residual),length(MetaData.dets_ged),:) # throw all runs/peaks together 


path_plot = "$(@__DIR__)/plots/p$partition/Residuals/"
if !ispath(path_plot)
    mkpath("$path_plot")
end 

# plot settings
Mode = :keV
fs = 14 # font size 
colors = [:dodgerblue, :darkorange, :lightseagreen, :crimson, :violet]
HistArg = Dict(:normalize => :none,
               :fill => :true,
               :foreground_color_legend => :silver,
               :background_color_legend => :white,
               :grid => :off,
               :xlims => (-10,10),
               :xguidefontsize => fs+3, :xtickfontsize => fs-2,
               :yguidefontsize => fs+2, :ytickfontsize => fs-2,
               :legendfontsize => fs-2)

# labels 
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
data_det_types = [:bege, :icpc, :ppc, :coax]
data_det = [data[findall(MetaData.dets_type .== data_det_types[1]),:], 
            data[findall(MetaData.dets_type .==  data_det_types[2]),:], 
            data[findall(MetaData.dets_type .==  data_det_types[3]),:], 
            data[findall(MetaData.dets_type .==  data_det_types[4]),:]]

# plot overview for all detectors : histograms below each other 
bin_center, bin_edges, _, counts = get_histstats(ustrip.(mvalue.(data)); nbins = nbins)
median_tot = median(filter(isfinite,reshape(data,:)))
fwhm_tot = get_fwhm(bin_center,counts)
pall = stephist(reshape(data,:), bins = bin_edges, 
                color = :darkgrey, fillalpha = 0.5, 
                framestyle = :box, 
                legend = :topleft,
                ylims = (0,:auto),
                label = "All detectors, \nall peaks",
                ylabel = "Occurrence",
                xlabel  = xlbl 
                ; HistArg...)

stat_str_all = "median = " * @sprintf("%.2f%s\n",ustrip(mvalue(median_tot)),xunit) *  @sprintf("fwhm = %.2f%s",fwhm_tot,xunit)
annotate!(stat_x, ylims()[2]*0.9, text(stat_str_all, :right, fs-2, :grey))
vline!([0],label = :none, color = :black, linestyle = :dash, linewidth = 2)
medians = fill(0.0, npeaks)
fwhm = fill(0.0, npeaks)
for i in eachindex(data_det)
    bin_center, counts = nothing, nothing
    bin_center, _ , _, counts = get_histstats(mvalue.(ustrip.(data_det[i])); nbins = 300)
    medians[i] = median(filter(isfinite, mvalue.(ustrip.(reshape(data_det[i],:)))))
    fwhm[i] = get_fwhm(bin_center,counts)
end
p = Vector(undef, length(data_det_types))
for i in eachindex(data_det_types)
    stat_str = "median = " * @sprintf("%.2f%s\n",medians[i],xunit) *  @sprintf("fwhm = %.2f%s",fwhm[i],xunit)
    lbl = uppercase("$(data_det_types[i])")
    p[i] = stephist(reshape(data_det[i],:), bins = bin_edges, 
            color = colors[i], 
            fillalpha = 0.8,  
            legend = :left,
            yaxis=false, 
            ylims = (0,:auto),
            ylabel = " ",
            label=lbl;
            HistArg...)  
    annotate!(stat_x, ylims()[2]*0.5, text(stat_str, :right, fs-2, colors[i]))
    vline!([0],label = :none, color = :black, linestyle = :dash, linewidth = 2)
    if i == length(data_det_types)
        plot!(p[i], xlabel  = xlbl )
    else
        # myxticks = xticks(p[i])
        # Set the labels to empty strings
        #plot!(p[i],bottom_margin = -7mm, xformatter=_->"")#, xticks = (myxticks, fill(" ", length(myxticks))))
    end
end
l = @layout([a{0.25h}; grid(5, 1)])
ptot = plot(pall,p...,layout =l, size = (650,1200), xlims = xl, left_margin = 10mm, right_margin = 5mm)

if cal == true
    fname = path_plot * "Ecal_Residuals_DetTypes_$(Mode)_part$(partition)_$(e_type).png"
elseif cal == false
    fname = path_plot * "Ecal_Residuals_DetTypes_$(Mode)_part$(partition)_$(e_type)_simplecal.png"
end
savefig(ptot,fname)
@info "save plot to $fname"



