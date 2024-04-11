#= comparison µ peak fits vs literature values of calibration peaks 
plot 2:
- residual histogram: µ - literature . for all peaks together
plot 1:
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
using StatsPlots
using Colors
include("../SanityPlots/utils.jl")

path_plot = "PackageDevScripts/DevelopSpecFits/EnergyCalibration/plots"
path_rel_plot = relpath(path_plot,pwd())
path_abs_plot = pwd() * "/" * path_rel_plot * "/"


l200 = LegendData(:l200)

#select data and dsp output 
partition = 1
e_type = :e_cusp_ctc

# open data
partinfo = partitioninfo(l200)[DataPartition(partition)]
filekey = start_filekey(l200, (partinfo[1].period, partinfo[1].run, :cal)) 
chinfo = Table(channelinfo(l200, filekey; system=:geds, only_processable=true))
dets_ged = chinfo.detector
energy_config = dataprod_config(l200).energy(filekey)

# dets_icpc = dets_ged[chinfo.det_type .== :icpc]

# load all ProbDicts (for all period-run combination in selected partition). speed up load probdict for partition 
pd_ecal_p1 = [l200.par.rpars.ecal[entr.period, entr.run] for entr in partinfo] # takes a while 
nruns = length(pd_ecal_p1)
th228_literature = sort(pd_ecal_p1[1][Symbol(dets_ged[1])][e_type].cal.peaks) # get literature values mit denen gefittet wurde. interpolation st Qbb
th228_names = Symbol.(keys(pd_ecal_p1[1][Symbol(dets_ged[1])][e_type].fit))
th228_vals = [pd_ecal_p1[1][Symbol(dets_ged[1])][e_type].fit[th228_names[i]].µ for i in eachindex(th228_names)]
th228_val_idx = sortperm(th228_vals)
th228_vals = th228_vals[th228_val_idx]
th228_names  = th228_names[th228_val_idx]
th228_names = filter(x-> x != :Tl208DEP && x != :Tl208SEP,th228_names)

Qbb = 2039.0u"keV"
npeaks = length(th228_literature)

# all calibration fit results for a specific detector   
µ_all = ones(length(dets_ged), nruns, npeaks) .* NaN * 1u"keV" .± NaN * 1u"keV"
pvalues = zeros((length(dets_ged), nruns))
for i = 1:nruns
    for (d, det) in enumerate(dets_ged)
        energy_config_ch = merge(energy_config.default, get(energy_config, det, PropDict()))
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
skew_frac = ones(length(dets_ged), nruns, npeaks) .* NaN  .± NaN  
for i = 1:nruns
    for (d, det) in enumerate(dets_ged)  
        for (p, pname) in enumerate(th228_names)
            if haskey(pd_ecal_p1[i][det][e_type].fit,pname)
                skew_frac[d,i,p] = pd_ecal_p1[i][det][e_type].fit[pname].skew_fraction
            else
                skew_frac[d,i,p] = NaN ± NaN
            end
        end
    end
end  


# prepare data for plot 
residual = permutedims(permutedims(µ_all,(3,1,2)) .- th228_literature,(2,3,1))
residual_norm = mvalue.(residual) ./ muncert.(residual)
colors = [:dodgerblue, :darkorange, :lightseagreen, :crimson, :violet]

i=1 
d=1
# common plot arguments 
fs = 14
colors = [:dodgerblue, :darkorange, :lightseagreen, :crimson, :violet]
HistArg = Dict(:foreground_color_legend => :silver,
               :background_color_legend => :white,
               :legend => :left,
               :grid => :on,
               :xguidefontsize => fs, :xtickfontsize => fs-2,
               :yguidefontsize => fs, :ytickfontsize => fs-2,
               :legendfontsize => fs-2,
               :alpha => 0.1,
               :markerstrokewidth => 0)
Mode = :keV #:norm#
if Mode == :norm
    data = residual_norm
    stat_x = 14
    xl = (-15,15)
    xlbl = L"\textit{E}_\textrm{cal.} - \textit{E}_\textrm{lit.} (σ)"
    xunit = "σ"
elseif Mode == :keV
    data = ustrip.(mvalue.(residual))
    stat_x = 0.9
    xl = (-1,1)
    xlbl = L"\textit{E}_\textrm{cal.} - \textit{E}_\textrm{lit.}" * "(keV)"
    xunit = "keV"
end


pall = scatter(mvalue.(reshape(data[:,:,:],:)),mvalue.(reshape(skew_frac[:,:,:],:)),#,nbins=(200, 20),
            color= :grey, # cgrad(:greys), colorbar = false, colorbar_title = "Occurrence";
            label = "All peaks",
            xlabel = xlbl ,
            ylabel = "Fraction in low-energy tail";
            HistArg...)
vline!([0],label = :none, color = :red, linestyle = :dash, linewidth = 3)
p = Vector(undef,npeaks)
for i = 1:npeaks
    rgb_color = parse(Colorant, colors[i])
    lbl = "E = $(round(typeof(1u"keV"), th228_literature[i], digits = 0))" 
    p[i] = scatter(mvalue.(reshape(data[:,:,i],:)),mvalue.(reshape(skew_frac[:,:,i],:)),#,nbins=(200, 20), 
            color=colors[i],#cgrad([rgb_color + 0.5*RGB(1,1,1), rgb_color, rgb_color + 0.3*RGB(0,0,0)]),
            label = lbl,
            yaxis= true, yticks = ([0,0.1,0.2]),#colorbar = false,  ylims = (0,:auto);
            ;HistArg...)
            vline!([0],label = :none, color = :red, linestyle = :solid, linewidth = 3)
    if i==npeaks
        plot!(p[i], xlabel  = xlbl )
    else
        # myxticks = xticks(p[i])
        # Set the labels to empty strings
        plot!(p[i],bottom_margin = -5mm, xformatter=_->"")#, xticks = (myxticks, fill(" ", length(myxticks))))
    end
    vline!([0],label = :none, color = :red, linestyle = :dash, linewidth = 3)
end
l = @layout([a{0.35h}; grid(5, 1)])
ptot = plot(pall,p...,layout =l, size = (650,1000),
     xlims = xl, left_margin = 10mm, right_margin = 5mm)

fname = path_abs_plot * "Ecal_Residuals_$(Mode)_LWfrac_part$(partition)_$(e_type).png"
savefig(ptot,fname)

