#= 
fit paramater of peaks fit (calibration data)
distribution and correlation overview 
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
using StatsPlots
using Colors
include("../SanityPlots/utils.jl")

path_plot = "PackageDevScripts/DevelopSpecFits/EnergyCalibration/plots/Ecal_PeakFits_Corr/"
path_rel_plot = relpath(path_plot,pwd())
path_abs_plot = pwd() * "/" * path_rel_plot * "/"

#select data and dsp output
l200 = LegendData(:l200) 
partition = 1
e_type = :e_cusp_ctc

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
for i = 1:nruns
    for (d, det) in enumerate(dets_ged)  
        for (p, pname) in enumerate(th228_names)
            if haskey(pd_ecal_p1[i][det][e_type].fit,pname)
                skew_frac[d,i,p] = pd_ecal_p1[i][det][e_type].fit[pname].skew_fraction
                skew_width[d,i,p] = pd_ecal_p1[i][det][e_type].fit[pname].skew_width
                µ[d,i,p] = pd_ecal_p1[i][det][e_type].fit[pname].µ
                σ[d,i,p] = pd_ecal_p1[i][det][e_type].fit[pname].σ
                fwhm[d,i,p] = pd_ecal_p1[i][det][e_type].fit[pname].fwhm
                background[d,i,p] = pd_ecal_p1[i][det][e_type].fit[pname].background
                n[d,i,p] = pd_ecal_p1[i][det][e_type].fit[pname].n
                step_amplitude[d,i,p] = pd_ecal_p1[i][det][e_type].fit[pname].step_amplitude

            else
                skew_frac[d,i,p] = NaN ± NaN
                skew_width[d,i,p] = NaN ± NaN
                µ[d,i,p] = NaN* u"keV"  .± NaN  * u"keV"
                σ[d,i,p] = NaN* u"keV"  .± NaN  * u"keV"
                fwhm[d,i,p] = NaN* u"keV"  .± NaN  * u"keV"
                background[d,i,p] = NaN ± NaN
                n[d,i,p] = NaN ± NaN
                step_amplitude[d,i,p] = NaN ± NaN
    
            end
        end
    end
end  


FitPars = Dict()
for i in eachindex(th228_names)
FitPars[th228_names[i]] = Table(µ = mvalue.(ustrip.(reshape(µ[:,:,i],:))), 
                σ = mvalue.(ustrip.(reshape(σ[:,:,i],:))),
                fwhm = mvalue.(ustrip.(reshape(fwhm[:,:,i],:))),
                n = mvalue.(reshape(n[:,:,i],:)),
                background =  mvalue.(reshape(background[:,:,i],:)),
                step_amplitude =  mvalue.(reshape(step_amplitude[:,:,i],:)),
                skew_frac =  mvalue.(reshape(skew_frac[:,:,i],:)),
                skew_width =  mvalue.(reshape(skew_width[:,:,i],:)))
end
cols = [:µ :σ :fwhm :n :background :step_amplitude :skew_frac :skew_width]

fs = 30
PltArg = Dict(:foreground_color_legend => :silver,
               :background_color_legend => :white,
               :grid => :off,
               :xguidefontsize => fs, :xtickfontsize => fs-2,
               :xlabelfontsize => fs + 20, :ylabelfontsize => fs + 20,
               :yguidefontsize => fs, :ytickfontsize => fs-2,
               :legendfontsize => fs + 20,
               :markerstrokewidth => 0)


for p in eachindex(th228_names)
    cols = columnnames(FitPars[th228_names[p]])
    parlims = Vector{Tuple}(undef,length(cols))
    plt = Matrix{Plots.Plot}(undef,length(cols),length(cols))
    for i in eachindex(cols)
        for j in eachindex(cols)
            x = getproperty(FitPars[th228_names[p]], cols[i])
            y = getproperty(FitPars[th228_names[p]], cols[j])
        if i==j 
                plt[i,i] = stephist(x,fill = true, nbins = 100, yaxis = false, legend = false,
                                    color = :black , # xlabel = cols[i]
                                    ; PltArg...)
                    parlims[i] = xlims()         
        elseif i>j
                plt[i,j] = scatter(x,y, 
                                label = "ρ = $(round(cor(x,y), digits = 2))", legend = :false,
                                alpha = 0.1, color = :palevioletred1, ms = 7,
                                xformatter=_->"", yformatter=_->""
                                ; PltArg...)
        else 
                corcoef = round(cor(filter(!isnan,x),filter(!isnan,y)),digits = 2)
                if abs(corcoef) < 0.01
                    corcoef = 0.0 
                end
                plt[i,j] = heatmap(fill(abs(corcoef),10,10), colorbar = false, c = (RGB(1,1,1) - abs(corcoef)*RGB(1,1,1)), framestyle = :box, xformatter=_->"", yformatter=_->"") # axis = false
                annotate!(5.5,5.5,text("ρ = $corcoef", :palevioletred1, :center, fs + 30, "Helvetica Bold"))
                
        end
            if i==1
                plot!(plt[1,j], ylabel = cols[j],  ylabelfontsize = fs + 20)
            end
            if j == length(cols)
                plot!(plt[i,j], xlabel = cols[i],  xlabelfontsize = fs + 20)
            end
        end
    end 
    for i in eachindex(cols)
        for j in eachindex(cols)
            if i>j
                plot!(plt[i,j],xlims = parlims[i], ylims = parlims[j])
            end
            if i<length(cols)
                plot!(plt[i,j], right_margin = -30mm)
            else
                plot!(plt[i,j], left_margin = -5mm)
            end
            if j>1
                plot!(plt[i,j], top_margin = -25mm)
            end
            if i==1
                plot!(plt[i,j], left_margin = 25mm)
            end
        
        end
    end 
   ptot =  plot(plt..., layout = (length(cols),length(cols)), size = (5000,4700))#,
                        #   bottom_margin = 5mm, left_margin = 20mm, top_margin = 5mm)

    fname = path_abs_plot * "Ecal_FitParCorr_$(th228_names[p])_part$(partition)_$(e_type).png"
    savefig(ptot,fname)
    @info "save plot to $fname"
end
