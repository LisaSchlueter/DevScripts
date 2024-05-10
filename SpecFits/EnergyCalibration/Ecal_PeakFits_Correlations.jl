#= 
fit paramater of peaks fit (calibration data)
distribution and correlation overview of actual fit parameters 
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

path_plot = "$(@__DIR__)/plots/p$partition/FitPar/"
if !ispath(path_plot)
    mkdir("$path_plot")
end 

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
# pd_ecal_p1 = [l200.par.rpars.ecal[entr.period, entr.run] for entr in partinfo] # takes a while 
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

#####################################################################################################################
# do big correlation plot 
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
               :xguidefontsize => fs, #:xtickfontsize => fs-2,
               :xlabelfontsize => fs + 20, :ylabelfontsize => fs + 20,
               :yguidefontsize => fs, #:ytickfontsize => fs-2,
               :legendfontsize => fs + 20,
               :markerstrokewidth => 0)
default(fontfamily="Helvetica")
th228_literature = sort(pd_ecal_p1[1][Symbol(dets_ged[1])][e_type].cal.peaks) 
for p in eachindex(th228_names)
    cols = columnnames(FitPars[th228_names[p]])
    par_names  =  String.(collect(cols))
    par_names[par_names .== "step_amplitude"] .= "bkg step" #L"\mathrm{bkg}_\mathrm{step}"
    par_names[par_names .== "background"] .= "bkg"
    par_names[par_names .== "skew_width"] .= "tail skew"#L"\textrm{tail}_\textrm{skew}"
    par_names[par_names .== "skew_frac"] .= "tail frac"#L"\textrm{tail}_\textrm{frac}"

    parlims = Vector{Tuple}(undef,length(cols))
    plt = Matrix{Plots.Plot}(undef,length(cols),length(cols))
    for i in eachindex(cols)
        for j in eachindex(cols)
            x = getproperty(FitPars[th228_names[p]], cols[i])
            y = getproperty(FitPars[th228_names[p]], cols[j])
            corcoef = round(cor(filter(!isnan,x),filter(!isnan,y)),digits = 2)
        if i==j 
                plt[i,i] = stephist(x,fill = true, nbins = 100, yaxis = false, legend = false,
                                    color = get(ColorSchemes.RdBu_9, 1, (-1,1)) ,  xformatter=_->"", # xlabel = cols[i]
                                    ; PltArg...)
                    parlims[i] = xlims()         
        elseif i>j
            malpha = 0.3
            if corcoef<0.05
                clr = :silver
            elseif corcoef<0.2
                clr = get(ColorSchemes.RdBu_9, corcoef, (-1,1)) -  0.2*RGB(1,1,1)
            else
                clr = get(ColorSchemes.RdBu_9, corcoef, (-1,1))
            end

                plt[i,j] = scatter(x,y, 
                                label = "ρ = $(round(cor(x,y), digits = 2))", legend = :false,
                                alpha = 0.3, color = clr, ms = 8,
                                xformatter=_->"", yformatter=_->""
                                ; PltArg...)
        else 
                if abs(corcoef) < 0.01
                    corcoef = 0.0 
                end   #(RGB(1,1,1) - abs(corcoef)*RGB(1,1,1))
                plt[i,j] = heatmap(fill(abs(corcoef),10,10), colorbar = false, c = get(ColorSchemes.RdBu_9, corcoef, (-1,1)), framestyle = :box, xformatter=_->"", yformatter=_->"") # axis = false
                annotate!(5,5.5, text("ρ = $corcoef", :black, :center, fs + 30, "Helvetica"))
        end
            if i==1
                plot!(plt[1,j], ylabel = par_names[j],  ylabelfontsize = fs + 40)
            end
            if j == length(cols)
                plot!(plt[i,j], xlabel = par_names[i],  xlabelfontsize = fs + 40)
            end
        end
    end 
    for i in eachindex(cols)
        for j in eachindex(cols)
            if i>j
                plot!(plt[i,j], xlims = parlims[i], ylims = parlims[j])
            end
            if i<length(cols)# everything but last column 
                plot!(plt[i,j], right_margin = -40mm)
            else
                plot!(plt[i,j], left_margin = -6mm)
            end
            if j>1
                plot!(plt[i,j], top_margin = -20mm) #-25
            end
            if i==1
                plot!(plt[i,j], left_margin = 25mm)
            end
            if j==length(cols)
                plot!(plt[i,j], bottom_margin = 7mm)
            end
        end
    end 
   ptot =  plot(plt..., layout = (length(cols),length(cols)), size = (5000,4700), 
                plot_title = "Partition $partition, $(e_type), $(th228_names[p])", plot_titlefontsize = 60, dpi = 300)

    fname = path_plot * "Ecal_FitParCorr_$(th228_names[p])_part$(partition)_$(e_type).png"
     savefig(ptot,fname)
    @info "save plot to $fname"
 end
