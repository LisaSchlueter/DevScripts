#= 
fit paramater of peaks fit (calibration data)
distribution and correlation overview of actual fit parameters 
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
partition = 3
e_type = :e_cusp_ctc
FitPars, MetaData = get_peakfit_rpars_partition(DataPartition(partition); reload = false, e_type = e_type);

path_plot = "$(@__DIR__)/plots/p$partition/FitPar/"
if !ispath(path_plot)
    mkdir("$path_plot")
end 
#####################################################################################################################
# do big correlation plot 
PltPars = Dict()
for i in eachindex(MetaData.th228_names)
    PltPars[MetaData.th228_names[i]] = Table(µ = mvalue.(ustrip.(reshape(FitPars.µ[:,:,i],:))), 
                σ = mvalue.(ustrip.(reshape(FitPars.σ[:,:,i],:))),
                n = mvalue.(reshape(FitPars.n[:,:,i],:)),
                background =  mvalue.(reshape(FitPars.background[:,:,i],:)),
                step_amplitude =  mvalue.(reshape(FitPars.step_amplitude[:,:,i],:)),
                skew_frac =  mvalue.(reshape(FitPars.skew_frac[:,:,i],:)),
                skew_width =  mvalue.(reshape(FitPars.skew_width[:,:,i],:)))
end
cols = [:µ :σ :n :background :step_amplitude :skew_frac :skew_width]
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
# th228_literature = sort(pd_ecal_p1[1][Symbol(dets_ged[1])][e_type].cal.peaks) 
for p in eachindex(MetaData.th228_names)
   # cols = columnnames(PltPars[MetaData.th228_names[p]])
    par_names  =  String.(collect(cols))
    par_names[par_names .== "step_amplitude"] .= "bkg step" #L"\mathrm{bkg}_\mathrm{step}"
    par_names[par_names .== "background"] .= "bkg"
    par_names[par_names .== "skew_width"] .= "tail skew"#L"\textrm{tail}_\textrm{skew}"
    par_names[par_names .== "skew_frac"] .= "tail frac"#L"\textrm{tail}_\textrm{frac}"

    parlims = Vector{Tuple}(undef,length(cols))
    plt = Matrix{Plots.Plot}(undef,length(cols),length(cols))
    for i in eachindex(cols)
        for j in eachindex(cols)
            x = getproperty(PltPars[MetaData.th228_names[p]], cols[i])
            y = getproperty(PltPars[MetaData.th228_names[p]], cols[j])
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
                plot_title = "Partition $partition, $(e_type), $(MetaData.th228_names[p])", plot_titlefontsize = 60, dpi = 300)

    fname = path_plot * "Ecal_FitParCorr_$(MetaData.th228_names[p])_part$(partition)_$(e_type).png"
     savefig(ptot,fname)
    @info "save plot to $fname"
 end
