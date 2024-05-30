#= comparison µ peak fits vs literature values of calibration peaks 
plot :
- histogram fit selected fit parameter:
- over view for all peaks and breakdown  
- for all peaks and breakdown of single peaks
- µ, \sigma  and fwhm are (re-)calibrated to keV (not just simple calibrated)

--> for all dets or only ICPC, Coax, Bege,....
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
cal = true
FitPars, MetaData = get_peakfit_rpars_partition(DataPartition(partition); reload = false, e_type = e_type, cal = cal);

# plotting path 
path_plot = "$(@__DIR__)/plots/p$partition/FitPar/"
if !ispath(path_plot)
    mkdir("$path_plot")
end 

# plot, common args and labels 
Tbl_lbl = ["E = $(round(typeof(1u"keV"), MetaData.th228_literature[i], digits = 0)) ($(MetaData.th228_names[i]))" for i=1:MetaData.npeaks]
colors = [:dodgerblue, :darkorange, :lightseagreen, :midnightblue, :crimson , :olive ,:violet]
fs = 14
HistArg = Dict(:normalize => :probability,
               :fill => :true,
               :foreground_color_legend => :silver,
               :background_color_legend => :white,
               :grid => :off,
               :xguidefontsize => fs+1, :xtickfontsize => fs-2,
               :yguidefontsize => fs+1, :ytickfontsize => fs-2,
               :legendfontsize => fs-2)

modes = [:fwhm, :σ]#, :n, :background, :step_amplitude, :skew_frac, :skew_width, :p_value]
det_types_list = [:all, unique(MetaData.dets_type)...]
det_mode = :all 

for Mode in modes
    nbins = 400
    pltLine = true 
    xl = (nothing ,nothing)

    if Mode == :fwhm
        data = ustrip.(mvalue.(FitPars.fwhm))
        xl = (0, maximum(filter(!isnan,reshape(data,:))) )
        maximum(xl) > 8 ? xl = (0,8) : nothing
        xlbl = "Peak fwhm"
        xunit = " (keV)"
    elseif Mode ==:σ
        data = ustrip.(mvalue.(FitPars.σ))
        xlbl = "Peak width σ"
        xl = (0, maximum(filter(!isnan,reshape(data,:))) )
        maximum(xl) > 4 ? xl = (0,4) : nothing
        xunit = " (keV)"
    elseif Mode ==:n
        data = ustrip.(mvalue.(FitPars.n))
        xlbl = "Peak amplitude n"
        xl = (0,round(median(filter(isfinite,reshape(data,:))), digits = 0)*5)
        xunit = ""
    elseif Mode ==:background
        data = ustrip.(mvalue.(FitPars.background))
        xlbl = "Background"
        nbins = 100
        xunit = " (per bin)"
    elseif Mode ==:step_amplitude
        data = ustrip.(mvalue.(FitPars.step_amplitude))
        xlbl = "Background step"
        xunit = ""
    elseif Mode ==:skew_frac
        data = ustrip.(mvalue.(FitPars.skew_frac))
        xlbl = "Low-energy tail fraction"
        xl = (0,0.1)
        nbins = 100
        xunit = ""
    elseif Mode ==:skew_width
        data = ustrip.(mvalue.(FitPars.skew_width))
        xlbl = "Low-energy tail width scaling"
        xl = (0, quantile(filter(isfinite,reshape(data,:)),0.95))
        xunit = ""
        nbins = 100
    elseif Mode == :p_value
        data = ustrip.(mvalue.(FitPars.pvalue))
        xl = (0,1)
        xlbl = "p-value"
        xunit = ""
        nbins = 25
    end
    all(isnothing.(xl)) ? xl = (minimum(filter(!isnan,reshape(data,:))), maximum(filter(!isnan,reshape(data,:))) )  : nothing 

   #for det_mode in det_types_list
        if det_mode == :all
            data_det = data
            nbins_det = nbins 
            path_plot_det = path_plot
        elseif Mode !== :p_value
            det_idx = findall(MetaData.dets_type .== det_mode)
            data_det = data[det_idx,:,:]
            nbins_det = round(typeof(1),nbins * length(det_idx)/length(MetaData.dets_ged) )
            @info "det_mode: $det_mode, nbins: $nbins_det"
            path_plot_det = path_plot * "$(det_mode)/"
            if !ispath(path_plot_det)
                mkdir("$path_plot_det")
            end
        end
        if Mode == :skew_width
            nbins_det = minimum(xl):0.0001:maximum(xl)
        end
        pall = stephist(reshape(data_det,:), nbins = nbins_det, 
                    color = :darkgrey, fillalpha = 0.5, 
                    framestyle = :box, 
                    legend = false,
                    ylims = (0,:auto),
                    label = false,
                    ylabel = "\n Occurrence (a.u.)",
                    xlabel  = xlbl * xunit,
                    xlims = xl,
                    ; HistArg...)
        plot!(pall,  yformatter=_->"")
        median_all = median(filter(isfinite,reshape(data,:)))
        lbl_str = "All peaks \n" * "median = " * @sprintf("%.2f ",median_all) * xunit
        lbl_str_xy = [maximum(xl) .- 0.05*diff([xl...]), maximum(ylims()) .- 0.15*diff([ylims()...])]
        annotate!(lbl_str_xy[1],lbl_str_xy[2], text(lbl_str, :right, fs-2, :grey))

        medians = [median(filter(isfinite, reshape(data[:,:,i],:))) for i = 1:MetaData.npeaks]
        if Mode == :p_value
            hline!([0.05], color = :black, label = "Uniform (expection)", linestyle = :dash, linewidth = 1.5)
        end

        p = Vector(undef,MetaData.npeaks)
        for i = 1:MetaData.npeaks
            lbl_str_p = "$(MetaData.th228_names[i]) ($(round(typeof(1u"keV"), MetaData.th228_literature[i], digits = 0)))\n"  * "median = " * @sprintf("%.2f ",medians[i]) * xunit
        
            p[i] = stephist(reshape(data_det[:,:,i],:), nbins = nbins_det, 
                    color = colors[i], 
                    fillalpha = 0.6,  
                    legend = false,
                    yaxis = false , 
                    ylims = (0,:auto);
                    HistArg...)  

            xlims!(xl)   
            annotate!(lbl_str_xy[1], maximum(ylims()) .- 0.25*diff([ylims()...]), text(lbl_str_p, :right, fs-2, colors[i])) 
  
            if i==MetaData.npeaks
                plot!(p[i], xlabel  =  xlbl * xunit)
            else
                plot!(p[i],bottom_margin = -7mm,ylabel = " \n ", yformatter=_->"") #xformatter=_->""
            end
        end
        
        l = @layout([a{0.15h}; grid(MetaData.npeaks+1, 1)])
        ptot = plot(pall,p..., layout = l,
                        plot_title = "Calibration peak fits \npartition $(partition), $(e_type), $(det_mode) dets",
                        plot_titlefontsize = fs-2,
                        size = (450,1500), left_margin = 3mm, right_margin = 5mm, top_margin = 7mm)
            
        if cal == true
            fname = path_plot_det * "Ecal_PeakFit_$(det_mode)_$(Mode)_part$(partition)_$(e_type).png"
        elseif cal == false
            fname = path_plot_det * "Ecal_PeakFit_$(det_mode)_$(Mode)_part$(partition)_$(e_type)_simplecal.png"
        end
        savefig(ptot,fname)
        @info "save plot to $fname"
    #end
end