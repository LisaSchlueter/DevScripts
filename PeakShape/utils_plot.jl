
function _get_def_peakcolors()
    def_colors = Dict(:Tl208a => :dodgerblue, 
            :Bi212a => :darkorange, 
            :Tl208b => :lightseagreen,
            :Tl208DEP => :midnightblue,
            :Bi212FEP => :crimson,
            :Tl208SEP => :olive,
            :Tl208FEP => :violet)
    return def_colors
end

function _get_def_histargs(; fs = 14)
    HistArg = Dict(:normalize => :probability,
               :fill => :true,
               :foreground_color_legend => :silver,
               :background_color_legend => :white,
               :grid => :off,
               :xguidefontsize => fs+1, :xtickfontsize => fs-2,
               :yguidefontsize => fs+1, :ytickfontsize => fs-2,
               :legendfontsize => fs-2)
    return HistArg
end

function _def_plot(; fs::Integer = 14)
    default(foreground_color_legend = :silver,
    background_color_legend = :white,
    grid = :off,
    framestyle = :semi,
    xtickfontsize = fs-2,
    ytickfontsize = fs-2,
    legendfontsize = fs-2,
    xlabelfontsize = fs + 2,
    ylabelfontsize = fs + 2)
end

# define plotting function 
function plothist_diff(FitPars1, FitPars2, MetaData1, MetaData2; par_mode::Symbol = :fwhm, det_mode::Symbol = :all, exclPeaks::Vector = [], saveplot::Bool = false, fontsize ::Int = 14, xlim = nothing, title::String = "")
    # select plot parameter and some plot-parameter specific settings/labels 
    HistArg = _get_def_histargs(;fs = fontsize)
    colors = _get_def_peakcolors()

    ParArg = Dict(:nbins => 400,
                :pltLine => true,
                :xl => (nothing ,nothing))
    if par_mode == :fwhm
        ParArg[:data] = ustrip.(mvalue.(FitPars1.fwhm)) .- ustrip.(mvalue.(FitPars2.fwhm))
        xl = (minimum(filter(!isnan,reshape(ParArg[:data],:))), maximum(filter(!isnan,reshape(ParArg[:data],:))) )
        maximum(xl) > 1 ? xl = (minimum(xl),1) : nothing
        minimum(xl) < -1 ? xl = (-1,maximum(xl)) : nothing

        ParArg[:xl] = isnothing(xlim) ? xl : xlim
        ParArg[:xlbl] = "Δ FWHM "
        ParArg[:xunit] = " (keV)"
    elseif par_mode == :µ
        ParArg[:data] = ustrip.(mvalue.(FitPars1.µ)) .- ustrip.(mvalue.(FitPars2.µ))
        xl = (minimum(filter(!isnan,reshape(ParArg[:data],:))), maximum(filter(!isnan,reshape(ParArg[:data],:))) )
        ParArg[:xlbl] = "Δ Peak position µ "
        ParArg[:xunit] = " (keV)"
    elseif par_mode ==:σ
        ParArg[:data] = ustrip.(mvalue.(FitPars1.σ)) .- ustrip.(mvalue.(FitPars2.σ))
        ParArg[:xlbl] = "Δ Width σ"
        xl = (0, maximum(filter(!isnan,reshape(ParArg[:data],:))) )
        maximum(xl) > 4 ? xl = (0,4) : nothing
        ParArg[:xl] =  isnothing(xlim) ? xl : xlim
        ParArg[:xunit]  = " (keV)"
    elseif par_mode ==:n
        ParArg[:data] = ustrip.(mvalue.(FitPars1.n)) .- ustrip.(mvalue.(FitPars2.n))
        ParArg[:xlbl] = "Δ Amplitude n"
        xl = (0,round(median(filter(isfinite,reshape(ParArg[:data],:))), digits = 0)*5) 
        ParArg[:xl]  = isnothing(xlim) ? xl : xlim
        ParArg[:xunit]  = ""
    elseif par_mode ==:background
        ParArg[:data] = ustrip.(mvalue.(FitPars1.background)) .- ustrip.(mvalue.(FitPars2.background))
        ParArg[:xlbl] = "Δ Background"
        ParArg[:nbins] = 100
        ParArg[:xunit]  = " (per bin)"
    elseif par_mode ==:step_amplitude
        ParArg[:data] = ustrip.(mvalue.(FitPars1.step_amplitude)) .- ustrip.(mvalue.(FitPars2.step_amplitude))
        ParArg[:xlbl] = "Δ Background step"
        ParArg[:xunit]  = ""
    elseif par_mode ==:skew_frac
        ParArg[:data] = ustrip.(mvalue.(FitPars1.skew_frac)) .- ustrip.(mvalue.(FitPars2.skew_frac))
        ParArg[:xlbl] = "Δ Low-energy tail fraction"
        xl = (0,0.1)
        ParArg[:xl] = isnothing(xlim) ? xl : xlim
        ParArg[:nbins] = 100
        ParArg[:xunit]  = ""
    elseif par_mode ==:skew_width
        ParArg[:data] = ustrip.(mvalue.(FitPars1.skew_width)) .- ustrip.(mvalue.(FitPars2.skew_width))
        ParArg[:xlbl] = "Δ Low-energy tail width scaling"
        xl  = (0, quantile(filter(isfinite,reshape(ParArg[:data],:)),0.95)) 
        ParArg[:xl] = isnothing(xlim) ? xl : xlim
        ParArg[:xunit]  = ""
        ParArg[:nbins] = 100
    elseif par_mode == :p_value
        ParArg[:data] = ustrip.(mvalue.(FitPars1.pvalue)) .- ustrip.(mvalue.(FitPars2.pvalue))
        ParArg[:xl] = isnothing(xlim) ? (0,1) : xlim
        ParArg[:xlbl] = "Δ p-value"
        ParArg[:xunit]  = ""
        ParArg[:nbins]  = 100
    elseif par_mode == :chi2
        ParArg[:data] = ustrip.(mvalue.(FitPars1.chi2)) .- ustrip.(mvalue.(FitPars2.chi2))
        # ParArg[:xl] = isnothing(xlim) ? (-100,100) : xlim
        ParArg[:xlbl] = L"Δ $\chi^2$"
        ParArg[:xunit]  = ""
        ParArg[:nbins]  = 100
    end
    if isnothing(xlim)
        all(isnothing.(ParArg[:xl])) ? ParArg[:xl] = (minimum(filter(!isnan,reshape(ParArg[:data],:))), maximum(filter(!isnan,reshape(ParArg[:data],:))) )  : nothing 
    else
        ParArg[:xl] = xlim
    end 
    th228_names = Symbol.(MetaData1.ecal_config.th228_names)
    th228_literature = MetaData1.ecal_config.th228_lines
    # peak and detector-type selection 
    CalPeakIdx = findall(x -> !(x in exclPeaks), th228_names)
    if det_mode == :all
        ParArg[:data_det] = ParArg[:data][:, :, CalPeakIdx]
        ParArg[:nbins_det] = ParArg[:nbins] 
    else
        local det_idx
        det_idx = findall(MetaData.chinfo.det_type .== det_mode)
        ParArg[:data_det] = ParArg[:data][det_idx, :, CalPeakIdx]
        if par_mode == :p_value
            ParArg[:nbins_det] = ParArg[:nbins]  
        else
            ParArg[:nbins_det] = round(typeof(1),ParArg[:nbins]  * length(det_idx)/length(MetaData.dets_ged) )
        end 
        @info "det_mode: $det_mode, nbins: $([:nbins_det]))"
    end
    if par_mode == :skew_width
        ParArg[:nbins_det] = minimum(ParArg[:xl]):0.0001:maximum(ParArg[:xl])
    end

    # plot overview : histograms below each other
    plt_data = filter(x-> minimum(ParArg[:xl]) <= x <= maximum(ParArg[:xl]),reshape(ParArg[:data_det],:))
    pall = stephist(plt_data, nbins = ParArg[:nbins_det], 
                color = :darkgrey, fillalpha = 0.5, 
                framestyle = :box, 
                legend = false,
                ylims = (0,:auto),
                label = false,
                ylabel = "\n Occurrence (a.u.)",
                xlabel  = ParArg[:xlbl] * ParArg[:xunit],
                xlims = ParArg[:xl],
                ; HistArg...)
    plot!(pall,  yformatter=_->"")
    median_all = median(filter(isfinite,reshape(ParArg[:data_det],:)))
    lbl_str = "All peaks \n" * "median = " * @sprintf("%.2g ",median_all) * ParArg[:xunit][3:end-1]
    lbl_str_xy = [maximum(ParArg[:xl]) .- 0.05*diff([ParArg[:xl]...]), maximum(ylims()) .- 0.15*diff([ylims()...])]

    annotate!(lbl_str_xy[1],lbl_str_xy[2], text(lbl_str, :right, fontsize-2, :grey))
    medians = [median(filter(isfinite, reshape(ParArg[:data_det][:,:,i],:))) for i in eachindex(CalPeakIdx)]
    # if par_mode == :p_value
    #     hline!([0.05], color = :black, label = "Uniform (expection)", linestyle = :dash, linewidth = 1.5)
    #     annotate!(0.03, 0.05 + 0.05*diff([ylims()...])[1], text("Uniform expection (5%)", :left, fontsize-2, :black))
    # end

    p = Vector(undef,length(CalPeakIdx))
    for i in eachindex(CalPeakIdx)
        lbl_str_p = "$(th228_names[CalPeakIdx[i]]) ($(round(typeof(1u"keV"), th228_literature[CalPeakIdx[i]], digits = 0)))\n"  * "median = " * @sprintf("%.2g ",medians[i]) * ParArg[:xunit][3:end-1]
        local plt_data = filter(x-> minimum(ParArg[:xl]) <= x <= maximum(ParArg[:xl]),reshape(ParArg[:data_det][:,:,i],:))
        p[i] = stephist(plt_data, nbins = ParArg[:nbins_det], 
                color = colors[th228_names[CalPeakIdx[i]]], 
                fillalpha = 0.6,  
                legend = false,
                yaxis = false , 
                ylims = (0,:auto);
                HistArg...)  

        xlims!(ParArg[:xl])   
        annotate!(lbl_str_xy[1], maximum(ylims()) .- 0.25*diff([ylims()...]), text(lbl_str_p, :right, fontsize-2, colors[th228_names[CalPeakIdx[i]]])) 

        if i==length(CalPeakIdx)
            plot!(p[i], xlabel  =  ParArg[:xlbl] * ParArg[:xunit])
        else
            plot!(p[i],bottom_margin = -7mm,ylabel = " \n ", yformatter=_->"") #xformatter=_->""
        end
    end
    if length(CalPeakIdx) <=5
        l = @layout([a{0.35 * h}; grid(length(CalPeakIdx), 1)])
    else
        l = @layout([a{0.25 * h}; grid(length(CalPeakIdx), 1)])
    end
    ptot = plot(pall, p..., layout = l,
                    size = (650, 180*(length(CalPeakIdx)+1)), left_margin = 3mm, right_margin = 5mm, top_margin = 7mm,
                    plot_title = "Peakshapes: $(MetaData1.cal_type) - $(MetaData2.cal_type) ",#"Calibration peak fits \nperiods $(periods), $(e_type), $(det_mode) dets",
                    plot_titlefontsize = fontsize-2)

    if saveplot == true
        path_plot = "$(@__DIR__)/plots/period$(join(string.(periods)))/FitParDiff/"
        if !ispath(path_plot)
            mkpath("$path_plot")
        end 
        if det_mode != :all
            path_plot = path_plot * "$(det_mode)/"
            if !ispath(path_plot)
                mkpath("$path_plot")
            end
        end        
        
        fname = path_plot * "Ecal_PeakFitDiff_$(det_mode)_$(par_mode)_period$(join(string.(periods)))_$(e_type)_$(MetaData1.cal_type)minus$(MetaData2.cal_type).png"
        savefig(ptot, fname)
        @info "save plot to $fname"
    end
    display(ptot)
end   