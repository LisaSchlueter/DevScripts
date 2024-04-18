#= comparison µ peak fits vs literature values of calibration peaks 
plot :
- histogram fit selected fit parameter:
- over view for all peaks and breakdown  
- for all peaks and breakdown of single peaks
- µ, \sigma  and fwhm are (re-)calibrated to keV (not just simple calibrated)
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
include("../SanityPlots/utils.jl")

l200 = LegendData(:l200)

#select data and dsp output 
partition = 1
e_type = :e_cusp_ctc

# open data
partinfo = partitioninfo(l200)[DataPartition(partition)]
filekey = start_filekey(l200, (partinfo[1].period, partinfo[1].run, :cal)) 
chinfo = Table(channelinfo(l200, filekey; system=:geds, only_processable=true))
dets_ged = chinfo.detector
dets_type = chinfo.det_type

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
# prepare for plotting 
th228_literature = sort([pd_ecal_p1[1][Symbol(dets_ged[1])][e_type].cal.peaks..., 1592u"keV", 2103u"keV"]) # get literature values mit denen gefittet wurde. interpolation st Qbb
th228_names = Symbol.(keys(pd_ecal_p1[1][Symbol(dets_ged[1])][e_type].fit))
IdxSort = sortperm(µ[:1,1,:])
th228_names  = th228_names[IdxSort]

# sort according to energy 
µ_sort= µ[:, :, IdxSort]
σ_sort= σ[:, :, IdxSort]
skew_frac_sort = skew_frac[:, :, IdxSort]
skew_width_sort = skew_width[:, :, IdxSort]
n_sort = n[:, :, IdxSort]
background_sort = background[:, :, IdxSort]
step_amplitude_sort = step_amplitude[:, :, IdxSort]
fwhm_sort= fwhm[:, :, IdxSort]
pval_sort = pvalues_peakfit[:, :, IdxSort]
pval = reshape(pval_sort,:)

# plot, common args and labels 
Tbl_lbl = ["E = $(round(typeof(1u"keV"), th228_literature[i], digits = 0)) ($(th228_names[i]))" for i=1:npeaks]
colors = [:dodgerblue, :darkorange, :lightseagreen, :midnightblue, :crimson , :olive ,:violet]
fs = 14
HistArg = Dict(:normalize => :probability,
               :fill => :true,
               :foreground_color_legend => :transparent,
               :background_color_legend => :transparent,
               :grid => :off,
               :xguidefontsize => fs, :xtickfontsize => fs-2,
               :yguidefontsize => fs, :ytickfontsize => fs-2,
               :legendfontsize => fs-2)

modes = [:µ, :fwhm, :σ, :n, :background, :step_amplitude, :skew_frac, :skew_width]
for Mode in modes
#Mode = :fwhm
nbins = 400
pltLine = true 
if Mode == :µ
    data = ustrip.(mvalue.(µ_sort))
    stat_x = 14
    #xl = (-15,15)
    xlbl = "Peak position µ" * " (keV)"
    xunit = "keV"
    nbins = 300
elseif Mode == :fwhm
    data = ustrip.(mvalue.(σ_sort))
    # stat_x = 0.9
    xl = (-1,1)
    xlbl = "Peak fwhm" * " (keV)"
    xunit = "keV"
elseif Mode ==:σ
    data = ustrip.(mvalue.(σ_sort))
    xlbl = "Peak width σ" * " (keV)"
elseif Mode ==:n
    data = ustrip.(mvalue.(n_sort))
    xlbl = "Peak amplitude n"
elseif Mode ==:background
    data = ustrip.(mvalue.(background_sort))
    xlbl = "Background / x keV"
    nbins = 100
elseif Mode ==:step_amplitude
    data = ustrip.(mvalue.(step_amplitude_sort))
    xlbl = "Background step / x keV"
    pltLine = false 
elseif Mode ==:skew_frac
    data = ustrip.(mvalue.(skew_frac_sort))
    xlbl = "Low-energy tail fraction"
    nbins = 100
    pltLine = false 
elseif Mode ==:skew_width
    data = ustrip.(mvalue.(skew_width_sort))
    xlbl = "Low-energy tail width"
    nbins = 100
    pltLine = false 
end

det_types_list = [:all, unique(dets_type)...]
for det_mode in det_types_list

    if det_mode == :all
        data_det = data
        nbins_det = nbins 
    else
        det_idx = findall(dets_type .== det_mode)
        data_det = data[det_idx,:,:]
        nbins_det = round(typeof(1),nbins * length(det_idx)/length(dets_ged) )
        # nbins_min = round(typeof(1),nbins*0.75)
        # nbins_det = ifelse(nbins_det<nbins_min,nbins_min,nbins_det)
        @info "det_mode: $det_mode, nbins: $nbins_det"
    end

    p = Vector(undef,npeaks)
    for i = 1:npeaks
        lbl = "E = $(round(typeof(1u"keV"), th228_literature[i], digits = 0)) \n $(th228_names[i])" 
        p[i] = stephist(reshape(data_det[:,:,i],:), bins = nbins_det, 
                color = colors[i], 
                fillalpha = 0.8,  
                legend = :left,
                yaxis = false , 
                ylims = (0,:auto),
                xlabel = xlbl,
                label=false;
                HistArg...)  
        median_par = median(filter(isfinite, reshape(data_det[:,:,i],:)))
        if pltLine==true
            vline!([median_par],label = :none, color = :black, linestyle = :dash, linewidth = 2)
        end 
    if Mode == :n
        xlims!(0, round(median_par, digits = 0)*5)
        xl = xlims()
        x_txt = xl[2] - 0.05 * (xl[2] - xl[1])
        annotate!(x_txt, ylims()[2]*0.5, text(lbl, :right, fs-2, colors[i]))
    elseif Mode == :background 
        xmin = minimum(filter(isfinite,reshape(data_det[:,:,i],:)))-10
        xlims!(ifelse(xmin>0,xmin,0), maximum(filter(isfinite,reshape(data_det[:,:,i],:)))  )
        xl = xlims()
        x_txt = xl[2] - 0.05 * (xl[2] - xl[1])
        annotate!(x_txt, ylims()[2]*0.5, text(lbl, :right, fs-2, colors[i]))
    elseif Mode == :step_amplitude 
        xlims!(0, round(median_par, digits = 2) *5  )
        xl = xlims()
        x_txt = xl[2] - 0.05 * (xl[2] - xl[1])
        annotate!(x_txt, ylims()[2]*0.5, text(lbl, :right, fs-2, colors[i]))
    elseif Mode == :skew_frac
        xlims!(0, 0.26)
        xl = xlims()
        x_txt = xl[2] - 0.1 * (xl[2] - xl[1])
        annotate!(x_txt, ylims()[2]*0.5, text(lbl, :right, fs-2, colors[i]))
    elseif Mode == :skew_width
        xlims!(0, 1)
        xl = xlims()
        x_txt = xl[2] - 0.1 * (xl[2] - xl[1])
        annotate!(x_txt, ylims()[2]*0.5, text(lbl, :right, fs-2, colors[i]))
    else 
        xlims!(round(median_par, digits = 2) - 1, round(median_par, digits = 2) + 1)
        xl = xlims()
        x_txt = xl[1] + 0.05 * (xl[2] - xl[1])
        annotate!(x_txt, ylims()[2]*0.5, text(lbl, :left, fs-2, colors[i]))
    end
    end
    ptot = plot(p..., layout = (:, 1),
                    plot_title = "Calibration peak fits \npartition $(partition), $(e_type), $(det_mode) dets",
                    plot_titlefontsize = fs-2,
                size = (500,1200), left_margin = 10mm, right_margin = 5mm, top_margin = 0mm )

    path_plot = "$(@__DIR__)/plots/Ecal_FitPar_Overview/"
    fname = path_plot * "Ecal_Residuals_$(det_mode)_$(Mode)_part$(partition)_$(e_type).png"
    savefig(ptot,fname)
    @info "save plot to $fname"
end
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

