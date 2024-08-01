"""
1. plot background model including linear slope with fit result
2. plot corresponding pseudo_prior for background slope for each peak 
"""
using Plots, LaTeXStrings
using Distributions
using LegendSpecFits
using LegendDataManagement
using LegendDataManagement.LDMUtils
using LegendHDF5IO
using TypedTables
using Unitful, Measures
using Measurements: value as mvalue, uncertainty as muncert
using StructArrays, PropDicts
using Revise
using Plots
using Distributions 
using ValueShapes

saveplot = true 

# data selection 
l200 = LegendData(:l200)
runsel = (DataPeriod(3), DataRun(0), :cal)
chinfo = channelinfo(l200, runsel; system=:geds, only_processable=true)
det = chinfo.detector[10]

# get energy spectrum in ADC
"""
    get_energy_ctc(data::LegendData, runsel::Tuple{DataPeriod, DataRun, Symbol}, detector::DetectorId; e_type::Symbol=:e_cusp)
# open hitchannel file and apply charge-trapping correction
"""
function get_energy_ctc(data::LegendData, runsel::Tuple{DataPeriod, DataRun, Symbol}, detector::DetectorId; e_type::Symbol=:e_cusp)
    channel = detector2channel(l200, runsel, detector)
    (period, run, category) = runsel
    data_ch = lh5open(data.tier[:jlhitch, category, period, run, channel])[channel, :dataQC][:]#.e_cusp
    ecusp_ctc_uncal = ljl_propfunc(data.par.rpars.ctc[period, run, detector][e_type].func).(data_ch)
    return ecusp_ctc_uncal
end
e_cusp_ctc_ADC = get_energy_ctc(l200, runsel, det)

# get config for calibration fits
energy_config = dataprod_config(l200).energy(start_filekey(l200, runsel))
energy_config_ch = merge(energy_config.default, get(energy_config, det, PropDict()))

quantile_perc = if energy_config_ch.quantile_perc isa String parse(Float64, energy_config_ch.quantile_perc) else energy_config_ch.quantile_perc end
th228_names = Symbol.(energy_config_ch.th228_names)
th228_lines_dict = Dict(th228_names .=> energy_config_ch.th228_lines)

# simple calibration 
result_simple, report_simple = simple_calibration(e_cusp_ctc_ADC, energy_config_ch.th228_lines, energy_config_ch.left_window_sizes, energy_config_ch.right_window_sizes,; calib_type=:th228, n_bins=energy_config_ch.n_bins, quantile_perc=quantile_perc, binning_peak_window=energy_config_ch.binning_peak_window)

# fit 1 peak with background slope  
peaks = collect(1:7)

result_fit_slope, report_fit_slope = fit_peaks(result_simple.peakhists[peaks], result_simple.peakstats[peaks], th228_names[peaks]; e_unit=result_simple.unit, calib_type=:th228, 
                                            fit_func = :f_fit_WithBkgSlope);#, pseudo_prior = pseudo_prior_add)


# get background model 
function background_model(x, background, background_step, background_slope, peak_pos)
    if background_slope == 0 
        return background_peakshape.(x, peak_pos,  2.0, background_step, background)
    else
        return background_peakshape.(x, peak_pos,  2.0, background_step, background, background_slope)
    end
end

function max_background_slope(h, ps)
    window_left = ps.peak_pos - minimum(h.edges[1])
    window_right = maximum(h.edges[1]) - ps.peak_pos
    pseudo_prior_std = ps.mean_background_std  / (window_left + window_right)
    pseudo_prior = LegendSpecFits.get_standard_pseudo_prior(h, ps, :f_fit_WithBkgSlope)
    x_bkg = range(-pseudo_prior_std, pseudo_prior_std, length = 1000)
    bkg_slope_min = x_bkg[findfirst(x -> x>0 , pdf.(pseudo_prior.background_slope, x_bkg))]
    bkg_slope_min = ifelse(bkg_slope_min > -pseudo_prior_std, bkg_slope_min, -pseudo_prior_std)
    bkg_slope_max = x_bkg[findlast(x -> x>0 , pdf.(pseudo_prior.background_slope, x_bkg))]
    bkg_slope_max = ifelse(bkg_slope_max < pseudo_prior_std, bkg_slope_max, pseudo_prior_std)    
    return [bkg_slope_min, bkg_slope_max]
end

#result_simple.peakhists[peak]
default(grid = false, framestyle = :semi,
xguidefontsize = 14, xtickfontsize = 10,
xlabelfontsize = 16, ylabelfontsize = 16,
yguidefontsize = 14, ytickfontsize = 10,
legendfontsize = 10, legendtitlefontsize = 10,
titlefontsize = 11,
foreground_color_legend = :silver,
background_color_legend = :white,
tick_direction = :in,
)

plts = []
for peak in eachindex(th228_names)
    local ps = result_simple.peakstats[peak]
    local h = result_simple.peakhists[peak]
    bkg_slope_fit = mvalue(result_fit_slope[th228_names[peak]].background_slope)
    bkg_slope_ps = max_background_slope(h, ps)

    # plot
    x_bkg =  collect(range(minimum(h.edges[1]),maximum(h.edges[1]), length = 1000))
    x_bkg_l = filter(x -> x< ps.peak_pos, x_bkg)
    x_bkg_r = filter(x -> x> ps.peak_pos, x_bkg)
    ptmp = plot(x_bkg_l, background_model(x_bkg_l, ps.mean_background, ps.mean_background_step,  bkg_slope_ps[2], ps.peak_pos), fillrange = background_model(x_bkg_l, ps.mean_background, ps.mean_background_step,  bkg_slope_ps[1], ps.peak_pos), 
            alpha = 0.5, fillcolor = :dodgerblue, color = :white, label = "Background-slope pseudoprior (1σ)", linewidth = 0.1)
    plot!(x_bkg_r, background_model(x_bkg_r, ps.mean_background, ps.mean_background_step,  bkg_slope_ps[1], ps.peak_pos), fillrange = background_model(x_bkg_r, ps.mean_background, ps.mean_background_step,  bkg_slope_ps[2], ps.peak_pos), 
            alpha = 0.5, fillcolor = :dodgerblue, color = :white, label = false, linewidth = 0.1)
    plot!(x_bkg,  background_model(x_bkg, ps.mean_background, ps.mean_background_step, 0, ps.peak_pos), label = "Background without slope", color = :black, lienwidth = 2)
    plot!(x_bkg, background_model(x_bkg, ps.mean_background, ps.mean_background_step, bkg_slope_fit, ps.peak_pos), label = "Background with fit slope", color = :darkorange, linewidth = 3, linestyle = :dash)
    xlims!(minimum(x_bkg), maximum(x_bkg))
    push!(plts, ptmp)
end

plot(plts..., layout = (length(plts), 1), size = (700, 450*length(plts)), legend = :topright, left_margin = 25mm, right_margin = 3mm, plot_title = "$(runsel[1]), $(runsel[2]), $det", titlefontsize = 12)
ylabel!("Counts")
xlabel!("Energy")
if saveplot 
    plt_path = "$(@__DIR__)/plots/"
    if !isdir(plt_path)
        mkpath(plt_path)
    end
    plt_name = plt_path * "background_slope_model_$(runsel[1])_$(runsel[2])_$det.png"
    savefig(plt_name)
end


# plot also pseudo prior distributions 
plts_slope = []
for peak in eachindex(th228_names)
   local ps = result_simple.peakstats[peak]
   local h = result_simple.peakhists[peak]
   window_left = ps.peak_pos - minimum(h.edges[1])
   window_right = maximum(h.edges[1]) - ps.peak_pos
   pseudo_prior_std = ps.mean_background_std  / (window_left + window_right)

    dist = LegendSpecFits.get_standard_pseudo_prior(h, ps, :f_fit_WithBkgSlope).background_slope
    x = range(-3*pseudo_prior_std, 3*pseudo_prior_std, length = 1000)
    ptmp = plot(x, pdf.(dist, x)/maximum(pdf.(dist, x)), label = "$(th228_names[peak]) \nbackground = $(round(ps.mean_background, digits = 0))", linewidth = 3, legend = :topleft)
    push!(plts_slope, ptmp)
end
plot(plts_slope..., layout = (length(plts_slope), 1), size = (700, 450*length(plts_slope)), legend = :bottom, left_margin = 25mm, right_margin = 3mm, plot_title = "$(runsel[1]), $(runsel[2]), $det", titlefontsize = 12)
ylabel!("Probability density (a.u.)")
xlabel!(L"Background slope (keV$^{-1}$)")
ylims!(0,1.1)

if saveplot 
    plt_path = "$(@__DIR__)/plots/"
    if !isdir(plt_path)
        mkpath(plt_path)
    end
    plt_name_slope = plt_path * "background_slope_pseudoprior_$(runsel[1])_$(runsel[2])_$det.png"
    savefig(plt_name_slope)
end