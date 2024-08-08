"""
1. plot background model including linear slope with fit result
2. plot corresponding pseudo_prior for background slope for each peak 
"""

using LaTeXStrings, Printf
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
include("./utils_plot.jl")
colors = _get_def_peakcolors()
saveplot = true 

# data selection 
l200 = LegendData(:l200)
runsel = (DataPeriod(3), DataRun(0), :cal)
chinfo = channelinfo(l200, runsel; system=:geds, only_processable=true)
det = chinfo.detector[22]
fit_func = :f_fit_bckExp

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

# fit peaks with background slope  
peaks = collect(1:7)
result_fit, report_fit = fit_peaks(result_simple.peakhists[peaks], result_simple.peakstats[peaks], th228_names[peaks]; e_unit=result_simple.unit, calib_type=:th228, fit_func = fit_func);

function strip_units(named_tuple)
    fnames = fieldnames(named_tuple)
    vals = map(mvalue,map(ustrip, values(named_tuple)[1:end-1]))
    return NamedTuple{fnames[1:end-1]}(vals)
end

# get background model 
function background_model(x, result_fit; background_exp = 0.0, background_center = 0.0 )
    v = strip_units(result_fit)
    if background_exp == 0 
        return background_peakshape(x, v.μ, v.σ, v.step_amplitude, v.background)
    else
        return background_peakshape(x, v.μ, v.σ, v.step_amplitude, v.background;  background_exp =  background_exp, background_center = background_center)
    end
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
for (i, peak) in enumerate(th228_names)
    local ps = result_simple.peakstats[i]
    local h = result_simple.peakhists[i]

    local _, _, bin_centers = LegendSpecFits._prepare_data(h)
    local background_center = ps.peak_pos 
    # plot
    local x =  collect(range(bin_centers[1], bin_centers[end], length = 1000))
    local xl = filter(x-> x < ps.peak_pos, x)
    local xr = filter(x-> x > ps.peak_pos, x)
    
    background_def(E) = background_model.(E, Ref(result_fit[peak])) 
    background_expfit(E) = background_model.(E, Ref(result_fit[peak]); background_exp = mvalue(result_fit[peak].background_exp), background_center = ps.peak_pos)
    background_exp(E, bexp) = background_model.(E, Ref(result_fit[peak]); background_exp = bexp, background_center = ps.peak_pos)
    
    # get pseudo prior for background slope
    local bexp = [0, exp10.(range(-7, -1.2, length = 5000))...]
    local prior_cdf = cdf.(LegendSpecFits.get_pseudo_prior(h, ps, fit_func).background_exp, bexp)
    stdcl =  diff(cdf.(Normal(0,1),[-1,1]))[1]
    bexp_std = round(bexp[findfirst(prior_cdf .> stdcl) ], digits = 3)

    local ptmp = plot(xl, background_def(xl), fillrange = background_exp(xl, bexp_std), 
            alpha = 0.5, fillcolor = :dodgerblue, color = :white, label =  @sprintf("Background decay constant pseudoprior %.2f%%" , 100*stdcl), linewidth = 0.1)
    plot!(xr,background_def(xr), fillrange = background_exp(xr, bexp_std),
            alpha = 0.5, fillcolor = :dodgerblue, color = :white, label = false, linewidth = 0.1)        
    plot!(x, background_def(x) , label = "Background model: τ = 0.0 ", color = :black, lienwidth = 2)
    plot!(x, background_expfit(x), label = @sprintf("Background model: τ = %.1g (%s fit) ",mvalue(result_fit[peak].background_exp), peak), color = :darkorange, linewidth = 3, linestyle = :dash)
    xlims!(minimum(x), maximum(x))
    push!(plts, ptmp)
end

plot(plts..., layout = (length(plts), 1), size = (700, 450*length(plts)),
 legend = :topright, left_margin = 25mm, right_margin = 3mm, 
 plot_title = "$(runsel[1]), $(runsel[2]), $det \n" * L"Background = step + B $\cdot \exp(-τ \cdot (E-E_\textrm{center})$)", titlefontsize = 12)
ylabel!("Counts")
xlabel!("Energy")
if saveplot 
    plt_path = "$(@__DIR__)/plots/"
    if !isdir(plt_path)
        mkpath(plt_path)
    end
    plt_name = plt_path * "background_exp_model_$(runsel[1])_$(runsel[2])_$det.png"
    savefig(plt_name)
end


# # plot also pseudo prior distributions 
# for (i, peak) in enumerate(th228_names)
i = 5; peak = th228_names[i]
ps = result_simple.peakstats[i]
h = result_simple.peakhists[i]
window_left = ps.peak_pos - minimum(h.edges[1])
window_right = maximum(h.edges[1]) - ps.peak_pos
# get pseudo prior for background slope
bexp = exp10.(range(-7, -1.2, length = 5000)) #range(0.0, 0.5, length = 1000)#
prior_cdf = cdf.(LegendSpecFits.get_pseudo_prior(h, ps, fit_func).background_exp, bexp)
bexp_median = round(bexp[findfirst(prior_cdf .> 0.5)], digits = 3)
stdcl =  diff(cdf.(Normal(0,1),[-1,1]))[1]
bexp_std = round(bexp[findfirst(prior_cdf .> stdcl) ], digits = 3)
prior_pdf = pdf.(LegendSpecFits.get_pseudo_prior(h, ps, fit_func).background_exp, bexp)

p1 = plot(bexp, prior_cdf, linewidth = 2, 
        label = false, color = colors[peak],
        legend = :bottomright,
        xlabel = "Background decay constant τ (1 / keV)", 
        ylabel = "Cum. probability (c.d.f.)",
        yscale = :log10,
        ylims = (1e-2 , 1),
        title =  "Pseudo-prior distribution for exponential background \n Weibull distribution with τ(50%) = $bexp_median,  τ($((round(100*stdcl, digits = 2)))%) = $bexp_std \n",
        grid = :xy,
        y_foreground_color_axis = colors[peak],
        y_foreground_color_text = colors[peak],
        y_foreground_color_border = colors[peak],
        y_guidefontcolor = colors[peak]
        )
p2 = twinx(p1)
plot!(p2, bexp, 100 .* (exp.(bexp*window_left) - exp.(-bexp*window_right)), 
     linewidth = 2, 
     linestyle = :dashdot, 
     color = :midnightblue,
     ylabel = "ΔB/B fit window (%)", 
     label = false,
     ylim = (1e-2, 1000),#(9e-4, 1000),
     yscale = :log10,
     y_foreground_color_axis = :midnightblue,
     y_foreground_color_text = :midnightblue,
     y_foreground_color_border = :midnightblue,
     y_guidefontcolor = :midnightblue
     )
hline!(p1,[0.50], color = :silver, label = false, linestyle = :solid )
xlims!(-1e-3, 0.051)

if saveplot 
    plt_path = "$(@__DIR__)/plots/"
    if !isdir(plt_path)
        mkpath(plt_path)
    end
    plt_name_slope = plt_path * "background_exp_pseudoprior.png"
    savefig(plt_name_slope)
end

