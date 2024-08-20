
using LegendDataManagement
using LegendDataManagement.LDMUtils
using LegendSpecFits
using LegendHDF5IO
using PropDicts
using Plots
using StatsPlots  
using Unitful, Measurements, Measures
using StatsBase, Optim
# include("../utils/utils_plot.jl") 

# channel 1, 
l200 = LegendData(:l200)
runsel = (DataPeriod(3), DataRun(0), :cal)
chinfo = channelinfo(l200, runsel; system=:geds, only_processable=true)

# load data
detIdx = 5
(period, run, category) = runsel
data_ch = read_ldata(l200, :jldsp, category, period, run, chinfo.channel[detIdx])

# get configs 
energy_config = dataprod_config(l200).energy(start_filekey(l200, runsel))
energy_config_ch = merge(energy_config.default, get(energy_config, chinfo.detector[detIdx], PropDict()))
ctc_config_ch = merge(energy_config.ctc.default, get(energy_config.ctc, chinfo.detector[detIdx], PropDict()))
quantile_perc = if energy_config_ch.quantile_perc isa String parse(Float64, energy_config_ch.quantile_perc) else energy_config_ch.quantile_perc end

# remove NaNs. Where are they coming from?
IdxKeepE = findall(isfinite, data_ch.e_cusp)
energy_ADC = data_ch.e_cusp[IdxKeepE]
qdrift = data_ch.qdrift[IdxKeepE]
IdxKeepQ = findall(isfinite, qdrift)
qdrift = qdrift[IdxKeepQ]
energy_ADC = energy_ADC[IdxKeepQ]

# simple calibration 
result_simple, report_simple = simple_calibration(energy_ADC, energy_config_ch.th228_lines, energy_config_ch.left_window_sizes, energy_config_ch.right_window_sizes,; calib_type=:th228, n_bins=energy_config_ch.n_bins, quantile_perc=quantile_perc, binning_peak_window=energy_config_ch.binning_peak_window)
plot(report_simple)

# charge trapping correction: does not work because of demain error; negative expectation value in poisson log likelihood
result_ctc, report_ctc = ctc_energy(energy_ADC .* result_simple.c, qdrift, ctc_config_ch.peak, (ctc_config_ch.left_window_size, ctc_config_ch.right_window_size), result_simple.c; e_expression="e_cusp")
result_ctc.converged
result_ctc.fct
plot(report_ctc)
# plot(report_ctc, true)

#######################################
# attempt to make 2dhistogram plot less ugly
# _def_plot()
 begin
    pspec = stephist(report_ctc.e_peak, fill = true, color = :darkgrey, label = "Before correction", legend = :topleft)
    stephist!(report_ctc.e_ctc, fill = true, color = :purple, alpha = 0.5, label = "After correction")
    ylims!(0, maximum(ylims()))
    plot!(xlabel = "Energy (keV)")
    xlims!(2590, 2641)
    xl = xlims()
    
    dens_before = BAT.kde((ustrip(report_ctc.e_peak), report_ctc.qdrift_peak ./ maximum(report_ctc.qdrift_peak)))
    dens_ctc = BAT.kde((ustrip(report_ctc.e_ctc), report_ctc.qdrift_peak ./ maximum(report_ctc.qdrift_peak)))
    pdens = plot(dens_before, c = reverse(cgrad(:greys)), fill = true, label = "before", colorbar = :none)
    plot!(dens_ctc, c = :plasma, alpha = 0.9, label = "after", framestyle = :semi, fill = false)
    # ylims!(0, quantile(report_ctc.qdrift_peak, 0.999))
    plot!(yformatter = :plain, legend = :topright, xlabel = "Energy (keV)", ylabel = "Drift time (a.u.)")
    annotate!(xl[1]+1, 0.16*quantile(report_ctc.qdrift_peak, 0.999), text("Before correction", :black, :left))
    annotate!(xl[1]+1, 0.08*quantile(report_ctc.qdrift_peak, 0.999), text("After correction", :darkviolet, :left))
    # annotate!(xl[1]+1, 0.25*quantile(report_ctc.qdrift_peak, 0.999), text("Tl208 FEP:", :black, :left))
    xlims!(xl)

    plot(pspec, pdens, layout = (2, 1), size = (800, 1000), thickness_scaling = 1.1, 
        plot_title = "Charge-trapping correction - Tl208 FEP peak ", plot_titlefontsize = 14,
        left_margin = 3mm)

    # mkpath("$(@__DIR__)/plots")
    # savefig("$(@__DIR__)/plots/ctc_kde_$(period)_$(run)_$(chinfo.detector[detIdx]).png")
 end



##########################################################
# charge trapping step by step
##########################################################
# match variable names with names in ctc.jl
e = energy_ADC.*result_simple.c 
window = (ctc_config_ch.left_window_size, ctc_config_ch.right_window_size)
peak = ctc_config_ch.peak

# create cut window around peak
cut = peak - first(window) .< e .< peak + last(window)
e_cut, qdrift_cut = e[cut], qdrift[cut]
e_unit = u"keV"
# calculate optimal bin width
bin_width_window = 5.0u"keV"
bin_width        = LegendSpecFits.get_friedman_diaconis_bin_width(e[peak - bin_width_window .< e .< peak + bin_width_window])
bin_width_qdrift = LegendSpecFits.get_friedman_diaconis_bin_width(qdrift[peak - bin_width_window .< e .< peak + bin_width_window])

# get FWHM before correction
# fit peak
h_before = fit(Histogram, ustrip.(e_unit, e_cut), ustrip(e_unit, minimum(e_cut)):ustrip(e_unit, bin_width):ustrip(e_unit, maximum(e_cut)))
ps_before = estimate_single_peak_stats(h_before)
result_before, report_before = fit_single_peak_th228(h_before, ps_before; uncertainty=true, low_e_tail = true)

fwhm = LegendSpecFits.estimate_fwhm(result_before)
v = result_before
vml = NamedTuple{Tuple(keys(v)[1:7])}(Measurements.value.(values(v))[1:7])

v = vml
if v.skew_fraction < 0.5
    e_low = v.μ  - v.σ 
    e_high = v.μ  + v.σ 
else
    e_low = vml.μ * (1 - vml.skew_width) 
    e_high = vml.μ * (1 + vml.skew_width)
end 

f_sigWithTail = Base.Fix2(LegendSpecFits.get_th228_fit_functions().f_sigWithTail, vml)
half_max_sig = maximum(f_sigWithTail.(e_low:0.001:e_high))/2
roots_low = LegendSpecFits.find_zero(x -> f_sigWithTail(x) - half_max_sig, e_low, maxiter=100)
roots_high = LegendSpecFits.find_zero(x -> f_sigWithTail(x) - half_max_sig,  e_high, maxiter=100)
roots_high - roots_low

plot(report_before, legend = :topleft)
result_before.background
result_before.step_amplitude
result_before.skew_width
result_before.skew_fraction
result_before.fwhm
result_before.gof.chi2

# get value for this channel 
qdrift_median = median(qdrift_cut)
fct_range = [ustrip(e_unit, 1e-4u"keV") / qdrift_median, ustrip(e_unit, 10u"keV") / qdrift_median]
fct_start = [ustrip(e_unit, 0.1u"keV") / qdrift_median]
fct = minimum(fct_range)
fct_fit = l200.par.rpars.ctc[period, run, chinfo.detector[detIdx]].e_cusp.fct

# debug LegendSpecFits.f_optimize_ctc()
e=ustrip.(e_unit, e_cut) 
bin_width=ustrip(e_unit, bin_width)
e_ctc = e .+ fct .* qdrift_cut
any(e_ctc .< 0)
# maximum negative before e negative
ftcmax = - e./qdrift_cut # if ftcmax becomes this value or smaller, then e becomes negative
maximum(ftcmax)

# fit peak
h = fit(Histogram, e_ctc, minimum(e_ctc):bin_width:maximum(e_ctc))
ps = estimate_single_peak_stats(h)
result_peak, report_peak = fit_single_peak_th228(h, ps; uncertainty=false)






# @recipe function f(report_ctc::NamedTuple{(:peak, :window, :fct, :bin_width, :bin_width_qdrift, :e_peak, :e_ctc, :qdrift_peak, :h_before, :h_after, :fwhm_before, :fwhm_after, :report_before, :report_after)})
#     layout := (2, 1)
#     size := (600, 700)
#     framestyle := :semi
#     grid := false 
#     left_margin --> (5, :mm)
#     right_margin --> (5, :mm)
#     foreground_color_legend := :silver
#     background_color_legend := :white
#     xtickfontsize := 12
#     xlabelfontsize := 14
#     ylabelfontsize := 14
#     ytickfontsize := 12
#     legendfontsize := 12
#     xl = (2590, 2641)
#     hist_bins = range(xl[1], xl[2], 100)
#     @series begin
#         seriestype := :stephist 
#         fill := true
#         color := :darkgrey
#         label := "Before correction"
#         legend := :topleft
#         subplot := 1
#         bins := hist_bins
#         report_ctc.e_peak
#     end
#     @series begin
#         seriestype := :stephist 
#         fill := true
#         alpha := 0.5
#         color := :purple
#         label := "After correction"
#         legend := :topleft
#         subplot := 1
#         ylims := (0, :auto)
#         xlims := xl
#         xlabel := "Energy"
#         bins := hist_bins
#         report_ctc.e_ctc
#     end
#     @series begin
#         seriestype := :line
#         subplot := 2
#         c := :binary
#         colorbar := :none
#         fill := true
#         label := "Before correction"
#         BAT.kde((ustrip(report_ctc.e_peak), report_ctc.qdrift_peak ./ maximum(report_ctc.qdrift_peak)))
#     end
#     @series begin
#         seriestype := :line
#         subplot := 2
#         c := :plasma
#         colorbar := :none
#         fill := false
#         label := "After correction"
#         xlims := xl
#         ylims := (0, 1)
#         yformatter := :plain
#         xlabel := "Energy (keV)"
#         ylabel := "Drift time (a.u.)"
#         BAT.kde((ustrip(report_ctc.e_ctc), report_ctc.qdrift_peak ./ maximum(report_ctc.qdrift_peak)))
#     end
# end