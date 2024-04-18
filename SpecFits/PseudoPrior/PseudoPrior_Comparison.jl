# investiagete pseudo priors peak fits of calibration data 
# test pseudo prior on skew widths
using LegendDataManagement
using LegendDataManagement.LDMUtils
using LegendHDF5IO, HDF5
using LegendSpecFits
using Distributions
using Plots, ColorSchemes
using TypedTables
using Measures
using Printf

l200 = LegendData(:l200)
ch_geds     = channelinfo(l200,(:p03,:r000,:cal),system = :geds, only_processable = true).channel;
ch_dets     = channelinfo(l200,(:p03,:r000,:cal),system = :geds, only_processable = true).detector;
nChannel = length(ch_geds)
chIdx = 2
period = 3
run = 0 
include("../utils.jl")
data_det = load_hitch(l200, DataPeriod(period), DataRun(run); mode=:cal,channel=ch_geds[chIdx])
result_simple, report_simple, th228_lines, th228_names = do_simple_calibration([data_det.e_cusp...]); # simple calibration, get histograms around the calibration lines and peakstats


for pidx in eachindex(th228_names)
    ps = result_simple.peakstats[pidx]
    result_reg, report_reg = fit_single_peak_th228(result_simple.peakhists[pidx],ps)
    # pseudo_prior = LegendSpecFits.NamedTupleDist(skew_width = weibull_from_mx(0.001, 1e-2),)

    # fit with very wide pseudo-priors; effectively no prior
    wide_pseudo_prior = LegendSpecFits.NamedTupleDist(
        μ = Uniform(ps.peak_pos-1000, ps.peak_pos+1000),
        σ = weibull_from_mx(ps.peak_sigma, 100*ps.peak_sigma),
        n = weibull_from_mx(ps.peak_counts, 100*ps.peak_counts),
        step_amplitude = weibull_from_mx(ps.mean_background_step, ps.mean_background_step + 100*ps.mean_background_std),
        skew_fraction = truncated(weibull_from_mx(0.01, 1), 0.0, 0.25),
        skew_width = weibull_from_mx(0.001, 1),
        background = weibull_from_mx(ps.mean_background, ps.mean_background + 100*ps.mean_background_std),
    )
    result_wide, report_wide = fit_single_peak_th228(result_simple.peakhists[pidx],ps; pseudo_prior = wide_pseudo_prior)

    pcm_reg = plot(report_reg,:cormat, title = "Standard pseudo-prior \n Correlation Matrix")
    pcm_wide = plot(report_wide,:cormat, title = "Loose pseudo-prior \n Correlation Matrix")

    p_reg = plot(report_reg)
    p_wide = plot(report_wide)
    pall = plot(pcm_reg,pcm_wide,p_reg,p_wide,layout=(2,2),size=(950,850), plot_title = "period-$period, run-$run, e_cusp, $(ch_dets[chIdx]), $(th228_names[pidx])", plot_titlefontsize = 14, dpi = 300)

    plt_path = "$(@__DIR__)/plots/"
    plt_name = plt_path * "PeakFit_PseudoPrior_p$(period)_r$(run)_$(ch_dets[chIdx])_ecusp_$(th228_names[pidx]).png"
    savefig(pall, plt_name)
    @info "save plot to $plt_name"
end