# script to test changes in LegendSpecFits package
# test new implementation of FWHM uncertainties and p-value/gof

using LegendSpecFits # make sure you're in development mode: dev Packages/LegendSpecFits.jl
using LegendDataManagement
using LegendHDF5IO, HDF5
using LegendDataTypes: fast_flatten, readdata
using Statistics, StatsBase
using Plots, ColorSchemes
using Distributions 
using Calculus, LinearAlgebra, Zygote
using ProgressMeter 
using TypedTables
using Printf
using BAT # for specfit_combined 
using InverseFunctions # for specfit_combined
using Optim 

include("utils.jl")

# get data path and filename
l200 = LegendData(:l200)
ch_geds     = channelinfo(l200,(:p03,:r000,:cal),system = :geds, only_processable = true).channel;
data        = load_hitch(l200,DataPeriod(3),DataRun(0); channel=ch_geds[1])
data_energy = data.e_cusp 

# ------------ do "simple" fit --> input for actual fit  ------------ #
result_simple, report_simple, th228_lines, th228_names = do_simple_calibration(data_energy); # simple calibration, get histograms around the calibration lines and peakstats

# cross-check; do fit iterative fits of gamma peaks 
result, report = fit_peaks(result_simple.peakhists, result_simple.peakstats, th228_lines,; uncertainty=true,iterative_fit=false);

#= combined fit 
# input: 
# - result_simple.peakhists: vector of histograms around the calibration lines
# - result_simple.peakstats: vector of statistics for the calibration line fits
=#
fit_peaks_combined(result_simple.peakhists, result_simple.peakstats, th228_lines;calib_type=:th228)

# copy elements of fit_peaks combined here to understand the function:
peakstats = result_simple.peakstats
peakhists = result_simple.peakhists
th228_combined_fit_functions = LegendSpecFits.th228_combined_fit_functions

pseudo_prior = LegendSpecFits.NamedTupleDist(
        e = BAT.ConstValueDist(th228_lines),
        # centroid parameter share
        μ_slope         = Uniform(0.99, 1.01),
        μ_intercept     = Uniform(-2.0, 2.0),
        # low-E tail parameter share
        skew_fraction_intercept = LogUniform(1e-4, 1e-1),
        skew_fraction_slope = Biweight(0, 1e-4),
        # single parameters
        σ              = product_distribution(weibull_from_mx.(peakstats.peak_sigma, 2*peakstats.peak_sigma)),
        n              = product_distribution(weibull_from_mx.(peakstats.peak_counts, 2*peakstats.peak_counts)),
        step_amplitude = product_distribution(weibull_from_mx.(peakstats.mean_background, 2*peakstats.mean_background)),
        skew_width     = product_distribution(fill(LogUniform(0.001, 0.1), length(peakstats))),
        background     = product_distribution(weibull_from_mx.(peakstats.mean_background, 2*peakstats.mean_background)),
    )
v_init = mean(pseudo_prior)
# transform back to frequency space
f_trafo = BAT.DistributionTransform(Normal, pseudo_prior)

f_loglike = let f_fit=th228_combined_fit_functions.f_fit, hists=peakhists
        v -> sum(hist_loglike.(Base.Fix2.(f_fit, expand_vars(v)), hists))
end

opt_r = optimize((-) ∘ f_loglike ∘ inverse(f_trafo), f_trafo(v_init), Optim.Options(iterations=100000))
    v_ml = inverse(f_trafo)(Optim.minimizer(opt_r))


    counts, bin_widths, bin_centers= LegendSpecFits._prepare_data(result_simple.peakhists[Idx])
model_fun = v -> LegendSpecFits._get_model_counts(th228_fit_functions.f_sigWithTail, v, bin_centers, bin_widths)