# script to test changes in LegendSpecFits package
# test new implementation of FWHM uncertainties and p-value/gof

# plan: look at all channels in this run. do all fits and compare fit parameters and goodness of fit. 
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
using Interpolations
include("utils.jl")

mode = "data"

if mode=="data"
    # get data path and filename
    l200 = LegendData(:l200)
    ch_geds     = channelinfo(l200,(:p03,:r000,:cal),system = :geds, only_processable = true).channel;
    data        = load_hitch(l200,DataPeriod(3),DataRun(0); channel=ch_geds[3]);
    data_energy = data.e_cusp; 
elseif mode=="mc"
    # get mc path and filename
    include("../../Packages/LegendSpecFits.jl/test/test_utils.jl")
    data_energy, _ = generate_mc_spectrum(200000)
end

# ------------ do "simple" fit --> input for actual fit  ------------ #
result_simple, report_simple, th228_lines_all, th228_names = do_simple_calibration(data_energy); # simple calibration, get histograms around the calibration lines and peakstats
plot(result_simple.peakhists)
vline!(result_simple.peakstats.peak_pos,color=:red,linestyle=:dash,linewidth=1.5,label = "simple calibration")
# ---------------------------------------------------------------------------------------
# Start global fit dev
# ---------------------------------------------------------------------------------------

# --------- 1. combined fit with no shared fit parameters 
th228_lines_sel = [1,2,3,5,7] #remove escape peaks: SEP and  DEP
th228_lines_fit = th228_lines_all[th228_lines_sel]
peakstats = result_simple.peakstats[th228_lines_sel]
peakhists = result_simple.peakhists[th228_lines_sel]
fit_function  = :f_fit_uncorr
result_combi,report_combi = LegendSpecFits.fit_peaks_combined_th228(peakhists,peakstats,th228_lines_fit,;fixed_position=false, fit_fun=fit_function, prior_flag=2)

#compare fith for-loop fit routine 
result,report = LegendSpecFits.fit_peaks(peakhists,peakstats,th228_lines_fit,;uncertainty=true)
for par in fieldnames(result_combi)
    println("Combined fit \tdiff \t|diff|/err")
    for i = 1:length(th228_lines_fit)
        println("$(par)=$(round(result_combi[par][i],digits=2)) \t$(round(result_combi[par][i]-result[th228_lines_fit[i]][par],digits=4)) \t $(abs(round((result_combi[par][i]-result[th228_lines_fit[i]][par])/result[th228_lines_fit[i]].err[par],digits=5)))")
    end
    println("------------------------------------------------------------")
end
println("pval single peak fits:")
for i = 1:length(th228_lines_fit)
    println("p = $(round(result[th228_lines_fit[i]].pval,digits=2))")
end
# println("µ=$(round(result_combi.µ[i],digits=2)) \t $(round(result[th228_lines_fit[i]].µ,digits=2)) \t bias = $(result_combi.µ[i]-result[th228_lines_fit[i]].µ) \t |bias|/err = $(abs(round((result_combi.µ[i]-result[th228_lines_fit[i]].µ)/result[th228_lines_fit[i]].err.µ,digits=5)))")
  
# -------------- investigate why there is a difference: do fit locally step-by-step
# # --> start fit routine 
pseudo_prior = LegendSpecFits.NamedTupleDist(
            µ              = product_distribution(Uniform.(peakstats.peak_pos .-10, peakstats.peak_pos .+10)),
            σ              = product_distribution(weibull_from_mx.(peakstats.peak_sigma, 2*peakstats.peak_sigma)),
            n              = product_distribution(weibull_from_mx.(peakstats.peak_counts, 2*peakstats.peak_counts)),
            step_amplitude = product_distribution(weibull_from_mx.(peakstats.mean_background, 2*peakstats.mean_background)),
            skew_fraction  = product_distribution(fill(Uniform(0.005, 0.25), length(peakstats))),
            skew_width     = product_distribution(fill(LogUniform(0.001, 0.1), length(peakstats))),
            background     = product_distribution(weibull_from_mx.(peakstats.mean_background, 2*peakstats.mean_background)),
        )

 # transform back to frequency space
 f_trafo = BAT.DistributionTransform(Normal, pseudo_prior)

 # start values for MLE
 v_init = mean(pseudo_prior)

 # create loglikehood function
 fit_fun = fit_function #account for different naming in LegendSpecFit
 f_loglike = let f_fit=LegendSpecFits.th228_combined_fit_functions[fit_fun], hists=peakhists
     v -> sum(hist_loglike.(Base.Fix2.(f_fit, expand_vars(v)), hists))
 end

 # MLE
 opt_r = optimize((-) ∘ f_loglike ∘ inverse(f_trafo), f_trafo(v_init), Optim.Options(iterations=100000))

 # best fit results
 v_ml = inverse(f_trafo)(Optim.minimizer(opt_r))

# # --> end fit routine 

# compare likelihoods 
 loglikesum = f_loglike(v_ml) #loglike of combined fit : all peaks  
 f_loglike2 = let f_fit=LegendSpecFits.th228_combined_fit_functions[fit_fun], hists=peakhists
    v -> (hist_loglike.(Base.Fix2.(f_fit, expand_vars(v)), hists))
end
loglike = -f_loglike2(v_ml) #loglike of combined fit : individual peaks  

f_loglike_single = zeros(length(th228_lines_fit))
for i= 1:length(peakhists)
    f_loglike_fit_peaks = let f_fit=LegendSpecFits.th228_fit_functions[:f_fit], h=peakhists[i]
        v -> -hist_loglike(Base.Fix2(f_fit, v), h)
    end
    v_ml_single = result[th228_lines_fit[i]]
    f_loglike_single[i] = f_loglike_fit_peaks(v_ml_single )

    println("Combi fit: $(round(loglike[i],digits=2))\t single peak fit: $(round(f_loglike_single[i],digits=2))")
end
#f_loglike_single

