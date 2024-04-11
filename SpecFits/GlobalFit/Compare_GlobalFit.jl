# compare global fit with single peak fits
# no shared fit parameters
# look at likelihoods, compare minima 

# plan: look at all channels in this run. do all fits and compare fit parameters and goodness of fit. 
using LegendSpecFits # make sure you're in development mode: dev Packages/LegendSpecFits.jl
using LegendDataManagement
using LegendHDF5IO, HDF5
using LegendDataTypes: fast_flatten, readdata
using Statistics, StatsBase
using Plots, ColorSchemes
using Distributions
using ProgressBars
using Printf
include("../utils.jl")

l200 = LegendData(:l200)
ch_geds     = channelinfo(l200,(:p03,:r000,:cal),system = :geds, only_processable = true).channel;
nChannel    =  length(ch_geds)

iter = ProgressBar(1:nChannel) # 35, 41,84 dont work
path = "PackageDevScripts/DevelopSpecFits/GlobalFit/results"
g_tols = [1e-8,1e-9,1e-10,1e-11,1e-12]

fit_fail = []

for i in iter  # loop over channels 
    data        = load_hitch(l200,DataPeriod(3),DataRun(0); channel=ch_geds[i]);
    data_energy = data.e_cusp; 
    # ------------ do "simple" fit --> input for actual fit  ------------ #
    result_simple, report_simple, th228_lines, th228_names = do_simple_calibration(data_energy); # simple calibration, get histograms around the calibration lines and peakstats
    set_description(iter, string(@sprintf("Fit %s (%.0f out of %.0f)", string(ch_geds[i]),i,nChannel)))
  
    for g_tol in g_tols # looper of g_tol (convergence criteria for gradient-based minimizer)
        try
            result,report = LegendSpecFits.fit_peaks_th228(result_simple.peakhists,result_simple.peakstats,th228_lines,;g_tol = g_tol) # do to: put this
            # prepare fit results for saving
            par_names_tuple = fieldnames(result[th228_lines[1]]) # tuple of symbols 
            par_names =  [string.(par_names_tuple )...][1:12] # vector of strings
            par_tuple = map(x -> [x...][1:end-1],[values(result)...])#vector of touples. each element is result from 1 peak fit 
            par_all = [map(x-> x[i],par_tuple) for i=1:length(par_names)]
            err_tuple =  values(map(x-> x[end],[values(result)...])) #vector of touples. each element is result from 1 peak fit 
            err_names = "err_" .* [string.(fieldnames(err_tuple[1]))...]
            err_all  = [map(x-> x[i],err_tuple) for i=1:length(err_names)]
           fname = "single_peak_fit_p3_run0_$(string(ch_geds[i]))_gtol$(g_tol).h5"
            save_result(fname,[par_names...,err_names...,"th228_lines"],[par_all...,err_all...,th228_lines],;path=path)
    
        catch
            push!(fit_fail, [ch_geds[i],g_tol])
            @info "fit failed"
        end

        #=
        # global fit
        prior = 2
        fit_func = :f_fit_uncorr
        result_combi, _, gof = LegendSpecFits.fit_peaks_combined_th228(result_simple.peakhists,result_simple.peakstats,th228_lines,;fixed_position=false, fit_fun=fit_func, prior_flag=prior)
        par_names_combi =  [string.(fieldnames(result_combi) )...]
        par_all_combi = [map(x-> x[i],expand_vars(result_combi)) for i=1:length(par_names_combi)]
        gof_names_combi =  [string.(fieldnames(gof) )...]
        gof_combi  = [map(x-> x[i],expand_vars(gof)) for i=1:length(gof_names_combi)] 
        fname_combi = "combi_peak_fit_p3_run0_$(string(ch_geds[i]))_prior$(prior)_$(string(fit_func)).h5"
        save_result(fname_combi,[par_names_combi...,gof_names_combi...,"th228_lines"],[par_all_combi...,gof_combi...,th228_lines],;path=path)
        =#
    end
end
