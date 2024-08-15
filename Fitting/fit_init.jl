# investigate initivalizationin gamma peak fits 
using LegendDataManagement
using LegendDataManagement.LDMUtils
using LegendSpecFits
import LegendSpecFits: get_th228_fit_functions, hist_loglike, get_pseudo_prior, tuple_to_array
using LegendHDF5IO
using TypedTables
using Unitful, Measures
using Measurements: value as mvalue, uncertainty as muncert
using StructArrays, PropDicts
using Revise
using Plots
using Distributions 
using ValueShapes
using LaTeXStrings, Printf
using LegendDataTypes: fast_flatten
using BAT,  Optim, InverseFunctions, ForwardDiff
include("../utils/utils_ecal.jl")
include("../utils/utils_plot.jl")

# data selection 
l200 = LegendData(:l200)
runsel = (DataPeriod(3), DataRun(0), :cal)
chinfo = channelinfo(l200, runsel; system=:geds, only_processable=true)
dets_ged  = chinfo.detector
det = dets_ged[22] 
fit_func = :f_fit#

# simple calibration and automatic fit
result_simple, _, th228_names, th228_lines_dict = simple_cal(l200, runsel, det)
result, report = fit_peaks(result_simple.peakhists, result_simple.peakstats, th228_names; e_unit=result_simple.unit, calib_type=:th228, fit_func = fit_func, uncertainty = true);
p = plot(broadcast(k -> plot(report[k], left_margin=20mm, top_margin=-5mm, bottom_margin=-2mm, title=string(k), ms=2), th228_names)...,
         layout=(length(report), 1), size=(1000,710*length(report)) , thickness_scaling=1.8, titlefontsize = 10, legendfontsize = 8, yguidefontsize = 9, xguidefontsize=11)

### single-peak fit with different intial parameter 
quantiles_i = 0.05:0.05:0.95
pval = ones(length(quantiles_i), 7)
µ = ones(length(quantiles_i), 7)
σ = ones(length(quantiles_i), 7)
skew_frac = ones(length(quantiles_i), 7)
for peak in 1:7
    h1 = result_simple.peakhists[peak]
    ps1 = result_simple.peakstats[peak]
    pseudo_prior = get_pseudo_prior(h1, ps1, fit_func)
    f_trafo = BAT.DistributionTransform(Normal, pseudo_prior) # transform back to frequency space
    fit_function = get_th228_fit_functions(; background_center = ps1.peak_pos)[fit_func] # get fit function with background center
    f_loglike = let f_fit = fit_function, h = h1 # create loglikehood function: f_loglike(v) that can be evaluated for any set of v (fit parameter)
        v -> hist_loglike(Base.Fix2(f_fit, v), h)
    end

    # define initial guess for fit parameter and do minimization
    function get_v_init(q::Real = 0.50)
        pseudo_prior_medians = [quantile(pseudo_prior[par], q) for par in keys(pseudo_prior)] # start values in parameter space
        return f_trafo(NamedTuple{keys(pseudo_prior)}(pseudo_prior_medians)) # transform pseudo prior medians into frequency space
    end

    # # default v_init 
    # v_init_default = Vector(mean(f_trafo.target_dist))# start values for MLE
    # if all(v_init_default .== get_v_init(0.5)) @info "v_init default is the same as median of pseudo-prior after distribution transform" end
    # opt_r_def = optimize((-) ∘ f_loglike ∘ inverse(f_trafo), v_init_default, Optim.Options(time_limit = 60, iterations = 5000)) # MLE
    # v_ml_def = inverse(f_trafo)(Optim.minimizer(opt_r_def))# best fit results
    # pval, chi2, dof = p_value_poissonll(fit_function, h1, v_ml_def) # based on likelihood ratio 
    # v_ml_def = merge(v_ml_def, (gof = (pvalue = pval, chi2 = chi2, dof = dof),) )

    # test quantiles systemtically 
    for (idx, q) in enumerate(quantiles_i)
        try 
        local opt_r = optimize((-) ∘ f_loglike ∘ inverse(f_trafo), get_v_init(q), Optim.Options(time_limit = 60, iterations = 5000)) # MLE
        local v_ml = inverse(f_trafo)(Optim.minimizer(opt_r))# best fit results
        µ[idx, peak] = mvalue(v_ml.µ)
        σ[idx, peak] = mvalue(v_ml.σ)
        skew_frac[idx, peak] = mvalue(v_ml.skew_fraction)
        pval[idx, peak], _, _ = p_value_poissonll(fit_function, h1, v_ml) # based on likelihood ratio 
        catch 
            @info "Optimization failed for peak $(th228_names[peak]) quantile $q"
            µ[idx, peak] = NaN
            σ[idx, peak] = NaN
            pval[idx, peak] = NaN
            skew_frac[idx, peak] = NaN
        end
    end
end

lst = [:solid :dash :dot :dashdot :solid :dash :dashdot]
_def_plot()
begin
    qidx =  findfirst(x-> x == 0.5, quantiles_i)
    pµ = plot(quantiles_i,  µ .- (µ[qidx,:] .* ones(7,length(quantiles_i)))', ylabel = " Δ µ (keV)", legend = false, linewidth = 2.5, linestyle = lst)
    vline!([0.5], color = :silver, linewidth = 2, linestyle = :dash, label = false)
    pσ = plot(quantiles_i, σ, ylabel = "σ (keV)", legend=false, linewidth = 2.5, xticks = [], linestyle = lst,  ylims = (0, 5))
    vline!([0.5], color = :silver, linewidth = 2, linestyle = :dash)
    pskewfrac = plot(quantiles_i, skew_frac, xlabel = "Pseudo-prior quantile for v_init", ylabel = "low-energy tail fraction", label = hcat(String.(th228_names)...), legend= false, linewidth = 2.5, linestyle = lst, ylims = (0, 1))
    vline!([0.5], color = :silver, linewidth = 2, linestyle = :dash, label = false)
    ppval = plot(quantiles_i, pval, xlabel = "Pseudo-prior quantile for v_init", ylabel = "p-value", legend=false, linewidth = 2.5, linestyle = lst)
    vline!([0.5], color = :silver, linewidth = 2, linestyle = :dash)
    pleg = plot(ones(5, 7), ones(5, 7), legend = :topleft, linestyle = lst, label = hcat(String.(th228_names)...), framestyle = :none, legend_columns = 7 )
    plot(pµ, pσ, pskewfrac, ppval, pleg, layout=(3,2), size=(1500, 1200), bottom_margin = 0mm, plot_title = "$(runsel[1]), $(runsel[2]), $(det)", thickness_scaling = 1.5)
end
pltdir = "$(@__DIR__)/plots/"
pltname = pltdir * "fit_init_pseudopriorquantile_$(runsel[1])_$(runsel[2])_$(det).png"
savefig(pltname)


