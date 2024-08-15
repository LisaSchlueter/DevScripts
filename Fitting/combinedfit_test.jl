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
det_types = reverse(unique(detector_type.(Ref(l200), dets_ged)))
fit_func = :f_fit#
det = dets_ged[22] 

result_simple, _, th228_names, _ = simple_cal(l200, runsel, det)
result, report = fit_peaks(result_simple.peakhists, result_simple.peakstats, th228_names; e_unit=result_simple.unit, calib_type=:th228, fit_func = fit_func, uncertainty = true);
p = plot(broadcast(k -> plot(report[k], left_margin=20mm, top_margin=-5mm, bottom_margin=-2mm, title=string(k), ms=2), th228_names)...,
         layout=(length(report), 1), size=(1000,710*length(report)) , thickness_scaling=1.8, titlefontsize = 10, legendfontsize = 8, yguidefontsize = 9, xguidefontsize=11)

### manual single-peak fit: 
h1 = result_simple.peakhists[5]
ps1 = result_simple.peakstats[5]
pseudo_prior = get_pseudo_prior(h1, ps1, fit_func)
# standard single-peak fit for reference. 
f_trafo = BAT.DistributionTransform(Normal, pseudo_prior) # transform back to frequency space

fit_function = get_th228_fit_functions(; background_center = ps1.peak_pos)[fit_func] # get fit function with background center
f_loglike = let f_fit = fit_function, h = h1 # create loglikehood function: f_loglike(v) that can be evaluated for any set of v (fit parameter)
    v -> hist_loglike(Base.Fix2(f_fit, v), h)
end

v_init = Vector(mean(f_trafo.target_dist))# start values for MLE
opt_r = optimize((-) ∘ f_loglike ∘ inverse(f_trafo), v_init, Optim.Options(time_limit = 60, iterations = 5000)) # MLE
v_ml = inverse(f_trafo)(Optim.minimizer(opt_r))# best fit results

pseudo_prior_medians = [quantile(pseudo_prior[p], 0.55) for p in keys(pseudo_prior)] # start values in parameter space
v_init2 = f_trafo(NamedTuple{keys(v_ml)}(pseudo_prior_medians)) # transform pseudo prior medians into frequency space
opt_r2 = optimize((-) ∘ f_loglike ∘ inverse(f_trafo), v_init2, Optim.Options(time_limit = 60, iterations = 5000)) # MLE
v_ml2 = inverse(f_trafo)(Optim.minimizer(opt_r2))
pval, chi2, dof = p_value_poissonll(fit_function, h1, v_ml2) # based on likelihood ratio 
for par in keys(v_ml)
    @printf("Δ%s: %.3f\n", String(par), mvalue(v_ml[par]) - mvalue(v_ml2[par]))
end
@printf("Δchi2 / pval: %.3f", v_ml.gof.pvalue - pval)

### NEW manual combined fit: 
h = result_simple.peakhists[2:3]
ps = result_simple.peakstats[2:3]
fixed_position = false
low_e_tail = true
pseudo_prior_def = collect(get_pseudo_prior.(h, ps, fit_func))
pars_singlepeak = keys(pseudo_prior_def[1])
tuple_keys = Tuple([Symbol(pars_singlepeak[i], Symbol(peak)) for peak in eachindex(ps) for i in eachindex(pars_singlepeak)])
tuple_vals = Tuple(Base.Iterators.flatten(values.(pseudo_prior_def)))
pseudo_prior_combi = NamedTupleDist(NamedTuple{tuple_keys}(tuple_vals))
f_trafo = BAT.DistributionTransform(Normal, pseudo_prior_combi)
v_init = Vector(mean(f_trafo.target_dist))
fit_functions = [get_th228_fit_functions(; background_center = peak_pos)[fit_func] for peak_pos in ps.peak_pos]

f_loglike = let f_fit=fit_functions, h=h
    v -> begin
        v_peaks =  [(µ = v[Symbol(:µ, Symbol(i))],
            σ = v[Symbol(:σ, Symbol(i))],
            n = v[Symbol(:n, Symbol(i))],
            step_amplitude = v[Symbol(:step_amplitude, Symbol(i))],
            skew_fraction = v[Symbol(:skew_fraction, Symbol(i))],
            skew_width = v[Symbol(:skew_width, Symbol(i))],
            background = v[Symbol(:background, Symbol(i))])  for i in eachindex(ps)]

        sum(hist_loglike.(Base.Fix2.(f_fit, v_peaks), h))
    end
end

# MLE
opt_r = optimize((-) ∘ f_loglike ∘ inverse(f_trafo), v_init, LBFGS(), Optim.Options(time_limit = 60, iterations = 10000))
v_ml = inverse(f_trafo)(Optim.minimizer(opt_r))

# compare
ustrip(mvalue(result[th228_names[1]].µ)) -  v_ml[Symbol(:μ, 1)]
ustrip(mvalue(result[th228_names[2]].µ)) -  v_ml[Symbol(:μ, 2)]
ustrip(mvalue(result[th228_names[1]].σ)) -  v_ml[Symbol(:σ, 1)]
ustrip(mvalue(result[th228_names[2]].σ)) -  v_ml[Symbol(:σ, 2)]
result[th228_names[1]].gof.converged


# # try to pass gradient to optimizer. analytical would be much better. 
# loglikegrad(v) =  ForwardDiff.gradient((-) ∘ f_loglike ∘ inverse(f_trafo), v)

# # f_loglike_array = let v_keys = keys(pseudo_prior_combi)
# #     v ->  -f_loglike(NamedTuple{v_keys}(v))
# # end
# # Calculate the Hessian matrix using ForwardDiff
# # Grad(v) =  ForwardDiff.hessian(f_loglike_array, tuple_to_array(v))
# # Grad(v_ml)