using LegendDataManagement
using LegendDataManagement.LDMUtils
using LegendSpecFits
import LegendSpecFits: get_th228_fit_functions, hist_loglike, get_pseudo_prior
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
using BAT,  Optim, InverseFunctions
include("../utils/utils_ecal.jl")
include("../utils/utils_plot.jl")

# data selection 
l200 = LegendData(:l200)
runsel = (DataPeriod(3), DataRun(0), :cal)
chinfo = channelinfo(l200, runsel; system=:geds, only_processable=true)
dets_ged  = chinfo.detector
det_types = reverse(unique(detector_type.(Ref(l200), dets_ged)))
fit_func = :f_fit#Slope

det = dets_ged[22] 
result_simple, _, th228_names, _ = simple_cal(l200, runsel, det)
result, report = fit_peaks(result_simple.peakhists, result_simple.peakstats, th228_names; e_unit=result_simple.unit, calib_type=:th228, fit_func = peakshape, uncertainty = true);
p = plot(broadcast(k -> plot(report[k], left_margin=20mm, top_margin=-5mm, bottom_margin=-2mm, title=string(k), ms=2), th228_names)..., layout=(length(report), 1), size=(1000,710*length(report)) , thickness_scaling=1.8, titlefontsize = 10, legendfontsize = 8, yguidefontsize = 9, xguidefontsize=11)

### NEW: 
h = result_simple.peakhists[5]
ps = result_simple.peakstats[5]
pseudo_prior = get_pseudo_prior(h, ps, fit_func)
# standard single-peak fit for reference. 
f_trafo = BAT.DistributionTransform(Normal, pseudo_prior)# transform back to frequency space
v_init = Vector(mean(f_trafo.target_dist))# start values for MLE
fit_function = get_th228_fit_functions(; background_center = ps.peak_pos)[fit_func] # get fit function with background center
f_loglike = let f_fit = fit_function, h = h # create loglikehood function: f_loglike(v) that can be evaluated for any set of v (fit parameter)
    v -> hist_loglike(Base.Fix2(f_fit, v), h)
end

opt_r = optimize((-) ∘ f_loglike ∘ inverse(f_trafo), v_init, Optim.Options(time_limit = 60, iterations = 5000)) # MLE
Optim.converged(opt_r)
v_ml = inverse(f_trafo)(Optim.minimizer(opt_r))# best fit results

opt_r = optimize((-) ∘ f_loglike ∘ inverse(f_trafo), v_init, Optim.Options(time_limit = 60, iterations = 500)) # MLE
Optim.converged(opt_r)
v_ml500 = inverse(f_trafo)(Optim.minimizer(opt_r))# best fit results

for k in keys(v_ml)
    @printf("Δ%s: %.2g \n", String(k), mvalue(v_ml[k]) - mvalue(v_ml500[k]))
end
