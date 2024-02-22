# script to test changes in LegendSpecFits package
# test new implementation of FWHM uncertainties and p-value/gof
using LegendSpecFits # make sure you're in development mode: dev Packages/LegendSpecFits.jl
using LegendDataManagement
using LegendHDF5IO, HDF5
using LegendDataTypes: fast_flatten, readdata
using Statistics, StatsBase
using Plots, ColorSchemes
using Distributions # new
using Calculus, LinearAlgebra, Zygote # new
using ProgressMeter # new
## ------------ get input for fit: energie histogram  ------------ ##
# get file and channel name
l200 = LegendData(:l200)
period     = search_disk(DataPeriod,l200.tier[:jlhitch,:cal])[1]     # use first period for which files exists
run        = search_disk(DataRun,l200.tier[:jlhitch,:cal,period])[1] # use first run for which files exists
path       = l200.tier[:jlhitch,:cal,period,run] .* "/" # get data path
filenames  = path .* readdir(path);
pattern =  r"(?<=-ch)(.*)(?=-tier_jlhit\.lh5)";  # regular expression, look for channel name in file
channels = [match(pattern, filename).match for filename in filenames]; # extract channel number from filename

# open file and get energy spectrum
data = h5open(x -> readdata(x, "ch$(channels[1])/dataQC"), filenames[1])
data_energy = data.e_cusp

# ------------ do "simple" fit --> input for actual fit  ------------ #
result_simple, report_simple, th228_lines, th228_names = do_simple_calibration(data_energy); # simple calibration, get histograms around the calibration lines and peakstats

## ------------ do actual fit    ------------ ##
Idx =5 # choose peak 5 as example 
data_energy_hist = result_simple.peakhists[Idx] # 1 histogram for each th peak
result_simple.peakstats[Idx] # fit results of simple calibration
result_fit, report_fit = fit_single_peak_th228(data_energy_hist, result_simple.peakstats[Idx] ; uncertainty=true);
fit_par_names= propertynames(result_fit)
par_name_str = collect(map(x -> string(x),fit_par_names[1:7]))
fit_par   = result_fit[fit_par_names[1:7]]
fit_err = result_fit[:err]; fit_err = fit_err[fit_par_names[1:7]]

## ------------ look at fit results. plot data and model   ------------ ##
# look at the separate steps: fit of a histogram = binned MLW fit.
data_hist_bin_edges = first(data_energy_hist.edges)
data_hist_counts = data_energy_hist.weights
data_hist_bin_centers = (data_hist_bin_edges[begin:end-1] .+ data_hist_bin_edges[begin+1:end]) ./ 2
data_hist_bin_widths = data_hist_bin_edges[begin+1:end] .- data_hist_bin_edges[begin:end-1]

# get model function
ModelFun_fitpar  = Base.Fix2(th228_fit_functions.f_fit, fit_par) # fix the fit parameters to ML best-estimate
ModelFun = data_hist_bin_widths.*map(energy->ModelFun_fitpar(energy), data_hist_bin_centers) # evaluate model at bin center (= binned measured energies)  

# plot data and model
plot(data_energy_hist,facecolor=:"dodgerblue",linewidth=0.0,fillalpha=0.8,label="Data")
plot!(data_hist_bin_centers,ModelFun,linewidth=3,color="orange",label="Best fit")
xlabel!("Energy (keV)", fontsize = 16); ylabel!("Counts"); 

## ------------ understand loglikelihood calculation, some details   ------------ ##
fit_fun = :f_fit
h = data_energy_hist
f_loglike = let f_fit=th228_fit_functions[fit_fun], h=data_energy_hist
        v -> hist_loglike(Base.Fix2(f_fit, v), h)
end

function f_loglike_fun(v) # this is equivalent
        hist_loglike(Base.Fix2(th228_fit_functions.f_fit, v), h)
end
f_loglike_fun(fit_par) == f_loglike(fit_par)  #equivalent 
   
# --------- calculate hessian matrix with f_log_array ------------ #
# to get hessian, one needs the loglikelihood function in a form that is acceps an array and returns a scalar 
# new implementation that fit function doesn't have to be hard coded 
f_loglike_array = let f_fit=gamma_peakshape, h=h # old version
        v ->  - hist_loglike(    x -> f_fit(x, v...), h     ) 
end

f_loglike_array_new = let f_fit=th228_fit_functions[fit_fun], h=h, v_keys = keys(fit_par)  # new version, variable fit function 
    # input v is an array. convert back to touple in order to use general fit_fun 
    fit_par_tup =  NamedTuple{keys(fit_par)}(fit_par_vec)
    v ->  - hist_loglike(    x -> f_fit(x,NamedTuple{v_keys}(v)), h) 
end

# hessian needs array. 
H     = ForwardDiff.hessian(f_loglike_array, [fit_par...])
H_new = ForwardDiff.hessian(f_loglike_array_new, [fit_par...])
H==H_new

loglike_bf = -hist_loglike(ModelFun_fitpar,h) # loglikelihood for best-fit

result, report = fit_single_peak_th228(data_energy_hist, result_simple.peakstats[Idx] ; uncertainty=true);
result.pval
plot(report)

Idx_fit =2
result, report = fit_single_peak_th228(result_simple.peakhists[Idx_fit], result_simple.peakstats[Idx_fit] ; uncertainty=true,fit_fun=fit_fun);
plot(report)
fit_fun2=:f_sigWithBkg
result_2, report_2 = fit_single_peak_th228(result_simple.peakhists[Idx_fit], result_simple.peakstats[Idx_fit] ; uncertainty=true, fit_fun=fit_fun2);


# -----------------------------------------------------------------
# --------------- develop iterative fits ------------------------
# reminder:
#fit_peaks(peakhists::Array, peakstats::StructArray, th228_lines::Array,; calib_type::Symbol=:th228, uncertainty::Bool=true, low_e_tail::Bool=true)
#th228_lines =  [583.191,  727.330,  860.564,  1592.53,    1620.50,    2103.53,    2614.51]
th228_sym = [:583.191,  :727.330,  :860.564,  :1592.53,    :1620.50,    :2103.53,    :2614.51]
#th228_names =  ["Tl208a", "Bi212a", "Tl208b", "Tl208DEP", "Bi212FEP", "Tl208SEP", "Tl208FEP"]
result, report = fit_peaks(result_simple.peakhists, result_simple.peakstats, th228_lines,; uncertainty=true,iterative_fit=true);

# look at fit and data 
Idx_sel = 7
fit_par_names= propertynames(result[:2103.53])[1:7]
par_name_str = collect(map(x -> string(x),fit_par_names))
fit_par   = result[th228_sym[Idx_sel]]
counts, bin_widths, bin_centers = LegendSpecFits._prepare_data(result_simple.peakhists[Idx_sel]);
model_counts =LegendSpecFits._get_model_counts(th228_fit_functions.f_fit, fit_par, bin_centers,bin_widths)
bar(bin_centers,counts,label="data")
plot!(bin_centers,model_counts,label="model, p = $(round(result[th228_sym[Idx_sel]].pval,digits=2))",linewidth=3,color="orange")
xlims!(result[th228_sym[Idx_sel]].µ-8*result[th228_sym[Idx_sel]].σ,result[th228_sym[Idx_sel]].µ+8*result[th228_sym[Idx_sel]].σ)
xlabel!("Energy (eV)"); ylabel!("Counts"); 

#f_sig
fit_fun = :f_fit
th228_fit_functions[fit_fun] == th228_fit_functions.f_fit
