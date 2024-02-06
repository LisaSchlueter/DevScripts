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

# aux functions
include("utils.jl")

# load hitch data 
l200 = LegendData(:l200)
ch_geds     = channelinfo(l200,(:p03,:r000,:cal),system = :geds, only_processable = true).channel;
data        = load_hitch(l200,DataPeriod(3),DataRun(0); channel=ch_geds[1])
data_energy = data.e_cusp 

## ------------ do "simple" fit --> input needed for actual fit  ------------ ##
result_simple, report_simple, th228_lines, th228_names = do_simple_calibration(data_energy); # simple calibration, get histograms around the calibration lines and peakstats

# look whats inside fit results "result_simple"
Idx =5
data_energy_hist = result_simple.peakhists[Idx] # data (energies) in form of a histogram
result_simple.peakstats[Idx] # fit results of simple calibration

## ------------ do actual fit  --> developemtn happends HERE  ------------ ##
#= comments on fit routine
gamma_peakshape(x::Real, μ::Real, σ::Real, n::Real,step_amplitude::Real, skew_fraction::Real, skew_width::Real,background::Real) # model function. x = energy vector, for which which model is evaluated
th228_fit_functions.f_fit = (x, v) -> gamma_peakshape(x, v.μ, v.σ, v.n, v.step_amplitude, v.skew_fraction, v.skew_width, v.background) # why wrap around?
=#
# do fit and look at fit results
#fit_single_peak_th228(h::Histogram, ps::NamedTuple{(:peak_pos, :peak_fwhm, :peak_sigma, :peak_counts, :mean_background), NTuple{5, T}}; uncertainty::Bool=true, fixed_position::Bool=false) where T<:Real
result_fit, report_fit = fit_single_peak_th228(data_energy_hist, result_simple.peakstats[Idx] ; uncertainty=true);
fit_par_names= propertynames(result_fit)
par_name_str = collect(map(x -> string(x),fit_par_names[1:7]))
fit_par   = result_fit[fit_par_names[1:7]]
fit_err = result_fit[:err]; fit_err = fit_err[fit_par_names[1:7]]


## ------------ look at residuals    ------------ ##

counts, bin_widths, bin_centers = LegendSpecFits._prepare_data(data_energy_hist)
f_fit =  th228_fit_functions.f_fit 
v_ml = fit_par
# get peakshape of best-fit 
model_counts =LegendSpecFits._get_model_counts(f_fit, v_ml, bin_centers,bin_widths)

# calculate bin-wise residuals 
residuals    = model_counts[model_counts.>0]-counts[model_counts.>0]
sigma        = sqrt.(model_counts[model_counts.>0])
residuals_norm = residuals./sigma

# calculate something like a bin-wise p-value (in case that makes sense)
dist = Poisson.(model_counts) # each bin: poisson distributed 
cdf_value_low = cdf.(dist, model_counts.-abs.(residuals)) 
cdf_value_up  = 1 .-cdf.(dist, model_counts.+abs.(residuals))  
p_value_binwise = cdf_value_low .+ cdf_value_up # ~probability that residuals are as observed or higher


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

## ------------ p-value  ------------ ##
# ------------ option 1 (chi-squared test) ------------ 
# should be fine as long counts not too low (should be >5 in each bin)
chi2    = sum((ModelFun[ModelFun.>0]-data_hist_counts[ModelFun.>0]).^2 ./ ModelFun[ModelFun.>0])
npar    = length(fit_par)
ndof    = length(data_hist_counts[ModelFun.>0])-npar
chi2red = chi2/ndof
pval    = ccdf(Chisq(ndof),chi2)
if any(ModelFun.<=5)
    @warn "WARNING: bin with <=$(minimum(ModelFun)) counts -  chi2 test might be not valid"
else  
    @info "p-value = $(round(pval,digits=2))"
end
chi2_LogLike   = 2*sum(ModelFun.*log.(ModelFun./data_hist_counts)+ModelFun-data_hist_counts)
pval_logLike   = ccdf(Chisq(ndof),chi2_LogLike)
#------------ p-value option 2 (likelihood MC)  ------------ 
# compare max. likelikehood (best-fit) with likelihood distribution
# problem: need to now max. likelihood pdf -->  computational expensive
#  MC for likelihood: 1) randomize energy histogram: using best-fit model and poisson uncertainties
#                     2) fit randomized histogram and get max. likelihood   
#                     3) calculate max. likelihood pdf from samples                      

# draw random samples from (correlated) multivariate distribution
#likelihood function:   hist_loglike(f_fit::Base.Callable, h::Histogram{<:Real,1})  
loglike_bf = -hist_loglike(ModelFun_fitpar,data_energy_hist) # loglikelihood for best-fit

# randomize  data: counts in histogram 
n_samples = 1000
dists = Poisson.(ModelFun) # create poisson distributions for each bin
data_samples_vec = rand.(dists,n_samples) # randomize counts in histogram
data_samples = [ [] for _ in 1:n_samples ]
for i = 1:n_samples
    data_samples[i] = map(x -> x[i],data_samples_vec)
end
histogram(data_samples[rand(1:n_samples)],label = "samlple"); histogram!(data_energy_hist,label = "data")# sanity check

loglike_bf_sample = NaN.*ones(n_samples)
hist_sample = data_energy_hist
@showprogress for i=1:n_samples
    hist_sample.weights = data_samples[i]
    result_fit_sample, report = fit_single_peak_th228(hist_sample, result_simple.peakstats[Idx] ; uncertainty=true) 
    fit_par_sample   = result_fit_sample[fit_par_names[1:7]]
    ModelFun_fitpar_sample  = Base.Fix2(th228_fit_functions.f_fit, fit_par_sample) # fix the fit parameters to ML best-estimate
    loglike_bf_sample[i] = -hist_loglike(ModelFun_fitpar_sample,hist_sample) # loglikelihood for best-fit
end

histogram(loglike_bf_sample)

plot(vline!([loglike_bf],color="red",label="Best fit log-likelihood",linewidth=3))
pval_samples = sum(loglike_bf_sample.<=loglike_bf)./n_samples

#cdf_approx = cumsum(loglike_bf_sample./sum(loglike_bf_sample))
# get model function

#=
ModelFun_fitpar_samples  = [Base.Fix2(th228_fit_functions.f_fit, named_tuples[i]) for i=1:n_samples] # fix the fit parameters to ML best-estimate
loglike_samples = NaN.*ones(n_samples)
KeepLog = ones(n_samples)
for i=1:n_samples
try loglike_samples[i] = hist_loglike(ModelFun_fitpar_samples[i],data_energy_hist) 
catch
    KeepLog[i] = 0
    @warn "WARNING: loglikelihood calculation failed for sample $i"
end
end
loglike_samples = loglike_samples[KeepLog.==1]
histogram(-loglike_samples,bins=100, label="MC log-likelihoods")
vline!([-loglike_bf],color="red",label="Best fit log-likelihood",linewidth=3)
xlabel!("- Loglikelihood")
#loglike_samples = [try hist_loglike(ModelFun_fitpar_samples[i],data_energy_hist) catch end  for i=1:n_samples] # loglikelihood for all samples
sum(-loglike_samples.<-loglike_bf) # of course, this is the best fit! OMG
