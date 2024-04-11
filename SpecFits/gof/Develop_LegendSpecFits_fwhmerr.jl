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

# aux functions
include("utils.jl")

# load hitch data 
l200 = LegendData(:l200)
ch_geds     = channelinfo(l200,(:p03,:r000,:cal),system = :geds, only_processable = true).channel;
data        = load_hitch(l200,DataPeriod(3),DataRun(0); channel=ch_geds[1])
data_energy = data.e_cusp 

# ------------ do "simple" fit --> input for actual fit  ------------ #
result_simple, report_simple, th228_lines, th228_names = do_simple_calibration(data_energy); # simple calibration, get histograms around the calibration lines and peakstats

function estimate_fwhm_local(v::NamedTuple)
    # get FWHM
    try
        half_max_sig = maximum(Base.Fix2(th228_fit_functions.f_sigWithTail, v).(v.μ - v.σ:0.001:v.μ + v.σ))/2
        roots_low = find_zero(x -> Base.Fix2(th228_fit_functions.f_sigWithTail, v)(x) - half_max_sig, v.μ - v.σ, maxiter=100)
        roots_high = find_zero(x -> Base.Fix2(th228_fit_functions.f_sigWithTail,v)(x) - half_max_sig, v.μ + v.σ, maxiter=100)
        return roots_high - roots_low
    catch e
        return NaN
    end
end

# ------------ fit single peak  ------------ #
Idx = 5 # choose 1 peak as example
result_fit, report_fit = fit_single_peak_th228(result_simple.peakhists[Idx], result_simple.peakstats[Idx] ; uncertainty=true);
fwhm_fit = result_fit.fwhm
fwhm_err_corr = result_fit.err.fwhm#0.32 --> improvement by factor 2. without constrainted sampling: 0.48

# get uncorrelated values with nearest SPE : 
fit_par_names  = keys(result_fit)
v_ml_err = NamedTuple{fit_par_names[1:7]}(sqrt.(abs.(diag(result_fit.covmat))))
v_ml = result_fit[fit_par_names[1:7]]
_, fwhm_uncorr = get_peak_fwhm_th228(v_ml, v_ml_err)#0.53 

# get uncorrelated values without nearest SPE : 
v_ml_err = NamedTuple{fit_par_names[1:7]}(sqrt.(abs.(diag(result_fit.covmat_raw))))
_, fwhm_uncorr_raw = get_peak_fwhm_th228(v_ml, v_ml_err) # 0.665


#  --------------- fwhm uncertainty: step-by-step --------------------------- 

#generate MC (model) spectra from correlated fit paramters
v = v_ml
v_err = result_fit.covmat
n = 1000
function get_fwhm_err(v,v_err,n,sel)
    if !isposdef(v_err)
        v_err = nearestSPD(v_err)
        @debug "Covariance matrix not positive definite. Using nearestSPD"
    end
    v_err = v_err[1:6,1:6] #remove background, keep only relevant for sampling 
    v_fitpar = v[keys(v)[1:size(v_err,1)]] # only fit parameter
    dist = MvNormal([v_fitpar...], v_err) # multivariate distribution using covariance matrix)
    v_mc = rand(dist, n) # Draw samples

    # constain fit_par_samples to physical values. warning hardcoded. tbd 
    Idx_keep = findall((v_mc[3,:].>0) .*                #  positive amplitude 
                        (v_mc[5,:].<0.25).*             # skew fraction 
                        (v_mc[5,:].>0) .*    #skew fraction 
                        (v_mc[6,:].>0))                 # positive skew width
   
    if sel == true
        v_mc_sel = v_mc[:,Idx_keep];
    else 
        v_mc_sel = v_mc;
    end
    n_sel = size(v_mc_sel,2)
    v_mc_sel_tuple = [NamedTuple{keys(v)[1:size(v_err,1)]}(v_mc_sel[:,i]) for i=1:n_sel] # convert back to NamedTuple 

    fwhm_mc = estimate_fwhm.(v_mc_sel_tuple)
    frac_nan = sum(isnan.(fwhm_mc))/n_sel
    fwhm_err = std(fwhm_mc[isfinite.(fwhm_mc)])
    return fwhm_mc, fwhm_err, n_sel,v_mc, v_mc_sel,v_mc_sel_tuple,frac_nan
end

fwhm_mc, fwhm_err, n_sel,v_mc, v_mc_sel,v_mc_sel_tuple, frac_nan = get_fwhm_err(v_ml,result_fit.covmat,20000)

# histogram of  sampled parameters: all and those kept for fehm err estimation
h = [] # selected 
for i=1:size(v_mc_sel,1)
    bin_width = (maximum(v_mc[i,:])-minimum(v_mc[i,:]))/100
    bin_edges = minimum(v_mc[i,:]):bin_width:maximum(v_mc[i,:])
    tmp = histogram(v_mc[i,:],bins=bin_edges,alpha=0.7, linecolor =:silver,fillcolor=:silver,label = "All samples") # sanity check: plot distribution of samples
    histogram!(v_mc_sel[i,:],bins=bin_edges,alpha=0.5,linecolor =:red, fillcolor=:red,label = "Constrained samples" )# sanity check: plot distribution of samples
       push!(h,tmp)
    xlabel!("$(string(fit_par_names[i]))"); ylabel!("Occurence")
    if i==1
       # xticks!(round(v_ml.µ-2*v_ml.σ, digits=0):1:round(v_ml.µ+2*v_ml.σ, digits=0))
       µ = round(v_ml.µ, digits=1) 
       xticks!([µ-0.4,µ,µ+0.4])
    end
end
plot(h[1],h[2],h[3],h[5],h[6],layout=(3,2),legend=false,size = (600,600),grid=false)
n = size(v_mc,2)
n_eff = round(n_sel/n,digits=2)
plt_name = "./plots/fwhm_err_sampling_$(n)samples_$(n_eff)eff.pdf"
savefig(plt_name)# plot correlation matrix


# look at fwhm uncertainty as a function of sample 
n_test          = 100:100:20000
fwhm_err_sel   = zeros(length(n_test))
n_sel_test      = zeros(length(n_test))
fwhm_err_all   = zeros(length(n_test))
n_sel_all      = zeros(length(n_test))
for i in 1:length(n_test)
    _, fwhm_err_sel[i], n_sel_test[i], = get_fwhm_err(v_ml,result_fit.covmat,n_test[i],true)
    _, fwhm_err_all[i], n_sel_all[i], = get_fwhm_err(v_ml,result_fit.covmat,n_test[i],false)
end
formatter(x) = @sprintf("%.0f", x)
plot(n_test,fwhm_err_sel,size=(600,400),label = "constrained",legendtitle="Sampling",foreground_color_legend=:silver,legendfontsize=12,guidefontsize=14,tickfontsize=12,xformatter=formatter,linewidth=2)
plot!(n_test,fwhm_err_all,size=(600,400),label="unconstrained",guidefontsize=14,tickfontsize=12,xformatter=formatter,linewidth=2)
.-mean(fwhm_err_sel),.-mean(fwhm_err_all)
xlabel!("Number of samples"); 
ylabel!("FWHM uncertainty (eV)")
title!("period $(string(period)), run $(string(run)), detector $(string(det_good))",fontsize=14)
plt_name2 = "./plots/fwhm_err_sampling_ConvergenceTest_CorrelatedSampling.pdf"
savefig(plt_name2)# plot correlation matrix

plot(n_test,fwhm_err_sel.-mean(fwhm_err_sel),size=(600,400),label = "constrained",legendtitle="Sampling",foreground_color_legend=:silver,legendfontsize=12,guidefontsize=14,tickfontsize=12,xformatter=formatter,linewidth=2)
plot!(n_test,fwhm_err_all.-mean(fwhm_err_all),size=(600,400),label="unconstrained",guidefontsize=14,tickfontsize=12,xformatter=formatter,linewidth=2)
xlabel!("Number of samples"); 
ylabel!("FWHM uncertainty - mean (eV)")
title!("period $(string(period)), run $(string(run)), detector $(string(det_good))",fontsize=14)
plt_name3 = "./plots/fwhm_err_sampling_ConvergenceTest_CorrelatedSamplingDev.pdf"
savefig(plt_name3)# plot correlation matrix

#plot(n_test,fwhm_err_sel,size=(600,400),label = "constrained",legendtitle="Sampling",foreground_color_legend=:silver,legendfontsize=12,guidefontsize=14,tickfontsize=12,xformatter=formatter,linewidth=2)
plot(n_test,fwhm_err_all,size=(600,400),label="unconstrained",legendtitle="Sampling",foreground_color_legend=:silver,legendfontsize=12,guidefontsize=14,tickfontsize=12,xformatter=formatter,linewidth=2)
xlabel!("Number of samples"); 
ylabel!("FWHM uncertainty (eV)")
title!("period $(string(period)), run $(string(run)), detector $(string(det_good))",fontsize=14)
plt_name4 = "./plots/fwhm_err_sampling_ConvergenceTest_UnonstrainedCorrelatedSampling.pdf"
savefig(plt_name4)# plot correlation matrix


# investigation of stuff
nan_indices = findall(isnan, fwhm_samples)
good_indices = 1:n_samples#findall(!isnan, fwhm_samples) 
max_idx = findall(fwhm_samples[good_indices].==maximum(fwhm_samples[good_indices]))[1]
fwhm_err_mc = std(fwhm_samples[good_indices])

counts, bin_widths, bin_centers= LegendSpecFits._prepare_data(result_simple.peakhists[Idx])
model_fun = v -> LegendSpecFits._get_model_counts(th228_fit_functions.f_sigWithTail, v, bin_centers, bin_widths)
model_bf = model_fun(fit_par)

""" attempt via interpolation. doesn't work as is
counts, bin_widths, bin_centers = _prepare_data(data_energy_hist); # prepare data
fwhm_samples_inter = [get_fwhm(bin_centers,model_fun(named_tuples[i]),named_tuples[i]) for i=1:n_samples]
@info "In $(100*sum(isnan.(fwhm_samples_inter))/n_samples)% of the samples, fwhm calculation via interpolation fails "
fwhm_err_inter_mc = std(fwhm_samples_inter[findall(!isnan, fwhm_samples_inter)])
"""

# look why fwhm estimations fails sometimes. plot example
plot(bin_centers,model_bf,linewidth=3,color="orange",label="Model best-fit, FWHM = $(round(fwhm_fit,digits=2)) eV",legend=:topright,guidefontsize = 16)
#plot!(bin_centers,model_fun(named_tuples[nan_indices[1]]),linestyle=:dot,linewidth=3,marker=".",color="blue",label="Sample model NaN FWHM")
plot!(bin_centers,model_fun(named_tuples[good_indices[1]]),linestyle=:dash ,linewidth=3,marker=".",color="red",label="Sample model FWHM = $(round(fwhm_samples[good_indices[1]],digits=2)) eV")
plot!(bin_centers,model_fun(named_tuples[max_idx]),linestyle=:dash ,linewidth=3,marker=".",color="green",label="Sample model max. FWHM = $(round(fwhm_samples[max_idx],digits=2)) eV")
xlabel!("Energy (keV)", fontsize = 16); ylabel!("Counts"); 
xlims!(result_fit.µ-5*fwhm_fit,result_fit.µ+5*fwhm_fit)


# calculate std of non-nan fwhm samples
fwhm_samples_good = fwhm_samples[good_indices]
histogram(fwhm_samples_good,bins=100, label="fwhm")
fwhm_samples_err = std(fwhm_samples_good)

