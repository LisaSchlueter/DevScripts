import PhysicalConstants.CODATA2018: N_A, m_u # avogadro number, atomic mass constant
using LegendDataManagement
using Unitful
using Plots
using Distributions 
using Optim, ForwardDiff
using Interpolations, LinearAlgebra
using Measurements
using Printf
using JLD2

include("0nbb_likelihood.jl")
include("0nbb_simspec.jl")
include("0nbb_fit.jl")

nSamples = 100000
reCalc = false
v_input = (exposure = 200, efficiency = 0.95, fit_window = 240, Qbb = 2039.0, Qbb_sigma = 2.5, sig_exp = 1.0e-26) 
v_mc = (signal_strength = 0, background = 5.2e-4)
signal_strength_scan = collect(0:0.02:1)

loglikelihood_scan = zeros(nSamples, length(signal_strength_scan))
# generate mc spectrum 

# prepare for saving results 
fpath = "$(@__DIR__)/results/"
fname = fpath * @sprintf("mc0nbbfit_profile_samples%i_sig%.2f_bck%.1e.jld2", nSamples, v_mc.signal_strength, v_mc.background)
if isfile(fname) && reCalc == false 
    @info "file $fname exists already - skip compuation"
else
    @info "do fits for $fname"
    results = Vector{Any}(undef, nSamples)
   
    for i in 1:nSamples
        local events_mc = simspec(v_input, v_mc...; plotFlag = false)
        local result, report = fit_0nbb_mc(v_input, events_mc; uncertainty = false)
        for (j, sig) in enumerate(signal_strength_scan)
            local result_scan, _ = fit_0nbb_mc(v_input, events_mc; fixsignal = sig, uncertainty = false)
            loglikelihood_scan[i, j] = result_scan.gof.loglikelihood_ml
        end
        results[i] = result
    end

    # print info on convergence 
    if any([results[i].converged for i in 1:nSamples] .== 0)
        println("Some fits did not converge.")
    else
        println("All fits converged.")
    end

    # save 
    if !isdir(fpath)
        mkdir(fpath)
    end

    jldsave(fname, results = results, loglikelihood_scan = loglikelihood_scan, v_input = v_input, v_mc = v_mc, signal_strength_scan = signal_strength_scan)
    @info "fits done - saved to  $fname"
end

# loglikelihood_bf = [results[i].gof.loglikelihood_ml for i in 1:nSamples]
# teststatistic = 2 .* (loglikelihood_scan .- loglikelihood_bf )
# plot(signal_strength_scan, teststatistic, linewidth = 2)

