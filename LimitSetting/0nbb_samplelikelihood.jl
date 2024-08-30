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

nSamples = 100000 # 200000
reCalc = true
signal_strength  = collect(0.0:0.1:5.0) #1 ./ collect(0.1:0.2:4.1)[1] # in units of 1e-26 yr^-1
v_input = (exposure = 200, efficiency = 0.95, fit_window = 240, Qbb = 2039.0, Qbb_sigma = 2.5, sig_exp = 1.0e-26) 

for sig in signal_strength
    local v_mc = (signal_strength = sig, background = 5.2e-4)
   
    # prepare for saving results 
    fpath = "$(@__DIR__)/results/"
    fname = fpath * @sprintf("mc0nbbfit_samples%i_sig%.2f_bck%.1e.jld2", nSamples, v_mc.signal_strength, v_mc.background)
    if isfile(fname) && reCalc == false 
        @info "file $fname exists already - skip compuation"
    else
        @info "do fits for $fname"
        local results = Vector{Any}(undef, nSamples)
        local results_true = Vector{Any}(undef, nSamples)

        for i in 1:nSamples
            local events_mc = simspec(v_input, v_mc...; plotFlag = false)
            local result, report = fit_0nbb_mc(v_input, events_mc; uncertainty = false)
            local result_true, _ = fit_0nbb_mc(v_input, events_mc; fixsignal = v_mc.signal_strength, uncertainty = false)
            results[i] = result
            results_true[i] = result_true
        end

        # print info on convergence 
        if any([results[i].converged for i in 1:nSamples] .== 0) || any([results_true[i].converged for i in 1:nSamples] .== 0)
            println("Some fits did not converge.")
        else
            println("All fits converged.")
        end

        # save 
        if !isdir(fpath)
            mkdir(fpath)
        end
        jldsave(fname, results = results, results_true = results_true, v_input = v_input, v_mc = v_mc)
        @info "fits done - saved to  $fname"
    end

    results = nothing
    results_true = nothing
    v_mc = nothing 
end

# loglikelihoods_bf = [results[i].gof.loglikelihood_ml for i in 1:nSamples]
# loglikelihoods_true = [results_true[i].gof.loglikelihood_ml for i in 1:nSamples]
# teststatistic = 2 .* (loglikelihoods_true .- loglikelihoods_bf)
# t90 = quantile(teststatistic, 0.9)
# halflifes =  [results[i].halflife_0nbb for i in 1:nSamples]./1e26 
# stephist(teststatistic, nbins = 1000, fill = false, yscale = :log10,
#         label = "Test Statistic", xlabel = "Test Statistic", ylabel = "Occurrence", title = "Profile Likelihood Ratio Test")
# halflifes =  [results[i].halflife_0nbb for i in 1:nSamples]./1e26 
# # scatter(Measurements.value.(halflifes), -log.(likelihoods))
# stephist(Measurements.value.(halflifes), fill = true, nbins = 1000, xlabel = "Half-life (1e26 yr)", label = @sprintf("Fit results samples. median = %.3f",quantile(Measurements.value.(halflifes), 0.5)))
# xlims!(0, 2); ylims!(0, 4200)
# vline!([v_mc.halflife_0nbb], color = :red, label = "MC truth = $(v_mc.halflife_0nbb)e26 yr", linewidth = 2)


# # 
# # # # test:
i = 2
v_mc = (signal_strength = signal_strength[i], background = 5.2e-4)
events_mc = simspec(v_input, v_mc...; plotFlag = true)
result, report = fit_0nbb_mc(v_input, events_mc; uncertainty = false)
# plot_fit(report)