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

include("../utils/utils_plot.jl")

v_input = (exposure = 200, efficiency = 0.95, fit_window = 240, Qbb = 2039, Qbb_sigma = 2.5, sig_exp = 1.0e-26) 
fpath = "$(@__DIR__)/results/"
nSamples = 100000 
signal_mc = collect(0.0:0.1:5.0) 
t90 = zeros(length(signal_mc))

for (i, sig) in enumerate(signal_mc)
    v_mc = (signal_strength = sig, background = 5.2e-4)
    fname = fpath * @sprintf("mc0nbbfit_samples%i_sig%.2f_bck%.1e.jld2", nSamples, v_mc.signal_strength, v_mc.background)
    if isfile(fname) 
        @info "load $fname..."
        file = jldopen(fname)
        if v_input != file["v_input"] || v_mc != file["v_mc"]
            error("Input parameters do not match.")
        end
        results = file["results"]
        results_true =  file["results_true"]
    else
        @info "result files for half-life $sig doesnt exist - skip "
        continue
    end

    #remove not converged fits 
    convergedIdx  = unique(vcat(findall(x-> x == 1, [results[i].converged for i in 1:nSamples]), findall(x-> x == 1, [results_true[i].converged for i in 1:nSamples])))
    @printf("converged fits: %.2f %%", 100*length(convergedIdx)/nSamples)
    loglikelihoods_bf = [results[i].gof.loglikelihood_ml for i in 1:nSamples][convergedIdx]
    loglikelihoods_true = [results_true[i].gof.loglikelihood_ml for i in 1:nSamples][convergedIdx]
    teststatistic = 2 .* (loglikelihoods_true .- loglikelihoods_bf)
    t90[i] = quantile(teststatistic, 0.9)

    if i == length(signal_mc)
        fname = fpath * @sprintf("LikelihoodRatio_thresholds_samples%i.jld2", nSamples)
        jldsave(fname, signal_mc = signal_mc, t90 = t90)
    end
end

_def_plot()
hline(signal_mc, [2.7], linewidth = 1.5, label = "Expectation Wilk's theorem", color = :silver)
plot!(signal_mc, t90, linewidth = 2, label = "$nSamples samples each")
xlabel!(raw"$1/T_{1/2}^\mathrm{true}$" * @sprintf(" (%.0e yr)",v_input.sig_exp))
ylabel!("Likelihood ratio 90%-quantile")
ylims!(0, 3.5)
