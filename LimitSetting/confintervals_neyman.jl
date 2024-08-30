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

# include("0nbb_peakshapes.jl")
include("0nbb_simspec.jl")
include("0nbb_fit.jl")
include("../utils/utils_plot.jl")

v_input = (exposure = 200, efficiency = 0.95, fit_window = 240, Qbb = 2039, Qbb_sigma = 2.5, sig_exp = 1.0e-26) 
nSamples = 1000#200000
fpath = "$(@__DIR__)/results/"
signal_mc = collect(0.0:1:10.0)

t90 = zeros(length(signal_mc))
signal_fit = Vector{Vector{Float64}}(undef, length(signal_mc))
 fpath = "$(@__DIR__)/results/"

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
        # @info "...done"
    else
        @info "result files for signal strength $sig doesnt exist - skip "
        continue
    end
    #remove not converged fits 
    convergedIdx  = unique(vcat(findall(x-> x == 1, [results[i].converged for i in 1:nSamples]), findall(x-> x == 1, [results_true[i].converged for i in 1:nSamples])))
    @printf("converged fits: %.2f %%", 100*length(convergedIdx)/nSamples)

    loglikelihoods_bf = [results[s].gof.loglikelihood_ml for s in 1:nSamples][convergedIdx]
    loglikelihoods_true = [results_true[s].gof.loglikelihood_ml for s in 1:nSamples][convergedIdx]
    teststatistic = 2 .* (loglikelihoods_true .- loglikelihoods_bf)
    t90[i] = quantile(teststatistic, 0.9)

    # halflifes_fit[i] = filter(x -> x<100, [results[s].halflife_0nbb./1e26 for s in 1:nSamples][convergedIdx])
    signal_fit[i] =  [results[s].signal_strength for s in 1:nSamples][convergedIdx]
end
# halflifes_mc_inv = 1 ./ halflifes_mc
# halflifes_fit_inv = [1 ./ halflifes_fit[i] for i in 1:length(halflifes_fit)]


_def_plot()
plot(quantile.(signal_fit, 0.9), signal_mc, 
            linewidth = 2.5, color = :dodgerblue, 
            label = "One-sided confidence belt (upper limits) 90% C.L.",
            legend = :outertop)
plot!(quantile.(signal_fit, 0.05), signal_mc, linewidth = 2, linestyle = :dash, color = :red, label = "Two-sided confidence belt 90% C.L.")
plot!(quantile.(signal_fit, 0.95), signal_mc, linewidth = 2, linestyle = :dash, color = :red, label = false)
ylims!(minimum(signal_mc), maximum(signal_mc))
xlabel!("1/T measured (1e-26 yr)")
ylabel!("1/T true (1e-26 yr)")

# stephist(halflifes_fit[1])
# minimum(halflifes_fit[1])
# maximum(halflifes_fit[1])
