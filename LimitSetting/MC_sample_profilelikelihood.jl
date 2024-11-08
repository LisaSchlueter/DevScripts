import PhysicalConstants.CODATA2018: N_A, m_u # avogadro number, atomic mass constant
using LegendDataManagement
using Unitful
using CairoMakie
using Distributions 
using Optim, ForwardDiff
using Interpolations, LinearAlgebra
using Measurements, Measures
using Printf
using JLD2

include("funcs/0nbb_likelihood.jl")
include("funcs/0nbb_simspec.jl")
include("funcs/0nbb_fit.jl")
include("../utils/utils_plot.jl")
nSamples = 100000
reCalc = false
v_input = (exposure = 200, efficiency = 0.95, fit_window = 240, Qbb = 2039.0, Qbb_sigma = 2.5, sig_exp = 1.0e-26) 
v_mc = (signal_strength = 0, background = 5.2e-4) # MC truth

# prepare for saving results 
fpath = "$(@__DIR__)/results/"
fname = fpath * @sprintf("mc0nbbfit_profile_samples%i_sig%.2f_bck%.1e.jld2", nSamples, v_mc.signal_strength, v_mc.background)
if isfile(fname) && reCalc == false 
    @info "file $fname exists already - load results"
    file = jldopen(fname)
    signal_strength_scan = file["signal_strength_scan"]
    loglikelihood_profile = file["loglikelihood_profile"]
    results_bf = file["results"]
    loglikelihood_bf = file["loglikelihood_bf"]
    close(file)
else
    @info "do fits for $fname"
    results = Vector{Any}(undef, nSamples)
    signal_strength_scan = collect(0:0.02:1)
    loglikelihood_profile = zeros(nSamples, length(signal_strength_scan))

    for i in 1:nSamples
        local events_mc = simspec(v_input, v_mc...; plotFlag = false)
        local result, report = fit_0nbb_mc(v_input, events_mc; uncertainty = false)
        for (j, sig) in enumerate(signal_strength_scan)
            local result_scan, _ = fit_0nbb_mc(v_input, events_mc; fixsignal = sig, uncertainty = false)
            loglikelihood_profile[i, j] = result_scan.gof.loglikelihood_ml
        end
        results[i] = result
    end
    loglikelihood_bf = [results[i].gof.loglikelihood_ml for i in 1:nSamples]

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

    jldsave(fname, results = results, loglikelihood_profile = loglikelihood_profile, loglikelihood_bf = loglikelihood_bf,
        v_input = v_input, v_mc = v_mc, signal_strength_scan = signal_strength_scan)
    @info "fits done - saved to  $fname"
end

# Profile log likelihood: plot a few example samples
_plot_theme()
f = Figure(size = (600, 400), dpi = 150)
ax = Axis(f[1, 1], xlabel = latexstring("\\textrm{Signal strength } 1/T_{1/2} \\, (10^{$(round(Int,log10(v_input.sig_exp)))} \\, \\textrm{yr}^{-1})"), 
            ylabel = latexstring("2 \\cdot (-\\ln\\mathcal{L}_\\textrm{profile} + \\ln\\mathcal{L}_\\textrm{profile}^\\textrm{best fit} ) "),
            title = "Profile likelihood ratio for MC spectra w/o signal")
loglikelihood_profile_norm = 2 .* (loglikelihood_profile .- loglikelihood_bf )  # also our test statistics 
psample1 = scatterlines!(ax, signal_strength_scan, loglikelihood_profile_norm[1, :], color = :red, markersize = 2.5, linewidth = 2.5)
psample2 = scatterlines!(ax, signal_strength_scan, loglikelihood_profile_norm[3, :], color = :dodgerblue, markersize = 3, linewidth = 2.5)
xlims!(0, 1)
axislegend(ax, [psample1, psample2], ["MC sample 1", "MC sample 2"]; position = :lt)
f

# optional load threshold values for 90% C.L. from file
fname_threshold = fpath * @sprintf("LikelihoodRatio_thresholds_samples%i.jld2", nSamples)
file_t90= jldopen(fname_threshold)
t90_inter = file_t90["interlin"]
close(file_t90)
xinter = range(0, 1, length = 1000)
p90 = lines!(ax, xinter, t90_inter.(xinter), color = :black, linewidth = 2, linestyle = :dash)
axislegend(ax, [psample1, psample2, p90], ["MC sample 1", "MC sample 2", "90% C.L. cut threshold"]; position = :lt)
f
