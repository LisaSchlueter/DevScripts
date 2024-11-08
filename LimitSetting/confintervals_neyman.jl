# Calculate classical confidence belts for the signal strength using the Neyman construction
# not based on profile likelihood, but full likelihood

using Interpolations#, LinearAlgebra
using Statistics
using Measures
using Printf, LaTeXStrings
using JLD2
using CairoMakie 
# using Plots
# # include("0nbb_peakshapes.jl")
# include("funcs/0nbb_simspec.jl")
# include("funcs/0nbb_fit.jl")
# include("../utils/utils_plot.jl")

v_input = (exposure = 200, efficiency = 0.95, fit_window = 240, Qbb = 2039, Qbb_bias = 0, Qbb_sigma = 2.5, sig_exp = 1.0e-26) 
nSamples = 100000

fpath = "$(@__DIR__)/results/"
signal_mc = collect(0.0:0.1:5.0)
t90 = zeros(length(signal_mc))
signal_fit = Vector{Vector{Float64}}(undef, length(signal_mc))

for (i, sig) in enumerate(signal_mc)
    v_mc = (signal_strength = sig, background = 5.2e-4)
    fname = fpath * @sprintf("mc0nbbfit_samples%i_sig%.2f_bck%.1e.jld2", nSamples, v_mc.signal_strength, v_mc.background)
    if isfile(fname) 
        @info "load $fname..."
        file = jldopen(fname)
        # if v_input != file["v_input"] || v_mc != file["v_mc"]
        #     error("Input parameters do not match.")
        # end
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
    signal_fit[i] =  [results[s].signal_strength for s in 1:nSamples][convergedIdx]
end

fs = 20
PlotTheme = Theme(
    size = (600, 420),
    fontsize = fs,
    dpi = 300,
    margin = 3mm,
    Fonts = (
        regular="Helvetica", math="Helvetica"),
    Axis = (
        titlefont = :regular,
        xgridvisible = true,
        ygridvisible = true,
        xlabelsize= fs + 6,
        ylabelsize= fs + 6,
        xtickalign=1,
        ytickalign=1,
    ),
    Legend = (
        framecolor = :silver,
        labelsize = fs + 2,
        position = :lt,
    )
)
set_theme!(PlotTheme)

inter_l = LinearInterpolation(quantile.(signal_fit, 0.05), signal_mc, extrapolation_bc = Line())
inter_h = LinearInterpolation(quantile.(signal_fit, 0.95), signal_mc, extrapolation_bc = Line())
inter_x = range(minimum(quantile.(signal_fit, 0.05)), stop = maximum(quantile.(signal_fit, 0.95)), length = 100)
lim05 = inter_l.(inter_x)
lim95 = inter_h.(inter_x)

f = Figure(size = (800, 600))
ax = Axis(f[1, 1], xlabel = "Measured signal strength 1/T (1e-26 yr)", ylabel = "True (MC) signal strength 1/T (1e-26 yr)")
pone = band!(quantile.(signal_fit, 0.9), zeros(length(signal_mc)), signal_mc, 
            color = :dodgerblue, alpha = 0.3,
            label = "One-sided confidence belt (upper limits) 90% C.L.")
ptwo = band!(inter_x, lim95, lim05, 
            color = :red, alpha = 0.3,
            label = "Two-sided confidence belt 90% C.L.")


lines!(quantile.(signal_fit, 0.9), signal_mc, 
            color = :dodgerblue, alpha = 1, linewidth = 2)
lines!(quantile.(signal_fit, 0.05), signal_mc, 
        color = :red2, alpha = 1, linewidth = 2)
lines!(quantile.(signal_fit, 0.95), signal_mc, 
        color = :red2, alpha = 1, linewidth = 2)
ylims!(minimum(signal_mc), maximum(signal_mc))
xlims!(0, maximum(signal_mc))

axislegend(ax, [pone, ptwo], ["One-sided confidence belt (upper limit) 90% C.L.", "Two-sided confidence belt  90% C.L."], position = :lt)

f
