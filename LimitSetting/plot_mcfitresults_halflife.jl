import PhysicalConstants.CODATA2018: N_A, m_u # avogadro number, atomic mass constant
using LegendDataManagement
using Unitful
using Plots
using Distributions 
using Optim, ForwardDiff
using Interpolations, LinearAlgebra
using Measurements
using Printf, LaTeXStrings
using JLD2

# include("0nbb_peakshapes.jl")
include("0nbb_simspec.jl")
include("0nbb_fit.jl")

v_input = (exposure = 200, efficiency = 0.95, fit_window = 240, Qbb = 2039, Qbb_sigma = 2.5, sig_exp = 1e-26) 
nSamples = 1000 # 200000

signal_strength  = collect(0:1:10) 

signal_median = zeros(length(signal_strength))

for (i, sig) in enumerate(signal_strength)
    v_mc = (signal_strength = sig, background = 5.2e-4)
    fpath = "$(@__DIR__)/results/"
    fname = fpath * @sprintf("mc0nbbfit_samples%i_sig%.2f_bck%.1e.jld2", nSamples, v_mc.signal_strength, v_mc.background)
    if isfile(fname) 
        @info "load $fname..."
        file = jldopen(fname)
        if v_input != file["v_input"] || v_mc != file["v_mc"]
            error("Input parameters do not match.")
        end
        results = file["results"]
        close(file)
    else
        @info "files not found $fname "
        continue
    end

    #remove not converged fits 
    convergedIdx  = findall(x-> x == 1, [results[i].converged for i in 1:nSamples])
    @printf("converged fits: %.2f %%", 100*length(convergedIdx)/nSamples)
    
    signal_fit = [results[s].signal_strength for s in 1:nSamples][convergedIdx] 
    signal_fit = filter(x -> x < 100, signal_fit)

    fs = 14
    default(foreground_color_legend = :silver,
    background_color_legend = :white,
    grid = :off,
    framestyle = :semi,
    xtickfontsize = fs-2,
    ytickfontsize = fs-2,
    legendfontsize = fs-2,
    xlabelfontsize = fs + 2,
    ylabelfontsize = fs + 2)

    signal_median[i] = quantile(Measurements.value.(signal_fit), 0.5)
    stephist(signal_fit, fill = true, color = :silver,
            size = (520, 370),
            xlabel = @sprintf("Signal strength (%.0e 1/yr)", v_input.sig_exp),
            ylabel = "Occurrence",
            label = @sprintf("Fit results samples; median = %.2f", quantile(Measurements.value.(signal_fit), 0.5)) )
    ylims!(0, ylims()[2])
    vline!([v_mc.signal_strength], color = :red, label = "MC truth: T1/2 = $(v_mc.signal_strength)", linewidth = 2.5)
   
    ppath = "$(@__DIR__)/plots/mcfits/"
    if !isdir(ppath)
        mkdir(ppath)
    end
    pname = ppath * @sprintf("mc0nbbfit_samples%i_sig%.2f_bck%.1e.png", nSamples, v_mc.signal_strength, v_mc.background)
    @info "\n save plot to $pname"
    savefig(pname)
    signal_fit = nothing 
    file = nothing 
end


plot(signal_median, signal_median, linewidth = 2.5, color = :silver, label = "Median fit results = MC truth")
scatter!(signal_strength, signal_median, linewidth = 2.5, color = :black, label = @sprintf("Median of %.1e sample fits", nSamples))
ylabel!(raw"$\langle1/T_{1/2}^\mathrm{fit}\rangle$" * @sprintf(" (%.0e yr)",v_input.sig_exp))
xlabel!(raw"$1/T_{1/2}^\mathrm{true}$" * @sprintf(" (%.0e yr)",v_input.sig_exp))
ppath = "$(@__DIR__)/plots/mcfits/"
savefig(ppath * "mc0nbbfitresults_signal_mctruth_vs_median_samples$nSamples.png")

