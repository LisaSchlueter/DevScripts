# This script does:
# 1. Load fits results from  MC sample fits for different signals strength
# 2. Calculate the profile likelihood ratio test statistic for each signal strength
# 3. Calculate the 90% quantile of the test statistic
# 4. Compare the observed 90% quantile with the expected 90% quantile from the chi2 distribution --> Validate coverage 
# Needed to calculate our cutoff threshold for the profile likelihood ratio test
import PhysicalConstants.CODATA2018: N_A, m_u # avogadro number, atomic mass constant
using LegendDataManagement
using Unitful
using CairoMakie
using Distributions 
using Optim, ForwardDiff
using Interpolations, LinearAlgebra
using Measurements, Measures
using Printf, LaTeXStrings
using JLD2
using Interpolations
include("../utils/utils_plot.jl")
reCalc = true 

v_input = (exposure = 200, efficiency = 0.95, fit_window = 240, Qbb = 2039, Qbb_sigma = 2.5, sig_exp = 1.0e-26) 
fpath = "$(@__DIR__)/results/"
nSamples = 100000 
signal_mc = collect(0.0:0.1:5.0) 
t90 = zeros(length(signal_mc))

fname_threshold = fpath * @sprintf("LikelihoodRatio_thresholds_samples%i.jld2", nSamples)
ppath = "$(@__DIR__)/plots/LikelihoodRatio/"
if !ispath(ppath)
    mkpath(ppath)
end
_plot_theme()

if isfile(fname_threshold) .&& reCalc == false
    file = jldopen(fname_threshold)
    signal_mc = file["signal_mc"]
    t90 = file["t90"]
    @info "read threshold values t90 from file $fname_threshold"
    close(file)
else
    for (i, sig) in enumerate(signal_mc)
        v_mc = (signal_strength = sig, background = 5.2e-4)
        local fname = fpath * @sprintf("mc0nbbfit_samples%i_sig%.2f_bck%.1e.jld2", nSamples, v_mc.signal_strength, v_mc.background)
        if isfile(fname) 
            @info "load $fname..."
            local file = jldopen(fname)
            if v_input != file["v_input"] || v_mc != file["v_mc"]
                error("Input parameters do not match.")
            end
            results = file["results"]
            results_true =  file["results_true"]
            close(file)
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
        t90_expected = round(quantile(Chisq(1), 0.9), digits = 2)
  
        f = Figure(size = (600, 400), dpi = 150)
        ax = Axis(f[1, 1], 
            yscale = log10,
            xlabel = L"$2 \cdot (\textrm{neg.} \ln\mathcal{L}(1/T_\textrm{true}) - \textrm{neg.} \ln\mathcal{L}(1/T_\textrm{bf}))$",
            ylabel = "Occurrence", 
            title = latexstring("\\textrm{Test statistic for } 1/T_{1/2} = $sig \\cdot 10^{$(round(Int,v_input.sig_exp))} \\, \\textrm{yr}^{-1}"))

        h1 = hist!(ax, teststatistic, color = :silver, strokewidth = 0, bins = 1000 )
        ylims!(1, nothing)
        xlims!(ifelse(minimum(teststatistic) < 0, minimum(teststatistic), 0), quantile(teststatistic, 0.999))

        ltheo = vlines!([t90_expected], color = :red, label = "Expected 90%-quantile: $t90_expected", linewidth = 2, alpha = 0.5)
        ldata = vlines!([t90[i]], color = :navy, label = "Observed 90%-quantile: $(round(t90[i], digits = 2))", linewidth = 2, alpha = 0.5) 
        axislegend(ax,[h1, ltheo, ldata], ["MC $(nSamples) samples", "Expected 90%-quantile: $t90_expected", "Observed 90%-quantile: $(round(t90[i], digits = 2))" ]; position = :rt)
        save(ppath * @sprintf("LikelihoodRatio_samples%i_sig%.2f_bck%.1e.png", nSamples, v_mc.signal_strength, v_mc.background), f)
    end 
    interlin = LinearInterpolation(signal_mc, t90,  extrapolation_bc = Line())
    jldsave(fname_threshold, signal_mc = signal_mc, t90 = t90, interlin = interlin)
end 

# plot 
f = Figure(size = (600, 400), dpi = 150)
ax = Axis(f[1, 1], 
            xlabel = latexstring("1/T_{1/2}^\\mathrm{true} \\, (10^{$(round(Int,log10(v_input.sig_exp)))} \\, \\textrm{yr}^{-1})") ,#* @sprintf(" (%.0e yr)",v_input.sig_exp), 
            ylabel =  latexstring("t_{90} \\, \\textrm{(threshold for 90% C.L.)}"),
            title = "90% C.L. threshold for profile likelihood ratio test statistic")
pref = hlines!(ax, [2.7], linewidth = 1.5, label = "Expectation: 2.71 (90% C.L., 1 dof)", color = :red2)
pdata = scatter!(ax, signal_mc, t90, label = "$nSamples samples each", color = :dodgerblue, markersize = 10)

axislegend(ax, [pref, pdata], ["Wilk's theorem expect. 2.71 (90% C.L., 1 dof)", 
                        "MC sampled threshold, 1e$(round(Int,log10(nSamples))) samples each"], 
            position = :lt)
ylims!(1.5, 3.5)
save(ppath * @sprintf("LikelihoodRatio_t90_Wilks_samples%i.png", nSamples), f)
f
