
using Plots
using JLD2
using Printf
using Interpolations, Roots
include("../utils/utils_plot.jl")

fpath = "$(@__DIR__)/results/"

# load likelihood ratio thresholds for 90% C.L. from file   
fnameWilks = fpath * @sprintf("LikelihoodRatio_thresholds_samples%i.jld2", 100000)
# jldsave(fname, signal_mc = signal_mc, t90 = t90)
fileW = jldopen(fnameWilks)
wilks_t90 = fileW["t90"]
wilks_signal_mc = fileW["signal_mc"]
close(fileW)

# load profile likelihood for mc simulations w/o signal from file 
nSamples = 100000
v_mc = (signal_strength = 0, background = 5.2e-4)
fnameS = fpath * @sprintf("mc0nbbfit_profile_samples%i_sig%.2f_bck%.1e.jld2", nSamples, v_mc.signal_strength, v_mc.background)
fileS = jldopen(fnameS)
loglikelihood_scan = fileS["loglikelihood_scan"]
signal_strength_scan = fileS["signal_strength_scan"]
results = fileS["results"]
close(fileS)
loglikelihood_bf = [results[i].gof.loglikelihood_ml for i in 1:nSamples]
loglikeratio = 2 .* (loglikelihood_scan .- loglikelihood_bf )

# plot 1 example of profile likelihood (ratio)
_def_plot()
plot(wilks_signal_mc, wilks_t90, linewidth = 3, label = "90% C.L. threshold ", color = :black)
plot!(signal_strength_scan, loglikeratio[1,:], linewidth = 2, linestyle = :dash, label = "Simulated spectrum w/o signal I", color = :red)
plot!(signal_strength_scan, loglikeratio[3,:], linewidth = 2, linestyle = :dashdot, label = "Simulated spectrum w/o signal II", color = :green)
plot!(signal_strength_scan, loglikeratio[20,:], linewidth = 2, linestyle = :dot, label = "Simulated spectrum w/o signal III", color = :orange)
ylims!(0, 3.5)
xlims!(0, 1.5)
xlabel!(raw"$1/T_{1/2}$" * @sprintf(" (%.0e yr)", 1e-26))
ylabel!("Likelihood ratio")


interp_W =  linear_interpolation(wilks_signal_mc, wilks_t90)
upperlimit  = zeros(nSamples)
for i in 1:nSamples
    local interp_sample =  linear_interpolation(signal_strength_scan, loglikeratio[i,:])
    local lossfunc(x) =  interp_W.(x) .- interp_sample.(x)
    try 
    upperlimit[i] = find_zero(lossfunc, (0.0, 1.0))
    catch e
       # println("Error: ", e)
        upperlimit[i] = NaN
    end
end
sum(isnan.(upperlimit)) / nSamples

stephist(1 ./ upperlimit, label = "Simulated spectra without signal", fill = true, color = :silver, nbins = 1000)
xlabel!(raw"90% C.L. Upper limit $T_{1/2}$" * @sprintf(" (%.0e yr)", 1e26))
ylims!(0, ylims()[2])