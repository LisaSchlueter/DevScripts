
using CairoMakie
using JLD2
using Printf, LaTeXStrings
using Measures
using Interpolations, Roots
include("../utils/utils_plot.jl")

# load likelihood ratio thresholds for 90% C.L. from file   
fpath = "$(@__DIR__)/results/"
fnamet90 = fpath * @sprintf("LikelihoodRatio_thresholds_samples%i.jld2", 100000)
file_t90 = jldopen(fnamet90)
t90_inter = file_t90["interlin"]
# t90 = filet90["t90"]
# wilks_signal_mc = filet90["signal_mc"]
close(file_t90)

# load profile likelihoods for mc simulations w/o signal from file 
nSamples = 100000
v_mc = (signal_strength = 0, background = 5.2e-4)
fnameS = fpath * @sprintf("mc0nbbfit_profile_samples%i_sig%.2f_bck%.1e.jld2", nSamples, v_mc.signal_strength, v_mc.background)
filet90 = jldopen(fnameS)
loglikelihood_profile = filet90["loglikelihood_profile"]
loglikelihood_bf = filet90["loglikelihood_bf"]
signal_strength_scan = filet90["signal_strength_scan"]
results = filet90["results"]
close(filet90)

# calcualte test statistic for each sample 
teststatistic = 2 .* (loglikelihood_profile .- loglikelihood_bf )

# get confidence intervals for each sample
function conf_interval(sig::Vector{Float64}, t::Vector{Float64}, t90_inter::Interpolations.Extrapolation)
    interp_sample =  linear_interpolation(sig, t)
    lossfunc(x) =  t90_inter.(x) .- interp_sample.(x)
    try
        return find_zeros(lossfunc, minimum(sig), maximum(sig))    
    catch
        return NaN
    end
end 
ConfInter = [conf_interval(signal_strength_scan, teststatistic[i,:], t90_inter) for i in 1:nSamples] 
upperlims = vcat(filter(x -> length(x)==1, ConfInter)...) # 1 intersection
closed = filter(x -> length(x)>1, ConfInter) # 2 intersections
# sanity checks 
if any(length.(ConfInter) .> 2)
    println("Warning: some samples have more than 2 intersections")
end
if any(.!isfinite.(upperlims)) 
    println("Warning: some samples have no intersection")
end


#####  PLOTS 
_plot_theme()
ppath = "$(@__DIR__)/plots/Limits/"
if !ispath(ppath)
    mkpath(ppath)
end

# Plot 1: test statistics (profile likelihood ratio) for a few example samples
sig_exp = 1e-26
f = Figure(size = (600, 400), dpi = 150)
ax = Axis(f[1, 1], xlabel = latexstring("\\textrm{Signal strength } 1/T_{1/2} \\, (10^{$(round(Int,log10(sig_exp)))} \\, \\textrm{yr}^{-1})"), 
            ylabel = latexstring("2 \\cdot (-\\ln\\mathcal{L}_\\textrm{profile} + \\ln\\mathcal{L}_\\textrm{profile}^\\textrm{best fit} ) "),
            title = "Profile likelihood ratio for MC spectra w/o signal",
            xminorgridvisible = true,  xminorticks = IntervalsBetween(5), xticks = 0:0.5:1)

psample1 = scatterlines!(ax, signal_strength_scan, teststatistic[1, :], color = :red, markersize = 2.5, linewidth = 2.5)
# psample2 = scatterlines!(ax, signal_strength_scan, teststatistic[3, :], color = :dodgerblue, markersize = 3, linewidth = 2.5)
psample3 = scatterlines!(ax, signal_strength_scan, teststatistic[18, :], color = :green, markersize = 3, linewidth = 2.5)
xlims!(0, 1)
xinter = range(0, 1, length = 1000)
p90 = lines!(ax, xinter, t90_inter.(xinter), color = :black, linewidth = 2, linestyle = :dash)
axislegend(ax, [psample1,  psample3, p90], 
            [latexstring("1/T_{1/2} \\, < \\,$(round(ConfInter[1][1], digits = 2)) \\cdot 10^{$(round(Int,log10(sig_exp)))} \\, \\textrm{yr}^{-1}"), 
            latexstring("$(round(ConfInter[18][1], digits = 2)) \\, < \\, 1/T_{1/2} \\, < \\, $(round(ConfInter[18][2], digits = 2))  \\cdot 10^{$(round(Int,log10(sig_exp)))} \\, \\textrm{yr}^{-1}"),
            "90% C.L. cut threshold"]; position = :lt)
f
save(ppath * "MC_ProfileLikelihoods_MCexamples.png", f)



# Plot 2: distribution of upper limits (MC truth = no signal)
f2 = Figure(size = (600, 400), dpi = 150)
ax2 = Axis(f2[1, 1], xlabel = latexstring("\\textrm{Upper limit on signal strength } 1/T_{1/2} \\, (10^{$(round(Int,log10(sig_exp)))} \\, \\textrm{yr}^{-1})"), 
            title = "Upper limits derived from MC spectra w/o signal, 1e$(round(Int, log10(nSamples))) samples",
            xminorgridvisible = true,  xminorticks = IntervalsBetween(2), 
            xticks = 0:0.2:1, yticks = 1000:1000:10000, titlesize = 18)

hist!(ax2, upperlims, label = "Simulated spectra without signal", color = :orchid2, bins = 100, strokecolor = :orchid2 )
ylims!(0, nothing)
xlims!(0, nothing)
f2
save(ppath * "UpperLimits_MCdistribution.png", f2)


