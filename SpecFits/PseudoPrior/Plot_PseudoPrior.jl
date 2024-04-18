#= plot pseudo-priors from peaks fits

standard_pseudo_prior = NamedTupleDist(
    μ = ifelse(fixed_position, ConstValueDist(ps.peak_pos), Uniform(ps.peak_pos-10, ps.peak_pos+10)),
    σ = weibull_from_mx(ps.peak_sigma, 2*ps.peak_sigma),
    n = weibull_from_mx(ps.peak_counts, 2*ps.peak_counts),
    step_amplitude = weibull_from_mx(ps.mean_background_step, ps.mean_background_step + 5*ps.mean_background_std),
    skew_fraction = ifelse(low_e_tail, truncated(weibull_from_mx(0.01, 0.05), 0.0, 0.25), ConstValueDist(0.0)),
    skew_width = ifelse(low_e_tail, weibull_from_mx(0.001, 1e-2), ConstValueDist(1.0)),
    background = weibull_from_mx(ps.mean_background, ps.mean_background + 5*ps.mean_background_std),
)
=#
using Distributions
using Plots
using LegendSpecFits
# 1. Weibull: used on σ, n, step_amplitude, skew_widths, and background 
# play around with distribution 
weibull_median = 0.001   # median of distribution 
weibull_quantile = 1e-2  # quantile for which: integral of weibfull from 0 to weibull_quantile =  0.6827 --> "how streched distribtion is toward larger values" 
d = weibull_from_mx(weibull_median, weibull_quantile)
quantile(d, 0.5) # cross-check median
quantile(d, 0.6827) # cross-check quantile

function gety_weibull(median, quantile,x)
    d = weibull_from_mx(median, quantile)
    y = pdf.(d,x)
    return y
end

# plot different weibulls 
x = exp10.(range(-5, stop= -1, length=100))
weibull_quantiles = [0.002,0.003,0.005,0.01,0.02]#[0.002:0.002:0.01...,0.02:0.01:0.1...]
p = plot()
for i in eachindex(weibull_quantiles)
    plot!(x,gety_weibull.(Ref(weibull_median), weibull_quantiles[i], Ref(x)),yscale=:log10,legend=true, label = "q = $(weibull_quantiles[i])", ylabel = "pdf", xlabel = "Fit parameter", title = "Weibull distribution", lw = 2, legendfontsize = 10, xguidefontsize = 10, yguidefontsize = 10, xlabelfontsize = 10, ylabelfontsize = 10, grid = true)
end 
plot!(p,legendtitle = "median = $weibull_median") 
ylims!(1e-2,1e4)


# plot truncated weibull. used in skew_fraction
quantile_truncated = 0.05
median_truncated = 0.01
trunc = [0,0.25]
d_trunc = truncated(weibull_from_mx(median_truncated, quantile_truncated), trunc[1], trunc[2])
y = pdf.(d_trunc,x)
p_trunc = plot(x,y,yscale=:log10,legend=true, label = "q = $quantile_truncated, truncated [$(trunc[1]), $(trunc[2])]", ylabel = "pdf", xlabel = "Fit parameter", 
        title = "Truncated Weibull distribution", lw = 2, legendfontsize = 10, xguidefontsize = 10, yguidefontsize = 10, 
        xlabelfontsize = 10, ylabelfontsize = 10, grid = true)
plot!(p_trunc,legendtitle = "median = $median_truncated") 

