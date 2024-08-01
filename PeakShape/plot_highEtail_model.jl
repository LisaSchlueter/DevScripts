
using Plots
using Distributions
using LegendSpecFits


µ = 1600.0
σ = 2.0
x = collect(1590:0.1:1610)
n = 100


skew_width = 0.01
skew_frac = 0.01

skew = µ * 0.01
peakshape_lowE = skew -> n*skew_frac.*LegendSpecFits.ex_gauss_pdf.(-x, -μ, σ, skew)
peakshape_highE = skew -> n*skew_frac.*LegendSpecFits.ex_gauss_pdf.(x, μ, σ, skew)
signal = LegendSpecFits.signal_peakshape.(x, μ, σ, n*(1-skew_frac), 0)


plot(x, peakshape_lowE(skew), label="Low-energy tail = $skew", lw=2)
plot!(x, peakshape_highE(skew), label="High-energy tail = $skew", lw=2)
plot!(x, LegendSpecFits.signal_peakshape.(x, μ, σ, n, skew_frac), label="Signal", lw=2, yaxis = :log10)
# plot!(x, LegendSpecFits.lowEtail_peakshape.(x, μ, σ, n, skew_frac, skew_width), label="Low-energy tail ", lw=2, yaxis = :log10, linestyle = :dash)


# x::Real, μ::Real, σ::Real, n::Real,
# step_amplitude::Real, skew_fraction::Real, skew_width::Real,
# background::Real