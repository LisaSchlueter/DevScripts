import PhysicalConstants.CODATA2018: N_A, m_u # avogadro number, atomic mass constant
using Distributions 
using Interpolations, LinearAlgebra
using Plots
include("0nbb_likelihood.jl")

function simspec(v_input::NamedTuple{((:exposure, :efficiency, :fit_window, :Qbb, :Qbb_sigma, :sig_exp))}, signal_strength::Real, background::Real; plotFlag::Bool = false)

    # number of expected events in fit window
    Nevents_lambda = lambda_sig(v_input, signal_strength) + lambda_bck(v_input, background)
    Nevents = rand(Poisson(Nevents_lambda))

    # randomize energy values (spectrum): sample from peakshape 
    e_sample = ustrip.(collect(v_input.Qbb - v_input.fit_window/2:0.1:v_input.Qbb + v_input.fit_window/2))
    signal_pdf = peakshape_0nbb.(Ref(v_input), Ref(signal_strength), Ref(background), e_sample)
    signal_cdf = cumsum(signal_pdf)./maximum(cumsum(signal_pdf))
    interp_cdf_inv = linear_interpolation(signal_cdf, e_sample) # inverse cdf
    bandwidth = maximum(signal_cdf)-minimum(signal_cdf)
    rand_i = minimum(signal_cdf).+bandwidth.*rand(Nevents); 
    events_mc = interp_cdf_inv.(rand_i) 

    if plotFlag
        plt = stephist(events_mc, nbins = ifelse(Nevents > 100, 100, Nevents), fill = true, label = "Simulated spectrum")
        xlims!(ustrip(v_input.Qbb - v_input.fit_window/2), ustrip(v_input.Qbb + v_input.fit_window/2))
        display(plt)
    end 
    return events_mc
end

# v_input = (exposure = 200, efficiency = 0.95, fit_window = 240, Qbb = 2039, Qbb_sigma = 2.5)
# v_test = (signal_strength = 2, background = 5.2e-4)
# simspec(v_input, v_test...)