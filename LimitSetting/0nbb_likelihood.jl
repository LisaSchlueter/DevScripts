import PhysicalConstants.CODATA2018: N_A, m_u # avogadro number, atomic mass constant
using Unitful
using Distributions 

"""
    lambda_sig(v_input::NamedTuple{((:exposure, :efficiency, :fit_window, :Qbb, :Qbb_sigma))}, halflife_0nbb::Real = 1..0)
This function calculates the signal strength for a given 0nbb halflife, exposure, efficiency, fit_window, Qbb, and Qbb_sigma.
The signal strength is the expected number of 0nbb events. 
!IMPORTANT! unit of signal_strength is 1e-26 years.
"""
function lambda_sig(v_input::NamedTuple{((:exposure, :efficiency, :fit_window, :Qbb, :Qbb_sigma, :sig_exp))}, signal_strength::Real = 1.0)
    # AtomicMass_Ge76 = 75.921402726 *  m_u  # weight of 1 Ge76 atom in [kg]
    # Ge_MolarMass =  AtomicMass_Ge76 * N_A  # weight of 1 mol of Ge76 atom in [kg]
    return  signal_strength * v_input.sig_exp * ustrip(v_input.exposure * v_input.efficiency * log(2) / (75.921402726 *  m_u))
end

"""
    background_strength(v_input::NamedTuple{((:exposure, :efficiency, :fit_window, :Qbb, :Qbb_sigma))}, background::Real = 1e-3)
This function calculates the background strength for a given background index [counts / (kg year)], exposure, fit_window
The background strength is the expected number of background in the fit window
"""
function lambda_bck(v_input::NamedTuple{((:exposure, :efficiency, :fit_window, :Qbb, :Qbb_sigma, :sig_exp))}, background::Real = 1e-3)
    return ustrip(background * v_input.exposure * v_input.fit_window) 
end

"""
    peakshape_0nbb(v_input::NamedTuple{((:exposure, :efficiency, :fit_window, :Qbb, :Qbb_sigma))}, signal_strength::Real, background::Real, energy::Real)
This function calculates the peakshape for a given signal_strength (1/0nbb halflife [1e-26 yr]), background, and energy.
The peakshape is the probability density function of the 0nbb signal and background.
"""
function peakshape_0nbb(v_input::NamedTuple{((:exposure, :efficiency, :fit_window, :Qbb, :Qbb_sigma, :sig_exp))}, signal_strength::Real, background::Real, energy::Real)
    N_signal = lambda_sig(v_input, signal_strength)
    N_background = lambda_bck(v_input, background)
    return 1 / (N_signal + N_background) .* (N_signal .* pdf(Normal(ustrip(v_input.Qbb), ustrip(v_input.Qbb_sigma)), energy) .+  N_background   / ustrip(v_input.fit_window)) 
end

function likelihood_partition(signal_strength::Real , background::Real, v_input::NamedTuple{((:exposure, :efficiency, :fit_window, :Qbb, :Qbb_sigma, :sig_exp))}, Energy::Vector{<:Real} = 2039.0 .+ [-1.0,0.0, 1.0])
    if lambda_sig(v_input, signal_strength) < 0 || lambda_bck(v_input, background) < 0
        return 0
    else 
        poisson_term = pdf(Poisson(lambda_sig(v_input, signal_strength) + lambda_bck(v_input, background)), length(Energy))
        likelihood = prod(peakshape_0nbb.(Ref(v_input), signal_strength, background, Energy))
        likelihood_ext = poisson_term * likelihood
        return likelihood_ext
    end
end