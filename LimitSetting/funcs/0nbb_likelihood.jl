import PhysicalConstants.CODATA2018: N_A, m_u # avogadro number, atomic mass constant
using Unitful
using Distributions 
using Measurements: value as mvalue
"""
    lambda_sig(v_input::NamedTuple{((:exposure, :efficiency, :fit_window, :Qbb, :Qbb_sigma))}, halflife_0nbb::Real = 1..0)
This function calculates the signal strength for a given 0nbb halflife, exposure, efficiency, fit_window, Qbb, and Qbb_sigma.
The signal strength is the expected number of 0nbb events. 
!IMPORTANT! unit of signal_strength is 1e-26 years.
"""
function lambda_sig(v_input::NamedTuple{((:exposure, :efficiency, :fit_window, :Qbb, :Qbb_bias, :Qbb_sigma, :sig_exp))}, signal_strength::Union{Real, Unitful.AbstractQuantity} = 1.0)
    # AtomicMass_Ge76 = 75.921402726 *  m_u  # weight of 1 Ge76 atom in [kg]
    # Ge_MolarMass =  AtomicMass_Ge76 * N_A  # weight of 1 mol of Ge76 atom in [kg]
    if !isa(signal_strength, Unitful.AbstractQuantity)
        atomic_mass_u = ustrip(m_u)
    else
        atomic_mass_u = m_u
    end
    return  mvalue(signal_strength * v_input.sig_exp * v_input.exposure * v_input.efficiency * log(2) / (75.921402726 *  atomic_mass_u))
end

"""
    background_strength(v_input::NamedTuple{((:exposure, :efficiency, :fit_window, :Qbb, :Qbb_sigma))}, background::Real = 1e-3)
This function calculates the background strength for a given background index [counts / (kg year)], exposure, fit_window
The background strength is the expected number of background in the fit window
"""
function lambda_bck(v_input::NamedTuple{((:exposure, :efficiency, :fit_window, :Qbb, :Qbb_bias, :Qbb_sigma, :sig_exp))}, background::Union{Real, Unitful.AbstractQuantity} = 1e-3)
    return mvalue(background * v_input.exposure * v_input.fit_window)
end

"""
    peakshape_0nbb(v_input::NamedTuple{((:exposure, :efficiency, :fit_window, :Qbb, :Qbb_sigma))}, signal_strength::Real, background::Real, energy::Real)
This function calculates the peakshape for a given signal_strength (1/0nbb halflife [1e-26 yr]), background, and energy.
The peakshape is the probability density function of the 0nbb signal and background.
"""
function peakshape_0nbb(v_input::NamedTuple{((:exposure, :efficiency, :fit_window, :Qbb, :Qbb_bias, :Qbb_sigma, :sig_exp))}, signal_strength::Union{Real, Unitful.AbstractQuantity}, background::Union{Real, Unitful.AbstractQuantity}, energy::Union{Real, Unitful.AbstractQuantity})
    N_signal = lambda_sig(v_input, signal_strength)
    N_background = lambda_bck(v_input, background)
    return 1 / (N_signal + N_background) * (N_signal * pdf(Normal(ustrip(mvalue(v_input.Qbb)), ustrip(mvalue(v_input.Qbb_sigma + v_input.Qbb_bias))), ustrip(energy)) +  N_background   / ustrip(v_input.fit_window)) 
end

"""
    likelihood_partition(signal_strength::Real , background::Real, v_input::NamedTuple{((:exposure, :efficiency, :fit_window, :Qbb, :Qbb_sigma, :sig_exp))}, Energy::Vector{<:Real} = 2039.0 .+ [-1.0,0.0, 1.0])

This functions calculate the likelihood of the 0nbb signal and background for a given signal_strength, background, v_input, and energy for 1 partition.
"""
function likelihood_partition(signal_strength::Real , background::Real, v_input::NamedTuple{((:exposure, :efficiency, :fit_window, :Qbb, :Qbb_bias, :Qbb_sigma, :sig_exp))}, Energy::Vector{<:Real} = 2039.0 .+ [-1.0,0.0, 1.0])
    if lambda_sig(v_input, signal_strength) < 0 || lambda_bck(v_input, background) < 0 # force positive signal and background
        return 0
    else 
        poisson_term = pdf(Poisson(lambda_sig(v_input, signal_strength) + lambda_bck(v_input, background)), length(Energy))
        likelihood = prod(peakshape_0nbb.(Ref(v_input), signal_strength, background, Energy))
        likelihood_ext = poisson_term * likelihood
        return likelihood_ext
    end
end