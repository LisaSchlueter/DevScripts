import PhysicalConstants.CODATA2018: N_A, m_u # avogadro number, atomic mass constant
using Distributions 
using Interpolations, LinearAlgebra
# using Plots
using Measurements: value as mvalue
include("0nbb_likelihood.jl")
include("0nbb_get_vinput.jl")

function simspec(v_input::NamedTuple{((:exposure, :efficiency, :fit_window, :Qbb, :Qbb_bias, :Qbb_sigma, :sig_exp))}, signal_strength::Union{Real, Unitful.AbstractQuantity}, background::Union{Real, Unitful.AbstractQuantity}; plotFlag::Bool = false)

    # number of expected events in fit window
    Nevents_lambda = lambda_sig(v_input, signal_strength) + lambda_bck(v_input, background)
    Nevents = rand(Poisson(mvalue(Nevents_lambda)))

    # randomize energy values (spectrum): sample from peakshape 
    e_unit = unit(v_input.Qbb)
    e_sample = collect(v_input.Qbb + mvalue(v_input.Qbb_bias) - v_input.fit_window/2:0.1*e_unit:v_input.Qbb + mvalue(v_input.Qbb_bias)  + v_input.fit_window/2)
    signal_pdf = peakshape_0nbb.(Ref(v_input), Ref(signal_strength), Ref(background), e_sample)
    signal_cdf = cumsum(signal_pdf)./maximum(cumsum(signal_pdf))
    interp_cdf_inv = linear_interpolation(signal_cdf, e_sample) # inverse cdf
    bandwidth = maximum(signal_cdf)-minimum(signal_cdf)
    rand_i = minimum(signal_cdf).+bandwidth.*rand(Nevents); 
    events_mc = interp_cdf_inv.(rand_i)

    if plotFlag
        plt = stephist(events_mc, nbins = ifelse(Nevents < 100, Nevents, 100), fill = true, label = "Simulated spectrum")
        xlims!(ustrip(v_input.Qbb + mvalue(v_input.Qbb_bias) - v_input.fit_window/2), ustrip(v_input.Qbb + mvalue(v_input.Qbb_bias) + v_input.fit_window/2))
        display(plt)
    end 
    return events_mc
end

# # test the function
# data = LegendData(:l200)
# dets_ged = DetectorId.(collect(keys(data.metadata.datasets.ovbb_partitions_pars)))
# v_input = fit_input(data, dets_ged[1], DataPartition(5))
# v_input = merge(v_input, (exposure = v_input.exposure * 1000,)) #boost exposure
# v_test = (signal_strength = 2, background = 5.2e-4)
# energies = simspec(v_input, v_test...; plotFlag = true) 

