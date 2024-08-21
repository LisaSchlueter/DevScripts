
# using PhysicalConstants.CODATA2018
import PhysicalConstants.CODATA2018: N_A, m_u # avogadro number, atomic mass constant
using LegendDataManagement
using Unitful
using Plots 
# # load meta data: input for efficiency, exposure,
# l200 = LegendData(:l200)
# dets_ged = collect(keys(l200.metadata.datasets.ovbb_partitions_pars))

# function exposure(data::LegendData = l200, detector::DetectorId = DetectorId(:V07647B))
#    data.metadata.datasets.ovbb_partitions_pars[detector].default.exposure.val
# end

# l200.metadata.datasets.ovbb_partitions_pars[dets_ged[1]].default.ovbb_acceptance.active_volume.val

# l200.metadata.datasets.ovbb_partitions_pars[dets_ged[2]].default

# l200.metadata.datasets.ovbb_partitions_pars[dets_ged[2]].part0005

v_input = (exposure = 103.7u"kg*yr", efficiency = 0.99, fit_window = 240u"keV", Qbb = 2039u"keV", Qbb_sigma = 2.5u"keV")

# create dummy data set for testing 
halflife_0nbb = 2*1e26u"yr"
background = 5.2e-4u"1/(keV * kg * yr)" 
signal_strength(v_input)
background_strength(v_input, background)

function signal_strength(v_input::NamedTuple, halflife_0nbb::Quantity = 1e27u"yr")
    AtomicMass_Ge76 = 75.921402726 *  m_u  # weight of 1 Ge76 atom in [kg]
    Ge_MolarMass =  AtomicMass_Ge76 * N_A  # weight of 1 mol of Ge76 atom in [kg]
    return 1 / halflife_0nbb * N_A / Ge_MolarMass * v_input.exposure * v_input.efficiency * log(2) 
end

function background_strength(v_input::NamedTuple, background::Quantity = 1e-3)
    return background * v_input.exposure * v_input.fit_window
end


function peakshape_pdf(v_input, halflife_0nbb, background, Eevent::Quantity)
    lambda_sig = signal_strength(v_input, halflife_0nbb)
    lambda_bck = background_strength(v_input, background)
    return 1 / (lambda_sig + lambda_bck) .* (lambda_sig .* pdf(Normal(ustrip(v_input.Qbb), ustrip(v_input.Qbb_sigma)), ustrip(Eevent)) .+  lambda_bck / ustrip(v_input.fit_window))
end

Eevents =  collect(range(v_input.Qbb -v_input.fit_window/2, v_input.Qbb + v_input.fit_window/2, step = 1.0u"keV"))
lambda_sig = signal_strength(v_input, halflife_0nbb)
lambda_bck = background_strength(v_input, background)

peakshape_pdf.(Ref(v_input), Ref(halflife_0nbb), Ref(background), Eevents)
plot(Eevents, peakshape_pdf.(Ref(v_input), Ref(halflife_0nbb), Ref(background), Eevents))

# v_input should be 1 per partition

# likelihood depends on halflife_0nbb, background, (input) and Eevent
function likelihood_partition(v_input::NamedTuple, halflife_0nbb, background, Eevents::Vector{<:Real} = 2039.0 .+ [-1.0,0.0, 1.0])
    lambda_sig = signal_strength(v_input, halflife_0nbb)
    lambda_bck = background_strength(v_input, background)

    Nevent = length(Eevent)
    poisson_term(halflife_0nbb, background) =  pdf(Poisson(lambda_sig + lambda_bck), Nevent)

    # likelihood = prod( 1 / (lambda_sig + lambda_bck) .* (lambda_sig .* pdf.(Normal(v_input.Qbb, v_input.Qbb_sigma), Eevent) .+  lambda_bck / v_input.fit_window))
    likelihood = prod(peakshape_pdf.(Ref(v_input), halflife_0nbb, background, Eevents))
    likelihood_ext = poisson_term * likelihood
    return likelihood_ext
end





