using LegendDataManagement
using PropDicts
using Measurements: value as mvalue
using Unitful

"""
        stripunits_nt(nt::NamedTuple)
util function that strips the units of every field in a NamedTuple
"""
function stripunits_nt(nt::NamedTuple)
        return NamedTuple{keys(nt)}(([(ustrip(getproperty(nt, k))) for k in keys(nt)]))   
end

"""
        stripunits_uncert(nt::NamedTuple)
util function that strips the uncertainties of every field in a NamedTuple
"""
function stripuncert_nt(nt::NamedTuple)
        return NamedTuple{keys(nt)}(([(mvalue(getproperty(nt, k))) for k in keys(nt)]))   
end

"""
        fit_input(data::LegendData, detector::DetectorId, partition::DataPartition)
get the input for the fit function from the metadata
        returns a NamedTuple with the following fields:
        - exposure: exposure (kg.yr)
        - efficiency: efficiency based on quality cuts and psd
        - fit_window: energy window for the fit
        - Qbb: Q-value (keV) of the 0nbb decay
        - Qbb_bias: energy bias (keV) of the energy scale
        - Qbb_sigma: energy resolution (keV) of the energy scale
        - sig_exp: signal strength multiplication factor  (reminder: for practial reasons, the half-life fit parameter is expessed as 1/T*sig_exp in likelihood)
"""
function fit_input(data::LegendData, detector::DetectorId, partition::DataPartition; stripunits = true)
        # load  metadata
        metadata_hardware = data.metadata.hardware.detectors.germanium.diodes[detector].production
        metadata = data.metadata.datasets.ovbb_partitions_pars[detector]
        if !haskey(metadata, Symbol(replace(string(partition), "part" => "part00")))
                @warn("No metadata.datasets for this partition ($partition) for this detector ($detector)")
                return nothing
        end
        metadata_part = merge(metadata.default, get(metadata, Symbol(replace(string(partition), "part" => "part00")), PropDict()))

        # exposure calculation
        mass = metadata_hardware.mass_in_g ./ 1000.0 .* u"kg"
        enrichment = measurement(metadata_hardware.enrichment.val, metadata_hardware.enrichment.unc)
        active_volume = measurement(metadata_part.ovbb_acceptance.active_volume.val, metadata_part.ovbb_acceptance.active_volume.unc)
        containment = measurement(metadata_part.ovbb_acceptance.containment.val, metadata_part.ovbb_acceptance.containment.unc)
        livetime = metadata_part.livetime_in_s / (60*60*24*365) .* u"yr"
        exposure = mass * enrichment * active_volume * containment * livetime

        # efficiency calculation: quality cuts * psd 
        eff_quality = measurement(metadata_part.ovbb_acceptance.quality.val, metadata_part.ovbb_acceptance.quality.unc)
        eff_psd = measurement(metadata_part.ovbb_acceptance.psd.val, metadata_part.ovbb_acceptance.psd.unc)
        efficiency = eff_quality * eff_psd

        # energy-scale parameters
        Qbb_sigma = measurement(metadata_part.fwhm_in_keV.val, metadata_part.fwhm_in_keV.unc) / 2.355 .*u"keV"
        Qbb_bias  = measurement(metadata_part.energy_bias_in_keV.val, metadata_part.energy_bias_in_keV.unc)  .*u"keV"
        Qbb = 2039.0 .* u"keV" # tbd
        fit_window = 240 .* u"keV" #tbd

        v_input = (exposure = exposure, 
                efficiency = efficiency, 
                fit_window = fit_window, 
                Qbb = Qbb, 
                Qbb_bias = Qbb_bias,
                Qbb_sigma = Qbb_sigma, 
                sig_exp = 1.0e-26) 

        # summarize for fit input
        if stripunits == true 
                return stripunits_nt(v_input)
        else
                return v_input
        end
end


# # test the function
# data = LegendData(:l200)
# dets_ged = DetectorId.(collect(keys(data.metadata.datasets.ovbb_partitions_pars)))
# v_input = fit_input(data, dets_ged[1], DataPartition(5))
 
