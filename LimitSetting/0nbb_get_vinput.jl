using LegendDataManagement
using PropDicts
using Measurements
using Unitful

function fit_input(data::LegendData, detector::DetectorId, partition::DataPartition)
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

        # efficiency calculation 
        eff_quality = measurement(metadata_part.ovbb_acceptance.quality.val, metadata_part.ovbb_acceptance.quality.unc)
        eff_psd = measurement(metadata_part.ovbb_acceptance.psd.val, metadata_part.ovbb_acceptance.psd.unc)
        efficiency = eff_quality * eff_psd

        # energy-scale parameters
        Qbb_sigma = measurement(metadata_part.fwhm_in_keV.val, metadata_part.fwhm_in_keV.unc) / 2.355 .*u"keV"
        Qbb_bias  = measurement(metadata_part.energy_bias_in_keV.val, metadata_part.energy_bias_in_keV.unc)  .*u"keV"
        Qbb = 2039.0 .* u"keV" # tbd
        fit_window = 240 .* u"keV" #tbd

        # summarize for fit input
        return (exposure = exposure, 
                efficiency = efficiency, 
                fit_window = fit_window, 
                Qbb = Qbb, 
                Qbb_bias = Qbb_bias,
                Qbb_sigma = Qbb_sigma, 
                sig_exp = 1.0e-26) 
end

data = LegendData(:l200)
dets_ged = DetectorId.(collect(keys(data.metadata.datasets.ovbb_partitions_pars)))

v_input = fit_input(data, dets_ged[1], DataPartition(5))
 
