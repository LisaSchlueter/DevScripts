```
script to understand calibration fits and do some plots along the way
```
# understand simple_calibration and do some plots along the way
using LegendDataManagement
using LegendDataManagement.LDMUtils
using LegendSpecFits
using HDF5, LegendHDF5IO
using TypedTables
using Plots, Measures, Printf

# define data set 
period = 3
run = 0
l200 = LegendData(:l200)

# open data
filekey = start_filekey(l200, (DataPeriod(period), DataRun(run), :cal)) 
chinfo = Table(channelinfo(l200, filekey; system=:geds, only_processable=true))
dets_ged = chinfo.detector
ch_ged = chinfo.channel
chIdx = 1 # select a detector
# load data file
hitchfilename = l200.tier[:jlhitch, filekey, ch_ged[chIdx]]
if !isfile(hitchfilename)
    @error "Hit file $hitchfilename not found"
    throw(ErrorException("Hit file not found"))
else
    data_hit = LHDataStore(hitchfilename, "r");
    data_ch_after_qc = data_hit[ch_ged[chIdx]].dataQC[:];
    close(data_hit)
    @info "load data from $hitchfilename"
end

# load uncalibrated energies and do charge trapping correction
e_type =:e_cusp
pars_ctc = get_values(l200.par.rpars.ctc[DataPeriod(period), DataRun(run)])
e_uncal_func = pars_ctc[dets_ged[chIdx]][e_type].func
e_uncal = ljl_propfunc(e_uncal_func).(data_ch_after_qc)


# fit spectrum: "simple calibration"
energy_config = dataprod_config(l200).energy(filekey).default
th228_lines = energy_config.th228_lines
th228_names = energy_config.th228_names
result_simple, report_simple = simple_calibration(e_uncal, th228_lines, energy_config.left_window_sizes, energy_config.right_window_sizes,; calib_type = :th228, n_bins= energy_config.n_bins)
plot(report_simple)


# fit peak shape: actual calibration, all of them at once 
result_fitall, report_fitall = fit_peaks(result_simple.peakhists, result_simple.peakstats, th228_names; e_unit=result_simple.unit, calib_type=:th228)
p = plot(broadcast(k -> plot(report_fitall[k], left_margin=20mm, top_margin=-5mm, bottom_margin=-2mm, title=string(k), ms=2), keys(report_fit))..., layout=(length(report_fitall), 1), size=(1000,710*length(report_fitall)) , thickness_scaling=1.8, titlefontsize = 10, legendfontsize = 8, yguidefontsize = 9, xguidefontsize=11)
              
# fit peak shape: actual calibration, 1 peak at a time
peakIdx = 1
result_fit, report_fit =  fit_single_peak_th228(result_simple.peakhists[peakIdx], result_simple.peakstats[peakIdx]; uncertainty = true, low_e_tail = true)
plot(report_fit)