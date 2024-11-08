using LegendDataManagement
using LegendSpecFits
using LegendDataManagement.LDMUtils
using LegendHDF5IO
using TypedTables
using PropDicts
using Plots
using Measures

# data selection 
runsel = (DataPeriod(3), DataRun(0), :cal)
e_type = :e_cusp_ctc 
ChannelIdx = 1

# set up LegendDataManagement
l200 = LegendData(:l200)
chinfo = channelinfo(l200, runsel; system=:geds, only_processable=true)
(period, run, category) = runsel
ch = chinfo.channel[ChannelIdx]
det = chinfo.detector[ChannelIdx]
filekey = start_filekey(l200, runsel) # first filekey in the run

# load config 
energy_config = dataprod_config(l200).energy(filekey)
energy_config_ch = merge(energy_config.default, get(energy_config, det, PropDict()))

# load data table 
data_ch_after_qc = read_ldata(l200, :jlhit, category, period, run, chinfo.channel[ChannelIdx]).dataQC

# energy in ADC with or without charge trapping correction 
if contains(String(e_type), "ctc")
    # load charge trapping correction and apply 
    e_type_name = Symbol(split(string(e_type), "_ctc")[1])
    e_uncal_func = get_values(l200.par.rpars.ctc[period, run, det])[e_type_name].func
    e_uncal = ljl_propfunc(e_uncal_func).(data_ch_after_qc)
else
    # energy in ADC 
    e_uncal = getproperty(data_ch_after_qc, e_type)
end

# simple claibration 
quantile_perc = if energy_config_ch.quantile_perc isa String parse(Float64, energy_config_ch.quantile_perc) else energy_config_ch.quantile_perc end
th228_names = Symbol.(energy_config_ch.th228_names)
th228_lines_dict = Dict(th228_names .=> energy_config_ch.th228_lines)
result_simple, report_simple = simple_calibration(e_uncal, energy_config_ch.th228_lines, energy_config_ch.left_window_sizes, energy_config_ch.right_window_sizes,; calib_type=:th228, n_bins=energy_config_ch.n_bins, quantile_perc=quantile_perc, binning_peak_window=energy_config_ch.binning_peak_window)
plot(report_simple)

# actual fit 
fit_functions = [fill(:gamma_def, 3)..., :gamma_tails_bckFlat, :gamma_def, :gamma_tails_bckFlat, :gamma_def]

result_fit, report_fit = fit_peaks(result_simple.peakhists, result_simple.peakstats, th228_names; 
                                    e_unit=result_simple.unit, calib_type=:th228, fit_func = fit_functions)
p = plot(broadcast(k -> plot(report_fit[k], left_margin=20mm, top_margin=-5mm, bottom_margin=-2mm, title=string(k), ms=2), keys(report_fit))..., layout=(length(report_fit), 1), size=(1000,710*length(report_fit)) , thickness_scaling=1.8, titlefontsize = 10, legendfontsize = 8, yguidefontsize = 9, xguidefontsize=11)
plot!(p, plot_title=get_plottitle(filekey, det, "Peak Fits"; additiional_type=string(e_type)), plot_titlelocation=(0.5,0.2), plot_titlefontsize = 12)


