
using LegendDataManagement
using LegendDataManagement.LDMUtils
using LegendSpecFits
using LegendHDF5IO
using PropDicts
using Plots
using Unitful, Measurements, Measures
using StatsBase

# processing config
reprocess = true 

# failed periods
FailDict = PropDict()
FailDict[1] = (period = DataPeriod(3), run = DataRun(0), detector = DetectorId(:V05612B), e_type = :e_cusp)
FailDict[2] = (period = DataPeriod(4), run = DataRun(2), detector = DetectorId(:V08682A), e_type = :e_535)
FailDict[3] = (period = DataPeriod(4), run = DataRun(2), detector = DetectorId(:V05267A), e_type = :e_trap)
FailDict[4] = (period = DataPeriod(4), run = DataRun(4), detector = DetectorId(:P00573B), e_type = :e_cusp)
# FailDict[5] = (period = DataPeriod(4), run = DataRun(2), detector = DetectorId(:V01389A), e_type = :e_trap) ?
# select data 
idx = 1
period = FailDict[idx].period
run = FailDict[idx].run
detector = FailDict[idx].detector
e_type = FailDict[idx].e_type
category = :cal 
l200 = LegendData(:l200)

##### start processor 
begin # preparation: load configs etc. 
    @info "Energy calibration for period $period and run $run"

    filekey = start_filekey(l200, (period, run, :cal))
    @info "Found filekey $filekey"

    chinfo = channelinfo(l200, filekey; system=:geds, only_processable=true)
    @info "Loaded channel info with $(length(chinfo)) channels"

    energy_config = dataprod_config(l200).energy(filekey)
    @debug "Loaded energy config: $(energy_config)"

    pars_ctc = get_values(l200.par.rpars.ctc[period, run])
    @debug "Loaded CTC parameters"

    @debug "Create pars db"
    # mkpath(joinpath(data_path(l200.par.rpars.ecal), string(period)))
    pars_db = PropDict(l200.par.rpars.ecal[period, run])

    pars_db = ifelse(reprocess, PropDict(), pars_db)
    if reprocess @info "Reprocess all channels" end

    # create log line Tuple
    log_nt = NamedTuple{(:Channel, :Detector, :Status, Symbol("Filter Type"), Symbol("FWHM Qbb"), Symbol("FWHM FEP"), Symbol("Cal. Constant"), :Error)}


    # function ch_energy_calibration(chinfo_ch::NamedTuple)
    chIdx = findfirst(Symbol.(chinfo.detector) .== Symbol(detector))
    chinfo_ch = chinfo[chIdx]
    ch  = chinfo_ch.channel
    det = chinfo_ch.detector
end 

# load data after quality cuts 
begin 
    @debug "Processing channel $ch ($det)"
    hitchfilename = l200.tier[:jlhit, filekey, ch]
    # load data file
    if !isfile(hitchfilename)
        @error "Hit file $hitchfilename not found"
        throw(ErrorException("Hit file not found"))
    end

    result_dict    = Dict{Symbol, NamedTuple}()
    log_info_dict  = Dict{Symbol, NamedTuple}()
    processed_dict = Dict{Symbol, Bool}()

    energy_config_ch = merge(energy_config.default, get(energy_config, det, PropDict()))

    energy_types = Symbol.(energy_config_ch.energy_types)

    # load data file
    if !isfile(hitchfilename)
        @error "Hit file $hitchfilename not found"
        throw(ErrorException("Hit file not found"))
    end


    # get data
    data_ch_after_qc = nothing
    try
        @debug "Load hit file"
        data_hit = lh5open(hitchfilename, "r");
        data_ch_after_qc = data_hit[ch, :jlhit, :dataQC][:];
        close(data_hit)
    catch e
        @error "Error in loading data for channel $ch: $(truncate_string(string(e)))"
        throw(ErrorException("Error data loader"))
    end

    quantile_perc = if energy_config_ch.quantile_perc isa String parse(Float64, energy_config_ch.quantile_perc) else energy_config_ch.quantile_perc end
    th228_names = Symbol.(energy_config_ch.th228_names)
    th228_lines_dict = Dict(th228_names .=> energy_config_ch.th228_lines)
end 

# apply charge trapping correction, in case energy type has _ctc
begin 
    @debug "Calibrate $e_type"
    # get data
    e_uncal, e_uncal_func = nothing, nothing
    try
        @debug "Get $e_type data"
        # open hit data file
        e_type_name = Symbol(split(string(e_type), "_ctc")[1])
        e_uncal = getproperty(data_ch_after_qc, e_type_name)
        e_uncal_func = "$e_type_name"
        if endswith(string(e_type), "_ctc")
            @debug "Apply CT correction for $e_type"
            e_uncal_func = pars_ctc[det][e_type_name].func
            e_uncal = ljl_propfunc(e_uncal_func).(data_ch_after_qc)
        end
    catch e
        @error "Error in $e_type data extraction for channel $ch: $(truncate_string(string(e)))"
        throw(ErrorException("Error in $e_type data extraction"))
    end
    GC.gc()
end 

# fit peaks for calibration 
begin 
    result_simple, report_simple = nothing, nothing
    try
        @debug "Get $e_type simple calibration"
        result_simple, report_simple = simple_calibration(e_uncal, energy_config_ch.th228_lines, energy_config_ch.left_window_sizes, energy_config_ch.right_window_sizes,; calib_type=:th228, n_bins=energy_config_ch.n_bins, quantile_perc=quantile_perc, binning_peak_window=energy_config_ch.binning_peak_window)
    catch e
        @error "Error in $e_type simple calibration for channel $ch: $(truncate_string(string(e)))"
        throw(ErrorException("Error in $e_type simple calibration"))
    end
    GC.gc()

    # get simple calibration constant
    m_cal_simple = result_simple.c

    result_fit, report_fit = nothing, nothing
    try
        @debug "Fit all $e_type peaks"
        result_fit, report_fit = fit_peaks(result_simple.peakhists, result_simple.peakstats, th228_names; e_unit=result_simple.unit, calib_type=:th228)
    catch e
        @error "Error in $e_type peak fitting for channel $ch: $(truncate_string(string(e)))"
        throw(ErrorException("Error in $e_type peak fitting"))
    end
    GC.gc()

    p = plot(broadcast(k -> plot(report_fit[k], left_margin=20mm, top_margin=-5mm, bottom_margin=-2mm, title=string(k), ms=2), keys(report_fit))..., layout=(length(report_fit), 1), size=(1000,710*length(report_fit)) , thickness_scaling=1.8, titlefontsize = 10, legendfontsize = 8, yguidefontsize = 9, xguidefontsize=11)
    plot!(p, plot_title=get_plottitle(filekey, det, "Peak Fits"; additiional_type=string(e_type)), plot_titlelocation=(0.5,0.2), plot_titlefontsize = 12)
end 

# fit calibration curve
begin 
    @debug "Get $e_type calibration values"

    result_calib, report_calib = nothing, nothing
    try
        μ_fit =  [result_fit[p].μ for p in th228_names if !(p in Symbol.(energy_config_ch.cal_fit_excluded_peaks))] ./ m_cal_simple
        pp_fit = [th228_lines_dict[p] for p in th228_names if !(p in Symbol.(energy_config_ch.cal_fit_excluded_peaks))]
        result_calib, report_calib = fit_calibration(energy_config_ch.cal_pol_order, μ_fit, pp_fit; e_expression=e_uncal_func)
        @debug "Found $e_type calibration curve: $(result_calib.func)"
    catch e
        @error "Error in $e_type calibration curve fitting for channel $ch: $(truncate_string(string(e)))"
        throw(ErrorException("Error in $e_type calibration curve fitting"))
    end
    # add not-fitted peaks to plot 
    μ_notfit =  [result_fit[p].μ for p in Symbol.(energy_config_ch.cal_fit_excluded_peaks)] ./ m_cal_simple
    pp_notfit = [th228_lines_dict[p] for p in Symbol.(energy_config_ch.cal_fit_excluded_peaks)]
    p = plot(report_calib, xerrscaling=1, additional_pts=(μ = μ_notfit, peaks = pp_notfit), size = (1100, 800), legend = :topleft )
    plot!(plot_title=get_plottitle(filekey, det, "Calibration Curve"; additiional_type=string(e_type)), plot_titlelocation=(0.5,-0.3), plot_titlefontsize=12)
    display(p)
    f_cal(x) = report_calib.f_fit(x) .* report_calib.e_unit .- first(report_calib.par)

end

# fit fwhm curve 
fwhm_pol_order = 1
begin 
    result_fwhm, report_fwhm = nothing, nothing
    # try
        fwhm_fit = f_cal.([result_fit[p].fwhm for p in th228_names if !(p in Symbol.(energy_config_ch.:fwhm_fit_excluded_peaks))] ./ m_cal_simple)
        pp_fit = [th228_lines_dict[p] for p in th228_names if !(p in Symbol.(energy_config_ch.:fwhm_fit_excluded_peaks))]    
        result_fwhm, report_fwhm = fit_fwhm(pp_fit, fwhm_fit; pol_order=fwhm_pol_order, e_type_cal=Symbol("$(e_type)_cal"), e_expression=e_uncal_func, uncertainty=true)
        @debug "Found $e_type FWHM: $(round(u"keV", result_fwhm.qbb, digits=2))"
    # catch e
        # @error "Error in $e_type FWHM fitting for channel $ch: $(truncate_string(string(e)))"
        # throw(ErrorException("Error in $e_type FWHM fitting"))
    # end
    fwhm_notfit =  f_cal.([result_fit[p].fwhm for p in Symbol.(energy_config_ch.cal_fit_excluded_peaks)] ./ m_cal_simple)
    pp_notfit = [th228_lines_dict[p] for p in Symbol.(energy_config_ch.:fwhm_fit_excluded_peaks)]
    p = plot(report_fwhm, additional_pts=(peaks = pp_notfit, fwhm = fwhm_notfit),  size = (1100, 800), legend = :topleft)
    plot!(plot_title=get_plottitle(filekey, det, "FWHM"; additiional_type=string(e_type)), plot_titlelocation=(0.5,-0.3))
end 


