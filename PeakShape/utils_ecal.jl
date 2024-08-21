
using JLD2
using LegendDataManagement
using LegendDataManagement.LDMUtils
using LegendSpecFits
using Unitful
using Measurements
using Measurements: value as mvalue, uncertainty as muncert
using TypedTables
using StatsBase
using Dates

function get_chinfo(l200::LegendData; periods::Vector{<:Int} = [3, 4])
    chinfos = [Table(channelinfo(l200, (period, run, :cal) ; system=:geds, only_processable=true)) for period in DataPeriod.(periods) for run in search_disk(DataRun, l200.tier[:raw, :cal, DataPeriod(period)]) ]
    dets_ged = [chinfo.detector for chinfo in chinfos]
    dets_ged = unique(vcat(dets_ged...))   
    ch_ged = [chinfo.channel for chinfo in chinfos]
    ch_ged = unique(vcat(ch_ged...)) 
    if !all([all(chinfos[1].detector  == chinfo.detector) for chinfo in chinfos])
        @info "detectors are not the same for all period - use union"
    end
   
    chinfo = chinfos[1]
    filekey = start_filekey(l200, (DataPeriod(periods[1]), search_disk(DataRun,l200.tier[:raw, :cal, DataPeriod(periods[1])])[1] , :cal)) 
    ecal_config = l200.metadata.jldataprod.config.energy(filekey).default
    return chinfo, ecal_config, dets_ged, ch_ged 
end

""" 
    get_peakfitpars_run(pd_rpars::Vector{PropDict})
load peak fit parameter from rpars for a given run and all channels 
"""
function get_rpars_peakfit(l200::LegendData, periods::Vector{<:Int}; reload::Bool = true, e_type::Symbol = :e_cusp_ctc , cal_type::Symbol = :ecal, centroid = false) 
    if e_type == :e_cusp_ctc
        path_data = "$(@__DIR__)/EcalData/"
      
        file_data = path_data * "$(cal_type)_peakfit_rpars_$(join(string.(periods))).jld2"

        if isfile(file_data) && reload == false 
            @info "load peak fit parameter from file $file_data..."
            file = jldopen(file_data)
            return file["FitPars"], file["MetaData"]
        end
    end
    @info "load peak fit parameter from rpars for period: $periods"
    chinfo, ecal_config, dets_ged, ch_ged  = get_chinfo(l200; periods = periods)
    th228_names = Symbol.(ecal_config.th228_names)
    th228_literature = ecal_config.th228_lines
    npeaks     = length(th228_names)
  
    # load all ProbDicts (for all period-run combination in selected partition). speed up load probdict for partition 
    pd_rpars = [l200.par.rpars[cal_type, period, run] for period in DataPeriod.(periods) for run in search_disk(DataRun,l200.tier[:raw, :cal, DataPeriod(period)])] 
    nruns = length(pd_rpars)
    period_runs = [(period, run)   for period in DataPeriod.(periods) for run in search_disk(DataRun,l200.tier[:raw, :cal, DataPeriod(period)])]
    period_runs_lbl = ["p$(parse(Int,replace(string(period),"p" => ""))), r$(parse(Int,replace(string(run),"r" => "")))" for (period, run) in period_runs]
    
    # load fit parameter 
    skew_frac = ones(length(dets_ged), nruns, npeaks) .* NaN  .± NaN 
    skew_width = ones(length(dets_ged), nruns, npeaks) .* NaN  .± NaN   
    skew_frac_highE = ones(length(dets_ged), nruns, npeaks) .* NaN  .± NaN 
    skew_width_highE = ones(length(dets_ged), nruns, npeaks) .* NaN  .± NaN   
    µ  = ones(length(dets_ged), nruns, npeaks) .* NaN * u"keV"  .± NaN  * u"keV"
    µ_cen  = ones(length(dets_ged), nruns, 7) .* NaN * u"keV"  .± NaN  * u"keV"
    σ  = ones(length(dets_ged), nruns, npeaks) .* NaN * u"keV"  .± NaN  * u"keV"
    fwhm  = ones(length(dets_ged), nruns, npeaks) .* NaN * u"keV"  .± NaN  * u"keV"
    fwhm_qbb  = ones(length(dets_ged), nruns) .* NaN * u"keV"  .± NaN  * u"keV"
    background = ones(length(dets_ged), nruns, npeaks) .* NaN  .± NaN 
    background_slope = ones(length(dets_ged), nruns, npeaks) .* NaN  .± NaN 
    n = ones(length(dets_ged), nruns, npeaks) .* NaN  .± NaN 
    step_amplitude = ones(length(dets_ged), nruns, npeaks) .* NaN  .± NaN 
    pvalue = ones(length(dets_ged), nruns, npeaks) .* NaN
    chi2 = ones(length(dets_ged), nruns, npeaks) .* NaN
    pvalue_ls = ones(length(dets_ged), nruns, npeaks) .* NaN
    chi2_ls = ones(length(dets_ged), nruns, npeaks) .* NaN
    dof = ones(length(dets_ged), nruns, npeaks) .* NaN
    for i = 1:nruns
        for (d, det) in enumerate(dets_ged)  
            loadLog = true
            if isa(pd_rpars[i], AbstractDict)
                if !haskey(pd_rpars[i], det) 
                    loadLog = false 
                elseif haskey(pd_rpars[i],det) 
                    if !haskey(pd_rpars[i][det],e_type)
                        loadLog = false 
                    end
                end
            else 
                loadLog = false 
            end
            if loadLog == true
                fwhm_qbb[d,i] = pd_rpars[i][det][e_type].fwhm.qbb
                if centroid ==true
                    µ_cen_ADC = pd_rpars[i][det][e_type].cal_cen.µ
                    cal_func_cen = ljl_propfunc(pd_rpars[i][det][e_type].cal_cen.func)
                    Tbl_cen = Table(e_cusp = µ_cen_ADC, qdrift = zeros(length(µ_cen_ADC)))
                    µ_cen_keV = collect(cal_func_cen.(Tbl_cen))
                    for (p, peakval) in enumerate(pd_rpars[i][det][e_type].cal_cen.peaks)
                        pIdx = findfirst(peakval .== th228_literature)
                        µ_cen[d,i,pIdx] = µ_cen_keV[p]
                    end
                end
                for (p, pname) in enumerate(th228_names)
                    if haskey(pd_rpars[i][det][e_type].fit, pname)
                        skew_frac[d,i,p] = pd_rpars[i][det][e_type].fit[pname].skew_fraction
                        skew_width[d,i,p] = pd_rpars[i][det][e_type].fit[pname].skew_width   
                        µ_ADC  = pd_rpars[i][det][e_type].fit[pname].µ ./ pd_rpars[i][det][e_type].m_cal_simple
                        σ_ADC = pd_rpars[i][det][e_type].fit[pname].σ ./ pd_rpars[i][det][e_type].m_cal_simple
                        fwhm_ADC = pd_rpars[i][det][e_type].fit[pname].fwhm ./ pd_rpars[i][det][e_type].m_cal_simple
                        Tbl = Table(e_cusp = [µ_ADC, σ_ADC, fwhm_ADC], qdrift = [0,0,0])
                        cal_func = ljl_propfunc(pd_rpars[i][det][e_type].cal.func)
                        µ[d,i,p] = collect(cal_func.(Tbl))[1] # calibrated peak position for period given run and detector 
                        σ[d,i,p] = collect(cal_func.(Tbl))[2] .- pd_rpars[i][det][e_type].cal.par[1] # calibrated (w/o offset)
                        fwhm[d,i,p] = collect(cal_func.(Tbl))[3] .- pd_rpars[i][det][e_type].cal.par[1] # calibrated (w/o offset) 
                        background[d,i,p] = pd_rpars[i][det][e_type].fit[pname].background
                        n[d,i,p] = pd_rpars[i][det][e_type].fit[pname].n
                        step_amplitude[d,i,p] = pd_rpars[i][det][e_type].fit[pname].step_amplitude
                        pvalue[d,i,p] = pd_rpars[i][det][e_type].fit[pname].gof.pvalue 
                        chi2[d,i,p] = pd_rpars[i][det][e_type].fit[pname].gof.chi2 
                        pvalue_ls[d,i,p] = pd_rpars[i][det][e_type].fit[pname].gof.pvalue_ls
                        chi2_ls[d,i,p] = pd_rpars[i][det][e_type].fit[pname].gof.chi2_ls 
                        dof[d,i,p] = pd_rpars[i][det][e_type].fit[pname].gof.dof
                        if haskey(pd_rpars[i][det][e_type].fit[pname], :background_slope)
                            background_slope[d,i,p] = pd_rpars[i][det][e_type].fit[pname].background_slope
                        else
                            background_slope[d,i,p] = NaN ± NaN
                        end
                        if haskey(pd_rpars[i][det][e_type].fit[pname], :skew_fraction_highE)
                            skew_frac_highE[d,i,p] = pd_rpars[i][det][e_type].fit[pname].skew_fraction_highE
                            skew_width_highE[d,i,p] = pd_rpars[i][det][e_type].fit[pname].skew_width_highE 
                        else
                            skew_frac_highE[d,i,p] = NaN ± NaN
                            skew_width_highE[d,i,p] = NaN ± NaN
                        end
                    else
                        skew_frac[d,i,p] = NaN ± NaN
                        skew_width[d,i,p] = NaN ± NaN
                        skew_frac_highE[d,i,p] = NaN ± NaN
                        skew_width_highE[d,i,p] = NaN ± NaN
                        µ[d,i,p] = NaN* u"keV"  .± NaN  * u"keV"
                        σ[d,i,p] = NaN* u"keV"  .± NaN  * u"keV"
                        fwhm[d,i,p] = NaN* u"keV"  .± NaN  * u"keV"
                        background[d,i,p] = NaN ± NaN
                        n[d,i,p] = NaN ± NaN
                        step_amplitude[d,i,p] = NaN ± NaN
                        pvalue[d,i,p] = NaN
                        chi2[d,i,p] = NaN
                        pvalue_ls[d,i,p] = NaN
                        chi2_ls[d,i,p] = NaN
                        dof[d,i,p] = NaN
                        background_slope[d,i,p] = NaN ± NaN
                    end
                end
            else
                skew_frac[d,i,:] = fill(NaN ± NaN, npeaks)
                skew_width[d,i,:] = fill(NaN ± NaN, npeaks)
                skew_frac_highE[d,i,:] = fill(NaN ± NaN, npeaks)
                skew_width_highE[d,i,:] = fill(NaN ± NaN, npeaks)
                µ[d,i,:] = fill(NaN* u"keV"  .± NaN  * u"keV", npeaks)
                σ[d,i,:] = fill(NaN* u"keV"  .± NaN  * u"keV", npeaks)
                fwhm[d,i,:] = fill(NaN* u"keV"  .± NaN  * u"keV", npeaks)
                fwhm_qbb[d,i] = NaN* u"keV"  .± NaN  * u"keV"
                background[d,i,:] = fill(NaN ± NaN, npeaks)
                n[d,i,:] = fill(NaN ± NaN,npeaks)
                step_amplitude[d,i,:] = fill(NaN ± NaN, npeaks)
                pvalue[d,i,:] = fill(NaN, npeaks)
                chi2[d,i,:] = fill(NaN, npeaks)
                pvalue_ls[d,i,:] = fill(NaN, npeaks)
                chi2_ls[d,i,:] = fill(NaN, npeaks)
                dof[d,i,:] = fill(NaN, npeaks)
                background_slope[d,i,:] = fill(NaN ± NaN, npeaks)
                µ_cen[d,i,:] = fill(NaN* u"keV"  .± NaN  * u"keV", npeaks)
                @info "$(period_runs[i]), det $det ($d) missing"
            end
        end
    end

    residual = permutedims(permutedims(µ,(3,1,2)) .- th228_literature,(2,3,1))

    FitPars = (µ = µ, 
                σ = σ, 
                fwhm = fwhm, 
                background = background, 
                n = n, 
                step_amplitude = step_amplitude, 
                skew_frac = skew_frac,
                skew_width = skew_width,
                skew_frac_highE = skew_frac_highE,
                skew_width_highE = skew_width_highE,
                pvalue = pvalue,
                chi2 = chi2,
                pvalue_ls = pvalue_ls,
                chi2_ls = chi2_ls,
                dof = dof,
                residual = residual,
                fwhm_qbb = fwhm_qbb,
                µ_cen = µ_cen,
                background_slope = background_slope)

    MetaData = (chinfo = chinfo, 
              ecal_config = ecal_config, 
              dets_ged = dets_ged, 
              ch_ged = ch_ged,
              period_runs = period_runs,
              period_runs_lbl = period_runs_lbl,
              cal_type = cal_type,)

    if e_type == :e_cusp_ctc
        if !ispath(path_data)
            mkdir("$path_data")
        end
        jldsave(file_data, FitPars = FitPars, MetaData = MetaData)
        @info "saved ecal fit pars to $file_data"
    end
    return FitPars, MetaData
end


function get_rpars_calcurve(l200::LegendData, period::Int; reload::Bool = false, e_type::Symbol = :e_cusp_ctc )
    path_data = "$(@__DIR__)/EcalData/"
    file_data = path_data * "Ecal_calcurve_rpars_$(join(string.(period))).jld2"

    if isfile(file_data) && reload == false && e_type == :e_cusp_ctc
        @info "load calibration curve parameter from file $file_data..."
        file = jldopen(file_data)
        return file["CalPars"], file["MetaData"]
    end

    @info "load calibration curve parameter from rpars for period: $period"
    chinfo, ecal_config, dets_ged, ch_ged  = get_chinfo(l200; periods = [period])
    th228_names = Symbol.(ecal_config.th228_names)
    th228_literature = ecal_config.th228_lines .* u"keV"
    npeaks     = length(th228_names)
  
    # load all ProbDicts (for all period-run combination in selected period). speed up load probdict 
    pd_rpars = [l200.par.rpars.ecal[DataPeriod(period), run] for run in search_disk(DataRun,l200.tier[:raw, :cal, DataPeriod(period)])] 
    nruns = length(pd_rpars)
    period_runs = [(DataPeriod(period), run)  for run in search_disk(DataRun,l200.tier[:raw, :cal, DataPeriod(period)])]
    period_runs_lbl = ["p$(parse(Int,replace(string(period),"p" => ""))), r$(parse(Int,replace(string(run),"r" => "")))" for (period, run) in period_runs]

    # load peaks 
    # good_det = ones(length(dets_ged)) .* true

    # re-shuffle results for a specific detector
    pvalues = zeros((length(dets_ged), nruns)) .* NaN   
    µ_ADC   = zeros((length(dets_ged), nruns, npeaks))  .* NaN  .± NaN 
    peaks_literature_keV   = zeros((length(dets_ged), nruns, npeaks))  .* NaN   .*u"keV"
    residuals_norm  = zeros((length(dets_ged), nruns, npeaks)).* NaN  .± NaN 
    y_fit  = zeros((length(dets_ged), nruns, npeaks)).* NaN .* u"keV" .± NaN .* u"keV"
    calcurve_x =  zeros((length(dets_ged), nruns, 100)) .* NaN 
    calcurve_y =  zeros((length(dets_ged), nruns, 100)) .* NaN * u"keV"

    for i = 1:nruns
        @info "run $i"
        for (d, det) in enumerate(dets_ged)
            loadLog = true
            if isa(pd_rpars[i], AbstractDict)
                if !haskey(pd_rpars[i],det) 
                    loadLog = false 
                    @info "det $det ($d) missing"
                elseif haskey(pd_rpars[i],det) 
                    if !haskey(pd_rpars[i][det],e_type)
                        loadLog = false 
                        @info "$(e_type) of det $det ($d) missing"
                    end      
                end
            else
                loadLog = false
            end
            if loadLog == true
                #@info "load det $det ($d)"
                pvalues[d,i] = pd_rpars[i][dets_ged[d]][e_type].cal.gof.pvalue
                peaks_tmp =  pd_rpars[i][dets_ged[d]][e_type].cal.peaks
                if length(peaks_tmp) !== npeaks
                    peak_idx = reduce(vcat, [findall(peaks_tmp[peak] .== th228_literature) for peak in eachindex(peaks_tmp)])
                    peaks_literature_keV[d,i,peak_idx] = peaks_tmp 
                    µ_ADC[d,i,peak_idx] =  pd_rpars[i][dets_ged[d]][e_type].cal.µ
                    residuals_norm[d, i, peak_idx] = pd_rpars[i][dets_ged[d]][e_type].cal.gof.residuals_norm
                else
                    peaks_literature_keV[d,i,:] = peaks_tmp 
                    µ_ADC[d,i,:] =  pd_rpars[i][dets_ged[d]][e_type].cal.µ
                    residuals_norm[d,i,:] = pd_rpars[i][dets_ged[d]][e_type].cal.gof.residuals_norm
                end
                calfunc_tmp = ljl_propfunc(pd_rpars[i][dets_ged[d]][e_type].cal.func)
                calcurve_x[d,i,:] = range(mvalue(minimum(filter(!isnan,µ_ADC[d,i,:])))-100, stop = mvalue(maximum(filter(!isnan,µ_ADC[d,i,:])))+100, length = 100)
                Tbl = Table(e_cusp = calcurve_x[d,i,:], qdrift = fill(0, length(calcurve_x[d,i,:])))
                calcurve_y[d,i,:] =   calfunc_tmp.(Tbl)
                Tbl_y = Table(e_cusp = µ_ADC[d,i,:], qdrift = fill(0, length(µ_ADC[d,i,:])))
                y_fit[d,i,:] = calfunc_tmp.(Tbl_y) 
            else
                @info "skip run $i det $det ($d)"
                µ_ADC[d,i,:] = fill(NaN, npeaks)
                residuals_norm[d,i,:] = fill(NaN, npeaks)
                pvalues[d,i] = NaN
                calcurve_x[d,i,:] = fill(NaN, 100)
                peaks_literature_keV[d,i,:] = fill(NaN, npeaks) .* u"keV"
                calcurve_y[d,i,:] = fill(NaN, 100) .* u"keV"
                y_fit[d,i,:] = fill(NaN, npeaks) .* u"keV" 
                good_det[d] = false
            end
        end
    end
    CalPars = (µ_ADC = µ_ADC, 
              peaks_literature_keV = peaks_literature_keV,
              th228_literature = th228_literature,
              residuals_norm = residuals_norm,
              pvalue = pvalues,
              y_fit = y_fit,
              calcurve_x = calcurve_x,
              calcurve_y = calcurve_y)
    MetaData = (chinfo = chinfo, 
            ecal_config = ecal_config, 
            dets_ged = dets_ged, 
            ch_ged = ch_ged,
            period_runs = period_runs,
            period_runs_lbl = period_runs_lbl,)

    if !ispath(path_data)
        mkdir("$path_data")
    end
    jldsave(file_data, CalPars = CalPars, MetaData = MetaData)
    @info "saved ecal cal pars to $file_data"
    return CalPars, MetaData
end


function get_energy_ctc(data::LegendData, runsel::Tuple{DataPeriod, DataRun, Symbol}, detector::DetectorId; e_type::Symbol=:e_cusp)
    channel = detector2channel(data, runsel, detector)
    (period, run, category) = runsel
    data_ch = lh5open(data.tier[:jlhitch, category, period, run, channel])[channel, :dataQC][:]
    e_type_name = Symbol(split(string(e_type), "_ctc")[1])
    if endswith(string(e_type), "_ctc")
        @info "Apply CT correction for $e_type"
        e_uncal_func = pars_ctc[detector][e_type_name].func
        e_uncal = ljl_propfunc(e_uncal_func).(data_ch)
    else
        e_uncal = getproperty(data_ch, e_type_name)
    end
    return e_uncal
end


function simple_cal(data::LegendData, runsel::Tuple{DataPeriod, DataRun, Symbol}, detector::DetectorId; e_type::Symbol=:e_cusp)
    e_ADC = get_energy_ctc(data, runsel, detector; e_type = e_type)
    # get config for calibration fits
    energy_config = dataprod_config(data).energy(start_filekey(data, runsel))
    energy_config_ch = merge(energy_config.default, get(energy_config, detector, PropDict()))

    quantile_perc = if energy_config_ch.quantile_perc isa String parse(Float64, energy_config_ch.quantile_perc) else energy_config_ch.quantile_perc end
    th228_names = Symbol.(energy_config_ch.th228_names)
    th228_lines_dict = Dict(th228_names .=> energy_config_ch.th228_lines)

    # simple calibration 
    result_simple, report_simple = simple_calibration(e_ADC, energy_config_ch.th228_lines, energy_config_ch.left_window_sizes, energy_config_ch.right_window_sizes,; calib_type=:th228, n_bins=energy_config_ch.n_bins, quantile_perc=quantile_perc, binning_peak_window=energy_config_ch.binning_peak_window)
    @info "simple calibration done"
    return result_simple, report_simple, th228_names, th228_lines_dict
end
