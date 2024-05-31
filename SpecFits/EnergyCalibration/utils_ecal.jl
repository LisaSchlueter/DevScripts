
using JLD2
using LegendDataManagement
using LegendDataManagement.LDMUtils
using Unitful
using Measurements
using Measurements: value as mvalue, uncertainty as muncert
using TypedTables
using StatsBase
function get_histstats(data; nbins::Int = 300, bin_edges::Vector= [])
    if length(size(data)) > 1
        data = reshape(data,:)
    end
    if any(isnan.(data))
        data = filter(isfinite,data)
    end
    if isempty(bin_edges)
        h = fit(Histogram, data, nbins = nbins) 
    else
        h = fit(Histogram, data, edges)
    end 
    counts = h.weights
    edges = h.edges[1]
    width = diff(edges)[1]
    center = edges .+ width/2
    return center, edges, width, counts
end


function get_fwhm(centers::StepRangeLen, counts::Vector)
    max_idx = argmax(counts)
    max_val = counts[max_idx]
    half_max = max_val / 2
    left_idx = findfirst(counts[1:max_idx] .>= half_max)
    right_idx = findfirst(counts[max_idx:end] .<= half_max) + (max_idx-1)
    try
        interp_left = LinearInterpolation(counts[left_idx-1:left_idx+1],centers[left_idx-1:left_idx+1])
        interp_right = LinearInterpolation(counts[right_idx+1:-1:right_idx-1],centers[right_idx+1:-1:right_idx-1])
        fwhm = interp_right(half_max) - interp_left(half_max)
        return fwhm 
    catch
        fwhm = centers[right_idx] - centers[left_idx]
        return fwhm
    end
end

""" 
    get_peakfitpars_run(probdict_part::Vector{PropDict})
load peak fit parameter from rpars for a given run and all channels 
"""
function get_peakfit_rpars_partition(partition::DataPartition; reload::Bool = false, e_type::Symbol = :e_cusp_ctc , cal::Bool = true)
    path_data = "$(@__DIR__)/EcalData/"
    if cal == true
        file_data = path_data * "Ecal_peakfit_rpars_$(partition).jld2"
    elseif cal == false
        file_data = path_data * "Ecal_peakfit_rpars_$(partition)_simplecal.jld2"
    end

    if isfile(file_data) && reload == false && e_type == :e_cusp_ctc
        @info "load peak fit parameter from file $file_data..."
        file = jldopen(file_data)
        return file["FitPars"], file["MetaData"]
    end
    @info "load peak fit parameter from rpars for $partition"
    l200 = LegendData(:l200)
    partinfo = partitioninfo(l200)[partition]
    filekey = start_filekey(l200, (partinfo[1].period, partinfo[1].run, :cal)) 
    chinfo = Table(channelinfo(l200, filekey; system=:geds, only_processable=true))
    dets_ged = chinfo.detector
    dets_type = chinfo.det_type

    # load all ProbDicts (for all period-run combination in selected partition). speed up load probdict for partition 
    probdict_part = [l200.par.rpars.ecal[entr.period, entr.run] for entr in partinfo] 
    nruns = length(probdict_part)
    # load ALL peaks 
    th228_names = Symbol.(keys(probdict_part[1][Symbol(dets_ged[1])][e_type].fit))
    npeaks = length(th228_names)
  
    # load fit parameter 
    skew_frac = ones(length(dets_ged), nruns, npeaks) .* NaN  .± NaN 
    skew_width = ones(length(dets_ged), nruns, npeaks) .* NaN  .± NaN   
    µ  = ones(length(dets_ged), nruns, npeaks) .* NaN * u"keV"  .± NaN  * u"keV"
    µ_simple  = ones(length(dets_ged), nruns, npeaks) .* NaN * u"keV"  .± NaN  * u"keV" # simple calibrated
    σ  = ones(length(dets_ged), nruns, npeaks) .* NaN * u"keV"  .± NaN  * u"keV"
    fwhm  = ones(length(dets_ged), nruns, npeaks) .* NaN * u"keV"  .± NaN  * u"keV"
    background = ones(length(dets_ged), nruns, npeaks) .* NaN  .± NaN 
    n = ones(length(dets_ged), nruns, npeaks) .* NaN  .± NaN 
    step_amplitude = ones(length(dets_ged), nruns, npeaks) .* NaN  .± NaN 
    pvalue = ones(length(dets_ged), nruns, npeaks) .* NaN
    for i = 1:nruns
        for (d, det) in enumerate(dets_ged)  
            loadLog = true
            if !haskey(probdict_part[i],det) 
                loadLog = false 
            elseif haskey(probdict_part[i],det) 
                if !haskey(probdict_part[i][det],e_type)
                    loadLog = false 
                end
            end
            if loadLog == true
                for (p, pname) in enumerate(th228_names)
                    if haskey(probdict_part[i][det][e_type].fit,pname)
                        skew_frac[d,i,p] = probdict_part[i][det][e_type].fit[pname].skew_fraction
                        skew_width[d,i,p] = probdict_part[i][det][e_type].fit[pname].skew_width
                        # µ[d,i,p] = probdict_part[i][det][e_type].fit[pname].µ #they are only simply energy calibrated 
                        # σ[d,i,p] = probdict_part[i][det][e_type].fit[pname].σ
                        # fwhm[d,i,p] = probdict_part[i][det][e_type].fit[pname].fwhm
                        if cal == true
                            µ_ADC  = probdict_part[i][det][e_type].fit[pname].µ ./ probdict_part[i][det][e_type].m_cal_simple
                            σ_ADC = probdict_part[i][det][e_type].fit[pname].σ ./ probdict_part[i][det][e_type].m_cal_simple
                            fwhm_ADC = probdict_part[i][det][e_type].fit[pname].fwhm ./ probdict_part[i][det][e_type].m_cal_simple
                            Tbl = Table(e_cusp = [µ_ADC, σ_ADC, fwhm_ADC], qdrift = [0,0,0])
                            cal_func = ljl_propfunc(probdict_part[i][det][e_type].cal.func)
                            µ[d,i,p] = collect(cal_func.(Tbl))[1] # calibrated peak position for period given run and detector 
                            σ[d,i,p] = collect(cal_func.(Tbl))[2] .- probdict_part[i][det][e_type].cal.par[1] # calibrated (w/o offset)
                            fwhm[d,i,p] = collect(cal_func.(Tbl))[3] .- probdict_part[i][det][e_type].cal.par[1] # calibrated (w/o offset) 
                        elseif cal == false 
                            µ[d,i,p] = probdict_part[i][det][e_type].fit[pname].µ 
                            σ[d,i,p] = probdict_part[i][det][e_type].fit[pname].σ 
                            fwhm[d,i,p] = probdict_part[i][det][e_type].fit[pname].fwhm
                        end
                        background[d,i,p] = probdict_part[i][det][e_type].fit[pname].background
                        n[d,i,p] = probdict_part[i][det][e_type].fit[pname].n
                        step_amplitude[d,i,p] = probdict_part[i][det][e_type].fit[pname].step_amplitude
                        pvalue[d,i,p] = probdict_part[i][det][e_type].fit[pname].gof.pvalue 
                    else
                        skew_frac[d,i,p] = NaN ± NaN
                        skew_width[d,i,p] = NaN ± NaN
                        µ[d,i,p] = NaN* u"keV"  .± NaN  * u"keV"
                        σ[d,i,p] = NaN* u"keV"  .± NaN  * u"keV"
                        fwhm[d,i,p] = NaN* u"keV"  .± NaN  * u"keV"
                        background[d,i,p] = NaN ± NaN
                        n[d,i,p] = NaN ± NaN
                        step_amplitude[d,i,p] = NaN ± NaN
                        pvalue[d,i,p] = NaN
                    end
                end
            else
                skew_frac[d,i,:] = fill(NaN ± NaN, npeaks)
                skew_width[d,i,:] = fill(NaN ± NaN, npeaks)
                µ[d,i,:] = fill(NaN* u"keV"  .± NaN  * u"keV", npeaks)
                σ[d,i,:] = fill(NaN* u"keV"  .± NaN  * u"keV", npeaks)
                fwhm[d,i,:] = fill(NaN* u"keV"  .± NaN  * u"keV", npeaks)
                background[d,i,:] = fill(NaN ± NaN, npeaks)
                n[d,i,:] = fill(NaN ± NaN,npeaks)
                step_amplitude[d,i,:] = fill(NaN ± NaN, npeaks)
                pvalue[d,i,:] = fill(NaN, npeaks)
                @info "run $i, det $det ($d) missing"
            end
        end
    end
    # # sort according to energy 
    th228_literature = sort([probdict_part[1][Symbol(dets_ged[1])][e_type].cal.peaks..., 1592u"keV", 2103u"keV"]) # get literature values mit denen gefittet wurde. interpolation st Qbb
    IdxSort = sortperm(µ[:1,1,:])
    th228_names  = th228_names[IdxSort]

    # sort according to energy 
    µ_sort= µ[:, :, IdxSort]
    σ_sort= σ[:, :, IdxSort]
    skew_frac_sort = skew_frac[:, :, IdxSort]
    skew_width_sort = skew_width[:, :, IdxSort]
    n_sort = n[:, :, IdxSort]
    background_sort = background[:, :, IdxSort]
    step_amplitude_sort = step_amplitude[:, :, IdxSort]
    fwhm_sort= fwhm[:, :, IdxSort]
    pvalue_sort = pvalue[:, :, IdxSort]
    #pval = reshape(pval_sort,:)

    # prepare data for plot 
    residual = permutedims(permutedims(µ_sort,(3,1,2)) .- th228_literature,(2,3,1))

    FitPars = (µ = µ_sort, 
                σ = σ_sort, 
                fwhm = fwhm_sort, 
                background = background_sort, 
                n = n_sort, 
                step_amplitude = step_amplitude_sort, 
                skew_frac = skew_frac_sort,
                skew_width = skew_width_sort,
                pvalue = pvalue_sort,
                residual = residual)

    MetaData = (dets_ged = dets_ged, dets_type = dets_type, npeaks = npeaks, th228_names = th228_names, th228_literature = th228_literature, nruns = nruns, partinfo = partinfo)
    if !ispath(path_data)
        mkdir("$path_data")
    end
    jldsave(file_data, FitPars = FitPars, MetaData = MetaData)
    @info "saved ecal fit pars to $file_data"
    return FitPars, MetaData
end

function get_calcurvefit_rpars_partition(partition::DataPartition; reload::Bool = false, e_type::Symbol = :e_cusp_ctc )
    path_data = "$(@__DIR__)/EcalData/"
    file_data = path_data * "Ecal_calcurve_rpars_$(partition).jld2"

    if isfile(file_data) && reload == false && e_type == :e_cusp_ctc
        @info "load calibration fit parameter from file $file_data..."
        file = jldopen(file_data)
        return file["CalPars"], file["MetaData"]
    end
    @info "load calibration fit parameter from rpars for $partition"
    l200 = LegendData(:l200)
    partinfo = partitioninfo(l200)[partition]
    filekey = start_filekey(l200, (partinfo[1].period, partinfo[1].run, :cal)) 
    chinfo = Table(channelinfo(l200, filekey; system=:geds, only_processable=true))
    dets_ged = chinfo.detector
    dets_type = chinfo.det_type

    # load all ProbDicts (for all period-run combination in selected partition). speed up load probdict for partition 
    probdict_part = [l200.par.rpars.ecal[entr.period, entr.run] for entr in partinfo] 
    nruns = length(probdict_part)
    # load peaks 
    th228_literature = sort([probdict_part[1][Symbol(dets_ged[1])][e_type].cal.peaks...]) # get literature values
    npeaks = length(th228_literature)
    good_det = ones(length(dets_ged)) .* true

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
            if !haskey(probdict_part[i],det) 
                loadLog = false 
                @info "det $det ($d) missing"
            elseif haskey(probdict_part[i],det) 
                if !haskey(probdict_part[i][det],e_type)
                    loadLog = false 
                    @info "$(e_type) of det $det ($d) missing"
                end      
            end
            if loadLog == true
                #@info "load det $det ($d)"
                pvalues[d,i] = probdict_part[i][dets_ged[d]][e_type].cal.gof.pvalue
                peaks_tmp =  probdict_part[i][dets_ged[d]][e_type].cal.peaks
                if length(peaks_tmp) !== npeaks
                    peak_idx = reduce(vcat, [findall(peaks_tmp[peak] .== th228_literature) for peak in eachindex(peaks_tmp)])
                    peaks_literature_keV[d,i,peak_idx] = peaks_tmp 
                    µ_ADC[d,i,peak_idx] =  probdict_part[i][dets_ged[d]][e_type].cal.µ
                    residuals_norm[d, i, peak_idx] = probdict_part[i][dets_ged[d]][e_type].cal.gof.residuals_norm
                else
                    peaks_literature_keV[d,i,:] = peaks_tmp 
                    µ_ADC[d,i,:] =  probdict_part[i][dets_ged[d]][e_type].cal.µ
                    residuals_norm[d,i,:] = probdict_part[i][dets_ged[d]][e_type].cal.gof.residuals_norm
                end
                calfunc_tmp = ljl_propfunc(probdict_part[i][dets_ged[d]][e_type].cal.func)
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

    MetaData = (dets_ged = dets_ged, good_det = good_det, dets_type = dets_type, npeaks = npeaks, th228_literature = th228_literature, nruns = nruns, partinfo = partinfo)
    if !ispath(path_data)
        mkdir("$path_data")
    end
    jldsave(file_data, CalPars = CalPars, MetaData = MetaData)
    @info "saved ecal cal pars to $file_data"
    return CalPars, MetaData
end

function get_fwhmcurvefit_rpars_partition(partition::DataPartition; reload::Bool = false, e_type::Symbol = :e_cusp_ctc )
    path_data = "$(@__DIR__)/EcalData/"
    file_data = path_data * "Ecal_fwhmcurve_rpars_$(partition).jld2"

    if isfile(file_data) && reload == false && e_type == :e_cusp_ctc
        @info "load fwhm fit parameter from file $file_data..."
        file = jldopen(file_data)
        return file["FWHMPars"], file["MetaData"]
    end
    @info "load fwhm fit parameter from rpars for $partition"
    l200 = LegendData(:l200)
    partinfo = partitioninfo(l200)[partition]
    filekey = start_filekey(l200, (partinfo[1].period, partinfo[1].run, :cal)) 
    chinfo = Table(channelinfo(l200, filekey; system=:geds, only_processable=true))
    dets_ged = chinfo.detector
    dets_type = chinfo.det_type

    # load all ProbDicts (for all period-run combination in selected partition). speed up load probdict for partition 
    probdict_part = [l200.par.rpars.ecal[entr.period, entr.run] for entr in partinfo] 
    nruns = length(probdict_part)
    # load peaks 
    th228_literature = sort([probdict_part[1][Symbol(dets_ged[1])][e_type].cal.peaks...]) # get literature values
    npeaks = length(th228_literature)
    good_det = ones(length(dets_ged)) .* true

    # re-shuffle results for a specific detector
    pvalues = zeros((length(dets_ged), nruns)) .* NaN   
    fwhm_fit  = zeros((length(dets_ged), nruns, npeaks)).* NaN .* u"keV" .± NaN .* u"keV"
    fwhm_qbb  = zeros((length(dets_ged), nruns)).* NaN .* u"keV" .± NaN .* u"keV"
    µ   = zeros((length(dets_ged), nruns, npeaks))  .* NaN  .*u"keV"  .± NaN     .*u"keV"
    peaks_literature_keV   = zeros((length(dets_ged), nruns, npeaks))  .* NaN   .*u"keV"
    residuals_norm  = zeros((length(dets_ged), nruns, npeaks)).* NaN  
    

    for i = 1:nruns
        @info "run $i"
        for (d, det) in enumerate(dets_ged)
            loadLog = true
            if !haskey(probdict_part[i],det) 
                loadLog = false 
                @info "det $det ($d) missing"
            elseif haskey(probdict_part[i],det) 
                if !haskey(probdict_part[i][det],e_type)
                    loadLog = false 
                    @info "$(e_type) of det $det ($d) missing"
                end      
            end
            if loadLog == true
               
                pvalues[d,i] = probdict_part[i][dets_ged[d]][e_type].fwhm.gof.pvalue 
                fwhm_qbb[d,i] = probdict_part[i][dets_ged[d]][e_type].fwhm.qbb 
                peaks_tmp =  probdict_part[i][dets_ged[d]][e_type].fwhm.peaks

                µ_ADC = probdict_part[i][dets_ged[d]][e_type].cal.µ
                Tbl = Table(e_cusp = [µ_ADC...], qdrift = fill(0, length(µ_ADC)))
                cal_func = ljl_propfunc(probdict_part[i][det][e_type].cal.func) 

                if length(peaks_tmp) !== npeaks
                    @info "not all peaks "
                    peak_idx = reduce(vcat, [findall(peaks_tmp[peak] .== th228_literature) for peak in eachindex(peaks_tmp)])
                    peaks_literature_keV[d,i,peak_idx] = peaks_tmp 
                    fwhm_fit[d,i,peak_idx] =  probdict_part[i][dets_ged[d]][e_type].fwhm.fwhm
                    residuals_norm[d, i, peak_idx] = probdict_part[i][dets_ged[d]][e_type].fwhm.gof.residuals_norm
                    µ[d,i,peak_idx] = collect(cal_func.(Tbl))# calibrated peak position for period given run and detector
                else
                    peaks_literature_keV[d,i,:] = peaks_tmp 
                    fwhm_fit[d,i,:] =  probdict_part[i][dets_ged[d]][e_type].fwhm.fwhm  
                    residuals_norm[d,i,:] = probdict_part[i][dets_ged[d]][e_type].fwhm.gof.residuals_norm
                    µ[d,i,:] = collect(cal_func.(Tbl))# calibrated peak position for period given run and detector
                end
            else
                @info "skip run $i det $det ($d)"
                pvalues[d,i] = NaN
                fwhm_fit[d,i,:] = fill(NaN, npeaks) .* u"keV" .± NaN .* u"keV"
                fwhm_qbb[d,i] = NaN * u"keV" .± NaN * u"keV"    
                residuals_norm[d,i,:] = fill(NaN, npeaks)    
                µ[d,i,:]  = fill(NaN, npeaks) .* u"keV" .± NaN .* u"keV"
                good_det[d] = false
            end
        end
    end
    FWHMPars = (µ = µ, 
                fwhm_fit = fwhm_fit,
                fwhm_qbb = fwhm_qbb,
              peaks_literature_keV = peaks_literature_keV,
              th228_literature = th228_literature,
              residuals_norm = residuals_norm,
              pvalue = pvalues)

    MetaData = (dets_ged = dets_ged, good_det = good_det, dets_type = dets_type, npeaks = npeaks, th228_literature = th228_literature, nruns = nruns, partinfo = partinfo)
    if !ispath(path_data)
        mkdir("$path_data")
    end
    jldsave(file_data, FWHMPars = FWHMPars, MetaData = MetaData)
    @info "saved ecal cal pars to $file_data"
    return FWHMPars, MetaData
end
# function get_peakfit_residuals(partition::DataPartition; reload::Bool = false,  e_type::Symbol= :e_cusp_ctc )
#     path_data = "$(@__DIR__)/EcalData/"
#     file_data = path_data * "Ecal_peakfit_residuals_$(partition)_simplecal.jld2"

#     if isfile(file_data) && reload == false && e_type == :e_cusp_ctc
#         @info "load peak fit parameter from file $file_data..."
#         file = jldopen(file_data)
#         return file["residual"], file["residual_norm"], file["MetaData"]
#     end
#     @info "load peak fit parameter from rpars for $partition"
#     l200 = LegendData(:l200)
#     partinfo = partitioninfo(l200)[partition]
#     filekey = start_filekey(l200, (partinfo[1].period, partinfo[1].run, :cal)) 
#     chinfo = Table(channelinfo(l200, filekey; system=:geds, only_processable=true))
#     dets_ged = chinfo.detector

#     # load all ProbDicts (for all period-run combination in selected partition). speed up load probdict for partition 
#     probdict_part = [l200.par.rpars.ecal[entr.period, entr.run] for entr in partinfo];
#     nruns = length(probdict_part)

#     # get peak information 
#     th228_literature = sort(probdict_part[1][Symbol(dets_ged[1])][e_type].cal.peaks) # get literature values used in fit
#     npeaks = length(th228_literature)

#     # load fit parameter and calibrate 
#     µ_all = ones(length(dets_ged), nruns, npeaks) .* NaN * 1u"keV" .± NaN * 1u"keV"

#     for i = 1:nruns
#         for (d, det) in enumerate(dets_ged)
#             loadLog = true
#             if !haskey(probdict_part[i],det) 
#                 loadLog = false 
#             elseif haskey(probdict_part[i],det) 
#                 if !haskey(probdict_part[i][det],e_type)
#                     loadLog = false 
#                 end
#             end

#             if loadLog == true
#                 @info "run $i -  det $det ($d)"
#                 µ_ADC = probdict_part[i][det][e_type].cal.µ
#                 Tbl = Table(e_cusp = µ_ADC, qdrift = fill(0.0,length(µ_ADC)))
#                 cal_func = ljl_propfunc(probdict_part[i][det][e_type].cal.func)
#                 peakpos_cal = sort(collect(cal_func.(Tbl))) # calibrated peak position for period given run and detector 
#                 if length(peakpos_cal) !== npeaks
#                     # find missing peak and put NaN
#                     peakpos_cal_tmp = Vector{Quantity{Measurement{Float64}}}(undef,npeaks)
#                     for p = 1:npeaks
#                         idx = findfirst(abs.(peakpos_cal .- th228_literature[p]) .<= 10u"keV")
#                         if !isnothing(idx)
#                             @info "peak missing: run $i, det $det ($d), peak $(th228_literature[p]), replace with NaN"
#                             peakpos_cal_tmp[p] = peakpos_cal[idx]
#                         else
#                             peakpos_cal_tmp[p] =  NaN * 1u"keV" ±  NaN * 1u"keV"
#                         end
#                     end
#                     peakpos_cal = peakpos_cal_tmp
#                 end
#                 µ_all[d,i,:] = peakpos_cal
#             else
#                 @info "run $i, det $d $det missing"
#                 µ_all[d,i,:] = fill(NaN * 1u"keV" ±  NaN * 1u"keV", npeaks)
#             end
#         end
#     end
#     # prepare data for plot 
#     residual = permutedims(permutedims(µ_all,(3,1,2)) .- th228_literature,(2,3,1))
#     residual_norm = mvalue.(residual) ./ muncert.(residual)
   
#     MetaData = (dets_ged = dets_ged, npeaks = npeaks, nruns = nruns, partinfo = partinfo, th228_literature = th228_literature)
#     if !ispath(path_data)
#         mkdir("$path_data")
#     end
#     jldsave(file_data, residual = residual, residual_norm = residual_norm, MetaData = MetaData)
#     @info "saved ecal residual to $file_data"
#     return residual, residual_norm, MetaData
# end

