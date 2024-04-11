using LegendDataManagement
l200 = LegendData(:l200)

# get pvalues runwise: input period, run 
"""
    _get_pval_runwise(period::DataPeriod,run::DataRun,e_type::Symbol;th228_lines::Vector{Symbol}=Symbol[],mode::Symbol=:cal) 
get pvalues for all detectors and all fits, for a given period, run, mode and energy type 
"""
function _get_pval_runwise(period::DataPeriod,run::DataRun,e_type::Symbol;th228_lines::Vector{Symbol}=Symbol[],mode::Symbol=:cal) 
    filekey = start_filekey(l200, (period, run, mode))
    if isempty(th228_lines)
        energy_config = dataprod_config(l200).energy(filekey)
        th228_lines = Symbol.(energy_config.default.th228_names)
    end
    dets_ged = map(x-> x.label,channelinfo(l200,filekey,system = :geds, only_processable = true).detector)
   
    ecal_vals = l200.par.rpars.ecal[period, run]
    pvalues = Dict()
    for line in th228_lines
        pvalues[line] = [ecal_vals[det][e_type].fit[line].gof.pvalue for det in dets_ged if haskey(ecal_vals, det) && haskey(ecal_vals[det], e_type)]
    end

    return pvalues, th228_lines
end

"""
_get_pval_periodwise(period::DataPeriod,e_type::Symbol;th228_lines::Vector{Symbol}=Symbol[],mode::Symbol=:cal) 
get pvalues for all detectors and all fits, for a given period, run, mode and energy type 
"""
function _get_pval_periodwise(period::DataPeriod,e_type::Symbol;th228_lines::Vector{Symbol}=Symbol[],mode::Symbol=:cal, det = "") 
    runs_period = search_disk(DataRun,l200.tier[:jldsp,mode,period])
    pvalues =  Vector{Any}(undef, length(runs_period))

    for (i,run) in enumerate(runs_period)
        if haskey(l200.par.rpars.ecal[period],Symbol(run))
            @debug "period $(period), run $(run)"
            filekey = start_filekey(l200, (period, run, mode))
            if isempty(th228_lines)
                energy_config = dataprod_config(l200).energy(filekey)
                th228_lines = Symbol.(energy_config.default.th228_names)
            end
            dets_ged = map(x-> x.label,channelinfo(l200,filekey,system = :geds, only_processable = true).detector)
            ecal_vals = l200.par.rpars.ecal[period, run]
            pvalues[i] = Dict()
            for line in th228_lines
                if isa(det,DetectorId)
                    pvalues[i][line] = ecal_vals[det][e_type].fit[line].gof.pvalue 
                else
                  pvalues[i][line] = [ecal_vals[det][e_type].fit[line].gof.pvalue for det in dets_ged if haskey(ecal_vals, det) && haskey(ecal_vals[det], e_type)]
                end
            end
        else
            deleteat!(runs_period,i)
        end 
    end
    return pvalues, th228_lines, runs_period
end

#=
"""
    _get_pval_partwise(partition::Symbol,mode::Symbol,e_type::Symbol)
get pvalues for all detectors and all fits, for a given partition, mode and energy type
doesnt work at, because ppars empty
"""
function _get_pval_partwise(partition::Symbol,mode::Symbol,e_type::Symbol;th228_lines::Vector{Symbol}=Symbol[]) 
   filekey = start_filekey(l200, (period, run, mode))
    if isempty(th228_lines)
        energy_config = dataprod_config(l200).energy(filekey)
        th228_lines = Symbol.(energy_config.default.th228_names)
    end

    dets_ged = map(x-> x.label,channelinfo(l200,filekey,system = :geds, only_processable = true).detector)
    runs = search_disk(DataRun,l200.tier[:jldsp,mode,period])
    ecal_vals = l200.par.ppars.ecal[part]

    pvalues = Dict()
    for line in th228_lines
    line =th228_lines[1]
        pvalues[line] = [ecal_vals[det][e_type].fit[line].gof.pvalue for det in dets_ged if haskey(ecal_vals, det) && haskey(ecal_vals[det], e_type)]
    end

    return pvalues, th228_lines
end
=#


function _get_pvaldist_pltname(period::DataPeriod,run::DataRun,e_type::Symbol;format::String="pdf")
    path_plot = "PackageDevScripts/DevelopSpecFits/SanityPlots/plots/pvaldist_runwise/"
    path_rel_plot = relpath(path_plot,pwd())
    path_abs_plot = pwd() * "/" * path_rel_plot * "/"
    if !ispath(path_abs_plot)
        mkdir(path_abs_plot)
    end
    fname = path_abs_plot * "peakfit_pvaldist_runwise_$(e_type)_$(period)_$(run).$(format)"
    return fname 
end
function _get_pvaldist_pltname(period::DataPeriod,e_type::Symbol;format::String="png")
    path_plot = "PackageDevScripts/DevelopSpecFits/SanityPlots/plots/pvaldist_periodwise/"
    path_rel_plot = relpath(path_plot,pwd())
    path_abs_plot = pwd() * "/" * path_rel_plot * "/"
    if !ispath(path_abs_plot)
        mkdir(path_abs_plot)
    end
    fname = path_abs_plot * "peakfit_pvaldist_periodwise_$(e_type)_$(period).$(format)"
    return fname 
end
function _get_pvaldist_pltname(period::DataPeriod;format::String="png")
    path_plot = "PackageDevScripts/DevelopSpecFits/SanityPlots/plots/pvaldist_periodwise/"
    path_rel_plot = relpath(path_plot,pwd())
    path_abs_plot = pwd() * "/" * path_rel_plot * "/"
    if !ispath(path_abs_plot)
        mkdir(path_abs_plot)
    end
    fname = path_abs_plot * "peakfit_pvaldist_badfits_periodwise_$(period).$(format)"
    return fname 
end
function _get_pvaldist_rptname(period::DataPeriod;format::String="md")
    path_plot = "PackageDevScripts/DevelopSpecFits/SanityPlots/report/pvaldist_periodwise/"
    path_rel_plot = relpath(path_plot,pwd())
    path_abs_plot = pwd() * "/" * path_rel_plot * "/"
    if !ispath(path_abs_plot)
        mkdir(path_abs_plot)
    end
    fname = path_abs_plot * "peakfit_pvaldist_periodwise_$(period).$(format)"
    return fname 
end

function _get_bestfit_pltname(period::DataPeriod,run::DataRun,det::DetectorId;format::String="png")
    path_plot = "PackageDevScripts/DevelopSpecFits/SanityPlots/plots/specfit_res/"
    path_rel_plot = relpath(path_plot,pwd())
    path_abs_plot = pwd() * "/" * path_rel_plot * "/"
    if !ispath(path_abs_plot)
        mkdir(path_abs_plot)
    end
    fname = path_abs_plot * "specfit_res_$(det)_$(period)_$(run).$(format)"
    return fname 
end

function _get_bestfit_pltname(period::DataPeriod,run::DataRun,det::DetectorId,e_type::Symbol;format::String="png")
    path_plot = "PackageDevScripts/DevelopSpecFits/SanityPlots/plots/specfit_res/"
    path_rel_plot = relpath(path_plot,pwd())
    path_abs_plot = pwd() * "/" * path_rel_plot * "/"
    if !ispath(path_abs_plot)
        mkdir(path_abs_plot)
    end
    fname = path_abs_plot * "specfit_res_$(det)_$(period)_$(run)_$(e_type).$(format)"
    return fname 
end
function _get_calicurve_pltname(period::DataPeriod,run::DataRun,det::DetectorId,e_type::Symbol;format::String="png")
    path_plot = "PackageDevScripts/DevelopSpecFits/SanityPlots/plots/specfit_calicurve/"
    path_rel_plot = relpath(path_plot,pwd())
    path_abs_plot = pwd() * "/" * path_rel_plot * "/"
    if !ispath(path_abs_plot)
        mkdir(path_abs_plot)
    end
    fname = path_abs_plot * "specfit_calicurve_$(det)_$(period)_$(run)_$(e_type).$(format)"
    return fname 
end


