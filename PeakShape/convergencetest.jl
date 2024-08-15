using LegendDataManagement
using LegendDataManagement.LDMUtils
using LegendSpecFits
import LegendSpecFits: get_th228_fit_functions, hist_loglike, get_pseudo_prior
using LegendHDF5IO
using TypedTables
using Unitful, Measures
using Measurements: value as mvalue, uncertainty as muncert
using StructArrays, PropDicts
using Revise
using Plots
using Distributions 
using ValueShapes
using LaTeXStrings, Printf
using BAT,  Optim, InverseFunctions
include("../utils/utils_ecal.jl")
include("../utils/utils_plot.jl")
# data selection 

periods = [3,4]
path_data = "$(@__DIR__)/EcalData/"
file_data = path_data * "ecal_peakfit_rpars_$(join(string.(periods)))_convergence.jld2"

reload = false

if isfile(file_data) && reload == false 
    @info "load peak fit parameter from file $file_data..."
    file = jldopen(file_data)
    (period_runs, dets_ged, periods, niterations, nconverged, converged_frac) = file["period_runs"], file["dets_ged"], file["periods"], file["niterations"], file["nconverged"], file["converged_frac"]
else
    l200 = LegendData(:l200)
    chinfo = channelinfo(l200, runsel; system=:geds, only_processable=true)
    dets_ged  = chinfo.detector
    # load fit pars 
    pd_rpars = [l200.par.rpars.ecal[period, run] for period in DataPeriod.(periods) for run in search_disk(DataRun,l200.tier[:raw, :cal, DataPeriod(period)])][1:end-3]
    period_runs = [(period, run)   for period in DataPeriod.(periods) for run in search_disk(DataRun,l200.tier[:raw, :cal, DataPeriod(period)])][1:end-3]
    nruns = length(pd_rpars) 
    niterations = l200.par.rpars.ecal[period, run, det].e_cusp_ctc.fit[:converged].niterations  
    nconverged = zeros(nruns, length(dets_ged), length(niterations) )
    for r = 1:nruns 
        @info "run: $r"
        for (d, det) in enumerate(dets_ged)
            @info "det: $det"
            if !haskey(pd_rpars[r], det)
                @info "no data for det: $det"
                nconverged[r, d, :] = NaN .* ones(length(niterations))
                continue
            elseif !haskey(pd_rpars[r][det], :e_cusp_ctc)
                @info "no data e_cuscp_ctc for det: $det"
                nconverged[r, d, :] = NaN .* ones(length(niterations))
                continue
            end
            converged_tmp = hcat(pd_rpars[r][det].e_cusp_ctc.fit[:converged].converged...)
            nconverged[r, d, :] = vec(sum(converged_tmp, dims=2))
        end
    end

    # average over runs and detectors
    converged_frac = zeros(length(niterations))
    for i = 1:length(niterations)
        local nfits =  7*length(filter(isfinite, nconverged[:, : , i ]))
        converged_frac[i] = sum(filter(isfinite, nconverged[:, : , i ]))./nfits 
    end

    jldsave(file_data, period_runs = period_runs, dets_ged = dets_ged, periods = periods, niterations = niterations, nconverged = nconverged, converged_frac = converged_frac)
    @info "saved ecal fit pars to $file_data"
end 

_def_plot()
plot(niterations, 100 .*converged_frac, linewidth = 2, xlabel="Max. iterations (minimizer)", ylabel="Converged fits (%)", legend=false)
pltname = "$(@__DIR__)/plots/Peakfits_Convergence.png"
savefig(pltname)