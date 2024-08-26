#= 
fit paramater of peaks fit (calibration data)
look peak fits with drastic outliers (in terms of fit parameter values)
stability investigation 
=#
using Measures
using Plots, Printf, LaTeXStrings
using JLD2
using LegendDataManagement
using LegendDataManagement.LDMUtils
using Unitful
using Measurements
using Measurements: value as mvalue, uncertainty as muncert
using TypedTables
using StatsBase
using DataFrames

#select data and dsp output
periods = [3, 4, 6, 7, 8, 9]
e_type = :e_zac

# load data 
@info "load peak fit parameter from rpars for period $period"
l200 = LegendData(:l200)
period_runs = [(period, run)   for period in DataPeriod.(periods) for run in search_disk(DataRun,l200.tier[:raw, :cal, DataPeriod(period)])]

# get all detectors 
chinfos = [Table(channelinfo(l200, (period, run, :cal) ; system=:geds, only_processable=true)) for period in DataPeriod.(periods) for run in search_disk(DataRun, l200.tier[:raw, :cal, DataPeriod(period)]) ]
dets_ged = [chinfo.detector for chinfo in chinfos]
dets_ged = unique(vcat(dets_ged...))   
ch_ged = [chinfo.channel for chinfo in chinfos]
ch_ged = unique(vcat(ch_ged...)) 

# load all ProbDicts (for all period-run combination in selected partition). speed up load probdict for partition 
probdict= [l200.par.rpars.ecal[period_runs[i][1], period_runs[i][2]] for i in eachindex(period_runs)] 

# find missing dets 
missing_dets_mat = [[haskey(probdict[r], dets_ged[det]) for r in eachindex(period_runs)] for det in eachindex(dets_ged)]

# write report: 
report_path = "$(@__DIR__)/Report/"
if !ispath(report_path)
    mkdir("$report_path")
end 

report_name = report_path * "Ecal_MissingDets_periods$(join(string.(periods))).md"
Tbl_dets = DataFrame(detector = String[], periods_run = String[])
for det in eachindex(dets_ged)
    if any(missing_dets_mat[det] .== 0 )
        @info "missing det $(dets_ged[det]) ($det) \n"
        if all(missing_dets_mat[det] .== 0)
            periods_run_str = "all runs in periods $(join(string(periods)))"
        else
            periods_run_str = ""
            for missingIdx in findall(missing_dets_mat[det] .== 0)
                periods_run_str *= string(period_runs[missingIdx][1]) * "-" * string(period_runs[missingIdx][2]) * ", "
            end
        end
        push!(Tbl_dets, (detector = "$(dets_ged[det]) ($det)" , periods_run = periods_run_str))
    end
end
report = lreport()
lreport!(report, "# Periods $periods")
lreport!(report, "## Missing (energy calibration peak-fit) rpars for detectors $e_type:")
lreport!(report, Tbl_dets)
writelreport(report_name, report)

# # find missing ctc and write report
# missing_ctc_mat = zeros(length(dets_ged), nruns)
# for det in eachindex(dets_ged)  
#     for i in 1:nruns
#         if missing_dets_mat[det][i] == true # detector is NOT missing, data exists 
#             if haskey(probdict_part[i][dets_ged[det]],e_type) == false # but no charge trapping correction
#                 missing_ctc_mat[det,i] = false
#             else
#                 missing_ctc_mat[det,i] = true
#             end
#         else
#             missing_ctc_mat[det,i] = NaN
#         end
#     end
# end

# if any(missing_ctc_mat .== false) == false # all existing dets have ctc 
#    Tbl_ctc = "all existing dets have ctc in partition $partition"
# else
#     Tbl_ctc = DataFrame(detector = String[], periods_run = String[])
#     for det in eachindex(dets_ged)
#         if any(missing_ctc_mat[det,:] .== 0 )
#             @info "missing ctc for  $(dets_ged[det]) ($det) \n"
#             if all(missing_ctc_mat[det] .== 0)
#                 periods_run_str = "all runs in partition $partition"
#             else
#                 periods_run_str = join(string.(partinfo[findall(missing_ctc_mat[det,:] .== 0)].period) .* "-" .* string.(partinfo[findall(missing_ctc_mat[det,:] .== 0)].run) .* ", ") 
#             end
#             push!(Tbl_ctc, (detector = "$(dets_ged[det]) ($det)" , periods_run = periods_run_str))
#         end
#     end
# end

# lreport!(report, "## Existing detectors, but no charge trapping correction (energy calibration peak-fits):")
# lreport!(report, Tbl_ctc)
# writelreport(report_name, report)


