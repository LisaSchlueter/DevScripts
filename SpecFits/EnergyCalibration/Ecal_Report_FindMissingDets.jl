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
partition = 3
e_type = :e_cusp_ctc
@info "load peak fit parameter from rpars for $partition"
l200 = LegendData(:l200)
partinfo = partitioninfo(l200)[DataPartition(partition)]
filekey = start_filekey(l200, (partinfo[1].period, partinfo[1].run, :cal)) 
chinfo = Table(channelinfo(l200, filekey; system=:geds, only_processable=true))
dets_ged = chinfo.detector
dets_type = chinfo.det_type
# load all ProbDicts (for all period-run combination in selected partition). speed up load probdict for partition 
probdict_part = [l200.par.rpars.ecal[entr.period, entr.run] for entr in partinfo] 
nruns = length(probdict_part)
# report path 
report_path = "$(@__DIR__)/Report/"
if !ispath(report_path)
    mkdir("$report_path")
end 

# find missing dets and write report
missing_dets_mat = [[haskey(probdict_part[i], dets_ged[det]) for i in 1:nruns] for det in eachindex(dets_ged)]
# write report. 
report_name = report_path * "Ecal_MissingDets_part$partition.md"
Tbl_dets = DataFrame(detector = String[], periods_run = String[])
for det in eachindex(dets_ged)
    if any(missing_dets_mat[det] .== 0 )
        @info "missing det $(dets_ged[det]) ($det) \n"
        if all(missing_dets_mat[det] .== 0)
            periods_run_str = "all runs in partition $partition"
        else
            periods_run_str = join(string.(partinfo[findall(missing_dets_mat[det] .== 0)].period) .* "-" .* string.(partinfo[findall(missing_dets_mat[det] .== 0)].run) .* ", ") 
        end
        push!(Tbl_dets, (detector = "$(dets_ged[det]) ($det)" , periods_run = periods_run_str))
    end
end
report = lreport()
lreport!(report, "# Partition $partition")
lreport!(report, "## Missing detectors (energy calibration peak-fits):")
lreport!(report, Tbl_dets)
writelreport(report_name, report)

# find missing ctc and write report
missing_ctc_mat = zeros(length(dets_ged), nruns)
for det in eachindex(dets_ged)  
    for i in 1:nruns
        if missing_dets_mat[det][i] == true # detector is NOT missing, data exists 
            if haskey(probdict_part[i][dets_ged[det]],e_type) == false # but no charge trapping correction
                missing_ctc_mat[det,i] = false
            else
                missing_ctc_mat[det,i] = true
            end
        else
            missing_ctc_mat[det,i] = NaN
        end
    end
end

if any(missing_ctc_mat .== false) == false # all existing dets have ctc 
   Tbl_ctc = "all existing dets have ctc in partition $partition"
else
    Tbl_ctc = DataFrame(detector = String[], periods_run = String[])
    for det in eachindex(dets_ged)
        if any(missing_ctc_mat[det,:] .== 0 )
            @info "missing ctc for  $(dets_ged[det]) ($det) \n"
            if all(missing_ctc_mat[det] .== 0)
                periods_run_str = "all runs in partition $partition"
            else
                periods_run_str = join(string.(partinfo[findall(missing_ctc_mat[det,:] .== 0)].period) .* "-" .* string.(partinfo[findall(missing_ctc_mat[det,:] .== 0)].run) .* ", ") 
            end
            push!(Tbl_ctc, (detector = "$(dets_ged[det]) ($det)" , periods_run = periods_run_str))
        end
    end
end

lreport!(report, "## Existing detectors, but no charge trapping correction (energy calibration peak-fits):")
lreport!(report, Tbl_ctc)
writelreport(report_name, report)


