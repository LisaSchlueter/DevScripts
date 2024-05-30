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
include("$(@__DIR__)/utils_ecal.jl")

#select data and dsp output
partition = 1
e_type = :e_cusp_ctc
FitPars, MetaData = get_peakfit_rpars_partition(DataPartition(partition); reload = false, e_type = e_type);


# FIND OUTLIERS
# look only at detector that are NOT missing. --> missing dets are already reported in Ecal_FindMissingDets.jl
dets_Idx = findall([!all(isnan.(FitPars.µ[det,:,:])) for det in eachindex(MetaData.dets_ged)])
dets_ged = MetaData.dets_ged[dets_Idx]

µ = FitPars.µ[dets_Idx,:,:]
µ_med = [[median(filter(!isnan,µ[det,:,p])) for p in eachindex(MetaData.th228_literature)] for det in eachindex(dets_Idx)]
µ_outlier_prct = [round.(100 .* abs.(mvalue.((µ[det,:,:] .- µ_med'[det]) ./µ_med'[det])), digits = 3) for det in eachindex(dets_Idx)]

σ = FitPars.σ[dets_Idx,:,:]
σ_med = [[median(filter(!isnan,σ[det,:,p])) for p in eachindex(MetaData.th228_literature)] for det in eachindex(dets_Idx)]
σ_outlier_prct = [round.(100 .* abs.(mvalue.((σ[det,:,:] .- σ_med'[det]) ./σ_med'[det])), digits = 3) for det in eachindex(dets_Idx)]

skew_frac = FitPars.skew_frac[dets_Idx,:,:]
skew_frac_med = [[median(filter(!isnan,skew_frac[det,:,p])) for p in eachindex(MetaData.th228_literature)] for det in eachindex(dets_Idx)]
skew_frac_outlier_prct = [round.(100 .* abs.(mvalue.((skew_frac[det,:,:] .- skew_frac_med'[det]) ./skew_frac_med'[det])), digits = 3) for det in eachindex(dets_Idx)]

skew_width = FitPars.skew_width[dets_Idx,:,:]
skew_width_med = [[median(filter(!isnan,skew_width[det,:,p])) for p in eachindex(MetaData.th228_literature)] for det in eachindex(dets_Idx)]
skew_width_outlier_prct = [round.(100 .* abs.(mvalue.((skew_width[det,:,:] .- skew_width_med'[det]) ./skew_width_med'[det])), digits = 3) for det in eachindex(dets_Idx)]

background = FitPars.background[dets_Idx,:,:]
background_med = [[median(filter(!isnan,background[det,:,p])) for p in eachindex(MetaData.th228_literature)] for det in eachindex(dets_Idx)]
background_outlier_prct = [round.(100 .* abs.(mvalue.((background[det,:,:] .- background_med'[det]) ./background_med'[det])), digits = 3) for det in eachindex(dets_Idx)]

step_amplitude = FitPars.step_amplitude[dets_Idx,:,:]
step_amplitude_med = [[median(filter(!isnan,step_amplitude[det,:,p])) for p in eachindex(MetaData.th228_literature)] for det in eachindex(dets_Idx)]
step_amplitude_outlier_prct = [round.(100 .* abs.(mvalue.((step_amplitude[det,:,:] .- step_amplitude_med'[det]) ./step_amplitude_med'[det])), digits = 3) for det in eachindex(dets_Idx)]

n = FitPars.n[dets_Idx,:,:]
n_med = [[median(filter(!isnan,n[det,:,p])) for p in eachindex(MetaData.th228_literature)] for det in eachindex(dets_Idx)]
n_outlier_prct = [round.(100 .* abs.(mvalue.((n[det,:,:] .- n_med'[det]) ./n_med'[det])), digits = 3) for det in eachindex(dets_Idx)]

### -------
# report outlier. detector - by detector and by parameter
function _get_Tbl_outlier(outlier_prct::Vector{Matrix{Float64}},par_str::String, par::Array , par_median::Array; percent::Int64 = 100)
    Tbl_Outlier_par = DataFrame(detector = String[], period_run = String[], peak = String[], par = String[], median_runs = String[], diff2med = String[])
    for det in eachindex(dets_Idx)
        if any(outlier_prct[det] .> percent)
            outlier_idx = findall(outlier_prct[det] .> percent)
            outlier_prct[det][outlier_idx[1]]
            for i in eachindex(outlier_idx)
                @info  "$(MetaData.partinfo[outlier_idx[i][1]].period)-$(MetaData.partinfo[outlier_idx[i][1]].run)-$(MetaData.th228_names[outlier_idx[i][2]]) " 
                det_orig_idx = findfirst(string.(MetaData.dets_ged) .== string(dets_ged[det]))
                push!(Tbl_Outlier_par,( 
                                    detector = "$(dets_ged[det]) ($det_orig_idx)", 
                                    period_run = "$(MetaData.partinfo[outlier_idx[i][1]].period)-$(MetaData.partinfo[outlier_idx[i][1]].run)",
                                    peak = "$(MetaData.th228_names[outlier_idx[i][2]])" ,
                                    par = "$(par_str) = $(par[det,outlier_idx[i]])",
                                    median_runs = "< $(par_str)> = $(par_median[det][outlier_idx[i][2]])",
                                    diff2med = "$(round(outlier_prct[det][outlier_idx[i]], digits = 0))%"))               
            end
        end  
    end 
    if isempty(Tbl_Outlier_par)
        Tbl_Outlier_par = "no outliers found (threshold: $percent%)"
    end
    return Tbl_Outlier_par
end
Tbl_Outlier_µ = _get_Tbl_outlier(µ_outlier_prct,"µ",µ , µ_med; percent = 50)
Tbl_Outlier_σ = _get_Tbl_outlier(σ_outlier_prct,"σ",σ , σ_med;  percent = 100)
Tbl_Outlier_n = _get_Tbl_outlier(n_outlier_prct, "n", n, n_med;  percent = 200)
Tbl_Outlier_skew_frac = _get_Tbl_outlier(skew_frac_outlier_prct,"skew_frac", skew_frac, skew_frac_med;  percent = 100)
Tbl_Outlier_skew_width = _get_Tbl_outlier(skew_width_outlier_prct,"skew_width", skew_width, skew_width_med;  percent = 100)
Tbl_Outlier_background = _get_Tbl_outlier(background_outlier_prct,"background", background, background_med;  percent = 100)
Tbl_Outlier_step_amplitude = _get_Tbl_outlier(step_amplitude_outlier_prct,"step_amplitude", step_amplitude, step_amplitude_med;  percent = 100)

# write report 
# report path 
report_path = "$(@__DIR__)/Report/"
if !ispath(report_path)
    mkdir("$report_path")
end 
report_name = report_path * "FitPars_Outliers_part$partition.md"
report = lreport()
lreport!(report, "# Partition $partition - calibration peak fits: fit parameter outlier ")
lreport!(report, "## µ"); lreport!(report, Tbl_Outlier_µ)
lreport!(report, "## σ"); lreport!(report, Tbl_Outlier_σ)
lreport!(report, "## n"); lreport!(report, Tbl_Outlier_n)
lreport!(report, "## skew_frac"); lreport!(report, Tbl_Outlier_skew_frac)
lreport!(report, "## skew_width"); lreport!(report, Tbl_Outlier_skew_width)
lreport!(report, "## background"); lreport!(report, Tbl_Outlier_background)
lreport!(report, "## step_amplitude"); lreport!(report, Tbl_Outlier_step_amplitude)
writelreport(report_name, report)

