using LegendSpecFits
using LegendDataManagement
using LegendHDF5IO, HDF5
using TypedTables
using PropertyFunctions
using StatsBase, Distributions, ValueShapes, BAT, Optim, LinearAlgebra
using InverseFunctions, IntervalSets, ForwardDiff
using Measurements, Unitful, Measures
using Format
using Plots, TypedTables
using LegendDataManagement.LDMUtils
l200 = LegendData(:l200)
include("../../../Packages/legend-julia-dataflow/src/utils.jl")
period = DataPeriod(3)
run = DataRun(4)

filekey = start_filekey(l200, (period, run, :cal))
chinfo = Table(channelinfo(l200, filekey; system=:geds, only_processable=true))# |> filterby(@pf $low_aoe_status .== :valid)
dets = chinfo.detector
chs = chinfo.channel 

# sel_det = 82 # p3, run 0, 
dets[sel_det]

# load data aoe 
e_type = :e_cusp
hitchfilename = get_hitchfilename(l200, filekey, chs[sel_det])
data_hit = LHDataStore(hitchfilename, "r");
pd_dataQC = data_hit["$(chs[sel_det])/dataQC/"]
a = collect(data_hit["$(chs[sel_det])/dataQC/a"])# get a
a_sel = isfinite.(a)
e = collect(data_hit["$(chs[sel_det])/dataQC/$e_type"]) # get e
a = a[a_sel]
e = e[a_sel]

# get energy calibration and calculate aoe 
pb_ecal = get_values(l200.par.rpars.ecal[period, run])[dets[sel_det]]
 #energy calibration parameter
cal_func_str = pb_ecal[e_type].cal.func
cal_func = ljl_propfunc(cal_func_str)
e_cal = collect(cal_func.(pd_dataQC)) # calibrated energy 
e_cal = e_cal[a_sel]
aoe = ustrip.(a ./ e_cal); # get aoe without units 

# get compton bands and fit each band, determine µ and σ
psd_config = dataprod_config(l200).psd(filekey).default
compton_bands  = psd_config.compton_bands
compton_window = psd_config.compton_window
compton_band_peakhists = generate_aoe_compton_bands(aoe, e_cal, compton_bands, compton_window)
result_fit, report_fit = fit_aoe_compton(compton_band_peakhists.peakhists, compton_band_peakhists.peakstats, compton_bands,; uncertainty=true)
µ = [result_fit[band].µ for band in compton_bands]
σ = [result_fit[band].σ for band in compton_bands]

# fit aoe corrections
result, report =  LegendSpecFits.fit_aoe_corrections(compton_bands, μ, σ; e_expression = cal_func_str)

#scatter(result.µ_compton.x,result.µ_compton.y)
plot(report.report_µ)
plot(report.report_σ)

#####  usage  #####
# test µ (Energy) correction 
aoe_µ_corr = aoe .-ljl_propfunc(result.µ_compton.func).(pd_dataQC)[a_sel] # test µ correction function 
scatter(e_cal[1:1000], aoe[1:1000],label = "uncorrected",xlabel = "Energy", ylabel = "A/E (keV^-1)")
scatter!(e_cal[1:1000], aoe_µ_corr[1:1000],label = "µ (energy) corrected")

# test σ normalization  
aoe_σ_corr = ljl_propfunc(result.σ_compton.func).(pd_dataQC)[a_sel] # test µ correction function 
all(isfinite,aoe_σ_corr)
any(x->x==0,aoe_σ_corr)
scatter(e_cal[1:1000], aoe_µ_corr[1:1000],label = "uncorrected",xlabel = "Energy", ylabel = "A/E (a.u)",c = :silver,markerstrokewidth=0,alpha = 0.8)
scatter!(e_cal[1:1000], aoe_µ_corr[1:1000]./aoe_σ_corr[1:1000],label = "µ (energy) corrected",alpha=0.2,markerstrokewidth=0)

# test total correction/normalization  
aoe_tot_corr = ljl_propfunc(result.func).(pd_dataQC)[a_sel] 
scatter(e_cal[1:5000], aoe_µ_corr[1:5000], label = "uncorrected",xlabel = "Energy", ylabel = "A/E (σ)",c = :silver,markerstrokewidth=0,alpha = 0.8)
scatter!(e_cal[1:5000], aoe_tot_corr[1:5000], label = "µ (energy) corrected, normalized σ",alpha=0.2,markerstrokewidth=0)

# try with already calibrated data
tbl_e = Table(a = a, e_cal = e_cal)
aoe_tot_corr2 = ljl_propfunc(result.func_ecal).(tbl_e)
scatter(e_cal[1:5000], aoe_µ_corr[1:5000], label = "uncorrected",xlabel = "Energy", ylabel = "A/E (σ)",c = :silver,markerstrokewidth=0,alpha = 0.8)
scatter!(e_cal[1:5000], aoe_tot_corr2[1:5000], label = "µ (energy) corrected, normalized σ",alpha=0.2,markerstrokewidth=0)


