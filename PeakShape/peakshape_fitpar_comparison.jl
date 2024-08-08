#= comparison gamma peak fits with and without linear background slope. 
--> for all dets or only ICPC, Coax, Bege,....
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
include("$(@__DIR__)/utils_ecal.jl")
include("$(@__DIR__)/utils_plot.jl")

# plot settings 
par_mode_all = [:fwhm, :σ, :n, :background, :step_amplitude, :skew_frac, :skew_width, :p_value] # all possible plot parameters
det_mode_all= [:all, :ppc, :icpc, :bege, :coax] # all possible "detector modes" -> subset of detectors
fontsize = 16 

# data selection 
e_type = :e_cusp_ctc
periods = [3, 4]
pltMode = :norm # :norm for normalized residuals, :keV for residuals in keV
exclPeaks = []# [:Tl208DEP, :Tl208SEP] # exclude peaks from plot

# load fit results 
l200 = LegendData(:l200)
FitPars, MetaData = get_rpars_peakfit(l200, periods; reload = false, e_type = e_type, cal_type = :ecal, centroid = true);
FitPars_tails, MetaData_tails = get_rpars_peakfit(l200, periods; reload = false, e_type = e_type, cal_type = :ecal_tails);
FitPars_bslope, MetaData_bslope = get_rpars_peakfit(l200, periods; reload = false, e_type = e_type, cal_type = :ecal_bslope);
FitPars_exp, MetaData_exp = get_rpars_peakfit(l200, periods; reload = false, e_type = e_type, cal_type = :ecal_exp);

# # plot: default vs high-energy tail 
# plothist_diff(FitPars_tails, FitPars, MetaData_tails, MetaData; par_mode = :chi2, det_mode = :all, xlim = (-5,5), exclPeaks = exclPeaks, saveplot = true)
# plothist_diff(FitPars_tails, FitPars, MetaData_tails, MetaData; par_mode = :p_value, det_mode = :all, xlim = (-0.5,0.5), exclPeaks = exclPeaks, saveplot = true)
# plothist_diff(FitPars_tails, FitPars, MetaData_tails, MetaData; par_mode = :fwhm, det_mode = :all, xlim = (-0.15,0.15), exclPeaks = exclPeaks, saveplot = true)
# plothist_diff(FitPars_tails, FitPars, MetaData_tails, MetaData; par_mode = :µ, det_mode = :all, xlim = (-0.2,0.2), exclPeaks = exclPeaks, saveplot = true)
# plothist_diff(FitPars_tails, FitPars, MetaData_tails, MetaData; par_mode = :skew_frac, det_mode = :all, xlim = (-0.1,0.1), exclPeaks = exclPeaks, saveplot = true)
# plothist_diff(FitPars_tails, FitPars, MetaData_tails, MetaData; par_mode = :skew_width, det_mode = :all, xlim = (-0.005,0.005), exclPeaks = exclPeaks, saveplot = true)
# plothist_diff(FitPars_tails, FitPars, MetaData_tails, MetaData; par_mode = :background, det_mode = :all, xlim= (-2, 2), exclPeaks = exclPeaks, saveplot = true)
# plothist_diff(FitPars_tails, FitPars, MetaData_tails, MetaData; par_mode = :step_amplitude, det_mode = :all, xlim= (-0.2, 0.2), exclPeaks = exclPeaks, saveplot = true)


# # plot: default vs background slope
# plothist_diff(FitPars_bslope, FitPars, MetaData_bslope, MetaData; par_mode = :chi2, det_mode = :all, xlim = (-5,5), exclPeaks = exclPeaks, saveplot = true)
# plothist_diff(FitPars_bslope, FitPars, MetaData_bslope, MetaData; par_mode = :p_value, det_mode = :all, xlim = (-0.5,0.5), exclPeaks = exclPeaks, saveplot = true)
# plothist_diff(FitPars_bslope, FitPars, MetaData_bslope, MetaData; par_mode = :fwhm, det_mode = :all, xlim = (-0.15,0.15), exclPeaks = exclPeaks, saveplot = true)
# plothist_diff(FitPars_bslope, FitPars, MetaData_bslope, MetaData; par_mode = :µ, det_mode = :all, xlim = (-0.2,0.2), exclPeaks = exclPeaks, saveplot = true)
# plothist_diff(FitPars_bslope, FitPars, MetaData_bslope, MetaData; par_mode = :skew_frac, det_mode = :all, xlim = (-0.1,0.1), exclPeaks = exclPeaks, saveplot = true)
# plothist_diff(FitPars_bslope, FitPars, MetaData_bslope, MetaData; par_mode = :skew_width, det_mode = :all, xlim = (-0.005,0.005), exclPeaks = exclPeaks, saveplot = true)
# plothist_diff(FitPars_bslope, FitPars, MetaData_bslope, MetaData; par_mode = :background, det_mode = :all, xlim= (-2, 2), exclPeaks = exclPeaks, saveplot = true)
# plothist_diff(FitPars_bslope, FitPars, MetaData_bslope, MetaData; par_mode = :skew_width, det_mode = :all, xlim = (-0.005,0.005), exclPeaks = exclPeaks, saveplot = true)
# plothist_diff(FitPars_bslope, FitPars, MetaData_bslope, MetaData; par_mode = :step_amplitude, det_mode = :all, xlim= (-0.2, 0.2), exclPeaks = exclPeaks, saveplot = true)

# plot: default vs background exp
plothist_diff(FitPars_exp, FitPars, MetaData_exp, MetaData; par_mode = :chi2, det_mode = :all, xlim = (-5,5), exclPeaks = exclPeaks, saveplot = true)
plothist_diff(FitPars_exp, FitPars, MetaData_exp, MetaData; par_mode = :p_value, det_mode = :all, xlim = (-0.5,0.5), exclPeaks = exclPeaks, saveplot = true)
plothist_diff(FitPars_exp, FitPars, MetaData_exp, MetaData; par_mode = :fwhm, det_mode = :all, xlim = (-0.15,0.15), exclPeaks = exclPeaks, saveplot = true)
plothist_diff(FitPars_exp, FitPars, MetaData_exp, MetaData; par_mode = :µ, det_mode = :all, xlim = (-0.2,0.2), exclPeaks = exclPeaks, saveplot = true)
plothist_diff(FitPars_exp, FitPars, MetaData_exp, MetaData; par_mode = :skew_frac, det_mode = :all, xlim = (-0.1,0.1), exclPeaks = exclPeaks, saveplot = true)
plothist_diff(FitPars_exp, FitPars, MetaData_exp, MetaData; par_mode = :skew_width, det_mode = :all, xlim = (-0.005,0.005), exclPeaks = exclPeaks, saveplot = true)
plothist_diff(FitPars_exp, FitPars, MetaData_exp, MetaData; par_mode = :background, det_mode = :all, xlim= (-2, 2), exclPeaks = exclPeaks, saveplot = true)
plothist_diff(FitPars_exp, FitPars, MetaData_exp, MetaData; par_mode = :skew_width, det_mode = :all, xlim = (-0.005,0.005), exclPeaks = exclPeaks, saveplot = true)
plothist_diff(FitPars_exp, FitPars, MetaData_exp, MetaData; par_mode = :step_amplitude, det_mode = :all, xlim= (-0.2, 0.2), exclPeaks = exclPeaks, saveplot = true)


FitPars_cen = merge(FitPars, (µ = FitPars.µ_cen,))

plothist_diff(FitPars_cen, FitPars, MetaData, MetaData; par_mode = :µ, 
    det_mode = :all, xlim = (-0.2,0.2), exclPeaks = exclPeaks, saveplot = true)

# comparson chi2: from least squares and from loglikelihood
chi2diff = reshape(FitPars.chi2_ls,:) .- reshape(FitPars.chi2,:);
GoodIdx = findall((!isnan).(chi2diff))
chi2diff = chi2diff[GoodIdx]
chi2_ls  = reshape(FitPars.chi2_ls,:)[GoodIdx]
chi2     = reshape(FitPars.chi2,:)[GoodIdx]

chi2_ls[findall(chi2diff .>= 1000 )] - chi2[findall(chi2diff .>= 1000 )]
stephist(filter(x-> abs(x) <=100, chi2diff), fill = true,  bins = 100, label = "least squares - loglikelihood", xlabel = "Δ χ²", ylabel = "counts", alpha = 0.5)
xlims!(-50, 50)
maximum(chi2diff)