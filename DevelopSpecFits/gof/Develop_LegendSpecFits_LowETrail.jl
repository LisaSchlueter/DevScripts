# script to test changes in LegendSpecFits package
# test new implementation of FWHM uncertainties and p-value/gof

using LegendSpecFits # make sure you're in development mode: dev Packages/LegendSpecFits.jl
using LegendDataManagement
using LegendHDF5IO, HDF5
using LegendDataTypes: fast_flatten, readdata
using Statistics, StatsBase
using Plots, ColorSchemes
using Distributions 
using Calculus, LinearAlgebra, Zygote
using ProgressMeter 
using TypedTables
using Printf
using ValueShapes
using ForwardDiff
include("utils.jl")

# get data path and filename
l200 = LegendData(:l200)
ch_geds     = channelinfo(l200,(:p03,:r000,:cal),system = :geds, only_processable = true).channel;
data        = load_hitch(l200,DataPeriod(3),DataRun(0); channel=ch_geds[2]);
data_energy = data.e_cusp; 

# ------------ do "simple" fit --> input for actual fit  ------------ #
result_simple, report_simple, th228_lines, th228_names = do_simple_calibration(data_energy); # simple calibration, get histograms around the calibration lines and peakstats

# cross-check; do fit iterative fits of gamma peaks 
#result, report = fit_peaks(result_simple.peakhists, result_simple.peakstats, th228_lines,; uncertainty=true,iterative_fit=false);
Idx = 1
fixed_position = false
ps = result_simple.peakstats[Idx]
h = result_simple.peakhists[Idx]

result, _          = fit_single_peak_th228(h,ps, ; uncertainty=true, low_e_tail=true,fit_fun=:f_fit)
result_tailflase,_ = fit_single_peak_th228(h,ps, ; uncertainty=true, low_e_tail=false,fit_fun=:f_fit) # why does skew parameters still have uncertainty?
    
#result_tailflase2,_ = fit_single_peak_th228(h,ps, ; uncertainty=true, low_e_tail=false,fit_fun=:f_sigWithbackground)


