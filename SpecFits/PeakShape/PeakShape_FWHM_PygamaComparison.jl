# understand simple_calibration and do some plots along the way
using LegendDataManagement
using LegendDataManagement.LDMUtils
using LegendSpecFits
using LegendHDF5IO, HDF5
using Measures
using Plots
using PDFmerger, Printf, LaTeXStrings, ColorSchemes
using TypedTables
using Unitful
using Measurements
using Measurements: value as mvalue, uncertainty as muncert
using StatsBase
using Distributions, ValueShapes, BAT
include("../utils.jl")

path_plot = "$(@__DIR__)/plots/"
if !ispath(path_plot)
    mkdir("$path_plot")
end


# select and load data from file 
periods_comparison = [3,3,3,3,3,3]
runs_comparison = [0,1,1,1,2,4]
det_name_comparison = ["V04545A","V01240A", "V04545A", "V05267A", "V04545A", "V04545A"]
# compareIdx = 1
for compareIdx in 1:6
    period = periods_comparison[compareIdx]
    run = runs_comparison[compareIdx]
    ch_name = det_name_comparison[compareIdx]
    @info "Comparison $compareIdx: $period $run $ch_name"
    l200 = LegendData(:l200)
    ch_geds     = channelinfo(l200,(DataPeriod(period),DataRun(run),:cal),system = :geds, only_processable = true).channel;
    dets_geds     = channelinfo(l200,(DataPeriod(period),DataRun(run),:cal),system = :geds, only_processable = true).detector;
    filekey = start_filekey(l200, (DataPeriod(period), DataRun(run), :cal)) 
    nChannel = length(ch_geds)
    chIdx = findall(string.(dets_geds) .== ch_name)[1]
    #data_det = load_hitchfile(lh5open, l200, (DataPeriod(period), DataRun(run), :cal), ch_geds[chIdx]) # using LDMUtils doesnt work 
    data_det    = load_hitch(l200,DataPeriod(3),DataRun(0); channel= ch_geds[chIdx]);

    # load and apply  charge trapping correction 
    @debug "Loaded CTC parameters"
    pars_ctc = get_values(l200.par.rpars.ctc[DataPeriod(period), DataRun(run)])[dets_geds[chIdx]].e_cusp
    e_cusp_uncal_ctc = ljl_propfunc(pars_ctc.func).(data_det) # uncalibrated, but charge trapping corrected energies in ADC

    # simple calibration 
    result_simple, report_simple, th228_lines, th228_names = do_simple_calibration(e_cusp_uncal_ctc); # simple calibration, get histograms around the calibration lines and peakstats

    # fit 1 peaks and plot result. 
    peakIdx = 7
    ps = result_simple.peakstats[peakIdx]
    @info "fitting peak $(th228_lines[peakIdx]) ($(th228_names[peakIdx]))"

    #use standard prior (May 30, 2024)
    standard_pseudo_prior = LegendSpecFits.NamedTupleDist(
        μ = Uniform(ps.peak_pos-10, ps.peak_pos+10),
        σ = weibull_from_mx(ps.peak_sigma, 2*ps.peak_sigma),
        n = weibull_from_mx(ps.peak_counts, 2*ps.peak_counts),
        step_amplitude = weibull_from_mx(ps.mean_background_step, ps.mean_background_step + 5*ps.mean_background_std),
        skew_fraction = truncated(weibull_from_mx(0.01, 0.05), 0.0, 0.1),
        skew_width = weibull_from_mx(0.001, 1e-2),
        background = weibull_from_mx(ps.mean_background, ps.mean_background + 5*ps.mean_background_std),)
    result_st, report_st = LegendSpecFits.fit_single_peak_th228(result_simple.peakhists[peakIdx], ps; pseudo_prior = standard_pseudo_prior, uncertainty = true) 
    plt_tmp = plot(report_st, title = "standard pseudo-prior,  p$period-r$run, $(dets_geds[chIdx])", titlefontsize = 12, guidefontsize = 14, bottom_margin = 1mm)
    plt_reg = plot(plt_tmp[1], legend = :outerright, size = (800,350),left_margin = 5mm, top_margin = 3mm)
    plot!(1,1,color = :white, label = "µ = $(result_st.µ) keV")
    plot!(1,1,color = :white, label = "fwhm = $(result_st.fwhm) keV")
    plot!(1,1,color = :white, label = "skew_frac = $(result_st.skew_fraction)")
    plot!(1,1,color = :white, label = "skew_width = $(round(result_st.skew_width,digits = 3))")   
    savefig(plt_reg, path_plot .* "tmp.pdf")
    fname = path_plot * "singlepeakfit__e_cusp_ctc_p$(period)r$(run)_det$(string(dets_geds[chIdx])).pdf"
    append_pdf!(fname, path_plot .* "tmp.pdf", cleanup = true)

    #use different priors
    pseudo_prior1 = LegendSpecFits.NamedTupleDist(
            μ = Uniform(ps.peak_pos-ps.peak_sigma, ps.peak_pos+ps.peak_sigma),
            σ = weibull_from_mx(ps.peak_sigma, 2*ps.peak_sigma),
            n = weibull_from_mx(ps.peak_counts, 2*ps.peak_counts),
            step_amplitude = weibull_from_mx(ps.mean_background_step, ps.mean_background_step + 5*ps.mean_background_std),
            skew_fraction = truncated(weibull_from_mx(0.005, 0.02), 0.0, 0.1),
            skew_width = weibull_from_mx(ps.peak_sigma/ps.peak_pos, 2*ps.peak_sigma/ps.peak_pos),
            background = weibull_from_mx(ps.mean_background, ps.mean_background + 5*ps.mean_background_std))

    @info "fitting peak $(th228_lines[peakIdx]) ($(th228_names[peakIdx]))"
    result_1, report_1 = LegendSpecFits.fit_single_peak_th228(result_simple.peakhists[peakIdx], ps; pseudo_prior = pseudo_prior1, uncertainty = true) 
    plt_tmp2 = plot(report_1, title = "new pseudo-prior,  p$period-r$run, $(dets_geds[chIdx])", titlefontsize = 12, guidefontsize = 14, bottom_margin = 1mm)
    plt_1 = plot(plt_tmp2[1], legend = :outerright, size = (800,350),left_margin = 5mm, top_margin = 3mm)
    plot!(1,1,color = :white, label = "µ = $(result_1.µ) keV")
    plot!(1,1,color = :white, label = "fwhm = $(result_1.fwhm) keV")
    plot!(1,1,color = :white, label = "skew_frac = $(result_1.skew_fraction)")
    plot!(1,1,color = :white, label = "skew_width = $(round(result_1.skew_width,digits = 5))")                         
    savefig(plt_1, path_plot .* "tmp.pdf")
    fname = path_plot * "singlepeakfit__e_cusp_ctc_p$(period)r$(run)_det$(string(dets_geds[chIdx])).pdf"
    append_pdf!(fname, path_plot .* "tmp.pdf", cleanup = true)    

end