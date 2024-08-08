using LegendDataManagement
using LegendDataManagement.LDMUtils
using LegendSpecFits
using LegendHDF5IO
using TypedTables
using Unitful, Measures
using Measurements: value as mvalue, uncertainty as muncert
using StructArrays, PropDicts
using Revise
using Plots
using Distributions 
using ValueShapes
using LaTeXStrings
include("./utils_ecal.jl")
# data selection 
l200 = LegendData(:l200)
runsel = (DataPeriod(3), DataRun(0), :cal)
chinfo = channelinfo(l200, runsel; system=:geds, only_processable=true)
dets_ged  = chinfo.detector
det_types = reverse(unique(detector_type.(Ref(l200), dets_ged)))
peakshape = :f_fit_bckExp#Slope


th228_names = []
for det_type in det_types
    dets_plt = findall(x->x==det_type, detector_type.(Ref(l200), dets_ged))[1:5]

    for det in dets_ged[dets_plt]
        println("Type: ", detector_type(l200, det), ", detector: ", det)
        ## --------------------------------- preparation for fit start --------------------------------- ##
        # get energy spectrum in ADC
        """
            get_energy_ctc(data::LegendData, runsel::Tuple{DataPeriod, DataRun, Symbol}, detector::DetectorId; e_type::Symbol=:e_cusp)
        # open hitchannel file and apply charge-trapping correction
        """


        try
            local result_simple, _, th228_names, _ = simple_cal(l200, runsel, det)
         catch e
             @info "hitch file not found or other problem. $e"
             continue 
         end
         if !@isdefined(result_simple)
            continue
         end 
        ## --------------------------------- preparation for fit end--------------------------------- ##   
        # fit 1: baselinefit with default peakshape
        result, report = fit_peaks(result_simple.peakhists, result_simple.peakstats, th228_names; e_unit=result_simple.unit, calib_type=:th228, fit_func = :f_fit);

        # fit 2: choose different peakshape 
        # pseudo_prior_add = NamedTupleDist(background_slope = Normal(0.0, 0.01))
        result_alt, report_alt = fit_peaks(result_simple.peakhists, result_simple.peakstats, th228_names; e_unit=result_simple.unit, calib_type=:th228, fit_func = peakshape);#, pseudo_prior = pseudo_prior_add)

        # plots 
        p_def = plot(broadcast(k -> plot(report[k], left_margin=20mm, top_margin=-5mm, bottom_margin=-2mm, title=string(k), ms=2), th228_names)..., layout=(length(report), 1), size=(1000,710*length(report)) , thickness_scaling=1.8, titlefontsize = 10, legendfontsize = 8, yguidefontsize = 9, xguidefontsize=11)
        p_alt = plot(broadcast(k -> plot(report_alt[k], left_margin=20mm, top_margin=-5mm, bottom_margin=-2mm, title=string(k), ms=2), th228_names)..., layout=(length(report_alt), 1), size=(1000,710*length(report_alt)) , thickness_scaling=1.8, titlefontsize = 10, legendfontsize = 8, yguidefontsize = 9, xguidefontsize=11)
        pall = plot(p_def, p_alt, layout=(1,2), size=(2500, 760*length(report_alt)), thickness_scaling=2, titlefontsize = 12, 
                legendfontsize = 8, yguidefontsize = 9, xguidefontsize=12, dpi = 300, 
                plot_title = " $(runsel[1]), $(runsel[2]) $det ($(detector_type(l200,det))) \n left: hypermet  | right: hypermet + $(replace(string(peakshape),"f_fit_" => ""))")

        # plot(report_fit[th228_names[4]], legend = :topleft, thickness_scaling=1.1, title =  "$(runsel[1]), $(runsel[2]), $det: " * string(th228_names[4]), titlefontsize = 12, ms = 3)
        pltpath = "$(@__DIR__)/plots/specfit/"
        if !isdir(pltpath)
            mkdir(pltpath)
        end
        figname = pltpath *  "specfit_hypermet_vs_$(replace(string(peakshape),"f_fit_" => ""))_$(runsel[1])_$(runsel[2])_$(det_type)_$(det).png"
        savefig(figname)
    end
end

# compare fit parameter for these fits 
function compare_fitpar(result1, result2)
    Table(
         µ = [result1.μ, result2.μ, mvalue(result1.μ -  result2.μ)], 
         σ = [result1.σ, result2.σ, mvalue(result1.σ -  result2.σ)],
         n = [result1.n, result2.n, mvalue(result1.n -  result2.n)], 
         step_amplitude = [result1.step_amplitude,result2.step_amplitude, mvalue(result1.step_amplitude -  result2.step_amplitude)], 
         skew_fraction = [result1.skew_fraction, result2.skew_fraction, mvalue(result1.skew_fraction -  result2.skew_fraction)],
         skew_fraction_highE = [0.0, result2.skew_fraction_highE, NaN], 
         skew_width = [result1.skew_width, result2.skew_width, mvalue(result1.skew_width -  result2.skew_width)],
         skew_width_highE = [0.0, result2.skew_width_highE, NaN], 
         background = [result1.background, result2.background, mvalue(result1.background -  result2.background)],
        # background_slope = [0.0, result2.background_slope, NaN],
         chi2 = [result1.gof.chi2, result2.gof.chi2, result1.gof.chi2 - result2.gof.chi2],
         pval = [result1.gof.pvalue, result2.gof.pvalue, result1.gof.pvalue - result2.gof.pvalue]
         )
 end

 for peak in eachindex(th228_names)
    println("\n\nPeak: ", th228_names[peak])
    local tbl = compare_fitpar(result[th228_names[peak]], result_alt[th228_names[peak]])
    display(tbl)
end



