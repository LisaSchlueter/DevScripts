# write report for "bad fits" pvalues 
using LegendDataManagement
using Dates
using TypedTables
using SplitApplyCombine: product, innerjoin,merge
l200 = LegendData(:l200)
include("utils.jl")
include("plot_pvaluedist.jl")
include("../../PrettyPlotArg.jl")
function _get_pval_report(period,e_type)
    # gather pvalues
    @info "load pvalues: period $(period), $(e_type)"
    pvalues, th228_lines,runs_period  = _get_pval_periodwise(period,e_type;mode=:cal);
    pvalues_all = Dict()
    for line in th228_lines
    pvalues_all[line] = filter(x->isa(x,Float64),reduce(vcat,[pvalues[run][line] for run=1:length(runs_period)] )) # sometime "PropDicts.MissingProperty" --> filter these out 
    end

    # make structarray for report input 
    peaks = ["global",string.(th228_lines)...]
    period_summary = [sum(pvalues_all[line] .<0.05)/length(pvalues_all[line]) for line in th228_lines]
    period_summary = [sum(period_summary)/7,period_summary...]
    run_summary = Vector{Any}(undef, length(runs_period))
    for i=1:length(runs_period)
        run_summary[i] = [sum(pvalues[i][line][isa.(pvalues[i][line],Float64)] .<0.05)/length(pvalues[i][line][isa.(pvalues[i][line],Float64)]) for line in th228_lines]
        run_summary[i] = [sum(run_summary[i])/length(run_summary[i]),run_summary[i]...]
    end
    pval_report = StructArray((peaks = peaks, period = round.(period_summary, digits=2), (Symbol("run", i) => round.(run_summary[i],digits=2) for i in 1:length(run_summary))...))
return pval_report, th228_lines
end

periods_all = search_disk(DataPeriod,l200.tier[:jldsp,:cal])

for period in periods_all[1:end]
#period = DataPeriod(4)
    e_types = collect(keys(l200.par.rpars.ecal.p03.r000.P00664A))

    #period = periods_all[1];


    report = lreport()
    lreport!(report, "# Calibration fits")
    lreport!(report,"## Goodness-of-fit  : Period $(parse(Int,match(r"\d+",string(period)).match))")
    lreport!(report, "Date of processing: $(now())")

    pval_reports = Vector{Any}(undef, length(e_types))
    pval_reports[i], th228_lines = _get_pval_report(period,e_types[i])
    for i = 1:length(e_types)
        try
            pval_reports[i], th228_lines = _get_pval_report(period,e_types[i])
            lreport!(report,"### $(e_types[i]): Fraction of fits with p<0.05:",pval_reports[i])
        catch 
            pval_reports[i] = 1
            deleteat!(pval_reports,i)
            deleteat!(e_types,i)
        end
    end

    # plot 
    bw =  0.7/length(e_types)
    p = bar(size=(900,500),framestyle = :box,grid=true)
    for i=1:length(e_types)
    bin_center = (1:8) .- 0.35 .+ (i-1).*bw 
    bar!(p,bin_center,pval_reports[i].period,label = string(e_types[i]),
                    ylabel="Fraction with  p-value < 0.05",
                    xlabel= "Peak fit", xticks=(1:8,["global",string.(th228_lines)...]),
                bar_width = bw*0.9,linecolor=:transparent,fillalpha = 0.7)
    end
    ylims!(0,1)
    fs = 15
    plot!(p,legend=:topright,legendfontsize=fs,legendtitle="Energy filter",legendtitlefontsize=20,
    left_margin=7mm,bottom_margin=7mm;
        PlotFontSize(fontsize=15)...,PrettyLeg(fontsize=15,location=:top)...)
    pltname =  _get_pvaldist_pltname(period)
    savefig(pltname)
    lreport!(report,p)
    # write report 


    #  @info "Write log report"
    report_name = _get_pvaldist_rptname(period)
    writelreport(report_name, report)
end