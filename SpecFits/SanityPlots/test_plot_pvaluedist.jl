using LegendDataManagement

l200 = LegendData(:l200)
include("utils.jl")
include("plot_pvaluedist.jl")


period = DataPeriod(3)
run = DataRun(2)
e_types = collect(keys(l200.par.rpars.ecal.p03.r000.P00664A))
e_type = e_types[7]
# ----------- runwise ----------------- #
runs = search_disk(DataRun,l200.tier[:jldsp,:cal,period])
for run in runs
pvalues, th228_lines  = _get_pval_runwise(period,run,e_type;mode=:cal)
plt_pval = plot_pvaluedist(pvalues,th228_lines;fontsize = fontsize=37)
plot!(plt_pval[end],legendtitle= "period $(parse(Int,(match(r"\d+",string(period)).match))), run  $(parse(Int,(match(r"\d+",string(run)).match))), \n $(string(e_type))",legendtitlefontsize=37)
fname = _get_pvaldist_pltname(period,run,e_type;format="png")
savefig(plt_pval,fname)
end
# ----------- periodwise ----------------- #
periods_all = search_disk(DataPeriod,l200.tier[:jldsp,:cal])
for  e_type in e_types
    for period in periods_all[1:end]
        @info "process period $(period)"
    
        try
            pvalues, th228_lines,runs_period  = _get_pval_periodwise(period,e_type;mode=:cal);
        catch
            continue
        end

        pvalues_all = Dict()
        for line in th228_lines
            pvalues_all[line] = filter(x->isa(x,Float64),reduce(vcat,[pvalues[run][line] for run=1:length(runs_period)] )) # sometime "PropDicts.MissingProperty" --> filter these out 
        end

        plt_pval = plot_pvaluedist(pvalues_all,th228_lines;fontsize = fontsize=45)
        plot!(plt_pval[end],legendtitle= "Period $(parse(Int,(match(r"\d+",string(period)).match))), $(string(e_type))",legendtitlefontsize=37)

        fname = _get_pvaldist_pltname(period,e_type;format="png")
        savefig(plt_pval,fname)
        @info "save fig to $(fname)"
    end
end 


