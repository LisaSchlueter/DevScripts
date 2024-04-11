
using LegendDataManagement
using JLD2 
using Plots

path = "PackageDevScripts/DevelopSpecFits/BasicPlots/results"
path_rel = relpath(path,pwd())
path_abs = pwd() * "/" * path_rel * "/"

path_plot = "PackageDevScripts/DevelopSpecFits/BasicPlots/plots"
path_rel_plot = relpath(path_plot,pwd())
path_abs_plot = pwd() * "/" * path_rel_plot * "/"

l200 = LegendData(:l200)
ch_geds     = channelinfo(l200,(:p03,:r000,:cal),system = :geds, only_processable = true).channel;
ch_dets     = channelinfo(l200,(:p03,:r000,:cal),system = :geds, only_processable = true).detector;
@load "$(path_abs)singlepeakfit_p3r0_faildets.jld2" failIdx
deleteat!(ch_geds, failIdx);deleteat!(ch_dets, failIdx)
nChannel = length(ch_geds)
fnames = path_abs .* map(x->"singlepeakfit_p3r0_det$(string(x)).jld2",ch_dets)
result_all,report_all,th228_lines = [load(fname,"result","report","th228_lines") for fname in fnames] #fname result report th228_lines th228_names result_simple report_simple 

# plot all single fits with model 
th228_lines = map(x->load(x,"th228_lines"),fnames)[1] #fname result report th228_lines th228_names result_simple report_simple 
result = map(x->load(x,"result"),fnames) #fname result report th228_lines th228_names result_simple report_simple 
pvals = [map(x->x[line].pval,result) for line in th228_lines]

p = []
for i=1:length(th228_lines)
    ptmp = histogram(pvals[i],bins=20,xlims=(0,1),fillcolor=:silver,color=:silver,normalize=:probability,
            label="$(nChannel) fits") #(string(ch_dets[chIdx])) \n (th228_lines[i]) keV \n
    if (sum(pvals[i].<0.05)/length(pvals[i]))>0.4
        ylims!(0,0.55) 
    else
        ylims!(0,0.25)
    end
    hline!([0.05],color=:red,linestyle=:dash,linewidth=2,label="Expectation: 5%")
    yl = ylims(); xlims!(0,1)
    if i==1
        #annotate!(0.05,yl[2]-0.2*(yl[2]-yl[1]),text("$(string(ch_dets[chIdx])) \n $(nChannel) fits",16,:left))
        plot!(legend = :topright,legendfontsize=16,legendtitle="$(th228_lines[i]) keV",legendtitlefontsize=16)
    else
        plot!(legend = false)
        annotate!(0.95,yl[2]-0.15*(yl[2]-yl[1]), text(@sprintf("%.2f keV \n ",th228_lines[i]),16,:right))
    end
      push!(p,ptmp)
end
plt_pval = plot(p...,layout=(2,4), size=(2200,900),bottom_margin=10mm, left_margin=12mm, top_margin=7mm,right_margin=7mm;PlotBasic()...,PlotFontSize(fontsize=20)...)
xlabel!("p-value\n")
ylabel!("Probability")
fname = path_abs_plot * "singlepeakfit_p3r0_pvaldist.pdf"
savefig(plt_pval,fname)

