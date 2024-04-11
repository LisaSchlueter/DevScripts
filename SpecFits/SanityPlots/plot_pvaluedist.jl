using LegendDataManagement
using Printf
using Plots
using Measures

l200 = LegendData(:l200)

function plot_pvaluedist(pvalues::Dict,th228_lines::Vector{Symbol};fontsize::Int=37)
    p = []
    for (i,line) in enumerate(th228_lines)
        ptmp = histogram(pvalues[line],bins=20,xlims=(0,1),
                        fillcolor=:silver,color=:silver,normalize=:probability,
                        framestyle = :box,legend = false) #(string(ch_dets[chIdx])) \n (th228_lines[i]) keV \n
        if (sum(pvalues[line].<0.05)/length(pvalues[line]))>0.4 && (sum(pvalues[line].<0.05)/length(pvalues[line]))<0.55
            ylims!(0,0.55) 
        elseif (sum(pvalues[line].<0.05)/length(pvalues[line]))>0.4 && (sum(pvalues[line].<0.05)/length(pvalues[line]))<0.25
            ylims!(0,0.25)
        else 
          ylims!(0,maximum(0.05+sum(pvalues[line].<0.05)/length(pvalues[line])))
        end
        hline!([0.05],color=:red,linestyle=:dash,linewidth=3)
        yl = ylims(); xlims!(0,1)
      
       
        annotate!(0.95,yl[2]-0.1*(yl[2]-yl[1]), text("$(string(line))",fontsize-9,:right))

          push!(p,ptmp)
          if line==th228_lines[end]
            ptmp2 = histogram(pvalues[line],bins=20,xlims=(0,1),fillcolor=:silver,normalize=:probability,
            label=" $(length(pvalues[line])) fits")
            hline!([0.05],color=:red,linestyle=:dash,linewidth=3,label=" Expectation: 5%")
            push!(p,ptmp2)
          end
    end

    plt_pval = plot(p...,
                    layout=(2,4),size=(4400,1500),
                    xlabel = "p-value\n",
                    ylabel = "\nProbability",
                    bottom_margin=15mm, left_margin=20mm, top_margin=7mm,right_margin=10mm,
                    grid = false,  
                    tickfontsize = fontsize-12,
                    xguidefont=font(halign=:center, pointsize=fontsize),
                    yguidefont=font(halign=:center, pointsize=fontsize), #family="monospace"
                    dpi = 300)

   # make legend nice 
    plot!(plt_pval[end],showaxis=false,xlabel="",ylabel="",
     legend =:left,legendfontsize=fontsize-6)
     plot!(plt_pval[end],range(0,2,10), ones(10), ribbon = ones(10),
      fillalpha=1,fillcolor=:white,linecolor=:white,label=false)
 
    return plt_pval 
end

#=
plt_pval = plot_pvaluedist(pvalues,th228_lines;fontsize = fontsize=45)

# save 
path_plot = "PackageDevScripts/DevelopSpecFits/SanityPlots/plots"
path_rel_plot = relpath(path_plot,pwd())
path_abs_plot = pwd() * "/" * path_rel_plot * "/"
fname = path_abs_plot * "singlepeakfit_p3r0_pvaldist.png"
savefig(plt_pval,fname)
@info "save fig to $(fname)"
=#




