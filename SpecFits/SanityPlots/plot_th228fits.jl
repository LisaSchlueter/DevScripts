using LegendSpecFits
using LegendDataManagement
using Plots
using Measures
using Printf

"""
    plot_th228fits(result::Dict,report::Dict,th228_names::Vector{String};residuals::Bool=true)
plot all fits of the th228 lines with the model and residuals
"""
function plot_th228fits(result::Dict,report::Dict,th228_names::Vector{String};residuals::Bool=true,fontsize::Int=37)
    p = []
    th228_lines = collect(keys(result))
    for (i,line) in enumerate(th228_lines)
        ptmp = plot(report[line],linewidth=2,legend = false,framestyle = :box,grid=false)
        yl = ylims(); xl = xlims()
        if isnan(result[line].gof.pvalue)
            nbins = length(result[line].bin_centers)
            annotate!(xl[2]-0.05*(xl[2]-xl[1]),yl[2]-0.5*(yl[2]-yl[1]), text(@sprintf("p = N/A, because \n%.0f bins and %.0f fit parameter ",nbins,nbins-result[line].gof.dof),fontsize-2,:right))
        else
             annotate!(xl[2]-0.05*(xl[2]-xl[1]),yl[2]-0.5*(yl[2]-yl[1]), text(@sprintf("p = %.2g",result[line].gof.pvalue),fontsize-2,:right))
        end
        annotate!(xl[1]+0.05*(xl[2]-xl[1]),yl[2]-0.5*(yl[2]-yl[1]),text("$(th228_names[i])",fontsize-2,:left))
       
        if residuals==true
            ylabel!("\n\nCounts")
            xlabel!("")
            gof = result[line].gof
            pres = plot(range(xl[1],xl[2],10),zeros(10),ribbon = 3 .*ones(10),fillalpha=0.5,fillcolor=:lightgrey,linecolor=:darkgrey,legend=false)
            plot!(range(xl[1],xl[2],10),zeros(10),ribbon = ones(10),fillalpha=0.5,fillcolor=:grey,linecolor=:darkgrey,legend=false)
            plot!(gof.bin_centers,gof.residuals_norm,size=(500,100),marker = :dot,markersize=2,color=:red,markeredgecolor = :transparent,linecolor = :transparent,legend=false)
           
           if maximum(abs.(gof.residuals_norm))<5
                ylims!(-5,5)
                yticks!([-3,0,3])
           else
                yext = maximum(abs.(gof.residuals_norm))
                ylims!(-yext,yext)
                if yext<10
                    yticks!([-5,0,5])
                else 
                  yticks!([-floor(Int,yext)+2,0,floor(Int,yext)-2])
                end
            end
            xlabel!("Energy (keV)")
            ylabel!("\n\nResiduals (σ)")
            layout = @layout([a; b{0.3h}])
            ptmp2 = plot(ptmp,pres,
                    layout=layout,
                    framestyle = :box,grid=false)
            push!(p,ptmp2)
            if line==th228_lines[end]
                pleg = plot(report[line])
                push!(p,pleg)
            end
        else
            push!(p,ptmp)
            if line==th228_lines[end]
                pleg = plot(report[line])
                push!(p,pleg)
            end  
        end
    end
   
    if residuals==true
        plt_size = (4400, 2300)
    else
        plt_size = (4400, 1900)
    end

    plt_bf = plot(p...,
                layout=(2,4), 
                size=plt_size,
                left_margin=30mm, right_margin=10mm,
                bottom_margin=20mm, top_margin=20mm, 
                xguidefont=font(halign=:center, pointsize=fontsize),#family="monospace"
                yguidefont=font(halign=:center, pointsize=fontsize),#family="monospace",
                xtickfontsize=fontsize-9,
                ytickfontsize=fontsize-9,
                dpi=300)

    if residuals==true

        
        #residual plots, space between subplots 
        for i=1:length(plt_bf)
            if mod(i,2)==0 
            plot!(plt_bf[i],top_margin=-10mm)
            end
         end
    else 
        xlabel!("Energy (keV)\n")
        ylabel!("Counts")     
    end

    # make legend nice 
    plot!(plt_bf[end],showaxis=false,xlabel="",ylabel="",legend =:left,legendfontsize=fontsize)
    xl = xlims(plt_bf[end]); ylmax = maximum(abs.(ylims(plt_bf[end])))
    plot!(plt_bf[end],range(xl[1],xl[2],10),(ylmax+1) .*ones(10),ribbon = (ylmax) .*ones(10),
        fillalpha=1,fillcolor=:white,linecolor=:white,label=false)

    # 2600 line, less xticks 
    xt = xticks(plt_bf[end-1])
    @info "$(length(xt[1]))"
    plot!(plt_bf[end-1],xticks=(xt[1][1:2:end]))
    plot!(plt_bf[end-2],xticks=(xt[1][1:2:end]))
    return plt_bf 
end


