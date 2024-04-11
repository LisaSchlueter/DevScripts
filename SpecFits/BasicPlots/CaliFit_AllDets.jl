# understand simple_calibration and do some plots along the way
using LegendDataManagement
using LegendSpecFits
using HDF5, LegendHDF5IO
using TypedTables
using Plots
using Measures
using JLD2 # saving files 
using Printf
include("../../PrettyPlotArg.jl")
include("../utils.jl")

path = "PackageDevScripts/DevelopSpecFits/BasicPlots/results"
path_rel = relpath(path,pwd())
path_abs = pwd() * "/" * path_rel * "/"

path_plt = "PackageDevScripts/DevelopSpecFits/BasicPlots/plots"
path_rel_plt = relpath(path_plt,pwd())
path_abs_plt = pwd() * "/" * path_rel_plt * "/"

l200 = LegendData(:l200)
ch_geds     = channelinfo(l200,(:p03,:r000,:cal),system = :geds, only_processable = true).channel;
ch_dets     = channelinfo(l200,(:p03,:r000,:cal),system = :geds, only_processable = true).detector;
nChannel = length(ch_geds)

function PlotBf(result,report,chIdx;residuals::Bool=false)
    p = []
    for line in th228_lines
        ptmp = plot(report[line])
        yl = ylims(); xl = xlims()
        if isnan(result[line].gof.pvalue)
            nbins = length(result[line].bin_centers)
            annotate!(xl[1]+0.05*(xl[2]-xl[1]),yl[2]-0.5*(yl[2]-yl[1]), text(@sprintf("p = N/A, because \n%.0f bins and %.0f fit parameter ",nbins,nbins-result[line].gof.dof),16,:left))
        else
             annotate!(xl[1]+0.05*(xl[2]-xl[1]),yl[2]-0.45*(yl[2]-yl[1]), text(@sprintf("p = %.2g",result[line].gof.pvalue),16,:left))
        end
        if line==th228_lines[1]
            annotate!(xl[2]-0.05*(xl[2]-xl[1]),yl[2]-0.45*(yl[2]-yl[1]),text("$(string(ch_dets[chIdx]))",16,:right))
            plot!(legend = :bottomright,legendfontsize=10)
        else
            plot!(legend = false)
        end
    
        if residuals==true
            ylabel!("\nCounts")
            xlabel!("") # remove xlabel from upper plot
            gof = result[line].gof
            pres = plot(range(xl[1],xl[2],10),zeros(10),ribbon = 3 .*ones(10),fillalpha=0.5,fillcolor=:lightgrey,linecolor=:darkgrey,legend=false)
            plot!(range(xl[1],xl[2],10),zeros(10),ribbon = ones(10),fillalpha=1,fillcolor=:lightgrey,linecolor=:darkgrey,legend=false)
            plot!(gof.bin_centers,gof.residuals_norm,size=(500,100),marker = :dot,markersize=1.5,color=:red,markeredgecolor = :transparent,linecolor = :transparent,legend=false)
           # 
           if maximum(abs.(gof.residuals_norm))<5
                ylims!(-5,5)
                yticks!([-3,0,3])
           else 
             yext = maximum(abs.(gof.residuals_norm))
             ylims!(-yext,yext)
             #yticks!(-ceil(Int,yext):3:ceil(Int,yext))
           end
            xlabel!("Energy (keV)")
            ylabel!("Residuals (σ)")
            layout = @layout([a; b{0.3h}])
            ptmp2 = plot(ptmp,pres,layout=layout,size=(600,400),top_margin=-10mm;PlotBasic()...,PlotFontSize(fontsize=10)...)
            push!(p,ptmp2)
        else
            push!(p,ptmp)
        end
       
    end
   
     # if residuals==false
         if residuals==false
            plt_bf = plot(p...,layout=(2,4), size=(2200,900),bottom_margin=10mm, left_margin=12mm, top_margin=7mm,right_margin=7mm;PlotBasic()...,PlotFontSize(fontsize=18)...)
            xlabel!("Energy (keV)\n")
            ylabel!("\nCounts")
         else 
             plt_bf = plot(p...,layout=(2,4), size=(2200,1100),bottom_margin=9mm, left_margin=12mm, top_margin=3mm,right_margin=7mm;PlotBasic()...,PlotFontSize(fontsize=18)...)
        end
   # else
   #     plt_bf = plot(p...,layout=(2,4), size=(2200,900),bottom_margin=10mm, left_margin=12mm, top_margin=7mm,right_margin=7mm;PlotBasic()...,PlotFontSize(fontsize=18)...)
       
   # end

    return plt_bf 
end


# ---- test 
chIdx = 2
data_det        = load_hitch(l200,DataPeriod(3),DataRun(0); channel=ch_geds[chIdx]);
result_simple, report_simple, th228_lines, th228_names = do_simple_calibration([data_det.e_cusp...]); # simple calibration, get histograms around the calibration lines and peakstats
result,report = LegendSpecFits.fit_peaks_th228(result_simple.peakhists,result_simple.peakstats,th228_lines) # do to: put this
plt_bf = PlotBf(result,report,chIdx)

p = plot(report[th228_lines[1]])
xl = xlims()
xlabel!("")
gof = result[th228_lines[1]].gof
pband = plot(range(xl[1],xl[2],10),zeros(10),ribbon = 3 .*ones(10),fillalpha=0.5,fillcolor=:lightgrey,linecolor=:darkgrey,legend=false)
pband = plot!(range(xl[1],xl[2],10),zeros(10),ribbon = ones(10),fillalpha=1,fillcolor=:lightgrey,linecolor=:darkgrey,legend=false)
pres = plot!(gof.bin_centers,gof.residuals_norm,size=(500,100),marker = :dot,markersize=3,color=:red,markeredgecolor = :transparent,linecolor = :transparent,legend=false)
ylims!(-5,5)
xlabel!("Energy (keV)")
ylabel!("Residuals (σ)")

# Create a custom layout where the first subplot is twice as large as the second
p_all = plot(p,pres,layout=layout,size=(600,400),top_margin=-2mm;PlotBasic()...,PlotFontSize(fontsize=10)...)
p_all2 = plot(p,pres,layout=layout,size=(600,400),top_margin=-3mm;PlotBasic()...,PlotFontSize(fontsize=10)...)

xlims!(xl)
# ---test ende 

failIdx = [] 
#for chIdx=85#nChannel
chIdx=84#84
    @info "start idx $chIdx $(ch_dets[chIdx])"
  #  if !isfile(path_abs_plt * "singlepeakfit_p3r0_det$(string(ch_dets[chIdx])).pdf")
        #try
                data_det        = load_hitch(l200,DataPeriod(3),DataRun(0); channel=ch_geds[chIdx]);
                result_simple, report_simple, th228_lines, th228_names = do_simple_calibration([data_det.e_cusp...]); # simple calibration, get histograms around the calibration lines and peakstats
                result,report = LegendSpecFits.fit_peaks_th228(result_simple.peakhists,result_simple.peakstats,th228_lines) # do to: put this
                plt_bf = PlotBf(result,report,chIdx;residuals=true)
                fname_plt = path_abs_plt * "singlepeakfit_p3r0_det$(string(ch_dets[chIdx])).pdf"
                savefig(plt_bf,fname_plt)
                fname = path_abs * "singlepeakfit_p3r0_det$(string(ch_dets[chIdx])).jld2"
                @save fname result plt_bf th228_lines th228_names result_simple
       # catch err
       #     @info "$err"
      #      push!(failIdx,chIdx)
       # end
  #  end
#end
