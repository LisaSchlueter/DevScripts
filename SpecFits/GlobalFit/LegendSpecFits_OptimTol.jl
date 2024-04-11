# script to test changes in LegendSpecFits package
# test new implementation of FWHM uncertainties and p-value/gof

# plan: look at all channels in this run. do all fits and compare fit parameters and goodness of fit. 
using LegendSpecFits # make sure you're in development mode: dev Packages/LegendSpecFits.jl
using LegendDataManagement
using LegendHDF5IO, HDF5
using LegendDataTypes: fast_flatten, readdata
using Statistics, StatsBase
using Plots, ColorSchemes
using Unitful, Measures 
#using Distributions 
#using Calculus, LinearAlgebra, Zygote
#using ProgressMeter 
using LaTeXStrings
using TypedTables
using Printf
#using BAT # for specfit_combined 
#using InverseFunctions # for specfit_combined
#using Optim 
#using Interpolations
include("../utils.jl")
include("_LoadFitResults.jl")
include("../../PrettyPlotArg.jl")

# Load data 
l200 = LegendData(:l200)
ch_geds     = channelinfo(l200,(:p03,:r000,:cal),system = :geds, only_processable = true).channel;
dets        = channelinfo(l200,(:p03,:r000,:cal),system = :geds, only_processable = true).detector;
deleteat!(ch_geds, [35,41,84]) #35, 41,84 fits dont work
deleteat!(dets, [35,41,84]) #35, 41,84 fits dont work
nChannel    =  length(ch_geds)
Idx_ppc, Idx_icpc, Idx_coax, Idx_bege = get_det_types(dets) # do assignment of different detector types 
th228_lines =  [583.191,  727.330,  860.564,  1592.53,    1620.50,    2103.53,    2614.51]
g_tols = [1e-8,1e-9,1e-10,1e-11,1e-12]
pval = zeros(length(g_tols),7,nChannel)
for i=1:nChannel 
    try
        pval[:,:,i] = LoadPeakFits_g_tol(ch_geds[i],g_tols)
    catch
        pval[:,:,i] = NaN.*ones(length(g_tols),7)
        @info "file doesnt exist $(i)"
    end
end

# each channel separately. 
pval_i = repeat(reshape(pval[1,:,:],1,7,86),5,1,1) # pval with default 1e-08 tolerance 
pval_diff = pval .- pval_i # improvement of pval with decreasing g_tol

plt_arg = Dict(
    :xscale=>:log10,
    :xticks => g_tols,
)
path_plot = "PackageDevScripts/DevelopSpecFits/GlobalFit/plots"
path_rel_plot = relpath(path_plot,pwd())
path_abs_plot = pwd() * "/" * path_rel_plot * "/"

plts = []

for IdxPeak=1:length(th228_lines)
    pval_diff_p     = pval_diff[:,IdxPeak,:]
    pval_mean       = mean(pval_diff_p,dims=2)
    pval_std        = std(pval_diff_p,dims=2)
    pval_quantiles  = hcat([quantile(pval_diff_p[tol,:],[0.5,0.6827,0.95,0.99]) for tol=1:length(g_tols) ]...)


    plotlyjs()
    p = plot(g_tols,pval_quantiles[1,:],ribbon = pval_quantiles[3,:],fill = (0,:lightgray,0.8),label =:none,color=:silver;plt_arg...)
    plot!(g_tols,pval_quantiles[1,:],ribbon = pval_quantiles[2,:],fill = (0,:orange,0.5),label = :none,color=:orange;plt_arg...)
    plot!(g_tols,pval_quantiles[1,:],color=:black,linewidth=2,linestyle = :dash,label = "Median";plt_arg...)
    # make nice legend 
    plot!(g_tols,pval_quantiles[1,:], fillrange = 0, fill = (0, :orange, 0.5), line = (:dash, :transparent), label = "Quantile: 68.3% ")
    plot!(g_tols,pval_quantiles[1,:], fillrange = 0, fill = (0, :lightgray, 0.8), line = (:dash, :transparent), label = "Quantile: 95% ")
    ylabel!("∆ p-value"); xlabel!("Fit tolerance")
    plot!(p,margin=2mm,legendtitle="Th228 peak: $(round(th228_lines[IdxPeak],digits=1)) keV ";PlotFontSize()...,PlotBasic()...,PrettyLeg(;location=:topright)...)
    pname = path_abs_plot * "OptimTol_p3_r0_allchannel_peak$(string.(th228_lines[IdxPeak])[1:3])keV.pdf"
    @info "saving plot to $pname"
    savefig(pname)

    push!(plts,p)
end

plot(plts...,layout=(7,1),size=(500,1500),margin=2mm,legend=false)
#L"^{228}Th peak %$(th228_lines[IdxPeak])"

# Add a dummy series with a filled area for the custom legend symbol


# ------------ ------------ #
#  all peaks have a large fraction of p-vaulues, take 2614 keV as example (highest fraction of bad fits )
# --------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------

# --------- come information about the global vs single peak fits
LogLikeDiff = negloglike_ml_combi .-negloglike_ml
nChannel    = size(LogLikeDiff,2)
t = Table(
Th228Line       = ["all peaks",th228_lines...],
FracCombiBetter =  [round(100*sum(LogLikeDiff.<0)/(length(LogLikeDiff)),digits=1),round.(100 .*sum(LogLikeDiff.<0,dims=2)./nChannel,digits=1) ...],
FracCombiWorse  =  [round(100*sum(LogLikeDiff.>0)/(length(LogLikeDiff)),digits=1),round.(100 .*sum(LogLikeDiff.>0,dims=2)./nChannel,digits=1) ...]
)

pvaldiff = vec(pval_combi-pval)
FracCombiBetter = sum(pvaldiff.>0)./length(pvaldiff)
sum(pvaldiff.<0.01)./length(pvaldiff)

# look at fits for which global fit performed at least 10% better than single peak fit, and both fits are "good" (pval>0.05)
#select channel and peak to test g_tol
Idx_sel = findall( pval_combi-pval.>0.05 .&& pval.>0.05)
ch_geds_sel = Idx_sel[1][2]# 
peak_sel  =  Idx_sel[1][1]
pval_combi[Idx_sel[1]]
pval[Idx_sel[1]]

# load data and do simple calibration
data        = load_hitch(l200,DataPeriod(3),DataRun(0); channel=ch_geds[ch_geds_sel]);
data_energy = data.e_cusp; 
result_simple, report_simple, th228_lines, th228_names = do_simple_calibration(data_energy); # simple calibration, get histograms around the calibration lines and peakstats
peakstats = result_simple.peakstats[peak_sel]
peakhists = result_simple.peakhists[peak_sel]
result_i,report_i = LegendSpecFits.fit_single_peak_th228(peakhists,peakstats,;uncertainty=true)
@info "cross check: p = $(pval[Idx_sel[1]]) values from file, $(result_i.pval) values from live fit" # make sure its the correct peak/channel 

# save results 

result = Dict{Real, NamedTuple}()
report = Dict{Real, NamedTuple}()
for g_tol in g_tols
    result_tmp,report_tmp = LegendSpecFits.fit_single_peak_th228(peakhists,peakstats,;uncertainty=true,g_tol = g_tol)
    result[g_tol] = result_tmp
    report[g_tol] = report_tmp
end
map(x->result[x].pval,g_tols) #improvement from 0.68 to 0.75. but global fit still better (0.78)



# look at all data : single peak fits with different tols 
g_tols = [1e-8,1e-9,1e-10,1e-11,1e-12]
pval_tol = zeros(length(g_tols),7,nChannel)
for i=1:length(g_tols)
    try
    _, _ , _ ,_, _, _, _, _, pval_tol[i,:,:], _,_ , _ ,_, _, _, _ = LoadSinglePeakFits(ch_geds[20:30];g_tol = g_tols[i])
    catch
        pval_tol[i,:,:] = NaN.*ones(7,nChannel)
        @info "file doesnt exist $(i)"
    end
end



#pval_all2 = [pval_all[:,:,1:33]...,pval_all[:,:,45:end]...]

pvalimp = [ pval_all[:,i,ch] .- pval_all[1,i,ch] for i=1:7 for ch=1:nChannel]
plot(g_tols,pvalimp,xscale=:log10,ylabel="p-value improvement",xlabel="Fit tolerance",guidefontsize=14,legend=false,box=:on)#,labels=th228_lines,title="Improvement of p-value with decreasing g_tol",xlabel="g_tol",ylabel="pval - pval[1]")


plot(g_tols,pval_all[:,7,3],xscale=:log10,ylabel="p-value improvement",xlabel="Fit tolerance",guidefontsize=14,legend=false,box=:on)#,labels=th228_lines,title="Improvement of p-value with decreasing g_tol",xlabel="g_tol",ylabel="pval - pval[1]")



#plot(g_tols,mean(pvalimp),ribbon = std(),xscale=:log10,ylabel="p-value improvement",xlabel="Fit tolerance",guidefontsize=14,legend=false,box=:on)#,labels=th228_lines,title="Improvement of p-value with decreasing g_tol",xlabel="g_tol",ylabel="pval - pval[1]")

#pvals = LoadPeakFits_g_tol(ch_geds[1],g_tols)


pvalimp = [ pvals[:,i] .- pvals[1,i] for i=1:7]
plot(g_tols,pvalimp,xscale=:log10,ylabel="p-value improvement",xlabel="Fit tolerance",guidefontsize=14)#,legend=:topleft,labels=th228_lines,title="Improvement of p-value with decreasing g_tol",xlabel="g_tol",ylabel="pval - pval[1]")




path = "PackageDevScripts/DevelopSpecFits/GlobalFit/results"
fname = "single_peak_fit_p3_run0_allch_gtols.h5"
save_result(fname,["g_tols","pval_tol","ch_geds","th228_lines"],[g_tols, pval_tol,ch_geds,th228_lines],;path=path)


# read some fieldnames

