# compare global fit with single peak fits
# plot script to fits in "compare_global_fit"
# no shared fit parameters
# look at likelihoods, compare minima 

# plan: look at all channels in this run. do all fits and compare fit parameters and goodness of fit. 
using LegendSpecFits # make sure you're in development mode: dev Packages/LegendSpecFits.jl
using LegendDataManagement
using LegendHDF5IO, HDF5
using LegendDataTypes: fast_flatten, readdata
using Statistics, StatsBase
using Plots, ColorSchemes
using Distributions
using ProgressBars
using Printf
using StatsPlots
using LaTeXStrings
using TypedTables
include("../utils.jl")
include("_LoadFitResults.jl")
l200 = LegendData(:l200)
ch_geds     = channelinfo(l200,(:p03,:r000,:cal),system = :geds, only_processable = true).channel;
dets    = channelinfo(l200,(:p03,:r000,:cal),system = :geds, only_processable = true).detector;
deleteat!(ch_geds, [35,41,84]) #35, 41,84 fits dont work
deleteat!(dets, [35,41,84]) #35, 41,84 fits dont work

μ, σ , n ,step_amplitude, skew_fraction, skew_width, background, negloglike_ml, pval, µ_err,σ_err , n_err ,step_amplitude_err, skew_fraction_err, skew_width_err, background_err = LoadSinglePeakFits(ch_geds)
μ_combi,σ_combi, n_combi ,step_amplitude_combi, skew_fraction_combi, skew_width_combi, background_combi,negloglike_ml_combi, pval_combi  = LoadCombiPeakFits(ch_geds)

# ------------------------------------------------------------
# -------------------- PLOTS --------------------------------
path_plot = "PackageDevScripts/DevelopSpecFits/GlobalFit/plots"
path_rel_plot = relpath(path_plot,pwd())
path_abs_plot = pwd() * "/" * path_rel_plot * "/"
plotarg = (grid=false, box = :on)

# look how many fits are "bad" per th228 lines 
th228_lines = [583.191, 727.33, 860.564, 1592.53, 1620.5, 2103.53, 2614.51]
BadFitsFrac = [sum(pval[i,:].<0.05)/length(pval[i,:]) for i=1:7]
BadFitsFrac_Combi = [sum(pval_combi[i,:].<0.05)/length(pval_combi[i,:]) for i=1:7]
xticklabels = string.(th228_lines)
bar(0.75:6.75,BadFitsFrac,bar_width = 0.25,xticks=(0.86:6.86,xticklabels),label="Single peak fits",color=:silver,alpha = 1,guidefontsize=12,legendfontsize=10)
bar!(1:7,bar_width = 0.25,BadFitsFrac_Combi,xlabel = "Th228 peak(keV)",ylabel="Fraction of fits with p < 5%",label="Combined fit",color=:dodgerblue,alpha = 1,box = :on,grid=false,ylims=(0,0.61))
pname = path_abs_plot * "BadFitsFrac_Peaks.pdf"
savefig(pname)



# look if bad fits come from detector species per th228 lines
detIdx_ppc = map(x->x[1]=='P',string.(dets))
detIdx_icpc = map(x->x[1]=='V',string.(dets))
detIdx_coax = map(x->x[1]=='C',string.(dets))
detIdx_bege = map(x->x[1]=='B',string.(dets))
@sprintf("%.0f PPC, %.0f ICPC, %.0f Coax, %.0f BeGe = %.0f detectors",sum(detIdx_ppc),sum(detIdx_icpc),sum(detIdx_coax),sum(detIdx_bege),sum(detIdx_ppc)+sum(detIdx_icpc)+sum(detIdx_coax)+sum(detIdx_bege))

BadFracFun = x-> [sum(pval[i,x].<0.05)/length(pval[i,:]) for i=1:7]
BadFracFun_Combi = x-> [sum(pval_combi[i,x].<0.05)/length(pval_combi[i,:]) for i=1:7]
BadFrac = BadFracFun.([detIdx_ppc,detIdx_icpc,detIdx_coax,detIdx_bege])
BadFrac_Combi = BadFracFun_Combi.([detIdx_ppc,detIdx_icpc,detIdx_coax,detIdx_bege])

colors = [:gray :green :orange :red]
p1 = groupedbar(0.75:6.75,[BadFrac[1] BadFrac[2] BadFrac[3] BadFrac[4]],bar_position=:stack,bar_width = 0.22,
        xlabel = "Th228 peak(keV)",ylabel="Fraction of fits with p < 5%",xticks=(0.75:6.75,xticklabels),
        label=["PPC" "ICPC" "Coax" "BeGe"], alpha=1,
        color = colors, box = :on,grid=false,
        legend=:topleft,guidefontsize=12,legendfontsize=10, legendtitle = "Single peak fits",
        ylims=(0,0.61))
p2 = groupedbar(0.75:6.75,[BadFrac_Combi[1] BadFrac_Combi[2] BadFrac_Combi[3] BadFrac_Combi[4]],bar_position=:stack,bar_width = 0.22,
        xlabel = "Th228 peak(keV)",ylabel="Fraction of fits with p < 5%",xticks=(0.75:6.75,xticklabels),
        label=["PPC" "ICPC" "Coax" "BeGe"], alpha=1,
        color = colors, box = :on,grid=false,
        legend=:topleft,guidefontsize=12,legendfontsize=10, legendtitle = "Combined fits",
        ylims=(0,0.61))
plot(p1,p2,layout=(2,1),size=(500,700))
pname2 = path_abs_plot * "BadFitsFrac_Peaks_DetectorTypes.pdf"
savefig(pname2)

# bad fits per detector type - sum over all lines
BadFracDetFun = x-> sum(pval[:,x].<0.05)/length(pval[:,x])
BadFracDetFun_Combi = x-> sum(pval_combi[:,x].<0.05)/length(pval_combi[:,x])

bar(0.75:3.75,[BadFracDetFun(detIdx_ppc), BadFracDetFun(detIdx_icpc), BadFracDetFun(detIdx_coax), BadFracDetFun(detIdx_bege)],bar_width = 0.22,
    xticks=(0.86:3.86,["PPC" "ICPC" "Coax" "BeGe"]),label="Single peak fits",color=:silver,alpha = 1,guidefontsize=12,legendfontsize=10)
bar!(1:4,[BadFracDetFun_Combi(detIdx_ppc), BadFracDetFun_Combi(detIdx_icpc), BadFracDetFun_Combi(detIdx_coax), BadFracDetFun_Combi(detIdx_bege)],bar_width = 0.22,
   label="Combined peak fit",color=:dodgerblue,alpha = 1,guidefontsize=12,legendfontsize=10,
   xlabel="Detector type",ylabel="Fraction of fits with p < 5%",box = :on,grid=false,ylims=(0,0.61),xtickfontsize=10)
pname3 = path_abs_plot * "BadFitsFrac_DetectorTypes.pdf"
savefig(pname3)
   
# plot fit results, but remove "bad fits"
pval_single_mask = pval .> 0.05
pval_combi_mask = pval_combi .> 0.05
pval_mask = pval_single_mask .* pval_combi_mask

# plot results
h1 = histogram((μ_combi .- μ)./µ_err,color=:blue, alpha=0.5,xlabel="Norm. residuals: µ (σ)",legend=false)
h2 = histogram((σ_combi .- σ)./σ_err,color=:blue, alpha=0.5,xlabel="Norm. residuals: σ (σ)",legend=false)
h3 = histogram((n_combi .- n)./n_err,color=:blue, alpha=0.5,xlabel="Norm. residuals: Gauss amplitude (σ)",legend=false)
h4 = histogram((step_amplitud_combi .- step_amplitude)./step_amplitude_err,color=:blue, alpha=0.5,xlabel="Norm. residuals: step. amplitude (σ)",legend=false)
h5 = histogram((skew_fraction_combi .- skew_fraction)./skew_fraction_err,color=:blue, alpha=0.5,xlabel="Norm. residuals: skew fraction (σ)",legend=false)
h6 = histogram((skew_width_combi .- skew_width)./skew_width_err,color=:blue, alpha=0.5,xlabel="Norm. residuals: skew width (σ)",legend=false)
h7 = histogram((background_combi .- background)./background_err,color=:blue, alpha=0.5,xlabel="Norm. residuals: background (σ)",legend=false)
plot(h1,h2,h3,h4,h5,h6,h7,layout=(7,1),size=(600,2000),legend=false)

negloglikdiff = negloglike_ml_combi[pval_mask] .- negloglike_ml[pval_mask]
@sprintf("NegLogLikelihood differnce smaller than 0.01 in %.1f%% of the fits",100*sum(abs.(negloglikdiff).<1e-2)./length(negloglikdiff))

pvaldiff = pval_combi[pval_mask] .- pval[pval_mask]
#minimum(negloglikdiff)
#maximum(negloglikdiff)
diff = [(μ_combi[pval_mask] .- μ[pval_mask])./µ_err[pval_mask],
        (σ_combi[pval_mask] .- σ[pval_mask])./σ_err[pval_mask],
        (n_combi[pval_mask] .- n[pval_mask])./n_err[pval_mask],
        (step_amplitude_combi[pval_mask] .- step_amplitude[pval_mask])./step_amplitude_err[pval_mask],
        (skew_fraction_combi[pval_mask] .- skew_fraction[pval_mask])./skew_fraction_err[pval_mask],
        (skew_width_combi[pval_mask] .- skew_width[pval_mask])./skew_width_err[pval_mask],
        (background_combi[pval_mask] .- background[pval_mask])./background_err[pval_mask]]

ynames = [L"$(\mu_\textrm{combined} - \mu_\textrm{single}) / \sigma_\mu$ ",
        L"$(\sigma_\textrm{combined} - \sigma_\textrm{single}) / \sigma_\sigma$ ",     
        L"$(n_\textrm{combined} - n_\textrm{single}) / \sigma_n$ ",
        L"$(\textrm{step amplitude}_\textrm{combined} - \textrm{step amplitude}_\textrm{single}) / \sigma_{\textrm{step amplitude}}$ ",
        L"$(\textrm{skew fraction}_\textrm{combined} - \textrm{skew fraction}_\textrm{single}) / \sigma_{\textrm{skew fraction}}$ ",
        L"$(\textrm{skew width}_\textrm{combined} - \textrm{skew width}_\textrm{single}) / \sigma_{\textrm{skew width}}$ ",
        L"$(\textrm{background}_\textrm{combined} - \textrm{background}_\textrm{single}) / \sigma_{\textrm{background}}$ "]  
plt_label = ["μ","σ","n","step_amplitude","skew_fraction","skew_width","background"]

for i=1:7
    hline([-3,3],color = :gray,linestyle=:dash,linewidth=2)       
    hline!([0],color = :gray,linestyle=:solid,linewidth=2)
    vline!([0],color = :gray,linestyle=:solid,linewidth=2) 
    p1 = scatter!(negloglikdiff,diff[i],legend=false,
            xlabel=L"$-log(\mathcal{L}_\textrm{combined})+log(\mathcal{L}_\textrm{single})$",
            ylabel=ynames[i], color=:dodgerblue,alpha=0.5,markerstrokewidth=0,markersize=5,
            xlims=(-5,7),box=:on,grid=false,guidefontsize=12)
    pname4 = path_abs_plot * "CompareGlobalFit_Likelihood_$(plt_label[i]).pdf"
    savefig(pname4)                    
end

for i=1:7
#i=1
        hline([-3,3],color = :gray,linestyle=:dash,linewidth=2)       
        hline!([0],color = :gray,linestyle=:solid,linewidth=2)
        vline!([0],color = :gray,linestyle=:solid,linewidth=2) 
        p1 = scatter!(pvaldiff,diff[i],legend=false,
                xlabel=L"$p_\textrm{combined}-p_\textrm{single}$",
                ylabel=ynames[i], color=:dodgerblue,alpha=0.5,markerstrokewidth=0,markersize=5,
                xlims=(-0.6,0.6),box=:on,grid=false,guidefontsize=12)
        pname5 = path_abs_plot * "CompareGlobalFit_pVal_$(plt_label[i]).pdf"
        savefig(pname5)                    
end

histogram(pvaldiff,legend=false,xlabel=L"$p_\textrm{combined}-p_\textrm{single}$",ylabel="Occurence",color=:dodgerblue,alpha=0.5,box=:on,grid=false,guidefontsize=14)
pname6 = path_abs_plot * "CompareGlobalFit_pValDiffHist.pdf"
savefig(pname6) 

histogram(abs.(pvaldiff),legend=false,xlabel=L"$|p_\textrm{combined}-p_\textrm{single}|$",ylabel="Occurence",color=:dodgerblue,alpha=0.5,box=:on,grid=false,guidefontsize=14,xlims=(0,0.6))
pname7 = path_abs_plot * "CompareGlobalFit_pValAbsDiffHist.pdf"
savefig(pname7)

quantiles = map(x->quantile(abs.(pvaldiff),x),[0.5,0.9,0.95,0.99])

t = Table(
        CombiBetterSingle= [sum(pvaldiff.>0)/length(pvaldiff) ...],
        pDiffSmaller_1p= [sum(abs.(pvaldiff).<0.01)/length(pvaldiff) ...],
        pDiffSmaller_5p= [sum(abs.(pvaldiff).<0.05)/length(pvaldiff) ...],
        pDiffLarger_5p= [sum(abs.(pvaldiff).>0.05)/length(pvaldiff) ...]
)