# script to test changes in LegendSpecFits package
# test new implementation of FWHM uncertainties and p-value/gof
# focus covariance matrix 
using LegendSpecFits # make sure you're in development mode: dev Packages/LegendSpecFits.jl
using LegendDataManagement
using LegendHDF5IO, HDF5
using LegendDataTypes: fast_flatten, readdata
using Statistics, StatsBase
using Plots, ColorSchemes
using Distributions # new
using Calculus, LinearAlgebra, Zygote # new
using ProgressMeter # new
using Printf
## ------------ get input for fit: energie histogram  ------------ ##

include("../utils.jl")

# get data path and filename
l200 = LegendData(:l200)
ch_geds     = channelinfo(l200,(:p03,:r000,:cal),system = :geds, only_processable = true).channel;
data        = load_hitch(l200,DataPeriod(3),DataRun(0); channel=ch_geds[1])
data_energy = [data.e_cusp...] 

## ------------ do "simple" fit --> input needed for actual fit  ------------ ##
result_simple, report_simple, th228_lines, th228_names = do_simple_calibration(data_energy); # simple calibration, get histograms around the calibration lines and peakstats

# look whats inside fit results "result_simple"
Idx =5
data_energy_hist = result_simple.peakhists[Idx] # data (energies) in form of a histogram
result_simple.peakstats[Idx] # fit results of simple calibration

## ------------ do actual fit  --> developemtn happends HERE  ------------ ##
#= comments on fit routine
gamma_peakshape(x::Real, μ::Real, σ::Real, n::Real,step_amplitude::Real, skew_fraction::Real, skew_width::Real,background::Real) # model function. x = energy vector, for which which model is evaluated
th228_fit_functions.f_fit = (x, v) -> gamma_peakshape(x, v.μ, v.σ, v.n, v.step_amplitude, v.skew_fraction, v.skew_width, v.background) # why wrap around?
=#
# do fit and look at fit results
#fit_single_peak_th228(h::Histogram, ps::NamedTuple{(:peak_pos, :peak_fwhm, :peak_sigma, :peak_counts, :mean_background), NTuple{5, T}}; uncertainty::Bool=true, fixed_position::Bool=false) where T<:Real
result_fit, report_fit = fit_single_peak_th228(data_energy_hist, result_simple.peakstats[Idx] ; uncertainty=true);
fit_par_names= propertynames(result_fit)
par_name_str = collect(map(x -> string(x),fit_par_names[1:7]))
fit_par   = result_fit[fit_par_names[1:7]]
fit_err = result_fit[:err]; fit_err = fit_err[fit_par_names[1:7]]

f_loglike_array = let f_fit=gamma_peakshape, h=data_energy_hist
    v -> - hist_loglike(x -> f_fit(x, v...), h)
end
## ------------ covariance matrix   ------------ ##
# there is problem the covmat_fit. --> fix covariance matrix
# --> first look at covariance matrix provided by fit routine (at the moment using ForwardDiff)
covmat_fit = result_fit.gof.covmat # covariance matrix from fit
isposdef(covmat_fit) #something wrong with covmat :(

function PlotCM(cm,mode)
    if mode==:cov
        colormap = :viridis
        title_str = "covariance"
    elseif mode==:corr
        cm = cor(cm);
        colormap = :greys
        title_str = "correlation"
    end
    heatmap(cm,#colorbar_title="\n\n"*title_str, 
    aspect_ratio=:equal,
    yflip=true, 
    c=cgrad(colormap, rev=true), 
    colorbar_title_location = :right,
    guidefontsize=16,
    tickfontsize=14,
    size=(1000, 600))
    xlims!(0.5,7.5); ylims!(0.5,7.5)
    yticks!(1:7,par_name_str); 
    xticks!(1:7,[par_name_str[1],par_name_str[2],par_name_str[3],"amp.","s_f","s_w","bkg"])
    title!("Fit parameter "*title_str *" matrix")
    annotate!(9.4, 4, text(title_str, :center, 16, rotation=90))
end
PlotCM(covmat_fit,:corr) ; #savefig("./plots/corrmat_forwarddiff.png")# plot correlation matrix
PlotCM(covmat_fit,:cov) ;# savefig("./plots/covmat_forwarddiff.png")# plot covariance n matrix

# --> plot fit parameter and uncertainties
# relative uncertainties for many parameters seem very large and reach into negative regime. --> maybe model too complicated given this statistics? why not simplyfy model? 
function Print_FitResult(cm)
    for i = 1:length(fit_par)
        if cm[i,i]<0
              err = - sqrt(abs(cm[i,i]))
             else
               err = sqrt(cm[i,i])
        end
       # info_str = "$(first(par_name_str[i],6)) \t =  $(round(fit_par[i],digits=2))   \t err = $(round(err,digits=2)) \t rel err = $(round(100*err/fit_par[i],digits=1))% \n"
        @printf("%s = %.2f \t err = %.2g \t rel err = %.2g %%\n",first(par_name_str[i],6),fit_par[i],err,100*err/fit_par[i])#$(first(par_name_str[i],6)) \t =  $(round(fit_par[i],digits=2))   \t err = $(round(err,digits=2)) \t rel err = $(round(100*err/fit_par[i],digits=1))% \n")
    end 
end
Print_FitResult(covmat_fit)

# ---> find nearest semi positive definite covmat (based on ForwardDiff hessian from fit )
#https://www.sciencedirect.com/science/article/pii/0024379588902236
function nearestSPD(A)
    B = (A + A') / 2 # make sure matrix is symmetric
    _, s, V = svd(B) # singular value decomposition (SVD), s = singular values (~eigenvalues), V = right singular vector  (~eigenvector)
    H = V * diagm(0 => max.(s, 0)) * V' # symmetric polar factor of B
    B = (B + H) / 2 # calculate nearest positive definite matrix
    B = (B + B') / 2  # make sure matrix is symmetric
    return B
end 
covmat = nearestSPD(covmat_fit)
corrmat = cor(covmat)
isposdef(covmat) # now its positive definite
all(diag(covmat).>0) # and all diagonal entries positive 
Print_FitResult(covmat)
PlotCM(covmat,:corr) ; #savefig("./plots/corrmat_forwarddiff_SPE.png")# plot correlation matrix
PlotCM(covmat,:cov) ; #savefig("./plots/covmat_forwarddiff_SPE.png")# plot covariance n matrix

#covmat_fit_tmp = covmat_fit
#covmat_fit_tmp[3,:] = NaN*ones(7); covmat_fit_tmp[:,3] = NaN*ones(7);
# compare uncertainties
for i = 1:length(fit_par)
    @info "uncertainty $(par_name_str[i]): $(round(fit_err[i],digits=2)) (fit)  \t vs. $(round(sqrt(covmat[i,i]),digits = 2)) (nearestSPD) \n"
end

# ---> try different ways to calculate hessian and covmat
# with Calculus.hessian

fit_par_arr =  [fit_par[f] for f in fieldnames(fit_par)]
H_calc = Calculus.hessian(f_loglike_array, fit_par_arr)
covmat_calc = inv(H_calc)
isposdef(covmat_calc) # also not semi positive definite 
Print_FitResult(covmat_calc)
PlotCM(covmat_calc,:corr) ; savefig("./plots/corrmat_calc.png")# plot correlation matrix
PlotCM(covmat_calc,:cov) ; savefig("./plots/covmat_calc.png")# plot covariance n matrix


# with Zygote.hessian
H_zyg = Zygote.hessian(f_loglike_array, fit_par_arr)
covmat_zyg = inv(H_zyg)
isposdef(covmat_zyg) # also  not semi positive definite 
Print_FitResult(covmat_zyg)
PlotCM(covmat_zyg,:corr) ; savefig("./plots/corrmat_zyg.png")# plot correlation matrix
PlotCM(covmat_zyg,:cov) ; savefig("./plots/covmat_zyg.png")# plot covariance n matrix


