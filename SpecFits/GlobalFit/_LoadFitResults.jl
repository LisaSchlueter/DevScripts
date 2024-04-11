using  HDF5

path = "PackageDevScripts/DevelopSpecFits/GlobalFit/results"
path_rel = relpath(path,pwd())
path_abs = pwd() * "/" * path_rel * "/"

# init single peak fits 
function LoadSinglePeakFits(ch_ged::Vector; g_tol::Float64 = 1e-8)
    nChannel    =  length(ch_geds)
    μ = zeros(7,nChannel)
    σ = zeros(7,nChannel)
    n = zeros(7,nChannel)
    step_amplitude = zeros(7,nChannel)
    skew_fraction = zeros(7,nChannel)
    skew_width = zeros(7,nChannel)
    background = zeros(7,nChannel)
    µ_err = zeros(7,nChannel)
    σ_err = zeros(7,nChannel)
    n_err = zeros(7,nChannel)
    step_amplitude_err = zeros(7,nChannel)
    skew_fraction_err = zeros(7,nChannel)
    skew_width_err = zeros(7,nChannel)
    background_err = zeros(7,nChannel)
    negloglike_ml = zeros(7,nChannel)
    pval = zeros(7,nChannel)

    # read single peak fits 
    fnames = path_abs .* map(x->"single_peak_fit_p3_run0_$(string(x))_gtol$(g_tol).h5",ch_geds)
    for i=1:nChannel
            file = h5open(fnames[i],"r")
            μ[:,i], σ[:,i] , n[:,i] ,step_amplitude[:,i], skew_fraction[:,i], skew_width[:,i], background[:,i], negloglike_ml[:,i], pval[:,i], µ_err[:,i],σ_err[:,i] , n_err[:,i] ,step_amplitude_err[:,i], skew_fraction_err[:,i], skew_width_err[:,i], background_err[:,i] = 
            [read(file,parname) for parname in  ["μ","σ","n","step_amplitude","skew_fraction","skew_width","background","negloglike_ml","pval","err_μ","err_σ","err_n","err_step_amplitude","err_skew_fraction","err_skew_width","err_background"]]
            close(file)
    end

    # re-arrange single fit results, because in different order than combined fit results
    Idx_order = [5,3,2,6,4,1,7]
    μ = µ[Idx_order,:]
    σ = σ[Idx_order,:]
    n = n[Idx_order,:]
    step_amplitude = step_amplitude[Idx_order,:]
    skew_fraction = skew_fraction[Idx_order,:]
    skew_width = skew_width[Idx_order,:]
    background = background[Idx_order,:]
    µ_err = µ_err[Idx_order,:]
    σ_err = σ_err[Idx_order,:]
    n_err = n_err[Idx_order,:]
    step_amplitude_err = step_amplitude_err[Idx_order,:]
    skew_fraction_err = skew_fraction_err[Idx_order,:]
    skew_width_err = skew_width_err[Idx_order,:]
    background_err = background_err[Idx_order,:]
    negloglike_ml = negloglike_ml[Idx_order,:]
    pval = pval[Idx_order,:]

    return   μ, σ , n ,step_amplitude, skew_fraction, skew_width, background, negloglike_ml, pval, µ_err,σ_err , n_err ,step_amplitude_err, skew_fraction_err, skew_width_err, background_err
end

function LoadCombiPeakFits(ch_geds)
    nChannel    =  length(ch_geds)
    prior = 2
    fit_func = :f_fit_uncorr
    # init combi peak fits 
    μ_combi = zeros(7,nChannel)
    σ_combi = zeros(7,nChannel)
    n_combi = zeros(7,nChannel)
    step_amplitude_combi = zeros(7,nChannel)
    skew_fraction_combi = zeros(7,nChannel)
    skew_width_combi = zeros(7,nChannel)
    background_combi = zeros(7,nChannel)
    negloglike_ml_combi = zeros(7,nChannel)
    pval_combi = zeros(7,nChannel)
    chi2_combi = zeros(7,nChannel)

    fnames_combi = path_abs .* map(x->"combi_peak_fit_p3_run0_$(string(x))_prior$(prior)_$(string(fit_func)).h5",ch_geds)
    for i=1:nChannel
            file = h5open(fnames_combi[i],"r")
            μ_combi[:,i],σ_combi[:,i], n[:,i] ,step_amplitude_combi[:,i], skew_fraction_combi[:,i], skew_width_combi[:,i], background_combi[:,i],negloglike_ml_combi[:,i], pval_combi[:,i]  = [read(file,parname) for parname in ["μ","σ","n","step_amplitude","skew_fraction","skew_width","background","negLogLike_peakwise","pval_peakwise"]]
            close(file)
    end
    return  μ_combi,σ_combi, n,step_amplitude_combi, skew_fraction_combi, skew_width_combi, background_combi,negloglike_ml_combi, pval_combi 
end

function LoadPeakFits_g_tol(ch_ged,g_tols)

    fnames = path_abs .* map(x->"single_peak_fit_p3_run0_$(string(ch_ged))_gtol$(x).h5",g_tols)
    pval = zeros(length(g_tols),7)
    for i=1:length(g_tols)
        file = h5open(fnames[i],"r")
        pval[i,:] = read(file,"pval")
        close(file)
    end
    return pval 
end
