
using ForwardDiff
using Optim
using Printf, LaTeXStrings
include("0nbb_likelihood.jl")
include("0nbb_get_vinput.jl")
include("0nbb_simspec.jl")

function fit_0nbb_mc(v_input::NamedTuple{((:exposure, :efficiency, :fit_window, :Qbb, :Qbb_bias, :Qbb_sigma, :sig_exp))}, events::Union{Vector{<:Real}, Vector{<:Unitful.AbstractQuantity}}; fixsignal::Union{Nothing, Real} = nothing, uncertainty = true)
    if isnothing(fixsignal)
        f_opt = let vin = stripunits_nt(stripuncert_nt(v_input)), e = ustrip.(events), lh = likelihood_partition
            v -> -log(lh(v..., vin, e))
        end
        v_init = mvalue.(ustrip.([1, 1/v_input.exposure]))
    else
        f_opt = let vin = stripunits_nt(stripuncert_nt(v_input)), e = ustrip.(events), lh = likelihood_partition
            v -> -log(lh(fixsignal, v..., vin, e))
        end
        v_init = mvalue.(ustrip.([1/v_input.exposure]))
    end
   

    opt_r = optimize(f_opt, v_init, LBFGS(), Optim.Options(time_limit = 60, show_trace=false, iterations = 10000); autodiff = :forward)
    v_ml  = Optim.minimizer(opt_r)
    loglikelihood_ml = Optim.minimum(opt_r)
    #likelihood_ml = exp( - Optim.minimum(opt_r))
    converged = Optim.converged(opt_r)
    gof = (loglikelihood_ml = loglikelihood_ml, )
    
    if isa(v_input.exposure, Unitful.AbstractQuantity)
        sig_u = 1 / unit(v_input.exposure) *u"kg" 
        bck_u = 1 / (unit(v_input.exposure) * unit(v_input.fit_window))
    else
        sig_u = bck_u = unit(1) 
    end

    if uncertainty 
        covmat = inv(ForwardDiff.hessian(f_opt, v_ml))
        v_ml_err = sqrt.(diag(abs.(covmat)))
        if isnothing(fixsignal)
            par = (signal_strength = measurement.(v_ml[1], v_ml_err[1]) * sig_u, background = measurement.(v_ml[2], v_ml_err[2]) * bck_u)
        else
            par = (signal_strength = fixsignal, background = measurement.(v_ml[1], v_ml_err[1]) * bck_u)
        end
    else
        if isnothing(fixsignal)
            par = (signal_strength = v_ml[1] * sig_u, background = v_ml[2] * bck_u)
        else
            par = (signal_strength = fixsignal, background = v_ml[1] * bck_u)
        end
    end

    result = (par..., converged = converged, gof = gof) # fit function with optimized parameters
    report = (par..., events = events, fit_window = v_input.fit_window, sig_exp = v_input.sig_exp,
             f_fit = x -> peakshape_0nbb.(Ref(v_input), Measurements.value(par.signal_strength), Measurements.value(par.background), x), gof = gof)
    return result, report
end

function plot_fit(report)
    fs = 14
    default(foreground_color_legend = :silver,
    background_color_legend = :white,
    grid = :off,
    framestyle = :semi,
    xtickfontsize = fs-2,
    ytickfontsize = fs-2,
    legendfontsize = fs-2,
    xlabelfontsize = fs + 2,
    ylabelfontsize = fs + 2)

    e = ustrip.(report.events)

    nbins = ifelse(length(e) < 100, 100, length(e))
    stephist(e, nbins = nbins, fill = true, color = :silver, alpha = 0.5, label = "MC data")
    bin_width = (maximum(e) - minimum(e)) / nbins
    e_plt = collect(2039 - ustrip(report.fit_window)/2:0.1:2039 + ustrip(report.fit_window)/2)
    plot!(e_plt, report.f_fit(e_plt) .* length(e) .* bin_width, linewidth = 2, color = :red, label = @sprintf("Fit: S = %.1e %s (T1/2 = %.1e %s)\n B = %.1e %s", ustrip(report.signal_strength) * report.sig_exp, unit(report.signal_strength), 1/(ustrip(report.signal_strength) * report.sig_exp), 1/unit(report.signal_strength), ustrip(report.background), unit(report.background)), legend = :topleft)
    xlabel!("Energy (keV)")
    ylabel!("Counts")
    ylims!(0, ylims()[2])
end

# # ## test the function
# # prepare (fake) data
# data = LegendData(:l200)
# dets_ged = DetectorId.(collect(keys(data.metadata.datasets.ovbb_partitions_pars)))
# v_input = fit_input(data, dets_ged[1], DataPartition(5); stripunits = false)
# v_input = merge(v_input, (exposure = v_input.exposure * 4000,)) #boost exposure
# if isa(v_input.Qbb, Unitful.AbstractQuantity)
#     v_test = (signal_strength = 1/2 * 1/u"yr", background = 1e-3.*1/(unit(v_input.exposure)*unit(v_input.Qbb)))
# else
#     v_test = (signal_strength = 1/2, background = 1e-3)
# end
# events = simspec(v_input, v_test...; plotFlag = true) 
# #fit 
# result, report = fit_0nbb_mc(v_input, events)
# plot_fit(report)
 
