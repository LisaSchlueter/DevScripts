include("0nbb_likelihood.jl")

function fit_0nbb_mc(v_input::NamedTuple{((:exposure, :efficiency, :fit_window, :Qbb, :Qbb_sigma, :sig_exp))}, events::Vector{<:Real}; fixsignal::Union{Nothing, Real} = nothing, uncertainty = true)
    if isnothing(fixsignal)
        f_opt = let v_input = v_input, e = events, lh = likelihood_partition
            v -> -log(lh(v..., v_input, e))
        end
        v_init = [1, 1/v_input.exposure]
    else
        f_opt = let v_input = v_input, e = events, lh = likelihood_partition
            v -> -log(lh(fixsignal, v..., v_input, e))
        end
        v_init = [1/v_input.exposure]
    end
   
    opt_r = optimize(f_opt, v_init, LBFGS(), Optim.Options(time_limit = 60, show_trace=false, iterations = 10000); autodiff = :forward)
    v_ml  = Optim.minimizer(opt_r)
    loglikelihood_ml = Optim.minimum(opt_r)
    #likelihood_ml = exp( - Optim.minimum(opt_r))
    converged = Optim.converged(opt_r)
    gof = (loglikelihood_ml = loglikelihood_ml, )
    
    if uncertainty 
        covmat = inv(ForwardDiff.hessian(f_opt, v_ml))
        v_ml_err = sqrt.(diag(abs.(covmat)))
        if isnothing(fixsignal)
            par = (signal_strength = measurement.(v_ml[1], v_ml_err[1]), background = measurement.(v_ml[2], v_ml_err[2]))
        else
            par = (signal_strength = fixsignal, background = measurement.(v_ml[1], v_ml_err[1]))
        end
    else
        if isnothing(fixsignal)
            par = (signal_strength = v_ml[1], background = v_ml[2])
        else
            par = (signal_strength = fixsignal, background = v_ml[1])
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

    nbins = ifelse(length(report.events) < 100, 100, length(report.events))
    stephist(report.events, nbins = nbins, fill = true, color = :silver, alpha = 0.5, label = "MC data")
    bin_width = (maximum(report.events) - minimum(report.events)) / nbins
    e_plt = ustrip.(collect(2039 - report.fit_window/2:0.1:2039 + report.fit_window/2))
    plot!(e_plt, report.f_fit(e_plt) .* length(report.events) .* bin_width, linewidth = 2, color = :red, label = @sprintf("Fit: S = %.1e 1/yr, B = %.1e cts / keV", report.signal_strength * report.sig_exp, report.background), legend = :topleft)
    xlabel!("Energy (keV)")
    ylabel!("Counts")
    ylims!(0, ylims()[2])
end