# try to prodove gradiend based minimzer.
using BAT,  Optim, InverseFunctions, ForwardDiff
using Distributions


x = collect(1:0.1:10)
y = 5*x.^2 .+ 1
yerr = sqrt.(y) 
y += yerr .* randn(length(x)) # randomize 

f_fit(x, v) = v[1] .* x.^2 .+ v[2] # linear fit function
v_init = ones(2)
 # function that is minimized -> chi-squared
f_opt = let X_val = x, Y_val = y, Y_err = yerr, f_fit = f_fit
    v -> sum((Y_val - f_fit(X_val, v)).^2 ./ Y_err.^2)
    end

gradient_opt = let X_val = x, Y_val = y, Y_err = yerr, f_fit = f_fit
    v -> begin
    common_grad = 2 .* (Y_val - f_fit(X_val, v)) ./ Y_err.^2
    [ sum(common_grad .* 2 .* X_val *v[1] ) ,  sum( common_grad .* 1 )]
    end
end

f_opt(v_init)
gradient_opt(v_init)

# do the fit without analytical gradient
opt_r   = @time optimize(f_opt, gradient_opt, v_init,; inplace = false)
v_fit  = Optim.minimizer(opt_r)

summary(opt_r)