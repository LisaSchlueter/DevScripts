"""
alternative to root search algorithm, turned out to be less robust
using Interpolations
#  --------------- fwhm: alternative calculation  ---------------------------
#function get_fwhm(bin_centers,model_bf,fit_par_mc)

fit_par_mc = named_tuples[4]
model_bf = model_fun(named_tuples[4])

    # find maximum
   # try
    interp = CubicSplineInterpolation(bin_centers,model_bf);
   # catch
   #     fwhm = NaN
   #     return 
   # end
    x_inter = fit_par_mc.μ - fit_par_mc.σ:0.001:fit_par_mc.μ + fit_par_mc.σ
    AmpHalf = maximum(interp(x_inter))/2

    # find robust interpolation window 
    Idx_low1= findfirst(x -> x.>(AmpHalf/500), model_bf) # start interpolation at 1/10 of max 
    Idx_low2= findlast(x -> x.<(fit_par_mc.μ), bin_centers)-1
    Idx_up1 = findfirst(x -> x.>(fit_par_mc.μ), bin_centers)+1
    Idx_up2 = Idx_up1+findlast(x -> x.>(AmpHalf/500), model_bf[Idx_up1:end])
    # sanity check 
  #  all([model_bf[i] <= model_bf[i+1] for i in Idx_low1:Idx_low2-1]) # all increasing
   # all([model_bf[i] >= model_bf[i+1] for i in Idx_up1:Idx_up2-1]) # all decreasing
   
   #plot(bin_centers[Idx_low1:Idx_low2],model_bf[Idx_low1:Idx_low2],marker=".")
    
   # try
    interp_low = LinearInterpolation(model_bf[Idx_low1:Idx_low2],collect(bin_centers[Idx_low1:Idx_low2]));
    xhalf_low = interp_low.(AmpHalf)
    interp_up = LinearInterpolation(reverse(model_bf[Idx_up1:Idx_up2]),reverse(collect(bin_centers[Idx_up1:Idx_up2])));
    xhalf_up = interp_up.(AmpHalf)
    #catch
     #   fwhm = NaN
    #    return 
   # end
    fwhm = xhalf_up - xhalf_low
    return fwhm
#end