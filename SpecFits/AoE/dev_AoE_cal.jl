include("startup_AoEfits.jl")

plt_path = "./PackageDevScripts/DevelopSpecFits/AoE/plots/"

e_cal, aoe = nothing, nothing
data_hit = LHDataStore(hitchfilename, "r");
tab_data = data_hit["$(ch)/dataQC/"]
# get a
a = tab_data.a[:];
# get energy for best resolution
ecal_func_str = pars_energy[det][e_type].cal.func
e_cal = collect(ljl_propfunc(ecal_func_str).(tab_data))
# get aoe
aoe = ustrip.(a ./ e_cal);
            
p = histogram2d(e_cal, aoe, nbins=(0:0.5:3000, 0.2:5e-4:0.8), xlims=(0, 3000), ylims=(0.1, 0.9), size=(1200, 800), label = "A/E uncorrected",color=cgrad(:magma), colorbar_scale=:log10, legend=:topleft, xlabel="Energy", ylabel="A/E (a.u.)", margin=5mm)
plot!(p,guidefontsize=25,xguidefontsize = 25,yguidefontsize = 25,xtickfontsize = 18,ytickfontsize=18)
xticks!(p, 0:500:3000)
title!(p, "AoE uncorrected",titlefontsize=25)
savefig(p,plt_path*"Aoe_uncorrected.png")

# start corrections/normalization 
result_fit, report_fit, compton_band_peakhists = nothing, nothing, nothing
# get compton band peak histograms with generated peakstats
compton_band_peakhists = generate_aoe_compton_bands(aoe, e_cal, compton_bands, compton_window)
result_fit, report_fit = fit_aoe_compton(compton_band_peakhists.peakhists, compton_band_peakhists.peakstats, compton_bands,; uncertainty=true)

#compton_bands = [band for band in keys(result_fit) if result_fit[band].gof.pvalue >= p_value_cut] # end here 
μ = [result_fit[band].μ for band in compton_bands]
σ = [result_fit[band].σ for band in compton_bands]

# fit μ and σ with correction functions
# aoe_corrections = nothing
result, report = nothing, nothing
result, report = fit_aoe_corrections(compton_bands, μ, σ,; e_expression = ecal_func_str)
        
# plot μ correction functions  
p = scatter(report.report_µ.x*u"keV", report.report_µ.y, ms=3, color=:black, layout = @layout[grid(2, 1, heights=[0.8, 0.2])], size = (800,570),label="Compton band fits: Gaussian µ(A/E)", margin=10mm, framestyle=:box)
plot!(xlabel="Energy (keV)", ylabel="µ(A/E) (a.u.)", xticks = ustrip.(minimum(report.report_µ.x):200: maximum(report.report_µ.x)), xlims=ustrip.((minimum(report.report_µ.x)-50, maximum(report.report_µ.x)+50)), legend = :topright,legendfontsize=15)
plot!(ylims=(0.98*mvalue(median(report.report_µ.y)), 1.02*mvalue(median(report.report_µ.y))), subplot=1, xlabel="", xticks = :none, bottom_margin=0mm)
plot!(report.report_µ.x, report.report_µ.f_fit(report.report_µ.x), linealpha=1, label="Best Fit: $(mvalue(round(result.µ_compton.par[1], digits=2))) + Ecal * $(mvalue(round(ustrip(result.µ_compton.par[2]) * 1e6, digits=2)))1e-6", line_width=5, color=:red, subplot=1, xformatter=_->"")
hline!([0.0],c = :black, linestyle = :dash, label = "", subplot=2)
plot!(report.report_µ.x, mvalue.((report.report_µ.y-report.report_µ.f_fit(report.report_µ.x))./muncert.(report.report_µ.y)) , label="", ylabel="Residuals (σ) \n", color=:black, st=:scatter, ylims = (-7,7), yticks = [-5,0,5], markershape=:circle, subplot=2, framestyle=:box)
plot!(p,guidefontsize=18,xguidefontsize = 18,yguidefontsize = 18,xtickfontsize = 15,ytickfontsize=15)
#title!(get_plottitle(filekey, det, "A/E μ"), subplot=1)
savefig(p,plt_path*"ComptonBand_µ.png" )

p = scatter(report.report_σ.x*u"keV", report.report_σ.y, ms=5, color=:black, layout = @layout[grid(2, 1, heights=[0.8, 0.2])],  size = (800,570),label="Compton band fits: Gaussian σ(A/E)", margin=10mm, framestyle=:box)
plot!(xlabel="Energy (keV)", ylabel="σ(A/E) (a.u.)", xticks = ustrip.(minimum(report.report_σ.x):200: maximum(report.report_σ.x)), xlims=ustrip.((minimum(report.report_σ.x)-50, maximum(report.report_σ.x)+50)), legend = :topright,legendfontsize=15)
plot!(ylims=(0.2*mvalue(minimum(report.report_σ.y)), 1.2*mvalue(maximum(report.report_σ.y))), subplot=1, xlabel="", xticks = :none, bottom_margin=0mm)
idx_plt = sortperm(report.report_σ.x)
plot!(report.report_σ.x, mvalue.(report.report_σ.f_fit(report.report_σ.x)), linealpha=1, label="Best Fit: sqrt($(round(mvalue(result.σ_compton.par[1])*1e6, digits=1))e-6 + $(round(ustrip(mvalue(result.σ_compton.par[2])), digits=2)) / Ecal^2)", line_width=5, color=:red, subplot=1, xformatter=_->"")  
hline!([0.0],c = :black, linestyle = :dash, label = "", subplot=2)
plot!(report.report_σ.x, mvalue.((report.report_σ.y-report.report_σ.f_fit(report.report_σ.x))./muncert.(report.report_σ.y)) , label="", ylabel="Residuals (σ) \n", line_width=2, color=:black, st=:scatter, ylims = (-20, 20), yticks = [-15,0,15], markershape=:circle, subplot=2, framestyle=:box)
plot!(p,guidefontsize=18,xguidefontsize = 18,yguidefontsize = 18,xtickfontsize = 15,ytickfontsize=15)
savefig(p,plt_path*"ComptonBand_σ.png" )

# correct aoe
aoe_corr = ljl_propfunc(result.func).(tab_data)#[a_sel] 
p = histogram2d(e_cal, aoe_corr, nbins=(0:0.5:3000, -20:0.1:10), size=(1200, 800), label = "A/E uncorrected",color=cgrad(:magma), colorbar_scale=:log10, legend=:topleft, xlabel="Energy", ylabel="A/E (a.u.)", margin=5mm)
plot!(p, xlims=(0, 3000), ylims=(-21,11))
plot!(p,guidefontsize=25,xguidefontsize = 25,yguidefontsize = 25,xtickfontsize = 18,ytickfontsize=18)
xticks!(p, 0:500:3000)
title!(p, "AoE energy corrected, normalized by σ",titlefontsize=25)
savefig(p,plt_path*"Aoe_corrected.png")

####################################################################################################
# sanity check of correction function: µ correction:
# 1. with functions 
    aoe_µcorr_funcstr = result.µ_compton.func
    aoe_µcorr = ljl_propfunc(aoe_µcorr_funcstr).(tab_data) # correction for each aoe 
# 2. by hand 
    aoe_µcorr_man = mvalue.(result.µ_compton.par[1] .+ e_cal .* result.µ_compton.par[2])
# 3. compare 
    @info "Function and manual µ correction is the same: $(all(aoe_µcorr .== aoe_µcorr_man))"
    all(aoe_µcorr .== aoe_µcorr_man)
    
# sanity check of correction function: σ correction:
# 1. with functions 
aoe_σcorr_funcstr = result.σ_compton.func
aoe_σcorr = ljl_propfunc(aoe_σcorr_funcstr).(tab_data) # correction for each aoe 
# 2. by hand 
aoe_σcorr_man = mvalue.(sqrt.( abs.(result.σ_compton.par[1]) .+  abs.(result.σ_compton.par[2]) ./ e_cal.^2 ))          
# 3. compare 
@info "Function and manual σ correction is the same: $(all(aoe_σcorr  .== aoe_σcorr_man))"
all(aoe_σcorr  .== aoe_σcorr_man)

# sanity check of correction function: total correction:
# 1. with functions 
aoe_corr_funcstr = result.func
aoe_corr = ljl_propfunc(aoe_corr_funcstr).(tab_data) # correction for each aoe 
# 2. by hand 
aoe_corr_man = (aoe .- aoe_µcorr)./ aoe_σcorr   
# 3. compare 
@info "Function and manual σ correction is the same: $(all(aoe_corr  .== aoe_corr_man))"
all(aoe_corr  .== aoe_corr_man)



