using LegendDataManagement
using LegendSpecFits
using TypedTables
using Unitful, Measurements
using Plots

l200 = LegendData(:l200)
period = DataPeriod(3)
run = DataRun(0)
dets = channelinfo(l200,(period,run,:cal),system = :geds, only_processable = true).detector;

# get literature values for th228 lines in [keV]
filekey = start_filekey(l200, (period, run, :cal))
th228_lines = dataprod_config(l200).energy(filekey).default.th228_lines  
th228_names = Symbol.(dataprod_config(l200).energy(filekey).default.th228_names)
th228_lines_dict = Dict(th228_names .=> th228_lines)

# get fit values of peak positions in [ADC]
det = Symbol(dets[1]) # select detector
e_type = :e_zac_ctc
pars_db = l200.par.rpars.ecal[period,run];
µ    = [pars_db[det][e_type].fit[th228_names[x]].µ for x in 1:length(th228_names)]./pars_db[det][e_type].m_calib
fwhm_uncal = [pars_db[det][e_type].fit[th228_names[x]].fwhm for x in 1:length(th228_names)]./pars_db[det][e_type].m_calib

# fit 
f_lin(x,p1,p2)  = p1 + p2 * x
f_quad(x,p1,p2,p3)  = p1 + p2 .* x .+ p3 .* x.^2
result, report       = chi2fit(f_lin, µ, ustrip.(th228_lines); uncertainty=true) 
result_poly, _       = chi2fit(1, µ, ustrip.(th228_lines); uncertainty=true) 
µ_cal = f_lin(µ,result.par...)u"keV"

plot_fit(µ,th228_lines,result)
plot_fit(µ,th228_lines,result_poly)

function plot_fit(x,y,result)
    x_bf = range(minimum(x),stop=maximum(x),length=100)
    y_bf = result.f_fit_v(x_bf)
    plt = plot(x,y,marker = :dor, markersize=3, label = "data",color = :black)
    plot!(x_bf,y_bf,label = "fit",marker=false,xlabel="Energy [ADC]",ylabel="Literature energy [keV]",title="Calibration fit")
    return plt
end

# to do: write automatic calibraiton function as output 
# how to read functions that are stored as strings in json files. 
# filename = "blabla.json"
# props = readlprobs(filename)
# ljl_propfunc(props.h)(x)