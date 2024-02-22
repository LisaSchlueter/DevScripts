using LegendDataManagement
using LegendHDF5IO, HDF5
using Printf

#function load_hitch(l200::LegendData=LegendData(:l200), period::DataPeriod=DataPeriod(3),run::DataRun=DataRun(0);mode::Symbol=:cal,channel::ChannelId=nothing)
function load_hitch(l200::LegendData, period::DataPeriod,run::DataRun; mode::Symbol=:cal,channel::Union{ChannelId,Bool}=true)
        # get data path and filename, find channel 
    path       = l200.tier[:jlhitch,mode,period,run] .* "/" # get data path
    filenames  = path .* readdir(path);
    pattern =  r"(?<=-ch)(.*)(?=-tier_jlhit\.lh5)";  # regular expression, look for channel name in file
    ch_all = [match(pattern, filename).match for filename in filenames]; # extract all channel number from filename
   
    ch_geds = channelinfo(l200,(:p03,:r000,:cal),system = :geds, only_processable = true).channel
    dets = channelinfo(l200,(:p03,:r000,:cal),system = :geds, only_processable = true).detector

    if !isa(channel,ChannelId) #if not assigned, select first channel
        channel = ch_geds[1]
    end

    ch_sel_idx = findfirst(ch_all .== string(channel)[3:end])
    # read file, load energy spectrum
   # data = h5open(x -> readdata(x, "ch$(ch_all[ch_sel_idx])/dataQC"), filenames[ch_sel_idx])
    data = lh5open(filenames[ch_sel_idx])["ch$(ch_all[ch_sel_idx])/dataQC"]

    det_sel = dets[findfirst(map(x-> x == channel,ch_geds))]
    @info "Load hitch from channel $(string(channel)) / detector: $(string(det_sel))"
    return data
end

function do_simple_calibration(energy::Vector{Float64})
    th228_lines =  [583.191,  727.330,  860.564,  1592.53,    1620.50,    2103.53,    2614.51]
    th228_names =  ["Tl208a", "Bi212a", "Tl208b", "Tl208DEP", "Bi212FEP", "Tl208SEP", "Tl208FEP"]
    window_sizes =  vcat([(25.0,25.0) for _ in 1:6], (30.0,30.0))
    n_bins =  10000;
    quantile_perc=0.995;
    result_simple, report_simple = simple_calibration(energy, th228_lines, window_sizes, n_bins=n_bins,; calib_type=:th228, quantile_perc=quantile_perc)
    return result_simple, report_simple, th228_lines, th228_names
end

function save_result(fname::String, par_names::Vector{String}, par::Union{Vector{Any},Vector{Vector{Float64}}}; path::String="JuliaBasics/results")
   # find path 
    path_rel = relpath(path,pwd())
    path_abs = pwd() * "/" * path_rel * "/"
   
    file = h5open(path_abs * fname, "w")
    for (i,name) in enumerate(par_names)
        if isa(par[i],StepRangeLen)
            par[i] = collect.(par[i])
        end
        write(file, name, par[i])
    end
    close(file)

    @info "Results written to file $(path_abs * fname)"
end

function get_det_types(dets::Vector)
    detIdx_ppc = map(x->x[1]=='P',string.(dets))
    detIdx_icpc = map(x->x[1]=='V',string.(dets))
    detIdx_coax = map(x->x[1]=='C',string.(dets))
    detIdx_bege = map(x->x[1]=='B',string.(dets))
    @sprintf("%.0f PPC, %.0f ICPC, %.0f Coax, %.0f BeGe = %.0f detectors",sum(detIdx_ppc),sum(detIdx_icpc),sum(detIdx_coax),sum(detIdx_bege),sum(detIdx_ppc)+sum(detIdx_icpc)+sum(detIdx_coax)+sum(detIdx_bege))
    return detIdx_ppc, detIdx_icpc, detIdx_coax, detIdx_bege
end