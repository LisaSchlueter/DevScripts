#function process_psd_calibration(processing_config::PropDict, l200::LegendData, period::DataPeriod, run::DataRun,; reprocess::Bool=false, timeout::Int=300)
using LegendHDF5IO, LegendDSP, LegendSpecFits, LegendDataTypes, LegendDataManagement
using IntervalSets, PropertyFunctions, TypedTables, PropDicts, StatsBase
using Unitful, Formatting, LaTeXStrings, Printf, Measures, Dates, Measurements
using Measurements: value as mvalue
using Measurements: uncertainty as muncert
using Plots
using Distributed, ProgressMeter, TimerOutputs
using HDF5
using LegendDataTypes: fast_flatten, flatten_by_key, map_chunked
using Base.Iterators, StructArrays
using LegendSpecFits
using LegendDataManagement.LDMUtils
pwd()
#include("../../../Packages/legend-julia-dataflow/src/data_utils.jl")
l200 = LegendData(:l200)
period = DataPeriod(3)
run = DataRun(0)
reprocess = true
processing_config = PropDict()
filekey = start_filekey(l200, (period, run, :cal))
chinfo = Table(channelinfo(l200, filekey; system=:geds, only_processable=true)) |> filterby(@pf $low_aoe_status .== :valid)
psd_config = dataprod_config(l200).psd(filekey)
pars_ctc = get_values(l200.par.rpars.ctc[period, run])
pars_energy = get_values(l200.par.rpars.ecal[period, run])
pars_db = ifelse(l200.par.rpars.aoecal[period, run] isa LegendDataManagement.NoSuchPropsDBEntry, PropDict(), l200.par.rpars.aoecal[period, run])
pars_db = ifelse(reprocess, PropDict(), pars_db)
# create log line Tuple
log_nt = NamedTuple{(:Channel, :Detector, :Status, Symbol("Number of fitted Bands"), Symbol("μ Correction Slope"), Symbol("μ Correction Intercept"), :Error)}
chinfo_ch = chinfo[1]
ch  = chinfo_ch.channel
det = chinfo_ch.detector

if !reprocess && haskey(pars_db, det)
        @debug "Channel $(det) already processed, skip"
        pars_det = pars_db[det]
        log_ch = log_nt(ch, det, "Success", pars_det.n_compton, pars_det.μ_scs[2], pars_det.μ_scs[1], "-")
        return (processed = false, log = log_ch)
end
hitchfilename = get_hitchfilename(l200, filekey, ch)
        # load data file
    if !isfile(hitchfilename)
        @error "Hit file $hitchfilename not found"
        throw(ErrorException("Hit file not found"))
    end

psd_config_ch = merge(psd_config.default, get(psd_config, det, PropDict()) )

compton_bands  = psd_config_ch.compton_bands
compton_window = psd_config_ch.compton_window
p_value_cut    = psd_config_ch.p_value # what is this? p values threshold 
e_type         = Symbol(psd_config_ch.energy_type_aoe)

if !haskey(pars_energy, det) || !haskey(pars_energy[det], e_type)
        @error "Energy calibration for $(det) not found"
        throw(ErrorException("Energy calibration for $(det) not found"))
end