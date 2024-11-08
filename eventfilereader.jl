using LegendDataManagement
using LegendDataManagement.LDMUtils
using LegendHDF5IO
using TypedTables
using PropertyFunctions
using Plots
using Unitful
using Dates
using StructArrays
l200 = LegendData(:l200)
t = read_ldata(l200, :jlevt, :phy, DataPeriod(3), DataRun(0))
# evt_p3 = read_ldata(l200, :jlevt, :phy, DataPeriod(3))

#evt_p3_r0.geds.multiplicity , has to be 1
#evt_p3_r0.geds.max_e_cusp_ctc_cal , max energy of all triggered events
#evt_p3_r0.geds.aoe_sg_classifier_ds_cut , psd cut 

# build physics spectrum of events that have: passed qc, multiplicity = 1, LAr cut, PSD cut
#                      &&  use only good channels (from meta data) && 
qc = t.geds.is_valid_qc .&& t.geds.is_valid_hit .&& t.geds.multiplicity .== 1
t_valid = t[findall(qc)]
lar = .!t.ged_spm.lar_cut  
t_afterlar = t[findall(qc .& lar)]
psd_cut = @pf(all($aoe_sg_classifier_ds_cut[$trig_e_ch_idxs])).(t_afterlar.geds) # $ in pf is 1 column 
t_afterpsd = t_afterlar[findall(psd_cut)]

# plot 
plot(size=(1200, 800), xunit=u"keV")
stephist!(t_valid.geds.max_e_cusp_ctc_cal, bins=0:2:3000, label="After QC", yscale=:log10)
stephist!(ustrip.(t_afterlar.geds.max_e_cusp_ctc_cal), bins=0:2:3000, label="After QC and LAr cut", yscale=:log10)
stephist!(ustrip.(t_afterpsd.geds.max_e_cusp_ctc_cal), bins=0:2:3000, label="After QC, LAr and PSD cut", yscale=:log10)
plot!(xlabel="Energy (keV)", ylabel="Counts / 2 keV", title="Energy spectrum (p03-p04)")
plot!(framestyle=:box, legend=:topright, thickness_scaling=1.8, margins=(2, :mm))
plot!(xlims=(0, 3000), ylims=(1e0, ylims()[2]), xticks=0:250:3000, dpi=600)

# get fwhm for events in roi 
cut_roi = t_afterpsd.geds.max_e_cusp_ctc_cal .> 1930u"keV" .&& t_afterpsd.geds.max_e_cusp_ctc_cal .< 2190u"keV"
t_roi = t_afterpsd[findall(cut_roi)]
t_roi.geds
ch = ChannelId.(t_roi.geds.max_e_ch)
val_sel = ValiditySelection.(Timestamp.(unix2datetime.(ustrip.(t_roi.tstart))), Ref(DataCategory(:phy)) )
det = channel2detector.(Ref(l200), val_sel, ch)

pd_ecal = l200.par.rpars.ecal.(val_sel)

getproperty.(getproperty.(getproperty.(pd_ecal, Symbol.(det)), :e_cusp_ctc), :fwhm).func_cal


for i in 1:length(t_roi)
    funcs = pd_ecal[i][Symbol(det[i])].e_cusp_ctc.fwhm.func_cal
end

# pd_ecal
# l200.par.rpars[cal_type, period, run]