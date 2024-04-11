using LegendEventAnalysis
using LegendDataManagement
using LegendHDF5IO

# functions: 
l200 = LegendData(:l200)

search_disk(Filekey,l200.tier[:jlhit,:cal,DataPeriod(3),DataRun(0)]) #.* "/" # get data path