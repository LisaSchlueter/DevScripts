using LegendDataManagement
using LegendSpecFits
using TypedTables
using Unitful, Measurements
using Plots

f_lin(x,p1,p2)  = p1 + p2 * x 
f_quad(x,p1,p2,p3)  = p1 * x^2 .+ p2 * x + p3
par_true = [5,2]

x = [1,2,3,4,5,6,7,8,9,10] #.± ones(10)
y = f_lin.(x,par_true...) .+ 0.5.*randn(10)

# fit 
result, report       = chi2fit(f_lin, x, y; uncertainty=true) 
result_poly, _       = chi2fit(1, x, y; uncertainty=true) 

