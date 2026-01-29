#=
Load daily import rates.
This can be used as an imput to simtree and simforest to generate infection import times (timport).
=#

using CSV
using DataFrames
using NBPMscape

## Read in daily mean imports
mean_total_mean_imports_num_vec = CSV.read( joinpath( @__DIR__, "..", "data", "daily_mean_imports_df.csv" ), DataFrame )[:,:mean_imports]
const GLEAM_DAILY_IMPORT_RATES = mean_total_mean_imports_num_vec