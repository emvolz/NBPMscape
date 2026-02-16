#= Investigating sensitivity of outbreak simulation to different infection importation distributions
- t-distribution without truncation
- t-distribution with truncation inversely proportional to number of imports. Truncation at 0.1% =(1/1000=1/nimports) quantile.
- t-distribution with greater truncation at 1% quantile.
- Poisson sampling from mean of the mean daily imports across all possible outbreak origin countries as output from GLEAM
=# 

using Distributions
using NBPMscape
using Plots
using StatsBase
using Statistics

# Initialize default parameters for NBPMscape
initialize_parameters();
p = NBPMscape.P;


## Looking at time between 1st and 2nd importation

# Non-truncated t-distribution
for i in 1:10
	println( map(_->rand(TDist(p.import_t_df))*p.import_t_s, 1:p.nimports) |> sort |> diff |> first )
end
# 14.521899819473575
# 21.848369254244332
# 9.546490633168602
# 27.845490384386636
# 3.561009859749454
# 59.798212057651696
# 5.423592872986916
# 298.1229260777258
# 3.7375332281158222
# 120.67275666350142


# Truncated version
for i in 1:10
	timports = map(_->rand(TDist(p.import_t_df))*p.import_t_s, 1:p.nimports) ; timports = timports[ timports .> quantile(timports,.01) ] ; 
	println(timports |> sort |> diff |> first )
end
# 0.13689941377399606
# 0.46853349720706916
# 0.45432702396122693
# 1.8592060793580671
# 0.5032125075112148
# 0.6145582449799534
# 2.1310981032663037
# 0.3711058134267944
# 2.541472582797809
# 1.1611962860681686


## Comparing different timport methods

## Version extracted from core.jl (Master branch, commit 816b16af46de8c3dc9d1440ac4ef05540a345c45)
# This truncates the t distribution at 0.1% quantile (given nimports = 1000, truncation quantile = 1/nimports), thereby removing outliers
# Define params
p=NBPMscape.P; initialtime=0.0; maxtime=90.0; maxgenerations=100; initialcontact=:H; importmodel=:TDist; max_cases=100000000;
timport_reps = 1000
timport_1_2_diff_vec = []

# Generate vector of differences between 1st and 2nd importation time
for i in 1:timport_reps
    if importmodel == :Poisson # for testing  
		nimports = rand(Poisson((maxtime-initialtime)*p.importrate))
		nimports = max(1, nimports )
		timports = map(_->rand(Uniform(initialtime,maxtime)), 1:nimports)
	elseif importmodel == :TDist  # based on 1st wave covid 
		nimports = 0
		while nimports == 0 
			timports = map(_->rand(TDist(p.import_t_df))*p.import_t_s, 1:p.nimports)
			# use theoretical quantiles to prohibit outliers messing things up 
			tlb = quantile( TDist( p.import_t_df ), 1/p.nimports )*p.import_t_s 
			timports = timports[ (timports .> tlb) .& (timports .< (tlb+maxtime) ) ]
			nimports = length( timports )
		end
	end
	sort!(timports)
	timports .-= minimum( timports )

	# Compute time between 1st and 2nd importation time
	timport_1_2_diff = timports |> sort |> diff |> first
    push!( timport_1_2_diff_vec, timport_1_2_diff)

end

# Repeat using t distribution without truncation
tdist_v1 = []
for i in 1:timport_reps
    t_import_1_2_diff = map(_->rand(TDist(p.import_t_df))*p.import_t_s, 1:p.nimports) |> sort |> diff |> first
    push!( tdist_v1 , t_import_1_2_diff )
end

# Repeart using t distribution with truncation at 1% quantile, thereby removing larger proportion of tail.
tdist_v2 = []
for i in 1:timport_reps
    timports = map(_->rand(TDist(p.import_t_df))*p.import_t_s, 1:p.nimports)
    timports = timports[ timports .> quantile(timports,.01) ] 
    t_import_1_2_diff = timports |> sort |> diff |> first
    push!( tdist_v2 , t_import_1_2_diff )
end

# Repeat with timports generated using Poisson sampling from GLEAM daily mean imports. The mean is taken of the daily mean imports across all possible outbreak countries
gleam_diff = []
nimports = 0
for i in 1:timport_reps # i=1
    while nimports == 0 
        # Find indices for elements in vector that are not 'missing'
        days_with_import_rate = findall(!ismissing, GLEAM_DAILY_IMPORT_RATES)
        # Generate series of daily importation counts by Poisson sampling daily time series of mean daily importations
        daily_importation_counts = [ rand(Poisson(λ)) for λ in skipmissing(GLEAM_DAILY_IMPORT_RATES[1:Int(maxtime)]) ]
        #sum(daily_importation_counts)
        #println(daily_importation_counts)
        timports = []
        for i in 1:length( daily_importation_counts ) # i=25
            t_start = days_with_import_rate[ i ] - 1 # Subtract 1 because 1st element in vector is for day 0.
            t_end = days_with_import_rate[ i ]
            # Generate importation time for each import using a Uniform distribution
            daily_timports = map(_->rand( Uniform(t_start, t_end )), 1:daily_importation_counts[ i ] )
            append!(timports, daily_timports)
        end
        nimports = length( timports )
    end
    t_import_1_2_diff = timports |> sort |> diff |> first
    push!( gleam_diff , t_import_1_2_diff )
end

### Plot results of difference between 1st and 2nd importation times using the four different importation time generation methods

# Constant x-axis limits (for comparison across methods)
xlimit_vals = [0,1500]
ylimit_vals = [0,1000]
h1 = histogram( timport_1_2_diff_vec, alpha = 0.5, xlimit = xlimit_vals, ylimit=ylimit_vals, label = "current NBPMscape"
                , title = "Time between 1st and 2nd importation\n with different distributions" )
h2 = histogram( tdist_v1  , alpha = 0.5, xlimit = xlimit_vals, ylimit=ylimit_vals, label = "full t-dist")
h3 = histogram( tdist_v2  , alpha = 0.5, xlimit = xlimit_vals, ylimit=ylimit_vals, label = "truncated t-dist")
h4 = histogram( gleam_diff, alpha = 0.5, xlimit = xlimit_vals, ylimit=ylimit_vals, label = "Poisson sampling from GLEAM output (equal weight to outbreak countries)"
				,xlabel = "Days between 1st and 2nd importation")

plot( h1,h2,h3,h4, layout = (4, 1), size = (650,1000))

# Variable x-axis limits (to better view data for individual methods)
h1 = histogram( timport_1_2_diff_vec, alpha = 0.5, ylimit=ylimit_vals, label = "current NBPMscape"
                , title = "Time between 1st and 2nd importation\n with different distributions" )
h2 = histogram( tdist_v1  , alpha = 0.5, ylimit=ylimit_vals, label = "full t-dist")
h3 = histogram( tdist_v2  , alpha = 0.5, ylimit=ylimit_vals, label = "truncated t-dist")
h4 = histogram( gleam_diff, alpha = 0.5, ylimit=ylimit_vals, label = "Poisson sampling from GLEAM output (equal weight to outbreak countries)"
				,xlabel = "Days between 1st and 2nd importation")

plot( h1,h2,h3,h4, layout = (4, 1), size = (650,1000))

## Looking at the number of importations over time
## using the different distributions

"""
Description	Transforms vector of times to time series of cumulative counts	
Example		times = [1.0, 2.0, 2.0, 5.0, 7.0]
			t_axis, N = cumulative_event_series(times)
"""
function cumulative_event_series(times::AbstractVector)
    cm = countmap(times)
    ts = sort(collect(keys(cm)))
    counts = [cm[t] for t in ts]
    return ts, cumsum(counts)
end

# Compute importation times with different t-distributions
timports_1 = []; timports_2 = []; timports_3 = []; timports_4 = []

## Version extracted from core.jl (Master branch, commit 816b16af46de8c3dc9d1440ac4ef05540a345c45)
# This truncates the t distribution at 0.1% quantile (given nimports = 1000, truncation quantile = 1/nimports), thereby removing outliers
for i in 1:timport_reps
    if importmodel == :Poisson # for testing  
		nimports = rand(Poisson((maxtime-initialtime)*p.importrate))
		nimports = max(1, nimports )
		timports = map(_->rand(Uniform(initialtime,maxtime)), 1:nimports)
	elseif importmodel == :TDist  # based on 1st wave covid 
		nimports = 0
		while nimports == 0 
			timports = map(_->rand(TDist(p.import_t_df))*p.import_t_s, 1:p.nimports)
			# use theoretical quantiles to prohibit outliers messing things up 
			tlb = quantile( TDist( p.import_t_df ), 1/p.nimports )*p.import_t_s 
			timports = timports[ (timports .> tlb) .& (timports .< (tlb+maxtime) ) ]
			nimports = length( timports )
		end
	end
	sort!(timports)
	timports .-= minimum( timports )
    
	push!( timports_1, timports)
end


# t distribution without truncation
for i in 1:timport_reps
    timports = map(_->rand(TDist(p.import_t_df))*p.import_t_s, 1:p.nimports)
	sort!(timports)
	timports .-= minimum( timports )
    push!( timports_2 , timports )
end

# Repeart using t distribution with truncation at 1% quantile, thereby removing larger proportion of tail.
for i in 1:timport_reps
    timports = map(_->rand(TDist(p.import_t_df))*p.import_t_s, 1:p.nimports)
    timports = timports[ timports .> quantile(timports,.01) ] 
	sort!(timports)
	timports .-= minimum( timports )
    push!( timports_3 , timports )
end

# Repeat with timports generated using Poisson sampling from GLEAM daily mean imports. 
# The mean is taken of the daily mean imports across all possible outbreak countries
for i in 1:timport_reps
    nimports = 0
    while nimports == 0 
        # Find indices for elements in vector that are not 'missing'
        days_with_import_rate = findall(!ismissing, GLEAM_DAILY_IMPORT_RATES)
        # Generate series of daily importation counts by Poisson sampling daily time series of mean daily importations
        daily_importation_counts = [ rand(Poisson(λ)) for λ in skipmissing(GLEAM_DAILY_IMPORT_RATES[1:Int(maxtime)]) ]
        #sum(daily_importation_counts)
        #println(daily_importation_counts)
        timports = []
        for i in 1:length( daily_importation_counts ) # t=11
            t_start = days_with_import_rate[ i ] - 1 # Subtract 1 because 1st element in vector is for day 0.
            t_end = days_with_import_rate[ i ]
            # Generate importation time for each import using a Uniform distribution
            daily_timports = map(_->rand( Uniform(t_start, t_end )), 1:daily_importation_counts[ i ] )
            append!(timports, daily_timports)
        end
        nimports = length( timports )
    end

    sort!(timports)
    timports .-= minimum( timports )
    push!( timports_4 , timports )
end

## Plot the first 10 timport series
# With fixed y-axis limits, so can compare across the different methods
xlimit_vals = [0,90]
ylimit_vals = [0,300000]

p1 = plot( cumulative_event_series(timports_1[1]), xlimit = xlimit_vals, ylimit = ylimit_vals, xlabel = "Days since first importation", label=""
        , title = "Current NBPMscape timports")
for i in 2:10;	p1 = plot!(cumulative_event_series(timports_1[i]), label=""); end

p2 = plot( cumulative_event_series(timports_2[1]), xlimit = xlimit_vals, ylimit = ylimit_vals, xlabel = "Days since first importation", ylabel = "Cumulative number of importations", label=""
        , title = "No truncation on t-distribution")
for i in 2:10; p2 = plot!(cumulative_event_series(timports_2[i]), label=""); end

p3 = plot( cumulative_event_series(timports_3[1]), xlimit = xlimit_vals, ylimit = ylimit_vals, xlabel = "Days since first importation", ylabel = "Cumulative number of importations", label=""
        , title = "t-distribution left tail truncated at 1% quantile")
for i in 2:10; p3 = plot!(cumulative_event_series(timports_3[i]), label=""); end

p4 = plot( cumulative_event_series(timports_4[1]), xlimit = xlimit_vals, ylimit = ylimit_vals, xlabel = "Days since first importation", ylabel = "Cumulative number of importations", label=""
        , title = "Poisson sampling from GLEAM output\n (equal weight to outbreak countries)")
for i in 2:10; p4 = plot!(cumulative_event_series(timports_4[i]), label=""); end

plot(p1, p2, p3, p4, layout = (4, 1), size = (650,1000)) #default size is 600,400

plot( cumulative_event_series(timports), xlabel = "Days since first importation" #,ylabel = "Cumulative number of importations"
			, xlimit = [0,90], ylimit = [0,p.nimports]
			, label="", title = "t-distribution left tail truncated at 1% quantile")

# Plot the first 10 timport series
# With fixed y-axis limits, so can better see data from individual methods
p1 = plot( cumulative_event_series(timports_1[1])
			, xlabel = "days since first importation"#, ylabel = "Cumulative number of importations"
			, xlimit = xlimit_vals
			, label="", title = "Current NBPMscape timports")
for i in 2:10; p1 = plot!(cumulative_event_series(timports_1[i]), label=""); end

p2 = plot( cumulative_event_series(timports_2[1]), xlabel = "Days since first importation", ylabel = "Cumulative number of importations"
			, xlimit = xlimit_vals
			, label="", title = "No truncation on t-distribution")
for i in 2:10; p2 = plot!(cumulative_event_series(timports_2[i]), label=""); end

p3 = plot( cumulative_event_series(timports_3[1]), xlabel = "Days since first importation" #,ylabel = "Cumulative number of importations"
			, xlimit = xlimit_vals
			, label="", title = "t-distribution left tail truncated at 1% quantile")
for i in 2:10; p3 = plot!(cumulative_event_series(timports_3[i]), label=""); end

p4 = plot( cumulative_event_series(timports_4[1]), xlimit = xlimit_vals, xlabel = "Days since first importation", ylabel = "Cumulative number of importations", label=""
        , title = "Poisson sampling from GLEAM output\n (equal weight to outbreak countries)")
for i in 2:10; p4 = plot!(cumulative_event_series(timports_4[i]), label=""); end

plot(p1, p2, p3, p4, layout = (4, 1), size = (650,1000)) #default size is 600,400

################
## TEST Run with updated version in secondary_care_sampling branch commit f676383ae33166a2550814d7f011676311a7de52
sims_simid_tinf_df_1632923_1 = load("covidlike-1.4.3-sims-nrep1_sims_simid_tinf_df_1632923.1.jld2","sims_simid_tinf_df")
sims_simid_tinf_df_1632923_2 = load("covidlike-1.4.3-sims-nrep1_sims_simid_tinf_df_1632923.2.jld2","sims_simid_tinf_df")
sims_max_cases_df_1632923_1 = load("covidlike-1.4.3-sims-nrep1_max_cases_df_1632923.1.jld2","sims_max_cases_df")
sims_max_cases_df_1632923_2 = load("covidlike-1.4.3-sims-nrep1_max_cases_df_1632923.2.jld2","sims_max_cases_df")
sims_filtered_G_1632923_1 = load("covidlike-1.4.3-sims-nrep1_filtered_G_1632923.1.jld2","sims_G_filtered")
sims_filtered_G_1632923_2 = load("covidlike-1.4.3-sims-nrep1_filtered_G_1632923.2.jld2","sims_G_filtered")

################

## Run simulation with updated timport code in simforest()
"""
if importmodel == :Poisson # for testing  
		nimports = rand(Poisson((maxtime-initialtime)*p.importrate))
		nimports = max(1, nimports )
		timports = map(_->rand(Uniform(initialtime,maxtime)), 1:nimports)
	elseif importmodel == :TDist  # based on 1st wave covid 
		nimports = 0
		while nimports == 0 
			timports = map(_->rand(TDist(p.import_t_df))*p.import_t_s, 1:p.nimports) # sort(timports) # scatter(timports)
			# use theoretical quantiles to prohibit outliers messing things up 
			#tlb = quantile( TDist( p.import_t_df ), 1/p.nimports )*p.import_t_s 
			tlb = quantile( TDist( p.import_t_df ), 0.01 )*p.import_t_s 
			#tlb = quantile(timports,.01) # scatter(timports)
			timports = timports[ (timports .> tlb) .& (timports .< (tlb+maxtime) ) ] # sort(timports)# scatter!(timports)
			nimports = length( timports )
		end
	end
	sort!(timports)
	timports .-= minimum( timports )
	importregions = map(_->sampleimportregion(),1:nimports)
"""

using Pkg
#Pkg.add("JLD2")
#Pkg.add("DataFrames")
#Pkg.add(url="https://github.com/emvolz/NBPMscape", rev="secondary_care_sampling") # commit d8d6e7bc78928336d2aff67fdaccd053e749ee2f at 20:57 19 Dec 2025
using NBPMscape 
using JLD2 
using DataFrames

initialize_parameters();
initialize_parameters("config/outbreak_params_covid19_like_dummy_site_files.yaml");
initialize_parameters("config/outbreak_params_covid19_like.yaml");

NREPS = 1 #500 #1000 #5000
sims_timport_test = map( _-> simforest(NBPMscape.P; initialtime=0.0, maxtime=60.0, maxgenerations=100, max_cases=100000000), 1:NREPS )
G = sims_timport_test[1].G
sims_timport_test_G = [G[ isfinite.(G.tgp) .| isfinite.(G.ted) .| isfinite.(G.thospital) .| isfinite.(G.ticu), : ]]
#@save "covidlike-1.4.1-sims-nrep1.jld2" sims
println(sims_timport_test_G[1][1:10,:])

initialize_parameters("config/outbreak_params_covid19_like.yaml");
sc_tds_timport_test = secondary_care_td(p = NBPMscape.P, sims = sims_timport_test_G )
#1×9 DataFrame
# Row │ sim_n  SC_TD    SC_3TD   n_SC_cases  n_SC_cases_sampled  n_SC_cases_sampled_positive  SC_simid                           tinf_relating_to_SC_TD  last_tinf_relating_to_SC_3TD 
#     │ Any    Any      Any      Any         Any                 Any                          Any                                Any                     Any
#─────┼───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#   1 │ 1      35.7408  49.7408  4607        154                 136                          ["4925ca80-f234-11f0-ab30-13d80a…  29.5636                 40.2519



############################################################################################################
# script copied from 2025_12_HARISS_analysis.jl

#=
Investigating time to detection using metagenomic sampling from a network of hospitals based on existing HARISS network
December 2025/January 2026

## USING UPDATED timport DISTRIBUTION AS ABOVE, i.e. with truncation at 1% quantile ##

=#

using CSV 
using DataFrames
using DifferentialEquations
using Distributions
using GLM
using HypothesisTests
using Interpolations
using JLD2
using JumpProcesses 
using LinearAlgebra
using NamedArrays
using NBPMscape
using Optim
#using Pkg
using Plots
using Plots.PlotMeasures 
using QuadGK
using Random
using RData 
using Revise
import SpecialFunctions as SF 
import StatsBase
using Statistics
using StatsPlots
import UUIDs 
using XLSX

# Simulation run on Imperial College High Performance Computing (HPC) system
# using 'covidlike-1.4.2_w_G_filtered_nrep10_config.jl'
# Code copied below for reference (with some commented out code removed for clarity)
# Each simulation run had 10 replicates but a PBS array (covidlike-1.4.2_w_G_filtered_array100_nrep10.pbs) was used to repeat this up to 100 times
# Due to the size of the files generated, only 33 runs in the array were combined, generating a total of 
# 330 simulation replicates (10 sim replicates x 33 runs in array)

"""
using Pkg
Pkg.add("JLD2")
Pkg.add("DataFrames")
Pkg.add(url="https://github.com/emvolz/NBPMscape", rev="secondary_care_sampling") # commit f5f6638e9ce786e55f37567deeee26af69176b63 at 11:23 16 Jan 2026
using NBPMscape 
using JLD2 
using DataFrames

initialize_parameters("config/outbreak_params_covid19_like.yaml"); # Real site files are not required for simulating the outbreak. They can be replaced when simulating sampling
# dummy site file names also in core.jl as these are checked when installing NBPMscape

NREPS = 10 #50 #500 #1000 #5000
sims = map( _-> simforest(NBPMscape.P; initialtime=0.0, maxtime=90.0, maxgenerations=100, max_cases=100000000), 1:NREPS )
#@save "covidlike-1.4.1-sims-nrep1.jld2" sims

### Also save versions filtered for G and cases accessing healthcare, i.e. GP, ED, hospital and ICU (reduces file size and more likely to be able to reload)

# Create vectors to store filtered dataframes: G filtered for ICU and GP rows
#sims_G_icu_filter = [DataFrame() for _ in 1:length(sims)]
#sims_G_gp_filter = [DataFrame() for _ in 1:length(sims)]
# Saving as a single file avoids some overlap of cases, for example where the individual visits the GP before ICU
sims_G_filtered = [DataFrame() for _ in 1:length(sims)]
# Create df to save the max_cases_df output - useful for checking the impact on the number of infections
sims_max_cases_df = [DataFrame() for _ in 1:length(sims)]
# Create df to save the simid and tinf values for each simulation replicate
sims_simid_tinf_df = [DataFrame() for _ in 1:length(sims)]

# Loop through sim replicates
for s in 1:length(sims)	

	try
		fo = sims[s]
		mc = fo.max_cases_df
		
		G = fo.G

		if size(mc,1) > 0
			sims_max_cases_df[s] = mc
			sims_simid_tinf_df[s] = unique(G[:, [:simid, :tinf]])	
		else
			sims_max_cases_df[s] = missing
			sims_simid_tinf_df[s] = missing
		end

		if size(G,1) > 0
			sims_simid_tinf_df[s] = unique(G[:, [:simid, :tinf]])
		else
			sims_simid_tinf_df[s] = missing
		end
		
	catch err
		@warn "Error in iteration" exception = err
		println("Error in sim number: ", s)
		continue
	end

	try

		fo = sims[s]

		# Filter for G, which is df containing infection information
		# and only retain cases (rows) that progressed to healthcareICU (i.e. capable of detection under primary care, secondary care or ICU sampling methodology)
		G_filtered = fo.G[ isfinite.(fo.G.tgp) .| isfinite.(fo.G.ted) .| isfinite.(fo.G.thospital) .| isfinite.(fo.G.ticu), : ]
		# Trim columns to save on storage space
		G_filtered = G_filtered[:,["pid"
						,"tinf", "tgp", "ted", "thospital", "ticu", "tstepdown", "tdischarge", "trecovered", "tdeceased"
						, "severity", "fatal", "iscommuter"
						, "homeregion"
						#, "commuteregion"
						#, "generation"
						#, "F", "G", "H"
						, "infector_age", "infectee_age"
						, "importedinfection"
						, "simid"
						]]
		
		
		if size(G_filtered,1) > 0
			sims_G_filtered[s] = G_filtered
		else
			sims_G_filtered[s] = missing
		end
	catch err
		@warn "Error in iteration" exception = err
		println("Error in sim number: ", s)
		continue
	end
	
end

@save "covidlike-1.4.2-sims-nrep10_filtered_G.jld2" sims_G_filtered
@save "covidlike-1.4.2-sims-nrep10_max_cases_df.jld2" sims_max_cases_df
@save "covidlike-1.4.2-sims-nrep10_sims_simid_tinf_df.jld2" sims_simid_tinf_df

"""


## Load simulated outbreak data
sims = combine_sim_reps(;  # sim_input_folders in vector format where multiple paths
                          sim_input_folders = [ "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1540261/1540261_filtered_G/" ]
                        , sim_object_name = "sims_G_filtered" #"sims" # "sims_G_icu_filter" # "sims_G_gp_filter" # "sim_G_filtered"
                        , nrep = 10
                        )

# Save to file
@save "covidlike-1.4.2-sims-filtered_G_nrep330_1540261.jld2" sims

# Reload file if necessary
#sims = load("covidlike-1.4.2-sims-filtered_G_nrep330_1540261.jld2", "sims")

# Define path to output folder for analyses files
output_folder = "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1540261_nrep330_analysis/"

## Compute times to detection for secondary care sampling (via HARISS network)
# using different model parameters representing different scenarios
# See config_file_record.xlsx for summary descriptions of scenarios
# Analysis for each set of parameters takes between ~1 and ~3 hours
for i in [1,2,3,4,5,6] # i=1
    initialize_parameters("config/HARISS/covid19_like_params_HARISS_$(i).yaml"); #    initialize_parameters("config/covid19_like_params_HARISS_1.yaml")
    sc_tds = secondary_care_td(sims = sims)
    CSV.write( joinpath( output_folder,"sc_tds_$(i).csv" ), sc_tds )
end
# Started at 00:08 on 20 Jan. Only five done by 09:19 on 20 Jan (1-2h per scenario)

for i in [6,7,13,19,25] # i=1
    initialize_parameters("config/HARISS/covid19_like_params_HARISS_$(i).yaml"); #    initialize_parameters("config/covid19_like_params_HARISS_1.yaml")
    sc_tds = secondary_care_td(sims = sims)
    CSV.write( joinpath( output_folder,"sc_tds_$(i).csv" ), sc_tds )
end
# Started at 00:24 on 21 Jan. 5 completed by 10:02 on 21 Jan. ~2h per run.

for i in [31,37,43,49,55,61,67] # i=1
    initialize_parameters("config/HARISS/covid19_like_params_HARISS_$(i).yaml"); #    initialize_parameters("config/covid19_like_params_HARISS_1.yaml")
    sc_tds = secondary_care_td(sims = sims)
    CSV.write( joinpath( output_folder,"sc_tds_$(i).csv" ), sc_tds )
end
# Started 17:30 on 21 Jan. Completed by 01:56 on 22 Jan. 8.5 / 7 ~ 1h15m per run


# Run analysis for 12 different scenarios (2 seasons x 2 lists of ICU sites x 3 numbers of samples) using ICUs and GPs
for j in [1]#,13,25,37,49,61,73,74,75,76,77,78]
    initialize_parameters("config/HARISS/covid19_like_params_HARISS_$(j).yaml")
    icu_tds = icu_td(sims = sims)
    gp_tds = gp_td(sims = sims)
    CSV.write( joinpath( output_folder,"icu_tds_$(j).csv" ), icu_tds )
    CSV.write( joinpath( output_folder,"gp_tds_$(j).csv" ), gp_tds )
end

# Started at 14:16 finished at 17:54
# Took ~3.5/12 = ~0.30 hrs per run
for j in [13,25,37,49,61,73,74,75,76,77,78]
    initialize_parameters("config/HARISS/covid19_like_params_HARISS_$(j).yaml")
    icu_tds = icu_td(sims = sims)
    gp_tds = gp_td(sims = sims)
    CSV.write( joinpath( output_folder,"icu_tds_$(j).csv" ), icu_tds )
    CSV.write( joinpath( output_folder,"gp_tds_$(j).csv" ), gp_tds )
end
# Started 01:56 on 22 Jan. Completed by 03:28 on 22 Jan. 1.5/11 ~ 8 mins per scenario


## Read in sc_tds files and compute median TD values
using CSV
using DataFrames

# Directory and base file name
path_name = "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1540261_nrep330_analysis/"
base = "sc_tds"
scenario_numbers = [1,2,3,4,5,6,7,13,19,25,31,37,43,49,55,61,67]
# Prepare a container to keep the DataFrames
sc_tds_df_dict = Dict{Int, DataFrame}()
# Prepare a df to hold the median TD values
sc_tds_df = DataFrame(
                    Scenario_number = zeros(Int, length(scenario_numbers))
                    ,sc_median_TD = zeros(Float64, length(scenario_numbers))
                    ,sc_median_3TD = zeros(Float64, length(scenario_numbers))
                    )

# Loop through scenarios, read in files and compute median TD values
for i in 1:length(scenario_numbers) 
    k = scenario_numbers[i]
    # Construct the file name
    file = joinpath(path_name, "$(base)_$(string(k)).csv")

    # Read CSV into a DataFrame and store in the dictionary with the scenario number (numeric suffix) as key
    sc_tds_df_dict[k] = CSV.read(file, DataFrame)
                        
    # Compute times to detection and populate df
    sc_tds_df[i,"Scenario_number"] = k
    sc_tds_df[i,"sc_median_TD"] = median(sc_tds_df_dict[k][:,"SC_TD"])
    sc_tds_df[i,"sc_median_3TD"] = median(sc_tds_df_dict[k][:,"SC_3TD"])
end

# View
println(sc_tds_df)
println(sc_tds_df_dict[1][1:10,:])

# Save
CSV.write( joinpath( output_folder, "sc_median_tds.csv" ), sc_tds_df )


# Compute median TDs for ICU and GP sampling
scenario_numbers_2 = [1,13,25,37,49,61,73,74,75,76,77,78] # Note that scenarios 73-78 are with different ICU sites
# Prepare a container to keep the DataFrames
icu_tds_df_dict = Dict{Int, DataFrame}()
gp_tds_df_dict = Dict{Int, DataFrame}()
# Prepare a df to hold the median TD values
icu_gp_tds_df = DataFrame(
                    Scenario_number = zeros(Int, length(scenario_numbers_2))
                    ,icu_median_TD = zeros(Float64, length(scenario_numbers_2))
                    ,icu_median_3TD = zeros(Float64, length(scenario_numbers_2))
                    ,gp_median_TD = zeros(Float64, length(scenario_numbers_2))
                    ,gp_median_3TD = zeros(Float64, length(scenario_numbers_2))
                    )

# Loop through scenarios, read in files and compute median TD values
for i in 1:length(scenario_numbers_2) 
    k = scenario_numbers_2[i]
    # Construct the file name, pad with leading zeros if needed
    # e.g., sc_tds_01.csv for k=1
    #file = joinpath(path_name, "$(base)$(lpad(string(k), 2, '0'))$(ext)")
    file_icu = joinpath(path_name, "icu_tds_$(string(k)).csv")
    file_gp  = joinpath(path_name,  "gp_tds_$(string(k)).csv")

    # Read CSV into a DataFrame and store in the dictionary with the scenario number (numeric suffix) as key
    icu_tds_df_dict[k] = CSV.read(file_icu, DataFrame)
    gp_tds_df_dict[k] = CSV.read(file_gp, DataFrame)
                        
    # Compute times to detection and populate df
    icu_gp_tds_df[i,"Scenario_number"] = k
    icu_gp_tds_df[i,"icu_median_TD"] = median(icu_tds_df_dict[k][:,"ICU_TD"])
    icu_gp_tds_df[i,"icu_median_3TD"] = median(icu_tds_df_dict[k][:,"ICU_3TD"])
    icu_gp_tds_df[i,"gp_median_TD"] = median(gp_tds_df_dict[k][:,"GP_TD"])
    icu_gp_tds_df[i,"gp_median_3TD"] = median(gp_tds_df_dict[k][:,"GP_3TD"])
end

# View
println(icu_gp_tds_df)

println(icu_tds_df_dict[1][1:10,:])
println(gp_tds_df_dict[1][1:10,:])

# Save
CSV.write( joinpath( output_folder, "icu_gp_median_tds.csv" ), icu_gp_tds_df )



## Plot results

"""
Function        plot_histogram(;sc_data,icu_data,gp_data)
Description     Plots a histogram of the times to detection for 
                simulation replicates for primary care (RCGP),
                secondary care (HARISS) and tertiary care (ICU).
                The median values are also added to the plot.
Arguments       sc_data     Vector containing the times to detection for
                            the secondary care sampling
                icu_data    Vector containing the times to detection for
                            the tertiary care sampling
                gp_data     Vector containing the times to detection for
                            the primary care sampling
Returns         Histogram plot to screen which can later be 
                saved outside of the function
Example         
    
    # Dataframes containing the time to detection results are saved to a dictionary
    # with the key representing the scenario number, which also matches the suffix 
    # on the parameter configuration file
    plot_histogram( sc_data = sc_tds_df_dict[1][:,:SC_TD]
                  , icu_data = icu_tds_df_dict[1][:,:ICU_TD]
                  , gp_data = gp_tds_df_dict[1][:,:GP_TD])
    # Save plot
    Plots.savefig( joinpath( output_folder,"histogram.png") )
"""
function plot_histogram(;sc_data, icu_data_1, icu_data_2, gp_data
                        , xlim,ylim)
    sc_median = median(sc_data)
    icu_median_1 = median(icu_data_1)
    icu_median_2 = median(icu_data_2)
    gp_median = median(gp_data)
    
    p1 = histogram( sc_data, color = :blue, alpha= 0.5, label = "HARISS"
    ,xlimit = xlim
    ,ylimit = ylim
             #, xlabel = "Time since first import to first detection (days)"
             #, ylabel = "Number of simulation replicates"
                )
             #, normalize = :pdf)
    p1 = vline!( [sc_median], color=:blue
                , label="HARISS median"
                , linewidth = 3) 
    # Add tertiary care TD results
    p2 = histogram( icu_data_1, color = :green, alpha =0.5, label="ICU (8 sites)"
                , ylabel = "Number of simulation replicates"
                ,xlimit = xlim
                ,ylimit = ylim)
                #, normalize = :pdf)
    p2 = vline!( [icu_median_1], color=:green
                , label="ICU (8 sites) median"
                , linewidth = 3) 
    # Add tertiary care TD results
    p3 = histogram( icu_data_2, color = :lightgreen, alpha =0.5, label="ICU (29 sites)"
                #, ylabel = "Number of simulation replicates"
                ,xlimit = xlim
                ,ylimit = ylim)
                #, normalize = :pdf)
    p3 = vline!( [icu_median_2], color=:lightgreen
                , label="ICU (29 sites) median"
                , linewidth = 3) 
    
    # Add primary care TD results
    p4 = histogram( gp_data, color=:red, alpha =0.5, label = "RCGP"
                    , xlabel = "Time since first import to first detection (days)"
                    ,xlimit = xlim
                    ,ylimit = ylim)
                #, normalize = :pdf)
    p4 = vline!( [gp_median], color=:red
                , label="RCGP median"
                , linewidth = 3)
    
    
    plot(p1, p2, p3, p4, layout = (4, 1)
        , size = (600,800)) #default size is 600,400
end


# Results for scenario 1
# Winter, 300 samples, Equal weighting, no age targeting
plot_histogram( sc_data = sc_tds_df_dict[1][:,:SC_TD]
              , icu_data_1 = icu_tds_df_dict[1][:,:ICU_TD]
              , icu_data_2 = icu_tds_df_dict[73][:,:ICU_TD]
              , gp_data = gp_tds_df_dict[1][:,:GP_TD]
              ,xlim = [0,125]
              ,ylim = [0,300])
# Save plot
Plots.savefig( joinpath( output_folder, "histogram_1.png" ) )

## Comparing results from varying sample size

# Results for scenario 37
# Summer, 300 samples, Equal weighting, no age targeting
plot_histogram( sc_data = sc_tds_df_dict[37][:,:SC_TD]
              , icu_data_1 = icu_tds_df_dict[37][:,:ICU_TD]
              , icu_data_2 = icu_tds_df_dict[76][:,:ICU_TD]
              , gp_data = gp_tds_df_dict[37][:,:GP_TD]
              ,xlim = [0,125]
              ,ylim = [0,300])
Plots.savefig( joinpath( output_folder, "histogram_37.png" ) )

# Results for scenario 13
# Winter, 500 samples, Equal weighting, no age targeting
plot_histogram( sc_data = sc_tds_df_dict[13][:,:SC_TD]
              , icu_data_1 = icu_tds_df_dict[13][:,:ICU_TD]
              , icu_data_2 = icu_tds_df_dict[74][:,:ICU_TD]
              , gp_data = gp_tds_df_dict[13][:,:GP_TD]
              ,xlim = [0,125]
              ,ylim = [0,300])
Plots.savefig( joinpath( output_folder, "histogram_13.png" ) )

# Results for scenario 25
# Winter, 1000 samples, Equal weighting, no age targeting
plot_histogram( sc_data = sc_tds_df_dict[25][:,:SC_TD]
              , icu_data_1 = icu_tds_df_dict[25][:,:ICU_TD]
              , icu_data_2 = icu_tds_df_dict[75][:,:ICU_TD]
              , gp_data = gp_tds_df_dict[25][:,:GP_TD]
              ,xlim = [0,125]
              ,ylim = [0,300])
Plots.savefig( joinpath( output_folder, "histogram_25.png" ) )

## Comparing results from varying sample age target
# Define data
sc_data_1 = sc_tds_df_dict[1][:,:SC_TD]; sc_median_1 = median(sc_data_1) # Winter, 300 samples, Equal weighting, no age targeting
sc_data_2 = sc_tds_df_dict[2][:,:SC_TD]; sc_median_2 = median(sc_data_2) # Winter, 300 samples, Equal weighting, 100% adult : 0% child (if this split cannot be achieved then sample from the other age category are taken)
sc_data_3 = sc_tds_df_dict[3][:,:SC_TD]; sc_median_3 = median(sc_data_3) # Winter, 300 samples, Equal weighting, 75% adult : 25% child (if this split cannot be achieved then sample from the other age category are taken)
sc_data_4 = sc_tds_df_dict[4][:,:SC_TD]; sc_median_4 = median(sc_data_4) # Winter, 300 samples, Equal weighting, 50% adult : 50% child (if this split cannot be achieved then sample from the other age category are taken)
sc_data_5 = sc_tds_df_dict[5][:,:SC_TD]; sc_median_5 = median(sc_data_5) # Winter, 300 samples, Equal weighting, 25% adult : 75% child (if this split cannot be achieved then sample from the other age category are taken)
sc_data_6 = sc_tds_df_dict[6][:,:SC_TD]; sc_median_6 = median(sc_data_6) # Winter, 300 samples, Equal weighting, 0% adult : 25% child (if this split cannot be achieved then sample from the other age category are taken)
# Plot histograms
xlim = [0,125]
ylim = [0,300]
    
    p1 = histogram( sc_data_1, color = :blue, alpha= 0.5, label = "HARISS - no age target",xlimit = xlim,ylimit = ylim ) #, normalize = :pdf)
    p1 = vline!( [sc_median_1], color=:blue, label="median", linewidth = 3) 
    
    p2 = histogram( sc_data_2, color = :blue, alpha= 0.5, label = "HARISS - 100% adult", xlimit = xlim, ylimit = ylim) #, normalize = :pdf)
    p2 = vline!( [sc_median_2], color=:blue, label="median", linewidth = 3) 
    
    p3 = histogram( sc_data_3, color = :blue, alpha= 0.5, label = "HARISS - 75% adult", xlimit = xlim, ylimit = ylim) #, normalize = :pdf)
    p3 = vline!( [sc_median_3], color=:blue, label="median", linewidth = 3) 
    
    p4 = histogram( sc_data_4, color = :blue, alpha= 0.5, label = "HARISS - 50% adult", xlimit = xlim, ylimit = ylim
                  , ylabel = "Number of simulation replicates") #, normalize = :pdf)
    p4 = vline!( [sc_median_4], color=:blue, label="median", linewidth = 3) 
    
    p5 = histogram( sc_data_5, color = :blue, alpha= 0.5, label = "HARISS - 25% adult", xlimit = xlim, ylimit = ylim) #, normalize = :pdf)
    p5 = vline!( [sc_median_5], color=:blue, label="median", linewidth = 3) 
    
    p6 = histogram( sc_data_6, color = :blue, alpha= 0.5, label = "HARISS - 0% adult", xlimit = xlim, ylimit = ylim
                    , xlabel = "Time since first import to first detection (days)" ) #, normalize = :pdf)
    p6 = vline!( [sc_median_6], color=:blue, label="median", linewidth = 3) 
    
    plot(p1, p2, p3, p4, p5, p6, layout = (6, 1)
        , size = (600,800)) #default size is 600,400

Plots.savefig( joinpath( output_folder, "histogram_1to6.png" ) )

## Comparing results from varying method of allocating samples to HARISS network public health laboratories (PHLs)

# Define data
sc_data_1 = sc_tds_df_dict[1][:,:SC_TD]; sc_median_1 = median(sc_data_1) # Winter, 300 samples, Equal weighting, no age targeting
sc_data_7 = sc_tds_df_dict[7][:,:SC_TD]; sc_median_7 = median(sc_data_7) # Winter, 300 samples, Weighted by ED attendance, no age targeting

sc_data_13 = sc_tds_df_dict[13][:,:SC_TD]; sc_median_13 = median(sc_data_13) # Winter, 500 samples, Equal weighting, no age targeting
sc_data_19 = sc_tds_df_dict[19][:,:SC_TD]; sc_median_19 = median(sc_data_19) # Winter, 500 samples, Weighted by ED attendance, no age targeting

sc_data_25 = sc_tds_df_dict[25][:,:SC_TD]; sc_median_25 = median(sc_data_25) # Winter, 1000 samples, Equal weighting, no age targeting
sc_data_31 = sc_tds_df_dict[31][:,:SC_TD]; sc_median_31 = median(sc_data_31) # Winter, 1000 samples, Weighted by ED attendance, no age targeting
# Plot histograms
xlim = [0,125]
ylim = [0,300]
# 300 samples, equal weight    
p1 = histogram( sc_data_1, color = :blue, alpha= 0.5, label = "HARISS - 300 samples, equal",xlimit = xlim,ylimit = ylim ) #, normalize = :pdf)
p1 = vline!( [sc_median_1], color=:blue, label="median", linewidth = 3)
# 300 samples, weighted
p2 = histogram( sc_data_7, color = :orange, alpha= 0.5, label = "HARISS - 300 samples, weighted", xlimit = xlim, ylimit = ylim)
p2 = vline!( [sc_median_7], color=:orange, label="median", linewidth = 3)
# 500 samples, equal weight
p3 = histogram( sc_data_13, color = :blue, alpha= 0.5, label = "HARISS - 500 samples, equal", xlimit = xlim, ylimit = ylim) #, normalize = :pdf)
p3 = vline!( [sc_median_13], color=:blue, label="median", linewidth = 3) 
# 500 samples, weighted
p4 = histogram( sc_data_19, color = :orange, alpha= 0.5, label = "HARISS - 500 samples, weighted", xlimit = xlim, ylimit = ylim
              , ylabel = "Number of simulation replicates") #, normalize = :pdf)
p4 = vline!( [sc_median_19], color=:orange, label="median", linewidth = 3) 
# 1000 samples, equal weight
p5 = histogram( sc_data_25, color = :blue, alpha= 0.5, label = "HARISS - 1000 samples, equal", xlimit = xlim, ylimit = ylim) #, normalize = :pdf)
p5 = vline!( [sc_median_25], color=:blue, label="median", linewidth = 3) 
# 1000 samples, weighted
p6 = histogram( sc_data_31, color = :orange, alpha= 0.5, label = "HARISS - 1000 samples, weighted", xlimit = xlim, ylimit = ylim
                    , xlabel = "Time since first import to first detection (days)" ) #, normalize = :pdf)
p6 = vline!( [sc_median_31], color=:orange, label="median", linewidth = 3) 
    
plot(p1, p2, p3, p4, p5, p6, layout = (6, 1)
        , size = (600,800)) #default size is 600,400

Plots.savefig( joinpath( output_folder, "histogram_equal_v_weighted.png" ) )

## Plot all median TDs together

# TD (1 case)
# Compute median TDs and 95% CI using bootstrapping for median TDs
sc_td_medians_dict = Dict{Int,NamedTuple}()
for (k,v) in sc_tds_df_dict # v=sc_tds_df_dict[1]
    sc_td_medians_dict[k] = NBPMscape.median_ci_bootstrap( vec = v[:,:SC_TD], n_boot = 1000, alpha = 0.05 ) # [ median, lower, upper ] median_ci_bootstrap function is in misc_functions.jl
end # sc_td_medians_dict[1]

icu_td_medians_dict = Dict{Int,NamedTuple}()
for (k,v) in icu_tds_df_dict # v=sc_tds_df_dict[1]
    icu_td_medians_dict[k] = median_ci_bootstrap( vec = v[:,:ICU_TD], n_boot = 1000, alpha = 0.05 ) # [ median, lower, upper ]
end 

gp_td_medians_dict = Dict{Int,NamedTuple}()
for (k,v) in gp_tds_df_dict # v=sc_tds_df_dict[1]
    gp_td_medians_dict[k] = median_ci_bootstrap( vec = v[:,:GP_TD], n_boot = 1000, alpha = 0.05 ) # [ median, lower, upper ]
end

y_labels = ["HARISS: 300 samples, equal weight, winter, no age target"
           ,"HARISS: 300 samples, equal weight, winter, 100% adult   "
           ,"HARISS: 300 samples, equal weight, winter, 75% adult    "
           ,"HARISS: 300 samples, equal weight, winter, 50% adult    "
           ,"HARISS: 300 samples, equal weight, winter, 25% adult    "
           ,"HARISS: 300 samples, equal weight, winter, 0% adult     "
           ,"HARISS: 300 samples, ED weighted,  winter, no age target"
           ,"ICU 8 sites: 300 samples, winter"
           ,"ICU 29 sites: 300 samples, winter"
           ,"RCGP: 300 samples, winter"
           ,""
           ,"HARISS: 500 samples, equal weight, winter, no age target"
           ,"HARISS: 500 samples, ED weighted, winter, no age target"
           ,"ICU 8 sites: 500 samples, winter"
           ,"ICU 29 sites: 500 samples, winter"
           ,"RCGP: 500 samples, winter"
           ,""
           ,"HARISS: 1000 samples, equal weight, winter, no age target"
           ,"HARISS: 1000 samples, ED weighted, winter, no age target"
           ,"ICU 8 sites: 1000 samples, winter"
           ,"ICU 29 sites: 1000 samples, winter"
           ,"RCGP: 1000 samples, winter"
           ,""
           ,"HARISS: 300 samples, equal weight, summer, no age target"
           ,"HARISS: 300 samples, ED weighted, summer, no age target"
           ,"ICU 8 sites: 300 samples, summer"
           ,"ICU 29 sites: 300 samples, summer"
           ,"RCGP: 300 samples, summer"
           ,""
           ,"HARISS: 500 samples, equal weight, summer, no age target"
           ,"HARISS: 500 samples, ED weighted, summer, no age target"
           ,"ICU 8 sites: 500 samples, summer"
           ,"ICU 29 sites: 500 samples, summer"
           ,"RCGP: 500 samples, summer"
           ,""
           ,"HARISS: 1000 samples, equal weight, summer, no age target"
           ,"HARISS: 1000 samples, ED weighted, summer, no age target"
           ,"ICU 8 sites: 1000 samples, summer"
           ,"ICU 29 sites: 1000 samples, summer"
           ,"RCGP: 1000 samples, summer"
           ]
medians = [sc_td_medians_dict[1].median, sc_td_medians_dict[2].median, sc_td_medians_dict[3].median, sc_td_medians_dict[4].median, sc_td_medians_dict[5].median, sc_td_medians_dict[6].median, sc_td_medians_dict[7].median
          ,icu_td_medians_dict[1].median, icu_td_medians_dict[73].median, gp_td_medians_dict[1].median 
          ,0.0
          ,sc_td_medians_dict[13].median, sc_td_medians_dict[19].median
          ,icu_td_medians_dict[13].median, icu_td_medians_dict[74].median, gp_td_medians_dict[13].median 
          ,0.0
          ,sc_td_medians_dict[25].median, sc_td_medians_dict[31].median
          ,icu_td_medians_dict[25].median, icu_td_medians_dict[75].median, gp_td_medians_dict[25].median 
          ,0.0
          ,sc_td_medians_dict[37].median, sc_td_medians_dict[43].median
          ,icu_td_medians_dict[37].median, icu_td_medians_dict[76].median, gp_td_medians_dict[37].median 
          ,0.0
          ,sc_td_medians_dict[49].median, sc_td_medians_dict[55].median
          ,icu_td_medians_dict[49].median, icu_td_medians_dict[77].median, gp_td_medians_dict[49].median 
          ,0.0
          ,sc_td_medians_dict[61].median, sc_td_medians_dict[67].median
          ,icu_td_medians_dict[61].median, icu_td_medians_dict[78].median, gp_td_medians_dict[61].median 
          ]
lower_errors = medians .- [sc_td_medians_dict[1].lower, sc_td_medians_dict[2].lower, sc_td_medians_dict[3].lower, sc_td_medians_dict[4].lower, sc_td_medians_dict[5].lower, sc_td_medians_dict[6].lower, sc_td_medians_dict[7].lower
                          ,icu_td_medians_dict[1].lower, icu_td_medians_dict[73].lower, gp_td_medians_dict[1].lower 
                          ,0.0
                          ,sc_td_medians_dict[13].lower, sc_td_medians_dict[19].lower
                          ,icu_td_medians_dict[13].lower, icu_td_medians_dict[74].lower, gp_td_medians_dict[13].lower 
                          ,0.0
                          ,sc_td_medians_dict[25].lower, sc_td_medians_dict[31].lower
                          ,icu_td_medians_dict[25].lower, icu_td_medians_dict[75].lower, gp_td_medians_dict[25].lower 
                          ,0.0
                          ,sc_td_medians_dict[37].lower, sc_td_medians_dict[43].lower
                          ,icu_td_medians_dict[37].lower, icu_td_medians_dict[76].lower, gp_td_medians_dict[37].lower 
                          ,0.0
                          ,sc_td_medians_dict[49].lower, sc_td_medians_dict[55].lower
                          ,icu_td_medians_dict[49].lower, icu_td_medians_dict[77].lower, gp_td_medians_dict[49].lower 
                          ,0.0
                          ,sc_td_medians_dict[61].lower, sc_td_medians_dict[67].lower
                          ,icu_td_medians_dict[61].lower, icu_td_medians_dict[78].lower, gp_td_medians_dict[61].lower
                          ]
upper_errors = [sc_td_medians_dict[1].upper, sc_td_medians_dict[2].upper, sc_td_medians_dict[3].upper, sc_td_medians_dict[4].upper, sc_td_medians_dict[5].upper, sc_td_medians_dict[6].upper, sc_td_medians_dict[7].upper
               ,icu_td_medians_dict[1].upper, icu_td_medians_dict[73].upper, gp_td_medians_dict[1].upper 
               ,0.0
               ,sc_td_medians_dict[13].upper, sc_td_medians_dict[19].upper
               ,icu_td_medians_dict[13].upper, icu_td_medians_dict[74].upper, gp_td_medians_dict[13].upper
               ,0.0
               ,sc_td_medians_dict[25].upper, sc_td_medians_dict[31].upper
               ,icu_td_medians_dict[25].upper, icu_td_medians_dict[75].upper, gp_td_medians_dict[25].upper
               ,0.0
               ,sc_td_medians_dict[37].upper, sc_td_medians_dict[43].upper
               ,icu_td_medians_dict[37].upper, icu_td_medians_dict[76].upper, gp_td_medians_dict[37].upper
               ,0.0
               ,sc_td_medians_dict[49].upper, sc_td_medians_dict[55].upper
               ,icu_td_medians_dict[49].upper, icu_td_medians_dict[77].upper, gp_td_medians_dict[49].upper
               ,0.0
               ,sc_td_medians_dict[61].upper, sc_td_medians_dict[67].upper
               ,icu_td_medians_dict[61].upper, icu_td_medians_dict[78].upper, gp_td_medians_dict[61].upper 
               ] .- medians

colors = vcat(fill(:blue, 7), [:green, :lightgreen, :red]
            ,:black
            ,fill(:blue, 2), [:green, :lightgreen, :red]
            ,:black
            ,fill(:blue, 2), [:green, :lightgreen, :red]
            ,:black
            ,fill(:blue, 2), [:green, :lightgreen, :red]
            ,:black
            ,fill(:blue, 2), [:green, :lightgreen, :red]
            ,:black
            ,fill(:blue, 2), [:green, :lightgreen, :red])

# Plot with error bars
# (CI using samples from fitted distribution probably not that useful here)
p_scatter = scatter(
    medians
    ,1:length(y_labels)
    ,xerror = (lower_errors, upper_errors)
    #, markercolor = colors
    #, markerstrokecolor = colors
    , color = colors
    ,seriestype = :scatter
    ,yticks = (1:length(y_labels), y_labels)
    ,xlabel = "Median time to detection (days)\n (with 95% CI via 1000 bootstraps)"
    #,xlims = [40,90] #[0,110]
    ,xlims = [0,90] #[0,110]
    ,xticks = 0:10:100
    ,ylabel = ""
    #,title = "Time to detection\n with different site selection"
    ,legend = false
    ,yflip = true
)
plot(p_scatter, layout = (1, 1)
        , size = (900,600)) #default size is 600,400

Plots.savefig( joinpath( output_folder, "scatter_all_TD.png" ) )

# TD3 (3 cases)
# Compute median TDs and 95% CI using bootstrapping for median TDs
sc_td3_medians_dict = Dict{Int,NamedTuple}()
for (k,v) in sc_tds_df_dict # v=sc_tds_df_dict[1]
    sc_td3_medians_dict[k] = median_ci_bootstrap( vec = v[:,:SC_3TD], n_boot = 1000, alpha = 0.05 ) # [ median, lower, upper ] median_ci_bootstrap function is in misc_functions.jl
end # sc_td_medians_dict[1]

icu_td3_medians_dict = Dict{Int,NamedTuple}()
for (k,v) in icu_tds_df_dict # v=sc_tds_df_dict[1]
    icu_td3_medians_dict[k] = median_ci_bootstrap( vec = v[:,:ICU_3TD], n_boot = 1000, alpha = 0.05 ) # [ median, lower, upper ]
end 

gp_td3_medians_dict = Dict{Int,NamedTuple}()
for (k,v) in gp_tds_df_dict # v=sc_tds_df_dict[1]
    gp_td3_medians_dict[k] = median_ci_bootstrap( vec = v[:,:GP_3TD], n_boot = 1000, alpha = 0.05 ) # [ median, lower, upper ]
end

medians_td3 = [sc_td3_medians_dict[1].median, sc_td3_medians_dict[2].median, sc_td3_medians_dict[3].median, sc_td3_medians_dict[4].median, sc_td3_medians_dict[5].median, sc_td3_medians_dict[6].median, sc_td3_medians_dict[7].median
                ,icu_td3_medians_dict[1].median, icu_td3_medians_dict[73].median, gp_td3_medians_dict[1].median 
                ,0.0
                ,sc_td3_medians_dict[13].median, sc_td3_medians_dict[19].median
                ,icu_td3_medians_dict[13].median, icu_td3_medians_dict[74].median, gp_td3_medians_dict[13].median 
                ,0.0
                ,sc_td3_medians_dict[25].median, sc_td3_medians_dict[31].median
                ,icu_td3_medians_dict[25].median, icu_td3_medians_dict[75].median, gp_td3_medians_dict[25].median 
                ,0.0
                ,sc_td3_medians_dict[37].median, sc_td3_medians_dict[43].median
                ,icu_td3_medians_dict[37].median, icu_td3_medians_dict[76].median, gp_td3_medians_dict[37].median 
                ,0.0
                ,sc_td3_medians_dict[49].median, sc_td3_medians_dict[55].median
                ,icu_td3_medians_dict[49].median, icu_td3_medians_dict[77].median, gp_td3_medians_dict[49].median 
                ,0.0
                ,sc_td3_medians_dict[61].median, sc_td3_medians_dict[67].median
                ,icu_td3_medians_dict[61].median, icu_td3_medians_dict[78].median, gp_td3_medians_dict[61].median 
                ]
lower_errors_td3 = medians_td3 .- [sc_td3_medians_dict[1].lower, sc_td3_medians_dict[2].lower, sc_td3_medians_dict[3].lower, sc_td3_medians_dict[4].lower, sc_td3_medians_dict[5].lower, sc_td3_medians_dict[6].lower, sc_td3_medians_dict[7].lower
                                    ,icu_td3_medians_dict[1].lower, icu_td3_medians_dict[73].lower, gp_td3_medians_dict[1].lower 
                                    ,0.0
                                    ,sc_td3_medians_dict[13].lower, sc_td3_medians_dict[19].lower
                                    ,icu_td3_medians_dict[13].lower, icu_td3_medians_dict[74].lower, gp_td3_medians_dict[13].lower 
                                    ,0.0
                                    ,sc_td3_medians_dict[25].lower, sc_td3_medians_dict[31].lower
                                    ,icu_td3_medians_dict[25].lower, icu_td3_medians_dict[75].lower, gp_td3_medians_dict[25].lower 
                                    ,0.0
                                    ,sc_td3_medians_dict[37].lower, sc_td3_medians_dict[43].lower
                                    ,icu_td3_medians_dict[37].lower, icu_td3_medians_dict[76].lower, gp_td3_medians_dict[37].lower 
                                    ,0.0
                                    ,sc_td3_medians_dict[49].lower, sc_td3_medians_dict[55].lower
                                    ,icu_td3_medians_dict[49].lower, icu_td3_medians_dict[77].lower, gp_td3_medians_dict[49].lower 
                                    ,0.0
                                    ,sc_td3_medians_dict[61].lower, sc_td3_medians_dict[67].lower
                                    ,icu_td3_medians_dict[61].lower, icu_td3_medians_dict[78].lower, gp_td3_medians_dict[61].lower
                                    ]
upper_errors_td3 = [sc_td3_medians_dict[1].upper, sc_td3_medians_dict[2].upper, sc_td3_medians_dict[3].upper, sc_td3_medians_dict[4].upper, sc_td3_medians_dict[5].upper, sc_td3_medians_dict[6].upper, sc_td3_medians_dict[7].upper
                    ,icu_td3_medians_dict[1].upper, icu_td3_medians_dict[73].upper, gp_td3_medians_dict[1].upper 
                    ,0.0
                    ,sc_td3_medians_dict[13].upper, sc_td3_medians_dict[19].upper
                    ,icu_td3_medians_dict[13].upper, icu_td3_medians_dict[74].upper, gp_td3_medians_dict[13].upper
                    ,0.0
                    ,sc_td3_medians_dict[25].upper, sc_td3_medians_dict[31].upper
                    ,icu_td3_medians_dict[25].upper, icu_td3_medians_dict[75].upper, gp_td3_medians_dict[25].upper
                    ,0.0
                    ,sc_td3_medians_dict[37].upper, sc_td3_medians_dict[43].upper
                    ,icu_td3_medians_dict[37].upper, icu_td3_medians_dict[76].upper, gp_td3_medians_dict[37].upper
                    ,0.0
                    ,sc_td3_medians_dict[49].upper, sc_td3_medians_dict[55].upper
                    ,icu_td3_medians_dict[49].upper, icu_td3_medians_dict[77].upper, gp_td3_medians_dict[49].upper
                    ,0.0
                    ,sc_td3_medians_dict[61].upper, sc_td3_medians_dict[67].upper
                    ,icu_td3_medians_dict[61].upper, icu_td3_medians_dict[78].upper, gp_td3_medians_dict[61].upper 
                    ] .- medians_td3

# Plot with error bars
# (CI using samples from fitted distribution probably not that useful here)
p_scatter_td3 = scatter(
    medians_td3
    ,1:length(y_labels)
    ,xerror = (lower_errors_td3, upper_errors_td3)
    #, markercolor = colors
    #, markerstrokecolor = colors
    , color = colors
    ,seriestype = :scatter
    ,yticks = (1:length(y_labels), y_labels)
    ,xlabel = "Median time to detect 3 cases (days)\n (with 95% CI via 1000 bootstraps)"
    #,xlims = [40,90] #[0,110]
    ,xlims = [0,90] #[0,110]
    ,xticks = 0:10:100
    ,ylabel = ""
    #,title = "Time to detection\n with different site selection"
    ,legend = false
    ,yflip = true
)
plot(p_scatter_td3, layout = (1, 1)
        , size = (900,600)) #default size is 600,400

Plots.savefig( joinpath( output_folder, "scatter_all_TD3.png" ) )



###################################################################################################################

## Computing and plotting combined times to detection

# Read in files 

using CSV
using DataFrames

# Directory and base file name
path_name = "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1540261_nrep330_analysis/"

# Match scenario numbers (from config file suffixes) between sampling setting
scenario_numbers_sc =     [1,2,3,4,5,6,7,13,19,25,31,37,43,49,55,61,67]
scenario_numbers_icu_8 =  [1,1,1,1,1,1,1,13,13,25,25,37,37,49,49,61,61]
scenario_numbers_icu_29 = [73,73,73,73,73,73,73,74,74,75,75,76,76,77,77,78,78]
scenario_numbers_gp =     [1,1,1,1,1,1,1,13,13,25,25,37,37,49,49,61,61]

# Prepare containers to store the individual DataFrames
sc_tds_df_dict     = Dict{Int, DataFrame}()
icu_8_tds_df_dict  = Dict{Int, DataFrame}()
icu_29_tds_df_dict = Dict{Int, DataFrame}()
gp_tds_df_dict     = Dict{Int, DataFrame}()

# Prepare a df to hold the median TD values and upper and lower 95% CI from bootstrapping
default_int_col   = zeros(Int,     length(scenario_numbers_sc))
default_float_col = zeros(Float64, length(scenario_numbers_sc))
tds_df = DataFrame(
                    Scenario_number = default_int_col
                    ,Scenario_number_icu_29 = default_int_col # ICU sampling at 29 sites uses different config files because it requires a different parameter value for the ICU sites
                    ## Median
                    # TD
                    ,sc_median_TD = default_float_col
                    ,icu_8_median_TD = default_float_col
                    ,icu_29_median_TD = default_float_col
                    ,gp_median_TD = default_float_col
                    # 3TD
                    ,sc_median_3TD = default_float_col
                    ,icu_8_median_3TD = default_float_col
                    ,icu_29_median_3TD = default_float_col
                    ,gp_median_3TD = default_float_col
                    # Combined TD
                    # SC + ICU 8
                    ,sc_icu_8_comb_median_TD = default_float_col
                    # SC + ICU 29
                    ,sc_icu_29_comb_median_TD = default_float_col
                    # SC + GP
                    ,sc_gp_comb_median_TD = default_float_col
                    # SC + ICU 8 + GP
                    ,sc_icu_8_gp_comb_median_TD = default_float_col
                    # SC + ICU 29 + GP
                    ,sc_icu_29_gp_comb_median_TD = default_float_col
                    # Combined 3TD
                    # SC + ICU 8
                    ,sc_icu_8_comb_median_3TD = default_float_col
                    # SC + ICU 29
                    ,sc_icu_29_comb_median_3TD = default_float_col
                    # SC + GP
                    ,sc_gp_comb_median_3TD = default_float_col
                    # SC + ICU 8 + GP
                    ,sc_icu_8_gp_comb_median_3TD = default_float_col
                    # SC + ICU 29 + GP
                    ,sc_icu_29_gp_comb_median_3TD = default_float_col
                    ## Lower CI
                    # TD
                    ,sc_median_TD_lower = default_float_col
                    ,icu_8_median_TD_lower = default_float_col
                    ,icu_29_median_TD_lower = default_float_col
                    ,gp_median_TD_lower = default_float_col
                    # 3TD
                    ,sc_median_3TD_lower = default_float_col
                    ,icu_8_median_3TD_lower = default_float_col
                    ,icu_29_median_3TD_lower = default_float_col
                    ,gp_median_3TD_lower = default_float_col
                    # Combined TD
                    # SC + ICU 8
                    ,sc_icu_8_comb_median_TD_lower = default_float_col
                    # SC + ICU 29
                    ,sc_icu_29_comb_median_TD_lower = default_float_col
                    # SC + GP
                    ,sc_gp_comb_median_TD_lower = default_float_col
                    # SC + ICU 8 + GP
                    ,sc_icu_8_gp_comb_median_TD_lower = default_float_col
                    # SC + ICU 29 + GP
                    ,sc_icu_29_gp_comb_median_TD_lower = default_float_col
                    # Combined 3TD
                    # SC + ICU 8
                    ,sc_icu_8_comb_median_3TD_lower = default_float_col
                    # SC + ICU 29
                    ,sc_icu_29_comb_median_3TD_lower = default_float_col
                    # SC + GP
                    ,sc_gp_comb_median_3TD_lower = default_float_col
                    # SC + ICU 8 + GP
                    ,sc_icu_8_gp_comb_median_3TD_lower = default_float_col
                    # SC + ICU 29 + GP
                    ,sc_icu_29_gp_comb_median_3TD_lower = default_float_col
                    ## Upper CI
                    # TD
                    ,sc_median_TD_upper = default_float_col
                    ,icu_8_median_TD_upper = default_float_col
                    ,icu_29_median_TD_upper = default_float_col
                    ,gp_median_TD_upper = default_float_col
                    # 3TD
                    ,sc_median_3TD_upper = default_float_col
                    ,icu_8_median_3TD_upper = default_float_col
                    ,icu_29_median_3TD_upper = default_float_col
                    ,gp_median_3TD_upper = default_float_col
                    # Combined TD
                    # SC + ICU 8
                    ,sc_icu_8_comb_median_TD_upper = default_float_col
                    # SC + ICU 29
                    ,sc_icu_29_comb_median_TD_upper = default_float_col
                    # SC + GP
                    ,sc_gp_comb_median_TD_upper = default_float_col
                    # SC + ICU 8 + GP
                    ,sc_icu_8_gp_comb_median_TD_upper = default_float_col
                    # SC + ICU 29 + GP
                    ,sc_icu_29_gp_comb_median_TD_upper = default_float_col
                    # Combined 3TD
                    # SC + ICU 8
                    ,sc_icu_8_comb_median_3TD_upper = default_float_col
                    # SC + ICU 29
                    ,sc_icu_29_comb_median_3TD_upper = default_float_col
                    # SC + GP
                    ,sc_gp_comb_median_3TD_upper = default_float_col
                    # SC + ICU 8 + GP
                    ,sc_icu_8_gp_comb_median_3TD_upper = default_float_col
                    # SC + ICU 29 + GP
                    ,sc_icu_29_gp_comb_median_3TD_upper = default_float_col
                    )


# Loop through scenarios, read in files and compute median TD values
# Secondary care
for i in 1:length(scenario_numbers_sc) 
    j = scenario_numbers_sc[i]
    # Construct the file name, pad with leading zeros if needed
    # e.g., sc_tds_01.csv for k=1
    #file = joinpath(path_name, "$(base)$(lpad(string(k), 2, '0'))$(ext)")
    file = joinpath(path_name, "sc_tds_$(string(j)).csv")

    # Read CSV into a DataFrame and store in the dictionary with the scenario number (numeric suffix) as key
    sc_tds_df_dict[j] = CSV.read(file, DataFrame)
                        
    # Compute times to detection and populate df
    tds_df[i,"Scenario_number"] = j
    tds_df[i,"sc_median_TD"] = median(sc_tds_df_dict[j][:,"SC_TD"])
    tds_df[i,"sc_median_3TD"] = median(sc_tds_df_dict[j][:,"SC_3TD"])
end

# ICU 8 sites
for i in 1:length(scenario_numbers_icu_8) 
    k = scenario_numbers_icu_8[i]
    # Construct the file name, pad with leading zeros if needed
    # e.g., sc_tds_01.csv for k=1
    #file = joinpath(path_name, "$(base)$(lpad(string(k), 2, '0'))$(ext)")
    file = joinpath(path_name, "icu_tds_$(string(k)).csv")

    # Read CSV into a DataFrame and store in the dictionary with the scenario number (numeric suffix) as key
    icu_8_tds_df_dict[k] = CSV.read(file, DataFrame)
                        
    # Compute times to detection and populate df
    tds_df[i,"icu_8_median_TD"] = median(icu_8_tds_df_dict[k][:,"ICU_TD"])
    tds_df[i,"icu_8_median_3TD"] = median(icu_8_tds_df_dict[k][:,"ICU_3TD"])
end

# ICU 29 sites
for i in 1:length(scenario_numbers_icu_29) 
    m = scenario_numbers_icu_29[i]
    # Construct the file name, pad with leading zeros if needed
    # e.g., sc_tds_01.csv for k=1
    #file = joinpath(path_name, "$(base)$(lpad(string(k), 2, '0'))$(ext)")
    file = joinpath(path_name, "icu_tds_$(string(m)).csv")

    # Read CSV into a DataFrame and store in the dictionary with the scenario number (numeric suffix) as key
    icu_29_tds_df_dict[m] = CSV.read(file, DataFrame)
                        
    # Compute times to detection and populate df
    tds_df[i,"Scenario_number_icu_29"] = m
    tds_df[i,"icu_29_median_TD"] = median(icu_29_tds_df_dict[m][:,"ICU_TD"])
    tds_df[i,"icu_29_median_3TD"] = median(icu_29_tds_df_dict[m][:,"ICU_3TD"])
end

# GP
for i in 1:length(scenario_numbers_gp) 
    n = scenario_numbers_gp[i]
    # Construct the file name, pad with leading zeros if needed
    # e.g., sc_tds_01.csv for k=1
    #file = joinpath(path_name, "$(base)$(lpad(string(k), 2, '0'))$(ext)")
    file = joinpath(path_name, "gp_tds_$(string(n)).csv")

    # Read CSV into a DataFrame and store in the dictionary with the scenario number (numeric suffix) as key
    gp_tds_df_dict[n] = CSV.read(file, DataFrame)
                        
    # Compute times to detection and populate df
    tds_df[i,"gp_median_TD"]  = median(gp_tds_df_dict[n][:,"GP_TD"])
    tds_df[i,"gp_median_3TD"] = median(gp_tds_df_dict[n][:,"GP_3TD"])
end

#println(tds_df)

#println(sc_tds_df_dict[1][1:10,:])

# Compute median TDs combining different sampling settings
for i in 1:length(scenario_numbers_sc) # i =1
    j = scenario_numbers_sc[i]
    k = scenario_numbers_icu_8[i]
    m = scenario_numbers_icu_29[i]
    n = scenario_numbers_gp[i]

    # Populate scenario numbers
    tds_df[i,"Scenario_number"] = j # Assumes that j == k == n
    tds_df[i,"Scenario_number_icu_29"] = m

    # Times to detection (for each simulation replicate) for the different sample settings and for the specific parameter set
    # TD
    sc_tds     =     sc_tds_df_dict[ j ][:,"SC_TD" ]
    icu_8_tds  =  icu_8_tds_df_dict[ k ][:,"ICU_TD"]
    icu_29_tds = icu_29_tds_df_dict[ m ][:,"ICU_TD"]
    gp_tds     =     gp_tds_df_dict[ n ][:,"GP_TD" ]
    #println(sample_setting_1[1:10],sample_setting_2[1:10],combined_min_TD[1:10])
    # 3TD
    sc_tds_3     =     sc_tds_df_dict[j][:,"SC_3TD" ]
    icu_8_tds_3  =  icu_8_tds_df_dict[k][:,"ICU_3TD"]
    icu_29_tds_3 = icu_29_tds_df_dict[m][:,"ICU_3TD"]
    gp_tds_3     =     gp_tds_df_dict[n][:,"GP_3TD" ]

    # Compute median times to detection and 95% CI using bootstrap
    # SC
    (tds_df[i,"sc_median_TD"], tds_df[i,"sc_median_TD_lower"], tds_df[i,"sc_median_TD_upper"]) = NBPMscape.median_ci_bootstrap( vec = sc_tds, n_boot = 1000, alpha = 0.05)
    (tds_df[i,"sc_median_3TD"], tds_df[i,"sc_median_3TD_lower"], tds_df[i,"sc_median_3TD_upper"]) = NBPMscape.median_ci_bootstrap( vec = sc_tds_3, n_boot = 1000, alpha = 0.05)
    # ICU 8 sites
    (tds_df[i,"icu_8_median_TD"], tds_df[i,"icu_8_median_TD_lower"], tds_df[i,"icu_8_median_TD_upper"]) = NBPMscape.median_ci_bootstrap( vec = icu_8_tds, n_boot = 1000, alpha = 0.05)
    (tds_df[i,"icu_8_median_3TD"], tds_df[i,"icu_8_median_3TD_lower"], tds_df[i,"icu_8_median_3TD_upper"]) = NBPMscape.median_ci_bootstrap( vec = icu_8_tds_3, n_boot = 1000, alpha = 0.05)
    # ICU 29 sites
    (tds_df[i,"icu_29_median_TD"], tds_df[i,"icu_29_median_TD_lower"], tds_df[i,"icu_29_median_TD_upper"]) = NBPMscape.median_ci_bootstrap( vec = icu_29_tds, n_boot = 1000, alpha = 0.05)
    (tds_df[i,"icu_29_median_3TD"], tds_df[i,"icu_29_median_3TD_lower"], tds_df[i,"icu_29_median_3TD_upper"]) = NBPMscape.median_ci_bootstrap( vec = icu_29_tds_3, n_boot = 1000, alpha = 0.05)
    # GP
    (tds_df[i,"gp_median_TD"], tds_df[i,"gp_median_TD_lower"], tds_df[i,"gp_median_TD_upper"]) = NBPMscape.median_ci_bootstrap( vec = gp_tds, n_boot = 1000, alpha = 0.05)
    (tds_df[i,"gp_median_3TD"], tds_df[i,"gp_median_3TD_lower"], tds_df[i,"gp_median_3TD_upper"]) = NBPMscape.median_ci_bootstrap( vec = gp_tds_3, n_boot = 1000, alpha = 0.05)
    

    # Compute combined times to detection and populate df
    # SC and ICU 8 sites
    sc_icu_8_comb = min.(sc_tds,   icu_8_tds  )
    sc_icu_8_comb_3 = min.(sc_tds_3,   icu_8_tds_3  )

    (tds_df[i,"sc_icu_8_comb_median_TD"], tds_df[i,"sc_icu_8_comb_median_TD_lower"], tds_df[i,"sc_icu_8_comb_median_TD_upper"]) = NBPMscape.median_ci_bootstrap( vec = sc_icu_8_comb, n_boot = 1000, alpha = 0.05)
    (tds_df[i,"sc_icu_8_comb_median_3TD"], tds_df[i,"sc_icu_8_comb_median_3TD_lower"], tds_df[i,"sc_icu_8_comb_median_3TD_upper"]) = NBPMscape.median_ci_bootstrap( vec = sc_icu_8_comb_3, n_boot = 1000, alpha = 0.05)

    # SC and ICU 29 sites
    sc_icu_29_comb = min.(sc_tds,   icu_29_tds  )
    sc_icu_29_comb_3 = min.(sc_tds_3,   icu_29_tds_3  )

    (tds_df[i,"sc_icu_29_comb_median_TD"], tds_df[i,"sc_icu_29_comb_median_TD_lower"], tds_df[i,"sc_icu_29_comb_median_TD_upper"]) = NBPMscape.median_ci_bootstrap( vec = sc_icu_29_comb, n_boot = 1000, alpha = 0.05)
    (tds_df[i,"sc_icu_29_comb_median_3TD"], tds_df[i,"sc_icu_29_comb_median_3TD_lower"], tds_df[i,"sc_icu_29_comb_median_3TD_upper"]) = NBPMscape.median_ci_bootstrap( vec = sc_icu_29_comb_3, n_boot = 1000, alpha = 0.05)

    # SC and GP
    sc_gp_comb = min.(sc_tds,   gp_tds  )
    sc_gp_comb_3 = min.(sc_tds_3,   gp_tds_3  )

    (tds_df[i,"sc_gp_comb_median_TD"], tds_df[i,"sc_gp_comb_median_TD_lower"], tds_df[i,"sc_gp_comb_median_TD_upper"]) = NBPMscape.median_ci_bootstrap( vec = sc_gp_comb, n_boot = 1000, alpha = 0.05)
    (tds_df[i,"sc_gp_comb_median_3TD"], tds_df[i,"sc_gp_comb_median_3TD_lower"], tds_df[i,"sc_gp_comb_median_3TD_upper"]) = NBPMscape.median_ci_bootstrap( vec = sc_gp_comb_3, n_boot = 1000, alpha = 0.05)

    # SC and ICU 8 sites and GP
    sc_icu_8_gp_comb = min.(sc_tds,   icu_8_tds, gp_tds )
    sc_icu_8_gp_comb_3 = min.(sc_tds_3, icu_8_tds_3,   gp_tds_3  )

    (tds_df[i,"sc_icu_8_gp_comb_median_TD"], tds_df[i,"sc_icu_8_gp_comb_median_TD_lower"], tds_df[i,"sc_icu_8_gp_comb_median_TD_upper"]) = NBPMscape.median_ci_bootstrap( vec = sc_icu_8_gp_comb, n_boot = 1000, alpha = 0.05)
    (tds_df[i,"sc_icu_8_gp_comb_median_3TD"], tds_df[i,"sc_icu_8_gp_comb_median_3TD_lower"], tds_df[i,"sc_icu_8_gp_comb_median_3TD_upper"]) = NBPMscape.median_ci_bootstrap( vec = sc_icu_8_gp_comb_3, n_boot = 1000, alpha = 0.05)

    # SC and ICU 29 sites and GP
    sc_icu_29_gp_comb = min.(sc_tds,   icu_29_tds, gp_tds )
    sc_icu_29_gp_comb_3 = min.(sc_tds_3, icu_29_tds_3,   gp_tds_3  )

    (tds_df[i,"sc_icu_29_gp_comb_median_TD"], tds_df[i,"sc_icu_29_gp_comb_median_TD_lower"], tds_df[i,"sc_icu_29_gp_comb_median_TD_upper"]) = NBPMscape.median_ci_bootstrap( vec = sc_icu_29_gp_comb, n_boot = 1000, alpha = 0.05)
    (tds_df[i,"sc_icu_29_gp_comb_median_3TD"], tds_df[i,"sc_icu_29_gp_comb_median_3TD_lower"], tds_df[i,"sc_icu_29_gp_comb_median_3TD_upper"]) = NBPMscape.median_ci_bootstrap( vec = sc_icu_29_gp_comb_3, n_boot = 1000, alpha = 0.05)

end

#CSV.write( joinpath( output_folder, "tds_df_inc_comb.csv" ), tds_df )
CSV.write( joinpath( output_folder, "tds_df_inc_comb_inc_CI.csv") , tds_df ) 

## Plot all median TDs together

# TD (1 case)
# Create base labels and repeat 6 times (representing the 300, 500, 1000 samples in winter and summer), then drop the last "" (to match length)
y_label_base = [
                "HARISS"
                ,"HARISS & ICU 8 sites"
                ,"HARISS & ICU 29 sites"
                ,"HARISS & RCGP"
                ,"HARISS & ICU 8 sites & RCGP"
                ,"HARISS & ICU 29 sites & RCGP"
                ,""
            ]
y_labels = vcat(repeat(y_label_base, 6)...)[1:end-1]

medians = [ # Scenario / config file 1
            tds_df[1,:sc_median_TD], tds_df[1,:sc_icu_8_comb_median_TD], tds_df[1,:sc_icu_29_comb_median_TD], tds_df[1,:sc_gp_comb_median_TD], tds_df[1,:sc_icu_8_gp_comb_median_TD], tds_df[1,:sc_icu_29_gp_comb_median_TD]
          ,0.0
          # Scenario / config file 13
          ,tds_df[8,:sc_median_TD], tds_df[8,:sc_icu_8_comb_median_TD], tds_df[8,:sc_icu_29_comb_median_TD], tds_df[8,:sc_gp_comb_median_TD], tds_df[8,:sc_icu_8_gp_comb_median_TD], tds_df[8,:sc_icu_29_gp_comb_median_TD]
          ,0.0
          # Scenario / config file 25
          ,tds_df[10,:sc_median_TD], tds_df[10,:sc_icu_8_comb_median_TD], tds_df[10,:sc_icu_29_comb_median_TD], tds_df[10,:sc_gp_comb_median_TD], tds_df[10,:sc_icu_8_gp_comb_median_TD], tds_df[10,:sc_icu_29_gp_comb_median_TD]
          ,0.0
          # Scenario / config file 37
          ,tds_df[12,:sc_median_TD], tds_df[12,:sc_icu_8_comb_median_TD], tds_df[12,:sc_icu_29_comb_median_TD], tds_df[12,:sc_gp_comb_median_TD], tds_df[12,:sc_icu_8_gp_comb_median_TD], tds_df[12,:sc_icu_29_gp_comb_median_TD]
          ,0.0
          # Scenario / config file 49
          ,tds_df[14,:sc_median_TD], tds_df[14,:sc_icu_8_comb_median_TD], tds_df[14,:sc_icu_29_comb_median_TD], tds_df[14,:sc_gp_comb_median_TD], tds_df[14,:sc_icu_8_gp_comb_median_TD], tds_df[14,:sc_icu_29_gp_comb_median_TD]
          ,0.0
          # Scenario / config file 61
          ,tds_df[16,:sc_median_TD], tds_df[16,:sc_icu_8_comb_median_TD], tds_df[16,:sc_icu_29_comb_median_TD], tds_df[16,:sc_gp_comb_median_TD], tds_df[16,:sc_icu_8_gp_comb_median_TD], tds_df[16,:sc_icu_29_gp_comb_median_TD]
          ]
lower_errors = medians .- [ # Scenario / config file 1
                            tds_df[1,:sc_median_TD_lower], tds_df[1,:sc_icu_8_comb_median_TD_lower], tds_df[1,:sc_icu_29_comb_median_TD_lower], tds_df[1,:sc_gp_comb_median_TD_lower], tds_df[1,:sc_icu_8_gp_comb_median_TD_lower], tds_df[1,:sc_icu_29_gp_comb_median_TD_lower]
                            ,0.0
                            # Scenario / config file 13
                            ,tds_df[8,:sc_median_TD_lower], tds_df[8,:sc_icu_8_comb_median_TD_lower], tds_df[8,:sc_icu_29_comb_median_TD_lower], tds_df[8,:sc_gp_comb_median_TD_lower], tds_df[8,:sc_icu_8_gp_comb_median_TD_lower], tds_df[8,:sc_icu_29_gp_comb_median_TD_lower]
                            ,0.0
                            # Scenario / config file 25
                            ,tds_df[10,:sc_median_TD_lower], tds_df[10,:sc_icu_8_comb_median_TD_lower], tds_df[10,:sc_icu_29_comb_median_TD_lower], tds_df[10,:sc_gp_comb_median_TD_lower], tds_df[10,:sc_icu_8_gp_comb_median_TD_lower], tds_df[10,:sc_icu_29_gp_comb_median_TD_lower]
                            ,0.0
                            # Scenario / config file 37
                            ,tds_df[12,:sc_median_TD_lower], tds_df[12,:sc_icu_8_comb_median_TD_lower], tds_df[12,:sc_icu_29_comb_median_TD_lower], tds_df[12,:sc_gp_comb_median_TD_lower], tds_df[12,:sc_icu_8_gp_comb_median_TD_lower], tds_df[12,:sc_icu_29_gp_comb_median_TD_lower]
                            ,0.0
                            # Scenario / config file 49
                            ,tds_df[14,:sc_median_TD_lower], tds_df[14,:sc_icu_8_comb_median_TD_lower], tds_df[14,:sc_icu_29_comb_median_TD_lower], tds_df[14,:sc_gp_comb_median_TD_lower], tds_df[14,:sc_icu_8_gp_comb_median_TD_lower], tds_df[14,:sc_icu_29_gp_comb_median_TD_lower]
                            ,0.0
                            # Scenario / config file 61
                            ,tds_df[16,:sc_median_TD_lower], tds_df[16,:sc_icu_8_comb_median_TD_lower], tds_df[16,:sc_icu_29_comb_median_TD_lower], tds_df[16,:sc_gp_comb_median_TD_lower], tds_df[16,:sc_icu_8_gp_comb_median_TD_lower], tds_df[16,:sc_icu_29_gp_comb_median_TD_lower]
                            ]
upper_errors = [ # Scenario / config file 1
                            tds_df[1,:sc_median_TD_upper], tds_df[1,:sc_icu_8_comb_median_TD_upper], tds_df[1,:sc_icu_29_comb_median_TD_upper], tds_df[1,:sc_gp_comb_median_TD_upper], tds_df[1,:sc_icu_8_gp_comb_median_TD_upper], tds_df[1,:sc_icu_29_gp_comb_median_TD_upper]
                            ,0.0
                            # Scenario / config file 13
                            ,tds_df[8,:sc_median_TD_upper], tds_df[8,:sc_icu_8_comb_median_TD_upper], tds_df[8,:sc_icu_29_comb_median_TD_upper], tds_df[8,:sc_gp_comb_median_TD_upper], tds_df[8,:sc_icu_8_gp_comb_median_TD_upper], tds_df[8,:sc_icu_29_gp_comb_median_TD_upper]
                            ,0.0
                            # Scenario / config file 25
                            ,tds_df[10,:sc_median_TD_upper], tds_df[10,:sc_icu_8_comb_median_TD_upper], tds_df[10,:sc_icu_29_comb_median_TD_upper], tds_df[10,:sc_gp_comb_median_TD_upper], tds_df[10,:sc_icu_8_gp_comb_median_TD_upper], tds_df[10,:sc_icu_29_gp_comb_median_TD_upper]
                            ,0.0
                            # Scenario / config file 37
                            ,tds_df[12,:sc_median_TD_upper], tds_df[12,:sc_icu_8_comb_median_TD_upper], tds_df[12,:sc_icu_29_comb_median_TD_upper], tds_df[12,:sc_gp_comb_median_TD_upper], tds_df[12,:sc_icu_8_gp_comb_median_TD_upper], tds_df[12,:sc_icu_29_gp_comb_median_TD_upper]
                            ,0.0
                            # Scenario / config file 49
                            ,tds_df[14,:sc_median_TD_upper], tds_df[14,:sc_icu_8_comb_median_TD_upper], tds_df[14,:sc_icu_29_comb_median_TD_upper], tds_df[14,:sc_gp_comb_median_TD_upper], tds_df[14,:sc_icu_8_gp_comb_median_TD_upper], tds_df[14,:sc_icu_29_gp_comb_median_TD_upper]
                            ,0.0
                            # Scenario / config file 61
                            ,tds_df[16,:sc_median_TD_upper], tds_df[16,:sc_icu_8_comb_median_TD_upper], tds_df[16,:sc_icu_29_comb_median_TD_upper], tds_df[16,:sc_gp_comb_median_TD_upper], tds_df[16,:sc_icu_8_gp_comb_median_TD_upper], tds_df[16,:sc_icu_29_gp_comb_median_TD_upper]
                            ] .- medians

colors = vcat(fill(:blue, 7), [:green, :lightgreen, :red]
            ,:black
            ,fill(:blue, 2), [:green, :lightgreen, :red]
            ,:black
            ,fill(:blue, 2), [:green, :lightgreen, :red]
            ,:black
            ,fill(:blue, 2), [:green, :lightgreen, :red]
            ,:black
            ,fill(:blue, 2), [:green, :lightgreen, :red]
            ,:black
            ,fill(:blue, 2), [:green, :lightgreen, :red])

# Plot with error bars
# (CI using samples from fitted distribution probably not that useful here)
p_scatter = scatter(
    medians
    ,1:length(y_labels)
    ,xerror = (lower_errors, upper_errors)
    #, markercolor = colors
    #, markerstrokecolor = colors
    #, color = colors
    ,seriestype = :scatter
    ,yticks = (1:length(y_labels), y_labels)
    ,xlabel = "Median time to detection (days)\n (with 95% CI via 1000 bootstraps)"
    #,xlims = [40,90] #[0,110]
    ,xlims = [10,90] #[0,110]
    ,xticks = 0:10:100
    ,ylabel = ""
    #,title = "Time to detection\n with different site selection"
    ,legend = false
    ,yflip = true
)
plot(p_scatter, layout = (1, 1)
        , size = (900,600)) #default size is 600,400

Plots.savefig( joinpath( output_folder, "scatter_TD_combined_examples.png" ) )


# TD3 results scatter plot
medians_td3 = [ # Scenario / config file 1
            tds_df[1,:sc_median_3TD], tds_df[1,:sc_icu_8_comb_median_3TD], tds_df[1,:sc_icu_29_comb_median_3TD], tds_df[1,:sc_gp_comb_median_3TD], tds_df[1,:sc_icu_8_gp_comb_median_3TD], tds_df[1,:sc_icu_29_gp_comb_median_3TD]
          ,0.0
          # Scenario / config file 13
          ,tds_df[8,:sc_median_3TD], tds_df[8,:sc_icu_8_comb_median_3TD], tds_df[8,:sc_icu_29_comb_median_3TD], tds_df[8,:sc_gp_comb_median_3TD], tds_df[8,:sc_icu_8_gp_comb_median_3TD], tds_df[8,:sc_icu_29_gp_comb_median_3TD]
          ,0.0
          # Scenario / config file 25
          ,tds_df[10,:sc_median_3TD], tds_df[10,:sc_icu_8_comb_median_3TD], tds_df[10,:sc_icu_29_comb_median_3TD], tds_df[10,:sc_gp_comb_median_3TD], tds_df[10,:sc_icu_8_gp_comb_median_3TD], tds_df[10,:sc_icu_29_gp_comb_median_3TD]
          ,0.0
          # Scenario / config file 37
          ,tds_df[12,:sc_median_3TD], tds_df[12,:sc_icu_8_comb_median_3TD], tds_df[12,:sc_icu_29_comb_median_3TD], tds_df[12,:sc_gp_comb_median_3TD], tds_df[12,:sc_icu_8_gp_comb_median_3TD], tds_df[12,:sc_icu_29_gp_comb_median_3TD]
          ,0.0
          # Scenario / config file 49
          ,tds_df[14,:sc_median_3TD], tds_df[14,:sc_icu_8_comb_median_3TD], tds_df[14,:sc_icu_29_comb_median_3TD], tds_df[14,:sc_gp_comb_median_3TD], tds_df[14,:sc_icu_8_gp_comb_median_3TD], tds_df[14,:sc_icu_29_gp_comb_median_3TD]
          ,0.0
          # Scenario / config file 61
          ,tds_df[16,:sc_median_3TD], tds_df[16,:sc_icu_8_comb_median_3TD], tds_df[16,:sc_icu_29_comb_median_3TD], tds_df[16,:sc_gp_comb_median_3TD], tds_df[16,:sc_icu_8_gp_comb_median_3TD], tds_df[16,:sc_icu_29_gp_comb_median_3TD]
          ]
lower_errors_td3 = medians_td3 .- [ # Scenario / config file 1
                            tds_df[1,:sc_median_3TD_lower], tds_df[1,:sc_icu_8_comb_median_3TD_lower], tds_df[1,:sc_icu_29_comb_median_3TD_lower], tds_df[1,:sc_gp_comb_median_3TD_lower], tds_df[1,:sc_icu_8_gp_comb_median_3TD_lower], tds_df[1,:sc_icu_29_gp_comb_median_3TD_lower]
                            ,0.0
                            # Scenario / config file 13
                            ,tds_df[8,:sc_median_3TD_lower], tds_df[8,:sc_icu_8_comb_median_3TD_lower], tds_df[8,:sc_icu_29_comb_median_3TD_lower], tds_df[8,:sc_gp_comb_median_3TD_lower], tds_df[8,:sc_icu_8_gp_comb_median_3TD_lower], tds_df[8,:sc_icu_29_gp_comb_median_3TD_lower]
                            ,0.0
                            # Scenario / config file 25
                            ,tds_df[10,:sc_median_3TD_lower], tds_df[10,:sc_icu_8_comb_median_3TD_lower], tds_df[10,:sc_icu_29_comb_median_3TD_lower], tds_df[10,:sc_gp_comb_median_3TD_lower], tds_df[10,:sc_icu_8_gp_comb_median_3TD_lower], tds_df[10,:sc_icu_29_gp_comb_median_3TD_lower]
                            ,0.0
                            # Scenario / config file 37
                            ,tds_df[12,:sc_median_3TD_lower], tds_df[12,:sc_icu_8_comb_median_3TD_lower], tds_df[12,:sc_icu_29_comb_median_3TD_lower], tds_df[12,:sc_gp_comb_median_3TD_lower], tds_df[12,:sc_icu_8_gp_comb_median_3TD_lower], tds_df[12,:sc_icu_29_gp_comb_median_3TD_lower]
                            ,0.0
                            # Scenario / config file 49
                            ,tds_df[14,:sc_median_3TD_lower], tds_df[14,:sc_icu_8_comb_median_3TD_lower], tds_df[14,:sc_icu_29_comb_median_3TD_lower], tds_df[14,:sc_gp_comb_median_3TD_lower], tds_df[14,:sc_icu_8_gp_comb_median_3TD_lower], tds_df[14,:sc_icu_29_gp_comb_median_3TD_lower]
                            ,0.0
                            # Scenario / config file 61
                            ,tds_df[16,:sc_median_3TD_lower], tds_df[16,:sc_icu_8_comb_median_3TD_lower], tds_df[16,:sc_icu_29_comb_median_3TD_lower], tds_df[16,:sc_gp_comb_median_3TD_lower], tds_df[16,:sc_icu_8_gp_comb_median_3TD_lower], tds_df[16,:sc_icu_29_gp_comb_median_3TD_lower]
                            ]
upper_errors_td3 = [ # Scenario / config file 1
                            tds_df[1,:sc_median_3TD_upper], tds_df[1,:sc_icu_8_comb_median_3TD_upper], tds_df[1,:sc_icu_29_comb_median_3TD_upper], tds_df[1,:sc_gp_comb_median_3TD_upper], tds_df[1,:sc_icu_8_gp_comb_median_3TD_upper], tds_df[1,:sc_icu_29_gp_comb_median_3TD_upper]
                            ,0.0
                            # Scenario / config file 13
                            ,tds_df[8,:sc_median_3TD_upper], tds_df[8,:sc_icu_8_comb_median_3TD_upper], tds_df[8,:sc_icu_29_comb_median_3TD_upper], tds_df[8,:sc_gp_comb_median_3TD_upper], tds_df[8,:sc_icu_8_gp_comb_median_3TD_upper], tds_df[8,:sc_icu_29_gp_comb_median_3TD_upper]
                            ,0.0
                            # Scenario / config file 25
                            ,tds_df[10,:sc_median_3TD_upper], tds_df[10,:sc_icu_8_comb_median_3TD_upper], tds_df[10,:sc_icu_29_comb_median_3TD_upper], tds_df[10,:sc_gp_comb_median_3TD_upper], tds_df[10,:sc_icu_8_gp_comb_median_3TD_upper], tds_df[10,:sc_icu_29_gp_comb_median_3TD_upper]
                            ,0.0
                            # Scenario / config file 37
                            ,tds_df[12,:sc_median_3TD_upper], tds_df[12,:sc_icu_8_comb_median_3TD_upper], tds_df[12,:sc_icu_29_comb_median_3TD_upper], tds_df[12,:sc_gp_comb_median_3TD_upper], tds_df[12,:sc_icu_8_gp_comb_median_3TD_upper], tds_df[12,:sc_icu_29_gp_comb_median_3TD_upper]
                            ,0.0
                            # Scenario / config file 49
                            ,tds_df[14,:sc_median_3TD_upper], tds_df[14,:sc_icu_8_comb_median_3TD_upper], tds_df[14,:sc_icu_29_comb_median_3TD_upper], tds_df[14,:sc_gp_comb_median_3TD_upper], tds_df[14,:sc_icu_8_gp_comb_median_3TD_upper], tds_df[14,:sc_icu_29_gp_comb_median_3TD_upper]
                            ,0.0
                            # Scenario / config file 61
                            ,tds_df[16,:sc_median_3TD_upper], tds_df[16,:sc_icu_8_comb_median_3TD_upper], tds_df[16,:sc_icu_29_comb_median_3TD_upper], tds_df[16,:sc_gp_comb_median_3TD_upper], tds_df[16,:sc_icu_8_gp_comb_median_3TD_upper], tds_df[16,:sc_icu_29_gp_comb_median_3TD_upper]
                            ] .- medians_td3

colors = vcat(fill(:blue, 7), [:green, :lightgreen, :red]
            ,:black
            ,fill(:blue, 2), [:green, :lightgreen, :red]
            ,:black
            ,fill(:blue, 2), [:green, :lightgreen, :red]
            ,:black
            ,fill(:blue, 2), [:green, :lightgreen, :red]
            ,:black
            ,fill(:blue, 2), [:green, :lightgreen, :red]
            ,:black
            ,fill(:blue, 2), [:green, :lightgreen, :red])

# Plot with error bars
# (CI using samples from fitted distribution probably not that useful here)
p_scatter = scatter(
    medians_td3
    ,1:length(y_labels)
    ,xerror = (lower_errors_td3, upper_errors_td3)
    #, markercolor = colors
    #, markerstrokecolor = colors
    #, color = colors
    ,seriestype = :scatter
    ,yticks = (1:length(y_labels), y_labels)
    ,xlabel = "Median time to detection (days)\n (with 95% CI via 1000 bootstraps)"
    #,xlims = [40,90] #[0,110]
    ,xlims = [10,90] #[0,110]
    ,xticks = 0:10:100
    ,ylabel = ""
    #,title = "Time to detection\n with different site selection"
    ,legend = false
    ,yflip = true
)
plot(p_scatter, layout = (1, 1)
        , size = (900,600)) #default size is 600,400

Plots.savefig( joinpath( output_folder, "scatter_3TD_combined_examples.png" ) )


#############################################################################################

### Simulation using timports generated using Poisson sampling of GLEAM output ###

# Simulation run on Imperial College High Performance Computing (HPC) system
# using 'covidlike-1.4.3_w_G_filtered_nrep10.jl' and 'covidlike-1.4.3_w_G_filtered_array100_nrep10.pbs'
# Code copied below for reference (with some commented out code removed for clarity)
# Each simulation run had 10 replicates but a PBS array (covidlike-1.4.3_w_G_filtered_array100_nrep10.pbs) was used to repeat this up to 100 times
# Due to the size of the files generated, only 34 runs were combined, generating a total of 
# 340 simulation replicates (10 sim replicates x 34 runs in array)

"""
#= 
- Simulate multiple replicates using default covid-like parameters 
- Run to test use of GLEAM model output for importation times 
-- importmodel=:GLEAM_outbreak_country_not_specified in simforest()
-- maxtime reduced to 60 because the number of importations is much higher than in other methods
=#
using Pkg
Pkg.add("JLD2")
Pkg.add("DataFrames")
#Pkg.add(url="https://github.com/emvolz/NBPMscape", rev="secondary_care_sampling") # commit f676383ae33166a2550814d7f011676311a7de52 at 16:56 29 Jan 2026
Pkg.add(url="https://github.com/emvolz/NBPMscape", rev="f676383ae33166a2550814d7f011676311a7de52") # 'secondary_care_sampling' branch committed at 16:56 29 Jan 2026
using NBPMscape 
using JLD2 
using DataFrames

initialize_parameters("config/outbreak_params_covid19_like.yaml"); # Real site files are not required for simulating the outbreak. They can be replaced when simulating sampling
# dummy site file names also in core.jl as these are checked when installing NBPMscape

NREPS = 10 #50 #500 #1000 #5000
#sims = map( _-> simforest(NBPMscape.P; initialtime=0.0, maxtime=90.0, maxgenerations=100, max_cases=100000000), 1:NREPS )
sims = map( _-> simforest(NBPMscape.P; initialtime=0.0, maxtime=60.0, maxgenerations=100, max_cases=100000000, importmodel=:GLEAM_outbreak_country_not_specified), 1:NREPS )
#@save "covidlike-1.4.3-sims-nrep1.jld2" sims

### Also save versions filtered for G and cases accessing healthcare, i.e. GP, ED, hospital and ICU (reduces file size and more likely to be able to reload)

# Create vectors to store filtered dataframes: G filtered for ICU and GP rows
#sims_G_icu_filter = [DataFrame() for _ in 1:length(sims)]
#sims_G_gp_filter = [DataFrame() for _ in 1:length(sims)]
# Saving as a single file avoids some overlap of cases, for example where the individual visits the GP before ICU
sims_G_filtered = [DataFrame() for _ in 1:length(sims)]
# Create df to save the max_cases_df output - useful for checking the impact on the number of infections
sims_max_cases_df = [DataFrame() for _ in 1:length(sims)]
# Create df to save the simid and tinf values for each simulation replicate
sims_simid_tinf_df = [DataFrame() for _ in 1:length(sims)]

# Loop through sim replicates
for s in 1:length(sims)	

	try
		fo = sims[s]
		mc = fo.max_cases_df
		
		G = fo.G

		if size(mc,1) > 0
			sims_max_cases_df[s] = mc
			sims_simid_tinf_df[s] = unique(G[:, [:simid, :tinf]])	
		else
			sims_max_cases_df[s] = missing
			sims_simid_tinf_df[s] = missing
		end

		if size(G,1) > 0
			sims_simid_tinf_df[s] = unique(G[:, [:simid, :tinf]])
		else
			sims_simid_tinf_df[s] = missing
		end
		
	catch err
		@warn "Error in iteration" exception = err
		println("Error in sim number: ", s)
		continue
	end

	try

		fo = sims[s]

		# Filter for G, which is df containing infection information
		# and only retain cases (rows) that progressed to healthcareICU (i.e. capable of detection under primary care, secondary care or ICU sampling methodology)
		G_filtered = fo.G[ isfinite.(fo.G.tgp) .| isfinite.(fo.G.ted) .| isfinite.(fo.G.thospital) .| isfinite.(fo.G.ticu), : ]
		# Trim columns to save on storage space
		G_filtered = G_filtered[:,["pid"
						,"tinf", "tgp", "ted", "thospital", "ticu", "tstepdown", "tdischarge", "trecovered", "tdeceased"
						, "severity", "fatal", "iscommuter"
						, "homeregion"
						#, "commuteregion"
						#, "generation"
						#, "F", "G", "H"
						, "infector_age", "infectee_age"
						, "importedinfection"
						, "simid"
						]]
		
		
		if size(G_filtered,1) > 0
			sims_G_filtered[s] = G_filtered
		else
			sims_G_filtered[s] = missing
		end
	catch err
		@warn "Error in iteration" exception = err
		println("Error in sim number: ", s)
		continue
	end
	
end

#Save filtered data to files
@save "covidlike-1.4.3-sims-nrep10_filtered_G.jld2" sims_G_filtered
@save "covidlike-1.4.3-sims-nrep10_max_cases_df.jld2" sims_max_cases_df
@save "covidlike-1.4.3-sims-nrep10_sims_simid_tinf_df.jld2" sims_simid_tinf_df

"""

## Load simulated outbreak data
sims = combine_sim_reps(;  # sim_input_folders in vector format where multiple paths
                          sim_input_folders = [ "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1633695/1633695_sims_G_filtered/" ]
                        , sim_object_name = "sims_G_filtered" #"sims" # "sims_G_icu_filter" # "sims_G_gp_filter" # "sim_G_filtered"
                        , nrep = 10
                        )

# Save to file
@save "covidlike-1.4.3-sims-filtered_G_nrep340_1633695.jld2" sims

# Reload file if necessary
#sims = load("covidlike-1.4.3-sims-filtered_G_nrep340_1633695.jld2", "sims")

# Define path to output folder for analyses files
output_folder = "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1633695_nrep340_analysis/"

## Compute times to detection for secondary care sampling (via HARISS network)
# using different model parameters representing different scenarios
# See config_file_record.xlsx for summary descriptions of scenarios
# Analysis for each set of parameters takes between ~1 and ~3 hours
for i in [1] # i=1
    initialize_parameters("config/HARISS/covid19_like_params_HARISS_$(i).yaml"); #    initialize_parameters("config/covid19_like_params_HARISS_1.yaml")
    sc_tds = secondary_care_td(sims = sims)
    CSV.write( joinpath( output_folder,"sc_tds_$(i).csv" ), sc_tds )
end
# Stared at 20:42 29 Jan 2026

for i in [2,3,4,5,6,7,13,19,25,31,37,43,49,55,61,67] # i=1
    initialize_parameters("config/HARISS/covid19_like_params_HARISS_$(i).yaml"); #    initialize_parameters("config/covid19_like_params_HARISS_1.yaml")
    sc_tds = secondary_care_td(sims = sims)
    CSV.write( joinpath( output_folder,"sc_tds_$(i).csv" ), sc_tds )
end
# Started at 23:51 on 29 Jan. Finished at 00:55 on 30 Jan 2026

# Run analysis for 12 different scenarios (2 seasons x 2 lists of ICU sites x 3 numbers of samples) using ICUs and GPs
for j in [1,13,25,37,49,61,73,74,75,76,77,78]
    initialize_parameters("config/HARISS/covid19_like_params_HARISS_$(j).yaml")
    icu_tds = icu_td(sims = sims)
    gp_tds = gp_td(sims = sims)
    CSV.write( joinpath( output_folder,"icu_tds_$(j).csv" ), icu_tds )
    CSV.write( joinpath( output_folder,"gp_tds_$(j).csv" ), gp_tds )
end

## Read in sc_tds files and compute median TD values
using CSV
using DataFrames

# Directory and base file name
path_name = "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1633695_nrep340_analysis/"
base = "sc_tds"
scenario_numbers = [1,2,3,4,5,6,7,13,19,25,31,37,43,49,55,61,67]
# Prepare a container to keep the DataFrames
sc_tds_df_dict = Dict{Int, DataFrame}()
# Prepare a df to hold the median TD values
sc_tds_df = DataFrame(
                    Scenario_number = zeros(Int, length(scenario_numbers))
                    ,sc_median_TD = zeros(Float64, length(scenario_numbers))
                    ,sc_median_3TD = zeros(Float64, length(scenario_numbers))
                    )

# Loop through scenarios, read in files and compute median TD values
for i in 1:length(scenario_numbers) 
    println(i)
    k = scenario_numbers[i]
    # Construct the file name
    file = joinpath(path_name, "$(base)_$(string(k)).csv")

    # Read CSV into a DataFrame and store in the dictionary with the scenario number (numeric suffix) as key
    sc_tds_df_dict[k] = CSV.read(file, DataFrame)
                        
    # Compute times to detection and populate df
    sc_tds_df[i,"Scenario_number"] = k
    sc_tds_df[i,"sc_median_TD"] = median(sc_tds_df_dict[k][:,"SC_TD"])
    sc_tds_df[i,"sc_median_3TD"] = median(sc_tds_df_dict[k][:,"SC_3TD"])
end

# View
println(sc_tds_df)
println(sc_tds_df_dict[1][1:10,:])

# Save
CSV.write( joinpath( output_folder, "sc_median_tds.csv" ), sc_tds_df )


# Compute median TDs for ICU and GP sampling
scenario_numbers_2 = [1,13,25,37,49,61,73,74,75,76,77,78] # Note that scenarios 73-78 are with different ICU sites
# Prepare a container to keep the DataFrames
icu_tds_df_dict = Dict{Int, DataFrame}()
gp_tds_df_dict = Dict{Int, DataFrame}()
# Prepare a df to hold the median TD values
icu_gp_tds_df = DataFrame(
                    Scenario_number = zeros(Int, length(scenario_numbers_2))
                    ,icu_median_TD = zeros(Float64, length(scenario_numbers_2))
                    ,icu_median_3TD = zeros(Float64, length(scenario_numbers_2))
                    ,gp_median_TD = zeros(Float64, length(scenario_numbers_2))
                    ,gp_median_3TD = zeros(Float64, length(scenario_numbers_2))
                    )

# Loop through scenarios, read in files and compute median TD values
for i in 1:length(scenario_numbers_2) 
    k = scenario_numbers_2[i]
    # Construct the file name, pad with leading zeros if needed
    # e.g., sc_tds_01.csv for k=1
    #file = joinpath(path_name, "$(base)$(lpad(string(k), 2, '0'))$(ext)")
    file_icu = joinpath(path_name, "icu_tds_$(string(k)).csv")
    file_gp  = joinpath(path_name,  "gp_tds_$(string(k)).csv")

    # Read CSV into a DataFrame and store in the dictionary with the scenario number (numeric suffix) as key
    icu_tds_df_dict[k] = CSV.read(file_icu, DataFrame)
    gp_tds_df_dict[k] = CSV.read(file_gp, DataFrame)
                        
    # Compute times to detection and populate df
    icu_gp_tds_df[i,"Scenario_number"] = k
    icu_gp_tds_df[i,"icu_median_TD"] = median(icu_tds_df_dict[k][:,"ICU_TD"])
    icu_gp_tds_df[i,"icu_median_3TD"] = median(icu_tds_df_dict[k][:,"ICU_3TD"])
    icu_gp_tds_df[i,"gp_median_TD"] = median(gp_tds_df_dict[k][:,"GP_TD"])
    icu_gp_tds_df[i,"gp_median_3TD"] = median(gp_tds_df_dict[k][:,"GP_3TD"])
end

# View
println(icu_gp_tds_df)

println(icu_tds_df_dict[1][1:10,:])
println(gp_tds_df_dict[1][1:10,:])

# Save
CSV.write( joinpath( output_folder, "icu_gp_median_tds.csv" ), icu_gp_tds_df )



## Plot results

# Results for scenario 1
# Winter, 300 samples, Equal weighting, no age targeting
plot_histogram( sc_data = sc_tds_df_dict[1][:,:SC_TD]
              , icu_data_1 = icu_tds_df_dict[1][:,:ICU_TD]
              , icu_data_2 = icu_tds_df_dict[73][:,:ICU_TD]
              , gp_data = gp_tds_df_dict[1][:,:GP_TD]
              ,xlim = [0,125]
              ,ylim = [0,300])
# Save plot
Plots.savefig( joinpath( output_folder, "histogram_1.png" ) )

## Comparing results from varying sample size

# Results for scenario 37
# Summer, 300 samples, Equal weighting, no age targeting
plot_histogram( sc_data = sc_tds_df_dict[37][:,:SC_TD]
              , icu_data_1 = icu_tds_df_dict[37][:,:ICU_TD]
              , icu_data_2 = icu_tds_df_dict[76][:,:ICU_TD]
              , gp_data = gp_tds_df_dict[37][:,:GP_TD]
              ,xlim = [0,125]
              ,ylim = [0,300])
Plots.savefig( joinpath( output_folder, "histogram_37.png" ) )

# Results for scenario 13
# Winter, 500 samples, Equal weighting, no age targeting
plot_histogram( sc_data = sc_tds_df_dict[13][:,:SC_TD]
              , icu_data_1 = icu_tds_df_dict[13][:,:ICU_TD]
              , icu_data_2 = icu_tds_df_dict[74][:,:ICU_TD]
              , gp_data = gp_tds_df_dict[13][:,:GP_TD]
              ,xlim = [0,125]
              ,ylim = [0,300])
Plots.savefig( joinpath( output_folder, "histogram_13.png" ) )

# Results for scenario 25
# Winter, 1000 samples, Equal weighting, no age targeting
plot_histogram( sc_data = sc_tds_df_dict[25][:,:SC_TD]
              , icu_data_1 = icu_tds_df_dict[25][:,:ICU_TD]
              , icu_data_2 = icu_tds_df_dict[75][:,:ICU_TD]
              , gp_data = gp_tds_df_dict[25][:,:GP_TD]
              ,xlim = [0,125]
              ,ylim = [0,300])
Plots.savefig( joinpath( output_folder, "histogram_25.png" ) )

## Comparing results from varying sample age target
# Define data
sc_data_1 = sc_tds_df_dict[1][:,:SC_TD]; sc_median_1 = median(sc_data_1) # Winter, 300 samples, Equal weighting, no age targeting
sc_data_2 = sc_tds_df_dict[2][:,:SC_TD]; sc_median_2 = median(sc_data_2) # Winter, 300 samples, Equal weighting, 100% adult : 0% child (if this split cannot be achieved then sample from the other age category are taken)
sc_data_3 = sc_tds_df_dict[3][:,:SC_TD]; sc_median_3 = median(sc_data_3) # Winter, 300 samples, Equal weighting, 75% adult : 25% child (if this split cannot be achieved then sample from the other age category are taken)
sc_data_4 = sc_tds_df_dict[4][:,:SC_TD]; sc_median_4 = median(sc_data_4) # Winter, 300 samples, Equal weighting, 50% adult : 50% child (if this split cannot be achieved then sample from the other age category are taken)
sc_data_5 = sc_tds_df_dict[5][:,:SC_TD]; sc_median_5 = median(sc_data_5) # Winter, 300 samples, Equal weighting, 25% adult : 75% child (if this split cannot be achieved then sample from the other age category are taken)
sc_data_6 = sc_tds_df_dict[6][:,:SC_TD]; sc_median_6 = median(sc_data_6) # Winter, 300 samples, Equal weighting, 0% adult : 25% child (if this split cannot be achieved then sample from the other age category are taken)
# Plot histograms
xlim = [0,125]
ylim = [0,300]
    
    p1 = histogram( sc_data_1, color = :blue, alpha= 0.5, label = "HARISS - no age target",xlimit = xlim,ylimit = ylim ) #, normalize = :pdf)
    p1 = vline!( [sc_median_1], color=:blue, label="median", linewidth = 3) 
    
    p2 = histogram( sc_data_2, color = :blue, alpha= 0.5, label = "HARISS - 100% adult", xlimit = xlim, ylimit = ylim) #, normalize = :pdf)
    p2 = vline!( [sc_median_2], color=:blue, label="median", linewidth = 3) 
    
    p3 = histogram( sc_data_3, color = :blue, alpha= 0.5, label = "HARISS - 75% adult", xlimit = xlim, ylimit = ylim) #, normalize = :pdf)
    p3 = vline!( [sc_median_3], color=:blue, label="median", linewidth = 3) 
    
    p4 = histogram( sc_data_4, color = :blue, alpha= 0.5, label = "HARISS - 50% adult", xlimit = xlim, ylimit = ylim
                  , ylabel = "Number of simulation replicates") #, normalize = :pdf)
    p4 = vline!( [sc_median_4], color=:blue, label="median", linewidth = 3) 
    
    p5 = histogram( sc_data_5, color = :blue, alpha= 0.5, label = "HARISS - 25% adult", xlimit = xlim, ylimit = ylim) #, normalize = :pdf)
    p5 = vline!( [sc_median_5], color=:blue, label="median", linewidth = 3) 
    
    p6 = histogram( sc_data_6, color = :blue, alpha= 0.5, label = "HARISS - 0% adult", xlimit = xlim, ylimit = ylim
                    , xlabel = "Time since first import to first detection (days)" ) #, normalize = :pdf)
    p6 = vline!( [sc_median_6], color=:blue, label="median", linewidth = 3) 
    
    plot(p1, p2, p3, p4, p5, p6, layout = (6, 1)
        , size = (600,800)) #default size is 600,400

Plots.savefig( joinpath( output_folder, "histogram_1to6.png" ) )

## Comparing results from varying method of allocating samples to HARISS network public health laboratories (PHLs)

# Define data
sc_data_1 = sc_tds_df_dict[1][:,:SC_TD]; sc_median_1 = median(sc_data_1) # Winter, 300 samples, Equal weighting, no age targeting
sc_data_7 = sc_tds_df_dict[7][:,:SC_TD]; sc_median_7 = median(sc_data_7) # Winter, 300 samples, Weighted by ED attendance, no age targeting

sc_data_13 = sc_tds_df_dict[13][:,:SC_TD]; sc_median_13 = median(sc_data_13) # Winter, 500 samples, Equal weighting, no age targeting
sc_data_19 = sc_tds_df_dict[19][:,:SC_TD]; sc_median_19 = median(sc_data_19) # Winter, 500 samples, Weighted by ED attendance, no age targeting

sc_data_25 = sc_tds_df_dict[25][:,:SC_TD]; sc_median_25 = median(sc_data_25) # Winter, 1000 samples, Equal weighting, no age targeting
sc_data_31 = sc_tds_df_dict[31][:,:SC_TD]; sc_median_31 = median(sc_data_31) # Winter, 1000 samples, Weighted by ED attendance, no age targeting
# Plot histograms
xlim = [0,125]
ylim = [0,300]
# 300 samples, equal weight    
p1 = histogram( sc_data_1, color = :blue, alpha= 0.5, label = "HARISS - 300 samples, equal",xlimit = xlim,ylimit = ylim ) #, normalize = :pdf)
p1 = vline!( [sc_median_1], color=:blue, label="median", linewidth = 3)
# 300 samples, weighted
p2 = histogram( sc_data_7, color = :orange, alpha= 0.5, label = "HARISS - 300 samples, weighted", xlimit = xlim, ylimit = ylim)
p2 = vline!( [sc_median_7], color=:orange, label="median", linewidth = 3)
# 500 samples, equal weight
p3 = histogram( sc_data_13, color = :blue, alpha= 0.5, label = "HARISS - 500 samples, equal", xlimit = xlim, ylimit = ylim) #, normalize = :pdf)
p3 = vline!( [sc_median_13], color=:blue, label="median", linewidth = 3) 
# 500 samples, weighted
p4 = histogram( sc_data_19, color = :orange, alpha= 0.5, label = "HARISS - 500 samples, weighted", xlimit = xlim, ylimit = ylim
              , ylabel = "Number of simulation replicates") #, normalize = :pdf)
p4 = vline!( [sc_median_19], color=:orange, label="median", linewidth = 3) 
# 1000 samples, equal weight
p5 = histogram( sc_data_25, color = :blue, alpha= 0.5, label = "HARISS - 1000 samples, equal", xlimit = xlim, ylimit = ylim) #, normalize = :pdf)
p5 = vline!( [sc_median_25], color=:blue, label="median", linewidth = 3) 
# 1000 samples, weighted
p6 = histogram( sc_data_31, color = :orange, alpha= 0.5, label = "HARISS - 1000 samples, weighted", xlimit = xlim, ylimit = ylim
                    , xlabel = "Time since first import to first detection (days)" ) #, normalize = :pdf)
p6 = vline!( [sc_median_31], color=:orange, label="median", linewidth = 3) 
    
plot(p1, p2, p3, p4, p5, p6, layout = (6, 1)
        , size = (600,800)) #default size is 600,400

Plots.savefig( joinpath( output_folder, "histogram_equal_v_weighted.png" ) )

## Plot all median TDs together

# TD (1 case)
# Compute median TDs and 95% CI using bootstrapping for median TDs
sc_td_medians_dict = Dict{Int,NamedTuple}()
for (k,v) in sc_tds_df_dict # v=sc_tds_df_dict[1]
    sc_td_medians_dict[k] = NBPMscape.median_ci_bootstrap( vec = v[:,:SC_TD], n_boot = 1000, alpha = 0.05 ) # [ median, lower, upper ] median_ci_bootstrap function is in misc_functions.jl
end # sc_td_medians_dict[1]

icu_td_medians_dict = Dict{Int,NamedTuple}()
for (k,v) in icu_tds_df_dict # v=sc_tds_df_dict[1]
    icu_td_medians_dict[k] = median_ci_bootstrap( vec = v[:,:ICU_TD], n_boot = 1000, alpha = 0.05 ) # [ median, lower, upper ]
end 

gp_td_medians_dict = Dict{Int,NamedTuple}()
for (k,v) in gp_tds_df_dict # v=sc_tds_df_dict[1]
    gp_td_medians_dict[k] = median_ci_bootstrap( vec = v[:,:GP_TD], n_boot = 1000, alpha = 0.05 ) # [ median, lower, upper ]
end

y_labels = ["HARISS: 300 samples, equal weight, winter, no age target"
           ,"HARISS: 300 samples, equal weight, winter, 100% adult   "
           ,"HARISS: 300 samples, equal weight, winter, 75% adult    "
           ,"HARISS: 300 samples, equal weight, winter, 50% adult    "
           ,"HARISS: 300 samples, equal weight, winter, 25% adult    "
           ,"HARISS: 300 samples, equal weight, winter, 0% adult     "
           ,"HARISS: 300 samples, ED weighted,  winter, no age target"
           ,"ICU 8 sites: 300 samples, winter"
           ,"ICU 29 sites: 300 samples, winter"
           ,"RCGP: 300 samples, winter"
           ,""
           ,"HARISS: 500 samples, equal weight, winter, no age target"
           ,"HARISS: 500 samples, ED weighted, winter, no age target"
           ,"ICU 8 sites: 500 samples, winter"
           ,"ICU 29 sites: 500 samples, winter"
           ,"RCGP: 500 samples, winter"
           ,""
           ,"HARISS: 1000 samples, equal weight, winter, no age target"
           ,"HARISS: 1000 samples, ED weighted, winter, no age target"
           ,"ICU 8 sites: 1000 samples, winter"
           ,"ICU 29 sites: 1000 samples, winter"
           ,"RCGP: 1000 samples, winter"
           ,""
           ,"HARISS: 300 samples, equal weight, summer, no age target"
           ,"HARISS: 300 samples, ED weighted, summer, no age target"
           ,"ICU 8 sites: 300 samples, summer"
           ,"ICU 29 sites: 300 samples, summer"
           ,"RCGP: 300 samples, summer"
           ,""
           ,"HARISS: 500 samples, equal weight, summer, no age target"
           ,"HARISS: 500 samples, ED weighted, summer, no age target"
           ,"ICU 8 sites: 500 samples, summer"
           ,"ICU 29 sites: 500 samples, summer"
           ,"RCGP: 500 samples, summer"
           ,""
           ,"HARISS: 1000 samples, equal weight, summer, no age target"
           ,"HARISS: 1000 samples, ED weighted, summer, no age target"
           ,"ICU 8 sites: 1000 samples, summer"
           ,"ICU 29 sites: 1000 samples, summer"
           ,"RCGP: 1000 samples, summer"
           ]
medians = [sc_td_medians_dict[1].median, sc_td_medians_dict[2].median, sc_td_medians_dict[3].median, sc_td_medians_dict[4].median, sc_td_medians_dict[5].median, sc_td_medians_dict[6].median, sc_td_medians_dict[7].median
          ,icu_td_medians_dict[1].median, icu_td_medians_dict[73].median, gp_td_medians_dict[1].median 
          ,0.0
          ,sc_td_medians_dict[13].median, sc_td_medians_dict[19].median
          ,icu_td_medians_dict[13].median, icu_td_medians_dict[74].median, gp_td_medians_dict[13].median 
          ,0.0
          ,sc_td_medians_dict[25].median, sc_td_medians_dict[31].median
          ,icu_td_medians_dict[25].median, icu_td_medians_dict[75].median, gp_td_medians_dict[25].median 
          ,0.0
          ,sc_td_medians_dict[37].median, sc_td_medians_dict[43].median
          ,icu_td_medians_dict[37].median, icu_td_medians_dict[76].median, gp_td_medians_dict[37].median 
          ,0.0
          ,sc_td_medians_dict[49].median, sc_td_medians_dict[55].median
          ,icu_td_medians_dict[49].median, icu_td_medians_dict[77].median, gp_td_medians_dict[49].median 
          ,0.0
          ,sc_td_medians_dict[61].median, sc_td_medians_dict[67].median
          ,icu_td_medians_dict[61].median, icu_td_medians_dict[78].median, gp_td_medians_dict[61].median 
          ]
lower_errors = medians .- [sc_td_medians_dict[1].lower, sc_td_medians_dict[2].lower, sc_td_medians_dict[3].lower, sc_td_medians_dict[4].lower, sc_td_medians_dict[5].lower, sc_td_medians_dict[6].lower, sc_td_medians_dict[7].lower
                          ,icu_td_medians_dict[1].lower, icu_td_medians_dict[73].lower, gp_td_medians_dict[1].lower 
                          ,0.0
                          ,sc_td_medians_dict[13].lower, sc_td_medians_dict[19].lower
                          ,icu_td_medians_dict[13].lower, icu_td_medians_dict[74].lower, gp_td_medians_dict[13].lower 
                          ,0.0
                          ,sc_td_medians_dict[25].lower, sc_td_medians_dict[31].lower
                          ,icu_td_medians_dict[25].lower, icu_td_medians_dict[75].lower, gp_td_medians_dict[25].lower 
                          ,0.0
                          ,sc_td_medians_dict[37].lower, sc_td_medians_dict[43].lower
                          ,icu_td_medians_dict[37].lower, icu_td_medians_dict[76].lower, gp_td_medians_dict[37].lower 
                          ,0.0
                          ,sc_td_medians_dict[49].lower, sc_td_medians_dict[55].lower
                          ,icu_td_medians_dict[49].lower, icu_td_medians_dict[77].lower, gp_td_medians_dict[49].lower 
                          ,0.0
                          ,sc_td_medians_dict[61].lower, sc_td_medians_dict[67].lower
                          ,icu_td_medians_dict[61].lower, icu_td_medians_dict[78].lower, gp_td_medians_dict[61].lower
                          ]
upper_errors = [sc_td_medians_dict[1].upper, sc_td_medians_dict[2].upper, sc_td_medians_dict[3].upper, sc_td_medians_dict[4].upper, sc_td_medians_dict[5].upper, sc_td_medians_dict[6].upper, sc_td_medians_dict[7].upper
               ,icu_td_medians_dict[1].upper, icu_td_medians_dict[73].upper, gp_td_medians_dict[1].upper 
               ,0.0
               ,sc_td_medians_dict[13].upper, sc_td_medians_dict[19].upper
               ,icu_td_medians_dict[13].upper, icu_td_medians_dict[74].upper, gp_td_medians_dict[13].upper
               ,0.0
               ,sc_td_medians_dict[25].upper, sc_td_medians_dict[31].upper
               ,icu_td_medians_dict[25].upper, icu_td_medians_dict[75].upper, gp_td_medians_dict[25].upper
               ,0.0
               ,sc_td_medians_dict[37].upper, sc_td_medians_dict[43].upper
               ,icu_td_medians_dict[37].upper, icu_td_medians_dict[76].upper, gp_td_medians_dict[37].upper
               ,0.0
               ,sc_td_medians_dict[49].upper, sc_td_medians_dict[55].upper
               ,icu_td_medians_dict[49].upper, icu_td_medians_dict[77].upper, gp_td_medians_dict[49].upper
               ,0.0
               ,sc_td_medians_dict[61].upper, sc_td_medians_dict[67].upper
               ,icu_td_medians_dict[61].upper, icu_td_medians_dict[78].upper, gp_td_medians_dict[61].upper 
               ] .- medians

colors = vcat(fill(:blue, 7), [:green, :lightgreen, :red]
            ,:black
            ,fill(:blue, 2), [:green, :lightgreen, :red]
            ,:black
            ,fill(:blue, 2), [:green, :lightgreen, :red]
            ,:black
            ,fill(:blue, 2), [:green, :lightgreen, :red]
            ,:black
            ,fill(:blue, 2), [:green, :lightgreen, :red]
            ,:black
            ,fill(:blue, 2), [:green, :lightgreen, :red])

# Plot with error bars
# (CI using samples from fitted distribution probably not that useful here)
p_scatter = scatter(
    medians
    ,1:length(y_labels)
    ,xerror = (lower_errors, upper_errors)
    #, markercolor = colors
    #, markerstrokecolor = colors
    , color = colors
    ,seriestype = :scatter
    ,yticks = (1:length(y_labels), y_labels)
    ,xlabel = "Median time to detection (days)\n (with 95% CI via 1000 bootstraps)"
    #,xlims = [40,90] #[0,110]
    ,xlims = [0,90] #[0,110]
    ,xticks = 0:10:100
    ,ylabel = ""
    #,title = "Time to detection\n with different site selection"
    ,legend = false
    ,yflip = true
)
plot(p_scatter, layout = (1, 1)
        , size = (900,600)) #default size is 600,400

Plots.savefig( joinpath( output_folder, "scatter_all_TD.png" ) )

# TD3 (3 cases)
# Compute median TDs and 95% CI using bootstrapping for median TDs
sc_td3_medians_dict = Dict{Int,NamedTuple}()
for (k,v) in sc_tds_df_dict # v=sc_tds_df_dict[1]
    sc_td3_medians_dict[k] = median_ci_bootstrap( vec = v[:,:SC_3TD], n_boot = 1000, alpha = 0.05 ) # [ median, lower, upper ] median_ci_bootstrap function is in misc_functions.jl
end # sc_td_medians_dict[1]

icu_td3_medians_dict = Dict{Int,NamedTuple}()
for (k,v) in icu_tds_df_dict # v=sc_tds_df_dict[1]
    icu_td3_medians_dict[k] = median_ci_bootstrap( vec = v[:,:ICU_3TD], n_boot = 1000, alpha = 0.05 ) # [ median, lower, upper ]
end 

gp_td3_medians_dict = Dict{Int,NamedTuple}()
for (k,v) in gp_tds_df_dict # v=sc_tds_df_dict[1]
    gp_td3_medians_dict[k] = median_ci_bootstrap( vec = v[:,:GP_3TD], n_boot = 1000, alpha = 0.05 ) # [ median, lower, upper ]
end

medians_td3 = [sc_td3_medians_dict[1].median, sc_td3_medians_dict[2].median, sc_td3_medians_dict[3].median, sc_td3_medians_dict[4].median, sc_td3_medians_dict[5].median, sc_td3_medians_dict[6].median, sc_td3_medians_dict[7].median
                ,icu_td3_medians_dict[1].median, icu_td3_medians_dict[73].median, gp_td3_medians_dict[1].median 
                ,0.0
                ,sc_td3_medians_dict[13].median, sc_td3_medians_dict[19].median
                ,icu_td3_medians_dict[13].median, icu_td3_medians_dict[74].median, gp_td3_medians_dict[13].median 
                ,0.0
                ,sc_td3_medians_dict[25].median, sc_td3_medians_dict[31].median
                ,icu_td3_medians_dict[25].median, icu_td3_medians_dict[75].median, gp_td3_medians_dict[25].median 
                ,0.0
                ,sc_td3_medians_dict[37].median, sc_td3_medians_dict[43].median
                ,icu_td3_medians_dict[37].median, icu_td3_medians_dict[76].median, gp_td3_medians_dict[37].median 
                ,0.0
                ,sc_td3_medians_dict[49].median, sc_td3_medians_dict[55].median
                ,icu_td3_medians_dict[49].median, icu_td3_medians_dict[77].median, gp_td3_medians_dict[49].median 
                ,0.0
                ,sc_td3_medians_dict[61].median, sc_td3_medians_dict[67].median
                ,icu_td3_medians_dict[61].median, icu_td3_medians_dict[78].median, gp_td3_medians_dict[61].median 
                ]
lower_errors_td3 = medians_td3 .- [sc_td3_medians_dict[1].lower, sc_td3_medians_dict[2].lower, sc_td3_medians_dict[3].lower, sc_td3_medians_dict[4].lower, sc_td3_medians_dict[5].lower, sc_td3_medians_dict[6].lower, sc_td3_medians_dict[7].lower
                                    ,icu_td3_medians_dict[1].lower, icu_td3_medians_dict[73].lower, gp_td3_medians_dict[1].lower 
                                    ,0.0
                                    ,sc_td3_medians_dict[13].lower, sc_td3_medians_dict[19].lower
                                    ,icu_td3_medians_dict[13].lower, icu_td3_medians_dict[74].lower, gp_td3_medians_dict[13].lower 
                                    ,0.0
                                    ,sc_td3_medians_dict[25].lower, sc_td3_medians_dict[31].lower
                                    ,icu_td3_medians_dict[25].lower, icu_td3_medians_dict[75].lower, gp_td3_medians_dict[25].lower 
                                    ,0.0
                                    ,sc_td3_medians_dict[37].lower, sc_td3_medians_dict[43].lower
                                    ,icu_td3_medians_dict[37].lower, icu_td3_medians_dict[76].lower, gp_td3_medians_dict[37].lower 
                                    ,0.0
                                    ,sc_td3_medians_dict[49].lower, sc_td3_medians_dict[55].lower
                                    ,icu_td3_medians_dict[49].lower, icu_td3_medians_dict[77].lower, gp_td3_medians_dict[49].lower 
                                    ,0.0
                                    ,sc_td3_medians_dict[61].lower, sc_td3_medians_dict[67].lower
                                    ,icu_td3_medians_dict[61].lower, icu_td3_medians_dict[78].lower, gp_td3_medians_dict[61].lower
                                    ]
upper_errors_td3 = [sc_td3_medians_dict[1].upper, sc_td3_medians_dict[2].upper, sc_td3_medians_dict[3].upper, sc_td3_medians_dict[4].upper, sc_td3_medians_dict[5].upper, sc_td3_medians_dict[6].upper, sc_td3_medians_dict[7].upper
                    ,icu_td3_medians_dict[1].upper, icu_td3_medians_dict[73].upper, gp_td3_medians_dict[1].upper 
                    ,0.0
                    ,sc_td3_medians_dict[13].upper, sc_td3_medians_dict[19].upper
                    ,icu_td3_medians_dict[13].upper, icu_td3_medians_dict[74].upper, gp_td3_medians_dict[13].upper
                    ,0.0
                    ,sc_td3_medians_dict[25].upper, sc_td3_medians_dict[31].upper
                    ,icu_td3_medians_dict[25].upper, icu_td3_medians_dict[75].upper, gp_td3_medians_dict[25].upper
                    ,0.0
                    ,sc_td3_medians_dict[37].upper, sc_td3_medians_dict[43].upper
                    ,icu_td3_medians_dict[37].upper, icu_td3_medians_dict[76].upper, gp_td3_medians_dict[37].upper
                    ,0.0
                    ,sc_td3_medians_dict[49].upper, sc_td3_medians_dict[55].upper
                    ,icu_td3_medians_dict[49].upper, icu_td3_medians_dict[77].upper, gp_td3_medians_dict[49].upper
                    ,0.0
                    ,sc_td3_medians_dict[61].upper, sc_td3_medians_dict[67].upper
                    ,icu_td3_medians_dict[61].upper, icu_td3_medians_dict[78].upper, gp_td3_medians_dict[61].upper 
                    ] .- medians_td3

# Plot with error bars
# (CI using samples from fitted distribution probably not that useful here)
p_scatter_td3 = scatter(
    medians_td3
    ,1:length(y_labels)
    ,xerror = (lower_errors_td3, upper_errors_td3)
    #, markercolor = colors
    #, markerstrokecolor = colors
    , color = colors
    ,seriestype = :scatter
    ,yticks = (1:length(y_labels), y_labels)
    ,xlabel = "Median time to detect 3 cases (days)\n (with 95% CI via 1000 bootstraps)"
    #,xlims = [40,90] #[0,110]
    ,xlims = [0,90] #[0,110]
    ,xticks = 0:10:100
    ,ylabel = ""
    #,title = "Time to detection\n with different site selection"
    ,legend = false
    ,yflip = true
)
plot(p_scatter_td3, layout = (1, 1)
        , size = (900,600)) #default size is 600,400

Plots.savefig( joinpath( output_folder, "scatter_all_TD3.png" ) )



###################################################################################################################

## Computing and plotting combined times to detection

# Read in files 

using CSV
using DataFrames

# Directory and base file name
path_name = "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1633695_nrep340_analysis/"

# Match scenario numbers (from config file suffixes) between sampling setting
scenario_numbers_sc =     [1,2,3,4,5,6,7,13,19,25,31,37,43,49,55,61,67]
scenario_numbers_icu_8 =  [1,1,1,1,1,1,1,13,13,25,25,37,37,49,49,61,61]
scenario_numbers_icu_29 = [73,73,73,73,73,73,73,74,74,75,75,76,76,77,77,78,78]
scenario_numbers_gp =     [1,1,1,1,1,1,1,13,13,25,25,37,37,49,49,61,61]

# Prepare containers to store the individual DataFrames
sc_tds_df_dict     = Dict{Int, DataFrame}()
icu_8_tds_df_dict  = Dict{Int, DataFrame}()
icu_29_tds_df_dict = Dict{Int, DataFrame}()
gp_tds_df_dict     = Dict{Int, DataFrame}()

# Prepare a df to hold the median TD values and upper and lower 95% CI from bootstrapping
default_int_col   = zeros(Int,     length(scenario_numbers_sc))
default_float_col = zeros(Float64, length(scenario_numbers_sc))
tds_df = DataFrame(
                    Scenario_number = default_int_col
                    ,Scenario_number_icu_29 = default_int_col # ICU sampling at 29 sites uses different config files because it requires a different parameter value for the ICU sites
                    ## Median
                    # TD
                    ,sc_median_TD = default_float_col
                    ,icu_8_median_TD = default_float_col
                    ,icu_29_median_TD = default_float_col
                    ,gp_median_TD = default_float_col
                    # 3TD
                    ,sc_median_3TD = default_float_col
                    ,icu_8_median_3TD = default_float_col
                    ,icu_29_median_3TD = default_float_col
                    ,gp_median_3TD = default_float_col
                    # Combined TD
                    # SC + ICU 8
                    ,sc_icu_8_comb_median_TD = default_float_col
                    # SC + ICU 29
                    ,sc_icu_29_comb_median_TD = default_float_col
                    # SC + GP
                    ,sc_gp_comb_median_TD = default_float_col
                    # SC + ICU 8 + GP
                    ,sc_icu_8_gp_comb_median_TD = default_float_col
                    # SC + ICU 29 + GP
                    ,sc_icu_29_gp_comb_median_TD = default_float_col
                    # Combined 3TD
                    # SC + ICU 8
                    ,sc_icu_8_comb_median_3TD = default_float_col
                    # SC + ICU 29
                    ,sc_icu_29_comb_median_3TD = default_float_col
                    # SC + GP
                    ,sc_gp_comb_median_3TD = default_float_col
                    # SC + ICU 8 + GP
                    ,sc_icu_8_gp_comb_median_3TD = default_float_col
                    # SC + ICU 29 + GP
                    ,sc_icu_29_gp_comb_median_3TD = default_float_col
                    ## Lower CI
                    # TD
                    ,sc_median_TD_lower = default_float_col
                    ,icu_8_median_TD_lower = default_float_col
                    ,icu_29_median_TD_lower = default_float_col
                    ,gp_median_TD_lower = default_float_col
                    # 3TD
                    ,sc_median_3TD_lower = default_float_col
                    ,icu_8_median_3TD_lower = default_float_col
                    ,icu_29_median_3TD_lower = default_float_col
                    ,gp_median_3TD_lower = default_float_col
                    # Combined TD
                    # SC + ICU 8
                    ,sc_icu_8_comb_median_TD_lower = default_float_col
                    # SC + ICU 29
                    ,sc_icu_29_comb_median_TD_lower = default_float_col
                    # SC + GP
                    ,sc_gp_comb_median_TD_lower = default_float_col
                    # SC + ICU 8 + GP
                    ,sc_icu_8_gp_comb_median_TD_lower = default_float_col
                    # SC + ICU 29 + GP
                    ,sc_icu_29_gp_comb_median_TD_lower = default_float_col
                    # Combined 3TD
                    # SC + ICU 8
                    ,sc_icu_8_comb_median_3TD_lower = default_float_col
                    # SC + ICU 29
                    ,sc_icu_29_comb_median_3TD_lower = default_float_col
                    # SC + GP
                    ,sc_gp_comb_median_3TD_lower = default_float_col
                    # SC + ICU 8 + GP
                    ,sc_icu_8_gp_comb_median_3TD_lower = default_float_col
                    # SC + ICU 29 + GP
                    ,sc_icu_29_gp_comb_median_3TD_lower = default_float_col
                    ## Upper CI
                    # TD
                    ,sc_median_TD_upper = default_float_col
                    ,icu_8_median_TD_upper = default_float_col
                    ,icu_29_median_TD_upper = default_float_col
                    ,gp_median_TD_upper = default_float_col
                    # 3TD
                    ,sc_median_3TD_upper = default_float_col
                    ,icu_8_median_3TD_upper = default_float_col
                    ,icu_29_median_3TD_upper = default_float_col
                    ,gp_median_3TD_upper = default_float_col
                    # Combined TD
                    # SC + ICU 8
                    ,sc_icu_8_comb_median_TD_upper = default_float_col
                    # SC + ICU 29
                    ,sc_icu_29_comb_median_TD_upper = default_float_col
                    # SC + GP
                    ,sc_gp_comb_median_TD_upper = default_float_col
                    # SC + ICU 8 + GP
                    ,sc_icu_8_gp_comb_median_TD_upper = default_float_col
                    # SC + ICU 29 + GP
                    ,sc_icu_29_gp_comb_median_TD_upper = default_float_col
                    # Combined 3TD
                    # SC + ICU 8
                    ,sc_icu_8_comb_median_3TD_upper = default_float_col
                    # SC + ICU 29
                    ,sc_icu_29_comb_median_3TD_upper = default_float_col
                    # SC + GP
                    ,sc_gp_comb_median_3TD_upper = default_float_col
                    # SC + ICU 8 + GP
                    ,sc_icu_8_gp_comb_median_3TD_upper = default_float_col
                    # SC + ICU 29 + GP
                    ,sc_icu_29_gp_comb_median_3TD_upper = default_float_col
                    )


# Loop through scenarios, read in files and compute median TD values
# Secondary care
for i in 1:length(scenario_numbers_sc) 
    j = scenario_numbers_sc[i]
    # Construct the file name, pad with leading zeros if needed
    # e.g., sc_tds_01.csv for k=1
    #file = joinpath(path_name, "$(base)$(lpad(string(k), 2, '0'))$(ext)")
    file = joinpath(path_name, "sc_tds_$(string(j)).csv")

    # Read CSV into a DataFrame and store in the dictionary with the scenario number (numeric suffix) as key
    sc_tds_df_dict[j] = CSV.read(file, DataFrame)
                        
    # Compute times to detection and populate df
    tds_df[i,"Scenario_number"] = j
    tds_df[i,"sc_median_TD"] = median(sc_tds_df_dict[j][:,"SC_TD"])
    tds_df[i,"sc_median_3TD"] = median(sc_tds_df_dict[j][:,"SC_3TD"])
end

# ICU 8 sites
for i in 1:length(scenario_numbers_icu_8) 
    k = scenario_numbers_icu_8[i]
    # Construct the file name, pad with leading zeros if needed
    # e.g., sc_tds_01.csv for k=1
    #file = joinpath(path_name, "$(base)$(lpad(string(k), 2, '0'))$(ext)")
    file = joinpath(path_name, "icu_tds_$(string(k)).csv")

    # Read CSV into a DataFrame and store in the dictionary with the scenario number (numeric suffix) as key
    icu_8_tds_df_dict[k] = CSV.read(file, DataFrame)
                        
    # Compute times to detection and populate df
    tds_df[i,"icu_8_median_TD"] = median(icu_8_tds_df_dict[k][:,"ICU_TD"])
    tds_df[i,"icu_8_median_3TD"] = median(icu_8_tds_df_dict[k][:,"ICU_3TD"])
end

# ICU 29 sites
for i in 1:length(scenario_numbers_icu_29) 
    m = scenario_numbers_icu_29[i]
    # Construct the file name, pad with leading zeros if needed
    # e.g., sc_tds_01.csv for k=1
    #file = joinpath(path_name, "$(base)$(lpad(string(k), 2, '0'))$(ext)")
    file = joinpath(path_name, "icu_tds_$(string(m)).csv")

    # Read CSV into a DataFrame and store in the dictionary with the scenario number (numeric suffix) as key
    icu_29_tds_df_dict[m] = CSV.read(file, DataFrame)
                        
    # Compute times to detection and populate df
    tds_df[i,"Scenario_number_icu_29"] = m
    tds_df[i,"icu_29_median_TD"] = median(icu_29_tds_df_dict[m][:,"ICU_TD"])
    tds_df[i,"icu_29_median_3TD"] = median(icu_29_tds_df_dict[m][:,"ICU_3TD"])
end

# GP
for i in 1:length(scenario_numbers_gp) 
    n = scenario_numbers_gp[i]
    # Construct the file name, pad with leading zeros if needed
    # e.g., sc_tds_01.csv for k=1
    #file = joinpath(path_name, "$(base)$(lpad(string(k), 2, '0'))$(ext)")
    file = joinpath(path_name, "gp_tds_$(string(n)).csv")

    # Read CSV into a DataFrame and store in the dictionary with the scenario number (numeric suffix) as key
    gp_tds_df_dict[n] = CSV.read(file, DataFrame)
                        
    # Compute times to detection and populate df
    tds_df[i,"gp_median_TD"]  = median(gp_tds_df_dict[n][:,"GP_TD"])
    tds_df[i,"gp_median_3TD"] = median(gp_tds_df_dict[n][:,"GP_3TD"])
end

#println(tds_df)

#println(sc_tds_df_dict[1][1:10,:])

# Compute median TDs combining different sampling settings
for i in 1:length(scenario_numbers_sc) # i =1
    j = scenario_numbers_sc[i]
    k = scenario_numbers_icu_8[i]
    m = scenario_numbers_icu_29[i]
    n = scenario_numbers_gp[i]

    # Populate scenario numbers
    tds_df[i,"Scenario_number"] = j # Assumes that j == k == n
    tds_df[i,"Scenario_number_icu_29"] = m

    # Times to detection (for each simulation replicate) for the different sample settings and for the specific parameter set
    # TD
    sc_tds     =     sc_tds_df_dict[ j ][:,"SC_TD" ]
    icu_8_tds  =  icu_8_tds_df_dict[ k ][:,"ICU_TD"]
    icu_29_tds = icu_29_tds_df_dict[ m ][:,"ICU_TD"]
    gp_tds     =     gp_tds_df_dict[ n ][:,"GP_TD" ]
    #println(sample_setting_1[1:10],sample_setting_2[1:10],combined_min_TD[1:10])
    # 3TD
    sc_tds_3     =     sc_tds_df_dict[j][:,"SC_3TD" ]
    icu_8_tds_3  =  icu_8_tds_df_dict[k][:,"ICU_3TD"]
    icu_29_tds_3 = icu_29_tds_df_dict[m][:,"ICU_3TD"]
    gp_tds_3     =     gp_tds_df_dict[n][:,"GP_3TD" ]

    # Compute median times to detection and 95% CI using bootstrap
    # SC
    (tds_df[i,"sc_median_TD"], tds_df[i,"sc_median_TD_lower"], tds_df[i,"sc_median_TD_upper"]) = NBPMscape.median_ci_bootstrap( vec = sc_tds, n_boot = 1000, alpha = 0.05)
    (tds_df[i,"sc_median_3TD"], tds_df[i,"sc_median_3TD_lower"], tds_df[i,"sc_median_3TD_upper"]) = NBPMscape.median_ci_bootstrap( vec = sc_tds_3, n_boot = 1000, alpha = 0.05)
    # ICU 8 sites
    (tds_df[i,"icu_8_median_TD"], tds_df[i,"icu_8_median_TD_lower"], tds_df[i,"icu_8_median_TD_upper"]) = NBPMscape.median_ci_bootstrap( vec = icu_8_tds, n_boot = 1000, alpha = 0.05)
    (tds_df[i,"icu_8_median_3TD"], tds_df[i,"icu_8_median_3TD_lower"], tds_df[i,"icu_8_median_3TD_upper"]) = NBPMscape.median_ci_bootstrap( vec = icu_8_tds_3, n_boot = 1000, alpha = 0.05)
    # ICU 29 sites
    (tds_df[i,"icu_29_median_TD"], tds_df[i,"icu_29_median_TD_lower"], tds_df[i,"icu_29_median_TD_upper"]) = NBPMscape.median_ci_bootstrap( vec = icu_29_tds, n_boot = 1000, alpha = 0.05)
    (tds_df[i,"icu_29_median_3TD"], tds_df[i,"icu_29_median_3TD_lower"], tds_df[i,"icu_29_median_3TD_upper"]) = NBPMscape.median_ci_bootstrap( vec = icu_29_tds_3, n_boot = 1000, alpha = 0.05)
    # GP
    (tds_df[i,"gp_median_TD"], tds_df[i,"gp_median_TD_lower"], tds_df[i,"gp_median_TD_upper"]) = NBPMscape.median_ci_bootstrap( vec = gp_tds, n_boot = 1000, alpha = 0.05)
    (tds_df[i,"gp_median_3TD"], tds_df[i,"gp_median_3TD_lower"], tds_df[i,"gp_median_3TD_upper"]) = NBPMscape.median_ci_bootstrap( vec = gp_tds_3, n_boot = 1000, alpha = 0.05)
    

    # Compute combined times to detection and populate df
    # SC and ICU 8 sites
    sc_icu_8_comb = min.(sc_tds,   icu_8_tds  )
    sc_icu_8_comb_3 = min.(sc_tds_3,   icu_8_tds_3  )

    (tds_df[i,"sc_icu_8_comb_median_TD"], tds_df[i,"sc_icu_8_comb_median_TD_lower"], tds_df[i,"sc_icu_8_comb_median_TD_upper"]) = NBPMscape.median_ci_bootstrap( vec = sc_icu_8_comb, n_boot = 1000, alpha = 0.05)
    (tds_df[i,"sc_icu_8_comb_median_3TD"], tds_df[i,"sc_icu_8_comb_median_3TD_lower"], tds_df[i,"sc_icu_8_comb_median_3TD_upper"]) = NBPMscape.median_ci_bootstrap( vec = sc_icu_8_comb_3, n_boot = 1000, alpha = 0.05)

    # SC and ICU 29 sites
    sc_icu_29_comb = min.(sc_tds,   icu_29_tds  )
    sc_icu_29_comb_3 = min.(sc_tds_3,   icu_29_tds_3  )

    (tds_df[i,"sc_icu_29_comb_median_TD"], tds_df[i,"sc_icu_29_comb_median_TD_lower"], tds_df[i,"sc_icu_29_comb_median_TD_upper"]) = NBPMscape.median_ci_bootstrap( vec = sc_icu_29_comb, n_boot = 1000, alpha = 0.05)
    (tds_df[i,"sc_icu_29_comb_median_3TD"], tds_df[i,"sc_icu_29_comb_median_3TD_lower"], tds_df[i,"sc_icu_29_comb_median_3TD_upper"]) = NBPMscape.median_ci_bootstrap( vec = sc_icu_29_comb_3, n_boot = 1000, alpha = 0.05)

    # SC and GP
    sc_gp_comb = min.(sc_tds,   gp_tds  )
    sc_gp_comb_3 = min.(sc_tds_3,   gp_tds_3  )

    (tds_df[i,"sc_gp_comb_median_TD"], tds_df[i,"sc_gp_comb_median_TD_lower"], tds_df[i,"sc_gp_comb_median_TD_upper"]) = NBPMscape.median_ci_bootstrap( vec = sc_gp_comb, n_boot = 1000, alpha = 0.05)
    (tds_df[i,"sc_gp_comb_median_3TD"], tds_df[i,"sc_gp_comb_median_3TD_lower"], tds_df[i,"sc_gp_comb_median_3TD_upper"]) = NBPMscape.median_ci_bootstrap( vec = sc_gp_comb_3, n_boot = 1000, alpha = 0.05)

    # SC and ICU 8 sites and GP
    sc_icu_8_gp_comb = min.(sc_tds,   icu_8_tds, gp_tds )
    sc_icu_8_gp_comb_3 = min.(sc_tds_3, icu_8_tds_3,   gp_tds_3  )

    (tds_df[i,"sc_icu_8_gp_comb_median_TD"], tds_df[i,"sc_icu_8_gp_comb_median_TD_lower"], tds_df[i,"sc_icu_8_gp_comb_median_TD_upper"]) = NBPMscape.median_ci_bootstrap( vec = sc_icu_8_gp_comb, n_boot = 1000, alpha = 0.05)
    (tds_df[i,"sc_icu_8_gp_comb_median_3TD"], tds_df[i,"sc_icu_8_gp_comb_median_3TD_lower"], tds_df[i,"sc_icu_8_gp_comb_median_3TD_upper"]) = NBPMscape.median_ci_bootstrap( vec = sc_icu_8_gp_comb_3, n_boot = 1000, alpha = 0.05)

    # SC and ICU 29 sites and GP
    sc_icu_29_gp_comb = min.(sc_tds,   icu_29_tds, gp_tds )
    sc_icu_29_gp_comb_3 = min.(sc_tds_3, icu_29_tds_3,   gp_tds_3  )

    (tds_df[i,"sc_icu_29_gp_comb_median_TD"], tds_df[i,"sc_icu_29_gp_comb_median_TD_lower"], tds_df[i,"sc_icu_29_gp_comb_median_TD_upper"]) = NBPMscape.median_ci_bootstrap( vec = sc_icu_29_gp_comb, n_boot = 1000, alpha = 0.05)
    (tds_df[i,"sc_icu_29_gp_comb_median_3TD"], tds_df[i,"sc_icu_29_gp_comb_median_3TD_lower"], tds_df[i,"sc_icu_29_gp_comb_median_3TD_upper"]) = NBPMscape.median_ci_bootstrap( vec = sc_icu_29_gp_comb_3, n_boot = 1000, alpha = 0.05)

end

#CSV.write( joinpath( output_folder, "tds_df_inc_comb.csv" ), tds_df )
CSV.write( joinpath( output_folder, "tds_df_inc_comb_inc_CI.csv") , tds_df ) 

## Plot all median TDs together

# TD (1 case)
# Create base labels and repeat 6 times (representing the 300, 500, 1000 samples in winter and summer), then drop the last "" (to match length)
y_label_base = [
                "HARISS"
                ,"HARISS & ICU 8 sites"
                ,"HARISS & ICU 29 sites"
                ,"HARISS & RCGP"
                ,"HARISS & ICU 8 sites & RCGP"
                ,"HARISS & ICU 29 sites & RCGP"
                ,""
            ]
y_labels = vcat(repeat(y_label_base, 6)...)[1:end-1]

medians = [ # Scenario / config file 1
            tds_df[1,:sc_median_TD], tds_df[1,:sc_icu_8_comb_median_TD], tds_df[1,:sc_icu_29_comb_median_TD], tds_df[1,:sc_gp_comb_median_TD], tds_df[1,:sc_icu_8_gp_comb_median_TD], tds_df[1,:sc_icu_29_gp_comb_median_TD]
          ,0.0
          # Scenario / config file 13
          ,tds_df[8,:sc_median_TD], tds_df[8,:sc_icu_8_comb_median_TD], tds_df[8,:sc_icu_29_comb_median_TD], tds_df[8,:sc_gp_comb_median_TD], tds_df[8,:sc_icu_8_gp_comb_median_TD], tds_df[8,:sc_icu_29_gp_comb_median_TD]
          ,0.0
          # Scenario / config file 25
          ,tds_df[10,:sc_median_TD], tds_df[10,:sc_icu_8_comb_median_TD], tds_df[10,:sc_icu_29_comb_median_TD], tds_df[10,:sc_gp_comb_median_TD], tds_df[10,:sc_icu_8_gp_comb_median_TD], tds_df[10,:sc_icu_29_gp_comb_median_TD]
          ,0.0
          # Scenario / config file 37
          ,tds_df[12,:sc_median_TD], tds_df[12,:sc_icu_8_comb_median_TD], tds_df[12,:sc_icu_29_comb_median_TD], tds_df[12,:sc_gp_comb_median_TD], tds_df[12,:sc_icu_8_gp_comb_median_TD], tds_df[12,:sc_icu_29_gp_comb_median_TD]
          ,0.0
          # Scenario / config file 49
          ,tds_df[14,:sc_median_TD], tds_df[14,:sc_icu_8_comb_median_TD], tds_df[14,:sc_icu_29_comb_median_TD], tds_df[14,:sc_gp_comb_median_TD], tds_df[14,:sc_icu_8_gp_comb_median_TD], tds_df[14,:sc_icu_29_gp_comb_median_TD]
          ,0.0
          # Scenario / config file 61
          ,tds_df[16,:sc_median_TD], tds_df[16,:sc_icu_8_comb_median_TD], tds_df[16,:sc_icu_29_comb_median_TD], tds_df[16,:sc_gp_comb_median_TD], tds_df[16,:sc_icu_8_gp_comb_median_TD], tds_df[16,:sc_icu_29_gp_comb_median_TD]
          ]
lower_errors = medians .- [ # Scenario / config file 1
                            tds_df[1,:sc_median_TD_lower], tds_df[1,:sc_icu_8_comb_median_TD_lower], tds_df[1,:sc_icu_29_comb_median_TD_lower], tds_df[1,:sc_gp_comb_median_TD_lower], tds_df[1,:sc_icu_8_gp_comb_median_TD_lower], tds_df[1,:sc_icu_29_gp_comb_median_TD_lower]
                            ,0.0
                            # Scenario / config file 13
                            ,tds_df[8,:sc_median_TD_lower], tds_df[8,:sc_icu_8_comb_median_TD_lower], tds_df[8,:sc_icu_29_comb_median_TD_lower], tds_df[8,:sc_gp_comb_median_TD_lower], tds_df[8,:sc_icu_8_gp_comb_median_TD_lower], tds_df[8,:sc_icu_29_gp_comb_median_TD_lower]
                            ,0.0
                            # Scenario / config file 25
                            ,tds_df[10,:sc_median_TD_lower], tds_df[10,:sc_icu_8_comb_median_TD_lower], tds_df[10,:sc_icu_29_comb_median_TD_lower], tds_df[10,:sc_gp_comb_median_TD_lower], tds_df[10,:sc_icu_8_gp_comb_median_TD_lower], tds_df[10,:sc_icu_29_gp_comb_median_TD_lower]
                            ,0.0
                            # Scenario / config file 37
                            ,tds_df[12,:sc_median_TD_lower], tds_df[12,:sc_icu_8_comb_median_TD_lower], tds_df[12,:sc_icu_29_comb_median_TD_lower], tds_df[12,:sc_gp_comb_median_TD_lower], tds_df[12,:sc_icu_8_gp_comb_median_TD_lower], tds_df[12,:sc_icu_29_gp_comb_median_TD_lower]
                            ,0.0
                            # Scenario / config file 49
                            ,tds_df[14,:sc_median_TD_lower], tds_df[14,:sc_icu_8_comb_median_TD_lower], tds_df[14,:sc_icu_29_comb_median_TD_lower], tds_df[14,:sc_gp_comb_median_TD_lower], tds_df[14,:sc_icu_8_gp_comb_median_TD_lower], tds_df[14,:sc_icu_29_gp_comb_median_TD_lower]
                            ,0.0
                            # Scenario / config file 61
                            ,tds_df[16,:sc_median_TD_lower], tds_df[16,:sc_icu_8_comb_median_TD_lower], tds_df[16,:sc_icu_29_comb_median_TD_lower], tds_df[16,:sc_gp_comb_median_TD_lower], tds_df[16,:sc_icu_8_gp_comb_median_TD_lower], tds_df[16,:sc_icu_29_gp_comb_median_TD_lower]
                            ]
upper_errors = [ # Scenario / config file 1
                            tds_df[1,:sc_median_TD_upper], tds_df[1,:sc_icu_8_comb_median_TD_upper], tds_df[1,:sc_icu_29_comb_median_TD_upper], tds_df[1,:sc_gp_comb_median_TD_upper], tds_df[1,:sc_icu_8_gp_comb_median_TD_upper], tds_df[1,:sc_icu_29_gp_comb_median_TD_upper]
                            ,0.0
                            # Scenario / config file 13
                            ,tds_df[8,:sc_median_TD_upper], tds_df[8,:sc_icu_8_comb_median_TD_upper], tds_df[8,:sc_icu_29_comb_median_TD_upper], tds_df[8,:sc_gp_comb_median_TD_upper], tds_df[8,:sc_icu_8_gp_comb_median_TD_upper], tds_df[8,:sc_icu_29_gp_comb_median_TD_upper]
                            ,0.0
                            # Scenario / config file 25
                            ,tds_df[10,:sc_median_TD_upper], tds_df[10,:sc_icu_8_comb_median_TD_upper], tds_df[10,:sc_icu_29_comb_median_TD_upper], tds_df[10,:sc_gp_comb_median_TD_upper], tds_df[10,:sc_icu_8_gp_comb_median_TD_upper], tds_df[10,:sc_icu_29_gp_comb_median_TD_upper]
                            ,0.0
                            # Scenario / config file 37
                            ,tds_df[12,:sc_median_TD_upper], tds_df[12,:sc_icu_8_comb_median_TD_upper], tds_df[12,:sc_icu_29_comb_median_TD_upper], tds_df[12,:sc_gp_comb_median_TD_upper], tds_df[12,:sc_icu_8_gp_comb_median_TD_upper], tds_df[12,:sc_icu_29_gp_comb_median_TD_upper]
                            ,0.0
                            # Scenario / config file 49
                            ,tds_df[14,:sc_median_TD_upper], tds_df[14,:sc_icu_8_comb_median_TD_upper], tds_df[14,:sc_icu_29_comb_median_TD_upper], tds_df[14,:sc_gp_comb_median_TD_upper], tds_df[14,:sc_icu_8_gp_comb_median_TD_upper], tds_df[14,:sc_icu_29_gp_comb_median_TD_upper]
                            ,0.0
                            # Scenario / config file 61
                            ,tds_df[16,:sc_median_TD_upper], tds_df[16,:sc_icu_8_comb_median_TD_upper], tds_df[16,:sc_icu_29_comb_median_TD_upper], tds_df[16,:sc_gp_comb_median_TD_upper], tds_df[16,:sc_icu_8_gp_comb_median_TD_upper], tds_df[16,:sc_icu_29_gp_comb_median_TD_upper]
                            ] .- medians

colors = vcat(fill(:blue, 7), [:green, :lightgreen, :red]
            ,:black
            ,fill(:blue, 2), [:green, :lightgreen, :red]
            ,:black
            ,fill(:blue, 2), [:green, :lightgreen, :red]
            ,:black
            ,fill(:blue, 2), [:green, :lightgreen, :red]
            ,:black
            ,fill(:blue, 2), [:green, :lightgreen, :red]
            ,:black
            ,fill(:blue, 2), [:green, :lightgreen, :red])

# Plot with error bars
# (CI using samples from fitted distribution probably not that useful here)
p_scatter = scatter(
    medians
    ,1:length(y_labels)
    ,xerror = (lower_errors, upper_errors)
    #, markercolor = colors
    #, markerstrokecolor = colors
    #, color = colors
    ,seriestype = :scatter
    ,yticks = (1:length(y_labels), y_labels)
    ,xlabel = "Median time to detection (days)\n (with 95% CI via 1000 bootstraps)"
    #,xlims = [40,90] #[0,110]
    ,xlims = [10,60] #[0,110]
    ,xticks = 0:10:100
    ,ylabel = ""
    #,title = "Time to detection\n with different site selection"
    ,legend = false
    ,yflip = true
)
plot(p_scatter, layout = (1, 1)
        , size = (900,600)) #default size is 600,400

Plots.savefig( joinpath( output_folder, "scatter_TD_combined_examples.png" ) )


# TD3 results scatter plot
medians_td3 = [ # Scenario / config file 1
            tds_df[1,:sc_median_3TD], tds_df[1,:sc_icu_8_comb_median_3TD], tds_df[1,:sc_icu_29_comb_median_3TD], tds_df[1,:sc_gp_comb_median_3TD], tds_df[1,:sc_icu_8_gp_comb_median_3TD], tds_df[1,:sc_icu_29_gp_comb_median_3TD]
          ,0.0
          # Scenario / config file 13
          ,tds_df[8,:sc_median_3TD], tds_df[8,:sc_icu_8_comb_median_3TD], tds_df[8,:sc_icu_29_comb_median_3TD], tds_df[8,:sc_gp_comb_median_3TD], tds_df[8,:sc_icu_8_gp_comb_median_3TD], tds_df[8,:sc_icu_29_gp_comb_median_3TD]
          ,0.0
          # Scenario / config file 25
          ,tds_df[10,:sc_median_3TD], tds_df[10,:sc_icu_8_comb_median_3TD], tds_df[10,:sc_icu_29_comb_median_3TD], tds_df[10,:sc_gp_comb_median_3TD], tds_df[10,:sc_icu_8_gp_comb_median_3TD], tds_df[10,:sc_icu_29_gp_comb_median_3TD]
          ,0.0
          # Scenario / config file 37
          ,tds_df[12,:sc_median_3TD], tds_df[12,:sc_icu_8_comb_median_3TD], tds_df[12,:sc_icu_29_comb_median_3TD], tds_df[12,:sc_gp_comb_median_3TD], tds_df[12,:sc_icu_8_gp_comb_median_3TD], tds_df[12,:sc_icu_29_gp_comb_median_3TD]
          ,0.0
          # Scenario / config file 49
          ,tds_df[14,:sc_median_3TD], tds_df[14,:sc_icu_8_comb_median_3TD], tds_df[14,:sc_icu_29_comb_median_3TD], tds_df[14,:sc_gp_comb_median_3TD], tds_df[14,:sc_icu_8_gp_comb_median_3TD], tds_df[14,:sc_icu_29_gp_comb_median_3TD]
          ,0.0
          # Scenario / config file 61
          ,tds_df[16,:sc_median_3TD], tds_df[16,:sc_icu_8_comb_median_3TD], tds_df[16,:sc_icu_29_comb_median_3TD], tds_df[16,:sc_gp_comb_median_3TD], tds_df[16,:sc_icu_8_gp_comb_median_3TD], tds_df[16,:sc_icu_29_gp_comb_median_3TD]
          ]
lower_errors_td3 = medians_td3 .- [ # Scenario / config file 1
                            tds_df[1,:sc_median_3TD_lower], tds_df[1,:sc_icu_8_comb_median_3TD_lower], tds_df[1,:sc_icu_29_comb_median_3TD_lower], tds_df[1,:sc_gp_comb_median_3TD_lower], tds_df[1,:sc_icu_8_gp_comb_median_3TD_lower], tds_df[1,:sc_icu_29_gp_comb_median_3TD_lower]
                            ,0.0
                            # Scenario / config file 13
                            ,tds_df[8,:sc_median_3TD_lower], tds_df[8,:sc_icu_8_comb_median_3TD_lower], tds_df[8,:sc_icu_29_comb_median_3TD_lower], tds_df[8,:sc_gp_comb_median_3TD_lower], tds_df[8,:sc_icu_8_gp_comb_median_3TD_lower], tds_df[8,:sc_icu_29_gp_comb_median_3TD_lower]
                            ,0.0
                            # Scenario / config file 25
                            ,tds_df[10,:sc_median_3TD_lower], tds_df[10,:sc_icu_8_comb_median_3TD_lower], tds_df[10,:sc_icu_29_comb_median_3TD_lower], tds_df[10,:sc_gp_comb_median_3TD_lower], tds_df[10,:sc_icu_8_gp_comb_median_3TD_lower], tds_df[10,:sc_icu_29_gp_comb_median_3TD_lower]
                            ,0.0
                            # Scenario / config file 37
                            ,tds_df[12,:sc_median_3TD_lower], tds_df[12,:sc_icu_8_comb_median_3TD_lower], tds_df[12,:sc_icu_29_comb_median_3TD_lower], tds_df[12,:sc_gp_comb_median_3TD_lower], tds_df[12,:sc_icu_8_gp_comb_median_3TD_lower], tds_df[12,:sc_icu_29_gp_comb_median_3TD_lower]
                            ,0.0
                            # Scenario / config file 49
                            ,tds_df[14,:sc_median_3TD_lower], tds_df[14,:sc_icu_8_comb_median_3TD_lower], tds_df[14,:sc_icu_29_comb_median_3TD_lower], tds_df[14,:sc_gp_comb_median_3TD_lower], tds_df[14,:sc_icu_8_gp_comb_median_3TD_lower], tds_df[14,:sc_icu_29_gp_comb_median_3TD_lower]
                            ,0.0
                            # Scenario / config file 61
                            ,tds_df[16,:sc_median_3TD_lower], tds_df[16,:sc_icu_8_comb_median_3TD_lower], tds_df[16,:sc_icu_29_comb_median_3TD_lower], tds_df[16,:sc_gp_comb_median_3TD_lower], tds_df[16,:sc_icu_8_gp_comb_median_3TD_lower], tds_df[16,:sc_icu_29_gp_comb_median_3TD_lower]
                            ]
upper_errors_td3 = [ # Scenario / config file 1
                            tds_df[1,:sc_median_3TD_upper], tds_df[1,:sc_icu_8_comb_median_3TD_upper], tds_df[1,:sc_icu_29_comb_median_3TD_upper], tds_df[1,:sc_gp_comb_median_3TD_upper], tds_df[1,:sc_icu_8_gp_comb_median_3TD_upper], tds_df[1,:sc_icu_29_gp_comb_median_3TD_upper]
                            ,0.0
                            # Scenario / config file 13
                            ,tds_df[8,:sc_median_3TD_upper], tds_df[8,:sc_icu_8_comb_median_3TD_upper], tds_df[8,:sc_icu_29_comb_median_3TD_upper], tds_df[8,:sc_gp_comb_median_3TD_upper], tds_df[8,:sc_icu_8_gp_comb_median_3TD_upper], tds_df[8,:sc_icu_29_gp_comb_median_3TD_upper]
                            ,0.0
                            # Scenario / config file 25
                            ,tds_df[10,:sc_median_3TD_upper], tds_df[10,:sc_icu_8_comb_median_3TD_upper], tds_df[10,:sc_icu_29_comb_median_3TD_upper], tds_df[10,:sc_gp_comb_median_3TD_upper], tds_df[10,:sc_icu_8_gp_comb_median_3TD_upper], tds_df[10,:sc_icu_29_gp_comb_median_3TD_upper]
                            ,0.0
                            # Scenario / config file 37
                            ,tds_df[12,:sc_median_3TD_upper], tds_df[12,:sc_icu_8_comb_median_3TD_upper], tds_df[12,:sc_icu_29_comb_median_3TD_upper], tds_df[12,:sc_gp_comb_median_3TD_upper], tds_df[12,:sc_icu_8_gp_comb_median_3TD_upper], tds_df[12,:sc_icu_29_gp_comb_median_3TD_upper]
                            ,0.0
                            # Scenario / config file 49
                            ,tds_df[14,:sc_median_3TD_upper], tds_df[14,:sc_icu_8_comb_median_3TD_upper], tds_df[14,:sc_icu_29_comb_median_3TD_upper], tds_df[14,:sc_gp_comb_median_3TD_upper], tds_df[14,:sc_icu_8_gp_comb_median_3TD_upper], tds_df[14,:sc_icu_29_gp_comb_median_3TD_upper]
                            ,0.0
                            # Scenario / config file 61
                            ,tds_df[16,:sc_median_3TD_upper], tds_df[16,:sc_icu_8_comb_median_3TD_upper], tds_df[16,:sc_icu_29_comb_median_3TD_upper], tds_df[16,:sc_gp_comb_median_3TD_upper], tds_df[16,:sc_icu_8_gp_comb_median_3TD_upper], tds_df[16,:sc_icu_29_gp_comb_median_3TD_upper]
                            ] .- medians_td3

colors = vcat(fill(:blue, 7), [:green, :lightgreen, :red]
            ,:black
            ,fill(:blue, 2), [:green, :lightgreen, :red]
            ,:black
            ,fill(:blue, 2), [:green, :lightgreen, :red]
            ,:black
            ,fill(:blue, 2), [:green, :lightgreen, :red]
            ,:black
            ,fill(:blue, 2), [:green, :lightgreen, :red]
            ,:black
            ,fill(:blue, 2), [:green, :lightgreen, :red])

# Plot with error bars
# (CI using samples from fitted distribution probably not that useful here)
p_scatter = scatter(
    medians_td3
    ,1:length(y_labels)
    ,xerror = (lower_errors_td3, upper_errors_td3)
    #, markercolor = colors
    #, markerstrokecolor = colors
    #, color = colors
    ,seriestype = :scatter
    ,yticks = (1:length(y_labels), y_labels)
    ,xlabel = "Median time to detection (days)\n (with 95% CI via 1000 bootstraps)"
    #,xlims = [40,90] #[0,110]
    ,xlims = [10,60] #[0,110]
    ,xticks = 0:10:100
    ,ylabel = ""
    #,title = "Time to detection\n with different site selection"
    ,legend = false
    ,yflip = true
)
plot(p_scatter, layout = (1, 1)
        , size = (900,600)) #default size is 600,400

Plots.savefig( joinpath( output_folder, "scatter_3TD_combined_examples.png" ) )




#############################################################################################
## Scenario 1 
## - comparing original timport distribution, version using 1% quantile truncation, and Poisson sampling from GLEAM output
using CSV
using DataFrames
sc_tds_1_timport_original = CSV.read( "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1287570_nrep2000_analysis/sc_tds_1.csv", DataFrame)
sc_tds_1_timport_trunc_1pc = CSV.read( "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1540261_nrep330_analysis/sc_tds_1.csv", DataFrame)
sc_tds_1_gleam_poisson = CSV.read( "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1633695_nrep340_analysis/sc_tds_1.csv", DataFrame)
println(sc_tds_1_gleam_poisson)

histogram( sc_tds_1_timport_original[:,:SC_TD] , label = "timport original", normalize = :pdf, alpha = 0.5
			, xlabel = "Time since first importation (days)", ylabel = "PDF"
            ,ylimit=[0,0.06])

histogram!( sc_tds_1_timport_trunc_1pc[:,:SC_TD] , label = "timport t-dist truncated at 1% quantile", normalize = :pdf, alpha = 0.5)
histogram!( sc_tds_1_gleam_poisson[:,:SC_TD] , label = "Poisson sampling from GLEAM\n mean daily imports", normalize = :pdf, alpha = 0.5)
vline!( [median(sc_tds_1_timport_original[:,:SC_TD])], color=:blue
                , label="median = 51.4 days"
                , linewidth = 3) 
vline!( [median(sc_tds_1_timport_trunc_1pc[:,:SC_TD])], color=:red
                , label="median = 33.9 days"
                , linewidth = 3) 
vline!( [median(sc_tds_1_gleam_poisson[:,:SC_TD])], color=:green
                , label="median = 33.8 days"
                , linewidth = 3) 

Plots.savefig( joinpath( "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1287570_nrep2000_analysis/", "hist_HARISS_TD_scenario1_timport_method_compariso.png" ) )
Plots.savefig( joinpath( "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1540261_nrep330_analysis/",  "hist_HARISS_TD_scenario1_timport_method_compariso.png" ) )
Plots.savefig( joinpath( "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1633695_nrep340_analysis/",  "hist_HARISS_TD_scenario1_timport_method_compariso.png" ) )

infectivitytoR(NBPMscape.P.infectivity; nsims = 1000) # Reproduction number