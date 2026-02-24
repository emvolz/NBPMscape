#=
Investigating metagenomic sampling in different health care settings 
(primary (GPs via RCGP RSC), secondary (via HARISS network), and tertiary care (ICU) and combinations thereof)
February 2026

Repeat of 2025_12_HARISS_analysis.jl with updates following code review and also using the GLEAM output to model
importation times (more importations compared with previous method)

Also, timeline analysis as per 2026_01_HARISS_analysis_timelines.jl 

=#

using CSV 
using DataFrames
using Dates
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
# using 'covidlike-1.5.0_w_G_filtered_nrep10.jl' and 'covidlike-1.5.0_w_G_filtered_array100_nrep10.pbs'
# Code copied below for reference (with some commented out code removed for clarity)
# Each simulation run had 10 replicates but a PBS array (covidlike-1.5.0_w_G_filtered_array100_nrep10.pbs) was 
# used to repeat this 100 times across three array runs, generating a total of 
# 1000 simulation replicates (10 sim replicates x 100 runs in array)
# Using the HPC with 1 CPU and 32GB it took less than 20 minutes per sim (3.25hrs (195 mins) for the longest batch of 10 sims)

"""
#= 
- Simulate multiple replicates using default covid-like parameters 
- Run for final analysis for research article
-- after changes made following code review 
-- importmodel=:GLEAM_outbreak_country_not_specified in simforest()
-- covid-like parameters
-- maxtime increased to 70
=#
using Pkg
Pkg.add("JLD2")
Pkg.add("DataFrames")
Pkg.add(url="https://github.com/emvolz/NBPMscape", rev="84dc5e3419be9ba95afce9f434d569e9fbc87306") # 'master' branch committed at 14:58 16 Feb 2026 (merge pull request #8)
using NBPMscape 
using JLD2 
using DataFrames

initialize_parameters("config/outbreak_params_covid19_like.yaml"); # Real site files are not required for simulating the outbreak. They can be replaced when simulating sampling
# dummy site file names also in core.jl as these are checked when installing NBPMscape

NREPS = 10 #50 #500 #1000 #5000
sims = map( _-> simforest(NBPMscape.P; initialtime=0.0, maxtime=70.0, maxgenerations=100, max_cases=100000000, importmodel=:GLEAM_outbreak_country_not_specified), 1:NREPS )
@save "covidlike-1.5.0-sims-nrep10.jld2" sims

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
# Create df to save the tinf values for each simulation replicate
sims_tinf_df = [Vector() for _ in 1:length(sims)]

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
	
	# Save times of infection directly from G. 
	# sims_simid_tinf_df only includes unique combinations of simid and tinf so if more than one infection happens at the 
	# same time in the same simtree then only one will be recorded
	try
		fo = sims[s]
		tinfs = fo.G[:,:tinf]
		if length(tinfs) > 0 
			sims_tinf_df[s] = tinfs
		else
			sims_tinf_df[s] = missing
		end
	catch err
		@warn "Error in iteration" exception = err
		println("Error in sim number: ", s)
		continue
	end
	
end

#Save filtered data to files
@save "covidlike-1.5.0-sims-nrep10_filtered_G.jld2" sims_G_filtered
@save "covidlike-1.5.0-sims-nrep10_max_cases_df.jld2" sims_max_cases_df
@save "covidlike-1.5.0-sims-nrep10_sims_simid_tinf_df.jld2" sims_simid_tinf_df
@save "covidlike-1.5.0-sims-nrep10_sims_tinf_df.jld2" sims_tinf_df
"""


## Load simulated outbreak data
sims = combine_sim_reps(;  # sim_input_folders in vector format where multiple paths
                          sim_input_folders = [ "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2026_02_final_run_for_paper/1717144/1717144_sims_G_filtered/"
                                              , "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2026_02_final_run_for_paper/1719024/1719024_sims_G_filtered/"
                                              , "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2026_02_final_run_for_paper/1719029/1719029_sims_G_filtered/" ]
                        , sim_object_name = "sims_G_filtered" #"sims" # "sims_G_icu_filter" # "sims_G_gp_filter" # "sim_G_filtered"
                        , nrep = 10
                        )
# Check df that shows errors in the combination process
sims_mismatch_df = sims[2]
# df is empty

# Refine sims
sims = sims[1]

# Save to file
@save "covidlike-1.5.0-sims-filtered_G_nrep1000_1717144_1719024_1719029.jld2" sims

# Reload file if necessary
#sims = load("covidlike-1.5.0-sims-filtered_G_nrep1000_1717144_1719024_1719029.jld2", "sims")

# Define path to output folder for analyses files
output_folder = "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2026_02_final_run_for_paper/1717144_1719024_1719029_analysis/"

## Compute times to detection for secondary care sampling (via HARISS network)
# using different model parameters representing different scenarios
# See config_file_record.xlsx for summary descriptions of scenarios

# Run analysis for secondary care sampling (via HARISS network)
# Analysis for each set of parameters takes ~30-60 minutes
for i in [1,2,3,4,5,6,7,13,19,25,31,37,43,49,55,61,67] # i=1
    initialize_parameters("config/2026_02_18/covid19_like_params_$(i).yaml")
    sc_tds = secondary_care_td(sims = sims)
    CSV.write( joinpath( output_folder,"sc_tds_$(i).csv" ), sc_tds )
end

# Run analysis for 12 different scenarios (2 seasons x 2 lists of ICU sites x 3 numbers of samples) using ICUs and GPs
# ~1 minute per run on laptop
for j in [1,13,25,37,49,61,73,74,75,76,77,78]
    initialize_parameters("config/2026_02_18/covid19_like_params_$(j).yaml")
    icu_tds = icu_td(sims = sims)
    gp_tds = gp_td(sims = sims)
    CSV.write( joinpath( output_folder,"icu_tds_$(j).csv" ), icu_tds )
    CSV.write( joinpath( output_folder,"gp_tds_$(j).csv" ), gp_tds )
end



## Read in sc_tds files and compute median TD values
using CSV
using DataFrames

# Directory and base file name
path_name = "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2026_02_final_run_for_paper/1717144_1719024_1719029_analysis/"
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
    file_gp = joinpath(path_name, "gp_tds_$(string(k)).csv")

    # Read CSV into a DataFrame and store in the dictionary with the scenario number (numeric suffix) as key
    icu_tds_df_dict[k] = CSV.read(file_icu, DataFrame)
    gp_tds_df_dict[k] = CSV.read(file_gp, DataFrame)
                        
    # Compute times to detection and populate df
    icu_gp_tds_df[i,"Scenario_number"] = k
    icu_gp_tds_df[i,"icu_median_TD"]  = median(icu_tds_df_dict[k][:,"ICU_TD"])
    icu_gp_tds_df[i,"icu_median_3TD"] = median(icu_tds_df_dict[k][:,"ICU_3TD"])
    icu_gp_tds_df[i,"gp_median_TD"]   = median(gp_tds_df_dict[k][:,"GP_TD"])
    icu_gp_tds_df[i,"gp_median_3TD"]  = median(gp_tds_df_dict[k][:,"GP_3TD"])
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
                        , xlim, ylim, bins = 20, xticks = 0:5:100 )
    sc_median = median(sc_data)
    icu_median_1 = median(icu_data_1)
    icu_median_2 = median(icu_data_2)
    gp_median = median(gp_data)
    
    p1 = histogram( sc_data, color = :blue, alpha= 0.5, label = "HARISS"
    ,xlimit = xlim
    ,ylimit = ylim
    ,bins = bins
    ,xticks = xticks
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
                ,ylimit = ylim
                ,bins = bins
                ,xticks = xticks)
                #, normalize = :pdf)
    p2 = vline!( [icu_median_1], color=:green
                , label="ICU (8 sites) median"
                , linewidth = 3) 
    # Add tertiary care TD results
    p3 = histogram( icu_data_2, color = :lightgreen, alpha =0.5, label="ICU (29 sites)"
                #, ylabel = "Number of simulation replicates"
                ,xlimit = xlim
                ,ylimit = ylim
                ,bins = bins
                ,xticks = xticks)
                #, normalize = :pdf)
    p3 = vline!( [icu_median_2], color=:lightgreen
                , label="ICU (29 sites) median"
                , linewidth = 3) 
    
    # Add primary care TD results
    p4 = histogram( gp_data, color=:red, alpha =0.5, label = "RCGP"
                    , xlabel = "Time between first import to first detection (days)"
                    ,xlimit = xlim
                    ,ylimit = ylim
                    ,bins = bins
                    ,xticks = xticks)
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
              , xlim = [0,100]
              , ylim = [0,400]
              , bins = 20
              , xticks = 0:5:100)
# Save plot
Plots.savefig( joinpath( output_folder, "histogram_1.png" ) )

## Comparing results from varying sample size

# Results for scenario 37
# Summer, 300 samples, Equal weighting, no age targeting
plot_histogram( sc_data = sc_tds_df_dict[37][:,:SC_TD]
              , icu_data_1 = icu_tds_df_dict[37][:,:ICU_TD]
              , icu_data_2 = icu_tds_df_dict[76][:,:ICU_TD]
              , gp_data = gp_tds_df_dict[37][:,:GP_TD]
              ,xlim = [0,100]
              ,ylim = [0,400]
              ,bins=20
             , xticks = 0:5:100)
Plots.savefig( joinpath( output_folder, "histogram_37.png" ) )

# Results for scenario 13
# Winter, 500 samples, Equal weighting, no age targeting
plot_histogram( sc_data = sc_tds_df_dict[13][:,:SC_TD]
              , icu_data_1 = icu_tds_df_dict[13][:,:ICU_TD]
              , icu_data_2 = icu_tds_df_dict[74][:,:ICU_TD]
              , gp_data = gp_tds_df_dict[13][:,:GP_TD]
              ,xlim = [0,100]
              ,ylim = [0,400]
              ,bins=20
              , xticks = 0:5:100)
Plots.savefig( joinpath( output_folder, "histogram_13.png" ) )

# Results for scenario 25
# Winter, 1000 samples, Equal weighting, no age targeting
plot_histogram( sc_data = sc_tds_df_dict[25][:,:SC_TD]
              , icu_data_1 = icu_tds_df_dict[25][:,:ICU_TD]
              , icu_data_2 = icu_tds_df_dict[75][:,:ICU_TD]
              , gp_data = gp_tds_df_dict[25][:,:GP_TD]
              ,xlim = [0,100]
              ,ylim = [0,400]
              ,bins=20
              , xticks = 0:5:100)
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
xlim = [0,100]
ylim = [0,400]
    
    p1 = histogram( sc_data_1, color = :blue, alpha= 0.5, label = "HARISS - no age target", xlimit = xlim, ylimit = ylim, bins=20, xticks = 0:5:100) #, normalize = :pdf)
    p1 = vline!( [sc_median_1], color=:blue, label="median", linewidth = 3) 
    
    p2 = histogram( sc_data_2, color = :blue, alpha= 0.5, label = "HARISS - 100% adult", xlimit = xlim, ylimit = ylim, bins=20, xticks = 0:5:100) #, normalize = :pdf)
    p2 = vline!( [sc_median_2], color=:blue, label="median", linewidth = 3) 
    
    p3 = histogram( sc_data_3, color = :blue, alpha= 0.5, label = "HARISS - 75% adult", xlimit = xlim, ylimit = ylim, bins=20, xticks = 0:5:100) #, normalize = :pdf)
    p3 = vline!( [sc_median_3], color=:blue, label="median", linewidth = 3) 
    
    p4 = histogram( sc_data_4, color = :blue, alpha= 0.5, label = "HARISS - 50% adult", xlimit = xlim, ylimit = ylim, bins=20, xticks = 0:5:100
                  , ylabel = "Number of simulation replicates") #, normalize = :pdf)
    p4 = vline!( [sc_median_4], color=:blue, label="median", linewidth = 3) 
    
    p5 = histogram( sc_data_5, color = :blue, alpha= 0.5, label = "HARISS - 25% adult", xlimit = xlim, ylimit = ylim, bins=20, xticks = 0:5:100) #, normalize = :pdf)
    p5 = vline!( [sc_median_5], color=:blue, label="median", linewidth = 3) 
    
    p6 = histogram( sc_data_6, color = :blue, alpha= 0.5, label = "HARISS - 0% adult", xlimit = xlim, ylimit = ylim, bins=20, xticks = 0:5:100
                    , xlabel = "Time between first import to first detection (days)" ) #, normalize = :pdf)
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
xlim = [0,100]
ylim = [0,400]
# 300 samples, equal weight    
p1 = histogram( sc_data_1, color = :blue, alpha= 0.5, label = "HARISS - 300 samples, equal",xlimit = xlim,ylimit = ylim , bins=20, xticks = 0:5:100) #, normalize = :pdf)
p1 = vline!( [sc_median_1], color=:blue, label="median", linewidth = 3)
# 300 samples, weighted
p2 = histogram( sc_data_7, color = :orange, alpha= 0.5, label = "HARISS - 300 samples, weighted", xlimit = xlim, ylimit = ylim, bins=20, xticks = 0:5:100)
p2 = vline!( [sc_median_7], color=:orange, label="median", linewidth = 3)
# 500 samples, equal weight
p3 = histogram( sc_data_13, color = :blue, alpha= 0.5, label = "HARISS - 500 samples, equal", xlimit = xlim, ylimit = ylim, bins=20, xticks = 0:5:100) #, normalize = :pdf)
p3 = vline!( [sc_median_13], color=:blue, label="median", linewidth = 3) 
# 500 samples, weighted
p4 = histogram( sc_data_19, color = :orange, alpha= 0.5, label = "HARISS - 500 samples, weighted", xlimit = xlim, ylimit = ylim, bins=20, xticks = 0:5:100
              , ylabel = "Number of simulation replicates") #, normalize = :pdf)
p4 = vline!( [sc_median_19], color=:orange, label="median", linewidth = 3) 
# 1000 samples, equal weight
p5 = histogram( sc_data_25, color = :blue, alpha= 0.5, label = "HARISS - 1000 samples, equal", xlimit = xlim, ylimit = ylim, bins=20, xticks = 0:5:100) #, normalize = :pdf)
p5 = vline!( [sc_median_25], color=:blue, label="median", linewidth = 3) 
# 1000 samples, weighted
p6 = histogram( sc_data_31, color = :orange, alpha= 0.5, label = "HARISS - 1000 samples, weighted", xlimit = xlim, ylimit = ylim, bins=20, xticks = 0:5:100
                    , xlabel = "Time between first import to first detection (days)" ) #, normalize = :pdf)
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
    ,xlims = [20,60] #[0,110]
    ,xticks = 0:5:100 #0:10:100
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
    ,xlims = [20,60] #[0,110]
    ,xticks = 0:5:100
    ,ylabel = ""
    #,title = "Time to detection\n with different site selection"
    ,legend = false
    ,yflip = true
)
plot(p_scatter_td3, layout = (1, 1)
        , size = (900,600)) #default size is 600,400

Plots.savefig( joinpath( output_folder, "scatter_all_TD3.png" ) )


## Computing and plotting combined times to detection

# Read in files 

using CSV
using DataFrames

# Directory and base file name
path_name = "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2026_02_final_run_for_paper/1717144_1719024_1719029_analysis/"

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
    ,xlims = [20,60] #[0,110]
    ,xticks = 0:5:100
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
    ,xlims = [20,60] #[0,110]
    ,xticks = 0:5:100
    ,ylabel = ""
    #,title = "Time to detection\n with different site selection"
    ,legend = false
    ,yflip = true
)
plot(p_scatter, layout = (1, 1)
        , size = (900,600)) #default size is 600,400

Plots.savefig( joinpath( output_folder, "scatter_3TD_combined_examples.png" ) )


############################
## Timeline plots and tables
############################
# Adapted from code in '2026_01_HARISS_analysis_timelines.jl'

## Scenario 1: 300 samples per week, equal allocation to HARISS sites, winter, no age targeting (see config/HARISS/config file record.xlsx for details of scenarios)

# Secondary sampling care via HARISS (Poisson sampling from GLEAM mean daily import output)
poisson_gleam_timport_med_td = median_td_sim_reps(; tds_file = "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2026_02_final_run_for_paper/1717144_1719024_1719029_analysis/sc_tds_1.csv"
                                                , td_column_name = "SC_TD"
                                                , nreps = 10
                                                )
#(median_td = 34.37471793483482, median_td_lower = 34.371756928670074, median_td_upper = 34.37767894099957
#, median_td_lower_sim_rep_n = 385, median_td_upper_sim_rep_n = 395, median_td_sim_rep_n = missing
#, file_n_med_lower = 39, sim_rep_n_in_file_n_med_lower = 5
#, file_n_med_upper = 40, sim_rep_n_in_file_n_med_upper = 5
#, file_n_med = missing, sim_rep_n_in_file_n_med = missing)

# Check which raw files these sim reps relate to 
sim_input_folders = [ "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2026_02_final_run_for_paper/1717144/1717144_sims_G_filtered/"
                    , "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2026_02_final_run_for_paper/1719024/1719024_sims_G_filtered/"
                    , "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2026_02_final_run_for_paper/1719029/1719029_sims_G_filtered/"]
sim_files = []
for i in 1:length(sim_input_folders) 
    sim_files_temp = readdir(sim_input_folders[i]; join=true)
    sim_files = vcat(sim_files, sim_files_temp)
end

println( sim_files[poisson_gleam_timport_med_td.file_n_med_lower] ) 
# C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2026_02_final_run_for_paper/1719024/1719024_sims_G_filtered/covidlike-1.5.0-sims-nrep10_filtered_G_1719024.61.jld2
println( sim_files[poisson_gleam_timport_med_td.file_n_med_upper] ) 
# C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2026_02_final_run_for_paper/1719024/1719024_sims_G_filtered/covidlike-1.5.0-sims-nrep10_filtered_G_1719024.62.jld2

cases_before_td_v2(; sims_tinf_df_file = "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2026_02_final_run_for_paper/1719024/covidlike-1.5.0-sims-nrep10_sims_tinf_df_1719024.61.jld2"
            , sim_rep_n = poisson_gleam_timport_med_td.sim_rep_n_in_file_n_med_lower # 5
            , td_value  = poisson_gleam_timport_med_td.median_td_lower # 34.371756928670074
             )
# 886

cases_before_td_v2(; sims_tinf_df_file = "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2026_02_final_run_for_paper/1719024/covidlike-1.5.0-sims-nrep10_sims_tinf_df_1719024.62.jld2"
            , sim_rep_n = poisson_gleam_timport_med_td.sim_rep_n_in_file_n_med_upper # 5
            , td_value = poisson_gleam_timport_med_td.median_td_upper # 34.37767894099957
             )
# 3830

# Read files containing TDs
sc_tds_1_1717144_1719024_1719029  = CSV.read( "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2026_02_final_run_for_paper/1717144_1719024_1719029_analysis/sc_tds_1.csv", DataFrame)
gp_tds_1_1717144_1719024_1719029  = CSV.read( "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2026_02_final_run_for_paper/1717144_1719024_1719029_analysis/gp_tds_1.csv", DataFrame)
icu_tds_1_1717144_1719024_1719029 = CSV.read( "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2026_02_final_run_for_paper/1717144_1719024_1719029_analysis/icu_tds_1.csv", DataFrame)

## Lower median sim rep
tds_1_1717144_1719024_1719029 = DataFrame(
                                            sc_td   =  sc_tds_1_1717144_1719024_1719029[ poisson_gleam_timport_med_td.median_td_lower_sim_rep_n,:].SC_TD 
                                            , gp_td   =  gp_tds_1_1717144_1719024_1719029[ poisson_gleam_timport_med_td.median_td_lower_sim_rep_n,:].GP_TD 
                                            , icu_td  = icu_tds_1_1717144_1719024_1719029[ poisson_gleam_timport_med_td.median_td_lower_sim_rep_n,:].ICU_TD 
                                            , sc_td3  =  sc_tds_1_1717144_1719024_1719029[ poisson_gleam_timport_med_td.median_td_lower_sim_rep_n,:].SC_3TD 
                                            , gp_td3  =  gp_tds_1_1717144_1719024_1719029[ poisson_gleam_timport_med_td.median_td_lower_sim_rep_n,:].GP_3TD 
                                            , icu_td3 = icu_tds_1_1717144_1719024_1719029[ poisson_gleam_timport_med_td.median_td_lower_sim_rep_n,:].ICU_3TD
) 

## Load files
# load sims_simid_tinf_df
sims_tinf_df_1719024_61 = load( "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2026_02_final_run_for_paper/1719024/covidlike-1.5.0-sims-nrep10_sims_tinf_df_1719024.61.jld2"
                                    , "sims_tinf_df") # 39th file in folders
sims_tinf_df_1719024_61_5_385 = sims_tinf_df_1719024_61[5]

# load filtered G df
sims_1719024_61 = load( "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2026_02_final_run_for_paper/1719024/1719024_sims_G_filtered/covidlike-1.5.0-sims-nrep10_filtered_G_1719024.61.jld2"
                        , "sims_G_filtered" )
sims_1719024_61_5_385 = sims_1719024_61[5] # HPC run 61, sim rep 5 (in HPC run 61), sim rep 385 in combined .jld2 file, 39th file in folder used to combine

# Load max_cases_df
sims_max_cases_df_1719024_61 = load( "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2026_02_final_run_for_paper/1719024/1719024_max_cases_df/covidlike-1.5.0-sims-nrep10_max_cases_df_1719024.61.jld2", "sims_max_cases_df" )
sims_max_cases_df_1719024_61_5_385 = sims_max_cases_df_1719024_61[5] # HPC run 61, sim rep 5 (in HPC run 61), sim rep 385 in combined .jld2 file, 39th file in folder used to combine

# Run function
inf_t_analysis_1719024_61_5_385 = infections_time_analysis_v2(; sims_tinf_df_object = sims_tinf_df_1719024_61_5_385
                                                                , sims_G_filtered_object = sims_1719024_61_5_385
                                                                , sims_max_cases_df_object = sims_max_cases_df_1719024_61_5_385
                                                                #, sim_rep_n
                                                                , td_gp = tds_1_1717144_1719024_1719029.gp_td
                                                                , td_sc = tds_1_1717144_1719024_1719029.sc_td
                                                                , td_icu = tds_1_1717144_1719024_1719029.icu_td
                                                                , td3_gp = tds_1_1717144_1719024_1719029.gp_td3
                                                                , td3_sc = tds_1_1717144_1719024_1719029.sc_td3
                                                                , td3_icu = tds_1_1717144_1719024_1719029.icu_td3
                                                                , base_date = Date(2020, 1, 15)
                                                                , maxtime = 70 #days
                                                                )

# Import count range                                                            
minimum( inf_t_analysis_1719024_61_5_385.numbers_by_day_df.import_count )
maximum( inf_t_analysis_1719024_61_5_385.numbers_by_day_df.import_count )

# Count outputs
inf_t_analysis_1719024_61_5_385.number_at_td
#8×4 DataFrame
# Row │ Description                  TD_via_GP  TD_via_SC  TD_via_ICU 
#     │ String                       Int64      Int64      Int64
#─────┼───────────────────────────────────────────────────────────────
#   1 │ All_infections                   65429        886       11980
#   2 │ GP_visits                         3290         45         501
#   3 │ ED_visits                         1591         25         239
#   4 │ Hospital_admissions               1324         23         211
#   5 │ ICU_admissions                     156          4          24
#   6 │ Deaths                             126          7          24
#   7 │ Infections_leading_to_death        836         17         164
#   8 │ Imported_infections               6927        272        3331

inf_t_analysis_1719024_61_5_385.number_at_td3
#8×4 DataFrame
# Row │ Description                  TD3_via_GP  TD3_via_SC  TD3_via_ICU 
#     │ String                       Int64       Int64       Int64
#─────┼──────────────────────────────────────────────────────────────────
#   1 │ All infections                    76138        3574        62677
#   2 │ GP visits                          3901         157         3162
#   3 │ ED visits                          1882          79         1517
#   4 │ Hospital admissions                1549          70         1260
#   5 │ ICU admissions                      184           9          151
#   6 │ Deaths                              151          11          121
#   7 │ Infections_leading_to_death        1002          42          802
#   8 │ Imported_infections                6927        1080         6927

plot_outbreak_analysis_combined(; data = inf_t_analysis_1719024_61_5_385
                                    , time_period = "days", x_labels = "date", limit_x_axis = true
                                    , td = poisson_gleam_timport_med_td.median_td_lower
                                    , td3 = only( tds_1_1717144_1719024_1719029.sc_td3 ) 
                                    , add_import_markers = false)
Plots.savefig( "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2026_02_final_run_for_paper/1717144_1719024_1719029_analysis/sim_analysis_HPC1719024_SC_med_TD_lower_array61_simrep5_simrepcombined385.png" )

# Doubling time estimates
# Local (in time)
inf_t_analysis_1719024_61_5_385.tinf_local_doubling_times_vec
# Global (in time)
inf_t_analysis_1719024_61_5_385.tinf_global_doubling_times_df
#3×2 DataFrame
# Row │ method                          doubling_time_estimates 
#     │ String                          Float64
#─────┼─────────────────────────────────────────────────────────
#   1 │ global_log_linear                               3.60195
#   2 │ global_wtd_log_linear                           4.76801
#   3 │ global_nonlinear_least_squares                  3.60195

plot(inf_t_analysis_1719024_61_5_385.tinf_local_doubling_times_vec, xlabel = "Time since first importation (days)", ylabel = "Estimated doubling time (days)", label = "Local estimate")
hline!( repeat( [ inf_t_analysis_1719024_61_5_385.tinf_global_doubling_times_df[3,:doubling_time_estimates] ], length(inf_t_analysis_1719024_61_5_385.tinf_local_doubling_times_vec) ), label = "Global estimate using nonlinear least squares")
Plots.savefig( "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2026_02_final_run_for_paper/1717144_1719024_1719029_analysis/sim_analysis_HPC1719024_SC_med_TD_lower_array61_simrep5_simrepcombined385_doubling_time.png" )            

## Upper median sim rep
tds_1_1717144_1719024_1719029_median_td_upper = DataFrame(
                                                              sc_td   =  sc_tds_1_1717144_1719024_1719029[ poisson_gleam_timport_med_td.median_td_upper_sim_rep_n,:].SC_TD 
                                                            , gp_td   =  gp_tds_1_1717144_1719024_1719029[ poisson_gleam_timport_med_td.median_td_upper_sim_rep_n,:].GP_TD 
                                                            , icu_td  = icu_tds_1_1717144_1719024_1719029[ poisson_gleam_timport_med_td.median_td_upper_sim_rep_n,:].ICU_TD 
                                                            , sc_td3  =  sc_tds_1_1717144_1719024_1719029[ poisson_gleam_timport_med_td.median_td_upper_sim_rep_n,:].SC_3TD 
                                                            , gp_td3  =  gp_tds_1_1717144_1719024_1719029[ poisson_gleam_timport_med_td.median_td_upper_sim_rep_n,:].GP_3TD 
                                                            , icu_td3 = icu_tds_1_1717144_1719024_1719029[ poisson_gleam_timport_med_td.median_td_upper_sim_rep_n,:].ICU_3TD
) 

## Load files
# load sims_simid_tinf_df
sims_tinf_df_1719024_62 = load( "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2026_02_final_run_for_paper/1719024/covidlike-1.5.0-sims-nrep10_sims_tinf_df_1719024.62.jld2"
                                    , "sims_tinf_df") # 40th file in folder
sims_tinf_df_1719024_62_5_395 = sims_tinf_df_1719024_62[5]

# load filtered G df
sims_1719024_62 = load( "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2026_02_final_run_for_paper/1719024/1719024_sims_G_filtered/covidlike-1.5.0-sims-nrep10_filtered_G_1719024.62.jld2"
                        , "sims_G_filtered" )
sims_1719024_62_5_395 = sims_1719024_62[5] # HPC run 62, sim rep 5 (in HPC run 62), sim rep 395 in combined .jld2 file, 40th file in folder used to combine

# Load max_cases_df
sims_max_cases_df_1719024_62 = load( "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2026_02_final_run_for_paper/1719024/1719024_max_cases_df/covidlike-1.5.0-sims-nrep10_max_cases_df_1719024.62.jld2", "sims_max_cases_df" )
sims_max_cases_df_1719024_62_5_395 = sims_max_cases_df_1719024_62[5] # HPC run 62, sim rep 5 (in HPC run 62), sim rep 395 in combined .jld2 file, 40th file in folder used to combine

# Run function
inf_t_analysis_1719024_62_5_395 = infections_time_analysis_v2(; sims_tinf_df_object = sims_tinf_df_1719024_62_5_395
                                                                , sims_G_filtered_object = sims_1719024_62_5_395
                                                                , sims_max_cases_df_object = sims_max_cases_df_1719024_62_5_395
                                                                #, sim_rep_n
                                                                , td_gp = tds_1_1717144_1719024_1719029_median_td_upper.gp_td
                                                                , td_sc = tds_1_1717144_1719024_1719029_median_td_upper.sc_td
                                                                , td_icu = tds_1_1717144_1719024_1719029_median_td_upper.icu_td
                                                                , td3_gp = tds_1_1717144_1719024_1719029_median_td_upper.gp_td3
                                                                , td3_sc = tds_1_1717144_1719024_1719029_median_td_upper.sc_td3
                                                                , td3_icu = tds_1_1717144_1719024_1719029_median_td_upper.icu_td3
                                                                , base_date = Date(2020, 1, 15)
                                                                , maxtime = 70 #days
                                                                )

# Import count range                                                            
minimum( inf_t_analysis_1719024_62_5_395.numbers_by_day_df.import_count )
maximum( inf_t_analysis_1719024_62_5_395.numbers_by_day_df.import_count )

# Count outputs
inf_t_analysis_1719024_62_5_395.number_at_td
#8×4 DataFrame
# Row │ Description                  TD_via_GP  TD_via_SC  TD_via_ICU 
#     │ String                       Int64      Int64      Int64
#─────┼───────────────────────────────────────────────────────────────
#   1 │ All_infections                   15650       3830        3724
#   2 │ GP_visits                          726        161         155
#   3 │ ED_visits                          355         83          79
#   4 │ Hospital_admissions                261         61          61
#   5 │ ICU_admissions                      36          9           9
#   6 │ Deaths                              28          3           3
#   7 │ Infections_leading_to_death        189         53          53
#   8 │ Imported_infections               2830        656         645

inf_t_analysis_1719024_62_5_395.number_at_td3
#8×4 DataFrame
# Row │ Description                  TD3_via_GP  TD3_via_SC  TD3_via_ICU 
#     │ String                       Int64       Int64       Int64
#─────┼──────────────────────────────────────────────────────────────────
#   1 │ All infections                    38249        3830        30239
#   2 │ GP visits                          1717         161         1368
#   3 │ ED visits                           861          83          680
#   4 │ Hospital admissions                 622          61          502
#   5 │ ICU admissions                       82           9           64
#   6 │ Deaths                               71           3           59
#   7 │ Infections_leading_to_death         453          53          367
#   8 │ Imported_infections                6713         656         5422

plot_outbreak_analysis_combined(; data = inf_t_analysis_1719024_62_5_395
                                    , time_period = "days", x_labels = "date", limit_x_axis = true
                                    , td = poisson_gleam_timport_med_td.median_td_upper
                                    , td3 = only( tds_1_1717144_1719024_1719029.sc_td3 ) 
                                    , add_import_markers = false)
Plots.savefig( "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2026_02_final_run_for_paper/1717144_1719024_1719029_analysis/sim_analysis_HPC1719024_SC_med_TD_upper_array62_simrep5_simrepcombined395.png" )

# Doubling time estimates
# Local (in time)
inf_t_analysis_1719024_62_5_395.tinf_local_doubling_times_vec
# Global (in time)
inf_t_analysis_1719024_62_5_395.tinf_global_doubling_times_df
#3×2 DataFrame
# Row │ method                          doubling_time_estimates 
#     │ String                          Float64
#─────┼─────────────────────────────────────────────────────────
#   1 │ global_log_linear                               3.78959
#   2 │ global_wtd_log_linear                           5.10707
#   3 │ global_nonlinear_least_squares                  5.562

plot(inf_t_analysis_1719024_62_5_395.tinf_local_doubling_times_vec, xlabel = "Time since first importation (days)", ylabel = "Estimated doubling time (days)", label = "Local estimate")
hline!( repeat( [ inf_t_analysis_1719024_62_5_395.tinf_global_doubling_times_df[3,:doubling_time_estimates] ], length(inf_t_analysis_1719024_62_5_395.tinf_local_doubling_times_vec) ), label = "Global estimate using nonlinear least squares")
Plots.savefig( "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2026_02_final_run_for_paper/1717144_1719024_1719029_analysis/sim_analysis_HPC1719024_SC_med_TD_upper_array62_simrep5_simrepcombined395_doubling_time.png" )            