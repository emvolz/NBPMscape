#=
Part 1 of a series of scripts for the investigation of the effect on surveillance performance under
different selection strategies for sampling sites. Surveillance performance determined by times 
to detection for 1 case (TD) and 3 cases (3TD) and the cumulative number of infections at 
time of detection.

The different strategies investigated here are:
 - geographic optimisation of sentinels using a hierarchical clustering &
   geographic stratification (HCGS) strategy.
 - NHS Trusts in London only;
 - NHS Trusts in regions with the largest population sizes; and
 - mSCAPE pilot study sites.

Outline of the 'Geographic optimisation analysis':
-- Part 1-------------------------------------------------------------------------------------
(1) Description of outbreak simulation
(2) Manipulate individual simulation replicate files
(3) Compile outbreak metrics
(4) Sampling infections to simulate TDs using different sampling strategies:
    (4.1) HCGS
        (4.1.1) Compute TDs at regional level with 100% of infections sampled
        (4.1.2) Compute correlations between TD distributions (1000 simulation replicates) for each region
-- Part 1 & Part 2-------------------------------------------------------------------------------------        
        (4.1.3) Perform hierarchical clustering to group regions together based on TD distribution correlations
                and allocate a number of sites to each region. Regions with lower median TD are prioritised and
                sites are distributed across both groups and regions.
-- Part 2-------------------------------------------------------------------------------------                        
        (4.1.4) Select sampling sites
        (4.1.5) Simulate sampling and TD
        (4.1.6) Add outbreak metrics relating to TDs
        (4.1.7) Plot results
    (4.2) London only
        (4.2.1) Select sampling sites
        (4.2.2) Simulate sampling and TD
        (4.2.3) Add outbreak metrics relating to TDs
        (4.2.4) Plot results
    (4.3) Population size
        (4.3.1) Select sampling sites
        (4.3.2) Simulate sampling and TD
        (4.3.3) Add outbreak metrics relating to TDs
        (4.3.4) Plot results
    (4.4) mSCAPE pilot sites
        (4.4.1) Select sampling sites
        (4.4.2) Simulate sampling and TD
        (4.4.3) Add outbreak metrics relating to TDs
        (4.4.4) Plot results
(5) Compile TD statistics for all strategies and parameter combinations
(6) Comparison of best ICU sampling performers
    (6.1) Plot best ICU sampling performers: TD and cumulative infections at TD
    (6.2) Compute statistics for relative performance
    (6.3) Report 95th percentile of TD and cumulative infections at TD
    (6.4) mSCAPE sampling equivalent to 8 site HCGS
    (6.5) Uncertainty in variance
    (6.6) Additional comparisons
(7) Plot TD variance vs TD 
    (7.1) for all sampling scenarios and strategies
    (7.2) Function for plotting
    (7.3) Plot HCGS sampling scenarios only, but split by samples per week and season
(8) Look at cumulative GP attendances at TD re syndromic surveillance
(9) Ranking sampling strategies using mean-risk objective optimisation
    (Conditional value-at-risk (CVaR) and loss function using variance as risk measure)
=#

#-----------------------------------
# Activate environment for analysis
#-----------------------------------
# Unless otherwise stated, run on Laptop - 11th Gen Intel(R) Core(TM) i7-1185G7 @ 3.00GHz (1.80 GHz) with 32.0 GB RAM
# Using Julia v1.12.3 in VS Code 1.141.0

using Pkg
Pkg.activate(@__DIR__) # ...NBPMscape\\paper\\analysis\\2_sampling_analysis\\scripts

# Check packages in use
# Using NBPMscape version from commit "c4d186da30f208564e28f53922d9db22cfe9824d" as at 15:17 15 Oct 2026
Pkg.status()

# Check the environment location
Base.active_project()
#"...NBPMscape\\paper\\analysis\\2_sampling_analysis\\scripts\\Project.toml"

#-----------------------------------
# Load required packages
#-----------------------------------
Pkg.instantiate()
using CSV
using DataFrames
using Distributions
using GLM
using JLD2
using NBPMscape # Should be version from commit "bb582cce590d1daa342833887ed13503d0bb8c34" as at 18:53 12 Oct 2026
using Plots
using StatsBase
using StatsPlots

############################################################################################
# (1) Description of outbreak simulation
############################################################################################

#=
    Simulation run on Imperial College High Performance Computing (HPC) system
    using a Julia script and PBS file to manage the process:

    Scripts: Julia - 'covidlike-1.6.1_w_G_filtered_nrep10.jl'
             PBS   - 'covidlike-1.6.1_w_G_filtered_array100_nrep10.pbs'
            (see NBPMscape/paper/analysis/1_outbreak_simulations/covid_like/scripts)
    HPC run numbers: 4300448
    Simulation 'maxtime' = 70 days
    Simulation replicates per PBS array run: 10
    Runs in PBS array: 100
    Total sim reps completed: 1000
    CPUs per array run: 1
    RAM per array run: 64GB
    Longest batch run time: ~4.1hrs (<24 minutes per sim rep)
=#

############################################################################################
# (2) Manipulate individual simulation replicate files
############################################################################################

# Define root directory
#root_dir = root_folder() # Function defined in NBPMscape/src/misc_functions
dir = @__DIR__
while !isdir(joinpath(dir, "NBPMscape"))
    parent = dirname(dir)
    parent == dir && error("NBPmscape not found")
    dir = parent
end
root_dir = dir

# Combine simid_tinf_df from simulation replicates
sims_simid_tinf_df = combine_sim_reps(;  # sim_input_folders in vector format where multiple paths
                                          sim_input_folders = [ normpath( joinpath( root_dir, "NBPMscape/paper/analysis/1_outbreak_simulations/covid_like/output/4300448_sims_simid_tinf_df/") )]
                                        , sim_object_name = "sims_simid_tinf_df" #"sims_max_cases_df" #"sims_G_filtered" #"sims" # "sims_G_icu_filter" # "sims_G_gp_filter" # "sim_G_filtered"
                                        , nrep = 10
                                        )
# Save to file
sims_simid_tinf_df = sims_simid_tinf_df[1]
@save normpath( joinpath( root_dir, "NBPMscape/paper/analysis/1_outbreak_simulations/covid_like/output/covidlike-1.6.1-sims_simid_tinf_df_nrep1000_4300448.jld2") ) sims_simid_tinf_df

# Reload file if necessary
#sims_simid_tinf_df = load( normpath( joinpath( root_dir, "NBPMscape/paper/analysis/1_outbreak_simulations/covid_like/output/covidlike-1.6.1-sims_simid_tinf_df_nrep1000_4300448.jld2") )
#                            , "sims_simid_tinf_df")

# Extract time series of outbreak metrics
tinf_vals = [Float64[] for _ in 1:length(sims_simid_tinf_df)] 
for i in 1:length(sims_simid_tinf_df)
    temp = sort(sims_simid_tinf_df[i][:,:tinf])
    #push!(temp , tinf_vals)
    tinf_vals[i] = temp
end

@save normpath( joinpath( root_dir, "NBPMscape/paper/analysis/1_outbreak_simulations/covid_like/output/covidlike-1.6.1_nrep1000_4300448_tinf_vals.jld2") ) tinf_vals
# Reload file if necessary
#tinf_vals = load( normpath( joinpath( root_dir, "NBPMscap/paper/analysis/1_outbreak_simulations/covid_like/output/covidlike-1.6.1_nrep1000_4300448_tinf_vals.jld2")), "tinf_vals")

# Reduce Julia memory usage
sims_simid_tinf_df = nothing
GC.gc()

### Load and combine simulated outbreak data
sims = combine_sim_reps(;  # sim_input_folders in vector format where multiple paths
                        sim_input_folders = [ normpath( joinpath( root_dir
                          , "NBPMscape/paper/analysis/1_outbreak_simulations/covid_like/output/4300448_sims_G_filtered/") ) ]
                        , sim_object_name = "sims_G_filtered" #"sims" # "sims_G_icu_filter" # "sims_G_gp_filter" # "sim_G_filtered"
                        , nrep = 10
                        )
# Check df that shows errors in the combination process
sims_mismatch_df = sims[2]
# df is empty

# Refine sims
sims = sims[1]

# Save to file
@save normpath( joinpath( root_dir, "NBPMscape/paper/analysis/1_outbreak_simulations/covid_like/output/covidlike-1.6.1-sims_G_filtered_nrep1000_4300448.jld2") ) sims

# Reload file if necessary
#sims = load( normpath( joinpath( root_dir, "NBPMscape/paper/analysis/1_outbreak_simulations/covid_like/output/covidlike-1.6.1-sims_G_filtered_nrep1000_4300448.jld2") )
#            , "sims")

# Combine sims_max_cases_df from simulation replicates
sims_max_cases_df = combine_sim_reps(;  # sim_input_folders in vector format so can include multiple paths
                                          sim_input_folders = [ normpath( joinpath( root_dir, "NBPMscape/paper/analysis/1_outbreak_simulations/covid_like/output/4300448_sims_max_cases_df/") )]
                                        , sim_object_name = "sims_max_cases_df" #"sims_simid_tinf_df" #"sims_G_filtered" #"sims" # "sims_G_icu_filter" # "sims_G_gp_filter" # "sim_G_filtered"
                                        , nrep = 10
                                        )

# Save to file
sims_max_cases_df[2] # Empty dataframe, so no errors in combining
sims_max_cases_df = sims_max_cases_df[1]
@save normpath( joinpath( root_dir, "NBPMscape/paper/analysis/1_outbreak_simulations/covid_like/output/covidlike-1.6.1-sims_max_cases_df_nrep1000_4300448.jld2") ) sims_max_cases_df

# Reload file if necessary
#sims_max_cases_df = load( normpath( joinpath( root_dir, "NBPMscape/paper/analysis/1_outbreak_simulations/covid_like/output/covidlike-1.6.1-sims_max_cases_df_nrep1000_4300448.jld2") ), "sims_max_cases_df")


# Define path to output folder for analyses files
geo_optim_folder = normpath( joinpath( root_dir, "NBPMscape/paper/analysis/2_sampling_analysis/covid_like/geo_optimisation" ) )

##################################################################################################
# (3) Compile outbreak metrics
##################################################################################################

### Compute outbreak count metrics at TD for different sampling parameter sets

# Compute time series for each outbreak metric count

sims_G_filtered = sims

# Create objects to store time series
vals_format = [Float64[] for _ in 1:length(sims_G_filtered)]
tgp_vals        = copy(vals_format)
ted_vals        = copy(vals_format)
thosp_vals      = copy(vals_format)
ticu_vals       = copy(vals_format)
tdeceased_vals  = copy(vals_format)
tinf_fatal_vals = copy(vals_format)
timport_vals    = copy(vals_format)

for i in 1:length(sims_G_filtered) # sims_G_filtered = sims # i=1
    tgp_vals[i]        = sort( sims_G_filtered[i][ isfinite.( sims_G_filtered[i].tgp       ), : ][:,:tgp      ] )
    ted_vals[i]        = sort( sims_G_filtered[i][ isfinite.( sims_G_filtered[i].ted       ), : ][:,:ted      ] )
    thosp_vals[i]      = sort( sims_G_filtered[i][ isfinite.( sims_G_filtered[i].thospital ), : ][:,:thospital] )
    ticu_vals[i]       = sort( sims_G_filtered[i][ isfinite.( sims_G_filtered[i].ticu      ), : ][:,:ticu     ] )
    tdeceased_vals[i]  = sort( sims_G_filtered[i][ isfinite.( sims_G_filtered[i].tdeceased ), : ][:,:tdeceased] )
    tinf_fatal_vals[i] = sort( sims_G_filtered[i][ (sims_G_filtered[i].fatal) , : ][:,:tinf] ) # Note: this may miss some infections because these are from the df filtered for infections with healthcare setting times, i.e. there may be some fatal infections that not entered healthcare yet
    timport_vals[i]    = sort( sims_max_cases_df[i][ isfinite.( sims_max_cases_df[i].timport ), : ][:,:timport] )
end

outbreak_metrics_time_series = ( tinf_vals = tinf_vals
                                , tgp_vals = tgp_vals
                                , ted_vals = ted_vals
                                , thosp_vals = thosp_vals
                                , ticu_vals = ticu_vals
                                , tdeceased_vals = tdeceased_vals
                                , tinf_fatal_vals = tinf_fatal_vals
                                , timport_vals = timport_vals )

sum(outbreak_metrics_time_series.tinf_vals[1]) # 3.867453029824505e7
outbreak_metrics_time_series.tgp_vals[1]==outbreak_metrics_time_series.ted_vals[1] # false

@save normpath( joinpath( root_dir, "NBPMscape/paper/analysis/1_outbreak_simulations/covid_like/output/covidlike-1.6.1_nrep1000_4300448_outbreak_metrics_time_series.jld2") ) outbreak_metrics_time_series


###########################################################################################################
#=
(4) Sampling infections to simulate TDs using different sampling strategies:
    (4.1) HCGS
        (4.1.1) Compute TDs at regional level with 100% of infections sampled and
                compute correlations between TD distributions (1000 simulation replicates) for each region
=#
###########################################################################################################

# Compute time to detection via ICU for each ITL2 region assuming every true positive case is sampled, i.e. the
# optimum time to detection in each region using sampling in ICUs

# How many replicates are in this simulation
n_replicates = length(sims) 

# Initialize parameters using defaults
initialize_parameters("config/outbreak_params_covid19_like_20261008.yaml")

# ITL2 region list
REGKEY.code

# Sampling proportion within region
psampled = 1 # Using 1 here means that we will sample all infections reaching ICU in each region

#length(REGKEY.code)
#s=1

# Create dictionary to hold the simulated time of infection for each ITL2 region
tinf_by_region_dict         = Dict(name => Float64[] for name in REGKEY.code)
tsample_by_region_dict      = Dict(name => Float64[] for name in REGKEY.code)
tinf_tsample_by_region_dict = Dict(name => Float64[] for name in REGKEY.code)
tinf_treport_by_region_dict = Dict(name => Float64[] for name in REGKEY.code)
treport_by_region_dict      = Dict(name => Float64[] for name in REGKEY.code)

# Create df to hold the number of ICU cases by region in each simulation replicate
n_icu_by_region_by_simrep = DataFrame( zeros( length(sims) , length(REGKEY.code)), REGKEY.code )

# Create df to record the TD by simulation replicate and by region
td_by_simrep_by_region = DataFrame( [fill(Inf, length(sims)) for _ in 1:length(REGKEY.code)], REGKEY.code)


# Loop through each simulation replicate and compile information
for s in 1:length(sims) #s=1

    println("Processing simulation replicate $(s) of $(length(sims))")
            
    # Filter for ICU cases 
    fo = sims[s]
    G = fo[ isfinite.(fo.ticu), : ]
    
    #println( names(G) )
    
    # Check if G contains any data
    if size(G,1) == 0
        continue
    end

    # Record number of ICU cases by region and by simulation
    icu_count_by_region = countmap(G.homeregion)
    
    for (k, v) in icu_count_by_region
        n_icu_by_region_by_simrep[ s, k ] = Int(v)
    end
    #println(icu_count_by_region)
    #println(n_icu_by_region_by_simrep[1,:])

    # List regions with an ICU admission (i.e. possibility of getting a sample and detection)
    icu_regions = unique(G.homeregion)
    
    # Loop through each region and record time of infection stats
    for r in icu_regions #REGKEY.code
        #Test
        #r=icu_regions[1]
        
        g_region = G[ (G.homeregion).==r, :]
        
        # Sample size 
	    n = rand( Binomial( size(g_region,1), psampled) ) #p.psampled ) )
	    # Subsample of ICU cases within region
	    g_region_sub = g_region[sample( 1:size(g_region,1), n, replace=false ), :]

        # Check if there are any rows before continuing
        if size(g_region_sub,1) == 0 
            continue
        elseif size(g_region_sub,1) > 0 
            
            # Generate sample times
            if NBPMscape.P.icu_only_sample_before_death == true
                # Set upper limit on sample time to the minimum of time of death and ICU admission + max lag
                tsample = map( g -> rand( Uniform( g.ticu[1], min( g.tdeceased[1], g.ticu[1] + NBPMscape.P.icu_swab_lag_max ) ) ) , eachrow(g_region_sub) )
            else
                # Set upper limit on sample time to ICU admission + max lag
                tsample = map( g -> rand( Uniform( g.ticu[1], g.ticu[1] + NBPMscape.P.icu_swab_lag_max)) , eachrow(g_region_sub) )
            end
            g_region_sub.tsample = tsample 
            #println(g_region_sub[:,[2,3,4,5,6,16]])

            # Simulate reports times
            treport = (tsample .+ rand(Uniform( NBPMscape.P.turnaroundtime_icu[1], NBPMscape.P.turnaroundtime_icu[2])) )
            
            g_region_sub.treport = treport
            #println(g_region_sub[:,[2,3,4,5,6,16,17]])

            # Find infection (row) with minimum report (detection) time
            min_treport, min_treport_row_index = findmin( g_region_sub.treport )

            ## Add times to respective dictionaries
            # Times of infection for ICU admissions
            push!(tinf_by_region_dict[ r ], g_region_sub.tinf[min_treport_row_index] ) # println(tinf_by_region_dict)
            # Times of sampling for ICU admissions
            push!(tsample_by_region_dict[ r ], g_region_sub.tsample[min_treport_row_index] )
            # Times of report / detection for ICU admissions
            push!(treport_by_region_dict[ r ], g_region_sub.treport[min_treport_row_index] )
            # Time between infection and sampling for ICU admissions
            push!(tinf_tsample_by_region_dict[ r ], g_region_sub.tsample[min_treport_row_index] - g_region_sub.tinf[min_treport_row_index] )
            # Time between infection and reporting of sample results for ICU admissions
            push!(tinf_treport_by_region_dict[ r ], g_region_sub.treport[min_treport_row_index] - g_region_sub.tinf[min_treport_row_index] )
            
            # Record TD (= min(treport)) by simulation replicate and by region
            td_by_simrep_by_region[ s , r ] = min_treport
        end
                
    end

end

# Save 
CSV.write( normpath( joinpath( geo_optim_folder, "icu_td_regional/td_by_simrep_by_region_nrep$(n_replicates)_p$(Int(psampled*100))_tsample_b4_tdeceased.csv") ), td_by_simrep_by_region )

# Check how many tinf = 0 for each region (this is the tinf relating to the earliest treport)
tinf_zero_count_by_region = DataFrame( zeros( 1, length( tinf_by_region_dict ) ), REGKEY.code )
tinf_zero_perc_by_region  = DataFrame( zeros( 1, length( tinf_by_region_dict ) ), REGKEY.code )
for (k, v) in tinf_by_region_dict
    tinf_zero_count_by_region[ 1, k ] =       count(==(0), v)
    tinf_zero_perc_by_region[  1, k ] = 100 * count(==(0), v) / length(v)
end
println("Number of sim reps without a TD in each region")
for col in names(tinf_zero_count_by_region)
        println(rpad(string(col), 12), " | ", Int(tinf_zero_count_by_region[1, col]))
end
println("The minimum number of sim reps without a TD is $(Int(minimum(tinf_zero_count_by_region[1,:])))")
println("The maximum number of sim reps without a TD is $(Int(maximum(tinf_zero_count_by_region[1,:])))")
println("Percentage of sim reps without a TD in each region")
for col in names(tinf_zero_perc_by_region)
        println(rpad(string(col), 12), " | ", tinf_zero_perc_by_region[1, col],"%")
end
println("The minimum percentage of sim reps without a TD is $(minimum(tinf_zero_perc_by_region[1,:]))%")
println("The maximum percentage of sim reps without a TD is $(maximum(tinf_zero_perc_by_region[1,:]))%")

# Check how many treport = 0 for each region
treport_zero_count_by_region = DataFrame( zeros( 1, length( treport_by_region_dict ) ), REGKEY.code )
treport_zero_perc_by_region  = DataFrame( zeros( 1, length( treport_by_region_dict ) ), REGKEY.code )
for (k, v) in treport_by_region_dict
    treport_zero_count_by_region[ 1, k ] =       count(==(0), v)
    treport_zero_perc_by_region[ 1, k ]  = 100 * count(==(0), v) / length(v)
end
println("Number of sim reps with treport = 0 in each region")
for col in names(treport_zero_count_by_region)
        println(rpad(string(col), 12), " | ", Int(treport_zero_count_by_region[1, col]))
end
println("The minimum number of sim reps with treport = 0 is $(Int(minimum(treport_zero_count_by_region[1,:])))")
println("The maximum number of sim reps with treport = 0 is $(Int(maximum(treport_zero_count_by_region[1,:])))")
println("Percentage of sim reps with treport = 0 in each region")
for col in names(treport_zero_perc_by_region)
        println(rpad(string(col), 12), " | ", treport_zero_perc_by_region[1, col],"%")
end
println("The minimum percentage of sim reps with treport = 0 is $(minimum(treport_zero_perc_by_region[1,:]))%")
println("The maximum percentage of sim reps with treport = 0 is $(maximum(treport_zero_perc_by_region[1,:]))%")

# Check how many treport - tinf = 0 for each region
tinf_treport_zero_count_by_region = DataFrame( zeros( 1, length( tinf_treport_by_region_dict ) ), REGKEY.code )
tinf_treport_zero_perc_by_region  = DataFrame( zeros( 1, length( tinf_treport_by_region_dict ) ), REGKEY.code )
for (k, v) in tinf_treport_by_region_dict
    tinf_treport_zero_count_by_region[ 1, k ] =       count(==(0), v)
    tinf_treport_zero_perc_by_region[ 1, k ]  = 100 * count(==(0), v) / length(v)
end
println("Number of sim reps with treport - tinf = 0 in each region")
for col in names(tinf_treport_zero_count_by_region)
        println(rpad(string(col), 12), " | ", Int(tinf_treport_zero_count_by_region[1, col]))
end
println("The minimum number of sim reps with treport - tinf = 0 is $(Int(minimum(tinf_treport_zero_count_by_region[1,:])))")
println("The maximum number of sim reps with treport - tinf = 0 is $(Int(maximum(tinf_treport_zero_count_by_region[1,:])))")
println("Percentage of sim reps with treport - tinf = 0 in each region")
for col in names(tinf_treport_zero_perc_by_region)
        println(rpad(string(col), 12), " | ", tinf_treport_zero_perc_by_region[1, col],"%")
end
println("Percentage of sim reps with treport - tinf = 0 for each region:")
println("  minimum = $(minimum(tinf_treport_zero_perc_by_region[1,:]))%")
println("  maximum = $(maximum(tinf_treport_zero_perc_by_region[1,:]))%")

# Compute % of simulation replicates that include an ICU case by region
CSV.write( normpath( joinpath( geo_optim_folder, "icu_td_regional/n_icu_cases_by_region_by_simrep_psampled_$(Int(psampled*100))_tsample_b4_tdeceased.csv") ), n_icu_by_region_by_simrep)

perc_simreps_icu_by_region = DataFrame( zeros(1,39), REGKEY.code)
for i in 1:ncol( n_icu_by_region_by_simrep ) # i=1
    perc_simreps_icu_by_region[1,i] = 100 * sum( n_icu_by_region_by_simrep[:,i] .> 0 ) / nrow(n_icu_by_region_by_simrep)
end
println("Number of sim reps with at least 1 ICU case by region")
for col in names(perc_simreps_icu_by_region)
    println(rpad(string(col), 12), " | ", perc_simreps_icu_by_region[1, col],"%")
end
println("Percentage of sim reps with at least 1 ICU case by region:")
println("  minimum = $(minimum(perc_simreps_icu_by_region[1,:]))%")
println("  maximum = $(maximum(perc_simreps_icu_by_region[1,:]))%")

# Check that column names match
names(perc_simreps_icu_by_region) == names(tinf_zero_perc_by_region)
# Plot 
Plots.scatter(    Vector(perc_simreps_icu_by_region[1,:])
                , Vector(tinf_zero_perc_by_region[1,:])
                , xlabel="% of sim replicates with ICU cases", ylabel = "% of TD with tinf = 0"
                , title = "By ITL2 region"
                , label = "ITL2 regions")

# Save plot
Plots.savefig( normpath( joinpath( geo_optim_folder,"icu_td_regional/perc_tinf_0_vs_perc_reps_w_ICU_samples_psampled_$(Int(psampled*100)).png") ) )
Plots.savefig( normpath( joinpath( root_dir,"NBPMscape/paper/results/plots/icu_td_regional/perc_tinf_0_vs_perc_reps_w_ICU_samples_psampled_$(Int(psampled*100)).png") ) )

# Compute median TD by region across all simulation replicates...
td_median_by_region = DataFrame( zeros(1,39), REGKEY.code)
for i in 1:ncol( td_median_by_region )
    td_median_by_region[1,i] = median( filter( isfinite, td_by_simrep_by_region[:,i] ) )
end
# ... and median TD by simulation replicate across all regions
td_median_by_simrep = Vector{Union{Missing, Float64}}(missing, length(sims)) #DataFrame( zeros(1000,1), :auto )
td_min_by_simrep    = Vector{Union{Missing, Float64}}(missing, length(sims)) #DataFrame( zeros(1000,1), :auto )
for j in 1:length( td_median_by_simrep ) # j=1
    #td_median_by_simrep[j,1] = median( filter( isfinite, td_by_simrep_by_region[j,:] ) )
    data_inf_filtered = filter( isfinite, [ getindex( td_by_simrep_by_region, j, k) for k in 1:size( td_by_simrep_by_region, 2)] ) 
    td_median_by_simrep[j] = isempty( data_inf_filtered ) ? missing : median( data_inf_filtered )
    td_min_by_simrep[j]    = isempty( data_inf_filtered ) ? missing : minimum( data_inf_filtered )
end
#TEST # StatsPlots.histogram(data_inf_filtered, breaks = 20)

# Plot to compare the median TD values by region and by simulation replicate
median_td_by_region = [getindex( td_median_by_region, 1, i) for i in 1:size( td_median_by_region, 2)] 
plt_h = Plots.histogram( median_td_by_region; #bins=30
              normalize=:pdf
             , label="median TDs by region (all sim reps)"
             , color=:blue
             , alpha=0.5
             , legend=:topleft
             , title="Comparison of median TDs by region and by simulation replicate" #"$(r) $(REGKEY[REGKEY.code .== r, :name][1])" 
             , xlabel="Days", ylabel = "Density"
             , size = (800, 400)
             , xlimit = [0,100]
             , ylimit = [0,0.3]
             , bins = 0:2:100
)
#median_td_by_simrep = [getindex( td_median_by_simrep, 1, i) for i in 1:size( td_median_by_simrep, 1)] 
plt_h = histogram!( td_median_by_simrep; #bins=50
              normalize=:pdf
             , label="median TDs by sim rep (all regions)"
             , color=:red, alpha=0.5#, legend=:topleft
             , bins = 0:2:100
             )
plt_h = vline!( [median(median_td_by_region)], color=:blue
            , label="median of median TDs by region (all sim reps)"
            , linewidth = 3) 
plt_h = vline!( [median(skipmissing(td_median_by_simrep))], color=:red
        , label="median of median TDs by sim rep (all regions)", linewidth = 3  )

# Add minimum TD values by region for each simulation replicate
plt_h = histogram!( td_min_by_simrep; #bins=50
              normalize=:pdf
             , label="minimum of regional TDs by sim rep"
             , color=:green, alpha=0.5 #, legend=:topleft
             , bins = 0:2:100
             )
plt_h = vline!( [median(skipmissing(td_min_by_simrep))], color=:green
        , label="median of minimum of regional TDs by sim rep", linewidth = 3  )
#vline!( [quantile( skipmissing(td_min_by_simrep) ,0.025)], color=:green
#        , label="95% CI: minimum of regional TDs by sim rep"
#        , linewidth = 3  , linestyle = :dash)
#vline!( [quantile( skipmissing(td_min_by_simrep) ,0.975)], color=:green
#        , linewidth = 3  , linestyle = :dash)

#anns = [(2, y[2], "Low"), (8, y[8], text("High", :blue, :center, 16))]
#anns = [(80,0.15,"Median of the median TD by region = ", round(median(median_td_by_region), digits=1))
#        ,(80,0.13,"Median of the minimum of regional TDs by sim rep = ", round( median(skipmissing(td_min_by_simrep)), digits=1))
#        ,(80,0.11,"with 95% CI: ",round(quantile( skipmissing(td_min_by_simrep) ,0.025), digits = 1)
#        ," to ", round(quantile( skipmissing(td_min_by_simrep) ,0.975), digits=1))
#]
#annotate!(plt_h, anns)
display(plt_h)

println("Median of the median TD by region = ", round(median(median_td_by_region), digits=1))
println("Median of the minimum of regional TDs by sim rep = ", round( median(skipmissing(td_min_by_simrep)), digits=1))
println("with 95% CI: ",round(quantile( skipmissing(td_min_by_simrep) ,0.025), digits = 1)
        ," to ", round(quantile( skipmissing(td_min_by_simrep) ,0.975), digits=1))

# Save plot
Plots.savefig( normpath( joinpath( geo_optim_folder,"icu_td_regional/median_TD_by_region_and_sim_rep_psampled_$(Int(psampled*100))_tsample_b4_tdeceased.png") ) )
Plots.savefig( normpath( joinpath( root_dir,"NBPMscape/paper/results/plots/icu_td_regional/median_TD_by_region_and_sim_rep_psampled_$(Int(psampled*100))_tsample_b4_tdeceased.png") ) )

##########################################################################################################
# (4) Sampling infections to simulate TDs using different sampling strategies:
#    (4.1) HCGS
#        (4.1.2) Compute correlations between TD distributions (1000 simulation replicates) for each region
##########################################################################################################

n = size(td_by_simrep_by_region, 2)
corrmat     = Matrix{Union{Missing, Float64}}(undef, n, n)
corrmat_idx = Matrix{Union{Missing, Float64}}(undef, n, n)

for i in 1:n # i=1; j=2
    for j in i:n
        # Extract columns as vectors
        x = td_by_simrep_by_region[!, i]
        y = td_by_simrep_by_region[!, j]
        # Find indices where both are not missing
        #idx = .!ismissing.(x) .& .!ismissing.(y)
        idx = isfinite.(x) .& isfinite.(y)
        #TEST
        #sum(isfinite.(x))
        #sum(isfinite.(y))
        #sum( isfinite.(x) .& isfinite.(y) )
        # Compute correlation on complete cases
        #corr = cor(skipmissing(x[idx]), skipmissing(y[idx]))
        x_filtered = skipmissing(x[idx])[:]
        y_filtered = skipmissing(y[idx])[:]
        corr = sum(idx) <= 1 ? missing : cor(x_filtered, y_filtered)
        corrmat[i, j] = corr
        corrmat[j, i] = corr  # Symmetric
        corrmat_idx[i, j] = sum(idx)
        corrmat_idx[j, i] = sum(idx) # Symmetric
    end
end
#println(corrmat_idx)
#mean(corrmat_idx) # 995.3 for psampled = 1

corr_df = DataFrame(corrmat, Symbol.(names(td_by_simrep_by_region))) # println(corr_df)
#rename!(corr_df, Symbol.(names(td_by_simrep_by_region)))

# Max correlation
maximum(skipmissing(corrmat)) # 1
max_corr = maximum(x for x in skipmissing(corrmat) if x != 1.0) # 0.4787323142501508 #0.47624984993358827
max_corr_ix = findall(isequal(max_corr), corrmat) # 39,15 and 15,39 # 31,4 and 4,31
max_cor_itl2_1 = names(td_by_simrep_by_region)[max_corr_ix[1][1]] # TLL5 #TLJ4
max_cor_itl2_2 = names(td_by_simrep_by_region)[max_corr_ix[1][2]] # TLG1 #TLD3

# Min correlation
min_corr = minimum(skipmissing(corrmat)) # 0.2953387121900001 #0.2782233212174087
min_corr_ix = findall(isequal(min_corr), corrmat) # 21,3 and 3,21 # 14,3 and 3,14
min_cor_itl2_1 = names(td_by_simrep_by_region)[min_corr_ix[1][1]] # TLH5 # TLF3
min_cor_itl2_2 = names(td_by_simrep_by_region)[min_corr_ix[1][2]] # TLD1 # TLD1

println("Regions with highest correlation between TDs:")
max_cor_regions = filter( row -> row.code in [max_cor_itl2_1, max_cor_itl2_2], REGKEY)
#Regions with highest correlation between TDs:
#2×2 DataFrame
# Row │ code    name                              
#     │ String  String                            
#─────┼───────────────────────────────────────────
#   1 │ TLG1    Herefordshire, Worcestershire an…
#   2 │ TLL5    South East Wales

println("Regions with lowest correlation between TDs:")
min_cor_regions = filter( row -> row.code in [min_cor_itl2_1, min_cor_itl2_2], REGKEY)
#Regions with lowest correlation between TDs:
#2×2 DataFrame
# Row │ code    name         
#     │ String  String
#─────┼──────────────────────
#   1 │ TLD1    Cumbria
#   2 │ TLH5    Norfolk

# CHECK for regions with lowest correlation
# Scatter plots for the maximum and minimum correlation region pairs
# Only includes TDs for simulations in which both regions had at least 1 ICU case 
idx = isfinite.(td_by_simrep_by_region[:,Symbol(min_cor_regions.code[1])]) .& isfinite.(td_by_simrep_by_region[:,Symbol(min_cor_regions.code[2])])
Plots.scatter( td_by_simrep_by_region[:,Symbol(min_cor_regions.code[1])][idx]
             , td_by_simrep_by_region[:,Symbol(min_cor_regions.code[2])][idx]
             , xlabel = "$(min_cor_regions.name[1]) ($(min_cor_regions.code[1]))\nTD (days)"
             , ylabel = "$(min_cor_regions.name[2]) ($(min_cor_regions.code[2]))\nTD (days)", label = ""
             , xlimit = [0,100], ylimit = [0,100]
             , color = :red)
Plots.savefig( normpath( joinpath( geo_optim_folder, "icu_td_regional/TD_min_correlation_region_pair_$(n_replicates)_replicates_psampled_$(Int(psampled*100)).png") ) )
Plots.savefig( normpath( joinpath( root_dir, "NBPMscape/paper/results/plots/icu_td_regional/TD_min_correlation_region_pair_$(n_replicates)_replicates_psampled_$(Int(psampled*100)).png") ) )

Plots.plot( td_by_simrep_by_region[:,Symbol(min_cor_regions.code[1])][idx], xlabel = "simulation replicate", ylabel = "TD (days)", label = "$(min_cor_regions.name[1]) ($(min_cor_regions.code[1]))" )
Plots.plot!( td_by_simrep_by_region[:,Symbol(min_cor_regions.code[2])][idx] , label = "$(min_cor_regions.name[2]) ($(min_cor_regions.code[2]))" )

# CHECK for regions with highest correlation
# Only includes TDs for simulations in which both regions had at least 1 ICU case 
idx = isfinite.(td_by_simrep_by_region[:,Symbol(max_cor_regions.code[1])]) .& isfinite.(td_by_simrep_by_region[:,Symbol(max_cor_regions.code[2])])
Plots.scatter( td_by_simrep_by_region[:,Symbol(max_cor_regions.code[1])][idx]
             , td_by_simrep_by_region[:,Symbol(max_cor_regions.code[2])][idx]
             , xlabel = "$(max_cor_regions.name[1]) ($(max_cor_regions.code[1]))\nTD (days)"
             , ylabel = "$(max_cor_regions.name[2]) ($(max_cor_regions.code[2]))\nTD (days)", label = ""
             , xlimit = [0,100], ylimit = [0,100]
             , color = :blue)
Plots.savefig( normpath( joinpath( geo_optim_folder, "icu_td_regional/TD_max_correlation_region_pair_$(n_replicates)_replicates_psampled_$(Int(psampled*100)).png") ) )
Plots.savefig( normpath( joinpath( root_dir, "NBPMscape/paper/results/plots/icu_td_regional/TD_max_correlation_region_pair_$(n_replicates)_replicates_psampled_$(Int(psampled*100)).png") ) )

# Plot of TD against sim rep number
Plots.plot( td_by_simrep_by_region[:,Symbol(max_cor_regions.code[1])][idx], xlabel = "simulation replicate", ylabel = "TD (days)", label = "$(max_cor_regions.name[1]) ($(max_cor_regions.code[1]))" )
Plots.plot!( td_by_simrep_by_region[:,Symbol(max_cor_regions.code[2])][idx] , label = "$(max_cor_regions.name[2]) ($(max_cor_regions.code[2]))" )

# Plot correlations
# Create region codes/names for labelling heatmap
x_names = names(td_by_simrep_by_region)
ITL2_dict = Dict(zip(REGKEY.code, REGKEY.name))
ITL2_names = [ITL2_dict[k] for k in x_names]
# Merge pairs into strings (e.g., "key:value" = "code:names")
y_names = [string(k, ":", v) for (k, v) in zip(x_names, ITL2_names)]
#println(y_names)

# Create correlation values for annotating the heatmap (if required)
anns = []
for i in 1:n
    for j in 1:n
        val = round(corrmat[i, j], digits=2)
        # Note: y-axis is reversed in Plots.jl heatmap to match matrix layout
        push!(anns, (j, n - i + 1, string(val)))
    end
end
# Plot heatmap
Plots.heatmap(  x_names # x-tick labels (columns)
                , xmirror = true # move x-axis labels to the top
                , xrotation = 90 # Rotate x-axis labels
                #, xticks = (1:n, names(td_by_simrep_by_region))
                , xticks = (collect(1:n) .- 0.5, x_names) #names(td_by_simrep_by_region)) # centers of each cell
                , y_names #names(td_by_simrep_by_region) # y-tick labels (rows)
                #, yticks = (1:n, names(td_by_simrep_by_region))
                , yticks = (collect(1:n) .- 0.5, y_names) #names(td_by_simrep_by_region))# centers of each cell
                , corrmat;              # the correlation matrix
                color = :RdBu        # diverging color map, blue to red
                , clim = (-1, 1)       # color limits for correlation coefficients
                , xlabel = "ITL2 region"
                , ylabel = "ITL2 region"
                , title = "Correlation Matrix for regional TD across $(n_replicates) simulation replicates"
                , yflip = true         # flips y-axis to match matrix layout
                , size = (1100, 800)
                #, annotations = anns
)
# Save heatmap
Plots.savefig( normpath( joinpath( geo_optim_folder, "icu_td_regional/median_TD_corr_simreps_$(n_replicates)_replicates_psampled_$(Int(psampled*100)).png") ) )
Plots.savefig( normpath( joinpath( root_dir, "NBPMscape/paper/results/plots/icu_td_regional/median_TD_corr_simreps_$(n_replicates)_replicates_psampled_$(Int(psampled*100)).png") ) )

# heatmap without diagonals
corrmat_wo_diag = corrmat
for i in 1:size(corrmat_wo_diag, 1)
    corrmat_wo_diag[i,i] = NaN
end
Plots.heatmap(  x_names # x-tick labels (columns)
                , xmirror = true # move x-axis labels to the top
                , xrotation = 90 # Rotate x-axis labels
                #, xticks = (1:n, names(td_by_simrep_by_region))
                , xticks = (collect(1:n) .- 0.5, x_names) #names(td_by_simrep_by_region)) # centers of each cell
                , y_names #names(td_by_simrep_by_region) # y-tick labels (rows)
                #, yticks = (1:n, names(td_by_simrep_by_region))
                , yticks = (collect(1:n) .- 0.5, y_names) #names(td_by_simrep_by_region))# centers of each cell
                , corrmat_wo_diag;              # the correlation matrix
                color = :RdBu        # diverging color map, blue to red
                #, clim = (-1, 1)       # color limits for correlation coefficients
                , xlabel = "ITL2 region"
                , ylabel = "ITL2 region"
                , title = "Correlation Matrix for regional TD across $(n_replicates) simulation replicates"
                , yflip = true         # flips y-axis to match matrix layout
                , size = (1100, 800)
                #, annotations = anns
)
# Save heatmap
Plots.savefig( normpath( joinpath( geo_optim_folder, "icu_td_regional/median_TD_corr_simreps_$(n_replicates)_replicates_psampled_$(Int(psampled*100))_no_diag.png") ) )
Plots.savefig( normpath( joinpath( root_dir, "NBPMscape/paper/results/plots/icu_td_regional/median_TD_corr_simreps_$(n_replicates)_replicates_psampled_$(Int(psampled*100))_no_diag.png") ) )

# Save correlation matrix to a file
corr_df_w_names = insertcols!(corr_df, 1, :ITL2_region => y_names )
println(corr_df_w_names)
CSV.write( normpath( joinpath( geo_optim_folder, "icu_td_regional/regional_TD_correlation_matrix_simreps_$(n_replicates)_psampled_$(Int(psampled*100)).csv") ), corr_df_w_names)

## Plot a correlation matrix with a reduced number of ITL2 regions
#region_names = popfirst!( names(corr_df) ) 
region_names = names(corr_df)
subset_regions = ["TLI4","TLI3","TLI6","TLH4","TLG3","TLF2","TLF1","TLC2","TLE3","TLE4","TLD3","TLD7","TLJ3","TLJ1","TLK1","TLL2","TLM1","TLN0"]
subset_regions = sort(subset_regions)
#println(subset_regions)
#["TLC2", "TLD3", "TLD7", "TLE3", "TLE4", "TLF1", "TLF2", "TLG3", "TLH4"
#, "TLI3", "TLI4", "TLI6", "TLJ1", "TLJ3", "TLK1", "TLL2", "TLM1", "TLN0"]
# Subset regions not currently in region names - due to change of code/region
missing_regions = subset_regions[ findall(i -> isnothing(i), indexin(subset_regions, region_names)) ]
# TLC2 changed in 2025 - replace with TLC4
# TLK1 changed in 2025 - replace with TLK5
# TLL2 changed in 2025 - replace with TLL5
subset_regions = [["TLC4","TLK5","TLL5"];subset_regions]
subset_regions = sort(subset_regions) #println(subset_regions)
# Find indices of the subset regions in full list of region names
subset_idx = findall(x -> x in subset_regions, region_names)
#println( region_names )
#println(subset_regions)
#println( region_names[subset_idx] )
#println(subset_idx)
# Subset the correlation matrix
filtered_corr = corr_df_w_names[subset_idx.-1, [1;subset_idx]]
# Plot heatmap
subset_x_names = names(filtered_corr)[2:end]
subset_y_names = filtered_corr[:,1]
Plots.heatmap(  subset_x_names # x-tick labels (columns)
        , xmirror = true # move x-axis labels to the top
        , xrotation = 90 # Rotate x-axis labels
        , xticks = (collect(1:n) .- 0.5, subset_x_names) # centers of each cell
        , subset_y_names # y-tick labels (rows)
        , yticks = (collect(1:n) .- 0.5, subset_y_names) #names(td_by_simrep_by_region))# centers of each cell
        , Matrix( filtered_corr[:,2:end] ) ;#corrmat;              # the correlation matrix
          color = :RdBu        # diverging color map, blue to red
        , clim = (-1, 1)       # color limits for correlation coefficients
        , xlabel = "ITL2 region"
        , ylabel = "ITL2 region"
        , title = "Correlation Matrix for subset of regional TD across $(n_replicates) simulation replicates"
        , yflip = true         # flips y-axis to match matrix layout
        , size = (1100, 800)
        #, annotations = anns
)

# Save heatmap plot and data
Plots.savefig( normpath( joinpath( geo_optim_folder,"icu_td_regional/regional_TD_corr_simreps_$(n_replicates)_psampled_$(Int(psampled*100))_reduced.png") ) )
Plots.savefig( normpath( joinpath( root_dir,"NBPMscape/paper/results/plots/icu_td_regional/regional_TD_corr_simreps_$(n_replicates)_psampled_$(Int(psampled*100))_reduced.png") ) )
CSV.write( normpath( joinpath( geo_optim_folder, "icu_td_regional/regional_TD_correlation_matrix_simreps_$(n_replicates)_psampled_$(Int(psampled*100))_reduced.csv") ), filtered_corr)

#println(countmap(corr_df_w_names[:,2:size(corr_df_w_names,2)])) 
#sort( countmap( corrmat ) )

# Remove diagonal values from subset correlation matrix
corrmat_filtered_wo_diag = copy(filtered_corr)
for i in 1:size(corrmat_filtered_wo_diag, 1)
    corrmat_filtered_wo_diag[i,i+1] = NaN
end

# Plot heatmap
subset_x_names = names(corrmat_filtered_wo_diag)[2:end]
subset_y_names = corrmat_filtered_wo_diag[:,1]
Plots.heatmap(  subset_x_names # x-tick labels (columns)
        , xmirror = true # move x-axis labels to the top
        , xrotation = 90 # Rotate x-axis labels
        , xticks = (collect(1:n) .- 0.5, subset_x_names) # centers of each cell
        , subset_y_names # y-tick labels (rows)
        , yticks = (collect(1:n) .- 0.5, subset_y_names) #names(td_by_simrep_by_region))# centers of each cell
        , Matrix( corrmat_filtered_wo_diag[:,2:end] ) ;#corrmat;              # the correlation matrix
          color = :RdBu        # diverging color map, blue to red
        , clim = (-1, 1)       # color limits for correlation coefficients
        , xlabel = "ITL2 region"
        , ylabel = "ITL2 region"
        , title = "Correlation Matrix for subset of regional TD across $(n_replicates) simulation replicates"
        , yflip = true         # flips y-axis to match matrix layout
        , size = (1100, 800)
        #, annotations = anns
        , tickfontsize = 10
)

# Save heatmap plot and data
Plots.savefig( normpath( joinpath( geo_optim_folder, "icu_td_regional/regional_TD_corr_simreps_$(n_replicates)_psampled_$(Int(psampled*100))_reduced_no_diag.png") ) )
Plots.savefig( normpath( joinpath( root_dir, "NBPMscape/paper/results/plots/icu_td_regional/regional_TD_corr_simreps_$(n_replicates)_psampled_$(Int(psampled*100))_reduced_no_diag.png") ) )
CSV.write( normpath( joinpath(geo_optim_folder, "icu_td_regional/regional_TD_correlation_matrix_simreps_$(n_replicates)_psampled_$(Int(psampled*100))_reduced_no_diag.csv") ), corrmat_filtered_wo_diag )     

# Plot heatmap - repeat but with limits equal to min and max correlation rather than -1 to +1
Plots.heatmap(  subset_x_names # x-tick labels (columns)
        , xmirror = true # move x-axis labels to the top
        , xrotation = 90 # Rotate x-axis labels
        , xticks = (collect(1:n) .- 0.5, subset_x_names) # centers of each cell
        , subset_y_names # y-tick labels (rows)
        , yticks = (collect(1:n) .- 0.5, subset_y_names) #names(td_by_simrep_by_region))# centers of each cell
        , Matrix( corrmat_filtered_wo_diag[:,2:end] ) ;#corrmat;              # the correlation matrix
          color = :RdBu        # diverging color map, blue to red
        #, clim = (-1, 1)       # color limits for correlation coefficients
        , xlabel = "ITL2 region"
        , ylabel = "ITL2 region"
        , title = "Correlation Matrix for subset of regional TD across $(n_replicates) simulation replicates"
        , yflip = true         # flips y-axis to match matrix layout
        , size = (1100, 800)
        #, annotations = anns
        , tickfontsize = 10
)

# Save heatmap plot and data
Plots.savefig( normpath( joinpath( geo_optim_folder, "icu_td_regional/regional_TD_corr_simreps_$(n_replicates)_psampled_$(Int(psampled*100))_reduced_no_diag_short_range.png") ) )
Plots.savefig( normpath( joinpath( root_dir, "NBPMscape/paper/results/plots/icu_td_regional/regional_TD_corr_simreps_$(n_replicates)_psampled_$(Int(psampled*100))_reduced_no_diag_short_range.png") ) )

# List region pairs with a correlation above a specified threshold
#corr_matrix = cor(Matrix(df))
varnames = names(td_by_simrep_by_region)
n = length(varnames)
if psampled == 0.1
    threshold = 0.75
elseif psampled == 1.0
    threshold = 0.25
else 
    threshold = 0.25
end

results = DataFrame(Row = String[], Col = String[], Correlation = Float64[])

for i in 1:n-1
    for j in i+1:n
        val = corrmat[i, j]
        if ismissing(val)
            continue
        else
            if abs(val) > threshold
                push!(results, (varnames[i], varnames[j], val))
            end
        end
    end
end
# 39*39
sort!(results, :Correlation, rev=true)#, by=abs)
println(results)

# Build a Dict for fast lookup
region_lookup = Dict(REGKEY.code .=> REGKEY.name)

# Prepare results DataFrame with extra columns
results = DataFrame(Region_1 = String[], Region_2 = String[]
                    , Region_1_name = String[], Region_2_name = String[]
                    , Correlation = Float64[]
                    , sim_replicates_n = Int[]
                    #, Stat_Power = Float64[], required_n = Int[] 
                    )

for i in 1:n-1
    for j in i+1:n
        val = round( corrmat[i, j], digits = 2)
        if ismissing(val)
            continue
        else
            if abs(val) > threshold
                row_key = varnames[i]
                col_key = varnames[j]
                row_val = get(region_lookup, row_key, missing)
                col_val = get(region_lookup, col_key, missing)
                sim_replicates_n = Int( ceil( corrmat_idx[i,j] ) )
                push!(results, (row_key, col_key, row_val, col_val, val
                                , sim_replicates_n
                                #, stat_power, , req_n 
                                )
                     )
            end
        end
    end
end

# Sort results by absolute correlation, high to low
sort!(results, :Correlation, rev=true)#, by=abs)
println(results)

#TEST
# Frequency table of homeregions
#println(countmap(G.homeregion))


### For each region, fit gamma distribution to infection times (tinf) corresponding to the minimum report time (treport), aka TD

# Simulation maxtime
maxtime = 70.0

# Initialise df to record regions and the median tinf and treport (TD) times
median_times_by_region = DataFrame( 
      ITL2_code          = String[]
    , ITL2_name          = String[]
    , tinf_median        = Float64[]
    , tinf_median_fit    = Float64[]
    , treport_median     = Float64[]
    , treport_median_fit = Float64[]
)

# Loop through all ITL2 regions
for r in REGKEY.code # TEST: r = "TLI7"; r = icu_regions[2] ; r=REGKEY.code[1]
    
    # Times to fit truncated gamma distribution to
    times = tinf_by_region_dict[ r ]
    #Test
    #histogram(times)
    CSV.write( joinpath( root_dir, "NBPMscape/paper/analysis/2_sampling_analysis/covid_like/geo_optimisation/icu_td_regional/median_TD_by_region_plots/sim_regionentry/tinf_times_$(r).csv"), DataFrame(value = times))
        
    # Define truncation range
    lower = 0 #0.00001
    upper = maxtime

    # Initial guesses
    #init_gamma = [10.0, 5.0] # shape*scale = mean. Median is ~50 for the tinf data # [1.0, 1.0]
    init_gamma = [mean(times)^2 / var(times), var(times)/mean(times)] # shape*scale = median. Mean is ~50 for the tinf data # [1.0, 1.0]
    #init_gamma_mixed = [init_gamma; 0.05] # For use in a gamma mixed model with probability of zero
    init_weibull = init_gamma #[20.0, 80.0]
    init_nbinom = [20.0, 0.5]
    init_lognorm = [mean(log.(times.+1)), std(log.(times.+1))]
    init_normal = [mean(times), std(times)]

    fit_tinf_results = fit_multi_dist(; times
                                    , init_gamma, init_weibull #, init_nbinom, init_lognorm, , init_normal
                                    , lower, upper)
    
    # Convert to a df and print
    fit_tinf_results_df_rows = [
                                (
                                  dist
                                , v[1] 
                                , v[2]
                                , v[3][1]
                                , v[3][2]
                                , missing  # spare column, can be changed as needed
                                ) for (dist, v) in fit_tinf_results[1]
                                ]
    fit_tinf_results_df = DataFrame(fit_tinf_results_df_rows
                                    , [:Distribution, :Negative_log_lik, :AIC, :param1, :param2, :param3])
    #println(fit_tinf_results_df)

    ## Compute fit data for plotting
    # Gamma
    gamma_shape_fit, gamma_scale_fit = fit_tinf_results[1]["Gamma"][3]
    tinf_fitted_gamma = Gamma(gamma_shape_fit, gamma_scale_fit)
    tinf_fitted_gamma_trunc = truncated(tinf_fitted_gamma, 0.0, upper)
    # Overlay fitted (truncated) gamma PDF
    xs = range(0, upper; length=300)
    norm_gamma = cdf(tinf_fitted_gamma, upper) # normalization for truncation
    tinf_pdf_vals_gamma = pdf.(tinf_fitted_gamma, xs) ./ norm_gamma # truncated PDF
    tinf_pdf_vals_gamma_trunc = pdf.(tinf_fitted_gamma_trunc, xs) #./ norm_gamma # truncated PDF
    #CHECK
    #Plots.plot(tinf_pdf_vals_gamma)
    #Plots.plot!(tinf_pdf_vals_gamma_trunc, col="red")
    
    # Weibull
    weibull_shape_fit, weibull_scale_fit = fit_tinf_results[1]["Weibull"][3]
    tinf_fitted_weibull = Gamma(weibull_shape_fit, weibull_scale_fit)
    tinf_fitted_weibull_trunc = truncated(tinf_fitted_weibull, 0.0, upper)
    # Overlay fitted (truncated) gamma PDF
    norm_weibull = cdf(tinf_fitted_weibull, upper) # normalization for truncation
    tinf_pdf_vals_weibull = pdf.(tinf_fitted_weibull, xs) ./ norm_weibull # truncated PDF
    tinf_pdf_vals_weibull_trunc = pdf.(tinf_fitted_weibull_trunc, xs) #./ norm_weibull # truncated PDF
    #CHECK
    #Plots.plot(tinf_pdf_vals_weibull)
    #Plots.plot!(tinf_pdf_vals_weibull_trunc, col="red")
    
    # Plot histogram (normalized to probability density)
    p1 = Plots.histogram(times; bins=30
                    , normalize=:pdf
                    , label = "Time of infection"#label="tinf"
                    , color=:lightblue, alpha=0.5, legend=:topright
                    , title="$(r) $(REGKEY[REGKEY.code .== r, :name][1])" 
                    , xlabel="Days", ylabel = "Density"
                    , xlims = [0,200]
                    , ylims = [0,0.15]
                   )

    # Plot tinf fit
    p1 = Plots.plot!(xs, tinf_pdf_vals_gamma; label="Gamma fit", color=:blue, lw=2)
    p1 = Plots.plot!(xs, tinf_pdf_vals_gamma; label="tinf truncated Gamma fit", color=:blue, lw=2)
    p1 = Plots.plot!(xs, tinf_pdf_vals_weibull; label="tinf truncated Weibull fit", color=:purple, lw=2)

    # Add median for tinf
    tinf_median = quantile(times, 0.5)
    times_no_zeros = sort( filter(x -> x != 0, times) )
    tinf_median_no_zeros = quantile(times_no_zeros, 0.5)
    #tinf_median = median(times)
    p1 = vline!([tinf_median], color=:blue
            , label="tinf median")
    p1 = vline!([tinf_median_no_zeros], color=:purple
            , label="tinf (no zeros) median")

    # Add median for fitted tinf
    tinf_median_fit_gamma = quantile(tinf_fitted_gamma_trunc, 0.5)
    p1 = vline!([tinf_median_fit_gamma], color=:blue, linestyle=:dash
            , label="tinf median truncated Gamma fit")
    tinf_median_fit_weibull = quantile(tinf_fitted_weibull_trunc, 0.5)
    p1 = vline!([tinf_median_fit_weibull], color=:purple, linestyle=:dash
            , label="tinf median truncated Weibull fit")

    # Add histogram for treport
    p1 = Plots.histogram!(treport_by_region_dict[ r ]; bins=30
                , normalize=true
                , label="treport", color=:pink, alpha=0.5)

    # Add median for treport
    treport_median = quantile(treport_by_region_dict[ r ], 0.5)
    p1 = vline!([treport_median], color=:red, label="treport median")
    
    # Fit gamma distribution to (treport - tinf) from earliest treport
    times = tinf_treport_by_region_dict[ r ]
    
    # Initial guesses
    #init_params = [1.0, 1.0]
    init_gamma = [mean(times)^2 / var(times), var(times)/mean(times)] # shape*scale = median. Mean is ~50 for the tinf data # [1.0, 1.0]
    init_weibull = init_gamma
    init_nbinom = [20.0, 0.5]
    init_lognorm = [mean(log.(times.+1)), std(log.(times.+1))]
    init_normal = [ mean(times), std(times)]
    
    fit_tinf_treport_results = fit_multi_dist(; times
                                               , init_gamma, init_weibull #, init_normal, init_nbinom, init_lognorm
                                               , lower, upper)
    #println(fit_tinf_treport_results)

    # Gamma fit
    tinf_treport_gamma_shape_fit, tinf_treport_gamma_scale_fit = fit_tinf_treport_results[1]["Gamma"][3]
    tinf_treport_fitted_gamma = Gamma(tinf_treport_gamma_shape_fit, tinf_treport_gamma_scale_fit)
    tinf_treport_median_gamma_fit = median( tinf_treport_fitted_gamma ) # Note that this distribution is not truncated
    
    # Add median for fitted truncated gamma tinf + fitted gamma (treport-tinf) = fitted treport
    p1 = vline!([tinf_treport_median_gamma_fit + tinf_median_fit_gamma], color=:red, linestyle=:dash
            , label="treport median Gamma fit")

    # Add histogram for (treport - tinf)
    p2 = Plots.histogram(tinf_treport_by_region_dict[ r ]; bins=30, normalize=:pdf
                   , label="treport-tinf", color=:green, alpha=0.5, legend=:topright
                   , xlabel="Days", ylabel="Density"
                   , xlims = [0,100]
                   , ylims = [0,0.15]
                    )

    
    # Overlay fitted gamma PDF
    tinf_treport_pdf_vals_gamma = pdf.(tinf_treport_fitted_gamma, xs)
    # Plot (treport - tinf) fit
    p2 = Plots.plot!(xs, tinf_treport_pdf_vals_gamma; label="treport-tinf gamma fit"
            , color=:green, lw=2)

    # Add median for (treport - tinf)
    tinf_treport_median = quantile(tinf_treport_by_region_dict[ r ], 0.5) # Note does not include any zeros
    p2 = vline!([tinf_treport_median], color=:green
            , label="treport-tinf median", lw = 2)

    # Add median for fitted (treport - tinf)
    tinf_treport_median_fit_gamma = quantile(tinf_treport_fitted_gamma, 0.5)
    p2 = vline!([tinf_treport_median_fit_gamma], color=:green, linestyle=:dash
            , label="treport-tinf median gamma fit")

    Plots.plot( size=(800, 1200), p1, p2, layout = (2,1))
    
    Plots.savefig( joinpath( geo_optim_folder, "icu_td_regional/TD_fits_by_region_tsample_b4_tdeceased/$(r).png") )
    Plots.savefig( joinpath( root_dir, "NBPMscape/paper/results/plots/icu_td_regional/TD_fits_by_region_tsample_b4_tdeceased/$(r).png") )
    
    # Add row to df to record the median value by region
    new_row = ( ITL2_code = r
                , ITL2_name = REGKEY[REGKEY.code .== r, :name][1]
                , tinf_median = tinf_median
                , tinf_median_fit = tinf_median_fit_gamma
                , treport_median = treport_median
                , treport_median_fit = (tinf_treport_median_fit_gamma + tinf_median_fit_gamma)
                )
    push!( median_times_by_region, new_row )
     
end

println( median_times_by_region )
median_times_by_region_sorted = sort( median_times_by_region , :treport_median_fit )
println( median_times_by_region_sorted )

# Add the percent of simulation replicates that include ICU cases for each region
median_times_by_region_sorted[!,:perc_simreps_have_icu_cases] = zeros(nrow(median_times_by_region_sorted))
for i in 1:nrow( median_times_by_region_sorted )
    median_times_by_region_sorted.perc_simreps_have_icu_cases[i] = perc_simreps_icu_by_region[ 1 , median_times_by_region_sorted[i,"ITL2_code"] ]
end    
println(median_times_by_region_sorted)
minimum(median_times_by_region_sorted.perc_simreps_have_icu_cases )
maximum(median_times_by_region_sorted.perc_simreps_have_icu_cases )

# For each region find the top 3 regions with the highest correlations and add name and correlation value
# Add columns
for i in 1:3
    median_times_by_region_sorted[!, Symbol("high_cor_$(i)_region")] = fill("", nrow(median_times_by_region_sorted))
    median_times_by_region_sorted[!, Symbol("high_cor_$(i)_correlation")] = zeros(nrow(median_times_by_region_sorted))
end

# Fill columns
for i in 1:length(median_times_by_region_sorted.ITL2_code)
    ## Find regions with the top 3 highest correlations with region of interest
    # Region correlations with region i
    region_correls = DataFrame( ITL2_region = corr_df_w_names[:,1]
                                , correlation = corr_df_w_names[:,median_times_by_region_sorted.ITL2_code[i]]
                              )
    # Remove correlation with self
    region_i = median_times_by_region_sorted.ITL2_code[i]
    self_index = findfirst(j -> first(region_correls.ITL2_region[j], 4) == first(region_i, 4)
                            , 1:nrow(region_correls))
    # Remove self-index row
    if self_index !== nothing
        region_correls_ex_self = filter(row -> row !== self_index, eachindex(eachrow(region_correls))) .|> x -> region_correls[x, :]
        region_correls_ex_self = DataFrame(region_correls_ex_self)  # Reconstruct df2 from filtered rows
        #println(region_correls_ex_self)
    end
    #println(region_correls_ex_self.ITL2_region )
    # Filter out regions with missing correlation (i.e. no overlapping simulation replicates with region i)
    region_correls_filtered = region_correls_ex_self[ findall(!ismissing, region_correls_ex_self[:,2] ), : ]
    # Select top 3 correlations with region i
    top_3 = sort( region_correls_filtered , :correlation , rev=true)[1:3,1:2]
    # Fill data in median_times_by_region_sorted df
    for k in 1:3
        median_times_by_region_sorted[!, Symbol("high_cor_$(k)_region")][i]      = top_3[k,1]
        median_times_by_region_sorted[!, Symbol("high_cor_$(k)_correlation")][i] = top_3[k,2]
    end
       
end

println(median_times_by_region_sorted)

CSV.write( normpath( joinpath( geo_optim_folder, "icu_td_regional/TD_by_region_sorted_simreps_$(n_replicates)_psampled_$(Int(psampled*100)).csv" ) ), median_times_by_region_sorted )

###############################################################################################################
# (4) Sampling infections to simulate TDs using different sampling strategies:
#    (4.1) HCGS
#        (4.1.3) Perform hierarchical clustering to group regions together based on TD distribution correlations
#                and allocate a number of sites to each region. Regions with lower median TD are prioritised and
#                sites are distributed across both groups and regions.
###############################################################################################################

# Use '1b_site_selection_2026_10.R' to use block model and hierarchical clustering
# to allocation numbers of sites to ITL2 regions.
# Then continue analysis with '1c_geographic_optimisation_part_2.jl'
