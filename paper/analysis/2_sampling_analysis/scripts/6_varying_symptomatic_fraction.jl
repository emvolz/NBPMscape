#=

Script for the investigating impact on metagenomic surveillance sampling performance in COVID-like outbreaks with different
symptomatic fraction. Five different symptomatic proportions (which are age disaggregated) investigated:
- Base case assumption (from Knock et al (2021))
- Base case - 20 percentage points across all ages
- Base case + 20 percentage points across all ages
- Base case + 40 percentage points across all ages
- 100% symptomatic across all ages

This file contains:
(0) Define root directory
(1) Plot age disaggregated symptomatic profile for each pathogen outbreak simulated
(2) Base case - 20pp
    (2.1) Information on simulated outbreaks
    (2.2) Load simulated outbreak files
    (2.3) Sample infections in primary and tertiary settings
    (2.4) Compute median TDs
(3) Base case + 20pp
    (3.1) Information on simulated outbreaks
    (3.2) Load simulated outbreak files
    (3.3) Sample infections in primary and tertiary settings
    (3.4) Compute median TDs
(4) Base case + 40pp
    (4.1) Information on simulated outbreaks
    (4.2) Load simulated outbreak files
    (4.3) Sample infections in primary and tertiary settings
    (4.4) Compute median TDs
(5) 100% symptomatic
    (5.1) Information on simulated outbreaks
    (5.2) Load simulated outbreak files
    (5.3) Sample infections in primary and tertiary settings
    (5.4) Compute median TDs
(6) Combine individual simulation results into a single df and save .csv
    (6.1) Read in TD data dictionaries
    (6.2) Merge into single dataframe
(7) Plots and TD differences
(8) Statistical tests

=#

#-----------------------------------
# Activate environment for analysis
#-----------------------------------
# Unless otherwise stated, run on Laptop - 11th Gen Intel(R) Core(TM) i7-1185G7 @ 3.00GHz (1.80 GHz) with 32.0 GB RAM
# Using Julia v1.12.3 in VS Code 1.141.0

using Pkg
Pkg.activate(@__DIR__) # ...NBPMscape\\paper\\analysis\\2_sampling_analysis\\scripts

# Check packages in use
# Using NBPMscape version from commit "2762757a687c64146016e6bc742d258e20525cb8" as at 17:00 16 Oct 2026
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
using JLD2
using NBPMscape # Should be version from commit "2762757a687c64146016e6bc742d258e20525cb8" as at 17:00 16 Oct 2026
using Plots
using Statistics

############################################################################################
# (0) Define root directory
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

############################################################################################
# (1) Plot age disaggregated symptomatic profile for each pathogen outbreak simulated
############################################################################################
## Plots of symptomatic proportions by age used in different analyses

# Read in values
covid_symp_base       = CSV.read( "data/severity_care_probabilities/covid_knock_symptomatic_ihr_ifr_by_age.csv",              DataFrame )
covid_symp_minus_20pp = CSV.read( "data/severity_care_probabilities/covid_knock_symptomatic_lower_ihr_ifr_by_age.csv",        DataFrame )
covid_symp_plus_20pp  = CSV.read( "data/severity_care_probabilities/covid_knock_symptomatic_higher_ihr_ifr_by_age.csv",       DataFrame )
covid_symp_plus_40pp  = CSV.read( "data/severity_care_probabilities/covid_knock_symptomatic_base_plus_40_ihr_ifr_by_age.csv", DataFrame )
covid_symp_all100     = CSV.read( "data/severity_care_probabilities/covid_knock_symptomatic_100_ihr_ifr_by_age.csv",          DataFrame )
influenza_symp_mild   = CSV.read( "data/severity_care_probabilities/influenza_symptomatic_ihr_ifr_by_age_mild.csv",           DataFrame )
influenza_symp_severe = CSV.read( "data/severity_care_probabilities/influenza_symptomatic_ihr_ifr_by_age_severe.csv",         DataFrame )

Plots.plot(covid_symp_base.age_single_year,        covid_symp_base.symptomatic_prob,       label = "COVID-like: base symptomatic",        linewidth=2, color = :goldenrod1, xlabel = "Age of infectee", ylabel = "Symptomatic proportion", legend=true, ylim=[0,1], foreground_color_legend = nothing)
Plots.plot!(covid_symp_minus_20pp.age_single_year, covid_symp_minus_20pp.symptomatic_prob, label = "COVID-like: base symptomatic - 20pp", linewidth=2, color = :green)
Plots.plot!(covid_symp_plus_20pp.age_single_year,  covid_symp_plus_20pp.symptomatic_prob,  label = "COVID-like: base symptomatic + 20pp", linewidth=2, color = :blue)
Plots.plot!(covid_symp_plus_40pp.age_single_year,  covid_symp_plus_40pp.symptomatic_prob,  label = "COVID-like: base symptomatic + 40pp", linewidth=2, color = :darkorange1)
Plots.plot!(covid_symp_all100.age_single_year,     covid_symp_all100.symptomatic_prob,     label = "COVID-like: 100% symptomatic",        linewidth=3, color = :purple)
Plots.plot!(influenza_symp_mild.age_single_year,   influenza_symp_mild.symptomatic_prob,   label = "Influenza-like: mild",                linewidth=2, color = :red)
Plots.plot!(influenza_symp_severe.age_single_year, influenza_symp_severe.symptomatic_prob, label = "Influenza-like: severe",              linewidth=2, color = :black, linestyle=:dash)

Plots.savefig( joinpath( root_dir, "NBPMscape/paper/analysis/1_outbreak_simulations/covid_like_varied_sympt_prop/symp_prop.png") )
Plots.savefig( joinpath( root_dir, "NBPMscape/paper/results/figures_supplementary/raw/symp_prop.png") )

##########################################################################################################################
# (1) Define Primary care and Critical care sampling parameters that are common to all outbrake simulations in this script
##########################################################################################################################

# Set sampling parameters for different Primary care sampling scenarios
gp_sampling_params = (  scenario_01 = (scenario_number =  1, gp_ari_consults = 327, gp_ari_swabs =  747, gp_swabs_mg =  300)
                      , scenario_13 = (scenario_number = 13, gp_ari_consults = 327, gp_ari_swabs =  747, gp_swabs_mg =  500)
                      , scenario_25 = (scenario_number = 25, gp_ari_consults = 327, gp_ari_swabs = 1000, gp_swabs_mg = 1000)
                      #, scenario_37 = (scenario_number = 37, gp_ari_consults = 180, gp_ari_swabs =  319, gp_swabs_mg =  300)
                      #, scenario_49 = (scenario_number = 49, gp_ari_consults = 180, gp_ari_swabs =  500, gp_swabs_mg =  500)
                      #, scenario_61 = (scenario_number = 61, gp_ari_consults = 180, gp_ari_swabs = 1000, gp_swabs_mg = 1000)
                      )

# Set sampling parameters for different Critical Care mSCAPE sampling scenarios
mscape_icu_sampling_sites_8  = CSV.read( joinpath( root_dir, "NBPMscape/data/nhs_trust_site_sample_targets_24pilot.csv" ), DataFrame  )
mscape_icu_sampling_sites_26 = CSV.read( joinpath( root_dir, "NBPMscape/data/nhs_trust_site_sample_targets_24pilot_25cohort_longlist.csv" ), DataFrame  )

icu_sampling_params = ( # 8 sites
                          scenario_01 = (scenario_number =  1, n_icu_samples_per_week =  300, icu_ari_admissions =  1440, icu_nhs_trust_sampling_sites_file = mscape_icu_sampling_sites_8)
                        , scenario_13 = (scenario_number = 13, n_icu_samples_per_week =  500, icu_ari_admissions =  1440, icu_nhs_trust_sampling_sites_file = mscape_icu_sampling_sites_8)
                        , scenario_25 = (scenario_number = 25, n_icu_samples_per_week = 1000, icu_ari_admissions =  1440, icu_nhs_trust_sampling_sites_file = mscape_icu_sampling_sites_8)
                        #, scenario_37 = (scenario_number = 37, n_icu_samples_per_week =  300, icu_ari_admissions =   793, icu_nhs_trust_sampling_sites_file = mscape_icu_sampling_sites_8)
                        #, scenario_49 = (scenario_number = 49, n_icu_samples_per_week =  500, icu_ari_admissions =   793, icu_nhs_trust_sampling_sites_file = mscape_icu_sampling_sites_8)
                        #, scenario_61 = (scenario_number = 61, n_icu_samples_per_week = 1000, icu_ari_admissions =   793, icu_nhs_trust_sampling_sites_file = mscape_icu_sampling_sites_8)
                        # 26 sites
                        , scenario_73 = (scenario_number = 73, n_icu_samples_per_week =  300, icu_ari_admissions =  1440, icu_nhs_trust_sampling_sites_file = mscape_icu_sampling_sites_26)
                        , scenario_74 = (scenario_number = 74, n_icu_samples_per_week =  500, icu_ari_admissions =  1440, icu_nhs_trust_sampling_sites_file = mscape_icu_sampling_sites_26)
                        , scenario_75 = (scenario_number = 75, n_icu_samples_per_week = 1000, icu_ari_admissions =  1440, icu_nhs_trust_sampling_sites_file = mscape_icu_sampling_sites_26)
                        #, scenario_76 = (scenario_number = 76, n_icu_samples_per_week =  300, icu_ari_admissions =   793, icu_nhs_trust_sampling_sites_file = mscape_icu_sampling_sites_26)
                        #, scenario_77 = (scenario_number = 77, n_icu_samples_per_week =  500, icu_ari_admissions =   793, icu_nhs_trust_sampling_sites_file = mscape_icu_sampling_sites_26)
                        #, scenario_78 = (scenario_number = 78, n_icu_samples_per_week = 1000, icu_ari_admissions =   793, icu_nhs_trust_sampling_sites_file = mscape_icu_sampling_sites_26)
                    )


##################################################################################################################################
# (2) Base case - 20pp
#     (2.1) Information on simulated outbreaks
##################################################################################################################################
#=
    Simulation run on Imperial College High Performance Computing (HPC) system
    using a Julia script and PBS file to manage the process:

    Scripts: Julia - 'covidlike-bm20-1.6.1_w_G_filtered_nrep1.jl'
             PBS   - 'covidlike-bm20-1.6.1_w_G_filtered_array100_nrep1.pbs'
            (see NBPMscape\paper\analysis\1_outbreak_simulations\covid_like_varied_sympt_prop\base_minus_20pp\scripts)
    HPC run numbers: 4301080
    Simulation 'maxtime' = 70 days
    Simulation replicates per PBS array run: 10
    Runs in PBS array: 100
    Total sim reps completed: 1000
    CPUs per array run: 1
    RAM per array run: 64GB
    Longest batch run time: ~0.9hrs (~5 minutes per sim rep)
=#
##################################################################################################################################
# (2) Base case - 20pp
#     (2.2) Load simulated outbreak files
##################################################################################################################################

## Load simulated outbreak data with LOWER symptomatic proportion
main_folder = joinpath(root_dir,"NBPMscape/paper/analysis/1_outbreak_simulations/covid_like_varied_sympt_prop")
sims_bm20 = combine_sim_reps(;  # sim_input_folders in vector format where multiple paths
                          sim_input_folders = [ "$(main_folder)/base_minus_20pp/output/4301080_sims_G_filtered" ]
                        , sim_object_name = "sims_G_filtered" #"sims" # "sims_G_icu_filter" # "sims_G_gp_filter" # "sim_G_filtered"
                        , nrep = 10
                        )
# Check df that shows errors in the combination process
sims_lower_mismatch_df = sims_bm20[2]
# df is empty

# Refine sims
sims_bm20 = sims_bm20[1]

# Save to file
@save joinpath( main_folder
                , "base_minus_20pp/output"
                , "covidlike-1.6.1-sims_G_filtered_nrep1000_4301080_sims_bm20.jld2") sims_bm20

# Reload file if necessary
#sims_bm20 = load( joinpath( main_folder
#                            , "base_minus_20pp/output"
#                            , "covidlike-1.6.1-sims_G_filtered_nrep1000_4301080_sims_bm20.jld2")
#                    , "sims_bm20")


##################################################################################################################################
# (2) Base case - 20pp
#     (2.3) Sample infections in primary and tertiary settings
##################################################################################################################################

# Define path to output folder for analyses files
output_folder_bm20 = joinpath( root_dir, "NBPMscape/paper/analysis/2_sampling_analysis/covid_like_varied_sympt_prop/base_minus_20pp")

## Compute times to detection
# using different model parameters representing different scenarios
# See 'NBPMscape/config/HARISS/config_file_record.xlsx' for summary descriptions of scenarios

# Start with base parameters
initialize_parameters("config/outbreak_params_covid19_like_symptomatic_base_minus_20_20261008.yaml")


# Run analysis for 6 different scenarios (2 seasons x 3 numbers of samples) using GPs (Primary care)
for i in 1:length(gp_sampling_params) # in [1,13,25,37,49,61] # i=1
    gp_tds = gp_td(  sims = sims_bm20
                   , gp_ari_consults = gp_sampling_params[i].gp_ari_consults
                   , gp_ari_swabs    = gp_sampling_params[i].gp_ari_swabs
                   , gp_swabs_mg     = gp_sampling_params[i].gp_swabs_mg
                   )
    CSV.write( joinpath( root_dir, "NBPMscape/paper/analysis/2_sampling_analysis/covid_like_varied_sympt_prop/base_minus_20pp/gp_td"
                        ,"gp_tds_$(gp_sampling_params[i].scenario_number).csv" )
                , gp_tds )
    # Also, save file with 'GP_simid' column removed to reduce file size
    gp_tds_small = gp_tds[:,["sim_n","GP_TD","GP_3TD","n_GP_cases","n_GP_cases_sampled","n_GP_cases_sampled_positive","tinf_relating_to_GP_TD","last_tinf_relating_to_GP_3TD"]]
    CSV.write( joinpath( root_dir, "NBPMscape/paper/analysis/2_sampling_analysis/covid_like_varied_sympt_prop/base_minus_20pp/gp_td_small"
                        ,"small_gp_tds_$(gp_sampling_params[i].scenario_number).csv" )
                , gp_tds_small )
end


# Double check loaded parameters required for ICU sampling
NBPMscape.P.pathogen_type # "virus"
NBPMscape.P.icu_sample_type #"regional" 
NBPMscape.P.sample_icu_cases_version #"number"
NBPMscape.P.icu_only_sample_before_death #true            

# Run analysis for ICU sampling scenarios (2 seasons x 3 numbers of samples x 2 numbers of sampling sites)
for j in 1:length(icu_sampling_params) # [1,13,25,37,49,61,73,74,75,76,77,78] #j=1
    #initialize_parameters("config/influenza_v3/mild/outbreak_params_influenza_like_mild_$(j).yaml")
    
    icu_tds = icu_td( p = NBPMscape.P
                    , sims = sims_bm20
                    , n_icu_samples_per_week   = icu_sampling_params[j].n_icu_samples_per_week
                    , icu_ari_admissions       = icu_sampling_params[j].icu_ari_admissions
                    , nhs_trust_sampling_sites = icu_sampling_params[j].icu_nhs_trust_sampling_sites_file
                    )
    
    CSV.write( joinpath( root_dir, "NBPMscape/paper/analysis/2_sampling_analysis/covid_like_varied_sympt_prop/base_minus_20pp"
                        ,"icu_td/icu_tds_$(icu_sampling_params[j].scenario_number).csv" ), icu_tds )

    # Also, save file with 'ICU_simid' column removed to reduce file size
    icu_tds_small = icu_tds[:,["sim_n","ICU_TD","ICU_3TD","n_ICU_cases","n_ICU_cases_sampled","n_ICU_cases_sampled_positive","tinf_relating_to_ICU_TD","last_tinf_relating_to_ICU_3TD"]]
    CSV.write( joinpath( root_dir, "NBPMscape/paper/analysis/2_sampling_analysis/covid_like_varied_sympt_prop/base_minus_20pp"
                        ,"icu_td_small/small_icu_tds_$(icu_sampling_params[j].scenario_number).csv" ), icu_tds_small )
end



##################################################################################################################################
# (2) Base case - 20pp
#     (2.4) Compute median TDs
##################################################################################################################################

# Directory
path_name_bm20 = output_folder_bm20

### Compute median TDs for ICU and GP sampling
#scenario_numbers_2 = [1,13,25,37,49,61,73,74,75,76,77,78] # Note that scenarios 73-78 are with different ICU sites
scenario_numbers_2 = [1,13,25,73,74,75] # Note that scenarios 73-75 are with different ICU sites. 1,13,25 are with 8 ICU sites and 73-75 are with 26 sites.
# Prepare a container to keep the DataFrames
icu_tds_df_dict_bm20 = Dict{Int, DataFrame}()
 gp_tds_df_dict_bm20 = Dict{Int, DataFrame}()
# Prepare a df to hold the median TD values
icu_gp_tds_df_bm20 = DataFrame(
                    Scenario_number    = zeros(Int,     length(scenario_numbers_2))
                    ,icu_median_TD     = zeros(Float64, length(scenario_numbers_2))
                    ,icu_TD_iqr_lower  = zeros(Float64, length(scenario_numbers_2))
                    ,icu_TD_iqr_upper  = zeros(Float64, length(scenario_numbers_2))
                    ,icu_TD_std        = zeros(Float64, length(scenario_numbers_2))
                    ,icu_median_3TD    = zeros(Float64, length(scenario_numbers_2))
                    ,icu_3TD_iqr_lower = zeros(Float64, length(scenario_numbers_2))
                    ,icu_3TD_iqr_upper = zeros(Float64, length(scenario_numbers_2))
                    ,icu_3TD_std       = zeros(Float64, length(scenario_numbers_2))
                    ,gp_median_TD      = zeros(Float64, length(scenario_numbers_2))
                    ,gp_TD_iqr_lower   = zeros(Float64, length(scenario_numbers_2))
                    ,gp_TD_iqr_upper   = zeros(Float64, length(scenario_numbers_2))
                    ,gp_TD_std         = zeros(Float64, length(scenario_numbers_2))
                    ,gp_median_3TD     = zeros(Float64, length(scenario_numbers_2))
                    ,gp_3TD_iqr_lower  = zeros(Float64, length(scenario_numbers_2))
                    ,gp_3TD_iqr_upper  = zeros(Float64, length(scenario_numbers_2))
                    ,gp_3TD_std        = zeros(Float64, length(scenario_numbers_2))
                    )

# Loop through scenarios, read in files and compute median TD values
for i in 1:length(scenario_numbers_2) 
    k = scenario_numbers_2[i]
    file_icu = joinpath(path_name_bm20, "icu_td_small/small_icu_tds_$(string(k)).csv")
    
    # Read CSV into a DataFrame and store in the dictionary with the scenario number (numeric suffix) as key
    icu_tds_df_dict_bm20[k] = CSV.read(file_icu, DataFrame)
                            
    # Compute times to detection and populate df
    icu_gp_tds_df_bm20[i,"Scenario_number"] = k

    # ICU
    icu_gp_tds_df_bm20[i,"icu_median_TD"]  = median(icu_tds_df_dict_bm20[k][:,"ICU_TD"])
    icu_gp_tds_df_bm20[i,"icu_median_3TD"] = median(icu_tds_df_dict_bm20[k][:,"ICU_3TD"])
    
    icu_gp_tds_df_bm20[i,"icu_TD_iqr_lower"] = quantile(icu_tds_df_dict_bm20[k][:,"ICU_TD"], 0.25)
    icu_gp_tds_df_bm20[i,"icu_TD_iqr_upper"] = quantile(icu_tds_df_dict_bm20[k][:,"ICU_TD"], 0.75)
    
    icu_gp_tds_df_bm20[i,"icu_3TD_iqr_lower"] = quantile(icu_tds_df_dict_bm20[k][:,"ICU_3TD"], 0.25)
    icu_gp_tds_df_bm20[i,"icu_3TD_iqr_upper"] = quantile(icu_tds_df_dict_bm20[k][:,"ICU_3TD"], 0.75)
    
    icu_gp_tds_df_bm20[i,"icu_TD_std"] = std(icu_tds_df_dict_bm20[k][:,"ICU_TD"])
    icu_gp_tds_df_bm20[i,"icu_3TD_std"] = std(icu_tds_df_dict_bm20[k][:,"ICU_3TD"])

    # GP
    if k < 73 # Scenario numbers of 73 and above relate to the ICU sampling scenarios with 26 sites
        # Read CSV into a DataFrame and store in the dictionary with the scenario number (numeric suffix) as key
        file_gp  = joinpath(path_name_bm20,  "gp_td_small/small_gp_tds_$(string(k)).csv")
        gp_tds_df_dict_bm20[k]  = CSV.read(file_gp,  DataFrame)
        
        icu_gp_tds_df_bm20[i,"gp_median_TD"]   = median(gp_tds_df_dict_bm20[k][:,"GP_TD"])
        icu_gp_tds_df_bm20[i,"gp_median_3TD"]  = median(gp_tds_df_dict_bm20[k][:,"GP_3TD"])

        icu_gp_tds_df_bm20[i,"gp_TD_iqr_lower"] = quantile(gp_tds_df_dict_bm20[k][:,"GP_TD"], 0.25)
        icu_gp_tds_df_bm20[i,"gp_TD_iqr_upper"] = quantile(gp_tds_df_dict_bm20[k][:,"GP_TD"], 0.75)
        
        icu_gp_tds_df_bm20[i,"gp_3TD_iqr_lower"] = quantile(gp_tds_df_dict_bm20[k][:,"GP_3TD"], 0.25)
        icu_gp_tds_df_bm20[i,"gp_3TD_iqr_upper"] = quantile(gp_tds_df_dict_bm20[k][:,"GP_3TD"], 0.75)
        
        icu_gp_tds_df_bm20[i,"gp_TD_std"] = std(gp_tds_df_dict_bm20[k][:,"GP_TD"])
        icu_gp_tds_df_bm20[i,"gp_3TD_std"] = std(gp_tds_df_dict_bm20[k][:,"GP_3TD"])
    end

end

# View
println(icu_gp_tds_df_bm20)

println(icu_tds_df_dict_bm20[1][1:10,:])
println(gp_tds_df_dict_bm20[1][1:10,:])

# Save
CSV.write( joinpath( output_folder_bm20, "icu_gp_median_tds_bm20_symp.csv" ), icu_gp_tds_df_bm20 )
# Reload if necessary
#icu_gp_tds_df_bm20 = CSV.read( joinpath( output_folder_lower, "icu_gp_median_tds_bm20_symp.csv" ), DataFrame)

# Save dictionaries containing TD results for all sim reps
@save joinpath(output_folder_bm20, "covidlike-1.6.1-sims_G_filtered_nrep1000_4301080_gp_tds_df_dict_bm20.jld2") gp_tds_df_dict_bm20
@save joinpath(output_folder_bm20, "covidlike-1.6.1-sims_G_filtered_nrep1000_4301080_icu_tds_df_dict_bm20.jld2") icu_tds_df_dict_bm20
# Reload if necessary
#gp_tds_df_dict_lower  = load( joinpath(output_folder_bm20, "covidlike-1.6.1-sims_G_filtered_nrep1000_4301080_gp_tds_df_dict_bm20.jld2"),  "gp_tds_df_dict_bm20")
#icu_tds_df_dict_lower = load( joinpath(output_folder_bm20, "covidlike-1.6.1-sims_G_filtered_nrep1000_4301080_icu_tds_df_dict_bm20.jld2"), "icu_tds_df_dict_bm20")

# Reduce memory usage
sims_bm20 = nothing

##################################################################################################################################
# (3) Base case + 20pp
#     (3.1) Information on simulated outbreaks
##################################################################################################################################
#=
    Simulation run on Imperial College High Performance Computing (HPC) system
    using a Julia script and PBS file to manage the process:

    Scripts: Julia - 'covidlike-bp20-1.6.1_w_G_filtered_nrep1.jl'
             PBS   - 'covidlike-bp20-1.6.1_w_G_filtered_array100_nrep1.pbs'
            (see NBPMscape\paper\analysis\1_outbreak_simulations\covid_like_varied_sympt_prop\base_plus_20pp\scripts)
    HPC run numbers: 4321681
    Simulation 'maxtime' = 60 days
    Simulation replicates per PBS array run: 10
    Runs in PBS array: 100
    Total sim reps completed: 1000
    CPUs per array run: 1
    RAM per array run: 64GB
    Longest batch run time: ~2.6hrs (~15 minutes per sim rep)
=#
##################################################################################################################################
# (3) Base case + 20pp
#     (3.2) Load simulated outbreak files
##################################################################################################################################

## Load simulated outbreak data with HIGHER symptomatic proportion
sims_bp20 = combine_sim_reps(;  # sim_input_folders in vector format where multiple paths
                          sim_input_folders = [ "$(main_folder)/base_plus_20pp/output/4321681_sims_G_filtered/"]
                        , sim_object_name = "sims_G_filtered" #"sims" # "sims_G_icu_filter" # "sims_G_gp_filter" # "sim_G_filtered"
                        , nrep = 10
                        )
# Check df that shows errors in the combination process
sims_bp20_mismatch_df = sims_bp20[2]
# df is empty

# Refine sims
sims_bp20 = sims_bp20[1]

# Save to file
@save joinpath( main_folder
                , "base_plus_20pp"
                , "covidlike-1.6.1-sims_G_filtered_nrep1000_4321681_sims_bp20.jld2") sims_bp20

# Reload file if necessary
#sims_bp20 = load( joinpath( main_folder
#                            , "base_plus_20pp"
#                            , "covidlike-1.6.1-sims_G_filtered_nrep1000_4321681.jld2")
#                    , "sims_bp20")


##################################################################################################################################
# (3) Base case + 20pp
#     (3.3) Sample infections in primary and tertiary settings
##################################################################################################################################

# Define path to output folder for analyses files
output_folder_bp20 = joinpath( root_dir, "NBPMscape/paper/analysis/2_sampling_analysis/covid_like_varied_sympt_prop/base_plus_20pp" )

## Compute times to detection
# using different model parameters representing different scenarios
# See 'NBPMscape/config/HARISS/config_file_record.xlsx' for summary descriptions of scenarios

# Start with base parameters
initialize_parameters("config/outbreak_params_covid19_like_symptomatic_base_plus_20_20261008.yaml")

# Run analysis for 3 different scenarios (winter x 3 numbers of samples) using GPs (Primary care)
for i in 1:length(gp_sampling_params) # in [1,13,25,37,49,61] # i=1
    gp_tds = gp_td(  sims = sims_bp20
                   , gp_ari_consults = gp_sampling_params[i].gp_ari_consults
                   , gp_ari_swabs    = gp_sampling_params[i].gp_ari_swabs
                   , gp_swabs_mg     = gp_sampling_params[i].gp_swabs_mg
                   )
    CSV.write( joinpath( root_dir, "NBPMscape/paper/analysis/2_sampling_analysis/covid_like_varied_sympt_prop/base_plus_20pp/gp_td"
                        ,"gp_tds_$(gp_sampling_params[i].scenario_number).csv" )
                , gp_tds )
    # Also, save file with 'GP_simid' column removed to reduce file size
    gp_tds_small = gp_tds[:,["sim_n","GP_TD","GP_3TD","n_GP_cases","n_GP_cases_sampled","n_GP_cases_sampled_positive","tinf_relating_to_GP_TD","last_tinf_relating_to_GP_3TD"]]

    CSV.write( joinpath( root_dir, "NBPMscape/paper/analysis/2_sampling_analysis/covid_like_varied_sympt_prop/base_plus_20pp/gp_td_small"
                        ,"small_gp_tds_$(gp_sampling_params[i].scenario_number).csv" )
                , gp_tds_small )
end


# Double check loaded parameters required for ICU sampling
NBPMscape.P.pathogen_type # "virus"
NBPMscape.P.icu_sample_type #"regional" 
NBPMscape.P.sample_icu_cases_version #"number"
NBPMscape.P.icu_only_sample_before_death #true            

# Run analysis for ICU sampling scenarios (2 seasons x 3 numbers of samples x 2 numbers of sampling sites)
for j in 1:length(icu_sampling_params) # [1,13,25,37,49,61,73,74,75,76,77,78] #j=1
    #initialize_parameters("config/influenza_v3/mild/outbreak_params_influenza_like_mild_$(j).yaml")
    
    icu_tds = icu_td( p = NBPMscape.P
                    , sims = sims_bp20
                    , n_icu_samples_per_week   = icu_sampling_params[j].n_icu_samples_per_week
                    , icu_ari_admissions       = icu_sampling_params[j].icu_ari_admissions
                    , nhs_trust_sampling_sites = icu_sampling_params[j].icu_nhs_trust_sampling_sites_file
                    )
    
    CSV.write( joinpath( root_dir, "NBPMscape/paper/analysis/2_sampling_analysis/covid_like_varied_sympt_prop/base_plus_20pp"
                        ,"icu_td/icu_tds_$(icu_sampling_params[j].scenario_number).csv" ), icu_tds )

    # Also, save file with 'ICU_simid' column removed to reduce file size
    icu_tds_small = icu_tds[:,["sim_n","ICU_TD","ICU_3TD","n_ICU_cases","n_ICU_cases_sampled","n_ICU_cases_sampled_positive","tinf_relating_to_ICU_TD","last_tinf_relating_to_ICU_3TD"]]
    
    CSV.write( joinpath( root_dir, "NBPMscape/paper/analysis/2_sampling_analysis/covid_like_varied_sympt_prop/base_plus_20pp"
                        ,"icu_td_small/small_icu_tds_$(icu_sampling_params[j].scenario_number).csv" ), icu_tds_small )
end

##################################################################################################################################
# (3) Base case + 20pp
#     (3.4) Compute median TDs
##################################################################################################################################

# Directory
path_name_bp20 = output_folder_bp20

# Compute median TDs for ICU and GP sampling
#scenario_numbers_2 = [1,13,25,37,49,61,73,74,75,76,77,78] # Note that scenarios 73-78 are with different ICU sites
scenario_numbers_2 = [1,13,25,73,74,75] # Note that scenarios 73-75 are with different ICU sites. 1,13,25 are with 8 ICU sites and 73-75 are with 26 sites.
# Prepare a container to keep the DataFrames
icu_tds_df_dict_bp20 = Dict{Int, DataFrame}()
gp_tds_df_dict_bp20 = Dict{Int, DataFrame}()
# Prepare a df to hold the median TD values
icu_gp_tds_df_bp20 = DataFrame(
                    Scenario_number    = zeros(Int,     length(scenario_numbers_2))
                    ,icu_median_TD     = zeros(Float64, length(scenario_numbers_2))
                    ,icu_TD_iqr_lower  = zeros(Float64, length(scenario_numbers_2))
                    ,icu_TD_iqr_upper  = zeros(Float64, length(scenario_numbers_2))
                    ,icu_TD_std        = zeros(Float64, length(scenario_numbers_2))
                    ,icu_median_3TD    = zeros(Float64, length(scenario_numbers_2))
                    ,icu_3TD_iqr_lower = zeros(Float64, length(scenario_numbers_2))
                    ,icu_3TD_iqr_upper = zeros(Float64, length(scenario_numbers_2))
                    ,icu_3TD_std       = zeros(Float64, length(scenario_numbers_2))
                    ,gp_median_TD      = zeros(Float64, length(scenario_numbers_2))
                    ,gp_TD_iqr_lower   = zeros(Float64, length(scenario_numbers_2))
                    ,gp_TD_iqr_upper   = zeros(Float64, length(scenario_numbers_2))
                    ,gp_TD_std         = zeros(Float64, length(scenario_numbers_2))
                    ,gp_median_3TD     = zeros(Float64, length(scenario_numbers_2))
                    ,gp_3TD_iqr_lower  = zeros(Float64, length(scenario_numbers_2))
                    ,gp_3TD_iqr_upper  = zeros(Float64, length(scenario_numbers_2))
                    ,gp_3TD_std        = zeros(Float64, length(scenario_numbers_2))
                    )

# Loop through scenarios, read in files and compute median TD values
for i in 1:length(scenario_numbers_2) 
    k = scenario_numbers_2[i]
    file_icu = joinpath(path_name_bp20, "icu_td_small","small_icu_tds_$(string(k)).csv")

    # Read CSV into a DataFrame and store in the dictionary with the scenario number (numeric suffix) as key
    icu_tds_df_dict_bp20[k] = CSV.read(file_icu, DataFrame)
    
                        
    # Compute times to detection and populate df
    icu_gp_tds_df_bp20[i,"Scenario_number"] = k
    # ICU
    icu_gp_tds_df_bp20[i,"icu_median_TD"]  = median(icu_tds_df_dict_bp20[k][:,"ICU_TD"])
    icu_gp_tds_df_bp20[i,"icu_median_3TD"] = median(icu_tds_df_dict_bp20[k][:,"ICU_3TD"])
    
    icu_gp_tds_df_bp20[i,"icu_TD_iqr_lower"] = quantile(icu_tds_df_dict_bp20[k][:,"ICU_TD"], 0.25)
    icu_gp_tds_df_bp20[i,"icu_TD_iqr_upper"] = quantile(icu_tds_df_dict_bp20[k][:,"ICU_TD"], 0.75)
    
    icu_gp_tds_df_bp20[i,"icu_3TD_iqr_lower"] = quantile(icu_tds_df_dict_bp20[k][:,"ICU_3TD"], 0.25)
    icu_gp_tds_df_bp20[i,"icu_3TD_iqr_upper"] = quantile(icu_tds_df_dict_bp20[k][:,"ICU_3TD"], 0.75)
    
    icu_gp_tds_df_bp20[i,"icu_TD_std"]  = std(icu_tds_df_dict_bp20[k][:,"ICU_TD"])
    icu_gp_tds_df_bp20[i,"icu_3TD_std"] = std(icu_tds_df_dict_bp20[k][:,"ICU_3TD"])

    # GP
    if k < 73 # Scenario numbers of 73 and above relate to the ICU sampling scenarios with 26 sites
        # Read CSV into a DataFrame and store in the dictionary with the scenario number (numeric suffix) as key
        file_gp  = joinpath(path_name_bp20,  "gp_td_small", "small_gp_tds_$(string(k)).csv")
        gp_tds_df_dict_bp20[k]  = CSV.read(file_gp,  DataFrame)
        
        icu_gp_tds_df_bp20[i,"gp_median_TD"]   = median(gp_tds_df_dict_bp20[k][:,"GP_TD"])
        icu_gp_tds_df_bp20[i,"gp_median_3TD"]  = median(gp_tds_df_dict_bp20[k][:,"GP_3TD"])

        icu_gp_tds_df_bp20[i,"gp_TD_iqr_lower"] = quantile(gp_tds_df_dict_bp20[k][:,"GP_TD"], 0.25)
        icu_gp_tds_df_bp20[i,"gp_TD_iqr_upper"] = quantile(gp_tds_df_dict_bp20[k][:,"GP_TD"], 0.75)
        
        icu_gp_tds_df_bp20[i,"gp_3TD_iqr_lower"] = quantile(gp_tds_df_dict_bp20[k][:,"GP_3TD"], 0.25)
        icu_gp_tds_df_bp20[i,"gp_3TD_iqr_upper"] = quantile(gp_tds_df_dict_bp20[k][:,"GP_3TD"], 0.75)
        
        icu_gp_tds_df_bp20[i,"gp_TD_std"]  = std(gp_tds_df_dict_bp20[k][:,"GP_TD"])
        icu_gp_tds_df_bp20[i,"gp_3TD_std"] = std(gp_tds_df_dict_bp20[k][:,"GP_3TD"])
    end
end

# View
println(icu_gp_tds_df_bp20)

println(icu_tds_df_dict_bp20[1][1:10,:])
println(gp_tds_df_dict_bp20[1][1:10,:])

# Save df
CSV.write( joinpath( output_folder_bp20, "icu_gp_median_tds_bp20_symp.csv" ), icu_gp_tds_df_bp20 )
# Reload if necessary
#icu_gp_tds_df_bp20 = CSV.read( joinpath( output_folder_bp20, "icu_gp_median_tds_bp20_symp.csv" ), DataFrame)

# Save dictionaries containing TD results for all sim reps
@save joinpath(output_folder_bp20, "covidlike-1.6.1-sims_G_filtered_nrep1000_4321681_gp_tds_df_dict_bp20.jld2") gp_tds_df_dict_bp20
@save joinpath(output_folder_bp20, "covidlike-1.6.1-sims_G_filtered_nrep1000_4321681_icu_tds_df_dict_bp20.jld2") icu_tds_df_dict_bp20
# Reload if necessary
#gp_tds_df_dict_bp20  = load( joinpath(output_folder_bp20, "covidlike-1.6.1-sims_G_filtered_nrep1000_4321681_gp_tds_df_dict_bp20.jld2"),  "gp_tds_df_dict_bp20")
#icu_tds_df_dict_bp20 = load( joinpath(output_folder_bp20, "covidlike-1.6.1-sims_G_filtered_nrep1000_4321681_icu_tds_df_dict_bp20.jld2"), "icu_tds_df_dict_bp20")


##################################################################################################################################
# (4) Base case + 40pp
#     (4.1) Information on simulated outbreaks
##################################################################################################################################
#=
    Simulation run on Imperial College High Performance Computing (HPC) system
    using a Julia script and PBS file to manage the process:

    Scripts: Julia - 'covidlike-bp40-1.6.1_w_G_filtered_nrep1.jl'
             PBS   - 'covidlike-bp40-1.6.1_w_G_filtered_array100_nrep1.pbs'
            (see NBPMscape\paper\analysis\1_outbreak_simulations\covid_like_varied_sympt_prop\base_plus_40pp\scripts)
    HPC run numbers: 4321703
    Simulation 'maxtime' = 50 days
    Simulation replicates per PBS array run: 10
    Runs in PBS array: 100
    Total sim reps completed: 1000
    CPUs per array run: 1
    RAM per array run: 64GB
    Longest batch run time: ~1.5hrs (~9 minutes per sim rep)
=#
##################################################################################################################################
# (4) Base case + 40pp
#     (4.2) Load simulated outbreak files
##################################################################################################################################

## Load simulated outbreak data with base case plus 40 percentage points symptomatic proportion
sims_bp40 = combine_sim_reps(;  # sim_input_folders in vector format where multiple paths
                          sim_input_folders = [ "$(main_folder)/base_plus_40pp/output/4321703_sims_G_filtered/"]
                        , sim_object_name = "sims_G_filtered" #"sims" # "sims_G_icu_filter" # "sims_G_gp_filter" # "sim_G_filtered"
                        , nrep = 10
                        )
# Check df that shows errors in the combination process
sims_bp40_mismatch_df = sims_bp40[2]
# df is empty

# Refine sims
sims_bp40 = sims_bp40[1]

# Save to file
@save joinpath( main_folder
                ,"base_plus_40pp"
                ,"covidlike-1.6.1-sims_G_filtered_nrep1000_4321703_sims_bp40.jld2" ) sims_bp40

# Reload file if necessary
#sims_bp40 = load( joinpath( main_folder
#                            ,"base_plus_40pp"
#                            ,"covidlike-1.6.1-sims_G_filtered_nrep1000_4321703_sims_bp40.jld2" )
#                , "sims_bp40" )

##################################################################################################################################
# (4) Base case + 40pp
#     (4.3) Sample infections in primary and tertiary settings
##################################################################################################################################

# Define path to output folder for analyses files
output_folder_bp40 = joinpath( root_dir, "NBPMscape/paper/analysis/2_sampling_analysis/covid_like_varied_sympt_prop/base_plus_40pp" )

## Compute times to detection
# using different model parameters representing different scenarios
# See 'NBPMscape/config/HARISS/config_file_record.xlsx' for summary descriptions of scenarios

# Start with base parameters
initialize_parameters("config/outbreak_params_covid19_like_symptomatic_base_plus_40_20261008.yaml")

# Run analysis for 3 different scenarios (winter x 3 numbers of samples) using GPs (Primary care)
for i in 1:length(gp_sampling_params) # in [1,13,25,37,49,61] # i=1
    gp_tds = gp_td(  sims = sims_bp40
                   , gp_ari_consults = gp_sampling_params[i].gp_ari_consults
                   , gp_ari_swabs    = gp_sampling_params[i].gp_ari_swabs
                   , gp_swabs_mg     = gp_sampling_params[i].gp_swabs_mg
                   )
    CSV.write( joinpath( root_dir, "NBPMscape/paper/analysis/2_sampling_analysis/covid_like_varied_sympt_prop/base_plus_40pp/gp_td"
                        ,"gp_tds_$(gp_sampling_params[i].scenario_number).csv" )
                , gp_tds )
    # Also, save file with 'GP_simid' column removed to reduce file size
    gp_tds_small = gp_tds[:,["sim_n","GP_TD","GP_3TD","n_GP_cases","n_GP_cases_sampled","n_GP_cases_sampled_positive","tinf_relating_to_GP_TD","last_tinf_relating_to_GP_3TD"]]

    CSV.write( joinpath( root_dir, "NBPMscape/paper/analysis/2_sampling_analysis/covid_like_varied_sympt_prop/base_plus_40pp/gp_td_small"
                        ,"small_gp_tds_$(gp_sampling_params[i].scenario_number).csv" )
                , gp_tds_small )
end

# Double check loaded parameters required for ICU sampling
NBPMscape.P.pathogen_type # "virus"
NBPMscape.P.icu_sample_type #"regional" 
NBPMscape.P.sample_icu_cases_version #"number"
NBPMscape.P.icu_only_sample_before_death #true            

# Run analysis for ICU sampling scenarios (2 seasons x 3 numbers of samples x 2 numbers of sampling sites)
for j in 1:length(icu_sampling_params) # [1,13,25,37,49,61,73,74,75,76,77,78] #j=1
       
    icu_tds = icu_td( p = NBPMscape.P
                    , sims = sims_bp40
                    , n_icu_samples_per_week   = icu_sampling_params[j].n_icu_samples_per_week
                    , icu_ari_admissions       = icu_sampling_params[j].icu_ari_admissions
                    , nhs_trust_sampling_sites = icu_sampling_params[j].icu_nhs_trust_sampling_sites_file
                    )
    
    CSV.write( joinpath( root_dir, "NBPMscape/paper/analysis/2_sampling_analysis/covid_like_varied_sympt_prop/base_plus_40pp"
                        ,"icu_td/icu_tds_$(icu_sampling_params[j].scenario_number).csv" ), icu_tds )

    # Also, save file with 'ICU_simid' column removed to reduce file size
    icu_tds_small = icu_tds[:,["sim_n","ICU_TD","ICU_3TD","n_ICU_cases","n_ICU_cases_sampled","n_ICU_cases_sampled_positive","tinf_relating_to_ICU_TD","last_tinf_relating_to_ICU_3TD"]]
    
    CSV.write( joinpath( root_dir, "NBPMscape/paper/analysis/2_sampling_analysis/covid_like_varied_sympt_prop/base_plus_40pp"
                        ,"icu_td_small/small_icu_tds_$(icu_sampling_params[j].scenario_number).csv" ), icu_tds_small )
end


##################################################################################################################################
# (4) Base case + 40pp
#     (4.4) Compute median TDs
##################################################################################################################################

# Read in files and compute median TD values

# Directory
path_name_bp40 = output_folder_bp40

# Compute median TDs for ICU and GP sampling
#scenario_numbers_2 = [1,13,25,37,49,61,73,74,75,76,77,78] # Note that scenarios 73-78 are with different ICU sites
scenario_numbers_2 = [1,13,25,73,74,75] # Note that scenarios 73-75 are with different ICU sites

# Prepare a container to keep the DataFrames
icu_tds_df_dict_bp40 = Dict{Int, DataFrame}()
gp_tds_df_dict_bp40 = Dict{Int, DataFrame}()

# Prepare a df to hold the median TD values
icu_gp_tds_df_bp40 = DataFrame(
                    Scenario_number    = zeros(Int,     length(scenario_numbers_2))
                    ,icu_median_TD     = zeros(Float64, length(scenario_numbers_2))
                    ,icu_TD_iqr_lower  = zeros(Float64, length(scenario_numbers_2))
                    ,icu_TD_iqr_upper  = zeros(Float64, length(scenario_numbers_2))
                    ,icu_TD_std        = zeros(Float64, length(scenario_numbers_2))
                    ,icu_median_3TD    = zeros(Float64, length(scenario_numbers_2))
                    ,icu_3TD_iqr_lower = zeros(Float64, length(scenario_numbers_2))
                    ,icu_3TD_iqr_upper = zeros(Float64, length(scenario_numbers_2))
                    ,icu_3TD_std       = zeros(Float64, length(scenario_numbers_2))
                    ,gp_median_TD      = zeros(Float64, length(scenario_numbers_2))
                    ,gp_TD_iqr_lower   = zeros(Float64, length(scenario_numbers_2))
                    ,gp_TD_iqr_upper   = zeros(Float64, length(scenario_numbers_2))
                    ,gp_TD_std         = zeros(Float64, length(scenario_numbers_2))
                    ,gp_median_3TD     = zeros(Float64, length(scenario_numbers_2))
                    ,gp_3TD_iqr_lower  = zeros(Float64, length(scenario_numbers_2))
                    ,gp_3TD_iqr_upper  = zeros(Float64, length(scenario_numbers_2))
                    ,gp_3TD_std        = zeros(Float64, length(scenario_numbers_2))
                    )

# Loop through scenarios, read in files and compute median TD values
# Tertiary care (ICU)
for i in 1:length(scenario_numbers_2) 
    k = scenario_numbers_2[i]
    file_icu = joinpath(path_name_bp40, "icu_td_small","small_icu_tds_$(string(k)).csv")
    
    # Read CSV into a DataFrame and store in the dictionary with the scenario number (numeric suffix) as key
    icu_tds_df_dict_bp40[k] = CSV.read(file_icu, DataFrame)
        
    # Compute times to detection and populate df
    icu_gp_tds_df_bp40[i,"Scenario_number"] = k

    # ICU
    icu_gp_tds_df_bp40[i,"icu_median_TD"]  = median(icu_tds_df_dict_bp40[k][:,"ICU_TD"])
    icu_gp_tds_df_bp40[i,"icu_median_3TD"] = median(icu_tds_df_dict_bp40[k][:,"ICU_3TD"])
    
    icu_gp_tds_df_bp40[i,"icu_TD_iqr_lower"] = quantile(icu_tds_df_dict_bp40[k][:,"ICU_TD"], 0.25)
    icu_gp_tds_df_bp40[i,"icu_TD_iqr_upper"] = quantile(icu_tds_df_dict_bp40[k][:,"ICU_TD"], 0.75)
    
    icu_gp_tds_df_bp40[i,"icu_3TD_iqr_lower"] = quantile(icu_tds_df_dict_bp40[k][:,"ICU_3TD"], 0.25)
    icu_gp_tds_df_bp40[i,"icu_3TD_iqr_upper"] = quantile(icu_tds_df_dict_bp40[k][:,"ICU_3TD"], 0.75)
    
    icu_gp_tds_df_bp40[i,"icu_TD_std"]  = std(icu_tds_df_dict_bp40[k][:,"ICU_TD"])
    icu_gp_tds_df_bp40[i,"icu_3TD_std"] = std(icu_tds_df_dict_bp40[k][:,"ICU_3TD"])

end

# Only have primary care (GP) results for the first 3 scenarios in this list
for i in 1:length( scenario_numbers_2[1:3] ) 
    
    k = scenario_numbers_2[i]

    file_gp  = joinpath(path_name_bp40, "gp_td_small", "small_gp_tds_$(string(k)).csv")

    # Read CSV into a DataFrame and store in the dictionary with the scenario number (numeric suffix) as key
    gp_tds_df_dict_bp40[k] = CSV.read(file_gp,  DataFrame)
    
    # Compute times to detection and populate df
    icu_gp_tds_df_bp40[i,"Scenario_number"] = k

    # GP
    icu_gp_tds_df_bp40[ i, "gp_median_TD"]   = median(gp_tds_df_dict_bp40[k][:,"GP_TD"])
    icu_gp_tds_df_bp40[ i, "gp_median_3TD"]  = median(gp_tds_df_dict_bp40[k][:,"GP_3TD"])

    icu_gp_tds_df_bp40[ i, "gp_TD_iqr_lower"] = quantile(gp_tds_df_dict_bp40[k][:,"GP_TD"], 0.25)
    icu_gp_tds_df_bp40[ i, "gp_TD_iqr_upper"] = quantile(gp_tds_df_dict_bp40[k][:,"GP_TD"], 0.75)
    
    icu_gp_tds_df_bp40[ i, "gp_3TD_iqr_lower"] = quantile(gp_tds_df_dict_bp40[k][:,"GP_3TD"], 0.25)
    icu_gp_tds_df_bp40[ i, "gp_3TD_iqr_upper"] = quantile(gp_tds_df_dict_bp40[k][:,"GP_3TD"], 0.75)
    
    icu_gp_tds_df_bp40[ i, "gp_TD_std"]  = std(gp_tds_df_dict_bp40[k][:,"GP_TD"])
    icu_gp_tds_df_bp40[ i, "gp_3TD_std"] = std(gp_tds_df_dict_bp40[k][:,"GP_3TD"])

end

# View
println(icu_gp_tds_df_bp40)

println(icu_tds_df_dict_bp40[1][1:10,:])
println( gp_tds_df_dict_bp40[1][1:10,:])

# Save df
CSV.write( joinpath( output_folder_bp40, "icu_gp_median_tds_bp40_symp.csv" ), icu_gp_tds_df_bp40 )

# Reload if necessary
#icu_gp_tds_df_bp40 = CSV.read( joinpath( output_folder_bp40, "icu_gp_median_tds_bp40_symp.csv" ), DataFrame)

# Save dictionaries containing TD results for all sim reps
@save joinpath(output_folder_bp40, "covidlike-1.6.1-sims_G_filtered_nrep1000_4321703_gp_tds_df_dict_bp40.jld2") gp_tds_df_dict_bp40
@save joinpath(output_folder_bp40, "covidlike-1.6.1-sims_G_filtered_nrep1000_4321703_icu_tds_df_dict_bp40.jld2") icu_tds_df_dict_bp40

# Reload if necessary
#gp_tds_df_dict_bp40  = load( joinpath(output_folder_bp40, "covidlike-1.6.1-sims_G_filtered_nrep1000_4321703_gp_tds_df_dict_bp40.jld2"),   "gp_tds_df_dict_bp40")
#icu_tds_df_dict_bp40 = load( joinpath(output_folder_bp40, "covidlike-1.6.1-sims-filtered_G_nrep1000_4321703_icu_tds_df_dict_bp40.jld2"), "icu_tds_df_dict_bp40")

# Reduce memory
sims_bp40 = nothing


##################################################################################################################################
# (5) 100% symptomatic
#     (5.1) Information on simulated outbreaks
##################################################################################################################################
#=
    Simulation run on Imperial College High Performance Computing (HPC) system
    using a Julia script and PBS file to manage the process:

    Scripts: Julia - 'covidlike-100-1.6.1_w_G_filtered_nrep1.jl'
             PBS   - 'covidlike-100-1.6.1_w_G_filtered_array100_nrep1.pbs'
            (see NBPMscape\paper\analysis\1_outbreak_simulations\covid_like_varied_sympt_prop\100pc_symptomatic\scripts)
    HPC run numbers: 4321713
    Simulation 'maxtime' = 45 days
    Simulation replicates per PBS array run: 10
    Runs in PBS array: 100
    Total sim reps completed: 1000
    CPUs per array run: 1
    RAM per array run: 64GB
    Longest batch run time: ~1.8hrs (~11 minutes per sim rep)
=#
##################################################################################################################################
# (5) 100% symptomatic
#     (5.2) Load simulated outbreak files
##################################################################################################################################

## Load simulated outbreak data
sims_100 = combine_sim_reps(;  # sim_input_folders in vector format where multiple paths
                          sim_input_folders = [ "$(main_folder)/100pc_symptomatic/output/4321713_sims_G_filtered/"]
                        , sim_object_name = "sims_G_filtered" #"sims" # "sims_G_icu_filter" # "sims_G_gp_filter" # "sim_G_filtered"
                        , nrep = 10
                        )
# Check df that shows errors in the combination process
sims_100_mismatch_df = sims_100[2]
# df is empty

# Refine sims
sims_100 = sims_100[1]

# Save to file
@save joinpath( main_folder
                , "100pc_symptomatic"
                , "covidlike-1.6.1-sims_G_filtered_nrep1000_4321713_sims_100.jld2" ) sims_100

# Reload file if necessary
#sims_100 = load( joinpath( main_folder
#                        , "100pc_symptomatic"
#                        , "covidlike-1.6.1-sims_G_filtered_nrep1000_4321713_sims_100.jld2" )
#                , "sims_100")

##################################################################################################################################
# (5) 100% symptomatic
#     (5.3) Sample infections in primary, secondary and tertiary settings
##################################################################################################################################

# Define path to output folder for analyses files
output_folder_100 = joinpath( root_dir, "NBPMscape/paper/analysis/2_sampling_analysis/covid_like_varied_sympt_prop/100pc_symptomatic" )

## Compute times to detection
# using different model parameters representing different scenarios
# See config_file_record.xlsx for summary descriptions of scenarios

## Compute times to detection
# using different model parameters representing different scenarios
# See 'NBPMscape/config/HARISS/config_file_record.xlsx' for summary descriptions of scenarios

# Start with base parameters
initialize_parameters("config/outbreak_params_covid19_like_symptomatic_100_20261008.yaml")

# Run analysis for 3 different scenarios (winter x 3 numbers of samples) using GPs (Primary care)
for i in 1:length(gp_sampling_params) # in [1,13,25,37,49,61] # i=1
    gp_tds = gp_td(  sims = sims_100
                   , gp_ari_consults = gp_sampling_params[i].gp_ari_consults
                   , gp_ari_swabs    = gp_sampling_params[i].gp_ari_swabs
                   , gp_swabs_mg     = gp_sampling_params[i].gp_swabs_mg
                   )
    CSV.write( joinpath( root_dir, "NBPMscape/paper/analysis/2_sampling_analysis/covid_like_varied_sympt_prop/100pc_symptomatic/gp_td"
                        ,"gp_tds_$(gp_sampling_params[i].scenario_number).csv" )
                , gp_tds )
    
    # Also, save file with 'GP_simid' column removed to reduce file size
    gp_tds_small = gp_tds[:,["sim_n","GP_TD","GP_3TD","n_GP_cases","n_GP_cases_sampled","n_GP_cases_sampled_positive","tinf_relating_to_GP_TD","last_tinf_relating_to_GP_3TD"]]

    CSV.write( joinpath( root_dir, "NBPMscape/paper/analysis/2_sampling_analysis/covid_like_varied_sympt_prop/100pc_symptomatic/gp_td_small"
                        ,"small_gp_tds_$(gp_sampling_params[i].scenario_number).csv" )
                , gp_tds_small )
end

# Double check loaded parameters required for ICU sampling
NBPMscape.P.pathogen_type # "virus"
NBPMscape.P.icu_sample_type #"regional" 
NBPMscape.P.sample_icu_cases_version #"number"
NBPMscape.P.icu_only_sample_before_death #true            

# Run analysis for ICU sampling scenarios (2 seasons x 3 numbers of samples x 2 numbers of sampling sites)
for j in 1:length(icu_sampling_params) # [1,13,25,37,49,61,73,74,75,76,77,78] #j=1
       
    icu_tds = icu_td( p = NBPMscape.P
                    , sims = sims_100
                    , n_icu_samples_per_week   = icu_sampling_params[j].n_icu_samples_per_week
                    , icu_ari_admissions       = icu_sampling_params[j].icu_ari_admissions
                    , nhs_trust_sampling_sites = icu_sampling_params[j].icu_nhs_trust_sampling_sites_file
                    )
    
    CSV.write( joinpath( root_dir, "NBPMscape/paper/analysis/2_sampling_analysis/covid_like_varied_sympt_prop/100pc_symptomatic"
                        ,"icu_td/icu_tds_$(icu_sampling_params[j].scenario_number).csv" ), icu_tds )

    # Also, save file with 'ICU_simid' column removed to reduce file size
    icu_tds_small = icu_tds[:,["sim_n","ICU_TD","ICU_3TD","n_ICU_cases","n_ICU_cases_sampled","n_ICU_cases_sampled_positive","tinf_relating_to_ICU_TD","last_tinf_relating_to_ICU_3TD"]]
    
    CSV.write( joinpath( root_dir, "NBPMscape/paper/analysis/2_sampling_analysis/covid_like_varied_sympt_prop/100pc_symptomatic"
                        ,"icu_td_small/small_icu_tds_$(icu_sampling_params[j].scenario_number).csv" ), icu_tds_small )
end

##################################################################################################################################
# (5) 100% symptomatic
#     (5.4) Compute median TDs
##################################################################################################################################

# Read in files and compute median TD values

# Directory and base file name
path_name_100 = output_folder_100

# Compute median TDs for ICU and GP sampling
#scenario_numbers_2 = [1,13,25,37,49,61,73,74,75,76,77,78] # Note that scenarios 73-78 are with different ICU sites
scenario_numbers_2 = [1,13,25,73,74,75] # Note that scenarios 73-75 are with different ICU sites

# Prepare a container to keep the DataFrames
icu_tds_df_dict_100 = Dict{Int, DataFrame}()
gp_tds_df_dict_100 = Dict{Int, DataFrame}()

# Prepare a df to hold the median TD values
icu_gp_tds_df_100 = DataFrame(
                    Scenario_number    = zeros(Int,     length(scenario_numbers_2))
                    ,icu_median_TD     = zeros(Float64, length(scenario_numbers_2))
                    ,icu_TD_iqr_lower  = zeros(Float64, length(scenario_numbers_2))
                    ,icu_TD_iqr_upper  = zeros(Float64, length(scenario_numbers_2))
                    ,icu_TD_std        = zeros(Float64, length(scenario_numbers_2))
                    ,icu_median_3TD    = zeros(Float64, length(scenario_numbers_2))
                    ,icu_3TD_iqr_lower = zeros(Float64, length(scenario_numbers_2))
                    ,icu_3TD_iqr_upper = zeros(Float64, length(scenario_numbers_2))
                    ,icu_3TD_std       = zeros(Float64, length(scenario_numbers_2))
                    ,gp_median_TD      = zeros(Float64, length(scenario_numbers_2))
                    ,gp_TD_iqr_lower   = zeros(Float64, length(scenario_numbers_2))
                    ,gp_TD_iqr_upper   = zeros(Float64, length(scenario_numbers_2))
                    ,gp_TD_std         = zeros(Float64, length(scenario_numbers_2))
                    ,gp_median_3TD     = zeros(Float64, length(scenario_numbers_2))
                    ,gp_3TD_iqr_lower  = zeros(Float64, length(scenario_numbers_2))
                    ,gp_3TD_iqr_upper  = zeros(Float64, length(scenario_numbers_2))
                    ,gp_3TD_std        = zeros(Float64, length(scenario_numbers_2))
                    )

# Loop through scenarios, read in files and compute median TD values
for i in 1:length(scenario_numbers_2) 
    k = scenario_numbers_2[i]
    
    # Construct the file name
    file_icu = joinpath(path_name_100, "icu_td_small", "small_icu_tds_$(string(k)).csv")

    # Read CSV into a DataFrame and store in the dictionary with the scenario number (numeric suffix) as key
    icu_tds_df_dict_100[k] = CSV.read(file_icu, DataFrame)

    # Compute times to detection and populate df
    icu_gp_tds_df_100[i,"Scenario_number"] = k

    # ICU
    icu_gp_tds_df_100[i,"icu_median_TD"]  = median(icu_tds_df_dict_100[k][:,"ICU_TD"])
    icu_gp_tds_df_100[i,"icu_median_3TD"] = median(icu_tds_df_dict_100[k][:,"ICU_3TD"])
    
    icu_gp_tds_df_100[i,"icu_TD_iqr_lower"] = quantile(icu_tds_df_dict_100[k][:,"ICU_TD"], 0.25)
    icu_gp_tds_df_100[i,"icu_TD_iqr_upper"] = quantile(icu_tds_df_dict_100[k][:,"ICU_TD"], 0.75)
    
    icu_gp_tds_df_100[i,"icu_3TD_iqr_lower"] = quantile(icu_tds_df_dict_100[k][:,"ICU_3TD"], 0.25)
    icu_gp_tds_df_100[i,"icu_3TD_iqr_upper"] = quantile(icu_tds_df_dict_100[k][:,"ICU_3TD"], 0.75)
    
    icu_gp_tds_df_100[i,"icu_TD_std"]  = std(icu_tds_df_dict_100[k][:,"ICU_TD"])
    icu_gp_tds_df_100[i,"icu_3TD_std"] = std(icu_tds_df_dict_100[k][:,"ICU_3TD"])

end

# Only have GP results for the first 3 scenarios in this list
for i in 1:length( scenario_numbers_2[1:3] ) 
    
    k = scenario_numbers_2[i]

    # Construct the file name
    file_gp  = joinpath(path_name_100, "gp_td_small", "small_gp_tds_$(string(k)).csv")

    # Read CSV into a DataFrame and store in the dictionary with the scenario number (numeric suffix) as key
    gp_tds_df_dict_100[k] = CSV.read(file_gp,  DataFrame)
    
    # Compute times to detection and populate df
    icu_gp_tds_df_100[i,"Scenario_number"] = k

    # GP
    icu_gp_tds_df_100[ i, "gp_median_TD"]   = median(gp_tds_df_dict_100[k][:,"GP_TD"])
    icu_gp_tds_df_100[ i, "gp_median_3TD"]  = median(gp_tds_df_dict_100[k][:,"GP_3TD"])

    icu_gp_tds_df_100[ i, "gp_TD_iqr_lower"] = quantile(gp_tds_df_dict_100[k][:,"GP_TD"], 0.25)
    icu_gp_tds_df_100[ i, "gp_TD_iqr_upper"] = quantile(gp_tds_df_dict_100[k][:,"GP_TD"], 0.75)
    
    icu_gp_tds_df_100[ i, "gp_3TD_iqr_lower"] = quantile(gp_tds_df_dict_100[k][:,"GP_3TD"], 0.25)
    icu_gp_tds_df_100[ i, "gp_3TD_iqr_upper"] = quantile(gp_tds_df_dict_100[k][:,"GP_3TD"], 0.75)
    
    icu_gp_tds_df_100[ i, "gp_TD_std"]  = std(gp_tds_df_dict_100[k][:,"GP_TD"])
    icu_gp_tds_df_100[ i, "gp_3TD_std"] = std(gp_tds_df_dict_100[k][:,"GP_3TD"])

end

# View
println(icu_gp_tds_df_100)

println(icu_tds_df_dict_100[1][1:10,:])
println( gp_tds_df_dict_100[1][1:10,:])

# Save df
CSV.write( joinpath( output_folder_100, "icu_gp_median_tds_100_symp.csv" ), icu_gp_tds_df_100 )
# Reload if necessary
#icu_gp_tds_df_100 = CSV.read( joinpath( output_folder_100, "icu_gp_median_tds_100_symp.csv" ), DataFrame)

# Save dictionaries containing TD results for all sim reps
@save joinpath(output_folder_100, "covidlike-1.6.1-sims_G_filtered_nrep1000_4321713_gp_tds_df_dict_100.jld2") gp_tds_df_dict_100
@save joinpath(output_folder_100, "covidlike-1.6.1-sims_G_filtered_nrep1000_4321713_icu_tds_df_dict_100.jld2") icu_tds_df_dict_100
# Reload if necessary
#gp_tds_df_dict_100  = load( joinpath(output_folder_100, "covidlike-1.6.1-sims_G_filtered_nrep1000_4321713_gp_tds_df_dict_100.jld2"),   "gp_tds_df_dict_100")
#icu_tds_df_dict_100 = load( joinpath(output_folder_100, "covidlike-1.6.1-sims_G_filtered_nrep1000_4321713_icu_tds_df_dict_100.jld2"), "icu_tds_df_dict_100")

# Reduce memory
sims_100 = nothing

########################################################################
# (6) Combine individual simulation results into a single df and save .csv
#   (6.1) Read in TD data dictionaries
########################################################################

# Reload dictionaries containing individual simulation replicate TD results
base_folder = joinpath( root_dir, "NBPMscape/paper/analysis/2_sampling_analysis" )

# Covid-like base symptomatic proportion - 20pp across all ages
bm20_folder = "covid_like_varied_sympt_prop/base_minus_20pp"
gp_tds_df_dict_bm20  = load( joinpath(base_folder, bm20_folder, "covidlike-1.6.1-sims_G_filtered_nrep1000_4301080_gp_tds_df_dict_bm20.jld2"),   "gp_tds_df_dict_bm20")
icu_tds_df_dict_bm20 = load( joinpath(base_folder, bm20_folder, "covidlike-1.6.1-sims_G_filtered_nrep1000_4301080_icu_tds_df_dict_bm20.jld2"), "icu_tds_df_dict_bm20")

# Covid-like base symptomatic proportion = base case
bc_folder = "covid_like"#/geo_optimisation"
gp_tds_df_dict_bc  = load( joinpath(base_folder, bc_folder, "covidlike-1.6.1_nrep1000_4300448_gp_tds_df_dict.jld2"),          "gp_tds_df_dict")
icu_tds_df_dict_bc = load( joinpath(base_folder, bc_folder, "geo_optimisation/covidlike-1.6.1_nrep1000_4300448_icu_tds_df_dict.jld2"), "icu_tds_df_dict")
# Trim ICU TD results dictionary to mSCAPE only
#icu_tds_df_dict_bc_mscape = Dict()
#for key in keys(icu_tds_df_dict_bc) # key="icu_tds_seasonal_ari_793_1000_sampling_sites_5blocks_20sites"
#    if occursin(["mscape","1440"],key)
#        icu_tds_df_dict_bc_mscape[key] = icu_tds_df_dict_bc[key]
#        println(key)
#    end
#end

# Covid-like base symptomatic proportion + 20pp across all ages
bp20_folder    = "covid_like_varied_sympt_prop/base_plus_20pp"
gp_tds_df_dict_bp20  = load( joinpath(base_folder, bp20_folder, "covidlike-1.6.1-sims_G_filtered_nrep1000_4321681_gp_tds_df_dict_bp20.jld2"),    "gp_tds_df_dict_bp20")
icu_tds_df_dict_bp20 = load( joinpath(base_folder, bp20_folder, "covidlike-1.6.1-sims_G_filtered_nrep1000_4321681_icu_tds_df_dict_bp20.jld2"),  "icu_tds_df_dict_bp20")

# Covid-like base symptomatic proportion + 40pp across all ages
bp40_folder    = "covid_like_varied_sympt_prop/base_plus_40pp"
gp_tds_df_dict_bp40  = load( joinpath(base_folder, bp40_folder,    "covidlike-1.6.1-sims_G_filtered_nrep1000_4321703_gp_tds_df_dict_bp40.jld2"),    "gp_tds_df_dict_bp40")
icu_tds_df_dict_bp40 = load( joinpath(base_folder, bp40_folder,    "covidlike-1.6.1-sims_G_filtered_nrep1000_4321703_icu_tds_df_dict_bp40.jld2"),  "icu_tds_df_dict_bp40")

# Covid-like but 100% symptomatic proportion across all ages
bp100_folder    = "covid_like_varied_sympt_prop/100pc_symptomatic"
gp_tds_df_dict_100  = load( joinpath(base_folder, bp100_folder,    "covidlike-1.6.1-sims_G_filtered_nrep1000_4321713_gp_tds_df_dict_100.jld2"),    "gp_tds_df_dict_100")
icu_tds_df_dict_100 = load( joinpath(base_folder, bp100_folder,    "covidlike-1.6.1-sims_G_filtered_nrep1000_4321713_icu_tds_df_dict_100.jld2"),  "icu_tds_df_dict_100")


################################################################################################################################################
# (6) Combine individual simulation results into a single df and save .csv
#   (6.2) Merge into single dataframe
################################################################################################################################################
## Merge data into a single df

# Initialise df to store results
td_symp_df = DataFrame(   symptomatic_proportion = "dummy"
                        , sampling_method = "dummy"
                        , scenario_number = 1
                        , scenario_name = "text"
                        , sim_rep = 1
                        , TD = -1.0
                        )

td_symp_df_base_case = copy(td_symp_df)

# Helper function to add data from dictionaries containing TDs to df
function add_td_dict_data_to_df(; td_dict, td_col_name, symptomatic_proportion_name, sampling_method, results_df, base_case_symp = false)
    for key in keys(td_dict) # key=1
        for i in 1:nrow(td_dict[key]) #i=1
            sp = symptomatic_proportion_name #"base_minus_20pp"
            sm = sampling_method #"gp"
            scen = key
            if base_case_symp
                sr = i
                scen_number = -1
                scen_name = key
            else
                scen_number = key
                scen_name = "TBC"
                sr = td_dict[key][i,:].sim_n
            end
            td = td_dict[key][i,Symbol(td_col_name)]

            row_to_add = [sp, sm, scen_number, scen_name, sr, td]

            push!(results_df, row_to_add)

        end
    end
    return results_df
end

# Add data to df from LOWER/bm20 symptomatic proportion (base case COVID-like minus 20 percentage points across all ages)
td_symp_df = add_td_dict_data_to_df( td_dict = gp_tds_df_dict_bm20, td_col_name = "GP_TD", symptomatic_proportion_name = "base_minus_20pp", sampling_method = "RCGP", results_df = td_symp_df )

# Remove 1st dummy row
delete!(td_symp_df, [1])

td_symp_df = add_td_dict_data_to_df( td_dict = icu_tds_df_dict_bm20, td_col_name = "ICU_TD", symptomatic_proportion_name = "base_minus_20pp", sampling_method = "ICU",    results_df = td_symp_df )

# Add data to df from HIGHER/bp20 symptomatic proportion (base case COVID-like plus 20 percentage points across all ages)
td_symp_df = add_td_dict_data_to_df( td_dict = gp_tds_df_dict_bp20,  td_col_name =  "GP_TD", symptomatic_proportion_name = "base_plus_20pp", sampling_method = "RCGP",   results_df = td_symp_df )
td_symp_df = add_td_dict_data_to_df( td_dict = icu_tds_df_dict_bp20, td_col_name = "ICU_TD", symptomatic_proportion_name = "base_plus_20pp", sampling_method = "ICU",    results_df = td_symp_df )

# Add data to df from base plus 40 percentage points symptomatic proportion
td_symp_df = add_td_dict_data_to_df( td_dict = gp_tds_df_dict_bp40,  td_col_name =  "GP_TD", symptomatic_proportion_name = "base_plus_40pp", sampling_method = "RCGP",   results_df = td_symp_df )
td_symp_df = add_td_dict_data_to_df( td_dict = icu_tds_df_dict_bp40, td_col_name = "ICU_TD", symptomatic_proportion_name = "base_plus_40pp", sampling_method = "ICU",    results_df = td_symp_df )

# Add data to df from 100% symptomatic proportion
td_symp_df = add_td_dict_data_to_df( td_dict = gp_tds_df_dict_100,  td_col_name =  "GP_TD", symptomatic_proportion_name = "all_100", sampling_method = "RCGP",   results_df = td_symp_df )
td_symp_df = add_td_dict_data_to_df( td_dict = icu_tds_df_dict_100, td_col_name = "ICU_TD", symptomatic_proportion_name = "all_100", sampling_method = "ICU",    results_df = td_symp_df )

# Add data to df from BASE symptomatic proportion
# NOTE: the base case results have scenario names rather than scenario numbers
td_symp_df_base_case = add_td_dict_data_to_df( td_dict = gp_tds_df_dict_bc,  td_col_name =  "GP_TD", symptomatic_proportion_name = "base", sampling_method = "RCGP",   results_df = td_symp_df_base_case, base_case_symp = true)
td_symp_df_base_case = add_td_dict_data_to_df( td_dict = icu_tds_df_dict_bc, td_col_name = "ICU_TD", symptomatic_proportion_name = "base", sampling_method = "ICU",    results_df = td_symp_df_base_case, base_case_symp = true)

# Remove 1st dummy row
delete!(td_symp_df_base_case, [1])

## Filter for base case scenarios that match those in the symptomatic proportion variation
println(unique(td_symp_df.scenario_number)) # [13, 25, 73, 74, 75, 1]
# See 'config/2026_02_18/config file record.xlsx' for descriptions associated with scenario numbers. Summarised below
#= 
Scenario number     Samples per week    Season      Sample distribution     Age targeting       ICU sites
1                   300                 Winter      Equal                   No specification    8
13                  500                 Winter      Equal                   No specification    8
25                  1000                Winter      Equal                   No specification    8
73                  300                 Winter      Equal                   No specification    26
74                  500                 Winter      Equal                   No specification    26
75                  1000                Winter      Equal                   No specification    26
=#

# Add scenario numbers for base case data that only has scenario names
for i in unique(td_symp_df_base_case.scenario_name)
    println(i)
end

# Create dictionary of scenario names in base case results that match results for varied symptomatic proportion
scenario_name_number_dict = Dict(
                                    "gp_tds_winter_300samples" => 1
                                    ,"gp_tds_winter_500samples" => 13
                                    ,"gp_tds_winter_1000samples" => 25
                                    ,"icu_tds_seasonal_ari_1440_300_sampling_sites_mscape_8sites" => 1
                                    ,"icu_tds_seasonal_ari_1440_500_sampling_sites_mscape_8sites" => 13
                                    ,"icu_tds_seasonal_ari_1440_1000_sampling_sites_mscape_8sites" => 25
                                    ,"icu_tds_seasonal_ari_1440_300_sampling_sites_mscape_26sites" => 73
                                    ,"icu_tds_seasonal_ari_1440_500_sampling_sites_mscape_26sites" => 74
                                    ,"icu_tds_seasonal_ari_1440_1000_sampling_sites_mscape_26sites" => 75
)

# Filter 
td_symp_df_base_case_filtered = filter!( row -> row.scenario_name in keys(scenario_name_number_dict), td_symp_df_base_case)

# Add scenario numbers for base case results
for i in 1:nrow(td_symp_df_base_case_filtered) # i=1
    td_symp_df_base_case_filtered[ i, :scenario_number] = scenario_name_number_dict[td_symp_df_base_case_filtered[ i, :scenario_name]]
end

unique(td_symp_df_base_case_filtered.scenario_number) # 13,  1, 25, 73, 75, 74

unique(td_symp_df_base_case_filtered.scenario_name) 
#12-element Vector{String}:
# "gp_tds_winter_500samples"
# "gp_tds_winter_300samples"
# "gp_tds_winter_1000samples"
# "icu_tds_seasonal_ari_1440_300_sampling_sites_mscape_26sites"
# "icu_tds_seasonal_ari_1440_1000_sampling_sites_mscape_26sites"
# "icu_tds_seasonal_ari_1440_500_sampling_sites_mscape_26sites"
# "icu_tds_seasonal_ari_1440_300_sampling_sites_mscape_8sites"
# "icu_tds_seasonal_ari_1440_1000_sampling_sites_mscape_8sites"
# "icu_tds_seasonal_ari_1440_500_sampling_sites_mscape_8sites"

# Combine symptomatic proportion results with relevant base case Results
td_symp_all_df = vcat(td_symp_df_base_case_filtered, td_symp_df )

# Save df
CSV.write( joinpath( base_folder, "covid_like_varied_sympt_prop/td_symp_df_GP_ICU.csv" ), td_symp_all_df )

################################################################################################################################################
# (7) Plots and TD differences
################################################################################################################################################

# Note that this only plots the base case symptomatic proportion, base case -20pp and base case +20pp

symp_tds_dict = Dict{String, DataFrame}()

# Base symptomatic proportion
path_base    = joinpath( root_dir, "NBPMscape/paper/analysis/2_sampling_analysis/covid_like" )
symp_tds_dict["base_symp_icu8_300s_winter"]  = CSV.read( joinpath( path_base, "geo_optimisation/icu_td/mscape_small", "small_icu_tds_seasonal_ari_1440_300_sampling_sites_mscape_8sites.csv"), DataFrame) #"small_icu_tds_1.csv"), DataFrame)
symp_tds_dict["base_symp_icu26_300s_winter"] = CSV.read( joinpath( path_base, "geo_optimisation/icu_td/mscape_small", "small_icu_tds_seasonal_ari_1440_300_sampling_sites_mscape_26sites.csv"), DataFrame) #"small_icu_tds_73.csv"), DataFrame)
symp_tds_dict["base_symp_gp_300s_winter"]    = CSV.read( joinpath( path_base, "gp_td_small", "small_gp_tds_1.csv"), DataFrame)

# Base symptomatic proportion minus 20 percentage points across all ages ("Lower" symptomatic proportion)
path_bm20    = joinpath( root_dir, "NBPMscape/paper/analysis/2_sampling_analysis/covid_like_varied_sympt_prop/base_minus_20pp" )
symp_tds_dict["lower_symp_icu8_300s_winter"]  = CSV.read( joinpath(path_bm20, "icu_td_small", "small_icu_tds_1.csv"),    DataFrame)
symp_tds_dict["lower_symp_icu26_300s_winter"] = CSV.read( joinpath(path_bm20, "icu_td_small", "small_icu_tds_73.csv"),    DataFrame)
symp_tds_dict["lower_symp_gp_300s_winter"]    = CSV.read( joinpath(path_bm20, "gp_td_small",  "small_gp_tds_1.csv"),    DataFrame)

# Base symptomatic proportion plus 20 percentage points across all ages ("Higher" symptomatic proportion)
path_bp20   = joinpath( root_dir, "NBPMscape/paper/analysis/2_sampling_analysis/covid_like_varied_sympt_prop/base_plus_20pp" )
symp_tds_dict["higher_symp_icu8_300s_winter"]  = CSV.read( joinpath(path_bp20, "icu_td_small", "small_icu_tds_1.csv"),    DataFrame)
symp_tds_dict["higher_symp_icu26_300s_winter"] = CSV.read( joinpath(path_bp20, "icu_td_small", "small_icu_tds_73.csv"),    DataFrame)
symp_tds_dict["higher_symp_gp_300s_winter"]    = CSV.read( joinpath(path_bp20, "gp_td_small",  "small_gp_tds_1.csv"),    DataFrame)



# Plot jitter
# 300 samples in winter scenario so comparable with other figures in paper

# Define colours
#https://juliagraphics.github.io/Colors.jl/stable/namedcolors/
rcgp_color_td = :gold        ; rcgp_color_cum_inf = :gold
mscape_8_color_td = :mediumpurple3 ; mscape_8_color_cum_inf = :mediumpurple3
mscape_26_color_td = :royalblue1   ; mscape_26_color_cum_inf = :royalblue1

lower_symp_color   = :coral
base_symp_color    = :purple
higher_symp_color  = :midnightblue

p_covid_symp_base_gp        = kernel_box_jitter_plot(;x = collect(skipmissing(symp_tds_dict["base_symp_gp_300s_winter"   ][:,:GP_TD])),  plot_color = rcgp_color_td,samp_strategy="RCGP",n_samples=300,n_sites="",n_blocks="",season="Winter",x_label="",x_min=0,x_max=80, annot_x_pos="middle") #, main_title = catchpop_keys[1])
p_covid_symp_base_mscape_8  = kernel_box_jitter_plot(;x = collect(skipmissing(symp_tds_dict["base_symp_icu8_300s_winter" ][:,:ICU_TD])), plot_color = mscape_8_color_td,samp_strategy="mSCAPE",n_samples=300,n_sites="8",n_blocks="",season="Winter",x_label="",x_min=0,x_max=80, annot_x_pos="middle") #, main_title = catchpop_keys[1])
p_covid_symp_base_mscape_26 = kernel_box_jitter_plot(;x = collect(skipmissing(symp_tds_dict["base_symp_icu26_300s_winter"][:,:ICU_TD])), plot_color = mscape_26_color_td,samp_strategy="mSCAPE",n_samples=300,n_sites="26",n_blocks="",season="Winter",x_label="",x_min=0,x_max=80, annot_x_pos="middle") #, main_title = catchpop_keys[1])

#p_covid_symp_lower_gp        = kernel_box_jitter_plot(;x = collect( skipmissing(symp_tds_dict["lower_symp_gp_300s_winter"][:,:GP_TD])), plot_color = rcgp_color_td,samp_strategy="RCGP",n_samples=300,n_sites="",n_blocks="",season="Winter",x_label="",x_min=0,x_max=80, annot_x_pos="middle") #, main_title = catchpop_keys[1])
p_covid_symp_lower_gp        = kernel_box_jitter_plot(;x = collect( skipmissing( filter(x -> isfinite(x), symp_tds_dict["lower_symp_gp_300s_winter"][:,:GP_TD] ))), plot_color = rcgp_color_td,samp_strategy="RCGP",n_samples=300,n_sites="",n_blocks="",season="Winter",x_label="",x_min=0,x_max=80, annot_x_pos="middle") #, main_title = catchpop_keys[1])
p_covid_symp_lower_mscape_8  = kernel_box_jitter_plot(;x = collect(skipmissing(symp_tds_dict["lower_symp_icu8_300s_winter"][:,:ICU_TD])), plot_color = mscape_8_color_td,samp_strategy="mSCAPE",n_samples=300,n_sites="8",n_blocks="",season="Winter",x_label="",x_min=0,x_max=80, annot_x_pos="middle") #, main_title = catchpop_keys[1])
p_covid_symp_lower_mscape_26 = kernel_box_jitter_plot(;x = collect(skipmissing(symp_tds_dict["lower_symp_icu26_300s_winter"][:,:ICU_TD])), plot_color = mscape_26_color_td,samp_strategy="mSCAPE",n_samples=300,n_sites="26",n_blocks="",season="Winter",x_label="",x_min=0,x_max=80, annot_x_pos="middle") #, main_title = catchpop_keys[1])

#p_covid_symp_higher_gp        = kernel_box_jitter_plot(;x = collect(skipmissing(symp_tds_dict["higher_symp_gp_300s_winter"][:,:GP_TD])), plot_color = rcgp_color_td,samp_strategy="RCGP",n_samples=300,n_sites="",n_blocks="",season="Winter",x_label="",x_min=0,x_max=80, annot_x_pos="middle") #, main_title = catchpop_keys[1])
p_covid_symp_higher_gp        = kernel_box_jitter_plot(;x = collect(skipmissing( filter(x -> isfinite(x), symp_tds_dict["higher_symp_gp_300s_winter"][:,:GP_TD]))), plot_color = rcgp_color_td,samp_strategy="RCGP",n_samples=300,n_sites="",n_blocks="",season="Winter",x_label="",x_min=0,x_max=80, annot_x_pos="middle") #, main_title = catchpop_keys[1])
p_covid_symp_higher_mscape_8  = kernel_box_jitter_plot(;x = collect(skipmissing(symp_tds_dict["higher_symp_icu8_300s_winter"][:,:ICU_TD])), plot_color = mscape_8_color_td,samp_strategy="mSCAPE",n_samples=300,n_sites="8",n_blocks="",season="Winter",x_label="",x_min=0,x_max=80, annot_x_pos="middle") #, main_title = catchpop_keys[1])
p_covid_symp_higher_mscape_26 = kernel_box_jitter_plot(;x = collect(skipmissing(symp_tds_dict["higher_symp_icu26_300s_winter"][:,:ICU_TD])), plot_color = mscape_26_color_td,samp_strategy="mSCAPE",n_samples=300,n_sites="26",n_blocks="",season="Winter",x_label="",x_min=0,x_max=80, annot_x_pos="middle") #, main_title = catchpop_keys[1])

Plots.plot( p_covid_symp_lower_gp,  p_covid_symp_lower_mscape_8,  p_covid_symp_lower_mscape_26
           ,p_covid_symp_base_gp,   p_covid_symp_base_mscape_8,   p_covid_symp_base_mscape_26
           ,p_covid_symp_higher_gp, p_covid_symp_higher_mscape_8, p_covid_symp_higher_mscape_26
           , layout = (9, 1)
           , size = (1200,3600)) #default size is 600,400

Plots.plot( p_covid_symp_lower_gp, p_covid_symp_base_gp, p_covid_symp_higher_gp
            , p_covid_symp_lower_mscape_8, p_covid_symp_base_mscape_8, p_covid_symp_higher_mscape_8
            , p_covid_symp_lower_mscape_26, p_covid_symp_base_mscape_26, p_covid_symp_higher_mscape_26
           , layout = (9, 1)
           , size = (1200,3600)) #default size is 600,400

Plots.plot(  p_covid_symp_lower_gp,        p_covid_symp_base_gp,        p_covid_symp_higher_gp
            ,p_covid_symp_lower_mscape_8,  p_covid_symp_base_mscape_8,  p_covid_symp_higher_mscape_8
            ,p_covid_symp_lower_mscape_26, p_covid_symp_base_mscape_26, p_covid_symp_higher_mscape_26
           , layout = (3, 3)
           , size = (3600,1200)) #default size is 600,400

# GP median / ICU (26 sites) median                                                                                                       
median(symp_tds_dict["lower_symp_gp_300s_winter"][ :,:GP_TD]) / median(symp_tds_dict["lower_symp_icu26_300s_winter"][ :,:ICU_TD]) # 1.110465755234229
median(symp_tds_dict["base_symp_gp_300s_winter"][  :,:GP_TD]) / median(symp_tds_dict["base_symp_icu26_300s_winter"][  :,:ICU_TD]) # 1.1304666762082662
median(symp_tds_dict["higher_symp_gp_300s_winter"][:,:GP_TD]) / median(symp_tds_dict["higher_symp_icu26_300s_winter"][:,:ICU_TD]) # 1.0865666483810124
# GP median - ICU (26 sites) median
median(symp_tds_dict["lower_symp_gp_300s_winter"][ :,:GP_TD]) - median(symp_tds_dict["lower_symp_icu26_300s_winter"][ :,:ICU_TD]) # 5.592409059038793
median(symp_tds_dict["base_symp_gp_300s_winter"][  :,:GP_TD]) - median(symp_tds_dict["base_symp_icu26_300s_winter"][  :,:ICU_TD]) # 5.826077168858603
median(symp_tds_dict["higher_symp_gp_300s_winter"][:,:GP_TD]) - median(symp_tds_dict["higher_symp_icu26_300s_winter"][:,:ICU_TD]) # 3.5569017875736577
# median( GP / ICU (26 sites) )
median(symp_tds_dict["lower_symp_gp_300s_winter"][ :,:GP_TD] ./ symp_tds_dict["lower_symp_icu26_300s_winter"][ :,:ICU_TD])        # 1.1074242852618692
median(symp_tds_dict["base_symp_gp_300s_winter"][  :,:GP_TD] ./ symp_tds_dict["base_symp_icu26_300s_winter"][  :,:ICU_TD])        # 1.1171268107896037
median(symp_tds_dict["higher_symp_gp_300s_winter"][:,:GP_TD] ./ symp_tds_dict["higher_symp_icu26_300s_winter"][:,:ICU_TD])        # 1.0892461380289606
# median( GP - ICU (26 sites) )
median(symp_tds_dict["lower_symp_gp_300s_winter"][ :,:GP_TD] .- symp_tds_dict["lower_symp_icu26_300s_winter"][ :,:ICU_TD])        # 5.328355384956517
median(symp_tds_dict["base_symp_gp_300s_winter"][  :,:GP_TD] .- symp_tds_dict["base_symp_icu26_300s_winter"][  :,:ICU_TD])        # 5.411715461946187
median(symp_tds_dict["higher_symp_gp_300s_winter"][:,:GP_TD] .- symp_tds_dict["higher_symp_icu26_300s_winter"][:,:ICU_TD])        # 3.839837121891307

# GP median / ICU (8 sites) median                                                                                                       
median(symp_tds_dict["lower_symp_gp_300s_winter"][ :,:GP_TD]) / median(symp_tds_dict["lower_symp_icu8_300s_winter"][ :,:ICU_TD]) # 1.0407802241754118
median(symp_tds_dict["base_symp_gp_300s_winter"][  :,:GP_TD]) / median(symp_tds_dict["base_symp_icu8_300s_winter"][  :,:ICU_TD]) # 1.0392978199851857
median(symp_tds_dict["higher_symp_gp_300s_winter"][:,:GP_TD]) / median(symp_tds_dict["higher_symp_icu8_300s_winter"][:,:ICU_TD]) # 1.0241189903186663
# GP median - ICU (8 sites) median
median(symp_tds_dict["lower_symp_gp_300s_winter"][ :,:GP_TD]) - median(symp_tds_dict["lower_symp_icu8_300s_winter"][ :,:ICU_TD]) # 2.202759139709613
median(symp_tds_dict["base_symp_gp_300s_winter"][  :,:GP_TD]) - median(symp_tds_dict["base_symp_icu8_300s_winter"][  :,:ICU_TD]) # 1.9088106162280525
median(symp_tds_dict["higher_symp_gp_300s_winter"][:,:GP_TD]) - median(symp_tds_dict["higher_symp_icu8_300s_winter"][:,:ICU_TD]) # 1.051444478708433
# median( GP / ICU (8 sites) )
median(symp_tds_dict["lower_symp_gp_300s_winter"][ :,:GP_TD] ./ symp_tds_dict["lower_symp_icu8_300s_winter"][ :,:ICU_TD])        # 1.0441019855063822
median(symp_tds_dict["base_symp_gp_300s_winter"][  :,:GP_TD] ./ symp_tds_dict["base_symp_icu8_300s_winter"][  :,:ICU_TD])        # 1.0309998649267524
median(symp_tds_dict["higher_symp_gp_300s_winter"][:,:GP_TD] ./ symp_tds_dict["higher_symp_icu8_300s_winter"][:,:ICU_TD])        # 1.0237583575663427
# median( GP - ICU (8 sites) )
median(symp_tds_dict["lower_symp_gp_300s_winter"][ :,:GP_TD] .- symp_tds_dict["lower_symp_icu8_300s_winter"][ :,:ICU_TD])        # 2.344140085255308
median(symp_tds_dict["base_symp_gp_300s_winter"][  :,:GP_TD] .- symp_tds_dict["base_symp_icu8_300s_winter"][  :,:ICU_TD])        # 1.5385031196979675
median(symp_tds_dict["higher_symp_gp_300s_winter"][:,:GP_TD] .- symp_tds_dict["higher_symp_icu8_300s_winter"][:,:ICU_TD])        # 1.0823450179759764

# It appears that increasing the symptomatic proportion reduces the absolute advantage of ICU (26 sites) sampling over GP sampling, but the relative advantage
# is maintained. For ICU (8 sites), increasing the symptomatic proportion reduces its advantage over GP sampling in both absolute and relative terms. However,
# could this be due to the saturation effect, where the sampling level is above the number of patients? Probably not at 300 samples per week and changing the
# symptomatic proportion shouldn't affect the number of hospitalisations or critical care admissions.

# Compute ratios of detection for individual simulations and plot those
lower_symp_gp_300s_winter_gp_icu26_ratios  = symp_tds_dict["lower_symp_gp_300s_winter" ][:,:GP_TD] ./ symp_tds_dict["lower_symp_icu26_300s_winter" ][:,:ICU_TD]
base_symp_gp_300s_winter_gp_icu26_ratios   = symp_tds_dict["base_symp_gp_300s_winter"  ][:,:GP_TD] ./ symp_tds_dict["base_symp_icu26_300s_winter"  ][:,:ICU_TD]
higher_symp_gp_300s_winter_gp_icu26_ratios = symp_tds_dict["higher_symp_gp_300s_winter"][:,:GP_TD] ./ symp_tds_dict["higher_symp_icu26_300s_winter"][:,:ICU_TD]

Plots.scatter(lower_symp_gp_300s_winter_gp_icu26_ratios)
Plots.scatter!(base_symp_gp_300s_winter_gp_icu26_ratios)
Plots.scatter!(higher_symp_gp_300s_winter_gp_icu26_ratios)

p_lower_gp_icu26_ratio  = kernel_box_jitter_plot(;x = collect( skipmissing( filter(x -> isfinite(x), lower_symp_gp_300s_winter_gp_icu26_ratios ))), plot_color = lower_symp_color,samp_strategy="RCGP TD / ICU (26 sites) TD (lower symptomatic) ",n_samples=300,n_sites="",n_blocks="",season="Winter",x_label="",x_min=0,x_max=7, annot_x_pos="middle") #, main_title = catchpop_keys[1])
p_base_gp_icu26_ratio   = kernel_box_jitter_plot(;x = collect( skipmissing( filter(x -> isfinite(x), base_symp_gp_300s_winter_gp_icu26_ratios ))), plot_color = base_symp_color,samp_strategy="RCGP TD / ICU (26 sites) TD (base symptomatic)",n_samples=300,n_sites="",n_blocks="",season="Winter",x_label="",x_min=0,x_max=7, annot_x_pos="middle") #, main_title = catchpop_keys[1])
p_higher_gp_icu26_ratio = kernel_box_jitter_plot(;x = collect( skipmissing( filter(x -> isfinite(x), higher_symp_gp_300s_winter_gp_icu26_ratios ))), plot_color = higher_symp_color,samp_strategy="RCGP TD / ICU (26 sites) TD (higher symptomatic)",n_samples=300,n_sites="",n_blocks="",season="Winter",x_label="Ratio of time to detection for sampling via RCGP / ICU (26 sites)",x_min=0,x_max=7, annot_x_pos="middle") #, main_title = catchpop_keys[1])

Plots.plot(  p_lower_gp_icu26_ratio, p_base_gp_icu26_ratio, p_higher_gp_icu26_ratio
           , layout = (3, 1)
           , size = (1600,1200)) #default size is 600,400

# Compute log of ratios of detection for individual simulations and plot those
lower_symp_gp_300s_winter_gp_icu26_ratios_log  = log.(symp_tds_dict["lower_symp_gp_300s_winter" ][:,:GP_TD] ./ symp_tds_dict["lower_symp_icu26_300s_winter" ][:,:ICU_TD])
base_symp_gp_300s_winter_gp_icu26_ratios_log   = log.(symp_tds_dict["base_symp_gp_300s_winter"  ][:,:GP_TD] ./ symp_tds_dict["base_symp_icu26_300s_winter"  ][:,:ICU_TD])
higher_symp_gp_300s_winter_gp_icu26_ratios_log = log.(symp_tds_dict["higher_symp_gp_300s_winter"][:,:GP_TD] ./ symp_tds_dict["higher_symp_icu26_300s_winter"][:,:ICU_TD])

p_lower_gp_icu26_ratio_log  = kernel_box_jitter_plot(;x = collect( skipmissing( filter(x -> isfinite(x), lower_symp_gp_300s_winter_gp_icu26_ratios_log  ))), plot_color = lower_symp_color, samp_strategy="RCGP TD / ICU (26 sites) TD (lower symptomatic)", n_samples=300,n_sites="",n_blocks="",season="Winter",x_label="",x_min=-2,x_max=2, annot_x_pos="left", dens_zero_below_zero = false) #, main_title = catchpop_keys[1])
p_base_gp_icu26_ratio_log   = kernel_box_jitter_plot(;x = collect( skipmissing( filter(x -> isfinite(x), base_symp_gp_300s_winter_gp_icu26_ratios_log   ))), plot_color = base_symp_color,  samp_strategy="RCGP TD / ICU (26 sites) TD (base symptomatic)",  n_samples=300,n_sites="",n_blocks="",season="Winter",x_label="",x_min=-2,x_max=2, annot_x_pos="left", dens_zero_below_zero = false) #, main_title = catchpop_keys[1])
p_higher_gp_icu26_ratio_log = kernel_box_jitter_plot(;x = collect( skipmissing( filter(x -> isfinite(x), higher_symp_gp_300s_winter_gp_icu26_ratios_log ))), plot_color = higher_symp_color,samp_strategy="RCGP TD / ICU (26 sites) TD (higher symptomatic)",n_samples=300,n_sites="",n_blocks="",season="Winter",x_label="Log of ratio of time to detection for sampling via RCGP / ICU (26 sites)",x_min=-2,x_max=2, annot_x_pos="left", dens_zero_below_zero = false) #, main_title = catchpop_keys[1])

Plots.plot(  p_lower_gp_icu26_ratio_log, p_base_gp_icu26_ratio_log, p_higher_gp_icu26_ratio_log
           , layout = (3, 1)
           , size = (1600,1200)) #default size is 600,400

# Compute absolute difference of detection times for individual simulations and plot those
lower_symp_gp_300s_winter_gp_icu26_diff  = symp_tds_dict["lower_symp_gp_300s_winter" ][:,:GP_TD] .- symp_tds_dict["lower_symp_icu26_300s_winter" ][:,:ICU_TD]
base_symp_gp_300s_winter_gp_icu26_diff   = symp_tds_dict["base_symp_gp_300s_winter"  ][:,:GP_TD] .- symp_tds_dict["base_symp_icu26_300s_winter"  ][:,:ICU_TD]
higher_symp_gp_300s_winter_gp_icu26_diff = symp_tds_dict["higher_symp_gp_300s_winter"][:,:GP_TD] .- symp_tds_dict["higher_symp_icu26_300s_winter"][:,:ICU_TD]

p_lower_gp_icu26_diff  = kernel_box_jitter_plot(;x = collect( skipmissing( filter(x -> isfinite(x), lower_symp_gp_300s_winter_gp_icu26_diff ))),  plot_color = lower_symp_color,samp_strategy="RCGP TD - ICU (26 sites) TD (lower symptomatic) ",n_samples=300,n_sites="",n_blocks="",season="Winter",x_label="",x_min=-30,x_max=60, annot_x_pos="left", dens_zero_below_zero = false) #, main_title = catchpop_keys[1])
p_base_gp_icu26_diff   = kernel_box_jitter_plot(;x = collect( skipmissing( filter(x -> isfinite(x), base_symp_gp_300s_winter_gp_icu26_diff ))),   plot_color = base_symp_color,samp_strategy="RCGP TD - ICU (26 sites) TD (base symptomatic)",n_samples=300,n_sites="",n_blocks="",season="Winter",x_label="",x_min=-30,x_max=60, annot_x_pos="left", dens_zero_below_zero = false) #, main_title = catchpop_keys[1])
p_higher_gp_icu26_diff = kernel_box_jitter_plot(;x = collect( skipmissing( filter(x -> isfinite(x), higher_symp_gp_300s_winter_gp_icu26_diff ))), plot_color = higher_symp_color,samp_strategy="RCGP TD - ICU (26 sites) TD (higher symptomatic)",n_samples=300,n_sites="",n_blocks="",season="Winter",x_label="Difference between time to detection\n via RCGP sampling and ICU (26 sites) sampling (days)",x_min=-30,x_max=60, annot_x_pos="left", dens_zero_below_zero = false) #, main_title = catchpop_keys[1])

lower_symp_gp_300s_winter_gp_icu26_diff_below_0 = filter(x -> x <0, lower_symp_gp_300s_winter_gp_icu26_diff)
length(lower_symp_gp_300s_winter_gp_icu26_diff_below_0) / length(lower_symp_gp_300s_winter_gp_icu26_diff) #0.301

base_symp_gp_300s_winter_gp_icu26_diff_below_0 = filter(x -> x <0, base_symp_gp_300s_winter_gp_icu26_diff)
length(base_symp_gp_300s_winter_gp_icu26_diff_below_0) / length(base_symp_gp_300s_winter_gp_icu26_diff) # 0.273

higher_symp_gp_300s_winter_gp_icu26_diff_below_0 = filter(x -> x <0, higher_symp_gp_300s_winter_gp_icu26_diff)
length(higher_symp_gp_300s_winter_gp_icu26_diff_below_0) / length(higher_symp_gp_300s_winter_gp_icu26_diff) # 0.323


Plots.plot(  p_lower_gp_icu26_diff, p_base_gp_icu26_diff, p_higher_gp_icu26_diff
           , layout = (3, 1)
           , size = (800,1200)) #default size is 600,400

# Final combined plot
Plots.plot(  p_lower_gp_icu26_diff,  p_lower_gp_icu26_ratio
            ,p_base_gp_icu26_diff,   p_base_gp_icu26_ratio
            ,p_higher_gp_icu26_diff, p_higher_gp_icu26_ratio
           , layout = (3, 2)
           , size = (1600,1300)) #default size is 600,400

Plots.savefig( joinpath( root_dir
                        , "NBPMscape/paper/analysis/2_sampling_analysis/covid_like_varied_sympt_prop"
                        ,"gp_icu26_td_comparison_varying_symp_prop.png"))


################################################################################################################################################
# (8) Statistical tests
################################################################################################################################################

# Analysis in (7) is inconclusive so further analysis on differences in time to detection (absolute differences and ratios) 
# with different symptomatic proportion including statistical tests performed in 
# 'NBPMscape\scripts\paper\2_sampling_analysis\6a_symp_prop_stat_test_inc_log_v2.R'
# This includes all symptomatic proportions investigated (base case, bc-20pp, bc+20pp, bc+40pp and 100% symptomatic).
# This script is used to produce tables and a figure in the supplementary material.