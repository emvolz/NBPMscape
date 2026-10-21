#=

Script simulating metagenomic surveillance sampling infections in different healthcare settings.
The outbreaks simulated are two versions of an influenza-like pathogen (mild and severe).
Sampling is simulated at primary and critical care (tertiary) healthcare settings, based on
Oxford-RCGP RSC and NHS-mSCAPE respectively.

This file contains:
(0) Definition of root directory
(1) Influenza-like MILD
    (1.1) Description of outbreak simulation
    (1.2) Check estimated epidemiological parameters from transmission chain 
          (those not directly input into model) before running full simulation:
          (1.2.1) Reproduction number resulting from model input epidemiological parameters
          (1.2.2) Generation time resulting from model input epidemiological parameters
    (1.3) Maniputlate simulation files
    (1.4) Compute times to detection
(2) Influenza-like SEVERE
    (2.1) Description of outbreak simulation
    (2.2) Check estimated epidemiological parameters from transmission chain 
          (those not directly input into model) before running full simulation:
          (2.2.1) Reproduction number resulting from model input epidemiological parameters
          (2.2.2) Generation time resulting from model input epidemiological parameters
    (2.3) Maniputlate simulation files
    (2.4) Compute times to detection

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
# (1) Influenza-like MILD
#     (1.1) Description of outbreak simulation
############################################################################################
#=
    Simulation run on Imperial College High Performance Computing (HPC) system
    using a Julia script and PBS file to manage the process:

    Scripts: Julia - 'influenzalike-mild-1.6.1_w_G_filtered_nrep1.jl'
             PBS   - 'influenzalike-mild-1.6.1_w_G_filtered_array100_nrep1.pbs'
            (see NBPMscape/paper/analysis/1_outbreak_simulations/influenza_like_mild/scripts)
    HPC run numbers: 4422974
    Simulation 'maxtime' = 80 days
    Simulation replicates per PBS array run: 1
    Runs in PBS array: 100
    Total sim reps completed: 100
    CPUs per array run: 1
    RAM per array run: 64GB
    Longest batch run time: ~3.2hrs (~3.2hrs per sim rep)
=#
############################################################################################
# (1) Influenza-like MILD
#     (1.2) Check estimated epidemiological parameters from transmission chain 
#           (those not directly input into model) before running full simulation:
#           (1.2.1) Reproduction number resulting from model input epidemiological parameters
############################################################################################

# Checking R for influenza MILD parameters and setting the infectivity value

# Load parameters
initialize_parameters("config/outbreak_params_influenza_like_mild_20261008.yaml");

#using gleam importations
# function simforest(p; initialtime=0.0, maxtime=30.0, maxgenerations::Int64=100, initialcontact=:H, importmodel=:GLEAM_outbreak_country_not_specified, max_cases=50000)

# Check current infectivity value, which is effectively a scaling parameter
NBPMscape.P.infectivity # 0.9

infectivitytoR(NBPMscape.P.infectivity; nsims = 1000) # Reproduction number
# Change first argument if necessary - Aiming for R~1.5
# Aiming for R = 1.46 as per Biggerstaff et al (2014), Estimates of the reproduction number for seasonal, pandemic, and zoonotic influenza: a systematic review of the literature, BMC Infectious Diseases, 14(1),DOI:10.1186/1471-2334-14-480
# Note that this is higher than estimated values of R0~1 in 'Epidemiological report of pandemic (H1N1) 2009 in the UK' published by Health Protection Agency https://www.qresearch.org/Public_Documents/Pan%20flu%20report_final_8October2010%20covered%5b2%5d.pdf, but a higher value seems more useful for the purpose of testing the model.
#infectivitytoR(0.9; nsims = 1000) # Reproduction number

# 1
#NGM = [0.554 0.296 0.301; 0.199 0.925 0.215; 0.589 0.459 0.819]
#1.4368935445411055
# 2
#NGM = [0.529 0.385 0.306; 0.259 0.914 0.187; 0.533 0.539 0.828]
#1.465268301588458
# 3
#NGM = [0.563 0.323 0.295; 0.217 0.958 0.215; 0.49 0.506 0.848]
#1.4639196115714939
# 4
#NGM = [0.566 0.331 0.354; 0.276 0.992 0.224; 0.532 0.54 0.915]
#1.5718328246823137
# 5
#NGM = [0.536 0.321 0.289; 0.301 0.987 0.176; 0.498 0.487 0.826]
#1.4619375502607461

############################################################################################
# (1) Influenza-like MILD
#     (1.2) Check estimated epidemiological parameters from transmission chain 
#           (those not directly input into model) before running full simulation:
#           (1.2.2) Generation time resulting from model input epidemiological parameters
############################################################################################

# Compute generation time to check parameters working as expected
NREPS = 5 #50 #500 #1000 #5000
sims = map( _-> simforest(NBPMscape.P; initialtime=0.0, maxtime=40.0, maxgenerations=100, max_cases=100000000, importmodel=:GLEAM_outbreak_country_not_specified), 1:NREPS )
println(sims[1].G.generation)
names(sims[1].D)

# Initialise vector to store time between donor and recipient being infected
generation_times = []

# Loop through recipients, record their time of infection and find the time of infection of those they infect. The difference is the individual generation time.
for j in 1:NREPS
    for i in 1:nrow(sims[j].D) #i=1
        donor_id = sims[j].D.recipient[i]
        donor_temp_tinf = sims[j].D.timetransmission[i]
        
        filtered_df = filter(row -> !ismissing(row.donor) && row.donor == donor_id, sims[j].D)

        gt_temp = filtered_df.timetransmission .- donor_temp_tinf
        
        generation_times = vcat(generation_times, gt_temp)

    end
end
mean(generation_times)     # 4.205102136607828
median(generation_times)   # 3.8974844747879587


############################################################################################
# (1) Influenza-like MILD
#     (1.3) Maniputlate simulation files
############################################################################################

## Load simulated outbreak data
main_folder = joinpath( root_dir, "NBPMscape/paper/analysis/1_outbreak_simulations")
sims = combine_sim_reps(;  # sim_input_folders in vector format where multiple paths
                          sim_input_folders = [ "$(main_folder)/influenza_like_mild/output/4422974_sims_G_filtered/" ]
                        , sim_object_name = "sims_G_filtered" #"sims" # "sims_G_icu_filter" # "sims_G_gp_filter" # "sim_G_filtered"
                        , nrep = 1
                        )
# Check df that shows errors in the combination process
sims_mismatch_df = sims[2]
# df is empty

# Refine sims
sims = sims[1]

# Save to file
@save joinpath( main_folder, "influenza_like_mild/output/influenzalike-mild-1.6.1-sims_G_filtered_nrep100_4422974.jld2" ) sims

# Reload file if necessary
#sims = load( joinpath( main_folder, "influenza_like_mild/output/influenzalike-mild-1.6.1-sims_G_filtered_nrep100_4422974.jld2" )
#            , "sims")


############################################################################################
# (1) Influenza-like MILD
#     (1.4) Compute times to detection
############################################################################################

# Define path to output folder for analyses files
output_folder = joinpath( root_dir, "NBPMscape/paper/analysis/2_sampling_analysis/influenza_like_mild/" )

## Compute times to detection with sampling in different healthcare settings
# using different model parameters representing different scenarios
# See 'NBPMscape/config/2026_02_18/config_file_record.xlsx' for summary descriptions of scenarios (and indicated below by numbers in for loops)

# Start with base parameters
initialize_parameters("config/outbreak_params_influenza_like_mild_20261008.yaml")

# Set sampling parameters for different Primary care sampling scenarios
gp_sampling_params = (  scenario_01 = (scenario_number =  1, gp_ari_consults = 327, gp_ari_swabs =  747, gp_swabs_mg =  300)
                      , scenario_13 = (scenario_number = 13, gp_ari_consults = 327, gp_ari_swabs =  747, gp_swabs_mg =  500)
                      , scenario_25 = (scenario_number = 25, gp_ari_consults = 327, gp_ari_swabs = 1000, gp_swabs_mg = 1000)
                      , scenario_37 = (scenario_number = 37, gp_ari_consults = 180, gp_ari_swabs =  319, gp_swabs_mg =  300)
                      , scenario_49 = (scenario_number = 49, gp_ari_consults = 180, gp_ari_swabs =  500, gp_swabs_mg =  500)
                      , scenario_61 = (scenario_number = 61, gp_ari_consults = 180, gp_ari_swabs = 1000, gp_swabs_mg = 1000)
                      )

# Run analysis for 6 different scenarios (2 seasons x 3 numbers of samples) using GPs (Primary care)
for i in 1:length(gp_sampling_params) # in [1,13,25,37,49,61] # i=1
    gp_tds = gp_td(  sims = sims
                   , gp_ari_consults = gp_sampling_params[i].gp_ari_consults
                   , gp_ari_swabs    = gp_sampling_params[i].gp_ari_swabs
                   , gp_swabs_mg     = gp_sampling_params[i].gp_swabs_mg
                   )
    CSV.write( joinpath( root_dir, "NBPMscape/paper/analysis/2_sampling_analysis/influenza_like_mild/gp_td"
                        ,"gp_tds_$(gp_sampling_params[i].scenario_number).csv" )
                , gp_tds )
    # Also, save file with 'GP_simid' column removed to reduce file size
    gp_tds_small = gp_tds[:,["sim_n","GP_TD","GP_3TD","n_GP_cases","n_GP_cases_sampled","n_GP_cases_sampled_positive","tinf_relating_to_GP_TD","last_tinf_relating_to_GP_3TD"]]
    CSV.write( joinpath( root_dir, "NBPMscape/paper/analysis/2_sampling_analysis/influenza_like_mild/gp_td_small"
                        ,"small_gp_tds_$(gp_sampling_params[i].scenario_number).csv" )
                , gp_tds_small )
end


# Set sampling parameters for different Critical Care mSCAPE sampling scenarios
mscape_icu_sampling_sites_8  = CSV.read( joinpath( root_dir, "NBPMscape/data/nhs_trust_site_sample_targets_24pilot.csv" ), DataFrame  )
mscape_icu_sampling_sites_26 = CSV.read( joinpath( root_dir, "NBPMscape/data/nhs_trust_site_sample_targets_24pilot_25cohort_longlist.csv" ), DataFrame  )

icu_sampling_params = ( # 8 sites
                          scenario_01 = (scenario_number =  1, n_icu_samples_per_week =  300, icu_ari_admissions =  1440, icu_nhs_trust_sampling_sites_file = mscape_icu_sampling_sites_8)
                        , scenario_13 = (scenario_number = 13, n_icu_samples_per_week =  500, icu_ari_admissions =  1440, icu_nhs_trust_sampling_sites_file = mscape_icu_sampling_sites_8)
                        , scenario_25 = (scenario_number = 25, n_icu_samples_per_week = 1000, icu_ari_admissions =  1440, icu_nhs_trust_sampling_sites_file = mscape_icu_sampling_sites_8)
                        , scenario_37 = (scenario_number = 37, n_icu_samples_per_week =  300, icu_ari_admissions =   793, icu_nhs_trust_sampling_sites_file = mscape_icu_sampling_sites_8)
                        , scenario_49 = (scenario_number = 49, n_icu_samples_per_week =  500, icu_ari_admissions =   793, icu_nhs_trust_sampling_sites_file = mscape_icu_sampling_sites_8)
                        , scenario_61 = (scenario_number = 61, n_icu_samples_per_week = 1000, icu_ari_admissions =   793, icu_nhs_trust_sampling_sites_file = mscape_icu_sampling_sites_8)
                        # 26 sites
                        , scenario_73 = (scenario_number = 73, n_icu_samples_per_week =  300, icu_ari_admissions =  1440, icu_nhs_trust_sampling_sites_file = mscape_icu_sampling_sites_26)
                        , scenario_74 = (scenario_number = 74, n_icu_samples_per_week =  500, icu_ari_admissions =  1440, icu_nhs_trust_sampling_sites_file = mscape_icu_sampling_sites_26)
                        , scenario_75 = (scenario_number = 75, n_icu_samples_per_week = 1000, icu_ari_admissions =  1440, icu_nhs_trust_sampling_sites_file = mscape_icu_sampling_sites_26)
                        , scenario_76 = (scenario_number = 76, n_icu_samples_per_week =  300, icu_ari_admissions =   793, icu_nhs_trust_sampling_sites_file = mscape_icu_sampling_sites_26)
                        , scenario_77 = (scenario_number = 77, n_icu_samples_per_week =  500, icu_ari_admissions =   793, icu_nhs_trust_sampling_sites_file = mscape_icu_sampling_sites_26)
                        , scenario_78 = (scenario_number = 78, n_icu_samples_per_week = 1000, icu_ari_admissions =   793, icu_nhs_trust_sampling_sites_file = mscape_icu_sampling_sites_26)
                    )

# Double check loaded parameters required for ICU sampling
NBPMscape.P.pathogen_type # "virus"
NBPMscape.P.icu_sample_type #"regional" 
NBPMscape.P.sample_icu_cases_version #"number"
NBPMscape.P.icu_only_sample_before_death #true            

# Run analysis for ICU sampling scenarios (2 seasons x 3 numbers of samples x 2 numbers of sampling sites)
for j in 1:length(icu_sampling_params) # [1,13,25,37,49,61,73,74,75,76,77,78] #j=1
    #initialize_parameters("config/influenza_v3/mild/outbreak_params_influenza_like_mild_$(j).yaml")
    
    icu_tds = icu_td( p = NBPMscape.P
                    , sims = sims
                    , n_icu_samples_per_week   = icu_sampling_params[j].n_icu_samples_per_week
                    , icu_ari_admissions       = icu_sampling_params[j].icu_ari_admissions
                    , nhs_trust_sampling_sites = icu_sampling_params[j].icu_nhs_trust_sampling_sites_file
                    )
    
    CSV.write( joinpath( output_folder,"icu_td/icu_tds_$(icu_sampling_params[j].scenario_number).csv" ), icu_tds )

    # Also, save file with 'ICU_simid' column removed to reduce file size
    icu_tds_small = icu_tds[:,["sim_n","ICU_TD","ICU_3TD","n_ICU_cases","n_ICU_cases_sampled","n_ICU_cases_sampled_positive","tinf_relating_to_ICU_TD","last_tinf_relating_to_ICU_3TD"]]
    CSV.write( joinpath( output_folder,"icu_td_small/small_icu_tds_$(icu_sampling_params[j].scenario_number).csv" ), icu_tds_small )
end

# Check that generating TD for all (or nearly all) simulations replicates
icu_tds_1 = CSV.read( joinpath(root_dir,"NBPMscape/paper/analysis/2_sampling_analysis/influenza_like_mild/icu_td_small/small_icu_tds_1.csv"), DataFrame)
sum(isfinite.(icu_tds_1[:,:ICU_TD]))  #
sum(isfinite.(icu_tds_1[:,:ICU_3TD])) #

gp_tds_1 = CSV.read( joinpath(root_dir,"NBPMscape/paper/analysis/2_sampling_analysis/influenza_like_mild/gp_td_small/small_gp_tds_1.csv"), DataFrame)
sum(isfinite.(gp_tds_1[:,:GP_TD]))  #
sum(isfinite.(gp_tds_1[:,:GP_3TD])) #

############################################################################################
# (2) Influenza-like SEVERE
#     (2.1) Description of outbreak simulation
############################################################################################
#=
    Simulation run on Imperial College High Performance Computing (HPC) system
    using a Julia script and PBS file to manage the process:

    Scripts: Julia - 'influenzalike-severe-1.6.1_w_G_filtered_nrep1.jl'
             PBS   - 'influenzalike-severe-1.6.1_w_G_filtered_array100_nrep1.pbs'
            (see NBPMscape/paper/analysis/1_outbreak_simulations/influenza_like_severe/scripts)
    HPC run numbers: 4327529
    Simulation 'maxtime' = 60 days
    Simulation replicates per PBS array run: 1
    Runs in PBS array: 100
    Total sim reps completed: 100
    CPUs per array run: 1
    RAM per array run: 64GB
    Longest batch run time: ~2.5hrs (~2.5hrs per sim rep)
=#
############################################################################################
# (2) Influenza-like SEVERE
#     (2.2) Check estimated epidemiological parameters from transmission chain 
#           (those not directly input into model) before running full simulation:
#           (2.2.1) Reproduction number resulting from model input epidemiological parameters
############################################################################################

# Checking R for influenza SEVERE parameters and setting the infectivity value

# Load parameters
initialize_parameters("config/outbreak_params_influenza_like_severe_20261008.yaml");

#using gleam importations
# function simforest(p; initialtime=0.0, maxtime=30.0, maxgenerations::Int64=100, initialcontact=:H, importmodel=:GLEAM_outbreak_country_not_specified, max_cases=50000)

# Check current infectivity value, which is effectively a scaling parameter
NBPMscape.P.infectivity # 1.25

infectivitytoR(NBPMscape.P.infectivity; nsims = 1000) # Reproduction number
# Change first argument if necessary - Aiming for R~2
#infectivitytoR(1.25; nsims = 1000) # Reproduction number

# 1
#NGM = [0.773 0.403 0.458; 0.36 1.308 0.275; 0.698 0.678 1.134]
#2.015486624767621
# 2
#NGM = [0.754 0.448 0.396; 0.35 1.238 0.29; 0.69 0.648 1.201]
#1.9940001682166724
# 3
#NGM = [0.744 0.406 0.424; 0.369 1.271 0.258; 0.711 0.649 1.183]
#1.989640530940663
# 4
#NGM = [0.819 0.408 0.43; 0.375 1.208 0.26; 0.786 0.679 1.278]
#2.0563917389049977
# 5
#NGM = [0.736 0.397 0.432; 0.347 1.234 0.267; 0.771 0.695 1.122]
#1.9734585329621916

############################################################################################
# (2) Influenza-like SEVERE
#     (2.2) Check estimated epidemiological parameters from transmission chain 
#           (those not directly input into model) before running full simulation:
#           (2.2.2) Generation time resulting from model input epidemiological parameters
############################################################################################

# Compute generation time to check parameters working as expected

# Load parameters
initialize_parameters("config/outbreak_params_influenza_like_severe_20261008.yaml");

NREPS = 5 #50 #500 #1000 #5000
sims = map( _-> simforest(NBPMscape.P; initialtime=0.0, maxtime=40.0, maxgenerations=100, max_cases=100000000, importmodel=:GLEAM_outbreak_country_not_specified), 1:NREPS )
println(sims[1].G.generation)
names(sims[1].D)

# Initialise vector to store time between donor and recipient being infected
generation_times = []

# Loop through recipients, record their time of infection and find the time of infection of those they infect. The difference is the individual generation time.
for j in 1:NREPS
    for i in 1:nrow(sims[j].D) #i=1
        donor_id = sims[j].D.recipient[i]
        donor_temp_tinf = sims[j].D.timetransmission[i]
        
        filtered_df = filter(row -> !ismissing(row.donor) && row.donor == donor_id, sims[j].D)

        gt_temp = filtered_df.timetransmission .- donor_temp_tinf
        
        generation_times = vcat(generation_times, gt_temp)

    end
end
mean(generation_times)     # 3.9802684414968765
median(generation_times)   # 3.69142621484945


############################################################################################
# (2) Influenza-like SEVERE
#     (2.3) Maniputlate simulation files
############################################################################################

## Load simulated outbreak data
main_folder = joinpath( root_dir, "NBPMscape/paper/analysis/1_outbreak_simulations")
sims = combine_sim_reps(;  # sim_input_folders in vector format where multiple paths
                          sim_input_folders = [ "$(main_folder)/influenza_like_severe/output/4327529_sims_G_filtered/" ]
                        , sim_object_name = "sims_G_filtered" #"sims" # "sims_G_icu_filter" # "sims_G_gp_filter" # "sim_G_filtered"
                        , nrep = 1
                        )
# Check df that shows errors in the combination process
sims_mismatch_df = sims[2]
# df is empty

# Refine sims
sims = sims[1]

# Save to file
@save joinpath( main_folder, "influenza_like_severe/output/influenzalike-severe-1.6.1-sims_G_filtered_nrep100_4327529.jld2" ) sims

# Reload file if necessary
#sims = load( joinpath( main_folder, "influenza_like_severe/output/influenzalike-severe-1.6.1-sims_G_filtered_nrep100_4327529.jld2" )
#            , "sims")


############################################################################################
# (2) Influenza-like SEVERE
#     (2.4) Compute times to detection
############################################################################################

# Define path to output folder for analyses files
output_folder = joinpath( root_dir, "NBPMscape/paper/analysis/2_sampling_analysis/influenza_like_severe/" )

## Compute times to detection for with sampling in different healthcare settings
# using different model parameters representing different scenarios
# See 'NBPMscape/config/2026_02_18/config_file_record.xlsx' for summary descriptions of scenarios (and indicated below by numbers in for loops)

# Start with base parameters
initialize_parameters("config/outbreak_params_influenza_like_severe_20261008.yaml")

# Set sampling parameters for different Primary care sampling scenarios
gp_sampling_params = (  scenario_01 = (scenario_number =  1, gp_ari_consults = 327, gp_ari_swabs =  747, gp_swabs_mg =  300)
                      , scenario_13 = (scenario_number = 13, gp_ari_consults = 327, gp_ari_swabs =  747, gp_swabs_mg =  500)
                      , scenario_25 = (scenario_number = 25, gp_ari_consults = 327, gp_ari_swabs = 1000, gp_swabs_mg = 1000)
                      , scenario_37 = (scenario_number = 37, gp_ari_consults = 180, gp_ari_swabs =  319, gp_swabs_mg =  300)
                      , scenario_49 = (scenario_number = 49, gp_ari_consults = 180, gp_ari_swabs =  500, gp_swabs_mg =  500)
                      , scenario_61 = (scenario_number = 61, gp_ari_consults = 180, gp_ari_swabs = 1000, gp_swabs_mg = 1000)
                      )

# Run analysis for 6 different scenarios (2 seasons x 3 numbers of samples) using GPs (Primary care)
for i in 1:length(gp_sampling_params) # in [1,13,25,37,49,61] # i=1
    gp_tds = gp_td(  sims = sims
                   , gp_ari_consults = gp_sampling_params[i].gp_ari_consults
                   , gp_ari_swabs    = gp_sampling_params[i].gp_ari_swabs
                   , gp_swabs_mg     = gp_sampling_params[i].gp_swabs_mg
                   )
    CSV.write( joinpath( root_dir, "NBPMscape/paper/analysis/2_sampling_analysis/influenza_like_severe/gp_td"
                        ,"gp_tds_$(gp_sampling_params[i].scenario_number).csv")
                , gp_tds, bufsize = 16 * 1024 * 1024 )
    # Also, save file with 'GP_simid' column removed to reduce file size
    gp_tds_small = gp_tds[:,["sim_n","GP_TD","GP_3TD","n_GP_cases","n_GP_cases_sampled","n_GP_cases_sampled_positive","tinf_relating_to_GP_TD","last_tinf_relating_to_GP_3TD"]]
    CSV.write( joinpath( root_dir, "NBPMscape/paper/analysis/2_sampling_analysis/influenza_like_severe/gp_td_small"
                        ,"small_gp_tds_$(gp_sampling_params[i].scenario_number).csv" )
                , gp_tds_small )
end

# Set sampling parameters for different Critical Care mSCAPE sampling scenarios
mscape_icu_sampling_sites_8  = CSV.read( joinpath( root_dir, "NBPMscape/data/nhs_trust_site_sample_targets_24pilot.csv" ), DataFrame  )
mscape_icu_sampling_sites_26 = CSV.read( joinpath( root_dir, "NBPMscape/data/nhs_trust_site_sample_targets_24pilot_25cohort_longlist.csv" ), DataFrame  )

icu_sampling_params = ( # 8 sites
                          scenario_01 = (scenario_number =  1, n_icu_samples_per_week =  300, icu_ari_admissions =  1440, icu_nhs_trust_sampling_sites_file = mscape_icu_sampling_sites_8)
                        , scenario_13 = (scenario_number = 13, n_icu_samples_per_week =  500, icu_ari_admissions =  1440, icu_nhs_trust_sampling_sites_file = mscape_icu_sampling_sites_8)
                        , scenario_25 = (scenario_number = 25, n_icu_samples_per_week = 1000, icu_ari_admissions =  1440, icu_nhs_trust_sampling_sites_file = mscape_icu_sampling_sites_8)
                        , scenario_37 = (scenario_number = 37, n_icu_samples_per_week =  300, icu_ari_admissions =   793, icu_nhs_trust_sampling_sites_file = mscape_icu_sampling_sites_8)
                        , scenario_49 = (scenario_number = 49, n_icu_samples_per_week =  500, icu_ari_admissions =   793, icu_nhs_trust_sampling_sites_file = mscape_icu_sampling_sites_8)
                        , scenario_61 = (scenario_number = 61, n_icu_samples_per_week = 1000, icu_ari_admissions =   793, icu_nhs_trust_sampling_sites_file = mscape_icu_sampling_sites_8)
                        # 26 sites
                        , scenario_73 = (scenario_number = 73, n_icu_samples_per_week =  300, icu_ari_admissions =  1440, icu_nhs_trust_sampling_sites_file = mscape_icu_sampling_sites_26)
                        , scenario_74 = (scenario_number = 74, n_icu_samples_per_week =  500, icu_ari_admissions =  1440, icu_nhs_trust_sampling_sites_file = mscape_icu_sampling_sites_26)
                        , scenario_75 = (scenario_number = 75, n_icu_samples_per_week = 1000, icu_ari_admissions =  1440, icu_nhs_trust_sampling_sites_file = mscape_icu_sampling_sites_26)
                        , scenario_76 = (scenario_number = 76, n_icu_samples_per_week =  300, icu_ari_admissions =   793, icu_nhs_trust_sampling_sites_file = mscape_icu_sampling_sites_26)
                        , scenario_77 = (scenario_number = 77, n_icu_samples_per_week =  500, icu_ari_admissions =   793, icu_nhs_trust_sampling_sites_file = mscape_icu_sampling_sites_26)
                        , scenario_78 = (scenario_number = 78, n_icu_samples_per_week = 1000, icu_ari_admissions =   793, icu_nhs_trust_sampling_sites_file = mscape_icu_sampling_sites_26)
                    )

# Double check loaded parameters required for ICU sampling
NBPMscape.P.pathogen_type # "virus"
NBPMscape.P.icu_sample_type #"regional" 
NBPMscape.P.sample_icu_cases_version #"number"
NBPMscape.P.icu_only_sample_before_death #true            

# Run analysis for ICU sampling scenarios (2 seasons x 3 numbers of samples x 2 numbers of sampling sites)
for j in 1:length(icu_sampling_params) # [1,13,25,37,49,61,73,74,75,76,77,78] #j=1
    #initialize_parameters("config/influenza_v3/mild/outbreak_params_influenza_like_mild_$(j).yaml")
    
    icu_tds = icu_td( p = NBPMscape.P
                    , sims = sims
                    , n_icu_samples_per_week   = icu_sampling_params[j].n_icu_samples_per_week
                    , icu_ari_admissions       = icu_sampling_params[j].icu_ari_admissions
                    , nhs_trust_sampling_sites = icu_sampling_params[j].icu_nhs_trust_sampling_sites_file
                    )
    
    CSV.write( joinpath( root_dir, "NBPMscape/paper/analysis/2_sampling_analysis/influenza_like_severe"
                        ,"icu_td/icu_tds_$(icu_sampling_params[j].scenario_number).csv" ), icu_tds )

    # Also, save file with 'ICU_simid' column removed to reduce file size
    icu_tds_small = icu_tds[:,["sim_n","ICU_TD","ICU_3TD","n_ICU_cases","n_ICU_cases_sampled","n_ICU_cases_sampled_positive","tinf_relating_to_ICU_TD","last_tinf_relating_to_ICU_3TD"]]
    CSV.write( joinpath( root_dir, "NBPMscape/paper/analysis/2_sampling_analysis/influenza_like_severe"
                        ,"icu_td_small/small_icu_tds_$(icu_sampling_params[j].scenario_number).csv" ), icu_tds_small )
end

# Check that generating TD for all (or nearly all) simulations replicates
icu_tds_1 = CSV.read( joinpath(root_dir,"NBPMscape/paper/analysis/2_sampling_analysis/influenza_like_severe/icu_td_small/small_icu_tds_1.csv"), DataFrame)
sum(isfinite.(icu_tds_1[:,:ICU_TD]))  #100
sum(isfinite.(icu_tds_1[:,:ICU_3TD])) #100

gp_tds_1 = CSV.read( joinpath(root_dir,"NBPMscape/paper/analysis/2_sampling_analysis/influenza_like_severe/gp_td_small/small_gp_tds_1.csv"), DataFrame)
sum(isfinite.(gp_tds_1[:,:GP_TD]))  #100
sum(isfinite.(gp_tds_1[:,:GP_3TD])) #100
