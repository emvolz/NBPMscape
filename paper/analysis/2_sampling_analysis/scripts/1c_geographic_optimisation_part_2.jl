#=
Part 2 of a series of scripts for the investigation of the effect on surveillance performance under
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
using Distributions
using GLM
using JLD2
using NBPMscape # Should be version from commit "2762757a687c64146016e6bc742d258e20525cb8" as at 17:00 16 Oct 2026
using Plots
using StatsBase
using StatsPlots

###############################################################################################################
# (4) Sampling infections to simulate TDs using different sampling strategies:
#    (4.1) HCGS
#        (4.1.3) Perform hierarchical clustering to group regions together based on TD distribution correlations
#                and allocate a number of sites to each region. Regions with lower median TD are prioritised and
#                sites are distributed across both groups and regions.
###############################################################################################################

#### After running 'site_selection_2026_03.R' to use block model and hierarchical clustering
#### to allocate numbers of sites to ITL2 regions the code below can be run

# Define root directory
#root_dir = root_folder() # Function defined in NBPMscape/src/misc_functions
dir = @__DIR__
while !isdir(joinpath(dir, "NBPMscape"))
    parent = dirname(dir)
    parent == dir && error("NBPmscape not found")
    dir = parent
end
root_dir = dir


geo_optim_folder = joinpath( root_dir, "NBPMscape/paper/analysis/2_sampling_analysis/covid_like/geo_optimisation")
site_selection_folder = joinpath( geo_optim_folder, "site_select_outputs" )

# Read in site numbers for each region with different combinations of block numbers and total site numbers,

# Define numbers of blocks used in hierarchical clustering in 'site_selection_2026_03.R'
# and numbers of sites allocated across regions
n_blocks_vals = [5, 10, 20, 30]
n_sites_vals  = [8, 10, 20, 26, 30, 40, 50, 60]

# Build list of filenames from combinations of n_blocks_vals and n_sites_vals
region_site_numbers_files = [
    "site_select_output_$(nb)blocks_$(ns)sites.csv"
    for nb in n_blocks_vals, ns in n_sites_vals
] |> vec

# Dynamically load each CSV file into a Dict, keyed by the filename prefix (without extension)
region_site_numbers_dict = Dict(
    splitext(f)[1] => CSV.read(joinpath(site_selection_folder, f), DataFrame)
    for f in region_site_numbers_files
)

##################################################################################################
# (4) Sampling infections to simulate TDs using different sampling strategies:
#     (4.1) HCGS
#         (4.1.4) Select sampling sites
##################################################################################################

sampling_site_file_output_folder = joinpath( geo_optim_folder, "selected_sites" )
if !isdir(sampling_site_file_output_folder)
    mkpath(sampling_site_file_output_folder)
end

# Define list of output filenames corresponding to region_site_numbers_csv_files list
sampling_site_files = [replace(f, "site_select_output_" => "sampling_sites_") for f in region_site_numbers_files]

# Loop through list of files containing the number of sites per region as allocated using hierarchical clustering
# and select sampling sites
for i in 1:length(region_site_numbers_files)
    select_nhs_trust_sampling_sites(; prioritise_adult_sampling = true
                                    , n_sites_max = 127
                                    , site_selection_folder = site_selection_folder
                                    , site_selection_file = region_site_numbers_files[i] #i=1
                                    , output_folder = sampling_site_file_output_folder
                                    , output_file = sampling_site_files[i] )
end


##################################################################################################
#(4) Sampling infections to simulate TDs using different sampling strategies:
#    (4.1) HCGS
#        (4.1.5) Simulate sampling and TD
##################################################################################################

# 32 different sets of sampling sites. 
# Run ICU sampling with 300, 500 and 1000 samples per week in winter and summer. 
# Generate plots of TD distribution as well as cumulative count metrics: infections at TD, fatal infections at TD

# Load parameters used in simulation
initialize_parameters("config/outbreak_params_covid19_like_20261008.yaml")

# Double check loaded parameters required for ICU sampling
NBPMscape.P.pathogen_type # "virus"
NBPMscape.P.icu_sample_type #"regional" 
NBPMscape.P.sample_icu_cases_version #"number"
NBPMscape.P.icu_only_sample_before_death #true            
NBPMscape.P.hosp_ari_admissions #27894

# Define variable parameters (192 sets = 32 * 3 * 2)
# 32 sampling_site_files 
samples_per_week = [300,500,1000]
seasonal_icu_ari_admissions = [1440,793]


# Reload simulation data file if necessary
#sims = load( normpath( joinpath( root_dir, "NBPMscape/paper/analysis/1_outbreak_simulations/covid_like/output/covidlike-1.6.1-sims_G_filtered_nrep1000_4300448.jld2") )
#            , "sims")

# Loop through variables and run ICU sampling to compute times to detection
# ~1 to 15 mins per loop on laptop (10 sites are quicker than 30 sites)
counter = 0
for i in 1:length(sampling_site_files) #i=1
    for j in 1:length(seasonal_icu_ari_admissions) #j=1
        for k in 1:length(samples_per_week) # k=2
            
            # Read in NHS Trust sampling sites
            nhs_trust_sampling_sites = CSV.read( joinpath( geo_optim_folder, "selected_sites", sampling_site_files[i] ), DataFrame  )

            # Compute times to detection by sampling ICU cases
            temp_icu_td = icu_td( p = NBPMscape.P, sims = sims
                                , n_icu_samples_per_week = samples_per_week[k]
                                , icu_ari_admissions = seasonal_icu_ari_admissions[j]
                                , nhs_trust_sampling_sites = nhs_trust_sampling_sites
                                )
            
            # Save output to file
            output_file = "icu_tds_seasonal_ari_$(seasonal_icu_ari_admissions[j])_$(samples_per_week[k])_$(sampling_site_files[i])"
            CSV.write( joinpath( geo_optim_folder, "icu_td/hcgs", output_file ), temp_icu_td)
            #temp_icu_td_k1 = CSV.read( joinpath( geo_optim_folder, "icu_td", "icu_tds_seasonal_ari_$(seasonal_icu_ari_admissions[j])_$(samples_per_week[1])_$(sampling_site_files[i]).csv" ), DataFrame )
            
            # Also save a smaller file as easier to manipulate
            # (remove column "ICU_simid")
            temp_icu_td_small = temp_icu_td[:,["sim_n", "ICU_TD", "ICU_3TD", "n_ICU_cases", "n_ICU_cases_sampled", "n_ICU_cases_sampled_positive", "tinf_relating_to_ICU_TD", "last_tinf_relating_to_ICU_3TD"]]
            output_file_small = "small_icu_tds_seasonal_ari_$(seasonal_icu_ari_admissions[j])_$(samples_per_week[k])_$(sampling_site_files[i])"
            CSV.write( joinpath( geo_optim_folder, "icu_td/hcgs_small", output_file_small ), temp_icu_td_small)
            counter = counter + 1
            println("Successfully completed $(counter) of $(length(sampling_site_files)*length(seasonal_icu_ari_admissions)*length(samples_per_week)) runs")
        end
    end
end

# Reset parameters
initialize_parameters("config/outbreak_params_covid19_like_20261008.yaml")

### Compute median TDs for different sampling parameter sets
# Read in saved data
path_name = joinpath( geo_optim_folder, "icu_td" )

# Prepare a container to keep the DataFrames
icu_tds_df_dict = Dict{String, DataFrame}()

# Loop through scenarios, read in files and compute median TD values
for i in 1:length(sampling_site_files) #i=1
    for j in 1:length(seasonal_icu_ari_admissions) #j=1
        for k in 1:length(samples_per_week) # k=1

            # Define parameter name to be used in dictionary. 
            # This is also the same as the file prefix as used above
            sampling_sites_prefix = sampling_site_files[i][1:findfirst(==('.'), sampling_site_files[i])-1]
            params_name = "small_icu_tds_seasonal_ari_$(seasonal_icu_ari_admissions[j])_$(samples_per_week[k])_$(sampling_sites_prefix)"
            file_icu = joinpath(path_name, "$(params_name).csv")
    
            # Read CSV into a DataFrame and store in the dictionary with the scenario number (numeric suffix) as key
            icu_tds_df_dict[params_name] = CSV.read(file_icu, DataFrame)[:,[:ICU_TD,:ICU_3TD]]
        end
    end
end


##################################################################################################
# (4) Sampling infections to simulate TDs using different sampling strategies:
#     (4.1) HCGS
#         (4.1.6) Add outbreak metrics relating to TDs
##################################################################################################

# Reload if necessary (produced in (3) above)
# Reload outbreak metrics dictionary if necessary
outbreak_metrics_time_series = load( joinpath( root_dir, "NBPMscape/paper/analysis/1_outbreak_simulations/covid_like/output"
                                              , "covidlike-1.6.1_nrep1000_4300448_outbreak_metrics_time_series.jld2") 
                                    , "outbreak_metrics_time_series" )

# Loop through dictionary and add columns to each df in dictionary for the cumulative outbreak count metrics
for key in keys(icu_tds_df_dict)
    
    df = icu_tds_df_dict[key] # df = icu_tds_df_dict["icu_tds_seasonal_ari_1440_500_sampling_sites_30blocks_10sites"]
    
    # Add new columns
    empty_col_vec = Vector{Union{Missing, Int}}(missing, nrow(df))
    # TD
    df.cum_tinf_td       = copy(empty_col_vec)
    df.cum_tgp_td        = copy(empty_col_vec)
    df.cum_ted_td        = copy(empty_col_vec)
    df.cum_thosp_td      = copy(empty_col_vec)
    df.cum_ticu_td       = copy(empty_col_vec)
    df.cum_tdeceased_td  = copy(empty_col_vec)
    df.cum_tinf_fatal_td = copy(empty_col_vec)
    df.cum_timport_td    = copy(empty_col_vec)
    # TD3
    df.cum_tinf_td3       = copy(empty_col_vec)
    df.cum_tgp_td3        = copy(empty_col_vec)
    df.cum_ted_td3        = copy(empty_col_vec)
    df.cum_thosp_td3      = copy(empty_col_vec)
    df.cum_ticu_td3       = copy(empty_col_vec)
    df.cum_tdeceased_td3  = copy(empty_col_vec)
    df.cum_tinf_fatal_td3 = copy(empty_col_vec)
    df.cum_timport_td3    = copy(empty_col_vec)

    for i in 1:1000 #i=1
        td  = df[i,:ICU_TD]
        td3 = df[i,:ICU_3TD]
        # Counts at TD
        df.cum_tinf_td[i]       = sum(outbreak_metrics_time_series.tinf_vals[i]       .<= td )
        df.cum_tgp_td[i]        = sum(outbreak_metrics_time_series.tgp_vals[i]        .<= td )
        df.cum_ted_td[i]        = sum(outbreak_metrics_time_series.ted_vals[i]        .<= td )
        df.cum_thosp_td[i]      = sum(outbreak_metrics_time_series.thosp_vals[i]      .<= td )
        df.cum_ticu_td[i]       = sum(outbreak_metrics_time_series.ticu_vals[i]       .<= td )
        df.cum_tdeceased_td[i]  = sum(outbreak_metrics_time_series.tdeceased_vals[i]  .<= td )
        df.cum_tinf_fatal_td[i] = sum(outbreak_metrics_time_series.tinf_fatal_vals[i] .<= td )
        df.cum_timport_td[i]    = sum(outbreak_metrics_time_series.timport_vals[i]    .<= td )
        # Counts at TD3
        df.cum_tinf_td3[i]       = sum(outbreak_metrics_time_series.tinf_vals[i]       .<= td3 )
        df.cum_tgp_td3[i]        = sum(outbreak_metrics_time_series.tgp_vals[i]        .<= td3 )
        df.cum_ted_td3[i]        = sum(outbreak_metrics_time_series.ted_vals[i]        .<= td3 )
        df.cum_thosp_td3[i]      = sum(outbreak_metrics_time_series.thosp_vals[i]      .<= td3 )
        df.cum_ticu_td3[i]       = sum(outbreak_metrics_time_series.ticu_vals[i]       .<= td3 )
        df.cum_tdeceased_td3[i]  = sum(outbreak_metrics_time_series.tdeceased_vals[i]  .<= td3 )
        df.cum_tinf_fatal_td3[i] = sum(outbreak_metrics_time_series.tinf_fatal_vals[i] .<= td3 )
        df.cum_timport_td3[i]    = sum(outbreak_metrics_time_series.timport_vals[i]    .<= td3 )
    end
    icu_tds_df_dict[key] = df
end

# Save as file
@save joinpath( geo_optim_folder, "covidlike-1.6.1_nrep1000_4300448_icu_tds_df_dict.jld2") icu_tds_df_dict

# Reload if necessary
#icu_tds_df_dict = load( joinpath( geo_optim_folder, "covidlike-1.6.1_nrep1000_4300448_icu_tds_df_dict.jld2"), "icu_tds_df_dict" )

# Reduce Julia memory
sims = []
sims_G_filtered = []



##################################################################################################
#(4) Sampling infections to simulate TDs using different sampling strategies:
#    (4.1) HCGS
#        (4.1.7) Plot results
##################################################################################################

# See 'NBPMscape/paper/analysis/2_sampling_analysis/scripts/supp_mat_icu_plots.jl' for plotting

###########################################################################################################
# (4) Sampling infections to simulate TDs using different sampling strategies:
#     (4.2) London only
#         (4.2.1) Select sampling sites
###########################################################################################################

### Method of allocating samples/sites in London only
## NHS Trusts in London (source: 'All acutes' tab in 202506-July-2025-beds-sitrep-data-finalversion-v2.xlsx from https://www.england.nhs.uk/statistics/statistical-work-areas/uec-sitrep/)
# Except for RP6	Moorfields Eye Hospital NHS Foundation Trust, all have ICU beds (adult and/or paediatric)
# Region	ICB Name	                                    Org Code	Org Name
# LONDON	NHS NORTH CENTRAL LONDON INTEGRATED CARE BOARD	RAL	        Royal Free London NHS Foundation Trust
# LONDON	NHS NORTH CENTRAL LONDON INTEGRATED CARE BOARD	RAN	        Royal National Orthopaedic Hospital NHS Trust
# LONDON	NHS NORTH CENTRAL LONDON INTEGRATED CARE BOARD	RKE	        Whittington Health NHS Trust
# LONDON	NHS NORTH CENTRAL LONDON INTEGRATED CARE BOARD	RP4	        Great Ormond Street Hospital for Children NHS Foundation Trust
# LONDON	NHS NORTH CENTRAL LONDON INTEGRATED CARE BOARD	RP6	        Moorfields Eye Hospital NHS Foundation Trust
# LONDON	NHS NORTH CENTRAL LONDON INTEGRATED CARE BOARD	RRV	        University College London Hospitals NHS Foundation Trust
# LONDON	NHS NORTH EAST LONDON INTEGRATED CARE BOARD	    R1H	        Barts Health NHS Trust
# LONDON	NHS NORTH EAST LONDON INTEGRATED CARE BOARD	    RF4	        Barking, Havering and Redbridge University Hospitals NHS Trust
# LONDON	NHS NORTH EAST LONDON INTEGRATED CARE BOARD	    RQX	        Homerton Healthcare NHS Foundation Trust
# LONDON	NHS NORTH WEST LONDON INTEGRATED CARE BOARD	    R1K	        London North West University Healthcare NHS Trust
# LONDON	NHS NORTH WEST LONDON INTEGRATED CARE BOARD	    RAS	        The Hillingdon Hospitals NHS Foundation Trust
# LONDON	NHS NORTH WEST LONDON INTEGRATED CARE BOARD	    RQM	        Chelsea and Westminster Hospital NHS Foundation Trust
# LONDON	NHS NORTH WEST LONDON INTEGRATED CARE BOARD	    RYJ	        Imperial College Healthcare NHS Trust
# LONDON	NHS SOUTH EAST LONDON INTEGRATED CARE BOARD	    RJ1	        Guy's and St Thomas' NHS Foundation Trust
# LONDON	NHS SOUTH EAST LONDON INTEGRATED CARE BOARD	    RJ2	        Lewisham and Greenwich NHS Trust
# LONDON	NHS SOUTH EAST LONDON INTEGRATED CARE BOARD	    RJZ	        King's College Hospital NHS Foundation Trust
# LONDON	NHS SOUTH WEST LONDON INTEGRATED CARE BOARD	    RAX	        Kingston Hospital NHS Foundation Trust
# LONDON	NHS SOUTH WEST LONDON INTEGRATED CARE BOARD	    RJ6	        Croydon Health Services NHS Trust
# LONDON	NHS SOUTH WEST LONDON INTEGRATED CARE BOARD	    RJ7	        St George's University Hospitals NHS Foundation Trust
# LONDON	NHS SOUTH WEST LONDON INTEGRATED CARE BOARD	    RPY	        The Royal Marsden NHS Foundation Trust
# LONDON	NHS SOUTH WEST LONDON INTEGRATED CARE BOARD	    RVR	        Epsom and St Helier University Hospitals NHS Trust

nhs_trusts_london = ["RAL","RAN","RKE","RP4","RP6","RRV","R1H","RF4","RQX","R1K","RAS","RQM","RYJ","RJ1","RJ2","RJZ","RAX","RJ6","RJ7","RPY","RVR"]
# Remove NHS Trusts unlikely to admit ARI into ICU
# LONDON	NHS NORTH CENTRAL LONDON INTEGRATED CARE BOARD	RP6	        Moorfields Eye Hospital NHS Foundation Trust
# LONDON	NHS NORTH CENTRAL LONDON INTEGRATED CARE BOARD	RAN	        Royal National Orthopaedic Hospital NHS Trust
nhs_trusts_london_non_ari = ["RAN","RP6"]
nhs_trusts_london_ari = [x for x in nhs_trusts_london if x ∉ nhs_trusts_london_non_ari]

ari_cc_bed_london = filter( row -> row.NHS_Trust_code in nhs_trusts_london_ari, ARI_CC_BED_SITREP)
ari_cc_bed_london_sorted = sort(ari_cc_bed_london, :Adult_critical_care_beds_available, rev=true)

# Rank by number of adult ICU beds adult and add largest
# However, RP4 Great Ormond Street Hospital for Children NHS Foundation Trust has the largest number of children's ICU beds
# but no adult beds so an exception is made for this Trust

# 8 sites
sampling_sites_london_8sites = ari_cc_bed_london_sorted[1:7,:]
append!( sampling_sites_london_8sites, filter( row -> row.NHS_Trust_code == "RP4" , ari_cc_bed_london_sorted) )
sampling_sites_london_8sites = sampling_sites_london_8sites[:, [ :NHS_Trust_code,:NHS_Trust_name
                                                                ,:Adult_critical_care_beds_available,:Paediatric_intensive_care_beds_available] ]
# Check how many of the largest paediatric ICUs are including when rank by adult ICU size
sort(ari_cc_bed_london, :Paediatric_intensive_care_beds_available, rev=true)
# 8 sites selected include top 3 for paediatric ICU beds and 5 of the 7 that have paediatric ICU beds
println(sampling_sites_london_8sites)
# Run line below once created 'sampling_sites_london_19sites'
println(sort(sampling_sites_london_19sites, :Paediatric_intensive_care_beds_available, rev=true)[:,[1,2,4]])

# 10 sites
sampling_sites_london_10sites = ari_cc_bed_london_sorted[1:9,:]
append!( sampling_sites_london_10sites, filter( row -> row.NHS_Trust_code == "RP4" , ari_cc_bed_london_sorted) )
sampling_sites_london_10sites = sampling_sites_london_10sites[:, [:NHS_Trust_code,:NHS_Trust_name
                                 ,:Adult_critical_care_beds_available,:Paediatric_intensive_care_beds_available] ]
# Check how many of the largest paediatric ICUs are including when rank by adult ICU size
# Top 4 paediatric ICUs also included in top 9 adult ICUs, so including RP4 GOSH, the 10 sites includes top 9 adult and top 5 paediatric
println(sampling_sites_london_10sites)
# Run line below once created 'sampling_sites_london_19sites'
println(sort(sampling_sites_london_19sites, :Paediatric_intensive_care_beds_available, rev=true)[:,[1,2,4]])

# All sites 19
sampling_sites_london_19sites = ari_cc_bed_london_sorted[:, [:NHS_Trust_code,:NHS_Trust_name
                                                            ,:Adult_critical_care_beds_available,:Paediatric_intensive_care_beds_available] ]

# Adjust so matches the format required for icu_td function
# Target column names = NHS_Trust_code, ICU_site_name, ICU_full_name, ICU_beds_adult, ICU_beds_paediatric, sample_target_per_year,sample_target_per_week

# Define datasets with their labels
london_datasets = [
    ("8sites",  sampling_sites_london_8sites),
    ("10sites", sampling_sites_london_10sites),
    ("19sites", sampling_sites_london_19sites),
]

for (label, df) in london_datasets

    # London only: $(label)
    rename!(df, :NHS_Trust_name                           => :ICU_site_name
               ,:Adult_critical_care_beds_available       => :ICU_beds_adult
               ,:Paediatric_intensive_care_beds_available => :ICU_beds_paediatric)

    insertcols!(df, 3, :ICU_full_name => Vector{Union{Missing, String}}(missing, nrow(df)))

    insertcols!(df, 6, :sample_target_per_year => Vector{Union{Missing, Int}}(missing, nrow(df))
                      ,:sample_target_per_week => Vector{Union{Missing, String}}(missing, nrow(df)))

    # Save file
    CSV.write( joinpath( geo_optim_folder, "selected_sites/sampling_sites_london_$(label).csv"), df)

end


#############################################################################################################
# (4) Sampling infections to simulate TDs using different sampling strategies:
#     (4.2) London only
#         (4.2.2) Simulate sampling and TD
#############################################################################################################

### Compute ICU times to detection (TD and TD3) using the London only NHS Trusts
# Load parameters used in simulation
initialize_parameters("config/outbreak_params_covid19_like_20261008.yaml")
# Double check loaded parameters required for ICU sampling
NBPMscape.P.pathogen_type # "virus"
NBPMscape.P.icu_sample_type #"regional" 
NBPMscape.P.sample_icu_cases_version #"number"
NBPMscape.P.icu_only_sample_before_death #true            

# Define variable parameters (18 sets)
sampling_site_files=["sampling_sites_london_8sites.csv","sampling_sites_london_10sites.csv","sampling_sites_london_19sites.csv"]
samples_per_week = [300,500,1000]
seasonal_icu_ari_admissions = [1440,793]

# Reload simulation data file if necessary
sims = load( normpath( joinpath( root_dir, "NBPMscape/paper/analysis/1_outbreak_simulations/covid_like/output/covidlike-1.6.1-sims_G_filtered_nrep1000_4300448.jld2") )
            , "sims")

# Loop through variables and run ICU sampling to compute times to detection
# ~1 to 15 mins per loop on laptop (10 sites are quicker than 30 sites)
counter = 0
for i in 1:length(sampling_site_files) #i=1
    for j in 1:length(seasonal_icu_ari_admissions) #j=1
        for k in 1:length(samples_per_week) # k=1
            
            # Read in NHS Trust sampling sites
            nhs_trust_sampling_sites = CSV.read( joinpath( geo_optim_folder, "selected_sites", sampling_site_files[i] )
                                                , DataFrame  )

            # Compute times to detection by sampling ICU cases
            temp_icu_td = icu_td( p = NBPMscape.P, sims = sims
                                , n_icu_samples_per_week = samples_per_week[k]
                                , icu_ari_admissions = seasonal_icu_ari_admissions[j]
                                , nhs_trust_sampling_sites = nhs_trust_sampling_sites
                                )

            # Save output to file
            output_file = "icu_tds_seasonal_ari_$(seasonal_icu_ari_admissions[j])_$(samples_per_week[k])_$(sampling_site_files[i])"
            CSV.write( joinpath( geo_optim_folder, "icu_td/london_only", output_file ), temp_icu_td)
            
            # Also save a smaller file as easier to manipulate
            # (remove column "ICU_simid")
            temp_icu_td_small = temp_icu_td[:,["sim_n", "ICU_TD", "ICU_3TD", "n_ICU_cases", "n_ICU_cases_sampled", "n_ICU_cases_sampled_positive", "tinf_relating_to_ICU_TD", "last_tinf_relating_to_ICU_3TD"]]
            output_file_small = "small_icu_tds_seasonal_ari_$(seasonal_icu_ari_admissions[j])_$(samples_per_week[k])_$(sampling_site_files[i])"
            CSV.write( joinpath( geo_optim_folder, "icu_td/london_only_small", output_file_small ), temp_icu_td_small)
            
            # Print progress update
            counter = counter + 1
            println("Successfully completed $(counter) of $(length(sampling_site_files)*length(seasonal_icu_ari_admissions)*length(samples_per_week)) runs")
        end
    end
end

# Reset parameters
initialize_parameters("config/outbreak_params_covid19_like_20261008.yaml")

### Compute median TDs for different sampling parameter sets
# Read in saved data
path_name = joinpath( geo_optim_folder, "icu_td/london_only_small" )

# Loop through scenarios, read in files and add TD values to dict
for i in 1:length(sampling_site_files) #i=1
    for j in 1:length(seasonal_icu_ari_admissions) #j=1
        for k in 1:length(samples_per_week) # k=1

            # Define parameter name to be used in dictionary. 
            # This is also the same as the file prefix as used above
            sampling_sites_prefix = sampling_site_files[i][1:findfirst(==('.'), sampling_site_files[i])-1]
            params_name = "icu_tds_seasonal_ari_$(seasonal_icu_ari_admissions[j])_$(samples_per_week[k])_$(sampling_sites_prefix)"
            # Use smaller size file for speed (without simid column)
            file_icu = joinpath(path_name, "small_$(params_name).csv")
    
            # Read CSV into a DataFrame and store in the dictionary with the scenario number (numeric suffix) as key
            icu_tds_df_dict[params_name] = CSV.read(file_icu, DataFrame)[:,[:ICU_TD,:ICU_3TD]]
        end
    end
end


# Loop through dictionary and add columns to the new dfs in dictionary (London only sites) for
# the cumulative outbreak count metrics

london_keys = [k for k in keys(icu_tds_df_dict) if occursin("london", k)]

for key in london_keys # keys(icu_tds_df_dict)
    df = icu_tds_df_dict[key] # df = icu_tds_df_dict["icu_tds_seasonal_ari_1440_300_sampling_sites_london_8sites"]
    
    # Add new columns
    empty_col_vec = Vector{Union{Missing, Int}}(missing, nrow(df))
    # TD
    df.cum_tinf_td       = copy(empty_col_vec)
    df.cum_tgp_td        = copy(empty_col_vec)
    df.cum_ted_td        = copy(empty_col_vec)
    df.cum_thosp_td      = copy(empty_col_vec)
    df.cum_ticu_td       = copy(empty_col_vec)
    df.cum_tdeceased_td  = copy(empty_col_vec)
    df.cum_tinf_fatal_td = copy(empty_col_vec)
    df.cum_timport_td    = copy(empty_col_vec)
    # TD3
    df.cum_tinf_td3       = copy(empty_col_vec)
    df.cum_tgp_td3        = copy(empty_col_vec)
    df.cum_ted_td3        = copy(empty_col_vec)
    df.cum_thosp_td3      = copy(empty_col_vec)
    df.cum_ticu_td3       = copy(empty_col_vec)
    df.cum_tdeceased_td3  = copy(empty_col_vec)
    df.cum_tinf_fatal_td3 = copy(empty_col_vec)
    df.cum_timport_td3    = copy(empty_col_vec)

    for i in 1:1000 #i=1
        td = df[i,:ICU_TD]
        td3 = df[i,:ICU_3TD]
        # Counts at TD
        df.cum_tinf_td[i]       = sum(outbreak_metrics_time_series.tinf_vals[i]       .<= td )
        df.cum_tgp_td[i]        = sum(outbreak_metrics_time_series.tgp_vals[i]        .<= td )
        df.cum_ted_td[i]        = sum(outbreak_metrics_time_series.ted_vals[i]        .<= td )
        df.cum_thosp_td[i]      = sum(outbreak_metrics_time_series.thosp_vals[i]      .<= td )
        df.cum_ticu_td[i]       = sum(outbreak_metrics_time_series.ticu_vals[i]       .<= td )
        df.cum_tdeceased_td[i]  = sum(outbreak_metrics_time_series.tdeceased_vals[i]  .<= td )
        df.cum_tinf_fatal_td[i] = sum(outbreak_metrics_time_series.tinf_fatal_vals[i] .<= td )
        df.cum_timport_td[i]    = sum(outbreak_metrics_time_series.timport_vals[i]    .<= td )
        # Counts at TD3
        df.cum_tinf_td3[i]       = sum(outbreak_metrics_time_series.tinf_vals[i]       .<= td3 )
        df.cum_tgp_td3[i]        = sum(outbreak_metrics_time_series.tgp_vals[i]        .<= td3 )
        df.cum_ted_td3[i]        = sum(outbreak_metrics_time_series.ted_vals[i]        .<= td3 )
        df.cum_thosp_td3[i]      = sum(outbreak_metrics_time_series.thosp_vals[i]      .<= td3 )
        df.cum_ticu_td3[i]       = sum(outbreak_metrics_time_series.ticu_vals[i]       .<= td3 )
        df.cum_tdeceased_td3[i]  = sum(outbreak_metrics_time_series.tdeceased_vals[i]  .<= td3 )
        df.cum_tinf_fatal_td3[i] = sum(outbreak_metrics_time_series.tinf_fatal_vals[i] .<= td3 )
        df.cum_timport_td3[i]    = sum(outbreak_metrics_time_series.timport_vals[i]    .<= td3 )
    end
    icu_tds_df_dict[key] = df
end

# Check
names(icu_tds_df_dict["icu_tds_seasonal_ari_1440_300_sampling_sites_london_10sites"])
df = icu_tds_df_dict["icu_tds_seasonal_ari_1440_300_sampling_sites_london_10sites"]
td = df[i,:ICU_TD]
sum(outbreak_metrics_time_series.tinf_vals[i] .<= td ) == df.cum_tinf_td[i]

icu_tds_df_dict["icu_tds_seasonal_ari_1440_300_sampling_sites_london_8sites"]
df = icu_tds_df_dict["icu_tds_seasonal_ari_1440_300_sampling_sites_london_8sites"]
td = df[i,:ICU_TD]
sum(outbreak_metrics_time_series.tinf_vals[i] .<= td ) == df.cum_tinf_td[i]

# Save as file
@save joinpath( geo_optim_folder, "covidlike-1.6.1_nrep1000_4300448_icu_tds_df_dict.jld2") icu_tds_df_dict


#########################################################################################
# (4) Sampling infections to simulate TDs using different sampling strategies:
#     (4.2) London only
#         (4.2.4) Plot results
#########################################################################################

# See 'NBPMscape/scripts/paper/supp_mat_icu_plots.jl' for plotting


#########################################################################################
# (4) Sampling infections to simulate TDs using different sampling strategies:
#     (4.3) Population size
#        (4.3.1) Select sampling sites
#########################################################################################

#### Method of selecting sites based on population size
println(ITL2SIZE)
println(sort(ITL2SIZE,:total_population_2022, rev = true))
# Sample (with replacement) from ITL2 regions with sampling weighted by region population size.
# Once have list of regions, select NHS Trust most likely to visit based on catchment area. 
# If ITL2 region is in the list twice then select different NHS Trusts

# Remove ITL2 regions in Wales because metagenomic sampling cannot currently be done here in the model
itl2_in_wales = ["TLL3","TLL4","TLL5"]
itl2size_ex_wales = filter( row -> !(row.ITL225CD in itl2_in_wales), ITL2SIZE)
itl2size_ex_wales.total_population_2022_prop = itl2size_ex_wales.total_population_2022 ./ sum(itl2size_ex_wales.total_population_2022)

# Select regions where sampling sites will be located, based on population size
region_list_8sites = wsample(itl2size_ex_wales.ITL225CD, itl2size_ex_wales.total_population_2022_prop, 8, replace=true) # TEST wsample([1,10],[0.1,0.9],10) 
region_list_26sites = wsample(itl2size_ex_wales.ITL225CD, itl2size_ex_wales.total_population_2022, 26, replace=true)

# Check on 8 sites selected
region_list_8sites_df = filter( row -> row.ITL225CD in region_list_8sites, itl2size_ex_wales)
region_list_8sites_df = sort(region_list_8sites_df,:total_population_2022,rev=true)
# Add number of sites selected in each region
n_sites = []
itl2_8sites_cm = countmap(region_list_8sites)
for i in 1:size(region_list_8sites_df,1) #i=1
    push!(n_sites,itl2_8sites_cm[region_list_8sites_df.ITL225CD[i]])
end
region_list_8sites_df.n_sites = n_sites
# Add population rank
itl2size_ex_wales_sorted_by_pop = sort(itl2size_ex_wales,:total_population_2022,rev=true)
pop_rank=[]
for i in 1:size(region_list_8sites_df,1) #i=1
    push!(pop_rank,only(findall(x -> x == region_list_8sites_df[i,:ITL225CD], itl2size_ex_wales_sorted_by_pop[!, :ITL225CD])))
end
region_list_8sites_df.pop_rank = pop_rank
println(region_list_8sites_df)

#6×6 DataFrame
# Row │ ITL225CD  ITL225NM                           total_population_2022  total_population_2022_prop  n_sites  pop_rank 
#     │ String    String                             Float64                Float64                     Any      Any
#─────┼───────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#   1 │ TLJ2      Surrey, East and West Sussex                   2.93556e6                   0.0514051  2        2
#   2 │ TLJ1      Berkshire, Buckinghamshire and O…              2.54967e6                   0.0446477  1        4
#   3 │ TLE4      West Yorkshire                                 2.37815e6                   0.0416442  1        5
#   4 │ TLI4      Inner London - East                            2.35549e6                   0.0412473  2        6
#   5 │ TLK6      North Somerset, Somerset and Dor…              1.58117e6                   0.0276881  1        17
#   6 │ TLG1      Herefordshire, Worcestershire an…              1.40554e6                   0.0246126  1        20

# Select regions where sampling sites will be located, based on population size
region_list_26sites = wsample(itl2size_ex_wales.ITL225CD, itl2size_ex_wales.total_population_2022, 26, replace=true)
# Check on 26 sites selected
region_list_26sites_df = filter( row -> row.ITL225CD in region_list_26sites, itl2size_ex_wales)
region_list_26sites_df = sort(region_list_26sites_df,:total_population_2022,rev=true)
# Add number of sites selected in each region
n_sites = []
itl2_26sites_cm = countmap(region_list_26sites)
for i in 1:size(region_list_26sites_df,1) #i=1
    push!(n_sites,itl2_26sites_cm[region_list_26sites_df.ITL225CD[i]])
end
region_list_26sites_df.n_sites = n_sites
# Add population rank
itl2size_ex_wales_sorted_by_pop = sort(itl2size_ex_wales,:total_population_2022,rev=true)
pop_rank=[]
for i in 1:size(region_list_26sites_df,1) #i=1
    push!(pop_rank,only(findall(x -> x == region_list_26sites_df[i,:ITL225CD], itl2size_ex_wales_sorted_by_pop[!, :ITL225CD])))
end
region_list_26sites_df.pop_rank = pop_rank
println(region_list_26sites_df)

# Reformat df so can use generate_nhs_trust_sampling_sites function to select NHS Trusts and save file
rename!(region_list_8sites_df, :ITL225CD => :ITL2_code
                             , :ITL225NM => :ITL2_name
                             , :n_sites => :sites_allocated
    )
rename!(region_list_26sites_df, :ITL225CD => :ITL2_code
                              , :ITL225NM => :ITL2_name
                              , :n_sites => :sites_allocated
    )
# Save to file
CSV.write( joinpath( root_dir, "NBPMscape/paper/analysis/2_sampling_analysis/covid_like/geo_optimisation/pop_size_site_numbers"
                    , "site_select_popsize_8sites.csv")
                    , region_list_8sites_df)
CSV.write( joinpath( root_dir, "NBPMscape/paper/analysis/2_sampling_analysis/covid_like/geo_optimisation/pop_size_site_numbers"
                    ,"site_select_popsize_26sites.csv")
                    , region_list_26sites_df)


# Select NHS Trusts for sampling site location
site_selection_folder = joinpath( root_dir, "NBPMscape/paper/analysis/2_sampling_analysis/covid_like/geo_optimisation/pop_size_site_numbers" )
sampling_site_file_output_folder = joinpath( root_dir, "NBPMscape/paper/analysis/2_sampling_analysis/covid_like/geo_optimisation/selected_sites" )

# Note that 'prioritise_adult_sampling = true' means that probabilities relating to adults for ITL2 -> NHS Trust will be used to select the first
# NHS Trust and children for the 2nd trust, and repeats in that order until all NHS Trusts have been selected.
select_nhs_trust_sampling_sites(; prioritise_adult_sampling = true,n_sites_max = 127
                                    ,site_selection_folder = site_selection_folder, site_selection_file = "site_select_popsize_8sites.csv"
                                    ,output_folder = sampling_site_file_output_folder, output_file = "sampling_sites_popsize_8sites.csv")

select_nhs_trust_sampling_sites(; prioritise_adult_sampling = true,n_sites_max = 127
                                    ,site_selection_folder = site_selection_folder, site_selection_file = "site_select_popsize_26sites.csv"
                                    ,output_folder = sampling_site_file_output_folder, output_file = "sampling_sites_popsize_26sites.csv")

# CHECK NHS Trusts selected and ITL2 regions with site allocations
# Some slightly odd selections using this logic:
# - In 8 site selection TLJ2 (Surrey, East and West Sussex) has 2 sites allocated. Adult prioritorisation selects RYR (University Hospitals Sussex NHS Foundation Trust)
#   with 35% probability of adult emergency attendance. The second site is selected based on child probability of attendance, however, the top probability with 93% is RYR 
#   which has already been selected and so the second highest probability is selected (RJ7 St George's University Hospitals NHS Foundation Trust) but this is only 5% for
#   children in TLJ2 and only 2% for adults. If only probabilities for adults were used then the 2nd choice site would be RTP with a 13% probability.


#########################################################################################
# (4) Sampling infections to simulate TDs using different sampling strategies:
#     (4.3) Population size
#        (4.3.2) Simulate sampling and TD
#########################################################################################

# Have 8 different sets of sampling sites. Run ICU sampling with 300, 500 and 1000 samples per week in winter and summer. 

# Load parameters used in simulation
initialize_parameters("config/outbreak_params_covid19_like_20261008.yaml")
# Double check loaded parameters required for ICU sampling
NBPMscape.P.pathogen_type # "virus"
NBPMscape.P.icu_sample_type #"regional" 
NBPMscape.P.sample_icu_cases_version #"number"
NBPMscape.P.icu_only_sample_before_death #true            

# Define variable parameters (12 sets)
sampling_site_files=["sampling_sites_popsize_8sites.csv","sampling_sites_popsize_26sites.csv"]
samples_per_week = [300,500,1000]
seasonal_icu_ari_admissions = [1440,793]

# Loop through variables and run ICU sampling to compute times to detection
# ~1 to 15 mins per loop on laptop (10 sites are quicker than 30 sites)
counter = 0
for i in 1:length(sampling_site_files) #i=1
    for j in 1:length(seasonal_icu_ari_admissions) #j=1
        for k in 1:length(samples_per_week) # k=1
            
            # Read in NHS Trust sampling sites
            nhs_trust_sampling_sites = CSV.read( joinpath( sampling_site_file_output_folder
                                                        ,sampling_site_files[i] ), DataFrame  )

            # Compute times to detection by sampling ICU cases
            temp_icu_td = icu_td( p = NBPMscape.P, sims = sims
                                , n_icu_samples_per_week = samples_per_week[k]
                                , icu_ari_admissions = seasonal_icu_ari_admissions[j]
                                , nhs_trust_sampling_sites = nhs_trust_sampling_sites
                                )

            # Save output to file
            output_folder = joinpath( geo_optim_folder, "icu_td/popsize" )
            output_file = "icu_tds_seasonal_ari_$(seasonal_icu_ari_admissions[j])_$(samples_per_week[k])_$(sampling_site_files[i])"
            CSV.write( joinpath( output_folder, output_file ), temp_icu_td)

            # Also save a smaller file as easier to manipulate
            # (remove column "ICU_simid")
            temp_icu_td_small = temp_icu_td[:,["sim_n", "ICU_TD", "ICU_3TD", "n_ICU_cases", "n_ICU_cases_sampled", "n_ICU_cases_sampled_positive", "tinf_relating_to_ICU_TD", "last_tinf_relating_to_ICU_3TD"]]
            output_file_small = "small_icu_tds_seasonal_ari_$(seasonal_icu_ari_admissions[j])_$(samples_per_week[k])_$(sampling_site_files[i])"
            CSV.write( joinpath( geo_optim_folder, "icu_td/popsize_small", output_file_small ), temp_icu_td_small)
            
            # Print progress update
            counter = counter + 1
            println("Successfully completed $(counter) of $(length(sampling_site_files)*length(seasonal_icu_ari_admissions)*length(samples_per_week)) runs")
        end
    end
end

# Reset parameters
initialize_parameters("config/outbreak_params_covid19_like_20261008.yaml")

### Compute median TDs for different sampling parameter sets
# Read in saved data
path_name = joinpath( root_dir, "NBPMscape/paper/analysis/2_sampling_analysis/covid_like/geo_optimisation/icu_td" )

# Loop through scenarios, read in files and compute median TD values
for i in 1:length(sampling_site_files) #i=1
    for j in 1:length(seasonal_icu_ari_admissions) #j=1
        for k in 1:length(samples_per_week) # k=1

            # Define parameter name to be used in dictionary. 
            # This is also the same as the file prefix as used above
            sampling_sites_prefix = sampling_site_files[i][1:findfirst(==('.'), sampling_site_files[i])-1]
            params_name = "icu_tds_seasonal_ari_$(seasonal_icu_ari_admissions[j])_$(samples_per_week[k])_$(sampling_sites_prefix)"
            file_icu = joinpath(path_name, "popsize_small/small_$(params_name).csv")
    
            # Read CSV into a DataFrame and store in the dictionary with the scenario number (numeric suffix) as key
            icu_tds_df_dict[params_name] = CSV.read(file_icu, DataFrame)[:,[:ICU_TD,:ICU_3TD]]
        end
    end
end

#########################################################################################
# (4) Sampling infections to simulate TDs using different sampling strategies:
#     (4.3) Population size
#        (4.3.3) Add outbreak metrics relating to TDs
#########################################################################################

### Compute outbreak count metrics at TD for different sampling parameter sets

# Loop through dictionary and add columns to the new dfs in dictionary 
# for the cumulative outbreak count metrics
icu_tds_df_dict["icu_tds_seasonal_ari_1440_300_sampling_sites_popsize_8sites"]
popsize_keys=[
             # 8 sites 
              "icu_tds_seasonal_ari_793_300_sampling_sites_popsize_8sites"
             ,"icu_tds_seasonal_ari_793_500_sampling_sites_popsize_8sites"
             ,"icu_tds_seasonal_ari_793_1000_sampling_sites_popsize_8sites"
             ,"icu_tds_seasonal_ari_1440_300_sampling_sites_popsize_8sites"
             ,"icu_tds_seasonal_ari_1440_500_sampling_sites_popsize_8sites"
             ,"icu_tds_seasonal_ari_1440_1000_sampling_sites_popsize_8sites"
             # 26 sites 
             ,"icu_tds_seasonal_ari_793_300_sampling_sites_popsize_26sites"
             ,"icu_tds_seasonal_ari_793_500_sampling_sites_popsize_26sites"
             ,"icu_tds_seasonal_ari_793_1000_sampling_sites_popsize_26sites"
             ,"icu_tds_seasonal_ari_1440_300_sampling_sites_popsize_26sites"
             ,"icu_tds_seasonal_ari_1440_500_sampling_sites_popsize_26sites"
             ,"icu_tds_seasonal_ari_1440_1000_sampling_sites_popsize_26sites"
             ]

# Reload if necessary (produced in (3) above)
# Reload outbreak metrics dictionary if necessary
#outbreak_metrics_time_series = load( joinpath( root_dir, "NBPMscape/paper/analysis/1_outbreak_simulations/covid_like/output"
#                                              , "covidlike-1.6.1_nrep1000_4300448_outbreak_metrics_time_series.jld2") 
#                                    , "outbreak_metrics_time_series" )


for key in popsize_keys # keys(icu_tds_df_dict)
    df = icu_tds_df_dict[key] # df = icu_tds_df_dict["icu_tds_seasonal_ari_1440_300_sampling_sites_5blocks_8sites"]
    
    # Add new columns
    empty_col_vec = Vector{Union{Missing, Int}}(missing, nrow(df))
    # TD
    df.cum_tinf_td       = copy(empty_col_vec)
    df.cum_tgp_td        = copy(empty_col_vec)
    df.cum_ted_td        = copy(empty_col_vec)
    df.cum_thosp_td      = copy(empty_col_vec)
    df.cum_ticu_td       = copy(empty_col_vec)
    df.cum_tdeceased_td  = copy(empty_col_vec)
    df.cum_tinf_fatal_td = copy(empty_col_vec)
    df.cum_timport_td    = copy(empty_col_vec)
    # TD3
    df.cum_tinf_td3       = copy(empty_col_vec)
    df.cum_tgp_td3        = copy(empty_col_vec)
    df.cum_ted_td3        = copy(empty_col_vec)
    df.cum_thosp_td3      = copy(empty_col_vec)
    df.cum_ticu_td3       = copy(empty_col_vec)
    df.cum_tdeceased_td3  = copy(empty_col_vec)
    df.cum_tinf_fatal_td3 = copy(empty_col_vec)
    df.cum_timport_td3    = copy(empty_col_vec)

    for i in 1:1000 #i=1
        td = df[i,:ICU_TD]
        td3 = df[i,:ICU_3TD]
        # Counts at TD
        df.cum_tinf_td[i]       = sum(outbreak_metrics_time_series.tinf_vals[i]       .<= td )
        df.cum_tgp_td[i]        = sum(outbreak_metrics_time_series.tgp_vals[i]        .<= td )
        df.cum_ted_td[i]        = sum(outbreak_metrics_time_series.ted_vals[i]        .<= td )
        df.cum_thosp_td[i]      = sum(outbreak_metrics_time_series.thosp_vals[i]      .<= td )
        df.cum_ticu_td[i]       = sum(outbreak_metrics_time_series.ticu_vals[i]       .<= td )
        df.cum_tdeceased_td[i]  = sum(outbreak_metrics_time_series.tdeceased_vals[i]  .<= td )
        df.cum_tinf_fatal_td[i] = sum(outbreak_metrics_time_series.tinf_fatal_vals[i] .<= td )
        df.cum_timport_td[i]    = sum(outbreak_metrics_time_series.timport_vals[i]    .<= td )
        # Counts at TD3
        df.cum_tinf_td3[i]       = sum(outbreak_metrics_time_series.tinf_vals[i]       .<= td3 )
        df.cum_tgp_td3[i]        = sum(outbreak_metrics_time_series.tgp_vals[i]        .<= td3 )
        df.cum_ted_td3[i]        = sum(outbreak_metrics_time_series.ted_vals[i]        .<= td3 )
        df.cum_thosp_td3[i]      = sum(outbreak_metrics_time_series.thosp_vals[i]      .<= td3 )
        df.cum_ticu_td3[i]       = sum(outbreak_metrics_time_series.ticu_vals[i]       .<= td3 )
        df.cum_tdeceased_td3[i]  = sum(outbreak_metrics_time_series.tdeceased_vals[i]  .<= td3 )
        df.cum_tinf_fatal_td3[i] = sum(outbreak_metrics_time_series.tinf_fatal_vals[i] .<= td3 )
        df.cum_timport_td3[i]    = sum(outbreak_metrics_time_series.timport_vals[i]    .<= td3 )
    end
    icu_tds_df_dict[key] = df
end

# Save as file
@save joinpath( root_dir, "NBPMscape/paper/analysis/2_sampling_analysis/covid_like/geo_optimisation"
            ,"covidlike-1.6.1_nrep1000_4300448_icu_tds_df_dict.jld2") icu_tds_df_dict
# Reload if necessary
#icu_tds_df_dict = load(joinpath( root_dir, "scripts/paper/2_sampling_analysis/covid_like/1717144_1719024_1719029_analysis/geo_optimisation"
#                                ,"covidlike-1.5.0_nrep1000_1717144_1719024_1719029_icu_tds_df_dict.jld2")
#                       , "icu_tds_df_dict")

#########################################################################################
# (4) Sampling infections to simulate TDs using different sampling strategies:
#     (4.3) Population size
#        (4.3.4) Plot results
#########################################################################################

# See 'NBPMscape/scripts/paper/supp_mat_icu_plots.jl' for plotting



#########################################################################################
#(4) Sampling infections to simulate TDs using different sampling strategies:
#    (4.4) mSCAPE pilot sites
#        (4.4.1) Select sampling sites
#        (4.4.2) Simulate sampling and TD
#########################################################################################

# Two sets of mSCAPE sites used in sampling investigation:
# - 8 sites from pilot study
# - 26 sites on longlist for potential expansion of programme (all 26 sites in England and excluding longlist sites in Scotland, Wales and Northern Ireland)

# Load parameters used in simulation
initialize_parameters("config/outbreak_params_covid19_like_20261008.yaml")
# Double check loaded parameters required for ICU sampling
NBPMscape.P.pathogen_type # "virus"
NBPMscape.P.icu_sample_type #"regional" 
NBPMscape.P.sample_icu_cases_version #"number"
NBPMscape.P.icu_only_sample_before_death #true            

# Define variable parameters (12 sets)
                            # label for output file             # File name to be read and used for sampling
sampling_site_files = [  (input_name = "nhs_trust_site_sample_targets_24pilot.csv", output_name = "sampling_sites_mscape_8sites.csv")  
                        ,(input_name = "nhs_trust_site_sample_targets_24pilot_25cohort_longlist.csv", output_name = "sampling_sites_mscape_26sites.csv")
                        ]
samples_per_week = [300,500,1000]
seasonal_icu_ari_admissions = [1440,793]

# Loop through variables and run ICU sampling to compute times to detection
# ~1 to 15 mins per loop on laptop (10 sites are quicker than 30 sites)
counter = 0
for i in 1:length(sampling_site_files) #i=2
    for j in 1:length(seasonal_icu_ari_admissions) #j=1
        for k in 1:length(samples_per_week) # k=1
            
            # Read in NHS Trust sampling sites
            nhs_trust_sampling_sites = CSV.read( joinpath( root_dir, "NBPMscape/data" #sampling_site_file_output_folder
                                                        ,sampling_site_files[i].input_name ), DataFrame  )

            # Compute times to detection by sampling ICU cases
            temp_icu_td = icu_td( p = NBPMscape.P
                                , sims = sims
                                , n_icu_samples_per_week = samples_per_week[k]
                                , icu_ari_admissions = seasonal_icu_ari_admissions[j]
                                , nhs_trust_sampling_sites = nhs_trust_sampling_sites
                                )

            # Save output to file
            output_folder = joinpath( geo_optim_folder, "icu_td/mscape" )
            output_file = "icu_tds_seasonal_ari_$(seasonal_icu_ari_admissions[j])_$(samples_per_week[k])_$(sampling_site_files[i].output_name)"
            CSV.write( joinpath( output_folder, output_file ), temp_icu_td)

            # Also save a smaller file as easier to manipulate
            # (remove column "ICU_simid")
            temp_icu_td_small = temp_icu_td[:,["sim_n", "ICU_TD", "ICU_3TD", "n_ICU_cases", "n_ICU_cases_sampled", "n_ICU_cases_sampled_positive", "tinf_relating_to_ICU_TD", "last_tinf_relating_to_ICU_3TD"]]
            output_file_small = "small_icu_tds_seasonal_ari_$(seasonal_icu_ari_admissions[j])_$(samples_per_week[k])_$(sampling_site_files[i].output_name)"
            CSV.write( joinpath( geo_optim_folder, "icu_td/mscape_small", output_file_small ), temp_icu_td_small)
            
            # Print progress update
            counter = counter + 1
            println("Successfully completed $(counter) of $(length(sampling_site_files)*length(seasonal_icu_ari_admissions)*length(samples_per_week)) runs")
        end
    end
end

# Reset parameters
initialize_parameters("config/outbreak_params_covid19_like_20261008.yaml")

# Read in saved data
path_name = joinpath( root_dir, "NBPMscape/paper/analysis/2_sampling_analysis/covid_like/geo_optimisation/icu_td/mscape_small" )


# Loop through scenarios, read in files and compute median TD values
td_results_file_params = []
for i in 1:length(sampling_site_files) #i=1
    for j in 1:length(seasonal_icu_ari_admissions) #j=1
        for k in 1:length(samples_per_week) # k=1

            # Define parameter name to be used in dictionary. 
            # This is also the same as the file prefix as used above
            sampling_sites_prefix = sampling_site_files[i].output_name[1:findfirst(==('.'), sampling_site_files[i].output_name)-1]
            params_name = "icu_tds_seasonal_ari_$(seasonal_icu_ari_admissions[j])_$(samples_per_week[k])_$(sampling_sites_prefix)"
            file_icu = joinpath(path_name, "small_$(params_name).csv")
    
            # Read CSV into a DataFrame and store in the dictionary with the scenario number (numeric suffix) as key
            icu_tds_df_dict[params_name] = CSV.read(file_icu, DataFrame)[:,[:ICU_TD,:ICU_3TD]]
            
            push!(td_results_file_params,params_name)
        end
    end
end

# Loop through scenarios, read in files and compute median TD values
for i in 1:length(td_results_file_params) #i=1
    params_name = td_results_file_params[i]
    file_icu = joinpath(path_name, "small_$(params_name).csv")
    # Read CSV into a DataFrame and store in the dictionary with the scenario number (numeric suffix) as key
    icu_tds_df_dict[params_name] = CSV.read(file_icu, DataFrame)[:,[:ICU_TD,:ICU_3TD]]
end

# CHECK icu_tds_df_dict["icu_tds_seasonal_ari_1440_300_sampling_sites_mscape_8sites"]
# CHECK icu_tds_df_dict["icu_tds_seasonal_ari_1440_300_sampling_sites_mscape_26sites"]

#########################################################################################
#(4) Sampling infections to simulate TDs using different sampling strategies:
#    (4.4) mSCAPE pilot sites
#        (4.4.3) Add outbreak metrics relating to TDs
#########################################################################################

### Compute outbreak count metrics at TD for different sampling parameter sets

# Loop through dictionary and add columns to the new dfs in dictionary for the 
# cumulative outbreak count metrics
icu_tds_df_dict["icu_tds_seasonal_ari_1440_300_sampling_sites_mscape_8sites"]

# Reload if necessary (produced in (3) above)
# Reload outbreak metrics dictionary if necessary
#outbreak_metrics_time_series = load( joinpath( root_dir, "NBPMscape/paper/analysis/1_outbreak_simulations/covid_like/output"
#                                              , "covidlike-1.6.1_nrep1000_4300448_outbreak_metrics_time_series.jld2") 
#                                    , "outbreak_metrics_time_series" )

# Note that this for loop is slightly different to those used earlier
for j in 1:length(td_results_file_params) # keys(icu_tds_df_dict) # j=1
    key = td_results_file_params[j] #.params_name
    df = icu_tds_df_dict[key] # df = icu_tds_df_dict["icu_tds_seasonal_ari_1440_300_sampling_sites_5blocks_8sites"]
    
    # Add new columns
    empty_col_vec = Vector{Union{Missing, Int}}(missing, nrow(df))
    # TD
    df.cum_tinf_td       = copy(empty_col_vec)
    df.cum_tgp_td        = copy(empty_col_vec)
    df.cum_ted_td        = copy(empty_col_vec)
    df.cum_thosp_td      = copy(empty_col_vec)
    df.cum_ticu_td       = copy(empty_col_vec)
    df.cum_tdeceased_td  = copy(empty_col_vec)
    df.cum_tinf_fatal_td = copy(empty_col_vec)
    df.cum_timport_td    = copy(empty_col_vec)
    # TD3
    df.cum_tinf_td3       = copy(empty_col_vec)
    df.cum_tgp_td3        = copy(empty_col_vec)
    df.cum_ted_td3        = copy(empty_col_vec)
    df.cum_thosp_td3      = copy(empty_col_vec)
    df.cum_ticu_td3       = copy(empty_col_vec)
    df.cum_tdeceased_td3  = copy(empty_col_vec)
    df.cum_tinf_fatal_td3 = copy(empty_col_vec)
    df.cum_timport_td3    = copy(empty_col_vec)

    for i in 1:1000 #i=1
        td = df[i,:ICU_TD]
        td3 = df[i,:ICU_3TD]
        # Counts at TD
        df.cum_tinf_td[i]       = sum(outbreak_metrics_time_series.tinf_vals[i]       .<= td )
        df.cum_tgp_td[i]        = sum(outbreak_metrics_time_series.tgp_vals[i]        .<= td )
        df.cum_ted_td[i]        = sum(outbreak_metrics_time_series.ted_vals[i]        .<= td )
        df.cum_thosp_td[i]      = sum(outbreak_metrics_time_series.thosp_vals[i]      .<= td )
        df.cum_ticu_td[i]       = sum(outbreak_metrics_time_series.ticu_vals[i]       .<= td )
        df.cum_tdeceased_td[i]  = sum(outbreak_metrics_time_series.tdeceased_vals[i]  .<= td )
        df.cum_tinf_fatal_td[i] = sum(outbreak_metrics_time_series.tinf_fatal_vals[i] .<= td )
        df.cum_timport_td[i]    = sum(outbreak_metrics_time_series.timport_vals[i]    .<= td )
        # Counts at TD3
        df.cum_tinf_td3[i]       = sum(outbreak_metrics_time_series.tinf_vals[i]       .<= td3 )
        df.cum_tgp_td3[i]        = sum(outbreak_metrics_time_series.tgp_vals[i]        .<= td3 )
        df.cum_ted_td3[i]        = sum(outbreak_metrics_time_series.ted_vals[i]        .<= td3 )
        df.cum_thosp_td3[i]      = sum(outbreak_metrics_time_series.thosp_vals[i]      .<= td3 )
        df.cum_ticu_td3[i]       = sum(outbreak_metrics_time_series.ticu_vals[i]       .<= td3 )
        df.cum_tdeceased_td3[i]  = sum(outbreak_metrics_time_series.tdeceased_vals[i]  .<= td3 )
        df.cum_tinf_fatal_td3[i] = sum(outbreak_metrics_time_series.tinf_fatal_vals[i] .<= td3 )
        df.cum_timport_td3[i]    = sum(outbreak_metrics_time_series.timport_vals[i]    .<= td3 )
    end
    icu_tds_df_dict[key] = df
end

# Save to file
@save joinpath( root_dir, "NBPMscape/paper/analysis/2_sampling_analysis/covid_like/geo_optimisation"
                , "covidlike-1.6.1_nrep1000_4300448_icu_tds_df_dict.jld2") icu_tds_df_dict
# Reload if necessary
#icu_tds_df_dict = load( joinpath( root_dir, "NBPMscape/paper/analysis/2_sampling_analysis/covid_like/geo_optimisation"
#                                    , "covidlike-1.6.1_nrep1000_4300448_icu_tds_df_dict.jld2")
#                        , "icu_tds_df_dict")

#########################################################################################
#(4) Sampling infections to simulate TDs using different sampling strategies:
#    (4.4) mSCAPE pilot sites
#        (4.4.4) Plot results
#########################################################################################

# See 'NBPMscape/scripts/paper/supp_mat_icu_plots.jl' for plotting

#########################################################################################
# (5) Compile TD statistics for all strategies and parameter combinations
#########################################################################################

# Reload dictionary containing data if necessary
#icu_tds_df_dict = load(joinpath( root_dir, "NBPMscape/paper/analysis/2_sampling_analysis/covid_like/geo_optimisation"
#                                ,"covidlike-1.6.1_nrep1000_4300448_icu_tds_df_dict.jld2")
#                       , "icu_tds_df_dict")


### Build df containing the keys (col1), TD, Var, SD, U_0.1, U_0.5, U_1, TD3, Var3, SD3, U3_0.1, U3_0.5, U3_1
# U_X represents different values of lambda from the proposed loss function for a given sample design θ:
# J(θ) = E[TD|θ] + λVar[TD|θ]
# 

all_keys = keys(icu_tds_df_dict)

results_df_col_format = [-1.0 for _ in 1:length(all_keys)]
results_df = DataFrame( parameter_set_name = ["" for _ in 1:length(all_keys)]
                        ,TD_median = copy(results_df_col_format)
                        ,iqr_lower = copy(results_df_col_format)
                        ,iqr_upper = copy(results_df_col_format)
                        ,TD_mean   = copy(results_df_col_format)
                        ,Var   = copy(results_df_col_format)
                        ,SD    = copy(results_df_col_format)
                        ,U_01  = copy(results_df_col_format)
                        ,U_05  = copy(results_df_col_format)
                        ,U_1   = copy(results_df_col_format)
                        ,TD3_median = copy(results_df_col_format)
                        ,iqr3_lower = copy(results_df_col_format)
                        ,iqr3_upper = copy(results_df_col_format)
                        ,TD3_mean = copy(results_df_col_format)
                        ,Var3  = copy(results_df_col_format)
                        ,SD3   = copy(results_df_col_format)
                        ,U3_01 = copy(results_df_col_format)
                        ,U3_05 = copy(results_df_col_format)
                        ,U3_1  = copy(results_df_col_format)
)

k = 0
#for k in 1:length(all_keys) # k=1; key = "icu_tds_seasonal_ari_1440_300_sampling_sites_10blocks_10sites"
for key in all_keys
    k = k +1
    println(key,k)
    results_df.parameter_set_name[k] = key

    results_df.TD_median[k] = median(icu_tds_df_dict[key][:,:ICU_TD])
    results_df.iqr_lower[k] = quantile(icu_tds_df_dict[key][:,:ICU_TD],0.25)
    results_df.iqr_upper[k] = quantile(icu_tds_df_dict[key][:,:ICU_TD],0.75)
    results_df.TD_mean[k] = mean(icu_tds_df_dict[key][:,:ICU_TD])
    results_df.Var[k] = var(icu_tds_df_dict[key][:,:ICU_TD])
    results_df.SD[k] = std(icu_tds_df_dict[key][:,:ICU_TD])
    results_df.U_01[k] = results_df.TD_median[k] + 0.1*results_df.Var[k] 
    results_df.U_05[k] = results_df.TD_median[k] + 0.5*results_df.Var[k] 
    results_df.U_1[k] = results_df.TD_median[k] + 1*results_df.Var[k] 
    
    results_df.TD3_median[k] = median(icu_tds_df_dict[key][:,:ICU_3TD])
    results_df.iqr3_lower[k] = quantile(icu_tds_df_dict[key][:,:ICU_3TD],0.25)
    results_df.iqr3_upper[k] = quantile(icu_tds_df_dict[key][:,:ICU_3TD],0.75)
    results_df.TD3_mean[k] = mean(icu_tds_df_dict[key][:,:ICU_3TD])
    results_df.Var3[k] = var(icu_tds_df_dict[key][:,:ICU_3TD])
    results_df.SD3[k] = std(icu_tds_df_dict[key][:,:ICU_3TD])
    results_df.U3_01[k] = results_df.TD3_median[k] .+ 0.1*results_df.Var3[k] 
    results_df.U3_05[k] = results_df.TD3_median[k] .+ 0.5*results_df.Var3[k] 
    results_df.U3_1[k] = results_df.TD3_median[k] .+ 1*results_df.Var3[k] 
    
end

CSV.write( joinpath( root_dir, "NBPMscape/paper/analysis/2_sampling_analysis/covid_like/geo_optimisation"
                    ,"icu_td_results_df.csv" )
                    , results_df )

######################################################################################################
# (6) Comparison of best ICU sampling performers
# (6.1) Plot best ICU sampling performers
######################################################################################################
### Plot best performers (from inspection of icu_td_results_df.csv)

# Helper function for filtering parameter names for multiple substrings
function contains_all(s::AbstractString, subs::Vector{String})
    all(sub -> occursin(sub, s), subs)
end

# Results for winter (1440 ARI ICU per week), with 8 sites sampling and 300 samples per week
results_df_winter_8sites_300samples = filter(row -> contains_all(row.parameter_set_name, ["1440","8sites","300"]), results_df)
results_df_winter_8sites_300samples = sort(results_df_winter_8sites_300samples, "TD_median")
results_df_winter_8sites_300samples[:,2:size(results_df_winter_8sites_300samples,2)] = round.( results_df_winter_8sites_300samples[:,2:size(results_df_winter_8sites_300samples,2)], digits=1)
println(results_df_winter_8sites_300samples)
println(results_df_winter_8sites_300samples.parameter_set_name)

#7×19 DataFrame
# Row │ parameter_set_name                                                  TD_median  iqr_lower  iqr_upper  TD_mean  Var      SD       U_01     U_05     U_1      TD3_median  iqr3_lower  iqr3_upper  TD3_mean  Var3     SD3      U3_01    U3_05    U3_1    
#     │ String                                                              Float64    Float64    Float64    Float64  Float64  Float64  Float64  Float64  Float64  Float64     Float64     Float64     Float64   Float64  Float64  Float64  Float64  Float64 
#─────┼─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#   1 │ icu_tds_seasonal_ari_1440_300_sampling_sites_20blocks_8sites        43.1       37.3       49.4       43.2     84.4      9.2     51.5     85.3    127.5        50.6        46.9        55.7      51.4     47.1      6.9     55.3     74.1     97.7
#   2 │ icu_tds_seasonal_ari_1440_300_sampling_sites_30blocks_8sites        43.4       38.1       48.9       43.2     78.7      8.9     51.2     82.7    122.0        50.7        46.1        55.5      51.0     49.3      7.0     55.6     75.3    100.0
#   3 │ icu_tds_seasonal_ari_1440_300_sampling_sites_10blocks_8sites        44.9       39.3       50.6       44.7     78.5      8.9     52.7     84.1    123.4        52.3        48.1        57.1      52.8     46.4      6.8     56.9     75.5     98.7
#   4 │ icu_tds_seasonal_ari_1440_300_sampling_sites_popsize_8sites         45.6       39.7       51.3       45.4     84.2      9.2     54.0     87.7    129.7        52.7        48.4        57.8      53.4     48.6      7.0     57.6     77.0    101.3
#   5 │ icu_tds_seasonal_ari_1440_300_sampling_sites_5blocks_8sites         45.7       40.8       52.0       46.1     75.8      8.7     53.3     83.6    121.5        53.6        49.1        58.5      54.0     48.4      7.0     58.4     77.8    101.9
#   6 │ icu_tds_seasonal_ari_1440_300_sampling_sites_london_8sites          45.9       40.2       51.8       45.9     84.2      9.2     54.4     88.1    130.2        53.0        48.7        58.2      53.7     52.3      7.2     58.2     79.1    105.3
#   7 │ icu_tds_seasonal_ari_1440_300_sampling_sites_mscape_8sites          48.6       42.4       54.4       48.4     78.2      8.8     56.4     87.7    126.7        55.7        51.3        61.1      56.2     48.4      7.0     60.5     79.9    104.1

# Conclusion: In this scenario (winter sampling ICU at 300 samples per week across 8 sites/NHS Trusts) the best performing strategy is HCGS.
#             While 20 blocks has a marginally earlier median time to detect 1 case (TD), the 30 block strategy has a lower variance. 
#             This results in a lower combined metric value across differing variance weights (10% (U01), 50% (U05) and 100% (U1)).
#             This is reversed for the time to detect 3 cases (TD3).
#             Difficult to select between the two, but chosen 'icu_tds_seasonal_ari_1440_300_sampling_sites_30blocks_8sites'

#https://juliagraphics.github.io/Colors.jl/stable/namedcolors/
hcgs_color = :mediumpurple3
popsize_color = :royalblue1
london_color = :gold       
mscape_color = :orangered

# Best performers with 8 sites, in winter, with 300 samples per week
# TD
p1 = kernel_box_jitter_plot(;x = icu_tds_df_dict["icu_tds_seasonal_ari_1440_300_sampling_sites_30blocks_8sites"][:,:ICU_TD], plot_color = hcgs_color,    samp_strategy="HCGS",            n_samples=300, n_sites=8, n_blocks=30, season="Winter", x_label="")
p2 = kernel_box_jitter_plot(;x = icu_tds_df_dict["icu_tds_seasonal_ari_1440_300_sampling_sites_popsize_8sites"][:, :ICU_TD], plot_color = popsize_color, samp_strategy="Population size", n_samples=300, n_sites=8,              season="Winter", x_label="")
p3 = kernel_box_jitter_plot(;x = icu_tds_df_dict["icu_tds_seasonal_ari_1440_300_sampling_sites_london_8sites"][:,  :ICU_TD], plot_color = london_color,  samp_strategy="London only",     n_samples=300, n_sites=8,              season="Winter", x_label="")
p4 = kernel_box_jitter_plot(;x = icu_tds_df_dict["icu_tds_seasonal_ari_1440_300_sampling_sites_mscape_8sites"][:,  :ICU_TD], plot_color = mscape_color,  samp_strategy="NHS-mSCAPE",      n_samples=300, n_sites=8,              season="Winter")

Plots.plot(p1,p2,p3,p4, layout = (4, 1)
        , size = (800,1200)) #default size is 600,400

Plots.savefig( joinpath( root_dir, "NBPMscape/paper/analysis/2_sampling_analysis/covid_like/geo_optimisation"
                        , "icu_tds_seasonal_ari_1440_300_8sites_best_performers.png" ) )

# Cumuatlive infections (same colours as TD)
p5 = kernel_box_jitter_plot(;x = collect(skipmissing(icu_tds_df_dict["icu_tds_seasonal_ari_1440_300_sampling_sites_30blocks_8sites"][:,:cum_tinf_td])), plot_color = hcgs_color,    samp_strategy="HCGS",            n_samples=300, n_sites=8, n_blocks=30, season="Winter", x_label="", x_min=0, x_max=150000, annot_x_pos="middle")
p6 = kernel_box_jitter_plot(;x = collect(skipmissing(icu_tds_df_dict["icu_tds_seasonal_ari_1440_300_sampling_sites_popsize_8sites"][:, :cum_tinf_td])), plot_color = popsize_color, samp_strategy="Population size", n_samples=300, n_sites=8,              season="Winter", x_label="", x_min=0, x_max=150000, annot_x_pos="middle")
p7 = kernel_box_jitter_plot(;x = collect(skipmissing(icu_tds_df_dict["icu_tds_seasonal_ari_1440_300_sampling_sites_london_8sites"][:,  :cum_tinf_td])), plot_color = london_color,  samp_strategy="London only",     n_samples=300, n_sites=8,              season="Winter", x_label="", x_min=0, x_max=150000, annot_x_pos="middle")
p8 = kernel_box_jitter_plot(;x = collect(skipmissing(icu_tds_df_dict["icu_tds_seasonal_ari_1440_300_sampling_sites_mscape_8sites"][:,  :cum_tinf_td])), plot_color = mscape_color,  samp_strategy="NHS-mSCAPE",      n_samples=300, n_sites=8,              season="Winter", x_label="Cumulative infections between first importation and TD", x_min=0, x_max=150000, annot_x_pos="middle")

Plots.plot(p5,p6,p7,p8, layout = (4, 1)
        , size = (800,1200)) #default size is 600,400

Plots.savefig( joinpath( root_dir, "NBPMscape/paper/analysis/2_sampling_analysis/covid_like/geo_optimisation"
                , "icu_tds_seasonal_ari_1440_300_8sites_best_performers_cum_inf_td.png" ) )

#println(icu_tds_df_dict["icu_tds_seasonal_ari_1440_300_sampling_sites_30blocks_8sites"][:,:cum_tinf_td])

# Plot combining TD and cumulative infections
Plots.plot(p1,p5,p2,p6,p3,p7,p4,p8, layout = (4, 2)
        , size = (1200,1200)) #default size is 600,400

Plots.savefig( joinpath( root_dir, "NBPMscape/paper/analysis/2_sampling_analysis/covid_like/geo_optimisation"
                , "icu_tds_seasonal_ari_1440_300_8sites_best_performers_td_and_cum_inf_td.png" ) )
Plots.savefig( joinpath( root_dir, "NBPMscape/paper/results/figures_main/raw"
                , "icu_tds_seasonal_ari_1440_300_8sites_best_performers_td_and_cum_inf_td.png" ) )


# Plots for time to detect 3 cases (3TD) are made in 
# 'NBPMscape/scripts/paper/supp_mat_icu_plots.jl'

######################################################################################################
# (6) Comparison of best ICU sampling performers
# (6.2) Compute statistics for relative performance
######################################################################################################

# Reload ICU TD data if necessary
icu_tds_df_dict = load( joinpath( root_dir, "NBPMscape/paper/analysis/2_sampling_analysis/covid_like/geo_optimisation"
                                    , "covidlike-1.6.1_nrep1000_4300448_icu_tds_df_dict.jld2")
                        , "icu_tds_df_dict")

# Compare best TD performer in winter, 8 sites, 300 samples per week (which was selected as HCGS 30 blocks)
# against mscape strategy performance

# Compare TD for each simulation replicate and take the mean and median
td_diff_hcgs_mscape = icu_tds_df_dict["icu_tds_seasonal_ari_1440_300_sampling_sites_30blocks_8sites"][:,:ICU_TD] .- icu_tds_df_dict["icu_tds_seasonal_ari_1440_300_sampling_sites_mscape_8sites"][:,:ICU_TD]

# Use bootstrapping to estimate uncertainty
td_diff_hcgs_mscape_median = median_ci_bootstrap(vec = td_diff_hcgs_mscape, n_boot=1000, alpha=0.05)
# (median = -4.6887511817578424, lower = -5.226565532137467, upper = -4.0455307197597525)
td_diff_hcgs_mscape_mean = mean_ci_bootstrap(vec = td_diff_hcgs_mscape, n_boot=1000, alpha=0.05)
# (mean = -5.143472814563969, lower = -5.701240611787326, upper = -4.574869363883158)

# Plot TD differences
Plots.histogram(td_diff_hcgs_mscape, label = "TD difference (1000 sim reps)", color=:red
                , title = "HCGS 30 blocks vs mSCAPE \n (8 sites, winter, 300 samples per week)", legend=:topleft
                , ylab = "Number of simulation replicates (out of 1000)")
# Add median info
vline!([td_diff_hcgs_mscape_median[1]], label="median = $(round(td_diff_hcgs_mscape_median[1]; digits=1))", color = :black) 
vline!([td_diff_hcgs_mscape_median[2]], label="median lower 95% bootstrap estimate = $(round(td_diff_hcgs_mscape_median[2]; digits=1))", color = :black, line=(:dash, 3)) 
vline!([td_diff_hcgs_mscape_median[3]], label="median upper 95% bootstrap estimate = $(round(td_diff_hcgs_mscape_median[3]; digits=1))", color = :black, line=(:dash, 3)) 
# Add mean info
vline!([td_diff_hcgs_mscape_mean[1]], label="mean = $(round(td_diff_hcgs_mscape_mean[1]; digits=1))", color = :blue) 
vline!([td_diff_hcgs_mscape_mean[2]], label="mean lower 95% bootstrap estimate = $(round(td_diff_hcgs_mscape_mean[2]; digits=1))", color = :blue, line=(:dash, 3)) 
vline!([td_diff_hcgs_mscape_mean[3]], label="mean upper 95% bootstrap estimate = $(round(td_diff_hcgs_mscape_mean[3]; digits=1))", color = :blue, line=(:dash, 3)) 

Plots.savefig( joinpath( root_dir, "NBPMscape/paper/analysis/2_sampling_analysis/covid_like/geo_optimisation"
                , "icu_tds_seasonal_ari_1440_300_8sites_HCGS30blocks_v_mscape_td.png" ) )

# Look at the IQR
quantile(td_diff_hcgs_mscape, 0.25) # -10.858797407037143
quantile(td_diff_hcgs_mscape, 0.75) # +0.3235341326700407

# Look at 95% quantile
quantile(td_diff_hcgs_mscape, 0.025) # -24.44850474967557
quantile(td_diff_hcgs_mscape, 0.975) # +13.237627336914393

# Also want to know the mean and upper and lower values for the percentage reduction in TD
# HCGS / mscape
td_perc_diff_hcgs_mscape = 100 .* ( icu_tds_df_dict["icu_tds_seasonal_ari_1440_300_sampling_sites_30blocks_8sites"][:,:ICU_TD] ./ icu_tds_df_dict["icu_tds_seasonal_ari_1440_300_sampling_sites_mscape_8sites"][:,:ICU_TD] )

Plots.plot( td_perc_diff_hcgs_mscape )

Plots.scatter(  icu_tds_df_dict["icu_tds_seasonal_ari_1440_300_sampling_sites_mscape_8sites"][:,:ICU_TD]
              , icu_tds_df_dict["icu_tds_seasonal_ari_1440_300_sampling_sites_30blocks_8sites"][:,:ICU_TD]
              , xlabel = "TD - mSCAPE"
             , ylabel = "TD - HCGS"
             , size = (800,400))

# Use bootstrapping to estimate uncertainty
td_perc_diff_hcgs_mscape_median = median_ci_bootstrap(vec = td_perc_diff_hcgs_mscape, n_boot=1000, alpha=0.05)
# (median = 90.28244536349186, lower = 89.19949099313176, upper = 92.0875056411349)
# adjust to make X% less than mscape...
td_perc_diff_hcgs_mscape_median_v2_median = 100 * (1 - (td_perc_diff_hcgs_mscape_median[1] / 100)) #9.71755463650814
td_perc_diff_hcgs_mscape_median_v2_lower  = 100 * (1 - (td_perc_diff_hcgs_mscape_median[2] / 100)) #10.800509006868243
td_perc_diff_hcgs_mscape_median_v2_upper  = 100 * (1 - (td_perc_diff_hcgs_mscape_median[3] / 100)) #7.912494358865107

td_perc_diff_hcgs_mscape_mean = mean_ci_bootstrap(vec = td_perc_diff_hcgs_mscape, n_boot=1000, alpha=0.05)
# (mean = 91.15728516918591, lower = 89.93241958203059, upper = 92.40677284854178)


# Compare cumulative infections at TD for each simulation replicate and take the mean and median
cum_inf_diff_hcgs_mscape = icu_tds_df_dict["icu_tds_seasonal_ari_1440_300_sampling_sites_30blocks_8sites"][:,:cum_tinf_td] .- icu_tds_df_dict["icu_tds_seasonal_ari_1440_300_sampling_sites_mscape_8sites"][:,:cum_tinf_td]

# Use bootstrapping to estimate uncertainty
cum_inf_diff_hcgs_mscape_median = median_ci_bootstrap(vec = cum_inf_diff_hcgs_mscape, n_boot=1000, alpha=0.05) 
# (median = -7901.0, lower = -9351.0, upper = -6482.962500000001)
cum_inf_diff_hcgs_mscape_mean = mean_ci_bootstrap(vec = cum_inf_diff_hcgs_mscape, n_boot=1000, alpha=0.05) 
# (mean = -12965.198, lower = -14248.73765, upper = -11603.2193)

# Plot cumulative infections at TD differences
Plots.histogram(cum_inf_diff_hcgs_mscape, label = "Difference in cumulative infections at TD (1000 sim reps)", color=:blue
                , title = "HCGS 30 blocks vs mSCAPE \n (8 sites, winter, 300 samples per week)", legend=:topleft
                , ylab = "Number of simulation replicates (out of 1000)")
# Add median info
vline!([cum_inf_diff_hcgs_mscape_median[1]], label="median = $(round(cum_inf_diff_hcgs_mscape_median[1]; digits=1))", color = :black) 
vline!([cum_inf_diff_hcgs_mscape_median[2]], label="median lower 95% bootstrap estimate = $(round(cum_inf_diff_hcgs_mscape_median[2]; digits=1))", color = :black, line=(:dash, 3)) 
vline!([cum_inf_diff_hcgs_mscape_median[3]], label="median upper 95% bootstrap estimate = $(round(cum_inf_diff_hcgs_mscape_median[3]; digits=1))", color = :black, line=(:dash, 3)) 
# Add mean info
vline!([cum_inf_diff_hcgs_mscape_mean[1]], label="mean = $(round(cum_inf_diff_hcgs_mscape_mean[1]; digits=1))", color = :lightblue) 
vline!([cum_inf_diff_hcgs_mscape_mean[2]], label="mean lower 95% bootstrap estimate = $(round(cum_inf_diff_hcgs_mscape_mean[2]; digits=1))", color = :lightblue, line=(:dash, 3)) 
vline!([cum_inf_diff_hcgs_mscape_mean[3]], label="mean upper 95% bootstrap estimate = $(round(cum_inf_diff_hcgs_mscape_mean[3]; digits=1))", color = :lightblue, line=(:dash, 3)) 

Plots.savefig( joinpath( root_dir, "NBPMscape/paper/analysis/2_sampling_analysis/covid_like/geo_optimisation"
                , "icu_tds_seasonal_ari_1440_300_8sites_HCGS30blocks_v_mscape_cum_inf_td.png" ) )

# Look at the IQR
quantile(cum_inf_diff_hcgs_mscape, 0.25) # -24609.0
quantile(cum_inf_diff_hcgs_mscape, 0.75) # +386.75

# Look at 95% quantile
quantile(cum_inf_diff_hcgs_mscape, 0.025) # -62189.8
quantile(cum_inf_diff_hcgs_mscape, 0.975) # +21415.55

# Also want to know the mean and upper and lower values for the percentage reduction in cumulative infections
# HCGS / mscape
cum_inf_perc_diff_hcgs_mscape = 100 .* ( icu_tds_df_dict["icu_tds_seasonal_ari_1440_300_sampling_sites_30blocks_8sites"][:,:cum_tinf_td] ./ icu_tds_df_dict["icu_tds_seasonal_ari_1440_300_sampling_sites_mscape_8sites"][:,:cum_tinf_td] )

cum_inf_at_td_hcgs_v_mscape = DataFrame( icu_tds_seasonal_ari_1440_300_sampling_sites_30blocks_8sites_cum_tinf_td = icu_tds_df_dict["icu_tds_seasonal_ari_1440_300_sampling_sites_30blocks_8sites"][:,:cum_tinf_td]
                                       ,icu_tds_seasonal_ari_1440_300_sampling_sites_mscape_8sites_cum_tinf_td = icu_tds_df_dict["icu_tds_seasonal_ari_1440_300_sampling_sites_mscape_8sites"][:,:cum_tinf_td]
                                       )

CSV.write( joinpath( root_dir
                    , "NBPMscape/paper/analysis/2_sampling_analysis/covid_like/geo_optimisation"
                    , "cum_inf_at_td_hcgs_v_mscape.csv")
         , cum_inf_at_td_hcgs_v_mscape )

Plots.plot( cum_inf_perc_diff_hcgs_mscape )

Plots.scatter(  icu_tds_df_dict["icu_tds_seasonal_ari_1440_300_sampling_sites_mscape_8sites"][:,:cum_tinf_td]
              , icu_tds_df_dict["icu_tds_seasonal_ari_1440_300_sampling_sites_30blocks_8sites"][:,:cum_tinf_td]
              , xlabel = "cumulative infections at TD - mSCAPE"
             , ylabel = "cumulative infections at TD - HCGS"
             , size = (675,600))

# Use bootstrapping to estimate uncertainty
cum_inf_perc_diff_hcgs_mscape_median = median_ci_bootstrap(vec = cum_inf_perc_diff_hcgs_mscape, n_boot=1000, alpha=0.05)
# (median = 42.937629305445725, lower = 39.011966364812416, upper = 49.575519197362304)
# adjust to make X% less than mscape...
cum_inf_perc_diff_hcgs_mscape_median_v2_median = 100 * (1 - (cum_inf_perc_diff_hcgs_mscape_median[1] / 100)) #57.062370694554275
cum_inf_perc_diff_hcgs_mscape_median_v2_lower  = 100 * (1 - (cum_inf_perc_diff_hcgs_mscape_median[2] / 100)) #60.98803363518759
cum_inf_perc_diff_hcgs_mscape_median_v2_upper  = 100 * (1 - (cum_inf_perc_diff_hcgs_mscape_median[3] / 100)) #50.424480802637696

cum_inf_perc_diff_hcgs_mscape_mean = mean_ci_bootstrap(vec = cum_inf_perc_diff_hcgs_mscape, n_boot=1000, alpha=0.05)
# (mean = 180.48997853178494, lower = 134.80467584999064, upper = 230.87844848470223)


####################################################################################################
# (6) Comparison of best ICU sampling performers
#     (6.3) Report 95th percentile of TD and cumulative infections at TD
####################################################################################################

# Reload ICU TD data if necessary
icu_tds_df_dict = load( joinpath( root_dir, "NBPMscape/paper/analysis/2_sampling_analysis/covid_like/geo_optimisation"
                                    , "covidlike-1.6.1_nrep1000_4300448_icu_tds_df_dict.jld2")
                        , "icu_tds_df_dict")

# 8 sites in winter

sampling_param_labels = [ "icu_tds_seasonal_ari_1440_300_sampling_sites_30blocks_8sites"
                        , "icu_tds_seasonal_ari_1440_300_sampling_sites_mscape_8sites"
                        , "icu_tds_seasonal_ari_1440_300_sampling_sites_popsize_8sites"
                        , "icu_tds_seasonal_ari_1440_300_sampling_sites_london_8sites"
                        ]

quantile_vals = [0.025,0.05,0.50,0.95,0.975]

# Create df to store data
quantile_comparison_df = DataFrame(   sampling_param_labels     = sampling_param_labels
                                    , ICU_TD_quantile_2_5       = [-1.0,-1.0,-1.0,-1.0]
                                    , ICU_TD_quantile_5_0       = [-1.0,-1.0,-1.0,-1.0]
                                    , ICU_TD_quantile_50_0      = [-1.0,-1.0,-1.0,-1.0]
                                    , ICU_TD_quantile_95_0      = [-1.0,-1.0,-1.0,-1.0]
                                    , ICU_TD_quantile_97_5      = [-1.0,-1.0,-1.0,-1.0]
                                    , cum_tinf_td_quantile_2_5  = [-1.0,-1.0,-1.0,-1.0]
                                    , cum_tinf_td_quantile_5_0  = [-1.0,-1.0,-1.0,-1.0]
                                    , cum_tinf_td_quantile_50_0 = [-1.0,-1.0,-1.0,-1.0]
                                    , cum_tinf_td_quantile_95_0 = [-1.0,-1.0,-1.0,-1.0]
                                    , cum_tinf_td_quantile_97_5 = [-1.0,-1.0,-1.0,-1.0]
                                    )

# Loop through sampling data, compute quantile values and populate df
for i in 1:length(sampling_param_labels)

    samp_param = sampling_param_labels[i]

    for j in 1:length(quantile_vals)
        quant_val = quantile_vals[j]
        quantile_comparison_df[i,j+1] = quantile( icu_tds_df_dict[samp_param][:,:ICU_TD], quant_val)
        quantile_comparison_df[i,j+6] = quantile( icu_tds_df_dict[samp_param][:,:cum_tinf_td], quant_val)
    end

end

# View quantile data
println(quantile_comparison_df.sampling_param_labels)
#["icu_tds_seasonal_ari_1440_300_sampling_sites_30blocks_8sites", "icu_tds_seasonal_ari_1440_300_sampling_sites_mscape_8sites", "icu_tds_seasonal_ari_1440_300_sampling_sites_popsize_8sites", "icu_tds_seasonal_ari_1440_300_sampling_sites_london_8sites"]

println(quantile_comparison_df)
#4×11 DataFrame
# Row │ sampling_param_labels                                           ICU_TD_quantile_2_5  ICU_TD_quantile_5_0  ICU_TD_quantile_50_0  ICU_TD_quantile_95_0  ICU_TD_quantile_97_5  cum_tinf_td_quantile_2_5  cum_tinf_td_quantile_5_0  cum_tinf_td_quantile_50_0  cum_tinf_td_quantile_95_0  cum_tinf_td_quantile_97_5 
#     │ String                                                          Float64              Float64              Float64               Float64               Float64               Float64                   Float64                   Float64                    Float64                    Float64                   
#─────┼──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#   1 │ icu_tds_seasonal_ari_1440_300_sampling_sites_30blocks_8sites    25.2472              27.6434               43.3668               57.4363               59.8399              260.825                    545.0                      6677.0                    28594.0                    35334.2
#   2 │ icu_tds_seasonal_ari_1440_300_sampling_sites_mscape_8sites      30.7218              33.7587               48.5729               63.0353               66.0667              754.05                    1348.3                     16807.0                    63160.8                    75710.9
#   3 │ icu_tds_seasonal_ari_1440_300_sampling_sites_popsize_8sites     26.363               29.5652               45.565                60.3602               62.428               374.45                     726.65                    10257.0                    40277.1                    47259.2
#   4 │ icu_tds_seasonal_ari_1440_300_sampling_sites_london_8sites      25.0915              30.4947               45.9353               60.2061               64.0725              415.875                    799.8                     10878.0                    46352.8                    56366.1

# Save to file
CSV.write( joinpath( root_dir, "NBPMscape/paper/analysis/2_sampling_analysis/covid_like/geo_optimisation/icu_tds_seasonal_ari_1440_300_8sites_best_performer_quantiles.csv"), quantile_comparison_df )


####################################################################################################
# (6) Comparison of best ICU sampling performers
#     (6.4) mSCAPE sampling equivalent to 8 site HCGS
####################################################################################################

# Comparing times to detection (TDs) resulting from different sampling strategies.
# Specifically, trying to compute the number of samples in a sampling strategy using the mSCAPE sampling sites
# required to match the median TD achieved using sites selected with the HCGS site selection strategy.
# The number of samples required to achieve this is estimated using a simple linear interpolation
# between TD values from existing simulations

# Reload ICU TD data if necessary
icu_tds_df_dict = load( joinpath( root_dir, "NBPMscape/paper/analysis/2_sampling_analysis/covid_like/geo_optimisation"
                                    , "covidlike-1.6.1_nrep1000_4300448_icu_tds_df_dict.jld2")
                        , "icu_tds_df_dict")

# Create df containing scenario we want to match the TD to and the scenarios that we will use to 
# linearly interpolate between to get a number of samples that produce an equialent TD
lin_interp_data_df = DataFrame( icu_tds_seasonal_ari_1440_300_sampling_sites_30blocks_8sites = icu_tds_df_dict["icu_tds_seasonal_ari_1440_300_sampling_sites_30blocks_8sites"][:, :ICU_TD]
                               ,icu_tds_seasonal_ari_1440_300_sampling_sites_mscape_8sites   = icu_tds_df_dict["icu_tds_seasonal_ari_1440_300_sampling_sites_mscape_8sites"  ][:, :ICU_TD]
                               ,icu_tds_seasonal_ari_1440_500_sampling_sites_mscape_8sites   = icu_tds_df_dict["icu_tds_seasonal_ari_1440_500_sampling_sites_mscape_8sites"  ][:, :ICU_TD]
                               ,icu_tds_seasonal_ari_1440_1000_sampling_sites_mscape_8sites  = icu_tds_df_dict["icu_tds_seasonal_ari_1440_1000_sampling_sites_mscape_8sites" ][:, :ICU_TD]
                               ,icu_tds_seasonal_ari_1440_300_sampling_sites_mscape_26sites  = icu_tds_df_dict["icu_tds_seasonal_ari_1440_300_sampling_sites_mscape_26sites" ][:, :ICU_TD]
                               ,icu_tds_seasonal_ari_1440_500_sampling_sites_mscape_26sites  = icu_tds_df_dict["icu_tds_seasonal_ari_1440_500_sampling_sites_mscape_26sites" ][:, :ICU_TD]
                               ,icu_tds_seasonal_ari_1440_1000_sampling_sites_mscape_26sites = icu_tds_df_dict["icu_tds_seasonal_ari_1440_1000_sampling_sites_mscape_26sites"][:, :ICU_TD]
                               )


# Helper function to perform linear interpolation and 
# return number of samples required to match the TD
"""
Function   n_sample_lin_interp

Description     Function to compute a linear interpolation between time to detection (TD)
                values for a sampling strategy (strategy 2) and return the number of samples required
                to match the TD of another strategy (strategy 1). To work, the 'comparator' TD (from strategy 1) should lie
                between the TD values from strategy 2 (lin_interp_upper and lin_interp_lower).

Arguments   comparator::Float64       TD of strategy that we want to match with the other strategy
            lin_interp_lower          TD generated from the other strategy with the lower number of samples for linear interpolation
                                      Note that 'lower' here refers to the number of samples and not the TD
            lin_interp_upper          TD generated from the other strategy with the upper number of samples for linear interpolation
                                      Note that 'upper' here refers to the number of samples and not the TD 
            n_samples_lower           The lower number of samples, used to generate the TD value in lin_interp_lower
            n_samples_upper           The upper number of samples, used to generate the TD value in lin_interp_upper

Returns     Number of samples required for strategy 2 to match the TD of strategy 1
            
Examples    n_sample_lin_interp(  comparator_val       = 45.1
                                , lin_interp_lower_val = 51.0
                                , lin_interp_upper_val = 42.1
                                , n_samples_lower      = 300
                                , n_samples_upper      = 1000  )
"""
function n_sample_lin_interp(; comparator_val::Float64       = 45.1
                             , lin_interp_lower_val::Float64 = 51.0
                             , lin_interp_upper_val::Float64 = 42.1
                             , n_samples_lower::Int64    = 300
                             , n_samples_upper::Int64    = 1000  )
       
    # Compute the change in TD per sample
    td_reduction_per_lin_interp_sample = (lin_interp_upper_val - lin_interp_lower_val) / (n_samples_upper-n_samples_lower)

    # Change in TD required from strategies being linearly interpolated to match the TD of the comparator strategy
    required_td_change = comparator_val - lin_interp_lower_val

    # Compute number of additional samples required to match the TDs
    additional_samples = required_td_change / td_reduction_per_lin_interp_sample
    # Number of samples required to match the TD
    n_samples_required = n_samples_lower + additional_samples
    return( n_samples_required )
end
    


# Helper function to perform bootstrap of linear interpolation
"""
Function   n_sample_lin_interp_boot(;df, n_boot=1000, alpha=0.05)

Description     Function to wrap bootstrapping around the linear interpolation function n_sample_lin_interp

Arguments   df::DataFrame   Dataframe containing TD values
            n_boot::Int64   Number of bootstrap repeats required
            alpha::Float64  Quantile value, i.e. 0.05 for 95% CI
            , n_boot::Int64 = 1000, alpha::Float64 = 0.05
                                , comparator::String
                                , lin_interp_lower::String
                                , lin_interp_upper::String
                                , n_samples_lower::Int64 = 300
                                , n_samples_upper::Int64 = 1000

Returns     Named tuple contaiing: 
            - n_samples_median: median of number of samples required to match the TD from actual simulation replicates
            - n_samples_median_lower: lower quantile
            - n_samples_median_upper: upper quantile
            - n_samples_median_perc_change: median percentage change (as proportion of n_samples_lower - multiply by 100 to get percentage)
            - n_samples_median_lower_perc_change: lower quantile
            - n_samples_median_upper_perc_change: upper quantile
            - boot_equiv_n_samples: vector of bootstrap sample numbers
            
Examples    # Compute median and upper and lower values for 95% CI
            # for each TD results df in a dictionary containing different scenarios
            df_example = DataFrame(  a       = [45.0,46.0,47.0,50.5,61.2,62.8]
                                   , b_lower = [51.0,52.0,53.0,53.5,53.8,54.2]
                                   , b_upper = [40.0,41.0,42.0,42.1,42.3,42.8])
            n_sample_lin_interp_boot( df = df_example
                                    , n_boot = 1000
                                    , alpha = 0.05
                                    , comparator = "a"
                                    , lin_interp_lower = "b_lower"
                                    , lin_interp_upper = "b_upper"
                                    , n_samples_lower = 300
                                    , n_samples_upper = 1000
                                )
            
"""
function n_sample_lin_interp_boot(; df::DataFrame, n_boot::Int64 = 1000, alpha::Float64 = 0.05
                                , comparator::String
                                , lin_interp_lower::String
                                , lin_interp_upper::String
                                , n_samples_lower::Int64 = 300
                                , n_samples_upper::Int64 = 1000
                                )
    
    boot_equiv_n_samples = Vector{Float64}(undef, n_boot) #size(df,1))

    for i in 1:1000 #i=1

        # Draw df row indices for bootstrapping 
        # - using indices ensures that using the same group of sim reps for the linear interpolation
        bootstap_indices = rand(1:size(df,1), 1000) # rand()
        
        # Compute median TD for comparator (e.g. HCGS sampling TDs) from bootstrap samples
        boot_comparator_median = median( df[ bootstap_indices, comparator] )
        
        # Compute median TDs for upper and lower sampling strategies (e.g. mSCAPE sampling TDs) from bootstrap samples
        # for linear interpolation
        boot_lin_interp_lower_median = median( df[ bootstap_indices, lin_interp_lower] )
        boot_lin_interp_upper_median = median( df[ bootstap_indices, lin_interp_upper] )
        
        # Compute the actual number of samples for the linearly intperolated strategy (e.g. mSCAPE)
        # required to match the TD of the comparator strategy (e.g. HCGS)
        # and add to bootstrapping list 
        boot_equiv_n_samples[i] = n_sample_lin_interp(  comparator_val       = boot_comparator_median
                                                      , lin_interp_lower_val = boot_lin_interp_lower_median
                                                      , lin_interp_upper_val = boot_lin_interp_upper_median
                                                      , n_samples_lower      = n_samples_lower
                                                      , n_samples_upper      = n_samples_upper  )
    end

    # Compute median value from the actual simulation replicates rather than the bootstrapped values
    n_samples_median = n_sample_lin_interp(  comparator_val       = median( df[ :, comparator]  )
                                           , lin_interp_lower_val = median( df[ :, lin_interp_lower] )
                                           , lin_interp_upper_val = median( df[ :, lin_interp_upper] )
                                           , n_samples_lower      = n_samples_lower
                                           , n_samples_upper      = n_samples_upper  ) 

    # Compute median and bootstrap lower and upper estimates of median
    #n_samples_median = median(equiv_n_samples)
    n_samples_median_lower = quantile(boot_equiv_n_samples, alpha/2)
    n_samples_median_upper = quantile(boot_equiv_n_samples, 1 - alpha/2)
    
    # Compute percentage change required from lower number of samples in linear interpolation
    n_samples_median_perc_change       = n_samples_median       / n_samples_lower - 1 
    n_samples_median_lower_perc_change = n_samples_median_lower / n_samples_lower - 1
    n_samples_median_upper_perc_change = n_samples_median_upper / n_samples_lower - 1
    
    # Plot median values from bootstrapping
    #p1 = Plots.plot( boot_equiv_n_samples )
    #p2 = Plots.histogram( boot_equiv_n_samples )
    #p3 = Plots.plot( size=(800, 1200), p1, p2, layout = (2,1))

    return ( n_samples_median = n_samples_median
           , n_samples_median_lower = n_samples_median_lower
           , n_samples_median_upper = n_samples_median_upper
           , n_samples_median_perc_change = n_samples_median_perc_change
           , n_samples_median_lower_perc_change = n_samples_median_lower_perc_change
           , n_samples_median_upper_perc_change = n_samples_median_upper_perc_change
           , boot_equiv_n_samples )

end


# With 8 mSCAPE sites
# Compare and linear interpolation:
# - HCGS 30 blocks 8 sites in winter 300 samples
# against
# - mSCAPE 8 sites in winter 300 samples
# - mSCAPE 8 sites in winter 1000 samples
lin_interp_hcgs_v_mscape_8sites_300_1000 = n_sample_lin_interp_boot(  
                                                    df = lin_interp_data_df
                                                    , n_boot = 1000
                                                    , alpha = 0.05
                                                    , comparator = "icu_tds_seasonal_ari_1440_300_sampling_sites_30blocks_8sites"
                                                    , lin_interp_lower = "icu_tds_seasonal_ari_1440_300_sampling_sites_mscape_8sites"
                                                    , lin_interp_upper = "icu_tds_seasonal_ari_1440_1000_sampling_sites_mscape_8sites"
                                                    , n_samples_lower = 300
                                                    , n_samples_upper = 1000
                                                )

# ( n_samples_median                   =  6659.189745876947
# , n_samples_median_lower             =-46255.588067792465
# , n_samples_median_upper             = 32629.026738008833
# , n_samples_median_perc_change       =    21.197299152923158
# , n_samples_median_lower_perc_change =  -155.1852935593082
# , n_samples_median_upper_perc_change =   107.76342246002945
# , boot_equiv_n_samples = [9020.61905373689, 8191.078353889993, 19114.608973920484, 6343.24891568868, 3748.3742198949108, 2796.1661009064233, 3447.500678275232, 2948.577904734735, 4334.3100286026865, 6027.8180744490755  …  9677.194005166388, -166308.75595127261, 7207.889902023792, 6013.759971520789, 7131.655135016546, 14061.920256799718, 2723.8930724317574, 5284.218063056961, 10565.40951524039, 8865.752881013685])

# Check results
Plots.plot(lin_interp_hcgs_v_mscape_8sites_300_1000.boot_equiv_n_samples)
Plots.histogram(lin_interp_hcgs_v_mscape_8sites_300_1000.boot_equiv_n_samples)

# Conclusion: In most simulation replicates the TD achieved with HCGS (8 sites, 300 samples) cannot be achieved
#             by mSCAPE with 8 sites even when increasing the number of samples to 1000.
#             The results are therefore nonsensical.

# With 26 mSCAPE sites and max range of samples
# Compare and linear interpolation:
# - HCGS 30 blocks 8 sites in winter 300 samples
# against
# - mSCAPE 26 SITES in winter 300 samples
# - mSCAPE 26 SITES in winter 1000 samples
lin_interp_hcgs_v_mscape_26sites_300_1000 = n_sample_lin_interp_boot(  
                                          df = lin_interp_data_df
                                        , n_boot = 1000
                                        , alpha = 0.05
                                        , comparator = "icu_tds_seasonal_ari_1440_300_sampling_sites_30blocks_8sites"
                                        , lin_interp_lower = "icu_tds_seasonal_ari_1440_300_sampling_sites_mscape_26sites"
                                        , lin_interp_upper = "icu_tds_seasonal_ari_1440_1000_sampling_sites_mscape_26sites"
                                        , n_samples_lower = 300
                                        , n_samples_upper = 1000
                                    )
# ( n_samples_median                   = 501.24000695993936
# , n_samples_median_lower             = 428.16461537760506
# , n_samples_median_upper             = 583.7718762614953
# , n_samples_median_perc_change       =   0.6708000231997979
# , n_samples_median_lower_perc_change =   0.4272153845920168
# , n_samples_median_upper_perc_change =   0.9459062542049843
# , boot_equiv_n_samples = [458.0137484986292, 530.1727889923766, 479.04990069481096, 501.6686942662915, 503.2968532619641, 499.85773671283215, 511.4523510644532, 497.15842071844315, 522.1521333706016, 488.89100092770536  …  485.37877292700955, 524.4337164506245, 455.72330037147094, 489.69466917794614, 509.36307637813127, 459.3016800571397, 489.7712107887529, 493.1531515534009, 445.1441188973847, 437.3413528013116])

Plots.plot(lin_interp_hcgs_v_mscape_26sites_300_1000.boot_equiv_n_samples)
Plots.histogram(lin_interp_hcgs_v_mscape_26sites_300_1000.boot_equiv_n_samples)

# Conclusion: mSCAPE sampling sites would achieve an equivalent median TD to HCGS 30 blocks 8 sites 
#             by increasing to 26 sites and increasing the number of samples by 
#             67% (43% to 95%, 95% estimate using 1000 bootstraps from 1000 simulation replicates)


# With 26 mSCAPE sites but try linear interpolation between 300 and 500 samples rather than 300 and 1000
# Compare and linear interpolation:
# - HCGS 30 blocks 8 sites in winter 300 samples
# against
# - mSCAPE 26 SITES in winter 300 samples
# - mSCAPE 26 SITES in winter 500 samples
lin_interp_hcgs_v_mscape_26sites_300_500 = n_sample_lin_interp_boot(  
                                          df = lin_interp_data_df
                                        , n_boot = 1000
                                        , alpha = 0.05
                                        , comparator = "icu_tds_seasonal_ari_1440_300_sampling_sites_30blocks_8sites"
                                        , lin_interp_lower = "icu_tds_seasonal_ari_1440_300_sampling_sites_mscape_26sites"
                                        , lin_interp_upper = "icu_tds_seasonal_ari_1440_500_sampling_sites_mscape_26sites"
                                        , n_samples_lower = 300
                                        , n_samples_upper = 500
                                    )
# ( n_samples_median = 413.72321218205127
# , n_samples_median_lower = 377.5253338306468
# , n_samples_median_upper = 466.3528370870019
# , n_samples_median_perc_change = 0.37907737394017094
# , n_samples_median_lower_perc_change = 0.2584177794354894
# , n_samples_median_upper_perc_change = 0.5545094569566729
# , boot_equiv_n_samples = [421.8448063064184, 434.1433268718867, 429.7078461413568, 402.660763248125, 434.31930486480906, 434.2324991265002, 414.20690531212756, 447.25218790049587, 393.8692726421558, 414.0900118152633  …  419.8105485984143, 445.7703594757387, 419.5111477479917, 444.662436651001, 409.1596636382099, 412.25143433160883, 390.53531938159296, 469.0427483240585, 449.5109326308735, 376.3215584577712])

Plots.plot(lin_interp_hcgs_v_mscape_26sites_300_500.boot_equiv_n_samples)
Plots.histogram(lin_interp_hcgs_v_mscape_26sites_300_500.boot_equiv_n_samples)

# Conclusion: mSCAPE sampling sites would achieve an equivalent median TD to HCGS 30 blocks 8 sites 
#             by increasing to 26 sites and increasing the number of samples by 
#             38% (26% to 55%, 95% estimate using 1000 bootstraps from 1000 simulation replicates)
#             Using the difference in TD between 300 and 500 samples is probably better because at
#             1000 samples per week we're probably sampling everyone (or actually sampling less than 1000 ARI patients)

# Repeat but using the biggest TD difference per sample between 300 to 500 samples, 500 to 1000, and 300 to 1000
# This should be the most conservative estimate of the increase in the number of samples required
equiv_n_samples_cons = Vector{Float64}(undef, size(lin_interp_data_df,1))

for i in 1:1000 #i=1
    bootstap_indices = rand(1:size(lin_interp_data_df,1), 1000)
    
    median_td_hcgs_8_300 = median( lin_interp_data_df[bootstap_indices,"icu_tds_seasonal_ari_1440_300_sampling_sites_30blocks_8sites"]  )
    
    median_td_mscape_26_300 = median( lin_interp_data_df[bootstap_indices,"icu_tds_seasonal_ari_1440_300_sampling_sites_mscape_26sites"]  )
    median_td_mscape_26_500 = median( lin_interp_data_df[bootstap_indices,"icu_tds_seasonal_ari_1440_500_sampling_sites_mscape_26sites"]  )
    median_td_mscape_26_1000 = median( lin_interp_data_df[bootstap_indices,"icu_tds_seasonal_ari_1440_1000_sampling_sites_mscape_26sites"]  )
    
    td_reduction_per_mscape_sample_300_500 = (median_td_mscape_26_500 - median_td_mscape_26_300) / (500-300)
    td_reduction_per_mscape_sample_500_1000 = (median_td_mscape_26_1000 - median_td_mscape_26_500) / (1000-500)
    td_reduction_per_mscape_sample_300_1000 = (median_td_mscape_26_1000 - median_td_mscape_26_300) / (1000-300)
    
    td_reduction_per_mscape_sample = min( td_reduction_per_mscape_sample_300_500
                                        , td_reduction_per_mscape_sample_500_1000
                                        , td_reduction_per_mscape_sample_300_1000)

    required_td_change = median_td_hcgs_8_300 - median_td_mscape_26_300

    additional_samples = required_td_change / td_reduction_per_mscape_sample 

    equiv_n_samples_cons[i] = 300 + additional_samples 

end

median(equiv_n_samples_cons)            # 417.1578318871384
quantile(equiv_n_samples_cons, 0.025)   # 372.66640929575146
quantile(equiv_n_samples_cons, 0.975)   # 462.5666216257017

median(equiv_n_samples_cons) / 300 -1           # 0.39052610629046125
quantile(equiv_n_samples_cons, 0.025) / 300 -1  # 0.24222136431917152
quantile(equiv_n_samples_cons, 0.975) / 300 -1  # 0.541888738752339

Plots.plot(equiv_n_samples_cons)
Plots.histogram(equiv_n_samples_cons)
Plots.histogram!(lin_interp_hcgs_v_mscape_26sites_300_500.boot_equiv_n_samples, alpha=0.5)

# This result is basically the same as just using the difference in TD between 300 and 500 samples
# (which is reassuring)


####################################################################################################
# (6) Comparison of best ICU sampling performers
#     (6.5) Uncertainty in variance
####################################################################################################

# Reload ICU TD data if necessary
icu_tds_df_dict = load( joinpath( root_dir, "NBPMscape/paper/analysis/2_sampling_analysis/covid_like/geo_optimisation"
                                    , "covidlike-1.6.1_nrep1000_4300448_icu_tds_df_dict.jld2")
                        , "icu_tds_df_dict")

# HCGS variance                        
var_hcgs = NBPMscape.var_ci_bootstrap( vec = icu_tds_df_dict["icu_tds_seasonal_ari_1440_300_sampling_sites_30blocks_8sites"][:,:ICU_TD]
                            , n_boot = 1000
                            , alpha = 0.05
                            )
#(variance = 78.6815120712962, lower = 71.46367046017664, upper = 86.3039145461815)

# mSCAPE variance
var_mscape = NBPMscape.var_ci_bootstrap( vec = icu_tds_df_dict["icu_tds_seasonal_ari_1440_300_sampling_sites_mscape_8sites"][:,:ICU_TD]
                            , n_boot = 1000
                            , alpha = 0.05
                            )
# (variance = 78.16915082234677, lower = 71.94321516668384, upper = 85.26995040660658)


####################################################################################################
# (6) Comparison of best ICU sampling performers
#       (6.6) Additional comparisons
####################################################################################################

# Reload ICU TD data if necessary
icu_tds_df_dict = load( joinpath( root_dir, "NBPMscape/paper/analysis/2_sampling_analysis/covid_like/geo_optimisation"
                                    , "covidlike-1.6.1_nrep1000_4300448_icu_tds_df_dict.jld2")
                        , "icu_tds_df_dict")

# TD vectors
td_winter_8sites_300samples_hcgs30b = icu_tds_df_dict["icu_tds_seasonal_ari_1440_300_sampling_sites_30blocks_8sites"][:,   :ICU_TD]                        
td_winter_8sites_300samples_mscape  = icu_tds_df_dict["icu_tds_seasonal_ari_1440_300_sampling_sites_mscape_8sites"][:,   :ICU_TD]                        

# TD absolute difference
td_diff_8sites_winter_300 = td_winter_8sites_300samples_hcgs30b .- td_winter_8sites_300samples_mscape
median_ci_bootstrap( vec = td_diff_8sites_winter_300, n_boot = 1000, alpha = 0.05 )
#(median = -4.6887511817578424, lower = -5.228863560980379, upper = -4.019035933996997)

# TD percentage difference
td_perce_diff_8sites_winter_300 = (td_winter_8sites_300samples_hcgs30b ./ td_winter_8sites_300samples_mscape) .-1
median_ci_bootstrap( vec = td_perce_diff_8sites_winter_300, n_boot = 1000, alpha = 0.05 )
#(median = -0.09717554636508147, lower = -0.10782416704065527, upper = -0.08038188286753256)


# Cum Inf vectors
cum_inf_winter_8sites_300samples_hcgs30b = icu_tds_df_dict["icu_tds_seasonal_ari_1440_300_sampling_sites_30blocks_8sites"][:, :cum_tinf_td]                        
cum_inf_winter_8sites_300samples_mscape  = icu_tds_df_dict["icu_tds_seasonal_ari_1440_300_sampling_sites_mscape_8sites"][:,   :cum_tinf_td]                        

# Cum inf absolute difference
cum_inf_diff_8sites_winter_300 = cum_inf_winter_8sites_300samples_hcgs30b .- cum_inf_winter_8sites_300samples_mscape
median_ci_bootstrap( vec = cum_inf_diff_8sites_winter_300, n_boot = 1000, alpha = 0.05 )
#(median = -7901.0, lower = -9370.474999999999, upper = -6500.712500000001)

# Cum inf perc difference
cum_inf_perc_diff_8sites_winter_300 = (cum_inf_winter_8sites_300samples_hcgs30b ./ cum_inf_winter_8sites_300samples_mscape) .-1
median_ci_bootstrap( vec = cum_inf_perc_diff_8sites_winter_300, n_boot = 1000, alpha = 0.05 )
#(median = -0.5706237069455427, lower = -0.6098808315079745, upper = -0.5039042292127414)


# Extract from Discussion:
# "However, in a summer scenario with the same parameters, and in both winter and summer scenarios using 26 sites and 300 samples per week,
#  the advantage of the partially optimised system is seen across mean TD (4-day advantage, 2-day, 1-day, respectively),
#  mean cumulative infections (-49\%, -29\%, -20\%), and variance of TD (-4.8\%, -12.3\%, -3.5\%).


# TD vectors
td_summer_8sites_300samples_hcgs30b  = icu_tds_df_dict["icu_tds_seasonal_ari_793_300_sampling_sites_30blocks_8sites"][:,   :ICU_TD]
td_summer_26sites_300samples_hcgs20b = icu_tds_df_dict["icu_tds_seasonal_ari_793_300_sampling_sites_20blocks_26sites"][:,  :ICU_TD]
td_winter_26sites_300samples_hcgs20b = icu_tds_df_dict["icu_tds_seasonal_ari_1440_300_sampling_sites_20blocks_26sites"][:, :ICU_TD]

td_summer_8sites_300samples_mscape  = icu_tds_df_dict["icu_tds_seasonal_ari_793_300_sampling_sites_mscape_8sites"][:,   :ICU_TD]
td_summer_26sites_300samples_mscape = icu_tds_df_dict["icu_tds_seasonal_ari_793_300_sampling_sites_mscape_26sites"][:,  :ICU_TD]
td_winter_26sites_300samples_mscape = icu_tds_df_dict["icu_tds_seasonal_ari_1440_300_sampling_sites_mscape_26sites"][:, :ICU_TD]

# cumulative infections vectors
cuminf_summer_8sites_300samples_hcgs30b  = icu_tds_df_dict["icu_tds_seasonal_ari_793_300_sampling_sites_30blocks_8sites"][:,   :cum_tinf_td]
cuminf_summer_26sites_300samples_hcgs20b = icu_tds_df_dict["icu_tds_seasonal_ari_793_300_sampling_sites_20blocks_26sites"][:,  :cum_tinf_td]
cuminf_winter_26sites_300samples_hcgs20b = icu_tds_df_dict["icu_tds_seasonal_ari_1440_300_sampling_sites_20blocks_26sites"][:, :cum_tinf_td]

cuminf_summer_8sites_300samples_mscape  = icu_tds_df_dict["icu_tds_seasonal_ari_793_300_sampling_sites_mscape_8sites"][:,   :cum_tinf_td]
cuminf_summer_26sites_300samples_mscape = icu_tds_df_dict["icu_tds_seasonal_ari_793_300_sampling_sites_mscape_26sites"][:,  :cum_tinf_td]
cuminf_winter_26sites_300samples_mscape = icu_tds_df_dict["icu_tds_seasonal_ari_1440_300_sampling_sites_mscape_26sites"][:, :cum_tinf_td]

## 8 sites, summer, 300 samples per week - HCGS 30 blocks vs mSCAPE

# TD absolute difference
td_diff_8sites_summer_300 = td_summer_8sites_300samples_hcgs30b .- td_summer_8sites_300samples_mscape
median_ci_bootstrap( vec = td_diff_8sites_summer_300, n_boot = 1000, alpha = 0.05 )
#(median = -3.832682948758716, lower = -4.566220125624945, upper = -3.1450343680964243)

# TD variance
var_bs_summer_8sites_hcgs30b = NBPMscape.var_ci_bootstrap( vec = td_summer_8sites_300samples_hcgs30b, n_boot = 1000, alpha = 0.05 )
# (variance = 80.54403961367801, lower = 73.36421148436223, upper = 87.99920163722801)
var_bs_summer_8sites_mscape = NBPMscape.var_ci_bootstrap( vec = td_summer_8sites_300samples_mscape, n_boot = 1000, alpha = 0.05 )
#(variance = 84.61984770965977, lower = 77.29819195374255, upper = 93.58969061006256)
100*(var_bs_summer_8sites_hcgs30b[1] / var_bs_summer_8sites_mscape[1] -1) # -4.816610058158355

# Cumulative infections at TD difference ( hcgs - mscape )
cum_inf_diff_8sites_summer_300 = cuminf_summer_8sites_300samples_hcgs30b .- cuminf_summer_8sites_300samples_mscape
median_ci_bootstrap( vec = cum_inf_diff_8sites_summer_300, n_boot = 1000, alpha = 0.05 )
#(median = -5705.0, lower = -6999.162499999999, upper = -4538.0)

# Cumulative infections at TD percentage difference ( hcgs / mscape - 1)
cum_inf_perc_diff_8sites_summer_300 = ( cuminf_summer_8sites_300samples_hcgs30b ./ cuminf_summer_8sites_300samples_mscape) .- 1
median_ci_bootstrap( vec = cum_inf_perc_diff_8sites_summer_300, n_boot = 1000, alpha = 0.05 )
#(median = -0.49001662133592516, lower = -0.575554303901723, upper = -0.40810071590812236)

# FOR 26 sites we switch the HCGS strategy from 30 blocks to 20 blocks 
# as this produces better performance in terms of TD and cumulative infections at TD 
# when considering both winter and summer scenarios

results_df_winter_26sites_300samples = filter(row -> contains_all(row.parameter_set_name, ["1440","26sites","300"]), results_df)
results_df_winter_26sites_300samples = sort(results_df_winter_26sites_300samples, "TD_median")
results_df_winter_26sites_300samples[:,2:size(results_df_winter_26sites_300samples,2)] = round.( results_df_winter_26sites_300samples[:,2:size(results_df_winter_26sites_300samples,2)], digits=1)
println(results_df_winter_26sites_300samples.parameter_set_name)
#["icu_tds_seasonal_ari_1440_300_sampling_sites_20blocks_26sites", "icu_tds_seasonal_ari_1440_300_sampling_sites_10blocks_26sites", "icu_tds_seasonal_ari_1440_300_sampling_sites_5blocks_26sites", "icu_tds_seasonal_ari_1440_300_sampling_sites_30blocks_26sites", "icu_tds_seasonal_ari_1440_300_sampling_sites_popsize_26sites", "icu_tds_seasonal_ari_1440_300_sampling_sites_mscape_26sites"]
#["", "", "", "", "", ""]
println(results_df_winter_26sites_300samples)

#6×19 DataFrame
# Row │ parameter_set_name                                                  TD_median  iqr_lower  iqr_upper  TD_mean  Var      SD       U_01     U_05     U_1      TD3_median  iqr3_lower  iqr3_upper  TD3_mean  Var3     SD3      U3_01    U3_05    U3_1    
#     │ String                                                              Float64    Float64    Float64    Float64  Float64  Float64  Float64  Float64  Float64  Float64     Float64     Float64     Float64   Float64  Float64  Float64  Float64  Float64 
#─────┼──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#   1 │ icu_tds_seasonal_ari_1440_300_sampling_sites_20blocks_26sites       42.3       37.0       47.6       42.1     72.2      8.5     49.6     78.4     114.5    49.6        45.5        54.5        50.2     47.0      6.9      54.3     73.1     96.6
#   2 │ icu_tds_seasonal_ari_1440_300_sampling_sites_10blocks_26sites       42.9       37.6       49.5       43.1     76.8      8.8     50.5     81.3     119.7    50.7        46.0        55.6        51.1     44.9      6.7      55.2     73.2     95.6
#   3 │ icu_tds_seasonal_ari_1440_300_sampling_sites_5blocks_26sites        43.0       37.0       48.9       42.6     79.7      8.9     50.9     82.8     122.7    50.4        46.4        55.5        51.0     45.2      6.7      54.9     73.0     95.6
#   4 │ icu_tds_seasonal_ari_1440_300_sampling_sites_30blocks_26sites       43.0       37.6       48.5       42.9     73.8      8.6     50.4     79.9     116.8    50.3        45.8        55.1        50.8     45.0      6.7      54.8     72.8     95.2
#   5 │ icu_tds_seasonal_ari_1440_300_sampling_sites_popsize_26sites        43.3       37.5       48.8       43.0     74.0      8.6     50.7     80.3     117.3    50.8        46.4        55.6        51.2     44.5      6.7      55.3     73.1     95.3
#   6 │ icu_tds_seasonal_ari_1440_300_sampling_sites_mscape_26sites         44.7       39.5       50.6       44.7     82.3      9.1     52.9     85.8     127.0    52.3        48.2        57.2        52.9     46.5      6.8      56.9     75.5     98.8

# Conclusion: In this scenario (winter sampling ICU at 300 samples per week across 26 sites/NHS Trusts) the best performing strategy is still HCGS.
#             The 20 block strategy has the earliets median time to detect 1 case (TD) and 3 cases (TD3) as well as the lowest variance in TD, although
#             it has a high variance in TD3.

## 26 sites, summer, 300 samples per week - HCGS 20 blocks vs mSCAPE
# TD absolute difference
td_diff_26sites_summer_300 = td_summer_26sites_300samples_hcgs20b .- td_summer_26sites_300samples_mscape
median_ci_bootstrap( vec = td_diff_26sites_summer_300, n_boot = 1000, alpha = 0.05 )
#(median = -1.1551298713462792, lower = -1.536559971829277, upper = -0.7069740028819531)

# TD variance
var_bs_summer_26sites_hcgs20b = NBPMscape.var_ci_bootstrap( vec = td_summer_26sites_300samples_hcgs20b, n_boot = 1000, alpha = 0.05 )
# (variance = 76.42422575745216, lower = 68.91010810871823, upper = 83.44347611958679)
var_bs_summer_26sites_mscape = NBPMscape.var_ci_bootstrap( vec = td_summer_26sites_300samples_mscape, n_boot = 1000, alpha = 0.05 )
# (variance = 79.18364232938477, lower = 71.75249416190358, upper = 86.63962428677797)
100*(var_bs_summer_26sites_hcgs20b[1] / var_bs_summer_26sites_mscape[1] -1) # -3.4848315772771588

# Cumulative infections at TD percentage difference ( hcgs / mscape - 1)
cum_inf_perc_diff_26sites_summer_300 = ( cuminf_summer_26sites_300samples_hcgs20b ./ cuminf_summer_26sites_300samples_mscape) .- 1
median_ci_bootstrap( vec = cum_inf_perc_diff_26sites_summer_300, n_boot = 1000, alpha = 0.05 )
#(median = -0.20030384287753905, lower = -0.25402747566294037, upper = -0.13216930475942723)

## 26 sites, winter, 300 samples per week - HCGS 30 blocks vs mSCAPE
# TD absolute difference
td_diff_26sites_winter_300 = td_winter_26sites_300samples_hcgs20b .- td_winter_26sites_300samples_mscape
median_ci_bootstrap( vec = td_diff_26sites_winter_300, n_boot = 1000, alpha = 0.05 )
#(median = -1.861840917347262, lower = -2.2415725976647636, upper = -1.3305794229405465)

# TD variance
var_bs_winter_26sites_hcgs20b = NBPMscape.var_ci_bootstrap( vec = td_winter_26sites_300samples_hcgs20b, n_boot = 1000, alpha = 0.05 )
# (variance = 72.19348724077926, lower = 65.06570217789239, upper = 79.66585364595427)
var_bs_winter_26sites_mscape = NBPMscape.var_ci_bootstrap( vec = td_winter_26sites_300samples_mscape, n_boot = 1000, alpha = 0.05 )
#(variance = 82.30679151988159, lower = 74.88215522837405, upper = 89.60977609811084)
100*(var_bs_winter_26sites_hcgs20b[1] / var_bs_winter_26sites_mscape[1] -1) # -12.287326589154446

# Cumulative infections at TD percentage difference ( hcgs / mscape - 1)
cum_inf_perc_diff_26sites_winter_300 = ( cuminf_winter_26sites_300samples_hcgs20b ./ cuminf_winter_26sites_300samples_mscape) .- 1
median_ci_bootstrap( vec = cum_inf_perc_diff_26sites_winter_300, n_boot = 1000, alpha = 0.05 )
#(median = -0.28859638186822545, lower = -0.3509732020449274, upper = -0.2233913914519754)

####################################################################################################
# (6) Comparison of best ICU sampling performers
#     (6.6) Median cumulative infections
####################################################################################################

# Reload ICU TD data if necessary
icu_tds_df_dict = load( joinpath( root_dir, "NBPMscape/paper/analysis/2_sampling_analysis/covid_like/geo_optimisation"
                                    , "covidlike-1.6.1_nrep1000_4300448_icu_tds_df_dict.jld2")
                        , "icu_tds_df_dict")

# Winter, 300 samples per week, 8 sites #

# Filter ICU TD results dictionary for specific sampling scenarios
filter_params = ["1440","300","8sites"]
winter_300s_8sites =  Dict( k => v for (k, v) in icu_tds_df_dict
             if all(s -> occursin(s, k), filter_params)  # Check all substrings in key
)

for k in keys(winter_300s_8sites)
    println("$(k)   $(median_ci_bootstrap( vec = icu_tds_df_dict[k][:, :cum_tinf_td], n_boot = 1000, alpha = 0.05 ))")
end

#icu_tds_seasonal_ari_1440_300_sampling_sites_10blocks_8sites   (median =  8978.0, lower =  8227.2125, upper =  9798.0              )
#icu_tds_seasonal_ari_1440_300_sampling_sites_20blocks_8sites   (median =  6318.0, lower =  5844.125,  upper =  7087.099999999999   )
#icu_tds_seasonal_ari_1440_300_sampling_sites_london_8sites     (median = 10878.0, lower = 10083.0,    upper = 11831.25             )
#icu_tds_seasonal_ari_1440_300_sampling_sites_30blocks_8sites   (median =  6677.0, lower =  6129.0,    upper =  7129.537499999999   )
#icu_tds_seasonal_ari_1440_300_sampling_sites_mscape_8sites     (median = 16807.0, lower = 15247.925,  upper = 18448.5375           )
#icu_tds_seasonal_ari_1440_300_sampling_sites_5blocks_8sites    (median = 11524.0, lower = 10571.5,    upper = 12650.137499999999   )
#icu_tds_seasonal_ari_1440_300_sampling_sites_popsize_8sites    (median = 10257.0, lower =  9036.175,  upper = 11088.0              )

# Summer, 300 samples per week, 8 sites #

# Filter ICU TD results dictionary for specific sampling scenarios
filter_params = ["793","300","8sites"]
summer_300s_8sites =  Dict( k => v for (k, v) in icu_tds_df_dict
                            if all(s -> occursin(s, k), filter_params)  # Check all substrings in key
                        )

for k in keys( summer_300s_8sites )
    println("$(k)   $(median_ci_bootstrap( vec = icu_tds_df_dict[k][:, :cum_tinf_td], n_boot = 1000, alpha = 0.05 ))")
end
#icu_tds_seasonal_ari_793_300_sampling_sites_30blocks_8sites (median =  6626.5, lower =  6103.9875, upper =  7341.0     )
#icu_tds_seasonal_ari_793_300_sampling_sites_20blocks_8sites (median =  6736.0, lower =  6322.9375, upper =  7185.0     )
#icu_tds_seasonal_ari_793_300_sampling_sites_popsize_8sites  (median = 10736.5, lower =  9883.0,    upper = 11817.0     )
#icu_tds_seasonal_ari_793_300_sampling_sites_5blocks_8sites  (median = 10898.0, lower =  9725.15,   upper = 12043.2875  )
#icu_tds_seasonal_ari_793_300_sampling_sites_mscape_8sites   (median = 14799.5, lower = 13238.2375, upper = 16053.675   )
#icu_tds_seasonal_ari_793_300_sampling_sites_10blocks_8sites (median =  9871.0, lower =  8875.9125, upper = 10607.15    )
#icu_tds_seasonal_ari_793_300_sampling_sites_london_8sites   (median = 11625.5, lower = 10614.6,    upper = 12490.0     )


# Winter, 26 sites, 300 samples #                      
# Filter ICU TD results dictionary for specific sampling scenarios
filter_params = ["1440","300","26sites"]
winter_300s_26sites =  Dict( k => v for (k, v) in icu_tds_df_dict
                            if all(s -> occursin(s, k), filter_params)  # Check all substrings in key
                        )

for k in keys( winter_300s_26sites )
    println("$(k)   $(median_ci_bootstrap( vec = icu_tds_df_dict[k][:, :cum_tinf_td], n_boot = 1000, alpha = 0.05 ))")
end

# icu_tds_seasonal_ari_1440_300_sampling_sites_5blocks_26sites   (median = 5931.5, lower = 5395.0,    upper = 6462.0            )
# icu_tds_seasonal_ari_1440_300_sampling_sites_10blocks_26sites  (median = 6382.5, lower = 5957.9375, upper = 6934.2125         )
# icu_tds_seasonal_ari_1440_300_sampling_sites_30blocks_26sites  (median = 6289.0, lower = 5777.0,    upper = 7003.0            )
# icu_tds_seasonal_ari_1440_300_sampling_sites_20blocks_26sites  (median = 5288.0, lower = 4792.7625, upper = 5821.299999999999 )
# icu_tds_seasonal_ari_1440_300_sampling_sites_mscape_26sites    (median = 8743.5, lower = 8034.9875, upper = 9765.062499999998 )
# icu_tds_seasonal_ari_1440_300_sampling_sites_popsize_26sites   (median = 6511.5, lower = 5858.5,    upper = 7061.275          )


# Summer, 26 sites, 300 samples #
# Filter ICU TD results dictionary for specific sampling scenarios
filter_params = ["793","300","26sites"]
summer_300s_26sites =  Dict( k => v for (k, v) in icu_tds_df_dict
                            if all(s -> occursin(s, k), filter_params)  # Check all substrings in key
                        )

for k in keys( summer_300s_26sites )
    println("$(k)   $(median_ci_bootstrap( vec = icu_tds_df_dict[k][:, :cum_tinf_td], n_boot = 1000, alpha = 0.05 ))")
end

# icu_tds_seasonal_ari_793_300_sampling_sites_20blocks_26sites (median = 3529.5, lower = 3306.975,           upper = 3832.075           )
# icu_tds_seasonal_ari_793_300_sampling_sites_30blocks_26sites (median = 3458.0, lower = 3235.8999999999996, upper = 3790.0             )
# icu_tds_seasonal_ari_793_300_sampling_sites_popsize_26sites  (median = 3949.5, lower = 3498.875,           upper = 4279.5             )
# icu_tds_seasonal_ari_793_300_sampling_sites_mscape_26sites   (median = 5128.0, lower = 4669.474999999999,  upper = 5741.299999999999  )
# icu_tds_seasonal_ari_793_300_sampling_sites_10blocks_26sites (median = 3686.5, lower = 3342.0,             upper = 4142.6875          )
# icu_tds_seasonal_ari_793_300_sampling_sites_5blocks_26sites  (median = 3979.5, lower = 3591.0,             upper = 4265.5             )


####################################################################################################
# (7) Plot TD variance vs TD 
#   (7.1) for all sampling scenarios and strategies
####################################################################################################

# Plot all scenarios: Variance against TD

# Helper function for filtering parameter names for multiple substrings
function contains_all(s::AbstractString, subs::Vector{String})
    all(sub -> occursin(sub, s), subs)
end

# Reload results_df if necessary
results_df = CSV.read( joinpath( root_dir, "scripts/paper/2_sampling_analysis/covid_like/1717144_1719024_1719029_analysis/geo_optimisation"
                    ,"icu_td_results_df.csv" )
                    , DataFrame )


## Create separate datasets filtering results by sampling strategy, season and number of sites

# HCGS
# winter, 8 sites 
hcgs_winter_8sites       = filter(row -> contains_all(row.parameter_set_name, ["blocks","1440","8sites"]), results_df)
hcgs_winter_8sites_300s  = filter(row -> contains_all(row.parameter_set_name, ["blocks","1440","8sites","300"]),  results_df)
hcgs_winter_8sites_500s  = filter(row -> contains_all(row.parameter_set_name, ["blocks","1440","8sites","500"]),  results_df)
hcgs_winter_8sites_1000s = filter(row -> contains_all(row.parameter_set_name, ["blocks","1440","8sites","1000"]), results_df)
# summer, 8 sites
hcgs_summer_8sites       = filter(row -> contains_all(row.parameter_set_name, ["blocks","793","8sites"]), results_df)
hcgs_summer_8sites_300s  = filter(row -> contains_all(row.parameter_set_name, ["blocks","793","8sites","300"]), results_df)
hcgs_summer_8sites_500s  = filter(row -> contains_all(row.parameter_set_name, ["blocks","793","8sites","500"]), results_df)
hcgs_summer_8sites_1000s = filter(row -> contains_all(row.parameter_set_name, ["blocks","793","8sites","1000"]), results_df)
# winter, 26 sites
hcgs_winter_26sites       = filter(row -> contains_all(row.parameter_set_name, ["blocks","1440","26sites"]), results_df)
hcgs_winter_26sites_300s  = filter(row -> contains_all(row.parameter_set_name, ["blocks","1440","26sites","300"]), results_df)
hcgs_winter_26sites_500s  = filter(row -> contains_all(row.parameter_set_name, ["blocks","1440","26sites","500"]), results_df)
hcgs_winter_26sites_1000s = filter(row -> contains_all(row.parameter_set_name, ["blocks","1440","26sites","1000"]), results_df)
# summer, 26 sites
hcgs_summer_26sites       = filter(row -> contains_all(row.parameter_set_name, ["blocks","793","26sites"]), results_df)
hcgs_summer_26sites_300s  = filter(row -> contains_all(row.parameter_set_name, ["blocks","793","26sites","300"]), results_df)
hcgs_summer_26sites_500s  = filter(row -> contains_all(row.parameter_set_name, ["blocks","793","26sites","500"]), results_df)
hcgs_summer_26sites_1000s = filter(row -> contains_all(row.parameter_set_name, ["blocks","793","26sites","1000"]), results_df)

#hcgs_td_err = [0.1, 0.15, 0.2]   # X error bars
#hcgs_var_err = [0.2, 0.25, 0.3]   # Y error bars

# Population size
popsize_winter_8sites = filter(row -> contains_all(row.parameter_set_name, ["popsize","1440","8sites"]), results_df)
popsize_winter_8sites_300s = filter(row -> contains_all(row.parameter_set_name, ["popsize","1440","8sites","300"]), results_df)
popsize_winter_8sites_500s = filter(row -> contains_all(row.parameter_set_name, ["popsize","1440","8sites","500"]), results_df)
popsize_winter_8sites_1000s = filter(row -> contains_all(row.parameter_set_name, ["popsize","1440","8sites","1000"]), results_df)

popsize_summer_8sites = filter(row -> contains_all(row.parameter_set_name, ["popsize","793","8sites"]), results_df)
popsize_summer_8sites_300s = filter(row -> contains_all(row.parameter_set_name, ["popsize","793","8sites","300"]), results_df)
popsize_summer_8sites_500s = filter(row -> contains_all(row.parameter_set_name, ["popsize","793","8sites","500"]), results_df)
popsize_summer_8sites_1000s = filter(row -> contains_all(row.parameter_set_name, ["popsize","793","8sites","1000"]), results_df)

popsize_winter_26sites = filter(row -> contains_all(row.parameter_set_name, ["popsize","1440","26sites"]), results_df)
popsize_winter_26sites_300s = filter(row -> contains_all(row.parameter_set_name, ["popsize","1440","26sites","300"]), results_df)
popsize_winter_26sites_500s = filter(row -> contains_all(row.parameter_set_name, ["popsize","1440","26sites","500"]), results_df)
popsize_winter_26sites_1000s = filter(row -> contains_all(row.parameter_set_name, ["popsize","1440","26sites","1000"]), results_df)

popsize_summer_26sites = filter(row -> contains_all(row.parameter_set_name, ["popsize","793","26sites"]), results_df)
popsize_summer_26sites_300s = filter(row -> contains_all(row.parameter_set_name, ["popsize","793","26sites","300"]), results_df)
popsize_summer_26sites_500s = filter(row -> contains_all(row.parameter_set_name, ["popsize","793","26sites","500"]), results_df)
popsize_summer_26sites_1000s = filter(row -> contains_all(row.parameter_set_name, ["popsize","793","26sites","1000"]), results_df)

# mSCAPE
mscape_winter_8sites = filter(row -> contains_all(row.parameter_set_name, ["mscape","1440","8sites"]), results_df)
mscape_winter_8sites_300s = filter(row -> contains_all(row.parameter_set_name, ["mscape","1440","8sites","300"]), results_df)
mscape_winter_8sites_500s = filter(row -> contains_all(row.parameter_set_name, ["mscape","1440","8sites","500"]), results_df)
mscape_winter_8sites_1000s = filter(row -> contains_all(row.parameter_set_name, ["mscape","1440","8sites","1000"]), results_df)

mscape_summer_8sites = filter(row -> contains_all(row.parameter_set_name, ["mscape","793","8sites"]), results_df)
mscape_summer_8sites_300s = filter(row -> contains_all(row.parameter_set_name, ["mscape","793","8sites","300"]), results_df)
mscape_summer_8sites_500s = filter(row -> contains_all(row.parameter_set_name, ["mscape","793","8sites","500"]), results_df)
mscape_summer_8sites_1000s = filter(row -> contains_all(row.parameter_set_name, ["mscape","793","8sites","1000"]), results_df)

mscape_winter_26sites = filter(row -> contains_all(row.parameter_set_name, ["mscape","1440","26sites"]), results_df)
mscape_winter_26sites_300s = filter(row -> contains_all(row.parameter_set_name, ["mscape","1440","26sites","300"]), results_df)
mscape_winter_26sites_500s = filter(row -> contains_all(row.parameter_set_name, ["mscape","1440","26sites","500"]), results_df)
mscape_winter_26sites_1000s = filter(row -> contains_all(row.parameter_set_name, ["mscape","1440","26sites","1000"]), results_df)

mscape_summer_26sites = filter(row -> contains_all(row.parameter_set_name, ["mscape","793","26sites"]), results_df)
mscape_summer_26sites_300s = filter(row -> contains_all(row.parameter_set_name, ["mscape","793","26sites","300"]), results_df)
mscape_summer_26sites_500s = filter(row -> contains_all(row.parameter_set_name, ["mscape","793","26sites","500"]), results_df)
mscape_summer_26sites_1000s = filter(row -> contains_all(row.parameter_set_name, ["mscape","793","26sites","1000"]), results_df)

# London only
london_winter_8sites = filter(row -> contains_all(row.parameter_set_name, ["london","1440","8sites"]), results_df)
london_winter_8sites_300s = filter(row -> contains_all(row.parameter_set_name, ["london","1440","8sites","300"]), results_df)
london_winter_8sites_500s = filter(row -> contains_all(row.parameter_set_name, ["london","1440","8sites","500"]), results_df)
london_winter_8sites_1000s = filter(row -> contains_all(row.parameter_set_name, ["london","1440","8sites","1000"]), results_df)

london_summer_8sites = filter(row -> contains_all(row.parameter_set_name, ["london","793","8sites"]), results_df)
london_summer_8sites_300s = filter(row -> contains_all(row.parameter_set_name, ["london","793","8sites","300"]), results_df)
london_summer_8sites_500s = filter(row -> contains_all(row.parameter_set_name, ["london","793","8sites","500"]), results_df)
london_summer_8sites_1000s = filter(row -> contains_all(row.parameter_set_name, ["london","793","8sites","1000"]), results_df)

london_winter_19sites = filter(row -> contains_all(row.parameter_set_name, ["london","1440","19sites"]), results_df)
london_winter_19sites_300s = filter(row -> contains_all(row.parameter_set_name, ["london","1440","19sites","300"]), results_df)
london_winter_19sites_500s = filter(row -> contains_all(row.parameter_set_name, ["london","1440","19sites","500"]), results_df)
london_winter_19sites_1000s = filter(row -> contains_all(row.parameter_set_name, ["london","1440","19sites","1000"]), results_df)

london_summer_19sites = filter(row -> contains_all(row.parameter_set_name, ["london","793","19sites"]), results_df)
london_summer_19sites_300s = filter(row -> contains_all(row.parameter_set_name, ["london","793","19sites","300"]), results_df)
london_summer_19sites_500s = filter(row -> contains_all(row.parameter_set_name, ["london","793","19sites","500"]), results_df)
london_summer_19sites_1000s = filter(row -> contains_all(row.parameter_set_name, ["london","793","19sites","1000"]), results_df)

## Create the scatter plot with error bars

#https://juliagraphics.github.io/Colors.jl/stable/namedcolors/
hcgs_color_td    = :mediumpurple3
popsize_color_td = :royalblue1
london_color_td  = :gold       
mscape_color_td  = :orangered

# Version using:
# - marker size = number of samples
# - colour = strategy
# - shape = number of Sites
# - line width = season

# Marker sizes
shape_8_sites     = :circle
shape_19_sites    = :utriangle
shape_26_sites    = :diamond
size_300_samples  = Int(300/100)+2
size_500_samples  = Int(500/100)+2
size_1000_samples = Int(1000/100)+2
markerstrokewidth_winter = 4 # outline thickness of markers
markerstrokewidth_summer = 1 # outline thickness of markers

# Range of TD values
minimum( results_df.TD_median) # 36.43458306288525
maximum( results_df.TD_median) # 48.57293908282003

# Range of Variance values
minimum( results_df.Var) # 67.46343565231585
maximum( results_df.Var) # 86.9112815037503

# Create the scatter plot with error bars
p = Plots.scatter( [1], [1]
                ,label = "", color = :white, marker = :circle, markersize = size_500_samples, markerstrokewidth = markerstrokewidth_winter, markerstrokecolor = :black
                ,legend = :bottomright, legend_frame = false, foreground_color_legend=:transparent, background_color_legend=:transparent, legend_columns = 3
                , alpha = 0.6
                , xaxis = [38,50], yaxis=[65,90]
)

Plots.scatter!([1],[1],label=" ",color = :white, marker = :circle, markersize = 0)
Plots.scatter!([1],[1],label=" ",color = :white, marker = :circle, markersize = 0)
Plots.scatter!([1],[1], color = hcgs_color_td, marker = shape_8_sites, markersize = size_300_samples, label="HCGS", alpha=0.6 )

Plots.scatter!([1],[1],label=" ",color = :white, marker = :circle, markersize = 0)
Plots.scatter!([1],[1],label="8 sites",color = :white, marker = shape_8_sites, markersize = size_500_samples)
Plots.scatter!([1],[1],label="Population size", color = popsize_color_td, marker = shape_8_sites, markersize = size_300_samples, alpha=0.6 )

Plots.scatter!([1],[1],label = "Winter", color = :white, marker = :circle, markersize = size_500_samples, markerstrokewidth = markerstrokewidth_winter, markerstrokecolor = :black)
Plots.scatter!([1],[1],label="19 sites",color = :white, marker = shape_19_sites, markersize = size_500_samples)
Plots.scatter!([1],[1],label="London only", color = london_color_td, marker = shape_8_sites, markersize = size_300_samples, alpha=0.6 )

Plots.scatter!([1],[1],label="Summer",color = :white, marker = :circle, markersize = size_500_samples)
Plots.scatter!([1],[1],label="26 sites",color = :white, marker = shape_26_sites, markersize = size_500_samples)
Plots.scatter!([1],[1],label="NHS-mSCAPE", color = mscape_color_td, marker = shape_8_sites, markersize = size_300_samples, alpha=0.6 )

#legend_marker_scale

#Plots.scatter!([1],[1],label=" ",color = :white, marker = :circle, markersize = 0)
#Plots.scatter!([1],[1],label=" ",color = :white, marker = :circle, markersize = 0)
#Plots.scatter!([1],[1],label=" ",color = :white, marker = :circle, markersize = 0)

#Plots.scatter!([1],[1],label="300 samples per week",color = :white, marker = :circle, markersize = size_300_samples)
#Plots.scatter!([1],[1],label="500 samples per week",color = :white, marker = :circle, markersize = size_500_samples)
#Plots.scatter!([1],[1],label="1000 samples per week",color = :white, marker = :circle, markersize = size_1000_samples)
#Plots.scatter!([1],[1],label=" ",color = :white, marker = :circle, markersize = 0)

# Hierarchical clustering and geographic stratification (HCGS) strategy
Plots.scatter!( hcgs_winter_8sites_300s[:,:TD_median], hcgs_winter_8sites_300s[:,:Var], color = hcgs_color_td, marker = shape_8_sites, markersize = size_300_samples, label="", alpha=0.6, markerstrokewidth = markerstrokewidth_winter, markerstrokecolor = :black )
Plots.scatter!( hcgs_winter_8sites_500s[:,:TD_median], hcgs_winter_8sites_500s[:,:Var], color = hcgs_color_td, marker = shape_8_sites, markersize = size_500_samples, label="", alpha=0.6, markerstrokewidth = markerstrokewidth_winter, markerstrokecolor = :black )
Plots.scatter!( hcgs_winter_8sites_500s[:,:TD_median], hcgs_winter_8sites_1000s[:,:Var], color = hcgs_color_td, marker = shape_8_sites, markersize = size_1000_samples, label="", alpha=0.6, markerstrokewidth = markerstrokewidth_winter, markerstrokecolor = :black )

Plots.scatter!( hcgs_winter_26sites_300s[:,:TD_median], hcgs_winter_26sites_300s[:,:Var], color = hcgs_color_td, marker = shape_26_sites, markersize = size_300_samples, label="", alpha=0.6, markerstrokewidth = markerstrokewidth_winter, markerstrokecolor = :black )
Plots.scatter!( hcgs_winter_26sites_500s[:,:TD_median], hcgs_winter_26sites_500s[:,:Var], color = hcgs_color_td, marker = shape_26_sites, markersize = size_500_samples, label="", alpha=0.6, markerstrokewidth = markerstrokewidth_winter, markerstrokecolor = :black )
Plots.scatter!( hcgs_winter_26sites_1000s[:,:TD_median], hcgs_winter_26sites_1000s[:,:Var], color = hcgs_color_td, marker = shape_26_sites, markersize = size_1000_samples, label="", alpha=0.6, markerstrokewidth = markerstrokewidth_winter, markerstrokecolor = :black )

Plots.scatter!( hcgs_summer_8sites_300s[:,:TD_median], hcgs_summer_8sites_300s[:,:Var], color = hcgs_color_td, marker = shape_8_sites, markersize = size_300_samples, label="", alpha=0.6 )
Plots.scatter!( hcgs_summer_8sites_500s[:,:TD_median], hcgs_summer_8sites_500s[:,:Var], color = hcgs_color_td, marker = shape_8_sites, markersize = size_500_samples, label="", alpha=0.6 )
Plots.scatter!( hcgs_summer_8sites_1000s[:,:TD_median], hcgs_summer_8sites_1000s[:,:Var], color = hcgs_color_td, marker = shape_8_sites, markersize = size_1000_samples, label="", alpha=0.6 )

Plots.scatter!( hcgs_summer_26sites_300s[:,:TD_median], hcgs_summer_26sites_300s[:,:Var], color = hcgs_color_td, marker = shape_26_sites, markersize = size_300_samples, label="", alpha=0.6 )
Plots.scatter!( hcgs_summer_26sites_500s[:,:TD_median], hcgs_summer_26sites_500s[:,:Var], color = hcgs_color_td, marker = shape_26_sites, markersize = size_500_samples, label="", alpha=0.6 )
Plots.scatter!( hcgs_summer_26sites_1000s[:,:TD_median], hcgs_summer_26sites_1000s[:,:Var], color = hcgs_color_td, marker = shape_26_sites, markersize = size_1000_samples, label="", alpha=0.6 )

# Population size site selection strategy
Plots.scatter!( popsize_winter_8sites_300s[:, :TD_median], popsize_winter_8sites_300s[:,:Var], color = popsize_color_td, marker = shape_8_sites, markersize = size_300_samples, label="", alpha=0.6, markerstrokewidth = markerstrokewidth_winter, markerstrokecolor = :black )
Plots.scatter!( popsize_winter_8sites_500s[:, :TD_median], popsize_winter_8sites_500s[:,:Var], color = popsize_color_td, marker = shape_8_sites, markersize = size_500_samples, label="", alpha=0.6, markerstrokewidth = markerstrokewidth_winter, markerstrokecolor = :black )
Plots.scatter!( popsize_winter_8sites_1000s[:, :TD_median], popsize_winter_8sites_1000s[:,:Var], color = popsize_color_td, marker = shape_8_sites, markersize = size_1000_samples, label="", alpha=0.6, markerstrokewidth = markerstrokewidth_winter, markerstrokecolor = :black )

Plots.scatter!( popsize_winter_26sites_300s[:,:TD_median], popsize_winter_26sites_300s[:,:Var], color = popsize_color_td, marker = shape_26_sites, markersize = size_300_samples, label="", alpha=0.6, markerstrokewidth = markerstrokewidth_winter, markerstrokecolor = :black )
Plots.scatter!( popsize_winter_26sites_500s[:,:TD_median], popsize_winter_26sites_500s[:,:Var], color = popsize_color_td, marker = shape_26_sites, markersize = size_500_samples, label="", alpha=0.6, markerstrokewidth = markerstrokewidth_winter, markerstrokecolor = :black )
Plots.scatter!( popsize_winter_26sites_1000s[:,:TD_median], popsize_winter_26sites_1000s[:,:Var], color = popsize_color_td, marker = shape_26_sites, markersize = size_1000_samples, label="", alpha=0.6, markerstrokewidth = markerstrokewidth_winter, markerstrokecolor = :black )

Plots.scatter!( popsize_summer_8sites_300s[:, :TD_median], popsize_summer_8sites_300s[:,:Var], color = popsize_color_td, marker = shape_8_sites, markersize = size_300_samples, label="", alpha=0.6 )
Plots.scatter!( popsize_summer_8sites_500s[:, :TD_median], popsize_summer_8sites_500s[:,:Var], color = popsize_color_td, marker = shape_8_sites, markersize = size_500_samples, label="", alpha=0.6 )
Plots.scatter!( popsize_summer_8sites_1000s[:, :TD_median], popsize_summer_8sites_1000s[:,:Var], color = popsize_color_td, marker = shape_8_sites, markersize = size_1000_samples, label="", alpha=0.6 )

Plots.scatter!( popsize_summer_26sites_300s[:,:TD_median], popsize_summer_26sites_300s[:,:Var], color = popsize_color_td, marker = shape_26_sites, markersize = size_300_samples, label="", alpha=0.6 )
Plots.scatter!( popsize_summer_26sites_500s[:,:TD_median], popsize_summer_26sites_500s[:,:Var], color = popsize_color_td, marker = shape_26_sites, markersize = size_500_samples, label="", alpha=0.6 )
Plots.scatter!( popsize_summer_26sites_1000s[:,:TD_median], popsize_summer_26sites_1000s[:,:Var], color = popsize_color_td, marker = shape_26_sites, markersize = size_1000_samples, label="", alpha=0.6 )

# London only site selection strategy
Plots.scatter!( london_winter_8sites_300s[:, :TD_median], london_winter_8sites_300s[:,:Var], color = london_color_td, marker = shape_8_sites, markersize = size_300_samples, label="", alpha=0.6, markerstrokewidth = markerstrokewidth_winter, markerstrokecolor = :black ) 
Plots.scatter!( london_winter_8sites_500s[:, :TD_median], london_winter_8sites_500s[:,:Var], color = london_color_td, marker = shape_8_sites, markersize = size_500_samples, label="", alpha=0.6, markerstrokewidth = markerstrokewidth_winter, markerstrokecolor = :black ) 
Plots.scatter!( london_winter_8sites_1000s[:, :TD_median], london_winter_8sites_1000s[:,:Var], color = london_color_td, marker = shape_8_sites, markersize = size_1000_samples, label="", alpha=0.6, markerstrokewidth = markerstrokewidth_winter, markerstrokecolor = :black ) 

Plots.scatter!( london_winter_19sites_300s[:,:TD_median], london_winter_19sites_300s[:,:Var], color = london_color_td, marker = shape_19_sites, markersize = size_300_samples, label="", alpha=0.6, markerstrokewidth = markerstrokewidth_winter, markerstrokecolor = :black )
Plots.scatter!( london_winter_19sites_500s[:,:TD_median], london_winter_19sites_500s[:,:Var], color = london_color_td, marker = shape_19_sites, markersize = size_500_samples, label="", alpha=0.6, markerstrokewidth = markerstrokewidth_winter, markerstrokecolor = :black )
Plots.scatter!( london_winter_19sites_1000s[:,:TD_median], london_winter_19sites_1000s[:,:Var], color = london_color_td, marker = shape_19_sites, markersize = size_1000_samples, label="", alpha=0.6, markerstrokewidth = markerstrokewidth_winter, markerstrokecolor = :black )

Plots.scatter!( london_summer_8sites_300s[:, :TD_median], london_summer_8sites_300s[:,:Var], color = london_color_td, marker = shape_8_sites, markersize = size_300_samples, label="", alpha=0.6 )
Plots.scatter!( london_summer_8sites_500s[:, :TD_median], london_summer_8sites_500s[:,:Var], color = london_color_td, marker = shape_8_sites, markersize = size_500_samples, label="", alpha=0.6 )
Plots.scatter!( london_summer_8sites_1000s[:, :TD_median], london_summer_8sites_1000s[:,:Var], color = london_color_td, marker = shape_8_sites, markersize = size_1000_samples, label="", alpha=0.6 )

Plots.scatter!( london_summer_19sites_300s[:,:TD_median], london_summer_19sites_300s[:,:Var], color = london_color_td, marker = shape_19_sites, markersize = size_300_samples, label="", alpha=0.6 )
Plots.scatter!( london_summer_19sites_500s[:,:TD_median], london_summer_19sites_500s[:,:Var], color = london_color_td, marker = shape_19_sites, markersize = size_500_samples, label="", alpha=0.6 )
Plots.scatter!( london_summer_19sites_1000s[:,:TD_median], london_summer_19sites_1000s[:,:Var], color = london_color_td, marker = shape_19_sites, markersize = size_1000_samples, label="", alpha=0.6 )

# mSCAPE sites
Plots.scatter!( mscape_winter_8sites_300s[:, :TD_median], mscape_winter_8sites_300s[:,:Var], color = mscape_color_td, marker = shape_8_sites, markersize = size_300_samples, label="", alpha=0.6, markerstrokewidth = markerstrokewidth_winter, markerstrokecolor = :black ) 
Plots.scatter!( mscape_winter_8sites_500s[:, :TD_median], mscape_winter_8sites_500s[:,:Var], color = mscape_color_td, marker = shape_8_sites, markersize = size_500_samples, label="", alpha=0.6, markerstrokewidth = markerstrokewidth_winter, markerstrokecolor = :black ) 
Plots.scatter!( mscape_winter_8sites_1000s[:, :TD_median], mscape_winter_8sites_1000s[:,:Var], color = mscape_color_td, marker = shape_8_sites, markersize = size_1000_samples, label="", alpha=0.6, markerstrokewidth = markerstrokewidth_winter, markerstrokecolor = :black ) 

Plots.scatter!( mscape_winter_26sites_300s[:,:TD_median], mscape_winter_26sites_300s[:,:Var], color = mscape_color_td, marker = shape_26_sites, markersize = size_300_samples, label="", alpha=0.6, markerstrokewidth = markerstrokewidth_winter, markerstrokecolor = :black ) 
Plots.scatter!( mscape_winter_26sites_500s[:,:TD_median], mscape_winter_26sites_500s[:,:Var], color = mscape_color_td, marker = shape_26_sites, markersize = size_500_samples, label="", alpha=0.6, markerstrokewidth = markerstrokewidth_winter, markerstrokecolor = :black ) 
Plots.scatter!( mscape_winter_26sites_1000s[:,:TD_median], mscape_winter_26sites_1000s[:,:Var], color = mscape_color_td, marker = shape_26_sites, markersize = size_1000_samples, label="", alpha=0.6, markerstrokewidth = markerstrokewidth_winter, markerstrokecolor = :black ) 

Plots.scatter!( mscape_summer_8sites_300s[:, :TD_median], mscape_summer_8sites_300s[:,:Var], color = mscape_color_td, marker = shape_8_sites, markersize = size_300_samples, label="", alpha=0.6 )
Plots.scatter!( mscape_summer_8sites_500s[:, :TD_median], mscape_summer_8sites_500s[:,:Var], color = mscape_color_td, marker = shape_8_sites, markersize = size_500_samples, label="", alpha=0.6 )
Plots.scatter!( mscape_summer_8sites_1000s[:, :TD_median], mscape_summer_8sites_1000s[:,:Var], color = mscape_color_td, marker = shape_8_sites, markersize = size_1000_samples, label="", alpha=0.6 )

Plots.scatter!( mscape_summer_26sites_300s[:,:TD_median], mscape_summer_26sites_300s[:,:Var], color = mscape_color_td, marker = shape_26_sites, markersize = size_300_samples, label="", alpha=0.6)#,markerstrokecolor=:black, markerstrokewidth=6 )
Plots.scatter!( mscape_summer_26sites_500s[:,:TD_median], mscape_summer_26sites_500s[:,:Var], color = mscape_color_td, marker = shape_26_sites, markersize = size_500_samples, label="", alpha=0.6)#,markerstrokecolor=:black, markerstrokewidth=6 )
Plots.scatter!( mscape_summer_26sites_1000s[:,:TD_median], mscape_summer_26sites_1000s[:,:Var], color = mscape_color_td, marker = shape_26_sites, markersize = size_1000_samples, label="", alpha=0.6)#,markerstrokecolor=:black, markerstrokewidth=6 )

# Add titles and labels
xlabel!(p, "Median time to detection from first importation (days)")
ylabel!(p, "Variance of time to detection")# (days^2)")
#ylabel!(p, L"\mathrm{Variance of time to detection (days^2)}")
#title!(p, "Title")

# Display the plot
display(p)

Plots.savefig( joinpath( root_dir
                        , "NBPMscape/paper/analysis/2_sampling_analysis/covid_like/geo_optimisation"
                        , "icu_tds_var_8_and_26_sites.png" ) )
Plots.savefig( joinpath( root_dir
                        , "NBPMscape/paper/results/figures_main/raw"
                        , "icu_tds_var_8_and_26_sites.png" ) )


####################################################################################################
# (7) Plot TD variance vs TD 
#   (7.2) Function for plotting
####################################################################################################

function plot_var_v_td(; n_blocks = n_blocks
                        ,n_blocks_shape = n_blocks_shape
                        ,n_sites = n_sites
                        ,n_sites_color = n_sites_color
                        ,n_samples = n_samples
                        ,n_samples_size = n_samples_size
                        ,seasonal_ari = seasonal_ari
                        ,seasonal_ari_lwidth = seasonal_ari_lwidth
                        ,site_allocation_strategies = site_allocation_strategies
                        ,xaxis_minmax = [0,100]
                        ,yaxis_minmax = [0,100]
                        )

    # Initialise plot
    p = Plots.scatter( [1], [1],label = ""
                    #, color = :white, marker = :circle, markersize = size_500_samples, markerstrokewidth = markerstrokewidth_winter, markerstrokecolor = :black
                    , color = :white, marker = :circle, markersize = n_samples_size[2], markerstrokewidth = 4, markerstrokecolor = :black
                    ,legend = :bottomright, legend_frame = false, foreground_color_legend=:transparent, background_color_legend=:transparent
                    , legend_columns = 2
                    , alpha = 0.6
                    , xaxis = xaxis_minmax, yaxis=yaxis_minmax
    )

    # Add legend entries
    # Blocks
    for i in 1:length(n_blocks)
        Plots.scatter!([1],[1],label = n_blocks[i], color = :white, marker = n_blocks_shape[i], markersize = 5, markerstrokewidth = 0.5)
    end
    # Blank space
    Plots.scatter!([1],[1],label=" ",color = :white, marker = :circle, markersize = 0)
    Plots.scatter!([1],[1],label=" ",color = :white, marker = :circle, markersize = 0)
    # Sites
    for i in 1:length(n_sites)
        Plots.scatter!([1],[1],label = n_sites[i], color = n_sites_color[i], marker = :circle, markersize = 5, markerstrokewidth = 1)
    end
    
    display(p)
    
    # Loop through each parameter set
    for i in 1:length(n_blocks) #i=1
        for j in 1:length(n_sites) #j=1
            for k in 1:length(n_samples) #k=1
                for m in 1:1#length(seasonal_ari) #m=2
                    for n in 1:length(site_allocation_strategies) #n=1
                        
                        # Filter results that match the parameter_set_name
                        #filtered_results = filter(row -> contains_all(row.parameter_set_name, ["blocks","793","30sites"]), results_df)
                        filtered_results = filter(row -> contains_all(row.parameter_set_name, [n_blocks[i],n_sites[j],n_samples[k],seasonal_ari[m],site_allocation_strategies[n]])
                                                , results_df) #println(results_df)
                        
                        # Determing markers for parameter combination
                        filtered_shape = n_blocks_shape[i]
                        filtered_color = n_sites_color[j]
                        filtered_size = n_samples_size[k]
                        filtered_lwidth = 1#seasonal_ari_lwidth[m]
                                                                
                        # Plot filtered results
                        Plots.scatter!( filtered_results[:,:TD_median], filtered_results[:,:Var]
                                        , color = filtered_color
                                        , marker = filtered_shape
                                        , markersize = filtered_size
                                        , markerstrokewidth = filtered_lwidth
                                        , markerstrokecolor = :black 
                                        , label = ""
                                        , alpha = 0.6 
                                        )                    

                    end
                end
            end
        end
    end
# Add titles and labels
xlabel!(p, "Median time to detection from first importation (days)")
ylabel!(p, "Variance of time to detection")# (days^2)")
#ylabel!(p, L"\mathrm{Variance of time to detection (days^2)}")
#title!(p, "Title")
display(p)
return(p)
end


####################################################################################################
# (7) Plot TD variance vs TD 
#   (7.3) Plot HCGS sampling scenarios only, but split by samples per week and season
####################################################################################################

# Define parameter values and associated marker options
# - marker size = number of samples
# - colour = number of sites
# - shape = number of blocks for HCGS
# - line width = season

## Split by sample size and season
n_blocks = ["5blocks","10blocks","20blocks","30blocks"]
n_blocks_shape = [:utriangle,:circle,:square,:diamond]
#n_sites = ["8sites","10sites","19sites","20sites","26sites","30sites","40sites","50sites","60sites"]
#n_sites_color = [:yellow,:red,:orange,:blue,:green,:purple,:,:,:]
n_sites = ["8sites","10sites","20sites","26sites","30sites","40sites","50sites","60sites"]
#n_sites_color = [:yellow,:red,:blue,:green,:purple,:,:,:]
n_sites_color = [:grey100,:grey90,:grey80,:grey70,:grey60,:grey50,:grey40,:grey30]
n_sites_color = [:red,:orange,:yellow,:green,:blue,:indigo,:violet,:black]
#n_samples = ["300","500","1000"]
#n_samples_size = [Int(300/100)+2,Int(500/100)+2,Int(1000/100)+2]
#n_samples = ["1000"]
#n_samples_size = [Int(500/100)+2,Int(500/100)+2]
seasonal_ari = ["1440","793"]
seasonal_ari_lwidth = [4,1]
site_allocation_strategies = [""] #["hcgs"] HCGS not actually listed in parameter name because can be distinguished from inclusion of blocks in the name # Only 1 so no need to change any marker options in this plot

# Range of TD values
minimum( results_df.TD_median) # 36.43458306288525
maximum( results_df.TD_median) # 48.57293908282003

# Range of Variance values
minimum( results_df.Var) # 67.46343565231585
maximum( results_df.Var) # 86.9112815037503

p_winter_300 = plot_var_v_td(; n_blocks = n_blocks
                        ,n_blocks_shape = n_blocks_shape
                        ,n_sites = n_sites
                        ,n_sites_color = n_sites_color
                        ,n_samples = ["300"]#n_samples
                        ,n_samples_size = [Int(300/100)+2,Int(300/100)+2]#n_samples_size
                        ,seasonal_ari = ["1440"] #seasonal_ari
                        ,seasonal_ari_lwidth = 1#seasonal_ari_lwidth
                        ,site_allocation_strategies = site_allocation_strategies
                        ,xaxis_minmax = [35,50]
                        ,yaxis_minmax = [65,90]
                        )
p_winter_500 = plot_var_v_td(; n_blocks = n_blocks
                        ,n_blocks_shape = n_blocks_shape
                        ,n_sites = n_sites
                        ,n_sites_color = n_sites_color
                        ,n_samples = ["500"]#n_samples
                        ,n_samples_size = [Int(500/100)+2,Int(500/100)+2]#n_samples_size
                        ,seasonal_ari = ["1440"] #seasonal_ari
                        ,seasonal_ari_lwidth = 1#seasonal_ari_lwidth
                        ,site_allocation_strategies = site_allocation_strategies
                        ,xaxis_minmax = [35,50]
                        ,yaxis_minmax = [65,90]
                        )
p_winter_1000 = plot_var_v_td(; n_blocks = n_blocks
                        ,n_blocks_shape = n_blocks_shape
                        ,n_sites = n_sites
                        ,n_sites_color = n_sites_color
                        ,n_samples = ["1000"]#n_samples
                        ,n_samples_size = [Int(1000/100)+2,Int(1000/100)+2]#n_samples_size
                        ,seasonal_ari = ["1440"] #seasonal_ari
                        ,seasonal_ari_lwidth = 1#seasonal_ari_lwidth
                        ,site_allocation_strategies = site_allocation_strategies
                        ,xaxis_minmax = [35,50]
                        ,yaxis_minmax = [65,90]
                        )

p_summer_300 = plot_var_v_td(; n_blocks = n_blocks
                        ,n_blocks_shape = n_blocks_shape
                        ,n_sites = n_sites
                        ,n_sites_color = n_sites_color
                        ,n_samples = ["300"]#n_samples
                        ,n_samples_size = [Int(300/100)+2,Int(300/100)+2]#n_samples_size
                        ,seasonal_ari = ["793"] #seasonal_ari
                        ,seasonal_ari_lwidth = 1#seasonal_ari_lwidth
                        ,site_allocation_strategies = site_allocation_strategies
                        ,xaxis_minmax = [35,50]
                        ,yaxis_minmax = [65,90]
                        )
p_summer_500 = plot_var_v_td(; n_blocks = n_blocks
                        ,n_blocks_shape = n_blocks_shape
                        ,n_sites = n_sites
                        ,n_sites_color = n_sites_color
                        ,n_samples = ["500"]#n_samples
                        ,n_samples_size = [Int(500/100)+2,Int(500/100)+2]#n_samples_size
                        ,seasonal_ari = ["793"] #seasonal_ari
                        ,seasonal_ari_lwidth = 1#seasonal_ari_lwidth
                        ,site_allocation_strategies = site_allocation_strategies
                        ,xaxis_minmax = [35,50]
                        ,yaxis_minmax = [65,90]
                        )
p_summer_1000 = plot_var_v_td(; n_blocks = n_blocks
                        ,n_blocks_shape = n_blocks_shape
                        ,n_sites = n_sites
                        ,n_sites_color = n_sites_color
                        ,n_samples = ["1000"]#n_samples
                        ,n_samples_size = [Int(1000/100)+2,Int(1000/100)+2]#n_samples_size
                        ,seasonal_ari = ["793"] #seasonal_ari
                        ,seasonal_ari_lwidth = 1#seasonal_ari_lwidth
                        ,site_allocation_strategies = site_allocation_strategies
                        ,xaxis_minmax = [35,50]
                        ,yaxis_minmax = [65,90]
                        )


Plots.plot(p_winter_300,p_summer_300,p_winter_500,p_summer_500,p_winter_1000,p_summer_1000, layout = (3, 2)
        , size = (1000,1000)) #default size is 600,400)

Plots.savefig( joinpath( root_dir
                , "NBPMscape/paper/analysis/2_sampling_analysis/covid_like/geo_optimisation"
                , "icu_tds_var_HCGS_only_split_by_col_season_row_sample_size.png" ) )
Plots.savefig( joinpath( root_dir
                , "NBPMscape/paper/results/figures_main/raw"
                , "icu_tds_var_HCGS_only_split_by_col_season_row_sample_size.png" ) )


#####################################################################################################################                
# (8) Look at cumulative GP attendances at TD re syndromic surveillance
#####################################################################################################################

# Reload ICU TD data if necessary
icu_tds_df_dict = load( joinpath( root_dir, "NBPMscape/paper/analysis/2_sampling_analysis/covid_like/geo_optimisation"
                                    , "covidlike-1.6.1_nrep1000_4300448_icu_tds_df_dict.jld2")
                        , "icu_tds_df_dict")

Plots.histogram( icu_tds_df_dict["icu_tds_seasonal_ari_1440_500_sampling_sites_popsize_26sites"].cum_tgp_td )
median(icu_tds_df_dict["icu_tds_seasonal_ari_1440_500_sampling_sites_popsize_26sites"].cum_tgp_td)  #  150.5
median(icu_tds_df_dict["icu_tds_seasonal_ari_1440_500_sampling_sites_popsize_26sites"].cum_tinf_td) # 4024.5

# Collate for all sampling strategies investigated:
# - median TD
# - cum_tgp_td at median TD
# - median cum_tgp_td 
# - TD at median cum_tgp_td
# - 2.5% and 97.5% quantiles for cum_tgp_td 
# - 2.5% and 97.5% quantiles for TD
# - cum_tinf_td at median TD
# - median cum_tinf_td and TD at median cum_tgp_td
# - 2.5% and 97.5% quantiles for cum_tgp_td and TD

cum_tgp_df = DataFrame( sampling_scenario = []
                        , median_td = []
                        , median_cum_tgp_td = []
                        , cum_tgp_td_quantile_025 = []
                        , cum_tgp_td_quantile_975 = []
                        , median_cum_tinf_td = []
                        , cum_tinf_td_quantile_025 = []
                        , cum_tinf_td_quantile_975 = []
                        , median_cum_ted_td = []
                        , cum_ted_td_quantile_025 = []
                        , cum_ted_td_quantile_975 = []
                        , median_cum_thosp_td = []
                        , cum_thosp_td_quantile_025 = []
                        , cum_thosp_td_quantile_975 = []
                        , median_cum_ticu_td = []
                        , cum_ticu_td_quantile_025 = []
                        , cum_ticu_td_quantile_975 = []
)

for k in keys(icu_tds_df_dict) #  k = "icu_tds_seasonal_ari_1440_500_sampling_sites_popsize_26sites"
    
    median_td                 = median(icu_tds_df_dict[k].ICU_TD)
    median_cum_tgp_td         = median(icu_tds_df_dict[k].cum_tgp_td)
    cum_tgp_td_quantile_025   = quantile(icu_tds_df_dict[k].cum_tgp_td, 0.025)
    cum_tgp_td_quantile_975   = quantile(icu_tds_df_dict[k].cum_tgp_td, 0.975)
    median_cum_tinf_td        = median(icu_tds_df_dict[k].cum_tinf_td)
    cum_tinf_td_quantile_025  = quantile(icu_tds_df_dict[k].cum_tinf_td, 0.025)
    cum_tinf_td_quantile_975  = quantile(icu_tds_df_dict[k].cum_tinf_td, 0.975)
    median_cum_ted_td         = median(icu_tds_df_dict[k].cum_ted_td)
    cum_ted_td_quantile_025   = quantile(icu_tds_df_dict[k].cum_ted_td, 0.025)
    cum_ted_td_quantile_975   = quantile(icu_tds_df_dict[k].cum_ted_td, 0.975)
    median_cum_thosp_td       = median(icu_tds_df_dict[k].cum_thosp_td)
    cum_thosp_td_quantile_025 = quantile(icu_tds_df_dict[k].cum_thosp_td, 0.025)
    cum_thosp_td_quantile_975 = quantile(icu_tds_df_dict[k].cum_thosp_td, 0.975)
    median_cum_ticu_td        = median(icu_tds_df_dict[k].cum_ticu_td)
    cum_ticu_td_quantile_025  = quantile(icu_tds_df_dict[k].cum_ticu_td, 0.025)
    cum_ticu_td_quantile_975  = quantile(icu_tds_df_dict[k].cum_ticu_td, 0.975)

    new_row = [ k, median_td
                 , median_cum_tgp_td, cum_tgp_td_quantile_025, cum_tgp_td_quantile_975
                 , median_cum_tinf_td, cum_tinf_td_quantile_025, cum_tinf_td_quantile_975
                 , median_cum_ted_td, cum_ted_td_quantile_025, cum_ted_td_quantile_975
                 , median_cum_thosp_td, cum_thosp_td_quantile_025, cum_thosp_td_quantile_975
                 , median_cum_ticu_td, cum_ticu_td_quantile_025, cum_ticu_td_quantile_975
                 ]
    
    push!(cum_tgp_df, new_row)
end
println( cum_tgp_df )

Plots.scatter(  cum_tgp_df.median_cum_tinf_td, cum_tgp_df.median_cum_tgp_td   
              , xlabel = "median cumulative infections at TD"
              , ylabel = "median cumulative healthcare visits at TD"
              , title = "Number of healthcare visits against infections\nat ICU time of detection\nfor a COVID-like pathogen"
              , legend = :topleft, label = "GP visits")
Plots.scatter!( cum_tgp_df.median_cum_tinf_td, cum_tgp_df.median_cum_ted_td,   label = "ED attendances"  )
Plots.scatter!( cum_tgp_df.median_cum_tinf_td, cum_tgp_df.median_cum_thosp_td, label = "Hospital admissions" )
Plots.scatter!( cum_tgp_df.median_cum_tinf_td, cum_tgp_df.median_cum_ticu_td,  label = "ICU admissions" )

mean( cum_tgp_df.median_cum_tgp_td ./ cum_tgp_df.median_cum_tinf_td ) # 0.03663608548275423

filter!(row -> !occursin( "icu_tds_seasonal_ari_1440_300_sampling_sites_30blocks_8sites"
                        , String( row["sampling_scenario"] ) ), cum_tgp_df)

# Plot
Plots.scatter( cum_tgp_df.median_cum_tgp_td, cum_tgp_df.median_cum_tinf_td 
            , xlabel = "median cumulative GP attendances at TD"
            , ylabel = "median cumulative infections at TD"
            , title = "GP attendances and infections\nat ICU time of detection\nfor a COVID-like pathogen")

CSV.write( joinpath( root_dir, "NBPMscape/paper/analysis/1_outbreak_simulations/covid_like/output/cum_tgp_df.csv"), cum_tgp_df )


#########################################################################################################################
# (9) Ranking sampling strategies using mean-risk objective optimisation
#     (Conditional value-at-risk (CVaR) and loss function using variance as risk measure)
#########################################################################################################################

# Reload ICU TD data if necessary
icu_tds_df_dict = load( joinpath( root_dir, "NBPMscape/paper/analysis/2_sampling_analysis/covid_like/geo_optimisation"
                                    , "covidlike-1.6.1_nrep1000_4300448_icu_tds_df_dict.jld2")
                        , "icu_tds_df_dict")

# Using functions in 'NBPMscape/src/sample_strategy_ranking.jl'

# Manipulate to format for strategy_rank function
icu_tds      = results_dict_to_df(df_dict = icu_tds_df_dict, col_to_extract = :ICU_TD )
icu_3tds     = results_dict_to_df(df_dict = icu_tds_df_dict, col_to_extract = :ICU_3TD )
cum_inf_tds  = results_dict_to_df(df_dict = icu_tds_df_dict, col_to_extract = :cum_tinf_td )
cum_inf_3tds = results_dict_to_df(df_dict = icu_tds_df_dict, col_to_extract = :cum_tinf_td3 )


## Example ranking options
# using loss function with variance as risk measure
strat_rank_results_icu_tds_loss = strategy_rank( df = icu_tds, rank_by = "loss", lambda=1.0, phi=1.0, alpha=0.95, var_quantile_method = "exact")
show(strat_rank_results_icu_tds_loss.results_df_sorted, truncate = 0)
display( strat_rank_results_icu_tds_loss.mean_risk_plot )

# using CVaR
strat_rank_results_icu_tds_cvar = strategy_rank( df = icu_tds, rank_by = "cvar", lambda=1.0, phi=0.5, alpha=0.95, var_quantile_method = "exact")
show(strat_rank_results_icu_tds_cvar.results_df_sorted, truncate = 0)
display( strat_rank_results_icu_tds_cvar.mean_risk_plot )

### Only strategy results for winter, 300 samples per week, 8 sites
# Filter
icu_tds_1440_300_8sites = filter_strat_names( df = icu_tds, substrings = ["1440", "300", "8sites"])

# Rank using CVaR
strat_rank_results_icu_tds_1440_300_8_cvar = strategy_rank( df = icu_tds_1440_300_8sites, rank_by = "cvar", lambda=1.0, phi=0.5, alpha=0.95, var_quantile_method = "exact")
show(strat_rank_results_icu_tds_1440_300_8_cvar.results_df_sorted, truncate = 0)
display( strat_rank_results_icu_tds_1440_300_8_cvar.mean_risk_plot )

# Rank using loss function with variance as risk measure
strat_rank_results_icu_tds_1440_300_8_loss = strategy_rank( df = icu_tds_1440_300_8sites, rank_by = "loss", lambda=0.5, phi=0.5, alpha=0.95, var_quantile_method = "exact")
show(strat_rank_results_icu_tds_1440_300_8_loss.results_df_sorted, allcols=true, truncate=0)
show(strat_rank_results_icu_tds_1440_300_8_loss.results_df_sorted, truncate=0)
display( strat_rank_results_icu_tds_1440_300_8_loss.mean_risk_plot )