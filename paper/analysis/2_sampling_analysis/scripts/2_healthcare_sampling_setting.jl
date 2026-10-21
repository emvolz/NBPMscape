#=
Script for the investigation of the effect on surveillance performance by sampling in different 
healthcare settings:
- primary (general practitioners modelled based on RSC-RCGP)
- tertiary (ICUs in hospitals modelled based on mSCAPE)

Surveillance performance determined by times 
to detection for 1 case (TD) and 3 cases (3TD) and the cumulative number of infections at 
time of detection.

This file contains:
(0) Define root directory
(1) Load sampling results for primary and tertiary care sampling
(2) Add outbreak metrics at TD values
(3) Replicate the datasets but with scenario number labelling rather than scenario description
(4) Plot median time to detection values for all sampling scenarios
    (4.1) Time to detect 1 case (TD)
    (4.2) Time to detect 3 cases (3TD)
(5) Analysis of multi-modal sampling - sampling from combinations of healthcare settings
    (5.1) Computing combined times to detection
    (5.2) Plotting combined times to detection
(6) Specific values for article text
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
# (1) Generate sampling results for primary care sampling and load results for acute tertiary care sampling
#    (1.1) Primary care sampling   - modelled on Oxford-RCGP RSC
#    (1.2) Tertiary care sampling  - modelled on mSCAPE
############################################################################################

# File names and parameter names

# Primary care (GPs)
# File numbers correspond to the configuration parameter file numbers in 'NBPMscape/config/2026_02_18'
gp_td_results_file_params = [ # RCGP 
                             (file_name = "gp_tds_1.csv",  params_name = "gp_tds_winter_300samples")
                            ,(file_name = "gp_tds_13.csv", params_name = "gp_tds_winter_500samples")
                            ,(file_name = "gp_tds_25.csv", params_name = "gp_tds_winter_1000samples")
                            ,(file_name = "gp_tds_37.csv", params_name = "gp_tds_summer_300samples")
                            ,(file_name = "gp_tds_49.csv", params_name = "gp_tds_summer_500samples")
                            ,(file_name = "gp_tds_61.csv", params_name = "gp_tds_summer_1000samples")
]

# Reload sims file if necessary
sims = load( joinpath( root_dir, "NBPMscape/paper/analysis/1_outbreak_simulations/covid_like/output/covidlike-1.6.1-sims_G_filtered_nrep1000_4300448.jld2" )
            , "sims")

### Compute times to detection with sampling in primary care (GPs)
# using different model parameters representing different scenarios
# See 'NBPMscape/config/2026_02_18/config_file_record.xlsx' for summary descriptions of scenarios (and indicated below by numbers in for loops)

# Start with base parameters
initialize_parameters("config/outbreak_params_covid19_like_20261008.yaml")

# Set sampling parameters for different scenarios
gp_sampling_params = (  scenario_01 = (scenario_number =  1, gp_ari_consults = 327, gp_ari_swabs =  747, gp_swabs_mg =  300)
                      , scenario_13 = (scenario_number = 13, gp_ari_consults = 327, gp_ari_swabs =  747, gp_swabs_mg =  500)
                      , scenario_25 = (scenario_number = 25, gp_ari_consults = 327, gp_ari_swabs = 1000, gp_swabs_mg = 1000)
                      , scenario_37 = (scenario_number = 37, gp_ari_consults = 180, gp_ari_swabs =  319, gp_swabs_mg =  300)
                      , scenario_49 = (scenario_number = 49, gp_ari_consults = 180, gp_ari_swabs =  500, gp_swabs_mg =  500)
                      , scenario_61 = (scenario_number = 61, gp_ari_consults = 180, gp_ari_swabs = 1000, gp_swabs_mg = 1000)
                      )

# Run analysis for 6 different scenarios (2 seasons x 3 numbers of samples) using GPs
for i in 1:length(gp_sampling_params) # in [1,13,25,37,49,61] # i=1
    gp_tds = gp_td(  sims = sims
                   , gp_ari_consults = gp_sampling_params[i].gp_ari_consults
                   , gp_ari_swabs    = gp_sampling_params[i].gp_ari_swabs
                   , gp_swabs_mg     = gp_sampling_params[i].gp_swabs_mg
                   )
    CSV.write( joinpath( root_dir, "NBPMscape/paper/analysis/2_sampling_analysis/covid_like/gp_td"
                        ,"gp_tds_$(gp_sampling_params[i].scenario_number).csv" )
                , gp_tds )
    # Also, save file with 'GP_simid' column removed to reduce file size
    gp_tds_small = gp_tds[:,["sim_n","GP_TD","GP_3TD","n_GP_cases","n_GP_cases_sampled","n_GP_cases_sampled_positive","tinf_relating_to_GP_TD","last_tinf_relating_to_GP_3TD"]]
    CSV.write( joinpath( root_dir, "NBPMscape/paper/analysis/2_sampling_analysis/covid_like/gp_td_small"
                        ,"small_gp_tds_$(gp_sampling_params[i].scenario_number).csv" )
                , gp_tds_small )
end

# Loop through scenarios, read in files and compute median TD values

# Prepare a container to keep the DataFrames
gp_tds_df_dict = Dict{String, DataFrame}()

for i in 1:length(gp_td_results_file_params) #i=1
    params_name = gp_td_results_file_params[i].params_name
    file_gp = joinpath( root_dir, "NBPMscape/paper/analysis/2_sampling_analysis/covid_like/gp_td_small"
                        , "small_$(gp_td_results_file_params[i].file_name)")
    # Read CSV into a DataFrame and store in the dictionary with the scenario number (numeric suffix) as key
    gp_tds_df_dict[params_name] = CSV.read(file_gp, DataFrame)[:,[:GP_TD,:GP_3TD]]
end

# CHECK gp_tds_df_dict["gp_tds_winter_300samples"]

# Save as files
# Primary care (GP)
@save joinpath( root_dir, "NBPMscape/paper/analysis/2_sampling_analysis/covid_like/"
                , "covidlike-1.6.1_nrep1000_4300448_gp_tds_df_dict.jld2") gp_tds_df_dict


# Read in dictionary containing data for ICU sampling
# ICU sampling simulated using code in:
# - 'NBPMscape/paper/analysis/2_sampling_analysis/scripts/1a_geogpraphic_optimisation_part_1.jl'
# - 'NBPMscape/paper/analysis/2_sampling_analysis/scripts/1b_site_selection_2026_10.R'
# - 'NBPMscape/paper/analysis/2_sampling_analysis/scripts/1c_geogpraphic_optimisation_part_2.jl'

# Includes results from four sampling site selection strategies including mSCAPE sites
icu_tds_df_dict = load(joinpath( root_dir, "NBPMscape/paper/analysis/2_sampling_analysis/covid_like/geo_optimisation"
                                ,"covidlike-1.6.1_nrep1000_4300448_icu_tds_df_dict.jld2")
                       , "icu_tds_df_dict")


###############################################################################################
# (2) Add outbreak metrics at TD values                       
###############################################################################################

### Compute outbreak count metrics at TD for different sampling parameter sets

# Loop through dictionary and add columns to the new dfs in dictionary 
# for the cumulative outbreak count metrics

## Load outbreak metrics dictionary 
## Created using 'NBPMscape/scripts/paper/2_sampling_analysis/1_geographic_optimisation.jl'
outbreak_metrics_time_series = load( joinpath( root_dir, "NBPMscape/paper/analysis/1_outbreak_simulations/covid_like/output"
                                              , "covidlike-1.6.1_nrep1000_4300448_outbreak_metrics_time_series.jld2") 
                                    , "outbreak_metrics_time_series" )

# Loop through the primary care sampling scenarios and 
# add outbreak metrics at the TD in each simulation replicate
for j in 1:length(gp_td_results_file_params) # keys(icu_tds_df_dict) # j=1
    key = gp_td_results_file_params[j].params_name
    df = gp_tds_df_dict[key]
    
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
        td = df[i,:GP_TD]
        td3 = df[i,:GP_3TD]
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
    gp_tds_df_dict[key] = df
end

# Save as file
# Primary care (GP)
@save joinpath( root_dir, "NBPMscape/paper/analysis/2_sampling_analysis/covid_like/"
                , "covidlike-1.6.1_nrep1000_4300448_gp_tds_df_dict.jld2") gp_tds_df_dict
# Reload if necessary
#gp_tds_df_dict = load( joinpath( root_dir, "NBPMscape/paper/analysis/2_sampling_analysis/covid_like/"
#                               , "covidlike-1.6.1_nrep1000_4300448_gp_tds_df_dict.jld2")
#                       , "gp_tds_df_dict")


##########################################################################################################
# (3) Replicate the datasets but with scenario number labelling rather than scenario description
##########################################################################################################

# Helper function to replace dictionary keys based on a lookup dictionary
function replace_keys(;dict::Dict, lookup::Dict)
    # Create a new dictionary to store updated keys
    new_dict = Dict{Any, Any}()

    for (key, value) in dict
        # If key exists in lookup, replace it; otherwise keep original
        new_key = get(lookup, key, key)

        # Handle potential duplicate keys after replacement
        if haskey(new_dict, new_key)
            # Merge or overwrite — here we overwrite for simplicity
            @warn "Duplicate key after replacement: $new_key. Overwriting value."
        end

        new_dict[new_key] = value
    end

    return new_dict
end

# Lookup dictionary for GP TD sampling scenario names and numbers
gp_td_scenarios = Dict(   "gp_tds_winter_300samples"  => 1 #"1"
                        , "gp_tds_winter_500samples"  => 13 #"13"
                        , "gp_tds_winter_1000samples" => 25 #"25"
                        , "gp_tds_summer_300samples"  => 37 #"37"
                        , "gp_tds_summer_500samples"  => 49 #"49"
                        , "gp_tds_summer_1000samples" => 61 #"61"
                        )

# Create dictionary with scenario numbers as keys
gp_tds_df_dict_scenario_numbers = replace_keys(dict = gp_tds_df_dict, lookup = gp_td_scenarios)


# Lookup dictionary for ICU TD sampling scenario names and numbers
icu_td_scenarios = Dict(  
                        # 8 sampling sites
                          "icu_tds_seasonal_ari_1440_300_sampling_sites_mscape_8sites"   => 1 # "1"
                        , "icu_tds_seasonal_ari_1440_500_sampling_sites_mscape_8sites"   => 13 #"13"
                        , "icu_tds_seasonal_ari_1440_1000_sampling_sites_mscape_8sites"  => 25 #"25"
                        , "icu_tds_seasonal_ari_793_300_sampling_sites_mscape_8sites"    => 37 #"37"
                        , "icu_tds_seasonal_ari_793_500_sampling_sites_mscape_8sites"    => 49 #"49"
                        , "icu_tds_seasonal_ari_793_1000_sampling_sites_mscape_8sites"   => 61 #"61"
                        # 26 sampling sites
                        , "icu_tds_seasonal_ari_1440_300_sampling_sites_mscape_26sites"  => 73 #"73"
                        , "icu_tds_seasonal_ari_1440_500_sampling_sites_mscape_26sites"  => 74 #"74"
                        , "icu_tds_seasonal_ari_1440_1000_sampling_sites_mscape_26sites" => 75 #"75"
                        , "icu_tds_seasonal_ari_793_300_sampling_sites_mscape_26sites"   => 76 #"76"
                        , "icu_tds_seasonal_ari_793_500_sampling_sites_mscape_26sites"   => 77 #"77"
                        , "icu_tds_seasonal_ari_793_1000_sampling_sites_mscape_26sites"  => 78 #"78"
                        )

# Trim ICU TD results dictionary to mSCAPE only
icu_tds_df_dict_mscape = Dict()
for key in keys(icu_tds_df_dict) # key="icu_tds_seasonal_ari_793_1000_sampling_sites_5blocks_20sites"
    if occursin("mscape",key)
        icu_tds_df_dict_mscape[key] = icu_tds_df_dict[key]
        println(key)
    end
end

# Create dictionary with scenario numbers as keys
icu_tds_df_dict_mscape_scenario_numbers = replace_keys(dict = icu_tds_df_dict_mscape, lookup = icu_td_scenarios)

# CHECK # icu_tds_df_dict_mscape_scenario_numbers[1] == icu_tds_df_dict_mscape["icu_tds_seasonal_ari_1440_300_sampling_sites_mscape_8sites"]  # true
# CHECK # icu_tds_df_dict_mscape_scenario_numbers[76] == icu_tds_df_dict_mscape["icu_tds_seasonal_ari_793_300_sampling_sites_mscape_26sites"] # true

# Save as files
# Primary care (GP)
@save joinpath( root_dir, "NBPMscape/paper/analysis/2_sampling_analysis/covid_like/"
                , "covidlike-1.6.1_nrep1000_4300448_gp_tds_df_dict_scenario_numbers.jld2") gp_tds_df_dict_scenario_numbers
# ICU mscape only
@save joinpath( root_dir, "NBPMscape/paper/analysis/2_sampling_analysis/covid_like/"
                , "covidlike-1.6.1_nrep1000_4300448_icu_tds_df_dict_mscape_scenario_numbers.jld2") icu_tds_df_dict_mscape_scenario_numbers


##########################################################################################################
# (4) Compute median time to detection (TD) values for all sampling scenarios
##########################################################################################################

# TD (1 case)
# Compute median TDs and 95% CI using bootstrapping for median TDs
icu_td_medians_dict = Dict{Int,NamedTuple}()
for (k,v) in icu_tds_df_dict_mscape_scenario_numbers #icu_tds_df_dict_scenario_numbers # v=icu_tds_df_dict_mscape_scenario_numbers["1"][:,:ICU_TD]
    icu_td_medians_dict[k] = median_ci_bootstrap( vec = v[:,:ICU_TD], n_boot = 1000, alpha = 0.05 ) # [ median, lower, upper ]
end 

gp_td_medians_dict = Dict{Int,NamedTuple}()
for (k,v) in gp_tds_df_dict_scenario_numbers # v=sc_tds_df_dict[1]
    gp_td_medians_dict[k] = median_ci_bootstrap( vec = v[:,:GP_TD], n_boot = 1000, alpha = 0.05 ) # [ median, lower, upper ]
end


##########################################################################################################
# (4) Plot median time to detection values for all sampling scenarios
#   (4.1) Time to detect 1 case (TD)
##########################################################################################################

## Plot all median TDs together

y_labels = [ # 300 samples, winter
            "Critical (8 sites)"    #"Tertiary (8 sites)" #"ICU 8 sites:  300 samples, winter"
           ,"Critical (26 sites)"   #"Tertiary (26 sites)"#"ICU 29 sites: 300 samples, winter"
           ,"Primary"               #"RCGP:                               300 samples, winter"
           ,""
           # 500 samples, winter
           ,"Critical (8 sites)"    #"ICU 8 sites:  500 samples, winter"
           ,"Critical (26 sites)"   #"ICU 29 sites: 500 samples, winter"
           ,"Primary"               #"RCGP:         500 samples, winter"
           ,""
           # 1000 samples, winter
           ,"Critical (8 sites)" #"ICU 8 sites:  1000 samples, winter"
           ,"Critical (26 sites)"#"ICU 29 sites: 1000 samples, winter"
           ,"Primary"            #"RCGP:         1000 samples, winter"
           ,""
           # 300 samples, summer
           ,"Critical (8 sites)" #"ICU 8 sites:  300 samples, summer"
           ,"Critical (26 sites)"#"ICU 29 sites: 300 samples, summer"
           ,"Primary"            #"RCGP:         300 samples, summer"
           ,""
           # 500 samples, summer
           ,"Critical (8 sites)" #"ICU 8 sites:  500 samples, summer"
           ,"Critical (26 sites)"#"ICU 29 sites: 500 samples, summer"
           ,"Primary"#"RCGP: 500 samples, summer"
           ,""
           # 1000 samples, summer
           ,"Critical (8 sites)" #"ICU 8 sites:  1000 samples, summer"
           ,"Critical (26 sites)"#"ICU 29 sites: 1000 samples, summer"
           ,"Primary"#"RCGP: 1000 samples, summer"
           ]

           # Critical care (8 sites)       # Critical care (26 sites)      # Primary care
medians = [icu_td_medians_dict[ 1].median, icu_td_medians_dict[73].median, gp_td_medians_dict[ 1].median 
          ,0.0
          ,icu_td_medians_dict[13].median, icu_td_medians_dict[74].median, gp_td_medians_dict[13].median 
          ,0.0
          ,icu_td_medians_dict[25].median, icu_td_medians_dict[75].median, gp_td_medians_dict[25].median 
          ,0.0
          ,icu_td_medians_dict[37].median, icu_td_medians_dict[76].median, gp_td_medians_dict[37].median 
          ,0.0
          ,icu_td_medians_dict[49].median, icu_td_medians_dict[77].median, gp_td_medians_dict[49].median 
          ,0.0
          ,icu_td_medians_dict[61].median, icu_td_medians_dict[78].median, gp_td_medians_dict[61].median 
          ]
                           # Critical care (8 sites)      # Critical care (26 sites)     # Primary care
lower_errors = medians .- [icu_td_medians_dict[ 1].lower, icu_td_medians_dict[73].lower, gp_td_medians_dict[ 1].lower 
                          ,0.0
                          ,icu_td_medians_dict[13].lower, icu_td_medians_dict[74].lower, gp_td_medians_dict[13].lower 
                          ,0.0
                          ,icu_td_medians_dict[25].lower, icu_td_medians_dict[75].lower, gp_td_medians_dict[25].lower 
                          ,0.0
                          ,icu_td_medians_dict[37].lower, icu_td_medians_dict[76].lower, gp_td_medians_dict[37].lower 
                          ,0.0
                          ,icu_td_medians_dict[49].lower, icu_td_medians_dict[77].lower, gp_td_medians_dict[49].lower 
                          ,0.0
                          ,icu_td_medians_dict[61].lower, icu_td_medians_dict[78].lower, gp_td_medians_dict[61].lower
                          ]
                # Critical care (8 sites)      # Critical care (26 sites)     # Primary care
upper_errors = [icu_td_medians_dict[ 1].upper, icu_td_medians_dict[73].upper, gp_td_medians_dict[ 1].upper 
               ,0.0
               ,icu_td_medians_dict[13].upper, icu_td_medians_dict[74].upper, gp_td_medians_dict[13].upper
               ,0.0
               ,icu_td_medians_dict[25].upper, icu_td_medians_dict[75].upper, gp_td_medians_dict[25].upper
               ,0.0
               ,icu_td_medians_dict[37].upper, icu_td_medians_dict[76].upper, gp_td_medians_dict[37].upper
               ,0.0
               ,icu_td_medians_dict[49].upper, icu_td_medians_dict[77].upper, gp_td_medians_dict[49].upper
               ,0.0
               ,icu_td_medians_dict[61].upper, icu_td_medians_dict[78].upper, gp_td_medians_dict[61].upper 
               ] .- medians

# Critical care (8 sites)
critical_8sites_colour = :lightblue
critical_26sites_colour = :blue
primarY_colour = :red
plot_colours = [critical_8sites_colour, critical_26sites_colour, primarY_colour]

colors = vcat( plot_colours
                ,:black
                ,plot_colours
                ,:black
                ,plot_colours
                ,:black
                ,plot_colours
                ,:black
                ,plot_colours
                ,:black
                ,plot_colours
                )

# Plot with error bars
# (CI using samples from fitted distribution probably not that useful here)
p_scatter_td = Plots.scatter(
                        medians
                        ,1:length(y_labels)
                        ,xerror = (lower_errors, upper_errors)
                        #, markercolor = colors
                        #, markerstrokecolor = colors
                        , color = colors
                        ,seriestype = :scatter
                        ,yticks = (1:length(y_labels), y_labels)
                        ,xlabel = "Median time to detection (days)\n (with 95% CI via 1000 bootstraps)"
                        ,xlims = [30,60] #[0,110]
                        ,xticks = 0:2:100 #0:10:100
                        ,ylabel = ""
                        #,title = "Time to detection\n with different site selection"
                        ,legend = false
                        ,yflip = true
                    )
Plots.plot(p_scatter_td, layout = (1, 1)
        , size = (900,600)) #default size is 600,400

Plots.savefig( joinpath( root_dir
                        , "NBPMscape/paper/analysis/2_sampling_analysis/covid_like/healthcare_setting"
                        , "scatter_all_TD_generic_names_GP_ICU.png" ) )
# Narrow plot
Plots.plot(p_scatter_td, layout = (1, 1)
        , size = (600,600)) #default size is 600,400
Plots.savefig( joinpath( root_dir
                        , "NBPMscape/paper/analysis/2_sampling_analysis/covid_like/healthcare_setting"
                        , "scatter_all_TD_narrow_generic_names_GP_ICU.png" ) )


#####################################################################################################################################
# (4) Plot median time to detection values for all sampling scenarios
#     (4.2) Time to detect 3 cases (3TD)
#####################################################################################################################################

# TD3 (3 cases)
# Compute median TDs and 95% CI using bootstrapping for median TDs
icu_td3_medians_dict = Dict{Int,NamedTuple}()
for (k,v) in icu_tds_df_dict_mscape_scenario_numbers #icu_tds_df_dict_scenario_numbers # v=sc_tds_df_dict[1]
    icu_td3_medians_dict[k] = median_ci_bootstrap( vec = v[:,:ICU_3TD], n_boot = 1000, alpha = 0.05 ) # [ median, lower, upper ]
end 

gp_td3_medians_dict = Dict{Int,NamedTuple}()
for (k,v) in gp_tds_df_dict_scenario_numbers # v=sc_tds_df_dict[1]
    gp_td3_medians_dict[k] = median_ci_bootstrap( vec = v[:,:GP_3TD], n_boot = 1000, alpha = 0.05 ) # [ median, lower, upper ]
end

medians_td3 = [  icu_td3_medians_dict[ 1].median, icu_td3_medians_dict[73].median, gp_td3_medians_dict[ 1].median 
                ,0.0
                ,icu_td3_medians_dict[13].median, icu_td3_medians_dict[74].median, gp_td3_medians_dict[13].median 
                ,0.0
                ,icu_td3_medians_dict[25].median, icu_td3_medians_dict[75].median, gp_td3_medians_dict[25].median 
                ,0.0
                ,icu_td3_medians_dict[37].median, icu_td3_medians_dict[76].median, gp_td3_medians_dict[37].median 
                ,0.0
                ,icu_td3_medians_dict[49].median, icu_td3_medians_dict[77].median, gp_td3_medians_dict[49].median 
                ,0.0
                ,icu_td3_medians_dict[61].median, icu_td3_medians_dict[78].median, gp_td3_medians_dict[61].median 
                ]
lower_errors_td3 = medians_td3 .- [  icu_td3_medians_dict[ 1].lower, icu_td3_medians_dict[73].lower, gp_td3_medians_dict[ 1].lower 
                                    ,0.0
                                    ,icu_td3_medians_dict[13].lower, icu_td3_medians_dict[74].lower, gp_td3_medians_dict[13].lower 
                                    ,0.0
                                    ,icu_td3_medians_dict[25].lower, icu_td3_medians_dict[75].lower, gp_td3_medians_dict[25].lower 
                                    ,0.0
                                    ,icu_td3_medians_dict[37].lower, icu_td3_medians_dict[76].lower, gp_td3_medians_dict[37].lower 
                                    ,0.0
                                    ,icu_td3_medians_dict[49].lower, icu_td3_medians_dict[77].lower, gp_td3_medians_dict[49].lower 
                                    ,0.0
                                    ,icu_td3_medians_dict[61].lower, icu_td3_medians_dict[78].lower, gp_td3_medians_dict[61].lower
                                    ]
upper_errors_td3 = [ icu_td3_medians_dict[ 1].upper, icu_td3_medians_dict[73].upper, gp_td3_medians_dict[ 1].upper 
                    ,0.0
                    ,icu_td3_medians_dict[13].upper, icu_td3_medians_dict[74].upper, gp_td3_medians_dict[13].upper
                    ,0.0
                    ,icu_td3_medians_dict[25].upper, icu_td3_medians_dict[75].upper, gp_td3_medians_dict[25].upper
                    ,0.0
                    ,icu_td3_medians_dict[37].upper, icu_td3_medians_dict[76].upper, gp_td3_medians_dict[37].upper
                    ,0.0
                    ,icu_td3_medians_dict[49].upper, icu_td3_medians_dict[77].upper, gp_td3_medians_dict[49].upper
                    ,0.0
                    ,icu_td3_medians_dict[61].upper, icu_td3_medians_dict[78].upper, gp_td3_medians_dict[61].upper 
                    ] .- medians_td3

# Plot with error bars
# (CI using samples from fitted distribution probably not that useful here)
p_scatter_td3 = Plots.scatter(
    medians_td3
    ,1:length(y_labels)
    ,xerror = (lower_errors_td3, upper_errors_td3)
    #, markercolor = colors
    #, markerstrokecolor = colors
    , color = colors
    ,seriestype = :scatter
    ,yticks = (1:length(y_labels), y_labels)
    ,xlabel = "Median time to detect 3 cases (days)\n (with 95% CI via 1000 bootstraps)"
    ,xlims = [30,60] #[0,110]
    ,xticks = 0:2:100
    ,ylabel = ""
    #,title = "Time to detection\n with different site selection"
    ,legend = false
    ,yflip = true
)
Plots.plot(p_scatter_td3, layout = (1, 1)
        , size = (900,600)) #default size is 600,400

Plots.savefig( joinpath( root_dir
                        , "NBPMscape/paper/analysis/2_sampling_analysis/covid_like/healthcare_setting"
                        , "scatter_all_TD3_generic_names_GP_ICU.png" ) )
# Narrow plot
Plots.plot(p_scatter_td3, layout = (1, 1)
        , size = (600,600)) #default size is 600,400
Plots.savefig( joinpath( root_dir
                        , "NBPMscape/paper/analysis/2_sampling_analysis/covid_like/healthcare_setting"
                        , "scatter_all_TD3_narrow_generic_names_GP_ICU.png" ) )



##################################################################################################################################
# (5) Analysis of multi-modal sampling - sampling from combinations of healthcare settings 
#    (5.1) Computing combined times to detection
##################################################################################################################################

# Match scenario numbers (from config file suffixes) between sampling setting (see 'NBPMscape/config/)
scenario_numbers_icu_8 =  [ 1,13,25,37,49,61]
scenario_numbers_icu_29 = [73,74,75,76,77,78]
scenario_numbers_gp =     [ 1,13,25,37,49,61]

# Prepare containers to store the individual DataFrames
icu_8_tds_df_dict_v2  = Dict{Int, DataFrame}()
icu_29_tds_df_dict_v2 = Dict{Int, DataFrame}()
gp_tds_df_dict_v2     = Dict{Int, DataFrame}()

# Prepare a df to hold the median TD values and upper and lower 95% CI from bootstrapping
default_int_col   = zeros(Int,     length(scenario_numbers_gp))
default_float_col = zeros(Float64, length(scenario_numbers_gp))
tds_df = DataFrame(
                    Scenario_number         = copy(default_int_col)
                    ,Scenario_number_icu_29 = copy(default_int_col) # ICU sampling at 29 sites uses different config files because it requires a different parameter value for the ICU sites
                    ## Median
                    # TD
                    ,icu_8_median_TD        = copy(default_float_col)
                    ,icu_29_median_TD       = copy(default_float_col)
                    ,gp_median_TD           = copy(default_float_col)
                    # 3TD
                    ,icu_8_median_3TD       = copy(default_float_col)
                    ,icu_29_median_3TD      = copy(default_float_col)
                    ,gp_median_3TD          = copy(default_float_col)
                    # Combined TD
                    ,gp_icu_8_comb_median_TD     = copy(default_float_col)      # GP + ICU 8
                    ,gp_icu_29_comb_median_TD    = copy(default_float_col)      # GP + ICU 29
                    # Combined 3TD
                    ,gp_icu_8_comb_median_3TD     = copy(default_float_col)     # GP + ICU 8
                    ,gp_icu_29_comb_median_3TD    = copy(default_float_col)     # GP + ICU 29
                    ## Lower CI
                    # TD
                    ,icu_8_median_TD_lower  = copy(default_float_col)
                    ,icu_29_median_TD_lower = copy(default_float_col)
                    ,gp_median_TD_lower     = copy(default_float_col)
                    # 3TD
                    ,icu_8_median_3TD_lower  = copy(default_float_col)
                    ,icu_29_median_3TD_lower = copy(default_float_col)
                    ,gp_median_3TD_lower     = copy(default_float_col)
                    # Combined TD
                    ,gp_icu_8_comb_median_TD_lower     = copy(default_float_col) # GP + ICU 8
                    ,gp_icu_29_comb_median_TD_lower    = copy(default_float_col) # GP + ICU 29
                    # Combined 3TD
                    ,gp_icu_8_comb_median_3TD_lower     = copy(default_float_col) # GP + ICU 8
                    ,gp_icu_29_comb_median_3TD_lower    = copy(default_float_col) # GP + ICU 29
                    ## Upper CI
                    # TD
                    ,icu_8_median_TD_upper  = copy(default_float_col)
                    ,icu_29_median_TD_upper = copy(default_float_col)
                    ,gp_median_TD_upper     = copy(default_float_col)
                    # 3TD
                    ,icu_8_median_3TD_upper  = copy(default_float_col)
                    ,icu_29_median_3TD_upper = copy(default_float_col)
                    ,gp_median_3TD_upper     = copy(default_float_col)
                    # Combined TD
                    ,gp_icu_8_comb_median_TD_upper     = copy(default_float_col) # GP + ICU 8
                    ,gp_icu_29_comb_median_TD_upper    = copy(default_float_col) # GP + ICU 29
                    # Combined 3TD
                    ,gp_icu_8_comb_median_3TD_upper     = copy(default_float_col) # GP + ICU 8
                    ,gp_icu_29_comb_median_3TD_upper    = copy(default_float_col) # GP + ICU 29
                    ## Lower IQR
                    # TD
                    ,icu_8_TD_iqr_lower  = copy(default_float_col)
                    ,icu_29_TD_iqr_lower = copy(default_float_col)
                    ,gp_TD_iqr_lower     = copy(default_float_col)
                    # 3TD
                    ,icu_8_3TD_iqr_lower  = copy(default_float_col)
                    ,icu_29_3TD_iqr_lower = copy(default_float_col)
                    ,gp_3TD_iqr_lower     = copy(default_float_col)
                    # Combined TD
                    ,gp_icu_8_comb_TD_iqr_lower     = copy(default_float_col) # GP + ICU 8
                    ,gp_icu_29_comb_TD_iqr_lower    = copy(default_float_col) # GP + ICU 29
                    # Combined 3TD
                    ,gp_icu_8_comb_3TD_iqr_lower     = copy(default_float_col) # GP + ICU 8
                    ,gp_icu_29_comb_3TD_iqr_lower    = copy(default_float_col) # GP + ICU 29
                    ## Upper IQR
                    # TD
                    ,icu_8_TD_iqr_upper  = copy(default_float_col)
                    ,icu_29_TD_iqr_upper = copy(default_float_col)
                    ,gp_TD_iqr_upper     = copy(default_float_col)
                    # 3TD
                    ,icu_8_3TD_iqr_upper  = copy(default_float_col)
                    ,icu_29_3TD_iqr_upper = copy(default_float_col)
                    ,gp_3TD_iqr_upper     = copy(default_float_col)
                    # Combined TD
                    ,gp_icu_8_comb_TD_iqr_upper     = copy(default_float_col) # GP + ICU 8
                    ,gp_icu_29_comb_TD_iqr_upper    = copy(default_float_col) # GP + ICU 29
                    # Combined 3TD
                    ,gp_icu_8_comb_3TD_iqr_upper     = copy(default_float_col) # GP + ICU 8
                    ,gp_icu_29_comb_3TD_iqr_upper    = copy(default_float_col) # GP + ICU 29
                    ## Standard Deviation
                    # TD
                    ,icu_8_TD_std  = copy(default_float_col)
                    ,icu_29_TD_std = copy(default_float_col)
                    ,gp_TD_std     = copy(default_float_col)
                    # 3TD
                    ,icu_8_3TD_std  = copy(default_float_col)
                    ,icu_29_3TD_std = copy(default_float_col)
                    ,gp_3TD_std     = copy(default_float_col)
                    # Combined TD
                    ,gp_icu_8_comb_TD_std     = copy(default_float_col) # GP + ICU 8
                    ,gp_icu_29_comb_TD_std    = copy(default_float_col) # Gp + ICU 29
                    # Combined 3TD
                    ,gp_icu_8_comb_3TD_std     = copy(default_float_col) # GP + ICU 8
                    ,gp_icu_29_comb_3TD_std    = copy(default_float_col) # GP + ICU 29
                    )


## Loop through scenarios, read in files and compute median TD values

# ICU 8 sites
for i in 1:length(scenario_numbers_icu_8) # i=1
    k = scenario_numbers_icu_8[i]
    #file = joinpath( path_name_icu
    #                , "icu_tds_$(string(k)).csv")

    # Read CSV into a DataFrame and store in the dictionary with the scenario number (numeric suffix) as key
    #icu_8_tds_df_dict_v2[k] = CSV.read(file, DataFrame)
    
    for key in keys(icu_tds_df_dict_mscape_scenario_numbers) # key =1
        if in(key, scenario_numbers_icu_8)
            icu_8_tds_df_dict_v2[key] = icu_tds_df_dict_mscape_scenario_numbers[key]
        end
    end
                        
    # Compute times to detection and populate df
    tds_df[i,"icu_8_median_TD"]       = median(icu_8_tds_df_dict_v2[k][:,"ICU_TD"])
    tds_df[i,"icu_8_median_3TD"]      = median(icu_8_tds_df_dict_v2[k][:,"ICU_3TD"])
    tds_df[i,"icu_8_TD_iqr_lower"]  = quantile(icu_8_tds_df_dict_v2[k][:,"ICU_TD"],  0.25)
    tds_df[i,"icu_8_3TD_iqr_lower"] = quantile(icu_8_tds_df_dict_v2[k][:,"ICU_3TD"], 0.25)
    tds_df[i,"icu_8_TD_iqr_upper"]  = quantile(icu_8_tds_df_dict_v2[k][:,"ICU_TD"],  0.75)
    tds_df[i,"icu_8_3TD_iqr_upper"] = quantile(icu_8_tds_df_dict_v2[k][:,"ICU_3TD"], 0.75)
    tds_df[i,"icu_8_TD_std"]             = std(icu_8_tds_df_dict_v2[k][:,"ICU_TD"] )
    tds_df[i,"icu_8_3TD_std"]            = std(icu_8_tds_df_dict_v2[k][:,"ICU_3TD"] )

end

# ICU 29 sites
for i in 1:length(scenario_numbers_icu_29) 
    m = scenario_numbers_icu_29[i]
    #file = joinpath( path_name_icu
    #                , "icu_tds_$(string(m)).csv")

    # Read CSV into a DataFrame and store in the dictionary with the scenario number (numeric suffix) as key
    #icu_29_tds_df_dict_v2[m] = CSV.read(file, DataFrame)
    
    for key in keys(icu_tds_df_dict_mscape_scenario_numbers) # key =1
        if in(key, scenario_numbers_icu_29)
            icu_29_tds_df_dict_v2[key] = icu_tds_df_dict_mscape_scenario_numbers[key]
        end
    end

    # Compute times to detection and populate df
    tds_df[i,"Scenario_number_icu_29"] = m
    tds_df[i,"icu_29_median_TD"]       = median(icu_29_tds_df_dict_v2[m][:,"ICU_TD"])
    tds_df[i,"icu_29_median_3TD"]      = median(icu_29_tds_df_dict_v2[m][:,"ICU_3TD"])
    tds_df[i,"icu_29_TD_iqr_lower"]  = quantile(icu_29_tds_df_dict_v2[m][:,"ICU_TD"],  0.25)
    tds_df[i,"icu_29_3TD_iqr_lower"] = quantile(icu_29_tds_df_dict_v2[m][:,"ICU_3TD"], 0.25)
    tds_df[i,"icu_29_TD_iqr_upper"]  = quantile(icu_29_tds_df_dict_v2[m][:,"ICU_TD"],  0.75)
    tds_df[i,"icu_29_3TD_iqr_upper"] = quantile(icu_29_tds_df_dict_v2[m][:,"ICU_3TD"], 0.75)
    tds_df[i,"icu_29_TD_std"]             = std(icu_29_tds_df_dict_v2[m][:,"ICU_TD"] )
    tds_df[i,"icu_29_3TD_std"]            = std(icu_29_tds_df_dict_v2[m][:,"ICU_3TD"] )
end

# GP
for i in 1:length(scenario_numbers_gp) 
    n = scenario_numbers_gp[i]
    #file = joinpath( path_name_gp
    #                , "gp_tds_$(string(n)).csv")

    # Read CSV into a DataFrame and store in the dictionary with the scenario number (numeric suffix) as key
    #gp_tds_df_dict_v2[n] = CSV.read(file, DataFrame)
    
    gp_tds_df_dict_v2 = gp_tds_df_dict_scenario_numbers                    
    
    # Compute times to detection and populate df
    tds_df[i,"gp_median_TD"]       = median(gp_tds_df_dict_v2[n][:,"GP_TD"])
    tds_df[i,"gp_median_3TD"]      = median(gp_tds_df_dict_v2[n][:,"GP_3TD"])
    tds_df[i,"gp_TD_iqr_lower"]  = quantile(gp_tds_df_dict_v2[n][:,"GP_TD"], 0.25)
    tds_df[i,"gp_3TD_iqr_lower"] = quantile(gp_tds_df_dict_v2[n][:,"GP_3TD"], 0.25)
    tds_df[i,"gp_TD_iqr_upper"]  = quantile(gp_tds_df_dict_v2[n][:,"GP_TD"], 0.75)
    tds_df[i,"gp_3TD_iqr_upper"] = quantile(gp_tds_df_dict_v2[n][:,"GP_3TD"], 0.75)
    tds_df[i,"gp_TD_std"]             = std(gp_tds_df_dict_v2[n][:,"GP_TD"] )
    tds_df[i,"gp_3TD_std"]            = std(gp_tds_df_dict_v2[n][:,"GP_3TD"] )
end

#println(tds_df)

# Save dictionaries with scenario numbers for keys (already have dictionaries saved with string descriptions as keys)
dict_output = joinpath( root_dir
                        , "NBPMscape/paper/analysis/2_sampling_analysis/covid_like/healthcare_setting" )
#gp_tds_df_dict_scenario_numbers = gp_tds_df_dict_v2
icu_8_tds_df_dict_scenario_numbers = icu_8_tds_df_dict_v2
icu_26_tds_df_dict_scenario_numbers = icu_29_tds_df_dict_v2

@save joinpath( dict_output, "covidlike-1.6.1_nrep1000_4300448_gp_tds_df_dict_scenario_numbers.jld2")         gp_tds_df_dict_scenario_numbers
@save joinpath( dict_output, "covidlike-1.6.1_nrep1000_4300448_icu_8_tds_df_dict_scenario_numbers.jld2")   icu_8_tds_df_dict_scenario_numbers
@save joinpath( dict_output, "covidlike-1.6.1_nrep1000_4300448_icu_26_tds_df_dict_scenario_numbers.jld2") icu_26_tds_df_dict_scenario_numbers

# If necessary, read files back in 
#@load joinpath( dict_output, "covidlike-1.6.1_nrep1000_4300448_gp_tds_df_dict_scenario_numbers.jld2") gp_tds_df_dict_scenario_numbers
#@load joinpath( dict_output, "covidlike-1.6.1_nrep1000_4300448_icu_8_tds_df_dict_scenario_numbers.jld2") icu_8_tds_df_dict_scenario_numbers
#@load joinpath( dict_output, "covidlike-1.6.1_nrep1000_4300448_icu_26_tds_df_dict_scenario_numbers.jld2") icu_26_tds_df_dict_scenario_numbers

# Compute median TDs combining different sampling settings
for i in 1:length(scenario_numbers_gp) # i =1
    k = scenario_numbers_icu_8[i]
    m = scenario_numbers_icu_29[i]
    n = scenario_numbers_gp[i]

    # Populate scenario numbers
    tds_df[i,"Scenario_number"] = n # Assumes that k == n
    tds_df[i,"Scenario_number_icu_29"] = m

    # Times to detection (for each simulation replicate) for the different sample settings and for the specific parameter set
    # TD
    icu_8_tds  =  icu_8_tds_df_dict_scenario_numbers[ k ][:,"ICU_TD"]
    icu_29_tds = icu_26_tds_df_dict_scenario_numbers[ m ][:,"ICU_TD"]
    gp_tds     =     gp_tds_df_dict_scenario_numbers[ n ][:,"GP_TD" ]
    #println(sample_setting_1[1:10],sample_setting_2[1:10],combined_min_TD[1:10])
    # 3TD
    icu_8_tds_3  =  icu_8_tds_df_dict_scenario_numbers[k][:,"ICU_3TD"]
    icu_29_tds_3 = icu_26_tds_df_dict_scenario_numbers[m][:,"ICU_3TD"]
    gp_tds_3     =     gp_tds_df_dict_scenario_numbers[n][:,"GP_3TD" ]

    # Compute median times to detection and 95% CI using bootstrap
    # ICU 8 sites
    (tds_df[i,"icu_8_median_TD"],  tds_df[i,"icu_8_median_TD_lower"],  tds_df[i,"icu_8_median_TD_upper"])  = NBPMscape.median_ci_bootstrap( vec = icu_8_tds,   n_boot = 1000, alpha = 0.05)
    (tds_df[i,"icu_8_median_3TD"], tds_df[i,"icu_8_median_3TD_lower"], tds_df[i,"icu_8_median_3TD_upper"]) = NBPMscape.median_ci_bootstrap( vec = icu_8_tds_3, n_boot = 1000, alpha = 0.05)
    # ICU 29 sites
    (tds_df[i,"icu_29_median_TD"],  tds_df[i,"icu_29_median_TD_lower"],  tds_df[i,"icu_29_median_TD_upper"])  = NBPMscape.median_ci_bootstrap( vec = icu_29_tds,   n_boot = 1000, alpha = 0.05)
    (tds_df[i,"icu_29_median_3TD"], tds_df[i,"icu_29_median_3TD_lower"], tds_df[i,"icu_29_median_3TD_upper"]) = NBPMscape.median_ci_bootstrap( vec = icu_29_tds_3, n_boot = 1000, alpha = 0.05)
    # GP
    (tds_df[i,"gp_median_TD"],  tds_df[i,"gp_median_TD_lower"],  tds_df[i,"gp_median_TD_upper"])  = NBPMscape.median_ci_bootstrap( vec = gp_tds,   n_boot = 1000, alpha = 0.05)
    (tds_df[i,"gp_median_3TD"], tds_df[i,"gp_median_3TD_lower"], tds_df[i,"gp_median_3TD_upper"]) = NBPMscape.median_ci_bootstrap( vec = gp_tds_3, n_boot = 1000, alpha = 0.05)
    

    # Compute combined times to detection and populate df
    # GP and ICU 8 sites
    gp_icu_8_comb   = min.(gp_tds,   icu_8_tds  )
    gp_icu_8_comb_3 = min.(gp_tds_3, icu_8_tds_3  )

    (tds_df[i,"gp_icu_8_comb_median_TD"],  tds_df[i,"gp_icu_8_comb_median_TD_lower"],  tds_df[i,"gp_icu_8_comb_median_TD_upper"])  = NBPMscape.median_ci_bootstrap( vec = gp_icu_8_comb,   n_boot = 1000, alpha = 0.05)
    (tds_df[i,"gp_icu_8_comb_median_3TD"], tds_df[i,"gp_icu_8_comb_median_3TD_lower"], tds_df[i,"gp_icu_8_comb_median_3TD_upper"]) = NBPMscape.median_ci_bootstrap( vec = gp_icu_8_comb_3, n_boot = 1000, alpha = 0.05)
    
    tds_df[i,"gp_icu_8_comb_TD_iqr_lower"]  = quantile( gp_icu_8_comb, 0.25)
    tds_df[i,"gp_icu_8_comb_3TD_iqr_lower"] = quantile( gp_icu_8_comb_3, 0.25)
    tds_df[i,"gp_icu_8_comb_TD_iqr_upper"]  = quantile( gp_icu_8_comb, 0.75)
    tds_df[i,"gp_icu_8_comb_3TD_iqr_upper"] = quantile( gp_icu_8_comb_3, 0.75)
    tds_df[i,"gp_icu_8_comb_TD_std"]             = std( gp_icu_8_comb)
    tds_df[i,"gp_icu_8_comb_3TD_std"]            = std( gp_icu_8_comb_3)

    # GP and ICU 29 sites
    gp_icu_29_comb   = min.(gp_tds,   icu_29_tds  )
    gp_icu_29_comb_3 = min.(gp_tds_3, icu_29_tds_3  )

    (tds_df[i,"gp_icu_29_comb_median_TD"],  tds_df[i,"gp_icu_29_comb_median_TD_lower"],  tds_df[i,"gp_icu_29_comb_median_TD_upper"])  = NBPMscape.median_ci_bootstrap( vec = gp_icu_29_comb,   n_boot = 1000, alpha = 0.05)
    (tds_df[i,"gp_icu_29_comb_median_3TD"], tds_df[i,"gp_icu_29_comb_median_3TD_lower"], tds_df[i,"gp_icu_29_comb_median_3TD_upper"]) = NBPMscape.median_ci_bootstrap( vec = gp_icu_29_comb_3, n_boot = 1000, alpha = 0.05)
    
    tds_df[i,"gp_icu_29_comb_TD_iqr_lower"]  = quantile( gp_icu_29_comb, 0.25)
    tds_df[i,"gp_icu_29_comb_3TD_iqr_lower"] = quantile( gp_icu_29_comb_3, 0.25)
    tds_df[i,"gp_icu_29_comb_TD_iqr_upper"]  = quantile( gp_icu_29_comb, 0.75)
    tds_df[i,"gp_icu_29_comb_3TD_iqr_upper"] = quantile( gp_icu_29_comb_3, 0.75)
    tds_df[i,"gp_icu_29_comb_TD_std"]             = std( gp_icu_29_comb)
    tds_df[i,"gp_icu_29_comb_3TD_std"]            = std( gp_icu_29_comb_3)

end

CSV.write( joinpath( root_dir
                    , "NBPMscape/paper/analysis/2_sampling_analysis/covid_like/healthcare_setting"                    
                    , "tds_df_inc_comb_GP_ICU_inc_CI.csv") , tds_df ) 
# Reload if necessary
tds_df = CSV.read(joinpath( root_dir
                    , "NBPMscape/paper/analysis/2_sampling_analysis/covid_like/healthcare_setting"                    
                    , "tds_df_inc_comb_GP_ICU_inc_CI.csv")
                    , DataFrame)

##################################################################################################################################
# (5) Analysis of multi-modal sampling - sampling from combinations of healthcare settings 
#    (5.2) Plotting combined times to detection
##################################################################################################################################

## Plot all median TDs together

# TD (1 case)
# Create base labels and repeat 6 times (representing the 300, 500, 1000 samples in winter and summer), then drop the last "" (to match length)
y_label_base = [
                 "Primary"
                ,"Critical (8 sites)"
                ,"Critical (26 sites)"
                ,"Primary & Critical (8 sites)"
                ,"Primary & Critical (26 sites)"
                ,""
            ]
y_labels = vcat(repeat(y_label_base, 6)...)[1:end-1]

medians = [] # Scenario / config file 1
for i in 1:6 # 1 to 6 represents scenario numbers 1,13,25,37,48,61
    push!(medians, tds_df[i,:gp_median_TD])
    push!(medians, tds_df[i,:icu_8_median_TD])
    push!(medians, tds_df[i,:icu_29_median_TD])
    push!(medians, tds_df[i,:gp_icu_8_comb_median_TD])
    push!(medians, tds_df[i,:gp_icu_29_comb_median_TD])
    i < 6 ? push!(medians, 0.0) : missing
end

medians_lower = [] # Scenario / config file 1
for i in 1:6 # 1 to 6 represents scenario numbers 1,13,25,37,48,61
    push!(medians_lower, tds_df[i,:gp_median_TD_lower])
    push!(medians_lower, tds_df[i,:icu_8_median_TD_lower])
    push!(medians_lower, tds_df[i,:icu_29_median_TD_lower])
    push!(medians_lower, tds_df[i,:gp_icu_8_comb_median_TD_lower])
    push!(medians_lower, tds_df[i,:gp_icu_29_comb_median_TD_lower])
    i < 6 ? push!(medians_lower, 0.0) : missing
end

lower_errors = medians .- medians_lower

medians_upper = [] # Scenario / config file 1
for i in 1:6 # 1 to 6 represents scenario numbers 1,13,25,37,48,61
    push!(medians_upper, tds_df[i,:gp_median_TD_upper])
    push!(medians_upper, tds_df[i,:icu_8_median_TD_upper])
    push!(medians_upper, tds_df[i,:icu_29_median_TD_upper])
    push!(medians_upper, tds_df[i,:gp_icu_8_comb_median_TD_upper])
    push!(medians_upper, tds_df[i,:gp_icu_29_comb_median_TD_upper])
    i < 6 ? push!(medians_upper, 0.0) : missing
end

upper_errors = medians_upper .-  medians

# Define colours
pc     = "#E69F00" # "Primary"                       # orange 
cc8    = "#56B4E9" # "Critical (8 sites)"            # light blue 
cc26   = "#0072B2" # "Critical (26 sites)"           # dark blue 
cc8pc  = "#C7856A" # "Primary & Critical (8 sites)"  # muted salmon/brown (c2 + c3) # mix of light blue and orange
cc26pc = "#7A68A6" # "Primary & Critical (26 sites)" # muted purple (c1 + c3) # approximately blue + orange mixtures 

plot_colours = [pc, cc8, cc26, cc8pc, cc26pc]

colors = vcat(plot_colours
            ,:black
            ,plot_colours
            ,:black
            ,plot_colours
            ,:black
            ,plot_colours
            ,:black
            ,plot_colours
            ,:black
            ,plot_colours
            )
plot_markers = [:square, :triangle, :utriangle, :circle, :diamond]
markers = vcat(plot_markers
            ,:star
            ,plot_markers
            ,:star
            ,plot_markers
            ,:star
            ,plot_markers
            ,:star
            ,plot_markers
            ,:star
            ,plot_markers
            )

# Plot with error bars
# (CI using samples from fitted distribution probably not that useful here)
p_scatter_td = Plots.scatter(
                                medians
                                ,1:length(y_labels)
                                ,xerror = (lower_errors, upper_errors)
                                #, markercolor = colors
                                #, markerstrokecolor = colors
                                #,markershape = markers
                                , color = colors
                                ,seriestype = :scatter
                                ,yticks = (1:length(y_labels), y_labels)
                                ,xlabel = "Median time to detection (days)\n (with 95% CI via 1000 bootstraps)"
                                ,xlims = [30,60] #[0,110]
                                ,xticks = 0:2:100
                                ,ylabel = ""
                                #,title = "Time to detection\n with different site selection"
                                ,legend = false
                                ,yflip = true
                            )
Plots.plot(p_scatter_td, layout = (1, 1)
        , size = (900,600)) #default size is 600,400

Plots.savefig( joinpath( root_dir
                        , "NBPMscape/paper/analysis/2_sampling_analysis/covid_like/healthcare_setting"
                        , "scatter_TD_combined_examples_generic_names_GP_ICU.png" ) )
# Narrow plot
Plots.plot(p_scatter_td, layout = (1, 1)
        , size = (600,600)) #default size is 600,400
Plots.savefig( joinpath( root_dir
                        , "NBPMscape/paper/analysis/2_sampling_analysis/covid_like/healthcare_setting"
                        , "scatter_TD_combined_examples_narrow_generic_names_GP_ICU.png" ) )
Plots.savefig( joinpath( root_dir
                        , "NBPMscape/paper/results/figures_main/raw"
                        , "scatter_TD_combined_examples_narrow_generic_names_GP_ICU.png" ) )

# TD3 results scatter plot
medians_td3 = []
for i in 1:6 # 1 to 6 represents scenario numbers 1,13,25,37,48,61
    push!(medians_td3, tds_df[i,:gp_median_3TD])
    push!(medians_td3, tds_df[i,:icu_8_median_3TD])
    push!(medians_td3, tds_df[i,:icu_29_median_3TD])
    push!(medians_td3, tds_df[i,:gp_icu_8_comb_median_3TD])
    push!(medians_td3, tds_df[i,:gp_icu_29_comb_median_3TD])
    i < 6 ? push!(medians_td3, 0.0) : missing
end

medians_lower_td3 = []
for i in 1:6 # 1 to 6 represents scenario numbers 1,13,25,37,48,61
    push!(medians_lower_td3, tds_df[i,:gp_median_3TD_lower])
    push!(medians_lower_td3, tds_df[i,:icu_8_median_3TD_lower])
    push!(medians_lower_td3, tds_df[i,:icu_29_median_3TD_lower])
    push!(medians_lower_td3, tds_df[i,:gp_icu_8_comb_median_3TD_lower])
    push!(medians_lower_td3, tds_df[i,:gp_icu_29_comb_median_3TD_lower])
    i < 6 ? push!(medians_lower_td3, 0.0) : missing
end

lower_errors_td3 = medians_td3 .- medians_lower_td3

medians_upper_td3 = []
for i in 1:6 # 1 to 6 represents scenario numbers 1,13,25,37,48,61
    push!(medians_upper_td3, tds_df[i,:gp_median_3TD_upper])
    push!(medians_upper_td3, tds_df[i,:icu_8_median_3TD_upper])
    push!(medians_upper_td3, tds_df[i,:icu_29_median_3TD_upper])
    push!(medians_upper_td3, tds_df[i,:gp_icu_8_comb_median_3TD_upper])
    push!(medians_upper_td3, tds_df[i,:gp_icu_29_comb_median_3TD_upper])
    i < 6 ? push!(medians_upper_td3, 0.0) : missing
end

upper_errors_td3 = medians_upper_td3 .-  medians_td3

# Plot with error bars
# (CI using samples from fitted distribution probably not that useful here)
p_scatter_td3 = Plots.scatter(
                                medians_td3
                                ,1:length(y_labels)
                                ,xerror = (lower_errors_td3, upper_errors_td3)
                                #, markercolor = colors
                                #, markerstrokecolor = colors
                                , color = colors
                                ,seriestype = :scatter
                                ,yticks = (1:length(y_labels), y_labels)
                                ,xlabel = "Median time to detection (days)\n (with 95% CI via 1000 bootstraps)"
                                ,xlims = [30,60] #[0,110]
                                ,xticks = 0:2:100
                                ,ylabel = ""
                                #,title = "Time to detection\n with different site selection"
                                ,legend = false
                                ,yflip = true
                            )
Plots.plot(p_scatter_td3, layout = (1, 1)
        , size = (900,600)) #default size is 600,400

Plots.savefig( joinpath( root_dir
                        , "NBPMscape/paper/analysis/2_sampling_analysis/covid_like/healthcare_setting"
                        , "scatter_3TD_combined_examples_generic_names_GP_ICU.png" ) )
# Narrow plot
Plots.plot(p_scatter_td3, layout = (1, 1)
        , size = (600,600)) #default size is 600,400
Plots.savefig( joinpath( root_dir
                        , "NBPMscape/paper/analysis/2_sampling_analysis/covid_like/healthcare_setting"
                        , "scatter_3TD_combined_examples_narrow_generic_names_GP_ICU.png" ) )
Plots.savefig( joinpath( root_dir
                        , "NBPMscape/paper/results/figures_main/raw"
                        , "scatter_3TD_combined_examples_narrow_generic_names_GP_ICU.png" ) )


##########################################################################################################
# (6) Specific values for article text
##########################################################################################################

# Reload if necessary
tds_df = CSV.read(joinpath( root_dir
                    , "NBPMscape/paper/analysis/2_sampling_analysis/covid_like/healthcare_setting"                    
                    , "tds_df_inc_comb_GP_ICU_inc_CI.csv")
                    , DataFrame)

# mSCAPE 26 sites, winter, 300 samples per week
tds_df[1,:].icu_29_median_TD # 44.65567253018948
tds_df[1,:].icu_29_TD_iqr_lower # 39.489891806084536
tds_df[1,:].icu_29_TD_iqr_upper # 50.59601633405625

# mSCAPE 26 sites, summer, 300 samples per week
tds_df[4,:].icu_29_median_TD # 41.89225443586476
tds_df[4,:].icu_29_TD_iqr_lower # 36.43208223570267
tds_df[4,:].icu_29_TD_iqr_upper # 47.702858574512206

# mSCAPE 8 sites, 300 samples per week
tds_df[1,:].icu_8_median_TD # 48.57293908282003
tds_df[4,:].icu_8_median_TD # 47.8929973601776

# Primary care, 300 samples per week
tds_df[1,:].gp_median_TD # 50.48174969904808
tds_df[4,:].gp_median_TD # 46.54820975275941

# Difference: mSCAPE 26 vs mSCAPE 8
tds_df[1,:].icu_29_median_TD - tds_df[1,:].icu_8_median_TD # -3.917266552630551
tds_df[4,:].icu_29_median_TD - tds_df[4,:].icu_8_median_TD # -6.000742924312846

# Difference: mSCAPE 26 vs Primary care
tds_df[1,:].icu_29_median_TD - tds_df[1,:].gp_median_TD # -5.826077168858603
tds_df[4,:].icu_29_median_TD - tds_df[4,:].gp_median_TD # -4.655955316894655

# mSCAPE 26 sites, summer, 500 samples per week
tds_df[5,:].icu_29_median_TD    # 40.79436087628645
tds_df[5,:].icu_29_TD_iqr_lower # 34.35226552820704
tds_df[5,:].icu_29_TD_iqr_upper # 46.108976356047314

# mSCAPE 26 sites, summer, 500 samples per week
tds_df[6,:].icu_29_median_TD    # 40.27908714315577
tds_df[6,:].icu_29_TD_iqr_lower # 34.69237409717453
tds_df[6,:].icu_29_TD_iqr_upper # 45.69875518071501

# Difference
tds_df[6,:].icu_29_median_TD - tds_df[5,:].icu_29_median_TD # -0.5152737331306767

# mSCAPE 8 sites, all scenarios
tds_df.icu_8_median_TD 
#6-element Vector{Float64}:
# 48.57293908282003
# 48.2961045814483
# 47.999858927516826
# 47.8929973601776
# 47.73818030044728
# 47.71537087106285
minimum(tds_df.icu_8_median_TD) # 47.71537087106285
maximum(tds_df.icu_8_median_TD) # 48.57293908282003

# Relative advantage of mSCAPE v Primary care excluding sampling saturation (i.e. 1000 samples per week in summer)
# Absolute differences
tds_df[1:5,:].icu_29_median_TD - tds_df[1:5,:].gp_median_TD
#5-element Vector{Float64}:
# -5.826077168858603
# -4.909213046758985
# -3.913632885649953
# -4.655955316894655
# -3.7228992941189603

# Absolute differences
tds_df[1:5,:].icu_29_median_TD ./ tds_df[1:5,:].gp_median_TD
#5-element Vector{Float64}:
# 0.8845904271624629
# 0.8962070409597128
# 0.9112271237293075
# 0.8999756308217923
# 0.9163717784996592