#=
Investigating time to detection using metagenomic sampling from a network of hospitals based on existing HARISS network
December 2025
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
import SpecialFunctions as SF 
import StatsBase
using Statistics
using StatsPlots
import UUIDs 
using XLSX

## Load simulated outbreak data
#sims = combine_sim_reps(;  # sim_input_folders in vector format where multiple paths
#                          sim_input_folders = [ "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1240854/1240854_filtered_G/" ]
#                        , sim_object_name = "sims_G_filtered" #"sims" # "sims_G_icu_filter" # "sims_G_gp_filter" # "sim_G_filtered"
#                        , nrep = 50
#                        )
#@save "covidlike-1.4.1-sims-filtered_G_nrep1850_1240854.jld2" sims
sims = load("covidlike-1.4.1-sims-filtered_G_nrep1850_1240854.jld2", "sims")

# Test sample_hosp_cases_n on simulated data
df = sample_hosp_cases_n(; p = NBPMscape.P
                                , hosp_cases = sims[1] # hosp_cases = sims[1].G
                                , hosp_ari_admissions = 2000 * 7 # Estimate for winter based on Figure 3 in https://assets.publishing.service.gov.uk/media/691edba0ee5aaf516084d41f/UKHSA-emergency-department-syndromic-surveillance-bulletin-2025-week-46.pdf [Accessed: 21 Nov 2025]
                                , hosp_ari_admissions_adult_p = 0.52 # Proportion of ED ARI admissions that are adults (16y and over) - estimated from UKHSA Emergency Department Syndromic Surveillance System 
                                , hosp_ari_admissions_child_p = 0.48 # Proportion of ED ARI admissions that are children (<16y) - estimated from UKHSA Emergency Department Syndromic Surveillance System 
                                , n_hosp_samples_per_week = 300 # Total number of hospital samples to be taken per week for metagenomic testing
                                , sample_allocation = "equal" # "equal" or "weighted"
                                , sample_proportion_adult = "free" # "free" or numeric decimal, e.g. 0.75. Indicates split of sample target between adults and children.
                                , hariss_nhs_trust_sampling_sites = CSV.read("C:/Users/kdrake/Documents/mSCAPE/2_sampling/HARRIS/hariss_nhs_trust_sampling_sites.csv", DataFrame)
                                , nhs_trust_catchment_pop = NHS_TRUST_CATCHMENT_POP_ADULT_CHILD
                                #, nhs_trust_ae_12m = AE_12M
                                , weight_samples_by = "ae_mean" # or "catchment_pop"
                                , phl_collection_dow = [2,5] # [Mon, Thur]. Day(s) of the week that swab samples will be collected from public health labs. Day of week codes: Sunday = 1,... Saturday = 7.
                                , ed_discharge_limit = 0.25 # days. Assume that people attending the Emergency Department are discharged within this time limit.
                                # Summer
                                #, ed_ari_destinations_adult = DataFrame( destination = [:discharged,:short_stay,:longer_stay]
                                #                                       , proportion_of_attendances = [0.575,0.034,0.391])
                                #, ed_ari_destinations_child = DataFrame( destination = [:discharged,:short_stay,:longer_stay]
                                #                                       , proportion_of_attendances = [0.851,0.012,0.137])
                                # Winter
                                , ed_ari_destinations_adult = DataFrame( destination = [:discharged,:short_stay,:longer_stay]
                                                                       , proportion_of_attendances = [0.628,0.030,0.342])
                                , ed_ari_destinations_child = DataFrame( destination = [:discharged,:short_stay,:longer_stay]
                                                                       , proportion_of_attendances = [0.861,0.014,0.125])
                                , swab_time_mode = 0.25 # Assume swabbing peaks at 6hrs (=0.25 days) after attendance/admission at hospital
                                , swab_proportion_at_48h = 0.9 # Assume 90% of swabs are taken within 48hrs (=2 days) of attendance/admission at hospital
                                , proportion_hosp_swabbed = 0.9 # Assume 90% of ARI attendances are swabbed
                                , initial_dow = 1 # Day of the week that initial case imported. Day of week codes: Sunday =1, Monday = 2, ..., Saturday = 7 
                                , only_sample_before_death = true # There is a possibilty of swabbing time being drawn after death so 'true' here will constrain tswab to tdeceased
                                )

#TEST # hosp_cases = sims[173]

# Compute times to detection in secondary care setting (e.g. via HARISS network)
sc_tds = secondary_care_td(; p = NBPMscape.P
                            , sims = sims #[1025:1850]
                            , pathogen_type = "virus"
                            , initial_dow = 1
                            # Sampling parameters
                            , hariss_courier_to_analysis = 1.0 # Time between courier collection from PHL and beginning of analysis
                            , hariss_turnaround_time = NBPMscape.P.turnaroundtime_hariss #[2,4] # Time to process sample and report results / declare detection
                            , n_hosp_samples_per_week = 300 # Total number of hospital samples to be taken per week
                            , sample_allocation = "equal" # "equal" or "weighted"
                            , sample_proportion_adult = "free" # "free" or numeric decimal, e.g. 0.75. Indicates split of sample target 
                                                                    # between adults and children. "free" indicates that no split is specified
                            , hariss_nhs_trust_sampling_sites = CSV.read("C:/Users/kdrake/Documents/mSCAPE/2_sampling/HARRIS/hariss_nhs_trust_sampling_sites.csv", DataFrame) # List of NHS Trusts in HARISS sampling network 
                            , weight_samples_by = "ae_mean" # or "catchment_pop"
                            , phl_collection_dow = [2,5] # Day(s) of week that swab samples will be collected from public health labs. Day of week codes: Sunday = 1,... Saturday = 7.
                            , swab_time_mode = 0.25 # Assume swabbing peaks at 6hrs (=0.25 days) after attendance/admission at hospital
                            , swab_proportion_at_48h = 0.9 # Assume 90% of swabs are taken within 48hrs (=2 days) of attendance/admission at hospital
                            , proportion_hosp_swabbed = 0.9 # Assume 90% of ARI attendances are swabbed
                            #, initial_dow::Int64 = 1 # Day of the week that initial case imported. Day of week codes: Sunday =1, Monday = 2, ..., Saturday = 7 
                            , only_sample_before_death = true # There is a possibilty of swabbing time being drawn after death so 'true' here will constrain tswab to tdeceased
                            # Hospital parameters
                            , ed_discharge_limit = 0.25 # days. Assume that people attending the Emergency Department are discharged within this time limit.
                            , nhs_trust_catchment_pop = NHS_TRUST_CATCHMENT_POP_ADULT_CHILD
                            #, nhs_trust_ae_12m = AE_12M
                            # Seasonal values
                            # Winter
                            , hosp_ari_admissions = Int64( round( 79148 / ((31+31+29)/7), digits = 0 ) ) # for winter and Int64(45360 / ((30+31+31)/7)) for summer. Estimate of weekly hospital ARI admissions (excluding pathogen X being simulated) - using Dec 2023, Jan 2024, Feb 2025 data for winter and Jun, Jul, Aug 2024 for summer
                            , hosp_ari_admissions_adult_p = 0.52 # for winter and 0.58 for summer.Proportion of ED ARI admissions that are adults (16y and over)
                            , hosp_ari_admissions_child_p = 0.48 # for winter and 0.42 for summer. Proportion of ED ARI admissions that are children (<16y)
                            , ed_ari_destinations_adult = DataFrame( destination = [:discharged,:short_stay,:longer_stay]
                                                                   , proportion_of_attendances = [0.628,0.030,0.342])
                            , ed_ari_destinations_child = DataFrame( destination = [:discharged,:short_stay,:longer_stay]
                                                                   , proportion_of_attendances = [0.861,0.014,0.125])
                            # Summer
                            #, hosp_ari_admissions = Int64( round( 45360 / ((30+31+31)/7), digits = 0 ) ) # for summer and Int64(79148 / ((31+31+29)/7)) # for winter. Estimate of weekly hospital ARI admissions (excluding pathogen X being simulated) - using Dec 2023, Jan 2024, Feb 2025 data for winter and Jun, Jul, Aug 2024 for summer
                            #, hosp_ari_admissions_adult_p = 0.58 # for summer and 0.52 for winter. Proportion of ED ARI admissions that are adults (16y and over)
                            #, hosp_ari_admissions_child_p = 0.42 # for summer and 0.48 for winter. Proportion of ED ARI admissions that are children (<16y)
                            #, ed_ari_destinations_adult = DataFrame( destination = [:discharged,:short_stay,:longer_stay]
                            #                                       , proportion_of_attendances = [0.575,0.034,0.391])
                            #, ed_ari_destinations_child = DataFrame( destination = [:discharged,:short_stay,:longer_stay]
                            #                                       , proportion_of_attendances = [0.851,0.012,0.137])
                            )


sc_median = median( sc_tds[:,:SC_TD]) # 53.6 days
CSV.write( "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1240854_nrep1850_analysis/sc_tds.csv", sc_tds )

# Compute times to detection in tertiary care setting
icu_tds = NBPMscape.icu_td(;  p=NBPMscape.P
                            , sims = sims
                            , icu_sample_type = "regional" # "regional" or "fixed". If "fixed" then p_icu will be used, note that this doesn't take into account test sensitivity or practical sampling proportion, which "regional" does by using {sample_icu_cases} function
                            , pathogen_type = "virus"
                            , site_stage = "current" 
                            , sample_icu_cases_version = "number" # "proportion"
                            , p_icu = 0.15 # ICU sampling proportion
                            , n_icu_samples_per_week = 300
                            , icu_ari_admissions = 1440 # [793, 1440] # Weekly ICU admission numbers [summer,winter]
                            , icu_turnaround_time = [2,4] # Time to process sample and report results / declare detection
                            , icu_ari_admissions_adult_p = 0.76 # Proportion of ICU ARI admissions that are adults (16y and over)
                            , icu_ari_admissions_child_p = 0.24 # Proportion of ICU ARI admissions that are children (<16y)
                            , nhs_trust_sampling_sites = CSV.read( "C:/Users/kdrake/Documents/mSCAPE/2_sampling/ICU/nhs_trust_site_sample_targets.csv" , DataFrame )
                            )

icu_median = median( icu_tds[:,:ICU_TD]) # 71.6 days
CSV.write( "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1240854_nrep1850_analysis/icu_tds.csv", icu_tds )

# Compute times to detection in primary care setting (e.g. via RCGP network)
gp_tds = NBPMscape.gp_td(;  p=NBPMscape.P
                            , sims = sims 
                            # Parameters for existing Oxford-RCGP RSC primary care surveillance
                            , gp_practices_total = 6199 # Total number of GP practices in England at July 2025. Source: BMA analysis (https://www.bma.org.uk/advice-and-support/nhs-delivery-and-workforce/pressures/pressures-in-general-practice-data-analysis) of NHS Digital General Practice Workforce Statistics (https://digital.nhs.uk/data-and-information/publications/statistical/general-and-personal-medical-services) [Accessed 2 Sep 2025]  
                            , gp_practices_swab = 300 # Number of GP practices taking swabs for virology surveillance. Source: Data quality report: national flu and COVID-19 surveillance report (27 May 2025)
                            , gp_swabs_mg = 300 # [319, 747]#[25698,46685] # Assumed number of swabs that are metagenomic sequenced for investigating impact
                            , pop_eng = 5.7106398e7 # Population of England. Source: ONS mid-year 2022 UK population data disaggregated by various geo levels and age groups and gender. [Accessed 6 November 2024] Available at https://www.ons.gov.uk/peoplepopulationandcommunity/populationandmigration/populationestimates/datasets/populationestimatesforukenglandandwalesscotlandandnorthernireland 
                            , gp_ari_consults = 327 # 180, 327 # Number of ARI consultations per 100k of England population per week [mean summer 2024, mean winter 2024/25]. Source: Analysis of data extracted from RCGP Research & Surveillance Centre (RSC) Virology Dashboard [Accessed 29 Aug 2025]
                            , gp_ari_swabs = 747 #[319, 747] #[25698,46685]# Number of swabs taken from suspected ARI per week [mean summer 2024, mean winter 2024/25]. Source: Analysis of data extracted from RCGP Research & Surveillance Centre (RSC) Virology Dashboard [Accessed 29 Aug 2025]
                            , pc_swab_turnaround_time = [2,4] # [min,max] number of days between swab sample being taken and results received. Source: Data quality report: national flu and COVID-19 surveillance report (27 May 2025).  Assume same for metagenomic testing.
                            , pathogen_type = "virus"
                            )

gp_median = median( gp_tds[:,:GP_TD]) # 72.5 days
CSV.write( "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1240854_nrep1850_analysis/gp_tds.csv", gp_tds )


## Plot results
# Secondary care results
histogram( sc_tds[:,:SC_TD], alpha= 0.5, label = "HARISS", xlabel = "Time since first import to first detection (days)", ylabel = "Number of sim reps")
vline!( [sc_median], color=:blue
            , label="HARISS median"
            , linewidth = 3) 
# Add tertiary care TD results
histogram!( icu_tds[:,:ICU_TD], alpha =0.5, label="ICU")
vline!( [icu_median], color=:green
            , label="ICU median"
            , linewidth = 3) 
# Add primary care TD results
histogram!( gp_tds[:,:GP_TD], alpha =0.5, label = "RCGP")
vline!( [gp_median], color=:red
            , label="RCGP median"
            , linewidth = 3)
# Save plot
Plots.savefig("C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1240854_nrep1850_analysis/histogram.png")