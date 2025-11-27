module NBPMscape

using JumpProcesses 
using DifferentialEquations
using Random
using Distributions
using HypothesisTests
using DataFrames
import UUIDs 
import StatsBase
#using Statistics
using Interpolations
import SpecialFunctions as SF 
using Plots 
using LinearAlgebra
#using Pkg
#Pkg.add("NamedArrays")
using NamedArrays

using RData 
using CSV 
# using Pkg.Artifacts

using GLM
using StatsPlots
using JLD2

const COMMUTERPROBPATH = joinpath(@__DIR__, "..", "data",  "commuting_ITL2_prob_list.rds")
const COMMUTEPROB = load( COMMUTERPROBPATH )   
const COMMUTERINPROBPATH = joinpath(@__DIR__, "..", "data",  "commuting_ITL2_inprob_list.rds")
const COMMUTEINPROB = load( COMMUTERINPROBPATH )  
const COMMUTERMPATH = joinpath(@__DIR__, "..", "data",  "commuting_ITL2_list.rds")
const COMMUTERM = load( COMMUTERMPATH )  
const REGKEYPATH = joinpath( @__DIR__, "..", "data", "ITL2_key2.rds" )
const REGKEY = load( REGKEYPATH )
const CAAPATH = joinpath( @__DIR__, "..", "data", "CAA_pax_2024_ITL2.rds" )
const CAAIMPORTS = load(CAAPATH)
# itl2size = CSV.read(joinpath( @__DIR__, "..", "data", "itl2_regions_england.csv"), DataFrame)
# const ITL2SIZE = filter( r->r.Code in REGKEY.code, itl2size )
itl2size =  load( joinpath( @__DIR__, "..", "data", "itl2_population2022.rds" ) )
#= > head(d1) 
  ITL225CD                               ITL225NM total_population_2022
1     TLC3                            Tees Valley                688756 =#
const ITL2SIZE = filter( r->r.ITL225CD in REGKEY.code, itl2size )
#ITL2SIZE[37:39,:ITL225CD]
#ITL225CD_wales = ITL2SIZE[37:39,:ITL225CD]
#ITL2SIZE_eng = ITL2SIZE[.!in(ITL2SIZE.ITL225CD, ITL225CD_wales), :]
#pop_eng = sum(ITL2SIZE[1:36,3])

### Load international traveller age groups (note this includes returning UK residents and visiting overseas residents)
# ONS travel pac data from International Passenger Survey for 2019
const INT_TRAVELLERS_AGE_GROUP = load( joinpath( @__DIR__, "..", "data", "travelpac_2019_age_group_weights.rds" ) ) #int_travellers_age_group
const INT_TRAVELLERS_AGE_SINGLE_YR = load( joinpath( @__DIR__, "..", "data", "travelpac_2019_age_single_year_weights.rds" ) ) #int_travellers_age_single_yr

### Load contact distributions
# Generated using socialmixr R package and UK POLYMOD survey data
# Artificial split by single age
#contact_distributions_age_single_yr = load( joinpath( @__DIR__, "..", "data", "polymod_contact_age_setting_distribution.rds" ) )
contact_distributions_age_group = load( joinpath( @__DIR__, "..", "data", "contact_setting_age_group_distributions.rds" ) )
# Remove the 'all' age group as not required
const CONTACT_DISTRIBUTIONS = filter(:age_group => !=("all"), contact_distributions_age_group)

### Load contact matrices (disaggregated by age group)
# Generated using socialmixr R package and UK POLYMOD survey data
# TODO double check and remove any double counting of individual contacts made in multiple settings
const CONTACT_MATRIX_AGE_GROUPS = load( joinpath( @__DIR__, "..", "data", "polymod_contact_matrix_age_groups.rds" ) )
#contact_matrix_home_df = load( joinpath( @__DIR__, "..", "data", "polymod_contact_matrix_home.rds" ) , convert = true )
contact_matrix_home = load( joinpath( @__DIR__, "..", "data", "polymod_contact_matrix_home.rds" ) , convert = true )
const CONTACT_MATRIX_HOME = NamedArray( Matrix( contact_matrix_home ) 
                                      , (CONTACT_MATRIX_AGE_GROUPS, CONTACT_MATRIX_AGE_GROUPS))
contact_matrix_school_work = load( joinpath( @__DIR__, "..", "data", "polymod_contact_matrix_school_work.rds" ) )
const CONTACT_MATRIX_SCHOOL_WORK = NamedArray( Matrix( contact_matrix_school_work ) 
                                              , (CONTACT_MATRIX_AGE_GROUPS, CONTACT_MATRIX_AGE_GROUPS) )
contact_matrix_other = load( joinpath( @__DIR__, "..", "data", "polymod_contact_matrix_other.rds" ) )
const CONTACT_MATRIX_OTHER = NamedArray( Matrix( contact_matrix_other ) 
                                        , (CONTACT_MATRIX_AGE_GROUPS, CONTACT_MATRIX_AGE_GROUPS) )

### Load severity probabilities
# Extracted from Knock et al (2021), Key Epidemiological Drivers and Impact of Interventions in the 2020 SARS-CoV-2 Epidemic in England’.
# Science Translational Medicine, 16(602), eabg4262
symptomatic_prob_by_age = load( joinpath( @__DIR__, "..", "data", "knock_symptomatic_prob_by_age.rds" ) )
symptomatic_prob_by_age = filter(row -> !(row.age_single_year in ["Care_home_workers", "Care_home_residents"])
                                      , symptomatic_prob_by_age)
const SYMPTOMATIC_PROB_BY_AGE = symptomatic_prob_by_age[:,["age_single_year","symptomatic_prob"]]
ihr_ifr_by_age = load( joinpath( @__DIR__, "..", "data", "knock_ihr_ifr_by_age.rds" ) )
ihr_ifr_by_age = filter(row -> !(row.age_single_year in ["Care_home_residents"])
                      , ihr_ifr_by_age)
ihr_ifr_by_age[:,(3:8)] = ihr_ifr_by_age[:,(3:8)] ./ 100 # Convert from percentage to proportion
const IHR_BY_AGE = ihr_ifr_by_age[:,["age_single_year","IHR"]]
const IFR_BY_AGE = ihr_ifr_by_age[:,["age_single_year","IFR"]]

#const IFR_BY_AGE = parse.(Float64, ihr_ifr_by_age.age_single_year)","IFR"])
#care_pathway_prob_by_age
const CARE_PATHWAY_PROB_BY_AGE = load( joinpath( @__DIR__, "..", "data", "knock_care_pathway_prob_by_age.rds" ) )

### Load probabilities of visiting a particular NHS Trust given residency in a particular ITL2 region
# Computed from NHS Acute (Hospitals) Trust Catchment Populations compiled by Office for Health Improvement & Disparities 
# at https://app.powerbi.com/view?r=eyJrIjoiODZmNGQ0YzItZDAwZi00MzFiLWE4NzAtMzVmNTUwMThmMTVlIiwidCI6ImVlNGUxNDk5LTRhMzUtNGIyZS1hZDQ3LTVmM2NmOWRlODY2NiIsImMiOjh9
# using 'Hospital Episode Statistics Admitted Patient Care (HES APC)' data for 'Emergency' admissions
# National Statistics Postcode Lookup (NSPL) used to convert from MSOA to ITL2 regions.
itl2_to_nhs_trust_prob_adult = CSV.read( joinpath( @__DIR__, "..", "data", "NHS_Trust_prob_from_ITL2_regions_adult_2020.csv" ), DataFrame )
rename!(itl2_to_nhs_trust_prob_adult, :Column1 => :NHS_Trust_code)
const ITL2_TO_NHS_TRUST_PROB_ADULT = itl2_to_nhs_trust_prob_adult
itl2_to_nhs_trust_prob_child = CSV.read( joinpath( @__DIR__, "..", "data", "NHS_Trust_prob_from_ITL2_regions_child_2020.csv" ), DataFrame )
rename!(itl2_to_nhs_trust_prob_child, :Column1 => :NHS_Trust_code)
const ITL2_TO_NHS_TRUST_PROB_CHILD = itl2_to_nhs_trust_prob_child
#TEST wsample( ITL2_TO_NHS_TRUST_PROB_CHILD[:, :NHS_Trust_code], ITL2_TO_NHS_TRUST_PROB_CHILD[:, :TLI3], 1)

### Load probabilities of sampling at particular NHS Trust
nhs_trust_icu_sample_prob = CSV.read( joinpath( @__DIR__, "..", "data", "selected_ICU_prop_NHS_Trust.csv" ), DataFrame )
rename!(nhs_trust_icu_sample_prob, Symbol("Org.Code") => :NHS_Trust_code)
rename!(nhs_trust_icu_sample_prob, Symbol("Org.Name") => :NHS_Trust_name)
const NHS_TRUST_ICU_SAMPLE_PROB = nhs_trust_icu_sample_prob

### Load probabilities of being sampled in ICU given ITL2 home region
#itl2_icu_sample_prob = CSV.read( joinpath( @__DIR__, "..", "data", "ICU_prob_sample_by_ITL2_2019.csv" ), DataFrame )
#itl2_icu_sample_prob = CSV.read( joinpath( @__DIR__, "..", "data", "ICU_prob_sample_by_ITL2_2020.csv" ), DataFrame )
itl2_icu_sample_prob = CSV.read( joinpath( @__DIR__, "..", "data", "ICU_prob_sample_by_ITL2_2020_mscapeOct25_sitrepJan25.csv" ), DataFrame )
rename!(itl2_icu_sample_prob, Symbol("Column1") => :site_stage_age) 
const ITL2_ICU_SAMPLE_PROB = itl2_icu_sample_prob

### Load critical care bed data
# Downloaded from https://www.england.nhs.uk/statistics/statistical-work-areas/uec-sitrep/
cc_bed_sitrep = CSV.read( joinpath( @__DIR__, "..", "data", "202501-January-2025-beds-sitrep-data-FINAL-revised.csv" ), DataFrame )
# Filter for NHS Trusts 
cc_bed_sitrep_nt = filter(row -> !ismissing(row[:Level]) && row[:Level] == "Provider", cc_bed_sitrep)
# Filter for critical care beds
metric_filtered = [ "Adult critical care beds available" #"Adult critical care beds occupied","Adult critical care occupancy rate"
                   ,"Paediatric intensive care beds available"]#, "Paediatric intensive care beds occupied", "Paediatric intensive care occupancy rate"
cc_bed_sitrep_nt = filter(row -> !ismissing(row[:Metric]) && in( row[:Metric], metric_filtered), cc_bed_sitrep_nt)
#unique(cc_bed_sitrep_nt[:,:Metric])
# Trim df
cc_bed_sitrep_nt = cc_bed_sitrep_nt[:,3:end]
# Reformat so Adult and Paediatric bed numbers are in separate columns
rename!(cc_bed_sitrep_nt, "Org Code" => :NHS_Trust_code)
rename!(cc_bed_sitrep_nt, "Org Name" => :NHS_Trust_name)
cc_bed_sitrep_nt.Value = parse.(Int, cc_bed_sitrep_nt.Value)
cc_bed_sitrep_nt = unstack(cc_bed_sitrep_nt, [:Region,:ICB, :NHS_Trust_code, :NHS_Trust_name, :Type], :Metric, :Value)
rename!(cc_bed_sitrep_nt, "Adult critical care beds available" => :Adult_critical_care_beds_available)
rename!(cc_bed_sitrep_nt, "Paediatric intensive care beds available" => :Paediatric_intensive_care_beds_available)
# Remove specialist NHS Trusts that are unlikely to admit ARI to ICU
cc_bed_sitrep_nt = filter(row -> in(row[:NHS_Trust_code], ITL2_TO_NHS_TRUST_PROB_ADULT[:,:NHS_Trust_code]), cc_bed_sitrep_nt) # Note that ITL2_TO_NHS_TRUST_PROB_ADULT[:,:NHS_Trust_code] == ITL2_TO_NHS_TRUST_PROB_CHILD[:,:NHS_Trust_code]
## Compute proportion of total critical care beds at each NHS Trust and add to df
cc_beds_total_adult = sum(cc_bed_sitrep_nt.Adult_critical_care_beds_available) 
cc_beds_total_child = sum(cc_bed_sitrep_nt.Paediatric_intensive_care_beds_available) 
cc_bed_sitrep_nt.p_adult = cc_bed_sitrep_nt.Adult_critical_care_beds_available ./ cc_beds_total_adult
cc_bed_sitrep_nt.p_child = cc_bed_sitrep_nt.Paediatric_intensive_care_beds_available ./ cc_beds_total_child
const ARI_CC_BED_SITREP = cc_bed_sitrep_nt

export simtree, simforest, sampleforest, simgendist, Infection, infectivitytoR
export transmissionrate, sampdegree, REGKEY, COMMUTEPROB #TODO 
include("core.jl")

include("combine_small_sim_reps.jl")
export combine_sim_reps, n_sims_w_results, n_sims_w_results_icu_gp

include("sims_filter.jl")
export sims_filter

include("sampling_infections.jl")
export icu_v_pc_td, analyse_td_columns, plot_hist

include("icu_sample_prob_functions.jl")
export sample_nhs_trust, icu_sample_prob, icu_sample_prob_region, sample_icu_cases, sample_icu_cases_n

include("estimate_severity_age_weights.jl")
export inf_severity_estimate, inf_age_estimate

end 
