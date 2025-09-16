module NBPMscape

using JumpProcesses 
using DifferentialEquations
using Random
using Distributions
using DataFrames
import UUIDs 
import StatsBase
using Statistics
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

export simtree, simforest, sampleforest, simgendist, Infection, infectivitytoR
export transmissionrate, sampdegree, REGKEY, COMMUTEPROB #TODO 
include("core.jl")

include("sampling_infections.jl")
export icu_v_pc_td, analyse_td_columns, plot_hist

include("sims_filter.jl")
export sims_filter

include("estimate_severity_weights.jl")
export inf_severity_estimate

include("combine_small_sim_reps.jl")
export combine_sim_reps, n_sims_w_results, n_sims_w_results_icu_gp

end 
