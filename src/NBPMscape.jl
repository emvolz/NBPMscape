module NBPMscape

# using Pkg.Artifacts
using CSV 
using DataFrames
using DataFramesMeta
using Dates
using DifferentialEquations
using Distributions
using ForwardDiff
using GLM
using HypothesisTests
using Interpolations
using JLD2
using JumpProcesses 
using LinearAlgebra
using NamedArrays
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

const COMMUTERPROBPATH = joinpath(@__DIR__, "..", "data/commuting",  "commuting_ITL2_prob_list.rds")
const COMMUTEPROB = load( COMMUTERPROBPATH )   
const COMMUTERINPROBPATH = joinpath(@__DIR__, "..", "data/commuting",  "commuting_ITL2_inprob_list.rds")
const COMMUTEINPROB = load( COMMUTERINPROBPATH )  
const COMMUTERMPATH = joinpath(@__DIR__, "..", "data/commuting",  "commuting_ITL2_list.rds")
const COMMUTERM = load( COMMUTERMPATH )  
const REGKEYPATH = joinpath( @__DIR__, "..", "data/commuting", "ITL2_key2.rds" )
const REGKEY = load( REGKEYPATH )
const CAAPATH = joinpath( @__DIR__, "..", "data/air_traffic", "CAA_pax_2024_ITL2.rds" )
const CAAIMPORTS = load(CAAPATH)

### Load international traveller age groups (note this includes returning UK residents and visiting overseas residents)
# ONS travel pac data from International Passenger Survey for 2019
const INT_TRAVELLERS_AGE_GROUP = load( joinpath( @__DIR__, "..", "data/international_travel", "travelpac_2019_age_group_weights.rds" ) ) #int_travellers_age_group
const INT_TRAVELLERS_AGE_SINGLE_YR = load( joinpath( @__DIR__, "..", "data/international_travel", "travelpac_2019_age_single_year_weights.rds" ) ) #int_travellers_age_single_yr


## Constants
export REGKEY, COMMUTEPROB

include("contact_data.jl")
export CONTACT_DISTRIBUTIONS, CONTACT_MATRIX_AGE_GROUPS
export CONTACT_MATRIX_HOME, CONTACT_MATRIX_SCHOOL_WORK, CONTACT_MATRIX_OTHER

include("nhs_trust_data.jl")
export NHS_TRUST_CATCHMENT_POP, NHS_TRUST_CATCHMENT_POP_ADULT_CHILD
export ITL2_TO_NHS_TRUST_PROB_ADULT, ITL2_TO_NHS_TRUST_PROB_CHILD
export ARI_CC_BED_SITREP, AE_12M
#export NHS_TRUST_ICU_SAMPLE_PROB, ITL2_ICU_SAMPLE_PROB

include("NHS_Wales_ICU_data.jl")
export NHS_WALES_ICU_BEDS #NHS_WALES_ICU_BEDS_ADULT, NHS_WALES_ICU_BEDS_CHILD
#export ITL2_TO_NHS_TRUST_PROB_ADULT_EW, ITL2_TO_NHS_TRUST_PROB_CHILD_EW, ARI_CC_BED_SITREP_EW

include("population_data.jl")
export ITL2SIZE

include("severity_care_probabilities.jl")
export SYMPTOMATIC_PROB_BY_AGE, IHR_BY_AGE, IFR_BY_AGE, CARE_PATHWAY_PROB_BY_AGE

include("load_import_data.jl")
export GLEAM_DAILY_IMPORT_RATES

## Functions
include("core.jl")
export simtree, simforest, sampleforest, simgendist, Infection, infectivitytoR
export transmissionrate, sampdegree
export initialize_parameters, create_default_parameters

include("combine_small_sim_reps.jl")
export combine_sim_reps, n_sims_w_results, n_sims_w_results_icu_gp

include("sims_filter.jl")
export sims_filter

include("sampling_infections.jl")
export icu_td, gp_td, secondary_care_td, icu_v_pc_td, analyse_td_columns, plot_hist

include("icu_sample_prob_functions.jl")
export sample_nhs_trust, icu_sample_prob, icu_sample_prob_region, sample_icu_cases, sample_icu_cases_n

include("hosp_sampling_functions.jl")
export sample_hosp_cases_n, courier_collection_days

include("estimate_severity_age_weights.jl")
export inf_severity_estimate, inf_age_estimate

include("distribution_fitting.jl")
export fit_multi_dist, dist_fit_plot
export nll_trunc_gamma, discretize_gamma_pmf, assign_bins, nll_disc_gamma, nll_trunc_weibull, nll_trunc_normal
export erlang_truncated_means, gamma_params_from_mode_cdf

include("outbreak_timeline_analysis_functions.jl")
export median_td_sim_reps, cases_before_td, cases_before_td_v2, infections_time_analysis, infections_time_analysis_v2
export local_doubling_time, global_doubling_time, _fit_log_linear, _fit_weighted_log_linear, _growth_rate_to_doubling_time, _fit_nonlinear_optim, _estimate_initial_r
export plot_outbreak_analysis, plot_outbreak_analysis_combined

include("misc_functions.jl")
export median_ci_bootstrap, allocate_with_rounding, generation_time, severity_rolling_mean, tinf_by_age


## Configure parameter input values
include("config.jl")
export load_config, validate_config, update_configurable_parameters, print_changes, convert_params_to_dfs

initialize_parameters();

end 
