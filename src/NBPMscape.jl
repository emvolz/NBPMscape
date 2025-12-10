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
using Plots.PlotMeasures 
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

using XLSX
using QuadGK

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



## Load Emergency Department (ED) attendances for 12 months from April 2024 to March 2025 (most recent as at Nov 2025)
## Can be used to allocate estimated ED ARI cases to NHS Trusts
# Load data
#folder_name = "C:/Users/kdrake/OneDrive - Imperial College London/Documents/mSCAPE/2_sampling/HARRIS/NHS_England_AE/"
files_vec = ["Monthly-AE-April-2024-revised.csv","Monthly-AE-May-2024.csv","Monthly-AE-June-2024.csv","Monthly-AE-July-2024-revised.csv","Monthly-AE-August-2024.csv"
            ,"Monthly-AE-September-2024.csv","Monthly-AE-October-2024-revised.csv","Monthly-AE-November-2024-revised.csv","Monthly-AE-December-2024.csv","Monthly-AE-January-2025.csv"
            ,"Monthly-AE-February-2025-revised.csv","Monthly-AE-March-2025.csv"]
# Read all files into a vector of DataFrames
dfs = [CSV.read(joinpath( @__DIR__, "..", "data/NHS_England_AE", f ), DataFrame)[:, [:("Org Code"), :("A&E attendances Type 1")]] for f in files_vec]
# Rename columns dynamically
for (i, df) in enumerate(dfs)
    month_label = i <= 9 ? "2024_$(i+3)" : "2025_$(i-9)"  # Currently for data from Apr 2024 to Mar 2025 - Adjust as necessary
    rename!(df, Symbol("A&E attendances Type 1") => Symbol(month_label))
end
# Merge each month's data into a single df
ae_12m = reduce((x, y) -> outerjoin(x, y, on = :("Org Code")), dfs)
# Remove rows for Totals
ae_12m = filter(row -> lowercase(strip(row[:"Org Code"])) != "total", ae_12m)
# Adjust for Trust merger
nt_changes_dict = Dict{String, String}()
#nt_changes_dict[("RAP","North Middlesex University Hospital NHS Trust")] = ("RAL","Royal Free London NHS Foundation Trust") # Merged with Royal Free in Jan 2025 https://www.royalfree.nhs.uk/news/royal-free-london-and-north-middlesex-university-hospital-to-merge-1-january-2025 
nt_changes_dict[("RAP")] = "RAL" # Merged with Royal Free in Jan 2025 https://www.royalfree.nhs.uk/news/royal-free-london-and-north-middlesex-university-hospital-to-merge-1-january-2025 
for r in eachrow(ae_12m)
    key = (r[:("Org Code")])
    if haskey(nt_changes_dict, key)
        new1 = nt_changes_dict[key]
        #println(new1)
        r[:("Org Code")] = new1
    end
end
ae_12m = combine( groupby(ae_12m, :("Org Code")), names(ae_12m, 2:13) .=> (x -> sum(skipmissing(x))))
# CHECK println( filter(:("Org Code") => ==("RAL"), ae_12m) ) # println( filter(:("Org Code") => ==("RAP"), ae_12m) )
# Compute 12-month mean and proportion of total 12-month mean
ae_12m.mean_12m = mean.(eachrow(ae_12m[:, 2:end]))
mean_total = sum(skipmissing(ae_12m.mean_12m))
ae_12m.mean_12m_prop = ae_12m.mean_12m ./ mean_total
const AE_12M = ae_12m
#println(AE_12M)

export simtree, simforest, sampleforest, simgendist, Infection, infectivitytoR
export transmissionrate, sampdegree, REGKEY, COMMUTEPROB
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

include("distribution_fitting.jl")
export fit_multi_dist, dist_fit_plot
export nll_trunc_gamma, discretize_gamma_pmf, assign_bins, nll_disc_gamma, nll_trunc_weibull, nll_trunc_normal
export erlang_truncated_means
 
include("misc_functions.jl")
export allocate_with_rounding, generation_time, severity_rolling_mean, tinf_by_age


# Use configuration file(s) to load P
#using YAML

#cfg = YAML.load_file("config/NBPMscape.yaml")
#P = NamedTuple(Symbol.(keys(cfg["P"])) .=> values(cfg["P"]))

#include("config.jl")
#using .NBPMscapeConfig


end 
