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

### Load NHS Trust catchment area data
# Source: https://app.box.com/s/qh8gzpzeo1firv1ezfxx2e6c4tgtrudl [Accessed: 2 Oct 2025]
# Also see Map: https://www.arcgis.com/home/item.html?id=3edbec15dac64c40bb6f62bccd159270#overview
# Dashboard: https://app.powerbi.com/view?r=eyJrIjoiODZmNGQ0YzItZDAwZi00MzFiLWE4NzAtMzVmNTUwMThmMTVlIiwidCI6ImVlNGUxNDk5LTRhMzUtNGIyZS1hZDQ3LTVmM2NmOWRlODY2NiIsImMiOjh9
nhs_trust_catchment_pop = XLSX.readxlsx(joinpath( @__DIR__, "..", "data", "2022 Trust Catchment Populations Worksheet.xlsx" ))
nhs_trust_catchment_pop = nhs_trust_catchment_pop["Trust Analysis!A:I"]
headers = Symbol.(nhs_trust_catchment_pop[1, :])
data = nhs_trust_catchment_pop[2:end, :]
nhs_trust_catchment_pop = DataFrame(data, headers)
# Filter for 'Emergency' admissions and Catchment year = 2020
nhs_trust_catchment_pop = filter(row -> row[:AdmissionType] == "Emergency" && row[:CatchmentYear] == 2020, nhs_trust_catchment_pop)
# Remove specialist Trusts that are unlikely to admit ARI cases to ICU or have ED attendances for ARI
specialist_trusts = ["RAN" # RAN Royal National Orthopaedic Hospital NHS Trust
                    ,"RP6" # RP6 Moorfields Eye Hospital NHS Foundation Trust
                    ,"RRJ" # RRJ The Royal Orthopaedic Hospital NHS Foundation Trust
                    ,"RL1" # RL1 The Robert Jones and Agnes Hunt Orthopaedic Hospital NHS Foundation Trust
                    ,"RBV" # RBV The Christie NHS Foundation Trust (Cancer specialist trust https://www.christie.nhs.uk/about-us/about-the-christie/a-profile-of-the-christie)
                    ,"REN" # REN The Clatterbridge Cancer Centre NHS Foundation Trust
                    ,"REP" # LIVERPOOL WOMEN'S NHS FOUNDATION TRUST - based on departments and services offered, unlikely to offer services for ARI https://www.nhs.uk/services/hospital/liverpool-women-s-nhs-foundation-trust/REP01/departments-and-services
                    ,"RET" # RET The Walton Centre NHS Foundation Trust (neuro specialist https://www.thewaltoncentre.nhs.uk/)
                    ]
nhs_trust_catchment_pop = filter(row -> !(row[:TrustCode] in specialist_trusts), nhs_trust_catchment_pop)
# Update NHS Trusts that have changed since catchment area data compiled
trust_changes_dict = Dict{Tuple{String,String}, Tuple{String,String}}()
trust_changes_dict[("RA4","Yeovil District Hospital NHS Foundation Trust")] = ("RH5","Somerset NHS Foundation Trust") # Merger in Apr 2023 https://www.somersetft.nhs.uk/?news=two-somerset-nhs-trusts-merge-to-create-unique-nhs-trust
trust_changes_dict[("RAP","North Middlesex University Hospital NHS Trust")] = ("RAL","Royal Free London NHS Foundation Trust") # Merged with Royal Free in Jan 2025 https://www.royalfree.nhs.uk/news/royal-free-london-and-north-middlesex-university-hospital-to-merge-1-january-2025
trust_changes_dict[("RBZ","Northern Devon Healthcare NHS Trust")] = ("RH8","Royal Devon University Healthcare NHS Trust") # Merger in April 2022 https://www.royaldevon.nhs.uk/news/northern-devon-healthcare-nhs-trust-and-royal-devon-and-exeter-nhs-foundation-trust-merge-to-become-the-royal-devon-university-healthcare-nhs-foundation-trust/
trust_changes_dict[("RVY","Southport And Ormskirk Hospital NHS Trust")] = ("RBN","Mersey and West Lancashire Teaching Hospitals NHS Trust") # Transferred in Jul 2023 https://www.england.nhs.uk/publication/southport-and-ormskirk-hospital-nhs-trust/
for r in eachrow(nhs_trust_catchment_pop)
    key = (r.TrustCode, r.TrustName)
    if haskey(trust_changes_dict, key)
        new1, new2 = trust_changes_dict[key]
        r.TrustCode = new1
        r.TrustName = new2
    end
end
# Compute catchment pop as % of total pop
total_pop = sum( nhs_trust_catchment_pop[ : , :Catchment])
nhs_trust_catchment_pop.catchment_prop_of_total = nhs_trust_catchment_pop.Catchment ./ total_pop
nt_catchment_pop_grouped = groupby(nhs_trust_catchment_pop[:,[:TrustCode,:Age,:catchment_prop_of_total]]
                                  #, :TrustCode) 
                                  ,[ :TrustCode, :Age])
const NHS_TRUST_CATCHMENT_POP = combine(nt_catchment_pop_grouped, :catchment_prop_of_total => sum => :catchment_prop_of_total_sum)
# Define helper function to split age group into child/adult proportions
function split_age_group(age::String, value::Float64)
    if age == "90+"  # all adult
        return (child=0.0, adult=value)
    else
        # Parse age range
        parts = split(age, "-")
        if length(parts) == 2
            start_age = parse(Int, parts[1])
            end_age = parse(Int, parts[2])
            years = end_age - start_age + 1
            # Calculate overlap
            child_years = max(0, min(17, end_age) - start_age + 1)
            adult_years = years - child_years
            child_prop = (child_years / years) * value
            adult_prop = (adult_years / years) * value
            return (child=child_prop, adult=adult_prop)
        else
            error("Unexpected age format: $age")
        end
    end
end
# Apply splitting logic to each row
df_split = transform(NHS_TRUST_CATCHMENT_POP, [:Age, :catchment_prop_of_total_sum] => ByRow((age, val) -> split_age_group(age, val)) => [:child_prop, :adult_prop])
# Aggregate by TrustCode
const NHS_TRUST_CATCHMENT_POP_ADULT_CHILD = combine(groupby(df_split, :TrustCode),
                                                    :catchment_prop_of_total_sum => sum => :catchment_prop_of_total_sum,
                                                    :child_prop => sum => :prop_child,
                                                    :adult_prop => sum => :prop_adult
                                                    )


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

#= TODO 
    Can potentially replace (some of) the functions in distribution_fitting.jl with a Julia version of
    'fitting_censored_data_function.R'. This R function accounts for sparsity of samples and right 
    censoring due to max simulation time, whereas the functions below currently only account for sparsity
    of samples (by fitting truncated data).
=#
include("distribution_fitting.jl")
export fit_multi_dist, dist_fit_plot, nll_trunc_gamma, discretize_gamma_pmf, assign_bins, nll_disc_gamma, nll_trunc_weibull, nll_trunc_normal

include("misc_functions.jl")
export allocate_with_rounding, generation_time, severity_rolling_mean, tinf_by_age



end 
