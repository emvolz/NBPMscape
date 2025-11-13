#=
Functions included here:

- sample_nhs_trust      Samples an NHS Trust code based on the ITL2 region and the age of the infectee. 
                        The function distinguishes between adult and paediatric ICU admissions, using
                        different probability distributions depending on age.

- icu_sample_prob       Calculates the probability of sampling an ICU patient for metagenomic testing based
                        on admission to particular NHS Trust, patient age, pathogen type, and site engagement stage.

- icu_sample_prob_region    Returns probability of a positive metagenomic sample based on home region, age of infectee, pathogen type and
                            stage categorisation of sites to include for testing

- sample_icu_cases      Samples cases based on a proportion of cases, either nationally or in specific regions. The number of background ARI
                        admissions is therefore not explicitly taken into account, but is implicit.

- sample_icu_cases_n    Samples cases based on a number of samples rather than just a proportion of cases. The number of samples is converted to
                        a probability using an estimate of the background number of ICU ARI cases and the modelled pathogen X ICU ARI cases

=#

# Load packages
using Revise
using NBPMscape
using Plots
using DataFrames
using DataFramesMeta

"""
Function: sample_nhs_trust(; region, infectee_age)

Description:    Samples an NHS Trust code based on the ITL2 region and the age of the infectee. 
                The function distinguishes between adult and paediatric ICU admissions, using
                different probability distributions depending on age.
                NHS Trust admission probabilities are computed using catchment area data for 
                emergency admissions in 2020 - see NBPMscape.jl for futher details on source.

Arguments:  region::String (keyword argument):      The ITL2 region code used to determine the sampling distribution for NHS Trusts.
            infectee_age::Int (keyword argument):   The age of the infectee. Determines whether to use the adult or child ICU trust distribution.

Returns:    NHS_Trust_code::String: sampled NHS Trust code from the appropriate distribution (adult or child) based
                                    on the infectee's age.

Logic:      If infectee_age > 16, the function assumes the infectee would be admitted to an adult
            ICU and samples from ITL2_TO_NHS_TRUST_PROB_ADULT. Otherwise, it samples from 
            ITL2_TO_NHS_TRUST_PROB_CHILD, assuming admission to a paediatric ICU.

Dependencies:   wsample: A weighted sampling function.
                ITL2_TO_NHS_TRUST_PROB_ADULT, ITL2_TO_NHS_TRUST_PROB_CHILD: DataFrames containing NHS 
                Trust codes and region-specific sampling probabilities.

Examples:   nhs_trust_cd_child = sample_nhs_trust( region = "TLI3", infectee_age = 5)
            nhs_trust_cd_adult = sample_nhs_trust( region = "TLI3", infectee_age = 50)
"""
function sample_nhs_trust(; region::String, infectee_age::Int)
    @assert infectee_age ≥ 0 "Age must be non-negative"
    region_sym = Symbol(region)

    trust_data = infectee_age > 16 ? ITL2_TO_NHS_TRUST_PROB_ADULT : ITL2_TO_NHS_TRUST_PROB_CHILD
    @assert region_sym in Symbol.(names(trust_data)) "Region not found in trust data"

    return wsample(trust_data.NHS_Trust_code, trust_data[:, region_sym])
end


# TODO If want to do deeper analysis using NHS Trust then could add nhs_trust to Infection struct - as at 12 Nov 2025 added after simulation during ICU sampling
# TODO Could add this code to the jump when enter ICU or maybe hospital (as long as have to be in hospital before move to ICU)

"""
Function: icu_sample_prob(p=NBPMscape.P; nhs_trust_cd, infectee_age, pathogen_type="virus", site_stage="current")

Description:    Calculates the probability of sampling an ICU patient for metagenomic testing based
                on admission to particular NHS Trust, patient age, pathogen type, and site engagement stage.

Arguments:  p::NamedTuple (default: `NBPMscape.P`)  A parameter object containing sampling
                                                    and sensitivity values.
            nhs_trust_cd::String    Code identifying the NHS Trust containing the hospital site(s) that are sampling.
            infectee_age::Int       Age of the patient being considered for sampling.
            pathogen_type::String   (default: "virus")  Type of pathogen (`"virus"`, `"bacteria"`, or `"fungi"`) is used to 
                                                        determine the metagenomic test sensitivity and thereby the probability
                                                        of a true positive result given an infection is sampled.
            site_stage::String (default: "current") Stage of site engagement ("current", "engagement", or "longlist") based 
                                                    on list of potential sampling sites for surveillance.

Returns: `Float64` The probability of sampling the patient, computed as:
                   sample_target_prob × metagenomic_test_sensitivity × proportion_of_NT_beds_at_selected_ICU

Notes:  - The proportion of NT beds is selected based on age group and site stage.
        - Metagenomic test sensitivity is pathogen-specific.

Examples:   icu_sample_prob(; p = NBPMscape.P, nhs_trust_cd = "R0A", infectee_age = 50
                            , pathogen_type = "virus", site_stage = "current")
            # Returns
            Float64
"""

function icu_sample_prob(; p = NBPMscape.P, nhs_trust_cd, infectee_age, pathogen_type="virus", site_stage="current")
    
    # Ensure pathogen_type and site_stage are in the correct format
    pathogen_type = lowercase(strip(pathogen_type))
    site_stage = lowercase(strip(site_stage))

    # Validate site_stage
    valid_stages = ["current", "engagement", "longlist"]
    site_stage ∈ valid_stages || error("Invalid site_stage: $site_stage")

    # Validate pathogen_type
    sensitivities = Dict(
        :virus => p.sensitivity_mg_virus,
        :bacteria => p.sensitivity_mg_bacteria,
        :fungi => p.sensitivity_mg_fungi
    )
    # Determine metagenomic test sensitivity based on pathogen type
    if haskey(sensitivities, Symbol(pathogen_type))
        mg_test_sensitivity = sensitivities[ Symbol(pathogen_type) ]
    else
        error("Unknown pathogen type: $pathogen_type")
    end


    # Helper function to get NHS Trust (NT) sample proportion
    function get_nt_sample_prob(nhs_trust_cd::String, infectee_age::Int, site_stage::String)
        age_key = infectee_age > 16 ? "adult" : "child"
        stage_suffix = site_stage == "current" ? "" :
                       site_stage == "engagement" ? "_E" :
                       site_stage == "longlist" ? "_E_L" : error("Invalid site_stage")

        col_name = Symbol("prop_$(age_key)_NT_beds_at_selected_ICU_C$(stage_suffix)")
        return NHS_TRUST_ICU_SAMPLE_PROB[NHS_TRUST_ICU_SAMPLE_PROB.NHS_Trust_code .== nhs_trust_cd, col_name][1]
    end
    
    # Probability of being sampled at a particular NHS Trust
    p_nt_sample = get_nt_sample_prob(nhs_trust_cd, infectee_age, site_stage)

    # Proportion of target actually sampled x test sensitivity x probability of being sampled at a particular NHS Trust
    return p.sample_target_prob * mg_test_sensitivity * p_nt_sample
end


#TODO add function description
"""
Function: icu_sample_prob_region

Description:    Returns probability of a positive metagenomic sample based on home region, age of infectee, pathogen type and
                stage categorisation of sites to include for testing

Arguments:  p::NamedTuple           (default: `NBPMscape.P`): A parameter object containing sampling and sensitivity values.
            region::String          Code identifying the ITL2 home region for infectee
            infectee_age::Int64      Age of infected individual
            pathogen_type::String   (default: "virus")  Type of pathogen (`"virus"`, `"bacteria"`, or `"fungi"`) is used to 
                                                        determine the metagenomic test sensitivity and thereby the probability
                                                        of a true positive result given an infection is sampled.
            site_stage::String      (default: "current") Stage of site engagement ("current", "engagement", or "longlist") based 
                                    on list of potential sampling sites for surveillance.

Returns:    Float64 The probability of sampling the patient, computed as:
            Proportion of target actually sampled x test sensitivity x probability of being sampled at a particular NHS Trust
            p.sample_target_prob * mg_test_sensitivity * p_nt_sample

Notes:      Metagenomic test sensitivity is pathogen-specific.

Examples:   icu_sample_prob_region(; p = NBPMscape.P, region = "TLI3", infectee_age = 50
                                   , pathogen_type="virus", site_stage="current")
"""
function icu_sample_prob_region(;p=NBPMscape.P, region::String, infectee_age::Int64
                                , pathogen_type::String="virus", site_stage::String="current")
    # Ensure pathogen_type and site_stage are in the correct format
    pathogen_type = lowercase(strip(pathogen_type))
    site_stage = lowercase(strip(site_stage))

    # Validate site_stage
    valid_stages = ["current", "engagement", "longlist"]
    site_stage ∈ valid_stages || error("Invalid site_stage: $site_stage")

    # Validate pathogen_type
    sensitivities = Dict(
        :virus => p.sensitivity_mg_virus,
        :bacteria => p.sensitivity_mg_bacteria,
        :fungi => p.sensitivity_mg_fungi
    )
    # Determine metagenomic test sensitivity based on pathogen type
    if haskey(sensitivities, Symbol(pathogen_type))
        mg_test_sensitivity = sensitivities[ Symbol(pathogen_type) ]
    else
        error("Unknown pathogen type: $pathogen_type")
    end

    # Helper function to get sample proportion for a particular region, also taking account of 
    # infectee age, and sampling site stage
    function get_region_sample_prob(region::String, infectee_age::Int64, site_stage::String)
        age_key = infectee_age > 16 ? "adult" : "child"
        stage_prefix = site_stage == "current" ? "Cur" :
                       site_stage == "engagement" ? "Cur_Eng" :
                       site_stage == "longlist" ? "Cur_Eng_Long" : error("Invalid site_stage")

        #row_name = Symbol("prop_$(age_key)_NT_beds_at_selected_ICU$(stage_suffix)")
        row_name = "$(stage_prefix)_$(age_key)"
        return ITL2_ICU_SAMPLE_PROB[ITL2_ICU_SAMPLE_PROB.site_stage_age .== row_name, region][1] # TODO ITL2_ICU_SAMPLE_PROB might need to brought into the function via an argument, e.g. itl2_icu_sample_prob = NBPMscape.ITL2_ICU_SAMPLE_PROB
    end
    
    # Probability of being sampled at a particular NHS Trust given a particular home region
    p_nt_sample = get_region_sample_prob(region, infectee_age, site_stage)
    
    # Proportion of target actually sampled x test sensitivity x probability of being sampled at a particular NHS Trust
    return p.sample_target_prob * mg_test_sensitivity * p_nt_sample
end


"""
Function:   sample_icu_cases

Description:    Samples cases based on a proportion of cases, either nationally or in specific regions. 
                The number of background ARI admissions is therefore not explicitly taken into account, but is implicit.

Arguments:  p         (default: NBPMscape.P) List of paremeters
            icu_cases::DataFrame    Version of the G dataframe output from {simforest} or {simtree} and filtered ICU 
                                    cases, i.e. rows with a value for ticu. For example:
                                    
                                    Row │ pid                                tinf     tgp       thospital  ticu      tstepdown  tdischarge  trecovered  tdeceased  severity    fatal  iscommuter  homeregion  commuteregion  generation  F      G      H         infector_age  infectee_age  importedinfection  simid                             
                                        │ String                             Float64  Float64   Float64    Float64   Float64    Float64     Float64     Float64    Symbol      Bool   Bool        String      String         Int64       Int64  Int64  Float64   Int8?         Int8          Bool               String
                                    ─────┼─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
                                    1 │ 6f992f00-b8a9-11f0-32d8-758b99f1…  60.1983  Inf         67.3632   69.9923   Inf              Inf    Inf         70.5889  verysevere   true        true  TLK4        TLK4                    4      1      0   3.91699            84            91              false  6f992f00-b8a9-11f0-32d8-758b99f1…
                                    2 │ 6f992f00-b8a9-11f0-32d8-758b99f1…  67.9981  Inf         69.429    72.2313    73.7281         Inf     94.9738   Inf       verysevere  false        true  TLK4        TLK4                    5      0      0   7.03817            83            57              false  6f992f00-b8a9-11f0-32d8-758b99f1…
                                    3 │ 6f992f00-b8a9-11f0-32d8-758b99f1…  73.9862  Inf         79.3539   80.7018   Inf              Inf    Inf         89.8647  verysevere   true        true  TLK4        TLK4                    6      0      0  17.9888             98            96              false  6f992f00-b8a9-11f0-32d8-758b99f1…
                                    4 │ 6f992f00-b8a9-11f0-32d8-758b99f1…  75.8935  Inf         84.6874   86.1361    97.9741         Inf    Inf        Inf       verysevere   true       false  TLK3        TLK3                    7      9      0   1.54777            38            97              false  6f992f00-b8a9-11f0-32d8-758b99f1…
            
            icu_sample_type::String     Either "fixed" (e.g. 0.15 of all ICU cases are sampled)
                                        or "regional" (there is a probability of sampling based on the infected individual's home region)
            icu_fixed_sample_prop::Float64  Sample proportion used when icu_sample_type = "fixed"
            pathogen_type::String   (default: "virus")  Type of pathogen (`"virus"`, `"bacteria"`, or `"fungi"`) is used to 
                                                        determine the metagenomic test sensitivity and thereby the probability
                                                        of a true positive result given an infection is sampled.
            site_stage::String      (default: "current") Stage of site engagement ("current", "engagement", or "longlist") based 
                                    on list of potential sampling sites for surveillance.
            
Returns:    A dataframe with rows sampled from the icu_cases dataframe.

Examples:   icu_cases_sub = sample_icu_cases(; icu_cases = icu_cases, icu_sample_type = "regional")
            icu_cases_sub = sample_icu_cases(; icu_cases = icu_cases, icu_sample_type = "fixed", icu_fixed_sample_prop = 0.15 )

"""

function sample_icu_cases(; p=NBPMscape.P
                          , icu_cases::DataFrame
                          , icu_sample_type::String
                          , icu_fixed_sample_prop = 0.15
                          , pathogen_type = "virus"
                          , site_stage = "current"
                          )

    # Ensure variables are in the correct format
    pathogen_type = lowercase(strip(pathogen_type))
    site_stage = lowercase(strip(site_stage))
    icu_sample_type = lowercase(strip(icu_sample_type))
    
    # Ensure only ICU cases included
    icu_cases = icu_cases[ isfinite.(icu_cases.ticu), : ] 
    if size(icu_cases,1) == 0
        error("icu_cases dataframe contains no infections (rows)")
    end
    
    # Two options for sampling method: "regional" and "fixed"
    if icu_sample_type == "regional"
        icu_cases_sampled_ix = [] # Initialise vector of indices for ICU cases to be sampled
        
        wales_regions = ["TLL3", "TLL4", "TLL5"] #TODO investigate whether can get catchment data for ITL2 in Wales 
                
        # Loop through ICU cases and determine whether they will be sampled and return a positive result
        for i in 1:nrow(icu_cases)
            if in( icu_cases[i,:homeregion], wales_regions) 
                p_sample = 0
            else
            
                # Determine probability of particular ICU case being sampled AND returning a positive result
                # (icu_sample_prob_region function includes an adjustment for the metagenomic test sensitivity)
                p_sample = icu_sample_prob_region(;p=NBPMscape.P
                                                , region = icu_cases[i,:homeregion]
                                                , infectee_age = icu_cases[i,:infectee_age]
                                                , pathogen_type = pathogen_type
                                                , site_stage = site_stage)
            end

            # If random number is less than probability of sampling then add the icu_case df row index to the list of cases to sample
            if rand() < p_sample 
                push!(icu_cases_sampled_ix, i)
            end
            
        end

        # Sample from the ICU cases
        icu_cases_sub = icu_cases[ icu_cases_sampled_ix, : ]

    elseif icu_sample_type == "fixed"
        
        # Obtain metagenomic test sensitivity for pathogen_type
        sensitivities = Dict(
                            "virus" => p.sensitivity_mg_virus,
                            "bacteria" => p.sensitivity_mg_bacteria,
                            "fungi" => p.sensitivity_mg_fungi
                            )
        mg_test_sensitivity = get(sensitivities, pathogen_type) do
            error("Unknown pathogen type: $pathogen_type")
        end
        p_positive_sample = p.sample_target_prob * mg_test_sensitivity * icu_fixed_sample_prop 
        # Sample sizes
        n_icu =  rand( Binomial( size(icu_cases,1), p_positive_sample ) )
        
        # Sample from the ICU cases
        icu_cases_sub = icu_cases[sample( 1:size(icu_cases,1), n_icu, replace=false ), :]
        
    end

    # Return df containing the sampled ICU cases
    return(icu_cases_sub)
end


"""
Function:   sample_icu_cases_n

Description:    Samples cases based on a number of samples rather than just a proportion of cases (as in {sample_icu_cases} function).
                The number of samples is converted to a probability using an estimate of the background number of ICU ARI cases
                and the modelled pathogen X ICU ARI cases.
                Summary of method:
                    - Allocate total number of weekly samples across NHS Trust sampling sites pro rata to the number of beds at each Trust / site
                    - Determine which NHS Trust each ICU case is admitted to based on probabilities
                    - Remove ICU cases not admitted to a Trust that has a sampling site
                    - Compute the proportion of the NHS Trust critical care / ICU beds that are at the sampling site
                    - Estimate background ICU ARI admissions to NHS Trusts based on total admissions and the number of beds per NHS Trust (pro rata allocation)
                    - Estimate weekly total pathogen X admissions to each NHS Trust with a sampling site
                    - Compute how many pathogen X cases are sampled at each NHS Trust sampling site with the probability
                      of sampling based on: the number of pathogen X ARI cases, the estimated number of background ARI cases, and
                      the number of samples to be taken
                    - Sample pathogen X ICU cases at each NHS Trust with a sampling site and return them in a dataframe
                    - Note that the metagenomic test sensitivity and logistical/practical probability of completing target for samples 
                      are not applied in this function
                    - It is assumed here that the number of samples targeted is achieved, which is different from when the target is 100% of 
                      cases and we assume this is not achieved for practical reasons
                    

Arguments:  p                                       (default: NBPMscape.P)
            icu_cases::DataFrame    Version of the G dataframe output from {simforest} or {simtree} and filtered for ICU 
                                    cases, i.e. rows with a value for ticu. For example:
                                    
                                    Row │ pid                                tinf     tgp       thospital  ticu      tstepdown  tdischarge  trecovered  tdeceased  severity    fatal  iscommuter  homeregion  commuteregion  generation  F      G      H         infector_age  infectee_age  importedinfection  simid                             
                                        │ String                             Float64  Float64   Float64    Float64   Float64    Float64     Float64     Float64    Symbol      Bool   Bool        String      String         Int64       Int64  Int64  Float64   Int8?         Int8          Bool               String
                                    ─────┼─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
                                    1 │ 6f992f00-b8a9-11f0-32d8-758b99f1…  60.1983  Inf         67.3632   69.9923   Inf              Inf    Inf         70.5889  verysevere   true        true  TLK4        TLK4                    4      1      0   3.91699            84            91              false  6f992f00-b8a9-11f0-32d8-758b99f1…
                                    2 │ 6f992f00-b8a9-11f0-32d8-758b99f1…  67.9981  Inf         69.429    72.2313    73.7281         Inf     94.9738   Inf       verysevere  false        true  TLK4        TLK4                    5      0      0   7.03817            83            57              false  6f992f00-b8a9-11f0-32d8-758b99f1…
                                    3 │ 6f992f00-b8a9-11f0-32d8-758b99f1…  73.9862  Inf         79.3539   80.7018   Inf              Inf    Inf         89.8647  verysevere   true        true  TLK4        TLK4                    6      0      0  17.9888             98            96              false  6f992f00-b8a9-11f0-32d8-758b99f1…
                                    4 │ 6f992f00-b8a9-11f0-32d8-758b99f1…  75.8935  Inf         84.6874   86.1361    97.9741         Inf    Inf        Inf       verysevere   true       false  TLK3        TLK3                    7      9      0   1.54777            38            97              false  6f992f00-b8a9-11f0-32d8-758b99f1…
            
            icu_ari_admissions::Int                 Estimate of weekly ICU ARI admissions (excluding pathogen X being simulated)
            icu_ari_admissions_adult_p::Float64     Estimate of proportion of ICU ARI admissions that are adults (16y and over). 
                                                    Note that this is a fixed value in the model but in reality it varies throughout the year.
            icu_ari_admissions_child_p::Float64     Estimate of proportion of ICU ARI admissions that are children (<16y). Note that this varies throughout the year.
                                                    Note that this is a fixed value in the model but in reality it varies throughout the year.
            n_icu_samples_per_week::Int             Total number of ICU samples to be taken per week
            nhs_trust_sampling_sites::DataFrame     List of NHS Trusts with sampling sites. Formatted as:
                                                    
                                                    NHS_Trust_code	ICU_site_name	                    ICU_beds_adult	ICU_beds_paediatric
                                                    AAA	            XYZ NHS Foundation Trust - site A	100	            50
 

Returns:    DataFrame of the same format as icu_cases (see above) but only contains the rows representing sampled ICU cases for 
            the simulated pathogen X

Examples:   icu_cases_sub = sample_icu_cases_n(; p = NBPMscape.P
                                               , icu_cases = icu_cases
                                               , icu_ari_admissions = 1440 # winter mean (Dec-24, Jan-25, Feb-25)
                                               , icu_ari_admissions_adult_p = 0.74
                                               , icu_ari_admissions_child_p = 0.26
                                               , n_icu_samples_per_week # Total number of ICU samples to be taken per week
                                               , nhs_trust_sites
                                                )

"""
function sample_icu_cases_n(; p = NBPMscape.P
                          , icu_cases::DataFrame # filtered from simulation of infections (G df filtered for value in ticu column)
                          , icu_ari_admissions::Int # Estimate of weekly ICU ARI admissions (excluding pathogen X being simulated)
                          , icu_ari_admissions_adult_p::Float64 # Proportion of ICU ARI admissions that are adults (16y and over)
                          , icu_ari_admissions_child_p::Float64 # Proportion of ICU ARI admissions that are children (<16y)
                          , n_icu_samples_per_week::Int # Total number of ICU samples to be taken per week
                          , nhs_trust_sampling_sites::DataFrame # List of NHS Trusts with sampling sites 
                          )

    # Load data on the number of critical care and ICU beds per NHS Trust (some filtering to remove non-ARI specialist units)
    nhs_trust_ari_cc_beds = copy(ARI_CC_BED_SITREP)
  
    ### NHS Trust catchment data is curently only for England
    wales_regions = ["TLL3", "TLL4", "TLL5"] #TODO investigate whether can get catchment data for ITL2 in Wales 
    icu_cases_Eng = filter(row -> !in( row[:homeregion], wales_regions), icu_cases)
    # Remove sites with no NHS Trust Code TODO look at how to incorporate sites without NHS Trust codes, i.e. non-NHS England sites
    nhs_trust_site_sample_targets = copy(nhs_trust_sampling_sites[:,1:5]) #copy( NHS_TRUST_SITE_SAMPLES_TARGETS )
    nhs_trust_site_sample_targets = nhs_trust_site_sample_targets[nhs_trust_site_sample_targets[:, 1] .!= "NA", :]
    ## Allocate sample target across sites (pro rata to critical care bed numbers)
    total_cc_beds_at_sampling_sites = sum( nhs_trust_site_sample_targets.ICU_beds_adult ) + sum( nhs_trust_site_sample_targets.ICU_beds_paediatric )
    nhs_trust_site_sample_targets.sample_target_per_week = Int.( round.( n_icu_samples_per_week .* ( (nhs_trust_site_sample_targets.ICU_beds_adult .+ nhs_trust_site_sample_targets.ICU_beds_paediatric) ./ total_cc_beds_at_sampling_sites) ))
    #sum(nhs_trust_site_sample_targets.sample_target_per_week)

    ### Allocate NHS Trust to each simulated ICU pathogen X case
    # Allocation is determined using probability of admission to a particular trust given home region
    # Probabilities are derived from MSOA - NHS Trust catchment area analysis
    icu_nhs_trust_cd = Vector{Any}(undef, size(icu_cases_Eng,1))
    for i in 1:size(icu_cases_Eng,1) #i=1
        hr = icu_cases_Eng[i,:].homeregion
        inf_age = icu_cases_Eng[i,:].infectee_age
        # Determine which NHS Trust individual is admitte to. Differentiated for adult and paediatric probabilities.
        if inf_age < 16
            icu_nhs_trust_cd[i] = String( only( wsample( ITL2_TO_NHS_TRUST_PROB_CHILD[:, :NHS_Trust_code], ITL2_TO_NHS_TRUST_PROB_CHILD[:, Symbol(hr)], 1) ) )
        elseif inf_age >= 16
            icu_nhs_trust_cd[i] = String( only( wsample( ITL2_TO_NHS_TRUST_PROB_ADULT[:, :NHS_Trust_code], ITL2_TO_NHS_TRUST_PROB_ADULT[:, Symbol(hr)], 1) ) )
        end
    end
    insertcols!(icu_cases_Eng, :icu_nhs_trust_cd => icu_nhs_trust_cd)
    
    # countmap(icu_nhs_trust_cd)

    ### Filter ICU cases for those that are in an NHS Trust with sampling sites (zero probability of the other ICU cases being sampled)
    # NOTE: code assumes only one sampling site per NHS Trust
    icu_cases_at_sampling_sites = filter(row -> in( row[:icu_nhs_trust_cd], nhs_trust_site_sample_targets[:,:NHS_Trust_code]), icu_cases_Eng)
    # countmap(icu_cases_at_sampling_sites.icu_nhs_trust_cd)

    if ( size(icu_cases_at_sampling_sites,1) > 0 )

        ### Compute what proportion of total critical care beds in the Trust are at the sampling site within that Trust
        # Trim NHS Trust bed data to minimal columns
        nhs_trusts_icu_ari_beds_trim = select(nhs_trust_ari_cc_beds, [:NHS_Trust_code, :Adult_critical_care_beds_available, :Paediatric_intensive_care_beds_available])
        # Merge with sample site data
        nhs_trust_site_sample_targets = leftjoin(nhs_trust_site_sample_targets, nhs_trusts_icu_ari_beds_trim, on = :NHS_Trust_code)#, kind = :left)
        # Compute proportion of Trust beds at sampling site
        nhs_trust_site_sample_targets.sampling_site_proportion_of_nt_adult = nhs_trust_site_sample_targets.ICU_beds_adult      ./ nhs_trust_site_sample_targets.Adult_critical_care_beds_available
        nhs_trust_site_sample_targets.sampling_site_proportion_of_nt_child = nhs_trust_site_sample_targets.ICU_beds_paediatric ./ nhs_trust_site_sample_targets.Paediatric_intensive_care_beds_available
        # Adjust proportions
        # If >0.9 round to 1, if > 1 round to 1, if NaN set to zero, otherwise leave as is
        for i in 1: size(nhs_trust_site_sample_targets,1) #i=14
            #println(i)
            # If missing then set to 0
            nhs_trust_site_sample_targets.sampling_site_proportion_of_nt_adult[i] = ismissing(nhs_trust_site_sample_targets.sampling_site_proportion_of_nt_adult[i]) ? 0 : nhs_trust_site_sample_targets.sampling_site_proportion_of_nt_adult[i]
            nhs_trust_site_sample_targets.sampling_site_proportion_of_nt_child[i] = ismissing(nhs_trust_site_sample_targets.sampling_site_proportion_of_nt_child[i]) ? 0 : nhs_trust_site_sample_targets.sampling_site_proportion_of_nt_child[i]
            # If >0.9 round to 1, if > 1 round to 1
            nhs_trust_site_sample_targets.sampling_site_proportion_of_nt_adult[i] = nhs_trust_site_sample_targets.sampling_site_proportion_of_nt_adult[i] > 0.9 ? 1.0 : nhs_trust_site_sample_targets.sampling_site_proportion_of_nt_adult[i]
            nhs_trust_site_sample_targets.sampling_site_proportion_of_nt_child[i] = nhs_trust_site_sample_targets.sampling_site_proportion_of_nt_child[i] > 0.9 ? 1.0 : nhs_trust_site_sample_targets.sampling_site_proportion_of_nt_child[i]
            # If NaN set to zero, otherwise leave as is
            nhs_trust_site_sample_targets.sampling_site_proportion_of_nt_adult[i] = isnan(nhs_trust_site_sample_targets.sampling_site_proportion_of_nt_adult[i])  ? 0 : nhs_trust_site_sample_targets.sampling_site_proportion_of_nt_adult[i]
            nhs_trust_site_sample_targets.sampling_site_proportion_of_nt_child[i] = isnan(nhs_trust_site_sample_targets.sampling_site_proportion_of_nt_child[i]) ? 0 : nhs_trust_site_sample_targets.sampling_site_proportion_of_nt_child[i]
        end
        #println(nhs_trust_site_sample_targets)
        #CSV.write("nhs_trust_site_sample_targets.csv",nhs_trust_site_sample_targets)
        
        ## Estimate weekly ARI ICU admissions per NHS Trust    

        # Create dfs to store info
        est_weekly_icu_ari_admissions_adult = zeros(Int, size(nhs_trust_ari_cc_beds,1)) #Vector{Int}(undef, size(nhs_trust_ari_cc_beds,1)) 
        est_weekly_icu_ari_admissions_child = copy(est_weekly_icu_ari_admissions_adult)
        
        # Estimate adult and child ICU ARI admissions per week
        icu_ari_admissions_adult = Int( round( icu_ari_admissions * icu_ari_admissions_adult_p ) )
        icu_ari_admissions_child = Int( round( icu_ari_admissions * icu_ari_admissions_child_p ) )
            
        for i in 1:length(nhs_trust_ari_cc_beds[:,:Adult_critical_care_beds_available]) #i=1
            nt = nhs_trust_ari_cc_beds[i, :NHS_Trust_code ]
            # Allocate background ICU ARI admissions to NHS Trusts based on number of beds per NHS Trust
            # (use probability based estimate)
            #est_weekly_icu_ari_admissions_adult[i] = rand( Binomial( icu_ari_admissions_adult, nhs_trust_ari_cc_beds.p_adult[i] ) )
            #est_weekly_icu_ari_admissions_child[i] = rand( Binomial( icu_ari_admissions_child, nhs_trust_ari_cc_beds.p_child[i] ) )
            # OR deterministic
            est_weekly_icu_ari_admissions_adult[i] = Int( round( icu_ari_admissions_adult * nhs_trust_ari_cc_beds.p_adult[i] ) )
            est_weekly_icu_ari_admissions_child[i] = Int( round( icu_ari_admissions_child * nhs_trust_ari_cc_beds.p_child[i] ) )
        end
        # check # sum(est_weekly_icu_ari_admissions_adult) # sum(est_weekly_icu_ari_admissions_child)
        # Add columns to df
        insertcols!(nhs_trust_ari_cc_beds, :est_weekly_icu_ari_admissions_adult => est_weekly_icu_ari_admissions_adult)
        insertcols!(nhs_trust_ari_cc_beds, :est_weekly_icu_ari_admissions_child => est_weekly_icu_ari_admissions_child)
        #CSV.write("nhs_trusts_icu_ari_beds.csv",nhs_trusts_icu_ari_beds)
        
        # Check that totals for ICU ARI admissions approximately match - particularly if using
        # TODO could allocate or remove any difference randomly across the NHS Trusts
        #icu_ari_admissions == sum( nhs_trust_ari_cc_beds[:,:est_weekly_icu_ari_admissions_adult] ) + sum( nhs_trust_ari_cc_beds[:,:est_weekly_icu_ari_admissions_child])
        
        ## Determine weekly totals for pathogen X ICU ARI admissions per NHS Trust (specifically, NHS Trusts that have a site taking sampling)
        # Number of weeks in simulation
        
        # Add to the ICU cases df the sim week corresponding to each time of admission to ICU (ticu)
        icu_cases_at_sampling_sites.ticu_wk = Int.( ceil.( icu_cases_at_sampling_sites.ticu ./ 7 ) ) #Int( ceil( 53 / 7 ) ) # print(icu_cases_at_sampling_sites)
        # Temporary df containing info on NHS Trust, week, and infectee age for all ICU cases at sampling sites
        #icu_cases_nt_wk_age = DataFrame( nt_cd   = icu_cases_at_sampling_sites.icu_nhs_trust_cd
        #                               , ticu_wk = icu_cases_at_sampling_sites.ticu_wk
        #                               , inf_age = icu_cases_at_sampling_sites.infectee_age
        #                               )
        #println(icu_cases_nt_wk_age)

        # Simulation weeks with ICU cases at sampling sites
        #if ( size(icu_cases_at_sampling_sites.ticu_wk, 1 ) > 0 )
            icu_weeks = minimum(icu_cases_at_sampling_sites.ticu_wk):1:maximum(icu_cases_at_sampling_sites.ticu_wk) |> collect
        #else
        #    icu_weeks = missing
        #end
        
        # Build df to store pathogen X ICU cases
        nt_weekly_icu_cases_px_adult = DataFrame( sampling_NHS_Trust = nhs_trust_site_sample_targets[:,:NHS_Trust_code] )
        
        # Add a column for each simulation week
        for i in icu_weeks
            insertcols!(nt_weekly_icu_cases_px_adult, Symbol("week_$(Int(i))") => fill(0, nrow(nt_weekly_icu_cases_px_adult))) 
        end
            
        # Create df to store all types of weekly cases
        nt_weekly_icu_cases_px_child  = copy(nt_weekly_icu_cases_px_adult) # check sum of values in df nt_weekly_icu_cases_px_adult # println(combine(nt_weekly_icu_cases_px_adult, names(nt_weekly_icu_cases_px_adult, Number) .=> sum))
        nt_weekly_icu_cases_bkg_adult = copy(nt_weekly_icu_cases_px_adult)
        nt_weekly_icu_cases_bkg_child = copy(nt_weekly_icu_cases_px_adult)
        nt_weekly_icu_cases_all_adult = copy(nt_weekly_icu_cases_px_adult)
        nt_weekly_icu_cases_all_child = copy(nt_weekly_icu_cases_px_adult)

        # Sum pathogen X ICU cases by week and by NHS Trust
        
        # Base for icu_week column index in df
        w0 = minimum(icu_weeks) # also accounts for column 1 being the NHS Trust code
        # Fill dataframes by week (w) and NHS Trust (nt)
        for w in icu_weeks # w=12 ; nt="RGT"
            for nt in nt_weekly_icu_cases_px_adult[:,:sampling_NHS_Trust] 
                #println([w," ",nt])
                
                ## Convert nt into df row index
                nt_ix = findfirst(x -> x == nt, nt_weekly_icu_cases_px_adult.sampling_NHS_Trust) 

                ## Convert icu_week to column index
                w_ix = w - w0 +2 # +2 because NHS Trust codes are in column 1

                ## Fill weekly pathogen X ICU cases by NHS Trust for:
                # Adults
                nt_ix = findfirst(x -> x == nt, nt_weekly_icu_cases_px_adult.sampling_NHS_Trust) 
                nt_wk_age_infections_adult = filter(row -> row.icu_nhs_trust_cd == nt && row.ticu_wk == w && row.infectee_age >= 16, icu_cases_at_sampling_sites)
                nt_weekly_icu_cases_px_adult[ nt_ix, w_ix ] = size( nt_wk_age_infections_adult, 1 )
                # Children
                nt_wk_age_infections_child = filter(row -> row.icu_nhs_trust_cd == nt && row.ticu_wk == w && row.infectee_age < 16, icu_cases_at_sampling_sites)
                nt_weekly_icu_cases_px_child[ nt_ix, w_ix ] = size( nt_wk_age_infections_child , 1 )
                
                ## Background (non-simulation) weekly ICU ARI cases per NHS Trust
                # Includes adjustment factor because sampling sites potentially only represent part of an NHS Trust's total ICU beds
                # Adult
                nt_weekly_adult = only( filter(row -> row[:NHS_Trust_code] == nt, nhs_trust_ari_cc_beds).est_weekly_icu_ari_admissions_adult )
                sampling_site_prop_nt_adult = only( filter(row -> row[:NHS_Trust_code] == nt, nhs_trust_site_sample_targets).sampling_site_proportion_of_nt_adult )
                nt_weekly_icu_cases_bkg_adult[ nt_ix, w_ix] = Int( round( nt_weekly_adult * sampling_site_prop_nt_adult ) )
                # Child
                nt_weekly_child = only( filter(row -> row[:NHS_Trust_code] == nt, nhs_trust_ari_cc_beds).est_weekly_icu_ari_admissions_child )
                sampling_site_prop_nt_child = only( filter(row -> row[:NHS_Trust_code] == nt, nhs_trust_site_sample_targets).sampling_site_proportion_of_nt_child )
                nt_weekly_icu_cases_bkg_child[ nt_ix, w_ix] = Int( round( nt_weekly_child * sampling_site_prop_nt_child ) )
                
                ## Total weekly ICU ARI cases per NHS Trust
                nt_weekly_icu_cases_all_adult[ nt_ix, w_ix] = nt_weekly_icu_cases_px_adult[ nt_ix, w_ix] + nt_weekly_icu_cases_bkg_adult[ nt_ix, w_ix]
                nt_weekly_icu_cases_all_child[ nt_ix, w_ix] = nt_weekly_icu_cases_px_child[ nt_ix, w_ix] + nt_weekly_icu_cases_bkg_child[ nt_ix, w_ix]
            end
        end
        
        # Check
        # nt_weekly_icu_cases_px_adult
        # nt_weekly_icu_cases_px_child
        # nt_weekly_icu_cases_bkg_adult
        # nt_weekly_icu_cases_bkg_child
        # nt_weekly_icu_cases_all_adult
        # nt_weekly_icu_cases_all_child

        ## Determine number of pathogen X ICU cases that are sampled each week in each NHS Trust sampling site
        # Create dfs to store numbers of pathogen X samples and start at zero
        n_icu_px_samples_adult = copy( nt_weekly_icu_cases_px_adult )
        n_icu_px_samples_child = copy( nt_weekly_icu_cases_px_child )
        n_icu_px_samples_adult[:,2:end] .= 0 # Make sure all elements are zero apart from the NHS Trust codes
        n_icu_px_samples_child[:,2:end] .= 0

        # Create df to store sampled case information
        icu_cases_sub = empty( icu_cases_at_sampling_sites )

        ## Determine which ICU cases get sampled from the total (background non-pathogen X plus simulated pathogen X)
        # Loop through sim weeks
        for w in icu_weeks #TEST w=10 ; nt="RGT"
            # Loop through sampling NHS Trusts
            for nt in nt_weekly_icu_cases_all_adult[:,:sampling_NHS_Trust]
                #println(w,"_",nt)
                nt_ix = findfirst(x -> x == nt, nt_weekly_icu_cases_px_adult.sampling_NHS_Trust) 
                
                ## Convert icu_week to column index
                w_ix = w - w0 +2 # +2 because NHS Trust codes are in column 1

                # Total samples per NHS Trust (site)
                sampling_nt_info = @subset( nhs_trust_site_sample_targets, :NHS_Trust_code .== nt )
                nt_site_weekly_samples = only( sampling_nt_info.sample_target_per_week )

                # Total samples per NHS Trust (site)
                nt_site_weekly_cases_total = nt_weekly_icu_cases_all_adult[ nt_ix, w_ix] + nt_weekly_icu_cases_all_child[ nt_ix, w_ix]
                
                # If the number of samples is less than the number of cases, then need to determine which cases are sampled
                # Adult
                cases_px_adult  = nt_weekly_icu_cases_px_adult[ nt_ix, w_ix ]
                cases_bkg_adult = nt_weekly_icu_cases_bkg_adult[ nt_ix, w_ix ]
                # Child
                cases_px_child  = nt_weekly_icu_cases_px_child[ nt_ix, w_ix ]
                cases_bkg_child = nt_weekly_icu_cases_bkg_child[ nt_ix, w_ix ]
                # (1) split samples between adult and child ICU cases pro rata
                if nt_site_weekly_samples < nt_site_weekly_cases_total 
                    nt_site_weekly_samples_adult = nt_site_weekly_cases_total == 0 ? 0 : Int( round( nt_site_weekly_samples * ( nt_weekly_icu_cases_all_adult[ nt_ix, w_ix] / nt_site_weekly_cases_total ) ) )
                    nt_site_weekly_samples_child = nt_site_weekly_cases_total == 0 ? 0 : Int( round( nt_site_weekly_samples * ( nt_weekly_icu_cases_all_child[ nt_ix, w_ix] / nt_site_weekly_cases_total ) ) )
                    # Determine number of pathogen X cases that are sampled in each NHS Trust (site) and each sim week
                    # Adult
                    n_icu_px_samples_adult[nt_ix,w_ix] = rand( Hypergeometric( cases_px_adult , cases_bkg_adult , nt_site_weekly_samples_adult ) )
                    # Child
                    n_icu_px_samples_child[nt_ix,w_ix] = rand( Hypergeometric( cases_px_child , cases_bkg_child , nt_site_weekly_samples_child ) )
                # (2) Otherwise if total samples > total cases then all cases are sampled
                elseif nt_site_weekly_samples >= nt_site_weekly_cases_total 
                    n_icu_px_samples_adult[nt_ix,w_ix] = cases_px_adult
                    n_icu_px_samples_child[nt_ix,w_ix] = cases_px_child
                end

                # Determine which pathogen X cases are sampled
                # Filter ICU case dfs for matching ICU week, infectee age and NHS Trust
                #icu_cases_temp_adult = icu_cases_at_sampling_sites[ icu_cases_at_sampling_sites.ticu_wk == w && icu_cases_at_sampling_sites.icu_nhs_trust_cd == nt && icu_cases_at_sampling_sites.infectee_age >= 16, :]
                icu_cases_temp_adult = @subset(icu_cases_at_sampling_sites, :ticu_wk .== w, :icu_nhs_trust_cd .== nt, :infectee_age .>= 16 ) #println(icu_cases_temp_adult)
                #icu_cases_temp_child = icu_cases_at_sampling_sites[ icu_cases_at_sampling_sites.ticu_wk == w && icu_cases_at_sampling_sites.icu_nhs_trust_cd == nt && icu_cases_at_sampling_sites.infectee_age <  16, :]
                icu_cases_temp_child = @subset(icu_cases_at_sampling_sites, :ticu_wk .== w, :icu_nhs_trust_cd .== nt, :infectee_age .< 16 )

                # Adult
                icu_cases_temp_sub_adult = icu_cases_temp_adult[ sample( 1:size(icu_cases_temp_adult,1), n_icu_px_samples_adult[nt_ix,w_ix], replace=false ), :]
                # Child
                icu_cases_temp_sub_child = icu_cases_temp_child[ sample( 1:size(icu_cases_temp_child,1), n_icu_px_samples_child[nt_ix,w_ix], replace=false ), :]
                
                # Add to df containing sampled pathogen X ICU cases
                append!( icu_cases_sub , icu_cases_temp_sub_adult )
                append!( icu_cases_sub , icu_cases_temp_sub_child )

            end
        end

        # Check
        # n_icu_px_samples_adult
        # n_icu_px_samples_child
        # println(icu_cases_sub)
        # println(icu_cases)
        # minimum( icu_cases.ticu )
        # minimum( icu_cases_sub.ticu )
    else
        # return an empty dataframe if the size of icu_cases_at_sampling_sites was determined to be = 0 as per the start of this if statement
        icu_cases_sub = empty( icu_cases_at_sampling_sites )

    end

    return( icu_cases_sub )
    # NOTE that the metagenomic test sensitivity has not been applied in this function - This enables the number of samples
    # and the number of positive samples to be recorded separately.
    # The logistical/practical probability of completing target for samples is not applied in this function either.
    # We assume here that the number of samples targeted is achieved. This is different from when the target is 100% of 
    # ICU cases, which is believed to be practically unrealistic to achieve and so an adjustment must be made.
end