# Load packages
#using Revise
#using NBPMscape
#using GLM, Statistics, StatsBase, Distributions
#using DataFrames, CSV 
#using DataFramesMeta
#using Plots, StatsPlots
#using JLD2

#=  This file contains a number of functions used in the analysis of time to detection 
    under different scenarios. The functions are:
    - icu_td             - computes times to detection for simulated outbreaks with ICU sampling 
    - gp_td              - computes times to detection for simulated outbreaks with GP sampling 
    - secondary_care_td  - computes times to detection for simulated outbreaks with sampling from hospitals in the HARISS network
    - icu_v_pc_td        - used in previous analyses but superseded by icu_td and gp_td - computes the times to detection for simulated outbreaks
    - analyse_td_columns - computes median values across simulation replicates and % of sim reps that have a (time to) detection
    - clean_td_results   - used in previous analyses but superseded - removes some rows (simulation replicates) and converts type of Inf from string
    - plot_hist          - plots a basic histogram for results of a chosen column of data after performing some data cleaning tasks 
                          (similar to those in {clean_TD_results})
=#

"""
Function: icu_td

Description:    Computes times to detection based on metagenomic sampling in Intensive Care Units (ICUs).

Arguments:  p                           Set of parameters relating to the outbreak and sampling
            sims                        Formatted as a vector of dataframes. Simulated outbreak data in the G dataframe output from simtree or simforest functions
            icu_sample_type::String     Two options: "regional" or "fixed". The latter assumes a fixed proportion of ICU admissions are sampled
                                        at all ICUs nationally. The proportion is defined by p_icu. 
                                        The "regional" option uses proportions that are specific to ITL2 regions based on the sites selected for
                                        sampling (also see argument 'sample_icu_cases_version'). 
                                        Note that with the "fixed" option, no account is taken of test sensitivity or practical sampling proportion,
                                        which "regional" does by using the {sample_icu_cases} or {sample_icu_cases_n} functions.
            pathogen_type::String       "virus", "bacteria" or "fungi". Used to determine the metagenomic test sensitivity to apply and thereby
                                        change the probability of a positive sample.
            site_stage::String          "current", "engagement" or "longlist". Refers to the list of selected sites and the stage of onboarding 
                                        for sampling.
            sample_icu_cases_version::String    "number" or "proportion". This choice determines the sampling function used, either sample_icu_cases_n or 
                                                sample_icu_cases respectively. The difference is that the former samples ICU cases using a specific number 
                                                of samples per week and requires sampling sites to be specified, while the latter samples a proportion of
                                                all ICU cases. 
            n_icu_samples_per_week::Int     Total number of metagenomic samples to be taken per week across all sites
            only_sample_before_death::Bool  true or false. Modelled sampling period is a uniform distribution between time of admission to ICU (ticu) 
                                            and a variable number days after admission which is set in p.icu_swab_lag_max. It is possible for the time
                                            of death (tdeceased) to be before the sample time. Setting only_sample_before_death = true constrains the
                                            upper limit on tsample to be equal to tdeceased.
            p_icu::Float                ICU sampling proportion. Only relevant for sample_icu_cases_version = "proportion".
            icu_ari_admissions::Int     Estimate of weekly ICU ARI admissions (excluding pathogen X being simulated), e.g. 793 mean in summer 2024 and
                                        1440 in winter 2024/25 in England
            icu_turnaround_time         Formatted as a vector of length 2, giving the lower and upper limits of the time to process a sample and report
                                        results / declare detection, e.g. [2,4]
            icu_ari_admissions_adult_p::Float64     Proportion of ICU ARI admissions that are adults (16y and over). Note that this varies seasonally 
                                                    so is an estimate.
            icu_ari_admissions_child_p::Float64     Proportion of ICU ARI admissions that are children (<16y). Note that this varies seasonally 
                                                    so is an estimate.
            nhs_trust_sampling_sites::DataFrame     List of NHS Trusts containing the sampling sites including code, name and adult and paediatric
                                                    critical care/ICU bed numbers. Format example below:
                                                    
                                                    NHS_Trust_code	ICU_site_name	            ICU_beds_adult	ICU_beds_paediatric
                                                    AAA	            XYZ NHS Foundation Trust	43	            0                  

Returns:    Dataframe containing times to detection (TD and 3TD) for each simulation replicate as well as metadata (sim rep number, simid, the time of infection
            relating to the TD and the last time of infection relating to the time of detection of the 3rd case (3TD). See example df in 'Examples' below.

Examples:    icu_tds = NBPMscape.icu_td(;  p=NBPMscape.P
                                        , sims = sims
                                        , icu_sample_type = "regional" # "regional" or "fixed". If "fixed" then p_icu will be used, note that this doesn't take into account test sensitivity or practical sampling proportion, which "regional" does by using {sample_icu_cases} function
                                        , pathogen_type = "virus"
                                        , site_stage = "current" 
                                        , p_icu = 0.15 # ICU sampling proportion
                                        #, icu_ari_admissions = 793 # 1440 # Weekly ICU admission numbers [summer,winter]
                                        , icu_turnaround_time = [2,4] # Time to process sample and report results / declare detection
                                        )
            # Returns
            1200×8 DataFrame
            Row │ sim_n  ICU_TD    ICU_3TD   n_ICU_cases  n_ICU_cases_sampled  ICU_simid                          tinf_relating_to_ICU_TD  last_tinf_relating_to_ICU_3TD 
                │ Int64  Float64   Float64   Int64        Int64                String                             Float64                  Float64
            ──────┼────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
                1 │     1  Inf       Inf                14                    0  ["6f992f00-b8a9-11f0-32d8-758b99…                Inf                             Inf
                2 │     2   82.8091   93.0427           60                    9  ["71ea31dc-b8a9-11f0-22a2-51bf2b…                 59.7301                         83.2172
                3 │     3   68.4958   80.1025         2473                   46  ["7959a25e-b8a9-11f0-386b-77359c…                 58.8563                         61.3809
                4 │     4  102.036   Inf                24                    2  ["5346f5de-b8aa-11f0-2a4d-e1fc8c…                 88.9252                        Inf
                5 │     5   38.5468   67.2674         7182                  786  ["55ea4ffe-b8aa-11f0-35cc-3b70d9…                 29.4317                         57.3555
                6 │     6   68.0169   77.2858          526                  123  ["bb622798-b8ac-11f0-3a44-b9eaff…                 58.8111                         65.3002
                7 │     7   51.6203   54.1451         9970                 1987  ["edcffff4-b8ac-11f0-1ce3-9ba211…                 43.0713                         48.5687
                8 │     8   92.0159   95.4548           19                    4  ["99674cf2-b8b0-11f0-0ba2-bf1638…                 82.0763                         84.1551
            ⋮   │   ⋮       ⋮         ⋮           ⋮                ⋮                           ⋮                             ⋮                           ⋮
            1194 │  1194   81.6484   89.1694           79                   21  ["5ba90d2c-b8f3-11f0-25cb-4d44ad…                 72.5058                         77.4501
            1195 │  1195   92.2708  Inf                10                    1  ["61b38152-b8f3-11f0-01e0-39d3bc…                 79.489                         Inf
            1196 │  1196   73.2853   81.3157         2110                  102  ["63626752-b8f3-11f0-2101-4f08da…                 64.6652                         58.7474
            1197 │  1197   55.6959   63.9925         2746                   80  ["1b615b74-b8f4-11f0-24c7-d3c01a…                 48.2511                         55.8263
            1198 │  1198   65.8552   92.8191          800                   22  ["10ad8cba-b8f5-11f0-019c-b78e40…                 49.5414                         84.2504
            1199 │  1199   57.7104   69.4578         6380                  390  ["536f0c54-b8f5-11f0-2c2f-e9af04…                 46.6834                         54.46
            1200 │  1200   68.9914  101.335            30                    3  ["7579ea10-b8f7-11f0-2219-3dd9d2…                 53.4871                         89.8118
            

"""
function icu_td(; p = NBPMscape.P
                , sims
                # ICU parameters
                , icu_sample_type = NBPMscape.P.icu_sample_type #"regional" # "regional" or "fixed". If "fixed" then p_icu will be used, note that this doesn't take into account test sensitivity or practical sampling proportion, which "regional" does by using {sample_icu_cases} function
                , pathogen_type = NBPMscape.P.pathogen_type # "virus"
                , site_stage = NBPMscape.P.icu_site_stage #"current" 
                , sample_icu_cases_version = NBPMscape.P.sample_icu_cases_version #"number" # or "proportion"
                , n_icu_samples_per_week = NBPMscape.P.n_icu_samples_per_week #300
                , only_sample_before_death = NBPMscape.P.icu_only_sample_before_death #true
                , p_icu = NBPMscape.P.p_sampled_icu #0.15 # ICU sampling proportion
                , icu_turnaround_time = NBPMscape.P.turnaroundtime_icu #[2,4] # Time to process sample and report results / declare detection
                , icu_ari_admissions::Int = NBPMscape.P.icu_ari_admissions # 1440 # Estimate of weekly ICU ARI admissions (excluding pathogen X being simulated)
                , icu_ari_admissions_adult_p::Float64 = NBPMscape.P.icu_ari_admissions_adult_p #0.76 # Proportion of ICU ARI admissions that are adults (16y and over)
                , icu_ari_admissions_child_p::Float64 = NBPMscape.P.icu_ari_admissions_child_p #0.24 # Proportion of ICU ARI admissions that are children (<16y)
                , nhs_trust_sampling_sites::DataFrame = NBPMscape.P.icu_nhs_trust_sampling_sites #DataFrame()
                )

    n_replicates = length(sims)
    
    # Create df to store results
    col_any = Vector{Any}(undef, n_replicates)
    fill!(col_any, Inf)
    sim_tds_cols = [copy(col_any) for _ in 1:9]
    sim_tds = DataFrame( sim_tds_cols, ["sim_n"
                                        # Time to detection for 1 case
                                        ,"ICU_TD"     # ICU sampling only
                                        # Time to detection for 3 cases
                                        ,"ICU_3TD" # ICU sampling only
                                        # Number of cases and samples
                                        ,"n_ICU_cases","n_ICU_cases_sampled","n_ICU_cases_sampled_positive"
                                        ,"ICU_simid"
                                        , "tinf_relating_to_ICU_TD" # useful for investigating impact of early termination of simtree within simforest
                                        , "last_tinf_relating_to_ICU_3TD" # useful for investigating impact of early termination of simtree within simforest
                                        ])

    # Loop through simulation replicates, sampling infection cases (adapted from sampleforest() function in 'core.jl' and 'median_TD_by_region.jl'),
    # and computing times to detection
    for s in 1:n_replicates
        #TEST
        println("$(s) of $(n_replicates) replicates")

        # Add simulation number to results df
        #TEST s=1
        sim_tds[s,:sim_n] = s
        fo = sims[s]
                   
        # Filter for ICU cases
        icu_cases = size(fo,1) == 0 ? DataFrame() : fo[ isfinite.(fo.ticu), : ]
        
        # Add simids
        if size(icu_cases,1) > 0
            simids_icu = unique( icu_cases.simid )  #unique( [icu_cases[1,:simid]] ) 
            sim_tds[s,:ICU_simid] = simids_icu #unique( [icu_cases[1,:simid] ]) 
        end
    
        ### Record infection time stats
        if sample_icu_cases_version =="proportion"
            ## Sample from all ICU cases
            icu_cases_sub = sample_icu_cases(; icu_cases = icu_cases
                                            , icu_sample_type = icu_sample_type
                                            , icu_fixed_sample_prop = p_icu
                                            , pathogen_type = pathogen_type #"virus"
                                            , site_stage = site_stage )#"current"  )
        elseif sample_icu_cases_version == "number"
            # Note that metagenomic test sensitivity is applied in sample_icu_cases but not in sample_icu_cases_n
            icu_cases_sub = sample_icu_cases_n(; p=NBPMscape.P
                                            , icu_cases = icu_cases
                                            #, pathogen_type::String = "virus" # 
                                            #, nhs_trust_site_sample_targets # In df format same as NHS_TRUST_SITE_SAMPLES_TARGETS. List of NHS Trusts with sampling sites with number of samples per week
                                                                              # NHS Trusts (code, name) containing sites sampling. Includes Trust, site,
                                                                              # beds (adult, paediatric), beds as proportion of NHS Trust beds (must match nhs_trust_cc_beds) and weekly samples
                                            #, nhs_trust_ari_cc_beds::DataFrame # List of critical care beds for each NHS Trusts
                                            , icu_ari_admissions = icu_ari_admissions # Estimate of weekly ICU ARI admissions (excluding pathogen X being simulated)
                                            , icu_ari_admissions_adult_p = icu_ari_admissions_adult_p # Proportion of ICU ARI admissions that are adults (16y and over)
                                            , icu_ari_admissions_child_p = icu_ari_admissions_child_p # Proportion of ICU ARI admissions that are children (<16y)
                                            , nhs_trust_sampling_sites = nhs_trust_sampling_sites
                                            , n_icu_samples_per_week = n_icu_samples_per_week
                                            )
            
            # Record number of cases in the sim rep and the number that were tested AND (below) returned positive results
            sim_tds[s,:n_ICU_cases] = size(icu_cases,1)
            sim_tds[s,:n_ICU_cases_sampled] = size(icu_cases_sub,1)

            # Obtain metagenomic test sensitivity for pathogen_type
            # Additional sub-sampling to account for metagenomic test sensitivity
            sensitivities = Dict(
                            "virus" => p.sensitivity_mg_virus,
                            "bacteria" => p.sensitivity_mg_bacteria,
                            "fungi" => p.sensitivity_mg_fungi
                            )
            mg_test_sensitivity = get(sensitivities, pathogen_type) do
                error("Unknown pathogen type: $pathogen_type")
            end
            # Positive sample sizes
            n_icu =  rand( Binomial( size(icu_cases_sub,1), mg_test_sensitivity ) )
            # Subsample of ICU cases
            #icu_cases_sub = icu_cases[sample( 1:size(icu_cases,1), n_icu, replace=false ), :]
            icu_cases_sub = icu_cases_sub[sample( 1:size(icu_cases_sub,1), n_icu, replace=false ), :]
            
            # Record number of true positives returned from metagenomic testing
            sim_tds[s,:n_ICU_cases_sampled_positive] = size(icu_cases_sub,1)

        end
        
        if size(icu_cases_sub,1) == 0
                
            # Record zero cases and zero sampled
            #sim_tds[s,:n_ICU_cases] = 0
            #sim_tds[s,:n_ICU_cases_sampled] = 0

            # Record Inf as time to detection for earliest:
            # - 1 ICU case
            sim_tds[s,:ICU_TD] = Inf
            # - 3 ICU cases
            sim_tds[s,:ICU_3TD] = Inf
            # and for tinf
            sim_tds[s,:tinf_relating_to_ICU_TD] = Inf
            sim_tds[s,:last_tinf_relating_to_ICU_3TD] = Inf

            # Define empty vector for top 3 times to detection - to be used later
            icu_cases_sub_top3_td = []

        elseif size(icu_cases_sub,1) > 0 
                
            # Generate sample times
            if only_sample_before_death == true
                # Set upper limit on sample time to the minimum of time of death and ICU admission plus the user defined p.icu_swab_lag_max
                icu_tsample = map( g -> rand( Uniform( g.ticu[1], min( g.tdeceased[1], g.ticu[1] + p.icu_swab_lag_max ) ) ) , eachrow(icu_cases_sub) )
            else
                # Set upper limit on sample time to ICU admission plus the user defined p.icu_swab_lag_max
                icu_tsample = map( g -> rand( Uniform( g.ticu[1], g.ticu[1] + p.icu_swab_lag_max)) , eachrow(icu_cases_sub) )
            end
            icu_cases_sub.tsample = icu_tsample 

            # Simulate reports times
            icu_treport = (icu_tsample .+ rand(Uniform(icu_turnaround_time[1], icu_turnaround_time[2])) )
            icu_cases_sub.treport = icu_treport
            
            # Find cases with 3 shortest times to detection (TD)
            icu_cases_sub_top3_td = first(sort(icu_cases_sub, :treport), 3) #println(icu_cases_sub_top3_td.treport)

            ## Add to results df for ICU times to detection
            # Record the times to detection for earliest:
            # - 1 ICU case
            sim_tds[s,:ICU_TD] = icu_cases_sub_top3_td[1,:treport]
            # and for tinf
            sim_tds[s,:tinf_relating_to_ICU_TD] = icu_cases_sub_top3_td[1,:tinf]
            
            # - 3 ICU cases
            if size(icu_cases_sub_top3_td,1) >= 3
                sim_tds[s,:ICU_3TD] = icu_cases_sub_top3_td[3,:treport]
                # and for tinf
                sim_tds[s,:last_tinf_relating_to_ICU_3TD] = icu_cases_sub_top3_td[3,:tinf]
            else
                sim_tds[s,:ICU_3TD] = Inf
                sim_tds[s,:last_tinf_relating_to_ICU_3TD] = Inf
            end
        end
        
    end # End of loop through different simulations
    return sim_tds
end # End of icu_td function



"""
Function: gp_td
Description:    Computes times to detection based on metagenomic sampling in primary care setting, i.e. General Practioner (GP) appointments modelled 
                parameters for the existing Oxford-RCGP RSC primary care surveillance.

Arguments:  p       Set of parameters relating to the outbreak and sampling
            sims    Formatted as a vector of dataframes. Simulated outbreak data in the G dataframe output from simtree or simforest functions
            # Parameters for existing Oxford-RCGP RSC primary care surveillance
            gp_practices_total      Total number of GP practices in England. 6199 at July 2025. Source: BMA analysis (https://www.bma.org.uk/advice-and-support/nhs-delivery-and-workforce/pressures/pressures-in-general-practice-data-analysis) of NHS Digital General Practice Workforce Statistics (https://digital.nhs.uk/data-and-information/publications/statistical/general-and-personal-medical-services) [Accessed 2 Sep 2025]  
            gp_practices_swab       Number of GP practices taking swabs for virology surveillance. 300 Source: Data quality report: national flu and COVID-19 surveillance report (27 May 2025)
            gp_swabs_mg             Assumed number of swabs that are taken for metagenomic sequencing
            pop_eng                 Population of England. 5.7106398e7 used in outbreak model. This is lower than current population but is disaggregated by various parameters. 
                                    Source: ONS mid-year 2022 UK population data disaggregated by various geo levels and age groups and gender. [Accessed 6 November 2024] Available at https://www.ons.gov.uk/peoplepopulationandcommunity/populationandmigration/populationestimates/datasets/populationestimatesforukenglandandwalesscotlandandnorthernireland 
            gp_ari_consults         Number of ARI consultations per 100k of England population per week. mean summer 2024 = 180 and mean winter 2024/25 = 327. 
                                    Source: Analysis of data extracted from RCGP Research & Surveillance Centre (RSC) Virology Dashboard [Accessed 29 Aug 2025]
            gp_ari_swabs            Number of swabs taken from suspected ARI per week. Mean summer 2024 = 319 and mean winter 2024/25 = 747. 
                                    Source: Analysis of data extracted from RCGP Research & Surveillance Centre (RSC) Virology Dashboard [Accessed 29 Aug 2025]
            pc_swab_turnaround_time     Formatted as a vector of length 2, giving the lower and upper limits of the time to process a sample and report
                                        results / declare detection, e.g. [2,4]. Source: Data quality report: national flu and COVID-19 surveillance report (27 May 2025).
            pathogen_type::String   "virus", "bacteria" or "fungi". Used to determine the metagenomic test sensitivity to apply and thereby
                                    change the probability of a positive sample.

Returns:    Dataframe containing times to detection (TD and 3TD) for each simulation replicate as well as metadata (sim rep number, simid, the time of infection
            relating to the TD and the last time of infection relating to the time of detection of the 3rd case (3TD). See example df in 'Examples' below.

Examples:   gp_tds = NBPMscape.gp_td(; p = NBPMscape.P
                                     , sims = sims 
                                     # Parameters for existing Oxford-RCGP RSC primary care surveillance
                                     , gp_practices_total = 6199 
                                     , gp_practices_swab = 300 
                                     , gp_swabs_mg = i 
                                     , pop_eng = 5.7106398e7 
                                     , gp_ari_consults = 327 
                                     , gp_ari_swabs = 747 
                                     , pc_swab_turnaround_time = [2,4] 
                                     , pathogen_type = "virus"
                                    )

            1200×9 DataFrame
            Row │ sim_n  GP_TD     GP_3TD       n_GP_cases      n_GP_cases_sampled          n_GP_cases_sampled_positive    GP_simid                           tinf_relating_to_GP_TD  last_tinf_relating_to_GP_3TD 
                │ Int64  Float64   Float64          Int64              Int64                           Int64               String                             Float64                 Float64
            ──────┼─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
                1 │     1   83.4112  Inf              171                   1                            1              ["676d6880-b8a9-11f0-016a-edfefc…                 76.1573     Inf
                2 │     2   87.1256  Inf              705                   1                            1              ["71ea31dc-b8a9-11f0-22a2-51bf2b…                 67.0375     Inf
                3 │     3   76.2388   81.1789       29465                  56                           51              ["7959a25e-b8a9-11f0-386b-77359c…                 70.5452     73.1876
                4 │     4  Inf       Inf              207                   0                            0              ["5346f5de-b8aa-11f0-2a4d-e1fc8c…                Inf          Inf
                5 │     5   67.2104   72.3096       81084                 127                          116              ["55ea4ffe-b8aa-11f0-35cc-3b70d9…                 55.1878     58.1821
                6 │     6   82.5854   90.8711        6083                   7                            6              ["bb622798-b8ac-11f0-3a44-b9eaff…                 76.704      86.4411
                7 │     7   40.5837   69.4841      115934                 197                          169              ["edcffff4-b8ac-11f0-1ce3-9ba211…                 26.0654     62.9518
                8 │     8  Inf       Inf              215                   0                            0              ["99674cf2-b8b0-11f0-0ba2-bf1638…                Inf          Inf
              ⋮   │   ⋮       ⋮         ⋮          ⋮               ⋮                        ⋮                                            ⋮                                      ⋮           ⋮
             1194 │  1194   80.2762  101.425          814                   4                            4              ["5ba8d488-b8f3-11f0-36c3-31aa2b…                 73.5756     87.6789
             1195 │  1195  Inf       Inf              223                   0                            0              ["61b38152-b8f3-11f0-01e0-39d3bc…                Inf          Inf
             1196 │  1196   80.815    88.2466       24827                  32                           28              ["63626752-b8f3-11f0-2101-4f08da…                 76.1437     78.3087
             1197 │  1197   80.1165   84.0567       32851                  39                           34              ["1b615b74-b8f4-11f0-24c7-d3c01a…                 69.1657     76.3517
             1198 │  1198   78.6828   87.9033        9216                  18                           14              ["10ad8cba-b8f5-11f0-019c-b78e40…                 72.5925     80.8266
             1199 │  1199   66.3393   69.1236       73811                 109                           94              ["536f0c54-b8f5-11f0-2c2f-e9af04…                 57.2872     62.7019
             1200 │  1200  Inf       Inf              415                   0                            0              ["7579ea10-b8f7-11f0-2219-3dd9d2…                Inf          Inf

"""
function gp_td(;  p=NBPMscape.P
                  , sims 
                  # Parameters for existing Oxford-RCGP RSC primary care surveillance
                  , gp_practices_total = NBPMscape.P.gp_practices_total #6199 # Total number of GP practices in England at July 2025. Source: BMA analysis (https://www.bma.org.uk/advice-and-support/nhs-delivery-and-workforce/pressures/pressures-in-general-practice-data-analysis) of NHS Digital General Practice Workforce Statistics (https://digital.nhs.uk/data-and-information/publications/statistical/general-and-personal-medical-services) [Accessed 2 Sep 2025]  
                  , gp_practices_swab = NBPMscape.P.gp_practices_swab #300 # Number of GP practices taking swabs for virology surveillance. Source: Data quality report: national flu and COVID-19 surveillance report (27 May 2025)
                  , gp_swabs_mg = NBPMscape.P.gp_swabs_mg #100 # 200 [319, 747]#[25698,46685] # Assumed number of swabs that are metagenomic sequenced for investigating impact
                  , pop_eng = NBPMscape.P.pop_eng #5.7106398e7 # Population of England. Source: ONS mid-year 2022 UK population data disaggregated by various geo levels and age groups and gender. [Accessed 6 November 2024] Available at https://www.ons.gov.uk/peoplepopulationandcommunity/populationandmigration/populationestimates/datasets/populationestimatesforukenglandandwalesscotlandandnorthernireland 
                  , gp_ari_consults = NBPMscape.P.gp_ari_consults #180 # 327 # Number of ARI consultations per 100k of England population per week [mean summer 2024, mean winter 2024/25]. Source: Analysis of data extracted from RCGP Research & Surveillance Centre (RSC) Virology Dashboard [Accessed 29 Aug 2025]
                  , gp_ari_swabs = NBPMscape.P.gp_ari_swabs #319 # [319,747] #[25698,46685]# Number of swabs taken from suspected ARI per week [mean summer 2024, mean winter 2024/25]. Source: Analysis of data extracted from RCGP Research & Surveillance Centre (RSC) Virology Dashboard [Accessed 29 Aug 2025]
                  , pc_swab_turnaround_time = NBPMscape.P.turnaroundtime_rcgp #[2,4] # [min,max] number of days between swab sample being taken and results received. Source: Data quality report: national flu and COVID-19 surveillance report (27 May 2025).  Assume same for metagenomic testing.
                  , pathogen_type = NBPMscape.P.pathogen_type #"virus"
                  )
    
    ### Probability of primary care surveillance events
    # Probability of visiting a General Practice that takes surveillance swabs
    prob_visiting_gp_swab = minimum( [1, gp_practices_swab / gp_practices_total ])
    
    # Probability of being swabbed given visited a GP that takes suriveillance swabs
    prob_swabbed_g_visit_gp_swab = minimum([1, (1e5*gp_ari_swabs) / (pop_eng*gp_ari_consults) / prob_visiting_gp_swab])
    
    # Probability that swab is selected for metagenomic sequencing given that patient has been selected for a surveillance swab
    prob_mg_g_swabbed = minimum([1, gp_swabs_mg / gp_ari_swabs ])
    
    # Primary care sampling proportion - equivalent to probability a swab being taken and a metagenomic sequence being taken if attend a GP
    p_pc_mg_g_gp_visit = prob_visiting_gp_swab * prob_swabbed_g_visit_gp_swab * prob_mg_g_swabbed
    
    # Check if primary care sampling proportion, p_pc, was entered as an input to the function
    p_pc = p_pc_mg_g_gp_visit
    
    ## How many replicates are in this simulation
    n_replicates = length(sims)
    
    # Create df to store results
    col_any = Vector{Any}(undef, n_replicates)
    fill!(col_any, Inf)
    sim_tds_cols = [copy(col_any) for _ in 1:9]
    sim_tds = DataFrame( sim_tds_cols, ["sim_n"
                                        # Time to detection for 1 case
                                        ,"GP_TD"      # primary care sampling only
                                        # Time to detection for 3 cases
                                        ,"GP_3TD"  # primary care sampling only 
                                        # Number of cases, samples and positive samples
                                        ,"n_GP_cases","n_GP_cases_sampled","n_GP_cases_sampled_positive"
                                        # simids
                                        ,"GP_simid"
                                        , "tinf_relating_to_GP_TD" # useful for investigating impact of early termination of simtree within simforest
                                        , "last_tinf_relating_to_GP_3TD" # useful for investigating impact of early termination of simtree within simforest
                                        
                                        ])

    # Loop through simulation replicates, sampling infection cases (adapted from sampleforest() function in 'core.jl' and 'median_TD_by_region.jl'),
    # and computing times to detection
    for s in 1:n_replicates
        #TEST println("$(s) of $(n_replicates) replicates")

        # Add simulation number to results df
        #TEST s=1
        sim_tds[s,:sim_n] = s
        fo = sims[s]
                   
        # Filter for ICU cases
        gp_cases = size(fo,1) == 0 ? DataFrame() : fo[ isfinite.(fo.tgp), : ]
        
        # Add simids to df
        sim_tds[s,:GP_simid] = size(gp_cases,1) > 0 ? unique( gp_cases.simid ) : missing 
        #if size(gp_cases,1) > 0
        #    simids_gp = unique( gp_cases.simid )  
        #    sim_tds[s,:PC_simid] = simids_gp
        #end

        ### Record infection time stats
        # Record time stats for GP cases
        
        # Also apply metagenomic test sensitivity so that the subsampled GP cases are those with positive test results
        # Obtain metagenomic test sensitivity for pathogen_type
        sensitivities = Dict(
                            "virus" => p.sensitivity_mg_virus,
                            "bacteria" => p.sensitivity_mg_bacteria,
                            "fungi" => p.sensitivity_mg_fungi
                            )
        mg_test_sensitivity = get(sensitivities, pathogen_type) do
            error("Unknown pathogen type: $pathogen_type")
        end

        # Sample sizes based on probability of infected case having a metagenomic sample taken at the GP
        n_gp_samples = rand( Binomial( size(gp_cases,1), p_pc )) #* mg_test_sensitivity) )
        n_gp_samples_positive = rand( Binomial( n_gp_samples, mg_test_sensitivity) )
        #println(n_gp)
            
        # Subsample of GP cases
        gp_cases_sub = gp_cases[sample( 1:size(gp_cases,1), n_gp_samples_positive, replace=false ), :]
        
        # Record number of cases in the sim rep and the number that were tested
        sim_tds[s,:n_GP_cases] = size(gp_cases,1)
        sim_tds[s,:n_GP_cases_sampled] = n_gp_samples
        sim_tds[s,:n_GP_cases_sampled_positive] = n_gp_samples_positive

        # If there are NO GP cases in the sample
        if size(gp_cases_sub,1) == 0
                
            ## Add to results df for ICU times to detection
            # Record Inf times to detection if no GP cases sampled:
            #  1 primary care case for TD
            sim_tds[s, :GP_TD ] = Inf
            #  3 primary care cases for 3TD
            sim_tds[s, :GP_3TD ] = Inf
            # and for tinf
            sim_tds[s,:tinf_relating_to_GP_TD] = Inf
            sim_tds[s,:last_tinf_relating_to_GP_3TD] = Inf

        elseif size(gp_cases_sub,1) > 0 # If there ARE GP cases in the sample
            
            # Generate sample times
            gp_tsample = gp_cases_sub.tgp
            gp_cases_sub.tsample = gp_tsample # primary care sampling is done at the time of GP visit (ignores home sampling)

            # Simulate reports times
            gp_treport = (gp_tsample .+ rand(Uniform(pc_swab_turnaround_time[1], pc_swab_turnaround_time[2])) )
            gp_cases_sub.treport = gp_treport

            ## Find cases with 3 shortest times to detection (TD)
            # Primary care cases only
            gp_cases_sub_top3_td = first(sort(gp_cases_sub, :treport), 3)
            
            ## Add to results df for primary care times to detection
            # Record the times to detection for earliest:
            # Primary care cases only
            #  1 primary care case for TD
            sim_tds[s, :GP_TD ] = gp_cases_sub_top3_td[1,:treport]
            # tinf
            sim_tds[s,:tinf_relating_to_GP_TD] = gp_cases_sub_top3_td[1,:tinf]
            #  3 primary care cases for 3TD
            #if size(gp_cases_sub_top3_td,1) >= 3
            #    sim_tds[s, :GP_3TD ] = gp_cases_sub_top3_td[3,:treport]
            #else
            #    sim_tds[s, :GP_3TD ] = Inf 
            #end
            sim_tds[s, :GP_3TD ] = size(gp_cases_sub_top3_td,1) >= 3 ? gp_cases_sub_top3_td[3,:treport] : Inf
            sim_tds[s,:last_tinf_relating_to_GP_3TD] = size(gp_cases_sub_top3_td,1) >= 3 ? gp_cases_sub_top3_td[3,:tinf] : Inf

        end # End of if loop checking whether there are any GP cases sampled

    end # End of loop through different simulations
    return sim_tds
end # End of gp_td function


"""
Function: secondary_care_td

Description:    Function to sample infections from hospital admissions and return time to detection for a pathogen X. 
                Parameters modelled on expansion of the Hospital-based acute respiratory infection sentinel surveillance (HARISS) system.
                Description of existing HARISS network and surveillance:
                - https://www.gov.uk/government/publications/sources-of-surveillance-data-for-influenza-covid-19-and-other-respiratory-viruses/data-quality-report-national-flu-and-covid-19-surveillance-report#hospital-based-acute-respiratory-infection-sentinel-surveillance-hariss-system-additional-surveillance-system
                - Symes et al (2025), 'Estimating the disease burden of respiratory syncytial virus (RSV)
                  in older adults in England during the 2023/24 season: a new national hospital-based surveillance system', 
                  medrxiv, doi: https://doi.org/10.1101/2025.04.17.25325639

Arguments:  p::NamedTuple   Set of parameters relating to the outbreak and sampling. Parameters loaded from defaults in core.jl and
                            some possibly from the configuration .yaml file. NBPMscape.P
            sims            Formatted as a vector of dataframes. Simulated outbreak data in the G dataframe output from simtree or simforest functions.
                            Example format:
                            17113×18 DataFrame
                            Row │ pid                                tinf     tgp       ted       thospital  ticu     tstepdown  tdischarge  trecovered  tdeceased  severity               fatal  iscommuter  homeregion  infector_age  infectee_age  importedinfection  simid                             
                                │ String                             Float64  Float64   Float64   Float64    Float64  Float64    Float64     Float64     Float64    Symbol                 Bool   Bool        String      Int8?         Int8          Bool               String
                            ───────┼──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
                                1 │ fec10c2a-fd4d-11f0-3874-5958aa14…  18.5291   24.2016  Inf        Inf           Inf        Inf    Inf          52.4344        Inf  moderate_GP            false       false  TLH5                  64            67              false  fec10c2a-fd4d-11f0-3874-5958aa14…
                                2 │ fec10c2a-fd4d-11f0-3874-5958aa14…  24.4748  Inf       Inf         27.4966      Inf        Inf     39.2252     39.2803        Inf  severe_hosp_long_stay  false       false  TLH5                  44            84              false  fec10c2a-fd4d-11f0-3874-5958aa14…
                                3 │ fec10c2a-fd4d-11f0-3874-5958aa14…  24.79     28.016   Inf        Inf           Inf        Inf    Inf          54.7924        Inf  moderate_GP            false        true  TLH5                  44            41              false  fec10c2a-fd4d-11f0-3874-5958aa14…
                            ⋮   │                 ⋮                     ⋮        ⋮         ⋮          ⋮         ⋮         ⋮          ⋮           ⋮           ⋮                ⋮              ⋮        ⋮           ⋮            ⋮             ⋮                ⋮                          ⋮
                            17111 │ 63ff8da8-fd4e-11f0-2087-f7ff9459…  54.3506   57.7829  Inf        Inf           Inf        Inf    Inf          80.6498        Inf  moderate_GP            false       false  TLD3                  15            16              false  63ff8da8-fd4e-11f0-2087-f7ff9459…
                            17112 │ 63ff8da8-fd4e-11f0-2087-f7ff9459…  59.861    65.7573   65.8896   Inf           Inf        Inf     66.2212     86.8499        Inf  moderate_ED            false        true  TLD3                   9            36              false  63ff8da8-fd4e-11f0-2087-f7ff9459…
                            17113 │ 63ff8da8-fd4e-11f0-2087-f7ff9459…  58.9611   60.9347  Inf         61.1105      Inf        Inf     73.3345     79.2828        Inf  severe_hosp_long_stay  false       false  TLD3                  92            89              false  63ff8da8-fd4e-11f0-2087-f7ff9459…
                            
            pathogen_type::String       "virus", "bacteria" or "fungi". Used to determine the metagenomic test sensitivity to apply and thereby
                                        change the probability of a positive sample.
            initial_dow::Int64          Day of the week that initial case imported. Day of week codes: Sunday =1, Monday = 2, ..., Saturday = 7 
            hosp_cases::DataFrame       Filtered from simulation of infections (G df filtered for value in thosp column)
            # Sampling parameters #
            hariss_courier_to_analysis::Float64     Time between courier collection from PHL and beginning of analysis
            hariss_turnaround_time::Vector          Time to process sample and report results / declare detection
            n_hosp_samples_per_week::Int            Total number of hospital samples to be taken per week
            sample_allocation::String               Method of allocating the total number of samples required across the public health laboratories.
                                                    Options are either "equal" (i.e. equal weighting or no weighting) or "weighted" using the method in 
                                                    'weight_samples_by' argument.
            sample_proportion_adult::Any            Indicates split of sample target between adults and children. Weighting to samples from adults (aged 16 and over)
                                                    as a numeric decimal, e.g. 0.75 with the remaining 0.25 to children. 
                                                    Alternatively "free" indicates that no split is specified and the samples taken are the most recent prior to courier collection.
            hariss_nhs_trust_sampling_sites::DataFrame  List of NHS Trusts in HARISS sampling network 
            weight_samples_by::String               Options are "ae_mean" or "catchment_pop". This defines the method of weighting the total samples across
                                                    the public health labs (PHLs) in the HARISS netork. "ae_mean" uses the mean number of A&E attendances at NHS Trusts that 
                                                    feed into each PHL. "catchment_pop" uses the populations in the catchment areas of the NHS Trusts.
            phl_collection_dow::Vector{Int64}       Day(s) of the week that swab samples will be collected from public health labs. Day of week codes: Sunday = 1,... Saturday = 7.
            phl_collection_time::Float64            Decimal value time of day that swab samples are collected from public health labs. E.g. midday = 0.5 (=12/24), 9am = 0.375(=9/24), 5pm = 0.7083 (=17/24)
            hosp_to_phl_cutoff_time_relative::Real  The time period (in days) between the cut-off time for the last swab that can be taken at a hospital and 
                                                    still be transported to the local public health lab (PHL) to be made available for courier collection at the phl_collection_time. 
                                                    E.g. hosp_to_phl_cutoff_time_relative = 1 and phl_collection_time = 0.5, indicates midday collection from the PHL and cut-off swab time
                                                    at the hospital of midday the day before. As another example, midday (12:00) collection from PHL (phl_collection_time = 0.5) and
                                                    the cutoff for swabs at 5pm (17:00) the day before would be hosp_to_phl_cutoff_time_relative = 0.7917 (=(12+(24-17))/24)
            swab_time_mode::Real                    The peak time (in days) for swabbing after initial attendance/admission at hospital.
            swab_proportion_at_48h::Real            The proportion of acute respiratory infection (ARI) hospital attendees/admissions are swabbed within 48hrs (=2 days) of arrival.
            proportion_hosp_swabbed::Real           Proportion of ARI hospital attendees that are swabbed, as a decimal value
            only_sample_before_death::Bool          There is a possibilty of swabbing time being drawn after death so 'true' here will constrain tswab to tdeceased
            # Hospital parameters #
            ed_discharge_limit::Float64             Assume that people attending the Emergency Department only (i.e. not admitted) are discharged within this time limit (in days).
            hosp_short_stay_limit::Float64          Assume that people that are admitted to hospital but only have a 'short stay' (this relates to a specific infection severity) 
                                                    are discharged within this time limit (in days).
            nhs_trust_catchment_pop::DataFrame      Proportion of population (of England) split by NHS Trust. Also disaggregated by adults and children. Example below:
                                                    126×4 DataFrame
                                                    Row │ TrustCode  catchment_prop_of_total_sum  prop_child   prop_adult 
                                                        │ Any        Float64                      Float64      Float64
                                                    ─────┼─────────────────────────────────────────────────────────────────
                                                    1 │ R0A                         0.021291    0.00634217   0.0149488
                                                    2 │ R0B                         0.00772294  0.00146184   0.0062611
                                                    3 │ R0D                         0.0104117   0.0018323    0.00857945
                                                    ⋮  │     ⋮                   ⋮                    ⋮           ⋮
                                                    124 │ RXW                         0.00811911  0.00169838   0.00642073
                                                    125 │ RYJ                         0.0100881   0.00187858   0.00820951
                                                    126 │ RYR                         0.0179986   0.00352307   0.0144755
            nhs_trust_ae_12m::DataFrame             Data on A&E attendances at NHS Trusts including the 12-month mean and proportions of the total.
                                                    199×15 DataFrame
                                                    Row │ NHS_Trust_code  2024_4_function  2024_5_function  2024_6_function  2024_7_function  2024_8_function  2024_9_function  2024_10_function  2024_11_function  2024_12_function  2025_1_function  2025_2_function  2025_3_function  mean_12m  mean_12m_prop 
                                                        │ String7         Int64            Int64            Int64            Int64            Int64            Int64            Int64             Int64             Int64             Int64            Int64            Int64            Float64   Float64
                                                    ─────┼────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
                                                    ⋮  │       ⋮                ⋮                ⋮                ⋮                ⋮                ⋮                ⋮                ⋮                 ⋮                 ⋮                 ⋮                ⋮                ⋮            ⋮            ⋮
                                                    32 │ RGP                        6950             7644             7323             7844             7757             7222              7574              7287              7082             6535             6286             7416   7243.33     0.00530124
                                                    33 │ RJ6                       12926            13703            13337            13230            12064            12509             12732             13200             13668            12997            11931            13255  12962.7      0.00948709
                                                    34 │ RK9                        8342             9272             9105             9222             8907             8884              9604              9116              9586             8603             8072             9353   9005.5      0.00659093
                                                    35 │ RXF                       17513            17925            16962            16592            15391            15790             16821             16762             16641            15436            14812            17104  16479.1      0.0120607
                                                    36 │ RXR                       12318            13526            12826            13481            12519            13080             13289             13334             13318            12294            11973            13340  12941.5      0.0094716
                                                    ⋮  │       ⋮                ⋮                ⋮                ⋮                ⋮                ⋮                ⋮                ⋮                 ⋮                 ⋮                 ⋮                ⋮                ⋮            ⋮            ⋮
            # Seasonal values #
            hosp_ari_admissions::Int                Estimate of weekly hospital ARI admissions (excluding pathogen X being simulated)
            hosp_ari_admissions_adult_p::Float64    Proportion of Emergency Department ARI admissions that are adults (18y and over)
            hosp_ari_admissions_child_p::Float64    Proportion of Emergency Department ARI admissions that are children (<18y)
            ed_ari_destinations_adult::DataFrame    Destinations and proportions of adults arriving at hospital emergency departments.
                                                    For example: DataFrame( destination = [:discharged,:short_stay,:longer_stay]
                                                                          , proportion_of_attendances = [0.628,0.030,0.342])
            ed_ari_destinations_child::DataFrame    Destinations and proportions of children arriving at hospital emergency departments.
                                                    For example: DataFrame( destination = [:discharged,:short_stay,:longer_stay]
                                                                          , proportion_of_attendances = [0.861,0.014,0.125])
                                                    
Returns:    A dataframe containing the median times to detect one and three cases for each simulation replicate in sims.
            The output also contains the number of secondary care cases, the number that are sampled and the number that 
            return a positive result based on the assumed metagenomic test sensitivity. The simids are also recorded (a unique code for each simtree run within the simforest function).
            Information regarding the times of infection relating to the time of detection for one and three cases is also recorded.
            Output dataframe is formatted as example below.

            340×9 DataFrame
            Row │ sim_n  SC_TD    SC_3TD   n_SC_cases  n_SC_cases_sampled  n_SC_cases_sampled_positive  SC_simid                           tinf_relating_to_SC_TD  last_tinf_relating_to_SC_3TD 
                │ Any    Any      Any      Any         Any                 Any                          Any                                Any                     Any
            ─────┼───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
            1 │ 1      33.7966  36.7966  11115       376                 334                          ["fec10c2a-fd4d-11f0-3874-5958aa…  26.3112                 25.3633
            2 │ 2      21.8271  28.8271  6186        269                 239                          ["649b0062-fd4e-11f0-0308-e19372…  10.4537                 21.8793
            3 │ 3      42.681   53.681   1681        100                 91                           ["8e2cfe26-fd4e-11f0-2ebd-83a077…  30.0985                 42.6665
            ⋮  │   ⋮       ⋮        ⋮         ⋮               ⋮                        ⋮                               ⋮                            ⋮                          ⋮
            338 │ 338    40.3644  43.3644  6684        313                 283                          ["2b910e9e-fd4f-11f0-2596-a17c08…  32.0248                 32.3078
            339 │ 339    34.4027  37.4027  9849        410                 369                          ["507353ca-fd4f-11f0-0e65-a16a8e…  20.9622                 23.8261
            340 │ 340    43.5813  50.5813  2114        119                 103                          ["88caa1b2-fd4f-11f0-1320-6db9fb…  33.0952                 39.0238


Examples:
            sims = load("covidlike-1.4.3-sims-filtered_G_nrep340_1633695.jld2", "sims")
            results1 = secondary_care_td(; p = NBPMscape.P
                            , sims
                            , pathogen_type = NBPMscape.P.pathogen_type
                            , initial_dow = NBPMscape.P.initial_dow
                            # Sampling parameters
                            , hariss_courier_to_analysis = NBPMscape.P.hariss_courier_to_analysis #1.0 # Time between courier collection from PHL and beginning of analysis
                            , hariss_turnaround_time = NBPMscape.P.turnaroundtime_hariss #[2,4] # Time to process sample and report results / declare detection
                            , n_hosp_samples_per_week = NBPMscape.P.n_hosp_samples_per_week # Total number of hospital samples to be taken per week
                            , sample_allocation = NBPMscape.P.sample_allocation #"equal" # "equal" or "weighted"
                            , sample_proportion_adult = NBPMscape.P.sample_proportion_adult #"free" # "free" or numeric decimal, e.g. 0.75. Indicates split of sample target 
                                                                    # between adults and children. "free" indicates that no split is specified
                            , hariss_nhs_trust_sampling_sites = NBPMscape.P.hariss_nhs_trust_sampling_sites # List of NHS Trusts in HARISS sampling network 
                            , weight_samples_by = NBPMscape.P.weight_samples_by #"ae_mean" # or "catchment_pop"
                            , phl_collection_dow = NBPMscape.P.phl_collection_dow #[2,5] # Day(s) of week that swab samples will be collected from public health labs. Day of week codes: Sunday = 1,... Saturday = 7.
                            , phl_collection_time = NBPMscape.P.phl_collection_time # Decimal value time of day that swab samples are collected from public health labs. E.g. midday = 0.5 (=12/24), 9am = 0.375(=9/24), 5pm = 0.7083 (=17/24)
                            , hosp_to_phl_cutoff_time_relative = NBPMscape.P.hosp_to_phl_cutoff_time_relative # day(s). The cut-off time for the last swab that can be taken at a hospital and transported to the local public health lab (PHL) and made available for courier collection at the phl_collection_time. E.g. hosp_to_phl_cutoff_time_relative = 1 and phl_collection_time = 0.5, indicates midday collection from the PHL and cut-off swab time at the hospital of midday the day before. As another example, midday (12:00) collection from PHL (phl_collection_time = 0.5) and the cutoff for swabs at 5pm (17:00) the day before would be hosp_to_phl_cutoff_time_relative = 0.7917 (=(12+(24-17))/24)
                            , swab_time_mode = NBPMscape.P.swab_time_mode #0.25 # Assume swabbing peaks at 6hrs (=0.25 days) after attendance/admission at hospital
                            , swab_proportion_at_48h = NBPMscape.P.swab_proportion_at_48h #0.9 # Assume 90% of swabs are taken within 48hrs (=2 days) of attendance/admission at hospital
                            , proportion_hosp_swabbed = NBPMscape.P.proportion_hosp_swabbed #0.9 # Assume 90% of ARI attendances are swabbed
                            #, initial_dow = 1 # Day of the week that initial case imported. Day of week codes: Sunday =1, Monday = 2, ..., Saturday = 7 
                            , only_sample_before_death = NBPMscape.P.hariss_only_sample_before_death #true # There is a possibilty of swabbing time being drawn after death so 'true' here will constrain tswab to tdeceased
                            # Hospital parameters
                            , ed_discharge_limit = NBPMscape.P.tdischarge_ed_upper_limit #0.25 # days. Assume that people attending the Emergency Department are discharged within this time limit.
                            , hosp_short_stay_limit = NBPMscape.P.tdischarge_hosp_short_stay_upper_limit # days. Assume that people attending the Emergency Department are discharged within this time limit.
                            , nhs_trust_catchment_pop = NHS_TRUST_CATCHMENT_POP_ADULT_CHILD
                            #, nhs_trust_ae_12m = AE_12M
                            # Seasonal values
                            , hosp_ari_admissions = NBPMscape.P.hosp_ari_admissions # Estimate of weekly hospital ARI admissions (excluding pathogen X being simulated)
                            , hosp_ari_admissions_adult_p = NBPMscape.P.hosp_ari_admissions_adult_p #0.52# Proportion of ED ARI admissions that are adults (18y and over)
                            , hosp_ari_admissions_child_p = NBPMscape.P.hosp_ari_admissions_child_p #0.48 # Proportion of ED ARI admissions that are children (<18y)
                            , ed_ari_destinations_adult = NBPMscape.P.ed_ari_destinations_adult #DataFrame( destination = [:discharged,:short_stay,:longer_stay]
                                                                                                                    #, proportion_of_attendances = [0.628,0.030,0.342])
                            , ed_ari_destinations_child = NBPMscape.P.ed_ari_destinations_child #DataFrame( destination = [:discharged,:short_stay,:longer_stay]
                                                                                                                    #, proportion_of_attendances = [0.861,0.014,0.125])
                            )

"""
function secondary_care_td(; p = NBPMscape.P
                            , sims
                            , pathogen_type::String = NBPMscape.P.pathogen_type #"virus"
                            , initial_dow::Int64 = NBPMscape.P.initial_dow #1 # Day of the week that initial case imported. Day of week codes: Sunday =1, Monday = 2, ..., Saturday = 7 
                            #, hosp_cases::DataFrame # filtered from simulation of infections (G df filtered for value in thosp column)
                            # Sampling parameters
                            , hariss_courier_to_analysis::Float64 = NBPMscape.P.hariss_courier_to_analysis #1.0 # Time between courier collection from PHL and beginning of analysis
                            , hariss_turnaround_time::Vector = NBPMscape.P.turnaroundtime_hariss #[2,4] # Time to process sample and report results / declare detection
                            , n_hosp_samples_per_week::Int = NBPMscape.P.n_hosp_samples_per_week # Total number of hospital samples to be taken per week
                            , sample_allocation::String = NBPMscape.P.sample_allocation #"equal" # "equal" or "weighted"
                            , sample_proportion_adult::Any = NBPMscape.P.sample_proportion_adult #"free" # "free" or numeric decimal, e.g. 0.75. Indicates split of sample target 
                                                                    # between adults and children. "free" indicates that no split is specified
                            , hariss_nhs_trust_sampling_sites::DataFrame = NBPMscape.P.hariss_nhs_trust_sampling_sites # List of NHS Trusts in HARISS sampling network 
                            , weight_samples_by::String = NBPMscape.P.weight_samples_by #"ae_mean" # or "catchment_pop"
                            , phl_collection_dow::Vector{Int64} = NBPMscape.P.phl_collection_dow #[2,5] # Day(s) of week that swab samples will be collected from public health labs. Day of week codes: Sunday = 1,... Saturday = 7.
                            , phl_collection_time::Float64 = NBPMscape.P.phl_collection_time # Decimal value time of day that swab samples are collected from public health labs. E.g. midday = 0.5 (=12/24), 9am = 0.375(=9/24), 5pm = 0.7083 (=17/24)
                            , hosp_to_phl_cutoff_time_relative::Real = NBPMscape.P.hosp_to_phl_cutoff_time_relative # day(s). The cut-off time for the last swab that can be taken at a hospital and transported to the local public health lab (PHL) and made available for courier collection at the phl_collection_time. E.g. hosp_to_phl_cutoff_time_relative = 1 and phl_collection_time = 0.5, indicates midday collection from the PHL and cut-off swab time at the hospital of midday the day before. As another example, midday (12:00) collection from PHL (phl_collection_time = 0.5) and the cutoff for swabs at 5pm (17:00) the day before would be hosp_to_phl_cutoff_time_relative = 0.7917 (=(12+(24-17))/24)
                            , swab_time_mode::Real = NBPMscape.P.swab_time_mode #0.25 # Assume swabbing peaks at 6hrs (=0.25 days) after attendance/admission at hospital
                            , swab_proportion_at_48h::Real = NBPMscape.P.swab_proportion_at_48h #0.9 # Assume 90% of swabs are taken within 48hrs (=2 days) of attendance/admission at hospital
                            , proportion_hosp_swabbed::Real = NBPMscape.P.proportion_hosp_swabbed #0.9 # Assume 90% of ARI attendances are swabbed
                            , only_sample_before_death::Bool = NBPMscape.P.hariss_only_sample_before_death #true # There is a possibilty of swabbing time being drawn after death so 'true' here will constrain tswab to tdeceased
                            # Hospital parameters
                            , ed_discharge_limit::Float64 = NBPMscape.P.tdischarge_ed_upper_limit #0.25 # days. Assume that people attending the Emergency Department are discharged within this time limit.
                            , hosp_short_stay_limit::Float64 = NBPMscape.P.tdischarge_hosp_short_stay_upper_limit # days. Assume that people attending the Emergency Department are discharged within this time limit.
                            , nhs_trust_catchment_pop::DataFrame = NHS_TRUST_CATCHMENT_POP_ADULT_CHILD
                            #, nhs_trust_ae_12m::DataFrame = AE_12M
                            # Seasonal values
                            , hosp_ari_admissions::Int = NBPMscape.P.hosp_ari_admissions # Estimate of weekly hospital ARI admissions (excluding pathogen X being simulated)
                            , hosp_ari_admissions_adult_p::Float64 = NBPMscape.P.hosp_ari_admissions_adult_p #0.52# Proportion of ED ARI admissions that are adults (18y and over)
                            , hosp_ari_admissions_child_p::Float64 = NBPMscape.P.hosp_ari_admissions_child_p #0.48 # Proportion of ED ARI admissions that are children (<18y)
                            , ed_ari_destinations_adult::DataFrame = NBPMscape.P.ed_ari_destinations_adult #DataFrame( destination = [:discharged,:short_stay,:longer_stay]
                                                                                                                    #, proportion_of_attendances = [0.628,0.030,0.342])
                            , ed_ari_destinations_child::DataFrame = NBPMscape.P.ed_ari_destinations_child #DataFrame( destination = [:discharged,:short_stay,:longer_stay]
                                                                                                                    #, proportion_of_attendances = [0.861,0.014,0.125])
                            )

    n_replicates = length(sims)
    
    # Create df to store results
    col_any = Vector{Any}(undef, n_replicates)
    fill!(col_any, Inf)
    sim_tds_cols = [copy(col_any) for _ in 1:9]
    sim_tds = DataFrame( sim_tds_cols, ["sim_n"
                                        # Time to detection for 1 case
                                        ,"SC_TD"     # ICU sampling only
                                        # Time to detection for 3 cases
                                        ,"SC_3TD" # ICU sampling only
                                        # Number of cases and samples
                                        ,"n_SC_cases","n_SC_cases_sampled", "n_SC_cases_sampled_positive"
                                        ,"SC_simid"
                                        , "tinf_relating_to_SC_TD" # useful for adjusting for right censoring (due to simulation max time) and investigating impact of early termination of simtree within simforest
                                        , "last_tinf_relating_to_SC_3TD" # useful for adjusting for right censoring (due to simulation max time) and investigating impact of early termination of simtree within simforest
                                        ])

    # Loop through simulation replicates, sampling infection cases (adapted from sampleforest() function in 'core.jl' and 'median_TD_by_region.jl'),
    # and computing times to detection
    for s in 1:n_replicates # s=1
        #TEST
        println("$(s) of $(n_replicates) replicates")

        # Add simulation number to results df
        #TEST s=1
        sim_tds[s,:sim_n] = s
        fo = sims[s]
                   
        # Filter for simulated hospital admissions and emergency department (ED) only  cases
        hosp_cases = size(fo,1) == 0 ? DataFrame() : fo[ isfinite.(fo.ted) .| isfinite.(fo.thospital) , : ]
        
        # Add simids
        if size(hosp_cases,1) > 0
            simids_hosp = unique( hosp_cases.simid )  #unique( [hosp_cases[1,:simid]] ) 
            sim_tds[s,:SC_simid] = simids_hosp #unique( [hosp_cases[1,:simid] ]) 
        end
    
        ### Record infection time stats
        # Note that metagenomic test sensitivity is not applied in sample_hosp_cases_n # println(hosp_cases_sub)
        hosp_cases_sub = sample_hosp_cases_n(; p = NBPMscape.P
                                             , hosp_cases = hosp_cases
                                             # Hospital parameters
                                             , nhs_trust_catchment_pop = nhs_trust_catchment_pop
                                             , hosp_short_stay_limit = hosp_short_stay_limit
                                             #, nhs_trust_ae_12m = nhs_trust_ae_12m
                                             , ed_discharge_limit = ed_discharge_limit
                                             # Seasonal parameters
                                             , hosp_ari_admissions = hosp_ari_admissions
                                             , hosp_ari_admissions_adult_p = hosp_ari_admissions_adult_p
                                             , hosp_ari_admissions_child_p = hosp_ari_admissions_child_p
                                             , ed_ari_destinations_adult = ed_ari_destinations_adult
                                             , ed_ari_destinations_child = ed_ari_destinations_child
                                             # Sample parameters
                                             , n_hosp_samples_per_week = n_hosp_samples_per_week
                                             , sample_allocation = sample_allocation
                                             , sample_proportion_adult = sample_proportion_adult
                                             , hariss_nhs_trust_sampling_sites = hariss_nhs_trust_sampling_sites
                                             , weight_samples_by = weight_samples_by
                                             , phl_collection_dow = phl_collection_dow
                                             , phl_collection_time = phl_collection_time
                                             , hosp_to_phl_cutoff_time_relative = hosp_to_phl_cutoff_time_relative
                                             , swab_time_mode = swab_time_mode
                                             , swab_proportion_at_48h = swab_proportion_at_48h
                                             , proportion_hosp_swabbed = proportion_hosp_swabbed
                                             , initial_dow = initial_dow
                                             , only_sample_before_death = only_sample_before_death
                                            )

        # Record number of cases in the sim rep and the number that were tested (record number that returned positive results below)
        sim_tds[s,:n_SC_cases] = size(hosp_cases,1)
        sim_tds[s,:n_SC_cases_sampled] = size(hosp_cases_sub,1) # hosp_cases_sub = DataFrame()

        # Obtain metagenomic test sensitivity for pathogen_type
        # Additional sub-sampling to account for metagenomic test sensitivity
        sensitivities = Dict(
                        "virus" => p.sensitivity_mg_virus,
                        "bacteria" => p.sensitivity_mg_bacteria,
                        "fungi" => p.sensitivity_mg_fungi
                        )
        mg_test_sensitivity = get(sensitivities, pathogen_type) do
            error("Unknown pathogen type: $pathogen_type")
        end
        # Positive sample sizes
        n_hosp =  rand( Binomial( size(hosp_cases_sub,1), mg_test_sensitivity ) )
        # Record number of sampled positive cases that returned positive results
        sim_tds[s,:n_SC_cases_sampled_positive] = n_hosp
        # Subsample of hospital & ED cases
        hosp_cases_sub = hosp_cases_sub[sample( 1:size(hosp_cases_sub,1), n_hosp, replace=false ), :]
        
        if size(hosp_cases_sub,1) == 0
                
            # Record zero cases and zero sampled
            #sim_tds[s,:n_SC_cases] = 0
            #sim_tds[s,:n_SC_cases_sampled] = 0
            #sim_tds[s,:n_SC_cases_sampled_positive] = 0

            # Record Inf as time to detection for earliest:
            # - 1 Secondary Care (SC) case
            sim_tds[s,:SC_TD] = Inf
            # - 3 Secondary Care (SC)  cases
            sim_tds[s,:SC_3TD] = Inf
            # and for tinf
            #sim_tds[s,:tinf_relating_to_SC_TD] = Inf
            #sim_tds[s,:last_tinf_relating_to_SC_3TD] = Inf

            # Define empty vector for top 3 times to detection - to be used later
            hosp_cases_sub_top3_td = []

        elseif size(hosp_cases_sub,1) > 0 
                
            # Simulate reports times
            hosp_treport = (hosp_cases_sub.tcourier .+ hariss_courier_to_analysis .+ rand(Uniform(hariss_turnaround_time[1], hariss_turnaround_time[2])) )
            hosp_cases_sub.treport = hosp_treport
            
            # Find cases with 3 shortest times to detection (TD)
            hosp_cases_sub_top3_td = first(sort(hosp_cases_sub, :treport), 3) #println(icu_cases_sub_top3_td.treport)

            ## Add to results df for ICU times to detection
            # Record the times to detection for earliest:
            # - 1 ICU case
            sim_tds[s,:SC_TD] = hosp_cases_sub_top3_td[1,:treport]
            # and for tinf
            sim_tds[s,:tinf_relating_to_SC_TD] = hosp_cases_sub_top3_td[1,:tinf]
            
            # - 3 ICU cases
            if size(hosp_cases_sub_top3_td,1) >= 3
                sim_tds[s,:SC_3TD] = hosp_cases_sub_top3_td[3,:treport]
                # and for tinf
                sim_tds[s,:last_tinf_relating_to_SC_3TD] = hosp_cases_sub_top3_td[3,:tinf]
            else
                sim_tds[s,:SC_3TD] = Inf
                sim_tds[s,:last_tinf_relating_to_SC_3TD] = Inf
            end
        end
        
    end # End of loop through different simulations
    return sim_tds
end # End of secondary_care_td function





# icu_vs_pc_td was used in a previous analysis but has been superseded by icu_td and gp_td
"""
Function: icu_v_pc_td 

Description:    Function to compute times to detection for 1 case (TD) and 3 cases (3TD)
                by sampling from:
                - primary care
                - ICU

# Arguments:    'sim_object_type':  Description of the format of the underlying simulation data. 
                                    There are three options: 
                                    -   "full", which is the data as output from {simtree} or {simforest}. 
                                        This option requires an input for "sims_file" and no input for 
                                        "sims_G_gp_filter_file" and "sims_G_icu_filter_file"
                                    -   "G_filtered_icu_gp_combined", which is the "full" data already filtered for the G dataframe and
                                        GP and ICU cases, but kept in a single file. This option requires an input for "sims_file" and no input for 
                                        "sims_G_gp_filter_file" and "sims_G_icu_filter_file"
                                    -   "G_filtered_icu_gp_separate", which is the "full" data already filtered for GP and ICU
                                         cases and split into two files. This option requires inputs for 
                                         "sims_G_gp_filter_file" and "sims_G_icu_filter_file" and no input for "sims_file"
                'file_or_object'    A binary choice between "file" and "object" that will determine whether or not the function
                                    needs to load a file. Loading files can be time consuming so if running the function
                                    multiple times, i.e. with different parameters, it is beneficial to load the file once
                                    and then input the object from the file in subsequent runs.
                'sims_file':                .jld2 file containing simulation replicates as output by {simtree} or {simforest}
                'sims_G_gp_filter_file':    .jld2 file containing simulation replicates (as output by {simtree} or {simforest}) BUT
                                            filtered for GP cases in the G dataframe (infection information) using {sims_filter}
                'sims_G_icu_filter_file':   .jld2 file containing simulation replicates (as output by {simtree} or {simforest}) BUT
                                            filtered for ICU cases in the G dataframe (infection information) using {sims_filter}
                'sims_object_name':     Name of object containing simulation replicates as output by {simtree} or {simforest},
                                        or as loaded from .jld2 file
                'sims_G_gp_filter_object_name': Name of object containing simulation replicates (as output by {simtree} or {simforest}) BUT
                                                filtered for GP cases in the G dataframe (infection information) using {sims_filter},
                                                or as loaded from .jld2 file
                'sims_G_icu_filter_object_name':    Name of object containing simulation replicates (as output by {simtree} or {simforest}) BUT
                                                    filtered for ICU cases in the G dataframe (infection information) using {sims_filter},
                                                    or as loaded from .jld2 file
                'icu_sample_type':      "regional" or "fixed". If "fixed" then p_icu will be used, note that test sensitivity and the practical sampling proportion are applied to this value so the positive sample proportion will be lower than p_icu. "regional" uses the {sample_icu_cases} function and is dependent on home region
                'p_icu':                ICU sampling proportion
                only_sample_before_death::Bool  true or false. Modelled sampling period is a uniform distribution between time of admission to ICU (ticu) 
                                                and 3 days after admission (ticu+3). It is possible for the time of death (tdeceased) to be before the sample time.
                                                Setting only_sample_before_death = true constrains the upper limit on tsample to be equal to tdeceased.
                'icu_ari_admissions':   Weekly ICU admission numbers. Formatted in a 2-element vector: [summer,winter]
                'icu_turnaround_time':  Time to process sample and report results / declare detection. 
                                        Formatted in a 2-element vector: [min,max]
                'gp_practices_total':   Total number of GP practices
                'gp_practices_swab':    Number of GP practices taking swabs for virology surveillance
                'gp_swabs_mg':          Assumed number of swabs that are metagenomic sequenced for investigating impact
                'pop_eng':              Population of England
                'gp_ari_consults':      Number of ARI consultations per 100k of England population per week. 
                                        Formatted in a 2-element vector: [mean summer 2024, mean winter 2024/25]
                'gp_ari_swabs':         Number of swabs taken from suspected ARI per week. 
                                        Formatted in a 2-element vector: [mean summer 2024, mean winter 2024/25]
                'pc_swab_turnaround_time':  Number of days between swab sample being taken and results received. 
                                            Formatted in a 2-element vector: [min,max] 

# Returns:  'sim_tds'   A dataframe containing the times to detection of 1 case and 3 cases for each simulation replicate
                        assuming either:
                                - only sampling ICU admissions
                                - only sampling primary care (eg GP) cases
                                - sampling from both settings
                        The output is of this format:
                        sim_n	ICU_TD	    PC_TD	ICU_PC_TD	ICU_3TD	    PC_3TD	ICU_PC_3TD	n_ICU_cases	n_ICU_cases_sampled	n_GP_cases	n_GP_cases_sampled	ICU_simid	                                GP_simid
                        1	    53.83937265	Inf	    53.83937265	60.83018273	Inf	    60.83018273	45	        5	                621	        0	                ["f6f6b432-8e44-11f0-3306-2724694310b1"]	["f6f6b432-8e44-11f0-3306-2724694310b1", "02e49496-8e45-11f0-256a-598e6634981b"]
                        2	    Inf	        Inf	    Inf	        Inf	        Inf	    Inf	        0	        0	                0	        0	                Inf	                                        Inf
                        ...     ...         ...     ...         ...         ...     ...         ...         ...                 ...         ...                 ...                                         ...

# Examples: td_results_test_1 = icu_v_pc_td(; sim_object_type = "full", file_or_object = "file"
                                            , sims_file = "sims.jdl2"
                                            , sims_object_name = "sims"
                                            , gp_swabs_mg = 100 , gp_ari_swabs = 319, gp_ari_consults = 180 )

            td_results_test_2 = icu_v_pc_td(; sim_object_type = "full", file_or_object = "object"
                                            , sims_file = "sims.jdl2", sims_object_name = "sims"
                                            , gp_swabs_mg = 100 , gp_ari_swabs = 319, gp_ari_consults = 180 )
    
            td_results_test_3 = icu_v_pc_td(; sim_object_type = "filtered", file_or_object = "file"
                                            , sims_G_gp_filter_file = "sims_G_gp_filter.jld2", sims_G_gp_filter_object_name = "sims_G_gp_filter" 
                                            , sims_G_icu_filter_file = "sims_G_icu_filter.jld2", sims_G_icu_filter_object_name = "sims_G_icu_filter" 
                                            , gp_swabs_mg = 100 , gp_ari_swabs = 319, gp_ari_consults = 180 )
    
            td_results_test_4 = icu_v_pc_td(; sim_object_type = "filtered", file_or_object = "object"
                                            , sims_G_gp_filter_object_name = "sims_G_gp_filter" 
                                            , sims_G_icu_filter_object_name = "sims_G_icu_filter" 
                                            , gp_swabs_mg = 100 , gp_ari_swabs = 319, gp_ari_consults = 180 )
"""
# icu_vs_pc_td was used in a previous analysis but has been superseded by icu_td and gp_td
function icu_v_pc_td(  ;  p=NBPMscape.P
                        , sim_object_type # = "full" #"G_filtered_icu_gp_combined" "G_filtered_icu_gp_separate"
                        , file_or_object # "file" or "object"
                        , sims_file = "" # "sims.jdl2"
                        , sims_G_gp_filter_file = "" # e.g. "sims_G_gp_filter.jdl2"
                        , sims_G_icu_filter_file = "" # e.g. "sims_G_icu_filter.jdl2"
                        , sims_object_name = "" # "sims" 
                        , sims_G_gp_filter_object_name = "" # e.g. "sims_G_gp_filter" 
                        , sims_G_icu_filter_object_name = "" # e.g. "sims_G_icu_filter" 
                        # ICU parameters
                        , icu_sample_type = NBPMscape.P.icu_sample_type # "regional" # "regional" or "fixed". If "fixed" then p_icu will be used, note that this doesn't take into account test sensitivity or practical sampling proportion, which "regional" does by using {sample_icu_cases} function
                        , pathogen_type = NBPMscape.P.pathogen_type #"virus"
                        , site_stage = NBPMscape.P.icu_site_stage #"current" 
                        , p_icu = NBPMscape.P.p_sampled_icu #0.15 # ICU sampling proportion
                        , only_sample_before_death = NBPMscape.P.icu_only_sample_before_death #true
                        , icu_ari_admissions = NBPMscape.P.icu_ari_admissions # 793 # 1440 # Weekly ICU admission numbers [summer,winter]
                        , icu_turnaround_time = NBPMscape.P.turnaroundtime_icu #[2,4] # Time to process sample and report results / declare detection
                       # Parameters for existing Oxford-RCGP RSC primary care surveillance
                        , gp_practices_total = NBPMscape.P.gp_practices_total #6199 # Total number of GP practices in England at July 2025. Source: BMA analysis (https://www.bma.org.uk/advice-and-support/nhs-delivery-and-workforce/pressures/pressures-in-general-practice-data-analysis) of NHS Digital General Practice Workforce Statistics (https://digital.nhs.uk/data-and-information/publications/statistical/general-and-personal-medical-services) [Accessed 2 Sep 2025]  
                        , gp_practices_swab = NBPMscape.P.gp_practices_swab #300 # Number of GP practices taking swabs for virology surveillance. Source: Data quality report: national flu and COVID-19 surveillance report (27 May 2025)
                        , gp_swabs_mg = NBPMscape.P.gp_swabs_mg #100 # 200 [319, 747]#[25698,46685] # Assumed number of swabs that are metagenomic sequenced for investigating impact
                        , pop_eng = NBPMscape.P.pop_eng #5.7106398e7 # Population of England. Source: ONS mid-year 2022 UK population data disaggregated by various geo levels and age groups and gender. [Accessed 6 November 2024] Available at https://www.ons.gov.uk/peoplepopulationandcommunity/populationandmigration/populationestimates/datasets/populationestimatesforukenglandandwalesscotlandandnorthernireland 
                        , gp_ari_consults = NBPMscape.P.gp_ari_consults #180 # 327 # Number of ARI consultations per 100k of England population per week [mean summer 2024, mean winter 2024/25]. Source: Analysis of data extracted from RCGP Research & Surveillance Centre (RSC) Virology Dashboard [Accessed 29 Aug 2025]
                        , gp_ari_swabs = NBPMscape.P.gp_ari_swabs #319 # 747] #[25698,46685]# Number of swabs taken from suspected ARI per week [mean summer 2024, mean winter 2024/25]. Source: Analysis of data extracted from RCGP Research & Surveillance Centre (RSC) Virology Dashboard [Accessed 29 Aug 2025]
                        , pc_swab_turnaround_time = NBPMscape.P.turnaroundtime_rcgp # [2,4] # [min,max] number of days between swab sample being taken and results received. Source: Data quality report: national flu and COVID-19 surveillance report (27 May 2025).  Assume same for metagenomic testing.
                     )
    
    ### Probability of primary care surveillance events
    # Probability of visiting a General Practice that takes surveillance swabs
    prob_visiting_gp_swab = minimum( [1, gp_practices_swab / gp_practices_total ])
    # Probability of being swabbed given visited a GP that takes suriveillance swabs
    prob_swabbed_g_visit_gp_swab = minimum([1, (1e5*gp_ari_swabs) / (pop_eng*gp_ari_consults) / prob_visiting_gp_swab])
    # Probability that swab is selected for metagenomic sequencing given that patient has been selected for a surveillance swab
    prob_mg_g_swabbed = minimum([1, gp_swabs_mg / gp_ari_swabs ])
    # Primary care sampling proportion - equivalent to probability a swab being taken and a metagenomic sequence being taken if attend a GP
    p_pc_mg_g_gp_visit = prob_visiting_gp_swab * prob_swabbed_g_visit_gp_swab * prob_mg_g_swabbed
    
    # Check if primary care sampling proportion, p_pc, was entered as an input to the function
    p_pc = p_pc_mg_g_gp_visit
    
    ### Load file if required...
    if file_or_object == "file"
        if sim_object_type == "full" || sim_object_type == "G_filtered_icu_gp_combined"
            sims = load(sims_file, sims_object_name)
        elseif sim_object_type == "G_filtered_icu_gp_separate"
            sims_G_icu_filter = load(sims_G_icu_filter_file, sims_G_icu_filter_object_name)
            sims_G_gp_filter  = load(sims_G_gp_filter_file,  sims_G_gp_filter_object_name)
        end
    elseif file_or_object == "object" # ... or just assign objects to variable names if already loaded
        if sim_object_type == "full" || sim_object_type == "G_filtered_icu_gp_combined"
            isdefined(Main, :sims_object_name) ? sims = getfield(Main, Symbol(sims_object_name)) : nothing
        elseif sim_object_type == "G_filtered_icu_gp_separate"
            isdefined(Main, :sims_G_icu_filter_object_name) ? sims_G_icu_filter = getfield(Main, Symbol(sims_G_icu_filter_object_name)) : nothing
            isdefined(Main, :sims_G_gp_filter_object_name)  ? sims_G_gp_filter  = getfield(Main, Symbol(sims_G_gp_filter_object_name))  : nothing
        end
    end

    # How many replicates are in this simulation
    if sim_object_type == "full" || sim_object_type == "G_filtered_icu_gp_combined"
        n_replicates = length(sims)
    elseif sim_object_type == "G_filtered_icu_gp_separate"
        # Check number of replicates is the same for ICU and GP filtered datasets
        if length(sims_G_gp_filter) != length(sims_G_icu_filter)
            error("Different numbers of sim reps: sims_G_filtered_gp (=$(length(sims_G_gp_filter))) sims_G_filtered_icu (=$(length(sims_G_icu_filter)))")
        else
            n_replicates = length(sims_G_gp_filter)
        end
    end

    # Create df to store results
    col_any = Vector{Any}(undef, n_replicates)
    fill!(col_any, Inf)
    sim_tds_cols = [copy(col_any) for _ in 1:13]
    sim_tds = DataFrame( sim_tds_cols, ["sim_n"
                                        # Time to detection for 1 case
                                        ,"ICU_TD"     # ICU sampling only
                                        ,"PC_TD"      # primary care sampling only
                                        ,"ICU_PC_TD"  # ICU AND primary care sampling
                                        # Time to detection for 3 cases
                                        ,"ICU_3TD" # ICU sampling only
                                        ,"PC_3TD"  # primary care sampling only 
                                        ,"ICU_PC_3TD"  # ICU AND primary care sampling 
                                        ,"n_ICU_cases","n_ICU_cases_sampled"
                                        ,"n_GP_cases","n_GP_cases_sampled"
                                        ,"ICU_simid","GP_simid"
                                        ])

    # Loop through simulation replicates, sampling infection cases (adapted from sampleforest() function in 'core.jl' and 'median_TD_by_region.jl'),
    # and computing times to detection
    for s in 1:n_replicates
        #TEST
        #println("$(s) of $(n_replicates) replicates")

        # Add simulation number to results df
        sim_tds[s,:sim_n] = s
        
        #s=1
        # If sim_object_type is not pre-filtered then need to filter for ICU and GP cases
        if sim_object_type == "full"
            fo = sims[s]
            if size(fo,1) == 0
                icu_cases = DataFrame() #nothing
                gp_cases = DataFrame() #nothing
            else
                # Filter for ICU cases 
                icu_cases = fo.G[ isfinite.(fo.G.ticu), : ]
                # Filter for GP cases (some of these will also become ICU cases)
                gp_cases = fo.G[ isfinite.(fo.G.tgp), : ]
            end
        elseif sim_object_type == "G_filtered_icu_gp_combined"
            fo = sims[s]
            if size(fo,1) == 0
                icu_cases = DataFrame() #nothing
                gp_cases = DataFrame() #nothing
            else
                # Filter for ICU cases 
                icu_cases = fo[ isfinite.(fo.ticu), : ]
                # Filter for GP cases (some of these will also become ICU cases)
                gp_cases = fo[ isfinite.(fo.tgp), : ]
            end
        elseif sim_object_type == "G_filtered_icu_gp_separate"
            icu_cases = sims_G_icu_filter[s]
            gp_cases = sims_G_gp_filter[s]
            if size(icu_cases,1) == 0
                icu_cases = DataFrame() #nothing
            end
            if size(gp_cases,1) == 0
                gp_cases = DataFrame() #nothing
            end
            # Double check that icu_cases and gp_cases are from the same simulation.
            # Note that simid is generated in simtree() function and is unique to 
            # the cases originating from a particular imported case and so if a 
            # simulation has more than one import (i.e. using simforest() function) 
            # then a single outbreak simulation can contain multiple unique simid values.
            # However, it would be likely that the ICU and GP cases datasets would
            # have at least one simid value in common, although not certain.
            if size(icu_cases,1) > 0 && size(gp_cases,1) > 0
                simids_icu = unique( icu_cases.simid )  #unique( [icu_cases[1,:simid]] ) 
                simids_gp = unique( gp_cases.simid )
                sim_tds[s,:ICU_simid] = simids_icu #unique( [icu_cases[1,:simid] ]) 
                sim_tds[s,:GP_simid] = simids_gp #unique( [gp_cases[1,:simid]] )
                #if icu_cases[1,:simid] != gp_cases[1,:simid]
                if isempty( intersect( simids_icu, simids_gp) )
                    #error("Sim IDs don't match between sims_G_filtered_gp ($(gp_cases[1,:simid])) and sims_G_filtered_icu ($(icu_cases[1,:simid]))for sim  rep $(s)")
                    println("Sim IDs do not match between sims_G_filtered_gp and sims_G_filtered_icu for sim rep $(s)")
                    sim_tds[s,:sim_n] = s
                    sim_tds[s,2:11] = ["No Sim ID match" for _ in 2:11]
                    continue
                end
            end
        end

        ### Record infection time stats
        ## First for ICU cases
        icu_cases_sub = sample_icu_cases(; icu_cases = icu_cases
                                         , icu_sample_type = icu_sample_type
                                         , icu_fixed_sample_prop = p_icu
                                         , pathogen_type = "virus"
                                         , site_stage = "current"  )

        # Record number of cases in the sim rep and the number that were tested
        sim_tds[s,:n_ICU_cases] = size(icu_cases,1)
        sim_tds[s,:n_ICU_cases_sampled] = size(icu_cases_sub,1) #n_icu

        if size(icu_cases_sub,1) == 0
            
            # Record zero cases and zero sampled
            #sim_tds[s,:n_ICU_cases] = 0
            #sim_tds[s,:n_ICU_cases_sampled] = 0

            # Record Inf as time to detection for earliest:
            # - 1 ICU case
            sim_tds[s,:ICU_TD] = Inf
            # - 3 ICU cases
            sim_tds[s,:ICU_3TD] = Inf

            # Define empty vector for top 3 times to detection - to be used later
            icu_cases_sub_top3_td = []

        elseif size(icu_cases_sub,1) > 0 
            
            # Generate sample times
            if only_sample_before_death == true
                # Set upper limit on sample time to the minimum of time of death and ICU admission plus the user defined parameter p.icu_swab_lag_max
                icu_tsample = map( g -> rand( Uniform( g.ticu[1], min( g.tdeceased[1], g.ticu[1]+ p.icu_swab_lag_max ) ) ) , eachrow(icu_cases_sub) )
            else
                # Set upper limit on sample time to ICU admission plus the user defined parameter p.icu_swab_lag_max
                icu_tsample = map( g -> rand( Uniform( g.ticu[1], g.ticu[1]+ p.icu_swab_lag_max)) , eachrow(icu_cases_sub) )
            end
            icu_cases_sub.tsample = icu_tsample 

            # Simulate reports times
            #treport = (size(g_region,1)>0) ? (tsample .+ turnaroundtime) : []
            icu_treport = (icu_tsample .+ rand(Uniform(icu_turnaround_time[1], icu_turnaround_time[2])) )
            icu_cases_sub.treport = icu_treport
            #println(g_region_sub[:,[2,3,4,5,6,16,17]])

            # Find cases with 3 shortest times to detection (TD)
            icu_cases_sub_top3_td = first(sort(icu_cases_sub, :treport), 3)

            ## Add to results df for ICU times to detection
            # Record the times to detection for earliest:
            # - 1 ICU case
            sim_tds[s,:ICU_TD] = icu_cases_sub_top3_td[1,:treport]
            # - 3 ICU cases
            if size(icu_cases_sub_top3_td,1) >= 3
                sim_tds[s,:ICU_3TD] = icu_cases_sub_top3_td[3,:treport]
            else
                sim_tds[s,:ICU_3TD] = Inf
            end
        end

        # Record time stats for GP cases
        
        # Also apply metagenomic test sensitivity so that the subsampled GP cases are those with positive test results
        # Obtain metagenomic test sensitivity for pathogen_type
        sensitivities = Dict(
                            "virus" => p.sensitivity_mg_virus,
                            "bacteria" => p.sensitivity_mg_bacteria,
                            "fungi" => p.sensitivity_mg_fungi
                            )
        mg_test_sensitivity = get(sensitivities, pathogen_type) do
            error("Unknown pathogen type: $pathogen_type")
        end

        # Sample sizes based on probability of infected case having a metagenomic sample taken at the GP
        n_gp = rand( Binomial( size(gp_cases,1), p_pc * mg_test_sensitivity) )
        #println(n_gp)
            
        # Subsample of GP cases
        
        gp_cases_sub = gp_cases[sample( 1:size(gp_cases,1), n_gp, replace=false ), :]
        
        # Record number of cases in the sim rep and the number that were tested
        sim_tds[s,:n_GP_cases] = size(gp_cases,1)
        sim_tds[s,:n_GP_cases_sampled] = n_gp

        # If there are NO GP cases in the sample
        if size(gp_cases_sub,1) == 0
                
            # Record zero cases and zero sampled
            #sim_tds[s,:n_GP_cases] = 0
            #sim_tds[s,:n_GP_cases_sampled] = 0

            ## Add to results df for ICU times to detection
            # Record Inf times to detection if no GP cases sampled:
            #  1 primary care case for TD
            sim_tds[s, :PC_TD ] = Inf
            #  3 primary care cases for 3TD
            sim_tds[s, :PC_3TD ] = Inf

            # ICU and primary care combined
            # If there are no GP cases then the combination of primary care and ICU is just ICU only (which was already entered into df above)
            # - 1 case
            sim_tds[s, :ICU_PC_TD] = sim_tds[s,:ICU_TD]
            # - 3 case
            sim_tds[s, :ICU_PC_3TD] = sim_tds[s,:ICU_3TD]
                
            elseif size(gp_cases_sub,1) > 0 # If there ARE GP cases in the sample
            
                # Generate sample times
                gp_tsample = gp_cases_sub.tgp
                gp_cases_sub.tsample = gp_tsample # primary care sampling is done at the time of GP visit (ignores home sampling)

                # Simulate reports times
                gp_treport = (gp_tsample .+ rand(Uniform(pc_swab_turnaround_time[1], pc_swab_turnaround_time[2])) )
                gp_cases_sub.treport = gp_treport

                ## Find cases with 3 shortest times to detection (TD)
                # Primary care cases only
                gp_cases_sub_top3_td = first(sort(gp_cases_sub, :treport), 3)
                
                # ICU and primary care combined
                #if @isdefined icu_cases_sub_top3_td # Check if there are any ICU cases...
                if size(icu_cases_sub_top3_td, 1) > 0 # Check if there are any ICU cases...
                    icu_gp_cases_sub_top3_td = sort( append!(icu_cases_sub_top3_td, gp_cases_sub_top3_td), :treport)
                else # ... if not then the combination is just primary care (GP) cases only
                    icu_gp_cases_sub_top3_td = gp_cases_sub_top3_td
                end

                ## Add to results df for primary care times to detection
                # Record the times to detection for earliest:
                # Primary care cases only
                #  1 primary care case for TD
                sim_tds[s, :PC_TD ] = gp_cases_sub_top3_td[1,:treport]
                #  3 primary care cases for 3TD
                if size(gp_cases_sub_top3_td,1) >= 3
                    sim_tds[s, :PC_3TD ] = gp_cases_sub_top3_td[3,:treport]
                else
                    sim_tds[s, :PC_3TD ] = Inf 
                end
                
                # Primary care COMBINED with ICU cases
                # - 1 case, either ICU or primary care
                sim_tds[s, :ICU_PC_TD] = icu_gp_cases_sub_top3_td[1,:treport]
                # - 3 case, either ICU or primary care
                if size(icu_gp_cases_sub_top3_td,1) >= 3
                    sim_tds[s, :ICU_PC_3TD] = icu_gp_cases_sub_top3_td[3,:treport]
                else
                    sim_tds[s, :ICU_PC_3TD] = Inf
                end

            end # End of if loop checking whether there are any GP cases sampled

    end # End of loop through different simulations
    return sim_tds
end # End of icu_v_pc_td function



"""
Function: analyse_td_columns

Description:    To be used in conjunction with {icu_vs_pc_td}, namely to process its output
                Function to compute:
                - median values for each time to detection scenario across multiple simulation replicates
                - % of simulations that return a time to detection in each time to detection scenario

Arguments:    'df':   A dataframe as output from the {icu_v_pc_td} function containing the times to detection of
                        1 case and 3 cases for each simulation replicate assuming either:
                        - only sampling ICU admissions
                        - only sampling primary care (eg GP) cases
                        - sampling from both settings
                        The dataframe is of this format:
                          Row | sim_n  ICU_TD   PC_TD    ICU_PC_TD  ICU_3TD  PC_3TD  ICU_PC_3TD  n_ICU_cases  n_ICU_cases_sampled  n_GP_cases  n_GP_cases_sampled  ICU_simid  GP_simid 
                              │ Any    Any      Any      Any        Any      Any     Any         Any          Any                  Any         Any                 Any        Any
                        ──────┼────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
                            1 │ 1      Inf      Inf      Inf        Inf      Inf     Inf         0            0                    1           0                   Inf        Inf
                            2 │ 2      Inf      Inf      Inf        Inf      Inf     Inf         0            0                    8           0                   Inf        Inf
                          ⋮   │   ⋮       ⋮        ⋮         ⋮         ⋮       ⋮         ⋮            ⋮                ⋮               ⋮               ⋮               ⋮         ⋮
                          999 │ 999    49.082   69.3822  49.082     55.1244  Inf     55.1244     104          17                   1096        1                   Inf        Inf
                         1000 │ 1000   Inf      Inf      Inf        Inf      Inf     Inf         0            0                    0           0                   Inf        Inf 
            
                Ideally this should be trimmed to only include the columns related to time to detected, i.e. columns 2 to 11
            
                        Row │  ICU_TD   PC_TD    ICU_PC_TD  ICU_3TD  PC_3TD  ICU_PC_3TD  n_ICU_cases  n_ICU_cases_sampled  n_GP_cases  n_GP_cases_sampled 
                            │    Any      Any      Any        Any      Any     Any         Any          Any                  Any         Any             
                        ──────┼────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
                            1 │      Inf      Inf      Inf        Inf      Inf     Inf         0            0                    1           0                   
                            2 │      Inf      Inf      Inf        Inf      Inf     Inf         0            0                    8           0                   
                          ⋮   │       ⋮        ⋮         ⋮         ⋮       ⋮         ⋮            ⋮                ⋮               ⋮               ⋮          
                          999 │    49.082   69.3822  49.082     55.1244  Inf     55.1244     104          17                   1096        1                   
                         1000 │   Inf      Inf      Inf        Inf      Inf     Inf         0            0                    0           0                   

Returns:  'result':   A dataframe containing three columns: description of parameter (mostly relating to time to detection), 
                      median value (mostly time to detection), and the percentage of simulation replicates that have a (time to) detection.
                      The output is of this format:

                        Row │ TD_description       Median_TD  Percentage_with_a_TD 
                            │ String               Float64    Float64
                        ─────┼──────────────────────────────────────────────────────
                        1 │ ICU_TD                 54.3793                  62.1
                        2 │ PC_TD                  61.0425                  25.4
                        3 │ ICU_PC_TD              53.93                    63.2
                        4 │ ICU_3TD                61.0959                  42.6
                        5 │ PC_3TD                 64.1294                   3.9
                        6 │ ICU_PC_3TD             60.8523                  44.0
                        7 │ n_ICU_cases            12.0                    100.0
                        8 │ n_ICU_cases_sampled     2.0                    100.0
                        9 │ n_GP_cases            143.0                    100.0
                        10 │ n_GP_cases_sampled      0.0                    100.0

Examples:   td_results_test_1_analysis = analyse_td_columns(td_results_test_1[:, 2:11])
            println( td_results_test_1_analysis )
    
"""
function analyse_td_columns(df::DataFrame)
    result = DataFrame(TD_description=String[], Median_TD=Float64[], Percentage_with_a_TD=Float64[])
    # How many rows (sim reps) have a Sim ID mismatch
    #n_simid_mismatch = sum(col_data .=="No Sim ID match")
    #println("$(n_simid_mismatch) have a Sim ID mismatch")
    for col in names(df)
        col_data = df[!, col]
        
        # Remove rows (sim reps) with a sim id mismatch
        col_data_filtered = filter(row -> row != "No Sim ID match", col_data)
        
        # Convert strings "Inf" to actual Inf (if present as string)
        #cleaned = [x in ["Inf", "Sim ID mismatch"] ? Inf : x for x in col_data]
        #cleaned = [x == "Inf" ? Inf : x for x in col_data_filtered]
        
        # Keep only finite numbers
        finite_vals = filter(x -> isfinite(x), col_data_filtered) #cleaned)
        percent_finite_ex_simid_mismatch = 100 * length(finite_vals) / length(col_data_filtered)
        median_val = isempty(finite_vals) ? NaN : median(skipmissing(finite_vals))
        push!(result, (col, median_val, percent_finite_ex_simid_mismatch))
    end
    return result
end


# used in previous analyses but superseded
"""
Function: clean_td_results
    
Description:    Function to remove rows (simulation replicates) containing "No Sim ID match" and 
                convert string values of Inf to numeric values.

Arguments:  'df':   A dataframe as output from the {icu_v_pc_td} function containing the times to detection of
                    1 case and 3 cases for each simulation replicate assuming either:
                    - only sampling ICU admissions
                    - only sampling primary care (eg GP) cases
                    - sampling from both settings
                    The dataframe is of this format:
                    Row │ sim_n  ICU_TD   PC_TD    ICU_PC_TD  ICU_3TD  PC_3TD  ICU_PC_3TD  n_ICU_cases  n_ICU_cases_sampled  n_GP_cases  n_GP_cases_sampled  ICU_simid  GP_simid 
                        │ Any    Any      Any      Any        Any      Any     Any         Any          Any                  Any         Any                 Any        Any
                    ──────┼────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
                        1 │ 1      Inf      Inf      Inf        Inf      Inf     Inf         0            0                    1           0                   Inf        Inf
                        2 │ 2      Inf      Inf      Inf        Inf      Inf     Inf         0            0                    8           0                   Inf        Inf
                    ⋮   │   ⋮       ⋮        ⋮         ⋮         ⋮       ⋮         ⋮            ⋮                ⋮               ⋮               ⋮               ⋮         ⋮
                    999 │ 999    49.082   69.3822  49.082     55.1244  Inf     55.1244     104          17                   1096        1                   Inf        Inf
                    1000 │ 1000   Inf      Inf      Inf        Inf      Inf     Inf         0            0                    0           0                   Inf        Inf 
            
Examples:   td_results_test_1_clean = clean_TD_results(td_results_test_1)#[:, 2:11])
            println( td_results_test_1_analysis )
"""
# used in previous analyses but superseded

function clean_td_results(df::DataFrame)
    result = df
    for col in names(df)
        col_data = df[!, col]
        # Remove rows (sim reps) with a sim id mismatch
        col_data_filtered = filter(row -> row != "No Sim ID match", col_data)
        # Convert strings "Inf" to actual Inf (if present as string)
        cleaned = [x == "Inf" ? Inf : x for x in col_data_filtered]
        # Keep only finite numbers
        finite_vals = filter(x -> isfinite(x), cleaned)
        
    end
    return result
end



"""
Function: plot_hist

Description:    Function to plot a histogram of a particular column a dataframe.

Arguments:  'df':   A dataframe as output from the {icu_v_pc_td} function containing the times to detection of
                    1 case and 3 cases for each simulation replicate assuming either:
                    - only sampling ICU admissions
                    - only sampling primary care (eg GP) cases
                    - sampling from both settings
                    The dataframe is of this format:
                    Row │ sim_n  ICU_TD   PC_TD    ICU_PC_TD  ICU_3TD  PC_3TD  ICU_PC_3TD  n_ICU_cases  n_ICU_cases_sampled  n_GP_cases  n_GP_cases_sampled  ICU_simid  GP_simid 
                        │ Any    Any      Any      Any        Any      Any     Any         Any          Any                  Any         Any                 Any        Any
                    ──────┼────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
                        1 │ 1      Inf      Inf      Inf        Inf      Inf     Inf         0            0                    1           0                   Inf        Inf
                        2 │ 2      Inf      Inf      Inf        Inf      Inf     Inf         0            0                    8           0                   Inf        Inf
                      ⋮   │   ⋮       ⋮        ⋮         ⋮         ⋮       ⋮         ⋮            ⋮                ⋮               ⋮               ⋮               ⋮         ⋮
                      999 │ 999    49.082   69.3822  49.082     55.1244  Inf     55.1244     104          17                   1096        1                   Inf        Inf
                     1000 │ 1000   Inf      Inf      Inf        Inf      Inf     Inf         0            0                    0           0                   Inf        Inf 
                    
                    Or output from {clean_TD_results}

Returns:    Plots histogram to screen

Examples:   plot_hist( df = td_results_test_1, col = 2)
    
"""
function plot_hist(;df::DataFrame, col)
        col_data = df[!, col]
        # Remove rows (sim reps) with a sim id mismatch
        col_data_filtered = filter(row -> row != "No Sim ID match", col_data)
        
        # Convert strings "Inf" to actual Inf (if present as string)
        #cleaned = [x in ["Inf", "Sim ID mismatch"] ? Inf : x for x in col_data]
        cleaned = [x == "Inf" ? Inf : x for x in col_data_filtered]
        
        # Keep only finite numbers
        #finite_vals = filter(x -> isfinite(x), col_data_filtered) #cleaned)
        finite_vals = filter(x -> isfinite(x), cleaned)
        histogram([finite_vals],label = names(df)[col])
        #println(names(df)[col])
end