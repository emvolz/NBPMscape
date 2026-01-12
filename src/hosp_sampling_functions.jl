
#=
Functions included here:

- courier_collection_days   Generates vector of times when swabs will be collected from public health laboratories
- sample_hosp_cases_n       

=#

"""
Function        courier_collection_days
Description     Generates vector of times when swabs will be collected from public health laboratories

Arguments   initial_dow::Integer                Day of the week at the start of the simulated outbreak, 
                                                i.e. when the first infected individual arrived in the 
            n_weeks::Integer                    Number of weeks of outbreak to generate collection days
            collection_dow::Vector{Integer}     Day of the week that swabs will be collected by courier, e.g. [2,5] or [2,4,6].
                                                Day codes: Sunday=1, Monday=2, ..., Saturday=7
                                            
Returns     Vector of integer times (days)
Example     max_week = 18
            phl_courier_collection_days = courier_collection_times(; initial_dow = 1, n_weeks = max_week, collection_dow = [2,4,6])
"""
function courier_collection_times(; initial_dow::Integer, n_weeks::Integer, collection_dow::Vector{Int64})
    # @assert initial_dow >= 1 && initial_dow <= 7
    # @assert n_weeks > 0
    times = []
    for i in 1:n_weeks #i=3
        for j in 1:length(collection_dow)
            # Collection time = collection day (x week number) - initial day + 0.5 for midday collection
            times_temp = collection_dow[j] +(7*(i-1)) - initial_dow + 0.5
            push!(times,times_temp)
        end
    end

    # Only return positive values for collection times, i.e. ignore the collection days of the week that are before the initial day of the week
    return( filter(>(0), times) )
end



"""
Function:   sample_hosp_cases_n

Description:    Samples cases based on a number of total target number of samples split across a number of local labs which are fed by individual hospitals.
                The number of samples taken from each local can be distributed equally or weighted by either catchment population or A&E attendances.
                Samples at local labs are then transferred to a central lab for metagenomic testing.
                Sampling Criteria:
                - Only sample from swabs taken within 48h of admission/presentation to hospital
                - Assume Gamma distribution for swabbing times in hospital
                - Samples transferred from hospital to local public health laboratory 
                - Courier collections n times per week from local PHL and transferred to central metagenomic testing lab - most recent m samples 
                Estimate number of background ARI hospital cases:
                - by age (ED data)
                - by destination after ED (ED data) 
                - at each HARISS participating site / NHS Trust
                - allocate cases pro rata using hospital catchment area populations or average A&E attendances
                Summary of code:
                1 Run simulation and save all GP, hospital, ICU cases (this will be a larger file so may need to trim down to the data we actually need)
                2 For HARISS, filter for hospital cases (may need to add levels of severity to care pathways)
                3 Allocate each case to an NHS Trust based on home region and probability of emergency admission.
                4 Filter for cases at Trusts that are part of the sampling network (e.g. HARISS)
                5 Simulate hospital arrival times for background ARI cases using Uniform distribution
                6 Add Gamma distributed time between hospital arrival and swabbing to get tswab (for both ARI background and pathogen X cases).
                  Estimated distributed is truncated depending on length of hospital stay.
                7 Merge linelists for ARI background and pathogen X cases
                6 Filter for cases sampled:
                (a) within 48h of presentation to hospital
                (b) between 24h and 48h of presentation to hospital
                6 Define courier pick up times (e.g. Monday and Thursday) and swab cut-off times to account for transfer between hospital and public health lab
                7 Filter for X/(number of weekly courier pickups) most recent samples (merging background ARI and pathogen X samples) prior to the courier 
                  pick-up cut-off times (i.e. filter eligible samples, sort by time of swabbing and select most recent 
                  (weekly sample target)/(number of weekly courier pickups) samples for each pick up day)
                  - option to split by adults and children, e.g. A adult samples and C child samples (assume hospitals with an ED also able to treat children,
                    albeit they may be transferred if admitted, but this is not included)
                8 Compute probabilities of these samples being background or pathogen X

Arguments:  p::NamedTuple                       (default: NBPMscape.P)
            hosp_cases::DataFrame   Version of the G dataframe output from {simforest} or {simtree} and filtered for hospital  
                                    cases, i.e. rows with a value for ted and/or thospital. For example:
                                    
                                    194603×18 DataFrame
                                        Row │ pid                                tinf      tgp       ted       thospital  ticu     tstepdown  tdischarge  trecovered  tdeceased  severity                fatal  iscommuter  homeregion  infector_age  infectee_age  importedinfection  simid                             
                                            │ String                             Float64   Float64   Float64   Float64    Float64  Float64    Float64     Float64     Float64    Symbol                  Bool   Bool        String      Int8?         Int8          Bool               String
                                    ────────┼────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
                                          1 │ b7ceb040-d6a6-11f0-0449-1f31a343…   9.32135   11.2064   11.9129   Inf           Inf        Inf     12.2683     28.2669        Inf  moderate_ED             false       false  TLH2                   4             8              false  b7ceb040-d6a6-11f0-0449-1f31a343…
                                          2 │ b7ceb040-d6a6-11f0-0449-1f31a343…  11.2938    16.0533   16.4357   Inf           Inf        Inf     16.5044     40.7698        Inf  moderate_ED             false        true  TLH2                   8            25              false  b7ceb040-d6a6-11f0-0449-1f31a343…
                                          3 │ b7ceb040-d6a6-11f0-0449-1f31a343…  14.9047    19.5912   23.0399   Inf           Inf        Inf     23.1218     40.2354        Inf  moderate_ED             false       false  TLH2                   5            66              false  b7ceb040-d6a6-11f0-0449-1f31a343…
                                          4 │ b7ceb040-d6a6-11f0-0449-1f31a343…  15.4063    25.3535  Inf         26.0261      Inf        Inf     26.67       36.3814        Inf  severe_hosp_short_stay  false       false  TLH2                   7            78              false  b7ceb040-d6a6-11f0-0449-1f31a343…
                                    ⋮       │                 ⋮                     ⋮         ⋮         ⋮          ⋮         ⋮         ⋮          ⋮           ⋮           ⋮                ⋮               ⋮        ⋮           ⋮            ⋮             ⋮                ⋮                          ⋮
                                    194600 │ cce99d42-d6aa-11f0-2949-b383acab…  85.1219   Inf        88.7939   Inf           Inf        Inf     89.2342    104.08          Inf  moderate_ED             false        true  TLI4                  30            57              false  cce99d42-d6aa-11f0-2949-b383acab…
                                    194601 │ cce99d42-d6aa-11f0-2949-b383acab…  89.8193    92.1103   92.1211   Inf           Inf        Inf     92.3438    117.742         Inf  moderate_ED             false       false  TLI4                  50            84              false  cce99d42-d6aa-11f0-2949-b383acab…
                                    194602 │ cce99d42-d6aa-11f0-2949-b383acab…  88.5175   Inf        94.1549   Inf           Inf        Inf     94.4465    127.32          Inf  moderate_ED             false        true  TLI4                  33            52              false  cce99d42-d6aa-11f0-2949-b383acab…
                                    194603 │ cce99d42-d6aa-11f0-2949-b383acab…  88.5591    98.4067   98.6414   Inf           Inf        Inf     99.1132    127.733         Inf  moderate_ED             false       false  TLI4                  33            11              false  cce99d42-d6aa-11f0-2949-b383acab…
            initial_dow::Int64                      Day of the week that initial case imported. This is used as a reference to map time since first importation
                                                    to days of the week, which is required for the courier sample collection timetable. 
                                                    Day of week codes: Sunday =1, Monday = 2, ..., Saturday = 7
            ## Sampling parameters ##
            n_hosp_samples_per_week::Int            Total number of hospital samples to be taken per week for metagenomic testing
            sample_allocation::String               Options are "equal" or "weighted", which refer to whether the total number of samples is 
                                                    split equally between public health laboratories or weighted by either catchment area population
                                                    or A&E (aka ED) attendances averaged over 12 months.
            weight_samples_by::String               Options are "ae_mean" or "catchment_pop". If 'sample_allocation' (described above) = "weighted", then the 
                                                    number of samples taken from each public health laboratory (PHLs) will be computed based on either the 
                                                    12-month mean proportion of A&E attendances at NHS Trusts that feed into the PHLs or the proportion of the
                                                    population that feeds into NHS Trusts and in turn the PHLs.
            sample_proportion_adult::Any            "free" indicates no specific weighting between adults and children, while a numeric value between 0 and 1
                                                    indicates the proportion of samples that should be adults (aged 16y and over), 
                                                    with sample_proportion_child = 1-sample_proportion_adult
            hariss_nhs_trust_sampling_sites::DataFrame  DataFrame containing NHS Trusts and public health laboratories in the hospital sampling network, e.g. HARISS network.
                                                        Example format:
                                                        Row │ NHS_Trust_code  NHS_Trust_name                      public_health_laboratory 
                                                            │ String3         String                             String15
                                                        ────┼─────────────────────────────────────────────────────────────────────────────
                                                        1 │ AA1             NHS Trust A1                        PHL A
                                                        2 │ BB1             NHS Trust B1                        PHL B
                                                        3 │ BB2             NHS Trust B2                        PHL B
                                                        4 │ CC1             NHS Trust CC1                       PHL C
            phl_collection_dow::Vector{Int64}       Day(s) of the week that swab samples will be collected from public health labs. 
                                                    Day of week codes: Sunday = 1,... Saturday = 7. An example value is [2,5] (Monday and Thursday).
                                                    Note that the collection time assumed is midday.
            swab_time_mode::Real                    Timing of swabbing of ARI patients in hospital is assumed to follow a Gamma distribution, with a peak
                                                    shortly after arrival and a certain percentage of patients swabbed within 2 days (= 48hrs) of arrival.
                                                    This argument defines the timing of the peak in days, for example 0.25 days (=6hrs) after 
                                                    attendance/admission at hospital.
            swab_proportion_at_48h::Real            See description of 'swab_time_mode' above. 'swab_proportion_at_48h' defines the proportion of swabs that
                                                    are taken within 2 days (=48hrs) of attendance/admission at hospital.
            proportion_hosp_swabbed::Real           For practical reasons, all ARI attendances at hospital may not be sampled. This argument sets the 
                                                    proportion that are sampled using a value between 0 and 1.
            
            only_sample_before_death::Bool          There is a possibilty that the swabbing time (tswab) drawn from the Gamma distribution is after the time of 
                                                    death. So 'true' here will constrain tswab to tdeceased.
            ## Hospital parameters ##
            nhs_trust_catchment_pop::DataFrame      NHS Trust catchment area populations as a proportion of the total population. For example:
                                                    NHS_TRUST_CATCHMENT_POP_ADULT_CHILD
                                                    126×4 DataFrame
                                                    Row │ TrustCode  catchment_prop_of_total_sum  prop_child   prop_adult 
                                                        │ Any        Float64                      Float64      Float64
                                                    ────┼─────────────────────────────────────────────────────────────────
                                                      1 │ R0A                         0.021291    0.00634217   0.0149488
                                                      2 │ R0B                         0.00772294  0.00146184   0.0062611
                                                      3 │ R0D                         0.0104117   0.0018323    0.00857945
                                                      4 │ R1F                         0.00237302  0.000444012  0.001929
                                                     ⋮  │     ⋮                   ⋮                    ⋮           ⋮
                                                    123 │ RXR                         0.00872926  0.00204953   0.00667973
                                                    124 │ RXW                         0.00811911  0.00169838   0.00642073
                                                    125 │ RYJ                         0.0100881   0.00187858   0.00820951
                                                    126 │ RYR                         0.0179986   0.00352307   0.0144755
            nhs_trust_ae_12m                        12-month data for A&E attendances by NHS Trust, including the 12-month mean as a 
                                                    proportion of total A&E attendances. For example:
                                                    AE_12M
                                                    199×15 DataFrame
                                                    Row │ NHS_Trust_code  2024_4_function  2024_5_function  2024_6_function  2024_7_function  2024_8_function  2024_9_function  2024_10_function  2024_11_function  2024_12_function  2025_1_function  2025_2_function  2025_3_function  mean_12m  mean_12m_prop 
                                                        │ String7         Int64            Int64            Int64            Int64            Int64            Int64            Int64             Int64             Int64             Int64            Int64            Int64            Float64   Float64       
                                                    ────┼────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
                                                     12 │ NTV0W                         0                0                0                0                0                0                 0                 0                 0                0                0                0      0.0      0.0
                                                     13 │ RFF                        8335             8850             8780             9089             8490             8612              9029              9113              9184             8663             7870             9281   8774.67     0.00642199
                                                     14 │ RNU                           0                0                0                0                0                0                 0                 0                 0                0                0                0      0.0      0.0
                                                     15 │ RW1                           0                0                0                0                0                0                 0                 0                 0                0                0                0      0.0      0.0
                                                     16 │ RWY                       15028            16416            15680            15931            14746            14877             15686             15295             16032            14485            14096            16307  15381.6      0.0112574
                                                     17 │ RXL                        8434             9230             8724             9152             8668             8420              8970              8663              9075             8263             7854             8950   8700.25     0.00636752
            ed_discharge_limit::Float64             Upper limit on the time, in days, that people attending the Emergency Department (ED) spent there before being discharged.
                                                    Note that this only includes people discharged directly from ED, i.e. those not admitted to hospital.
                                                    Also, note that this should match the value used in {simtree} and/or {simforest}.
            hosp_short_stay_limit::Float64          Upper limit on the time, in days, that people admitted to hospital (i.e. not discharged directly from ED) spend
                                                    there before being discharged. Note that this does not includes people discharged directly from ED.
                                                    Also, note that this should match the value used in {simtree} and/or {simforest}.
            ## Seasonal values ##
            ed_ari_destinations_adult::DataFrame    Destinations for adulta (16y and over) after arrival at hospital and proportions. Note that can be seasonal. 
                                                    Example format:
                                                    DataFrame( destination = [:discharged,:short_stay,:longer_stay]
                                                             , proportion_of_attendances = [0.575,0.034,0.391])                                                    
            ed_ari_destinations_child::DataFrame    Destinations for children (<16y) after arrival at hospital and proportions. Note that can be seasonal.
                                                    Example format:
                                                    DataFrame( destination = [:discharged,:short_stay,:longer_stay]
                                                             , proportion_of_attendances = [0.851,0.012,0.137])
            hosp_ari_admissions::Int                Estimate of weekly hospital ARI admissions (excluding pathogen X being simulated)
            hosp_ari_admissions_adult_p::Float64    Estimate of proportion of ED ARI admissions that are adults (16y and over)
                                                    Note that this is a fixed value in the model but in reality it varies throughout the year.
            hosp_ari_admissions_child_p::Float64    Estimate of proportion of ED ARI admissions that are children (<16y)
                                                    Note that this is a fixed value in the model but in reality it varies throughout the year.
            
Returns:    DataFrame of the same format as hosp_cases (see above) but only contains the rows representing sampled hospital cases for 
            the simulated pathogen X

Examples:   
        
hosp_cases_sub = sample_hosp_cases_n(; p = NBPMscape.P
                                    , hosp_cases::DataFrame # filtered from simulation of infections (G df filtered for value in thosp column)
                                    , initial_dow::Int64 = 1 # Day of the week that initial case imported. Day of week codes: Sunday =1, Monday = 2, ..., Saturday = 7 
                                    # Sampling parameters
                                    , n_hosp_samples_per_week::Int # Total number of hospital samples to be taken per week
                                    , sample_allocation::String = "equal" # "equal" or "weighted"
                                    , sample_proportion_adult::Any = "free" # "free" or numeric decimal, e.g. 0.75. Indicates split of sample target 
                                                                            # between adults and children. "free" indicates that no split is specified
                                    , hariss_nhs_trust_sampling_sites::DataFrame # List of NHS Trusts in HARISS sampling network
                                                                    , swab_time_mode::Real = 0.25 # Assume swabbing peaks at 6hrs (=0.25 days) after attendance/admission at hospital
                                    , swab_proportion_at_48h::Real = 0.9 # Assume 90% of swabs are taken within 48hrs (=2 days) of attendance/admission at hospital
                                    , proportion_hosp_swabbed::Real = 0.9 # Assume 90% of ARI attendances are swabbed
                                    , only_sample_before_death::Bool = true # There is a possibilty of swabbing time being drawn after death so 'true' here will constrain tswab to tdeceased
                                    # Hospital parameters
                                    , ed_discharge_limit::Float64 = NBPMscape.P.tdischarge_ed_upper_limit # days. Assume that people attending the Emergency Department are discharged within this time limit.
                                    , hosp_short_stay_limit::Float64 = NBPMscape.P.tdischarge_hosp_short_stay_upper_limit # days. Assume that people attending the Emergency Department are discharged within this time limit.
                                    , nhs_trust_catchment_pop = NHS_TRUST_CATCHMENT_POP_ADULT_CHILD
                                    #, nhs_trust_ae_12m = AE_12M
                                    # Seasonal values
                                    , hosp_ari_admissions::Int # Estimate of weekly hospital ARI admissions (excluding pathogen X being simulated)
                                    , hosp_ari_admissions_adult_p::Float64 = 0.52# Proportion of ED ARI admissions that are adults (16y and over)
                                    , hosp_ari_admissions_child_p::Float64 = 0.48 # Proportion of ED ARI admissions that are children (<16y)
                                    , ed_ari_destinations_adult::DataFrame = DataFrame( destination = [:discharged,:short_stay,:longer_stay]
                                                                                        , proportion_of_attendances = [0.628,0.030,0.342]
                                                                                        )
                                    , ed_ari_destinations_child::DataFrame = DataFrame( destination = [:discharged,:short_stay,:longer_stay]
                                                                                        , proportion_of_attendances = [0.861,0.014,0.125]
                                                                                        )
                                    , weight_samples_by = "ae_mean" # or "catchment_pop"
                                    , phl_collection_dow::Vector{Int64} = [2,5] # Day(s) of week that swab samples will be collected from public health labs. Day of week codes: Sunday = 1,... Saturday = 7.
                                    )
   
    
                                    553×24 DataFrame
                                    Row │ pid                                tinf     tgp       ted       thospital  ticu     tstepdown  tdischarge  trecovered  tdeceased  severity                fatal  iscommuter  homeregion  infector_age  infectee_age  importedinfection  simid                              hosp_nhs_trust_cd  public_health_laboratory  ted_or_thospital  tswab     tswab_sub_48h  tcourier 
                                        │ String                             Float64  Float64   Float64   Float64    Float64  Float64    Float64     Float64     Float64    Symbol                  Bool   Bool        String      Int8?         Int8          Bool               String                             Any                String15                  Float64           Float64   Bool           Float64?
                                    ─────┼─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
                                    1 │ 6c780ade-d6b9-11f0-21b0-717e8d16…  17.7552   21.735   Inf         23.3375      Inf        Inf     24.0794     41.7767        Inf  severe_hosp_short_stay  false        true  TLI7                  30            63              false  6c780ade-d6b9-11f0-21b0-717e8d16…  RYJ                Imperial                           23.3375   23.6481           true      25.5
                                    2 │ 04c8bc98-d6ba-11f0-2861-0dc9995b…  31.3438   33.4322   34.6315   Inf           Inf        Inf     34.8906     45.3063        Inf  moderate_ED             false       false  TLJ1             missing            54               true  04c8bc98-d6ba-11f0-2861-0dc9995b…  RDU                Frimley                            34.6315   34.7582           true      36.5
                                    3 │ 04c8bc98-d6ba-11f0-2861-0dc9995b…  38.5534   41.9337  Inf         42.3331      Inf        Inf    Inf         Inf             Inf  severe_hosp_long_stay    true       false  TLJ1                  67            96              false  04c8bc98-d6ba-11f0-2861-0dc9995b…  RDU                Frimley                            42.3331   42.4869           true      43.5
                                    4 │ 6c780ade-d6b9-11f0-21b0-717e8d16…  40.6852   44.2124   44.7715   Inf           Inf        Inf     44.8347     66.4395        Inf  moderate_ED             false        true  TLI5                  30            30              false  6c780ade-d6b9-11f0-21b0-717e8d16…  RAL                Royal Free                         44.7715   44.98             true      46.5
                                    5 │ 04c8bc98-d6ba-11f0-2861-0dc9995b…  46.8122  Inf       Inf         56.0806      Inf        Inf     56.5193     79.8787        Inf  severe_hosp_short_stay  false       false  TLJ1                  43            78              false  04c8bc98-d6ba-11f0-2861-0dc9995b…  RTH                Oxford                             56.0806   56.2707           true      57.5
                                    ⋮  │                 ⋮                     ⋮        ⋮         ⋮          ⋮         ⋮         ⋮          ⋮           ⋮           ⋮                ⋮               ⋮        ⋮           ⋮            ⋮             ⋮                ⋮                          ⋮                          ⋮                     ⋮                     ⋮             ⋮            ⋮           ⋮
                                    550 │ 6c780ade-d6b9-11f0-21b0-717e8d16…  87.8574   95.8388  104.501    Inf           Inf        Inf    104.752     120.368         Inf  moderate_ED             false       false  TLJ3                  14            14              false  6c780ade-d6b9-11f0-21b0-717e8d16…  RHM                Southampton                       104.501   104.552            true     106.5
                                    551 │ 6c780ade-d6b9-11f0-21b0-717e8d16…  86.7316  102.831   105.283    Inf           Inf        Inf    105.286     127.861         Inf  moderate_ED             false       false  TLI7                   4             4              false  6c780ade-d6b9-11f0-21b0-717e8d16…  RYJ                Imperial                          105.283   105.431            true     106.5
                                    552 │ 6c780ade-d6b9-11f0-21b0-717e8d16…  87.3668   96.5047  108.095    Inf           Inf        Inf    108.234     111.909         Inf  moderate_ED             false        true  TLI7                  53            38              false  6c780ade-d6b9-11f0-21b0-717e8d16…  RYJ                Imperial                          108.095   108.291            true     109.5
                                    553 │ 6c780ade-d6b9-11f0-21b0-717e8d16…  88.2568  Inf       Inf        105.896       Inf        Inf    107.983     110.4           Inf  severe_hosp_long_stay   false        true  TLI7                  34            54              false  6c780ade-d6b9-11f0-21b0-717e8d16…  RYJ                Imperial                          105.896   107.89             true     109.5
"""
function sample_hosp_cases_n(; p = NBPMscape.P
                                , hosp_cases::DataFrame # filtered from simulation of infections (G df filtered for value in thosp column)
                                , initial_dow::Int64 = 1 # Day of the week that initial case imported. Day of week codes: Sunday =1, Monday = 2, ..., Saturday = 7 
                                # Sampling parameters
                                , n_hosp_samples_per_week::Int # Total number of hospital samples to be taken per week
                                , sample_allocation::String = "equal" # "equal" or "weighted"
                                , sample_proportion_adult::Any = "free" # "free" or numeric decimal, e.g. 0.75. Indicates split of sample target 
                                                                        # between adults and children. "free" indicates that no split is specified
                                , hariss_nhs_trust_sampling_sites::DataFrame # List of NHS Trusts in HARISS sampling network
                                , swab_time_mode::Real = 0.25 # Assume swabbing peaks at 6hrs (=0.25 days) after attendance/admission at hospital
                                , swab_proportion_at_48h::Real = 0.9 # Assume 90% of swabs are taken within 48hrs (=2 days) of attendance/admission at hospital
                                , proportion_hosp_swabbed::Real = 0.9 # Assume 90% of ARI attendances are swabbed
                                , only_sample_before_death::Bool = true # There is a possibilty of swabbing time being drawn after death so 'true' here will constrain tswab to tdeceased
                                # Hospital parameters
                                , ed_discharge_limit::Float64 = NBPMscape.P.tdischarge_ed_upper_limit # days. Assume that people attending the Emergency Department are discharged within this time limit.
                                , hosp_short_stay_limit::Float64 = NBPMscape.P.tdischarge_hosp_short_stay_upper_limit # days. Assume that people attending the Emergency Department are discharged within this time limit.
                                , nhs_trust_catchment_pop = NHS_TRUST_CATCHMENT_POP_ADULT_CHILD
                                #, nhs_trust_ae_12m = AE_12M
                                # Seasonal values
                                , hosp_ari_admissions::Int # Estimate of weekly hospital ARI admissions (excluding pathogen X being simulated)
                                , hosp_ari_admissions_adult_p::Float64 = 0.52# Proportion of ED ARI admissions that are adults (16y and over)
                                , hosp_ari_admissions_child_p::Float64 = 0.48 # Proportion of ED ARI admissions that are children (<16y)
                                , ed_ari_destinations_adult::DataFrame = DataFrame( destination = [:discharged,:short_stay,:longer_stay]
                                                                                  , proportion_of_attendances = [0.628,0.030,0.342]
                                                                                    )
                                , ed_ari_destinations_child::DataFrame = DataFrame( destination = [:discharged,:short_stay,:longer_stay]
                                                                                    , proportion_of_attendances = [0.861,0.014,0.125]
                                                                                    )
                                , weight_samples_by = "ae_mean" # or "catchment_pop"
                                , phl_collection_dow::Vector{Int64} = [2,5] # Day(s) of week that swab samples will be collected from public health labs. Day of week codes: Sunday = 1,... Saturday = 7.
                                )
    
    # Generate Gamma distribution for times between attendance/admission to hospital and swabbing, given the mode and proportion taken within 48hrs 
    swab_time_gamma_d = gamma_params_from_mode_cdf( mode_val = swab_time_mode, cdf_at_2 = swab_proportion_at_48h
                                                        , lower_shape = 1.0 + 1e-6, upper_shape = 10)

    # Ensure only hospital cases included
    hosp_cases = hosp_cases[ isfinite.(hosp_cases.ted) .| isfinite.(hosp_cases.thospital) , : ] 

    if size( hosp_cases, 1 ) == 0 
        return( DataFrame() )
    end

    # Find maximum time for ED attendance or hospital admission
    # First extract ted and thospital times
    ted_vals = hosp_cases[!, :ted]
    ted_finite_vals = ted_vals[isfinite.(ted_vals)] == Float64[] ? 0 : ted_vals[isfinite.(ted_vals)]
    thospital_vals = hosp_cases[!, :thospital]
    thospital_finite_vals = thospital_vals[isfinite.(thospital_vals)] == Float64[] ? 0 : thospital_vals[isfinite.(thospital_vals)]

    # Then find the maximum value
    #max_ted_thosp = maximum( [maximum( filter(row -> isfinite(row[:ted]), hosp_cases).ted ), maximum( filter(row -> isfinite(row[:thospital]), hosp_cases).thospital )])
    max_ted_thosp = maximum( [maximum( ted_finite_vals ), maximum( thospital_finite_vals )]) #max_ted_thosp =0
    max_week = Int( ceil( max_ted_thosp / 7 ) + 1 ) # Add an extra week to ensure that simulated enough background ARI swab times to cover the pathogen X swab times and the courier collection times etc
    
    # If there are no hospital cases then return an empty df, otherwise sample hospital cases
    if size(hosp_cases,1) == 0
        #error("hosp_cases dataframe contains no infections (rows)")
        return(hosp_cases_sub = empty( hosp_cases ))
    end

    ### NHS Trust catchment data is curently only for England
    wales_regions = ["TLL3", "TLL4", "TLL5"]
    hosp_cases_Eng = filter(row -> !in( row[:homeregion], wales_regions), hosp_cases)
    
    # Assume that only a proportion of hospital cases are swabbed, e.g. 90% swabbed and 10% not (on average) for practical reasons
    n_pathx_cases_swabbed =  rand( Binomial( size(hosp_cases_Eng,1), proportion_hosp_swabbed ) )
    # Sample from the ICU cases
    hosp_cases_Eng = hosp_cases_Eng[sample( 1:size(hosp_cases_Eng,1), n_pathx_cases_swabbed, replace=false ), :]
        
    # Remove sites with no NHS Trust Code 
    # Not relevant for current HARISS PHLs but if expanded would need to look at how to incorporate sites without NHS Trust codes, i.e. non-NHS England sites
    nhs_trust_site_sample_targets = copy(hariss_nhs_trust_sampling_sites)
    nhs_trust_site_sample_targets = nhs_trust_site_sample_targets[nhs_trust_site_sample_targets[:, 1] .!= "NA", :] # There are no non-NHS England sites in HARISS network currently but if there were in the future then this would change the allocation of samples (i.e. if in reality samples were allocated to Wales, Scotland or NI, here in the model they would be reallocated to England)
    
    ## Allocate sample target across sites
    
    if weight_samples_by == "catchment_pop"

        # Add catchment area population data 
        rename!(nhs_trust_catchment_pop, :TrustCode => :NHS_Trust_code )
        nhs_trust_site_sample_targets = innerjoin(nhs_trust_site_sample_targets, nhs_trust_catchment_pop, on = :NHS_Trust_code)
        # CHECK on NHS Trusts remove from nhs_trust_site_sample_targets by innerjoin - filter(row -> !(row.NHS_Trust_code in skipmissing(nhs_trust_catchment_pop.NHS_Trust_code)), hariss_nhs_trust_sampling_sites) 
        # Note that any specialist NHS Trusts that have been removed from the nhs_trust_catchment_pop already will also be removed from nhs_trust_site_sample_targets by the innerjoin
        
        # Create df for Public Health Laboratory catchment population and sample targets
        # Group NHS Trust catchment area populations by PHL
        nhs_trust_site_sample_targets_groupby_phl = groupby(nhs_trust_site_sample_targets, :public_health_laboratory)
        phl_sample_targets = combine(nhs_trust_site_sample_targets_groupby_phl, :public_health_laboratory, :catchment_prop_of_total_sum => sum => :catchment_prop_of_total_pop)
        phl_sample_targets = unique(phl_sample_targets)
        # Compute proportion of PHL catchment population at each PHL
        phl_sample_targets[!,:catchment_pop_prop_of_phls] .= phl_sample_targets.catchment_prop_of_total_pop / sum( phl_sample_targets.catchment_prop_of_total_pop )

        # Allocate samples across the public health laboratories (PHL) (not the sites/NHS Trusts)
        if sample_allocation == "equal"
            # Split the number of samples across PHL's equally
            phl_sample_targets[!,:sample_target_per_week] .= Int( round( n_hosp_samples_per_week / nrow(phl_sample_targets); digits = 0 ) )
        
        elseif sample_allocation == "weighted"
            # Split the number of samples across PHLs according to the total NHS Trust catchment population feeding into each PHL
            phl_sample_targets[!,:sample_target_per_week] = allocate_with_rounding(total = n_hosp_samples_per_week, weights = phl_sample_targets.catchment_pop_prop_of_phls)
            # Check sum(phl_sample_targets[!,:sample_target_per_week]) == n_hosp_samples_per_week
        end
    
    elseif weight_samples_by == "ae_mean"
        
        # Add A&E attendance data (aka Emergency Department (ED))
        #rename!(AE_12M, :("Org Code") => :NHS_Trust_code )
        nhs_trust_site_sample_targets = innerjoin(nhs_trust_site_sample_targets, AE_12M[:,[:NHS_Trust_code,:mean_12m_prop]], on = :NHS_Trust_code) # CHECK sum( AE_12M.mean_12m_prop )
                
        # Create df for Public Health Laboratory proportion of 12m A&E mean attendances and sample targets
        # Group NHS Trust A&E attendance data by PHL
        nhs_trust_site_sample_targets_groupby_phl = groupby(nhs_trust_site_sample_targets, :public_health_laboratory)
        phl_sample_targets = combine(nhs_trust_site_sample_targets_groupby_phl, :public_health_laboratory, :mean_12m_prop => sum => :phl_sum_mean_12m_prop)
        phl_sample_targets = unique(phl_sample_targets)
        # Compute proportion of PHL A&E attendances at each PHL
        phl_sample_targets[!,:ae_12m_mean_prop_of_phls] .= phl_sample_targets.phl_sum_mean_12m_prop / sum( phl_sample_targets.phl_sum_mean_12m_prop )

        # Allocate samples across the public health laboratories (PHL) (not the sites/NHS Trusts)
        if sample_allocation == "equal"
            # Split the number of samples across PHL's equally
            phl_sample_targets[!,:sample_target_per_week] .= Int( round( n_hosp_samples_per_week / nrow(phl_sample_targets); digits = 0 ) )
        
        elseif sample_allocation == "weighted"
            # Split the number of samples across PHLs according to the total NHS Trust 12m mean A&E attendance feeding into each PHL
            phl_sample_targets[!,:sample_target_per_week] = allocate_with_rounding(total = n_hosp_samples_per_week, weights = phl_sample_targets.ae_12m_mean_prop_of_phls)
            # Check sum(phl_sample_targets[!,:sample_target_per_week]) == n_hosp_samples_per_week
        end

    end

    ## Allocate site sample target across adults and children
    # If input is "free" then no allocation between adults and children is necessary
    # These two lines also serve to create the columns which can be filled in the
    # following for loop if required
    #insertcols!(phl_sample_targets, :sample_target_per_week_adult => Inf)
    insertcols!(phl_sample_targets, :sample_target_per_week_adult => -1)
    insertcols!(phl_sample_targets, :sample_target_per_week_child => -1)
    # Otherwise use the adult proportion to allocate the site target between adults and children
    if sample_proportion_adult != "free" # Test # sample_proportion_adult = 0.75
        for i in 1:nrow(phl_sample_targets) # i=1
            phl_sample_targets[ i ,[:sample_target_per_week_adult, :sample_target_per_week_child]] = allocate_with_rounding(  total = phl_sample_targets[!,:sample_target_per_week][i]
                                                                                                                            , weights = [sample_proportion_adult,  (1-sample_proportion_adult)])
        end
    end
    
    ### Allocate NHS Trust and PHL to each simulated pathogen X hospital case
    # Allocation is determined using probability of admission to a particular trust given home region
    # Probabilities are derived from MSOA - NHS Trust catchment area analysis
    hosp_nhs_trust_cd = Vector{Any}(undef, size(hosp_cases_Eng,1))
    #hosp_phl = Vector{Any}(undef, size(hosp_cases_Eng,1))
    for i in 1:size(hosp_cases_Eng,1) #i=1
        hr = hosp_cases_Eng[i,:].homeregion
        inf_age = hosp_cases_Eng[i,:].infectee_age
        # Determine which NHS Trust individual is admitted to. Differentiated for adult and paediatric probabilities.
        if inf_age < 16
            hosp_nhs_trust_cd[i] = String( only( wsample( ITL2_TO_NHS_TRUST_PROB_CHILD[:, :NHS_Trust_code], ITL2_TO_NHS_TRUST_PROB_CHILD[:, Symbol(hr)], 1) ) )
        elseif inf_age >= 16
            hosp_nhs_trust_cd[i] = String( only( wsample( ITL2_TO_NHS_TRUST_PROB_ADULT[:, :NHS_Trust_code], ITL2_TO_NHS_TRUST_PROB_ADULT[:, Symbol(hr)], 1) ) )
        end
    end
    insertcols!(hosp_cases_Eng, :hosp_nhs_trust_cd => hosp_nhs_trust_cd)
    
    # Temporary key for NHS Trust code and Public Health Laboratory
    temp_nt_phl = nhs_trust_site_sample_targets[:,[:NHS_Trust_code,:public_health_laboratory]]
    rename!(temp_nt_phl, :NHS_Trust_code => :("hosp_nhs_trust_cd"))
    
    # Add PHL to df for each infected individual 
    # This will also filter the hospital cases for those attending/admitted to NHS Trusts that feed into PHLs being sampled (zero probability of the other hospital cases being sampled)
    hosp_cases_at_sampling_sites = innerjoin(hosp_cases_Eng, temp_nt_phl, on = :hosp_nhs_trust_cd) # Note that this reduces the cases to those that attend an NHS Trust that feeds into a sampling PHL
    # CHECK # StatsBase.countmap(hosp_nhs_trust_cd)
    # CHECK # StatsBase.countmap(hosp_cases_Eng.hosp_nhs_trust_cd)
    # CHECK # StatsBase.countmap(hosp_cases_at_sampling_sites.public_health_laboratory) 

    # Proportion of hospital cases that are admitted to ICU
    # In England
    # nrow( filter(row -> isfinite( row.ticu ) , hosp_cases_Eng)) / nrow(hosp_cases_Eng)
    # At sampling sites
    # nrow( filter(row -> isfinite( row.ticu ) , hosp_cases_at_sampling_sites)) / nrow(hosp_cases_at_sampling_sites)
    
    
    if ( size(hosp_cases_at_sampling_sites,1) == 0 )
        # return an empty dataframe if the size of icu_cases_at_sampling_sites
        #hosp_cases_sub = empty( hosp_cases_at_sampling_sites )
        #pathx_samples = empty( hosp_cases_at_sampling_sites )
        return( hosp_cases_sub = empty( hosp_cases_at_sampling_sites ) )
    else
                
        ## Estimate weekly ARI A&E/ED attendances per PHL

        # Create dfs to store info
        #est_weekly_ed_ari_adult = zeros(Int, size(AE_12M,1)) #Vector{Int}(undef, size(nhs_trust_ari_cc_beds,1)) 
        #est_weekly_ed_ari_child = copy(est_weekly_ed_ari_adult)
        
        # Estimate adult and child ARI A&E/ED attendances and hospital admissions per PHL per week (as an estimate of background ARI)
        phl_sample_targets[!,:est_weekly_bkg_hosp_ari_admissions] = Int.( round.( hosp_ari_admissions * phl_sample_targets.phl_sum_mean_12m_prop) )
        phl_sample_targets[!,:est_weekly_bkg_hosp_ari_admissions_adult] = Int.( round.( hosp_ari_admissions_adult_p * phl_sample_targets.est_weekly_bkg_hosp_ari_admissions) )
        phl_sample_targets[!,:est_weekly_bkg_hosp_ari_admissions_child] = Int.( round.( hosp_ari_admissions_child_p * phl_sample_targets.est_weekly_bkg_hosp_ari_admissions) )
        # CHECK # sum(phl_sample_targets.est_weekly_bkg_hosp_ari_admissions_adult) + sum(phl_sample_targets.est_weekly_bkg_hosp_ari_admissions_child) == sum(phl_sample_targets.est_weekly_bkg_hosp_ari_admissions)
        
        ## Determine timetable for courier collection of samples from PHLs and arrival at central testing lab
        # Courier collections twice per week - assume Monday and Thursday at midday. 
        # Assume swabs taken at hospitals before midday the day before are available for collection at the PHL.

        phl_courier_collection_times = NBPMscape.courier_collection_times( initial_dow = initial_dow, n_weeks = max_week, collection_dow = Int64.(phl_collection_dow ))
        # TEST
        #courier_collection_times( initial_dow = 1, n_weeks = max_week, collection_dow = Int64.(phl_collection_dow ))
        #courier_collection_times( initial_dow = 2, n_weeks = max_week, collection_dow = Int64.(phl_collection_dow ))
        #courier_collection_times( initial_dow = 3, n_weeks = max_week, collection_dow = Int64.(phl_collection_dow ))
        #courier_collection_times( initial_dow = 4, n_weeks = max_week, collection_dow = Int64.(phl_collection_dow ))
        #courier_collection_times( initial_dow = 5, n_weeks = max_week, collection_dow = Int64.(phl_collection_dow ))
        #courier_collection_times( initial_dow = 6, n_weeks = max_week, collection_dow = Int64.(phl_collection_dow ))
        #courier_collection_times( initial_dow = 7, n_weeks = max_week, collection_dow = Int64.(phl_collection_dow ))

        ## Simulate swab times for background ARI ED attendance and hospital admissions
        
        # Define distributions for time between arrival at hospital and swab time depending on length of stay in hospital
        tswab_dists_bkg = Dict(:discharged => Truncated(swab_time_gamma_d,0.0,ed_discharge_limit),
                                :short_stay => Truncated(swab_time_gamma_d, 0.0, hosp_short_stay_limit),
                                :longer_stay => swab_time_gamma_d
                                )

        tswab_dists_pathx = Dict(:moderate_ED => Truncated(swab_time_gamma_d,0.0,ed_discharge_limit)
                                ,:severe_hosp_short_stay => Truncated(swab_time_gamma_d, 0.0, hosp_short_stay_limit)
                                ,:severe_hosp_long_stay => swab_time_gamma_d
                                ,:verysevere => swab_time_gamma_d
                                )

        # Empty vector of DataFrames
        ari_bkg_times_by_phl_df = DataFrame(  pid = String[]
                                            , public_health_laboratory = String[]
                                            , age_group = String[]
                                            , thosp = Float64[]
                                            , tswab = Float64[]
                                            , ED_destination = Symbol[]
                                           )
        
        for j in 1:nrow(phl_sample_targets) # j=1
            
            # Public Health Lab (PHL) name
            phl_name = String(phl_sample_targets[j,:public_health_laboratory])
            # Weekly numbers of swabs for individual PHLs
            phl_bkg_ari_swabs_n_adult = Int( round(phl_sample_targets[j,:est_weekly_bkg_hosp_ari_admissions_adult] * proportion_hosp_swabbed , digits = 0) )
            phl_bkg_ari_swabs_n_child = Int( round(phl_sample_targets[j,:est_weekly_bkg_hosp_ari_admissions_child] * proportion_hosp_swabbed, digits = 0 ) )
                
            for wn in 1:max_week # wn=15

                # Create df for each PHL and DataFrames and push them
                # Severity
                ED_destination_adult = wsample( ed_ari_destinations_adult[:,1] , ed_ari_destinations_adult[:,2], phl_bkg_ari_swabs_n_adult) # CHECK # StatsBase.countmap(ED_destination_adult)
                ED_destination_child = wsample( ed_ari_destinations_child[:,1] , ed_ari_destinations_child[:,2], phl_bkg_ari_swabs_n_child) # CHECK # StatsBase.countmap(ED_destination_child)
                
                # Simulate the time (in days) for hospital arrival, thosp = start time of week + random draw from 7 days
                thosp_adult = ((wn-1)*7) .+ rand(phl_bkg_ari_swabs_n_adult).*7 
                thosp_child = ((wn-1)*7) .+ rand(phl_bkg_ari_swabs_n_child).*7
                
                # Simulate the time (in days) that swab was taken, tswab = time of arrival in hospital + draw from Gamma distribution and dependent on length of hospital stay
                # samples different time distributions depending on duration of stay in hospital
                tswab_adult = thosp_adult .+ rand.(getindex.(Ref(tswab_dists_bkg), ED_destination_adult))
                tswab_child = thosp_child .+ rand.(getindex.(Ref(tswab_dists_bkg), ED_destination_child))

                # Adult
                df_temp_adult = DataFrame( pid = String("background_ari")
                                            , public_health_laboratory = repeat([phl_name], phl_bkg_ari_swabs_n_adult) 
                                            , age_group = repeat(["adult"], phl_bkg_ari_swabs_n_adult)
                                            , thosp = thosp_adult
                                            , tswab = tswab_adult
                                            , ED_destination = ED_destination_adult
                                            )
                
                # Child
                df_temp_child = DataFrame( pid = String("background_ari")
                                            , public_health_laboratory = repeat([phl_name], phl_bkg_ari_swabs_n_child) 
                                            , age_group = repeat(["child"], phl_bkg_ari_swabs_n_child)
                                            , thosp = thosp_child
                                            , tswab = tswab_child
                                            , ED_destination = ED_destination_child
                                            )
                                
                # Add to main df                    
                append!(ari_bkg_times_by_phl_df, df_temp_adult)
                append!(ari_bkg_times_by_phl_df, df_temp_child)
            end
        end
            
        # Does the swab time meet the criteria of being taken within 48h of arrival at hospital?
        tswab_diff = (ari_bkg_times_by_phl_df.tswab .- ari_bkg_times_by_phl_df.thosp)
        ari_bkg_times_by_phl_df.tswab_sub_48h = Bool.( tswab_diff .< (48.0/24.0) )
        # CHECK # StatsBase.countmap(ari_bkg_times_by_phl_df.tswab_sub_48h)
        
        # Filter ARI background cases for those that meet 48h swab time cut-off
        ari_bkg_times_by_phl_df = filter( row -> row.tswab_sub_48h == true, ari_bkg_times_by_phl_df ) # CHECK # StatsBase.countmap(ari_bkg_times_by_phl_df.tswab_sub_48h)

        ## Simulate swab times for pathogen X

        # Helper function to select ted or thospital, depending which is defined
        function get_finite_value(row, col1, col2)
            v1 = row[col1]
            v2 = row[col2]
            return isfinite(v1) ? v1 : v2  # Assumes v2 is finite if v1 is Inf
        end
        
        # Determine time of swab but drawing the time between attendance/admission and swabbing from a Gamma distribution with mode 0.25 and cdf = 90% at 2 days.
        #hosp_cases_at_sampling_sites.tswab = [get_finite_value(row, :ted, :thospital) for row in eachrow(hosp_cases_at_sampling_sites)] .+ rand(swab_time_gamma_d, nrow(hosp_cases_at_sampling_sites))
        hosp_cases_at_sampling_sites.ted_or_thospital = [get_finite_value(row, :ted, :thospital) for row in eachrow(hosp_cases_at_sampling_sites)]
        hosp_cases_at_sampling_sites.tswab = hosp_cases_at_sampling_sites.ted_or_thospital .+ rand.(getindex.(Ref(tswab_dists_pathx), hosp_cases_at_sampling_sites.severity))

        # Potentially constrain time of swabbing for deceased individuals to time of death, tdeceased
        if only_sample_before_death == true
            hosp_cases_at_sampling_sites.tswab .= ifelse.(hosp_cases_at_sampling_sites.tswab .> hosp_cases_at_sampling_sites.tdeceased
                                                         , hosp_cases_at_sampling_sites.tdeceased
                                                         , hosp_cases_at_sampling_sites.tswab
                                                         )
            # CHECK # sum( hosp_cases_at_sampling_sites.tswab .> hosp_cases_at_sampling_sites.tdeceased )
            # Alternative method removes swabs where tswab > tdeceased
            #filter(row -> row.tswab .> row.tdeceased, hosp_cases_at_sampling_sites)
        end

        # Does the swab time meet the within 48h criteria?
        tswab_diff = (hosp_cases_at_sampling_sites.tswab .- hosp_cases_at_sampling_sites.ted_or_thospital)
        hosp_cases_at_sampling_sites.tswab_sub_48h = Bool.( tswab_diff .< 48.0/24.0 )
        # CHECK # StatsBase.countmap(hosp_cases_at_sampling_sites.tswab_sub_48h)
        
        # CHECK # hosp_cases_at_sampling_sites[:,[:pid,:ted,:thospital,:tswab, :tswab_sub_48h]]
        # CHECK # filter( row -> row.tswab_sub_48h == false, hosp_cases_at_sampling_sites )
     
        # Assume only swabs taken before midday the day before at a hospital can be collected from the PHL
        swab_time_cut_offs = phl_courier_collection_times .- 1

        # Filter for swabs taken within 48h of attendance or admission
        hosp_cases_at_sampling_sites_sub48h = filter( row -> row.tswab_sub_48h == true, hosp_cases_at_sampling_sites)

        ## Identify top n samples (depending on PHL, total number and PHL allocation method) for adults and/or children
        
        # Build df containing samples collected from each site
        pathx_samples = copy(hosp_cases_at_sampling_sites_sub48h)
        pathx_samples = empty(pathx_samples)
        # Add column for sample courier times
        pathx_samples.tcourier = Vector{Union{Missing, Float64}}(missing, nrow(pathx_samples))

        if size(hosp_cases_at_sampling_sites_sub48h,1) == 0
          return( DataFrame() )
        
        else
          # Loop through cut-off times for courier collection
            for c in 1:length(swab_time_cut_offs) # c=2
                # Filter swabs for those meeting the required criteria and meeting the cut-off time (between the current cut-off time and previous cut-off time)
                path_x_cases_temp = filter( row -> (row.tswab <= swab_time_cut_offs[c] && row.tswab > ( c == 1 ? 0 : swab_time_cut_offs[c-1])),  hosp_cases_at_sampling_sites_sub48h ) #path_x_cases_temp = filter( row -> (row.tswab <= swab_time_cut_offs[c] && row.tswab > swab_time_cut_offs[c-1]),  hosp_cases_at_sampling_sites_sub48h )
                bkg_cases_temp    = filter( row -> (row.tswab <= swab_time_cut_offs[c] && row.tswab > ( c == 1 ? 0 : swab_time_cut_offs[c-1])),  ari_bkg_times_by_phl_df ) # bkg_cases_temp = filter( row -> (row.tswab <= swab_time_cut_offs[c] && row.tswab > max(try swab_time_cut_offs[c-1]),  ari_bkg_times_by_phl_df )

                # Loop through PHLs. Combined path_x and background cases and sort by tswab. Take the top n samples based on PHL sample targets.
                for i in 1:nrow(phl_sample_targets) # i=1
                    phl_name = phl_sample_targets[i,:public_health_laboratory]
                  
                    # Filter background and pathogen X cases by PHL
                    path_x_cases_temp_phl = filter( row -> row.public_health_laboratory == phl_name, path_x_cases_temp)
                    bkg_cases_temp_phl = filter( row -> row.public_health_laboratory == phl_name, bkg_cases_temp)

                    # Select samples based on age group target
                    if sample_proportion_adult == "free"
                  
                        # Combine pathx and bkg cases dfs
                        temp_df = append!( path_x_cases_temp_phl[:,[:pid,:tswab]] , bkg_cases_temp_phl[:,[:pid,:tswab]] )
                    
                        # Number of samples for this PHL
                        n_samples_per_week = only(filter( row -> row.public_health_laboratory == phl_name, phl_sample_targets)[:,:sample_target_per_week]) 
                        n_collections_per_week = length(phl_collection_dow)
                        n_samples_per_collection_vec = allocate_with_rounding(;total = n_samples_per_week, weights = fill(1/n_collections_per_week,n_collections_per_week))
                        collection_dow = Int64(phl_courier_collection_times[c] % 7 + initial_dow - 0.5) # Finds the day of the week of this collection (adjusted for 0.5 midday collection time/cut-off) 
                        collection_number_this_week = findfirst(==(collection_dow), phl_collection_dow) # Finds which number collection day this is, i.e. collection 2 out of 3 per week
                        n_samples_for_this_collection = n_samples_per_collection_vec[ collection_number_this_week ]
                        # Samples selected to send to PHL (i.e. the most recent n_samples or if less available then all of them)
                        samples_sent_to_phl = sort!( temp_df, :tswab , rev = true)[ 1 : min( n_samples_for_this_collection, nrow(temp_df) ) , : ]
                        pathx_sample_pids =  filter( row -> row.pid != "background_ari", samples_sent_to_phl )[:,:pid]
                        
                        # Add swab courier collection time
                        pathx_samples_temp = filter( row -> row.pid in pathx_sample_pids, hosp_cases_at_sampling_sites_sub48h)
                        pathx_samples_temp.tcourier .= phl_courier_collection_times[ c ]

                        # Add courier collected samples to df
                        append!( pathx_samples, pathx_samples_temp )
                  
                    else     
                        # When there is a split between adult and child specified:

                        # Split cases by adult and child cases
                        path_x_cases_temp_phl_adult = filter( row -> row.infectee_age > 16, path_x_cases_temp_phl)
                        path_x_cases_temp_phl_child = filter( row -> row.infectee_age <= 16, path_x_cases_temp_phl)
                        bkg_cases_temp_phl_adult = filter( row -> row.age_group == "adult", bkg_cases_temp_phl)
                        bkg_cases_temp_phl_child = filter( row -> row.age_group == "child", bkg_cases_temp_phl)

                        # Combine pathx and bkg cases dfs
                        temp_df_adult = append!( path_x_cases_temp_phl_adult[:,[:pid,:tswab]] , bkg_cases_temp_phl_adult[:,[:pid,:tswab]] )
                        temp_df_child = append!( path_x_cases_temp_phl_child[:,[:pid,:tswab]] , bkg_cases_temp_phl_child[:,[:pid,:tswab]] )
                                        
                        # Number of samples for this PHL, split by adult and child cases
                        n_samples_adult_per_week = Int(only(filter( row -> row.public_health_laboratory == phl_name, phl_sample_targets)[:,:sample_target_per_week_adult]))
                        n_samples_child_per_week = Int(only(filter( row -> row.public_health_laboratory == phl_name, phl_sample_targets)[:,:sample_target_per_week_child]))
                        
                        # Convert number of samples per week to number of samples for individual collection
                        n_collections_per_week = length(phl_collection_dow)
                        n_samples_adult_per_collection_vec = allocate_with_rounding(;total = n_samples_adult_per_week, weights = fill(1/n_collections_per_week,n_collections_per_week))
                        n_samples_child_per_collection_vec = allocate_with_rounding(;total = n_samples_child_per_week, weights = fill(1/n_collections_per_week,n_collections_per_week))
                        collection_dow = Int64(phl_courier_collection_times[c] % 7 + initial_dow - 0.5) # Finds the day of the week of this collection (adjusted for 0.5 midday collection time/cut-off) 
                        collection_number_this_week = findfirst(==(collection_dow), phl_collection_dow) # Finds which number collection day this is, i.e. collection 2 out of 3 per week
                        n_samples_adult_for_this_collection = n_samples_adult_per_collection_vec[ collection_number_this_week ]
                        n_samples_child_for_this_collection = n_samples_child_per_collection_vec[ collection_number_this_week ]
                        

                        # If not enough samples in either adult or child age groups then make up the sample target from the other age group
                        n_swabs_adult = nrow(temp_df_adult) # n_swabs_adult = 10
                        n_swabs_child = nrow(temp_df_child) # n_swabs_child = 10

                        # Shortfall in swabs to meet sample targets
                        n_shortfall_adult = -(min(0, n_swabs_adult - n_samples_adult_for_this_collection)) # n_swabs_adult = 10 # n_samples_adult = 12
                        n_shortfall_child = -(min(0, n_swabs_child - n_samples_child_for_this_collection))

                        # Adjust sample targets to account for shortfalls
                        n_samples_adult_adj = min( n_swabs_adult, n_samples_adult_for_this_collection + n_shortfall_child )
                        n_samples_child_adj = min( n_swabs_child, n_samples_child_for_this_collection + n_shortfall_adult )
                        
                        # Samples selected to send to PHL (i.e. the most recent n_samples or if less available then all of them)
                        samples_sent_to_phl_adult = sort!( temp_df_adult, :tswab , rev = true)[ 1 : n_samples_adult_adj , : ]
                        samples_sent_to_phl_child = sort!( temp_df_child, :tswab , rev = true)[ 1 : n_samples_child_adj , : ]
                        
                        pathx_sample_pids_adult =  filter( row -> row.pid != "background_ari", samples_sent_to_phl_adult )[:,:pid]
                        pathx_sample_pids_child =  filter( row -> row.pid != "background_ari", samples_sent_to_phl_child )[:,:pid]

                        # Add swab courier collection time
                        pathx_samples_temp_adult = filter( row -> row.pid in pathx_sample_pids_adult, hosp_cases_at_sampling_sites_sub48h)
                        pathx_samples_temp_adult.tcourier .= phl_courier_collection_times[ c ]
                        pathx_samples_temp_child = filter( row -> row.pid in pathx_sample_pids_child, hosp_cases_at_sampling_sites_sub48h)
                        pathx_samples_temp_child.tcourier .= phl_courier_collection_times[ c ]
                        
                        append!( pathx_samples, pathx_samples_temp_adult )
                        append!( pathx_samples, pathx_samples_temp_child )

                    end

                end
          
            end
        
        end
    
    end
    return( pathx_samples )
end