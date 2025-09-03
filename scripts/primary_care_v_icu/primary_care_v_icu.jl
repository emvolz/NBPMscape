#### Script to investigate the addition of primary care sampling
#= 
Description:
- Primary care sampling is modelled on the current Oxford-RCGP RSC surveillance
- Script compares ICU sampling vs Primary care sampling vs combination of both
Outline:
1) Load simulation
2) Extract sample of cases
    a) 25% of ICU cases
    b) X% of GP cases based on RCGP-RSC surveillance parameters and calculating probabilities of 
        i) going to a GP that takes surveillance swabs
        ii) individual is swabbed
        iii) metagenomics is performed on the swab
3) Compute report/detection times for sampled cases
4) Compute and plot results for median detection times across n simulation replicates
    a) ICU only
    b) primary care only
    c) combination of both
=#

# Load packages
using Revise
using NBPMscape
using GLM, Statistics
using Distributions
using DataFrames
using StatsBase
using CSV 

#using Pkg
#Pkg.add(PackageSpec(name="JLD2", version="0.5.15"))#version="1.11.2"))#version="0.4.54"))
using JLD2
# ] status JLD2
#version = Pkg.TOML.parsefile(joinpath(pkgdir(JLD2), "Project.toml"))["version"]
#println(version)

# Load simulation
sims = load("covidlike-1.1.1-sims.jld2", "sims")
#@load "covidlike-1.1.1-sims.jld2" sims

# TODO Or adapt to use (because smaller files have less issues with reloading)
#@load "covidlike-1.1-sims_filtered_G_icu.jld2" sims_filtered_G_icu
#@load "covidlike-1.1-sims_filtered_G_gp.jld2"  sims_filtered_G_gp

# How many replicates are in this simulation
n_replicates = length(sims)

# ICU sampling proportion
p_icu = 0.25 # Assumption TODO needs refining 

# Time to process sample and report results / declare detection
#icu_turnaround_time = 3
icu_turnaround_time = [2,4] # 'icu_v_pc()' function assumes this is equal to 'pc_swab_turnaround_time'

# Parameters for existing Oxford-RCGP RSC primary care surveillance
gp_practices_total = 6199 # Total number of GP practices in England at July 2025. Source: BMA analysis (https://www.bma.org.uk/advice-and-support/nhs-delivery-and-workforce/pressures/pressures-in-general-practice-data-analysis) of NHS Digital General Practice Workforce Statistics (https://digital.nhs.uk/data-and-information/publications/statistical/general-and-personal-medical-services) [Accessed 2 Sep 2025]  
gp_practices_swab = 300 # Number of GP practices taking swabs for virology surveillance. Source: Data quality report: national flu and COVID-19 surveillance report (27 May 2025)
gp_swabs_mg = [100,200] #[319, 747] # Assumed number of swabs that are metagenomic sequenced for investigating impact
#ITL225CD_wales = ITL2SIZE[37:39,:ITL225CD]
#ITL2SIZE_eng = ITL2SIZE[.!in(ITL2SIZE.ITL225CD, ITL225CD_wales), :]
#ITL2SIZE_eng = ITL2SIZE[ ITL2SIZE.ITL225CD .!= , :]
#pop_eng = sum(ITL2SIZE[1:36,3])
pop_eng = 5.7106398e7 # Population of England. Source: ONS mid-year 2022 UK population data disaggregated by various geo levels and age groups and gender. [Accessed 6 November 2024] Available at https://www.ons.gov.uk/peoplepopulationandcommunity/populationandmigration/populationestimates/datasets/populationestimatesforukenglandandwalesscotlandandnorthernireland 
gp_ari_consults = [180, 327] # Number of ARI consultations per 100k of England population per week [mean summer 2024, mean winter 2024/25]. Source: Analysis of data extracted from RCGP Research & Surveillance Centre (RSC) Virology Dashboard [Accessed 29 Aug 2025]
gp_ari_swabs = [319, 747] # Number of swabs taken from suspected ARI per week [mean summer 2024, mean winter 2024/25]. Source: Analysis of data extracted from RCGP Research & Surveillance Centre (RSC) Virology Dashboard [Accessed 29 Aug 2025]
pc_swab_turnaround_time = [2,4] # [min,max] number of days between swab sample being taken and results received. Source: Data quality report: national flu and COVID-19 surveillance report (27 May 2025).  Assume same for metagenomic testing.

function icu_v_pc_td(; p_icu = 0.25 # Assumption TODO needs refining 
                     #, icu_turnaround_time = [2,4] # Time to process sample and report results / declare detection
                       # Parameters for existing Oxford-RCGP RSC primary care surveillance
                     , gp_practices_total = 6199 # Total number of GP practices in England at July 2025. Source: BMA analysis (https://www.bma.org.uk/advice-and-support/nhs-delivery-and-workforce/pressures/pressures-in-general-practice-data-analysis) of NHS Digital General Practice Workforce Statistics (https://digital.nhs.uk/data-and-information/publications/statistical/general-and-personal-medical-services) [Accessed 2 Sep 2025]  
                     , gp_practices_swab = 300 # Number of GP practices taking swabs for virology surveillance. Source: Data quality report: national flu and COVID-19 surveillance report (27 May 2025)
                     , gp_swabs_mg = [100,200] #[319, 747] # Assumed number of swabs that are metagenomic sequenced for investigating impact
                     , pop_eng = 5.7106398e7 # Population of England. Source: ONS mid-year 2022 UK population data disaggregated by various geo levels and age groups and gender. [Accessed 6 November 2024] Available at https://www.ons.gov.uk/peoplepopulationandcommunity/populationandmigration/populationestimates/datasets/populationestimatesforukenglandandwalesscotlandandnorthernireland 
                     , gp_ari_consults = [180, 327] # Number of ARI consultations per 100k of England population per week [mean summer 2024, mean winter 2024/25]. Source: Analysis of data extracted from RCGP Research & Surveillance Centre (RSC) Virology Dashboard [Accessed 29 Aug 2025]
                     , gp_ari_swabs = [319, 747] # Number of swabs taken from suspected ARI per week [mean summer 2024, mean winter 2024/25]. Source: Analysis of data extracted from RCGP Research & Surveillance Centre (RSC) Virology Dashboard [Accessed 29 Aug 2025]
                     , pc_swab_turnaround_time = [2,4] # [min,max] number of days between swab sample being taken and results received. Source: Data quality report: national flu and COVID-19 surveillance report (27 May 2025).  Assume same for metagenomic testing.
                     #, p_pc # TODO Needs to be uncommented when want to set this directly
                    )
    # Probability of primary care surveillance events
    prob_visiting_gp_swab = minimum( [1, gp_practices_swab / gp_practices_total ])
    prob_swabbed_summer = minimum([1, (100000*gp_ari_swabs[1]) / (pop_eng*gp_ari_consults[1])])
    prob_swabbed_winter = minimum([1, (100000*gp_ari_swabs[2]) / (pop_eng*gp_ari_consults[2])])
    prob_mg_min_summer = minimum([1, gp_swabs_mg[1] / gp_ari_swabs[1] ])
    prob_mg_min_winter = minimum([1, gp_swabs_mg[1] / gp_ari_swabs[2] ])
    prob_mg_max_summer = minimum([1, gp_swabs_mg[2] / gp_ari_swabs[1] ])
    prob_mg_max_winter = minimum([1, gp_swabs_mg[2] / gp_ari_swabs[2] ])

    # Primary care sampling proportion - equivalent to probability a swab being taken and a metagenomic sequence being taken if attend a GP
    p_pc_mg_min_summer = prob_visiting_gp_swab * prob_swabbed_summer * prob_mg_min_summer
    p_pc_mg_min_winter = prob_visiting_gp_swab * prob_swabbed_winter * prob_mg_min_winter
    p_pc_mg_max_summer = prob_visiting_gp_swab * prob_swabbed_summer * prob_mg_max_summer
    p_pc_mg_max_winter = prob_visiting_gp_swab * prob_swabbed_winter * prob_mg_max_winter
    # Check if primary care sampling proportion, p_pc, was entered as an input to the function
    if @isdefined p_pc 
        #continue
    else # Define based on calculated probabilities
        p_pc = [ p_pc_mg_min_summer, p_pc_mg_max_summer, p_pc_mg_min_winter, p_pc_mg_max_winter ]
    end    
    println(p_pc)
    #p_pc = 0.25 # Test to compare TDs with ICU with 0.25 sampling proportion
    #swabs_at_25_summer
    #swabs_at_25_winter

    # Create df to store results
    sim_tds = DataFrame( [fill(Inf,length(sims)) for _ in 1:19], ["sim_n"
                                                                ,"ICU_TD","ICU_3TD"                              # time to detection results using ICU sampling only
                                                                ,"PC_TD_min_summer","PC_3TD_min_summer"          # time to detection results using primary care sampling only with LOWER  number of metagenomic swab tests in SUMMER
                                                                ,"PC_TD_max_summer","PC_3TD_max_summer"          # time to detection results using primary care sampling only with HIGHER number of metagenomic swab tests in SUMMER
                                                                ,"PC_TD_min_winter","PC_3TD_min_winter"          # time to detection results using primary care sampling only with LOWER  number of metagenomic swab tests in WINTER
                                                                ,"PC_TD_max_winter","PC_3TD_max_winter"          # time to detection results using primary care sampling only with HIGHER number of metagenomic swab tests in WINTER
                                                                ,"ICU_PC_TD_min_summer","ICU_PC_3TD_min_summer"  # time to detection results using ICU AND primary care sampling  with LOWER  number of metagenomic swab tests in SUMMER
                                                                ,"ICU_PC_TD_max_summer","ICU_PC_3TD_max_summer"  # time to detection results using ICU AND primary care sampling  with HIGHER number of metagenomic swab tests in SUMMER
                                                                ,"ICU_PC_TD_min_winter","ICU_PC_3TD_min_winter"  # time to detection results using ICU AND primary care sampling  with LOWER  number of metagenomic swab tests in WINTER
                                                                ,"ICU_PC_TD_max_winter","ICU_PC_3TD_max_winter"  # time to detection results using ICU AND primary care sampling  with HIGHER number of metagenomic swab tests in WINTER
                                                                ])

    # Adaptation of sampleforest() function in 'core.jl' and 'median_TD_by_region.jl'
    for s in 1:length(sims)

        # TEST
        #println(s)

        # Add simulation number to results df
        sim_tds[s,:sim_n] = s
        
        #s=1
        #fo = sims_test[s]#sims[s]
            
        # Filter for ICU cases 
        fo = sims[s]
        icu_cases = fo.G[ isfinite.(fo.G.ticu), : ]
        
        # Filter for GP cases (some of these will also become ICU cases)
        gp_cases = fo.G[ isfinite.(fo.G.tgp), : ]    

        ### Record infection time stats
        ## First for ICU cases
        # Sample sizes
        n_icu =  rand( Binomial( size(icu_cases,1), p_icu ) )
        # Subsample of ICU cases
        icu_cases_sub = icu_cases[sample( 1:size(icu_cases,1), n_icu, replace=false ), :]

        if size(icu_cases_sub,1) == 0
            # Record Inf as time to detection for earliest:
            # - 1 ICU case
            sim_tds[s,:ICU_TD] = Inf
            # - 3 ICU cases
            sim_tds[s,:ICU_3TD] = Inf

        elseif size(icu_cases_sub,1) > 0 
                
            #= Sample time has uniform distribution between time of admission to ICU and time of recovery
                TODO THIS MAY NEED TO BE UPDATED FOR LATER VERSIONS WITH MORE COMPLEX CARE PATHWAYS
                NOTE TIME BETWEEN ticu and trecovered can be large and so with uniform distribution of 
                sampling the tsample can be a long time after tinf. 
                Example I saw was tinf = 48.9, ticu = 52.6, trecovered = 91.1, and treport = 90.7
                If we're modelling 100% of ICU cases being sampled then more likely to sampled closer to ticu
            =#
            
            # Generate sample times
            #tsample = (size(g_region,1)>0) ? map( g -> rand( Uniform( g.ticu[1], g.trecovered[1])) , eachrow(g_region) ) : []
            icu_tsample = map( g -> rand( Uniform( g.ticu[1], g.trecovered[1])) , eachrow(icu_cases_sub) )
            # TODO Possible alternative
            #icu_tsample = map( g -> rand( Uniform( g.ticu[1], g.ticu[1]+3)) , eachrow(icu_cases_sub) )
            icu_cases_sub.tsample = icu_tsample 

            # Simulate reports times
            #treport = (size(g_region,1)>0) ? (tsample .+ turnaroundtime) : []
            icu_treport = (icu_tsample .+ rand(Uniform(pc_swab_turnaround_time[1], pc_swab_turnaround_time[2])) )
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

        # Record time stats for GP cases using different sampling parameters 
        # Sample sizes
        n_gp = Vector{Int64}(undef,length(p_pc))
        for i in 1:length(p_pc)
            
            # Sample sizes based on probability of infected case having a metagenomic sample taken at the GP
            n_gp = rand( Binomial( size(gp_cases,1), p_pc[i] ) )
            #println(n_gp)
            
            # Subsample of GP cases
            gp_cases_sub = gp_cases[sample( 1:size(gp_cases,1), n_gp, replace=false ), :]
        
            # If there are NO GP cases in the sample
            if size(gp_cases_sub,1) == 0
                
                #p_pc = [ p_pc_mg_min_summer, p_pc_mg_max_summer, p_pc_mg_min_winter, p_pc_mg_max_winter ]
                
                ## Add to results df for ICU times to detection
                # Record Inf times to detection if no GP cases sampled:
                #  1 primary care case for TD
                sim_tds[s, 2*i + 2 ] = Inf
                #  3 primary care cases for 3TD
                sim_tds[s, 2*i + 3 ] = Inf

                # ICU and primary care combined
                # If there are no GP cases then the combination of primary care and ICU is just ICU only (which was already entered into df above)
                # - 1 case
                sim_tds[s, 2*i + 10] = sim_tds[s,:ICU_TD]
                # - 3 case
                sim_tds[s, 2*i + 11] = sim_tds[s,:ICU_3TD]
                
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
                if @isdefined icu_cases_sub_top3_td # Check if there are any ICU cases...
                    icu_gp_cases_sub_top3_td = sort( append!(icu_cases_sub_top3_td, gp_cases_sub_top3_td), :treport)
                else # ... if not then the combination is just primary care (GP) cases only
                    icu_gp_cases_sub_top3_td = gp_cases_sub_top3_td
                end

                ## Add to results df for primary care times to detection
                # Record the times to detection for earliest:
                # Primary care cases only
                #  1 primary care case for TD
                sim_tds[s, 2*i + 2 ] = gp_cases_sub_top3_td[1,:treport]
                #  3 primary care cases for 3TD
                if size(gp_cases_sub_top3_td,1) >= 3
                    sim_tds[s, 2*i + 3 ] = gp_cases_sub_top3_td[3,:treport]
                else
                    sim_tds[s, 2*i + 3 ] = Inf 
                end
                
                # Primary care COMBINED with ICU cases
                # - 1 case, either ICU or primary care
                sim_tds[s, 2*i + 10] = icu_gp_cases_sub_top3_td[1,:treport]
                # - 3 case, either ICU or primary care
                if size(icu_gp_cases_sub_top3_td,1) >= 3
                    sim_tds[s, 2*i + 11] = icu_gp_cases_sub_top3_td[3,:treport]
                else
                    sim_tds[s, 2*i + 11] = Inf
                end

            end # End of if loop checking whether there are any GP cases sampled

        end # End of different primary care sampling parameters loop

    end # End of loop through different simulations
    return sim_tds
end # End of icu_v_pc_td function

#### Save as .csv file for inspection
## Versions Only recording a combined TD or 3TD when there is a value for primary care TD
#CSV.write("scripts/primary_care_v_icu/sim_TDs.csv", sim_tds) # using 100 and 200 metagenomic samples per week
#CSV.write("scripts/primary_care_v_icu/sim_TDs_all_gp_samples.csv", sim_tds) # using all samples for metagenomic sequencing each week
#CSV.write("scripts/primary_care_v_icu/sim_TDs_primary_care_p25.csv", sim_tds) # metagenomic sampling of 25% of all ARI consultations

## Versions recording a combined TD or 3TD even if there is no value for primary care TD (i.e. the combined value = ICU value)
#CSV.write("scripts/primary_care_v_icu/sim_TDs_v2.csv", sim_tds) # using 100 and 200 metagenomic samples per week
#CSV.write("scripts/primary_care_v_icu/sim_TDs_all_gp_samples_v2.csv", sim_tds) # using all samples for metagenomic sequencing each week
#CSV.write("scripts/primary_care_v_icu/sim_TDs_primary_care_p25_v2.csv", sim_tds) # metagenomic sampling of 25% of all ARI consultations

### Compute median values for each time to detection scenario
### Compute % of simulations that return a time to detection in each time to detection scenario

function analyse_columns(df::DataFrame)
    result = DataFrame(TD_description=String[], Median_TD=Float64[], Percentage_with_a_TD=Float64[])
    for col in names(df)
        col_data = df[!, col]
        # Convert strings "Inf" to actual Inf (if present as string)
        cleaned = [x == "Inf" ? Inf : x for x in col_data]
        # Keep only finite numbers
        finite_vals = filter(x -> isfinite(x), cleaned)
        percent_finite = 100 * length(finite_vals) / length(col_data)
        median_val = isempty(finite_vals) ? NaN : median(skipmissing(finite_vals))
        push!(result, (col, median_val, percent_finite))
    end
    return result
end


###### Analysis

# 1 # using 100 and 200 metagenomic samples per week
sim_tds_100_200_samples = icu_v_pc_td(; gp_swabs_mg = [100,200] ) #[319, 747] # Assumed number of swabs that are metagenomic sequenced for investigating impact
# Save as .csv file for inspection
CSV.write("scripts/primary_care_v_icu/sim_TDs_v2.csv", sim_tds_100_200_samples) # using 100 and 200 metagenomic samples per week
sim_tds_100_200_samples_analysis = analyse_columns(sim_tds_100_200_samples[:,2:19]) # Not essential but remove the column containing the simulation number
println(sim_tds_100_200_samples_analysis)

# 2 # using all samples for metagenomic sequencing each week
sim_tds_all_gp_samples = icu_v_pc_td(; gp_swabs_mg = [319, 747] ) #[100,200] # Assumed number of swabs that are metagenomic sequenced for investigating impact
# Save as .csv file for inspection
CSV.write("scripts/primary_care_v_icu/sim_TDs_all_gp_samples_v2.csv", sim_tds_all_gp_samples) # using 100 and 200 metagenomic samples per week
sim_tds_all_gp_samples_analysis = analyse_columns(sim_tds_all_gp_samples[:,2:19]) # Not essential but remove the column containing the simulation number
println(sim_tds_all_gp_samples_analysis)

# 3 # Sampling 25% of all ARI consultations
sim_tds_primary_care_p25 = icu_v_pc_td(; p_pc = 0.25 ) # Test to compare TDs with ICU with 0.25 sampling proportion # If p_pc not input then it is calculated from other inputs
# Save as .csv file for inspection
CSV.write("scripts/primary_care_v_icu/sim_TDs_primary_care_p25_v2.csv", sim_tds_primary_care_p25) # metagenomic sampling of 25% of all ARI consultations
sim_tds_primary_care_p25_analysis = analyse_columns(sim_tds_primary_care_p25[:,2:19]) # Not essential but remove the column containing the simulation number
println(sim_tds_primary_care_p25_analysis)

# 3(a) # Sampling 25% of all ARI consultations but using number of swabs instead of directly defining the proportion
# Same for all seasons - note that the differences between min and max are not meaningful in this case - it is just variance from random sampling
sim_tds_primary_care_25_summer = icu_v_pc_td(; gp_swabs_mg = [25698, 25698] #[100,200] # Assumed number of swabs that are metagenomic sequenced for investigating impact
                                             , gp_ari_swabs = [25698, 25698] #[319, 747] # Number of swabs taken from suspected ARI per week [mean summer 2024, mean winter 2024/25]. Source: Analysis of data extracted from RCGP Research & Surveillance Centre (RSC) Virology Dashboard [Accessed 29 Aug 2025]
                                            )
# Save as .csv file for inspection
CSV.write("scripts/primary_care_v_icu/sim_tds_primary_care_25_summer_v2.csv", sim_tds_primary_care_25_summer) # using 100 and 200 metagenomic samples per week
sim_tds_primary_care_25_summer_analysis = analyse_columns(sim_tds_primary_care_25_summer[:,2:19]) # Not essential but remove the column containing the simulation number
println(sim_tds_primary_care_25_summer_analysis)

# 4 # Assume hit the targetted values for the number of GP swabs (1000 per week or 20 swabs per practice x 300 practices = 6000) and using all samples for metagenomic sequencing each week
# Same for all seasons - note that the differences between min and max are not meaningful in this case - it is just variance from random sampling
sim_tds_higher_gp_swab_target_1000 = icu_v_pc_td(; gp_swabs_mg = [1000, 1000] #[100,200] # Assumed number of swabs that are metagenomic sequenced for investigating impact
                                             , gp_ari_swabs = [1000, 1000] #[319, 747] # Number of swabs taken from suspected ARI per week [mean summer 2024, mean winter 2024/25]. Source: Analysis of data extracted from RCGP Research & Surveillance Centre (RSC) Virology Dashboard [Accessed 29 Aug 2025]
                                            )
# Save as .csv file for inspection
CSV.write("scripts/primary_care_v_icu/sim_tds_higher_gp_swab_target_1000_v2.csv", sim_tds_higher_gp_swab_target_1000) # using 100 and 200 metagenomic samples per week
sim_tds_higher_gp_swab_target_1000_analysis = analyse_columns(sim_tds_higher_gp_swab_target_1000[:,2:19]) # Not essential but remove the column containing the simulation number
println(sim_tds_higher_gp_swab_target_1000_analysis)

# 5 # Assume hit the targetted values for the number of GP swabs (20 swabs per practice x 300 practices = 6000) and using all samples for metagenomic sequencing each week
# Same for all seasons - note that the differences between min and max are not meaningful in this case - it is just variance from random sampling
sim_tds_higher_gp_swab_target_6000 = icu_v_pc_td(; gp_swabs_mg = [6000, 6000] #[100,200] # Assumed number of swabs that are metagenomic sequenced for investigating impact
                                             , gp_ari_swabs = [6000, 6000] #[319, 747] # Number of swabs taken from suspected ARI per week [mean summer 2024, mean winter 2024/25]. Source: Analysis of data extracted from RCGP Research & Surveillance Centre (RSC) Virology Dashboard [Accessed 29 Aug 2025]
                                            )
# Save as .csv file for inspection
CSV.write("scripts/primary_care_v_icu/sim_tds_higher_gp_swab_target_6000_v2.csv", sim_tds_higher_gp_swab_target_6000) # using 100 and 200 metagenomic samples per week
sim_tds_higher_gp_swab_target_6000_analysis = analyse_columns(sim_tds_higher_gp_swab_target_6000[:,2:19]) # Not essential but remove the column containing the simulation number
println(sim_tds_higher_gp_swab_target_6000_analysis)

# Create compilation of data
# Times to detection
combined_df = DataFrame([zeros(Float64,7) for _ in 1:6], [:Number_of_PC_mg_samples
                                                        , :ICU_only
                                                        , :PC_only_summer
                                                        , :Combined_summer
                                                        , :PC_only_winter
                                                        , :Combined_winter])
#                    Number_of_PC_mg_samples  ICU_only                                          PC_only_summer                                      Combined_summer                                     PC_only_winter                                      Combined_winter
combined_df[1,:] = [ 100,                     sim_tds_100_200_samples_analysis[1,2],            sim_tds_100_200_samples_analysis[3,2],              sim_tds_100_200_samples_analysis[11,2],             sim_tds_100_200_samples_analysis[7,2],              sim_tds_100_200_samples_analysis[15,2]           ]
combined_df[2,:] = [ 200,                     sim_tds_100_200_samples_analysis[1,2],            sim_tds_100_200_samples_analysis[5,2],              sim_tds_100_200_samples_analysis[13,2],             sim_tds_100_200_samples_analysis[9,2],              sim_tds_100_200_samples_analysis[17,2]           ]                                
combined_df[3,:] = [ 319,                     sim_tds_all_gp_samples_analysis[1,2],             sim_tds_all_gp_samples_analysis[3,2],               sim_tds_all_gp_samples_analysis[11,2],              sim_tds_all_gp_samples_analysis[7,2],               sim_tds_all_gp_samples_analysis[15,2]            ]                     
combined_df[4,:] = [ 747,                     sim_tds_all_gp_samples_analysis[1,2],             sim_tds_all_gp_samples_analysis[5,2],               sim_tds_all_gp_samples_analysis[13,2],              sim_tds_all_gp_samples_analysis[9,2],               sim_tds_all_gp_samples_analysis[17,2]            ]
combined_df[5,:] = [ 1000,                    sim_tds_higher_gp_swab_target_1000_analysis[1,2], sim_tds_higher_gp_swab_target_1000_analysis[3,2],   sim_tds_higher_gp_swab_target_1000_analysis[11,2],  sim_tds_higher_gp_swab_target_1000_analysis[7,2],   sim_tds_higher_gp_swab_target_1000_analysis[15,2]]
combined_df[6,:] = [ 6000,                    sim_tds_higher_gp_swab_target_6000_analysis[1,2], sim_tds_higher_gp_swab_target_6000_analysis[3,2],   sim_tds_higher_gp_swab_target_6000_analysis[11,2],  sim_tds_higher_gp_swab_target_1000_analysis[7,2],   sim_tds_higher_gp_swab_target_6000_analysis[15,2]]
combined_df[7,:] = [ 0.25,                    sim_tds_primary_care_p25_analysis[1,2],           sim_tds_primary_care_p25_analysis[3,2],             sim_tds_primary_care_p25_analysis[11,2],            sim_tds_primary_care_p25_analysis[7,2],             sim_tds_primary_care_p25_analysis[15,2]          ]
combined_df[7,5] = combined_df[7,3]
combined_df[7,6] = combined_df[7,4]
println(combined_df)
# Difference between combined and ICU_only for summer
println( combined_df[:,4]-combined_df[:,2] )
# Difference between combined and ICU_only for winter
println( combined_df[:,6]-combined_df[:,2] )


# % of sim reps
combined_sim_rep_df = DataFrame([zeros(Float64,7) for _ in 1:6], [:Number_of_PC_mg_samples
                                                        , :ICU_only
                                                        , :PC_only_summer
                                                        , :Combined_summer
                                                        , :PC_only_winter
                                                        , :Combined_winter])
#                           Number_of_PC_mg_samples   ICU_only                                          PC_only_summer                                      Combined_summer                                     PC_only_winter                                      Combined_winter
combined_sim_rep_df[1,:] = [ 100,                     sim_tds_100_200_samples_analysis[1,3],            sim_tds_100_200_samples_analysis[3,3],              sim_tds_100_200_samples_analysis[11,3],             sim_tds_100_200_samples_analysis[7,3],              sim_tds_100_200_samples_analysis[15,3]           ]
combined_sim_rep_df[2,:] = [ 200,                     sim_tds_100_200_samples_analysis[1,3],            sim_tds_100_200_samples_analysis[5,3],              sim_tds_100_200_samples_analysis[13,3],             sim_tds_100_200_samples_analysis[9,3],              sim_tds_100_200_samples_analysis[17,3]           ]                                
combined_sim_rep_df[3,:] = [ 319,                     sim_tds_all_gp_samples_analysis[1,3],             sim_tds_all_gp_samples_analysis[3,3],               sim_tds_all_gp_samples_analysis[11,3],              sim_tds_all_gp_samples_analysis[7,3],               sim_tds_all_gp_samples_analysis[15,3]            ]                     
combined_sim_rep_df[4,:] = [ 747,                     sim_tds_all_gp_samples_analysis[1,3],             sim_tds_all_gp_samples_analysis[5,3],               sim_tds_all_gp_samples_analysis[13,3],              sim_tds_all_gp_samples_analysis[9,3],               sim_tds_all_gp_samples_analysis[17,3]            ]
combined_sim_rep_df[5,:] = [ 1000,                    sim_tds_higher_gp_swab_target_1000_analysis[1,3], sim_tds_higher_gp_swab_target_1000_analysis[3,3],   sim_tds_higher_gp_swab_target_1000_analysis[11,3],  sim_tds_higher_gp_swab_target_1000_analysis[7,3],   sim_tds_higher_gp_swab_target_1000_analysis[15,3]]
combined_sim_rep_df[6,:] = [ 6000,                    sim_tds_higher_gp_swab_target_6000_analysis[1,3], sim_tds_higher_gp_swab_target_6000_analysis[3,3],   sim_tds_higher_gp_swab_target_6000_analysis[11,3],  sim_tds_higher_gp_swab_target_1000_analysis[7,3],   sim_tds_higher_gp_swab_target_6000_analysis[15,3]]
combined_sim_rep_df[7,:] = [ 0.25,                    sim_tds_primary_care_p25_analysis[1,3],           sim_tds_primary_care_p25_analysis[3,3],             sim_tds_primary_care_p25_analysis[11,3],            sim_tds_primary_care_p25_analysis[7,3],             sim_tds_primary_care_p25_analysis[15,3]          ]
combined_sim_rep_df[7,5] = combined_sim_rep_df[7,3]
combined_sim_rep_df[7,6] = combined_sim_rep_df[7,4]
println(combined_sim_rep_df)



using DataFrames
using Pkg
Pkg.add("StatsPlots")
using StatsPlots  # This imports the @df macro for easier DataFrame plotting

# Plot the DataFrame with points and lines
# Summer
x_replace_summer = [100,200,319,747,1000,6000,25698]
# plot separately so coloured by y-series
p = plot(xlabel = "Number of primary care metagenomic samples"
        , ylabel = "Median time to detection of 1 case (TD) in days"
        , xscale=:log10
        , ylim=(0,70)
        , legend = :bottomleft
        )
@df combined_df scatter!(x_replace_summer
                        , :ICU_only
                        , color=:red
                        , marker=:circle
                        , markersize=4
                        , label="ICU only")
@df combined_df scatter!(x_replace_summer
                        , :PC_only_summer
                        , color=:blue
                        , marker=:circle
                        , markersize=4
                        , label="Primary care only (summer)")
@df combined_df scatter!(x_replace_summer
                        , :Combined_summer
                        , color=:green
                        , marker=:circle
                        , markersize=4
                        , label="Combined (summer)")
@df combined_df plot!(x_replace_summer, :ICU_only, color=:red, linewidth=2, label="")
@df combined_df plot!(x_replace_summer, :PC_only_summer, color=:blue, linewidth=2, label="")
@df combined_df plot!(x_replace_summer, :Combined_summer, color=:green, linewidth=2, label="")

# Winter
x_replace_winter = [100,200,319,747,1000,6000,46685]
# plot separately so coloured by y-series
#p = plot(xlabel = "Number of primary care metagenomic samples"
#        , ylabel = "Time to detection of 1 case (days)"
#        , xscale=:log10
#        , ylim=(0,70)
#        , legend = :topright
#        )
@df combined_df scatter!(x_replace_winter
                        , :ICU_only
                        , color=:red
                        , marker=:circle
                        , markersize=4
                        , label="")
@df combined_df scatter!(x_replace_winter
                        , :PC_only_winter
                        , color=:lightblue
                        , marker=:circle
                        , markersize=4
                        , label="Primary care only (winter)")
@df combined_df scatter!(x_replace_winter
                        , :Combined_winter
                        , color=:lightgreen
                        , marker=:circle
                        , markersize=4
                        , label="Combined (winter)")
@df combined_df plot!(x_replace_winter, :ICU_only, color=:red, linewidth=2, label="")
@df combined_df plot!(x_replace_winter, :PC_only_winter, color=:lightblue, linewidth=2, label="")
@df combined_df plot!(x_replace_winter, :Combined_winter, color=:lightgreen, linewidth=2, label="")

# Save to file
savefig("scripts/primary_care_v_icu/TD_vs_PC_sample_size.png")


# Analysis with altered tsample for ICU sampling
# Currently random sample from uniform distribution between ticu and trecovered but this can be quite large
# Try with random sample from uniform distribution between ticu and ticu+3

# 1 # using 100 and 200 metagenomic samples per week
sim_tds_100_200_samples_shorter_icu_tsample = icu_v_pc_td(; gp_swabs_mg = [100,200] ) #[319, 747] # Assumed number of swabs that are metagenomic sequenced for investigating impact
# Save as .csv file for inspection
CSV.write("scripts/primary_care_v_icu/sim_TDs_shorter_icu_tsample_v2.csv", sim_tds_100_200_samples_shorter_icu_tsample) # using 100 and 200 metagenomic samples per week
sim_tds_100_200_samples_shorter_icu_tsample_analysis = analyse_columns(sim_tds_100_200_samples_shorter_icu_tsample[:,2:19]) # Not essential but remove the column containing the simulation number
println(sim_tds_100_200_samples_analysis)


### looking at the number of ICU cases and GP cases per week
sims_G_gp_filter = load("covidlike-1.1.1-sims_filtered_G_gp.jld2", "sims_G_gp_filter")
sims_G_icu_filter = load("covidlike-1.1.1-sims_filtered_G_icu.jld2", "sims_G_icu_filter")

n_cases_df = DataFrame([zeros(Int,1000) for _ in 1:2], [:n_ICU_cases
                                                        , :n_GP_cases
                                                        ])

for s in 1:1000
    n_cases_df[s,1] = size(sims_G_icu_filter[s],1)
    n_cases_df[s,2] = size(sims_G_gp_filter[s],1)
end
println(n_cases_df)
median(n_cases_df[:,1])
median(n_cases_df[:,2])
n_cases_mat = hcat(n_cases_df[:,1],n_cases_df[:,2])
histogram(n_cases_mat, label=["ICU cases" "GP cases"], bins=100, alpha=0.6)