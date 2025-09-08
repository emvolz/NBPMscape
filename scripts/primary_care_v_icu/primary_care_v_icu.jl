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
#@load "covidlike-1.1.1-sims_filtered_G_icu_combined_nrep53000.jld2" sims_filtered_G_icu
#@load "covidlike-1.1.1-sims_filtered_G_gp_combined_nrep53000.jld2"  sims_filtered_G_gp

# Obtain population of England from constant in NBPMscape
# TODO not currently defined in Main
#ITL225CD_wales = ITL2SIZE[37:39,:ITL225CD]
#ITL2SIZE_eng = ITL2SIZE[.!in(ITL2SIZE.ITL225CD, ITL225CD_wales), :]
#ITL2SIZE_eng = ITL2SIZE[ ITL2SIZE.ITL225CD .!= , :]
#pop_eng = sum(ITL2SIZE[1:36,3])

function icu_v_pc_td(; p_icu = 0.15 # ICU sampling proportion TODO Assumption needs refining 
                     , icu_ari_admissions = [793,1440] # Weekly ICU admission numbers [summer,winter]
                     , icu_turnaround_time = [2,4] # Time to process sample and report results / declare detection
                       # Parameters for existing Oxford-RCGP RSC primary care surveillance
                     , gp_practices_total = 6199 # Total number of GP practices in England at July 2025. Source: BMA analysis (https://www.bma.org.uk/advice-and-support/nhs-delivery-and-workforce/pressures/pressures-in-general-practice-data-analysis) of NHS Digital General Practice Workforce Statistics (https://digital.nhs.uk/data-and-information/publications/statistical/general-and-personal-medical-services) [Accessed 2 Sep 2025]  
                     , gp_practices_swab = 300 # Number of GP practices taking swabs for virology surveillance. Source: Data quality report: national flu and COVID-19 surveillance report (27 May 2025)
                     , gp_swabs_mg = [100,200] #[319, 747]#[25698,46685] # Assumed number of swabs that are metagenomic sequenced for investigating impact
                     , pop_eng = 5.7106398e7 # Population of England. Source: ONS mid-year 2022 UK population data disaggregated by various geo levels and age groups and gender. [Accessed 6 November 2024] Available at https://www.ons.gov.uk/peoplepopulationandcommunity/populationandmigration/populationestimates/datasets/populationestimatesforukenglandandwalesscotlandandnorthernireland 
                     , gp_ari_consults = [180, 327] # Number of ARI consultations per 100k of England population per week [mean summer 2024, mean winter 2024/25]. Source: Analysis of data extracted from RCGP Research & Surveillance Centre (RSC) Virology Dashboard [Accessed 29 Aug 2025]
                     , gp_ari_swabs = [319, 747] #[25698,46685]# Number of swabs taken from suspected ARI per week [mean summer 2024, mean winter 2024/25]. Source: Analysis of data extracted from RCGP Research & Surveillance Centre (RSC) Virology Dashboard [Accessed 29 Aug 2025]
                     , pc_swab_turnaround_time = [2,4] # [min,max] number of days between swab sample being taken and results received. Source: Data quality report: national flu and COVID-19 surveillance report (27 May 2025).  Assume same for metagenomic testing.
                     #, p_pc # TODO Needs to be uncommented when want to set this directly
                     , sim_object # = "filtered" or  "full" 
                     )
    
    ### Probability of primary care surveillance events
    # Probability of visiting a General Practice that takes surveillance swabs
    prob_visiting_gp_swab = minimum( [1, gp_practices_swab / gp_practices_total ])
    # Probability of being swabbed given visited a GP that takes suriveillance swabs
    prob_swabbed_g_visit_gp_swab_summer = minimum([1, (1e5*gp_ari_swabs[1]) / (pop_eng*gp_ari_consults[1]) / prob_visiting_gp_swab])
    prob_swabbed_g_visit_gp_swab_winter = minimum([1, (1e5*gp_ari_swabs[2]) / (pop_eng*gp_ari_consults[2]) / prob_visiting_gp_swab])
    # Probability that swab is selected for metagenomic sequencing given that patient has been selected for a surveillance swab
    prob_mg_g_swabbed_min_summer = minimum([1, gp_swabs_mg[1] / gp_ari_swabs[1] ])
    prob_mg_g_swabbed_min_winter = minimum([1, gp_swabs_mg[1] / gp_ari_swabs[2] ])
    prob_mg_g_swabbed_max_summer = minimum([1, gp_swabs_mg[2] / gp_ari_swabs[1] ])
    prob_mg_g_swabbed_max_winter = minimum([1, gp_swabs_mg[2] / gp_ari_swabs[2] ])

    # Primary care sampling proportion - equivalent to probability a swab being taken and a metagenomic sequence being taken if attend a GP
    p_pc_mg_g_gp_visit_min_summer = prob_visiting_gp_swab * prob_swabbed_g_visit_gp_swab_summer * prob_mg_g_swabbed_min_summer
    p_pc_mg_g_gp_visit_min_winter = prob_visiting_gp_swab * prob_swabbed_g_visit_gp_swab_winter * prob_mg_g_swabbed_min_winter
    p_pc_mg_g_gp_visit_max_summer = prob_visiting_gp_swab * prob_swabbed_g_visit_gp_swab_summer * prob_mg_g_swabbed_max_summer
    p_pc_mg_g_gp_visit_max_winter = prob_visiting_gp_swab * prob_swabbed_g_visit_gp_swab_winter * prob_mg_g_swabbed_max_winter
    # Check if primary care sampling proportion, p_pc, was entered as an input to the function
    #if @isdefined p_pc 
        #continue
    #else # Define based on calculated probabilities
        p_pc = [ p_pc_mg_g_gp_visit_min_summer
                , p_pc_mg_g_gp_visit_max_summer
                , p_pc_mg_g_gp_visit_min_winter
                , p_pc_mg_g_gp_visit_max_winter ]
    #end    
    #println("p_pc = ",p_pc)
    #p_pc = 0.25 # Test to compare TDs with ICU with 0.25 sampling proportion
    #swabs_at_25_summer
    #swabs_at_25_winter

    # How many replicates are in this simulation
    if sim_object == "full"
        n_replicates = length(sims)
    elseif sim_object == "filtered"
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
    sim_tds_cols = [copy(col_any) for _ in 1:25]
    sim_tds = DataFrame( sim_tds_cols, ["sim_n"
                                        # Time to detection for 1 case
                                        ,"ICU_TD"                # ICU sampling only
                                        ,"PC_TD_min_summer"      # primary care sampling only with LOWER  number of metagenomic swab tests in SUMMER
                                        ,"PC_TD_max_summer"      # primary care sampling only with HIGHER number of metagenomic swab tests in SUMMER
                                        ,"PC_TD_min_winter"      # primary care sampling only with LOWER  number of metagenomic swab tests in WINTER
                                        ,"PC_TD_max_winter"      # primary care sampling only with HIGHER number of metagenomic swab tests in WINTER
                                        ,"ICU_PC_TD_min_summer"  # ICU AND primary care sampling  with LOWER  number of metagenomic swab tests in SUMMER
                                        ,"ICU_PC_TD_max_summer"  # ICU AND primary care sampling  with HIGHER number of metagenomic swab tests in SUMMER
                                        ,"ICU_PC_TD_min_winter"  # ICU AND primary care sampling  with LOWER  number of metagenomic swab tests in WINTER
                                        ,"ICU_PC_TD_max_winter"  # ICU AND primary care sampling  with HIGHER number of metagenomic swab tests in WINTER
                                        # Time to detection for 3 cases
                                        ,"ICU_3TD"                # ICU sampling only
                                        ,"PC_3TD_min_summer"      # primary care sampling only with LOWER  number of metagenomic swab tests in SUMMER
                                        ,"PC_3TD_max_summer"      # primary care sampling only with HIGHER number of metagenomic swab tests in SUMMER
                                        ,"PC_3TD_min_winter"      # primary care sampling only with LOWER  number of metagenomic swab tests in WINTER
                                        ,"PC_3TD_max_winter"      # primary care sampling only with HIGHER number of metagenomic swab tests in WINTER
                                        ,"ICU_PC_3TD_min_summer"  # ICU AND primary care sampling  with LOWER  number of metagenomic swab tests in SUMMER
                                        ,"ICU_PC_3TD_max_summer"  # ICU AND primary care sampling  with HIGHER number of metagenomic swab tests in SUMMER
                                        ,"ICU_PC_3TD_min_winter"  # ICU AND primary care sampling  with LOWER  number of metagenomic swab tests in WINTER
                                        ,"ICU_PC_3TD_max_winter"  # ICU AND primary care sampling  with HIGHER number of metagenomic swab tests in WINTER
                                        ,"n_ICU_cases","n_ICU_cases_sampled"
                                        ,"n_GP_cases","n_GP_cases_sampled"
                                        ,"ICU_simid","GP_simid"
                                        ])

    # Adaptation of sampleforest() function in 'core.jl' and 'median_TD_by_region.jl'
    for s in 1:n_replicates
        #TEST
        #println("$(s) of $(n_replicates) replicates")

        # Add simulation number to results df
        sim_tds[s,:sim_n] = s
        
        #s=1
        # If sim_object is not pre-filtered then need to filter for ICU and GP cases
        if sim_object == "full"
            # Filter for ICU cases 
            fo = sims[s]
            icu_cases = fo.G[ isfinite.(fo.G.ticu), : ]
            # Filter for GP cases (some of these will also become ICU cases)
            gp_cases = fo.G[ isfinite.(fo.G.tgp), : ]    
        elseif sim_object == "filtered"
            icu_cases = sims_G_icu_filter[s]
            gp_cases = sims_G_gp_filter[s]
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
                    sim_tds[s,2:23] = ["No Sim ID match" for _ in 2:23]
                    continue
                end
            end
        end

        ### Record infection time stats
        ## First for ICU cases
        # Sample sizes
        n_icu =  rand( Binomial( size(icu_cases,1), p_icu ) )
        # Subsample of ICU cases
        icu_cases_sub = icu_cases[sample( 1:size(icu_cases,1), n_icu, replace=false ), :]

        # Record number of cases in the sim rep and the number that were tested
        sim_tds[s,:n_ICU_cases] = size(icu_cases,1)
        sim_tds[s,:n_ICU_cases_sampled] = n_icu

        if size(icu_cases_sub,1) == 0
            
            # Record zero cases and zero sampled
            #sim_tds[s,:n_ICU_cases] = 0
            #sim_tds[s,:n_ICU_cases_sampled] = 0

            # Record Inf as time to detection for earliest:
            # - 1 ICU case
            sim_tds[s,:ICU_TD] = Inf
            # - 3 ICU cases
            sim_tds[s,:ICU_3TD] = Inf

        elseif size(icu_cases_sub,1) > 0 
            
            #= Sample time has uniform distribution between time of admission to ICU and time of recovery
                TODO THIS MAY NEED TO BE UPDATED FOR LATER VERSIONS WITH MORE COMPLEX CARE PATHWAYS
            =#
            # Generate sample times
            icu_tsample = map( g -> rand( Uniform( g.ticu[1], g.ticu[1]+3)) , eachrow(icu_cases_sub) )
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

        # Record time stats for GP cases using different sampling parameters 
        # Sample sizes
        n_gp = Vector{Int64}(undef,length(p_pc))
        for i in 1:length(p_pc)
            
            # Sample sizes based on probability of infected case having a metagenomic sample taken at the GP
            n_gp = rand( Binomial( size(gp_cases,1), p_pc[i] ) )
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

                #p_pc = [ p_pc_mg_min_summer, p_pc_mg_max_summer, p_pc_mg_min_winter, p_pc_mg_max_winter ]
                
                ## Add to results df for ICU times to detection
                # Record Inf times to detection if no GP cases sampled:
                #  1 primary care case for TD
                sim_tds[s, i + 2 ] = Inf
                #  3 primary care cases for 3TD
                sim_tds[s, i + 6 ] = Inf

                # ICU and primary care combined
                # If there are no GP cases then the combination of primary care and ICU is just ICU only (which was already entered into df above)
                # - 1 case
                sim_tds[s, i + 6] = sim_tds[s,:ICU_TD]
                # - 3 case
                sim_tds[s, 2*i + 15] = sim_tds[s,:ICU_3TD]
                
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
                sim_tds[s, i + 2 ] = gp_cases_sub_top3_td[1,:treport]
                #  3 primary care cases for 3TD
                if size(gp_cases_sub_top3_td,1) >= 3
                    sim_tds[s, i + 11 ] = gp_cases_sub_top3_td[3,:treport]
                else
                    sim_tds[s, i + 11 ] = Inf 
                end
                
                # Primary care COMBINED with ICU cases
                # - 1 case, either ICU or primary care
                sim_tds[s, i + 6] = icu_gp_cases_sub_top3_td[1,:treport]
                # - 3 case, either ICU or primary care
                if size(icu_gp_cases_sub_top3_td,1) >= 3
                    sim_tds[s, i + 15] = icu_gp_cases_sub_top3_td[3,:treport]
                else
                    sim_tds[s, i + 15] = Inf
                end

            end # End of if loop checking whether there are any GP cases sampled

        end # End of different primary care sampling parameters loop

    end # End of loop through different simulations
    return sim_tds
end # End of icu_v_pc_td function

function icu_v_pc_td_2(; p_icu = 0.15 # ICU sampling proportion TODO Assumption needs refining 
                        , icu_ari_admissions = 793 # 1440 # Weekly ICU admission numbers [summer,winter]
                        , icu_turnaround_time = [2,4] # Time to process sample and report results / declare detection
                       # Parameters for existing Oxford-RCGP RSC primary care surveillance
                        , gp_practices_total = 6199 # Total number of GP practices in England at July 2025. Source: BMA analysis (https://www.bma.org.uk/advice-and-support/nhs-delivery-and-workforce/pressures/pressures-in-general-practice-data-analysis) of NHS Digital General Practice Workforce Statistics (https://digital.nhs.uk/data-and-information/publications/statistical/general-and-personal-medical-services) [Accessed 2 Sep 2025]  
                        , gp_practices_swab = 300 # Number of GP practices taking swabs for virology surveillance. Source: Data quality report: national flu and COVID-19 surveillance report (27 May 2025)
                        , gp_swabs_mg = 100 # 200 [319, 747]#[25698,46685] # Assumed number of swabs that are metagenomic sequenced for investigating impact
                        , pop_eng = 5.7106398e7 # Population of England. Source: ONS mid-year 2022 UK population data disaggregated by various geo levels and age groups and gender. [Accessed 6 November 2024] Available at https://www.ons.gov.uk/peoplepopulationandcommunity/populationandmigration/populationestimates/datasets/populationestimatesforukenglandandwalesscotlandandnorthernireland 
                        , gp_ari_consults = 180 # 327 # Number of ARI consultations per 100k of England population per week [mean summer 2024, mean winter 2024/25]. Source: Analysis of data extracted from RCGP Research & Surveillance Centre (RSC) Virology Dashboard [Accessed 29 Aug 2025]
                        , gp_ari_swabs = 319 # 747] #[25698,46685]# Number of swabs taken from suspected ARI per week [mean summer 2024, mean winter 2024/25]. Source: Analysis of data extracted from RCGP Research & Surveillance Centre (RSC) Virology Dashboard [Accessed 29 Aug 2025]
                        , pc_swab_turnaround_time = [2,4] # [min,max] number of days between swab sample being taken and results received. Source: Data quality report: national flu and COVID-19 surveillance report (27 May 2025).  Assume same for metagenomic testing.
                        , sim_object # = "filtered" or  "full" 
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
    
    # How many replicates are in this simulation
    if sim_object == "full"
        n_replicates = length(sims)
    elseif sim_object == "filtered"
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

    # Adaptation of sampleforest() function in 'core.jl' and 'median_TD_by_region.jl'
    for s in 1:n_replicates
        #TEST
        #println("$(s) of $(n_replicates) replicates")

        # Add simulation number to results df
        sim_tds[s,:sim_n] = s
        
        #s=1
        # If sim_object is not pre-filtered then need to filter for ICU and GP cases
        if sim_object == "full"
            fo = sims[s]
            # Filter for ICU cases 
            icu_cases = fo.G[ isfinite.(fo.G.ticu), : ]
            # Filter for GP cases (some of these will also become ICU cases)
            gp_cases = fo.G[ isfinite.(fo.G.tgp), : ]    
        elseif sim_object == "filtered"
            icu_cases = sims_G_icu_filter[s]
            gp_cases = sims_G_gp_filter[s]
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
        # Sample sizes
        n_icu =  rand( Binomial( size(icu_cases,1), p_icu ) )
        # Subsample of ICU cases
        icu_cases_sub = icu_cases[sample( 1:size(icu_cases,1), n_icu, replace=false ), :]

        # Record number of cases in the sim rep and the number that were tested
        sim_tds[s,:n_ICU_cases] = size(icu_cases,1)
        sim_tds[s,:n_ICU_cases_sampled] = n_icu

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
            
            #= Sample time has uniform distribution between time of admission to ICU and time of recovery
                TODO THIS MAY NEED TO BE UPDATED FOR LATER VERSIONS WITH MORE COMPLEX CARE PATHWAYS
            =#
            # Generate sample times
            icu_tsample = map( g -> rand( Uniform( g.ticu[1], g.ticu[1]+3)) , eachrow(icu_cases_sub) )
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
        # Sample sizes
        # Sample sizes based on probability of infected case having a metagenomic sample taken at the GP
        n_gp = rand( Binomial( size(gp_cases,1), p_pc ) )
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

function clean_TD_results(df::DataFrame)
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

function analyse_columns(df::DataFrame)
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
println(result)


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
        histogram([finite_vals])    
        println(names(df)[col])
end

#### Simplified analysis with simplified function - 7 Sep 2025

# 100 mg samples - summer - current RCGP protocol
TDs_mg100_swab319_ari180 = icu_v_pc_td_2(; gp_swabs_mg = 100 , gp_ari_swabs = 319, gp_ari_consults = 180, sim_object = "filtered")
CSV.write("scripts/primary_care_v_icu/TDs_mg100_swab319_ari180.csv", TDs_mg100_swab319_ari180) 
TDs_mg100_swab319_ari180_analysis = analyse_columns(TDs_mg100_swab319_ari180[:,2:11]) # Not essential but remove the column containing the simulation number

# 100 mg samples - winter - current RCGP protocol
TDs_mg100_swab747_ari327 = icu_v_pc_td_2(; gp_swabs_mg = 100, gp_ari_swabs = 747, gp_ari_consults = 327, sim_object = "filtered")
CSV.write("scripts/primary_care_v_icu/TDs_mg100_swab747_ari327.csv", TDs_mg100_swab747_ari327) 
TDs_mg100_swab747_ari327_analysis = analyse_columns(TDs_mg100_swab747_ari327[:,2:11]) # Not essential but remove the column containing the simulation number

# 200 mg samples - summer - current RCGP protocol
TDs_mg200_swab319_ari180 = icu_v_pc_td_2(; gp_swabs_mg = 200 , gp_ari_swabs = 319, gp_ari_consults = 180, sim_object = "filtered")
CSV.write("scripts/primary_care_v_icu/TDs_mg200_swab319_ari180.csv", TDs_mg200_swab319_ari180) 
TDs_mg200_swab319_ari180_analysis = analyse_columns(TDs_mg200_swab319_ari180[:,2:11]) # Not essential but remove the column containing the simulation number

# 200 mg samples - winter - current RCGP protocol
TDs_mg200_swab747_ari327 = icu_v_pc_td_2(; gp_swabs_mg = 200, gp_ari_swabs = 747, gp_ari_consults = 327, sim_object = "filtered")
CSV.write("scripts/primary_care_v_icu/TDs_mg200_swab747_ari327.csv", TDs_mg200_swab747_ari327) 
TDs_mg200_swab747_ari327_analysis = analyse_columns(TDs_mg200_swab747_ari327[:,2:11]) # Not essential but remove the column containing the simulation number

# 319 mg samples - summer - current RCGP protocol
TDs_mg319_swab319_ari180 = icu_v_pc_td_2(; gp_swabs_mg = 319 , gp_ari_swabs = 319, gp_ari_consults = 180, sim_object = "filtered")
CSV.write("scripts/primary_care_v_icu/TDs_mg319_swab319_ari180.csv", TDs_mg319_swab319_ari180) 
TDs_mg319_swab319_ari180_analysis = analyse_columns(TDs_mg319_swab319_ari180[:,2:11]) # Not essential but remove the column containing the simulation number

# 319 mg samples - winter - current RCGP protocol
TDs_mg319_swab747_ari327 = icu_v_pc_td_2(; gp_swabs_mg = 319, gp_ari_swabs = 747, gp_ari_consults = 327, sim_object = "filtered")
CSV.write("scripts/primary_care_v_icu/TDs_mg319_swab747_ari327.csv", TDs_mg319_swab747_ari327) 
TDs_mg319_swab747_ari327_analysis = analyse_columns(TDs_mg319_swab747_ari327[:,2:11]) # Not essential but remove the column containing the simulation number

# 747 mg samples - summer - More swabs than current RCGP 
TDs_mg747_swab747_ari180 = icu_v_pc_td_2(; gp_swabs_mg = 747, gp_ari_swabs = 747, gp_ari_consults = 180, sim_object = "filtered")
CSV.write("scripts/primary_care_v_icu/TDs_mg747_swab747_ari180.csv", TDs_mg747_swab747_ari180) 
TDs_mg747_swab747_ari180_analysis = analyse_columns(TDs_mg747_swab747_ari180[:,2:11]) # Not essential but remove the column containing the simulation number

# 747 mg samples - winter - current RCGP 
TDs_mg747_swab747_ari327 = icu_v_pc_td_2(; gp_swabs_mg = 747, gp_ari_swabs = 747, gp_ari_consults = 327, sim_object = "filtered")
CSV.write("scripts/primary_care_v_icu/TDs_mg747_swab747_ari327.csv", TDs_mg747_swab747_ari327) 
TDs_mg747_swab747_ari327_analysis = analyse_columns(TDs_mg747_swab747_ari327[:,2:11]) # Not essential but remove the column containing the simulation number

# 1000 mg samples - summer - More swabs than current RCGP - Gu et al (2024) target
TDs_mg1000_swab1000_ari180 = icu_v_pc_td_2(; gp_swabs_mg = 1000, gp_ari_swabs = 1000, gp_ari_consults = 180, sim_object = "filtered")
CSV.write("scripts/primary_care_v_icu/TDs_mg1000_swab1000_ari180.csv", TDs_mg1000_swab1000_ari180) 
TDs_mg1000_swab1000_ari180_analysis = analyse_columns(TDs_mg1000_swab1000_ari180[:,2:11]) # Not essential but remove the column containing the simulation number

# 1000 mg samples - winter - More swabs than current RCGP - Gu et al (2024) target
TDs_mg1000_swab1000_ari327 = icu_v_pc_td_2(; gp_swabs_mg = 1000, gp_ari_swabs = 1000, gp_ari_consults = 327, sim_object = "filtered")
CSV.write("scripts/primary_care_v_icu/TDs_mg1000_swab1000_ari327.csv", TDs_mg1000_swab1000_ari327) 
TDs_mg1000_swab1000_ari327_analysis = analyse_columns(TDs_mg1000_swab1000_ari327[:,2:11]) # Not essential but remove the column containing the simulation number

# 6000 mg samples - summer - More swabs than current RCGP - Leston et al (2022) target
TDs_mg6000_swab6000_ari180 = icu_v_pc_td_2(; gp_swabs_mg = 6000, gp_ari_swabs = 6000, gp_ari_consults = 180, sim_object = "filtered")
CSV.write("scripts/primary_care_v_icu/TDs_mg6000_swab6000_ari180.csv", TDs_mg6000_swab6000_ari180) 
TDs_mg6000_swab6000_ari180_analysis = analyse_columns(TDs_mg6000_swab6000_ari180[:,2:11]) # Not essential but remove the column containing the simulation number

# 6000 mg samples - winter - More swabs than current RCGP - Leston et al (2022) target
TDs_mg6000_swab6000_ari327 = icu_v_pc_td_2(; gp_swabs_mg = 6000, gp_ari_swabs = 6000, gp_ari_consults = 327, sim_object = "filtered")
CSV.write("scripts/primary_care_v_icu/TDs_mg6000_swab6000_ari327.csv", TDs_mg6000_swab6000_ari327) 
TDs_mg6000_swab6000_ari327_analysis = analyse_columns(TDs_mg6000_swab6000_ari327[:,2:11]) # Not essential but remove the column containing the simulation number

# 15419 mg samples - summer - More swabs than current RCGP - 15% of summer ARI GP consultations
TDs_mg15419_swab15419_ari180 = icu_v_pc_td_2(; gp_swabs_mg = 15419, gp_ari_swabs = 15419, gp_ari_consults = 180, sim_object = "filtered")
CSV.write("scripts/primary_care_v_icu/TDs_mg15419_swab15419_ari180.csv", TDs_mg15419_swab15419_ari180) 
TDs_mg15419_swab15419_ari180_analysis = analyse_columns(TDs_mg15419_swab15419_ari180[:,2:11]) # Not essential but remove the column containing the simulation number

# 15419 mg samples - winter - More swabs than current RCGP - 15% of summer ARI GP consultations
TDs_mg15419_swab15419_ari327 = icu_v_pc_td_2(; gp_swabs_mg = 15419, gp_ari_swabs = 15419, gp_ari_consults = 327, sim_object = "filtered")
CSV.write("scripts/primary_care_v_icu/TDs_mg15419_swab15419_ari327.csv", TDs_mg15419_swab15419_ari327) 
TDs_mg15419_swab15419_ari327_analysis = analyse_columns(TDs_mg15419_swab15419_ari327[:,2:11]) # Not essential but remove the column containing the simulation number

# 28011 mg samples - summer - More swabs than current RCGP - 15% of winter ARI GP consultations
TDs_mg28011_swab28011_ari180 = icu_v_pc_td_2(; gp_swabs_mg = 28011, gp_ari_swabs = 28011, gp_ari_consults = 180, sim_object = "filtered")
CSV.write("scripts/primary_care_v_icu/TDs_mg28011_swab28011_ari180.csv", TDs_mg28011_swab28011_ari180) 
TDs_mg28011_swab28011_ari180_analysis = analyse_columns(TDs_mg28011_swab28011_ari180[:,2:11]) # Not essential but remove the column containing the simulation number

# 28011 mg samples - winter - More swabs than current RCGP - 15% of winter ARI GP consultations
TDs_mg28011_swab28011_ari327 = icu_v_pc_td_2(; gp_swabs_mg = 28011, gp_ari_swabs = 28011, gp_ari_consults = 327, sim_object = "filtered")
CSV.write("scripts/primary_care_v_icu/TDs_mg28011_swab28011_ari327.csv", TDs_mg28011_swab28011_ari327) 
TDs_mg28011_swab28011_ari327_analysis = analyse_columns(TDs_mg28011_swab28011_ari327[:,2:11]) # Not essential but remove the column containing the simulation number

# 186738 mg samples - summer - More swabs than current RCGP - 100% of winter ARI GP consultations
TDs_mg186738_swab186738_ari180 = icu_v_pc_td_2(; gp_swabs_mg = 186738, gp_ari_swabs = 186738, gp_ari_consults = 180, sim_object = "filtered")
CSV.write("scripts/primary_care_v_icu/TDs_mg186738_swab186738_ari180.csv", TDs_mg186738_swab186738_ari180) 
TDs_mg186738_swab186738_ari180_analysis = analyse_columns(TDs_mg186738_swab186738_ari180[:,2:11]) # Not essential but remove the column containing the simulation number

# 186738 mg samples - winter - More swabs than current RCGP - 100% of winter ARI GP consultations
TDs_mg186738_swab186738_ari327 = icu_v_pc_td_2(; gp_swabs_mg = 186738, gp_ari_swabs = 186738, gp_ari_consults = 327, sim_object = "filtered")
CSV.write("scripts/primary_care_v_icu/TDs_mg186738_swab186738_ari327.csv", TDs_mg186738_swab186738_ari327) 
TDs_mg186738_swab186738_ari327_analysis = analyse_columns(TDs_mg186738_swab186738_ari327[:,2:11]) # Not essential but remove the column containing the simulation number

plot_hist(df = TDs_mg6000_swab6000_ari180, col = 6)
plot_hist(df = TDs_mg186738_swab186738_ari180, col = 3)
median(TDs_mg186738_swab186738_ari327[:,3])

#### Note that there is convergence in the TD if mg samples and swabs are increased but the number of GPs swabbing is not increased. Prob of sampling is capped at 5%, because you are capturing all the cases at 5% of GPs but none at the 95% of GPs that are not swabbing.
#### Increase gp_practices_swab to match % of ARI consultations being mg sequenced
gp_practices_swab = 300

# 15419 mg samples - summer - More swabs than current RCGP - 15% of summer ARI GP consultations - number of swabbing GPs increased
TDs_mg15419_swab15419_ari180_gp929 = icu_v_pc_td_2(; gp_swabs_mg = 15419, gp_ari_swabs = 15419, gp_ari_consults = 180, gp_practices_swab = floor(0.15*6199), sim_object = "filtered")
CSV.write("scripts/primary_care_v_icu/TDs_mg15419_swab15419_ari180_gp929.csv", TDs_mg15419_swab15419_ari180_gp929) 
TDs_mg15419_swab15419_ari180_gp929_analysis = analyse_columns(TDs_mg15419_swab15419_ari180_gp929[:,2:11]) # Not essential but remove the column containing the simulation number

# 15419 mg samples - winter - More swabs than current RCGP - 15% of summer ARI GP consultations - number of swabbing GPs increased
TDs_mg15419_swab15419_ari327_gp929 = icu_v_pc_td_2(; gp_swabs_mg = 15419, gp_ari_swabs = 15419, gp_ari_consults = 327, gp_practices_swab = floor(0.15*6199), sim_object = "filtered")
CSV.write("scripts/primary_care_v_icu/TDs_mg15419_swab15419_ari327_gp929.csv", TDs_mg15419_swab15419_ari327_gp929) 
TDs_mg15419_swab15419_ari327_gp929_analysis = analyse_columns(TDs_mg15419_swab15419_ari327_gp929[:,2:11]) # Not essential but remove the column containing the simulation number

# 28011 mg samples - summer - More swabs than current RCGP - 15% of winter ARI GP consultations - number of swabbing GPs increased
TDs_mg28011_swab28011_ari180_gp929 = icu_v_pc_td_2(; gp_swabs_mg = 28011, gp_ari_swabs = 28011, gp_ari_consults = 180, gp_practices_swab = floor(0.15*6199), sim_object = "filtered")
CSV.write("scripts/primary_care_v_icu/TDs_mg28011_swab28011_ari180_gp929.csv", TDs_mg28011_swab28011_ari180_gp929) 
TDs_mg28011_swab28011_ari180_gp929_analysis = analyse_columns(TDs_mg28011_swab28011_ari180_gp929[:,2:11]) # Not essential but remove the column containing the simulation number

# 28011 mg samples - winter - More swabs than current RCGP - 15% of winter ARI GP consultations - number of swabbing GPs increased
TDs_mg28011_swab28011_ari327_gp929 = icu_v_pc_td_2(; gp_swabs_mg = 28011, gp_ari_swabs = 28011, gp_ari_consults = 327, gp_practices_swab = floor(0.15*6199), sim_object = "filtered")
CSV.write("scripts/primary_care_v_icu/TDs_mg28011_swab28011_ari327_gp929.csv", TDs_mg28011_swab28011_ari327_gp929) 
TDs_mg28011_swab28011_ari327_gp929_analysis = analyse_columns(TDs_mg28011_swab28011_ari327_gp929[:,2:11]) # Not essential but remove the column containing the simulation number

# 102792 mg samples - summer - More swabs than current RCGP - 100% of summer ARI GP consultations - number of swabbing GPs kept at 300
TDs_mg102792_swab102792_ari180_gp300 = icu_v_pc_td_2(; gp_swabs_mg = 102792, gp_ari_swabs = 102792, gp_ari_consults = 180, gp_practices_swab = 300, sim_object = "filtered")
CSV.write("scripts/primary_care_v_icu/TDs_mg102792_swab102792_ari180_gp300.csv", TDs_mg102792_swab102792_ari180_gp300) 
TDs_mg102792_swab102792_ari180_gp300_analysis = analyse_columns(TDs_mg102792_swab102792_ari180_gp300[:,2:11]) # Not essential but remove the column containing the simulation number

# 102792 mg samples - summer - More swabs than current RCGP - 100% of summer ARI GP consultations - number of swabbing GPs increased
TDs_mg102792_swab102792_ari180_gp6199 = icu_v_pc_td_2(; gp_swabs_mg = 102792, gp_ari_swabs = 102792, gp_ari_consults = 180, gp_practices_swab = 6199, sim_object = "filtered")
CSV.write("scripts/primary_care_v_icu/TDs_mg102792_swab102792_ari180_gp6199.csv", TDs_mg102792_swab102792_ari180_gp6199) 
TDs_mg102792_swab102792_ari180_gp6199_analysis = analyse_columns(TDs_mg102792_swab102792_ari180_gp6199[:,2:11]) # Not essential but remove the column containing the simulation number

# 102792 mg samples - winter - More swabs than current RCGP - 100% of summer ARI GP consultations - number of swabbing GPs kept at 300
TDs_mg102792_swab102792_ari180_gp300 = icu_v_pc_td_2(; gp_swabs_mg = 102792, gp_ari_swabs = 102792, gp_ari_consults = 327, gp_practices_swab = 300, sim_object = "filtered")
CSV.write("scripts/primary_care_v_icu/TDs_mg102792_swab102792_ari180_gp300.csv", TDs_mg102792_swab102792_ari180_gp300) 
TDs_mg102792_swab102792_ari180_gp300_analysis = analyse_columns(TDs_mg102792_swab102792_ari180_gp300[:,2:11]) # Not essential but remove the column containing the simulation number

# 102792 mg samples - winter - More swabs than current RCGP - 100% of summer ARI GP consultations - number of swabbing GPs kept at 300
TDs_mg102792_swab102792_ari327_gp300 = icu_v_pc_td_2(; gp_swabs_mg = 102792, gp_ari_swabs = 102792, gp_ari_consults = 327, gp_practices_swab = 300, sim_object = "filtered")
CSV.write("scripts/primary_care_v_icu/TDs_mg102792_swab102792_ari327_gp300.csv", TDs_mg102792_swab102792_ari327_gp300) 
TDs_mg102792_swab102792_ari327_gp300_analysis = analyse_columns(TDs_mg102792_swab102792_ari327_gp300[:,2:11]) # Not essential but remove the column containing the simulation number

# 102792 mg samples - winter - More swabs than current RCGP - 100% of summer ARI GP consultations - number of swabbing GPs increased
TDs_mg102792_swab102792_ari327_gp6199 = icu_v_pc_td_2(; gp_swabs_mg = 102792, gp_ari_swabs = 102792, gp_ari_consults = 327, gp_practices_swab = 6199, sim_object = "filtered")
CSV.write("scripts/primary_care_v_icu/TDs_mg102792_swab102792_ari327_gp6199.csv", TDs_mg102792_swab102792_ari327_gp6199) 
TDs_mg102792_swab102792_ari327_gp6199_analysis = analyse_columns(TDs_mg102792_swab102792_ari327_gp6199[:,2:11]) # Not essential but remove the column containing the simulation number

# 186738 mg samples - winter - More swabs than current RCGP - 100% of winter ARI GP consultations - number of swabbing GPs increased
TDs_mg186738_swab186738_ari327_gp6199 = icu_v_pc_td_2(; gp_swabs_mg = 186738, gp_ari_swabs = 186738, gp_ari_consults = 327, gp_practices_swab = 6199, sim_object = "filtered")
CSV.write("scripts/primary_care_v_icu/TDs_mg186738_swab186738_ari327_gp6199.csv", TDs_mg186738_swab186738_ari327_gp6199) 
TDs_mg186738_swab186738_ari327_gp6199_analysis = analyse_columns(TDs_mg186738_swab186738_ari327_gp6199[:,2:11]) # Not essential but remove the column containing the simulation number



# Create compilation of data
# Times to detection 1 case
simple_combined_df = DataFrame([zeros(Float64,10) for _ in 1:8], [:Number_of_PC_mg_samples
                                                        , :ICU_only
                                                        , :PC_only_summer
                                                        , :Combined_summer
                                                        , :Improvement_summer
                                                        , :PC_only_winter
                                                        , :Combined_winter
                                                        , :Improvement_winter])
#                         Number of
#                         PC mg samples  ICU_only                                               PC_only_summer                                      Combined_summer                                     Improvement summer                                                                                      PC_only_winter                                          Combined_winter                                         Improvement winter
simple_combined_df[1,:] = [ 100,         TDs_mg100_swab319_ari180_analysis[1,2],                TDs_mg100_swab319_ari180_analysis[2,2],             TDs_mg100_swab319_ari180_analysis[3,2],             TDs_mg100_swab319_ari180_analysis[3,2] - TDs_mg100_swab319_ari180_analysis[1,2],                        TDs_mg100_swab747_ari327_analysis[2,2],                 TDs_mg100_swab747_ari327_analysis[3,2],                 TDs_mg100_swab319_ari180_analysis[3,2] - TDs_mg100_swab747_ari327_analysis[1,2]           ]
simple_combined_df[2,:] = [ 200,         TDs_mg200_swab319_ari180_analysis[1,2],                TDs_mg200_swab319_ari180_analysis[2,2],             TDs_mg200_swab319_ari180_analysis[3,2],             TDs_mg200_swab319_ari180_analysis[3,2] - TDs_mg200_swab319_ari180_analysis[1,2],                        TDs_mg200_swab747_ari327_analysis[2,2],                 TDs_mg200_swab747_ari327_analysis[3,2],                 TDs_mg200_swab319_ari180_analysis[3,2] - TDs_mg200_swab747_ari327_analysis[1,2]           ]
simple_combined_df[3,:] = [ 319,         TDs_mg319_swab319_ari180_analysis[1,2],                TDs_mg319_swab319_ari180_analysis[2,2],             TDs_mg319_swab319_ari180_analysis[3,2],             TDs_mg319_swab319_ari180_analysis[3,2] - TDs_mg319_swab319_ari180_analysis[1,2],                        TDs_mg319_swab747_ari327_analysis[2,2],                 TDs_mg319_swab747_ari327_analysis[3,2],                 TDs_mg319_swab319_ari180_analysis[3,2] - TDs_mg319_swab747_ari327_analysis[1,2]           ]
simple_combined_df[4,:] = [ 747,         TDs_mg747_swab747_ari180_analysis[1,2],                TDs_mg747_swab747_ari180_analysis[2,2],             TDs_mg747_swab747_ari180_analysis[3,2],             TDs_mg747_swab747_ari180_analysis[3,2] - TDs_mg747_swab747_ari180_analysis[1,2],                        TDs_mg747_swab747_ari327_analysis[2,2],                 TDs_mg747_swab747_ari327_analysis[3,2],                 TDs_mg747_swab747_ari327_analysis[3,2] - TDs_mg747_swab747_ari327_analysis[1,2]           ]
simple_combined_df[5,:] = [ 1000,        TDs_mg1000_swab1000_ari180_analysis[1,2],              TDs_mg1000_swab1000_ari180_analysis[2,2],           TDs_mg1000_swab1000_ari180_analysis[3,2],           TDs_mg1000_swab1000_ari180_analysis[3,2] - TDs_mg1000_swab1000_ari180_analysis[1,2],                    TDs_mg1000_swab1000_ari327_analysis[2,2],                 TDs_mg1000_swab1000_ari327_analysis[3,2],             TDs_mg1000_swab1000_ari327_analysis[3,2] - TDs_mg1000_swab1000_ari327_analysis[1,2]           ]
simple_combined_df[6,:] = [ 6000,        TDs_mg6000_swab6000_ari180_analysis[1,2],              TDs_mg6000_swab6000_ari180_analysis[2,2],           TDs_mg6000_swab6000_ari180_analysis[3,2],           TDs_mg6000_swab6000_ari180_analysis[3,2] - TDs_mg6000_swab6000_ari180_analysis[1,2],                    TDs_mg6000_swab6000_ari327_analysis[2,2],                TDs_mg6000_swab6000_ari327_analysis[3,2],              TDs_mg6000_swab6000_ari327_analysis[3,2] - TDs_mg6000_swab6000_ari327_analysis[1,2]           ]
simple_combined_df[7,:] = [ 15419,       TDs_mg15419_swab15419_ari180_gp929_analysis[1,2],      TDs_mg15419_swab15419_ari180_gp929_analysis[2,2],   TDs_mg15419_swab15419_ari180_gp929_analysis[3,2],   TDs_mg15419_swab15419_ari180_gp929_analysis[3,2] - TDs_mg15419_swab15419_ari180_gp929_analysis[1,2],    TDs_mg15419_swab15419_ari327_gp929_analysis[2,2],        TDs_mg15419_swab15419_ari327_gp929_analysis[3,2],      TDs_mg15419_swab15419_ari327_gp929_analysis[3,2] - TDs_mg15419_swab15419_ari327_gp929_analysis[1,2]           ]
simple_combined_df[8,:] = [ 28011,       TDs_mg28011_swab28011_ari180_gp929_analysis[1,2],      TDs_mg28011_swab28011_ari180_gp929_analysis[2,2],   TDs_mg28011_swab28011_ari180_gp929_analysis[3,2],   TDs_mg28011_swab28011_ari180_gp929_analysis[3,2] - TDs_mg28011_swab28011_ari180_gp929_analysis[1,2],    TDs_mg28011_swab28011_ari327_gp929_analysis[2,2],       TDs_mg28011_swab28011_ari327_gp929_analysis[3,2],       TDs_mg28011_swab28011_ari327_gp929_analysis[3,2] - TDs_mg28011_swab28011_ari327_gp929_analysis[1,2]           ]
simple_combined_df[9,:] = [ 102792,     TDs_mg102792_swab102792_ari180_gp6199_analysis[1,2],   TDs_mg102792_swab102792_ari180_gp6199_analysis[2,2],TDs_mg102792_swab102792_ari180_gp6199_analysis[3,2] , TDs_mg102792_swab102792_ari180_gp6199_analysis[3,2] - TDs_mg102792_swab102792_ari180_gp6199_analysis[1,2],  TDs_mg102792_swab102792_ari327_gp6199_analysis[2,2],    TDs_mg102792_swab102792_ari327_gp6199_analysis[3,2],    TDs_mg102792_swab102792_ari327_gp6199_analysis[3,2] - TDs_mg102792_swab102792_ari327_gp6199_analysis[1,2]           ]
simple_combined_df[10,:] = [ 186738,    TDs_mg186738_swab186738_ari327_gp6199_analysis[1,2],   -1,-1 , -1,  TDs_mg186738_swab186738_ari327_gp6199_analysis[2,2],    TDs_mg186738_swab186738_ari327_gp6199_analysis[3,2],    TDs_mg186738_swab186738_ari327_gp6199_analysis[3,2] - TDs_mg186738_swab186738_ari327_gp6199_analysis[1,2]           ]

println(simple_combined_df)

# Times to detection 3 cases
simple_combined_3TD_df = DataFrame([zeros(Float64,10) for _ in 1:8], [:Number_of_PC_mg_samples
                                                        , :ICU_only
                                                        , :PC_only_summer
                                                        , :Combined_summer
                                                        , :Improvement_summer
                                                        , :PC_only_winter
                                                        , :Combined_winter
                                                        , :Improvement_winter])
#                         Number of
#                         PC mg samples  ICU_only                                               PC_only_summer                                      Combined_summer                                     Improvement summer                                                                                          PC_only_winter                                          Combined_winter                                         Improvement winter
simple_combined_3TD_df[1,:] = [ 100,         TDs_mg100_swab319_ari180_analysis[4,2],                TDs_mg100_swab319_ari180_analysis[5,2],             TDs_mg100_swab319_ari180_analysis[6,2],             TDs_mg100_swab319_ari180_analysis[6,2] - TDs_mg100_swab319_ari180_analysis[4,2],                            TDs_mg100_swab747_ari327_analysis[5,2],                 TDs_mg100_swab747_ari327_analysis[6,2],                 TDs_mg100_swab319_ari180_analysis[6,2] - TDs_mg100_swab747_ari327_analysis[4,2]           ]
simple_combined_3TD_df[2,:] = [ 200,         TDs_mg200_swab319_ari180_analysis[4,2],                TDs_mg200_swab319_ari180_analysis[5,2],             TDs_mg200_swab319_ari180_analysis[6,2],             TDs_mg200_swab319_ari180_analysis[6,2] - TDs_mg200_swab319_ari180_analysis[4,2],                            TDs_mg200_swab747_ari327_analysis[5,2],                 TDs_mg200_swab747_ari327_analysis[6,2],                 TDs_mg200_swab319_ari180_analysis[6,2] - TDs_mg200_swab747_ari327_analysis[4,2]           ]
simple_combined_3TD_df[3,:] = [ 319,         TDs_mg319_swab319_ari180_analysis[4,2],                TDs_mg319_swab319_ari180_analysis[5,2],             TDs_mg319_swab319_ari180_analysis[6,2],             TDs_mg319_swab319_ari180_analysis[6,2] - TDs_mg319_swab319_ari180_analysis[4,2],                            TDs_mg319_swab747_ari327_analysis[5,2],                 TDs_mg319_swab747_ari327_analysis[6,2],                 TDs_mg319_swab319_ari180_analysis[6,2] - TDs_mg319_swab747_ari327_analysis[4,2]           ]
simple_combined_3TD_df[4,:] = [ 747,         TDs_mg747_swab747_ari180_analysis[4,2],                TDs_mg747_swab747_ari180_analysis[5,2],             TDs_mg747_swab747_ari180_analysis[6,2],             TDs_mg747_swab747_ari180_analysis[6,2] - TDs_mg747_swab747_ari180_analysis[4,2],                            TDs_mg747_swab747_ari327_analysis[5,2],                 TDs_mg747_swab747_ari327_analysis[6,2],                 TDs_mg747_swab747_ari327_analysis[6,2] - TDs_mg747_swab747_ari327_analysis[4,2]           ]
simple_combined_3TD_df[5,:] = [ 1000,        TDs_mg1000_swab1000_ari180_analysis[4,2],              TDs_mg1000_swab1000_ari180_analysis[5,2],           TDs_mg1000_swab1000_ari180_analysis[6,2],           TDs_mg1000_swab1000_ari180_analysis[6,2] - TDs_mg1000_swab1000_ari180_analysis[4,2],                        TDs_mg1000_swab1000_ari327_analysis[5,2],                 TDs_mg1000_swab1000_ari327_analysis[6,2],             TDs_mg1000_swab1000_ari327_analysis[6,2] - TDs_mg1000_swab1000_ari327_analysis[4,2]           ]
simple_combined_3TD_df[6,:] = [ 6000,        TDs_mg6000_swab6000_ari180_analysis[4,2],              TDs_mg6000_swab6000_ari180_analysis[5,2],           TDs_mg6000_swab6000_ari180_analysis[6,2],           TDs_mg6000_swab6000_ari180_analysis[6,2] - TDs_mg6000_swab6000_ari180_analysis[4,2],                        TDs_mg6000_swab6000_ari327_analysis[5,2],                TDs_mg6000_swab6000_ari327_analysis[6,2],              TDs_mg6000_swab6000_ari327_analysis[6,2] - TDs_mg6000_swab6000_ari327_analysis[4,2]           ]
simple_combined_3TD_df[7,:] = [ 15419,       TDs_mg15419_swab15419_ari180_gp929_analysis[4,2],      TDs_mg15419_swab15419_ari180_gp929_analysis[5,2],   TDs_mg15419_swab15419_ari180_gp929_analysis[6,2],   TDs_mg15419_swab15419_ari180_gp929_analysis[6,2] - TDs_mg15419_swab15419_ari180_gp929_analysis[4,2],        TDs_mg15419_swab15419_ari327_gp929_analysis[5,2],        TDs_mg15419_swab15419_ari327_gp929_analysis[6,2],      TDs_mg15419_swab15419_ari327_gp929_analysis[6,2] - TDs_mg15419_swab15419_ari327_gp929_analysis[4,2]           ]
simple_combined_3TD_df[8,:] = [ 28011,       TDs_mg28011_swab28011_ari180_gp929_analysis[4,2],      TDs_mg28011_swab28011_ari180_gp929_analysis[5,2],   TDs_mg28011_swab28011_ari180_gp929_analysis[6,2],   TDs_mg28011_swab28011_ari180_gp929_analysis[6,2] - TDs_mg28011_swab28011_ari180_gp929_analysis[4,2],        TDs_mg28011_swab28011_ari327_gp929_analysis[5,2],       TDs_mg28011_swab28011_ari327_gp929_analysis[6,2],       TDs_mg28011_swab28011_ari327_gp929_analysis[6,2] - TDs_mg28011_swab28011_ari327_gp929_analysis[4,2]           ]
simple_combined_3TD_df[9,:] = [ 102792,     TDs_mg102792_swab102792_ari180_gp6199_analysis[4,2],   TDs_mg102792_swab102792_ari180_gp6199_analysis[5,2], TDs_mg102792_swab102792_ari180_gp6199_analysis[6,2] , TDs_mg102792_swab102792_ari180_gp6199_analysis[6,2] - TDs_mg102792_swab102792_ari180_gp6199_analysis[4,2],  TDs_mg102792_swab102792_ari327_gp6199_analysis[5,2],    TDs_mg102792_swab102792_ari327_gp6199_analysis[6,2],    TDs_mg102792_swab102792_ari327_gp6199_analysis[6,2] - TDs_mg102792_swab102792_ari327_gp6199_analysis[4,2]           ]
simple_combined_3TD_df[10,:] = [ 186738,    TDs_mg186738_swab186738_ari327_gp6199_analysis[4,2],   -1,                                                  -1 ,                                                -1,                                                                                                         TDs_mg186738_swab186738_ari327_gp6199_analysis[5,2],    TDs_mg186738_swab186738_ari327_gp6199_analysis[6,2],    TDs_mg186738_swab186738_ari327_gp6199_analysis[6,2] - TDs_mg186738_swab186738_ari327_gp6199_analysis[4,2]           ]

println(simple_combined_3TD_df)

# % of replicates with times to detection 1 case
simple_combined_simrep_perc_df = DataFrame([zeros(Float64,10) for _ in 1:8], [:Number_of_PC_mg_samples
                                                        , :ICU_only
                                                        , :PC_only_summer
                                                        , :Combined_summer
                                                        , :Improvement_summer
                                                        , :PC_only_winter
                                                        , :Combined_winter
                                                        , :Improvement_winter])
#                         Number of
#                         PC mg samples  ICU_only                                               PC_only_summer                                      Combined_summer                                                 Improvement summer                                                                                      PC_only_winter                                          Combined_winter                                         Improvement winter
simple_combined_simrep_perc_df[1,:] = [ 100,         TDs_mg100_swab319_ari180_analysis[1,3],                TDs_mg100_swab319_ari180_analysis[2,3],             TDs_mg100_swab319_ari180_analysis[3,3],             TDs_mg100_swab319_ari180_analysis[3,3] - TDs_mg100_swab319_ari180_analysis[1,3],                        TDs_mg100_swab747_ari327_analysis[2,3],                 TDs_mg100_swab747_ari327_analysis[3,3],                 TDs_mg100_swab319_ari180_analysis[3,3] - TDs_mg100_swab747_ari327_analysis[1,3]           ]
simple_combined_simrep_perc_df[2,:] = [ 200,         TDs_mg200_swab319_ari180_analysis[1,3],                TDs_mg200_swab319_ari180_analysis[2,3],             TDs_mg200_swab319_ari180_analysis[3,3],             TDs_mg200_swab319_ari180_analysis[3,3] - TDs_mg200_swab319_ari180_analysis[1,3],                        TDs_mg200_swab747_ari327_analysis[2,3],                 TDs_mg200_swab747_ari327_analysis[3,3],                 TDs_mg200_swab319_ari180_analysis[3,3] - TDs_mg200_swab747_ari327_analysis[1,3]           ]
simple_combined_simrep_perc_df[3,:] = [ 319,         TDs_mg319_swab319_ari180_analysis[1,3],                TDs_mg319_swab319_ari180_analysis[2,3],             TDs_mg319_swab319_ari180_analysis[3,3],             TDs_mg319_swab319_ari180_analysis[3,3] - TDs_mg319_swab319_ari180_analysis[1,3],                        TDs_mg319_swab747_ari327_analysis[2,3],                 TDs_mg319_swab747_ari327_analysis[3,3],                 TDs_mg319_swab319_ari180_analysis[3,3] - TDs_mg319_swab747_ari327_analysis[1,3]           ]
simple_combined_simrep_perc_df[4,:] = [ 747,         TDs_mg747_swab747_ari180_analysis[1,3],                TDs_mg747_swab747_ari180_analysis[2,3],             TDs_mg747_swab747_ari180_analysis[3,3],             TDs_mg747_swab747_ari180_analysis[3,3] - TDs_mg747_swab747_ari180_analysis[1,3],                        TDs_mg747_swab747_ari327_analysis[2,3],                 TDs_mg747_swab747_ari327_analysis[3,3],                 TDs_mg747_swab747_ari327_analysis[3,3] - TDs_mg747_swab747_ari327_analysis[1,3]           ]
simple_combined_simrep_perc_df[5,:] = [ 1000,        TDs_mg1000_swab1000_ari180_analysis[1,3],              TDs_mg1000_swab1000_ari180_analysis[2,3],           TDs_mg1000_swab1000_ari180_analysis[3,3],           TDs_mg1000_swab1000_ari180_analysis[3,3] - TDs_mg1000_swab1000_ari180_analysis[1,3],                    TDs_mg1000_swab1000_ari327_analysis[2,3],                 TDs_mg1000_swab1000_ari327_analysis[3,3],             TDs_mg1000_swab1000_ari327_analysis[3,3] - TDs_mg1000_swab1000_ari327_analysis[1,3]           ]
simple_combined_simrep_perc_df[6,:] = [ 6000,        TDs_mg6000_swab6000_ari180_analysis[1,3],              TDs_mg6000_swab6000_ari180_analysis[2,3],           TDs_mg6000_swab6000_ari180_analysis[3,3],           TDs_mg6000_swab6000_ari180_analysis[3,3] - TDs_mg6000_swab6000_ari180_analysis[1,3],                    TDs_mg6000_swab6000_ari327_analysis[2,3],                TDs_mg6000_swab6000_ari327_analysis[3,3],              TDs_mg6000_swab6000_ari327_analysis[3,3] - TDs_mg6000_swab6000_ari327_analysis[1,3]           ]
simple_combined_simrep_perc_df[7,:] = [ 15419,       TDs_mg15419_swab15419_ari180_gp929_analysis[1,3],      TDs_mg15419_swab15419_ari180_gp929_analysis[2,3],   TDs_mg15419_swab15419_ari180_gp929_analysis[3,3],   TDs_mg15419_swab15419_ari180_gp929_analysis[3,3] - TDs_mg15419_swab15419_ari180_gp929_analysis[1,3],    TDs_mg15419_swab15419_ari327_gp929_analysis[2,3],        TDs_mg15419_swab15419_ari327_gp929_analysis[3,3],      TDs_mg15419_swab15419_ari327_gp929_analysis[3,3] - TDs_mg15419_swab15419_ari327_gp929_analysis[1,3]           ]
simple_combined_simrep_perc_df[8,:] = [ 28011,       TDs_mg28011_swab28011_ari180_gp929_analysis[1,3],      TDs_mg28011_swab28011_ari180_gp929_analysis[2,3],   TDs_mg28011_swab28011_ari180_gp929_analysis[3,3],   TDs_mg28011_swab28011_ari180_gp929_analysis[3,3] - TDs_mg28011_swab28011_ari180_gp929_analysis[1,3],    TDs_mg28011_swab28011_ari327_gp929_analysis[2,3],       TDs_mg28011_swab28011_ari327_gp929_analysis[3,3],       TDs_mg28011_swab28011_ari327_gp929_analysis[3,3] - TDs_mg28011_swab28011_ari327_gp929_analysis[1,3]           ]
simple_combined_simrep_perc_df[9,:] = [ 102792,     TDs_mg102792_swab102792_ari180_gp6199_analysis[1,3],   TDs_mg102792_swab102792_ari180_gp6199_analysis[2,3],TDs_mg102792_swab102792_ari180_gp6199_analysis[3,3] , TDs_mg102792_swab102792_ari180_gp6199_analysis[3,3] - TDs_mg102792_swab102792_ari180_gp6199_analysis[1,3],  TDs_mg102792_swab102792_ari327_gp6199_analysis[2,3],    TDs_mg102792_swab102792_ari327_gp6199_analysis[3,3],    TDs_mg102792_swab102792_ari327_gp6199_analysis[3,3] - TDs_mg102792_swab102792_ari327_gp6199_analysis[1,3]           ]
simple_combined_simrep_perc_df[10,:] = [ 186738,    TDs_mg186738_swab186738_ari327_gp6199_analysis[1,3],   -1,-1 , -1,  TDs_mg186738_swab186738_ari327_gp6199_analysis[2,3],    TDs_mg186738_swab186738_ari327_gp6199_analysis[3,3],    TDs_mg186738_swab186738_ari327_gp6199_analysis[3,3] - TDs_mg186738_swab186738_ari327_gp6199_analysis[1,3]           ]

println(simple_combined_simrep_perc_df)

# % of replicates with times to detection 3 cases
simple_combined_3TD_simrep_perc_df = DataFrame([zeros(Float64,10) for _ in 1:8], [:Number_of_PC_mg_samples
                                                        , :ICU_only
                                                        , :PC_only_summer
                                                        , :Combined_summer
                                                        , :Improvement_summer
                                                        , :PC_only_winter
                                                        , :Combined_winter
                                                        , :Improvement_winter])
#                         Number of
#                         PC mg samples  ICU_only                                               PC_only_summer                                      Combined_summer                                     Improvement summer                                                                                          PC_only_winter                                          Combined_winter                                         Improvement winter
simple_combined_3TD_simrep_perc_df[1,:] = [ 100,         TDs_mg100_swab319_ari180_analysis[4,3],                TDs_mg100_swab319_ari180_analysis[5,3],             TDs_mg100_swab319_ari180_analysis[6,3],             TDs_mg100_swab319_ari180_analysis[6,3] - TDs_mg100_swab319_ari180_analysis[4,3],                            TDs_mg100_swab747_ari327_analysis[5,3],                 TDs_mg100_swab747_ari327_analysis[6,3],                 TDs_mg100_swab319_ari180_analysis[6,3] - TDs_mg100_swab747_ari327_analysis[4,3]           ]
simple_combined_3TD_simrep_perc_df[2,:] = [ 200,         TDs_mg200_swab319_ari180_analysis[4,3],                TDs_mg200_swab319_ari180_analysis[5,3],             TDs_mg200_swab319_ari180_analysis[6,3],             TDs_mg200_swab319_ari180_analysis[6,3] - TDs_mg200_swab319_ari180_analysis[4,3],                            TDs_mg200_swab747_ari327_analysis[5,3],                 TDs_mg200_swab747_ari327_analysis[6,3],                 TDs_mg200_swab319_ari180_analysis[6,3] - TDs_mg200_swab747_ari327_analysis[4,3]           ]
simple_combined_3TD_simrep_perc_df[3,:] = [ 319,         TDs_mg319_swab319_ari180_analysis[4,3],                TDs_mg319_swab319_ari180_analysis[5,3],             TDs_mg319_swab319_ari180_analysis[6,3],             TDs_mg319_swab319_ari180_analysis[6,3] - TDs_mg319_swab319_ari180_analysis[4,3],                            TDs_mg319_swab747_ari327_analysis[5,3],                 TDs_mg319_swab747_ari327_analysis[6,3],                 TDs_mg319_swab319_ari180_analysis[6,3] - TDs_mg319_swab747_ari327_analysis[4,3]           ]
simple_combined_3TD_simrep_perc_df[4,:] = [ 747,         TDs_mg747_swab747_ari180_analysis[4,3],                TDs_mg747_swab747_ari180_analysis[5,3],             TDs_mg747_swab747_ari180_analysis[6,3],             TDs_mg747_swab747_ari180_analysis[6,3] - TDs_mg747_swab747_ari180_analysis[4,3],                            TDs_mg747_swab747_ari327_analysis[5,3],                 TDs_mg747_swab747_ari327_analysis[6,3],                 TDs_mg747_swab747_ari327_analysis[6,3] - TDs_mg747_swab747_ari327_analysis[4,3]           ]
simple_combined_3TD_simrep_perc_df[5,:] = [ 1000,        TDs_mg1000_swab1000_ari180_analysis[4,3],              TDs_mg1000_swab1000_ari180_analysis[5,3],           TDs_mg1000_swab1000_ari180_analysis[6,3],           TDs_mg1000_swab1000_ari180_analysis[6,3] - TDs_mg1000_swab1000_ari180_analysis[4,3],                        TDs_mg1000_swab1000_ari327_analysis[5,3],                 TDs_mg1000_swab1000_ari327_analysis[6,3],             TDs_mg1000_swab1000_ari327_analysis[6,3] - TDs_mg1000_swab1000_ari327_analysis[4,3]           ]
simple_combined_3TD_simrep_perc_df[6,:] = [ 6000,        TDs_mg6000_swab6000_ari180_analysis[4,3],              TDs_mg6000_swab6000_ari180_analysis[5,3],           TDs_mg6000_swab6000_ari180_analysis[6,3],           TDs_mg6000_swab6000_ari180_analysis[6,3] - TDs_mg6000_swab6000_ari180_analysis[4,3],                        TDs_mg6000_swab6000_ari327_analysis[5,3],                TDs_mg6000_swab6000_ari327_analysis[6,3],              TDs_mg6000_swab6000_ari327_analysis[6,3] - TDs_mg6000_swab6000_ari327_analysis[4,3]           ]
simple_combined_3TD_simrep_perc_df[7,:] = [ 15419,       TDs_mg15419_swab15419_ari180_gp929_analysis[4,3],      TDs_mg15419_swab15419_ari180_gp929_analysis[5,3],   TDs_mg15419_swab15419_ari180_gp929_analysis[6,3],   TDs_mg15419_swab15419_ari180_gp929_analysis[6,3] - TDs_mg15419_swab15419_ari180_gp929_analysis[4,3],        TDs_mg15419_swab15419_ari327_gp929_analysis[5,3],        TDs_mg15419_swab15419_ari327_gp929_analysis[6,3],      TDs_mg15419_swab15419_ari327_gp929_analysis[6,3] - TDs_mg15419_swab15419_ari327_gp929_analysis[4,3]           ]
simple_combined_3TD_simrep_perc_df[8,:] = [ 28011,       TDs_mg28011_swab28011_ari180_gp929_analysis[4,3],      TDs_mg28011_swab28011_ari180_gp929_analysis[5,3],   TDs_mg28011_swab28011_ari180_gp929_analysis[6,3],   TDs_mg28011_swab28011_ari180_gp929_analysis[6,3] - TDs_mg28011_swab28011_ari180_gp929_analysis[4,3],        TDs_mg28011_swab28011_ari327_gp929_analysis[5,3],       TDs_mg28011_swab28011_ari327_gp929_analysis[6,3],       TDs_mg28011_swab28011_ari327_gp929_analysis[6,3] - TDs_mg28011_swab28011_ari327_gp929_analysis[4,3]           ]
simple_combined_3TD_simrep_perc_df[9,:] = [ 102792,     TDs_mg102792_swab102792_ari180_gp6199_analysis[4,3],   TDs_mg102792_swab102792_ari180_gp6199_analysis[5,3], TDs_mg102792_swab102792_ari180_gp6199_analysis[6,3] , TDs_mg102792_swab102792_ari180_gp6199_analysis[6,3] - TDs_mg102792_swab102792_ari180_gp6199_analysis[4,3],  TDs_mg102792_swab102792_ari327_gp6199_analysis[5,3],    TDs_mg102792_swab102792_ari327_gp6199_analysis[6,3],    TDs_mg102792_swab102792_ari327_gp6199_analysis[6,3] - TDs_mg102792_swab102792_ari327_gp6199_analysis[4,3]           ]
simple_combined_3TD_simrep_perc_df[10,:] = [ 186738,    TDs_mg186738_swab186738_ari327_gp6199_analysis[4,3],   -1,                                                  -1 ,                                                -1,                                                                                                         TDs_mg186738_swab186738_ari327_gp6199_analysis[5,3],    TDs_mg186738_swab186738_ari327_gp6199_analysis[6,3],    TDs_mg186738_swab186738_ari327_gp6199_analysis[6,3] - TDs_mg186738_swab186738_ari327_gp6199_analysis[4,3]           ]

println(simple_combined_3TD_simrep_perc_df)


# Plot the DataFrame with points and lines
##### 1 case ####
# Summer
simple_combined_df_ex100winter = simple_combined_df[1:9,:]
x_replace_summer = simple_combined_df_ex100winter[:,1]#[100,200,319,747,1000,6000,15419]
# plot separately so coloured by y-series
p = plot(xlabel = "Number of primary care metagenomic samples"
        , ylabel = "Median time to detection of 1 case (TD) in days"
        , xscale=:log10
        , ylim=(0,70)
        , legend = :bottomleft
        #, size =(1600,1200)
        , xticks = [100, 1000, 10000, 100000, 1000000]
        )
@df simple_combined_df_ex100winter scatter!(x_replace_summer
                        , :ICU_only
                        , color=:red
                        , marker=:circle
                        , markersize=4
                        , label="ICU only")
@df simple_combined_df_ex100winter scatter!(x_replace_summer
                        , :PC_only_summer
                        , color=:blue
                        , marker=:circle
                        , markersize=4
                        , label="Primary care only (summer)")
@df simple_combined_df_ex100winter scatter!(x_replace_summer
                        , :Combined_summer
                        , color=:green
                        , marker=:circle
                        , markersize=4
                        , label="Combined (summer)")
@df simple_combined_df_ex100winter plot!(x_replace_summer, :ICU_only, color=:red, linewidth=2, label="")
@df simple_combined_df_ex100winter plot!(x_replace_summer, :PC_only_summer, color=:blue, linewidth=2, label="")
@df simple_combined_df_ex100winter plot!(x_replace_summer, :Combined_summer, color=:green, linewidth=2, label="")

# Winter
x_replace_winter = simple_combined_df[:,1] #[100,200,319,747,1000,6000,28011]
# plot separately so coloured by y-series
#p = plot(xlabel = "Number of primary care metagenomic samples"
#        , ylabel = "Time to detection of 1 case (days)"
#        , xscale=:log10
#        , ylim=(0,70)
#        , legend = :topright
#        )
@df simple_combined_df scatter!(x_replace_winter
                        , :ICU_only
                        , color=:red
                        , marker=:circle
                        , markersize=4
                        , label="")
@df simple_combined_df scatter!(x_replace_winter
                        , :PC_only_winter
                        , color=:lightblue
                        , marker=:circle
                        , markersize=4
                        , label="Primary care only (winter)")
@df simple_combined_df scatter!(x_replace_winter
                        , :Combined_winter
                        , color=:lightgreen
                        , marker=:circle
                        , markersize=4
                        , label="Combined (winter)")
@df simple_combined_df plot!(x_replace_winter, :ICU_only, color=:red, linewidth=2, label="")
@df simple_combined_df plot!(x_replace_winter, :PC_only_winter, color=:lightblue, linewidth=2, label="")
@df simple_combined_df plot!(x_replace_winter, :Combined_winter, color=:lightgreen, linewidth=2, label="")

# Save to file
savefig("scripts/primary_care_v_icu/TD_vs_PC_sample_size.png")

##### 3 cases ####
# Summer
simple_combined_3TD_df_ex100winter = simple_combined_3TD_df[1:9,:]
x_replace_summer = simple_combined_3TD_df_ex100winter[:,1]#[100,200,319,747,1000,6000,15419]
# plot separately so coloured by y-series
p = plot(xlabel = "Number of primary care metagenomic samples"
        , ylabel = "Median time to detection\n for 3 cases (3TD) in days"
        , xscale=:log10
        , ylim=(0,70)
        , legend = :bottomleft
        #, size =(1600,1200)
        , xticks = [100, 1000, 10000, 100000, 1000000]
        )
@df simple_combined_3TD_df_ex100winter scatter!(x_replace_summer
                        , :ICU_only
                        , color=:red
                        , marker=:circle
                        , markersize=4
                        , label="ICU only")
@df simple_combined_3TD_df_ex100winter scatter!(x_replace_summer
                        , :PC_only_summer
                        , color=:blue
                        , marker=:circle
                        , markersize=4
                        , label="Primary care only (summer)")
@df simple_combined_3TD_df_ex100winter scatter!(x_replace_summer
                        , :Combined_summer
                        , color=:green
                        , marker=:circle
                        , markersize=4
                        , label="Combined (summer)")
@df simple_combined_3TD_df_ex100winter plot!(x_replace_summer, :ICU_only, color=:red, linewidth=2, label="")
@df simple_combined_3TD_df_ex100winter plot!(x_replace_summer, :PC_only_summer, color=:blue, linewidth=2, label="")
@df simple_combined_3TD_df_ex100winter plot!(x_replace_summer, :Combined_summer, color=:green, linewidth=2, label="")

# Winter
x_replace_winter = simple_combined_3TD_df[:,1] #[100,200,319,747,1000,6000,28011]
# plot separately so coloured by y-series
#p = plot(xlabel = "Number of primary care metagenomic samples"
#        , ylabel = "Time to detection of 1 case (days)"
#        , xscale=:log10
#        , ylim=(0,70)
#        , legend = :topright
#        )
@df simple_combined_3TD_df scatter!(x_replace_winter
                        , :ICU_only
                        , color=:red
                        , marker=:circle
                        , markersize=4
                        , label="")
@df simple_combined_3TD_df scatter!(x_replace_winter
                        , :PC_only_winter
                        , color=:lightblue
                        , marker=:circle
                        , markersize=4
                        , label="Primary care only (winter)")
@df simple_combined_3TD_df scatter!(x_replace_winter
                        , :Combined_winter
                        , color=:lightgreen
                        , marker=:circle
                        , markersize=4
                        , label="Combined (winter)")
@df simple_combined_3TD_df plot!(x_replace_winter, :ICU_only, color=:red, linewidth=2, label="")
@df simple_combined_3TD_df plot!(x_replace_winter, :PC_only_winter, color=:lightblue, linewidth=2, label="")
@df simple_combined_3TD_df plot!(x_replace_winter, :Combined_winter, color=:lightgreen, linewidth=2, label="")

# Save to file
savefig("scripts/primary_care_v_icu/3TD_vs_PC_sample_size.png")




###### Analysis 7 Sep 2025 - updated parameters
# 1 # using 100 and 200 metagenomic samples per week
sim_tds_100_200_samples = icu_v_pc_td(; gp_swabs_mg = [100,200] #[319, 747] # Assumed number of swabs that are metagenomic sequenced for investigating impact
                                      , sim_object = "filtered" ) 
# Save as .csv file for inspection
CSV.write("scripts/primary_care_v_icu/sim_TDs_2025_09_07.csv", sim_tds_100_200_samples) # using 100 and 200 metagenomic samples per week
sim_tds_100_200_samples_analysis = analyse_columns(sim_tds_100_200_samples[:,2:23]) # Not essential but remove the column containing the simulation number
println(sim_tds_100_200_samples_analysis)

# 2 # using all samples for metagenomic sequencing each week
sim_tds_all_gp_samples = icu_v_pc_td(; gp_swabs_mg = [319, 747] #[100,200] # Assumed number of swabs that are metagenomic sequenced for investigating impact
                                    , sim_object = "filtered" )
# Save as .csv file for inspection
CSV.write("scripts/primary_care_v_icu/sim_TDs_all_gp_samples_2025_09_07.csv", sim_tds_all_gp_samples) # using 100 and 200 metagenomic samples per week
sim_tds_all_gp_samples_analysis = analyse_columns(sim_tds_all_gp_samples[:,2:23]) # Not essential but remove the column containing the simulation number and simid cols
println(sim_tds_all_gp_samples_analysis)

# 3 # Assume hit the targetted values for the number of GP swabs (1000 per week or 20 swabs per practice x 300 practices = 6000) and using all samples for metagenomic sequencing each week
# Same for all seasons - note that the differences between min and max are not meaningful in this case - it is just variance from random sampling
sim_tds_higher_gp_swab_target_1000 = icu_v_pc_td(; gp_swabs_mg = [1000, 1000] #[100,200] # Assumed number of swabs that are metagenomic sequenced for investigating impact
                                             , gp_ari_swabs = [1000, 1000] #[319, 747] # Number of swabs taken from suspected ARI per week [mean summer 2024, mean winter 2024/25]. Source: Analysis of data extracted from RCGP Research & Surveillance Centre (RSC) Virology Dashboard [Accessed 29 Aug 2025]
                                            , sim_object = "filtered"
                                            )
# Save as .csv file for inspection
CSV.write("scripts/primary_care_v_icu/sim_tds_higher_gp_swab_target_1000_2025_09_07.csv", sim_tds_higher_gp_swab_target_1000) # using 100 and 200 metagenomic samples per week
sim_tds_higher_gp_swab_target_1000_analysis = analyse_columns(sim_tds_higher_gp_swab_target_1000[:,2:23]) # Not essential but remove the column containing the simulation number
println(sim_tds_higher_gp_swab_target_1000_analysis)

# 4 # Assume hit the targetted values for the number of GP swabs (20 swabs per practice x 300 practices = 6000) and using all samples for metagenomic sequencing each week
# Same for all seasons - note that the differences between min and max are not meaningful in this case - it is just variance from random sampling
sim_tds_higher_gp_swab_target_6000 = icu_v_pc_td(; gp_swabs_mg = [1000, 6000] #[100,200] # Assumed number of swabs that are metagenomic sequenced for investigating impact
                                             , gp_ari_swabs = [1000, 6000] #[319, 747] # Number of swabs taken from suspected ARI per week [mean summer 2024, mean winter 2024/25]. Source: Analysis of data extracted from RCGP Research & Surveillance Centre (RSC) Virology Dashboard [Accessed 29 Aug 2025]
                                             , sim_object = "filtered"
                                             )
# Save as .csv file for inspection
CSV.write("scripts/primary_care_v_icu/sim_tds_higher_gp_swab_target_6000_2025_09_07.csv", sim_tds_higher_gp_swab_target_6000) # using 100 and 200 metagenomic samples per week
sim_tds_higher_gp_swab_target_6000_analysis = analyse_columns(sim_tds_higher_gp_swab_target_6000[:,2:23]) # Not essential but remove the column containing the simulation number
println(sim_tds_higher_gp_swab_target_6000_analysis)
Histogram(sim_tds_higher_gp_swab_target_6000[:,2],breaks=100)
using StatsPlots
plot_hist(df = sim_tds_higher_gp_swab_target_6000, col = 18)
plot_hist(df = sim_tds_higher_gp_swab_target_6000, col = 14)

# 4(a)
test_results = icu_v_pc_td_2(; gp_swabs_mg = 1000 #[100,200] # Assumed number of swabs that are metagenomic sequenced for investigating impact
                             , gp_ari_swabs = 1000 #[319, 747] # Number of swabs taken from suspected ARI per week [mean summer 2024, mean winter 2024/25]. Source: Analysis of data extracted from RCGP Research & Surveillance Centre (RSC) Virology Dashboard [Accessed 29 Aug 2025]
                             , sim_object = "filtered"
                             )


# 5 # Sampling 15% of all ARI consultations
#sim_tds_primary_care_p15 = icu_v_pc_td(; p_pc = 0.15 # Test to compare TDs with ICU with 0.15 sampling proportion # If p_pc not input then it is calculated from other inputs
#                                    , sim_object = "filtered" )
# Save as .csv file for inspection
#CSV.write("scripts/primary_care_v_icu/sim_TDs_primary_care_p15_2025_09_07.csv", sim_tds_primary_care_p15) # metagenomic sampling of 25% of all ARI consultations
#sim_tds_primary_care_p15_analysis = analyse_columns(sim_tds_primary_care_p15[:,2:23]) # Not essential but remove the column containing the simulation number
#println(sim_tds_primary_care_p15_analysis)

# 6 # Sampling 15% of all ARI consultations but using number of swabs instead of directly defining the proportion
# Same for all seasons - note that the differences between min and max are not meaningful in this case - it is just variance from random sampling
sim_tds_primary_care_15_summer = icu_v_pc_td(; gp_swabs_mg = [15419, 15419] #[100,200] # Assumed number of swabs that are metagenomic sequenced for investigating impact
                                             , gp_ari_swabs = [15419, 15419] #[319, 747] # Number of swabs taken from suspected ARI per week [mean summer 2024, mean winter 2024/25]. Source: Analysis of data extracted from RCGP Research & Surveillance Centre (RSC) Virology Dashboard [Accessed 29 Aug 2025]
                                            , sim_object = "filtered"
                                            )
# Save as .csv file for inspection
CSV.write("scripts/primary_care_v_icu/sim_tds_primary_care_15_summer_2025_09_07.csv", sim_tds_primary_care_15_summer) # using 100 and 200 metagenomic samples per week
sim_tds_primary_care_15_summer_analysis = analyse_columns(sim_tds_primary_care_15_summer[:,2:23]) # Not essential but remove the column containing the simulation number
println(sim_tds_primary_care_15_summer_analysis)

# 7 # Sampling 15% of all ARI consultations but using number of swabs instead of directly defining the proportion
# Same for all seasons - note that the differences between min and max are not meaningful in this case - it is just variance from random sampling
sim_tds_primary_care_15_winter = icu_v_pc_td(; gp_swabs_mg = [28011, 28011] #[100,200] # Assumed number of swabs that are metagenomic sequenced for investigating impact
                                             , gp_ari_swabs = [28011, 28011] #[319, 747] # Number of swabs taken from suspected ARI per week [mean summer 2024, mean winter 2024/25]. Source: Analysis of data extracted from RCGP Research & Surveillance Centre (RSC) Virology Dashboard [Accessed 29 Aug 2025]
                                            , sim_object = "filtered"
                                            )
# Save as .csv file for inspection
CSV.write("scripts/primary_care_v_icu/sim_tds_primary_care_p15_winter_2025_09_07.csv", sim_tds_primary_care_15_winter) # using 100 and 200 metagenomic samples per week
sim_tds_primary_care_15_winter_analysis = analyse_columns(sim_tds_primary_care_15_winter[:,2:23]) # Not essential but remove the column containing the simulation number
println(sim_tds_primary_care_15_winter_analysis)


##### Compilation of results

# Create compilation of data
# Times to detection
combined_df = DataFrame([zeros(Float64,9) for _ in 1:6], [:Number_of_PC_mg_samples
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
#combined_df[7,:] = [ 0.15,                    sim_tds_primary_care_p15_analysis[1,2],           sim_tds_primary_care_p15_analysis[3,2],             sim_tds_primary_care_p15_analysis[11,2],            sim_tds_primary_care_p15_analysis[7,2],             sim_tds_primary_care_p15_analysis[15,2]          ]
#combined_df[7,5] = combined_df[7,3]
#combined_df[7,6] = combined_df[7,4]
combined_df[8,:] = [ 15419,                   sim_tds_primary_care_15_summer_analysis[1,2],     sim_tds_primary_care_15_summer_analysis[3,2],       sim_tds_primary_care_15_summer_analysis[11,2],      sim_tds_primary_care_15_summer_analysis[7,2],       sim_tds_primary_care_15_summer_analysis[15,2]    ]
combined_df[9,:] = [ 28011,                   sim_tds_primary_care_15_winter_analysis[1,2],     sim_tds_primary_care_15_winter_analysis[3,2],       sim_tds_primary_care_15_winter_analysis[11,2],      sim_tds_primary_care_15_winter_analysis[7,2],       sim_tds_primary_care_15_winter_analysis[15,2]    ]


println(combined_df)
# Difference between combined and ICU_only for summer
println( combined_df[:,4]-combined_df[:,2] )
# Difference between combined and ICU_only for winter
println( combined_df[:,6]-combined_df[:,2] )


# % of sim reps TODO
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
combined_sim_rep_df[7,:] = [ 0.15,                    sim_tds_primary_care_p15_analysis[1,3],           sim_tds_primary_care_p15_analysis[3,3],             sim_tds_primary_care_p15_analysis[11,3],            sim_tds_primary_care_p15_analysis[7,3],             sim_tds_primary_care_p15_analysis[15,3]          ]
combined_sim_rep_df[7,5] = combined_sim_rep_df[7,3]
combined_sim_rep_df[7,6] = combined_sim_rep_df[7,4]
println(combined_sim_rep_df)

##### TODO PLOTS

###### Analysis

# 1 # using 100 and 200 metagenomic samples per week
sim_tds_100_200_samples = icu_v_pc_td(; gp_swabs_mg = [100,200] #[319, 747] # Assumed number of swabs that are metagenomic sequenced for investigating impact
                                    
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

# 3.0 # Sampling 15% of all ARI consultations
sim_tds_primary_care_p15 = icu_v_pc_td(; p_pc = 0.15 ) # Test to compare TDs with ICU with 0.15 sampling proportion # If p_pc not input then it is calculated from other inputs
# Save as .csv file for inspection
CSV.write("scripts/primary_care_v_icu/sim_TDs_primary_care_p15_v2.csv", sim_tds_primary_care_p15) # metagenomic sampling of 25% of all ARI consultations
sim_tds_primary_care_p15_analysis = analyse_columns(sim_tds_primary_care_p15[:,2:19]) # Not essential but remove the column containing the simulation number
println(sim_tds_primary_care_p15_analysis)

# 3.0 # Sampling 10% of all ARI consultations
sim_tds_primary_care_p10 = icu_v_pc_td(; p_pc = 0.10 ) # Test to compare TDs with ICU with 0.15 sampling proportion # If p_pc not input then it is calculated from other inputs
# Save as .csv file for inspection
CSV.write("scripts/primary_care_v_icu/sim_TDs_primary_care_p10_v2.csv", sim_tds_primary_care_p10) # metagenomic sampling of 25% of all ARI consultations
sim_tds_primary_care_p10_analysis = analyse_columns(sim_tds_primary_care_p10[:,2:19]) # Not essential but remove the column containing the simulation number
println(sim_tds_primary_care_p10_analysis)

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
combined_df[7,:] = [ 0.15,                    sim_tds_primary_care_p15_analysis[1,2],           sim_tds_primary_care_p15_analysis[3,2],             sim_tds_primary_care_p15_analysis[11,2],            sim_tds_primary_care_p15_analysis[7,2],             sim_tds_primary_care_p15_analysis[15,2]          ]
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
combined_sim_rep_df[7,:] = [ 0.15,                    sim_tds_primary_care_p15_analysis[1,3],           sim_tds_primary_care_p15_analysis[3,3],             sim_tds_primary_care_p15_analysis[11,3],            sim_tds_primary_care_p15_analysis[7,3],             sim_tds_primary_care_p15_analysis[15,3]          ]
combined_sim_rep_df[7,5] = combined_sim_rep_df[7,3]
combined_sim_rep_df[7,6] = combined_sim_rep_df[7,4]
println(combined_sim_rep_df)



using DataFrames
using Pkg
Pkg.add("StatsPlots")
using StatsPlots  # This imports the @df macro for easier DataFrame plotting

# Plot the DataFrame with points and lines
# Summer
x_replace_summer = [100,200,319,747,1000,6000,15419]
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
x_replace_winter = [100,200,319,747,1000,6000,28011]
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

# Analysis with lower p_icu for ICU sampling
# 1 # using 100 and 200 metagenomic samples per week
sim_tds_100_200_samples_p_icu_10 = icu_v_pc_td(; p_icu = 0.10 # 0.25
                                                , gp_swabs_mg = [100,200] ) #[319, 747] # Assumed number of swabs that are metagenomic sequenced for investigating impact
# Save as .csv file for inspection
CSV.write("scripts/primary_care_v_icu/sim_TDs_p_icu_10_v2.csv", sim_tds_100_200_samples_p_icu_10) # using 100 and 200 metagenomic samples per week
sim_tds_100_200_samples_p_icu_10_analysis = analyse_columns(sim_tds_100_200_samples_p_icu_10[:,2:19]) # Not essential but remove the column containing the simulation number
println(sim_tds_100_200_samples_p_icu_10_analysis)


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