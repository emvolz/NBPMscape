# Load packages
using Revise
using NBPMscape
using GLM, Statistics, StatsBase, Distributions
using DataFrames, CSV 
using Plots, StatsPlots
using JLD2

#=  This file contains a number of functions used in the analysis of time to detection 
    under different scenarios. The functions are:
    - icu_v_pc_td - computes the times to detection for simulated outbreaks
    - analyse_td_columns - computes median values across simulation replicates and % of sim reps that have a (time to) detection
    - clean_td_results TODO NOT YET WORKING - removes some rows (simulation replicates) and converts type of Inf from string
    - plot_hist - plots a basic histogram for results of a chosen column of data after performing some data cleaning tasks 
                  (similar to those in {clean_TD_results})
=#


"""
Function: icu_v_pc_td
    Function to compute times to detection for 1 case (TD) an 3 cases (3TD)
    by sampling from:
    - primary care
    - ICU

# Arguments
    'sim_object_type':          Dscription of the format of the underlying simulation data. 
                                There are two options: 
                                -   "full", which is the data as output from {simtree} or {simforest}. 
                                    This option requires an input for "sims_file" and no input for 
                                    "sims_G_gp_filter_file" and "sims_G_icu_filter_file"
                                -   "filtered", which is the "full" data already filtered for GP and ICU
                                     cases and split into two files. This option requires inputs for 
                                     "sims_G_gp_filter_file" and "sims_G_icu_filter_file" and no input for "sims_file"
    'file_or_object'            A binary choice between "file" and "object" that will determine whether or not the function
                                needs to load a file. Loading files can be time consuming so if running the function
                                multiple times, i.e. with different parameters, it is beneficial to load the file once
                                and then input the object from the file in subsequent runs.
    'sims_file':                .jld2 file containing simulation replicates as output by {simtree} or {simforest}
    'sims_G_gp_filter_file':    .jld2 file containing simulation replicates (as output by {simtree} or {simforest}) BUT
                                filtered for GP cases in the G dataframe (infection information) using {sims_filter}
    'sims_G_icu_filter_file':   .jld2 file containing simulation replicates (as output by {simtree} or {simforest}) BUT
                                filtered for ICU cases in the G dataframe (infection information) using {sims_filter}
    'sims_object_name':         Name of object containing simulation replicates as output by {simtree} or {simforest},
                                or as loaded from .jld2 file
    'sims_G_gp_filter_object_name': Name of object containing simulation replicates (as output by {simtree} or {simforest}) BUT
                                    filtered for GP cases in the G dataframe (infection information) using {sims_filter},
                                    or as loaded from .jld2 file
    'sims_G_icu_filter_object_name':Name of object containing simulation replicates (as output by {simtree} or {simforest}) BUT
                                    filtered for ICU cases in the G dataframe (infection information) using {sims_filter},
                                    or as loaded from .jld2 file
    'p_icu':                    ICU sampling proportion
    'icu_ari_admissions':       Weekly ICU admission numbers. Formatted in a 2-element vector: [summer,winter]
    'icu_turnaround_time':      Time to process sample and report results / declare detection. 
                                Formatted in a 2-element vector: [min,max]
    'gp_practices_total':       Total number of GP practices
    'gp_practices_swab':        Number of GP practices taking swabs for virology surveillance
    'gp_swabs_mg':              Assumed number of swabs that are metagenomic sequenced for investigating impact
    'pop_eng':                  Population of England
    'gp_ari_consults':          Number of ARI consultations per 100k of England population per week. 
                                Formatted in a 2-element vector: [mean summer 2024, mean winter 2024/25]
    'gp_ari_swabs':             Number of swabs taken from suspected ARI per week. 
                                Formatted in a 2-element vector: [mean summer 2024, mean winter 2024/25]
    'pc_swab_turnaround_time':  Number of days between swab sample being taken and results received. 
                                Formatted in a 2-element vector: [min,max] 

# Returns
    'sim_tds'                   A dataframe containing the times to detection of 1 case and 3 cases for each simulation replicate
                                assuming either:
                                - only sampling ICU admissions
                                - only sampling primary care (eg GP) cases
                                - sampling from both settings
                                The output is of this format:
                                sim_n	ICU_TD	    PC_TD	ICU_PC_TD	ICU_3TD	    PC_3TD	ICU_PC_3TD	n_ICU_cases	n_ICU_cases_sampled	n_GP_cases	n_GP_cases_sampled	ICU_simid	                                GP_simid
                                1	    53.83937265	Inf	    53.83937265	60.83018273	Inf	    60.83018273	45	        5	                621	        0	                ["f6f6b432-8e44-11f0-3306-2724694310b1"]	["f6f6b432-8e44-11f0-3306-2724694310b1", "02e49496-8e45-11f0-256a-598e6634981b"]
                                2	    Inf	        Inf	    Inf	        Inf	        Inf	    Inf	        0	        0	                0	        0	                Inf	                                        Inf
                                ...     ...         ...     ...         ...         ...     ...         ...         ...                 ...         ...                 ...                                         ...

# Examples
    td_results_test_1 = icu_v_pc_td(; sim_object_type = "full", file_or_object = "file"
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

function icu_v_pc_td(;    sim_object_type # = "filtered" or "full"
                        , file_or_object # "file" or "object"
                        , sims_file = "" # "sims.jdl2"
                        , sims_G_gp_filter_file = "" # e.g. "sims_G_gp_filter.jdl2"
                        , sims_G_icu_filter_file = "" # e.g. "sims_G_icu_filter.jdl2"
                        , sims_object_name = "" # "sims" 
                        , sims_G_gp_filter_object_name = "" # e.g. "sims_G_gp_filter" 
                        , sims_G_icu_filter_object_name = "" # e.g. "sims_G_icu_filter" 
                        # ICU parameters
                        , p_icu = 0.15 # ICU sampling proportion TODO Assumption needs refining 
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
        if sim_object_type == "full"
            sims = load(sims_file, sims_object_name)
        elseif sim_object_type == "filtered"
            sims_G_icu_filter = load(sims_G_icu_filter_file, sims_G_icu_filter_object_name)
            sims_G_gp_filter  = load(sims_G_gp_filter_file,  sims_G_gp_filter_object_name)
        end
    elseif file_or_object == "object" # ... or just assign objects to variable names if already loaded
        if sim_object_type == "full"
            isdefined(Main, :sims_object_name) ? sims = getfield(Main, Symbol(sims_object_name)) : nothing
        elseif sim_object_type == "filtered"
            isdefined(Main, :sims_G_icu_filter_object_name) ? sims_G_icu_filter = getfield(Main, Symbol(sims_G_icu_filter_object_name)) : nothing
            isdefined(Main, :sims_G_gp_filter_object_name)  ? sims_G_gp_filter  = getfield(Main, Symbol(sims_G_gp_filter_object_name))  : nothing
        end
    end

    # How many replicates are in this simulation
    if sim_object_type == "full"
        n_replicates = length(sims)
    elseif sim_object_type == "filtered"
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
            # Filter for ICU cases 
            icu_cases = fo.G[ isfinite.(fo.G.ticu), : ]
            # Filter for GP cases (some of these will also become ICU cases)
            gp_cases = fo.G[ isfinite.(fo.G.tgp), : ]    
        elseif sim_object_type == "filtered"
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



"""
Function: analyse_td_columns
    Function to compute:
        - median values for each time to detection scenario across multiple simulation replicates
        - % of simulations that return a time to detection in each time to detection scenario

# Arguments
    'df':   A dataframe as output from the {icu_v_pc_td} function containing the times to detection of
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
            
            Ideally this should be trimmed to only include the columns related to time to detected, i.e. columns 2 to 11
            
              Row │  ICU_TD   PC_TD    ICU_PC_TD  ICU_3TD  PC_3TD  ICU_PC_3TD  n_ICU_cases  n_ICU_cases_sampled  n_GP_cases  n_GP_cases_sampled 
                  │    Any      Any      Any        Any      Any     Any         Any          Any                  Any         Any             
            ──────┼────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
                1 │      Inf      Inf      Inf        Inf      Inf     Inf         0            0                    1           0                   
                2 │      Inf      Inf      Inf        Inf      Inf     Inf         0            0                    8           0                   
              ⋮   │       ⋮        ⋮         ⋮         ⋮       ⋮         ⋮            ⋮                ⋮               ⋮               ⋮          
              999 │    49.082   69.3822  49.082     55.1244  Inf     55.1244     104          17                   1096        1                   
             1000 │   Inf      Inf      Inf        Inf      Inf     Inf         0            0                    0           0                   
# Returns
    'result':   A dataframe containing three columns: description of parameter (mostly relating to time to detection), 
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
# Examples
    td_results_test_1_analysis = analyse_td_columns(td_results_test_1[:, 2:11])
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


# TODO NOT YET COMPLETED
"""
Function: clean_td_results
    Function to remove rows (simulation replicates) containing "No Sim ID match" and 
    convert string values of Inf to numeric values.

    # Arguments
    'df':   A dataframe as output from the {icu_v_pc_td} function containing the times to detection of
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
            
# Examples
    td_results_test_1_clean = clean_TD_results(td_results_test_1)#[:, 2:11])
    println( td_results_test_1_analysis )
"""
# TODO NOT YET COMPLETED

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
    Function to plot a histogram of a particular column a dataframe.

    # Arguments
    'df':   A dataframe as output from the {icu_v_pc_td} function containing the times to detection of
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

# Returns
    Plots histogram to screen

# Examples
    plot_hist( df = td_results_test_1, col = 2)
    
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