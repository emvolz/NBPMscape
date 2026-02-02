#=
Functions used in analysing the timeline of the simulated outbreak

- infections_time_analysis
        Analysis of simulation output files to create time series for outbreak, including 
        the number of infections, import, numbers at each healthcare stage. Also computes the
        counts of infections, imports and healthcare setting visits/admissions at times of
        detection. In addition, estimates the doubling time of the outbreak.

- plot_outbreak_analysis
        Plots counts and cumulative values for 8 metrics relating to infections and healthcare seeking

- plot_outbreak_analysis_combined
        Plots counts and cumulative values for 8 metrics relating to infections and healthcare seeking
        across 5 plots on a 3x2 grid
        TODO MAY CURRENTLY ONLY WORK WITH time_period = "days" and x_labels = "date" options

- median_td_sim_reps                
        Function to return the median time to detection from a file containing 
        times to detection for many simulation replicates

- cases_before_td
        Counts the number of infected cases before the first case is detected.

- local_doubling_time
        Returns a vector of local doubling times assuming exponential growth between observations.

- global_doubling_time  
        Computes the doubling time of a time series (cumulative infections) 
        in the context of the simulated outbreak.

sub-functions used in the global_doubling_time function

- _fit_log_linear                   Returns doubling time using log linear fitting method
- _fit_weighted_log_linear          Returns doubling time using weighted log linear fitting method
- _growth_rate_to_doubling_time     Coverts growth rate to doubling time
- _fit_nonlinear_optim              Returns doubling time using least squares fit method via Optim package
- _estimate_initial_r               Estimates growth rate, r, assuming exponential growth

=#


"""
Function        median_td_sim_reps

Description     Function to return the median time to detection from a file containing 
                times to detection for many simulation replicates. Also identifies the
                simulation replicate number and the file number (large simulations can be 
                run on an HPC in arrays and the multiple files containing smaller numbers
                of simulation replicates can be combined for analysis)

Arguments       tds_file::String        Path and filename for file containing the times to detection data
                td_column_name::String  Name of column containing the time to detection information
                nreps                   Number of simulation replicates per simulation file (assuming that the TDs in the file input here are compiled from multiple simulation files)

Returns         median_td           Median time to detection from all simulation replicates
                median_td_lower     TD value below median when total number of sim reps is even
                median_td_upper     TD value above median when total number of sim reps is even
                median_td_lower_sim_rep_n   Sim rep number relating to median_td_lower
                median_td_upper_sim_rep_n   Sim rep number relating to median_td_upper
                median_td_sim_rep_n         Sim rep number relating to median_td when total number of sim reps is odd
                file_n_med_lower                File number containing the sim rep with the lower median value 
                                                (Note that if the combine_sim_reps() function was used to combine
                                                sim reps from multiple files then the order may not be numerical, 
                                                i.e. file.11.jld2 comes before file.2.jld2)
                sim_rep_n_in_file_n_med_lower   The number of the sim rep within the smaller file prior to combination,
                                                for the sim rep producing the TD immediatley below the TD 
                                                (when the number of sim reps is even and the median TD is therefore the midpoint)
                file_n_med_upper                File number containing the sim rep with the upper median value
                sim_rep_n_in_file_n_med_upper   The number of the sim rep within the smaller file prior to combination,
                                                for the sim rep producing the TD immediatley above the TD 
                                                (when the number of sim reps is even and the median TD is therefore the midpoint)
                file_n_med                      File number containing the sim rep with the median value 
                                                (when the number of sim reps is odd)
                sim_rep_n_in_file_n_med         The number of the sim rep within the smaller file prior to combination,
                                                for the sim rep producing the median TD (when the number of sim reps is
                                                odd and the TD is therefore exact)

Examples
                median_td_sim_reps(; tds_file = "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1287570_nrep2000_analysis/sc_tds_1.csv"
                                   , td_column_name = "SC_TD"
                                   , nreps = 50 )
                # Returns
                (median_td = 51.384143237678416, median_td_lower = 51.38412541748787, median_td_upper = 51.384161057868965
                , median_td_lower_sim_rep_n = 1004, median_td_upper_sim_rep_n = 52, median_td_sim_rep_n = missing)

"""
function median_td_sim_reps(; tds_file::String = "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1287570_nrep2000_analysis/sc_tds_1.csv"
                            , td_column_name::String = "SC_TD"
                            , nreps = 50
                            )
    # Read in file containing times to detection for multiple simulations
    tds_df = CSV.read( tds_file, DataFrame)     #file = joinpath( td_file, "sc_tds_$(string(scenario_number)).csv")

    # Number of simulation replicates
    sim_reps = nrow(tds_df)

    # Extract median time to detection (TD)
    median_td = median( tds_df[:, Symbol(td_column_name) ] )

    # Identify TDs above and below the median (if the number of sim reps is even)
    if iseven( sim_reps )

        # Times to detection either side of median value
        median_td_lower = sort( tds_df, Symbol(td_column_name) )[ Int(sim_reps / 2)     , Symbol(td_column_name) ]
        median_td_upper = sort( tds_df, Symbol(td_column_name) )[ Int(sim_reps / 2) + 1 , Symbol(td_column_name) ]
        # Simulation replicate numbers for TDs closest to median value
        median_td_lower_sim_rep_n = only( findall(==(median_td_lower), tds_df[:, Symbol(td_column_name) ] ) )
        median_td_upper_sim_rep_n = only( findall(==(median_td_upper), tds_df[:, Symbol(td_column_name) ] ) )
        # File numbers and sim rep within file that relate to median TD value(s)
        file_n_med_lower = median_td_lower_sim_rep_n % nreps == 0 ? Int( median_td_lower_sim_rep_n / nreps ) : Int( floor( median_td_lower_sim_rep_n / nreps ) ) + 1 
        sim_rep_n_in_file_n_med_lower = Int( median_td_lower_sim_rep_n - ((file_n_med_lower-1) * nreps) )
        file_n_med_upper = median_td_upper_sim_rep_n % nreps == 0 ? Int( median_td_upper_sim_rep_n / nreps ) : Int( floor( median_td_upper_sim_rep_n / nreps ) ) + 1 
        sim_rep_n_in_file_n_med_upper = Int( median_td_upper_sim_rep_n - ((file_n_med_upper-1) * nreps) )
        
        # If sim_reps is even then this has no value
        median_td_sim_rep_n = missing
        # File numbers and sim rep within file that relate to median TD value
        file_n_med = missing
        sim_rep_n_in_file_n_med = missing
        
    elseif isodd( sim_reps )
        
        # Median time to detection
        median_td_sim_rep_n = only( findall(==(median_td), tds_df[:, Symbol(td_column_name) ] ) )
        # File numbers and sim rep within file that relate to median TD value
        file_n_med = median_td_sim_rep_n % nreps == 0 ? Int( median_td_sim_rep_n / nreps ) : Int( floor( median_td_sim_rep_n / nreps ) ) + 1 
        sim_rep_n_in_file_n_med = Int( median_td_sim_rep_n - ((file_n_med-1) * nreps) )
        
        # If sim_reps is odd then these have no value
        median_td_lower_sim_rep_n = missing
        median_td_upper_sim_rep_n = missing
        file_n_med_lower = missing
        sim_rep_n_in_file_n_med_lower = missing
        file_n_med_upper = missing
        sim_rep_n_in_file_n_med_upper = missing
        
    end

    return( ( median_td = median_td
            , median_td_lower = median_td_lower
            , median_td_upper = median_td_upper
            , median_td_lower_sim_rep_n = median_td_lower_sim_rep_n
            , median_td_upper_sim_rep_n = median_td_upper_sim_rep_n
            , median_td_sim_rep_n = median_td_sim_rep_n
            , file_n_med_lower = file_n_med_lower
            , sim_rep_n_in_file_n_med_lower = sim_rep_n_in_file_n_med_lower
            , file_n_med_upper = file_n_med_upper
            , sim_rep_n_in_file_n_med_upper = sim_rep_n_in_file_n_med_upper
            , file_n_med = file_n_med
            , sim_rep_n_in_file_n_med = sim_rep_n_in_file_n_med             
            ) )
end


"""
Function        cases_before_td

Description     Counts the number of infected cases before the first case is detected. 
                Requires the time to detection as an input.

Arguments       sims_simid_tinf_df_file     Name of .jld2 file containing vecotr of dataframes, each containing
                                            simids and times of infection, as output from simulation.
                sim_rep_n                   Number of the particular simulation (element within vector) within file
                                            that the TD value relates to.
                td_value                    The number of cases of infection before or equal to this time will be counted

Returns         Number of cases of infection before the td_value

Examples        
            cases_before_td(; sims_simid_tinf_df_file = "covidlike-1.4.1-sims-nrep50_sims_simid_tinf_df_1287570.29.jld2" # 21st file in folder ordered by name
                            , sim_rep_n = sc_original_timport_med_td.sim_rep_n_in_file_n_med_lower # 4
                            , td_value = sc_original_timport_med_td.median_td_lower #51.38412541748787
                            )
            # 1151

"""
function cases_before_td(; sims_simid_tinf_df_file = "covidlike-1.4.1-sims-nrep50_sims_simid_tinf_df_1287570.10.jld2"
            , sim_rep_n
            , td_value
             )
    # Load sims_simid_tinf_df for the individual sim reps relating to the median TD
    sims_simid_tinf_df = load( sims_simid_tinf_df_file, "sims_simid_tinf_df")
    
    # simid and tinf for simulation replicates that correspond to median TD
    simid_tinf_df = sims_simid_tinf_df[ sim_rep_n ]
    #println( sort(simid_tinf_df, :tinf)[1:51,:] )

    # Count number of cases before detection
    n_cases_before_td = nrow( filter(row -> row.tinf <= td_value, sort(simid_tinf_df, :tinf) ) )
    
    return( n_cases_before_td )
end


"""
Function        infections_time_analysis

Description     Analysis of simulation output files to create time series for outbreak, including 
                the number of infections, import, numbers at each healthcare stage. 
                Also computes the counts of infections, imports and healthcare setting visits/admissions
                at times of detection. In addition, estimates the doubling time of the outbreak.

Arguments       sims_simid_tinf_df_object   Name of object containing the sims_simid_tinf_df data preloaded from a single simulation output file (see examples below)
                sims_G_filtered_object      Name of object containing the sims_G_filtered data preloaded from a single simulation output file (see examples below)
                sims_max_cases_df_object    Name of object containing the sims_max_cases_df data preloaded from a single simulation output file (see examples below)
                td_gp                       Time to detect 1st case using GP sampling
                td_sc                       Time to detect 1st case using Secondary Care sampling, e.g. via HARISS network
                td_icu                      Time to detect 1st case using ICU sampling
                td3_gp                      Time to detect 3rd case using GP sampling
                td3_sc                      Time to detect 3rd case using Secondary Care sampling, e.g. via HARISS network
                td3_icu                     Time to detect 3rd case using ICU sampling
                base_date                   Example date of first import into UK, e.g. Date(2020, 1, 15)
                maxtime                     Parameter used when simtree or simforest run used to create the input data, e.g. 90 #days
                
                Note: td values used here should be for the same sim rep otherwise results (of number of m at TD) 
                are not really comparable

Returns         Returns a tuple containing:
                - Values for tinf, tgp, ted, thosp, ticu, tdeceased, tinf for fatal infections, timport
                - A df containing counts and cumulative values for these times by day and date (given input base date)
                - A df containing counts and cumulative values for these times by week (given input base date)
                - Counts at time of 1st detection (TD) for infections, timports and healthcare settings
                - Counts at time of 3rd detection (3TD) for infections, timports and healthcare settings
                - Estimated infection doubling times (local time series and df of global single values)
              
Examples    ## Example 1    
            ## Load files
            # load sims_simid_tinf_df
            sims_simid_tinf_df_1287570_29 = load( "covidlike-1.4.1-sims-nrep50_sims_simid_tinf_df_1287570.29.jld2", "sims_simid_tinf_df") # 21st file in folder
            sims_simid_tinf_df_1287570_29_4_1004 = sims_simid_tinf_df_1287570_29[4]
            # Format
            257052×2 DataFrame
                Row │ simid                              tinf     
                    │ String                             Float64
            ────────┼─────────────────────────────────────────────
                1 │ c2464d9c-de46-11f0-3e1c-39eb3c3e…   0.0
                2 │ c2464d9c-de46-11f0-3e1c-39eb3c3e…   9.12721
                3 │ c2464d9c-de46-11f0-3e1c-39eb3c3e…   9.47821
                4 │ c2464d9c-de46-11f0-3e1c-39eb3c3e…  11.8267
                5 │ c2464d9c-de46-11f0-3e1c-39eb3c3e…  13.497
            ⋮    │                 ⋮                     ⋮
            257049 │ 441c2d96-de47-11f0-22fd-1933a5e6…  45.7849
            257050 │ 441c2d96-de47-11f0-22fd-1933a5e6…  51.7974
            257051 │ 441c8816-de47-11f0-0046-23333003…  48.2385
            257052 │ 441cb0d6-de47-11f0-3ef4-a3f7f86f…  48.348
            
            # load filtered G df
            #sims_1287570 = load("covidlike-1.4.1-sims-filtered_G_nrep2000_1287570.jld2", "sims") # sim used in HARISS analysis Jan 2026
            #sims_1287570_29_4_1004 = sims_1287570[1004] # HPC run 29, sim rep 4 (in HPC run 29), sim rep 1004 in combined .jld2 file, 21st file in folder used to combine
            # or
            sims_1287570_29 = load( "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1287570/1287570_sims_G_filtered_40/covidlike-1.4.1-sims-nrep50_filtered_G_1287570.29.jld2", "sims_G_filtered" )
            sims_1287570_29_4_1004 = sims_1287570_29[4] # HPC run 29, sim rep 4 (in HPC run 29), sim rep 1004 in combined .jld2 file, 21st file in folder used to combine
            # Format
            36857×18 DataFrame
            Row │ pid                                tinf      tgp       ted       thospital  ticu      tstepdown  tdischarge  trecovered  tdeceased  severity                fatal  iscommuter  homeregion  infector_age  infectee_age  importedinfection  simid                             
                │ String                             Float64   Float64   Float64   Float64    Float64   Float64    Float64     Float64     Float64    Symbol                  Bool   Bool        String      Int8?         Int8          Bool               String
            ───────┼─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
                1 │ c2464d9c-de46-11f0-3e1c-39eb3c3e…   9.12721  Inf        22.2307   Inf       Inf             Inf     22.6371     33.8434   Inf       moderate_ED             false       false  TLI7                  35             9              false  c2464d9c-de46-11f0-3e1c-39eb3c3e…
                2 │ c2464d9c-de46-11f0-3e1c-39eb3c3e…  19.692     28.5249  Inf         28.5495  Inf             Inf     28.8786     52.4668   Inf       severe_hosp_short_stay  false       false  TLG3                  11            38              false  c2464d9c-de46-11f0-3e1c-39eb3c3e…
                3 │ c2464d9c-de46-11f0-3e1c-39eb3c3e…  23.2961    28.7786  Inf         29.4848  Inf             Inf     29.8335     47.6694   Inf       severe_hosp_short_stay  false       false  TLG3                  74            75              false  c2464d9c-de46-11f0-3e1c-39eb3c3e…
                4 │ c2464d9c-de46-11f0-3e1c-39eb3c3e…  23.8932    27.7191  Inf        Inf       Inf             Inf    Inf          46.3213   Inf       moderate_GP             false        true  TLG3                  57            31              false  c2464d9c-de46-11f0-3e1c-39eb3c3e…
                5 │ c2464d9c-de46-11f0-3e1c-39eb3c3e…  28.5098    32.3764  Inf         32.5412  Inf             Inf    Inf          36.8546   Inf       severe_hosp_long_stay   false       false  TLG3                  64            79              false  c2464d9c-de46-11f0-3e1c-39eb3c3e…
            ⋮   │                 ⋮                     ⋮         ⋮         ⋮          ⋮         ⋮          ⋮          ⋮           ⋮           ⋮                ⋮               ⋮        ⋮           ⋮            ⋮             ⋮                ⋮                          ⋮
            36854 │ 43ce16e2-de47-11f0-2f43-cfccec59…  88.0894   100.312   Inf        Inf       Inf             Inf    Inf         131.415    Inf       moderate_GP             false        true  TLG1                  58            55              false  43ce16e2-de47-11f0-2f43-cfccec59…
            36855 │ 43ce16e2-de47-11f0-2f43-cfccec59…  86.6421    99.9856  Inf        Inf       Inf             Inf    Inf         112.335    Inf       moderate_GP             false       false  TLG1                  50            50              false  43ce16e2-de47-11f0-2f43-cfccec59…
            36856 │ 43ce16e2-de47-11f0-2f43-cfccec59…  87.7424    92.5081   93.8507   Inf       Inf             Inf     94.0891    112.133    Inf       moderate_ED             false        true  TLG1                  50            31              false  43ce16e2-de47-11f0-2f43-cfccec59…
            36857 │ 43ce16e2-de47-11f0-2f43-cfccec59…  88.0336    90.2208  Inf         90.7357  Inf             Inf    106.769     106.866    Inf       severe_hosp_long_stay   false       false  TLG1                  79            98              false  43ce16e2-de47-11f0-2f43-cfccec59…

            # Load max_cases_df
            sims_max_cases_df_29 = load( "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1287570/1287570_max_cases_df/covidlike-1.4.1-sims-nrep50_max_cases_df_1287570.29.jld2", "sims_max_cases_df" )
            sims_max_cases_df_29_4_1004 = sims_max_cases_df_29[4] # HPC run 29, sim rep 4 (in HPC run 29), sim rep 1004 in combined .jld2 file, 21st file in folder used to combine
            
            # Run function
            using NBPMscape
            result_1 = infections_time_analysis(; sims_simid_tinf_df_object = sims_simid_tinf_df_1287570_29_4_1004
                                    , sims_G_filtered_object = sims_1287570_29_4_1004
                                    , sims_max_cases_df_object = sims_max_cases_df_29_4_1004
                                    #, sim_rep_n
                                    , td_gp = Inf
                                    , td_sc = 51.38412541748787 # sc_original_timport_med_td.median_td_lower
                                    , td_icu = Inf
                                    , td3_gp = Inf
                                    , td3_sc = 51.38412541748787 # sc_original_timport_med_td.median_td_lower
                                    , td3_icu = 62.4 #Inf
                                    , base_date = Date(2020, 1, 15)
                                    , maxtime = 90 #days
                                    )
            plot(result_1.tinf_local_doubling_times_vec, xlabel = "Time since first importation (days)", ylabel = "Estimated doubling time (days)", label = "Local estimate")
            hline!( repeat( [ result_1.tinf_global_doubling_times_df[3,:doubling_time_estimates] ], length(result_1.tinf_local_doubling_times_vec) ), label = "Global estimate using nonlinear least squares")
            
            ### Example 2
            ## Load files
            # load sims_simid_tinf_df
            sims_simid_tinf_df_1287570_10 = load( "covidlike-1.4.1-sims-nrep50_sims_simid_tinf_df_1287570.10.jld2", "sims_simid_tinf_df")
            sims_simid_tinf_df_1287570_10_2_52 = sims_simid_tinf_df_1287570_10[2]
            
            # load filtered G df
            sims_1287570_10 = load( "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1287570/1287570_sims_G_filtered_40/covidlike-1.4.1-sims-nrep50_filtered_G_1287570.10.jld2", "sims_G_filtered" )
            sims_1287570_10_2_52 = sims_1287570_10[2] # HPC run 10, sim rep 2 (in HPC run 10), sim rep 52 in combined .jld2 file, 2nd file in folder used to combine
            
            # Run function
            using NBPMscape
            result_2 = infections_time_analysis(; sims_simid_tinf_df_object = sims_simid_tinf_df_1287570_10_2_52
                                    , sims_G_filtered_object = sims_1287570_10_2_52
                                    , sims_max_cases_df_object
                                    #, sim_rep_n
                                    , td_gp = Inf
                                    , td_sc = 51.384161057868965 # sc_original_timport_med_td.median_td_upper
                                    , td_icu = Inf
                                    , base_date = Date(2020, 1, 15)
                                    , maxtime = 90 #days
                                    )
            plot(result_2.tinf_doubling_time)

            ### Example 3
            ## Load files
            # load sims_simid_tinf_df
            sims_simid_tinf_df_1287570_10 = load( "covidlike-1.4.1-sims-nrep50_sims_simid_tinf_df_1287570.10.jld2", "sims_simid_tinf_df")
            sims_simid_tinf_df_1287570_10_48_98 = sims_simid_tinf_df_1287570_10[48]
            
            # load filtered G df
            sims_1287570_10 = load( "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1287570/1287570_sims_G_filtered_40/covidlike-1.4.1-sims-nrep50_filtered_G_1287570.10.jld2", "sims_G_filtered" )
            sims_1287570_10_48_98 = sims_1287570_10[48] # HPC run 10, sim rep 48 (in HPC run 10), sim rep 98 in combined .jld2 file, 2nd file in folder used to combine
            
            # Run function
            using NBPMscape
            result_3 = infections_time_analysis(; sims_simid_tinf_df_object = sims_simid_tinf_df_1287570_10_48_98
                                    , sims_G_filtered_object = sims_1287570_10_48_98
                                    , sims_max_cases_df_object
                                    #, sim_rep_n
                                    , td_gp = Inf
                                    , td_sc = Inf
                                    , td_icu = 75.09736315548189 # icu_original_timport_med_td.median_td_lower
                                    , base_date = Date(2020, 1, 15)
                                    , maxtime = 90 #days
                                    )
            plot(result_3.tinf_doubling_time)
            println(result_3.numbers_by_day_df )


            ### Example 4
            ## Load files
            # load sims_simid_tinf_df
            sims_simid_tinf_df_1287570_25 = load( "covidlike-1.4.1-sims-nrep50_sims_simid_tinf_df_1287570.25.jld2", "sims_simid_tinf_df")
            sims_simid_tinf_df_1287570_25_46_846 = sims_simid_tinf_df_1287570_25[46]
            
            # load filtered G df
            sims_1287570_25 = load( "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1287570/1287570_sims_G_filtered_40/covidlike-1.4.1-sims-nrep50_filtered_G_1287570.25.jld2", "sims_G_filtered" )
            sims_1287570_25_46_846 = sims_1287570_25[46] # HPC run 25, sim rep 46 (in HPC run 25), sim rep 846 in combined .jld2 file, 17th file in folder used to combine
            
            # Run function
            using NBPMscape
            result_4 = infections_time_analysis(; sims_simid_tinf_df_object = sims_simid_tinf_df_1287570_25_46_846
                                    , sims_G_filtered_object = sims_1287570_25_46_846
                                    , sims_max_cases_df_object
                                    #, sim_rep_n
                                    , td_gp = Inf
                                    , td_sc = Inf
                                    , td_icu = 75.10994255126862 # icu_original_timport_med_td.median_td_lower
                                    , base_date = Date(2020, 1, 15)
                                    , maxtime = 90 #days
                                    )
            plot(result_4.tinf_doubling_time)
            println(result_4.numbers_by_day_df )

            # Example output

            (tinf_vals = [0.0, 0.03478013489480247, 2.4248266293955467, 4.464891428195476, 4.48118165654884, 6.201697887758393, 6.3317832676838925, 6.716844124695573, 8.285459037710579, 8.313550324889103  …  59.99970312954884, 59.999834827986724, 59.99984243122344, 59.999905594742074, 59.999908714564626, 59.99992738727209, 59.99995031226779, 59.999952971419944, 59.999957641587116, 59.99997127157053]
            ,tgp_vals = [11.141408139404426, 11.29260291034947, 13.182479878990822, 17.288673186845028, 19.206619838971037, 19.33464181126501, 24.059684339546457, 24.10565668216914, 24.202279118254065, 24.388865623065083  …  83.70799852949868, 83.76187911541113, 84.22609496726362, 84.37036285859429, 84.81783942075201, 84.8222572841406, 85.26541917912618, 85.70437968943983, 86.27659241401565, 91.85744683366991]
            , ted_vals = [14.690887432529896, 19.789268125849805, 24.27541977085941, 24.53853788944042, 24.868801943087526, 25.658732987813917, 25.974258778383003, 26.063578735462865, 26.15034707510045, 26.316121548065798  …  77.9145477599536, 78.01662623645498, 78.52976789731275, 78.59431411255724, 80.14880150760266, 80.16865182784166, 80.20633334334485, 80.70248942629836, 82.19281531052067, 84.07649284819529]
            , thosp_vals = [18.044942483438767, 24.613624471074356, 24.76422132533206, 26.44744996267905, 26.71853758423846, 26.719262985866838, 27.164600490914548, 27.285677386971766, 27.90663266525602, 28.309723051171424  …  74.54064006070756, 74.91121803418768, 75.20131714206731, 75.25530705731367, 75.45401618373673, 75.72512232963587, 76.18898415304399, 76.32467091302583, 77.31633513291749, 77.73888766834655]
            , ticu_vals = [28.790088511717705, 29.188825111167226, 31.613686583907942, 33.300596585450755, 33.875422336131955, 34.466277965885936, 35.068054108843604, 36.23236992124786, 36.37795120407955, 36.479524160244026  …  72.5637975838133, 72.60952657833035, 73.33750830118223, 74.1158684139547, 74.55209275219406, 74.95585947848632, 75.82455849161357, 76.24669717778605, 79.39632426340467, 81.02297108297965]
            , tdeceased_vals = [32.59998191505391, 33.15956471226622, 34.80748568243381, 37.46874901393653, 38.49439101121762, 38.85768052446824, 39.23563400459135, 39.56773028680076, 39.96550263966643, 42.28483780336335  …  88.1879009200781, 88.30569140211442, 90.44866567419487, 90.70926785749214, 90.70982282730074, 90.84008515740447, 91.12248912081404, 92.28497663009773, 94.1484490951511, 96.22986933740899]
            , tinf_fatal_vals = [20.50654477771215, 22.35431753659192, 23.70692495869025, 24.47857437784323, 25.73701928154305, 26.91840242975509, 27.02352613278702, 27.126895125749762, 27.138158569399714, 27.596627056684838  …  59.961263552085995, 59.96539713708319, 59.9671215149394, 59.98057573618309, 59.9814576687087, 59.985049623540505, 59.986678118464546, 59.989500639430695, 59.999905594742074, 59.999908714564626]
            , timport_vals = [0.0, 2.340837262689952, 3.2215609566280463, 4.2241649886594885, 5.729299305623627, 8.120519883269566, 8.77531308836599, 9.0709250905687, 9.604593066130029, 9.903098253698744  …  80.07679190570906, 80.42503833501675, 82.40192848627491, 82.93511361090285, 83.87843020949067, 83.95777541403224, 83.9934719940582, 85.62213473581303, 86.12425427855253, 87.19774784232376]
            
            ,numbers_by_day_df = 61×18 DataFrame
            Row │ day    date        all_count  all_cumulative  gp_count  gp_cumulative  ed_count  ed_cumulative  hosp_count  hosp_cumulative  icu_count  icu_ ⋯
                │ Int64  Date        Int64      Int64           Int64     Int64          Int64     Int64          Int64       Int64            Int64      Int6 ⋯
            ─────┼───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
            1 │     1  2020-01-15          2               2         0              0         0              0           0                0          0       ⋯
            2 │     2  2020-01-16          0               2         0              0         0              0           0                0          0        
            3 │     3  2020-01-17          1               3         0              0         0              0           0                0          0        
            4 │     4  2020-01-18          0               3         0              0         0              0           0                0          0        
            ⋮  │   ⋮        ⋮           ⋮            ⋮            ⋮            ⋮           ⋮            ⋮            ⋮              ⋮             ⋮           ⋱
            59 │    59  2020-03-13      11088           84910       539           4236       281           2102         221             1692         19       ⋯
            60 │    60  2020-03-14      12452           97362       605           4841       316           2418         263             1955         16        
            61 │    61  2020-03-15          0           97362       738           5579       339           2757         309             2264         25        
            7 columns and 54 rows omitted

            , numbers_by_week_df = 9×17 DataFrame
            Row │ week     all_count  all_cumulative  gp_count  gp_cumulative  ed_count  ed_cumulative  hosp_count  hosp_cumulative  icu_count  icu_cumulative ⋯
                │ String   Int64      Int64           Int64     Int64          Int64     Int64          Int64       Int64            Int64      Int64          ⋯
            ─────┼───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
            1 │ 2020-3           5               5         0              0         0              0           0                0          0               0 ⋯
            2 │ 2020-4          12              17         2              2         0              0           0                0          0               0  
            3 │ 2020-5          77              94         2              4         1              1           1                1          0               0  
            4 │ 2020-6         334             428        11             15         6              7           2                3          0               0  
            ⋮  │    ⋮         ⋮            ⋮            ⋮            ⋮           ⋮            ⋮            ⋮              ⋮             ⋮            ⋮        ⋱
            7 │ 2020-9       10041           15722       563            819       272            385         235              325         23              41 ⋯
            8 │ 2020-10      26508           42230      1320           2139       652           1037         502              827         58              99  
            9 │ 2020-11      55132           97362      3440           5579      1720           2757        1437             2264        119             218  
            6 columns and 2 rows omitted
            
            , number_at_td = 8×4 DataFrame
            Row │ Description                  TD_via_GP  TD_via_SC  TD_via_ICU 
                │ String                       Int64      Int64      Int64
            ─────┼───────────────────────────────────────────────────────────────
            1 │ All_infections                   16084       1966        8874
            2 │ GP_visits                          831         70         441
            3 │ ED_visits                          398         32         207
            4 │ Hospital_admissions                331         27         159
            5 │ ICU_admissions                      42          4          29
            6 │ Deaths                              35          2          13
            7 │ Infections_leading_to_death        218         29         128
            8 │ Imported_infections                574        154         422
            
            , number_at_td3 = 8×4 DataFrame
            Row │ Description                  TD3_via_GP  TD3_via_SC  TD3_via_ICU 
                │ String                       Int64       Int64       Int64
            ─────┼──────────────────────────────────────────────────────────────────
            1 │ All infections                    26536       17543        23504
            2 │ GP visits                          1378         892         1202
            3 │ ED visits                           673         440          593
            4 │ Hospital admissions                 522         356          457
            5 │ ICU admissions                       65          45           57
            6 │ Deaths                               61          42           53
            7 │ Infections_leading_to_death         336         237          308
            8 │ Imported_infections                 692         595          671
            
            , tinf_local_doubling_times_vec = Union{Missing, Float64}[missing, Inf, 1.7095112913514547, Inf, 1.3569154488567239, Inf, 1.4747698473569484, Inf, 3.10628371950539, Inf  …  4.9516773125821025, 4.9512434774593945, 4.79308971453852, 4.990747898273267, 4.800606980753248, 5.042320264207685, 5.030600707165202, 4.9533481626461295, 5.065231789053088, Inf]
            
            , tinf_global_doubling_times_df = 3×2 DataFrame
            Row │ method                          doubling_time_estimates 
                │ String                          Float64
            ─────┼─────────────────────────────────────────────────────────
            1 │ global_log_linear                               3.72596
            2 │ global_wtd_log_linear                           4.83769
            3 │ global_nonlinear_least_squares                  5.19894)
"""
function infections_time_analysis(; sims_simid_tinf_df_object
                                  , sims_G_filtered_object
                                  , sims_max_cases_df_object
                                  #, sim_rep_n
                                  , td_gp = Inf
                                  , td_sc = Inf
                                  , td_icu = Inf
                                  , td3_gp = Inf
                                  , td3_sc = Inf
                                  , td3_icu = Inf
                                  , base_date = Date(2020, 1, 15)
                                  , maxtime = 90 #days
                                 )
    
    ## Time series of infections
    tinf_vals = sort( sims_simid_tinf_df_object[:,:tinf] ) 
    tgp_vals       = sort( sims_G_filtered_object[ isfinite.( sims_G_filtered_object.tgp       ), : ][:,:tgp      ] )
    ted_vals       = sort( sims_G_filtered_object[ isfinite.( sims_G_filtered_object.ted       ), : ][:,:ted      ] )
    thosp_vals     = sort( sims_G_filtered_object[ isfinite.( sims_G_filtered_object.thospital ), : ][:,:thospital] )
    ticu_vals      = sort( sims_G_filtered_object[ isfinite.( sims_G_filtered_object.ticu      ), : ][:,:ticu     ] )
    tdeceased_vals = sort( sims_G_filtered_object[ isfinite.( sims_G_filtered_object.tdeceased ), : ][:,:tdeceased] ) # TODO May be lower because from df filtered for infections with healthcare setting times
    tinf_fatal_vals = sort( sims_G_filtered_object[ (sims_G_filtered_object.fatal) , : ][:,:tinf] ) # TODO May be lower because from df filtered for infections with healthcare setting times# TODO May be lower because from df filtered for infections with healthcare setting times
    timport_vals = sort( sims_max_cases_df_object[ isfinite.( sims_max_cases_df_object.timport ), : ][:,:timport] )

    ## Count number of infections before detection
    # One detection (TD)
    number_at_td = DataFrame(  Description = ["All_infections","GP_visits","ED_visits","Hospital_admissions","ICU_admissions","Deaths","Infections_leading_to_death","Imported_infections"]
                            , TD_via_GP  = [ sum(tinf_vals .<= td_gp), sum(tgp_vals .<= td_gp), sum(ted_vals .<= td_gp), sum(thosp_vals .<= td_gp), sum(ticu_vals .<= td_gp), sum(tdeceased_vals .<= td_gp), sum(tinf_fatal_vals .<= td_gp), sum(timport_vals .<= td_gp) ]
                            , TD_via_SC  = [ sum(tinf_vals .<= td_sc), sum(tgp_vals .<= td_sc), sum(ted_vals .<= td_sc), sum(thosp_vals .<= td_sc), sum(ticu_vals .<= td_sc), sum(tdeceased_vals .<= td_sc), sum(tinf_fatal_vals .<= td_sc), sum(timport_vals .<= td_sc) ]
                            , TD_via_ICU = [ sum(tinf_vals .<= td_icu), sum(tgp_vals .<= td_icu), sum(ted_vals .<= td_icu), sum(thosp_vals .<= td_icu), sum(ticu_vals .<= td_icu), sum(tdeceased_vals .<= td_icu), sum(tinf_fatal_vals .<= td_icu), sum(timport_vals .<= td_icu) ]
                            )
    # Three detections (TD)
    number_at_td3 = DataFrame(  Description = ["All infections","GP visits","ED visits","Hospital admissions","ICU admissions","Deaths","Infections_leading_to_death","Imported_infections"]
                            , TD3_via_GP  = [ sum(tinf_vals .<= td3_gp), sum(tgp_vals .<= td3_gp), sum(ted_vals .<= td3_gp), sum(thosp_vals .<= td3_gp), sum(ticu_vals .<= td3_gp), sum(tdeceased_vals .<= td3_gp), sum(tinf_fatal_vals .<= td3_gp), sum(timport_vals .<= td3_gp) ]
                            , TD3_via_SC  = [ sum(tinf_vals .<= td3_sc), sum(tgp_vals .<= td3_sc), sum(ted_vals .<= td3_sc), sum(thosp_vals .<= td3_sc), sum(ticu_vals .<= td3_sc), sum(tdeceased_vals .<= td3_sc), sum(tinf_fatal_vals .<= td3_sc), sum(timport_vals .<= td3_sc) ]
                            , TD3_via_ICU = [ sum(tinf_vals .<= td3_icu), sum(tgp_vals .<= td3_icu), sum(ted_vals .<= td3_icu), sum(thosp_vals .<= td3_icu), sum(ticu_vals .<= td3_icu), sum(tdeceased_vals .<= td3_icu), sum(tinf_fatal_vals .<= td3_icu), sum(timport_vals .<= td3_icu) ]
                            )
    
    ## Convert infection times to time series of infections by day and week
    
    # Helper function to convert infection times to time series of infection by DAY
    function convert_tinf_to_time_series_days(; times, base_date) base_date = Date(2020,1,15)
    
        # Create df with times of infection and conversion to days
        times_days_df = DataFrame(tinf = times)
        times_days_df.day = Int.( floor.(times) .+ 1 )

        # Group infections by day number and date
        day_infs_df = DataFrame( days = collect(1:maxtime+1))
        day_infs_gdf = groupby( times_days_df, :day)
        day_infs_df_final = combine(day_infs_gdf, nrow => :count)
        day_infs_df_final.cum_cases = cumsum( day_infs_df_final.count )
        day_infs_df_final.date = base_date .+ Day.(day_infs_df_final.day) .- Day(1)

        return( day_infs_df_final )
    end

    # Helper function to convert infection times to time series of infection by WEEK
    function convert_tinf_to_time_series_weeks(; times, base_date)
        
        # Convert times to Dates (integer + fractional days)
        dates = base_date + Dates.Day.(floor.(times)) #+ Dates.Millisecond.(round.(mod.(times, 1) * 86400000))

        # Create DataFrame and count by year-week
        week_infs_df = DataFrame(date = dates)
        week_infs_df.week_groups = string.(year.(dates), '-', week.(dates))

        # Group and count
        week_infs_gdf = groupby(week_infs_df, :week_groups)
        week_infs_df_final = combine(week_infs_gdf, nrow => :count)
        week_infs_df_final.cum_cases = cumsum( week_infs_df_final.count )

        return( week_infs_df_final )
    end

    # Compile daily time series
    tinf_n_by_day       = convert_tinf_to_time_series_days(; times = tinf_vals,       base_date = base_date) # times = [0.1,1.5,2.8]
    tgp_n_by_day        = convert_tinf_to_time_series_days(; times = tgp_vals,        base_date = base_date)
    ted_n_by_day        = convert_tinf_to_time_series_days(; times = ted_vals,        base_date = base_date)
    thosp_n_by_day      = convert_tinf_to_time_series_days(; times = thosp_vals,      base_date = base_date)
    ticu_n_by_day       = convert_tinf_to_time_series_days(; times = ticu_vals,       base_date = base_date)
    tdeceased_n_by_day  = convert_tinf_to_time_series_days(; times = tdeceased_vals,  base_date = base_date)
    tinf_fatal_n_by_day = convert_tinf_to_time_series_days(; times = tinf_fatal_vals, base_date = base_date)
    timport_n_by_day    = convert_tinf_to_time_series_days(; times = timport_vals,    base_date = base_date)


    # Create base df with day and date columns
    days = collect(1:maxtime+1)
    day_date_df = DataFrame(
        day = days,
        date = base_date.+ Day.(days) .-Day(1) # Subtract 1 because day 1  and NOT day 0 is the base date
    )

    # Initialise the df to contain the daily count and cumulative values
    numbers_by_day_df = copy(day_date_df)
    count_types = ["all", "gp", "ed", "hosp", "icu", "deceased", "fatal", "import"]

    for count_type in count_types
        numbers_by_day_df[!, Symbol("$(count_type)_count")] = zeros(Int, maxtime+1)
        numbers_by_day_df[!, Symbol("$(count_type)_cumulative")] = zeros(Int, maxtime+1)
    end

    # Process each data type using a helper function
    function process_count_data!(df, source_df, prefix)
        temp_df = sort(leftjoin(day_date_df, source_df[:, [:day, :count]], on = :day), :day)
        replace!(temp_df.count, missing => 0)
        temp_df.cumulative = cumsum(temp_df.count)
        
        df[!, Symbol("$(prefix)_count")] = Int.(temp_df.count)
        df[!, Symbol("$(prefix)_cumulative")] = Int.(temp_df.cumulative)
    end

    # Apply to all data sources and fill df
    data_sources_day = [
        (tinf_n_by_day, "all")
        ,(tgp_n_by_day, "gp") 
        ,(ted_n_by_day, "ed")
        ,(thosp_n_by_day, "hosp")
        ,(ticu_n_by_day, "icu")
        ,(tdeceased_n_by_day, "deceased")
        ,(tinf_fatal_n_by_day, "fatal")
        ,(timport_n_by_day, "import")
    ]

    for (source_df, prefix) in data_sources_day
        process_count_data!(numbers_by_day_df, source_df, prefix)
    end
    # CSV.write("test.csv", numbers_by_day_df)

    # Produce weekly df by converting daily df
    # Base on daily df
    numbers_by_week_df = copy(numbers_by_day_df)
    # Add weeks
    numbers_by_week_df.week = string.(year.(numbers_by_week_df.date), '-', week.(numbers_by_week_df.date))
    # Remove day and date columns
    select!(numbers_by_week_df, Not([:day,:date]))
    # Sum data by week and reduce to one row per week
    numbers_by_week_df = combine(
                                  groupby(numbers_by_week_df, :week),
                                  names(numbers_by_week_df, Not(:week)) .=> sum .=> names(numbers_by_week_df, Not(:week))
    )
    # Cumulative values will be incorrect after sum, so should be recalculated
    # Define the column pairs
    cols = [:all, :gp, :ed, :hosp, :icu, :deceased, :fatal, :import]
    # Apply cumsum to all count columns and assign to cumulative columns
    for col in cols
        numbers_by_week_df[!, Symbol(col, :_cumulative)] = cumsum(numbers_by_week_df[!, Symbol(col, :_count)])
    end
    
    ## Compute estimates of doubling times using different methods
    tinf_local_doubling_times_vec  = local_doubling_time(t = numbers_by_day_df.day, x = numbers_by_day_df.all_cumulative)
    tinf_global_doubling_times_df = DataFrame( method = ["global_log_linear", "global_wtd_log_linear","global_nonlinear_least_squares"]
                                , doubling_time_estimates = [ global_doubling_time(t = numbers_by_day_df.day, x = numbers_by_day_df.all_cumulative, method=:log_linear)
                                                            , global_doubling_time(t = numbers_by_day_df.day, x = numbers_by_day_df.all_cumulative, method=:weighted_log_linear)
                                                            , global_doubling_time(t = numbers_by_day_df.day, x = numbers_by_day_df.all_cumulative, method=:nonlinear_optim)
                                                            ]
    )


    # Compile output
    output = (  # Simulated times of infection
                tinf_vals  = tinf_vals
              , tgp_vals   = tgp_vals
              , ted_vals   = ted_vals
              , thosp_vals = thosp_vals
              , ticu_vals  = ticu_vals
              , tdeceased_vals  = tdeceased_vals
              , tinf_fatal_vals  = tinf_fatal_vals
              , timport_vals  = timport_vals
              # Daily time series
              , numbers_by_day_df = numbers_by_day_df
              # Weekly time series
              , numbers_by_week_df = numbers_by_week_df
              # Number of infections at first time to detection
              , number_at_td
              # Number of infections at third time to detection
              , number_at_td3
              # Estimated infection doubling times
              , tinf_local_doubling_times_vec
              , tinf_global_doubling_times_df
              )

    return( output ) # output.numbers_by_week_df # output.ticu_vals
end

"""
# old code from function above


    day_date_df = DataFrame( day = collect(1:maxtime+1)
                               , date = base_date .+ Day.(collect(1:maxtime+1))
    )
    
    numbers_by_day_df = DataFrame(  day = collect(1:maxtime+1)
                               , date = base_date .+ Day.(collect(1:maxtime+1))
                               , all_count       = zeros(maxtime+1)
                               , gp_count        = zeros(maxtime+1)
                               , ed_count        = zeros(maxtime+1)
                               , hosp_count      = zeros(maxtime+1)
                               , icu_count       = zeros(maxtime+1)
                               , all_cumulative  = zeros(maxtime+1)
                               , gp_cumulative   = zeros(maxtime+1)
                               , ed_cumulative   = zeros(maxtime+1)
                               , hosp_cumulative = zeros(maxtime+1)                            
                               , icu_cumulative  = zeros(maxtime+1)
    )

    # Add to final df
    # All infections
    tinf_n_by_day_full = sort( leftjoin( day_date_df, tinf_n_by_day[:,[:day,:count]], on = :day ), :day)
    replace!( tinf_n_by_day_full.count, missing => 0)
    tinf_n_by_day_full.cumulative = cumsum( tinf_n_by_day_full.count )
    numbers_by_day_df[:,:all_count] = Int.(tinf_n_by_day_full.count)
    numbers_by_day_df[:,:all_cumulative] = Int.(tinf_n_by_day_full.cumulative)
    # GP visits
    tgp_n_by_day_full = sort( leftjoin( day_date_df, tgp_n_by_day[:,[:day,:count]], on = :day ), :day)
    replace!( tgp_n_by_day_full.count, missing => 0)
    tgp_n_by_day_full.cumulative = cumsum( tgp_n_by_day_full.count )
    numbers_by_day_df[:,:gp_count] = Int.(tgp_n_by_day_full.count)
    numbers_by_day_df[:,:gp_cumulative] = Int.(tgp_n_by_day_full.cumulative)
    # ED visits
    ted_n_by_day_full = sort( leftjoin( day_date_df, ted_n_by_day[:,[:day,:count]], on = :day ), :day)
    replace!( ted_n_by_day_full.count, missing => 0)
    ted_n_by_day_full.cumulative = cumsum( ted_n_by_day_full.count )
    numbers_by_day_df[:,:ed_count] = Int.(ted_n_by_day_full.count)
    numbers_by_day_df[:,:ed_cumulative] = Int.(ted_n_by_day_full.cumulative)
    # Hospital visits
    thosp_n_by_day_full = sort( leftjoin( day_date_df, thosp_n_by_day[:,[:day,:count]], on = :day ), :day)
    replace!( thosp_n_by_day_full.count, missing => 0)
    thosp_n_by_day_full.cumulative = cumsum( thosp_n_by_day_full.count )
    numbers_by_day_df[:,:hosp_count] = Int.(thosp_n_by_day_full.count)
    numbers_by_day_df[:,:hosp_cumulative] = Int.(thosp_n_by_day_full.cumulative)
    # ICU visits
    ticu_n_by_day_full = sort( leftjoin( day_date_df, ticu_n_by_day[:,[:day,:count]], on = :day ), :day)
    replace!( ticu_n_by_day_full.count, missing => 0)
    ticu_n_by_day_full.cumulative = cumsum( ticu_n_by_day_full.count )
    numbers_by_day_df[:,:icu_count] = Int.(ticu_n_by_day_full.count)
    numbers_by_day_df[:,:icu_cumulative] = Int.(ticu_n_by_day_full.cumulative)
    #CSV.write( "test.csv", numbers_by_day_df)

"""

"""
Function        local_doubling_time(t, x)

Description     Return a vector of local doubling times (same length as `x`, with `missing`
                at the first index) assuming exponential growth between observations.
                `t` and `x` must be same length vectors of times and cumulative values.

Arguments       t   Vector of times
                x   Vector of cumulative values

Returns         Vector of doubling time estimate at each time point

Example
                result = local_doubling_time(t = [0,5,10,15,20], x = [1,2,4,8,16])
                plot(result)
                result_1 = local_doubling_time(t = tinf_n_by_day_full.day, x = tinf_n_by_day_full.cumulative)
                mean(result_1)
                plot!(result_1)
                result_2 = local_doubling_time(t = tgp_n_by_day_full.day, x = tgp_n_by_day_full.cumulative)
                plot!(result_2)

"""
function local_doubling_time(;t::AbstractVector, x::AbstractVector)
    length(t) == length(x) || throw(ArgumentError("t and x must have same length"))
    n = length(t)
    Td = Vector{Union{Missing, Float64}}(undef, n)
    Td[1] = missing
    for i in 2:n
        dt = t[i] - t[i-1]
        @assert dt > 0 "times must be strictly increasing"
        ratio = x[i] / x[i-1]
        @assert ratio > 0 "values must be positive"
        # local exponential growth rate r between i-1 and i
        r = log(ratio) / dt
        Td[i] = log(2) / r   # local doubling time
    end
    return Td
end



using LinearAlgebra, Statistics

"""
Function        global_doubling_time(t::AbstractVector, x::AbstractVector; method=:log_linear) -> Float64

Description     Calculates global doubling time, i.e. single value for time series, by
                fitting exponential growth model to entire dataset.

Arguments       `t`: Time vector (must be strictly increasing)
                `x`: Value vector (must be positive)
                `method`: Fitting method (`:log_linear`, `:weighted_log_linear`, or ':nonlinear_optim')

Returns         Global doubling time as Float64, or `Inf` if no growth, `NaN` if decay.

Mathematical Model  Fits x(t) = x₀ * exp(r*t), where doubling time Td = ln(2)/r. 
                    t is time and r is the growth rate.

Examples
                ##Examples
                # Perfect exponential growth
                t1 = [0, 1, 2, 3, 4]
                x1 = [1, 2, 4, 8, 16]
                Td1 = global_doubling_time(t1, x1)  # ≈ 1.0

                # Noisy exponential growth
                t2 = [0, 1, 2, 3, 4]
                x2 = [1.1, 1.9, 4.2, 7.8, 15.9]
                Td2 = global_doubling_time(t2, x2, method=:weighted_log_linear)

                # Compare with local method
                local_times = doubling_time(t1, x1)
                global_time = global_doubling_time(t1, x1)

                println("Local doubling times: ", local_times)
                println("Global doubling time: ", global_time)

                local_times = doubling_time(t2, x2)
                global_time = global_doubling_time(t2, x2)

                println("Local doubling times: ", local_times)
                println("Global doubling time: ", global_time)

                t = tinf_n_by_day_full.day
                x = tinf_n_by_day_full.cumulative
                Td = global_doubling_time(t, x, method = :nonlinear_optim)  # ≈ 1.0

"""
function global_doubling_time(;t::AbstractVector, x::AbstractVector, method=:log_linear) 
    # Input validation 
    length(t) == length(x) || throw(ArgumentError("t and x must have same length")) 
    n = length(t)
    n < 2 && throw(ArgumentError("need at least 2 data points"))

    # Check for positive values
    all(x.> 0) || throw(ArgumentError("all values must be positive for exponential fitting"))

    # Check for monotonic time
    all(diff(t).> 0) || throw(ArgumentError("times must be strictly increasing"))

    if method == :log_linear
        return _fit_log_linear(t, x)
    elseif method == :weighted_log_linear
        return _fit_weighted_log_linear(t, x)
    elseif method == :nonlinear_optim
        return _fit_nonlinear_optim(t, x)
    else
        throw(ArgumentError("method must be :log_linear, :weighted_log_linear, or :nonlinear_optim"))
    end
end

function _fit_log_linear(t::AbstractVector, x::AbstractVector) 
    # Log of cumulative values
    log_x = log.(x)
    # Linear regression: log_x = log(x0) +r*t = a + b*t
    n = length(t)
    t_mean = mean(t)
    log_x_mean = mean(log_x)

    # Calculate slope (growth rate r)
    numerator = sum((t.- t_mean).* (log_x.- log_x_mean))
    denominator = sum((t.- t_mean).^2)

    if abs(denominator) < eps(Float64)
        return Inf  # No time variation
    end

    r = numerator / denominator

    # Convert growth rate to doubling time
    return _growth_rate_to_doubling_time(r)
end

function _fit_weighted_log_linear(t::AbstractVector, x::AbstractVector) 
    # Log of cumulative values
    log_x = log.(x) 
    weights = x # Weight by magnitude of cumulative values

    # Weighted linear regression
    w_sum = sum(weights)
    t_weighted_mean = sum(weights.* t) / w_sum
    log_x_weighted_mean = sum(weights.* log_x) / w_sum

    numerator = sum(weights.* (t.- t_weighted_mean).* (log_x.- log_x_weighted_mean))
    denominator = sum(weights.* (t.- t_weighted_mean).^2)

    # If the value of the denominator is less than the precision of Float64 type then 
    # return doubling time as infinite
    if abs(denominator) < eps(Float64)
        return Inf
    end

    r = numerator / denominator

    # Convert growth rate to doubling time
    return _growth_rate_to_doubling_time(r)
end

function _growth_rate_to_doubling_time(r::Float64) 
    if abs(r) < eps(Float64)
         return Inf # No growth 
    elseif r > 0 
        # Doubling time = log(2) / growth rate
        return log(2) / r # Positive growth 
    else 
        return NaN # Decay (negative growth) 
    end 
end

using Optim
function _fit_nonlinear_optim(t::AbstractVector, x::AbstractVector)
    # Initial parameter guess from log-linear fit
    r_initial = _estimate_initial_r(t, x)
    x0_initial = x[1] * exp(-r_initial * t[1])
    
    # Define objective function (sum of squared residuals)
    function objective(params)
        x0, r = params
        if x0 <= 0
            return Inf  # Invalid parameter
        end
        
        predicted = x0.* exp.(r.* t)
        return sum((x.- predicted).^2)
    end
    
    # Optimize
    initial_params = [x0_initial, r_initial]
    result = optimize(objective, initial_params, BFGS())
    
    if Optim.converged(result)
        _, r_fitted = Optim.minimizer(result)
        return _growth_rate_to_doubling_time(r_fitted)
    else
        @warn "Nonlinear optimization did not converge. Using log-linear fallback."
        return _fit_log_linear(t, x)
    end
end

function _estimate_initial_r(t::AbstractVector, x::AbstractVector)
    # Simple estimate: r ≈ log(x_end/x_start) / (t_end - t_start)
    return log(x[end] / x[1]) / (t[end] - t[1])
end


"""
Function        plot_outbreak_analysis(; data = sim_analysis_1287570_29_4_1004
                                , time_period = "days" # or "weeks"
                                , x_labels = "date" #"number" # or "date")

Description     Plots counts and cumulative values for 8 metrics relating to infections and healthcare seeking

Arguments       data::NamedTuple        Output from infections_time_analysis() function
                time_period::String     Which time variable should be plotted on the x-axis, "days" or "weeks"
                x_labels::String        Which time variable should be labelled on the x-axis, "date" or "number"
                                        Only applies when time_period = "days"

Returns         4x2 grid of plots of counts (LHS y-axis) and cumulative counts (RHS y-axis) against time (x-axis).
                Plots are: Infections, GP Visits, Emergency Department Visits, Hospital admissions, ICU admissions,
                Deaths, Infections leading to death, Infections imported to UK 

Examples

        plot_outbreak_analysis(; data = sim_analysis_1287570_29_4_1004
                                , time_period = "weeks"#"days" # or "weeks"
                                , x_labels = "number" #"date" #"number" # or "date"
                                )

"""
function plot_outbreak_analysis(; data = sim_analysis_1287570_29_4_1004
                                , time_period = "days" # or "weeks"
                                , x_labels = "date" #"number" # or "date"
    )

    @assert time_period in ("days","weeks")
    @assert x_labels in ("number","date")

    # Select df for plotting - daily or weekly data - and x-labels
    if time_period == "days"
        plot_data = data.numbers_by_day_df 
            # Choose either numbers or dates for x vales
        x = x_labels == "number" ? data.numbers_by_day_df.day : data.numbers_by_day_df.date
        xlabel = x_labels == "number" ? "Days since first infection imported to UK" : "First infection imported to UK \n on $(minimum(data.numbers_by_day_df.date))"
    elseif time_period == "weeks"
        plot_data = data.numbers_by_week_df
        # X will be week numbers
        x = plot_data.week
        xlabel = "First infection imported to UK\n on $(minimum(data.numbers_by_day_df.date))"
    end

    # Define data types for plotting
    data_types = ["all","gp","ed","hosp","icu","deceased","fatal","import"]
    plot_titles = ["Infections", "GP visits", "ED visits", "Hospital admissions","ICU admissions","Deaths","Infections leading to death","Infections imported to UK"]


    # Select plot backend
    gr()

    # Store base plots (not twin plots)
    base_plots = []

    for i in 1:length(data_types)
        count_bars = plot_data[:, Symbol("$(data_types[i])_count")]
        cum_line   = plot_data[:, Symbol("$(data_types[i])_cumulative")]
        
        # Create base bar plot
        p_base = bar(x, count_bars
                    ,label = "Daily"
                    ,color = :blue
                    ,alpha = 0.6
                    ,xlabel = xlabel
                    ,ylabel = "Count"
                    ,title = plot_titles[i]
                    ,xrotation = 45
                    ,ylims = (0, maximum(count_bars) * 1.1)  # Start from 0
                    ,bottom_margin = 15Plots.mm
                    ,left_margin = 30Plots.mm
                    ,right_margin = 30Plots.mm
                    ,xtick_direction = :out
                    ,ytick_direction = :out
                )
        
        # Add twin axis - but don't store the twin, store the base
        p_twin = twinx(p_base)
        plot!(p_twin, x, cum_line
            ,label = "Cumulative"
            ,color = :red
            ,linewidth = 2
            ,ylabel = "Cumulative infections"
            ,ylims = (0, maximum(cum_line) * 1.1)  # Start from 0
            ,ytick_direction = :out
            )
        
        # Push the BASE plot, not the twin
        push!(base_plots, p_base)
    end

    # Now combine the base plots
    combined = plot(base_plots..., layout = (2, 4), size = (1400, 800))
    display(combined)

end


"""
Function        plot_outbreak_analysis_combined(; data, time_period, x_labels, td, td3, limit_x_axis )

Description     Plots counts and cumulative values for 8 metrics relating to infections and healthcare seeking
                across 5 plots on a 3x2 grid

Arguments       data::NamedTuple        Output from infections_time_analysis() function (see above)
                time_period::String     Which time variable should be plotted on the x-axis, "days" or "weeks"
                x_labels::String        Which time variable should be labelled on the x-axis, "date" or "number"
                                        (Only applies when time_period = "days")
                td::Float64             Time to 1st detection
                td3::Float64            Time to 3rd detection
                limit_x_axis::Bool      true or false. If true then the x-axis maximum will be td3 + 10 days. This 
                                        makes it easier to see values around the times of detection because values
                                        increase exponentially as the outbreak progresses.

Returns         4x2 grid of plots of counts (LHS y-axis) and cumulative counts (RHS y-axis) against time (x-axis).
                Plots are: Infections, GP Visits, Emergency Department Visits, Hospital admissions, ICU admissions,
                Deaths, Infections leading to death, Infections imported to UK 

Example

        plot_outbreak_analysis_combined(; data = inf_t_analysis_1633695_5_8_188
                                    , time_period = "days", x_labels = "date", limit_x_axis = true
                                    , td = poisson_gleam_timport_med_td.median_td_upper
                                    , td3 = only( tds_1_1633695_median_td_upper.sc_td3 ) )

"""
function plot_outbreak_analysis_combined(; data = sim_analysis_1287570_29_4_1004
                                        , time_period = "days"#"weeks" # or "days"
                                        , x_labels = "date" #"number" #
                                        , td = Inf # td = 50
                                        , td3 = Inf # td3 = 60
                                        , limit_x_axis = true )

    @assert time_period in ("days","weeks")
    @assert x_labels in ("number","date")
    @assert td > 0
    @assert td3 > 0

    # Select df for plotting - daily or weekly data - and x-labels
    if time_period == "days"
        
        # If option to limit x-axis (time period) then filter data to TD3 + 10 days
        plot_data = limit_x_axis == true ? filter( row -> row.day <= td3 + 10, data.numbers_by_day_df ) : data.numbers_by_day_df 
        y1_label = "Count per day"
        
        # Choose either numbers or dates for x values
        #x = x_labels == "number" ? data.numbers_by_day_df.day : data.numbers_by_day_df.date
        #xlabel = x_labels == "number" ? "Days since first infection imported to UK" : "First infection imported to UK \n on $(minimum(data.numbers_by_day_df.date))"
        if x_labels == "number"
            x = plot_data.day #data.numbers_by_day_df.day
            xlabel = "Days since first infection imported to UK"
            td_value = td
            td3_value = td3
            # Need to define xticks as format is particular
            x_labels_n = maximum(x) - minimum(x) +1
            x_labels_step = Int(floor(x_labels_n/10))
            x_tick_positions = maximum( collect(0:x_labels_step:x_labels_n) ) == x_labels_n ? collect(0:x_labels_step:(x_labels_n - 1)) : collect(0:x_labels_step:x_labels_n) #range(1,n,length =k)  #positions = range(1, x_labels_n, length=x_labels_step)  # [1.0, 11.0, ..., 72.0]
            x_tick_labels = x[x_tick_positions.+1] #range(1,n,length =k) 
            bar_xticks = ( x_tick_positions .- 0.5, x_tick_labels)
            # Create importation df by removing days/weeks with zero importations and selecting appropriate day, date, week column for use with x-axis
            importations = filter( row -> row.import_count > 0, plot_data[:,[:day,:import_count]] ) #println(plot_data)

        elseif x_labels == "date"
            x = plot_data.date #data.numbers_by_day_df.date
            xlabel = "First infection imported to UK \n on $(minimum(data.numbers_by_day_df.date))"
            # x_ticks as date need to be formatted correctly
            x_labels_n = Dates.value(maximum(x) - minimum(x)) + 1  # date range
            x_labels_step = Int(floor(x_labels_n/10))   # number of labels
            x_tick_positions = maximum( collect(0:x_labels_step:x_labels_n) ) == x_labels_n ? collect(0:x_labels_step:(x_labels_n - 1)) : collect(0:x_labels_step:x_labels_n) #range(1,n,length =k)  #positions = range(1, x_labels_n, length=x_labels_step)  # [1.0, 11.0, ..., 72.0]
            x_tick_labels = x[x_tick_positions.+1] # collect( minimum(x):Day(7):maximum(x) )
            bar_xticks = ( x_tick_labels, x_tick_labels)  # Map positions to labels
            #x_tick_labels = ( collect( minimum(x):Day(7):maximum(x) ) , Dates.format.( collect( minimum(x):Day(7):maximum(x) ) , dateformat"yyyy-mm-dd") )
            #x_tick_labels = Dates.format.( collect( minimum(x):Day(7):maximum(x) ) , dateformat"yyyy-mm-dd") #( collect( minimum(x):Day(7):maximum(x) ) , Dates.format.( collect( minimum(x):Day(7):maximum(x) ) , dateformat"yyyy-mm-dd") )
            # Create importation df by removing days/weeks with zero importations and selecting appropriate day, date, week column for use with x-axis
            importations = filter( row -> row.import_count > 0, plot_data[:,[:date,:import_count]] ) #println(plot_data)

            if isfinite(td)
                td_value  = minimum(data.numbers_by_day_df.date) + Day( ceil( td ) )
            end
            if isfinite(td3)
                td3_value = minimum(data.numbers_by_day_df.date) + Day( ceil( td3 ) )
            end
        end
    
    elseif time_period == "weeks"
        
        y1_label = "Count per week"

        # If option to limit x-axis (time period) then filter data to TD3 + 10 days
        date_limit = minimum(data.numbers_by_day_df.date) + Day( ceil(td3) + 10 )
        week_limit = string.(year.(date_limit), '-', week.(date_limit)) 
        row_idx = findfirst(==(week_limit), data.numbers_by_week_df.week)

        # Add filtered df to data NamedTuple
        if row_idx === nothing
            data = merge( data, (numbers_by_week_df_sub = data.numbers_by_week_df[1:0, :], ) )  # empty dataframe
        else
            data = merge( data, (numbers_by_week_df_sub = data.numbers_by_week_df[1:row_idx, :], ) )  # rows 1 to row_idx, all columns
        end
        
        plot_data = limit_x_axis == true ? data.numbers_by_week_df_sub : data.numbers_by_week_df 
        #plot_data = data.numbers_by_week_df
        
        # X will be week numbers
        x = plot_data.week
        # x_ticks as date need to be formatted correctly
        x_labels_n = length(x)
        x_labels_step = 1 #Int(floor(x_labels_n/10))   # number of labels
        x_tick_positions = collect(1:x_labels_step:x_labels_n) #range(1,n,length =k)  #positions = range(1, x_labels_n, length=x_labels_step)  # [1.0, 11.0, ..., 72.0]
        x_tick_labels = x[x_tick_positions] # collect( minimum(x):Day(7):maximum(x) )
        bar_xticks = ( x_tick_positions .- 0.5, x_tick_labels)  # Map positions to labels
        xlabel = "Year and week\nFirst infection imported to UK\n on $(minimum(data.numbers_by_day_df.date))"
        # Create importation df by removing days/weeks with zero importations and selecting appropriate day, date, week column for use with x-axis
        importations = filter( row -> row.import_count > 0, plot_data[:,[:week,:import_count]] ) #println(plot_data)

        # Convert times to detection from day numbers to yyyy-ww format
        # TD: time to 1st detection
        if isfinite(td)
            td_date = minimum(data.numbers_by_day_df.date) + Day( ceil(td) )
            td_week = string.(year.(td_date), '-', week.(td_date))
            td_value = findfirst(==(td_week), x)
        end
        # TD3: time to 3rd detection
        if isfinite(td3)
            td3_date = minimum(data.numbers_by_day_df.date) + Day( ceil(td3) )
            td3_week = string.( year.( td3_date ), '-', week.( td3_date ))
            td3_value = findfirst(==(td3_week), x)
        end

    end
    
    # Select plot backend
    gr()

    # Store base plots (not twin plots)
    base_plots = []

    
    ## Plot 1: Importations
    plot_title = "Infections imported to UK"
    count_bars_import = plot_data[:, :import_count     ]
    cum_line_import   = plot_data[:, :import_cumulative]
    importations_y = maximum( count_bars_import ) * 0.1

    # Create base bar plot
    p_import = bar(   x
                    , count_bars_import
                    , label = time_period == "days" ? "Daily (LHS)" : "Weekly (LHS)"
                    , color = :blue, alpha = 0.6
                    #, xlabel = xlabel
                    , ylabel = y1_label
                    , xrotation = 45, xtick_direction = :out, ytick_direction = :out, bottom_margin = 15Plots.mm
                    #, xticks = x_labels == "number" ? xticks : (x_ticks, Dates.format.(x_ticks, dateformat"yyyy-mm-dd"))
                    #, xticks = ( collect(1:length(x_tick_labels)), x_tick_labels)  # Map positions to labels
                    , xticks = bar_xticks  # Map positions to labels
                    , title = plot_title
                    , ylims = (0, maximum(count_bars_import) * 1.1)
                )

    # Add dummy data so can add to legend
    plot!( x[1:2], [0,0], label = "Cumulative (RHS)", color = :red, linewidth = 2)
    
    # Add times to detection
    vline!( [ td_value  ], label = "1st case detected (TD)" )
    vline!( [ td3_value ], label = "3rd case detected (TD3)" )

    # Add import dates
    scatter!(importations[:,1], repeat( [importations_y] , length(importations[:,1]) )
            , marker = :diamond, markersize = sqrt.(importations.import_count) #* 5 # 8
            , markerstrokewidth = 1.5, color = :red, label = "Imported infection(s)")

    # Add twin axis - but don't store the twin, store the base
    p_twin = twinx(p_import)
    
    # Plot cumulative values on secondary y-axis
    plot!(  p_twin
          , x
          , cum_line_import
          , label = "" # "Cumulative (RHS)"
          , color = :red, linewidth = 2
          , ylabel = "Cumulative", ytick_direction = :out
          , ylims = (0, maximum(cum_line_import) * 1.1)
          )
    
    # Position legend
    plot!(legend=(0.1,0.9))
        
    # Push the BASE plot, not the twin
    push!(base_plots, p_import)

    ## Plot 2: Infections
    plot_title = "Infections"
    count_bars_inf = plot_data[:, :all_count     ]
    cum_line_inf   = plot_data[:, :all_cumulative]
    
    importations_y = maximum( count_bars_inf ) * 0.1

    # Create base bar plot
    p_inf = bar(  x
                , count_bars_inf
                , label = time_period == "days" ? "Daily (LHS)" : "Weekly (LHS)"
                , color = :blue, alpha = 0.6
                #, xlabel = xlabel
                , ylabel = y1_label
                , xrotation = 45, xtick_direction = :out, ytick_direction = :out, bottom_margin = 15Plots.mm
                #, xticks = x_labels == "number" ? xticks : (x_ticks, Dates.format.(x_ticks, dateformat"yyyy-mm-dd"))
                , xticks = bar_xticks #(1:length(x_tick_labels), x_tick_labels)  # Map positions to labels
                , title = plot_title
                , ylims = (0, maximum(count_bars_inf) * 1.1)
                )
    
    # Add dummy data so can add to legend
    plot!( x[1:2], [0,0], label = "Cumulative (RHS)", color = :red, linewidth = 2)                
    
    # Add times to detection
    vline!( [ td_value  ], label = "1st case detected (TD)" )
    vline!( [ td3_value ], label = "3rd case detected (TD3)" )

    # Add import dates
    scatter!(importations[:,1], repeat( [importations_y] , length(importations[:,1]) )
            , marker = :diamond, markersize = sqrt.(importations.import_count) #* 5 # 8
            , markerstrokewidth = 1.5, color = :red, label = "Imported infection(s)")

    # Add twin axis - but don't store the twin, store the base
    p_twin = twinx(p_inf)
    
    # Plot cumulative values on secondary y-axis
    plot!( p_twin
         , x
         , cum_line_inf
         , label = "" #"Cumulative (RHS)"
         , color = :red, linewidth = 2
         , ylabel = "Cumulative", ytick_direction = :out
         , ylims = (0, maximum(cum_line_inf) * 1.1)
        )
    
    # Position legend
    plot!(legend=(0.1,0.9))

    # Push the BASE plot, not the twin
    push!(base_plots, p_inf)
    
    ## Plot 3: Infections (log10)
    plot_title = "Infections (log10)"
    count_bars_inf_log10 = log10.(plot_data[:, :all_count     ])
    cum_line_inf_log10   = log10.(plot_data[:, :all_cumulative])
    importations_y = maximum( count_bars_inf_log10 ) * 0.1 # y values for markers of importations

    # Create base bar plot
    p_inf_log10 = bar( x
                     , count_bars_inf_log10
                     , label = time_period == "days" ? "Daily (10^y) (LHS)" : "Weekly (10^y) (LHS)"
                     , color = :blue, alpha = 0.6
                     #, xlabel = xlabel
                     , ylabel = y1_label * " (10^y)"
                     , xrotation = 45, xtick_direction = :out, ytick_direction = :out, bottom_margin = 15Plots.mm
                     #, xticks = x_labels == "number" ? xticks : (x_ticks, Dates.format.(x_ticks, dateformat"yyyy-mm-dd"))
                     , xticks = bar_xticks #(0.5:length(x_tick_labels)-0.5, x_tick_labels)  # Map positions to labels
                     , title = plot_title
                     , ylims = (0, maximum(count_bars_inf_log10) * 1.1)
                    #, yscale = :log10
                     )
    
    # Add dummy data so can add to legend
    plot!( x[1:2], [0,0], label = "Cumulative (RHS)", color = :red, linewidth = 2)

    # Add times to detection
    vline!( [ td_value  ], label = "1st case detected (TD)", linewidth = 2 )
    vline!( [ td3_value ], label = "3rd case detected (TD3)", linewidth = 2 )

    # Add import dates
    scatter!(importations[:,1], repeat( [importations_y] , length(importations[:,1]) )
            , marker = :diamond, markersize = sqrt.(importations.import_count) #* 5 # 8
            , markerstrokewidth = 1.5, color = :red, label = "Imported infection(s)")

    # Add twin axis - but don't store the twin, store the base
    p_twin = twinx(p_inf_log10)
    
    plot!( p_twin
         , x
         , cum_line_inf_log10
         , label = "" #"Cumulative (10^y) (RHS)"
         , color = :red, linewidth = 2
         , ylabel = "Cumulative (10^y)", ytick_direction = :out
         , ylims = (0, maximum(cum_line_inf_log10) * 1.1)
         #, yscale = :log10
        )
    
    # Position legend
    plot!(legend=(0.1,0.9))

    # Push the BASE plot, not the twin
    push!(base_plots, p_inf_log10)

    
    ## Plot 4: GP visits, ED visits, Hospital admissions
    plot_title = "GP visits, ED visits, Hospital admissions"
    
    # GP data
    count_bars_gp = plot_data[:, :gp_count ] #[plot_data[:, :gp_count ], plot_data[:, :ed_count ], plot_data[:, :hosp_count ]]
    cum_line_gp   = plot_data[:, :gp_cumulative]
    
    # Compute y-axis limits
    gp_ed_hosp_max_count = maximum( [ maximum(plot_data[:, :gp_count ]), maximum(plot_data[:, :ed_count ]), maximum(plot_data[:, :hosp_count ]) ] )
    gp_ed_hosp_max_cum = maximum( [ maximum(plot_data[:, :gp_cumulative ]), maximum(plot_data[:, :ed_cumulative ]), maximum(plot_data[:, :hosp_cumulative ]) ] )
    y_limits_count_bars = (0, gp_ed_hosp_max_count * 1.1)  # Start from 0
    y_limits_cum_line   = (0, gp_ed_hosp_max_cum * 1.1)  # Start from 0

    # y-axis positioning for importation markers
    importations_y = y_limits_count_bars[2] * 0.1

    # Create base bar plot
    # using GP visit count data
    p_gp_ed_hosp = bar(   x, count_bars_gp
                        , label = "GP visits" * (time_period == "days" ? " - Daily (LHS)" : " - Weekly (LHS)")
                        , color = :darkred, alpha = 0.6
                        #, xlabel = xlabel
                        , ylabel = y1_label
                        , xrotation = 45, xtick_direction = :out, ytick_direction = :out, bottom_margin = 15Plots.mm
                        #, xticks = x_labels == "number" ? xticks : (x_ticks, Dates.format.(x_ticks, dateformat"yyyy-mm-dd"))
                        , xticks = bar_xticks #(1:length(x_tick_labels), x_tick_labels)  # Map positions to labels
                        #, group = [Symbol(":GP visits (LHS)"), Symbol("GP visits (LHS)"), Symbol("(Hosp visits (LHS)")]
                        #, bar_position = :dodge
                        , title = plot_title
                        , ylims = y_limits_count_bars
                        )
    # Add ED visit count data
    count_bars_ed = plot_data[:, :ed_count     ]
    cum_line_ed   = plot_data[:, :ed_cumulative]
    p_gp_ed_hosp = bar!(  x, count_bars_ed
                        , label = "ED visits" * (time_period == "days" ? " - Daily (LHS)" : " - Weekly (LHS)")
                        , color = :blue2, alpha = 0.6
                        )
    
    # Add hospital admission count data
    count_bars_hosp = plot_data[:, :hosp_count     ]
    cum_line_hosp   = plot_data[:, :hosp_cumulative]
    p_gp_ed_hosp = bar!(  x, count_bars_hosp
                        , label = "Hospital admissions" * (time_period == "days" ? " - Daily (LHS)" : " - Weekly (LHS)")
                        , color = :cyan, alpha = 0.6
                        )
    # Add dummy cumulative data so can add entries to legend
    plot!( x[1:2], [0,0], label = "GP visits - Cumulative (RHS)",           color = :darkred, linewidth = 2)
    plot!( x[1:2], [0,0], label = "ED visits - Cumulative (RHS)",           color = :blue2,   linewidth = 2)
    plot!( x[1:2], [0,0], label = "Hospital admissions - Cumulative (RHS)", color = :cyan,    linewidth = 2)

    # Add times to detection
    vline!( [ td_value  ], label = "1st case detected (TD)", linewidth = 2 )
    vline!( [ td3_value ], label = "3rd case detected (TD3)", linewidth = 2 )

    # Add import dates
    scatter!(importations[:,1], repeat( [importations_y] , length(importations[:,1]) )
            , marker = :diamond, markersize = sqrt.(importations.import_count) #* 5 # 8
            , markerstrokewidth = 1.5, color = :red, label = "Imported infection(s)")

    # Add twin axis - but don't store the twin, store the base
    p_twin = twinx(p_gp_ed_hosp)
    
    # Plot cumulative values on secondary y-axis
    # GP visits
    plot!( p_twin
         , x, cum_line_gp
         , label = "" # "GP visits - Cumulative (RHS)"
         , color = :darkred , linewidth = 2
         , ylabel = "Cumulative", ytick_direction = :out
         , ylims = y_limits_cum_line
        )

    # ED visits
    #p_twin = twinx(p_gp_ed_hosp)
    plot!( p_twin
         , x, cum_line_ed
         , label = "" # "ED visits - Cumulative (RHS)"
         , linewidth = 2
         #, ylabel = "Cumulative", ytick_direction = :out
         , color = :blue2
         , ylims = y_limits_cum_line
        )
    
    # Hospital admissions
    #p_twin = twinx(p_gp_ed_hosp)
    plot!( p_twin
         , x, cum_line_hosp
         , label = "" # "Hospital admissions - Cumulative (RHS)"
         , color = :cyan, linewidth = 2
         #, ylabel = "Cumulative", ytick_direction = :out
         , ylims = y_limits_cum_line
        )

    # Position legend
    plot!(legend=(0.1,0.9))

    # Push the BASE plot, not the twin
    push!(base_plots, p_gp_ed_hosp)


    ## Plot 5: ICU admissions, Deaths, Infections leading to death
    plot_title = "ICU admissions, Deaths, Infections leading to death"
    
    # Define y-axis limits
    icu_deceased_fatal_max_count = maximum( [ maximum(plot_data[:, :icu_count ]), maximum(plot_data[:, :deceased_count ]), maximum(plot_data[:, :fatal_count ]) ] )
    icu_deceased_fatal_max_cum = maximum( [ maximum(plot_data[:, :icu_cumulative ]), maximum(plot_data[:, :deceased_cumulative ]), maximum(plot_data[:, :fatal_cumulative ]) ] )
    y_limits_count_bars = (0, icu_deceased_fatal_max_count * 1.1)  # Start from 0
    y_limits_cum_line   = (0, icu_deceased_fatal_max_cum * 1.1)  # Start from 0
    importations_y = y_limits_count_bars[2] * 0.1
    
    # Create base bar plot
    
    # Fatal infections plot data
    count_bars_fatal = plot_data[:, :fatal_count     ]
    cum_line_fatal   = plot_data[:, :fatal_cumulative]
    
    p_icu_deceased_fatal = bar( x
                                , count_bars_fatal
                                , label = "Infections leading to death" * (time_period == "days" ? " - Daily (LHS)" : " - Weekly (LHS)")
                                , color = :blue, alpha = 0.6
                                #, xlabel = xlabel
                                , ylabel = y1_label
                                #, xticks = x_labels == "number" ? xticks : (x_ticks, Dates.format.(x_ticks, dateformat"yyyy-mm-dd"))
                                , xticks = bar_xticks #(1:length(x_tick_labels), x_tick_labels)  # Map positions to labels
                                , xrotation = 45, bottom_margin = 15Plots.mm, xtick_direction = :out, ytick_direction = :out
                                , title = plot_title
                                , ylims = y_limits_count_bars )

    # ICU plot data
    count_bars_icu = plot_data[:, :icu_count ]
    cum_line_icu   = plot_data[:, :icu_cumulative]
    
    p_icu_deceased_fatal = bar!( x, count_bars_icu
                                , label = "ICU admissions" * (time_period == "days" ? " - Daily (LHS)" : " - Weekly (LHS)")
                                , color = :red, alpha = 0.6 )
    
    # Deceased plot data
    count_bars_deceased = plot_data[:, :deceased_count     ]
    cum_line_deceased   = plot_data[:, :deceased_cumulative]
    
    p_icu_deceased_fatal = bar!( x, count_bars_deceased
                                , label = "Deaths" * (time_period == "days" ? " - Daily (LHS)" : " - Weekly (LHS)")
                                , color = :green, alpha = 0.6 )

    # Add dummy data so can add to legend
    plot!( x[1:2], [0,0], label = "Infections leading to death - Cumulative (RHS)", color = :blue, linewidth = 2)    
    plot!( x[1:2], [0,0], label = "ICU admissions - Cumulative (RHS)", color = :red, linewidth = 2)    
    plot!( x[1:2], [0,0], label = "Deaths - Cumulative (RHS)", color = :lightgreen, linewidth = 2)    
    
    # Add times to detection
    vline!( [ td_value  ], label = "1st case detected (TD)", linewidth = 2 )
    vline!( [ td3_value ], label = "3rd case detected (TD3)", linewidth = 2 )

    # Add import dates
    scatter!(importations[:,1], repeat( [importations_y] , length(importations[:,1]) )
            , marker = :diamond, markersize = sqrt.(importations.import_count) #* 5 # 8
            , markerstrokewidth = 1.5, color = :red, label = "Imported infection(s)")

    # Add twin axis - but don't store the twin, store the base
    p_twin = twinx(p_icu_deceased_fatal)
    
    plot!( p_twin
         , x, cum_line_fatal
         , label = "" #, label = "Infections leading to death - Cumulative (RHS)"
         , color = :blue, linewidth = 2
         , ylabel = "Cumulative", ytick_direction = :out
         , ylims = y_limits_cum_line )

    #p_twin = twinx(p_icu_deceased_fatal)
    plot!( p_twin
         , x, cum_line_icu
         , label = "" #"ICU admissions - Cumulative (RHS)"
         , color = :red, linewidth = 2
         , ylims = y_limits_cum_line )
    
    #p_twin = twinx(p_icu_deceased_fatal)
    plot!( p_twin
         , x, cum_line_deceased
         , label = "" #"Deaths - Cumulative (RHS)"
         , color = :lightgreen, linewidth = 2
         , ytick_direction = :out
         , ylims = y_limits_cum_line )

    # Position legend
    plot!(legend=(0.1,0.9))
    
    # Push the BASE plot, not the twin
    push!(base_plots, p_icu_deceased_fatal)

    ## Combine the base plots
    combined = plot(base_plots..., layout = (2, 3), size = (1800, 800), bottom_margin = 10Plots.mm, left_margin = 10Plots.mm, right_margin = 10Plots.mm) #,top_margin = 10Plots.mm)
    display(combined)

end