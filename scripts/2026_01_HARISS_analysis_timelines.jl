#=
Follow up from '2025_12_HARISS_analysis.jl'

Timeline analysis objectives

1) Show case per week since time 0;
2) re-code time zero as e..g. Jan 15 2020-- just to reinforce the comparison; 
3) Show cumulative imports over time in this simulation, or in addition a simple time of event plot on one axis; 
4) show care-seeking, e.g. cumulative case in gp, hosp, icu , and then 
5) when the actual TD is , which depends on when a case presents to a participating sentinel 

For specific scenarios in previous simulation/analysis compute various outbreak metrics:
- Daily (or weekly) and cumulative number of:
    -- imports, infections, GP visits, ED visits, Hospital admissions, ICU admissions, deaths, infections leading to death
- Values for these metrics at the time of 1st and 3rd detection
- Estimate doubling time for infections
- Combined plot on 3x2 grid showing:
    (1) Infections imported to UK
    (2) Infections
    (3) Infections (log10 scale)
    (4) GP visits, ED visits, Hospital admissions
    (5) ICU admissions, deaths, infections leading to death
    (All with imports per unit time marked on plot)

Summary of steps, some of which are performed by functions in 'src/outbreak_timeline_analysis_functions.jl':
    (A) Load sim output G df filtered for GP, hospital (multiple types) and ICU cases
    (B) Load / compute median TDs for required sampling scenario
    (C) Identify simulation replicate(s) relating to median TD
    (D) Load sims_simid_tinf_df for sim rep relating to median TD
    (E) Compute cumulative cases, total and split by healthcare setting, for particular sim rep 
        (potentially also look at computing the median cumulative cases - probably need to do on HPC as files are large)
    (F) Visualise data

Specific scenario investigated here is scenario 1 from '2025_12_HARISS_analysis.jl': 300 samples per week, equal allocation to HARISS sites, winter, no age targeting.
Results are those for the simulation replicate producing the time to 1st detection (TD) via secondary care HARISS sampling.

See config/HARISS/config file record.xlsx for details of parameters used in different outbreak/sampling scenarios.

This could be extended to look at other sampling scenarios, such as 
    - base case ICU only
    - base case HARISS + ICU combined 
    - shortest all combined (GP, HARISS and ICU)

Timelines for three different methods of timports are also considered:
# (1) Original timport with t-distribution truncated at quantile 0.1% (=1/1000=1/nimports)
# (2) t-distribution truncated at quantile 1%
# (3) Poisson sampling of Mean by outbreak country of mean daily imports from GLEAM

=#

# Load packages
using CSV
using DataFrames
using Dates
using JLD2
using NBPMscape
using Plots


# Comparing timelines for three different methods of timports:
# (1) Original timport with t-distribution truncated at quantile 0.1% (=1/1000=1/nimports) - HPC simulation ref 1287570
# (2) t-distribution truncated at quantile 1% - HPC simulation ref 1540261
# (3) Poisson sampling of Mean by outbreak country of mean daily imports from GLEAM - HPC simulation ref 1633695

#########################################################################################
# (1) Original timport with t-distribution truncated at quantile 0.1% (=1/1000=1/nimports)
#########################################################################################

## Sampling Scenario 1: 300 samples per week, equal allocation to HARISS sites, winter, no age targeting (see config/HARISS/config file record.xlsx for details of scenarios)

# Secondary sampling care via HARISS (original timport distribution)
sc_original_timport_med_td = median_td_sim_reps(; tds_file = "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1287570_nrep2000_analysis/sc_tds_1.csv"
                                                , td_column_name = "SC_TD"
                                                , nreps = 50
                                                )
#(median_td = 51.384143237678416, median_td_lower = 51.38412541748787, median_td_upper = 51.384161057868965
#, median_td_lower_sim_rep_n = 1004, median_td_upper_sim_rep_n = 52, median_td_sim_rep_n = missing
#, file_n_med_lower = 21, sim_rep_n_in_file_n_med_lower = 4
#, file_n_med_upper = 2, sim_rep_n_in_file_n_med_upper = 2
#, file_n_med = missing, sim_rep_n_in_file_n_med = missing)

# Check which raw files these sim reps relate to 
sim_input_folders = [ "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1287570/1287570_sims_G_filtered_40/" ]
sim_files = []
    for i in 1:length(sim_input_folders) #sim_input_folders = ["C:/Users/kdrake/OneDrive - Imperial College London/Documents/mSCAPE/3_results/from_hpc/2025_10_maxtime_maxcases/1096588/1096588_sims_simid_tinf_df_nrep5000/"]
        sim_files_temp = readdir(sim_input_folders[i]; join=true)
        sim_files = vcat(sim_files, sim_files_temp)
    end
println(sim_files)
println(sim_files[sc_original_timport_med_td.file_n_med_lower]) # C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1287570/1287570_sims_G_filtered_40/covidlike-1.4.1-sims-nrep50_filtered_G_1287570.29.jld2

# Files would have been read in by combine_sim_reps in order of first digit
# so sim rep number 52 in the combined file would be sim rep number 2 in covidlike-1.4.1-sims-nrep50_filtered_G_1287570.10.jld2 because this is the 2nd file and each file contains 50 sim reps
# sim rep number 1004 in the combined file would be sim rep number 4 in the 21st file which is covidlike-1.4.1-sims-nrep50_filtered_G_1287570.29.jld2

cases_before_td(; sims_simid_tinf_df_file = "covidlike-1.4.1-sims-nrep50_sims_simid_tinf_df_1287570.29.jld2" # 21st file in folder
                , sim_rep_n = sc_original_timport_med_td.sim_rep_n_in_file_n_med_lower # 4
                , td_value = sc_original_timport_med_td.median_td_lower #51.38412541748787
                )
# 1151

println(sim_files[sc_original_timport_med_td.file_n_med_upper]) # C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1287570/1287570_sims_G_filtered_40/covidlike-1.4.1-sims-nrep50_filtered_G_1287570.10.jld2
cases_before_td(; sims_simid_tinf_df_file = "covidlike-1.4.1-sims-nrep50_sims_simid_tinf_df_1287570.10.jld2"
            , sim_rep_n = sc_original_timport_med_td.sim_rep_n_in_file_n_med_upper # 2
            , td_value = sc_original_timport_med_td.median_td_upper # 51.384161057868965
             )
# 1097

# Read files containing TDs
sc_tds_1 = CSV.read( "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1287570_nrep2000_analysis/sc_tds_1.csv", DataFrame)
gp_tds_1 = CSV.read( "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1287570_nrep2000_analysis/gp_tds_1.csv", DataFrame)
icu_tds_1 = CSV.read( "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1287570_nrep2000_analysis/icu_tds_1.csv", DataFrame)

## Lower median sim rep
sc_td   = sc_tds_1[sc_original_timport_med_td.median_td_lower_sim_rep_n,:].SC_TD # 51.38412541748787
sc_td3  = sc_tds_1[sc_original_timport_med_td.median_td_lower_sim_rep_n,:].SC_3TD # 62.38412541748787
gp_td   = gp_tds_1[sc_original_timport_med_td.median_td_lower_sim_rep_n,:].GP_TD # 68.43986734324332
gp_td3  = gp_tds_1[sc_original_timport_med_td.median_td_lower_sim_rep_n,:].GP_3TD # 74.94585336783484
icu_td  = icu_tds_1[sc_original_timport_med_td.median_td_lower_sim_rep_n,:].ICU_TD # 75.99272546525181
icu_td3 = icu_tds_1[sc_original_timport_med_td.median_td_lower_sim_rep_n,:].ICU_3TD # 88.49696050514001

## Load simulation files
# load sims_simid_tinf_df
sims_simid_tinf_df_1287570_29 = load( "covidlike-1.4.1-sims-nrep50_sims_simid_tinf_df_1287570.29.jld2", "sims_simid_tinf_df") # 21st file in folder
sims_simid_tinf_df_1287570_29_4_1004 = sims_simid_tinf_df_1287570_29[4]

# load filtered G df
#sims_1287570 = load("covidlike-1.4.1-sims-filtered_G_nrep2000_1287570.jld2", "sims") # sim used in HARISS analysis Jan 2026
# or
sims_1287570_29 = load( "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1287570/1287570_sims_G_filtered_40/covidlike-1.4.1-sims-nrep50_filtered_G_1287570.29.jld2", "sims_G_filtered" )
sims_1287570_29_4_1004 = sims_1287570_29[4] # HPC run 29, sim rep 4 (in HPC run 29), sim rep 1004 in combined .jld2 file, 21st file in folder used to combine

# Load max_cases_df
sims_max_cases_df_1287570_29 = load( "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1287570/1287570_max_cases_df/covidlike-1.4.1-sims-nrep50_max_cases_df_1287570.29.jld2", "sims_max_cases_df" )
sims_max_cases_df_1287570_29_4_1004 = sims_max_cases_df_1287570_29[4] # HPC run 29, sim rep 4 (in HPC run 29), sim rep 1004 in combined .jld2 file, 21st file in folder used to combine

# Run function
inf_t_analysis_1287570_29_4_1004 = infections_time_analysis(; sims_simid_tinf_df_object = sims_simid_tinf_df_1287570_29_4_1004
                                                            , sims_G_filtered_object = sims_1287570_29_4_1004
                                                            , sims_max_cases_df_object = sims_max_cases_df_1287570_29_4_1004
                                                            #, sim_rep_n
                                                            , td_gp = 68.43986734324332 # gp_td
                                                            , td_sc = 51.38412541748787 # sc_td and sc_original_timport_med_td.median_td_lower
                                                            , td_icu = 75.99272546525181 # icu_td
                                                            , td3_gp = 74.94585336783484 # gp_td3
                                                            , td3_sc = 62.38412541748787 # sc_td3
                                                            , td3_icu = 88.49696050514001 # icu_td3
                                                            , base_date = Date(2020, 1, 15)
                                                            , maxtime = 90 #days
                                                            )

# Plot
plot_outbreak_analysis_combined(; data = inf_t_analysis_1287570_29_4_1004 #sim_analysis_1287570_29_4_1004
                                , time_period = "weeks", x_labels = "number"
                                , limit_x_axis = true
                                , td = sc_original_timport_med_td.median_td_lower
                                , td3 = sc_td3 )
Plots.savefig( "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1287570_nrep2000_analysis/sim_analysis_HPC1287570_SC_med_TD_lower_array29_simrep4_simrepcombined1004.png" )

plot(inf_t_analysis_1287570_29_4_1004.tinf_local_doubling_times_vec, xlabel = "Time since first importation (days)", ylabel = "Estimated doubling time (days)", label = "Local estimate")
hline!( repeat( [ inf_t_analysis_1287570_29_4_1004.tinf_global_doubling_times_df[3,:doubling_time_estimates] ], length(inf_t_analysis_1287570_29_4_1004.tinf_local_doubling_times_vec) ), label = "Global estimate using nonlinear least squares")
            

# Upper median sim rep
plot_outbreak_analysis_combined(; data = sim_analysis_1287570_29_4_1004
                                , time_period = "days", x_labels = "date", limit_x_axis = true
                                , td = sc_original_timport_med_td.median_td_upper
                                , td3 = Inf)

## Upper median sim rep
## Load files
# load sims_simid_tinf_df
sims_simid_tinf_df_1287570_10 = load( "covidlike-1.4.1-sims-nrep50_sims_simid_tinf_df_1287570.10.jld2", "sims_simid_tinf_df") # 2nd file in folder
sims_simid_tinf_df_1287570_10_2_52 = sims_simid_tinf_df_1287570_10[2]

# load filtered G df
#sims_1287570 = load("covidlike-1.4.1-sims-filtered_G_nrep2000_1287570.jld2", "sims") # sim used in HARISS analysis Jan 2026
#sims_1287570_29_4_1004 = sims_1287570[1004] # HPC run 29, sim rep 4 (in HPC run 29), sim rep 1004 in combined .jld2 file, 21st file in folder used to combine
# or
sims_1287570_10 = load( "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1287570/1287570_sims_G_filtered_40/covidlike-1.4.1-sims-nrep50_filtered_G_1287570.10.jld2", "sims_G_filtered" )
sims_1287570_10_2_52 = sims_1287570_10[2] # HPC run 10, sim rep 2 (in HPC run 10), sim rep 52 in combined .jld2 file, 2nd file in folder used to combine

# Load max_cases_df
sims_max_cases_df_1287570_10 = load( "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1287570/1287570_max_cases_df/covidlike-1.4.1-sims-nrep50_max_cases_df_1287570.10.jld2", "sims_max_cases_df" )
sims_max_cases_df_1287570_10_2_52 = sims_max_cases_df_1287570_29[4] # HPC run 29, sim rep 4 (in HPC run 29), sim rep 1004 in combined .jld2 file, 21st file in folder used to combine

sc_td   = sc_tds_1[sc_original_timport_med_td.median_td_upper_sim_rep_n,:].SC_TD # 51.384161057868965
sc_td3  = sc_tds_1[sc_original_timport_med_td.median_td_upper_sim_rep_n,:].SC_3TD # 65.38416105786897
gp_td   = gp_tds_1[sc_original_timport_med_td.median_td_upper_sim_rep_n,:].GP_TD # 74.70880023426385
gp_td3  = gp_tds_1[sc_original_timport_med_td.median_td_upper_sim_rep_n,:].GP_3TD # 79.69644029376208
icu_td  = icu_tds_1[sc_original_timport_med_td.median_td_upper_sim_rep_n,:].ICU_TD # 77.2867388635035
icu_td3 = icu_tds_1[sc_original_timport_med_td.median_td_upper_sim_rep_n,:].ICU_3TD # 88.37056293484767

# Run function
inf_t_analysis_1287570_10_2_52 = infections_time_analysis(; sims_simid_tinf_df_object = sims_simid_tinf_df_1287570_10_2_52
                                                            , sims_G_filtered_object = sims_1287570_10_2_52
                                                            , sims_max_cases_df_object = sims_max_cases_df_1287570_10_2_52
                                                            #, sim_rep_n
                                                            , td_gp = 68.43986734324332 # gp_td
                                                            , td_sc = 51.38412541748787 # sc_td and sc_original_timport_med_td.median_td_lower
                                                            , td_icu = 75.99272546525181 # icu_td
                                                            , td3_gp = 74.94585336783484 # gp_td3
                                                            , td3_sc = 62.38412541748787 # sc_td3
                                                            , td3_icu = 88.49696050514001 # icu_td3
                                                            , base_date = Date(2020, 1, 15)
                                                            , maxtime = 90 #days
                                                            )
plot(inf_t_analysis_1287570_10_2_52.tinf_local_doubling_times_vec, xlabel = "Time since first importation (days)", ylabel = "Estimated doubling time (days)", label = "Local estimate")
hline!( repeat( [ inf_t_analysis_1287570_10_2_52.tinf_global_doubling_times_df[3,:doubling_time_estimates] ], length(inf_t_analysis_1287570_10_2_52.tinf_local_doubling_times_vec) ), label = "Global estimate using nonlinear least squares")

using NBPMscape
plot_outbreak_analysis_combined(; data = inf_t_analysis_1287570_10_2_52
                                    , time_period = "days", x_labels = "date", limit_x_axis = true
                                    , td = sc_original_timport_med_td.median_td_upper
                                    , td3 = sc_td3 )
Plots.savefig( "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1287570_nrep2000_analysis/sim_analysis_HPC1287570_SC_med_TD_upper_array10_simrep2_simrepcombined52.png" )


### ICU sampling (original timport distribution)
icu_original_timport_med_td = median_td_sim_reps(; tds_file = "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1287570_nrep2000_analysis/icu_tds_1.csv"
                                            , td_column_name = "ICU_TD"
                                            , nreps = 50
                                        )
#(median_td = 75.10365285337525, median_td_lower = 75.09736315548189, median_td_upper = 75.10994255126862
#, median_td_lower_sim_rep_n = 98, median_td_upper_sim_rep_n = 846, median_td_sim_rep_n = missing
#, file_n_med_lower = 2, sim_rep_n_in_file_n_med_lower = 48
#, file_n_med_upper = 17, sim_rep_n_in_file_n_med_upper = 46
#, file_n_med = missing, sim_rep_n_in_file_n_med = missing)

# Find file suffix to use with simid_tinf file
println(sim_files[icu_original_timport_med_td.file_n_med_lower]) # C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1287570/1287570_sims_G_filtered_40/covidlike-1.4.1-sims-nrep50_filtered_G_1287570.10.jld2
cases_before_td(; sims_simid_tinf_df_file = "covidlike-1.4.1-sims-nrep50_sims_simid_tinf_df_1287570.10.jld2"
            , sim_rep_n = icu_original_timport_med_td.sim_rep_n_in_file_n_med_lower
            , td_value = icu_original_timport_med_td.median_td_lower
             )
# 7164

println(sim_files[icu_original_timport_med_td.file_n_med_upper]) # C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1287570/1287570_sims_G_filtered_40/covidlike-1.4.1-sims-nrep50_filtered_G_1287570.25.jld2
cases_before_td(; sims_simid_tinf_df_file = "covidlike-1.4.1-sims-nrep50_sims_simid_tinf_df_1287570.25.jld2"
            , sim_rep_n = icu_original_timport_med_td.sim_rep_n_in_file_n_med_upper
            , td_value = icu_original_timport_med_td.median_td_upper
             )
# 12353


#########################################################################################
# (2) t-distribution truncated at quantile 1%
#########################################################################################
## Scenario 1: 300 samples per week, equal allocation to HARISS sites, winter, no age targeting (see config/HARISS/config file record.xlsx for details of scenarios)

# Secondary sampling care via HARISS (timport distribution truncated at 1% quantile)
sc_trunc_timport_med_td = median_td_sim_reps(; tds_file = "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1540261_nrep330_analysis/sc_tds_1.csv"
                                                , td_column_name = "SC_TD"
                                                , nreps = 10
                                                )
#(median_td = 33.89031011238485, median_td_lower = 33.88279861583523, median_td_upper = 33.89782160893448
#, median_td_lower_sim_rep_n = 42, median_td_upper_sim_rep_n = 152, median_td_sim_rep_n = missing
#, file_n_med_lower = 5, sim_rep_n_in_file_n_med_lower = 2
#, file_n_med_upper = 16, sim_rep_n_in_file_n_med_upper = 2
#, file_n_med = missing, sim_rep_n_in_file_n_med = missing)

# Check which raw files these sim reps relate to 
sim_input_folders = [ "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1540261/1540261_max_cases_df/" ]
sim_files = []
for i in 1:length(sim_input_folders) #sim_input_folders = ["C:/Users/kdrake/OneDrive - Imperial College London/Documents/mSCAPE/3_results/from_hpc/2025_10_maxtime_maxcases/1096588/1096588_sims_simid_tinf_df_nrep5000/"]
    sim_files_temp = readdir(sim_input_folders[i]; join=true)
    sim_files = vcat(sim_files, sim_files_temp)
end

println( sim_files[sc_trunc_timport_med_td.file_n_med_lower] ) # C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1540261/1540261_max_cases_df/covidlike-1.4.2-sims-nrep10_max_cases_df_1540261.2.jld2

cases_before_td(; sims_simid_tinf_df_file = "covidlike-1.4.2-sims-nrep10_sims_simid_tinf_df_1540261.2.jld2"
            , sim_rep_n = sc_trunc_timport_med_td.sim_rep_n_in_file_n_med_lower # 2
            , td_value  = sc_trunc_timport_med_td.median_td_lower #33.88279861583523
             )
# 1036

println(sim_files[sc_trunc_timport_med_td.file_n_med_upper]) # C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1540261/1540261_max_cases_df/covidlike-1.4.2-sims-nrep10_max_cases_df_1540261.41.jld2
cases_before_td(; sims_simid_tinf_df_file = "covidlike-1.4.2-sims-nrep10_sims_simid_tinf_df_1540261.41.jld2"
            , sim_rep_n = sc_trunc_timport_med_td.sim_rep_n_in_file_n_med_upper # 2
            , td_value = sc_trunc_timport_med_td.median_td_upper # 33.89782160893448
             )
# 982

# Read files containing TDs
sc_tds_1_1540261 = CSV.read( "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1540261_nrep330_analysis/sc_tds_1.csv", DataFrame)
gp_tds_1_1540261 = CSV.read( "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1540261_nrep330_analysis/gp_tds_1.csv", DataFrame)
icu_tds_1_1540261 = CSV.read( "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1540261_nrep330_analysis/icu_tds_1.csv", DataFrame)

## Lower median sim rep
sc_td   = sc_tds_1_1540261[ sc_trunc_timport_med_td.median_td_lower_sim_rep_n,:].SC_TD # 33.88279861583523
sc_td3  = sc_tds_1_1540261[ sc_trunc_timport_med_td.median_td_lower_sim_rep_n,:].SC_3TD # 40.88279861583523
gp_td   = gp_tds_1_1540261[ sc_trunc_timport_med_td.median_td_lower_sim_rep_n,:].GP_TD # 40.153293776690376
gp_td3  = gp_tds_1_1540261[ sc_trunc_timport_med_td.median_td_lower_sim_rep_n,:].GP_3TD # 60.07154460821236
icu_td  = icu_tds_1_1540261[sc_trunc_timport_med_td.median_td_lower_sim_rep_n,:].ICU_TD # 29.084285548358007
icu_td3 = icu_tds_1_1540261[sc_trunc_timport_med_td.median_td_lower_sim_rep_n,:].ICU_3TD # 56.4389912636898

## Load files
# load sims_simid_tinf_df
sims_simid_tinf_df_1540261_2 = load( "covidlike-1.4.2-sims-nrep10_sims_simid_tinf_df_1540261.2.jld2", "sims_simid_tinf_df") # 5th file in folder
sims_simid_tinf_df_1540261_2_2_42 = sims_simid_tinf_df_1540261_2[2]

# load filtered G df
#sims_1287570 = load("covidlike-1.4.1-sims-filtered_G_nrep2000_1287570.jld2", "sims") # sim used in HARISS analysis Jan 2026
#sims_1287570_29_4_1004 = sims_1287570[1004] # HPC run 29, sim rep 4 (in HPC run 29), sim rep 1004 in combined .jld2 file, 21st file in folder used to combine
# or
#sims_1540261_2 = load( "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1540261/1540261_sims_G_filtered/covidlike-1.4.1-sims-nrep50_filtered_G_1287570.29.jld2", "sims_G_filtered" )
sims_1540261_2 = load( "covidlike-1.4.2-sims-nrep10_filtered_G_1540261.2.jld2", "sims_G_filtered" )
sims_1540261_2_2_42 = sims_1540261_2[2] # HPC run 2, sim rep 2 (in HPC run 2), sim rep 42 in combined .jld2 file, 5th file in folder used to combine

# Load max_cases_df
sims_max_cases_df_1540261_2 = load( "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1540261/1540261_max_cases_df/covidlike-1.4.2-sims-nrep10_max_cases_df_1540261.2.jld2", "sims_max_cases_df" )
sims_max_cases_df_1540261_2_2_42 = sims_max_cases_df_1540261_2[2] # # HPC run 2, sim rep 2 (in HPC run 2), sim rep 42 in combined .jld2 file, 5th file in folder used to combine

# Run function
inf_t_analysis_1540261_2_2_42 = infections_time_analysis(; sims_simid_tinf_df_object = sims_simid_tinf_df_1540261_2_2_42
                                                            , sims_G_filtered_object = sims_1540261_2_2_42
                                                            , sims_max_cases_df_object = sims_max_cases_df_1540261_2_2_42
                                                            #, sim_rep_n
                                                            , td_gp = 40.153293776690376 # gp_td
                                                            , td_sc = 33.88279861583523 # sc_td and sc_trunc_timport_med_td.median_td_lower
                                                            , td_icu = 29.084285548358007 # icu_td
                                                            , td3_gp = 60.07154460821236 # gp_td3
                                                            , td3_sc = 40.88279861583523 # sc_td3
                                                            , td3_icu = 56.4389912636898 # icu_td3
                                                            , base_date = Date(2020, 1, 15)
                                                            , maxtime = 90 #days
                                                            )

using NBPMscape

minimum( inf_t_analysis_1540261_2_2_42.numbers_by_day_df.import_count )
maximum( inf_t_analysis_1540261_2_2_42.numbers_by_day_df.import_count )

p = plot_outbreak_analysis_combined(; data = inf_t_analysis_1540261_2_2_42
                                    , time_period = "days", x_labels = "date", limit_x_axis = true
                                    , td = sc_trunc_timport_med_td.median_td_lower
                                    , td3 = sc_td3 )
Plots.savefig( "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1540261_nrep330_analysis/sim_analysis_HPC1540261_SC_med_TD_lower_array2_simrep2_simrepcombined42.png" )

plot(inf_t_analysis_1540261_2_2_42.tinf_local_doubling_times_vec, xlabel = "Time since first importation (days)", ylabel = "Estimated doubling time (days)", label = "Local estimate")
hline!( repeat( [ inf_t_analysis_1540261_2_2_42.tinf_global_doubling_times_df[3,:doubling_time_estimates] ], length(inf_t_analysis_1540261_2_2_42.tinf_local_doubling_times_vec) ), label = "Global estimate using nonlinear least squares")
            

## Upper median sim rep
sc_td   = sc_tds_1_1540261[ sc_trunc_timport_med_td.median_td_upper_sim_rep_n,:].SC_TD   # 33.89782160893448
sc_td3  = sc_tds_1_1540261[ sc_trunc_timport_med_td.median_td_upper_sim_rep_n,:].SC_3TD  # 36.89782160893448
gp_td   = gp_tds_1_1540261[ sc_trunc_timport_med_td.median_td_upper_sim_rep_n,:].GP_TD   # 62.52443915611895
gp_td3  = gp_tds_1_1540261[ sc_trunc_timport_med_td.median_td_upper_sim_rep_n,:].GP_3TD  # 65.7794322387617
icu_td  = icu_tds_1_1540261[sc_trunc_timport_med_td.median_td_upper_sim_rep_n,:].ICU_TD  # 33.59618110972124
icu_td3 = icu_tds_1_1540261[sc_trunc_timport_med_td.median_td_upper_sim_rep_n,:].ICU_3TD # 42.19153883577215

## Load files
# load sims_simid_tinf_df
sims_simid_tinf_df_1540261_41 = load( "covidlike-1.4.2-sims-nrep10_sims_simid_tinf_df_1540261.41.jld2", "sims_simid_tinf_df") # 16th file in folder
sims_simid_tinf_df_1540261_41_2_152 = sims_simid_tinf_df_1540261_41[2]

# load filtered G df
sims_1540261_41 = load( "covidlike-1.4.2-sims-nrep10_filtered_G_1540261.41.jld2", "sims_G_filtered" )
sims_1540261_41_2_152 = sims_1540261_41[2] # HPC run 41, sim rep 2 (in HPC run 41), sim rep 152 in combined .jld2 file, 16th file in folder used to combine

# Load max_cases_df
sims_max_cases_df_1540261_41 = load( "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1540261/1540261_max_cases_df/covidlike-1.4.2-sims-nrep10_max_cases_df_1540261.41.jld2", "sims_max_cases_df" )
sims_max_cases_df_1540261_41_2_152 = sims_max_cases_df_1540261_41[2] # HPC run 41, sim rep 2 (in HPC run 41), sim rep 152 in combined .jld2 file, 16th file in folder used to combine

# Run function
inf_t_analysis_1540261_41_2_152 = infections_time_analysis(; sims_simid_tinf_df_object = sims_simid_tinf_df_1540261_41_2_152
                                                            , sims_G_filtered_object   = sims_1540261_41_2_152
                                                            , sims_max_cases_df_object = sims_max_cases_df_1540261_41_2_152
                                                            #, sim_rep_n
                                                            , td_gp = 62.52443915611895 # gp_td
                                                            , td_sc = 33.89782160893448 # sc_td and sc_trunc_timport_med_td.median_td_lower
                                                            , td_icu = 33.59618110972124 # icu_td
                                                            , td3_gp = 65.7794322387617 # gp_td3
                                                            , td3_sc = 36.89782160893448 # sc_td3
                                                            , td3_icu = 42.19153883577215 # icu_td3
                                                            , base_date = Date(2020, 1, 15)
                                                            , maxtime = 90 #days
                                                            )

minimum( inf_t_analysis_1540261_41_2_152.numbers_by_day_df.import_count )
maximum( inf_t_analysis_1540261_41_2_152.numbers_by_day_df.import_count )

p = plot_outbreak_analysis_combined(; data = inf_t_analysis_1540261_2_2_42
                                    , time_period = "days", x_labels = "date", limit_x_axis = true
                                    , td = sc_trunc_timport_med_td.median_td_lower
                                    , td3 = sc_td3 )
Plots.savefig( "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1540261_nrep330_analysis/sim_analysis_HPC1540261_SC_med_TD_upper_array41_simrep2_simrepcombined152.png" )

plot(inf_t_analysis_1540261_41_2_152.tinf_local_doubling_times_vec, xlabel = "Time since first importation (days)", ylabel = "Estimated doubling time (days)", label = "Local estimate")
hline!( repeat( [ inf_t_analysis_1540261_41_2_152.tinf_global_doubling_times_df[3,:doubling_time_estimates] ], length(inf_t_analysis_1540261_41_2_152.tinf_local_doubling_times_vec) ), label = "Global estimate using nonlinear least squares")


#### ICU sampling (timport distribution truncated at 1% quantile)
icu_trunc_timport_med_td = median_td_sim_reps(; tds_file = "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1540261_nrep330_analysis/icu_tds_1.csv"
                                              , td_column_name = "ICU_TD"
                                              , nreps = 10
                                              )
#(median_td = 49.4174091646332, median_td_lower = 49.39978473912122, median_td_upper = 49.43503359014517
#, median_td_lower_sim_rep_n = 161, median_td_upper_sim_rep_n = 277, median_td_sim_rep_n = missing
#, file_n_med_lower = 17, sim_rep_n_in_file_n_med_lower = 1
#, file_n_med_upper = 28, sim_rep_n_in_file_n_med_upper = 7
#, file_n_med = missing, sim_rep_n_in_file_n_med = missing)

# Check which raw files these sim reps relate to 
sim_input_folders = [ "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1540261/1540261_max_cases_df/" ]
sim_files = []
for i in 1:length(sim_input_folders) #sim_input_folders = ["C:/Users/kdrake/OneDrive - Imperial College London/Documents/mSCAPE/3_results/from_hpc/2025_10_maxtime_maxcases/1096588/1096588_sims_simid_tinf_df_nrep5000/"]
    sim_files_temp = readdir(sim_input_folders[i]; join=true)
    sim_files = vcat(sim_files, sim_files_temp)
end

# Find file suffix to use with simid_tinf file
println(sim_files[icu_trunc_timport_med_td.file_n_med_lower]) # C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1540261/1540261_max_cases_df/covidlike-1.4.2-sims-nrep10_max_cases_df_1540261.43.jld2
cases_before_td(; sims_simid_tinf_df_file = "covidlike-1.4.2-sims-nrep10_sims_simid_tinf_df_1540261.43.jld2"
                , sim_rep_n = icu_trunc_timport_med_td.sim_rep_n_in_file_n_med_lower # 1
                , td_value = icu_trunc_timport_med_td.median_td_lower # 49.39978473912122
                )
# 11,163

println(sim_files[icu_trunc_timport_med_td.file_n_med_upper]) # C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1540261/1540261_max_cases_df/covidlike-1.4.2-sims-nrep10_max_cases_df_1540261.82.jld2
cases_before_td(; sims_simid_tinf_df_file = "covidlike-1.4.2-sims-nrep10_sims_simid_tinf_df_1540261.82.jld2"
            , sim_rep_n = icu_trunc_timport_med_td.sim_rep_n_in_file_n_med_upper # 7
            , td_value = icu_trunc_timport_med_td.median_td_upper # 49.43503359014517
             )
# 16,789



#########################################################################################
# (3) Poisson sampling of Mean by outbreak country of mean daily imports from GLEAM
#########################################################################################
# Note that this simulation was run with maxtime = 60 days instead of 90 days for sims (1) and (2) above

## Scenario 1: 300 samples per week, equal allocation to HARISS sites, winter, no age targeting (see config/HARISS/config file record.xlsx for details of scenarios)

# Secondary sampling care via HARISS (Poisson sampling from GLEAM mean daily import output)
poisson_gleam_timport_med_td = median_td_sim_reps(; tds_file = "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1633695_nrep340_analysis/sc_tds_1.csv"
                                                , td_column_name = "SC_TD"
                                                , nreps = 10
                                                )
#(median_td = 33.810776984947765, median_td_lower = 33.78833208911515, median_td_upper = 33.83322188078037, 
# median_td_lower_sim_rep_n = 170, median_td_upper_sim_rep_n = 188, median_td_sim_rep_n = missing, 
# file_n_med_lower = 17, sim_rep_n_in_file_n_med_lower = 10, 
# file_n_med_upper = 19, sim_rep_n_in_file_n_med_upper = 8, 
# file_n_med = missing, sim_rep_n_in_file_n_med = missing)

# Check which raw files these sim reps relate to 
sim_input_folders = [ "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1633695/1633695_max_cases_df/" ]
sim_files = []
for i in 1:length(sim_input_folders) 
    sim_files_temp = readdir(sim_input_folders[i]; join=true)
    sim_files = vcat(sim_files, sim_files_temp)
end

println( sim_files[poisson_gleam_timport_med_td.file_n_med_lower] ) 
# C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1633695/1633695_max_cases_df/covidlike-1.4.3-sims-nrep10_max_cases_df_1633695.38.jld2
println( sim_files[poisson_gleam_timport_med_td.file_n_med_upper] ) 
# C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1633695/1633695_max_cases_df/covidlike-1.4.3-sims-nrep10_max_cases_df_1633695.5.jld2

cases_before_td(; sims_simid_tinf_df_file = "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1633695/1633695_sims_simid_tinf_df/covidlike-1.4.3-sims-nrep10_sims_simid_tinf_df_1633695.38.jld2"
            , sim_rep_n = poisson_gleam_timport_med_td.sim_rep_n_in_file_n_med_lower # 10
            , td_value  = poisson_gleam_timport_med_td.median_td_lower # 33.78833208911515
             )
# 1966

cases_before_td(; sims_simid_tinf_df_file = "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1633695/1633695_sims_simid_tinf_df/covidlike-1.4.3-sims-nrep10_sims_simid_tinf_df_1633695.5.jld2"
            , sim_rep_n = poisson_gleam_timport_med_td.sim_rep_n_in_file_n_med_upper # 8
            , td_value = poisson_gleam_timport_med_td.median_td_upper # 33.83322188078037
             )
# 3842

# Read files containing TDs
sc_tds_1_1633695  = CSV.read( "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1633695_nrep340_analysis/sc_tds_1.csv", DataFrame)
gp_tds_1_1633695  = CSV.read( "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1633695_nrep340_analysis/gp_tds_1.csv", DataFrame)
icu_tds_1_1633695 = CSV.read( "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1633695_nrep340_analysis/icu_tds_1.csv", DataFrame)

## Lower median sim rep
tds_1_1633695 = DataFrame(
                              sc_td   =  sc_tds_1_1633695[ poisson_gleam_timport_med_td.median_td_lower_sim_rep_n,:].SC_TD 
                            , gp_td   =  gp_tds_1_1633695[ poisson_gleam_timport_med_td.median_td_lower_sim_rep_n,:].GP_TD 
                            , icu_td  = icu_tds_1_1633695[ poisson_gleam_timport_med_td.median_td_lower_sim_rep_n,:].ICU_TD 
                            , sc_td3  =  sc_tds_1_1633695[ poisson_gleam_timport_med_td.median_td_lower_sim_rep_n,:].SC_3TD 
                            , gp_td3  =  gp_tds_1_1633695[ poisson_gleam_timport_med_td.median_td_lower_sim_rep_n,:].GP_3TD 
                            , icu_td3 = icu_tds_1_1633695[ poisson_gleam_timport_med_td.median_td_lower_sim_rep_n,:].ICU_3TD
) 

## Load files
# load sims_simid_tinf_df
sims_simid_tinf_df_1633695_38 = load( "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1633695/1633695_sims_simid_tinf_df/covidlike-1.4.3-sims-nrep10_sims_simid_tinf_df_1633695.38.jld2"
                                    , "sims_simid_tinf_df") # 17th file in folder
sims_simid_tinf_df_1633695_38_10_170 = sims_simid_tinf_df_1633695_38[10]

# load filtered G df
sims_1633695_38 = load( "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1633695/1633695_sims_G_filtered/covidlike-1.4.3-sims-nrep10_filtered_G_1633695.38.jld2"
                        , "sims_G_filtered" )
sims_1633695_38_10_170 = sims_1633695_38[10] # HPC run 38, sim rep 10 (in HPC run 38), sim rep 170 in combined .jld2 file, 17th file in folder used to combine

# Load max_cases_df
sims_max_cases_df_1633695_38 = load( "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1540261/1540261_max_cases_df/covidlike-1.4.2-sims-nrep10_max_cases_df_1540261.2.jld2", "sims_max_cases_df" )
sims_max_cases_df_1633695_38_10_170 = sims_max_cases_df_1633695_38[10] # # HPC run 38, sim rep 10 (in HPC run 38), sim rep 170 in combined .jld2 file, 17th file in folder used to combine

# Run function
inf_t_analysis_1633695_38_10_170 = infections_time_analysis(; sims_simid_tinf_df_object = sims_simid_tinf_df_1633695_38_10_170
                                                            , sims_G_filtered_object = sims_1633695_38_10_170
                                                            , sims_max_cases_df_object = sims_max_cases_df_1633695_38_10_170
                                                            #, sim_rep_n
                                                            , td_gp = tds_1_1633695.gp_td
                                                            , td_sc = tds_1_1633695.sc_td
                                                            , td_icu = tds_1_1633695.icu_td
                                                            , td3_gp = tds_1_1633695.gp_td3
                                                            , td3_sc = tds_1_1633695.sc_td3
                                                            , td3_icu = tds_1_1633695.icu_td3
                                                            , base_date = Date(2020, 1, 15)
                                                            , maxtime = 60 #days
                                                            )

# Import count range                                                            
minimum( inf_t_analysis_1633695_38_10_170.numbers_by_day_df.import_count )
maximum( inf_t_analysis_1633695_38_10_170.numbers_by_day_df.import_count )

# Count outputs
inf_t_analysis_1633695_38_10_170.number_at_td
#8×4 DataFrame
# Row │ Description                  TD_via_GP  TD_via_SC  TD_via_ICU 
#     │ String                       Int64      Int64      Int64
#─────┼───────────────────────────────────────────────────────────────
#   1 │ All_infections                   16084       1966        8874
#   2 │ GP_visits                          831         70         441
#   3 │ ED_visits                          398         32         207
#   4 │ Hospital_admissions                331         27         159
#   5 │ ICU_admissions                      42          4          29
#   6 │ Deaths                              35          2          13
#   7 │ Infections_leading_to_death        218         29         128
#   8 │ Imported_infections                574        154         422

inf_t_analysis_1633695_38_10_170.number_at_td3
#8×4 DataFrame
# Row │ Description                  TD3_via_GP  TD3_via_SC  TD3_via_ICU 
#     │ String                       Int64       Int64       Int64
#─────┼──────────────────────────────────────────────────────────────────
#   1 │ All infections                    26536       17543        23504
#   2 │ GP visits                          1378         892         1202
#   3 │ ED visits                           673         440          593
#   4 │ Hospital admissions                 522         356          457
#   5 │ ICU admissions                       65          45           57
#   6 │ Deaths                               61          42           53
#   7 │ Infections_leading_to_death         336         237          308
#   8 │ Imported_infections                 692         595          671

plot_outbreak_analysis_combined(; data = inf_t_analysis_1633695_38_10_170
                                    , time_period = "days", x_labels = "date", limit_x_axis = true
                                    , td = poisson_gleam_timport_med_td.median_td_lower
                                    , td3 = only( tds_1_1633695.sc_td3 ) )
Plots.savefig( "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1633695_nrep340_analysis/sim_analysis_HPC1633695_SC_med_TD_lower_array38_simrep10_simrepcombined170.png" )

# Doubling time estimates
# Local (in time)
inf_t_analysis_1633695_38_10_170.tinf_local_doubling_times_vec
# Global (in time)
inf_t_analysis_1633695_38_10_170.tinf_global_doubling_times_df
#3×2 DataFrame
# Row │ method                          doubling_time_estimates 
#     │ String                          Float64
#─────┼─────────────────────────────────────────────────────────
#   1 │ global_log_linear                               3.72596
#   2 │ global_wtd_log_linear                           4.83769
#   3 │ global_nonlinear_least_squares                  5.19894

plot(inf_t_analysis_1633695_38_10_170.tinf_local_doubling_times_vec, xlabel = "Time since first importation (days)", ylabel = "Estimated doubling time (days)", label = "Local estimate")
hline!( repeat( [ inf_t_analysis_1633695_38_10_170.tinf_global_doubling_times_df[3,:doubling_time_estimates] ], length(inf_t_analysis_1633695_38_10_170.tinf_local_doubling_times_vec) ), label = "Global estimate using nonlinear least squares")
            

## Upper median sim rep
tds_1_1633695_median_td_upper = DataFrame(
                                              sc_td   =  sc_tds_1_1633695[ poisson_gleam_timport_med_td.median_td_upper_sim_rep_n,:].SC_TD 
                                            , gp_td   =  gp_tds_1_1633695[ poisson_gleam_timport_med_td.median_td_upper_sim_rep_n,:].GP_TD 
                                            , icu_td  = icu_tds_1_1633695[ poisson_gleam_timport_med_td.median_td_upper_sim_rep_n,:].ICU_TD 
                                            , sc_td3  =  sc_tds_1_1633695[ poisson_gleam_timport_med_td.median_td_upper_sim_rep_n,:].SC_3TD 
                                            , gp_td3  =  gp_tds_1_1633695[ poisson_gleam_timport_med_td.median_td_upper_sim_rep_n,:].GP_3TD 
                                            , icu_td3 = icu_tds_1_1633695[ poisson_gleam_timport_med_td.median_td_upper_sim_rep_n,:].ICU_3TD
) 

## Load files
# load sims_simid_tinf_df
sims_simid_tinf_df_1633695_5 = load( "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1633695/1633695_sims_simid_tinf_df/covidlike-1.4.3-sims-nrep10_sims_simid_tinf_df_1633695.5.jld2"
                                    , "sims_simid_tinf_df") # 19th file in folder
sims_simid_tinf_df_1633695_5_8_188 = sims_simid_tinf_df_1633695_5[8]

# load filtered G df
sims_1633695_5 = load( "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1633695/1633695_sims_G_filtered/covidlike-1.4.3-sims-nrep10_filtered_G_1633695.5.jld2"
                        , "sims_G_filtered" )
sims_1633695_5_8_188 = sims_1633695_5[8] # HPC run 5, sim rep 8 (in HPC run 8), sim rep 188 in combined .jld2 file, 19th file in folder used to combine

# Load max_cases_df
sims_max_cases_df_1633695_5 = load( "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1540261/1540261_max_cases_df/covidlike-1.4.2-sims-nrep10_max_cases_df_1540261.2.jld2", "sims_max_cases_df" )
sims_max_cases_df_1633695_5_8_188 = sims_max_cases_df_1633695_5[8] # # HPC run 5, sim rep 8 (in HPC run 5), sim rep 188 in combined .jld2 file, 17th file in folder used to combine

# Run function
inf_t_analysis_1633695_5_8_188 = infections_time_analysis(; sims_simid_tinf_df_object = sims_simid_tinf_df_1633695_5_8_188
                                                            , sims_G_filtered_object = sims_1633695_5_8_188
                                                            , sims_max_cases_df_object = sims_max_cases_df_1633695_5_8_188
                                                            #, sim_rep_n
                                                            , td_gp = tds_1_1633695_median_td_upper.gp_td
                                                            , td_sc = tds_1_1633695_median_td_upper.sc_td
                                                            , td_icu = tds_1_1633695_median_td_upper.icu_td
                                                            , td3_gp = tds_1_1633695_median_td_upper.gp_td3
                                                            , td3_sc = tds_1_1633695_median_td_upper.sc_td3
                                                            , td3_icu = tds_1_1633695_median_td_upper.icu_td3
                                                            , base_date = Date(2020, 1, 15)
                                                            , maxtime = 60 #days
                                                            )

# Imports
minimum( inf_t_analysis_1633695_38_10_170.numbers_by_day_df.import_count )
maximum( inf_t_analysis_1633695_38_10_170.numbers_by_day_df.import_count )

plot_outbreak_analysis_combined(; data = inf_t_analysis_1633695_5_8_188
                                    , time_period = "days", x_labels = "date", limit_x_axis = true
                                    , td = poisson_gleam_timport_med_td.median_td_upper
                                    , td3 = only( tds_1_1633695_median_td_upper.sc_td3 ) )
Plots.savefig( "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1633695_nrep340_analysis/sim_analysis_HPC1633695_SC_med_TD_upper_array5_simrep8_simrepcombined188.png" )

# Doubling time estimates
plot(inf_t_analysis_1633695_5_8_188.tinf_local_doubling_times_vec, xlabel = "Time since first importation (days)", ylabel = "Estimated doubling time (days)", label = "Local estimate")
hline!( repeat( [ inf_t_analysis_1633695_5_8_188.tinf_global_doubling_times_df[3,:doubling_time_estimates] ], length(inf_t_analysis_1633695_5_8_188.tinf_local_doubling_times_vec) ), label = "Global estimate using nonlinear least squares")