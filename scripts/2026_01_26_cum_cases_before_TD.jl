#=
21 Jan 2026
Computing number of cases prior to first detection for:
    - HARISS sampling (300 samples per week, winter, equal sample weight to HARISS labs, no age targeting)
    - ICU (8 sites in 2024 pilot) sampling (300 samples per week, winter)
Also, look at the same but with slightly different method of generating imported infections
=#

using CSV 
using DataFrames
using DifferentialEquations
using Distributions
using GLM
using HypothesisTests
using Interpolations
using JLD2
using JumpProcesses 
using LinearAlgebra
using NamedArrays
using NBPMscape
using Optim
#using Pkg
using Plots
using Plots.PlotMeasures 
using QuadGK
using Random
using RData 
using Revise
import SpecialFunctions as SF 
import StatsBase
using Statistics
using StatsPlots
import UUIDs 
using XLSX


#####################################

## Identify simulation replicate relating to median TD and load relevant data

# Load G df
using NBPMscape
initialize_parameters();
#initialize_parameters("config/outbreak_params_covid19_like.yaml");
sims_1287570 = load("covidlike-1.4.1-sims-filtered_G_nrep2000_1287570.jld2", "sims") # sim used in HARISS analysis Jan 2026


## Scenario 1: 300 samples per week, equal allocation to HARISS sites, winter, no age targeting

# Secondary sampling care via HARISS (original timport distribution)
sc_original_timport_med_td = NBPMscape.median_td_sim_reps(; tds_file = "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1287570_nrep2000_analysis/sc_tds_1.csv"
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

NBPMscape.cases_before_td(; sims_simid_tinf_df_file = "covidlike-1.4.1-sims-nrep50_sims_simid_tinf_df_1287570.29.jld2" # 21st file in folder
            , sim_rep_n = sc_original_timport_med_td.sim_rep_n_in_file_n_med_lower # 4
            , td_value = sc_original_timport_med_td.median_td_lower #51.38412541748787
             )
# 1151

println(sim_files[sc_original_timport_med_td.file_n_med_upper]) # C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1287570/1287570_sims_G_filtered_40/covidlike-1.4.1-sims-nrep50_filtered_G_1287570.10.jld2
NBPMscape.cases_before_td(; sims_simid_tinf_df_file = "covidlike-1.4.1-sims-nrep50_sims_simid_tinf_df_1287570.10.jld2"
            , sim_rep_n = sc_original_timport_med_td.sim_rep_n_in_file_n_med_upper # 2
            , td_value = sc_original_timport_med_td.median_td_upper # 51.384161057868965
             )
# 1097


# ICU sampling (original timport distribution)
icu_original_timport_med_td = NBPMscape.median_td_sim_reps(; tds_file = "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1287570_nrep2000_analysis/icu_tds_1.csv"
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
NBPMscape.cases_before_td(; sims_simid_tinf_df_file = "covidlike-1.4.1-sims-nrep50_sims_simid_tinf_df_1287570.10.jld2"
            , sim_rep_n = icu_original_timport_med_td.sim_rep_n_in_file_n_med_lower
            , td_value = icu_original_timport_med_td.median_td_lower # 75.09736315548189
             )
# 7164

println(sim_files[icu_original_timport_med_td.file_n_med_upper]) # C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1287570/1287570_sims_G_filtered_40/covidlike-1.4.1-sims-nrep50_filtered_G_1287570.25.jld2
NBPMscape.cases_before_td(; sims_simid_tinf_df_file = "covidlike-1.4.1-sims-nrep50_sims_simid_tinf_df_1287570.25.jld2"
            , sim_rep_n = icu_original_timport_med_td.sim_rep_n_in_file_n_med_upper
            , td_value = icu_original_timport_med_td.median_td_upper # 75.10994255126862
             )
# 12353


# Secondary sampling care via HARISS (timport distribution truncated at 1% quantile)
sc_trunc_timport_med_td = NBPMscape.median_td_sim_reps(; tds_file = "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1540261_nrep330_analysis/sc_tds_1.csv"
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

NBPMscape.cases_before_td(; sims_simid_tinf_df_file = "covidlike-1.4.2-sims-nrep10_sims_simid_tinf_df_1540261.2.jld2"
            , sim_rep_n = sc_trunc_timport_med_td.sim_rep_n_in_file_n_med_lower # 2
            , td_value = sc_trunc_timport_med_td.median_td_lower #33.88279861583523
             )
# 1036

println(sim_files[sc_trunc_timport_med_td.file_n_med_upper]) # C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1540261/1540261_max_cases_df/covidlike-1.4.2-sims-nrep10_max_cases_df_1540261.41.jld2
NBPMscape.cases_before_td(; sims_simid_tinf_df_file = "covidlike-1.4.2-sims-nrep10_sims_simid_tinf_df_1540261.41.jld2"
            , sim_rep_n = sc_trunc_timport_med_td.sim_rep_n_in_file_n_med_upper # 2
            , td_value = sc_trunc_timport_med_td.median_td_upper # 33.89782160893448
             )
# 982


# ICU sampling (timport distribution truncated at 1% quantile)
icu_trunc_timport_med_td = NBPMscape.median_td_sim_reps(; tds_file = "C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1540261_nrep330_analysis/icu_tds_1.csv"
                                              , td_column_name = "ICU_TD"
                                              , nreps = 10
                                              )
#(median_td = 49.4174091646332, median_td_lower = 49.39978473912122, median_td_upper = 49.43503359014517
#, median_td_lower_sim_rep_n = 161, median_td_upper_sim_rep_n = 277, median_td_sim_rep_n = missing
#, file_n_med_lower = 17, sim_rep_n_in_file_n_med_lower = 1
#, file_n_med_upper = 28, sim_rep_n_in_file_n_med_upper = 7
#, file_n_med = missing, sim_rep_n_in_file_n_med = missing)

# Find file suffix to use with simid_tinf file
println(sim_files[icu_trunc_timport_med_td.file_n_med_lower]) # C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1540261/1540261_max_cases_df/covidlike-1.4.2-sims-nrep10_max_cases_df_1540261.43.jld2
NBPMscape.cases_before_td(; sims_simid_tinf_df_file = "covidlike-1.4.2-sims-nrep10_sims_simid_tinf_df_1540261.43.jld2"
                , sim_rep_n = icu_trunc_timport_med_td.sim_rep_n_in_file_n_med_lower # 1
                , td_value = icu_trunc_timport_med_td.median_td_lower # 49.39978473912122
                )
# 11,163

println(sim_files[icu_trunc_timport_med_td.file_n_med_upper]) # C:/Users/kdrake/Documents/mSCAPE/3_results/from_hpc/2025_12_secondary_care/1540261/1540261_max_cases_df/covidlike-1.4.2-sims-nrep10_max_cases_df_1540261.82.jld2
NBPMscape.cases_before_td(; sims_simid_tinf_df_file = "covidlike-1.4.2-sims-nrep10_sims_simid_tinf_df_1540261.82.jld2"
            , sim_rep_n = icu_trunc_timport_med_td.sim_rep_n_in_file_n_med_upper # 7
            , td_value = icu_trunc_timport_med_td.median_td_upper # 49.43503359014517
             )
# 16,789