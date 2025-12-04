#### Script to investigate the addition of primary care sampling to ICU sampling
### September 2025
#= 
Description:
- Primary care sampling is modelled on the current Oxford-RCGP RSC surveillance
- Script compares ICU sampling vs Primary care sampling vs combination of both
Outline:
1) Load simulation
2) Extract sample of cases
    a) 15% of ICU cases
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
using GLM, Statistics, Distributions, StatsBase
using DataFrames
using CSV 
using Plots, StatsPlots

#using Pkg
#Pkg.add(PackageSpec(name="JLD2", version="0.5.15"))#version="1.11.2"))#version="0.4.54"))
using JLD2
# ] status JLD2
#version = Pkg.TOML.parsefile(joinpath(pkgdir(JLD2), "Project.toml"))["version"]
#println(version)

# Check R for current parameters in NBPMscape
infectivitytoR(2, nsims=10000)
# nsim = 1,000: 2.14 2.19 2.01 2.13 2.16
# nsim = 10,000: 2.07 2.11

# Load simulation
sims = load("covidlike-1.1.1-sims.jld2", "sims")
#@load "covidlike-1.1.1-sims.jld2" sims

## Combine and save multiple files each containing 1000 simulation replicates that have been filtered
## to retain ICU cases only (G dataframe)
sims_G_icu_filter = combine_sim_reps( sim_input_folders = [ "C:/Users/kdrake/OneDrive - Imperial College London/Documents/mSCAPE/3_results/from_hpc/primary_care_Sep_2025/955898/G_filtered_icu"]
                                       #                     ,"C:/Users/kdrake/OneDrive - Imperial College London/Documents/mSCAPE/3_results/from_hpc/primary_care_Sep_2025/955898/G_filtered_gp"]
                                    , sim_object_name = "sims_G_icu_filter" #"sims_G_gp_filter"# "sims" #  
                                    , nrep = 1000
                                    )
@save "covidlike-1.3.1-sims_filtered_G_icu_combined_nrep64000_955898.jld2" sims_G_icu_filter #@save "covidlike-1.3.1-sims_filtered_G_icu_combined_nrep$(nrep_sims_G_icu_filter)_955898.jld2" sims_G_icu_filter

## Combine and save multiple files each containing 1000 simulation replicates that have been filtered 
## to retain GP cases only (G dataframe)
sims_G_gp_filter = combine_sim_reps( sim_input_folders = [ "C:/Users/kdrake/OneDrive - Imperial College London/Documents/mSCAPE/3_results/from_hpc/primary_care_Sep_2025/955898/G_filtered_gp"]
                                    , sim_object_name = "sims_G_gp_filter" #"sims_G_icu_filter" ## "sims" #  
                                    , nrep = 1000
                                    )
@save "covidlike-1.3.1-sims_filtered_G_gp_combined_nrep64000_955898.jld2" sims_G_gp_filter #@save "covidlike-1.3.1-sims_filtered_G_icu_combined_nrep$(nrep_sims_G_gp_filter)_955898.jld2" sims_G_gp_filter

### Check whether can reload files - large files can be difficult to reload
sims_G_icu_filter = load("covidlike-1.3.1-sims_filtered_G_icu_combined_nrep64000_955898.jld2", "sims_G_icu_filter") #load("covidlike-1.3.1-sims_filtered_G_icu_combined_nrep$(nrep_sims_G_icu_filter).jld2", "sims_G_icu_filter")
sims_G_gp_filter = load("covidlike-1.3.1-sims_filtered_G_gp_combined_nrep64000_955898.jld2", "sims_G_gp_filter") #sims_G_gp_filter = load("covidlike-1.3.1-sims_filtered_G_gp_combined_nrep$(nrep_sims_G_gp_filter).jld2", "sims_G_gp_filter")

# If using files created from smaller files, which have less issues with reloading
# Versions BEFORE addition of age disaggregation of severity
#@load "covidlike-1.1.1-sims_filtered_G_icu_combined_nrep53000.jld2" sims_filtered_G_icu
#@load "covidlike-1.1.1-sims_filtered_G_gp_combined_nrep53000.jld2"  sims_filtered_G_gp
# Versions AFTER addition of age disaggregation of severity
sims_G_icu_filter = load("covidlike-1.3.1-sims_filtered_G_icu_combined_nrep64000_955898.jld2", "sims_G_icu_filter")
sims_G_gp_filter = load("covidlike-1.3.1-sims_filtered_G_gp_combined_nrep64000_955898.jld2", "sims_G_gp_filter")

# Obtain population of England from constant in NBPMscape
# (not currently defined in Main)
#ITL225CD_wales = ITL2SIZE[37:39,:ITL225CD]
#ITL2SIZE_eng = ITL2SIZE[.!in(ITL2SIZE.ITL225CD, ITL225CD_wales), :]
#ITL2SIZE_eng = ITL2SIZE[ ITL2SIZE.ITL225CD .!= , :]
#pop_eng = sum(ITL2SIZE[1:36,3])

### Compute median values for each time to detection scenario
### Compute % of simulations that return a time to detection in each time to detection scenario

#### Simplified analysis with simplified function - 7 Sep 2025

# 100 mg samples - summer - current RCGP protocol
TDs_mg100_swab319_ari180 = icu_v_pc_td(; gp_swabs_mg = 100 , gp_ari_swabs = 319, gp_ari_consults = 180, sim_object = "filtered")
CSV.write("scripts/primary_care_v_icu/TDs_mg100_swab319_ari180.csv", TDs_mg100_swab319_ari180) 
TDs_mg100_swab319_ari180 = CSV.read("scripts/primary_care_v_icu/TDs_mg100_swab319_ari180.csv", DataFrame)
TDs_mg100_swab319_ari180_analysis = analyse_td_columns(TDs_mg100_swab319_ari180[:,2:11]) # Not essential but remove the column containing the simulation number

# 100 mg samples - winter - current RCGP protocol
TDs_mg100_swab747_ari327 = icu_v_pc_td(; gp_swabs_mg = 100, gp_ari_swabs = 747, gp_ari_consults = 327, sim_object = "filtered")
CSV.write("scripts/primary_care_v_icu/TDs_mg100_swab747_ari327.csv", TDs_mg100_swab747_ari327) 
TDs_mg100_swab747_ari327_analysis = analyse_td_columns(TDs_mg100_swab747_ari327[:,2:11]) # Not essential but remove the column containing the simulation number

# 200 mg samples - summer - current RCGP protocol
TDs_mg200_swab319_ari180 = icu_v_pc_td(; gp_swabs_mg = 200 , gp_ari_swabs = 319, gp_ari_consults = 180, sim_object = "filtered")
CSV.write("scripts/primary_care_v_icu/TDs_mg200_swab319_ari180.csv", TDs_mg200_swab319_ari180) 
TDs_mg200_swab319_ari180_analysis = analyse_td_columns(TDs_mg200_swab319_ari180[:,2:11]) # Not essential but remove the column containing the simulation number

# 200 mg samples - winter - current RCGP protocol
TDs_mg200_swab747_ari327 = icu_v_pc_td(; gp_swabs_mg = 200, gp_ari_swabs = 747, gp_ari_consults = 327, sim_object = "filtered")
CSV.write("scripts/primary_care_v_icu/TDs_mg200_swab747_ari327.csv", TDs_mg200_swab747_ari327) 
TDs_mg200_swab747_ari327_analysis = analyse_td_columns(TDs_mg200_swab747_ari327[:,2:11]) # Not essential but remove the column containing the simulation number

# 319 mg samples - summer - current RCGP protocol
TDs_mg319_swab319_ari180 = icu_v_pc_td(; gp_swabs_mg = 319 , gp_ari_swabs = 319, gp_ari_consults = 180, sim_object = "filtered")
CSV.write("scripts/primary_care_v_icu/TDs_mg319_swab319_ari180.csv", TDs_mg319_swab319_ari180) 
TDs_mg319_swab319_ari180_analysis = analyse_td_columns(TDs_mg319_swab319_ari180[:,2:11]) # Not essential but remove the column containing the simulation number

# 319 mg samples - winter - current RCGP protocol
TDs_mg319_swab747_ari327 = icu_v_pc_td(; gp_swabs_mg = 319, gp_ari_swabs = 747, gp_ari_consults = 327, sim_object = "filtered")
CSV.write("scripts/primary_care_v_icu/TDs_mg319_swab747_ari327.csv", TDs_mg319_swab747_ari327) 
TDs_mg319_swab747_ari327_analysis = analyse_td_columns(TDs_mg319_swab747_ari327[:,2:11]) # Not essential but remove the column containing the simulation number

# 747 mg samples - summer - More swabs than current RCGP 
TDs_mg747_swab747_ari180 = icu_v_pc_td(; gp_swabs_mg = 747, gp_ari_swabs = 747, gp_ari_consults = 180, sim_object = "filtered")
CSV.write("scripts/primary_care_v_icu/TDs_mg747_swab747_ari180.csv", TDs_mg747_swab747_ari180) 
TDs_mg747_swab747_ari180_analysis = analyse_td_columns(TDs_mg747_swab747_ari180[:,2:11]) # Not essential but remove the column containing the simulation number

# 747 mg samples - winter - current RCGP 
TDs_mg747_swab747_ari327 = icu_v_pc_td(; gp_swabs_mg = 747, gp_ari_swabs = 747, gp_ari_consults = 327, sim_object = "filtered")
CSV.write("scripts/primary_care_v_icu/TDs_mg747_swab747_ari327.csv", TDs_mg747_swab747_ari327) 
TDs_mg747_swab747_ari327_analysis = analyse_td_columns(TDs_mg747_swab747_ari327[:,2:11]) # Not essential but remove the column containing the simulation number

# 1000 mg samples - summer - More swabs than current RCGP - Gu et al (2024) target
TDs_mg1000_swab1000_ari180 = icu_v_pc_td(; gp_swabs_mg = 1000, gp_ari_swabs = 1000, gp_ari_consults = 180, sim_object = "filtered")
CSV.write("scripts/primary_care_v_icu/TDs_mg1000_swab1000_ari180.csv", TDs_mg1000_swab1000_ari180) 
TDs_mg1000_swab1000_ari180_analysis = analyse_td_columns(TDs_mg1000_swab1000_ari180[:,2:11]) # Not essential but remove the column containing the simulation number

# 1000 mg samples - winter - More swabs than current RCGP - Gu et al (2024) target
TDs_mg1000_swab1000_ari327 = icu_v_pc_td(; gp_swabs_mg = 1000, gp_ari_swabs = 1000, gp_ari_consults = 327, sim_object = "filtered")
CSV.write("scripts/primary_care_v_icu/TDs_mg1000_swab1000_ari327.csv", TDs_mg1000_swab1000_ari327) 
TDs_mg1000_swab1000_ari327_analysis = analyse_td_columns(TDs_mg1000_swab1000_ari327[:,2:11]) # Not essential but remove the column containing the simulation number

# 6000 mg samples - summer - More swabs than current RCGP - Leston et al (2022) target
TDs_mg6000_swab6000_ari180 = icu_v_pc_td(; gp_swabs_mg = 6000, gp_ari_swabs = 6000, gp_ari_consults = 180, sim_object = "filtered")
CSV.write("scripts/primary_care_v_icu/TDs_mg6000_swab6000_ari180.csv", TDs_mg6000_swab6000_ari180) 
TDs_mg6000_swab6000_ari180_analysis = analyse_td_columns(TDs_mg6000_swab6000_ari180[:,2:11]) # Not essential but remove the column containing the simulation number

# 6000 mg samples - winter - More swabs than current RCGP - Leston et al (2022) target
TDs_mg6000_swab6000_ari327 = icu_v_pc_td(; gp_swabs_mg = 6000, gp_ari_swabs = 6000, gp_ari_consults = 327, sim_object = "filtered")
CSV.write("scripts/primary_care_v_icu/TDs_mg6000_swab6000_ari327.csv", TDs_mg6000_swab6000_ari327) 
TDs_mg6000_swab6000_ari327_analysis = analyse_td_columns(TDs_mg6000_swab6000_ari327[:,2:11]) # Not essential but remove the column containing the simulation number

# 15419 mg samples - summer - More swabs than current RCGP - 15% of summer ARI GP consultations
TDs_mg15419_swab15419_ari180 = icu_v_pc_td(; gp_swabs_mg = 15419, gp_ari_swabs = 15419, gp_ari_consults = 180, sim_object = "filtered")
CSV.write("scripts/primary_care_v_icu/TDs_mg15419_swab15419_ari180.csv", TDs_mg15419_swab15419_ari180) 
TDs_mg15419_swab15419_ari180_analysis = analyse_td_columns(TDs_mg15419_swab15419_ari180[:,2:11]) # Not essential but remove the column containing the simulation number

# 15419 mg samples - winter - More swabs than current RCGP - 15% of summer ARI GP consultations
TDs_mg15419_swab15419_ari327 = icu_v_pc_td(; gp_swabs_mg = 15419, gp_ari_swabs = 15419, gp_ari_consults = 327, sim_object = "filtered")
CSV.write("scripts/primary_care_v_icu/TDs_mg15419_swab15419_ari327.csv", TDs_mg15419_swab15419_ari327) 
TDs_mg15419_swab15419_ari327_analysis = analyse_td_columns(TDs_mg15419_swab15419_ari327[:,2:11]) # Not essential but remove the column containing the simulation number

# 28011 mg samples - summer - More swabs than current RCGP - 15% of winter ARI GP consultations
TDs_mg28011_swab28011_ari180 = icu_v_pc_td(; gp_swabs_mg = 28011, gp_ari_swabs = 28011, gp_ari_consults = 180, sim_object = "filtered")
CSV.write("scripts/primary_care_v_icu/TDs_mg28011_swab28011_ari180.csv", TDs_mg28011_swab28011_ari180) 
TDs_mg28011_swab28011_ari180_analysis = analyse_td_columns(TDs_mg28011_swab28011_ari180[:,2:11]) # Not essential but remove the column containing the simulation number

# 28011 mg samples - winter - More swabs than current RCGP - 15% of winter ARI GP consultations
TDs_mg28011_swab28011_ari327 = icu_v_pc_td(; gp_swabs_mg = 28011, gp_ari_swabs = 28011, gp_ari_consults = 327, sim_object = "filtered")
CSV.write("scripts/primary_care_v_icu/TDs_mg28011_swab28011_ari327.csv", TDs_mg28011_swab28011_ari327) 
TDs_mg28011_swab28011_ari327_analysis = analyse_td_columns(TDs_mg28011_swab28011_ari327[:,2:11]) # Not essential but remove the column containing the simulation number

# 186738 mg samples - summer - More swabs than current RCGP - 100% of winter ARI GP consultations
TDs_mg186738_swab186738_ari180 = icu_v_pc_td(; gp_swabs_mg = 186738, gp_ari_swabs = 186738, gp_ari_consults = 180, sim_object = "filtered")
CSV.write("scripts/primary_care_v_icu/TDs_mg186738_swab186738_ari180.csv", TDs_mg186738_swab186738_ari180) 
TDs_mg186738_swab186738_ari180_analysis = analyse_td_columns(TDs_mg186738_swab186738_ari180[:,2:11]) # Not essential but remove the column containing the simulation number

# 186738 mg samples - winter - More swabs than current RCGP - 100% of winter ARI GP consultations
TDs_mg186738_swab186738_ari327 = icu_v_pc_td(; gp_swabs_mg = 186738, gp_ari_swabs = 186738, gp_ari_consults = 327, sim_object = "filtered")
CSV.write("scripts/primary_care_v_icu/TDs_mg186738_swab186738_ari327.csv", TDs_mg186738_swab186738_ari327) 
TDs_mg186738_swab186738_ari327_analysis = analyse_td_columns(TDs_mg186738_swab186738_ari327[:,2:11]) # Not essential but remove the column containing the simulation number

plot_hist(df = TDs_mg6000_swab6000_ari180, col = 6)
plot_hist(df = TDs_mg186738_swab186738_ari180, col = 3)
median(TDs_mg186738_swab186738_ari327[:,3])

#### Note that there is convergence in the TD if mg samples and swabs are increased but the number of GPs swabbing is not increased. Prob of sampling is capped at 5%, because you are capturing all the cases at 5% of GPs but none at the 95% of GPs that are not swabbing.
#### Increase gp_practices_swab to match % of ARI consultations being mg sequenced
gp_practices_swab = 300

# 15419 mg samples - summer - More swabs than current RCGP - 15% of summer ARI GP consultations - number of swabbing GPs increased
TDs_mg15419_swab15419_ari180_gp929 = icu_v_pc_td(; gp_swabs_mg = 15419, gp_ari_swabs = 15419, gp_ari_consults = 180, gp_practices_swab = floor(0.15*6199), sim_object = "filtered")
CSV.write("scripts/primary_care_v_icu/TDs_mg15419_swab15419_ari180_gp929.csv", TDs_mg15419_swab15419_ari180_gp929) 
TDs_mg15419_swab15419_ari180_gp929_analysis = analyse_td_columns(TDs_mg15419_swab15419_ari180_gp929[:,2:11]) # Not essential but remove the column containing the simulation number

# 15419 mg samples - winter - More swabs than current RCGP - 15% of summer ARI GP consultations - number of swabbing GPs increased
TDs_mg15419_swab15419_ari327_gp929 = icu_v_pc_td(; gp_swabs_mg = 15419, gp_ari_swabs = 15419, gp_ari_consults = 327, gp_practices_swab = floor(0.15*6199), sim_object = "filtered")
CSV.write("scripts/primary_care_v_icu/TDs_mg15419_swab15419_ari327_gp929.csv", TDs_mg15419_swab15419_ari327_gp929) 
TDs_mg15419_swab15419_ari327_gp929_analysis = analyse_td_columns(TDs_mg15419_swab15419_ari327_gp929[:,2:11]) # Not essential but remove the column containing the simulation number

# 28011 mg samples - summer - More swabs than current RCGP - 15% of winter ARI GP consultations - number of swabbing GPs increased
TDs_mg28011_swab28011_ari180_gp929 = icu_v_pc_td(; gp_swabs_mg = 28011, gp_ari_swabs = 28011, gp_ari_consults = 180, gp_practices_swab = floor(0.15*6199), sim_object = "filtered")
CSV.write("scripts/primary_care_v_icu/TDs_mg28011_swab28011_ari180_gp929.csv", TDs_mg28011_swab28011_ari180_gp929) 
TDs_mg28011_swab28011_ari180_gp929_analysis = analyse_td_columns(TDs_mg28011_swab28011_ari180_gp929[:,2:11]) # Not essential but remove the column containing the simulation number

# 28011 mg samples - winter - More swabs than current RCGP - 15% of winter ARI GP consultations - number of swabbing GPs increased
TDs_mg28011_swab28011_ari327_gp929 = icu_v_pc_td(; gp_swabs_mg = 28011, gp_ari_swabs = 28011, gp_ari_consults = 327, gp_practices_swab = floor(0.15*6199), sim_object = "filtered")
CSV.write("scripts/primary_care_v_icu/TDs_mg28011_swab28011_ari327_gp929.csv", TDs_mg28011_swab28011_ari327_gp929) 
TDs_mg28011_swab28011_ari327_gp929_analysis = analyse_td_columns(TDs_mg28011_swab28011_ari327_gp929[:,2:11]) # Not essential but remove the column containing the simulation number

# 102792 mg samples - summer - More swabs than current RCGP - 100% of summer ARI GP consultations - number of swabbing GPs kept at 300
TDs_mg102792_swab102792_ari180_gp300 = icu_v_pc_td(; gp_swabs_mg = 102792, gp_ari_swabs = 102792, gp_ari_consults = 180, gp_practices_swab = 300, sim_object = "filtered")
CSV.write("scripts/primary_care_v_icu/TDs_mg102792_swab102792_ari180_gp300.csv", TDs_mg102792_swab102792_ari180_gp300) 
TDs_mg102792_swab102792_ari180_gp300_analysis = analyse_td_columns(TDs_mg102792_swab102792_ari180_gp300[:,2:11]) # Not essential but remove the column containing the simulation number

# 102792 mg samples - summer - More swabs than current RCGP - 100% of summer ARI GP consultations - number of swabbing GPs increased
TDs_mg102792_swab102792_ari180_gp6199 = icu_v_pc_td(; gp_swabs_mg = 102792, gp_ari_swabs = 102792, gp_ari_consults = 180, gp_practices_swab = 6199, sim_object = "filtered")
CSV.write("scripts/primary_care_v_icu/TDs_mg102792_swab102792_ari180_gp6199.csv", TDs_mg102792_swab102792_ari180_gp6199) 
TDs_mg102792_swab102792_ari180_gp6199_analysis = analyse_td_columns(TDs_mg102792_swab102792_ari180_gp6199[:,2:11]) # Not essential but remove the column containing the simulation number

# 102792 mg samples - winter - More swabs than current RCGP - 100% of summer ARI GP consultations - number of swabbing GPs kept at 300
TDs_mg102792_swab102792_ari180_gp300 = icu_v_pc_td(; gp_swabs_mg = 102792, gp_ari_swabs = 102792, gp_ari_consults = 327, gp_practices_swab = 300, sim_object = "filtered")
CSV.write("scripts/primary_care_v_icu/TDs_mg102792_swab102792_ari180_gp300.csv", TDs_mg102792_swab102792_ari180_gp300) 
TDs_mg102792_swab102792_ari180_gp300_analysis = analyse_td_columns(TDs_mg102792_swab102792_ari180_gp300[:,2:11]) # Not essential but remove the column containing the simulation number

# 102792 mg samples - winter - More swabs than current RCGP - 100% of summer ARI GP consultations - number of swabbing GPs kept at 300
TDs_mg102792_swab102792_ari327_gp300 = icu_v_pc_td(; gp_swabs_mg = 102792, gp_ari_swabs = 102792, gp_ari_consults = 327, gp_practices_swab = 300, sim_object = "filtered")
CSV.write("scripts/primary_care_v_icu/TDs_mg102792_swab102792_ari327_gp300.csv", TDs_mg102792_swab102792_ari327_gp300) 
TDs_mg102792_swab102792_ari327_gp300_analysis = analyse_td_columns(TDs_mg102792_swab102792_ari327_gp300[:,2:11]) # Not essential but remove the column containing the simulation number

# 102792 mg samples - winter - More swabs than current RCGP - 100% of summer ARI GP consultations - number of swabbing GPs increased
TDs_mg102792_swab102792_ari327_gp6199 = icu_v_pc_td(; gp_swabs_mg = 102792, gp_ari_swabs = 102792, gp_ari_consults = 327, gp_practices_swab = 6199, sim_object = "filtered")
CSV.write("scripts/primary_care_v_icu/TDs_mg102792_swab102792_ari327_gp6199.csv", TDs_mg102792_swab102792_ari327_gp6199) 
TDs_mg102792_swab102792_ari327_gp6199_analysis = analyse_td_columns(TDs_mg102792_swab102792_ari327_gp6199[:,2:11]) # Not essential but remove the column containing the simulation number

# 186738 mg samples - winter - More swabs than current RCGP - 100% of winter ARI GP consultations - number of swabbing GPs increased
TDs_mg186738_swab186738_ari327_gp6199 = icu_v_pc_td(; gp_swabs_mg = 186738, gp_ari_swabs = 186738, gp_ari_consults = 327, gp_practices_swab = 6199, sim_object = "filtered")
CSV.write("scripts/primary_care_v_icu/TDs_mg186738_swab186738_ari327_gp6199.csv", TDs_mg186738_swab186738_ari327_gp6199) 
TDs_mg186738_swab186738_ari327_gp6199_analysis = analyse_td_columns(TDs_mg186738_swab186738_ari327_gp6199[:,2:11]) # Not essential but remove the column containing the simulation number

# Create compilation of data
# Vector of dataframes containing results for summer
TD_results_dfs_summer = [TDs_mg100_swab319_ari180_analysis[1:6,:]
                        , TDs_mg200_swab319_ari180_analysis[1:6,:]
                        , TDs_mg319_swab319_ari180_analysis[1:6,:]
                        , TDs_mg747_swab747_ari180_analysis[1:6,:]
                        , TDs_mg1000_swab1000_ari180_analysis[1:6,:]
                        , TDs_mg6000_swab6000_ari180_analysis[1:6,:]
                        , TDs_mg15419_swab15419_ari180_gp929_analysis[1:6,:]
                        , TDs_mg28011_swab28011_ari180_gp929_analysis[1:6,:]
                        , TDs_mg102792_swab102792_ari180_gp6199_analysis[1:6,:]
                        , TDs_mg186738_swab186738_ari327_gp6199_analysis[1:6,:] # placeholder values to be replaced
                        ] 

# Vector of dataframes containing results for summer
TD_results_dfs_winter = [ TDs_mg100_swab747_ari327_analysis[1:6,:]
                        , TDs_mg200_swab747_ari327_analysis[1:6,:]
                        , TDs_mg319_swab747_ari327_analysis[1:6,:]
                        , TDs_mg747_swab747_ari327_analysis[1:6,:]
                        , TDs_mg1000_swab1000_ari327_analysis[1:6,:]
                        , TDs_mg6000_swab6000_ari327_analysis[1:6,:]
                        , TDs_mg15419_swab15419_ari327_gp929_analysis[1:6,:]
                        , TDs_mg28011_swab28011_ari327_gp929_analysis[1:6,:]
                        , TDs_mg102792_swab102792_ari327_gp6199_analysis[1:6,:]
                        , TDs_mg186738_swab186738_ari327_gp6199_analysis[1:6,:]
                        ] 

# Create a new dataframe with results for times to detection of 1 case
results_df = DataFrame( Number_of_PC_mg_samples = [100,200,319,747,1000,6000,15419,28011,102792,186738])
results_df.ICU_only           = [df[1, :Median_TD] for df in TD_results_dfs_summer]
results_df.PC_only_summer     = [df[2, :Median_TD] for df in TD_results_dfs_summer]
results_df.Combined_summer    = [df[3, :Median_TD] for df in TD_results_dfs_summer]
results_df.Improvement_summer = [df[3, :Median_TD] for df in TD_results_dfs_summer] - [df[1, :Median_TD] for df in TD_results_dfs_summer]
results_df.PC_only_winter     = [df[2, :Median_TD] for df in TD_results_dfs_winter]
results_df.Combined_winter    = [df[3, :Median_TD] for df in TD_results_dfs_winter]
results_df.Improvement_winter = [df[3, :Median_TD] for df in TD_results_dfs_winter] - [df[1, :Median_TD] for df in TD_results_dfs_winter]
results_df[10,3:5] = [-1,-1,-1]
println(results_df)
CSV.write("scripts/primary_care_v_icu/TDs_1case_w_gp_adj.csv", results_df) 


## Results when number of GPs kept at 300
# Vector of dataframes containing results for summer
TD_results_dfs_summer_gp300 = [ TDs_mg15419_swab15419_ari180_analysis[1:6,:]
                            , TDs_mg28011_swab28011_ari180_analysis[1:6,:]
                            , TDs_mg102792_swab102792_ari180_gp300_analysis[1:6,:]
                            , TDs_mg186738_swab186738_ari327_analysis[1:6,:] # placeholder values to be replaced
                        ] 

# Vector of dataframes containing results for summer
TD_results_dfs_winter_gp300 = [ TDs_mg15419_swab15419_ari327_analysis[1:6,:]
                        , TDs_mg28011_swab28011_ari327_analysis[1:6,:]
                        , TDs_mg102792_swab102792_ari327_gp300_analysis[1:6,:]
                        , TDs_mg186738_swab186738_ari327_analysis[1:6,:]
                        ] 

# Create a new dataframe with results for times to detection of 1 case
results_df_gp300 = DataFrame( Number_of_PC_mg_samples = [15419,28011,102792,186738])
results_df_gp300.ICU_only           = [df[1, :Median_TD] for df in TD_results_dfs_summer_gp300]
results_df_gp300.PC_only_summer     = [df[2, :Median_TD] for df in TD_results_dfs_summer_gp300]
results_df_gp300.Combined_summer    = [df[3, :Median_TD] for df in TD_results_dfs_summer_gp300]
results_df_gp300.Improvement_summer = [df[3, :Median_TD] for df in TD_results_dfs_summer_gp300] - [df[1, :Median_TD] for df in TD_results_dfs_summer_gp300]
results_df_gp300.PC_only_winter     = [df[2, :Median_TD] for df in TD_results_dfs_winter_gp300]
results_df_gp300.Combined_winter    = [df[3, :Median_TD] for df in TD_results_dfs_winter_gp300]
results_df_gp300.Improvement_winter = [df[3, :Median_TD] for df in TD_results_dfs_winter_gp300] - [df[1, :Median_TD] for df in TD_results_dfs_winter_gp300]
results_df_gp300[4,3:5] = [-1,-1,-1]
println(results_df_gp300)
CSV.write("scripts/primary_care_v_icu/TDs_1case_w_gp_300.csv", results_df_gp300) 

# Repeat for time to detection of 3 cases (3TD)
# Create a new dataframe with results for times to detection of 3 case
results_3td_df = DataFrame( Number_of_PC_mg_samples = [100,200,319,747,1000,6000,15419,28011,102792,186738])
results_3td_df.ICU_only           = [df[4, :Median_TD] for df in TD_results_dfs_summer]
results_3td_df.PC_only_summer     = [df[5, :Median_TD] for df in TD_results_dfs_summer]
results_3td_df.Combined_summer    = [df[6, :Median_TD] for df in TD_results_dfs_summer]
results_3td_df.Improvement_summer = [df[6, :Median_TD] for df in TD_results_dfs_summer] - [df[4, :Median_TD] for df in TD_results_dfs_summer]
results_3td_df.PC_only_winter     = [df[5, :Median_TD] for df in TD_results_dfs_winter]
results_3td_df.Combined_winter    = [df[6, :Median_TD] for df in TD_results_dfs_winter]
results_3td_df.Improvement_winter = [df[6, :Median_TD] for df in TD_results_dfs_winter] - [df[4, :Median_TD] for df in TD_results_dfs_winter]
results_3td_df[10,3:5] = [-1,-1,-1]
println(results_3td_df)
CSV.write("scripts/primary_care_v_icu/TDs_3cases_w_gp_adj.csv", results_3td_df) 

## Results when number of GPs kept at 300
# Vector of dataframes containing results for summer
TD_results_dfs_summer_gp300 = [ TDs_mg15419_swab15419_ari180_analysis[1:6,:]
                            , TDs_mg28011_swab28011_ari180_analysis[1:6,:]
                            , TDs_mg102792_swab102792_ari180_gp300_analysis[1:6,:]
                            , TDs_mg186738_swab186738_ari327_analysis[1:6,:] # placeholder values to be replaced
                        ] 

# Vector of dataframes containing results for summer
TD_results_dfs_winter_gp300 = [ TDs_mg15419_swab15419_ari327_analysis[1:6,:]
                        , TDs_mg28011_swab28011_ari327_analysis[1:6,:]
                        , TDs_mg102792_swab102792_ari327_gp300_analysis[1:6,:]
                        , TDs_mg186738_swab186738_ari327_analysis[1:6,:]
                        ] 

# Create a new dataframe with results for times to detection of 1 case
results_df_gp300 = DataFrame( Number_of_PC_mg_samples = [15419,28011,102792,186738])
results_df_gp300.ICU_only           = [df[1, :Median_TD] for df in TD_results_dfs_summer_gp300]
results_df_gp300.PC_only_summer     = [df[2, :Median_TD] for df in TD_results_dfs_summer_gp300]
results_df_gp300.Combined_summer    = [df[3, :Median_TD] for df in TD_results_dfs_summer_gp300]
results_df_gp300.Improvement_summer = [df[3, :Median_TD] for df in TD_results_dfs_summer_gp300] - [df[1, :Median_TD] for df in TD_results_dfs_summer_gp300]
results_df_gp300.PC_only_winter     = [df[2, :Median_TD] for df in TD_results_dfs_winter_gp300]
results_df_gp300.Combined_winter    = [df[3, :Median_TD] for df in TD_results_dfs_winter_gp300]
results_df_gp300.Improvement_winter = [df[3, :Median_TD] for df in TD_results_dfs_winter_gp300] - [df[1, :Median_TD] for df in TD_results_dfs_winter_gp300]
results_df_gp300[4,3:5] = [-1,-1,-1]
# println(results_df_gp300)
# CSV.write("scripts/primary_care_v_icu/TDs_3cases_w_gp_300.csv", results_df_gp300) 

# % of replicates with times to detection 1 case

# Create a new dataframe with results for times to detection of 1 case
simrep_perc_df = DataFrame( Number_of_PC_mg_samples = [100,200,319,747,1000,6000,15419,28011,102792,186738])
simrep_perc_df.ICU_only           = [df[1, :Percentage_with_a_TD] for df in TD_results_dfs_summer]
simrep_perc_df.PC_only_summer     = [df[2, :Percentage_with_a_TD] for df in TD_results_dfs_summer]
simrep_perc_df.Combined_summer    = [df[3, :Percentage_with_a_TD] for df in TD_results_dfs_summer]
simrep_perc_df.Improvement_summer = [df[3, :Percentage_with_a_TD] for df in TD_results_dfs_summer] - [df[1, :Percentage_with_a_TD] for df in TD_results_dfs_summer]
simrep_perc_df.PC_only_winter     = [df[2, :Percentage_with_a_TD] for df in TD_results_dfs_winter]
simrep_perc_df.Combined_winter    = [df[3, :Percentage_with_a_TD] for df in TD_results_dfs_winter]
simrep_perc_df.Improvement_winter = [df[3, :Percentage_with_a_TD] for df in TD_results_dfs_winter] - [df[1, :Percentage_with_a_TD] for df in TD_results_dfs_winter]
simrep_perc_df[10,3:5] = [-1,-1,-1]
println(simrep_perc_df)
CSV.write("scripts/primary_care_v_icu/simrep_perc_1case_w_gp_adj.csv", simrep_perc_df) 

# % of replicates with times to detection 3 cases
# Create a new dataframe with results for times to detection of 1 case
simrep_perc_3td_df = DataFrame( Number_of_PC_mg_samples = [100,200,319,747,1000,6000,15419,28011,102792,186738])
simrep_perc_3td_df.ICU_only           = [df[4, :Percentage_with_a_TD] for df in TD_results_dfs_summer]
simrep_perc_3td_df.PC_only_summer     = [df[5, :Percentage_with_a_TD] for df in TD_results_dfs_summer]
simrep_perc_3td_df.Combined_summer    = [df[6, :Percentage_with_a_TD] for df in TD_results_dfs_summer]
simrep_perc_3td_df.Improvement_summer = [df[6, :Percentage_with_a_TD] for df in TD_results_dfs_summer] - [df[4, :Percentage_with_a_TD] for df in TD_results_dfs_summer]
simrep_perc_3td_df.PC_only_winter     = [df[5, :Percentage_with_a_TD] for df in TD_results_dfs_winter]
simrep_perc_3td_df.Combined_winter    = [df[6, :Percentage_with_a_TD] for df in TD_results_dfs_winter]
simrep_perc_3td_df.Improvement_winter = [df[6, :Percentage_with_a_TD] for df in TD_results_dfs_winter] - [df[4, :Percentage_with_a_TD] for df in TD_results_dfs_winter]
simrep_perc_3td_df[10,3:5] = [-1,-1,-1]
println(simrep_perc_3td_df)
CSV.write("scripts/primary_care_v_icu/simrep_perc_3cases_w_gp_adj.csv", simrep_perc_3td_df) 

# Plot the DataFrame with points and lines
##### 1 case ####
# Summer
results_df_ex100winter = results_df[1:9,:]
x_replace_summer = results_df_ex100winter[:,1]#[100,200,319,747,1000,6000,15419]
# plot separately so coloured by y-series
p = plot(xlabel = "Number of primary care metagenomic samples"
        , ylabel = "Median time to detection\n of 1 case (TD) in days"
        , xscale=:log10
        , ylim=(0,70)
        , legend = :bottomleft
        #, size =(1600,1200)
        , xticks = [100, 1000, 10000, 100000, 1000000]
        )
@df results_df_ex100winter scatter!(x_replace_summer
                        , :ICU_only
                        , color=:red
                        , marker=:circle
                        , markersize=4
                        , label="ICU only")
@df results_df_ex100winter scatter!(x_replace_summer
                        , :PC_only_summer
                        , color=:blue
                        , marker=:circle
                        , markersize=4
                        , label="Primary care only (summer)")
@df results_df_ex100winter scatter!(x_replace_summer
                        , :Combined_summer
                        , color=:green
                        , marker=:circle
                        , markersize=4
                        , label="Combined (summer)")
@df results_df_ex100winter plot!(x_replace_summer, :ICU_only, color=:red, linewidth=2, label="")
@df results_df_ex100winter plot!(x_replace_summer, :PC_only_summer, color=:blue, linewidth=2, label="")
@df results_df_ex100winter plot!(x_replace_summer, :Combined_summer, color=:green, linewidth=2, label="")

# Winter
x_replace_winter = results_df[:,1] #[100,200,319,747,1000,6000,28011]
# plot separately so coloured by y-series
#p = plot(xlabel = "Number of primary care metagenomic samples"
#        , ylabel = "Time to detection of 1 case (days)"
#        , xscale=:log10
#        , ylim=(0,70)
#        , legend = :topright
#        )
@df results_df scatter!(x_replace_winter
                        , :ICU_only
                        , color=:red
                        , marker=:circle
                        , markersize=4
                        , label="")
@df results_df scatter!(x_replace_winter
                        , :PC_only_winter
                        , color=:lightblue
                        , marker=:circle
                        , markersize=4
                        , label="Primary care only (winter)")
@df results_df scatter!(x_replace_winter
                        , :Combined_winter
                        , color=:lightgreen
                        , marker=:circle
                        , markersize=4
                        , label="Combined (winter)")
@df results_df plot!(x_replace_winter, :ICU_only, color=:red, linewidth=2, label="")
@df results_df plot!(x_replace_winter, :PC_only_winter, color=:lightblue, linewidth=2, label="")
@df results_df plot!(x_replace_winter, :Combined_winter, color=:lightgreen, linewidth=2, label="")

# Save to file
savefig("scripts/primary_care_v_icu/TD_vs_PC_sample_size.png")

##### 3 cases ####
# Summer
results_3td_df_ex100winter = results_3td_df[1:9,:]
x_replace_summer = results_3td_df_ex100winter[:,1]#[100,200,319,747,1000,6000,15419]
# plot separately so coloured by y-series
p = plot(xlabel = "Number of primary care metagenomic samples"
        , ylabel = "Median time to detection\n for 3 cases (3TD) in days"
        , xscale=:log10
        , ylim=(0,70)
        , legend = :bottomleft
        #, size =(1600,1200)
        , xticks = [100, 1000, 10000, 100000, 1000000]
        )
@df results_3td_df_ex100winter scatter!(x_replace_summer
                        , :ICU_only
                        , color=:red
                        , marker=:circle
                        , markersize=4
                        , label="ICU only")
@df results_3td_df_ex100winter scatter!(x_replace_summer
                        , :PC_only_summer
                        , color=:blue
                        , marker=:circle
                        , markersize=4
                        , label="Primary care only (summer)")
@df results_3td_df_ex100winter scatter!(x_replace_summer
                        , :Combined_summer
                        , color=:green
                        , marker=:circle
                        , markersize=4
                        , label="Combined (summer)")
@df results_3td_df_ex100winter plot!(x_replace_summer, :ICU_only, color=:red, linewidth=2, label="")
@df results_3td_df_ex100winter plot!(x_replace_summer, :PC_only_summer, color=:blue, linewidth=2, label="")
@df results_3td_df_ex100winter plot!(x_replace_summer, :Combined_summer, color=:green, linewidth=2, label="")

# Winter
x_replace_winter = results_3td_df[:,1] #[100,200,319,747,1000,6000,28011]
# plot separately so coloured by y-series
#p = plot(xlabel = "Number of primary care metagenomic samples"
#        , ylabel = "Time to detection of 1 case (days)"
#        , xscale=:log10
#        , ylim=(0,70)
#        , legend = :topright
#        )
@df results_3td_df scatter!(x_replace_winter
                        , :ICU_only
                        , color=:red
                        , marker=:circle
                        , markersize=4
                        , label="")
@df results_3td_df scatter!(x_replace_winter
                        , :PC_only_winter
                        , color=:lightblue
                        , marker=:circle
                        , markersize=4
                        , label="Primary care only (winter)")
@df results_3td_df scatter!(x_replace_winter
                        , :Combined_winter
                        , color=:lightgreen
                        , marker=:circle
                        , markersize=4
                        , label="Combined (winter)")
@df results_3td_df plot!(x_replace_winter, :ICU_only, color=:red, linewidth=2, label="")
@df results_3td_df plot!(x_replace_winter, :PC_only_winter, color=:lightblue, linewidth=2, label="")
@df results_3td_df plot!(x_replace_winter, :Combined_winter, color=:lightgreen, linewidth=2, label="")

# Save to file
savefig("scripts/primary_care_v_icu/3TD_vs_PC_sample_size.png")



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



#### Age profile of ICU and GP cases

# Initialise df to store ages
ages_icu_df = DataFrame( sim_rep_n = []
                        ,ICU_ages = [] )
ages_gp_df = DataFrame( sim_rep_n = []
                        ,GP_ages = [] )
                    
# Gather ages
for s in 1:10000 #length(sims_G_gp_filter)
    
    if size(sims_G_icu_filter[s],1) >0
        ## ICU
        # ICU ages
        icu_case_ages = sims_G_icu_filter[s].infectee_age
        # ICU simulation replicate number
        icu_sim_rep_temp = [ s for _ in 1:length(icu_case_ages)]
        # Add to temp df
        ages_icu_df_temp = DataFrame( sim_rep_n = icu_sim_rep_temp
                                , ICU_ages = icu_case_ages )
        # Append to output df
        ages_icu_df = append!( ages_icu_df, ages_icu_df_temp)
    end
    
    if size(sims_G_gp_filter[s],1) >0
        ## GP
        # GP ages
        gp_case_ages = sims_G_gp_filter[s].infectee_age
        # ICU simulation replicate number
        gp_sim_rep_temp = [ s for _ in 1:length(gp_case_ages)]
        # Add to temp df
        ages_gp_df_temp = DataFrame( sim_rep_n = gp_sim_rep_temp
                                , GP_ages = gp_case_ages )
        # Append to output df
        ages_gp_df = append!( ages_gp_df, ages_gp_df_temp)
    end
end

# Plot distributions
using StatsBase, Plots

# Sample data
#n1 = 1
#n2 = 30
#x1 = filter( row -> row.sim_rep_n in 1:30 , ages_icu_df)[:,2]
x1 = ages_icu_df[:,2]
x2 = ages_gp_df[:,2]

mean_icu_age = mean(x1)
mean_gp_age = mean(x2)

median_icu_age = median(x1)
median_gp_age = median(x2)

# Create histogram objects using the same bin edges
histogram([x1, x2]
    ,normalize = :pdf
    ,alpha = 0.5
    ,bins = 20 #100
    ,label = ["ICU case ages" "GP case ages"]
    ,xlabel = "Age in years"
    , ylabel = "Density"
    , colors = [:blue,:red]
    , linewidth = 2)


## Check how many simulations contain results
n_sims_w_results_icu_gp( sims_file_icu = "covidlike-1.1.1-sims_filtered_G_icu_combined_nrep53000.jld2" 
                , sims_file_gp = "covidlike-1.1.1-sims_filtered_G_gp_combined_nrep53000.jld2") 
# ICU: 90% and GP: 90% 

n_sims_w_results_icu_gp( sims_file_icu = "covidlike-1.3.1-sims_filtered_G_icu_combined_nrep84000.jld2" 
                , sims_file_gp = "covidlike-1.3.1-sims_filtered_G_gp_combined_nrep84000.jld2") 
# Combination of HPC run 851401+854085 
# ICU: 35% and GP: 34% (because R was too low at ~1)

n_sims_w_results_icu_gp( sims_file_icu = "covidlike-1.3.1-sims_filtered_G_icu_combined_nrep64000_955898.jld2" 
                , sims_file_gp = "covidlike-1.3.1-sims_filtered_G_gp_combined_nrep64000_955898.jld2") 
# ICU: 80% and GP: 80%
