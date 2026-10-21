#Script for the comparing the sampling performance in three different 
#simulated pathogen outbreaks.
#
#This file contains:
#(0) Define root directory
#(1) Load sampling and detection data for different pathogens:
#  (1.1) COVID-like outbreak
#  (1.2) Influenza-like mild outbreak
#  (1.3) Influenza-like severe outbreak
#(2) Plot comparing the time to detection performance in different healthcare settings across different pathogens
#(3) Statistics for article text
#

#-----------------------------------
# Activate environment for analysis
#-----------------------------------
# Unless otherwise stated, run on Laptop - 11th Gen Intel(R) Core(TM) i7-1185G7 @ 3.00GHz (1.80 GHz) with 32.0 GB RAM
# Using Julia v1.12.3 in VS Code 1.141.0

using Pkg
Pkg.activate(@__DIR__) # ...NBPMscape\\paper\\analysis\\2_sampling_analysis\\scripts

# Check packages in use
# Using NBPMscape version from commit "2762757a687c64146016e6bc742d258e20525cb8" as at 17:00 16 Oct 2026
Pkg.status()

# Check the environment location
Base.active_project()
#"...NBPMscape\\paper\\analysis\\2_sampling_analysis\\scripts\\Project.toml"

#-----------------------------------
# Load required packages
#-----------------------------------
Pkg.instantiate()

using CSV
using DataFrames
using JLD2
using NBPMscape # Should be version from commit "2762757a687c64146016e6bc742d258e20525cb8" as at 17:00 16 Oct 2026
using Plots
using Statistics

############################################################################################
# (0) Define root directory
############################################################################################

# Define root directory
#root_dir = root_folder() # Function defined in NBPMscape/src/misc_functions
dir = @__DIR__
while !isdir(joinpath(dir, "NBPMscape"))
    parent = dirname(dir)
    parent == dir && error("NBPmscape not found")
    dir = parent
end
root_dir = dir

############################################################################################
# (1) Load sampling and detection data for different pathogens:
#   (1.1) COVID-like outbreak
############################################################################################

# Load COVID-19-like dictionaries containing time to detection (TD) data
# (Dictionaries created in 'NBPMscape/paper/analysis/2_sampling_analysis/scripts/2_healthcare_sampling_setting.jl')
input_dict_covid = joinpath( root_dir, "NBPMscape/paper/analysis/2_sampling_analysis/covid_like/healthcare_setting")
covid_gp_tds_df_dict     = load( joinpath( input_dict_covid, "covidlike-1.6.1_nrep1000_4300448_gp_tds_df_dict_scenario_numbers.jld2"),         "gp_tds_df_dict_scenario_numbers" )
covid_icu_8_tds_df_dict  = load( joinpath( input_dict_covid, "covidlike-1.6.1_nrep1000_4300448_icu_8_tds_df_dict_scenario_numbers.jld2"),   "icu_8_tds_df_dict_scenario_numbers" )
covid_icu_26_tds_df_dict = load( joinpath( input_dict_covid, "covidlike-1.6.1_nrep1000_4300448_icu_26_tds_df_dict_scenario_numbers.jld2"), "icu_26_tds_df_dict_scenario_numbers" )

############################################################################################
# (1) Load sampling and detection data for different pathogens:
#   (1.2) Influenza-like mild outbreak
############################################################################################

# Load MILD Influenza-like dictionaries containing TD data
# See 'NBPMscape/paper/analysis/2_sampling_analysis/scripts/4_influenza_like_outbreaks.jl' for how they were produced
folder_mild_flu = joinpath( root_dir, "NBPMscape/paper/analysis/2_sampling_analysis/influenza_like_mild" )

mild_flu_icu_8_tds_df_dict  = Dict{Int, DataFrame}()
mild_flu_icu_26_tds_df_dict = Dict{Int, DataFrame}()
mild_flu_gp_tds_df_dict     = Dict{Int, DataFrame}()

# See 'NBPMscape/config/2026_02_18/config_file_record.xlsx' for summary descriptions of scenarios
#(Column alignment indicates comparable scenarios)
mild_flu_scenario_numbers_gp     = [ 1, 13, 25, 37, 49, 61]
mild_flu_scenario_numbers_icu_8  = [ 1, 13, 25, 37, 49, 61]
mild_flu_scenario_numbers_icu_26 = [73, 74, 75, 76, 77, 78]

# Loop through scenarios, read in files
# Primary care (GPs modelled on Oxford-RCGP RSC)
for i in 1:length(mild_flu_scenario_numbers_gp) 
    k = mild_flu_scenario_numbers_gp[i]
    # Construct the file name, pad with leading zeros if needed
    file_gp = joinpath(folder_mild_flu, "gp_td_small/small_gp_tds_$(string(k)).csv")
    # Read CSV into a DataFrame and store in the dictionary with the scenario number (numeric suffix) as key
    mild_flu_gp_tds_df_dict[k] = CSV.read(file_gp, DataFrame)
end

# ICU 8 sites (modelled on mSCAPE pilot study)
for i in 1:length(mild_flu_scenario_numbers_icu_8) 
    k = mild_flu_scenario_numbers_icu_8[i]
    # Construct the file name, pad with leading zeros if needed
    file_icu_8 = joinpath(folder_mild_flu, "icu_td_small/small_icu_tds_$(string(k)).csv")
    # Read CSV into a DataFrame and store in the dictionary with the scenario number (numeric suffix) as key
    mild_flu_icu_8_tds_df_dict[k] = CSV.read(file_icu_8, DataFrame)
end

# ICU 26 sites (modelled on mSCAPE longlist)
for i in 1:length(mild_flu_scenario_numbers_icu_26) 
    k = mild_flu_scenario_numbers_icu_26[i]
    # Construct the file name, pad with leading zeros if needed
    file_icu_26 = joinpath(folder_mild_flu, "icu_td_small/small_icu_tds_$(string(k)).csv")
    # Read CSV into a DataFrame and store in the dictionary with the scenario number (numeric suffix) as key
    mild_flu_icu_26_tds_df_dict[k] = CSV.read(file_icu_26, DataFrame)
end


############################################################################################
# (1) Load sampling and detection data for different pathogens:
#   (1.3) Influenza-like severe outbreak
############################################################################################

# Load SEVERE Influenza-like dictionaries containing TD data
# See 'NBPMscape/paper/analysis/2_sampling_analysis/scripts/4_influenza_like_outbreaks.jl' for how they were produced
folder_severe_flu = joinpath( root_dir, "NBPMscape/paper/analysis/2_sampling_analysis/influenza_like_severe" )

severe_flu_icu_8_tds_df_dict  = Dict{Int, DataFrame}()
severe_flu_icu_26_tds_df_dict = Dict{Int, DataFrame}()
severe_flu_gp_tds_df_dict     = Dict{Int, DataFrame}()

# See 'NBPMscape/config/2026_02_18/config_file_record.xlsx' for summary descriptions of scenarios
#(Column alignment indicates comparable scenarios)
severe_flu_scenario_numbers_gp     = [ 1, 13, 25, 37, 49, 61]
severe_flu_scenario_numbers_icu_8  = [ 1, 13, 25, 37, 49, 61]
severe_flu_scenario_numbers_icu_26 = [73, 74, 75, 76, 77, 78]

# Loop through scenarios, read in files
# Primary care (GPs modelled on Oxford-RCGP RSC)
for i in 1:length(severe_flu_scenario_numbers_gp) 
    k = severe_flu_scenario_numbers_gp[i]
    # Construct the file name, pad with leading zeros if needed
    file_gp = joinpath(folder_severe_flu, "gp_td_small/small_gp_tds_$(string(k)).csv")
    # Read CSV into a DataFrame and store in the dictionary with the scenario number (numeric suffix) as key
    severe_flu_gp_tds_df_dict[k] = CSV.read(file_gp, DataFrame)          
end

# ICU 8 sites (modelled on mSCAPE pilot study)
for i in 1:length(severe_flu_scenario_numbers_icu_8) 
    k = severe_flu_scenario_numbers_icu_8[i]
    # Construct the file name, pad with leading zeros if needed
    file_icu_8 = joinpath(folder_severe_flu, "icu_td_small/small_icu_tds_$(string(k)).csv")
    # Read CSV into a DataFrame and store in the dictionary with the scenario number (numeric suffix) as key
    severe_flu_icu_8_tds_df_dict[k] = CSV.read(file_icu_8, DataFrame)
end

# ICU 26 sites (modelled on mSCAPE longlist)
for i in 1:length(severe_flu_scenario_numbers_icu_26) 
    k = severe_flu_scenario_numbers_icu_26[i]
    # Construct the file name, pad with leading zeros if needed
    file_icu_26 = joinpath(folder_severe_flu, "icu_td_small/small_icu_tds_$(string(k)).csv")
    # Read CSV into a DataFrame and store in the dictionary with the scenario number (numeric suffix) as key
    severe_flu_icu_26_tds_df_dict[k] = CSV.read(file_icu_26, DataFrame)
end


####################################################################################################################################################
# (2) Plot comparing the time to detection performance in different healthcare settings across different pathogens
####################################################################################################################################################

# Plots for figure comparing time to detection across different pathogens:
# COVID-like, mild Influenza-like, severe Influenza-like


### Plot
# (dictionary keys follow the scenario naming convention in 'NBPMscape\config\HARISS\config file record.xlsx')
# Build individual plots
x_maximum = 100

# Define colours
pc = :red
cc8 = :lightblue
cc26 = :darkblue

# COVID-like TD
p_covid_gp_td     = kernel_box_jitter_plot(;x = covid_gp_tds_df_dict[1][:,:GP_TD],       plot_color = pc,   samp_strategy="Primary care",             n_samples=300,n_sites="",n_blocks="",season="Winter", x_max = x_maximum, x_ticks=false, x_label = "")
p_covid_icu_8_td  = kernel_box_jitter_plot(;x = covid_icu_8_tds_df_dict[1][:,:ICU_TD],   plot_color = cc8,  samp_strategy="Critical care (8 sites)",  n_samples=300,n_sites="",n_blocks="",season="Winter", x_max = x_maximum, x_ticks=false, x_label = "")
p_covid_icu_26_td = kernel_box_jitter_plot(;x = covid_icu_26_tds_df_dict[73][:,:ICU_TD], plot_color = cc26, samp_strategy="Critical care (26 sites)", n_samples=300,n_sites="",n_blocks="",season="Winter", x_max = x_maximum, x_ticks=false, x_label = "")

# COVID-like 3TD
# Note: covid_gp_tds_df_dict[1][:,:GP_3TD] contains a single Inf value (i.e. there was no 3TD in 1 out of 100 simulations)
p_covid_gp_3td     = kernel_box_jitter_plot(;x = filter(isfinite, covid_gp_tds_df_dict[1][:,:GP_3TD]), plot_color = pc,   samp_strategy="Primary care",             n_samples=300,n_sites="",n_blocks="",season="Winter", x_max = x_maximum, x_label = "Time since first imported infection (days)")
p_covid_icu_8_3td  = kernel_box_jitter_plot(;x = covid_icu_8_tds_df_dict[1][:,:ICU_3TD],               plot_color = cc8,  samp_strategy="Critical care (8 sites)",  n_samples=300,n_sites="",n_blocks="",season="Winter", x_max = x_maximum, x_ticks=false, x_label = "")
p_covid_icu_26_3td = kernel_box_jitter_plot(;x = covid_icu_26_tds_df_dict[73][:,:ICU_3TD],             plot_color = cc26, samp_strategy="Critical care (26 sites)", n_samples=300,n_sites="",n_blocks="",season="Winter", x_max = x_maximum, x_ticks=false, x_label = "")

# Mild influenza TD
p_mild_flu_gp_td     = kernel_box_jitter_plot(;x = mild_flu_gp_tds_df_dict[1][:,:GP_TD],       plot_color = pc,   samp_strategy="Primary care",             n_samples=300,n_sites="",n_blocks="",season="Winter", x_max = x_maximum, x_ticks=false, x_label = "")
p_mild_flu_icu_8_td  = kernel_box_jitter_plot(;x = mild_flu_icu_8_tds_df_dict[1][:,:ICU_TD],   plot_color = cc8,  samp_strategy="Critical care (8 sites)",  n_samples=300,n_sites="",n_blocks="",season="Winter", x_max = x_maximum, x_ticks=false, x_label = "")
p_mild_flu_icu_26_td = kernel_box_jitter_plot(;x = mild_flu_icu_26_tds_df_dict[73][:,:ICU_TD], plot_color = cc26, samp_strategy="Critical care (26 sites)", n_samples=300,n_sites="",n_blocks="",season="Winter", x_max = x_maximum, x_ticks=false, x_label = "")

# Mild influenza 3TD
p_mild_flu_gp_3td     = kernel_box_jitter_plot(;x = mild_flu_gp_tds_df_dict[1][:,:GP_3TD],       plot_color = pc,   samp_strategy="Primary care",            n_samples=300,n_sites="",n_blocks="",season="Winter", x_max = x_maximum, x_label = "Time since first imported infection (days)")
p_mild_flu_icu_8_3td  = kernel_box_jitter_plot(;x = mild_flu_icu_8_tds_df_dict[1][:,:ICU_3TD],   plot_color = cc8,  samp_strategy="Critical care (8 sites)", n_samples=300,n_sites="",n_blocks="",season="Winter", x_max = x_maximum, x_ticks=false, x_label = "")
p_mild_flu_icu_26_3td = kernel_box_jitter_plot(;x = mild_flu_icu_26_tds_df_dict[73][:,:ICU_3TD], plot_color = cc26, samp_strategy="Critical care (26 sites)",n_samples=300,n_sites="",n_blocks="",season="Winter", x_max = x_maximum, x_ticks=false, x_label = "")

# Severe influenza TD
p_severe_flu_gp_td     = kernel_box_jitter_plot(;x = severe_flu_gp_tds_df_dict[1][:,:GP_TD],       plot_color = pc,   samp_strategy="Primary care",             n_samples=300,n_sites="",n_blocks="",season="Winter", x_max = x_maximum, x_ticks=false, x_label = "")
p_severe_flu_icu_8_td  = kernel_box_jitter_plot(;x = severe_flu_icu_8_tds_df_dict[1][:,:ICU_TD],   plot_color = cc8,  samp_strategy="Critical care (8 sites)",  n_samples=300,n_sites="",n_blocks="",season="Winter", x_max = x_maximum, x_ticks=false, x_label = "")
p_severe_flu_icu_26_td = kernel_box_jitter_plot(;x = severe_flu_icu_26_tds_df_dict[73][:,:ICU_TD], plot_color = cc26, samp_strategy="Critical care (26 sites)", n_samples=300,n_sites="",n_blocks="",season="Winter", x_max = x_maximum, x_ticks=false, x_label = "")

# Severe influenza 3TD
p_severe_flu_gp_3td     = kernel_box_jitter_plot(;x = severe_flu_gp_tds_df_dict[1][:,:GP_3TD],       plot_color = pc,        samp_strategy="Primary care",             n_samples=300,n_sites="",n_blocks="",season="Winter", x_max = x_maximum, x_label = "Time since first imported infection (days)")
# Note: severe_flu_icu_8_tds_df_dict[1][:,:ICU_3TD] contains a single Inf value (i.e. there was no 3TD in 1 out of 100 simulations)
p_severe_flu_icu_8_3td  = kernel_box_jitter_plot(;x = filter(isfinite,severe_flu_icu_8_tds_df_dict[1][:,:ICU_3TD]),   plot_color = cc8,      samp_strategy="Critical care (8 sites)",  n_samples=300,n_sites="",n_blocks="",season="Winter", x_max = x_maximum, x_ticks=false, x_label = "")
p_severe_flu_icu_26_3td = kernel_box_jitter_plot(;x = severe_flu_icu_26_tds_df_dict[73][:,:ICU_3TD], plot_color = cc26, samp_strategy="Critical care (26 sites)", n_samples=300,n_sites="",n_blocks="",season="Winter", x_max = x_maximum, x_ticks=false, x_label = "")


# Ordered by care setting - primary at the top and tertiary at the bottom
Plots.plot(  p_covid_gp_td,      p_mild_flu_gp_td,      p_severe_flu_gp_td
           , p_covid_icu_8_td,   p_mild_flu_icu_8_td,   p_severe_flu_icu_8_td
           , p_covid_icu_26_td,  p_mild_flu_icu_26_td,  p_severe_flu_icu_26_td
           , p_covid_gp_3td,     p_mild_flu_gp_3td,     p_severe_flu_gp_3td
           , p_covid_icu_8_3td,  p_mild_flu_icu_8_3td,  p_severe_flu_icu_8_3td
           , p_covid_icu_26_3td, p_mild_flu_icu_26_3td, p_severe_flu_icu_26_3td
           , layout = (6, 3)
        #, size = (1200,1800)) #default size is 600,400
        , size = (1800,1800)) #default size is 600,400

# Ordered by time to detection performance from shortest at the top to longest at the bottom
Plots.plot(  p_covid_icu_26_td,  p_mild_flu_icu_26_td,  p_severe_flu_icu_26_td
           , p_covid_icu_8_td,   p_mild_flu_icu_8_td,   p_severe_flu_icu_8_td
           , p_covid_gp_td,      p_mild_flu_gp_td,      p_severe_flu_gp_td
           , p_covid_icu_26_3td, p_mild_flu_icu_26_3td, p_severe_flu_icu_26_3td
           , p_covid_icu_8_3td,  p_mild_flu_icu_8_3td,  p_severe_flu_icu_8_3td
           , p_covid_gp_3td,     p_mild_flu_gp_3td,     p_severe_flu_gp_3td
           , layout = (6, 3)
        #, size = (1200,1800)) #default size is 600,400
        , size = (2100,2100)) #default size is 600,400

Plots.savefig( joinpath( root_dir, "NBPMscape/paper/analysis/2_sampling_analysis/jitter_TD_3TD_COVID_mild_severe_flu_GP_ICU.png" ) )
Plots.savefig( joinpath( root_dir, "NBPMscape/paper/results/figures_main/raw/jitter_TD_3TD_COVID_mild_severe_flu_GP_ICU.png" ) )



####################################################################################################
# (3) Statistics for article text
####################################################################################################

# Calculating ratios of TDs (primary care / acute secondary care) for different pathogens to be used in 
# section '2.4 Optimal pathogen agnostic sampling strategy may differ from optimal sampling strategies for specific pathogens'

# Focus on winter, 300 samples per week (scenario 1 - see 'NBPMscape/config/2026_02_18/config file record.xlsx')

# Initialise df to store results
stats_df = DataFrame( pathogen = ["covid-like", "mild_influenzs_like","severe_influenzs_like"]
                    , primary_care_td = [-1.0,-1.0,-1.0]
                    , critical_care_26_td = [-1.0,-1.0,-1.0]
                    , critical_care_8_td = [-1.0,-1.0,-1.0]
                    , pc_icu26_td_ratio = [-1.0,-1.0,-1.0] )

stats_df[1,"primary_care_td"] = median(covid_gp_tds_df_dict[1][:,:GP_TD])
stats_df[2,"primary_care_td"] = median(mild_flu_gp_tds_df_dict[1][:,:GP_TD])
stats_df[3,"primary_care_td"] = median(severe_flu_gp_tds_df_dict[1][:,:GP_TD])

stats_df[1,"critical_care_26_td"] = median(covid_icu_26_tds_df_dict[73][:,:ICU_TD])
stats_df[2,"critical_care_26_td"] = median(mild_flu_icu_26_tds_df_dict[73][:,:ICU_TD])
stats_df[3,"critical_care_26_td"] = median(severe_flu_icu_26_tds_df_dict[73][:,:ICU_TD])

stats_df[1,"critical_care_8_td"] = median(covid_icu_8_tds_df_dict[1][:,:ICU_TD])
stats_df[2,"critical_care_8_td"] = median(mild_flu_icu_8_tds_df_dict[1][:,:ICU_TD])
stats_df[3,"critical_care_8_td"] = median(severe_flu_icu_8_tds_df_dict[1][:,:ICU_TD])

stats_df[1,"pc_icu26_td_ratio"] = stats_df[1,"primary_care_td"] / stats_df[1,"critical_care_26_td"]
stats_df[2,"pc_icu26_td_ratio"] = stats_df[2,"primary_care_td"] / stats_df[2,"critical_care_26_td"]
stats_df[3,"pc_icu26_td_ratio"] = stats_df[3,"primary_care_td"] / stats_df[3,"critical_care_26_td"]

println(stats_df)

#3×5 DataFrame
# Row │ pathogen               primary_care_td  critical_care_26_td  critical_care_8_td  pc_icu26_td_ratio 
#     │ String                 Float64          Float64              Float64             Float64           
#─────┼────────────────────────────────────────────────────────────────────────────────────────────────────
#   1 │ covid-like                     50.4817              44.6557             48.5729           1.13047
#   2 │ mild_influenzs_like            45.1746              56.3412             59.2789           0.801804
#   3 │ severe_influenzs_like          39.7573              36.7628             40.5825           1.08146

# Helper functions to:
# Compute median, lower IQR and upper IQR
iqr_stats(dict, key, col) = [median(dict[key][:,col]), quantile(dict[key][:,col], 0.25), quantile(dict[key][:,col], 0.75)]
# Build column with IQR stats for GP sampling, ICU (8 sites) and ICU (26 sites)
build_col(gp, icu8, icu26) = vcat(
    iqr_stats(gp,     1,  :GP_TD),
    iqr_stats(icu8,   1, :ICU_TD),
    iqr_stats(icu26, 73, :ICU_TD),
)
# Define labels for statistics to add to the df
metrics = [
    "primary_care_td",     "primary_care_td_iqr_lower",     "primary_care_td_iqr_upper",
    "critical_care_8_td",  "critical_care_8_td_iqr_lower",  "critical_care_8_td_iqr_upper",
    "critical_care_26_td", "critical_care_26_td_iqr_lower", "critical_care_26_td_iqr_upper",
]

stats_df = DataFrame(
    metric                = metrics,
    covid_like            = build_col(covid_gp_tds_df_dict,       covid_icu_8_tds_df_dict,       covid_icu_26_tds_df_dict),
    mild_influenza_like   = build_col(mild_flu_gp_tds_df_dict,    mild_flu_icu_8_tds_df_dict,    mild_flu_icu_26_tds_df_dict),
    severe_influenza_like = build_col(severe_flu_gp_tds_df_dict,  severe_flu_icu_8_tds_df_dict,  severe_flu_icu_26_tds_df_dict),
)

#9×4 DataFrame
# Row │ metric                         covid_like  mild_influenza_like  severe_influenza_like 
#     │ String                         Float64     Float64              Float64               
#─────┼───────────────────────────────────────────────────────────────────────────────────────
#   1 │ primary_care_td                   50.4817              45.1746                39.7573
#   2 │ primary_care_td_iqr_lower         44.2949              36.8502                35.3593
#   3 │ primary_care_td_iqr_upper         55.9022              51.0068                43.9126
#   4 │ critical_care_8_td                48.5729              59.2789                40.5825
#   5 │ critical_care_8_td_iqr_lower      42.3974              54.844                 34.6103
#   6 │ critical_care_8_td_iqr_upper      54.3711              65.3266                45.0778
#   7 │ critical_care_26_td               44.6557              56.3412                36.7628
#   8 │ critical_care_26_td_iqr_lower     39.4899              51.2638                32.194
#   9 │ critical_care_26_td_iqr_upper     50.596               62.9721                42.7772