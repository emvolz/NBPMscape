#=

Script for the investigation of the optimal number of sampling sites
for surveillance performance, as determined by times to detection (TD)
and the variance of the time of detection.

This file contains:
(0) Define root directory
(1) Load data computed in '1_geographic_optimisation.jl'
(2) Define helper functions
(3) Plot TD and variance of TD against number of sampling sites.
    (only HCGS strategy selected sites)
    Separate plots for different:
    - number of hierarchical clustering blocks
    - number of samples per week
    - season 
    (used in supplementary material)
(4) Plot TD and variance of TD separately against number of sampling sites.
   Only one sampling scenario. Used in main text.
   (4.1) TD
   (4.2) Variance of TD

=#

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
#using CMPlot
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
# (1) Load data computed in '1_geographic_optimisation.jl'
############################################################################################

# Load results created using "NBPMscape/paper/analysis/2_sampling_analysis/scripts/1c_geographic_optimisation_part_2.jl"
results_df = CSV.read( joinpath( root_dir, "NBPMscape/paper/analysis/2_sampling_analysis/covid_like/geo_optimisation"
                    ,"icu_td_results_df.csv" )
                    , DataFrame )

icu_tds_df_dict = load( joinpath( root_dir 
                                , "NBPMscape/paper/analysis/2_sampling_analysis/covid_like/geo_optimisation"
                                , "covidlike-1.6.1_nrep1000_4300448_icu_tds_df_dict.jld2")
                        , "icu_tds_df_dict")


############################################################################################
# (2) Define helper functions                        
############################################################################################

# Helper function to compute the median, variance, and 95% estimate using bootstrap
function median_var_boot(; vec::Vector, n_boot::Int64=1000, alpha::Float64=0.05)
    med = median(vec)
    med_boot_samples = [median(rand(vec, length(vec))) for _ in 1:n_boot]
    med_lower = quantile(med_boot_samples, alpha/2)
    med_upper = quantile(med_boot_samples, 1 - alpha/2)

    vec_var = var(vec)
    vec_var_boot_samples = [var(rand(vec, length(vec))) for _ in 1:n_boot]
    vec_var_lower = quantile(vec_var_boot_samples, alpha/2)
    vec_var_upper = quantile(vec_var_boot_samples, 1 - alpha/2)

    return (median = med, median_lower = med_lower, median_upper = med_upper, variance = vec_var, variance_lower = vec_var_lower, variance_upper = vec_var_upper)
end


# Function to generate data for plotting from TD values from multiple simulations
function td_plot_data(; number_of_sites_vec::Vector{Int64} = [8,10,20,26,30,40,50,60]
                      , plot_scenario::String
                      )
    
    # Initialise vectors
    td_median = fill(-1.0,length(number_of_sites_vec))
    td_median_error_lower = copy(td_median)
    td_median_error_upper = copy(td_median)
    td_var = copy(td_median)
    td_var_error_lower = copy(td_median)
    td_var_error_upper = copy(td_median)

    for i in 1:length(number_of_sites_vec)
        td_vec = icu_tds_df_dict[plot_scenario * "$(number_of_sites_vec[i])" * "sites"][:,:ICU_TD]
        
        median_var_boot_results = median_var_boot(; vec = td_vec, n_boot = 1000, alpha = 0.05)
        
        td_median[i] = median_var_boot_results.median
        td_median_error_lower[i] = median_var_boot_results.median - median_var_boot_results.median_lower
        td_median_error_upper[i] = median_var_boot_results.median_upper - median_var_boot_results.median

        td_var[i] = median_var_boot_results.variance
        td_var_error_lower[i] = median_var_boot_results.variance - median_var_boot_results.variance_lower
        td_var_error_upper[i] = median_var_boot_results.variance_upper - median_var_boot_results.variance
    end

    return (td_median = td_median, td_median_error_lower = td_median_error_lower, td_median_error_upper = td_median_error_upper
            , td_var = td_var, td_var_error_lower = td_var_error_lower, td_var_error_upper = td_var_error_upper)
end

function plot_med_var(; x::Vector{Int64}
                      , y1::Vector{Float64}, y1_lower_errors::Vector{Float64}, y1_upper_errors::Vector{Float64}
                      , y2::Vector{Float64}, y2_lower_errors::Vector{Float64}, y2_upper_errors::Vector{Float64}
                      , title::String )

    # Median TD on the left y-axis
    p = Plots.plot(#scatter(
        [x ; x ]
        , [y1; y1]
        , yerror = (y1_lower_errors, y1_upper_errors)
        , seriestype = :line
        , color = :blue
        , linecolor = :blue
        , errorbarcolor = :blue 
        , marker = (:circle, 6, :blue, 0.5)
        #, label  = "Median time to detection (LHS)"
        #, ylabel = "Median time to detect 1st case\n (days since first imported infection)"
        #, xlabel = "Number of ICU sampling sites"
        , label  = ""
        , ylabel = ""
        , xlabel = ""
        , title  = ""
        , legend = :topleft
        , alpha = 0.5
        , ylimit = [35,50]
    )
    # Variance in TD on the right y-axis
    Plots.plot!(#scatter(
        twinx(p)
        , [x ; x ]
        , [y2; y2]
        , yerror = (y2_lower_errors, y2_upper_errors)
        , seriestype = :line
        , color = :red
        , linecolor = :red
        , errorbarcolor = :red
        , marker = (:square, 6, :red, 0.5)
        #, label = "Variance (RHS)"
        #, ylabel = "Variance in time to detection"
        #,xlabel = "Number of ICU sampling sites"
        , title  = title # eg "Winter, 300 samples, 10 blocks"
        , label = ""
        , ylabel = ""
        , legend = :topright
        , alpha = 0.5
        , ylimit = [65,85]
    )

    display(p)
    return (p)
end


# Helper function: builds scenario prefix, runs pipeline, returns plot
function make_combo_plot(number_of_sites_vec, season_code, n_samples, n_blocks)
    prefix = "icu_tds_seasonal_ari_$(season_code)_$(n_samples)_sampling_sites_$(n_blocks)blocks_"
    plot_data_results = td_plot_data(; number_of_sites_vec = number_of_sites_vec, plot_scenario = prefix)
    return plot_med_var(;
        x                 = number_of_sites_vec,
        y1                = plot_data_results.td_median,
        y1_lower_errors   = plot_data_results.td_median_error_lower,
        y1_upper_errors   = plot_data_results.td_median_error_upper,
        y2                = plot_data_results.td_var,
        y2_lower_errors   = plot_data_results.td_var_error_lower,
        y2_upper_errors   = plot_data_results.td_var_error_upper,
        title             = ""
    )
end

############################################################################################
# (3) Plot TD and variance of TD against number of sampling sites.
#     (only HCGS strategy selected sites)
#     Separate plots for different:
#     - number of hierarchical clustering blocks
#     - number of samples per week
#     - season 
#     (used in supplementary material)
############################################################################################

# Numbers of sites
number_of_sites_vec = [8,10,20,26,30,40,50,60]

# Numbers of sites
number_of_sites_vec = [8, 10, 20, 26, 30, 40, 50, 60]

# Season code: Winter = 1440, Summer = 793
season_codes = Dict("w" => 1440, "s" => 793)
sample_counts = [300, 500, 1000]
block_counts  = [10, 20, 30]

# Generate all plots and store in a dictionary, e.g. plots["w_300s_10b"]
plots = Dict(
    "$(season)_$(n_samples)s_$(n_blocks)b" =>
        make_combo_plot(number_of_sites_vec, code, n_samples, n_blocks)
    for (season, code) in season_codes
    for n_samples in sample_counts
    for n_blocks  in block_counts
)

# Plot layout is 
#                          10 block, 20, 30, 10 blocks, 20, 30
# winter      300 samples
#             500 samples
#            1000 samples
# summer      300 samples
#             500 samples
#            1000 samples

# Winter
p_winter = Plots.plot(  plots["w_300s_10b"],  plots["w_300s_20b"],  plots["w_300s_30b"]
                      , plots["w_500s_10b"],  plots["w_500s_20b"],  plots["w_500s_30b"]
                      , plots["w_1000s_10b"], plots["w_1000s_20b"], plots["w_1000s_30b"]
        , layout = (3, 3)
        , size = (1200,600)) #default size is 600,400

# Summer
p_summer = Plots.plot(  plots["s_300s_10b"],  plots["s_300s_20b"],  plots["s_300s_30b"]
                      , plots["s_500s_10b"],  plots["s_500s_20b"],  plots["s_500s_30b"]
                      , plots["s_1000s_10b"], plots["s_1000s_20b"], plots["s_1000s_30b"]
        , layout = (3, 3)
        , size = (1200,600)) #default size is 600,400

Plots.plot( p_winter, p_summer, layout = (2,1), size = (1000, 1000))
Plots.savefig( joinpath( root_dir
                        , "NBPMscape/paper/analysis/2_sampling_analysis/covid_like/optimal_number_of_sites"
                        , "ICU_TD_median_and_variance_plot_winter_top_summer_bottom_rows_samples_cols_blocks.png" ) )
Plots.savefig( joinpath( root_dir
                        , "NBPMscape/paper/results/figures_supplementary/raw"
                        , "ICU_TD_median_and_variance_plot_winter_top_summer_bottom_rows_samples_cols_blocks.png" ) )                        

# Blank plot with just the legend so can add in Inkscape for figure
p = Plots.plot( [0 ; 0 ], [0; 0]
        , seriestype = :line
        , color = :blue
        , linecolor = :blue
        , label  = "Median time to detection (LHS)"
        , marker = (:circle, 6, :blue, 0.5)
        , legend = :topleft
        , alpha = 0.5
    )
Plots.plot!(twinx(p), [0 ; 0 ], [0; 0]
        , seriestype = :line
        , color = :red
        , linecolor = :red
        , label = "Variance (RHS)"
        , marker = (:square, 6, :red, 0.5)
        , ylabel = "Variance in time to detection"
        , legend = :topright
        , alpha = 0.5
            )
Plots.savefig( joinpath( root_dir
                        , "NBPMscape/paper/analysis/2_sampling_analysis/covid_like/optimal_number_of_sites"
                        , "ICU_TD_median_and_variance_plot_winter_top_summer_bottom_rows_samples_cols_blocks_legend.png" ) )
Plots.savefig( joinpath( root_dir
                        , "NBPMscape/paper/results/figures_supplementary/raw"
                        , "ICU_TD_median_and_variance_plot_winter_top_summer_bottom_rows_samples_cols_blocks_legend.png" ) )


#############################################################################################
# (4) Plot TD and variance of TD separately against number of sampling sites.
#    Only one sampling scenario. Used in main text.
#    (4.1) TD
#############################################################################################

### Only for winter and 20 blocks, but multiple samples per week

# Define plot data
plot_data_w_300s_20b  = td_plot_data(; number_of_sites_vec = number_of_sites_vec, plot_scenario = "icu_tds_seasonal_ari_1440_300_sampling_sites_20blocks_")
plot_data_w_500s_20b  = td_plot_data(; number_of_sites_vec = number_of_sites_vec, plot_scenario = "icu_tds_seasonal_ari_1440_500_sampling_sites_20blocks_")
plot_data_w_1000s_20b = td_plot_data(; number_of_sites_vec = number_of_sites_vec, plot_scenario = "icu_tds_seasonal_ari_1440_1000_sampling_sites_20blocks_")

## Time to detection (TD) vs number of sampling sites
# Plot median TDs
x = number_of_sites_vec

y1 = plot_data_w_300s_20b.td_median
y1_lower_errors = plot_data_w_300s_20b.td_median_error_lower
y1_upper_errors = plot_data_w_300s_20b.td_median_error_upper

y2 = plot_data_w_500s_20b.td_median
y2_lower_errors = plot_data_w_500s_20b.td_median_error_lower
y2_upper_errors = plot_data_w_500s_20b.td_median_error_upper

y3 = plot_data_w_1000s_20b.td_median
y3_lower_errors = plot_data_w_1000s_20b.td_median_error_lower
y3_upper_errors = plot_data_w_1000s_20b.td_median_error_upper

color1 = :skyblue1
color2 = :blue
color3 = :navyblue

p_median = Plots.plot( [x ; x ], [y1; y1]
              , yerror = (y1_lower_errors, y1_upper_errors)
            , seriestype = :line, color = color1, linecolor = color1, errorbarcolor = color1
            , marker = (:circle, 6, color1, 0.5)
            , label  = "300 samples per week"
            , ylabel = "Median time to detect 1st case\n (days since first imported infection)"
            , xlabel = "Number of ICU sampling sites"
            , legend = :topleft
            , background_color_legend = nothing
            , foreground_color_legend = nothing
            , alpha = 0.5
            , ylimit = [35,50]
                )

Plots.plot!( [x ; x ], [y2; y2]
              , yerror = (y2_lower_errors, y2_upper_errors)
            , seriestype = :line, color = color2, linecolor = color2, errorbarcolor = color2
            , marker = (:square, 6, color2, 0.5)
            , label  = "500 samples per week"
            , alpha = 0.5
            )

Plots.plot!( [x ; x ], [y3; y3]
            , yerror = (y3_lower_errors, y3_upper_errors)
            , seriestype = :line, color = color3, linecolor = color3, errorbarcolor = color3
            , marker = (:diamond, 6, color3, 0.5)
            , label  = "1000 samples per week"
            , alpha = 0.5
            )

Plots.savefig( joinpath( root_dir
                        , "NBPMscape/paper/analysis/2_sampling_analysis/covid_like/optimal_number_of_sites"
                        , "ICU_TD_median_winter_HCGS_20_blocks.png" ) )
                        
Plots.savefig( joinpath( root_dir
                        , "NBPMscape/paper/results/figures_main/raw"
                        , "ICU_TD_median_winter_HCGS_20_blocks.png" ) )

#############################################################################################
# (4) Plot TD and variance of TD separately against number of sampling sites.
#    Only one sampling scenario. Used in main text.
#    (4.2) Variance of TD
#############################################################################################

## Variance of time to detection (TD) vs number of sampling sites
x = number_of_sites_vec

y1 = plot_data_w_300s_20b.td_var
y1_lower_errors = plot_data_w_300s_20b.td_var_error_lower
y1_upper_errors = plot_data_w_300s_20b.td_var_error_upper

y2 = plot_data_w_500s_20b.td_var
y2_lower_errors = plot_data_w_500s_20b.td_var_error_lower
y2_upper_errors = plot_data_w_500s_20b.td_var_error_upper

y3 = plot_data_w_1000s_20b.td_var
y3_lower_errors = plot_data_w_1000s_20b.td_var_error_lower
y3_upper_errors = plot_data_w_1000s_20b.td_var_error_upper

color1 = :lightsalmon
color2 = :red #:indianred1
color3 = :red4 #:red

p_var = Plots.plot( [x ; x ], [y1; y1]
              , yerror = (y1_lower_errors, y1_upper_errors)
            , seriestype = :line, color = color1, linecolor = color1, errorbarcolor = color1
            , marker = (:circle, 6, color1, 0.5)
            , label  = "300 samples per week"
            , ylabel = "Variance of time to detection"
            , xlabel = "Number of ICU sampling sites"
            , legend = :bottomright
            , foreground_color_legend = nothing
            , background_color_legend = nothing
            , alpha = 0.5
            , ylimit = [60,90]
                )

Plots.plot!( [x ; x ], [y2; y2]
              , yerror = (y2_lower_errors, y2_upper_errors)
            , seriestype = :line, color = color2, linecolor = color2, errorbarcolor = color2
            , marker = (:square, 6, color2, 0.5)
            , label  = "500 samples per week"
            , alpha = 0.5
            )

Plots.plot!( [x ; x ], [y3; y3]
            , yerror = (y3_lower_errors, y3_upper_errors)
            , seriestype = :line, color = color3, linecolor = color3, errorbarcolor = color3
            , marker = (:diamond, 6, color3, 0.5)
            , label  = "1000 samples per week"
            , alpha = 0.5
            )

Plots.savefig( joinpath( root_dir
                        , "NBPMscape/paper/analysis/2_sampling_analysis/covid_like/optimal_number_of_sites"
                        , "ICU_TD_variance_winter_HCGS_20_blocks.png" ) )

p_combined = Plots.plot( p_median, p_var, layout = (2,1), size = (700,700) )
Plots.savefig( joinpath( root_dir
                        , "NBPMscape/paper/analysis/2_sampling_analysis/covid_like/optimal_number_of_sites"
                        , "ICU_TD_median_and_variance_winter_HCGS_20_blocks.png" ) )
Plots.savefig( joinpath( root_dir
                        , "NBPMscape/paper/results/figures_main/raw"
                        , "ICU_TD_median_and_variance_winter_HCGS_20_blocks.png" ) )

