#= Miscellaneous functions

- median_ci_bootstrap:    Computes median and (1-alpha)% bootstrap estimate

- mean_ci_bootstrap:      Computes mean and (1-alpha)% bootstrap estimate

- var_ci_bootstrap:       Computes variance and (1-alpha)% bootstrap estimate

- allocate_with_rounding:   allocates a number across a number of categories based on weights
                            ensuring integer values are allocated and the sum of allocations
                            is equal to the original total, e.g. total number of samples
                            allocated across NHS Trusts but the allocations must be integer values
                            and the sum must be equal to the total

- generation_time:  Computes generation times from df containing data on infector (:donor), infectee (:recipient) and time
                    of infection (:timetransmission) for multiple simulation replicates.
                    A Gamma distribution is fitted to the generation times and plotted.
                    Mean and median generation times are computed from the fitted distribution and
                    the raw results.

- severity_rolling_mean     Produces a line plot of the rolling mean age of infected individuals disaggregated by
                            infection severity. This is done by:
                            - combining data in the G dataframe from multiple simulation replicates in an object
                              named 'sims'
                            - disaggregating by infection severity
                            - computing the rolling mean age between time of importation to the UK and the 
                              maximum time of infection in the simulation (maxtime)

- tinf_by_age   Generate three plots:
                (1) boxplots of time of infection vs age group for individual simulation replicates
                (2) boxplots of time of infection vs age group for individual simulation replicates combined
                (3) boxplots of time of infection vs age group disaggregated by infection severity

- kernel_box_jitter_plot    Generates plot of times to detection (TD) displaying:
                            (1) line of distribution
                            (2) boxplot
                            (3) jitter with points
                            Information will also be annotated on the plot. 
                            This includes statistics computed from the data (e.g. mean, variance, etc)
                            and text supplied as arguments to the function.

- root_folder     Function to define root directory/folder

=#

"""
Function    median_ci_bootstrap(;vec, n_boot=1000, alpha=0.05)

Description     Function to compute median and (1-alpha)% CI via bootstrap

Arguments   vec::Vector     Vector of values for median and CI to be computed on
            n_boot::Int64   Number of bootstrap repeats required
            alpha::Float64  Quantile value, i.e. 0.05 for 95% CI

Returns     Median of vector and the lower and upper quantiles or median
            values generated using bootstrapping. Returned as a NamedTuple.

Examples    # Compute median and upper and lower values for 95% CI
            # for each TD results df in a dictionary containing different scenarios
            v = rand(100000)
            median_ci_bootstrap( vec = v, n_boot = 1000, alpha = 0.05 ) 
"""
function median_ci_bootstrap(; vec::Vector, n_boot::Int64=1000, alpha::Float64=0.05)
    med = median(vec)
    boot_samples = [median(rand(vec, length(vec))) for _ in 1:n_boot]
    lower = quantile(boot_samples, alpha/2)
    upper = quantile(boot_samples, 1 - alpha/2)
    return (median = med, lower = lower, upper = upper)
end

"""
Function    mean_ci_bootstrap(vec; n_boot=2000, alpha=0.05)

Description     Function to compute mean and (1-alpha)% bootstrap estimate

Arguments   vec::Vector     Vector of values for median and bootstrap estimates to be computed on
            n_boot::Int64   Number of bootstrap repeats required
            alpha::Float64  Quantile value, i.e. 0.05 for 95% boostrap estimate

Returns     Mean of vector and the lower and upper quantiles or mean
            values generated using bootstrapping. Returned as a NamedTuple.

Examples    # Compute mean and upper and lower values for 95% bootstrap estimate
            # for each TD results df in a dictionary containing different scenarios
            v = rand(100000)
            mean_ci_bootstrap( vec = v, n_boot = 1000, alpha = 0.05 ) 
"""
function mean_ci_bootstrap(; vec::Vector, n_boot::Int64=1000, alpha::Float64=0.05)
    μ = mean(vec)
    boot_samples = [mean(rand(vec, length(vec))) for _ in 1:n_boot]
    lower = quantile(boot_samples, alpha/2)
    upper = quantile(boot_samples, 1 - alpha/2)
    return (mean = μ, lower = lower, upper = upper)
end

"""
Function    var_ci_bootstrap(vec; n_boot=1000, alpha=0.05)

Description     Function to compute variance and (1-alpha)% bootstrap estimate

Arguments   vec::Vector     Vector of values for variance and bootstrap estimates to be computed on
            n_boot::Int64   Number of bootstrap repeats required
            alpha::Float64  Quantile value, i.e. 0.05 for 95% boostrap estimate

Returns     Variance of vector and the lower and upper quantiles 
            generated using bootstrapping. Returned as a NamedTuple.

Examples    # Compute variance and upper and lower values for 95% bootstrap estimate
            v = rand(100000)
            var_ci_bootstrap( vec = v, n_boot = 1000, alpha = 0.05 ) 
"""
function var_ci_bootstrap(; vec::Vector, n_boot::Int64=1000, alpha::Float64=0.05)
    vec_var = var(vec)
    boot_samples = [ var( rand(vec, length(vec) ) ) for _ in 1:n_boot]
    lower = quantile(boot_samples, alpha/2)
    upper = quantile(boot_samples, 1 - alpha/2)
    return (variance = vec_var, lower = lower, upper = upper)
end

"""
Function:       allocate_with_rounding

Description:    Allocates a number across a number of categories/groups based on weights
                ensuring integer values are allocated and the sum of allocations
                is equal to the original total, e.g. total number of samples
                allocated across NHS Trusts but the allocations must be integer values
                and the sum must be equal to the total.

Arguments:      total::Int          Number to be allocated across categories/groups
                weights::Vector     Weightings for each category/group

Returns:        Vector of integer values

Examples:       alloc = allocate_with_rounding( total = 10, weights = [0.15, 0.15, 0.3, 0.4])
                # Checks
                println(alloc)  
                sum(alloc)
                alloc / sum(alloc)
"""
function allocate_with_rounding(;total, weights)
    weightsum = sum(weights)
    # Ideal unrounded allocations
    exact = total .* (weights ./ weightsum)
    # Integer part (floor)
    allocation = floor.(Int, exact)
    # Compute total remainder after integer (floor) allocation
    remainder = total - sum(allocation)
    # Fractional remainders
    fractional = exact .- allocation
    # Find indices of categories/groups with the largest fractional remainders (for distributing leftover units)
    idx = partialsortperm(fractional, rev=true, 1:remainder)
    # Add 1 to the allocations of categories/groups with the largest fractional parts
    allocation[idx] .+= 1
    return Int64.(allocation)
end

using DataFrames
using Distributions
using StatsBase
using Plots

"""
Function:       generation_time

Description:    Computes generation times from df containing data on infector (:donor), infectee (:recipient) and time
                of infection (:timetransmission) for multiple simulation replicates.
                A Gamma distribution is fitted to the generation times and plotted.
                Mean and median generation times are computed from the fitted distribution and
                the raw results.

Arguments:      tinf_df::DataFrame  Dataframe with three columns: :donor, :recipient, :timetransmission
                
Returns:        Plot to screen and dataframe with mean and median values for generation time
                from data and fitted Gamma distribution

Examples:       # Load transmission data
                sims_simid_tinf_df = load("C:/Users/kdrake/OneDrive - Imperial College London/Documents/mSCAPE/3_results/from_hpc/2025_10_maxtime_maxcases/1101621/covidlike-1.3.14-sims_simid_tinf_nrep1200_1101621.jld2", "sims_simid_tinf_df")
                # Run generation time function
                Tg_results = generation_time( tinf_df = sims_simid_tinf_df )
                # Returns dataframe

                # Save returned plot
                savefig("examples/generation_time.png")

"""

function generation_time(; tinf_df )
    
# Compute generation times
    # For each row, find infector's (donor's) infection time and subtract from the time of transmission to the infectee (recipient)
    infector_times = Dict(row.recipient => row.timetransmission for row in eachrow(tinf_df))
    
    # Filter out imports which have no donor ID
    tinf_df_wo_imports = filter( row -> !ismissing(row[:donor]), tinf_df)

    generation_times = [row.timetransmission - infector_times[row.donor] for row in eachrow(tinf_df_wo_imports)]
    
    # Fit Gamma distribution
    fit_gamma = fit(Gamma, generation_times)
    println("Fitted Gamma parameters: shape = $(fit_gamma.α), scale = $(fit_gamma.θ)")

    # Compute mean and median
    mean_gt = mean(fit_gamma)
    median_gt = quantile(fit_gamma, 0.5)
    println("Mean generation time: ", mean_gt)
    println("Median generation time: ", median_gt)

    # Plot histogram and fitted density
    histogram(generation_times, normalize=true, alpha=0.5, label="Observed", xlabel="Generation Time", ylabel="Density")
    x = range(minimum(generation_times), stop=maximum(generation_times), length=200)

    # Add mean and median generation times for the simulated data and Gamma fit to a df
    df = DataFrame( statistic = [ "mean" , "median" ]
                    , simulated_data = [ mean(generation_times), median(generation_times)]
                    , Gamma_fit = [ mean_gt, median_gt ] 
                    )
    return(df)

end

"""
Function        severity_rolling_mean

Description     Produces a line plot of the rolling mean age of infected individuals disaggregated by infection severity.
                This is done by:
                - combining data in the G dataframe from multiple simulation replicates in an object named 'sims'
                - disaggregating by infection severity
                - computing the rolling mean age between time of importation to the UK and the maximum time of infection
                    in the simulation (maxtime)

Arguments   sims            object containing simulation data output from simtree or simforest, including G dataframes
            format_G        Determine whether the 'sims' object holds the data in dataframes named 'G' or not
            rolling_window  Number of days to include in the rolling window
            maxtime         maxtime used when running simulation to create sims
            plot_save_name  path and filename for output plot .png file

Returns     Line plot of the rolling mean age of infected individuals by time since importation to the UK.
            Plot is saved to .png file

Examples
            # Load file
            sims = load("covidlike-1.3.6-sims-nreps10.jld2", "sims")
            # Run function
            severity_rolling_mean(; sims = sims, format_G = true, rolling_window = 3,  maxtime = 90, plot_save_name = "examples/tinf_age_rolling_mean_3d.png")    
            severity_rolling_mean(; sims = sims, format_G = true, rolling_window = 10, maxtime = 90, plot_save_name = "examples/tinf_age_rolling_mean_10d.png")    

"""    

function severity_rolling_mean(; sims, rolling_window = 3, maxtime = 90, format_G::Bool = false, plot_save_name)    

    ### Combine data from simulation replicates
    nreps = length(sims)
    tinf_age_severity_dfs = Vector{DataFrame}(undef, nreps) 
    for i in 1:length(sims)
        if format_G
            tinf_age_severity_dfs[i] = sims[i].G[:,[:tinf,:infectee_age,:severity]]
        else
            tinf_age_severity_dfs[i] = sims[i][:,[:tinf,:infectee_age,:severity]]
        end
    end
    combined_df = vcat(tinf_age_severity_dfs...)

    ### Split out by severity
    #groups = groupby(combined_df, :severity)
    severity_dict = Dict(unique(combined_df.severity) .=> [g for g in groupby(combined_df, :severity)])
    # Create vector of dfs split out by severity
    severity_dfs_vec = collect(values(severity_dict))

    ### Compute rolling averages and plot
    for j in 1:length(severity_dfs_vec) #j=1
        
        # Create temporary df
        df = severity_dfs_vec[j]

        # Create df to be filled
        rolling_avg_age = DataFrame( time = collect(1:1:maxtime) # maxtime = 90
                                    ,rolling_10d_mean_age = Vector{Float64}(undef, maxtime) 
                                    )

        # Calculate rolling mean age based on time of infection
        for i in 1:maxtime # i=89
            # Define start of rolling window
            start_window = i - rolling_window +1 # +1 because includeing the current day
            # Select ages where time of infection (tinf) is within the window
            ages_in_window = df.infectee_age[(df.tinf .>= start_window) .& (df.tinf .<= i)]
            # Compute mean for rolling window
            rolling_avg_age[i,:rolling_10d_mean_age] = mean(ages_in_window)
        end
        
        # Plot rolling means
        if j == 1 # First plot includes labels and titles
            Plots.plot(rolling_avg_age[:,:time], rolling_avg_age[:,:rolling_10d_mean_age], xlimit=[0,maxtime], ylimit=[0,100]
                                ,label = string(only(unique( df[:,:severity] )))
                                ,linewidth=2
                                , xlabel = "Time since importation into the UK (days)"
                                , ylabel = "$(rolling_window)-day rolling mean age of \n infected individuals (years)"
                            )
        else
            Plots.plot!(rolling_avg_age[:,:time], rolling_avg_age[:,:rolling_10d_mean_age], xlimit=[0,maxtime], ylimit=[0,100]
                             ,label = string(only(unique( df[:,:severity] )))
                            ,linewidth = 2
                            )
        end

    end

    # Save plot to file
    Plots.savefig( plot_save_name )

end

"""
Function        tinf_by_age

Description     Generate three plots:
                (1) boxplots of time of infection vs age group for individual simulation replicates
                (2) boxplots of time of infection vs age group for individual simulation replicates combined
                (3) boxplots of time of infection vs age group disaggregated by infection severity
                
Arguments   - sims                      object containing simulation data output from simtree or simforest, including G dataframes
            - format_G::Bool            Determine whether the 'sims' object holds the data in dataframes named 'G' or not
            - age_group_width::Integer  bin width (e.g., 5 years)
            - min_age::Integer          minimum age
            - max_age::Integer          maximum age
            - plot_file_prefix          path and filename prefix to save plots    

Returns     Three plots in .png files as described above, and a vector of dataframes. 
            Each df contains data relating to a different infection severity, with information on
            the time of infection, infectee age, and age group. Example below for 'moderate_ED' severity.
            Outputting the data as a vector of dfs allows more flexibility to change how the data is plotted.

                  Row │ tinf     infectee_age  severity     age_group 
                      │ Float64  Int8          Symbol       Cat…      
            ──────────┼───────────────────────────────────────────────
                    1 │ 27.2723            65  moderate_ED  65-69
                    2 │ 30.5866            44  moderate_ED  40-44
                    3 │ 38.0046            80  moderate_ED  80-84
                    4 │ 38.35              27  moderate_ED  25-29

Example     
            # Load file
            sims = load("covidlike-1.3.6-sims-nreps10.jld2", "sims")
            # Run function to generate plots
            tinf_by_age(; sims = sims, format_G::Bool = true
                        ,age_group_width = 5, min_age = 0, max_age = 100, format_G
                        , plot_file_prefix = "examples/test_prefix"
                        )
"""
function tinf_by_age(; sims, format_G::Bool = true
                     , age_group_width::Integer = 5, min_age::Integer = 0, max_age::Integer = 100
                     , plot_file_prefix
                     )
    
    # Collate G dataframes from each simulation replicate in sims and store in a vector
    # Also trim columns to only those required: tinf, infectee_age and severity
    
    tinf_age_severity_dfs = Vector{DataFrame}(undef, length(sims)) 
    for i in 1:length(sims)
        if format_G
            tinf_age_severity_dfs[i] = sims[i].G[:,[:tinf,:infectee_age,:severity]]
        else
            tinf_age_severity_dfs[i] = sims[i][:,[:tinf,:infectee_age,:severity]]
        end
    end
    
    # Compute age bins and labels
    age_bins = collect( min_age : age_group_width : (max_age + age_group_width))
    labels = ["$(age_bins[i])-$(age_bins[i+1]-1)" for i in 1:length(age_bins)-1]

    # Assign age groups
    for i in 1:length(sims)
        #tinf_age_dfs[i][!, :age_group] = CategoricalArrays.cut(tinf_age_dfs[i].infectee_age, age_bins; labels=labels)
        tinf_age_severity_dfs[i][!, :age_group] = CategoricalArrays.cut(tinf_age_severity_dfs[i].infectee_age, age_bins; labels=labels)
    end
    
    ## Plot infection time against age group for each individual simulation replicate

    # Violin plots
    #plots = [@df tinf_age_severity_dfs[i] StatsPlots.violin(:age_group, :tinf, legend=false) for i in 1:length(tinf_age_severity_dfs)]
    #Plots.plot(plots..., layout=(10,1), size=(1000, 1000))
    
    # Boxplots
    #using Plots.PlotMeasures
    #plots = [@df tinf_age_dfs[i] StatsPlots.boxplot(:age_group, :tinf
    #                                                #, xlabel="Age Group", ylabel="Infection Time (days)"
    #                                                #, title="Distribution of Infection Times by Age Group"
    #                                                , size=(1000,400), legend = false, left_margin=10mm) for i in 1:length(tinf_age_dfs)]

    # Only plot separate simulation replicates if there are 10 or less
    if length(sims) <= 10
        plots = [@df tinf_age_severity_dfs[i] StatsPlots.boxplot(:age_group, :tinf
                                                        #, xlabel="Age Group", ylabel="Infection Time (days)"
                                                        #, title="Distribution of Infection Times by Age Group"
                                                        , size=(1000,400), legend = false, left_margin=10mm) for i in 1:length(tinf_age_severity_dfs)]
        #Plots.plot(plots..., layout=(10,1), size=(1200, 2000), left_margin = 10mm)
        Plots.plot(plots..., layout=(length(sims),1), size=(1200, 2000), left_margin = 10mm)
        
        # Save plot to file
        Plots.savefig("$(plot_file_prefix)_tinf_age.png")
    end

    # Plot data combined from all sim reps
    combined_df = vcat(tinf_age_severity_dfs...)

    @df combined_df StatsPlots.boxplot(:age_group, :tinf
                                                    #, xlabel="Age Group", ylabel="Infection Time (days)"
                                                    #, title="Distribution of Infection Times by Age Group"
                                                    , size=(1000,400), legend = false, left_margin=10mm)
    # Save plot to file
    Plots.savefig("$(plot_file_prefix)_tinf_age_nrep$(length(sims)).png")

    ## Plot infection time against age group disggregated by severity but with simulation replicates combined

    # Group by severity
    severity_dict = Dict(unique(combined_df.severity) .=> [g for g in groupby(combined_df, :severity)])
    # Add separate dataframes for each infection severity level to a vector of dataframes
    severity_dfs_vec = collect(values(severity_dict))

    # Create vector to store severity types
    severity_types = Vector{Symbol}(undef, length(severity_dfs_vec)) 

    for i in 1:1:length(severity_dfs_vec)
        severity_types[i] = only(unique( severity_dfs_vec[i][:,:severity] ))
    end
    
    plots = [@df severity_dfs_vec[i] StatsPlots.boxplot(:age_group, :tinf
                                                        #, xlabel="Age Group", ylabel="Infection Time (days)"
                                                        , title = severity_types[i]#string(only(unique( severity_dfs_vec[i][:,:severity] )))
                                                        , size=(1000,400), legend = false, left_margin=10mm, color = palette(:default)[i]) for i in 1:length(severity_dfs_vec)]
    
    Plots.plot(plots..., layout=(length(severity_dfs_vec),1), size=(1200, 2000), left_margin = 10mm)
    
    # Save plot to file
    Plots.savefig("$(plot_file_prefix)_tinf_age_severity_nrep$(length(sims)).png")

    return( severity_dfs_vec )
end


"""
TODO function description
Example
        convert_t_to_day_of_week(initial_dow = 1, t = 6.0)
"""
function convert_t_to_day_of_week(;t, initial_dow)
    d = Int( floor( (initial_dow-1) + t ) % 7  ) + 1
    d
end


"""
Function        kernel_box_jitter_plot

Description     Generates plot of times to detection (TD) displaying:
                (1) line of distribution
                (2) boxplot
                (3) jitter with points
                Information will also be annotated on the plot. 
                This includes statistics computed from the data (e.g. mean, variance, etc)
                and text supplied as arguments to the function.
                Requires at least: StatsPlots, KernelDensity, Statistics, Distributions, Random

Arguments   
            # Data
            - x
            - plot_color
            # Data for label/text on plot
            - samp_strategy
            - n_samples
            - n_sites
            - n_blocks=""
             - season
            # Axis information
            - x_label = "Time since first imported infection (days)"
            - x_min = 0
            - x_max = 80
            - x_ticks = true
            - annot_x_pos = "left"
            - dens_zero_below_zero = true

Returns     Plot as described above

Examples    
            # Example 1
            # Define time to detection dataset
            df = icu_tds_df_dict["icu_tds_seasonal_ari_1440_500_sampling_sites_30blocks_10sites"]
            x = df[:,:ICU_TD]
            # Define colour
            hcgs_color = :mediumpurple3 #(i.e. :royalblue1, :orangered,:gold, :mediumpurple3)
            # Run function to generate plots
            p1 = kernel_box_jitter_plot(;x = icu_tds_df_dict["icu_tds_seasonal_ari_1440_300_sampling_sites_10blocks_10sites"][:,:ICU_TD]
                                        , plot_color = hcgs_color
                                        ,samp_strategy="HCGS"
                                        ,n_samples=300
                                        ,n_sites=10
                                        ,n_blocks=10
                                        ,season="Winter")
            
            # Example 2
            kernel_box_jitter_plot(;x = rand(Gamma(2.0,3.0), 1000) .+30
                        , plot_color = :red
                        ,samp_strategy="HCGS"
                        ,n_samples=300
                        ,n_sites=10
                        ,n_blocks=10
                        ,season="Winter")

            
"""
function kernel_box_jitter_plot(;x, plot_color, samp_strategy, n_samples, n_sites, n_blocks="", season
                                , x_label = "Time since first imported infection (days)", x_min = 0, x_max = 80, x_ticks = true, annot_x_pos = "left"
                                , dens_zero_below_zero = true )#, main_title)
    
    # Available named colors
    # https://juliagraphics.github.io/Colors.jl/stable/namedcolors/
    # For color blindness, four good colors are orange, blue, purple and yellow (i.e. :royalblue1, :orangered,:gold, :mediumpurple3)

    # Summary stats (mean + 95% CI)
    n   = length(x)
    μ   = mean(x)
    #σ   = std(x)
    #se  = σ / sqrt(n)
    #α   = 0.05
    #tval = quantile(TDist(n - 1), 1 - α/2)
    #ci  = (μ - tval * se, μ + tval * se)
    # Bootstrap estimate of the median TD, computed using n resamples with replacement and taking the 2.5% and 97.5% quantiles
    bootstrap_estimate = NBPMscape.mean_ci_bootstrap( vec = x, n_boot = 1000, alpha = 0.05 ) # [ mean, lower, upper ] mean_ci_bootstrap function is in misc_functions.jl

    # Kernel density
    kd   = kde(x)
    #kd   = kde(x; boundary=(0, x_max)) # Specify boundaries for density
    xs   = kd.x
    dens = kd.density ./ maximum(kd.density)  # normalize to 0–1
    # Force density to 0 for x < 0
    if dens_zero_below_zero 
        dens[xs.< 0].= 0
    end

    # Kernel density with boundary correction (values cannot be below 0)
    # For information on KDE bias corrction by reflection, see https://search.r-project.org/CRAN/refmans/evmix/html/bckden.html and Boneva, L.I., Kendall, D.G. and Stefanov, I. (1971). Spline transformations: Three new diagnostic aids for the statistical data analyst (with discussion). Journal of the Royal Statistical Society B, 33, 1-70.
    #x_reflected = vcat(x, -x)  # Reflect across x=0
    #kd   = kde(x_reflected)
    #xs   = kd.x
    #dens = kd.density
    ## Keep only x >= 0 and double the density (due to reflection)
    #mask = xs.>= 0
    #xs = xs[mask]
    #dens = 2 .* dens[mask]
    #dens = dens./ maximum(dens)  # normalize to 0–1

    # Visual placement parameters
    yloc           = 1.0           # vertical "track" for this raincloud
    jitter         = 0.06          # vertical jitter for points
    density_gap    = 0.35          # distance between box midline and density start
    density_height = 0.45          # max density height above its start

    default(legend=false, background_color=:white, grid=false, framestyle=:box)
    p = Plots.plot(; size=(900, 280))#, main = main_title)

    # 1) BOX: use a categorical position at y=1 via group=fill(1, n)
    g = fill(1, n)
    #if xticks == () : 
    StatsPlots.boxplot!(p, x;#, g;
        orientation = :horizontal
        ,fillalpha   = 0.14
        ,fillcolor   = plot_color #:orangered
        ,linecolor   = :black
        ,whiskercolor = :black
        ,xlabel = x_label #"Time since first imported infection (days)"
        ,xticks = x_ticks
        ,guidefontsize = 14     # Axis labels font size
        ,tickfontsize = 12      # Tick labels font size
        ,outliers    = false # These will be plotted in 2) scatter below
    )

    # 2) SCATTER (overlays the boxplot): draw AFTER the box
    scatter!(p, x, yloc .+ randn(n) .* jitter;
        ms = 5, mc = plot_color #:orangered
        , ma = 0.40, msw = 0
        ,alpha = 0.5
    )

    # 3) DENSITY (above): draw at y = yloc + density_gap + scaled_height
    ys = yloc .+ density_gap .+ dens .* density_height
    plot!(p, xs, ys; lw=3, c=plot_color) #:orangered)

    # 4) Mean tick and annotation
    plot!(p, [μ, μ], [yloc - 0.18, yloc + 0.18]; c=:black, lw=3)
    #annot_str = "Mean: $(round(μ, digits=1))\n95% CI: [$(round(ci[1], digits=1)), $(round(ci[2], digits=1))]"
    #annot_str = "Mean: $(round(μ, digits=1))\n[95% bootstrap est: $(round(bootstrap_estimate[1], digits=1)), $(round(bootstrap_estimate[2], digits=1))]\nVariance: $(round(var(x),digits=1))\nStandard deviation: $(round(std(x),digits=1))"
    annot_str = "Mean: $(round(μ, digits=1))\n[95%: $(round(bootstrap_estimate[2], digits=1)), $(round(bootstrap_estimate[3], digits=1))]\nVar: $(round(var(x),digits=1))\nSt.Dev.: $(round(std(x),digits=1))"
    xmin, xmax = (x_min,x_max) #(0,80) #extrema(x)
    xpad       = 0.06 * (xmax - xmin)
    annot_x    = xmin + 0.5 * (xmax - xmin) #xmin + 0.75 * (xmax - xmin)
    annot_y    = yloc + density_gap + density_height + 0.2 #0.1
    #annotate!(p, annot_x, annot_y, text(annot_str, 12, :right, :black))
    annotate!(p, xmax, annot_y, text(annot_str, 12, :right, :black))
    #scatter!(p, [annot_x - 0.02*(xmax-xmin)], [annot_y]; marker=:rect, ms=10, mc=plot_color #:orangered
    #            , ma=0.7, msw=0)

    # Add text describing data
    if n_blocks == ""
        blocks_text = ""
    else
        blocks_text = "blocks in HC"
    end

    annot_data_desc = "$(samp_strategy)\n$(season)\n$(n_samples) samples\n$(n_sites) sites\n$(n_blocks) $(blocks_text)"
    
    if annot_x_pos == "left"
        annot_x = xmin
    elseif annot_x_pos == "middle"
        annot_x = xmin + 0.25 * (xmax - xmin) #xmin + 0.75 * (xmax - xmin)
    end
    
    annotate!(p, annot_x, annot_y, text(annot_data_desc, 12, :left, :black))
    
    samp_strategy, n_samples, n_sites, season

    # Cosmetics
    xlims!(p, xmin - xpad, xmax + xpad)
    yticks!(p, ([1.0], [""]))
    ylims!(p, 0.3, 2.5) #ylims!(p, 0.3, 2.5)

    #display(p)
    return p
end



"""
Function        root_folder

Description     Function to define root directory/folder

Arguments   None required
            

Returns     the path to the Project.toml file, which is can then be used as the root directory
            from which to locate other files and folders

Examples    root_folder()
            # Example return
            "C:/Users/user/Documents/GitHub/NBPMscape"
            # This can then be combined using joinpath()
            results_df = CSV.read( joinpath( root_dir, "scripts/paper/2_sampling_analysis/covid_like/1717144_1719024_1719029_analysis/geo_optimisation"
                                            ,"icu_td_results_df.csv" )
                                    , DataFrame )
                        
"""
function root_folder()
    dir = @__DIR__
    while !isfile(joinpath(dir, "Project.toml"))
        parent = dirname(dir)
        parent == dir && error("Project.toml not found")
        dir = parent
    end
    return dir
end