#= Miscellaneous functions

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

=#

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
    return allocation
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
            rolling_window  Number of days to include in the rolling window
            maxtime         maxtime used when running simulation to create sims
            plot_save_name  path and filename for output plot .png file

Returns     Line plot of the rolling mean age of infected individuals by time since importation to the UK.
            Plot is saved to .png file

Examples
            # Load file
            sims = load("covidlike-1.3.6-sims-nreps10.jld2", "sims")
            # Run function
            severity_rolling_mean(; sims = sims, rolling_window = 3, maxtime = 90, plot_save_name = "examples/tinf_age_rolling_mean_3d.png")    
            severity_rolling_mean(; sims = sims, rolling_window = 10, maxtime = 90, plot_save_name = "examples/tinf_age_rolling_mean_10d.png")    

"""    

function severity_rolling_mean(; sims, rolling_window = 3, maxtime = 90, plot_save_name)    

    ### Combine data from simulation replicates
    nreps = length(sims)
    tinf_age_severity_dfs = Vector{DataFrame}(undef, nreps) 
    for i in 1:length(sims)
        tinf_age_severity_dfs[i] = sims[i].G[:,[:tinf,:infectee_age,:severity]]
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
            - age_group_width::Integer  bin width (e.g., 5 years)
            - min_age::Integer          minimum age
            - max_age::Integer          maximum age
            - plot_file_prefix          path and filename prefix to save plots    

Returns     Three plots in .png files as described above

Example     
            # Load file
            sims = load("covidlike-1.3.6-sims-nreps10.jld2", "sims")
            # Run function to generate plots
            tinf_by_age(; sims = sims, age_group_width = 5, min_age = 0, max_age = 100
                     , plot_file_prefix = "examples/test_prefix"
                     )
"""
function tinf_by_age(; sims
                     , age_group_width::Integer = 5, min_age::Integer = 0, max_age::Integer = 100
                     , plot_file_prefix
                     )
    
    # Collate G dataframes from each simulation replicate in sims and store in a vector
    # Also trim columns to only those required: tinf, infectee_age and severity
    
    tinf_age_severity_dfs = Vector{DataFrame}(undef, length(sims)) 
    for i in 1:length(sims)
        tinf_age_severity_dfs[i] = sims[i].G[:,[:tinf,:infectee_age,:severity]]
    end
    
    # Compute age bins and labels
    age_bins = collect( min_age : age_group_width : (max_age + age_group_width))
    labels = ["$(age_bins[i])-$(age_bins[i+1]-1)" for i in 1:length(age_bins)-1]

    # Assign age groups
    for i in 1:length(sims)
        tinf_age_dfs[i][!, :age_group] = CategoricalArrays.cut(tinf_age_dfs[i].infectee_age, age_bins; labels=labels)
        tinf_age_severity_dfs[i][!, :age_group] = CategoricalArrays.cut(tinf_age_severity_dfs[i].infectee_age, age_bins; labels=labels)
    end
    
    ## Plot infection time against age group for each individual simulation replicate

    # Violin plots
    #plots = [@df tinf_age_severity_dfs[i] StatsPlots.violin(:age_group, :tinf, legend=false) for i in 1:length(tinf_age_severity_dfs)]
    #Plots.plot(plots..., layout=(10,1), size=(1000, 1000))
    
    # Boxplots
    #using Plots.PlotMeasures
    plots = [@df tinf_age_dfs[i] StatsPlots.boxplot(:age_group, :tinf
                                                    #, xlabel="Age Group", ylabel="Infection Time (days)"
                                                    #, title="Distribution of Infection Times by Age Group"
                                                    , size=(1000,400), legend = false, left_margin=10mm) for i in 1:length(tinf_age_dfs)]
    Plots.plot(plots..., layout=(10,1), size=(1200, 2000), left_margin = 10mm)
    
    # Save plot to file
    Plots.savefig("$(plot_file_prefix)_tinf_age.png")

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

    for i in 1:1:length(severity_dfs_vec)
        severity_types[i] = only(unique( severity_dfs_vec[i][:,:severity] ))
    end
    
    plots = [@df severity_dfs_vec[i] StatsPlots.boxplot(:age_group, :tinf
                                                        #, xlabel="Age Group", ylabel="Infection Time (days)"
                                                        , title = severity_types[i]#string(only(unique( severity_dfs_vec[i][:,:severity] )))
                                                        , size=(1000,400), legend = false, left_margin=10mm, color = :lightgreen) for i in 1:length(severity_dfs_vec)]
    
    Plots.plot(plots..., layout=(length(severity_dfs_vec),1), size=(1200, 2000), left_margin = 10mm)
    
    # Save plot to file
    Plots.savefig("$(plot_file_prefix)_tinf_age_severity_nrep$(length(sims)).png")

end