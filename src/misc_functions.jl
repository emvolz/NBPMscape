#= Miscellaneous functions

- allocate_with_rounding:   allocates a number across a number of categories based on weights
                            ensuring integer values are allocated and the sum of allocations
                            is equal to the original total, e.g. total number of samples
                            allocated across NHS Trusts but the allocations must be integer values
                            and the sum must be equal to the total

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