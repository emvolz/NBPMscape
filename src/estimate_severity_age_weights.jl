#=
Functions to estimate post simulation infection severity and age distributions
(1) Function to estimate frequency/percentage of infection in each severity category
    - inf_severity_estimate
(2) Age distribution
    - inf_age_estimate
    - group_ages
=#

#using JLD2
#using DataFrames
#using StatsPlots

"""
Function: inf_severity_estimate

This function takes a .jld2 file containing a single simulation or multiple simulation replicates and
computes the percentage of infections that are of each severity category in each simulation replicate.
It plots the distribution of these percentages across all simulation repliacetes and computes the median
percentage value for the each infection severity category.

# Arguments
    'sims_file':            .jld2 file containing simulation replicates as output by {simtree} or {simforest}
                            For example, "sims", or "sims_G_gp_filter" or "sims_G_icu_filter" if pre-filtered for certain cases.
    'sim_object_name':      The name of the object saved to the .jld2 file. 
                            For example, "sims", "sims_G_icu_filter", or "sims_G_gp_filter" for files filtered using {sims_filter}.
                            Without this name the data cannot be loaded. 
    'infection_severity':   Vector of infection severity categories (formatted as symbols), which must match those in the G dataframe 
                            in the simulation results file.
# Returns
    Infection severity distribution for each simulation in file is shown in a violin plot. 
    This can then be saved to file as shown below.
    Median percentages also printed to screen and also returned from function.

# Example
    # Estimate infection severity weights
    inf_severity_estimate( sims_file = "covidlike-1.3.1-sims-nrep1000_955898_2.jld2"
                        , sims_object_name = "sims" 
                        , infection_severity = [:asymptomatic, :mild, :moderate, :severe, :verysevere ]
                        )
     Row │ infection_severity  median_percentage_across_simreps 
         │ Symbol              Float64
    ─────┼──────────────────────────────────────────────────────
       1 │ asymptomatic                               52.2536
       2 │ mild                                       38.9568
       3 │ moderate                                    4.70901
       4 │ severe                                      3.40471
       5 │ verysevere                                  0.541461

    # Save figure
    savefig("examples/infection_severity_weights.png")

"""

function inf_severity_estimate(;  sims_file
                                , sims_object_name = "sims"
                                , infection_severity = [:asymptomatic, :mild, :moderate, :severe, :verysevere ]
                                )

    # Load simulation file
    sims = load(sims_file, sims_object_name)

    # Column to search in
    column_name = :severity 
    # Possible values for infection severity
    #infection_severity = [:asymptomatic, :mild, :moderate, :severe, :verysevere ]  

    # Initialise vectors to store results
    inf_severity_counts = [] ; inf_severity_percs = [] ; total_infections = []

    # Loop through the simulation replicates and gather frequency and percentages for each severity type
    for s in sims
        counts = [count(==(val), s.G[!, column_name]) for val in infection_severity]
        push!(inf_severity_counts, counts)

        total = size(s.G,1)
        push!(total_infections, total)

        percentages = counts ./ total #[count(==(val), s[!, column_name]) / total for val in infection_severity]
        push!(inf_severity_percs, percentages)
    end

    # Change format of data so it can be plotted
    inf_severity_percs_separated = [Float64[] for _ in 1:length(infection_severity)]

    for i in 1:length(infection_severity)
        for s in inf_severity_percs
            push!(inf_severity_percs_separated[i], s[i])
        end
    end

    # Series names
    #infection_severity_names = ["asymptomatic" "mild" "moderate" "severe" "verysevere"]
    infection_severity_names = string.(infection_severity) #["asymptomatic" "mild" "moderate" "severe" "verysevere"]
    infection_severity_names = reshape(infection_severity_names,1,5)
    
    # Create the violin plot
    colors = [:red, :blue, :green, :orange, :purple]
    violin( infection_severity_names
          , inf_severity_percs_separated*100
          , legend=false
          , palette = colors
          , xlabel = "Infection severity"
          , ylabel = "% of cases\n in 1000 simulations")

    severity_median_percentage = DataFrame( infection_severity = infection_severity
                                           , median_percentage_across_simreps = [
                                                                                  median(inf_severity_percs_separated[1])*100
                                                                                , median(inf_severity_percs_separated[2])*100
                                                                                , median(inf_severity_percs_separated[3])*100
                                                                                , median(inf_severity_percs_separated[4])*100
                                                                                , median(inf_severity_percs_separated[5])*100
                                                                               ]
                                           )
    println(severity_median_percentage)
    return(severity_median_percentage)
end

# Save figure
#savefig("scripts/primary_care_v_icu/infection_severity_weights.png")


#violin(total_infections, legend=false, palette = colors
#        , ylabel = "total infections per simulation"
#        , yscale = :log10)


"""
Function: inf_age_estimate

Description:    This function takes a .jld2 file containing a single simulation or multiple simulation replicates and
                plots the age distribution disaggregated by ICU and GP care stages.

Arguments       'sims_file':        .jld2 file containing simulation replicates as output by {simtree} or {simforest}
                                    For example, "sims", or "sims_G_gp_filter" or "sims_G_icu_filter" if pre-filtered for certain cases.
                
                'sim_object_name':  The name of the object saved to the .jld2 file. 
                                    For example, "sims", "sims_G_icu_filter", or "sims_G_gp_filter" for files filtered using {sims_filter}.
                                    Without this name the data cannot be loaded. 
                
                'sims_file_type::String':   Options are either "full" output from {simtree} or {simforest} or filtered for the G dataframe, "filtered_G". 
    
Returns   Two histogram plots on screen, which can then be saved using the command in the example below
    
Example
    # Estimate infection age distribution disaggregated by ICU and GP
    inf_age_estimate( sims_file = "covidlike-1.3.1-sims-nrep1000_955898_2.jld2"
                        , sims_object_name = "sims" 
                        , sims_file_type = "filtered_G"
                        )
    # Save figure
    savefig("examples/age_distribution_icu_gp.png")


    gp_icu_age_groups = inf_age_estimate( sims_file = "C:/Users/kdrake/OneDrive - Imperial College London/Documents/mSCAPE/3_results/from_hpc/2025_10_maxtime_maxcases/1101621/covidlike-1.3.14-sims_filtered_G_icu_gp_nrep1200_1101621.jld2"
                                        , sims_object_name = "sims")
    println(gp_icu_age_groups)
    CSV.write("age_groups_gp.csv", gp_icu_age_groups[1])
    CSV.write("age_groups_icu.csv", gp_icu_age_groups[2])

"""

function inf_age_estimate(; sims_file
                            , sims_object_name = "sims"
                            , sims_file_type = "full"
                            )

    @assert sims_file_type in ["full","filtered_G"]

    # Load simulation file
    sims = load(sims_file, sims_object_name)
    
    # Initialise vectors to store results
    infectee_age_gp = [] ; infectee_age_icu = []

    # Loop through the simulation replicates and gather infectee ages by care pathway (note that individuals can visit a GP AND ICU but might be admitted to ICU without GP visit)
    for s in sims #s=sims[1]

        if sims_file_type == "filtered_G"
            age_gp = s[ isfinite.(s.tgp), :infectee_age ]
            infectee_age_gp = vcat(infectee_age_gp, age_gp)

            age_icu = s[ isfinite.(s.ticu), :infectee_age ]
            infectee_age_icu = vcat(infectee_age_icu, age_icu)
        
        elseif sims_file_type == "full"
            
            age_gp = s.G[ isfinite.(s.G.tgp), :infectee_age ]
            infectee_age_gp = vcat(infectee_age_gp, age_gp)

            age_icu = s.G[ isfinite.(s.G.ticu), :infectee_age ]
            infectee_age_icu = vcat(infectee_age_icu, age_icu)           
        end
    end

    # Plot infectee age distributions in a histogram
    plot1 = Plots.histogram(infectee_age_gp
                #, bins=30
                , title=""
                , xlabel="Infectee age", ylabel="Number of infected individuals ($(length(sims)) sim reps)"
                #, normalize = :pdf
                , alpha = 0.3 
                , label = "GP"
                )
    Plots.histogram!(infectee_age_icu
                #, bins=30
                #, title=""
                #, xlabel="Infectee age", ylabel="Proportion"
                #, normalize = :pdf
                , alpha = 0.3 
                , label = "ICU"
                )

    plot2 = Plots.histogram(infectee_age_icu
                , bins=21
                , title=""
                , xlabel="Infectee age", ylabel="Proportion"
                , normalize = :pdf
                , alpha = 0.3 
                , label = "ICU"
                )
    Plots.histogram!(infectee_age_gp
                , bins=21
                #, title=""
                #, xlabel="Infectee age", ylabel="Proportion"
                , normalize = :pdf
                , alpha = 0.3 
                , label = "GP"
                )
    Plots.plot(plot1, plot2, layout = (1, 2), legend = false)
    
    # Group ages and return in df
    age_groups_gp = group_ages(ages = infectee_age_gp,  width=5, start=0, top=100)
    age_groups_icu = group_ages(ages = infectee_age_icu,  width=5, start=0, top=100)
    
    return( [age_groups_gp, age_groups_icu] )
end



"""
Function        group_ages(; ages, width=5, start=0, top=100)

Description     Transform a vector of ages into a table (DataFrame) that counts people in
                `width`-year age groups, starting at `start`, with an open-ended top bin `top+`.

Arguments
                - ages::AbstractVector         : numeric ages (Int or Real); `missing` allowed
                - width::Integer = 5           : bin width (e.g., 5 years)
                - start::Integer = 0           : lower bound of the first bin (e.g., 0)
                - top::Integer = 100           : open-ended top bin threshold (e.g., 100+)

Returns
                - DataFrame with columns: `AgeGroup` (String), `Count` (Int)

Example usage
                ages = [12, 18, 23, 27, 34, 36, 42, 45, 49, 52, 58, 61, 67, 72, 75, 80, 100, 101, 105]#, missing, 4.9]
                df = group_ages(ages; width=5, start=0, top=100)
                println(df)

"""
function group_ages(; ages::Vector, width::Integer=5, start::Integer=0, top::Integer=100)
    # Filter out missings and invalid ages
    #clean = [a for a in ages if !ismissing(a) && isnumeric(a) && a ≥ start]

    # Floor non-integer ages to their integer year (e.g., 42.7 -> 42)
    clean_int = floor.(Int, ages)#clean)

    # Build labels for finite bins up to (top-1)
    max_finite_upper = top - 1
    finite_edges = collect(start:width:max_finite_upper)  # lower bounds
    labels = ["$(lb)-$(min(lb + width - 1, max_finite_upper))" for lb in finite_edges]

    # Add open-ended top bin label
    top_label = "$(top)+"
    push!(labels, top_label)

    # Initialize counts
    counts = Dict(label => 0 for label in labels)

    # Count ages into bins
    for age in clean_int
        if age ≥ top
            counts[top_label] += 1
        else
            # Compute the finite bin index
            # Example: start=0, width=5 -> 0-4, 5-9, ...
            bin_index = (age - start) ÷ width
            lb = start + bin_index * width
            ub = min(lb + width - 1, max_finite_upper)
            label = "$(lb)-$(ub)"
            counts[label] += 1
        end
    end

    # Build DataFrame in the desired order
    DataFrame(
        AgeGroup = labels,
        Count = [counts[lbl] for lbl in labels],
    )
end