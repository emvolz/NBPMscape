#### Functions to estimate post simulation infection severity and age distributions
#(1) Function to estimate frequency/percentage of infection in each severity category
#(2) Age distribution

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
    savefig("scripts/primary_care_v_icu/infection_severity_weights.png")

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

This function takes a .jld2 file containing a single simulation or multiple simulation replicates and
plots the age distribution disaggregated by ICU and GP care stages.

# Arguments
    'sims_file':            .jld2 file containing simulation replicates as output by {simtree} or {simforest}
                            For example, "sims", or "sims_G_gp_filter" or "sims_G_icu_filter" if pre-filtered for certain cases.
    'sim_object_name':      The name of the object saved to the .jld2 file. 
                            For example, "sims", "sims_G_icu_filter", or "sims_G_gp_filter" for files filtered using {sims_filter}.
                            Without this name the data cannot be loaded. 
    
# Returns
    Two histogram plots on screen, which can then be saved using the command in the example below

# Example
    # Estimate infection age distribution disaggregated by ICU and GP
    inf_age_estimate( sims_file = "covidlike-1.3.1-sims-nrep1000_955898_2.jld2"
                        , sims_object_name = "sims" 
                        )
    # Save figure
    savefig("examples/age_distribution_icu_gp.png")

"""

function inf_age_estimate(;  sims_file
                                , sims_object_name = "sims"
                                )

    # Load simulation file
    sims = load(sims_file, sims_object_name)
    
    # Initialise vectors to store results
    infectee_age_gp = [] ; infectee_age_icu = []

    # Loop through the simulation replicates and gather infectee ages by care pathway (note that individuals can visit a GP AND ICU but might be admitted to ICU without GP visit)
    for s in sims #s=sims[2]
        age_gp = s[ isfinite.(s.tgp), :infectee_age ]
        infectee_age_gp = vcat(infectee_age_gp, age_gp)

        age_icu = s[ isfinite.(s.ticu), :infectee_age ]
        infectee_age_icu = vcat(infectee_age_icu, age_icu)
    end

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

end