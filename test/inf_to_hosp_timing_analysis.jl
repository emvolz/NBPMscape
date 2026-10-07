# Script to check correspondence between time delay inputs (implemented via jump process rates)
# and the resulting time delays in the outbreak simulation

#------------------------------------------------------------------------------------
# Initial draft created using Claude Sonnet 5.0 via dAIsy on 7 Oct 2026
# with prompt:
#    Please build out this pseudocode outline into full Julia code:
#    Create object to store data
#
#   for i in 1:nsims
#
#   filter for severity == moderate_ED
#
#   t_inf_to_gp_moderate = tgp .- tinf
#   t_gp_to_hosp_moderate = ted .- tgp
#   t_inf_to_hosp_moderate = ted .- tinf
#
#   filter for severity in [severe, very_severe]
#
#   t_inf_to_hosp_sev_v_sev = thosp .- tinf
#   t_inf_to_gp_sev_v_sev = tgp .- tinf
#   t_gp_to_hosp_sev_v_sev = thosp .- tgp
#
#   # Add values to data storage object
#
#   # Plot distributions for each variable and label with the mean and median values
#------------------------------------------------------------------------------------

using DataFrames
using Random
using Statistics
using Plots
using Distributions
using JLD2
using NBPMscape

# =================================================================
# 1. Load simulation data
# =================================================================
root_dir = root_folder() # Function defined in NBPMscape/src/misc_functions
sims = load( joinpath( root_dir, "scripts/paper/1_outbreak_simulations/covid_like/1717144_1719024_1719029_analysis/sim_files/covidlike-1.5.0-sims-filtered_G_nrep1000_1717144_1719024_1719029.jld2" )
            , "sims")

            nsims = length(sims)   # number of simulation replicates
unique(sims[1].severity)
#5-element Vector{Symbol}:
 #:moderate_GP
 #:verysevere
 #:moderate_ED
 #:severe_hosp_long_stay
 #:severe_hosp_short_stay


# =================================================================
# 2. Object to store data
#    A Dict of vectors (one per timing variable), pooled across
#    all simulations.
# =================================================================
results = Dict{Symbol, Vector{Float64}}(
    :t_inf_to_gp_moderate    => Float64[],
    :t_gp_to_hosp_moderate   => Float64[],
    :t_inf_to_hosp_moderate  => Float64[],
    :t_inf_to_hosp_sev_v_sev => Float64[],
    :t_inf_to_gp_sev_v_sev   => Float64[],
    :t_gp_to_hosp_sev_v_sev  => Float64[]
)

# =================================================================
# 3. Helper to strip Inf / -Inf values.
#    (Default values for times (e.g. tgp, ted etc) is Inf,
#    so if the infected individual has not reach that healthcare 
#    stage when the simulation ends then the vlaue will be Inf)
# =================================================================
function remove_inf(x::AbstractVector{<:Real})
    n_before = length(x)
    cleaned  = filter(v -> !isinf(v), x)
    #n_removed = n_before - length(cleaned)
    #if n_removed > 0
    #    println("  Removed $n_removed Inf/-Inf value(s)")
    #end
    return cleaned
end

# =================================================================
# 4. Main simulation loop
#       :tinf     -> time of infection
#       :tgp      -> time of GP presentation
#       :ted      -> time of ED presentation
#       :thospital-> time of hospital admission (note that for severe and very severe it is assumed that admission is via ED but no time is spent in ED and a ted is not recorded)
# =================================================================
for i in 1:nsims

    # --- Get this simulation's individual-level data ---
    sim_data = sims[i] #simulate_one_run(n_people_per_sim)

    # -------------------------------------------------------
    # Filter for severity == moderate_ED
    # -------------------------------------------------------
    df_mod = filter(row -> row.severity == :moderate_ED, sim_data)
    
    t_inf_to_gp_moderate   = df_mod.tgp .- df_mod.tinf
    t_gp_to_hosp_moderate  = df_mod.ted .- df_mod.tgp
    t_inf_to_hosp_moderate = df_mod.ted .- df_mod.tinf

    # -------------------------------------------------------
    # Filter for severity in [severe, very_severe]
    # -------------------------------------------------------
    df_sev = filter(row -> row.severity in (:severe_hosp_short_stay,:severe_hosp_long_stay,:very_severe), sim_data)

    t_inf_to_hosp_sev_v_sev = df_sev.thospital .- df_sev.tinf
    t_inf_to_gp_sev_v_sev   = df_sev.tgp       .- df_sev.tinf
    t_gp_to_hosp_sev_v_sev  = df_sev.thospital .- df_sev.tgp

    # -------------------------------------------------------
    # Remove Inf / -Inf values before storing
    # -------------------------------------------------------
    t_inf_to_gp_moderate    = remove_inf(t_inf_to_gp_moderate)
    t_gp_to_hosp_moderate   = remove_inf(t_gp_to_hosp_moderate)
    t_inf_to_hosp_moderate  = remove_inf(t_inf_to_hosp_moderate)
    t_inf_to_hosp_sev_v_sev = remove_inf(t_inf_to_hosp_sev_v_sev)
    t_inf_to_gp_sev_v_sev   = remove_inf(t_inf_to_gp_sev_v_sev)
    t_gp_to_hosp_sev_v_sev  = remove_inf(t_gp_to_hosp_sev_v_sev)

    # -------------------------------------------------------
    # Add values to data storage object
    # -------------------------------------------------------
    append!(results[:t_inf_to_gp_moderate],    t_inf_to_gp_moderate)
    append!(results[:t_gp_to_hosp_moderate],   t_gp_to_hosp_moderate)
    append!(results[:t_inf_to_hosp_moderate],  t_inf_to_hosp_moderate)
    append!(results[:t_inf_to_hosp_sev_v_sev], t_inf_to_hosp_sev_v_sev)
    append!(results[:t_inf_to_gp_sev_v_sev],   t_inf_to_gp_sev_v_sev)
    append!(results[:t_gp_to_hosp_sev_v_sev],  t_gp_to_hosp_sev_v_sev)
end

# =================================================================
# 5. Time delays input in configuration file that generate the data
#    in the simulation via jump processes (core.jl)
#    (time delays are actually converted to rates = 1 / time delay)
#    These are overlaid on each plot for visual comparison against the
#    simulated distribution.
# =================================================================
input_values = Dict(
    :t_inf_to_gp_moderate    => 3.0,
    :t_gp_to_hosp_moderate   => 2.0,
    :t_inf_to_hosp_moderate  => 5.0,
    :t_inf_to_hosp_sev_v_sev => 4.0,
    :t_inf_to_gp_sev_v_sev   => 3.0,
    :t_gp_to_hosp_sev_v_sev  => 1.0
)

# =================================================================
# 5. Plot distributions for each variable, labeled with mean/median
# =================================================================
function plot_with_stats(data::Vector{Float64}, title_str::String, input_value::Real;
                          bins::Int = 30, color = :blue)#:steelblue)

    # Safety net: ensure no Inf/-Inf values slip through to plotting
    data = remove_inf(data)

    mu  = mean(data)
    med = median(data)

    p = histogram(data, bins = bins, normalize = :pdf, color = color, alpha = 0.6,
                  label = "Distribution", legend = :topright,
                  xlabel = "Time (days)", ylabel = "Density",
                  title = title_str, titlefontsize = 10
                  , xlim = [0,30])

    vline!(p, [mu],  color = :red,   linewidth = 2, linestyle = :dash,
           label = "Mean = $(round(mu, digits = 2))")
    vline!(p, [med], color = :black, linewidth = 2, linestyle = :dot,
           label = "Median = $(round(med, digits = 2))")
    vline!(p, [input_value], color = :green, linewidth = 2, linestyle = :solid,
           label = "Input value = $(round(input_value, digits = 2))")

    return p
end

plot_titles = Dict(
    :t_inf_to_gp_moderate    => "Infection → GP (Moderate/ED)",
    :t_gp_to_hosp_moderate   => "GP → ED (Moderate/ED)",
    :t_inf_to_hosp_moderate  => "Infection → ED (Moderate/ED)",
    :t_inf_to_hosp_sev_v_sev => "Infection → Hospital (Severe/Very Severe)",
    :t_inf_to_gp_sev_v_sev   => "Infection → GP (Severe/Very Severe)",
    :t_gp_to_hosp_sev_v_sev  => "GP → Hospital (Severe/Very Severe)"
)

# Preserve a sensible, consistent ordering for the plot grid
ordered_keys = [
    :t_inf_to_gp_moderate,
    :t_gp_to_hosp_moderate,
    :t_inf_to_hosp_moderate,
    :t_inf_to_gp_sev_v_sev,
    :t_gp_to_hosp_sev_v_sev,
    :t_inf_to_hosp_sev_v_sev
]

plots_list = [plot_with_stats(results[key], plot_titles[key], input_values[key]) for key in ordered_keys]

final_plot = plot(plots_list..., layout = (2, 3), size = (1200, 1000))
display(final_plot)
savefig(final_plot, "test/inf_to_hosp_timing_distributions.png")

println("Done. Plot saved to inf_to_hosp_timing_distributions.png")
