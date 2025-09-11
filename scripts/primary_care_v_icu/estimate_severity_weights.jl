#### Counting severity to estimate % in each category

using JLD2
using DataFrames
using StatsPlots


# Load simulation file
sims = load("C:/Users/kdrake/OneDrive - Imperial College London/Documents/mSCAPE/3_results/from_hpc/primary_care_Sep_2025/955898/covidlike-1.3.1-sims-nrep1000_955898_2.jld2", "sims")

# Column to search in
column_name = :severity 
# Possible values for infection severity
infection_severity = [:asymptomatic, :mild, :moderate, :severe, :verysevere ]  

# Initialise vectors to store results
inf_severity_counts = []
inf_severity_percs = []
total_infections = []

# Loop through the simulation and gather frequency and percentages for each severity type
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
infection_severity_names = ["asymptomatic" "mild" "moderate" "severe" "verysevere"]

# Create the violin plot
colors = [:red, :blue, :green, :orange, :purple]
violin(infection_severity_names, inf_severity_percs_separated*100, legend=false, palette = colors, ylabel = "% of cases\n in 1000 simulations")

median(inf_severity_percs_separated[1])*100
median(inf_severity_percs_separated[2])*100
median(inf_severity_percs_separated[3])*100
median(inf_severity_percs_separated[4])*100
median(inf_severity_percs_separated[5])*100


# Save figure
savefig("scripts/primary_care_v_icu/infection_severity_weights.png")


violin(total_infections, legend=false, palette = colors
        , ylabel = "total infections per simulation"
        , yscale = :log10)