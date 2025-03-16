# Examples

## Basic Simulation

```julia
using NBPMscape
using Plots

# Run a simulation using default parameters
results = simbp(NBPMscape.P; 
    initialtime=2000.0, 
    maxtime=2010.0, 
    maxgenerations=10, 
    initialcontact=:G
)

# Print summary statistics
println("Total infections: ", nrow(results.G))
println("Transmission events: ", nrow(results.D))

# Plot the infection timeline
scatter(results.G.timeinfected, results.G.generation, 
    xlabel="Year", ylabel="Generation", 
    title="Infection Timeline", legend=false)
```

## Custom Parameters

You can modify the default parameters:

```julia
# Create modified parameters
my_params = merge(NBPMscape.P, (
    τ₀ = 0.002,  # Higher initial transmission probability
    μ = 0.002,   # Faster molecular clock
    ω = 0.7      # Higher variance for genetic distances
))

# Run simulation with custom parameters
results = simbp(my_params; 
    initialtime=2000.0, 
    maxtime=2020.0, 
    maxgenerations=15
)
```

## Analyzing Results

```julia
using NBPMscape
using DataFrames
using Plots

# Run simulation
results = simbp(NBPMscape.P)

# Get information about infections
infections_df = results.G

# Calculate statistics by generation
gen_stats = combine(groupby(infections_df, :generation),
    nrow => :count,
    :timeinfected => mean => :mean_infection_time
)

# Plot number of infections by generation
bar(gen_stats.generation, gen_stats.count,
    xlabel="Generation", ylabel="Number of infections",
    title="Infections by generation", legend=false)
```
