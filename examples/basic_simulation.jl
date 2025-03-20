using NBPMscape
using Plots
using DataFrames

# Run a simulation using default parameters
results = simtree(NBPMscape.P
	, initialtime=2000.0
	, maxtime=2010.0
	, maxgenerations=10
	, initialcontact=:G
)

# Print summary statistics
println("Total infections: ", nrow(results.G))
println("Transmission events: ", nrow(results.D))

# Plot the infection timeline
scatter(results.G.timeinfected, results.G.generation
	, xlabel="Year"
	, ylabel="Generation"
	, title="Infection Timeline"
	, legend=false
)

# Optional: Save the plot
# savefig("infection_timeline.png")

# Access specific infection information
if nrow(results.G) > 0
	println("\nSample of infected individuals:")
	for i in 1:min(5, nrow(results.G))
		println("ID: ", results.G.pid[i]
			, ", Generation: ", results.G.generation[i]
			, ", Infected: ", results.G.timeinfected[i]
			, ", Diagnosed: ", results.G.timediagnosed[i]
		)
	end
end
