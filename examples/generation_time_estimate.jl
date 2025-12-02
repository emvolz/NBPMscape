#= 
- Version description
- Simulate multiple replicates using default covid-like parameters
- Includes age diaggregation of:
                                - contact numbers and assortativity
                                - commuting (only people aged 16 years and over)
                                - infection severity
- Plus more detailed care pathways and rates
- Plus death COMPARTMENT
- Plus recording whether infection was imported
- Change to icurate to ONLY include moving from :admittedhospital to :admittedicu
- gprate changed from 1/3 to 1/5
- changes to carestage and infstages in jumps and trecovered now defined using distribution and jumps
- maxtime increased from 60 to 90
- Prep for sub-sampling
- Serialise 
=#

using Pkg
using NBPMscape 
using JLD2 

# Run simulation
NREPS = 10 #10 #1000
sims = map( _-> simforest(NBPMscape.P; initialtime=0.0, maxtime=90.0, max_cases = 10000000, maxgenerations=100), 1:NREPS )
@save "covidlike-1.3.6-sims-nreps10.jld2" sims

# Reload if necessary
sims = load("covidlike-1.3.6-sims-nreps10.jld2", "sims")

# View transmission pair information
sims[1].D
CSV.write("covidlike-1.3.6-sims-nreps10_D.csv", sims[1].D)

# Compute generation times
# For single simulation
Tg_results = generation_time(tinf_df = sims[3].D[:,[:donor,:recipient,:timetransmission]])

# For multiple simulations
nreps = length(sims)                           
Tg_results_df = DataFrame(
    sim_n = Vector{Int}(undef, nreps),
    Gamma_fit_median = Vector{Float64}(undef, nreps),
    Gamma_fit_mean = Vector{Float64}(undef, nreps)
)

for i in 1:nreps #i=1
    Tg_results = generation_time(tinf_df = sims[i].D[:,[:donor,:recipient,:timetransmission]])
    Tg_results_df[i,:sim_n] = i
    Tg_results_df[i,:Gamma_fit_median] = only(filter( row -> (row[:statistic] == "median"), Tg_results)[:,:Gamma_fit])
    Tg_results_df[i,:Gamma_fit_mean] = only(filter( row -> (row[:statistic] == "mean"), Tg_results)[:,:Gamma_fit])
end