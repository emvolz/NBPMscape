# Check if simulations are independent
# If they are then they can be summed to created simulations with larger numbers of replicates
# Depends on whether they are initiated with a seed

using JLD2

# Load simulations
sims_nrep1000 = load("covidlike-1.0-sims-regionentry.jld2" , "sims") # sim with 1000 replicates. preliminary simulation re-run with regionentry adjusted
sims_nrep5000 = load("covidlike-1.0-sims-regionentry-nrep5000.jld2" , "sims") # sim with 5000 replicates. preliminary simulation re-run with regionentry adjusted

sims_nrep1000[1:1000] == sims_nrep5000[1:1000]
sims_nrep1000[1].G
sims_nrep5000[1].G

# Conclusion - the different simulation runs look different and so could be summed
