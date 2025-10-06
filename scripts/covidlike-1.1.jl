#= 
- Simulate multiple replicates using default covid-like parameters 
- Plus age disaggregation of contact numbers and assortativity
- Plus age disaggregation of commuting (only people aged 16 years and above)
- Prep for sub-sampling
- Serialise 
=#

using Pkg
Pkg.instantiate()
Pkg.resolve()

using NBPMscape 
using JLD2 

NREPS = 1000
sims = map( _-> simforest(NBPMscape.P; initialtime=0.0, maxtime=60.0, maxgenerations=100), 1:NREPS )
@save "covidlike-1.1.1-sims.jld2" sims 


