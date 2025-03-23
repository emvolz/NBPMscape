#= 
- Simulate multiple replicates using default covid-like parameters 
- Prep for sub-sampling
- Serialise 
=#

using NBPMscape 
using JLD2 

NREPS = 1000
sims = map( _-> simforest(NBPMscape.P; initialtime=0.0, maxtime=60.0, maxgenerations=100), 1:NREPS )
@save "covidlike-1.0-sims.jld2" sims 


