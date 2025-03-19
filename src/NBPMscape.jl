module NBPMscape

using JumpProcesses 
using DifferentialEquations
using Random
using Distributions
using DataFrames
import UUIDs 
import StatsBase 
using Interpolations
import SpecialFunctions as SF 
using Plots 

using RData 
using CSV 
# using Pkg.Artifacts

using Revise
using Debugger 

const COMMUTERPROBPATH = joinpath(@__DIR__, "..", "data",  "commuting_ITL2_prob_list.rds")
const COMMUTEPROB = load( COMMUTERPROBPATH )   
const COMMUTERINPROBPATH = joinpath(@__DIR__, "..", "data",  "commuting_ITL2_inprob_list.rds")
const COMMUTEINPROB = load( COMMUTERINPROBPATH )  
const COMMUTERMPATH = joinpath(@__DIR__, "..", "data",  "commuting_ITL2_list.rds")
const COMMUTERM = load( COMMUTERMPATH )  
const REGKEYPATH = joinpath( @__DIR__, "..", "data", "ITL2_key2.rds" )
const REGKEY = load( REGKEYPATH )

export simbp, simgeneration, simgendist, Infection
export transmissionrate, sampdegree #TODO 
include("core.jl")

end 
