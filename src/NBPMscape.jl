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
using LinearAlgebra

using RData 
using CSV 
# using Pkg.Artifacts

# TODO 
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
const CAAPATH = joinpath( @__DIR__, "..", "data", "CAA_pax_2024_ITL2.rds" )
const CAAIMPORTS = load(CAAPATH)
# itl2size = CSV.read(joinpath( @__DIR__, "..", "data", "itl2_regions_england.csv"), DataFrame)
# const ITL2SIZE = filter( r->r.Code in REGKEY.code, itl2size )
itl2size =  load( joinpath( @__DIR__, "..", "data", "itl2_population2022.rds" ) )
#= > head(d1) 
  ITL225CD                               ITL225NM total_population_2022
1     TLC3                            Tees Valley                688756 =#
const ITL2SIZE = filter( r->r.ITL225CD in REGKEY.code, itl2size )

export simtree, simforest, sampleforest, simgendist, Infection, infectivitytoR
export transmissionrate, sampdegree, REGKEY, COMMUTEPROB #TODO 
include("core.jl")

end 
