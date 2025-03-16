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
using Debugger 

export simbp, simgeneration, simgendist, Infection

include("core.jl")

end 
