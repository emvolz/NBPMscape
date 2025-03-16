# NBPMscape.jl

```@meta
CurrentModule = NBPMscape
```

Documentation for [NBPMscape](https://github.com/emvolz/NBPMscape.jl).

NBPMscape is a Julia package for simulating disease transmission and cluster growth.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/emvolz/NBPMscape.jl")
```

## Quick Start

```julia
using NBPMscape

# Run a simulation
results = simbp(NBPMscape.P; 
    initialtime=1990.0, 
    maxtime=2020.0, 
    maxgenerations=100, 
    initialcontact=:G
)

# Access simulation results
results.G  # DataFrame with infection information
results.D  # DataFrame with transmission information
results.infections  # Array of Infection objects
results.H  # Transmission history
```
