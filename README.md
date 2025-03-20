# NBPMscape

A Julia package for simulating disease transmission and cluster growth.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/emvolz/NBPMscape.jl")
```

## Usage

```julia
using NBPMscape

# Run a simulation
results = simtree(NBPMscape.P; 
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

## Parameters

The package includes default parameters in `NBPMscape.P`. You can modify these parameters:

```julia
my_params = merge(NBPMscape.P, (μ = 0.002, ω = 0.7))
results = simtree(my_params)
```

## License

MIT
