# API Reference

```@meta
CurrentModule = NBPMscape
```

## Functions

```@docs
simtree
simgeneration
simgendist
```

## Default Parameters

NBPMscape comes with a set of default parameters stored in `NBPMscape.P`. These parameters can be customized by creating a modified version:

```julia
my_params = merge(NBPMscape.P, (μ = 0.002, ω = 0.7))
```

Key parameters include:

- `τ₀`, `τ₁`: Initial and final transmission probabilities per act
- `ν`: Decay rate of transmission probability
- `δ`: Diagnosis rate
- `κ`: Rate of viral suppression after diagnosis
- `μ`: Mean clock rate for genetic distance
- `ω`: Variance inflation for genetic distance
