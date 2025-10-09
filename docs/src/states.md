# Systems and States

## Systems

```@docs
System
length(::System)
sim(::System)
⊗(::System, ::System)
```

## States

```@docs
Limits
State
length(::State)
maxlinkdim(::State)
mix
truncate(::State)
trace(::State)
trace2
norm(::State)
hermitianize
hermitianity
RandomState
partial_trace
```

## Simulations

```@docs
Simulation
get_sim_file
```