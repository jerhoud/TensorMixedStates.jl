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
mix(::State)
truncate(::State)
trace(::State)
trace2(::State)
norm(::State)
dag(::State)
hermitianize(::State)
hermitianity(::State)
random_state
partial_trace(::State, ::Vector{Int})
```

## Simulations

```@docs
Simulation
get_sim_file
```