# High level interface

A simulation is described by a `SimData`, which lists its phases, and run by `runTMS`, which
returns the `Simulation` it ends with. The [Manual](@ref) shows it at work.

## Running a simulation

```@docs
runTMS
stopped
SimData
Simulation
get_sim_file
Data
data_to_frame
```

## Phases

```@docs
AbstractPhase
CreateState
LoadState
SaveState
ToMixed
Evolve
Algo
Tdvp
ApproxW
Gates
GroundState
SteadyState
Thermalize
PartialTrace
Weaken
```

## Phases, algorithms and representations of one's own

How to use them is explained in [Extending TMS](@ref).

```@docs
TensorMixedStates.run_phase
TensorMixedStates.creates_state
TensorMixedStates.phase_system
run_steps
resume_step
resume_time
committed_time
close_sim_files
TensorMixedStates.evolve
TensorMixedStates.write_state
TensorMixedStates.read_state
```
