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
DataToFrame
```

## Phases

```@docs
Phases
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
Dmrg
SteadyState
PartialTrace
Weaken
```

## Phases, algorithms and representations of one's own

How to use them is explained in [Extending TMS](@ref).

```@docs
TensorMixedStates.run_phase
run_steps
resume_step
close_sim_files
TensorMixedStates.evolve
TensorMixedStates.write_state
TensorMixedStates.read_state
```
