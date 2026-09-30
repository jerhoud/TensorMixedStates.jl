# High level interface

A simulation is described by a `SimData`, which lists its phases, and run by `runTMS`, which
returns the `Simulation` it ends with. The [Manual](@ref) shows it at work.

## Running a simulation

```@docs
runTMS
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

## Phases of one's own

```@docs
TensorMixedStates.run_phase
run_steps
resume_step
close_sim_files
```
