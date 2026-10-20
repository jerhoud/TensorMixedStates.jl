# High level interface

A simulation is described by a `SimData`, which lists its phases, and run by `runTMS`, which
returns the `Simulation` it ends with. The [High Level Interface](@ref) page shows it at work.

## Running a simulation

```@docs
runTMS
stopped
cleanTMS
SimData
Simulation
Data
data_to_frame
```

## Phases

```@docs
TensorMixedStates.AbstractPhase
CreateState
LoadState
SaveState
CopyState
Update
ToMixed
Weaken
PartialTrace
Gates
Evolve
Thermalize
GroundState
SteadyState
```

## Phases, methods, destinations and representations of one's own

How to use them is explained in [Extending TMS](@ref).

```@docs
TensorMixedStates.prepare_phase
TensorMixedStates.run_phase
TensorMixedStates.step_done
TensorMixedStates.provided_values
TensorMixedStates.evolve_step
TensorMixedStates.search
TensorMixedStates.write_row
TensorMixedStates.check_list
TensorMixedStates.resume_mark
TensorMixedStates.write_state
TensorMixedStates.read_state
TensorMixedStates.expect_value
```
