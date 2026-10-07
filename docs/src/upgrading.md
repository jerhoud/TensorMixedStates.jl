# Upgrading from 1.x

Version 2.0 renames and removes a few names, so that the whole interface follows the same
rules, and changes some defaults. A program written for 1.x may need the edits below, most of
them mechanical: a name it still uses is refused with an `UndefVarError`, a `MethodError` or
an error on an unsupported keyword, rather than being taken for something else.

The files written by `save_state` or `SaveState` with 1.x are read by 2.0. The checkpoints are
not: a run interrupted under 1.x cannot be resumed under 2.0. Finish it before upgrading, or
start it again with `runTMS(sim_data; restart = true)`.

## Renamed

| 1.x | 2.0 |
|:----|:----|
| `measures`, `final_measures`, `measures_period` (fields of the phases) | `measurements`, `final_measurements`, `measurements_period` |
| `final_measures` (field of `SimData`) | `final_measurements` |
| `tolerance` (`GroundState`, `SteadyState`) | `tol` |
| `n_expand`, `n_hermitianize` (`Tdvp`, `ApproxW`, `tdvp`, `approx_W`, `thermal_state`) | `expand_period`, `hermitianize_period` |
| `tdvp(...; nsweeps, first_sweep)`, `approx_W(...; nsweeps, first_sweep)` | `nsteps`, `first_step` |
| `variance(h, state)` | `variance(state, h)` |
| `partial_trace(state, positions; keepers = true)` | `partial_trace(state, positions; keep = true)` |
| `PartialTrace(trace_positions = p)` | `PartialTrace(positions = p)` |
| `PartialTrace(keep_positions = p)` | `PartialTrace(positions = p, keep = true)` |
| `has_fermionic` | `hasfermionic` |
| `log_msg` | `log_message` |

`dmrg`, `steady_state`, `GroundState` and `SteadyState` keep `nsweeps`: they count sweeps,
where `tdvp` and `approx_W` count time steps. The phase `Evolve` is given its `duration` and
its `time_step` as before.

For instance, a phase written for 1.x as

```julia
Evolve(duration = 2., time_step = 0.1, algo = Tdvp(n_expand = 5), evolver = -im * H,
       measures = "data" => [X, Z], measures_period = 2)
```

is written for 2.0 as

```julia
Evolve(duration = 2., time_step = 0.1, algo = Tdvp(expand_period = 5), evolver = -im * H,
       measurements = "data" => [X, Z], measurements_period = 2)
```

## Removed

These names were deprecated since 1.2.0 or 1.3.0:

| 1.x | 2.0 |
|:----|:----|
| `EE` | `EntanglementEntropy` |
| `Linkdim` | `MaxLinkdim` |
| `Mutual_Info_Renyi2` | `MutualInfoRenyi2` |
| `DataToFrame` | `data_to_frame` |
| `Dmrg` | `GroundState` |
| `steady_state(...; alg)` | `steady_state(...; mpo_algo)` |

`tdvp`, `approx_W`, `dmrg` and `steady_state` no longer pass the options they do not know on
to ITensorMPS: `cutoff`, `maxdim` and `mindim` are given through
`limits = Limits(cutoff = ..., maxdim = ..., mindim = ...)`.

## Phases of your own

A phase of your own is now a subtype of [`TensorMixedStates.AbstractPhase`](@ref), and its
field of final measurements is `final_measurements`:

```julia
using TensorMixedStates: AbstractPhase

Base.@kwdef struct Kicks <: AbstractPhase
    name::String = "kicks"
    time_start = nothing
    final_measurements = []
    nkicks::Int = 4
end
```

`SimData` refuses a struct that is not a subtype of `AbstractPhase`. The union `Phases` of the
phase types of the library is gone, `AbstractPhase` taking its place. See
[Extending TMS](@ref) for what a phase of your own can do in 2.0.

`run_steps`, `resume_step` and `Algo` are no longer exported, the interfaces of one's own
being experimental: import them, as with `using TensorMixedStates: run_steps`, or write them
in full.

## Defaults that changed

Some results change without any edit of the program:

- **The default truncation.** `Limits()` has the cutoff `eps()`, about 2.2e-16, the default
  of the solvers of ITensorMPS, rather than 0. It discards the singular values below about
  1.5e-8 of the norm, which removes the noise of rounding that made the bond dimension grow;
  the bond dimensions of a run with the default limits are smaller. `Limits(cutoff = 0)`
  keeps everything. A serious simulation chooses its own limits in any case.
- **`steady_state`** and `SteadyState` solve each local step with `Krylov(dim = 8,
  maxiter = 3)` by default, which reaches steady states the former default stalled before:
  the results of a steady state search change, for the better.
- **`tdvp`** tests the convergence of each local exponentiation after every Krylov vector: it
  is several times faster, the states agreeing with those of 1.x to about 1e-13.
- **The MPOs** of the operators are compacted, with fewer channels for the same operator, and
  the approximations WI and WII of `approx_W` change at the order of their error.
- **The time functions** of a time dependent evolver take real values: a complex function is
  written as two terms, its real and its imaginary parts.
- **`Gate`** refuses an operator already placed on sites: `Gate(X)(1)` rather than
  `Gate(X(1))`.
- **`data_to_frame`** orders its columns by the names of the measurements.

The changelog of the package lists every change, the additions of 2.0 included.
