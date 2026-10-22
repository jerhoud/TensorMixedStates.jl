# Upgrading from 1.x

TensorMixedStates 2.0 is a rewrite of the package. What it computes is the same, but its
interface was made smaller and uniform, so that a program written for 1.x needs edits, most of
them mechanical, which this page lists from version 1.6. A name or a keyword that 2.0 removed,
renamed or moved is refused with an error rather than taken for something else, and for most of
them the error Julia prints ends with a line beginning with `TensorMixedStates 2.0:` that says
what to write instead.

## Saved states and checkpoints

The states saved by `save_state` or `SaveState` with 1.x are read by 2.0. The checkpoints are
not: 2.0 refuses to resume a simulation stopped under 1.x. Finish it with 1.x, or start it over,
by putting a file named `restart` in its directory or by calling `cleanTMS(sim_data)` before
`runTMS`.

## A program converted

These lines of the manual of 1.6 compute a ground state, three evolutions and a steady state:

```julia
energy, ground = dmrg(hamiltonian, RandomState{Pure}(mysystem, 8); nsweeps = 10,
                      limits = Limits(maxdim = 20))
myevolved = tdvp(-im * sum(Z(i) for i in 1:4), 0.5, State{Pure}(System(4, Qubit()), "+");
                 nsweeps = 10)
mydecayed = tdvp(sum(Dissipator(Sm)(i) for i in 1:4), 1.0, myrho; nsweeps = 10)
mydecayed = approx_W(sum(Dissipator(Sm)(i) for i in 1:4), 1.0, myrho; order = 4, w = 2,
                     nsweeps = 10)
value, mysteady = steady_state(sum(Dissipator(Sm)(i) + 0.5Dissipator(Sp)(i) for i in 1:4),
                               myrho; nsweeps = 10)
```

With 2.0, they are written

```julia
energy, ground = ground_state(hamiltonian, random_state(Pure, mysystem, 8);
                              algo = Dmrg(nsweeps = 10), limits = Limits(maxdim = 20))
myevolved = evolve(-im * sum(Z(i) for i in 1:4), 0.5, State{Pure}(System(4, Qubit()), "+");
                   nsteps = 10)
mydecayed = evolve(sum(Dissipator(Sm)(i) for i in 1:4), 1.0, myrho; nsteps = 10)
mydecayed = evolve(sum(Dissipator(Sm)(i) for i in 1:4), 1.0, myrho;
                   algo = ApproxW(order = 4, w = 2), nsteps = 10)
value, mysteady = steady_state(sum(Dissipator(Sm)(i) + 0.5Dissipator(Sp)(i) for i in 1:4),
                               myrho; algo = Dmrg(nsweeps = 10))
```

and they give the same results: the same ground state energy, -12.381489999654791, and the
same expectation values to 1e-14.

## Solvers

The solvers are four functions, `evolve`, `ground_state`, `steady_state` and `thermalize`,
each given its method as `algo`, see [Algorithms](algorithms.md). The options of a method are
the fields of `algo`, and no option is passed on to ITensorMPS any more: the truncation is given
by `limits = Limits(cutoff = …, maxdim = …, mindim = …)`, and a Krylov method by
`krylov = Krylov(…)`.

| 1.6 | 2.0 |
|:----|:----|
| `tdvp(evolver, t, state; nsweeps = n)` | `evolve(evolver, t, state; nsteps = n)`, or `time_step = dt`, one of them required where `tdvp` made one sweep by default, `algo = Tdvp()` by default |
| `approx_W(evolver, t, state; order, w, nsweeps = n)` | `evolve(evolver, t, state; algo = ApproxW(; order, w), nsteps = n)` |
| `Tdvp(n_expand, n_hermitianize)`, `ApproxW(n_hermitianize)` | `Tdvp(; expand_period, hermitianize_period)`, `ApproxW(; order, hermitianize_period)` |
| `first_sweep = k` | `step = k`, `time` the time the evolution starts at, `duration` and `nsteps` those of the steps left |
| `coefs` of an evolver depending on time, `hs => coefs` of `Evolve` | a vector of terms, `[A, B => g]` for ``A + g(t) B``, each with its real function of the time if it has one |
| `time_start` | `time` |
| `dmrg(H, state; nsweeps, noise)` | `ground_state(H, state; algo = Dmrg(; nsweeps, noise, tol))` |
| `steady_state(L, state; nsweeps, alg)` | `steady_state(L, state; algo = Dmrg(; nsweeps))` |
| `mpo_algo` of `steady_state` and `SteadyState` | none: the MPO of ``L^\dagger L`` is built from that of ``L`` with as few channels as its rank, `mpo_limits` truncating it if given |
| `RandomState{Pure}(system, d)` | `random_state(Pure, system, d)`, `random_state(Mixed, system, d)` or `random_state(state, d)` |
| `TdvpObserver`, `ApproxWObserver`, `DmrgObserver`, `observer!` | call the solver step by step and measure between the calls, or measure in a phase |
| `make_approx_W1`, `make_approx_W2` | not exported: `TensorMixedStates.make_approx_W1` |

`Dmrg` has no default number of sweeps, which `dmrg` took as 1. `ground_state` no longer takes
an MPO: give the operator, or its `PreMPO(state, operator)`, prepared once for several calls. A
mixed state is refused by `ground_state`: search the pure state and mix it. `ground_state`,
`steady_state` and `thermalize` return named tuples, `(energy, state)`, `(residual, state)` and
`(log_trace, state)`, which destructure as before.

## Measurements

| 1.6 | 2.0 |
|:----|:----|
| `measure(state, list, t)` | `measure(state, list; time = t)` |
| `StateFunc("name", f)` | `"name" => c -> f(c.state)` |
| `TimeFunc("name", f)` | `"name" => c -> f(c.time)` |
| `RealValue`, `ImaginaryValue` | `"name" => op`, `"name" => c -> imag(expect(c.state, op; complex = true))` |
| `ComplexValue` | removed: the value of an operator that is not Hermitian is complex |
| `Check("name", obs1, obs2, tol)` | `"name" => c -> Check(a, b; tol)`, the values `a` and `b` computed from `c`, as `expect(c.state, op)` |
| `trace2`, `Trace2` | `purity`, `Purity` |
| `TraceError`, `HermiticityError` | `TraceDeviation`, `HermiticityDeviation` |
| `EE(cut)` | `EntanglementEntropy(cut)` |
| `EntanglementEntropy(cut, n)`, the entropy then the spectrum | `EntanglementEntropy(cut)`, `EntanglementSpectrum(cut, n)` |
| `Linkdim` | `MaxLinkdim` |
| `SubRenyi2(k)`, `renyi2(state, k)` | `SubRenyi2(1:k)`, `renyi2(state, 1:k)` |
| `MutualInfoRenyi2(k)`, `Mutual_Info_Renyi2(k)`, `mutual_info_renyi2(state, k)` | `MutualInfoRenyi2(1:k)`, `mutual_info_renyi2(state, 1:k)` |
| `variance(h, state)` | `variance(state, h)` |
| `ee, spectrum = entanglement_entropy(state, cut)` | `entanglement_entropy(state, cut)`, `entanglement_spectrum(state, cut)` |
| `entanglement_by_sector(state, cut)[QN("N", 2)]` | `entanglement_by_sector(state, cut)[N => 2]`, the charges as `fully_mixed` takes them |
| `partial_trace(state, positions)` | `partial_trace(state; trace_out = positions)` |
| `partial_trace(state, positions; keepers = true)` | `partial_trace(state; keep = positions)` |
| `:sweep` in a list of measurements | `:step` |

The value of a Hermitian operator is real, its imaginary part dropped with a warning when it is
more than rounding, and that of an operator that is not Hermitian is complex, whatever the state:
the columns of a file do not depend on it. For the same reason, `Trace` is complex. `expect`,
`expect1` and `expect2` follow that rule, where 1.x gave the complex numbers computed, which
`complex = true` gives.

A list of measurements given to `measure` is one list: `measure(state, [Measure(…),
Measure(…)])` measures the two as one, where 1.x returned a group for each. A measurement
alone gives its value, where 1.x gave a vector of one pair: `last(only(measure(state, m)))`
is now `measure(state, m)`, and `measure(state, [m])` gives the pair.

## Simulations

| 1.6 | 2.0 |
|:----|:----|
| `measures`, `final_measures`, `measures_period` | `measurements`, `final_measurements`, `measurements_period` |
| `GroundState(nsweeps, noise, tolerance)` | `GroundState(algo = Dmrg(; nsweeps, noise, tol))`, `algo` being required |
| `SteadyState(nsweeps, tolerance)` | `SteadyState(algo = Dmrg(; nsweeps, tol))`, `algo` being required |
| `Dmrg(…)`, the phase | `GroundState(…)`: `Dmrg` is now the method of a search |
| `PartialTrace(trace_positions = p)` | `PartialTrace(trace_out = p)` |
| `PartialTrace(keep_positions = p)` | `PartialTrace(keep = p)` |
| `CreateState(seed = s)` | `SimData(seed = s)`, the draws of every phase using the generator of the run |
| `SimData(time_format, data_format)` | `SimData(; format)`, for every number written |
| `SimData(checkpoint_interval, max_time, threading)` | `runTMS(sim_data; checkpoint_interval, max_time, threading)`, or the variables `TMS_CHECKPOINT_INTERVAL`, `TMS_MAX_TIME`, `TMS_THREADING` |
| `runTMS(sim_data; restart = true)` | a file `restart` in the directory, or `cleanTMS(sim_data)` before `runTMS` |
| `runTMS(sim_data; clean = true)` | `cleanTMS(sim_data)` |
| `SetState(state)` | `Reset(state)`, a channel: at a rate `γ` in an evolver, `γ * Dissipator(Reset(state))` |
| `log_msg` | `log_message` |
| `output(sim, list; energy = e)` | `output(Simulation(sim, sim.state; energy = e), list)` |
| `get_sim_file`, `close_sim_files` | `output(sim, destination => list)`, or a `TensorMixedStates.Destination` of one's own; `runTMS` closes the files |
| `DataToFrame` | `data_to_frame` |
| `Simulation(state; …)`, a function given a `Simulation` | the functions take the state, `sim.state`; a `Simulation` is made by `runTMS` |
| `create_graph_state` | `CreateState(graph_state(graph))` |
| `CreateState(type = Pure(), …)` | `CreateState(type = Pure, …)`, the type as in `State{Pure}`, by default that of a `State` given |

The data files are tables: a line naming the columns, then a line per time measured.
`sim.data` holds the rows of each `Data` destination, a row reading a value by its name,
`row["X(1)"]`, and the rows a table with the columns of the data files, which `data_to_frame`
and `DataFrame` take.

A phase of one's own is a subtype of `TensorMixedStates.AbstractPhase`, which replaces the
union `Phases`, with the fields `name`, `time_start` and `final_measurements`, and it runs by
`run_phase(phase, sim, done)`, or step by step by `prepare_phase`, `step_count` and `run_step`,
see [Extending TMS](extending.md). `run_steps`, `resume_step` and `Algo` are gone,
`TensorMixedStates.EvolutionAlgorithm` and `TensorMixedStates.SearchAlgorithm` taking the place
of the last.

## Operators, sites and states

| 1.6 | 2.0 |
|:----|:----|
| `@def_operators(site, [type => [name = matrix, …]])` | `@def_operators T [name = matrix, …]`, `T` the site type |
| `Operator{N}(name, definition, type)` | `Operator{N}(name, definition)` |
| `Op`, `GenericOp`, `IndexedOp`, `SimpleOp` | not exported: `TensorMixedStates.GenericOp` |
| `OpType`, `plain_op`, `fermionic_op`, `selfadjoint_op`, `involution_op` | removed: the parity of an operator comes from `fermion_parities` of its site |
| `isfermionic`, `has_fermionic` | removed, as above |
| `named(…; type)`, `controlled(op; type)` | without `type` |
| `Left(A)`, `Right(A)` | the channel `Gate(A)`, ``\rho \to A \rho A^\dagger``, or the generator `Evolver(K)`, ``\rho \to K \rho + \rho K^\dagger`` |
| `identity_operator` | `Id`, and `matrix(Id, site)` for its matrix |
| `matrix("X", site)` | `matrix(X, site)` |
| `tensor(op, site…)`, the ITensor of an operator | `matrix(op, site…)`, its matrix, or `matrix(op, system)` on a whole system: `tensor` is `⊗` |
| `Fup`, `Fdn` of `Electron` and `Tj` | `parity(Nup)`, `parity(Ndn)` |
| `Zd`, `Xd`, `Hd`, `Sumd` of `Qudit` | `Z`, `X` of `Qudits`, `H`, `Sum` of `Circuits` |
| `H`, `S`, `T`, `Phase` of `Qubits` | in the module `Circuits`: `using .Circuits` |
| `Qubits.Swap`, `Qubits.controlled` | `Swap` and `controlled`, exported by TensorMixedStates for sites of any type |
| a local state given as a number `n`, as `Proj(1)` | `Level => n`, the basis state `n` counted from 0, as `Proj(Level => 1)` |
| `conserve_string` | `conserving(T; conserve)`, in the keyword constructor of a site type |
| `state(site, name)` | `TensorMixedStates.local_state(site, name)` |
| `Index(site)`, `SysIndex`, `system.pure_indices` | `siteind(Pure, system, i)`, `siteinds(Pure, system)` |
| `State(system, state)` | removed: systems with the same sites are equal |
| `load_state(file, name; system)` | `load_state(file, name)` |
| `simplify`, `compact` | removed: an operator is lowered on the system it acts on, and its MPO compacted when it is built |

## Defaults that changed

Some results change without any edit beyond those above:

- **The truncation.** A field of `Limits` left out takes the default of the function the limits
  are given to, almost everywhere the cutoff `1e-26`, which discards the singular values below
  1e-13 of the norm, what rounding leaves, where 1.x kept everything. The bond dimensions are
  smaller, the results differing at the level of rounding. `Limits(maxdim = 50)` truncates at
  `1e-26` as well, and `Limits(cutoff = 0)` keeps everything.
- **`Tdvp`** enlarges the bonds of a state whose bond dimension is below both 5 and `maxdim`
  before a step, which makes an evolution from a product state much more precise, and
  `expand_period` adds such expansions periodically. It tests the convergence of each local
  exponentiation after every Krylov vector, where 1.x waited for the Krylov space to be full.
- **`ApproxW`** multiplies by its MPOs exactly before truncating, `apply_algo = "naive"`, where
  1.x used `"densitymatrix"`, and its approximation WII, `w = 2`, is built after Zaletel et al.:
  its results change at the order of the error of the method, with smaller bond dimensions.
  With an evolver depending on time, the orders 3 and 4 use the commutator-free Magnus
  integrator of order 4 of Blanes and Moan, where 1.x read the time functions at the middle of
  each step.
- **`steady_state`** and `SteadyState` solve each local problem with `Krylov(dim = 8,
  maxiter = 3)`, which the spectrum of ``L^\dagger L``, crowded near zero, needs, where 1.x
  took those of ITensorMPS, 3 and 1: they converge where the searches of 1.x could stall.
  `Krylov(dim = 3, maxiter = 1)` with `Limits(cutoff = 0)` gives the results of 1.x again.
- **The values measured.** Whether an operator is Hermitian is decided numerically: an operator
  1.x could not prove Hermitian, as `Sp(1)Sm(2) + Sm(1)Sp(2)`, has a real value, and an
  anti-Hermitian one its complex value, where 1.x gave its imaginary part as `Im(name)`. A pair
  of operators is named `(X, Y)`, where 1.x named it `XY`.
- **The files.** Every number is written with `SimData(format = "%.8g")`, the time included,
  which 1.x wrote with four digits. A `.json` file is the array of an object per time, where
  1.x wrote an object per measurement, and a `.jsonl` file the same objects as JSON Lines,
  written as they are measured.
- **`runTMS`** no longer changes the working directory: a relative path in the code of a phase
  is taken from the directory `runTMS` was called in, where 1.x ran in that of the simulation.
- **`data_to_frame`** gives the columns of the data files, a vector or a matrix spread over
  several, as `X(1)`, `X(2)`, in the order the measurements were first made, where 1.x gave a
  column per measurement, a vector or a matrix in each cell.
- **The random draws** of a run use its own generator, seeded by `SimData(; seed)`, or by a seed
  drawn and written in `stamp`, and kept by its checkpoints.
