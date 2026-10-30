# Changelog

Notable changes to TensorMixedStates are recorded here, in the format of
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

This file starts with the release in preparation. For what came before, see the
[tags](https://github.com/jerhoud/TensorMixedStates.jl/tags); version 1.2.8 is the one
archived as the [codebase release](https://doi.org/10.21468/SciPostPhysCodeb.72-r1.28) of
the reference article.

## [Unreleased]

### Added

- `apply_algo = "zipup"` for `apply` of an MPO and `ApproxW`: faster when truncating much, less accurate.
- `random_state(Pure, system, states, linkdims)`, a random pure state in the sector of `states`.
- `Proj(N(3) => 1)`, `Dephase(N(1))` and `Reset(N(3) => 1)` are placed where the operator is.
- The log notes each launch, the bond dimension of the MPOs, the checkpoints and the peak memory, and for each step the largest bond dimension and the duration.
- The `stamp` notes the machine and the versions of ITensors, ITensorMPS, NDTensors, KrylovKit and of the packages of the script.

### Changed

- TensorMixedStates needs ITensorMPS 0.3.45 and ITensors 0.8, for zipup.
- `expect`, `expect1`, `expect2` and `measure` take an operator of no definite charge by its part of no flux, without copying the state.
- A dense operator on many sites is split into terms faster: 15 s rather than 85 s on 10 qubits.

## [2.0.1] - 2026-10-28

### Fixed

- The log writes times in the format of the text tables: `0.3`, not `0.30000000000000004`.
- `mpo_limits` and `PreMPO(…; limits)` work on a system that conserves a charge, where they raised a `FluxError`.
- A truncated MPO keeps or drops channels of equal weight together, so a truncated hopping stays hermitian.
- `steady_state` warns when the state found has a purity above one, as for a combination whose traces nearly cancel.
- `Proj(v)` projects on the unit vector of `v`: `Proj([1, 1])` was twice a projector.
- `Reset` takes a vector for its unit vector and a matrix divided by its trace: `Reset([1, 1])` doubled the trace.
- The noise dissipators refuse a negative rate, `Limits` an empty vector, and `Dmrg` an empty or negative noise.
- `save_state` keeps the state it replaces until the new one is written whole, so that a write that fails leaves it.
- A checkpoint survives a power cut: its files are synced to the disk before it replaces the previous one, but on Windows.
- From Julia 1.11 on, the hints for the names of 1.x show only in a module that sees TensorMixedStates.

## [2.0.0] - 2026-10-24

A rewrite of the package: the same computations, through a smaller and uniform interface. A
program written for 1.6 needs edits, most of them mechanical, which the page "Upgrading from
1.x" of the documentation lists with the defaults that changed; an error on a name of 1.x says
what to write instead. The states saved by 1.x are read, its checkpoints are not resumed. The
interfaces of one's own, on the page "Extending TMS", are experimental and may change in a
version 2.x.

### Added

- `?X` on an operator of a site type: what it is on each site type that has it, read off their
  helps, and the modules exporting it; the help of `AbstractSite` lists the operators a site
  type gets from those it declares, as `Sm` from `Sp`.
- `thermalize` and the phase `Thermalize`, the thermal state of a Hamiltonian, or its canonical
  state from `fully_mixed` of a sector.
- Channels and their dissipators: `depolarizing_channel`, `dephasing_channel`, `reset_channel`,
  and for qubits `amplitude_damping_channel` and `thermal_relaxation_channel`; `Dephase`, a
  measurement whose result is not read.
- Prepared states: `superposition`, `mixture`, `ghz_state`, `w_state`, `dicke_state`,
  `dimer_state`, `aklt_state`, `mps_state`, `dense_state`, `fully_mixed`, `slater_state` and
  `fermi_sea`.
- `probabilities`, `sample(state, op)` and `collapse`, the measurement of an operator placed on
  sites, its result drawn or given.
- `reduced_density_matrix`, `von_neumann_entropy` and `log_negativity`, with their measurements.
- `Proj(A => λ)`, the projector on an eigenspace, and `Level => n`, the basis state `n` of any
  site.
- `"Up"` and `"Dn"`, or `"↑"` and `"↓"`, the states of a `Spin` of largest and smallest `Sz`, as
  on a qubit, and the names of a spin one half for the free spins of `aklt_state`.
- The phases `Update`, which changes the simulation by a function, and `CopyState`, which keeps
  a state under a name; the field `on` of the phases, which act on a named state, measured by a
  list `:φ => [...]` or by `Overlap(:φ, op)`.
- `cleanTMS`, `stopped`, and `Krylov`, the parameters of the Krylov method of each local step.
- The algebra of `PreMPO`s: `pre'`, `a * b`, `a ^ p`, `a + b`, `a - b`, `c * a` and `sum(v)`, each
  with as few channels as its rank.
- `limits` of `PreMPO`, and `mpo_limits` of `ground_state`, `evolve`, `thermalize` and their phases,
  which truncate the MPO in its triangular form, the one `ApproxW` needs.
- The module `Circuits`: the gates `H`, `S`, `T`, `Phase` and `Sum`, which `Qubits` and `Qudits` no
  longer export, the rotations `Rx` to `Rzz`, and `CX`, `CNOT`, `CY`, `CZ`, `CCX`, `Toffoli`.
- `tr`, the trace of a state under the name of `LinearAlgebra`.
- A system iterates over its sites, `for site in system`, and `symmetries(state)` gives what the
  system of a state conserves.
- `limits` of `mix`, which truncate the density matrix as `ToMixed` does, neither truncating with
  the default `Limits()`.
- `variance` and `Variance` of a density matrix ρ, tr(ρH†H) - |tr(ρH)|².
- `matrix(op, system)` and `Vector(ψ)`, the dense matrix of an operator on a small system and
  the vector of a pure state, to compare with an exact computation.
- `a ⊗ b`, the product of two states on the product of their systems, and `permute_sites`, which
  moves the sites of a state.
- `dag` and `'` of an operator placed on sites, and `'` of a generic one: `H = -(h + h')`.
- The rows of a `Data` read a value by its name, `row["X(1)"]`, and form a table of `Tables.jl`
  with the columns of a text table, which `DataFrame` and `CSV.write` take.
- An error on a keyword of ITensors, as `maxdim`, `cutoff` or `dt`, or on one a phase does not
  have, ends with what to write; one on an unknown state names those the site knows, one on an
  operator the site types that have it.
- A quick reference, pages for users of QuTiP or ITensors and on choosing and checking the
  parameters, and a skill for AI agents in `skills/tensormixedstates`.

### Changed

- The solvers are `ground_state`, `evolve`, `steady_state` and `thermalize`, given their method
  as `algo = Dmrg(…)`, `Tdvp(…)` or `ApproxW(…)` and their truncation as `limits`.
- A truncation left unspecified discards the singular values below 1e-13 of the norm, what
  rounding leaves, where 1.x kept everything; `Limits(cutoff = 0)` keeps everything.
- `measure`, the outputs, `expect`, `expect1` and `expect2` give the value of a Hermitian operator
  as real, and that of any other as complex, whatever the state; `complex = true` gives `expect`,
  `expect1` and `expect2` the complex numbers computed, as 1.x did.
- A system keeps prepared the operators of several terms measured twice on its states, and the
  systems it weakens to: `expect` of such an operator no longer lowers and compacts it again; a
  phase prepares its lists of measurements once, as the MPO of `Variance`.
- The phases, `SimData`, `runTMS`, the measurements and the operators take new names and forms,
  listed in "Upgrading from 1.x".
- `evolve` and `thermalize` take their number of steps, `nsteps`, or their step, `time_step` or
  `beta_step`, with no default, where 1.x made one step; `Evolve` and `Thermalize` take either.
- An evolver depending on time is a vector of terms, `[A, B => g]` for A + g(t) B, where 1.x
  paired the vectors of the terms and of their functions; its `PreMPO` keeps the functions.
- `CreateState` takes the type of its state, `type = Pure`, where 1.x took `Pure()`, by default
  that of a `State` given, as `CreateState(ψ)`, and `CreateState{Pure}(system, state)`.
- `ground_state`, `steady_state` and `thermalize` return named tuples, `(energy, state)`,
  `(residual, state)` and `(log_trace, state)`, which destructure as the tuples of 1.x did.
- `entanglement_entropy(state, cut)` and `EntanglementEntropy(cut)` give the entropy, and the
  new `entanglement_spectrum` and `EntanglementSpectrum(cut, n)` the spectrum, where 1.x gave
  both.
- A zero entanglement entropy, as that of a product state, is 0.0, where 1.x gave -0.0,
  written -0 in a file.
- A system whose conserved quantities take the four components ITensors allows, not all on each
  site conserving something, is refused, where ITensors failed during a computation.
- `SubRenyi2(positions)` became `Renyi2(positions)`, the measurement of `renyi2(state, positions)`,
  as `Renyi2` is that of `renyi2(state)`.
- `entanglement_by_sector` reads a sector by its charges as `fully_mixed` takes them,
  `sectors[N => 2]`, where 1.x indexed it by a `QN`.
- `SteadyState` provides the eigenvalue of ``L^\dagger L`` as `:residual`, where 1.x named it
  `:energy`, which its lists of measurements cannot read, it being that of another phase.
- `graph_state` takes the system first, `graph_state(system, graph; others)`, as the other
  prepared states, its sites in no edge taking `others`, where 1.x made a system of qubits.
- `partial_trace` and `PartialTrace` take the sites to keep, `keep = …`, or those to trace out,
  `trace_out = …`, where 1.x took the positions traced out unless `keepers = true`.
- `inner`, `fidelity` and `hs_fidelity` take two states of the same sites that conserve different
  things on what both conserve, as `Fidelity` and `Overlap` then do.
- A system is shown as Julia builds it, `System(100, Qubit())`, where each of its sites was
  written, and an error on two systems says how they differ.
- `Limits` refuses a `maxdim` below 1 and a negative `cutoff` or `mindim`, and the methods and the
  phases a negative period, where `maxdim = 0` failed at the first truncation.
- A site given twice is refused by the functions and measurements taking positions, and the sites
  kept by `partial_trace` and `reduced_density_matrix` are given in increasing order.
- `measure` of one measurement returns its value, where 1.x returned a vector of one pair; a
  list returns its pairs.
- A `.json` file is the array of an object per time, where 1.x wrote an object per measurement
  once the run ended, and a `.jsonl` file the same objects as JSON Lines, written as measured.
- A text table writes every number with thirteen significant digits, `SimData(format = "%.13g")`,
  where 1.x wrote the values with eight and the time with four.
- `steady_state` builds the MPO of ``L^\dagger L`` from that of ``L`` with as few channels as its
  rank, where the product of the MPOs had the square of their bond dimension, and its
  `mpo_limits` truncate it in its triangular form.
- A product of sums that would develop into more terms than its factors have, as a moment `H^4`,
  and on a density matrix a gate or a dissipator of an operator of many terms, are built from the
  MPOs of their factors, at a cost linear in the number of sites.
- `Op`, `GenericOp`, `IndexedOp` and `SimpleOp` are no longer exported.
- `Gates` is made of the steps of `apply`, a batch of gates or an MPO on some sites, between
  which a run can write a checkpoint or stop, where 1.x applied the whole product in one go.

### Removed

- `dmrg`, `tdvp`, `approx_W` and their observers: call the solver step by step and measure
  between the calls, or measure in a phase.
- `StateFunc`, `TimeFunc`, `RealValue`, `ImaginaryValue` and `ComplexValue`: a list of
  measurements takes `"name" => c -> f(c.state)`.
- `simplify`, `compact` and the operator types of `OpType`: an operator is lowered on the system
  it acts on, its parity given by its site.
- `Left` and `Right`, which `Gate(A)` and `Evolver(K)` replace.
- `get_sim_file` and `close_sim_files`: `output` writes to a destination, `runTMS` closes the
  files.
- The names deprecated since 1.2.0 and 1.3.0: `EE`, `Linkdim`, `Mutual_Info_Renyi2` and
  `DataToFrame`.
- `mpo_algo` of `steady_state` and `SteadyState`, the method of the product of the MPOs.

## [1.6.2] - 2026-10-20

### Fixed

- `tdvp`, and `Evolve` with `Tdvp()`, lost precision with the norm of the state, as for a
  density matrix near the fully mixed state on many sites: each step is now made at norm one.

## [1.6.1] - 2026-10-09

### Changed

- A float coefficient of integer value prints as that integer, `2.0X` as `2X`, which renames a
  measurement written with one: an integer coefficient is now stored as a float, see Fixed.

### Fixed

- `partial_trace`, the `PartialTrace` phase, and through it `renyi2`, `mutual_info_renyi2`,
  `SubRenyi2` and `MutualInfoRenyi2`, keep the fermionic signs of the sites traced out.
- `sample` of a mixed state drew every site after about the 1074th of a chain in its last
  state, a probability falling below the smallest float.
- `state` and `matrix` give a copy of the array a site declares, which changing changed it for
  everyone.
- A state can be measured from several threads at once, where its caches, filled with no
  lock, made measurements fail and could crash Julia.
- A projector on a fermionic site refuses a state of no definite parity, which was taken as
  even when named and placed without string when given as a vector.
- A non integer power of a superoperator holding a fermionic operator, as
  `Gate(C + dag(C))^0.5`, is refused on any site but the first, where an MPO placed it wrong.
- A tensor product with a fermionic factor of several sites, as `named(C ⊗ N, "CN")`, was
  taken as even: its adjoint, its matrix and its dissipator were wrong.
- A gate of fermionic operators and of an operator of several sites defined by a matrix, as
  `C(1) * C(4) * M(1, 2)`, gave a wrong state, a string being moved past part of `M`.
- `named` refuses a matrix that does not satisfy its `type`, as a matrix declared an involution
  that is not self adjoint, which `simplify` squared to the identity.
- A `SteadyState` checkpointed on its last sweep hands on a state of trace one once resumed,
  where it handed on the eigenvector dmrg gives, of trace 1.39 in a test.
- An `Evolve` whose duration and time step have opposite signs evolves over its duration,
  where it made no step; `tdvp` and `approx_W` refuse `nsweeps` below one.
- `dmrg` with no sweep left, its `first_sweep` past `nsweeps`, gives the energy of the state
  it is given, where it gave 0.
- On a mixed state the time functions of an evolver take real values: a complex one gave a
  state of complex trace, and is written as its real and imaginary parts, each its own term.
- The docstrings of `approx_W` and `ApproxW` say that a time dependent evolver, its functions
  taken at the middle of each step, limits the order to 2: orders 3 and 4 converge as order 2.
- An integer coefficient is stored as a float, whose products do not wrap around:
  `prod(2Sz(i) for i in 1:63)` measured -1 on sites up, and `(2X)^64` was `0Id`.
- A checkpoint gives back an integer measurement other than an `Int` with its type, where it
  came back as another type, or with JSON 0.21 as a wrong number above `typemax(Int64)`.
- The values of one call of `output` sent to one `Data` by several pairs make one event, one
  row of `data_to_frame`, rather than one per pair.
- `SetState`, on a fermionic site with fermionic sites on its left, changed the sign of the
  correlations across it, or erased them: it lacked their Jordan-Wigner strings.

## [1.6.0] - 2026-09-30

This release adds exported names, among them `set_threading`, `run_steps` and
`close_sim_files`, which is why it is a minor version rather than a patch.

### Added

- `set_threading` and `threading_settings`, with `ThreadingState`, set and read how ITensors
  threads the contractions, in a dense mode or in a block sparse one for systems that conserve
  something, and `SimData` has a `threading` field: `:dense`, the default, `:blocks`, `:auto`,
  which chooses before each phase, or `nothing`.

- The `stamp` file of a simulation records how the run is threaded.

- Phases of one's own are supported: `run_steps` runs one written as a loop of steps, which is
  checkpointed, stopped and resumed between two steps, and `resume_step` gives one driving a
  solver the sweeps it has already done. The manual describes them.

- `close_sim_files` is exported, so that a `Simulation` built by hand writes its json files.

### Changed

- `runTMS` runs a simulation in the dense mode unless its `SimData` asks otherwise, and puts
  the threading it found back when it returns.

- The quantities a site conserves are kept sorted by name, as ITensors sorts the components of
  a charge: two sites conserving the same quantities are equal whatever the order of their
  declaration. A checkpoint written by an earlier version for a site declaring several
  quantities in another order is refused as belonging to another simulation.

- `resume_sweeps!(sim.checkpoint)`, which was not exported, is replaced by `resume_step(sim)`.

- The identity on density matrices prints as `Gate(Id)`, and `Limits` as the call that builds
  it.

- The documentation is reviewed: statements the code contradicted are corrected, pitfalls and
  defaults are documented, and the reference of the high level interface is reorganized.

### Fixed

- A simulation stopped, resumed and completed does nothing when run again, where without
  periodic checkpoints it resumed from the checkpoint of the stop.

- The log of a resumed simulation gives the simulation time it resumes from.

- `mutual_info_renyi2`, `renyi2`, `MutualInfoRenyi2` and `SubRenyi2` take a cut from 0 to the
  number of sites, `SubRenyi2(k)` standing for the sites `1:k`, where it failed when measured,
  and give 0 for empty positions, given in a vector of any element type.

- `weaken` gives a state or a system back unchanged when the target names what it conserves in
  another order, where it built a new system with new indices.

- A simulation whose phases hold an anonymous function resumes when its program is included
  again in the same Julia session, where its checkpoint was refused. A checkpoint written by an
  earlier version for such phases is refused once.

## [1.5.0] - 2026-09-29

This release adds exported names, among them `RealValue`, `ImaginaryValue` and
`ComplexValue`, which is why it is a minor version rather than a patch.

### Added

- `Electron` and `Tj` have the mixed state `"MixedSpin"`, or `"↑|↓"`, one electron of fully
  mixed spin, which a strong conservation of `Ntot` allows where `"FullyMixed"` is refused.

- `named` defines an operator from a matrix or a function of the sites as well as from an
  expression, its type read off the matrix, and builds it on the sites it is given: `named(m,
  "MySwap", Qubit())` acts on two qubits.

- An operator of several sites defined by a matrix, a function or an expression `simplify`
  cannot develop, such as `exp(X ⊗ X)`, can be given its sites, `Operator{2}("P2", m,
  selfadjoint_op, Spin(1))`: split into one site factors, it goes into a hamiltonian and can
  be measured, where it could only be applied as a gate (#14).

- `RealValue`, `ImaginaryValue` and `ComplexValue` declare the values of a measurement real,
  purely imaginary or complex (#19); without them `measure` finds the kind of an operator by
  itself.

### Changed

- The checkpoint and output machinery is rebuilt, a checkpoint recording the state, the
  counts and the outputs of one and the same moment. What a program sees does not change.

- The state of a checkpoint is written to `checkpoint-1.h5` or `checkpoint-2.h5` instead of
  `checkpoint.h5`.

- `CreateState` keeps the time of the simulation unless given its `time_start`. It set it to
  0.

- An `Evolve` phase covers the duration it is given, its time step being adjusted. The
  duration was adjusted instead: a duration of 1 in steps of 0.3 stopped at 0.9.

- `Zd` on a `Qudit` is `mod(N, d)`, so that it is conserved modulo `d`. On `Qudit(2)` it was
  conserved as an integer.

- With periodic checkpoints, a checkpoint is written after the last phase, so that running a
  completed simulation again does nothing.

- `EntanglementEntropy(pos, n)` always writes `n` eigenvalues, padded with zeros.

- `CreateState` refuses by a message to randomize a `State` object into a mixed state, which
  has no purification to draw from.

- `measure` gives each value the kind of its measurement: real, imaginary under the name
  `Im(name)`, or complex. `expect`, `expect1` and `expect2` still give the values as
  computed.

- A `Check` writes each of its two values according to its kind, and still compares them as
  computed.

- The warning about a large dropped part goes through `@warn`, also reports a real part
  dropped from an imaginary value, and is relative to the modulus.

- `simplify` takes the adjoint of a projector to be the projector, so that a projector is
  measured as real.

- A json destination writes a complex number as `{"re": …, "im": …}`.

- The checkpoint file is at version 3, and one of an earlier version is refused.

- A product of more than six factors is abbreviated in the name of a measurement, as a long
  sum already was: `Z(1)*Z(2)*Z(3)*...*Z(9)*Z(10)`.

- The docstrings are rewritten, shorter and clearer.

### Removed

- `output(sim, io, header, data)`, which only the tests used.

### Fixed

- Powers of operators follow two rules: an integer power of zero or above is a product, any
  other a function of the operator, taken as `exp` is. This fixes non integer powers that
  were refused or wrong, complex exponents, `^` of a placed operator and powers of the
  identity left unreduced.

- There is a single identity of each kind. `Id ⊗ Id`, `Left(Id)`, `Id(3)` and the others were
  each recognised in some places only, so that `Id(2) - X(2) * X(2)` did not cancel.

- `expect` checks the sites of an operator as written, as `make_mpo` and `apply` do.

- The coefficient of a time dependent term of several sites was silently raised to the power
  of its number of sites, in `make_mpo`, `make_approx_W1` and `make_approx_W2`.

- `apply(A*B, state)` computed `B*A`. A product of gates is the operator it denotes, its
  rightmost factor acting first: a script writing a circuit from left to right must reverse
  it.

- `dag(C ⊗ C)` lost its sign, so that a hamiltonian written with the adjoint of a hopping
  term was not self adjoint.

- The matrix and the tensor of a tensor product left out the Jordan-Wigner strings between
  its factors, which made some operators of several sites silently wrong.

- An operator of several sites placed on a repeated site, `Swap(1, 1)`, is refused.

- `exp` and `mod` of an operator of one site that is not even lost their Jordan-Wigner string
  after a fermionic site.

- A Jordan-Wigner string was moved across a projector on a vector as if it had a definite
  parity.

- The terms of an odd sum on one site, `(C + dag(C))(3)`, share their string, so that `apply`
  takes it. Products that differ only by their last factor are gathered, which also shrinks
  MPOs: a sum of `Dissipator(C + dag(C))` on four sites went from 40 terms to 13.

- `RandomState{Mixed}` no longer truncates what it draws, which left a trace off one and
  negative eigenvalues, and comes back on the system asked for. `RandomState(state,
  linkdims)` no longer overwrites a state of one site.

- `State{Pure}(system, 1)` was taken for the amplitudes of a single site.

- `sample(state, pos)` failed on a system conserving something.

- A mixed state saved from a partial trace of a charged system could not be used once loaded.

- `SetState` under a strong symmetry gave a state of the wrong trace. It is refused.

- The tensor product of two systems could put one index on several sites, and a gate then
  failed or never returned.

- `weaken` failed on a site whose `conserve` field is not a string.

- A quantity taken modulo, declared with `@def_operators` and conserved, was conserved as an
  integer.

- The example of `@def_operators` declared a fermionic operator as a `plain_op`.

- A `Measure` given among other measurements was taken for a single one. It stands for its
  measurements.

- `dest => "text"` in `final_measures` was written to the log, and failed on a `Data` or json
  destination.

- `data_to_frame` joined values on their time, pairing values never measured together. It
  gives one row per measurement set, which each value records under `"events"`.

- A matrix is written in a json file as the list of its rows, not of its columns.

- A resumed simulation writes what the uninterrupted one writes, which it did not in several
  cases, among them an interrupt during a dmrg sweep, a tolerance reached on the first
  resumed sweep and the final measurements of a stopped phase.

- A value that is not finite no longer makes the checkpoint fail or leaves a json file empty.

- The fingerprint of the phases failed on a phase holding a dictionary or a set.

- `dmrg` with a hamiltonian on a mixed state, and so `GroundState`, is refused. It minimised
  the superoperator `ρ ↦ Hρ + ρH`.

- `steady_state`, and so `SteadyState`, gives its state a trace of one.

- The error messages on the charge of an MPO printed the opposite flux.

- `RandomState{Mixed}` on a system conserving something strongly sent the user to another
  form, which refused as well.

- `ToMixed` ignored its limits when the state was already mixed.

- Docstrings and pages of the manual that the code contradicted are corrected.

- Well defined operations that were refused are accepted, among them `Gate` of a placed sum,
  `controlled(C)`, a renamed operator of several sites, a range as a measurement and an
  integer `noise`.

- `apply` of an `Evolver`, and `Fidelity` or `Overlap` of a reference conserving less than
  the measured state, are refused by a message.

- Complex measurements lost their imaginary part when written by `output` (#19). They take
  two columns.

- The part of a `Check` made on a vector observable is written number by number.

- A vector of measurements inside a set, `"data" => [[Z(1), X(2)]]`, stands for its
  measurements.

- A `Data` destination resumed from a checkpoint gives its matrices back as matrices.

- A complex simulation time no longer makes JSON.jl 1.5 fail.

- A term whose factor vanishes on its site, `C(1) * C(1)`, is left out of the MPO. It failed
  on a charged system (#13).

- A matrix whose size does not fit its sites is refused by a message naming the operator and
  the sites.

- A simulation whose first phase is neither `CreateState` nor `LoadState` is refused when its
  `SimData` is built.

- `Swap` and `controlled(Swap)` go into an MPO and `expect` on qubits that conserve
  something, and the MPO of `Swap` is real.

- An operator computed through an eigendecomposition, as `exp` or a non integer power, is no
  longer refused on charged sites for the rounding it holds between charges.

- A power of a superoperator, `(Gate(X)^2)(1)`, goes into an MPO.

- A resumed `GroundState` or `SteadyState` numbers its sweeps as the phase does.

- A term of coefficient zero is left out of a sum, so that `0 * C + dag(C)` is fermionic.

- A state whose sites are of a type defined in a nested module is loaded back, and its
  simulation resumes.

- `Right` of an operator of several sites carrying no definite charge is refused by a
  message, as `Left` is.

- `expect` refuses by a message an operator placed on no site, and a superoperator.

- A measurement on a charged pure state no longer fails once its orthogonality centre has
  moved past its first site.

- `isfermionic` extends the function of ITensors instead of shadowing it.

- `tensor` of a matrix given one site for several identical ones lays it on as many indices
  of that site, charges included.

- A conserved quantity whose name holds `:`, `;` or `%`, or ends in `*`, is refused.

- Sites conserving a quantity of the same name modulo different numbers are refused on one
  system.

- `@def_operators` refuses a name that stands for an operator with a definition of its own.

- Applying a fermionic pure operator to a mixed state checks its positions as written.

- A time dependent evolver is refused unless it is given one time function per term.

- `partial_trace` keeps the trace of the state. It gave a trace of 1, and NaN for a traceless
  state.

- `data_to_frame` renames a measurement named `time` or `event` to `time_1` or `event_1`.

- `Xd` on `Qudit(1)` is the identity. It was zero.

- The type of an operator is checked against its matrix wherever it is placed, and a wrong
  one is refused: `Operator{1}("Bad", [0 1; 0 0], involution_op)` squared to `Id`.

- No operator stores a signed zero, which kept `simplify` from merging equal terms.

- `Proj` and `SetState` are ordered by the values of their state, so that equal ones merge.

- `apply` takes a null gate, which makes the state null.

- An operator prints so that it reads back as itself: `X ⊗ (Y*Z)` printed `X⊗Y*Z`.

- `weaken(system, symmetries(system))` is the system itself.

- `MutualInfoRenyi2(k)` is named `MutualInfoRenyi2(1:k)`, after the sites it stands for.

- A checkpoint is replaced whole or not at all, even when the process is killed while writing
  it.

- An interrupt no longer leaves a checkpoint whose outputs are ahead of its state.

- A resumed dmrg writes the same log and measurements as the uninterrupted one.

- A complex simulation time with no imaginary part stays complex through a checkpoint.

- The fingerprint of a simulation no longer changes with the version of Julia, which had its
  checkpoint refused after an upgrade.

- Without a directory, an interrupt of `runTMS` goes on to the caller, instead of returning
  as if the run had completed.

- The `error` marker of a failed run is removed when the simulation is run again.

- A phase of one's own driving a solver no longer runs all its sweeps again when resumed: it
  is resumed from its start, unless it reads its resume point with `resume_sweeps!`.

- A dmrg search stopped for a checkpoint no longer logs `Done, dmrg final energy`.

- The fingerprint of the phases no longer depends on what the program imports.

- A state file reads back a site field of a floating point type other than `Float64`, and an
  integer beyond `Int`.

- A state a site declares is the one its name gives, before its generic forms and
  `"FullyMixed"`.

- A local state of the wrong size is refused by a message naming the site.

- Whether a tensor carries a definite charge is computed beforehand, rather than read from an
  error of ITensors.

- A destination named as a file `runTMS` writes in the simulation directory, as `log`, `stop`
  or `checkpoint.json`, is refused.

- A declaration of kind given a range or a view, `RealValue(1:3)`, holds for each of its
  elements. It failed on a `MethodError`.

- The docstrings of `Fermion`, `Electron` and `Tj` list `F`.

## [1.4.0] - 2026-09-24

This release adds conserved quantities and the exported names that go with them, which is
why it is a minor version rather than a patch.

### Added

- Conserved quantities. A site declares one with `conserve`, naming it by one of its own
  operators: `Fermion(conserve = N)`, `Electron(conserve = (Ntot, 2Sz))`,
  `Boson(4, conserve = parity(N))`. Its indices then carry that charge and its tensors become
  block sparse, throughout states, operators, MPOs, evolution, ground and steady states and
  the state files. A state is confined to the sector it was built in, and an operator carrying
  no definite charge is refused by a message naming it.

- `strong`, which declares a conserved quantity a strong symmetry rather than the weak one
  assumed otherwise. Weak asks that the density matrix commute with the charge, which any jump
  operator of definite charge preserves, particle loss included; strong asks that every jump
  commute with it, which keeps the charge of the ket apart from that of the bra and cuts the
  blocks finer, in exchange for a state living in a single sector. Dephasing is the usual
  strong case, and `partial_trace` is unavailable there.

- `weaken`, which takes a state, a system or a simulation down to a lower level of
  conservation, strong to weak or weak to none, or to exactly the quantities a target names,
  building the system it lands on. A phase may thus evolve under a strong symmetry and the
  next one under a weak one, where a jump moving the charge becomes possible; there is no way
  back. `symmetries` reports what a system conserves, and how, in the form `weaken` takes.

- `Weaken`, the phase doing the same within a simulation.

- `RandomState{Mixed}(system, states, linkdims)`, which draws a random density matrix from
  the states its purification starts from. A system that conserves something has no sector to
  draw one in otherwise, and what tracing half of the purification leaves is a mixture over
  the sectors around the one named. `CreateState` uses it for a mixed state given with a
  `state` and `randomize`, which it used to refuse.

- `entanglement_by_sector`, the entanglement across a cut of a pure state resolved by the
  charge the left part carries: for each charge, its probability and the entropy and spectrum
  held inside it, which add up to the entanglement entropy with the number entropy on top.

- `flux`, the charge an operator carries on a site, `parity` and `mod`, which reduce a
  quantity modulo an integer, and `named`, which renames an operator so that two quantities
  conserved separately do not go under one name.

- `N`, the number of excitations, on `Qubit` and on `Spin`, where the other site types
  already had it. On a qubit it is `1/2 - Sz`, that is the projector on `"Dn"`, and on a spin
  it is `s - Sz`, the counting of the Holstein-Primakoff mapping. Its eigenvalues are integers
  for every spin, half integer ones included. Note that `Sm` is then what raises it, `"Up"`
  being the empty state.

### Changed

- The state file format is now version 2. The fields of a site are written as a string each,
  together with the kind of value they hold, so that a site whose field is not a number can
  be saved at all — a symbol naming what the site conserves, for instance. Version 1 files
  are still read, which matters beyond old files: a checkpoint left by an earlier version
  holds a state in that format, so resuming one depends on it.

### Removed

- **`A` on `Fermion`, and `Aup` and `Adn` on `Electron` and `Tj`.** They were the bare local
  operators, without their Jordan-Wigner string, which serve to write the strings by hand as
  one does with ITensor. `simplify` inserts the strings itself, and being declared as plain
  operators they were taken to commute with `F`, so a string crossing one of them on its site
  gave the wrong sign: `C(3) * A(1)` came out as the opposite of `C(3) * C(1)`. `A` remains
  the destruction operator of `Boson` and `Qboson`.

### Fixed

- `save_state` no longer destroys the state already saved under a name when it cannot write
  the new one. It deleted the old group before reading the fields of the sites, so a state it
  was going to refuse took the previous one with it and put nothing in its place. The fields
  are read before the file is opened now, and a refusal names the site and the field it
  cannot carry instead of surfacing as a `MethodError` raised by `convert` inside HDF5.

- `matrix` and `tensor` accept a single site for a superoperator acting on several identical
  ones, as their documentation says and as they already did for other operators:
  `matrix(Left(Swap), Qubit())` raised a `DimensionMismatch`.

- **A sum of fermionic operators on one site had the wrong sign when the Jordan-Wigner string
  of a later operator crossed it.** `simplify` gathers the terms of one site into a single
  factor and took that factor to commute with `F`: `C(3) * (C + dag(C))(1)` came out as the
  opposite of `C(3) * C(1) + C(3) * dag(C)(1)`, and a hamiltonian written with the later site
  first, `sum(γ(i+1) * γ(i))` with `γ = C + dag(C)`, had its first bond with the wrong sign.
  The parity of every factor is now read, and a factor that has none stops the string.

- **`dag` of a non integer power was the power of `dag`**, which only holds when the operator
  has no negative eigenvalue: `dag(sqrt(X))` simplified to `sqrt(X)`, so `expect` returned
  the conjugate of the right value and the MPO of `Dissipator(sqrt(X))` was wrong.

- **A fermionic operator inside a superoperator or a tensor product was applied as a gate
  without its Jordan-Wigner string.** `apply` only looked for one on a one site factor, so
  `Gate(C)(2)`, `Left(C)(3)` and `(C ⊗ dag(C))(1, 2)` were built from bare matrices, the last
  one coming out as the opposite of `C(1) * dag(C)(2)`, which it is by definition. The whole
  operator is now looked into, and what cannot be placed as a gate is refused: a dissipator
  of a fermionic operator, which turns into a sum, and a function of a fermionic operator
  acting on several sites, such as its exponential, which leaves no room for a string.

- `partial_trace` refuses a position the state does not have, and names it. Tracing out a
  site beyond the last one left the state whole without a word, and keeping one raised a
  `BoundsError`. `renyi2` and `mutual_info_renyi2` given positions, and the `PartialTrace`
  phase, go through it.

- `entanglement_entropy` returned `NaN` when a singular value was exactly zero, which a
  `mindim` above the Schmidt rank keeps: a zero now adds nothing to the entropy.

- `expect` of a superoperator, `Left(X)(1)` for instance, says that it takes an observable,
  where it raised a `MethodError` on `length`.

- A sum mixing a `Proj` given by a number with one given by a name, `Proj(1) + Proj("Up")`,
  raised a `MethodError` on `isless` when simplified. `Proj` is ordered by how its state
  prints, as `SetState` already was.

- `Limits(cutoff = 0)` is accepted, the cutoff being converted to a float, where a field of
  union type refused an integer with a `MethodError` on `convert`.

- Under Julia 1.10, a non integer power of an operator whose matrix is real, diagonal and
  has a negative entry raised a `DomainError`, so that `Sz^0.5` could not be used at all.
  Its matrix is computed complex, as later versions of Julia do.

- A gate of several sites taking a state to zero, and the sum or difference of two zero
  states, raised a `BoundsError` inside ITensors. `mindim` defaulted to 0, which ITensors
  does not expect and which made it truncate a spectrum of zeros past its first value. It
  defaults to 1, the least a bond can have, and a smaller value is taken as 1; a checkpoint
  written with the former default is still resumed.

- Asking for the `F` of a site type that has none, the identity, wrote it into the operator
  library, so that declaring `F` for that type afterwards was refused as a redefinition.

- A non integer power of an operator with a complex coefficient, `(im * X)^0.5` for
  instance, took the coefficient out whole and was another determination of the power than
  the principal one its matrix has. Only the modulus of the coefficient comes out now, its
  phase staying inside the power, as the sign of a negative one already did.

- An operator acting on several sites at once with no expression to be replaced by, one
  defined by a matrix or a function of one such as `exp(X ⊗ X)`, failed inside `expect` and
  `make_mpo` on a `MethodError` or a `BoundsError`. It is refused by a message saying that
  it can only be applied as a gate, and the manual says so where it shows one.

## [1.3.0] - 2026-09-21

This release carries a user visible change of behaviour, the removal of MKL, which is why
it is a minor version rather than a patch.

### Added

- `Qudit` site type, and richer operator sets on several existing site types.
- Checkpointing: a simulation can be stopped and restarted from where it stopped, driven by
  `checkpoint_interval` and `max_time`.
- `SaveState` and `LoadState` phases, to write a state to disk and read it back.
- `sample`, to draw measurement outcomes from a state.
- `SetState`, and an integer parameter to `Proj` and `state`, which `controlled` now uses.
- `square_lattice`, alongside the existing graph helpers.
- Multi site function operators.
- `matrix` and `tensor` for `Dissipator`.
- `Evolve` accepts vectors in `Limits`.
- `mindim` in `Limits`, the bond dimension a truncation is not allowed to go below. It is
  passed on to every ITensor call that takes one, and defaults to zero, which is no
  minimum. `maxdim` keeps the last word when the two ask for opposite things.
- `inner` and `dot` between two states, and the fidelities that follow: `fidelity` for two
  pure representations and for a pure one against a mixed one, `hs_fidelity` for the
  normalised Hilbert-Schmidt overlap of two mixed ones. The Uhlmann fidelity of two mixed
  states is deliberately absent, needing the spectrum of a density operator, and asking
  for it says so and points at `hs_fidelity`.
- `State(system, state)`, which puts a state on another system of the same sites, and the
  `system` keyword of `load_state`, which reads one straight onto an existing system. A
  `System` carries ITensor indices of its own, so this is what makes two states built
  apart comparable at all.
- `variance(hamiltonian, state)` and the `Variance(hamiltonian)` measurement, the
  convergence check of a ground state search and what gives it an error bar.
- `Fidelity(ref)` and `Overlap(ref)`, to follow either against a reference state while a
  simulation runs.
- `CITATION.cff`, `CONTRIBUTING.md` and this changelog.

### Changed

- **Eleven names are no longer exported**, being machinery no program writes: the types the
  operator algebra builds (`Identity`, `AtIndex`, `JW`, `JW_F`, `Multi_F`), the type
  parameters (`PM`, `GI`, `Generic`, `Indexed`), `sim`, which is `System(system.sites)`
  under another name and was called once in the whole package, and `removeMulti`. All of
  them remain reachable as `TensorMixedStates.name`, and none is a name a program writes.
- **The `alg` keyword of `steady_state` is now `mpo_algo`**, the name the `SteadyState`
  phase already gave the same thing. The old one goes on working for this cycle and warns.
- **MKL is no longer a dependency.** It was used for the single purpose of switching the
  BLAS backend of the whole Julia session, which is a decision that belongs to the
  application rather than to a library, and `MKL_jll` ships `x86_64` Linux and Windows
  binaries only, so the package could not be installed on macOS or on any ARM machine.
  OpenBLAS, which Julia ships, is now used as it comes. Users on `x86_64` who want the
  performance of MKL load it themselves before TMS; the Installation section of the manual
  says how.
- `@def_operators` binds an operator name only the first time it sees it. A name already in
  scope is registered for the new site and checked against what it already stands for,
  rather than bound again. This makes the sharing of a name between site types explicit and
  enforced instead of accidental, refuses a declaration whose `OpType` disagrees, and
  removes an error that made the package unusable from Julia 1.10 and 1.11 as soon as a
  user declared an operator whose name a loaded site module exported.
- **`approx_W` now takes `order` with no default and defaults `w` to 2**, which aligns it
  with the `ApproxW` phase whose docstring already described that behaviour. Code calling
  `approx_W` directly without `order` no longer runs, and code that left `w` out now gets
  WII where it got WI. Phases are unaffected: `ApproxW` passes both explicitly.
- **A `SimData` can no longer be used as a phase of another simulation.** It has the shape
  of a phase only because `runTMS` runs the top level one through the same machinery. Use
  nested vectors to build a list of phases in pieces, which is the documented way and is
  flattened on construction.
- `julia = "1.10.5"` in `[compat]`, which reads as `[1.10.5, 2.0.0)`. The package no longer
  becomes uninstallable on the day a new Julia minor version is released.
- `simplify` always expands multi site operators, which makes its result predictable.
- Global identifiers are `const`.
- Tensors are kept real instead of complex whenever possible.
- Measurement names, the display of phases and the were made more regular.
- The reference article is published: it is
  [SciPost Phys. Codebases 72 (2026)](https://doi.org/10.21468/SciPostPhysCodeb.72), and the
  README and manual cite it instead of the preprint.

### Deprecated

- **`EE` is now `EntanglementEntropy` and `Linkdim` is now `MaxLinkdim`**, each named after
  the function it measures, as the other state functions are after theirs. The labels
  written to the output files follow the new names, so a column that read `EE(3)` now reads
  `EntanglementEntropy(3)`. `Linkdim` is a binding rather than a function, and a deprecated
  binding only warns on a qualified access, so this line is its notice.
- **`DataToFrame` is now `data_to_frame`**, spelled like the other functions of the package.
- **`Mutual_Info_Renyi2` is renamed `MutualInfoRenyi2`**, which matches the naming of every
  other measurement. The old spelling goes on working and forwards to the new one, but it is
  marked deprecated and the label written to the output files is the new name, so a script
  reading a column by its header has to follow.
- **`Dmrg` is marked deprecated**, use `GroundState`. It has carried the word in its
  docstring for several versions without telling anyone; it is now a deprecated binding.

Both warnings only show with `--depwarn=yes`, which is what running a test suite does. An
ordinary run stays silent, so this file and the docstrings are the notice.

### Fixed

Many bugs have been fixed, in particular:
- **`expect`, `expect1` and `expect2` returned unnormalised values on a pure state in certain corner cases.**
- **Operator equality was structural for some forms and identity based for others**, so
  `2X(1)*Y(2)`, `dag(X*Y)`, `Left(X*Y)`, `Phase(0.3)` and `controlled(Z)` never compared
  equal to themselves. `==` and `hash` are now read from the type, once, for the whole
  hierarchy. No result was wrong, the two being false together, but a measurement shared
  by two observables was computed twice and terms that could have merged did not.
- `renyi2`, `mutual_info_renyi2` and `partial_trace` take any vector of integers, a range
  included, where they demanded a `Vector{Int}` and refused `1:3` with a `MethodError`.
  `partial_trace` on a pure representation now says what to do instead of raising one.
- `renyi2` and `mutual_info_renyi2`.
- `dag` and `expect2` on fermionic operators, and `expect` in the presence of `Multi_F`.
- Time dependent evolution, which was broken.
- Files and data sets are properly terminated.
- Applying gates is safer, and `expect` no longer trips on an unsimplified expression.

### Documentation

- The examples of the manual and of the reference pages are executed at every documentation
  build, so the build fails when an example stops working. Converting them found and fixed
  five errors that indented code blocks had been hiding.
- Versioned documentation is published: `stable` now points at the latest release instead of
  the development branch.
- A new section on how the MPO is built and what it costs: TMS gives each term of a sum a
  channel of its own and does not compress, so a long ranged operator gets a much larger MPO
  than `ITensorMPS`' `OpSum` would give, while a short ranged one costs exactly the same.
  The reason is that the form obtained is what WI and WII need and what makes time dependent
  evolvers free.
- New sections on the choice of the BLAS backend and on reusing an operator name across
  site types, and a word asking authors who use TMS to acknowledge it and to get in touch.
- A word on the lowercase names the package brings into scope. The manual advised keeping
  one's own identifiers lowercase, which avoids the capitalised operator names but points
  straight at the fifty lowercase function names of the package, `state`, `sim`, `output`
  and the like. Assigning to one of them shadows it, and Julia 1.10 and 1.11 refuse the
  assignment outright when the name has already been used. The manual now says so, and
  neither the manual nor the examples shadow one any more.

### Development

- Continuous integration runs on the development branch, not only on `main`, and the matrix
  covers the minimum Julia version that `[compat]` promises, the current stable one and the
  upcoming release, on Linux and on Apple Silicon.
- Coverage is measured on every CI run and published to Codecov, with a badge in the
  README. The tests already ran instrumented, so this only gathers and sends what was
  produced and thrown away before.
- The test suite went from about 190 assertions to more than 600, organised in groups that
  can be run selectively.
- The precompilation workload covers `dmrg` and `apply`, which it did not, so a first call
  to either takes about a second instead of five to eight. Measured on the whole of a first
  session: 100 seconds with no workload, 44 with the old one, 28 with this one, for six
  seconds more of precompilation.
- The examples of `examples/high_level` are run by a step of their own in CI, on one job
  of the matrix, so that one of them breaking is noticed. They were shortened to make that
  affordable, three and a half minutes for the five. The examples of `examples/article`
  stay at the sizes the article published and are not run.
- The Codecov patch status is informational. It marks a commit red when the lines it
  changes are less covered than the project as a whole, which counts a rewritten error
  message as new untested code although it adds none. The figure is still reported on the
  commit and in pull requests, it just no longer fails.
