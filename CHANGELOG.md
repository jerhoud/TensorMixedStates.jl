# Changelog

Notable changes to TensorMixedStates are recorded here, in the format of
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/). The package should follows
[semantic versioning](https://semver.org/spec/v2.0.0.html).

This file starts with the release in preparation. For what came before, see the
[tags](https://github.com/jerhoud/TensorMixedStates.jl/tags); version 1.2.8 is the one
archived as the [codebase release](https://doi.org/10.21468/SciPostPhysCodeb.72-r1.28) of
the reference article.

## [Unreleased]

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
  the sectors around the one named.

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

### Fixed

- `save_state` no longer destroys the state already saved under a name when it cannot write
  the new one. It deleted the old group before reading the fields of the sites, so a state it
  was going to refuse took the previous one with it and put nothing in its place. The fields
  are read before the file is opened now, and a refusal names the site and the field it
  cannot carry instead of surfacing as a `MethodError` raised by `convert` inside HDF5.

- `matrix` and `tensor` accept a single site for a superoperator acting on several identical
  ones, as their documentation says and as they already did for other operators:
  `matrix(Left(Swap), Qubit())` raised a `DimensionMismatch`.

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
