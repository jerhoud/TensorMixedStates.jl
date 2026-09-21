# Changelog

Notable changes to TensorMixedStates are recorded here, in the format of
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/). The package should follows
[semantic versioning](https://semver.org/spec/v2.0.0.html).

This file starts with the release in preparation. For what came before, see the
[tags](https://github.com/jerhoud/TensorMixedStates.jl/tags); version 1.2.8 is the one
archived as the [codebase release](https://doi.org/10.21468/SciPostPhysCodeb.72-r1.28) of
the reference article.

## [Unreleased]

The next release carries a user visible change of behaviour, the removal of MKL, so it is
expected to be 1.3.0 rather than a patch.

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
- `inner` and `dot` between two states, and the fidelities that follow: `fidelity` for two
  pure representations and for a pure one against a mixed one, `hs_fidelity` for the
  normalised Hilbert-Schmidt overlap of two mixed ones. The Uhlmann fidelity of two mixed
  states is deliberately absent, needing the spectrum of a density operator.
- `State(system, state)`, which puts a state on another system of the same sites, and the
  `system` keyword of `load_state`, which reads one straight onto an existing system. A
  `System` carries ITensor indices of its own, so this is what makes two states built
  apart comparable at all.
- `CITATION.cff`, `CONTRIBUTING.md` and this changelog.

### Changed

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
- `simplify` now merges the factors of a generic product until nothing moves, instead of
  making a single pass. `X*Y*Y*X` gives `Id` where it gave `X*X`, and the result no longer
  depends on how many times `simplify` is called. Indexed operators, and therefore every
  MPO, were already reduced fully and are unchanged.
- Global identifiers are `const`.
- Tensors are kept real instead of complex whenever possible.
- Error messages name what was not found and what was expected.
- Measurement names, the display of phases and the were made more regular.
- The reference article is published: it is
  [SciPost Phys. Codebases 72 (2026)](https://doi.org/10.21468/SciPostPhysCodeb.72), and the
  README and manual cite it instead of the preprint.

### Deprecated

- **`Mutual_Info_Renyi2` is renamed `MutualInfoRenyi2`**, which matches the naming of every
  other measurement. The old spelling goes on working and forwards to the new one, but it is
  marked deprecated and the label written to the output files is the new name, so a script
  reading a column by its header has to follow.
- **`Dmrg` is marked deprecated**, use `GroundState`. It has carried the word in its
  docstring for several versions without telling anyone; it is now a deprecated binding.

Both warnings only show with `--depwarn=yes`, which is what running a test suite does. An
ordinary run stays silent, so this file and the docstrings are the notice.

### Fixed

- **`expect`, `expect1` and `expect2` returned unnormalised values on a pure state in certain corner cases.**
- **Operator equality was structural for some forms and identity based for others**, so
  `2X(1)*Y(2)`, `dag(X*Y)`, `Left(X*Y)`, `Phase(0.3)` and `controlled(Z)` never compared
  equal to themselves. `==` and `hash` are now read from the type, once, for the whole
  hierarchy. No result was wrong, the two being false together, but a measurement shared
  by two observables was computed twice and terms that could have merged did not.
- **A period of zero behaved differently in each of the four places the library has one.**
  `measures_period = 0` raised a division by zero on the first sweep, a negative
  `n_expand` or `n_hermitianize` was read as one sweep out of two, and a negative
  `checkpoint_interval` put the next checkpoint in the past and kept it there, writing the
  whole state to disk on every sweep. A period below one now means never, everywhere:
  `sweep_due` carries the rule for the three sweep counters and `checkpoint_due` applies
  it to the interval in seconds.
- **Site index tags no longer depend on what the user imported.** The type name was
  printed with its module prefix when the site module was not in scope, and ITensors cuts
  a tag at 16 characters, so every site type came out tagged `TensorMixedState`.
- `renyi2`, `mutual_info_renyi2` and `partial_trace` take any vector of integers, a range
  included, where they demanded a `Vector{Int}` and refused `1:3` with a `MethodError`.
  `partial_trace` on a pure representation now says what to do instead of raising one.
- **A `SimData` nested inside the phases of another silently skipped its inner phases.**
  The loop it opened shared the phase counter of the loop around it, so every inner phase
  whose index was below the outer one was passed over, on the very first run and with no
  checkpoint involved. It is now refused outright.
- An object that is not a phase, or one that has the fields of a phase but no `run_phase`
  method, now says so instead of surfacing as a `MethodError` or a `FieldError` from the
  middle of a run.
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
- The free boson testset no longer measures the arithmetic of the machine it runs on. It
  evolved at a bond dimension where the truncation is unstable, so the result followed the
  BLAS thread count of the runner: the same code gave a green Windows job one hour and a
  red one the next, missing its tolerance by one percent. Raising the bond dimension from
  10 to 16 removes the dependence, and the recorded references hold to 8.6e-7.
