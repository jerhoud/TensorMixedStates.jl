# Remaining work

What is left after a full review of the library, and after the correctness pass that
followed it. Every item below was checked against the source, not copied from the review:
where a claim is quoted, it is one that still reproduces.

The bug list of that review (B1–B25) and its wrong documentation examples (D1–D12) are
done. Two of the reported bugs, B17 and B19, turned out not to be bugs; what was real
about them was fixed instead. The test suite went from about 190 assertions to 601, all
passing.

Items are ordered by what they cost against what they buy.

---

## 1. Documentation

The largest remaining gap between what the library does and what a reader is told.

### 1.1 Make the examples executable — **done**

Thirty-two `@example` blocks are now run at every documentation build, spread over
`manual.md` (17), `operators.md` (8), `measurements.md` (5) and `others.md` (2). The
remaining blocks are genuine pseudo-code (`options...`, `[...]`, undefined placeholders)
or the forty site simulation, and were turned into fenced `julia` blocks so that they are
at least highlighted. `index.md` and `sites.md` only contained non-runnable blocks and got
the same treatment.

`@example` was preferred to `jldoctest` everywhere: it fails the build when the code
raises, which is the guarantee wanted, without comparing output text that contains random
`ITensor` index ids and version dependent float formatting. The single existing
`jldoctest` of `src/Sites.jl` is kept, and `docs/make.jl` now sets a `DocTestSetup` so
that future doctests in docstrings have the site modules in scope.

Verified by breaking one example on purpose: the build stops with
`makedocs encountered an error [:example_block]`. The documentation CI is now blocking.

The conversion found and fixed five errors that the indented blocks had been hiding:

- the manual was not followable — after `using .Qubits` it used `Boson(4)`, `Spin(3/2)`
  and `Fermion()`, three `UndefVarError` for whoever copied it;
- `state = (state1 + state2) / 2` added a pure state of `system1` to a mixed state of
  `system2`, which the sentence just above it says is impossible;
- `mix(purestate)` used a variable defined nowhere;
- `Swap = Operator{2}("Swap", ...)` shadowed the predefined `Swap` of `Qubits`, right
  after the manual advises keeping one's own identifiers lowercase;
- the three measurement destination examples were indented under a bullet, which in
  CommonMark is running text, not code.

Two things were added in passing: a `@contents` block at the top of `manual.md` and
`operators.md` (part of 1.2, free here), and `TensorMixedStates` declared in
`docs/Project.toml` with `[sources] path = ".."`, so that `julia --project=docs
docs/make.jl` works without the `Pkg.develop` step that only CI performed.

### 1.2 Short fixes — **done**, except one deferred point

- The site type count is fixed: `manual.md` says eight and lists `Qudit`, and the
  parametric constructors now list the four of them, `Qudit(dim)`, `Boson(dim)`,
  `Spin(s)` and `Qboson(q, dim)`, with what each argument means.
- **`SimpleOp` has its own docstring** (`src/Operators.jl`). The published page used to
  show the docstring of `Op` — wrong signature, wrong description, and the `ig` for `is`
  typo of that docstring on display. The replacement is two lines, matching `Op`,
  `GenericOp` and `IndexedOp` next to it: the type is `GenericOp{Pure, 1}` and it is
  abstract. Nothing more, because no user ever writes `SimpleOp`: it appears only in
  internal signatures and never in an error message. The knowledge that reading a `name`
  field off a `SimpleOp` is wrong already lives where it is needed, in the comment at
  `src/Measure.jl:124`.
- **The `modules` point of the review does not reproduce.** Adding the eight site
  submodules to `modules = [...]` was tried and changes nothing: those modules carry no
  docstring of their own, they only re-export constants defined in `TensorMixedStates`,
  so `checkdocs` has nothing new to look at. `make.jl` was left alone.

  What does surface the gap is `checkdocs = :all`, run once locally: 45 docstrings appear
  in no `@docs` block. Of those, only `Swap` (`src/Qubits.jl`) and `Sumd`
  (`src/Qudits.jl`) are things a user writes, and they now have a reference entry in
  `sites.md`. The rest — the concrete operator subtypes, the `Base` overloads, the site
  library lookups, the whole checkpointing machinery — is implementation and stays out
  deliberately. `checkdocs = :exports` is the expression of that choice and stays.
- **`measures_period` was already published** through the docstrings of `Evolve`,
  `GroundState` and `SteadyState`; what was missing was a mention in the manual's
  narrative, which it now has. The `:sweep` and `:energy` symbol measurements are
  documented in `measurements.md`, including the trap: a symbol the running algorithm does
  not provide yields an empty value rather than an error, so `:energy` in an `Evolve`
  phase writes nothing.
- `@contents` was added to `manual.md` and `operators.md` in 1.1, and `index.md` now ends
  with an `@index`. Note that the index only lists `@docs` entries, so `X`, `Sp`, `N` and
  the other site operator constants are still anchorless: they are described as bullets in
  the docstring of `Qubit` and the like. Giving them real docstrings is a different job.
- **Spelling** fixed everywhere in the docstrings and the markdown: `measurments`,
  `alorithm`, `beyween`, `reprensent`, `tranformed` (twice), `substracted`, `simultation`,
  `parametred` (three times), `ig` for `is`, and `features requests`.

  `descritpion` is deliberately left as it is. It is a *field name* of `SimData`, so the
  docstring spells it the way the code does; correcting the word would mean renaming the
  field, and the API is not changed for a spelling mistake.

**Deferred**: documenting `TMS_SKIP_PRECOMPILE_WORKLOAD`. It is a development shortcut,
not a user facing option, so its place is the `CONTRIBUTING.md` of 2.4 rather than the
published manual. To be written there.

---

## 2. Packaging and release

### 2.1 Versioned documentation has never been published

`origin/gh-pages` contains only `dev/`: `versions.js` reads `DOC_VERSIONS = ["dev"]` and
`DOCUMENTER_STABLE = "dev"`, despite tags up to `v1.2.9`. The README's Documentation link
therefore points at the development docs. `documentation.yml` does trigger on tags, so the
cause has to be looked for in the Actions logs.

### 2.2 Continuous integration

- **1.10 is never tested.** The matrix is `['1.11', 'pre']` while `[compat]` promises
  `1.10.5` and `docs/src/index.md` promises it to the user.
- **CI only triggers on pushes to `main`** (`CI.yml:4-6`), while development happens on
  `dev`. Fifty commits have just landed there without CI running once.
- **No coverage measurement and no badge.** For a package this size it is the cheapest
  thing left to add.
- Neither macOS nor Windows. `test/checkpoint.jl` changes the working directory and
  manipulates files, which is precisely what breaks on Windows.
- **`test/Project.toml` and `docs/Project.toml` have no `[compat]` section**, so a future
  release of Aqua or Documenter can redden CI for reasons unrelated to the package.
- `documentation.yml` uses `julia-actions/cache@v2` where `CI.yml` uses `@v3`, has no
  `concurrency` block and no `timeout-minutes`, and still carries the PkgTemplates
  boilerplate comment.

### 2.3 Two decisions

- **MKL is a hard dependency**, used once, in `src/TensorMixedStates.jl`, to switch the
  BLAS backend. `MKL_jll` ships x86 binaries only, so the package is **not installable on
  Apple Silicon**, and this is also what keeps macOS out of the CI matrix. Moving it to an
  extension behind a platform test, or at least documenting the fallback to
  `LinearAlgebra.BLAS`, would be a real portability gain. The Installation section of
  `docs/src/index.md` mentions no platform restriction.
- **`julia = "1.10.5 - 1.13"` is inclusive**, so `[1.10.5, 1.14.0)`. The package becomes
  uninstallable the day 1.14 is released, until a new release is cut, and the bound has
  already had to be raised once. Unless a 1.14 incompatibility is known, `julia = "1.10.5"`
  is the convention and removes the recurring maintenance.

### 2.4 What a published package is missing

By value: `CITATION.cff` (the README and `index.md` both ask for the HAL preprint to be
cited, but without a machine readable file GitHub's "Cite this repository" button never
appears); `CHANGELOG.md`; coverage badge; `CONTRIBUTING.md` — non-trivial process
knowledge already exists in the wrong place, namely the test group convention in
`test/runtests.jl:6-13`, the selective run via `Pkg.test(test_args = [...])`, the
"what belongs here" file headers and the `TMS_SKIP_PRECOMPILE_WORKLOAD=true` shortcut;
a fuller README (a ten line runnable quick start, the installation line, **the licence** —
GPL-3.0 appears nowhere in the README, and it is a strong copyleft choice next to ITensors
under Apache-2.0, which deserves to be visible — and the logo that already exists);
issue and PR templates.

---

## 3. Tests

- **No RNG seed is fixed anywhere** in `test/`. `test/observables.jl:7-40` samples 2000 and
  4000 times with `atol = 0.05` and `0.03`, about 4.5σ and 4.7σ, so a false failure roughly
  once in 10^5 runs. `sample` already accepts an `rng` that the tests never pass, so this
  is free to fix. `test/checkpoint.jl` is already deliberately deterministic, triggering on
  a measurement counter rather than a clock; the same reflex belongs on the RNG.
- **Restoring `exit_on_sigint` has no test.** It can only be checked by sending a real
  SIGINT from a child process, which was judged not worth its cost. It is the one fix of
  the correctness pass with no coverage.
- **Recorded reference values** in `test/algorithms.jl:81-87`, `:103-106` and `:122-124`
  are compared to 11–15 significant digits with no comment saying where the numbers come
  from. If they come from a previous run, they freeze the implementation's behaviour at
  that date rather than checking it. The steady state values nearby are visibly exact
  fractions and do not have this problem.

---

## 4. API and design

These need decisions, not patches.

- **The implicit invariant on site operator names is the one to worry about.** Every site
  module defines its operator constants in the same module (`src/Sites.jl`), so `Sp`, `Sz`,
  `N`, `A`, `S` are defined two to four times. It works only because the returned
  `Operator{1}(name, nothing, type)` happens to be identical. The day a new site type
  registers one of those names with a different `OpType`, the last definition silently wins
  for every site type. `add_operator` only checks for duplicates within one site type.
- **`Simulation` is immutable on the surface only**: the copy constructor shares `files`,
  `data` and `checkpoint` by reference. State is threaded functionally, but the whole
  checkpoint bookkeeping is a mutable channel shared between every copy. That was the root
  cause of the interrupted-phase checkpoint bug.
- **`runTMS` is neither reentrant nor usable in parallel**: `cd`, `Base.exit_on_sigint`, a
  global `Random.seed!`. Defensible for a script driver, but written nowhere.
- **`Dmrg` is declared deprecated** (`src/PhaseTypes.jl`) yet is still exported and
  documented without a deprecation warning.
- **Divergent defaults** between phase structs and solver signatures: `ApproxW.w = 2`
  against `approx_W(...; w = 1)`; `order` with no default on the phase side and `= 1` on the
  solver side. Divergent naming: `alg` on the function, `mpo_algo` on the phase.
- **`run_phase` has no fallback method**, so an unsupported object in `phases` gives a raw
  `MethodError`. It is also a clean user extension point — define a struct with
  `name`/`time_start`/`final_measures` plus a `run_phase` method — but this is documented
  nowhere and `run_phase` has no docstring.
- **`const Phases` does not include `SimData`**, although `SimData` *is* a phase, and
  `flatten_phases` does not flatten a nested `SimData`, which would therefore be accepted
  and would silently break the phase numbering a checkpoint records.
- Minor: duplicated exports (`Limits`, `mix`, `dag`, `tensor`); `removeMulti` in camelCase
  among `make_mpo`, `partial_trace`, `graph_base_size`; `trace2` and `Purity` are the same
  function under two names; no `show` for `JW`, `Multi_F`, `Proj`, `SetState`, and
  `test/operators.jl` depends on the default struct display.

---

## 5. Performance

Only two targets are worth the effort, and both are on the contraction side rather than in
Julia-level micro-optimisation.

1. **Memoise the local tensors in `expect2`** (`src/Observables.jl`). For each pair `(i,j)`
   it rebuilds the one site tensors from scratch: library lookup, the matrix product
   `o1*F`, a new `Index`, an `ITensor`, a `combiner` and two contractions. That is
   O(n²·|ops|) tensor constructions that depend only on `(site, operator)`. This is the
   most profitable single change in the file.
2. **The MPO is never compressed.** `PreMPO!` (`src/Mpo.jl`) allocates one private channel
   per term, so the MPO bond dimension is `2 + #(terms crossing the link)`, where
   ITensorMPS' `OpSum` construction applies an SVD compression. For a long range
   Hamiltonian with O(n²) terms this gives a bond dimension in O(n²) instead of O(n). It is
   *the* major algorithmic difference with ITensorMPS and it is not documented. The
   counterpart is real: the resulting triangular structure is exactly what WI/WII need, and
   the coefficients can be changed from sweep to sweep without rebuilding the term list.
   Whatever is decided, the cost belongs in the documentation.

Lesser points, all measurable but small: scalar `setindex` filling of the MPO tensors,
`Matrix(undef, n, n)` giving a `Matrix{Any}` for the `expect2` result, single pass
simplification of generic products (`simplify(X*Y*Y*X)` gives `ProdOp([X,X])` rather than
`Id`, which costs MPO bond dimension), abstract field types in the hot symbolic layer, and
`stop_requested` doing an `isfile` on every sweep.

---

## 6. Examples

The six scripts of `examples/article/` and the five of `examples/high_level/` are all up to
date with the current API and abundantly commented, but:

- **none is run by CI**, and none uses the recent features (`checkpoint_interval`,
  `max_time`, `square_lattice`, `Qudit`, time dependent evolvers);
- **`examples/high_level/ising_quench.jl` writes a description that lies**: `limits` sets
  `maxdim = 50` and the description text says `maxdim = 100`. That text is written to the
  run's `description` file, so the simulation records a false parameter;
- `examples/article/4_Free_fermions_with source.jl` **has a space in its filename**, unlike
  its five siblings; simulation names with spaces appear in three more scripts, and `name`
  becomes a directory name;
- `examples/article/1_Fermion_chain_with_dephasing.jl` has a comment saying the evolver
  must be `-im * hamiltonian + dissipators` immediately above a line writing `im*sum(...)`.
  The code is right and the comment contradicts it;
- there is no README in `examples/high_level/`, and `examples/article/README` should be
  `README.md` to be rendered by GitHub;
- no indication of running time anywhere, and these are publication sizes;
- **`docs/src/index.md` says "Working examples are presented in the folder `examples`"
  with no hyperlink**, no mention of the two subfolders and no word on what each shows. A
  reader of the published site has no path to the eleven scripts.
