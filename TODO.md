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

The one deferred point, documenting `TMS_SKIP_PRECOMPILE_WORKLOAD`, is **done**: it is a
development shortcut rather than a user facing option, so it went into the `CONTRIBUTING.md`
of 2.4, with the warning not to rely on it anywhere a user could end up.

---

## 2. Packaging and release

### 2.1 Versioned documentation — **done**

The cause was in the Actions logs: no run had ever been triggered by a tag. TagBot pushed
every tag with the default `GITHUB_TOKEN`, and GitHub does not trigger workflows on pushes
made with that token, so the `tags: '*'` of `documentation.yml` never fired. The repository
had no secret and no deploy key at all, so the `ssh: ${{ secrets.DOCUMENTER_KEY }}` of
`TagBot.yml` was empty and silently fell back to the token.

A `Documenter` deploy key with write access and the `DOCUMENTER_KEY` secret are now in
place, so future tags deploy on their own. `v1.2.9` itself was built from the tag and
assembled with Documenter's own `expand_versions`, `generate_version_file`,
`generate_redirect_file` and `rm_and_add_symlink`, because the workflow file at that commit
carried no `workflow_dispatch` trigger to dispatch on. `gh-pages` now holds `v1.2.9/` with
the `stable`, `v1.2` and `v1` symlinks, and the root redirects to `stable` instead of
`dev`. `workflow_dispatch` was added to `documentation.yml` so any later tag can be rebuilt
by hand with `gh workflow run documentation.yml --ref vX.Y.Z`.

### 2.2 Continuous integration — **done except coverage**

- CI and the documentation build now trigger on `dev` as well as `main`, so work stops
  landing untested.
- The matrix is `['1.10', '1', 'pre']`: the minimum `[compat]` promises, the current
  stable, and the upcoming release. The pinned `'1.11'` was dropped, being neither.
- `test/Project.toml` and `docs/Project.toml` have `[compat]` sections, so a future release
  of Aqua, DataFrames or Documenter can no longer redden CI on its own.
- `documentation.yml` was aligned with `CI.yml`: `cache@v3`, a `concurrency` block,
  `timeout-minutes`, and the PkgTemplates boilerplate comment removed.

- macOS is in the matrix, on Apple Silicon, now that 2.3 has removed the MKL dependency
  that made the package unavailable there. The matrix is written as an explicit `include`
  list, so each row pairs its own architecture.

- **Coverage is collected and published.** The tests already ran instrumented, since
  `coverage` defaults to true in `julia-runtest`, so all that was missing was gathering the
  `.cov` files and sending them: `julia-actions/julia-processcoverage` then
  `codecov/codecov-action`, on each of the four matrix jobs, which Codecov merges. The badge
  is in the README. Read the list of uncovered lines rather than the percentage: a covered
  line was executed, which does not mean its result was checked.

**Windows is in the matrix now, to find out.** The claim that `test/checkpoint.jl` would
break there was never tried, and the two Julia 1.10 failures showed that a job nobody runs
is a job that hides things. The group is the only one touching the disk, so if anything
goes wrong it will be there, and what actually fails will be fixed rather than what was
imagined. The old wording, kept for the record: `test/checkpoint.jl` changes the working
directory and manipulates files, which is precisely what would break there, so adding it
means fixing the test first.

### 2.3 Two decisions — **both taken**

- **MKL was dropped as a dependency.** Its whole use was the `using MKL` of
  `src/TensorMixedStates.jl:10`, which switches the BLAS backend of the entire Julia
  session. The situation was worse than the review stated: the `MKL = "0.7 - 0.9"` bound
  resolves to MKL.jl 0.9, which requires MKL_jll 2025, and that generation ships `x86_64`
  Linux and Windows only. Intel dropped macOS from oneMKL, so the package was installable
  on no macOS at all, Intel included, and on no ARM machine.

  The weak dependency route of the review was examined and rejected: an extension only
  fires when the user loads MKL themselves, and at that point MKL has already switched the
  backend on its own, so the extension would be empty.

  What replaces it is a `BLAS backend` subsection in the Installation section of
  `docs/src/index.md`: OpenBLAS is used as it comes, and an `x86_64` user who wants MKL
  loads it before TMS. This is a user visible performance change for those users and
  belongs in the CHANGELOG of 2.4. Verified after the removal: the package loads, BLAS
  reports `libopenblas64_`, and `expect1` returns exact values on a two qubit state.

- **`julia = "1.10.5"`**, which Pkg reads as `[1.10.5, 2.0.0)`. The package no longer
  becomes uninstallable the day 1.14 ships, and the bound no longer has to be raised at
  every Julia release.

### 2.4 What a published package was missing — **done, except the coverage badge**

- **`CITATION.cff`** carries the SciPost DOI rather than the HAL preprint, with the codebase
  release `72-r1.28` as a second reference, so that GitHub's "Cite this repository" button
  appears and offers the right thing. Grégoire Misguich's ORCID comes from the Crossref
  record of the article; Jérôme Houdayer has none registered there, so the field is absent
  rather than invented. No email address, by choice.
- **`CHANGELOG.md`** starts at the release in preparation, which collects the fifty six
  commits made since v1.2.9 under Added, Changed, Fixed, Documentation and Development.
  Earlier versions are left to the tags: their commit messages were not written to be turned
  into a changelog, and reconstructing them would have produced plausible fiction.
- **`CONTRIBUTING.md`** gathers the process knowledge that lived in the wrong places: the
  test group convention and the selective run, the "what belongs here" file headers, the
  `TMS_SKIP_PRECOMPILE_WORKLOAD=true` shortcut deferred from 1.2, the documentation build
  and why a broken `@example` turns CI red, the `checkdocs = :exports` choice, the rule that
  site operators go through `@def_operators` and never `add_operator`, and two paragraphs
  saying where performance actually lies — the contractions, and the term count of the MPO
  rather than the time simplification takes.
- **The README** has the logo that already existed, a quick start that was run before being
  pasted, with its real output, the installation line, links to the two example folders, the
  published article, and the licence: GPL-3.0-or-later, which was the deliberate choice
  between the two readings, next to ITensor under Apache 2.0.
- **Issue forms** for bugs and for feature requests, links to the documentation and to
  `CONTRIBUTING.md`, and a pull request template whose checklist names the four things a
  change has to have done.

The coverage badge, which was waiting on the decision of 2.2, is there too.

---

## 3. Tests

- **The RNG is pinned — done.** `test/runtests.jl` seeds the global generator once, which
  covers `RandomState`, and the sampling tests of `test/observables.jl` pass a generator of
  their own to `sample`. That second part matters as much as the first: leaning on the
  global generator would have made the frequencies depend on how much randomness the groups
  before them drew, so running one group alone gave a different answer from running the
  suite. The tolerances of 4.5σ and 4.7σ were left as they are, since with a fixed seed the
  outcome is no longer a lottery. Note that a seed pins the stream for one Julia version,
  not across versions.
- **Restoring `exit_on_sigint` has no test.** It can only be checked by sending a real
  SIGINT from a child process, which was judged not worth its cost. It is the one fix of
  the correctness pass with no coverage.
- **Tolerances set at the floating point floor.** `test/evolve.jl` compared a `Tdvp` and an
  `ApproxW(order = 4, w = 2)` evolution to `cos(2t)` and `sin(2t)` with `1e-14`, which the
  Julia 1.10 job missed by `1.3e-14` — rounding noise that differs with the BLAS build, not
  an error of the algorithm. Both are now `1e-13`, which is still far tighter than any real
  error would be. Worth a look at the other bounds before adding one: a tolerance should say
  how accurate the method is, not how the machine happened to round that day.
- **Recorded reference values — two of the three replaced by exact ones.** They were
  compared to 11–15 significant digits with no word on where they came from, which froze the
  behaviour of the day they were recorded instead of checking it.

  `test/reference/` now holds two scripts that compute the references from first principles,
  in `LinearAlgebra` alone, sharing no code with what they check: `ising_ed.jl` diagonalizes
  the 64 dimensional Hilbert space of the 6 qubit ring, and `fermion_lindblad.jl` does the
  dense Lindblad evolution of the 32 dimensional Fock space of the 5 site chain, by hand
  Jordan-Wigner then the exponential of the vectorized Liouvillian. Neither is run by the
  suite; they exist so the numbers can be regenerated and argued with. The old recorded
  values turned out to be right, agreeing with the exact ones to 2e-8 and 8e-8, which is the
  error of the algorithms themselves and sits well inside the 1e-7 and 1e-6 tolerances. The
  testsets now carry the exact values and a comment saying where they come from, so the
  tolerance means what it should: the accuracy of the method, not the reproducibility of a
  past run.

  **Left open**: the free boson testset. Four sites of dimension 7 give 2401 Fock states and
  a vectorized Liouvillian of 2401², so a dense reference is out of reach. A reference built
  from the Gaussian moments of that quadratic Lindbladian would close it; until then the
  testset says in a comment that it is a regression check, which is at least honest.

---

## 4. API and design

These need decisions, not patches.

- **The implicit invariant on site operator names — fixed.** Every site module declared its
  operator constants in the same module, so `Sp`, `Sz`, `N`, `A`, `S` were bound two to four
  times, and it worked only because the returned `Operator{1}(name, nothing, type)` happened
  to be identical each time. A new site type registering one of those names with a different
  `OpType` would have silently won for every site type.

  The Julia 1.10 job of the new CI matrix showed the other face of the same design:
  `@def_operators` expanded to a `const NAME = add_operator(...)` in the calling module, so
  a user declaring an operator whose name a loaded site module exports hit
  `cannot assign a value to imported variable`, a hard error of the language up to Julia
  1.11. Julia 1.12 changed the binding rules and hid the problem on recent versions only.

  `@def_operators` now binds a name once. The decision is taken at expansion time, from
  `isdefined(__module__, sym)`: a name already in scope yields no binding at all, only a
  registration for the new site preceded by `check_shared_operator`, which refuses a name
  bound to something else or declared with another `OpType`. The check runs before
  `add_operator`, so a refused declaration leaves the operator library untouched. Both
  faces go away: the invariant is enforced instead of hoped for, and nothing is ever
  rebound, so the Julia 1.10 error cannot occur. Verified by loading the package —
  `Bosons.N === Fermions.N === Qudits.N === Qbosons.N` — and by `test/sites.jl`, which
  declares `N` for a site of its own and checks that the mismatching declaration is
  refused. Documented in the `Defining new site types` section of `sites.md`.
- **`Simulation` sharing — documented.** The copy constructor hands the new object the very
  `files`, `data` and `checkpoint` of the old one. That is deliberate: those are the parts
  that must not fork, and sharing them is what lets a copy made inside a phase advance the
  same checkpoint. It was also the root cause of the interrupted-phase bug, which is reason
  enough to say it out loud, and the docstring of `Simulation` now does.
- **`runTMS` is neither reentrant nor usable in parallel — documented.** It changes the
  working directory, sets `Base.exit_on_sigint`, and reseeds the global generator for a
  `CreateState` carrying a `seed`. Defensible for a script driver, but it was written
  nowhere; the docstring now says to use separate processes or `output`, and that
  parallelism inside one simulation is a different matter and works as usual.
- **`Dmrg` is really deprecated now**, through `Base.@deprecate_binding`, which is Julia's
  tool for a binding as opposed to `Base.@deprecate` for a function. `Mutual_Info_Renyi2`
  was renamed `MutualInfoRenyi2` in the same movement and its old spelling deprecated with
  `Base.@deprecate`, the label it writes to the output files following the new name.

  Three things worth remembering, all measured rather than assumed. Both macros export the
  old name themselves unless told not to, so `false` is passed and the export lists stay the
  single source of truth. A docstring cannot sit above either macro, since they expand to a
  toplevel block, so it is attached with `@doc` afterwards. And the warnings are invisible to
  ordinary users: `--depwarn` defaults to `no`, and only a test run turns it on.

  **Left open**: a deprecated *binding* warns only on a qualified access,
  `TensorMixedStates.Dmrg`. After `using TensorMixedStates`, which is how everyone writes it,
  `Dmrg` is silent even with `--depwarn=yes`. The function deprecation of
  `Mutual_Info_Renyi2` does not have this problem and warns as expected. Making `Dmrg` warn
  on use would mean turning it from a type alias into a function forwarding to the
  `GroundState` constructor, which would break `x isa Dmrg` and any use in type position.
- **Divergent defaults — aligned on the phase.** `approx_W` takes `order` with no default
  and `w = 2`, which is what the docstring of `ApproxW` described all along, so the
  documentation stops contradicting the code. The call in `Precompile.jl` and one test
  called `approx_W` without `order` and were fixed; no phase is affected, since `ApproxW`
  passes both explicitly. The change is in the CHANGELOG, being a break for direct callers.

  The divergent naming, `alg` on the function against `mpo_algo` on the phase, is left as it
  is on purpose: renaming a keyword is an API change, and the API is not changed for a
  naming blemish.
- **`run_phase` — done.** It has a docstring saying that it is the extension point, and two
  failures now say what is wrong instead of surfacing as a `MethodError` or a `FieldError`
  from the middle of a run: `check_is_phase` in `log_phase` catches an object that has none
  of the three fields a phase is read through, and the fallback `run_phase` catches one that
  has them but no method. Both messages name what to define.
- **A `SimData` is not a phase — refused.** It has the shape of one, with `name`,
  `time_start` and `final_measures`, only because `runTMS` runs the top level one through
  `log_phase` like any phase, which is where the first line of the log comes from. That made
  `run_phase(::Simulation, ::SimData)` reachable for a `SimData` sitting inside `phases`, and
  there it was worse than the numbering problem the review suspected: the loop it opened
  shared the phase counter of the loop around it, so it skipped every inner phase whose index
  was below the outer one. Reproduced without any checkpoint involved — the same two
  evolutions wrote twelve lines flat and six lines nested, the log showing one `Time
  evolution` where there should have been two. `flatten_phases` now refuses it and points at
  vector nesting, which is the documented way and works. `const Phases` was right all along.
- Minor: duplicated exports are gone, each of `Limits`, `mix`, `dag` and `tensor` being
  exported once now, from where it is defined or from the file included first. Left alone on
  purpose, because the API is not changed for a naming blemish: `removeMulti` in camelCase
  among `make_mpo`, `partial_trace`, `graph_base_size`; `trace2` and `Purity` are the same
  function under two names; no `show` for `JW`, `Multi_F`, `Proj`, `SetState`, and
  `test/operators.jl` depends on the default struct display.

---

## 5. Performance

One target is left, on the contraction side, where the time actually goes.

1. ~~**Memoise the local tensors in `expect2`**~~ — **set aside for now.** The item was
   ambiguous and the ambiguity is what mattered. Memoising the *environments*, or the
   transfer matrices per site, which is the usual way to speed up all-pairs correlators,
   costs χ⁴ per site for a transfer matrix carrying four link indices: out of the question
   at publication sizes. What the item actually pointed at is narrower, the
   `tensor_obs(state, o(i))` calls inside the `map(ops)` of the inner loop. Those tensors
   are built by `tensor(system, ::AtIndex)` in `src/Systems.jl`, whose indices are the site
   indices and their primes — no link index enters, so a `(site, operator)` cache would hold
   n·|ops| tensors of d² entries, a few hundred bytes on a chain of forty qubits. It would
   remove O(n²·|ops|) constructions, but those are dictionary lookups and small matrix work
   sitting next to ITensor contractions in the same loop, so the gain would not show in a
   measurement. Left undone deliberately rather than forgotten.
2. **The MPO is never compressed — kept as it is, and now documented.** `PreMPO!` gives
   each term a channel of its own, so the bond dimension on a link is
   `2 + #(terms crossing it)`. Measured on `sum(Z(i)Z(j))` over all pairs: 27, 102 and 402
   for 10, 20 and 40 sites, against 3 in every case for `ITensorMPS`' `OpSum`, which
   compresses with an SVD. On a short ranged operator there is no difference at all, the
   Ising chain giving 3 on both sides whatever the length.

   Compressing was set aside: the triangular form the construction leaves is exactly what
   WI and WII need, and keeping each term separate is what lets a coefficient change from
   sweep to sweep, which is what makes time dependent evolvers free. What was missing was
   that none of this was written anywhere, so a user coming from `ITensorMPS` met the
   difference without an explanation. The `MPO` section of `others.md` now carries it, with
   the formula, the measured table and the reason.

The `Matrix{Any}` of `expect2` was looked at and closed. It never reaches the caller:
`expect2` returns `Matrix{Float64}` for one pair of operators and `Vector{Matrix{Float64}}`
for several, because `unroll` rebuilds a concretely typed array from the values. The
container is an internal intermediate whose cost is n² boxed stores against n² ITensor
contractions in the same loop, and typing it ahead would mean computing a cell first or
going through `promote_op`, for no measurable gain. All that was changed is that
`Matrix{Any}(undef, n, n)` now says what it is.

The single pass simplification of generic products is **done**, and the review had it
wrong about why it mattered. `simplify(X*Y*Y*X)` did give `X*X`: `simplify_core_prod` kept
one current base and pushed each finished run into its result, so a run collapsing to the
identity left its two neighbours adjacent without anyone looking back. Each call peeled one
layer, `X*Y*Z*Z*Y*X` giving `X*Y*Y*X`.

The claimed cost in MPO bond dimension does **not** reproduce, though. An MPO is built from
indexed operators, and the indexed path already merged to a fixed point: it compares each
factor with the last one kept, in a `while change` loop, which is the bubble sort its
comments mention. Measured on four operators, including one spanning three sites,
`maxlinkdim` was the minimum in every case. What was real was that `simplify` was not
idempotent on generic expressions, its result depending on how many times it was called.

The fix takes the shape of the indexed loop, minus the reordering, which has no meaning
between two factors on one site since they do not commute.

Lesser points, all measurable but small: scalar `setindex` filling of the MPO tensors,
abstract field types in the hot symbolic layer, and `stop_requested` doing an `isfile` on
every sweep.

---

## 6. Examples

The six scripts of `examples/article/` and the five of `examples/high_level/` are all up to
date with the current API and abundantly commented, but:

- **none is run by CI**, and none uses the recent features (`checkpoint_interval`,
  `max_time`, `square_lattice`, `Qudit`, time dependent evolvers);
- ~~`ising_quench.jl` writes a description that lies~~ — **done**. It said `maxdim = 100`
  where `limits` set 50, and that text goes to the run's `description` file, so the
  simulation recorded a false parameter. Rather than correct the number, the description is
  now built from the values themselves, `$(limits.maxdim)` and the like, with the algorithm
  and the time step hoisted into variables so they have one definition. It cannot drift
  again, and it gained what it used to omit: the rendered algorithm reads
  `ApproxW(order = 4, w = 2, n_hermitianize = 5)`;
- ~~a space in a filename, and simulation names with spaces~~ — **done**.
  `4_Free_fermions_with_source.jl` matches its five siblings, and the three `SimData` names
  carrying spaces, in `gates.jl`, `complete_graph_tdvp.jl` and `precession.jl`, no longer
  do, since that name becomes a directory. The `name` of a phase is only a log label and
  was left as it reads;
- ~~a comment contradicting the evolver~~ — **done**, by explaining rather than by changing
  the code, which was right: the hamiltonian of a tight binding chain is minus its hopping
  sum, so `-im * H` is `+im` times that sum, which is what the line writes;
- there is no README in `examples/high_level/`, and `examples/article/README` should be
  `README.md` to be rendered by GitHub;
- no indication of running time anywhere, and these are publication sizes;
- **`docs/src/index.md` says "Working examples are presented in the folder `examples`"
  with no hyperlink**, no mention of the two subfolders and no word on what each shows. A
  reader of the published site has no path to the eleven scripts.
