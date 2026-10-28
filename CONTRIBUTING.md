# Contributing to TensorMixedStates

Contributions are welcome, from a one line typo fix to a new site type or a new algorithm.
This file gathers the process knowledge that is otherwise scattered in the repository.

## Getting set up

```
git clone https://github.com/jerhoud/TensorMixedStates.jl
cd TensorMixedStates.jl
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Fixes of the released version are made on `main`, everything else on `dev`, where the next
version is prepared. Which branch a given change belongs on, and what else it needs, is the
subject of the next section.

## Where a change goes

`main` carries the released series: the last release, and the fixes that will make its next
patch release. `dev` is where the next minor version is prepared: new features, changes of
behaviour, and every fix of `main`, merged from there. Three facts about this repository
decide the rest.

The documentation is published from `dev` and from tags. A push to `dev` builds it, runs its
doctests and `@example` blocks, which is the check that matters, and publishes it as the `dev`
documentation; a tag publishes `vX.Y.Z` and moves `stable` onto it; a push to `main` builds
it and publishes nothing. A documentation fix is therefore visible in the `dev` documentation
as soon as it reaches `dev`, and in `stable` at the next release.

Continuous integration runs on `main`, on `dev`, on tags and on pull requests, and otherwise
only when dispatched by hand. A topic branch pushed on its own runs nothing, so open a pull
request, a draft one if the work is unfinished, to have the tests run on it: against `main`
for a fix of the released version, against `dev` for anything else.

Changes flow one way, from `main` into `dev`. Every fix that lands on `main` is merged into
`dev` at once, so that `dev` always contains `main`. Nothing goes from `dev` to `main` but a
minor or a major release, replayed on `main`, whose history is kept linear, see *Releasing*.

A fix makes the package do what its documentation says, or what it plainly should, without
adding to its interface. A new function, keyword or option, a change of documented behaviour,
a raised lower bound of compat, a rewrite or an optimisation are not fixes, whatever bug they
also cure, and go on `dev`. Widening a compat bound, as the pull requests of CompatHelper do
against `main`, is a fix. When in doubt, choose `dev`: a fix released with the next minor
version does less harm than a patch release that breaks something.

| Kind of change | Branch | Changelog | Released in |
|---|---|---|---|
| A typo or a docstring | `main`, or `dev` for what only `dev` has | no | — |
| A bug fix | `main`, or a topic branch off it | yes, *Fixed* | the next patch release |
| An urgent fix of a serious bug | `main` | yes, *Fixed* | a patch release, at once |
| A bug of what only `dev` has | `dev` | no | — |
| A documentation rework | topic branch off `dev` | yes, *Changed* | the next minor release |
| A small improvement, its interface settled | `dev` | yes, *Added* or *Changed* | the next minor release |
| A large feature, or an interface that may still change | topic branch off `dev` | yes, *Added*, written at the end | the next minor release |

Documentation alone does not earn a changelog entry, which is what the "user visible" of the
[pull request template](.github/pull_request_template.md) comes to in practice, and neither
does a bug of what no release has had. Without push rights, every line of that table is a
pull request against the branch it names, and that template lists what to check before
opening one.

Some of those lines need more than a row.

**A fix** is made on `main` and merged into `dev` right away:

```
git switch main            # or a topic branch off it, merged by a pull request
# the fix, the regression test that would have caught it, its line under Fixed
git switch dev && git merge main
```

The last line is the one that is easy to forget. Without it `dev` loses the fix and stops
containing `main`, which complicates every merge afterwards. The merge conflicts in
`CHANGELOG.md` whenever `dev` has entries of its own: keep both, the fix under *Fixed* in the
`Unreleased` section of `dev`. A second working tree on `main`, made once by
`git worktree add ../TensorMixedStates-main main`, saves stashing the work in progress on
`dev` to make a fix.

**An urgent fix** needs nothing more: `main` holding only releases and fixes, it can be
released at any moment, as below.

**A long lived branch stays off `dev`**, whether it carries a documentation rework or a
feature. Cut it from `dev` and merge it back there once what it brings could be published as
it is: tested, documented, and with an interface you are ready to keep until the next major
version, since a published name cannot be taken back before then. `dev` has to remain
releasable as the next minor version at any moment; a half finished rewrite sitting there
takes that away.

Commit often on it; a commit is already a state you can return to. A branch of its own is
for what a commit does not protect: `git branch before-rebase` ahead of a rebase or a
reset, which leave the old commits reachable only from the reflog, and a second branch when
two approaches have to live side by side. Rebase on `dev` now and then so that the final
merge stays small, bearing in mind that rebasing a branch already pushed rewrites it and
needs `--force-with-lease`.

### When to release

A patch release costs users nothing; a minor release commits the package, every name and
keyword it publishes staying until the next major version.

A fix is released according to what it fixes. A wrong result given without warning, or an
error that blocks a common use, is released within days, since someone may publish a wrong
computation in the meantime. Other fixes, of an error message, a corner case or the
documentation, wait for the next urgent one or until a few have gathered, but not for more
than a month or two: a fix left unreleased brings reports of bugs already fixed.

A minor release gathers a coherent set of features rather than one each. It is made when they
are finished, tested, documented and used for a while in real computations, which is where a
wrong name or signature shows up while it can still be changed, or when someone needs one of
them. A few months between two is a reasonable pace.

A major release is for incompatible changes only, gathered and rare, announced where possible
by deprecation warnings in a minor release before.

A release that is due is not held back by something on `dev` that is not ready. What is
unfinished is kept out of the public interface for that release, neither exported nor
documented, or documented as experimental, since only those names commit the package. A
change that cannot be hidden, such as a change of behaviour, is reverted on `dev` before the
version commit and restored right after, by reverting the revert. A feature is never released
alone from `main` by cherry picking it: `main` would then carry a feature, and the two
branches would diverge.

### Releasing

This is the maintainer's part. A patch release is made from `main`, a minor or a major one
from what `dev` holds, by the same version commit: it sets the version in `Project.toml`, in
`CITATION.cff` with its `date-released`, and in the placeholder of
`.github/ISSUE_TEMPLATE/bug_report.yml`, and closes the changelog, its `Unreleased` section
becoming that of the version, dated like `CITATION.cff`, the day of the registration.

A patch release:

1. On `main`, the version commit, pushed.
2. A comment `@JuliaRegistrator register` on that commit, followed by the release notes,
   registers it in the General registry.
3. TagBot creates the `vX.Y.Z` tag once the registry pull request is merged.
4. The tag deploys `vX.Y.Z` and moves `stable` onto it.
5. `main` is merged into `dev`, keeping the version of `dev` in `Project.toml`. The entries
   the release took leave the `Unreleased` section of `dev` for the section of the release,
   which the merge brings.

A minor or a major release goes through a branch of its own, `prepare-X.Y.0`: the rule
"Protect main" of the repository keeps the history of `main` linear, and `dev` holds the
merges of `main`, so its commits are replayed on `main` rather than fast forwarded.

1. After a last merge of `main` into `dev` if a fix is missing there,
   `git switch -c prepare-X.Y.0 dev`, then `git rebase main`, which replays the commits of
   `dev` on `main` and leaves the merges out. The conflicts the merges resolved come back,
   and `git rerere`, enabled once by `git config --global rerere.enabled true`, replays the
   resolutions it recorded then.
2. `git diff dev prepare-X.Y.0` is empty: the branch holds what `dev` holds, as tested.
3. On the branch, the version commit. A pull request against `main` has CI run on it; it is
   not merged from GitHub, whose buttons would add a merge commit or rewrite the commits,
   and closes as merged at the next step.
4. `main` is fast forwarded onto the branch and pushed, then registered, tagged and
   documented as for a patch release, steps 2 to 4. The branch is deleted.
5. `dev` restarts from `main`: `git switch dev && git reset --hard main`, a commit giving it
   the version of the next series, and `git push --force-with-lease origin dev`. This is the
   only time `dev` is rewritten, when it holds nothing that `main` does not.

Between two releases `dev` carries the version it prepares with the suffix `-DEV`, as
`2.1.0-DEV`, so that a checkout of it is never taken for a release.

A fix of an older series, once `main` has moved on to a newer one, goes on a branch
`release-X.Y` cut from the last tag of that series, the fix and then the version commit as
above. CI runs there only when dispatched by hand, the release is registered by
`@JuliaRegistrator register branch=release-X.Y`, and the changelog of `dev` takes the
section of that release word for word, as it did for 1.6.1 and 1.6.2.

Tags and `gh-pages` do not refresh on their own in a clone, so `git fetch --tags` before
judging what is released.

## Running the tests

The whole suite:

```
julia --project=. -e 'using Pkg; Pkg.test()'
```

It is split into groups, listed in `GROUPS` at the top of `test/runtests.jl`. Each group
lives in the file of the same name, and the header comment of that file says what belongs
there. Read that comment before adding a test, and put your test where it belongs rather
than at the end of whichever file you happen to have open. To add a group, create the file
and add its name to `GROUPS`.

The whole suite takes a while, so while you work run only the groups your change touches:

```
julia --project=. -e 'using Pkg; Pkg.test(test_args = ["observables", "sites"])'
```

Set `TMS_SKIP_PRECOMPILE_WORKLOAD=true` in the environment to skip the precompilation
workload of `src/Precompile.jl`. That workload exists to make the first call fast for users,
and it costs a great deal of time when you are loading or testing the package over and over.
It is a development shortcut only: do not set it in CI and do not rely on it anywhere a user
could end up.

Coverage is measured on every CI run and sent to Codecov, which shows the lines of `src/`
that no test reaches. A pull request that adds code without tests will show up there. Note
that a covered line only means it was executed, not that its result was checked, so read the
list of uncovered lines rather than the percentage.

Compare against an exact value whenever an exact value exists. A number recorded from a
previous run freezes the behaviour of the day it was recorded rather than checking anything.
When no closed form is available, an independent computation often still is:
`test/reference/` holds scripts that produce the reference values of `evolve.jl` by exact
diagonalization, by dense Lindblad evolution and from Gaussian moments, in `LinearAlgebra`
alone, and the files of earlier versions that `states_io.jl` reads back. They share no code
with what they check, which is the whole point, and the suite does not run them — they are
there so the numbers in the tests can be regenerated and argued with. If you must
record a number, say in a comment where it comes from and that it is a regression check.

Randomness is pinned: `runtests.jl` seeds the global generator, and a test that draws should
pass a generator of its own rather than lean on that seed, so that running one group alone
gives the same answer as running the whole suite.

## Building the documentation

```
julia --project=docs docs/make.jl
```

Two things to know. The `@example` blocks of the manual and of the reference pages are
executed at every build and the build fails if one of them raises, so a broken example turns
CI red. Keep your examples runnable and use a plain fenced `julia` block for genuine
pseudo-code. And `checkdocs = :exports` is a deliberate choice: the published reference
documents what a user calls. The concrete operator subtypes, the `Base` overloads and the
machinery behind the phases have no `@docs` entry on purpose.

## Conventions

Code, comments, docstrings and commit messages are written in English.

Site operators are declared with `@def_operators`, never by adding a method to
`operator_definition` yourself. The macro is what binds the operator name, checks it against a
name already in scope and keeps the library consistent; a method written by hand goes around
all of that. The same holds for `@def_states` and `@create_site_module`.

Global identifiers are `const`.

Fermions rest on one invariant, the Jordan-Wigner convention: a fermionic operator `A` placed
on site `i` is its matrix there with `F` on every site before it, so that `A(i) * B(j)`, with
`i < j`, is `(A * F)(i) * F(i+1) * … * F(j-1) * B(j)` on the matrices of the sites. The
descent keeps the strings of a term apart, as the sites where they change (`Term`), and these
places lay them and must agree: `laid_strings`, for the dense matrices of terms and the dense
blocks of `apply`; `string_factors` and `direct_pass`, for the MPOs, whose channels carry the
strings over the sites a term crosses; `strung_kron`, for the matrix of a tensor product;
`string_signs`, for a matrix of several sites split into terms; `expect2`, along its
correlation matrix; and `trace_signs`, for `partial_trace`. The testset "A random fermionic
oracle" of `test/observables.jl` checks them against dense matrices written from the
convention alone. A new place laying strings goes through `laid_strings`, or is added to this
list and to that testset.

Control flow is written with `if`, never with the short-circuit operators. `cond && return`,
`x isa T || error(…)` and `flag && do_something()` are out, and so is every other use of
`&&` or `||` for its side effect. Inside a boolean expression, as in `if a && b`, they are
ordinary and welcome.

A public name takes the form most of the package has. Version 2.0 brought the names published
before into line with the rules below, but for the exception noted; from then on a published
name is kept as it is, since renaming it would break programs for no gain:

- words are written out: `measurements`, `log_message`, the tolerance `tol` excepted;
- a predicate is one word, as in Base: `stopped`;
- a measurement takes the state first, as `expect` and `variance`, and a solver the operator
  first, as in ITensorMPS: `evolve(evolver, duration, state)`, `ground_state(hamiltonian,
  state)`;
- a number of time steps is `nsteps`, a number of dmrg sweeps `nsweeps`;
- something done every `k` steps is given by `<what>_period`, as `measurements_period`;
- a field does not repeat the name of its struct: `Krylov(dim = 8)`, not `krylovdim`;
- the suffix `Error` is kept for an exception, a measurement of how far a state is from a
  property being a `Deviation`, as `TraceDeviation`;
- a constructor returns an object of its own type, a function making one of another type
  being written in lower case, as `random_state`.

## Performance work

Two things are worth knowing before optimising anything.

The running time is dominated by the ITensor contractions, not by the Julia code around
them. Rewriting an inner loop, shaving an allocation or chasing a type instability in the
symbolic layer will not show up in a measurement. The contraction side will.

When a change touches the descent of operators or the compaction into MPOs, the measure of
success is the number of channels of the resulting MPO, its bond dimension, which
`test/compact.jl` compares with the least a triangular MPO can have, not the time the descent
takes. A descent that runs twice as fast but leaves more channels behind is a loss, because
the bond dimension is what the whole computation is paid on.

One design point that looks like an oversight and is not: `apply` never develops a product.
It takes it factor by factor from the right, each as a dense block on its sites when small
enough and through its MPO otherwise, see `pieces` and `piece_steps!`, so that a product of
`Swap`s stays a sequence of gates rather than the sum it would develop into.

## Reporting a problem

Open an issue with the smallest piece of code that reproduces what you see, the version of
TMS and the version of Julia. If what you have is a question rather than a bug, or if you
want to say what you are doing with TMS and what you would need from it, the contact address
is given in the [reference article](https://doi.org/10.21468/SciPostPhysCodeb.72).

## Licence

TensorMixedStates is distributed under the GNU General Public License, version 3 or later.
By contributing you agree that your contribution is distributed under the same terms.
