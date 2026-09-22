# Contributing to TensorMixedStates

Contributions are welcome, from a one line typo fix to a new site type or a new algorithm.
This file gathers the process knowledge that is otherwise scattered in the repository.

## Getting set up

```
git clone https://github.com/jerhoud/TensorMixedStates.jl
cd TensorMixedStates.jl
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Development happens on the `dev` branch and `main` carries the released state. Which branch
a given change belongs on, and what else it needs, is the subject of the next section.

## Where a change goes

`dev` is where work lands and `main` carries the released state, together with the tag that
names it. Three facts about this repository decide the rest.

The documentation is published from `main` and from tags only. A push to `dev` builds it and
runs its doctests and `@example` blocks, which is the check that matters, but publishes
nothing: `main` deploys the `dev` documentation and a tag deploys `vX.Y.Z` and moves
`stable` onto it. A documentation fix is therefore verified as soon as it reaches `dev`, and
visible only once it reaches `main`.

Continuous integration runs on `main`, on `dev`, on tags and on pull requests, and otherwise
only when dispatched by hand. A topic branch pushed on its own runs nothing, so open a pull
request against `dev`, a draft one if the work is unfinished, to have the tests run on it.

`main` carrying the tag means that any commit put on it moves the branch past the released
version. Only a release does that.

| Kind of change | Branch | Changelog | Version |
|---|---|---|---|
| A typo or a docstring | `dev` | no | — |
| A documentation rework | topic branch off `dev` | yes, *Changed* | patch |
| A bug fix | `dev` | yes, *Fixed* | patch |
| An urgent fix of a serious bug | see below | yes, *Fixed* | patch, released at once |
| A small improvement | `dev` | yes, *Added* or *Changed* | patch, minor if it adds an exported name |
| A large feature | topic branch off `dev` | yes, *Added*, written at the end | minor |

Documentation alone does not earn a changelog entry, which is what the "user visible" of the
[pull request template](.github/pull_request_template.md) comes to in practice. Without push
rights, every line of that table is a pull request against `dev`, and that template lists
what to check before opening one.

Some of those lines need more than a row.

**A long lived branch stays off `dev`**, whether it carries a documentation rework or a
feature. Cut it from `dev` and merge it back there. Only an urgent fix ever branches from
the tag. `dev` has to remain mergeable into `main` at any moment, because that is what
makes an urgent release possible; a half finished rewrite sitting there takes it away.

Commit often on it; a commit is already a state you can return to. A branch of its own is
for what a commit does not protect: `git branch before-rebase` ahead of a rebase or a
reset, which leave the old commits reachable only from the reflog, and a second branch when
two approaches have to live side by side. Rebase on `dev` now and then so that the final
merge stays small, bearing in mind that rebasing a branch already pushed rewrites it and
needs `--force-with-lease`.

**An urgent fix** depends on what `dev` holds. If it holds nothing you would refuse to
publish, fix it there and fast forward `main` onto it. If it holds unfinished work, branch
from the tag instead, merge into `main`, release, and merge `main` back into `dev`:

```
git switch -c hotfix/short-name vX.Y.Z     # the released tag
# the fix, plus the regression test that would have caught it
git switch main && git merge hotfix/short-name
# bump the version in Project.toml, commit, push, release
git switch dev && git merge main
```

The last line is the one that is easy to forget. Without it `dev` loses the fix and stops
being a fast forward of `main`, which complicates every merge afterwards.

### Releasing

This is the maintainer's part, and it is the same for every version.

1. `main` is fast forwarded onto `dev`.
2. The `version` field of `Project.toml` is bumped on `main`, committed and pushed, which
   publishes the `dev` documentation.
3. That commit is registered in the General registry through JuliaRegistrator.
4. TagBot creates the `vX.Y.Z` tag once the registry pull request is merged.
5. The tag deploys `vX.Y.Z` and moves `stable` onto it.

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
`test/reference/` holds scripts that produce the reference values of `algorithms.jl` by
exact diagonalization and by dense Lindblad evolution, in `LinearAlgebra` alone. They share
no code with what they check, which is the whole point, and the suite does not run them —
they are there so the numbers in the tests can be regenerated and argued with. If you must
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

Site operators are declared with `@def_operators`, never by calling `add_operator` directly.
The macro is what binds the operator name, checks it against a name already in scope and
keeps the library consistent; calling `add_operator` yourself goes around all of that. The
same holds for `@def_states` and `@create_site_module`.

Global identifiers are `const`.

Control flow is written with `if`, never with the short-circuit operators. `cond && return`,
`x isa T || error(…)` and `flag && do_something()` are out, and so is every other use of
`&&` or `||` for its side effect. Inside a boolean expression, as in `if a && b`, they are
ordinary and welcome.

## Performance work

Two things are worth knowing before optimising anything.

The running time is dominated by the ITensor contractions, not by the Julia code around
them. Rewriting an inner loop, shaving an allocation or chasing a type instability in the
symbolic layer will not show up in a measurement. The contraction side will.

When a change touches the symbolic layer or the construction of the MPO, the measure of
success is the number of terms in the resulting MPO, not the time simplification takes. A
simplification that runs twice as fast but leaves more terms behind is a loss, because the
term count is what the bond dimension is paid on, and the bond dimension is what the whole
computation is paid on.

One design point that looks like an oversight and is not: applying gates deliberately
bypasses `simplify`, and only the MPO path simplifies.

## Reporting a problem

Open an issue with the smallest piece of code that reproduces what you see, the version of
TMS and the version of Julia. If what you have is a question rather than a bug, or if you
want to say what you are doing with TMS and what you would need from it, the contact address
is given in the [reference article](https://doi.org/10.21468/SciPostPhysCodeb.72).

## Licence

TensorMixedStates is distributed under the GNU General Public License, version 3 or later.
By contributing you agree that your contribution is distributed under the same terms.
