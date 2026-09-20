# Contributing to TensorMixedStates

Contributions are welcome, from a one line typo fix to a new site type or a new algorithm.
This file gathers the process knowledge that is otherwise scattered in the repository.

## Getting set up

```
git clone https://github.com/jerhoud/TensorMixedStates.jl
cd TensorMixedStates.jl
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Development happens on the `dev` branch, `main` carries the released state, and continuous
integration runs on both.

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
