# Threads and performance

Most of the running time goes into the tensor contractions of ITensors. How they use the cores
of the machine comes down to three choices.

**Starting Julia.** Start it with as many threads as the machine has:

```sh
julia --threads=auto my_script.jl
```

or `julia -t auto my_script.jl` for short; `--threads=4`, or `-t 4`, gives exactly four. The
garbage collector then runs on as many threads from Julia 1.12 on, and on half as many
before, which speeds up the runs, and they are what the `:blocks` mode below runs on. Started
without them, Julia still runs the products of dense matrices on several cores, BLAS having
threads of its own, but everything else on a single one. The threads of Julia are fixed when
it starts: they cannot be added from within a program.

**The BLAS library.** On `x86_64` machines running Linux or Windows, MKL, which replaces
OpenBLAS for the products of dense matrices, is often faster, see [BLAS backend](@ref). This
choice is independent of the two others.

**The mode.** The `threading` keyword of `runTMS`, or the environment variable
`TMS_THREADING`, or `set_threading`, chooses how the contractions use the threads of Julia:

- `:dense`: BLAS runs each product of matrices on several threads, and Strided, which ITensors
  uses for the permutations of dense tensors, runs on a single one, as ITensors recommends.
  Julia starts Strided on as many threads as it has itself, where they compete with those of
  BLAS. `runTMS` applies this mode by default; a program calling the functions of TMS directly
  should start with `set_threading(:dense)`;
- `:blocks` asks ITensors to run the products of the blocks of block sparse tensors in parallel
  on the threads of Julia, BLAS running each of them on a single thread. The tensors of a
  system are block sparse when it conserves something, see [Conserving a quantity](@ref), so
  this mode is meant for such systems: it is worth trying on them, with OpenBLAS above all;
- `:auto`, for `runTMS` only, chooses before each phase: `:blocks` when the system of the
  state conserves something and Julia has several threads, `:dense` otherwise and until there
  is a state.
  `set_threading(mysystem)` makes the same choice once.

```julia
using MKL                      # if MKL is installed, on x86_64 Linux or Windows: before TMS
using TensorMixedStates
set_threading(:dense)          # first thing in a program calling the functions of TMS directly

old = set_threading(:blocks)   # set_threading returns the threading it replaces,
set_threading(old)             # which it takes back

runTMS(sim_data; threading = :blocks)
```

The threading `set_threading` returns is a `ThreadingState`, which gives the threads of BLAS and
of Strided and whether block sparse multithreading is on. One can be built by hand to try other
settings: `set_threading(ThreadingState(blas = 2))` puts BLAS on two threads and leaves the rest
as it is.

The `stamp` file of a simulation records all these settings, so that the running times of two
runs can be compared.

On a laptop with four cores, an Intel i7-8650U, for two sweeps of `Dmrg` and two steps of
`Tdvp` on a Heisenberg chain of 20 sites, at bond dimensions 64, 256 and 768, the running times
changed as follows, each setting run once, two runs of the same setting differing by up to 6 %:

| | dense tensors | block sparse tensors |
|:---|:---|:---|
| `--threads=auto` with `:dense`, rather than `julia` alone | −22 to +4 % | −25 to −5 % |
| MKL rather than OpenBLAS | −18 to −2 % | −6 to +6 % |
| Strided left on the threads of Julia, rather than on one as `:dense` puts it | −9 to +6 % | −7 to +8 % |
| `:blocks` rather than `:dense`, with OpenBLAS | −7 to +60 % | −9 to +47 % |
| `:blocks` rather than `:dense`, with MKL | +43 to +71 % | +19 to +48 % |

The threads of Julia gain from the bond dimension 256 on, MKL most at the smallest, 64, and
`:blocks` costs most on dense tensors; leaving Strided on the threads of Julia changes nothing
beyond the scatter of the runs.

The documentations of Julia and of ITensors agree that the way to find the best settings is to
try them on a few sweeps of your own calculation.

## BLAS backend

Most of the running time of TMS is spent in the tensor contractions of ITensors, which
themselves call BLAS. Julia ships with OpenBLAS and TMS uses it as it comes: switching the
BLAS backend affects the whole Julia session, so that choice is left to you rather than
made by a library you load.

On `x86_64` machines running Linux or Windows, Intel's MKL is often noticeably faster on
these contractions. To use it, add `MKL` to your project and load it before TMS:

```julia
using MKL
using TensorMixedStates
```

MKL is not distributed for macOS nor for ARM machines, where OpenBLAS is the only option.

## Long sums

A sum written with a generator, `sum(Z(i) * Z(j) for (i, j) in g)`, is built term by term,
each term copying the sum so far, so that its cost grows as the square of the number of terms.
Written with brackets, `sum([Z(i) * Z(j) for (i, j) in g])`, it is the sum of a vector, which
Julia adds by halves: the 79800 terms of the pairs of 400 sites take 0.6 s rather than 18 s.
Below a few thousand terms, the two are alike.
