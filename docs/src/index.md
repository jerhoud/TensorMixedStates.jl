# Home

TensorMixedStates (TMS) is a Julia library to make simulations of closed or open quantum systems using Matrix Product States representations.

## Features

TMS uses Matrix Product State representations for the density matrix of the system. It proposes a large set of features:
manipulations of systems and states, a rich set of sites and operators easily extensible by the user, powerful algorithms: computation of ground states using DMRG, Hamiltonian and Lindbladian evolution with TDVP (and others), applications of gates (including noisy gates).

Being based on ITensor, TMS delivers high performance computations and naturally runs in parallel.

The interface is user friendly. In particular, it features a very expressive syntax for operators allowing easy definitions of operators such as observables, gates, Hamiltonians or Lindbladians. Moreover the optional high level interface allows the writing of simple simulations in a few lines of code.

## References

TMS is described in the following article, published in SciPost Physics Codebases. Please
cite both the article and the codebase release it documents, which is the convention
SciPost asks for:

Jérôme Houdayer and Grégoire Misguich, *TensorMixedStates: A Julia library for simulating
pure and mixed quantum states using matrix product states*,
[SciPost Phys. Codebases **72** (2026)](https://doi.org/10.21468/SciPostPhysCodeb.72).

Jérôme Houdayer and Grégoire Misguich, *Codebase release 1.28 for TensorMixedStates*,
[SciPost Phys. Codebases **72-r1.28** (2026)](https://doi.org/10.21468/SciPostPhysCodeb.72-r1.28).

In BibTeX:

```bibtex
@article{TensorMixedStates,
  title     = {{TensorMixedStates}: A {Julia} library for simulating pure and mixed
               quantum states using matrix product states},
  author    = {Houdayer, J{\'e}r{\^o}me and Misguich, Gr{\'e}goire},
  journal   = {SciPost Phys. Codebases},
  pages     = {72},
  year      = {2026},
  publisher = {SciPost},
  doi       = {10.21468/SciPostPhysCodeb.72},
}

@article{TensorMixedStatesRelease,
  title     = {{Codebase} release 1.28 for {TensorMixedStates}},
  author    = {Houdayer, J{\'e}r{\^o}me and Misguich, Gr{\'e}goire},
  journal   = {SciPost Phys. Codebases},
  pages     = {72-r1.28},
  year      = {2026},
  publisher = {SciPost},
  doi       = {10.21468/SciPostPhysCodeb.72-r1.28},
}
```

### Acknowledgements and feedback

If TMS helped produce the data of a publication, a word about it in your acknowledgements
would be warmly appreciated. It costs you a single line, and it is what makes it possible
to keep the library developed and maintained.

We would also be delighted to hear from you directly, and you should not hesitate for a
moment. What are you doing with TMS? What works well, what gets in your way, what do you
need that is not there yet? The contact address is the one given in the article. Such
messages are read with care, and they are what shapes what comes next.

## Installation

To use TMS, you need to have Julia installed on your system. Installing julia is usually easy and fast, see [The Julia Programming Language](https://julialang.org/) for instructions. TMS requires at least Julia version 1.10.5 to run.

To install TMS in Julia, launch the Julia interface (by typing 'julia' on the command line) and type

```
]add TensorMixedStates
```

Note the "]" required to enter julia package management system.

After downloading TMS, Julia will automatically compile and install it. This process usually takes a couple of minutes and does not require interactions on the user part.

### BLAS backend

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

## Using TMS

To use TMS in your code you need to write

```julia
using TensorMixedStates
```

Julia script names are usually written with a .jl extension. Once you have written your script, you can execute it with

```sh
julia my_script.jl
```

to run it on a single processor, or

```sh
julia --threads=4 my_script.jl
julia --threads=auto my_script.jl
```

to use multi-threading (see julia documentation for more details on multi-threading).

You can also use TMS in the interactive julia interpreter.

## Bug reports and feature requests

TMS is still in development and certainly contains bugs. If you think you have found one please report it on the
[Github page](https://github.com/jerhoud/TensorMixedStates.jl).

Feature requests may also be sent on the [Github page](https://github.com/jerhoud/TensorMixedStates.jl).

## Documentation

You are currently reading it!

You can have access to inline documentation on TMS at the julia prompt simply by typing "?" followed by the function name or type name you are interested in. For example

```
?runTMS
```

For this to work you must have first imported TMS with

```julia
using TensorMixedStates
```

## Examples

Working examples are presented in the folder `examples` in the repository. 

## Module

```@docs
TensorMixedStates
```

## Index

```@index
```
