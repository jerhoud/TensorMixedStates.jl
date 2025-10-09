# Home

TensorMixedStates (TMS) is a Julia library to make simulations of closed or open quantum systems using Matrix Product States representations.

## Features

TMS uses Matrix Product State representations for the density matrix of the system. It proposes a large set of features:
manipulations of systems and states, a rich set of sites and operators easily extensible by the user, powerful algorithms: computation of ground states using DMRG, Hamiltonian and Lindbladian evolution with TDVP (and others), applications of gates (including noisy gates).

Being based on ITensor, TMS delivers high performance computations and naturally runs in parallel.

The interface is user friendly. In particular, it features a very expressive syntax for operators allowing easy definitions of operators such as observables, gates, Hamiltonians or Lindbladians. Moreover the optional high level interface allows the writing of simple simulations in a few lines of code.

## References

To cite this software, please cite the following reference article

[TensorMixedStates: A Julia library for simulating pure and mixed quantum states using matrix product states](https://hal.science/hal-04945872)

## Installation

To use TMS, you need to have Julia installed on your system. Installing julia is usually esay and fast, see [The Julia Programming Language](https://julialang.org/) for instructions. TMS requires at least Julia version 1.10.5 to run.

To install TMS in Julia, launch the Julia interface (by typing 'julia' on the command line) and type

    ]add TensorMixedStates

Note the "]" required to enter julia package management system.

After downloading TMS, Julia will automatically compile and install it. This process usually takes a couple of minutes and does not require interactions on the user part.

## Using TMS

To use TMS in your code you need to write

    using TensorMixedStates

Julia script names are usually written with a .jl extension. Once you have written your script, you can execute it with

    julia my_script.jl

to run it on a single processor, or
    
    julia --threads=4 my_script.jl
    julia --threads=auto my_script.jl

to use multi-threading (see julia documentation for more details on multi-threading).

You can also use TMS in the interactive julia interpreter.

## Bug reports and features requests

TMS is still in development and certainly contains bugs. If you think you have found one please report it on the
[Github page](https://github.com/jerhoud/TensorMixedStates.jl).

Features requests may also be sent on the [Github page](https://github.com/jerhoud/TensorMixedStates.jl).

## Documentation

Your are currently reading it!

You can have access to inline documentation on TMS at the julia prompt simply by typing "?" followed by the function name or type name you are interested in. For example

    ?runTMS

For this to work you must have first imported TMS with

    using TensorMixedStates

## Examples

Working examples are presented in the folder `examples` in the repository. 
