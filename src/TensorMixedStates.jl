"""
    A module to make numerical simulations of closed or open quantum systems using Matrix Product States 
"""
module TensorMixedStates

import Base: *, +, -, /, ^, exp, sqrt, mod, show, length, getindex, isless, ==, hash
import ITensors: matrix, truncate, dim, Index, dag, norm, sim, inner, flux
import LinearAlgebra: dot
import ITensorMPS: maxlinkdim, apply, state, expect, normalize, measure!, checkdone!, tdvp, dmrg, sample

using ITensors, ITensorMPS, Printf, Dates, JSON, Random, LinearAlgebra, HDF5

# Core
include("Operators.jl")
include("Sites.jl")
include("Systems.jl")
include("Mixer.jl")
include("States.jl")
include("Simplify.jl")

# Low Level interface
include("Gates.jl")
include("Mpo.jl")
include("Observables.jl")
include("Measure.jl")
include("RandomState.jl")
include("Io.jl")

# High level interface
include("Checkpoint.jl")
include("Simulation.jl")
include("Output.jl")
include("Observers.jl")
include("Solvers.jl")
include("PhaseTypes.jl")
include("Phases.jl")
include("Run.jl")

# Utilities
include("Graphs.jl")

# Sites
include("Qubits.jl")
include("Fermions.jl")
include("Bosons.jl")
include("Spins.jl")
include("Electrons.jl")
include("Tjs.jl")
include("Qbosons.jl")
include("Qudits.jl")

# Precompilation
include("Precompile.jl")

end
