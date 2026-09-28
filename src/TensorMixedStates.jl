"""
    A module to make numerical simulations of closed or open quantum systems using Matrix Product States 
"""
module TensorMixedStates

import Base: *, +, -, /, ^, exp, sqrt, mod, show, length, getindex, isless, ==, hash
import ITensors: matrix, truncate, dim, Index, dag, norm, sim, inner, flux, isfermionic
import LinearAlgebra: dot
import ITensorMPS: maxlinkdim, apply, state, expect, normalize, checkdone!, tdvp, dmrg, sample

using ITensors, ITensorMPS, Printf, Dates, JSON, Random, LinearAlgebra, HDF5
import Logging

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
include("Destinations.jl")
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
include("sites/Qubits.jl")
include("sites/Fermions.jl")
include("sites/Bosons.jl")
include("sites/Spins.jl")
include("sites/Electrons.jl")
include("sites/Tjs.jl")
include("sites/Qbosons.jl")
include("sites/Qudits.jl")

# Precompilation
include("Precompile.jl")

end
