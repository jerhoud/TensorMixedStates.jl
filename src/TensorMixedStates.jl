"""
    A module to make numerical simulations of closed or open quantum systems using Matrix Product States 
"""
module TensorMixedStates

import Base: *, +, -, /, ^, exp, sqrt, show, length, getindex, isless, ==
import ITensors: matrix, truncate, dim, Index, dag, norm, sim
import ITensorMPS: maxlinkdim, apply, state, expect, normalize, measure!, checkdone!, tdvp, dmrg

using MKL, ITensors, ITensorMPS, Printf, Dates, JSON, Random

# Core
include("Operators.jl")
include("Sites.jl")
include("Mixer.jl")
include("Systems.jl")
include("States.jl")
include("Simplify.jl")

# Low Level interface
include("Gates.jl")
include("Mpo.jl")
include("Observables.jl")
include("Solvers.jl")
include("Measure.jl")
include("RandomState.jl")

# High level interface
include("Simulation.jl")
include("Observers.jl")
include("Output.jl")
include("Run.jl")
include("PhaseTypes.jl")
include("Phases.jl")

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

# Precompilation
#include("Precompile.jl")

end
