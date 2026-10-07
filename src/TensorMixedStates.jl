# The module TensorMixedStates: what it imports from its dependencies, and its source files in
# the order they are read.

"""
    TensorMixedStates

a library for simulating closed and open quantum systems with matrix product states. A state
is either pure, its wave function as an MPS, or mixed, its density matrix as an MPS, on which
Lindbladian evolution and noisy gates act.
"""
module TensorMixedStates

import Base: *, +, -, /, ^, exp, sqrt, mod, show, length, getindex, isless, ==, hash, isapprox
import ITensors: matrix, truncate, dim, Index, dag, norm, sim, inner, flux, isfermionic
import LinearAlgebra: dot
import ITensorMPS: maxlinkdim, apply, state, expect, normalize, checkdone!, tdvp, dmrg, sample

using ITensors, ITensorMPS, Printf, Dates, JSON, Random, LinearAlgebra, HDF5
import Logging

# The names of the experimental interfaces of one's own, phases, algorithms and
# representations, are not exported, see docs/src/extending.md

# Core
include("Operators.jl")
include("Sites.jl")
include("Systems.jl")
include("Mixer.jl")
include("Definitions.jl")
include("States.jl")
include("Simplify.jl")
include("Compact.jl")

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
include("Threading.jl")
include("Output.jl")
include("Observers.jl")
include("Solvers.jl")
include("Phases.jl")
include("Run.jl")

# The phases, each written through the interface of Phases.jl
include("phases/CreateState.jl")
include("phases/LoadState.jl")
include("phases/SaveState.jl")
include("phases/ToMixed.jl")
include("phases/Evolve.jl")
include("phases/Gates.jl")
include("phases/search.jl")
include("phases/GroundState.jl")
include("phases/SteadyState.jl")
include("phases/Thermalize.jl")
include("phases/PartialTrace.jl")
include("phases/Weaken.jl")

# Utilities
include("Graphs.jl")
include("FreeFermions.jl")
include("PreparedStates.jl")

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
