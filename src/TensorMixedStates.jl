# The module TensorMixedStates: what it imports from its dependencies, and its source files in
# the order they are read. Version 2 is written from its lowest layer up, the code of 1.x
# staying on the branch dev.

"""
    TensorMixedStates

a library for simulating closed and open quantum systems with matrix product states. A state
is either pure, its wave function as an MPS, or mixed, its density matrix as an MPS, on which
Lindbladian evolution and noisy gates act.
"""
module TensorMixedStates

import Base: *, +, -, /, ^, exp, sqrt, mod, show, isless, ==, hash, length, getindex, truncate
import LinearAlgebra: norm, normalize, dot
import ITensors: dag, dim, matrix, flux, inner, siteind, siteinds
import ITensorMPS: maxlinkdim

using ITensors, ITensorMPS, LinearAlgebra, Random, HDF5

# Sites, systems and states
include("Operators.jl")
include("Sites.jl")
include("Matrices.jl")
include("Conservation.jl")
include("Systems.jl")
include("States.jl")
include("PreparedStates.jl")
include("RandomStates.jl")
include("Io.jl")

# Operators on a system
include("Split.jl")
include("Lowering.jl")
include("Compact.jl")
include("Mpo.jl")

# Site types
include("sites/Qubits.jl")
include("sites/Fermions.jl")
include("sites/Bosons.jl")
include("sites/Spins.jl")
include("sites/Electrons.jl")
include("sites/Tjs.jl")
include("sites/Qbosons.jl")
include("sites/Qudits.jl")

end
