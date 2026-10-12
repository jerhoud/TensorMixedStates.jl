# The module TensorMixedStates: what it imports from its dependencies, and its source files in
# the order they are read, a directory for each layer: core/ for sites, systems and states,
# engine/ for the operators on a system, measurements/ for what is measured and written, then
# the site types in sites/. Version 2 is written from its lowest layer up, the code of 1.x
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
import ITensorMPS: maxlinkdim, expect, apply, sample

using ITensors, ITensorMPS, LinearAlgebra, Random, HDF5

# Sites, systems and states
include("core/Operators.jl")
include("core/Sites.jl")
include("core/Matrices.jl")
include("core/Noise.jl")
include("core/Conservation.jl")
include("core/Systems.jl")
include("core/States.jl")
include("core/PreparedStates.jl")
include("core/RandomStates.jl")
include("core/Io.jl")
include("core/Graphs.jl")

# Operators on a system
include("engine/Split.jl")
include("engine/Lowering.jl")
include("engine/Compact.jl")
include("engine/Mpo.jl")
include("engine/Expect.jl")
include("engine/Apply.jl")
include("engine/PartialTrace.jl")
include("engine/GateStates.jl")

# Measurements and outputs
include("measurements/Measurements.jl")
include("measurements/Sampling.jl")

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
