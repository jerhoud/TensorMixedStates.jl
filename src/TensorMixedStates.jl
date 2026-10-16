# The module TensorMixedStates: what it imports from its dependencies, and its source files in
# the order they are read, a directory for each layer: core/ for sites, systems and states,
# engine/ for the operators on a system, measurements/ for what is measured and written,
# solvers/ for the evolutions and the searches, runner/ for the simulations runTMS runs, the
# site types in sites/, the hints Julia prints after an error, then the precompilation
# workload. Version 2 is written from its lowest layer up, the code of 1.x staying on the
# branch dev.

"""
    TensorMixedStates

a library for simulating closed and open quantum systems with matrix product states. A state
is either pure, its wave function as an MPS, or mixed, its density matrix as an MPS, on which
Lindbladian evolution and noisy gates act.

    using TensorMixedStates, .Qubits
    sys = System(10, Qubit())
    ham = -sum(Z(i) * Z(i + 1) for i in 1:9) - 0.5 * sum(X(i) for i in 1:10)
    E, ψ = ground_state(ham, random_state(Pure, sys, 8); algo = Dmrg(nsweeps = 10))
    ρ = evolve(-im * ham + sum(Dissipator(sqrt(0.1) * Sm)(i) for i in 1:10), 1.0, mix(ψ);
               nsteps = 20, limits = Limits(maxdim = 64))
    measure(ρ, [Z, X(1) * X(10), Purity])

The evolver is `-im * ham` plus dissipators, `"Up"` is ``|0\rangle`` on a qubit, and a value
is real for a Hermitian operator. The page "Quick reference" of the documentation, shipped with
the package as `joinpath(pkgdir(TensorMixedStates), "docs", "src", "quickref.md")`, gives the
whole library on one page, and "Upgrading from 1.x" what changed since version 1.
"""
module TensorMixedStates

import Base: *, +, -, /, ^, exp, sqrt, mod, show, showerror, ==, hash, length, getindex,
             getproperty, setproperty!, propertynames, truncate, open, close
import LinearAlgebra: norm, normalize, dot
import ITensors: dag, dim, matrix, flux, inner, siteind, siteinds
import ITensorMPS: maxlinkdim, expect, apply, sample

using ITensors, ITensorMPS, LinearAlgebra, Random, HDF5
using Serialization: serialize, deserialize
import Printf, JSON, Dates, Logging

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
include("engine/Ordering.jl")
include("engine/Lowering.jl")
include("engine/Compact.jl")
include("engine/Algebra.jl")
include("engine/Mpo.jl")
include("engine/Expect.jl")
include("engine/Apply.jl")
include("engine/PartialTrace.jl")
include("engine/GateStates.jl")

# Measurements and outputs
include("measurements/Measurements.jl")
include("measurements/Sampling.jl")
include("measurements/Measure.jl")
include("measurements/Outputs.jl")

# Solvers
include("solvers/Algorithms.jl")
include("solvers/ApproxW.jl")
include("solvers/Evolve.jl")
include("solvers/Search.jl")
include("solvers/Thermalize.jl")

# Simulations
include("runner/Threading.jl")
include("runner/Checkpoint.jl")
include("runner/Simulation.jl")
include("runner/Phases.jl")
include("runner/Run.jl")
include("runner/phases/CreateState.jl")
include("runner/phases/LoadState.jl")
include("runner/phases/Update.jl")
include("runner/phases/Transforms.jl")
include("runner/phases/Solvers.jl")

# Site types
include("sites/Qubits.jl")
include("sites/Fermions.jl")
include("sites/Bosons.jl")
include("sites/Spins.jl")
include("sites/Electrons.jl")
include("sites/Tjs.jl")
include("sites/Qbosons.jl")
include("sites/Qudits.jl")

# The hints of the errors of a program written for 1.x
include("Hints.jl")

# What precompilation runs, so that a first script does not wait for it to compile
include("Precompile.jl")

end
