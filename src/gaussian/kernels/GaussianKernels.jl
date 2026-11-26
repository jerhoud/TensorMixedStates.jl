# The numerical kernels of the Gaussian backend, on plain matrices and using no type of the
# library, so that they can become a package of their own: pure fermionic Gaussian states with
# their phase, their overlaps, evolutions and measurements. Conventions of Bravyi and Gosset
# (Commun. Math. Phys. 356, 451 (2017), arXiv:1609.00735): γ_{2j-1} = c_j + c_j†,
# γ_{2j} = i(c_j† - c_j), covariance M_ab = (-i/2)⟨[γ_a, γ_b]⟩, so that the vacuum has
# M_{2j-1,2j} = 1; the Jordan-Wigner string of mode j runs over the modes before it.
module GaussianKernels

using LinearAlgebra
using SkewLinearAlgebra: logabspfaffian

export LogComplex, logpfaffian, fock_covariance, majorana_on_fock, PureGaussian, nmodes,
       find_support, describe, apply_monomial, overlap, cross_term, expect_monomial,
       transition_matrix

include("Numbers.jl")
include("Fock.jl")
include("States.jl")
include("Overlaps.jl")

end
