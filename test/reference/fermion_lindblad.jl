# Reference values for the "Free fermions with source" testset of `algorithms.jl`, by dense
# Lindblad evolution of the 5 site chain: Jordan-Wigner by hand, then the exponential of the
# vectorized Liouvillian on the 32 dimensional Fock space.
#
# This file is not part of the test suite and nothing runs it automatically. It is here so
# that the numbers written in `algorithms.jl` can be checked and regenerated rather than
# taken on faith: run it with `julia test/reference/fermion_lindblad.jl`. It uses
# LinearAlgebra only, on purpose — a reference produced by the code under test would check
# nothing.
#
# The values it prints agree with what the evolution of the test produces to about 8e-8,
# which is the error of `Tdvp` at `time_step = 0.05` with `maxdim = 16`, and the tolerance
# of the test is 1e-6. `trace` is printed as a check on the Liouvillian itself: a Lindblad
# evolution preserves it, so anything but 1 means the reference is wrong, not the library.
using LinearAlgebra

const I2 = ComplexF64[1 0; 0 1]
const Cm = ComplexF64[0 1; 0 0]     # annihilation, same basis order as Fermions.jl
const Fm = ComplexF64[1 0; 0 -1]    # Jordan-Wigner string

n = 5
c(i) = foldl(kron, [j < i ? Fm : j == i ? Cm : I2 for j in 1:n])

H = sum(c(i)' * c(i + 1) + c(i + 1)' * c(i) for i in 1:n-1)
L = sqrt(2 * 0.2) * c(3)'           # the Dissipator of the test

d = 2^n
Id = Matrix{ComplexF64}(I, d, d)
# dρ/dt = -i[H,ρ] + LρL† - ½{L†L, ρ}, vectorized column major: vec(AρB) = (Bᵀ⊗A) vec(ρ)
liou = -im * (kron(Id, H) - kron(transpose(H), Id)) +
       kron(transpose(L'), L) -
       0.5 * (kron(Id, L' * L) + kron(transpose(L' * L), Id))

rho0 = zeros(ComplexF64, d, d)
rho0[1, 1] = 1.0                    # every site in state "0"
rho = reshape(exp(liou * 1.0) * vec(rho0), d, d)

println("trace    : ", round(real(tr(rho)); digits = 12))
println("N        : ", [round(real(tr(rho * c(i)' * c(i))); digits = 13) for i in 1:n])
println("C3†Ci    : ", [round(tr(rho * c(3)' * c(i)); digits = 13) for i in 1:n])
println("Purity   : ", round(real(tr(rho * rho)); digits = 13))
