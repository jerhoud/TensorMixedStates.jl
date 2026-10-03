# Reference values for the "Free bosons with source" testset of `algorithms.jl`, from the
# Gaussian moments of its Lindbladian, which is quadratic: a hopping hamiltonian
# H = Σ h_ab a_a† a_b on 4 sites and a gain L = √g a_2†, g = (2√0.1)² = 0.4.
#
# This file is not part of the test suite and nothing runs it automatically. It is here so
# that the numbers written in `algorithms.jl` can be checked and regenerated rather than
# taken on faith: run it with `julia test/reference/boson_gaussian.jl`. It uses LinearAlgebra
# only, on purpose — a reference produced by the code under test would check nothing.
#
# The correlations G_ij = <a_i† a_j> follow a closed linear equation, from
# d<O>/dt = i<[H, O]> + g <a_2 O a_2† - ½{a_2 a_2†, O}>:
#
#     dG/dt = i(hG - Gh) + g/2 (PG + GP) + g P,    P = e_2 e_2ᵀ,
#
# solved from the vacuum, G = 0, by the exponential of the matrix of the equation with its
# constant term appended. The state stays Gaussian, of zero mean and with no <a_i a_j>, so
# its purity is 1 / det(I + 2G). The test cuts each site at 6 bosons, which with at most 0.12
# boson per site changes nothing at the tolerance of the test.
using LinearAlgebra

n = 4
g = (2 * sqrt(0.1))^2
k = 2
h = zeros(n, n)
for i in 1:n-1
    h[i, i+1] = h[i+1, i] = 1
end
P = zeros(n, n)
P[k, k] = 1
Idn = Matrix{Float64}(I, n, n)
# vectorized column major: vec(AGB) = (Bᵀ⊗A) vec(G)
M = im * (kron(Idn, h) - kron(transpose(h), Idn)) + g / 2 * (kron(Idn, P) + kron(transpose(P), Idn))
augmented = [M vec(g * P); zeros(1, n * n) 0]
G = reshape((exp(augmented * 0.3) * [zeros(n * n); 1])[1:n*n], n, n)

println("N        : ", [round(real(G[i, i]); digits = 13) for i in 1:n])
println("A2†Ai    : ", [round(G[2, i]; digits = 13) for i in 1:n])
println("Purity   : ", round(real(1 / det(Idn + 2G)); digits = 13))
