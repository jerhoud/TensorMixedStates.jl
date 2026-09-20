# Reference values for the "Ising chain" testset of `algorithms.jl`, by exact
# diagonalization of the 6 qubit ring in its full 64 dimensional Hilbert space.
#
# This file is not part of the test suite and nothing runs it automatically. It is here so
# that the numbers written in `algorithms.jl` can be checked and regenerated rather than
# taken on faith: run it with `julia test/reference/ising_ed.jl`. It uses LinearAlgebra
# only, on purpose — a reference produced by the code under test would check nothing.
#
# The values it prints agree with what the evolution of the test produces to about 2e-8,
# which is the error of `ApproxW(order = 4, w = 2)` at `time_step = 0.02` with
# `maxdim = 8`, and the tolerance of the test is 1e-7.

using LinearAlgebra

const I2 = ComplexF64[1 0; 0 1]
const Xm = ComplexF64[0 1; 1 0]
const Ym = ComplexF64[0 -im; im 0]
const Zm = ComplexF64[1 0; 0 -1]

n = 6
at(o, i) = foldl(kron, [j == i ? o : I2 for j in 1:n])
prodat(o, sites) = foldl(*, [at(o, i) for i in sites])

# H = sum_i Z_i Z_{i+1} (periodic) - sum_i X_i, matching the evolver -im*H of the test
H = sum(at(Zm, i) * at(Zm, i + 1) for i in 1:n-1) + at(Zm, n) * at(Zm, 1) -
    sum(at(Xm, i) for i in 1:n)

plus = ComplexF64[1, 1] / sqrt(2)          # "X+"
psi0 = foldl(kron, [plus for _ in 1:n])
psi = exp(-im * Matrix(H) * 1.0) * psi0

ev(o) = real(psi' * o * psi)

println("X on each site : ", [round(ev(at(Xm, i)); digits = 11) for i in 1:n])
println("Y on each site : ", [round(ev(at(Ym, i)); digits = 11) for i in 1:n])
println("Z on each site : ", [round(ev(at(Zm, i)); digits = 11) for i in 1:n])
println("ZZ 12,23,16    : ", [round(ev(at(Zm,a)*at(Zm,b)); digits=11) for (a,b) in ((1,2),(2,3),(1,6))])
println("YY 12,23,16    : ", [round(ev(at(Ym,a)*at(Ym,b)); digits=11) for (a,b) in ((1,2),(2,3),(1,6))])
println("XX 12,23,16    : ", [round(ev(at(Xm,a)*at(Xm,b)); digits=11) for (a,b) in ((1,2),(2,3),(1,6))])
println("XXXX           : ", [round(ev(prodat(Xm, s)); digits=11) for s in ([1,2,3,4],[2,3,4,5],[4,5,6,1])])
println("XXXXXX         : ", round(ev(prodat(Xm, 1:6)); digits = 11))

m = reshape(psi, 8, 8)                      # cut between sites 3 and 4
p = svdvals(m) .^ 2
p = p[p .> 1e-30]
println("EE(3)          : ", round(-sum(p .* log.(p)); digits = 11))
