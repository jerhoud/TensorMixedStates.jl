# The thermal state of a transverse field Ising chain, by an evolution in imaginary time from the
# infinite temperature state. thermalize also returns the logarithm of the partition function, up
# to the constant n log 2, whose derivative in β is minus the energy: two ways of computing the
# energy, which agree when the step and the truncation are good, both checked here against the
# exact thermal state, the chain being small enough to be written as a dense matrix.

using TensorMixedStates, .Qubits, LinearAlgebra

n = 8
β = 1.0
ising = -sum(Z(i) * Z(i + 1) for i in 1:n-1) - 0.7 * sum(X(i) for i in 1:n)
hot = State{Mixed}(System(n, Qubit()), "FullyMixed")
lim = Limits(maxdim = 64)

logtrace, thermal = thermalize(ising, β, hot; nsteps = 40, limits = lim)
energy = real(expect(thermal, ising))
# the derivative of log Z, by two other temperatures close by
δ = 0.01
plus, _ = thermalize(ising, β + δ, hot; nsteps = 40, limits = lim)
minus, _ = thermalize(ising, β - δ, hot; nsteps = 40, limits = lim)
energy_from_z = -(plus - minus) / (2δ)

# exactly
x, z = matrix(X, Qubit()), matrix(Z, Qubit())
on(m, i) = kron([ j == i ? m : Matrix(1.0I, 2, 2) for j in 1:n ]...)
dense = -sum(on(z, i) * on(z, i + 1) for i in 1:n-1) - 0.7 * sum(on(x, i) for i in 1:n)
λ = eigvals(Hermitian(dense))
exact_logz = log(sum(exp.(-β * λ)))
exact_energy = sum(λ .* exp.(-β * λ)) / sum(exp.(-β * λ))

println("log Z  ", logtrace + n * log(2), "   exact ", exact_logz)
println("energy ", energy, ", from log Z ", energy_from_z, "   exact ", exact_energy)
@assert abs(logtrace + n * log(2) - exact_logz) < 1e-6
@assert abs(energy - exact_energy) < 1e-6
@assert abs(energy_from_z - exact_energy) < 1e-3
