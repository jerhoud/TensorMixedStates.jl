# The thermal states of an XXZ chain over a range of temperatures, by a Thermalize phase cooling
# the infinite temperature state step by step: each step gives its inverse temperature :beta and
# the logarithm of the partition function, :log_trace, up to the constant n log 2, written with
# the energy as JSON Lines, one object per temperature. The last one is checked against the
# exact thermal state, the chain being small enough to be written as a dense matrix.

using TensorMixedStates, .Qubits, LinearAlgebra

n = 6
Δ = 2.0
xxz = sum(X(i) * X(i + 1) + Y(i) * Y(i + 1) + Δ * Z(i) * Z(i + 1) for i in 1:n-1)
energy = "energy" => xxz

sim = runTMS(SimData(
    name = "XXZ_thermal",
    description = "thermal states of the XXZ chain, Δ = $Δ, from β = 0 to 2",
    phases = [
        CreateState{Mixed}(n, Qubit(), "FullyMixed"),
        Thermalize(hamiltonian = xxz, beta = 2.0, beta_step = 0.05, limits = Limits(maxdim = 64),
                   measurements = "thermal.jsonl" => [:beta, :log_trace, energy],
                   measurements_period = 4),
    ],
))

x, y, z = (matrix(op, Qubit()) for op in (X, Y, Z))
on(m, i) = kron([ j == i ? m : Matrix(1.0I, 2, 2) for j in 1:n ]...)
dense = sum(on(x, i) * on(x, i + 1) + on(y, i) * on(y, i + 1) + Δ * on(z, i) * on(z, i + 1)
            for i in 1:n-1)
λ = eigvals(Hermitian(dense))
w = exp.(-2.0 * λ)
exact_energy = sum(λ .* w) / sum(w)
println("energy at β = 2: ", expect(sim.state, xxz), "   exact ", exact_energy)
println("log Z at β = 2: ", sim.log_trace + n * log(2), "   exact ", log(sum(w)))
@assert abs(expect(sim.state, xxz) - exact_energy) < 1e-4
@assert abs(sim.log_trace + n * log(2) - log(sum(w))) < 1e-4
