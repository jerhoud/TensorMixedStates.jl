# Quantum trajectories: a chain of qubits in a transverse field, pumped from down to up by a
# dissipator, unravelled into trajectories. CreateState{Trajectory} is the branch point: each
# sample is one trajectory, a pure state evolving under the effective Hamiltonian between jumps
# drawn at random, and the files and the Data hold the averages over the trajectories, with
# their errors. The value :jumps counts the jumps of each trajectory.
#
# The same evolution of the density matrix, State{Mixed}, gives the reference: the average of
# the trajectories is checked against it, and the mean number of jumps against the integral of
# the rate of the jumps, γ Σ_i ⟨Sp† Sp⟩(i), over the evolution.

using TensorMixedStates, .Qubits

n = 4
γ = 0.2
ham = -sum(Z(i) * Z(i + 1) for i in 1:n-1) - 0.8 * sum(X(i) for i in 1:n)
pump = sum(Dissipator(sqrt(γ) * Sp)(i) for i in 1:n)         # "Dn" to "Up" at the rate γ
rate = γ * sum((dag(Sp) * Sp)(i) for i in 1:n)
time_step = 0.02

measured = [Z, "rate" => rate]
# ApproxW rather than the default Tdvp: on 4 qubits, a step of Tdvp costs about seven times more
evolution(period, measured) = Evolve(evolver = -im * ham + pump, duration = 2.0,
                                     time_step = time_step, algo = ApproxW(order = 2),
                                     limits = Limits(maxdim = 16),
                                     measurements_period = period,
                                     measurements = Data("z") => measured)

mixed = runTMS(SimData(
    name = "damped_mixed",
    description = "the density matrix of a pumped chain of qubits",
    phases = [CreateState{Mixed}(System(n, Qubit()), "Dn";
                                 final_measurements = Data("z") => measured),
              evolution(1, measured)],
))

trajectories = runTMS(SimData(
    name = "damped_trajectories",
    description = "the trajectories of a pumped chain of qubits",
    seed = 1,
    samples = 200,
    phases = [CreateState{Trajectory}(System(n, Qubit()), "Dn"), evolution(10, [measured; :jumps])],
    final_measurements = EachSample(Data("jumps")) => :jumps,
))

# the averages of the trajectories, every 10 steps, against the density matrix at the same time,
# within 4 errors and the bias of the time step
reference = Dict(round(row["time"]; digits = 6) => row for row in mixed.data["z"])
for row in trajectories.data["z"], i in 1:n
    exact = reference[round(row["time"]; digits = 6)]["Z"][i]     # the vector, without samples
    @assert abs(row["Z($i)"] - exact) < 4 * row["err(Z($i))"] + 1e-3
end
final = trajectories.data["z"][end]
println("at t = ", final["time"], ", ⟨Z(1)⟩ = ", round(final["Z(1)"]; digits = 3), " ± ",
        round(final["err(Z(1))"]; digits = 3), " over ", final["N"], " trajectories, ",
        round(reference[round(final["time"]; digits = 6)]["Z"][1]; digits = 3),
        " for the density matrix")

# the mean number of jumps against the integral of the rate, by the trapezoidal rule
times = [row["time"] for row in mixed.data["z"]]
rates = [row["rate"] for row in mixed.data["z"]]
expected = sum((times[k+1] - times[k]) * (rates[k] + rates[k+1]) / 2 for k in 1:length(times)-1)
println("jumps per trajectory ", round(final["jumps"]; digits = 3), " ± ",
        round(final["err(jumps)"]; digits = 3), ", integral of the rate ",
        round(expected; digits = 3))
@assert abs(final["jumps"] - expected) < 4 * final["err(jumps)"] + 1e-3
counts = [Int(row["jumps"]) for row in trajectories.data["jumps"]]
println("trajectories with 0, 1, 2, … jumps: ", [count(==(k), counts) for k in 0:maximum(counts)])
