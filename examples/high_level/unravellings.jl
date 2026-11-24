# Two unravellings of one dephasing: an Ising chain in a transverse field, each qubit dephasing
# in its Z basis at the rate γ. dephasing_dissipator(γ) jumps to an eigenstate of Z, a
# measurement of the qubit; dephasing_dissipator(γ; jumps = :phases) applies Z at the rate γ / 2,
# which changes no population. Both give the same density matrix, so the same averages, but not
# the same trajectories: the spread of a correlation from one trajectory to another, std in the
# files, differs, and with it the number of trajectories an error bar needs.
#
# The averages of both are checked against the density matrix, State{Mixed}.

using TensorMixedStates, .Qubits

n = 6
γ = 0.1
ham = -sum(Z(i) * Z(i + 1) for i in 1:n-1) - sum(X(i) for i in 1:n)
correlation = Z(n ÷ 2) * Z(n ÷ 2 + 1)

dephasing(jumps) = sum(dephasing_dissipator(γ; jumps)(i) for i in 1:n)
# maxdim 4^3, the whole space of the density matrix, which the reference needs
evolution(jumps) = Evolve(evolver = -im * ham + dephasing(jumps),
                          duration = 2.0, time_step = 0.05, algo = ApproxW(order = 2),
                          limits = Limits(maxdim = 64),
                          final_measurements = Data("end") => ["ZZ" => correlation, X(n ÷ 2)])

simulate(type, jumps, samples) = runTMS(SimData(
    name = "unravelling_$(type)_$jumps",
    description = "an Ising chain dephasing, as a $type state, jumps $jumps",
    seed = 1,
    samples = samples,
    phases = [CreateState{type}(System(n, Qubit()), "Up"), evolution(jumps)],
)).data["end"] |> only

mixed = simulate(Mixed, :projectors, nothing)
spread = Dict{Symbol, Float64}()
for jumps in (:projectors, :phases)
    row = simulate(Trajectory, jumps, 200)
    spread[jumps] = row["std(ZZ)"]
    println("jumps $jumps: ⟨Z($(n ÷ 2))Z($(n ÷ 2 + 1))⟩ = ", round(row["ZZ"]; digits = 3),
            " ± ", round(row["err(ZZ)"]; digits = 3), ", spread ", round(row["std(ZZ)"]; digits = 3),
            "; ⟨X($(n ÷ 2))⟩ = ", round(row["X($(n ÷ 2))"]; digits = 3), " ± ",
            round(row["err(X($(n ÷ 2)))"]; digits = 3))
    @assert abs(row["ZZ"] - mixed["ZZ"]) < 4 * row["err(ZZ)"] + 1e-3
    @assert abs(row["X($(n ÷ 2))"] - mixed["X($(n ÷ 2))"]) < 4 * row["err(X($(n ÷ 2)))"] + 1e-3
end
println("density matrix: ⟨Z($(n ÷ 2))Z($(n ÷ 2 + 1))⟩ = ", round(mixed["ZZ"]; digits = 3),
        "; ⟨X($(n ÷ 2))⟩ = ", round(mixed["X($(n ÷ 2))"]; digits = 3))
# the phases, which leave the populations as they are, spread the correlation less
@assert spread[:phases] < spread[:projectors]
