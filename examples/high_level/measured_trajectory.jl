# A quantum trajectory: a chain of qubits kicked, its first qubit measured after each kick, the
# state collapsing on the result drawn. An Update phase runs a function of the context at each
# of its steps, which reads and changes the state and provides the result as a value, :result,
# measured with the rest. Its draws use the generator of the run, seeded, so that running the
# program again gives the same trajectory, a resumed run included.

using TensorMixedStates, .Qubits

n = 6
# a kick: the coupling of each pair of neighbours, then a rotation of each qubit
kick = prod(exp(-0.4im * X)(i) for i in 1:n) * prod(exp(-0.3im * Z ⊗ Z)(i, i + 1) for i in 1:n-1)
limits = Limits(maxdim = 32)

function step!(c)
    c.result, _, c.state = collapse(apply(kick, c.state; limits), Z(1); rng = c.rng)
end

sim = runTMS(SimData(
    name = "trajectory",
    description = "a chain of qubits kicked, its first qubit measured after each kick",
    seed = 7,
    phases = [
        CreateState{Pure}(n, Qubit(), "Up"),
        Update(step!; name = "Kicked and measured", nsteps = 30, provides = [:result],
               measurements = [ "trajectory" => [:step, :result, Z(1), Z(2)],
                                Data("trajectory") => [:result, Z(1), EntanglementEntropy(3)] ]),
    ],
))

rows = sim.data["trajectory"]
results = [ Dict(row)["result"] for row in rows ]
println("results of the measurements: ", Int.(results))
# after a measurement, the first qubit is in the state of the result
@assert all(Dict(row)["result"] ≈ Dict(row)["Z(1)"] for row in rows)
