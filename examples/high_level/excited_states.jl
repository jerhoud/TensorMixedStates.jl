# The ground state and the first excited state of the frustrated spin chain at its
# Majumdar-Ghosh point, J2 = J1/2, where the ground state of an open chain of an even number of
# sites is exactly a product of singlets on the pairs (1, 2), (3, 4), ...: the search is checked
# against it. The ground state found is kept under a name by CopyState, and the next search is
# given it as the state to stay orthogonal to.

using TensorMixedStates, .Spins

n = 12
s(i, j) = Sx(i) * Sx(j) + Sy(i) * Sy(j) + Sz(i) * Sz(j)
hamiltonian = sum(s(i, i + 1) for i in 1:n-1) + 0.5 * sum(s(i, i + 2) for i in 1:n-2)
singlets = dimer_state(System(n, Spin(1/2)), [ (i, i + 1) for i in 1:2:n-1 ], "1/2", "-1/2")
search = Dmrg(nsweeps = 20, tol = 1e-12)
limits = Limits(maxdim = 64)

sim = runTMS(SimData(
    name = "Majumdar_Ghosh_$n",
    description = "the ground state and the first excited state of the J1-J2 chain, J2 = J1/2",
    seed = 1,
    phases = [
        CreateState(type = Pure, system = System(n, Spin(1/2)), randomize = 8),
        GroundState(hamiltonian = hamiltonian, algo = search, limits = limits,
                    measurements = "ground" => [:step, :energy],
                    final_measurements = "found" => [:energy, Fidelity(singlets)]),
        CopyState(:state => :ground),
        CreateState(type = Pure, system = System(n, Spin(1/2)), randomize = 8),
        GroundState(hamiltonian = hamiltonian, algo = search, limits = limits,
                    orthogonal_to = [:ground],
                    measurements = "excited" => [:step, :energy, Fidelity(:ground)]),
    ],
))

e0 = expect(sim.ground, hamiltonian)
e1 = expect(sim.state, hamiltonian)
println("ground state energy ", e0, ", that of the singlets ", expect(singlets, hamiltonian))
println("fidelity of the ground state to the singlets ", fidelity(sim.ground, singlets))
println("gap ", e1 - e0, ", overlap of the two states ", fidelity(sim.state, sim.ground))
@assert fidelity(sim.ground, singlets) > 1 - 1e-8
@assert fidelity(sim.state, sim.ground) < 1e-8
