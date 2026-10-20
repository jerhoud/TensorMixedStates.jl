# Free fermions hopping on a chain, at half filling: the Fermi sea TMS builds directly, and the
# ground state a search finds in the sector of the number of particles, compared with each other
# and with the exact solution, whose single particle states are sine waves.

using TensorMixedStates, .Fermions, LinearAlgebra

n = 20
nf = n ÷ 2
hamiltonian = -sum(dag(C)(i) * C(i + 1) + dag(C)(i + 1) * C(i) for i in 1:n-1)
sys = System(n, Fermion(conserve = N))

# the Fermi sea, exact as an MPS, and the ground state searched from a state of nf fermions
sea = fermi_sea(sys, hamiltonian, nf)
e, ground = ground_state(hamiltonian, State{Pure}(sys, [isodd(i) ? "Occ" : "Emp" for i in 1:n]);
                         algo = Dmrg(nsweeps = 10), limits = Limits(maxdim = 100))

# the exact solution: the energies -2cos(kπ/(n+1)) of the sine waves, the nf lowest filled
modes = [ sqrt(2 / (n + 1)) * sin(k * π * i / (n + 1)) for i in 1:n, k in 1:nf ]
exact_energy = sum(-2cos(k * π / (n + 1)) for k in 1:nf)
exact_correlations = modes * modes'                        # ⟨c†_i c_j⟩

correlations = expect2(ground, (dag(C), C))                # the matrix of ⟨c†_i c_j⟩
println("energy of the search ", e, "   of the Fermi sea ", expect(sea, hamiltonian),
        "   exact ", exact_energy)
println("overlap of the two states ", abs(inner(sea, ground)))
println("largest error on ⟨c†_i c_j⟩ ", maximum(abs.(correlations - exact_correlations)))
@assert abs(e - exact_energy) < 1e-8
@assert abs(expect(sea, hamiltonian) - exact_energy) < 1e-10
@assert maximum(abs.(correlations - exact_correlations)) < 1e-6
