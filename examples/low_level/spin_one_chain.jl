# The ground state of the spin one Heisenberg chain and its first excitation, found twice: as
# the lowest state orthogonal to the ground state, and as the ground state of the sector of
# total Sz = 1, the excitation being a triplet. Both are checked against the exact
# diagonalization of the chain, small enough here to be written as a dense matrix.

using TensorMixedStates, .Spins, LinearAlgebra

n = 6
heisenberg = sum(Sz(i) * Sz(i + 1) + (Sp(i) * Sm(i + 1) + Sm(i) * Sp(i + 1)) / 2 for i in 1:n-1)
lim = Limits(maxdim = 64)

# the ground state, then the lowest state orthogonal to it
sys = System(n, Spin(1))
e0, ground = ground_state(heisenberg, random_state(Pure, sys, 8); algo = Dmrg(nsweeps = 10),
                          limits = lim)
e1, excited = ground_state(heisenberg, random_state(Pure, sys, 8); algo = Dmrg(nsweeps = 20),
                           limits = lim, orthogonal_to = [ground])

# the same excitation as the ground state of the sector Sz = 1, its spin being conserved
conserving = System(n, Spin(1, conserve = 2Sz))
up = ["1", "0", "1", "-1", "0", "0"]                       # a state of total Sz = 1
e1_sector, _ = ground_state(heisenberg, State{Pure}(conserving, up); algo = Dmrg(nsweeps = 10),
                            limits = lim)

# the exact spectrum, from the dense matrix of the chain
exact = eigvals(Hermitian(matrix(heisenberg, sys)))

println("ground state energy  ", e0, "   exact ", exact[1])
println("first excitation     ", e1 - e0, "   in the sector Sz = 1 ", e1_sector - e0,
        "   exact ", exact[2] - exact[1])
@assert abs(e0 - exact[1]) < 1e-8
@assert abs(e1 - exact[2]) < 1e-8
@assert abs(e1_sector - exact[2]) < 1e-8
