# The fermion chain with dephasing of the article, examples/article/1_Fermion_chain_with_dephasing.jl,
# made to conserve its number of particles and then continued with a loss.
#
# Dephasing commutes with the number of particles, so the chain may conserve it strongly: every
# jump operator commutes with it, the state stays in its sector, and the tensors are cut into
# the finest blocks. A loss moves that number, which a strong symmetry forbids, so the second
# phase weakens it first. The weak symmetry that is left still has the density matrix commute
# with the number, but lets the state spread over several sectors.
#
# The sizes are small, the two evolutions taking about half a minute together once compiled.
# The article used n = 40, maxdim = 200 and a duration of 4.

using TensorMixedStates, .Fermions

limits = Limits(
    cutoff = 1e-30,
    maxdim = 32,
)

n = 8
gamma = 0.75        # dephasing
kappa = 0.2         # loss, in the second phase only
time_step = 0.05

output(n) = [
    "density.dat" => N,
    # constant while the symmetry is strong, then decaying as exp(-kappa t)
    "particles.dat" => sum(N(i) for i in 1:n),
    "OSEE.dat" => EntanglementEntropy(n ÷ 2, 4),
    "purity.dat" => [Purity, Trace, MaxLinkdim],
]

# the hamiltonian of a tight binding chain is minus its hopping sum, so -im * H is +im times it
hopping(n) = im * sum(dag(C)(i)C(i+1) + dag(C)(i+1)C(i) for i in 1:n-1)
dephasing(n) = sum(Dissipator(sqrt(4 * gamma) * N)(i) for i in 1:n)
loss(n) = sum(Dissipator(sqrt(kappa) * C)(i) for i in 1:n)

evolve(evolver) = Evolve(
    algo = Tdvp(),
    limits = limits,
    duration = 0.25,
    time_step = time_step,
    evolver = evolver,
    measures = output(n),
    measures_period = 2,
)

sim_data(n) = SimData(
    name = "Fermion_chain_conserved_$n",
    description = """
        Spinless fermion chain with $n sites, open boundary conditions and dephasing
        gamma = $gamma, conserving its number of particles strongly, then weakly once a
        loss kappa = $kappa is added
        cutoff = $(limits.cutoff)
        maxdim = $(limits.maxdim)
        time_step = $time_step
    """,
    phases = [
        CreateState(
            name = "Initialization in the alternating state",
            type = Mixed(),
            system = System(n, Fermion(conserve = strong(N))),
            state = [i % 2 == 0 ? "Occ" : "Emp" for i in 1:n],
            final_measures = output(n),
        ),
        evolve(hopping(n) + dephasing(n)),
        Weaken(),
        evolve(hopping(n) + dephasing(n) + loss(n)),
    ],
)

runTMS(sim_data(n); restart = true)
