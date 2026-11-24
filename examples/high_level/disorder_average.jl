# A disorder average: free fermions on a chain with random on-site energies, a domain wall
# melting. Each sample draws its energies in Realizations, with the generator c.rng of the
# sample, and evolves the wall under its own Hamiltonian; the files and the Data hold the
# averages over the samples, with their errors, and EachSample a row for each sample. The
# samples stop once the errors on the density are below the tolerance.
#
# Each sample is checked against the exact solution of its free fermions, and the average
# against an average of the exact solution over independent draws of the disorder.

using TensorMixedStates, .Fermions, LinearAlgebra, Random

n = 10
W = 2.0                 # on-site energies uniform in [-W, W]
duration = 2.0
hopping = -sum(dag(C)(i) * C(i + 1) + dag(C)(i + 1) * C(i) for i in 1:n-1)
wall = [i <= n ÷ 2 ? "Occ" : "Emp" for i in 1:n]
energies(rng) = W * (2 * rand(rng, n) .- 1)

sim = runTMS(SimData(
    name = "disorder_average",
    description = "a domain wall of free fermions melting, averaged over random on-site energies",
    seed = 1,
    samples = 1000,                      # the most that run
    tol = "density" => 0.03,             # stop once every err(N(i)) is below 0.03
    phases = [
        CreateState{Pure}(System(n, Fermion(conserve = N)), wall),
        Realizations() do c
            ε = energies(c.rng)
            [Update(d -> d.energies = ε; provides = [:energies]),
             Evolve(evolver = -im * (hopping + sum(ε[i] * N(i) for i in 1:n)), duration = duration,
                    time_step = 0.1, algo = ApproxW(order = 4), limits = Limits(maxdim = 32))]
        end,
    ],
    final_measurements = ["density" => N, Data("density") => N,
                          EachSample(Data("samples")) => [:energies, N]],
))

# the exact density of free fermions: ⟨n_i(t)⟩ = Σ_j |U_ij|² ⟨n_j(0)⟩ with U = exp(-iht)
function exact_density(ε)
    h = SymTridiagonal(ε, -ones(n - 1))
    U = exp(-im * duration * Matrix(h))
    return abs2.(U) * [i <= n ÷ 2 ? 1.0 : 0.0 for i in 1:n]
end

samples = sim.data["samples"]
# a row of EachSample, as a Data without samples, holds the vector of N; an averaged one, N(i)
worst = maximum(maximum(abs.(row["N"] - exact_density(row["energies"]))) for row in samples)
println(length(samples), " samples, largest error of a sample against its exact density ", worst)
@assert worst < 1e-4

# the average of the exact density over independent draws, against that of the samples
draws = 20000
rng = Xoshiro(0)
exact = [exact_density(energies(rng)) for _ in 1:draws]
mean_exact = sum(exact) / draws
err_exact = sqrt.(sum(x -> (x - mean_exact) .^ 2, exact) / (draws - 1) / draws)
average = only(sim.data["density"])
for i in 1:n
    m, e = average["N($i)"], average["err(N($i))"]
    println("site $i: ", round(m; digits = 4), " ± ", round(e; digits = 4),
            "   exact average ", round(mean_exact[i]; digits = 4))
    @assert abs(m - mean_exact[i]) < 4 * sqrt(e^2 + err_exact[i]^2)
end
