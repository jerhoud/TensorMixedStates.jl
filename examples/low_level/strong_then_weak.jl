# Fermions hopping on a chain under dephasing conserve their number strongly: every jump
# operator commutes with it, the density matrix stays in its sector, and its tensors are cut in
# the finest blocks. A loss changes the number, which a strong symmetry refuses: the state is
# weakened first, which keeps the number of particles a good quantum number of the density
# matrix without confining it to one sector. With a loss at the same rate on every site, the
# number of particles decays exactly as exp(-κt).

using TensorMixedStates, .Fermions

n = 6
γ = 0.5          # dephasing
κ = 0.3          # loss
hopping = im * sum(dag(C)(i) * C(i + 1) + dag(C)(i + 1) * C(i) for i in 1:n-1)   # -im * H
dephasing = sum(Dissipator(sqrt(γ) * N)(i) for i in 1:n)
loss = sum(Dissipator(sqrt(κ) * C)(i) for i in 1:n)
particles = sum(N(i) for i in 1:n)
# the whole bond dimension of a density matrix of 6 sites: with less, the number of particles
# below is off by the error of the truncation
lim = Limits(maxdim = 64)

ρ = State{Mixed}(System(n, Fermion(conserve = strong(N))), [isodd(i) ? "Occ" : "Emp" for i in 1:n])
ρ = evolve(hopping + dephasing, 1.0, ρ; nsteps = 20, limits = lim)
println("after the dephasing, ", real(expect(ρ, particles)), " particles, ",
        symmetries(ρ.system), " conserved")

# the loss is refused under the strong symmetry
try
    evolve(hopping + dephasing + loss, 0.1, ρ)
catch e
    println("refused: ", first(split(sprint(showerror, e), "\n")))
end

ρ = weaken(ρ)                            # one level down: N conserved weakly
t = 2.0
ρ = evolve(hopping + dephasing + loss, t, ρ; nsteps = 40, limits = lim)
println("after the loss, ", real(expect(ρ, particles)), " particles, exactly ", 3exp(-κ * t),
        ", ", symmetries(ρ.system), " conserved")
@assert abs(real(expect(ρ, particles)) - 3exp(-κ * t)) < 1e-8
