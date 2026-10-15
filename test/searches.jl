# The searches of a state: ground_state, the lowest eigenvector of a Hamiltonian, and
# steady_state, that of L†L for a Lindbladian L, by Dmrg, against exact diagonalization and
# steady states known exactly; the tolerance that stops the sweeps, a search resumed after some
# of its sweeps, the noise and the limits of each sweep, the parameters of the Krylov method and
# their defaults field by field, the MPOs a PreMPO keeps, and what the searches refuse. The
# searches of one's own are in extension.jl, the phases GroundState and SteadyState resumed in
# runner.jl.

@testset "Dmrg" begin
    # from a random state, the ground state of a field, every spin up
    e, st = ground_state(-sum(Z(i) for i in 1:5), random_state(Pure, System(5, Qubit()), 10);
                         algo = Dmrg(nsweeps = 2), limits = Limits(maxdim = 10))
    @test e ≈ -5 atol = 1e-7
    @test real(expect1(st, X)) ≈ zeros(5) atol = 1e-7
    @test real(expect1(st, Y)) ≈ zeros(5) atol = 1e-7
    @test real(expect1(st, Z)) ≈ ones(5) atol = 1e-7
    @test norm(st) ≈ 1
end

@testset "Ground state against exact diagonalization" begin
    # a transverse field Ising chain, whose ground state is the lowest eigenvector of its dense
    # Hamiltonian, with a variance the energy alone cannot show
    LA = TMS.LinearAlgebra
    n = 6
    sys = System(n, Qubit())
    x, z = [0. 1. ; 1. 0.], [1. 0. ; 0. -1.]
    on(m, i) = kron([ k == i ? m : [1. 0. ; 0. 1.] for k in 1:n ]...)
    h = -sum(Z(i) * Z(i + 1) for i in 1:n-1) - sum(X(i) for i in 1:n)
    dense = -sum(on(z, i) * on(z, i + 1) for i in 1:n-1) - sum(on(x, i) for i in 1:n)
    f = LA.eigen(LA.Symmetric(dense))
    r = random_state(Pure, sys, 16)
    e, gs = ground_state(h, r; algo = Dmrg(nsweeps = 12), limits = Limits(cutoff = 1e-14, maxdim = 32))
    @test e ≈ f.values[1] atol = 1e-10
    @test e ≈ real(expect(gs, h))
    @test abs(f.vectors[:, 1]' * dense_vector(gs)) ≈ 1 atol = 1e-10
    @test norm(gs) ≈ 1
    @test variance(gs, h) < 1e-8
    # and the limits of each sweep, the bond dimension growing to its last value
    e, gs = ground_state(h, r; algo = Dmrg(nsweeps = 8),
                         limits = Limits(cutoff = 1e-14, maxdim = [2, 4, 8]))
    @test maxlinkdim(gs) ≤ 8
    @test e ≈ f.values[1] atol = 1e-8
    # on fermions, the two modes of negative energy of a chain of four filled, -2cos(π/5) -
    # 2cos(2π/5), from a state of two fermions, the number being conserved by the hopping
    h = -sum(dag(C)(i) * C(i+1) + dag(C)(i+1) * C(i) for i in 1:3)
    st = State{Pure}(System(4, Fermion()), ["Occ", "Emp", "Occ", "Emp"])
    @test first(ground_state(h, st; algo = Dmrg(nsweeps = 3), limits = Limits(maxdim = 16))) ≈ -sqrt(5)
end

@testset "The tolerance stops the sweeps" begin
    # the search stops once the energy changes by less than `tol` from one sweep to the next,
    # which the end of each sweep, where the phases will hook, shows
    n = 6
    h = -sum(Z(i) * Z(i + 1) for i in 1:n-1) - sum(X(i) for i in 1:n)
    r = random_state(Pure, System(n, Qubit()), 16)
    mpo = make_mpo(r, h)
    lim = Limits(cutoff = 1e-14, maxdim = 32)
    function sweeps(algo; kwargs...)
        seen = []
        e, st = TMS.search(algo, mpo, r, lim;
                           sweep_done = (k, e, st, stop) -> push!(seen, (k, e, st, stop)),
                           kwargs...)
        return e, st, seen
    end
    e, st, seen = sweeps(Dmrg(nsweeps = 20, tol = 1e-10))
    k = length(seen)
    @test 3 ≤ k < 20
    @test first.(seen) == 1:k
    @test abs(seen[k][2] - seen[k-1][2]) < 1e-10
    @test abs(seen[k-1][2] - seen[k-2][2]) ≥ 1e-10
    @test seen[k][2] ≈ e
    # the last sweep says that the search ends there, which a phase commits as its end
    @test last.(seen) == [ falses(k - 1) ; true ]
    # the state of each sweep is normalized, as the one returned
    @test norm(seen[k][3]) ≈ 1
    @test norm(seen[k][3] - st) < 1e-12
    # with no tolerance, every sweep is made, the last one ending the search
    @test first.(last(sweeps(Dmrg(nsweeps = 5)))) == 1:5
    @test last.(last(sweeps(Dmrg(nsweeps = 5)))) == [ falses(4) ; true ]
    @test first(ground_state(h, r; algo = Dmrg(nsweeps = 20, tol = 1e-10), limits = lim)) ≈ e
end

@testset "A search resumed after some of its sweeps" begin
    # a search resumed after its last sweep is asked for none, on which ITensorMPS gave an
    # energy of 0: it is that of its last sweep, as the search it resumes gave it, or else that
    # of the state given
    h = -sum(Z(i) for i in 1:2) - 0.5 * X(1)
    st = State{Pure}(System(2, Qubit()), "Up")
    e, s = TMS.search(Dmrg(nsweeps = 2), make_mpo(st, h), st, Limits(); done = 2)
    @test e ≈ real(expect(st, h))
    @test norm(s - st) < 1e-14
    @test first(TMS.search(Dmrg(nsweeps = 2), make_mpo(st, h), st, Limits(); done = 2,
                           energy = -1.5)) == -1.5
    # the sweeps are counted after those made, and take the rest of the schedules: the bond
    # dimension of the last sweep, 4, where a search of one sweep keeps the first, 1
    n = 4
    hx = sum(X(i) * X(i + 1) for i in 1:n - 1)
    up = State{Pure}(System(n, Qubit()), "Up")
    mpo = make_mpo(up, hx)
    lim = Limits(maxdim = [1, 1, 4])
    for done in 0:2
        seen = Int[]
        _, s = TMS.search(Dmrg(nsweeps = 3), mpo, up, lim; done,
                          sweep_done = (k, _, _, _) -> push!(seen, k))
        @test seen == done+1:3
        @test maxlinkdim(s) > 1
    end
    @test maxlinkdim(last(TMS.search(Dmrg(nsweeps = 1), mpo, up, lim))) == 1
    later = TMS.later_steps
    @test later([1, 2, 3], 1) == [2, 3]
    @test later([1, 2, 3], 5) == [3]
    @test later(0.5, 3) == 0.5
    @test later(Limits(cutoff = 1e-10, maxdim = [1, 2, 3]), 2) == Limits(cutoff = 1e-10, maxdim = [3])
    # the energy of the last sweep made is compared with that of the first one resumed: from
    # the ground state, one sweep is enough to stop, where two are needed without it
    e0, gs = ground_state(hx, up; algo = Dmrg(nsweeps = 4), limits = Limits(maxdim = 8))
    for (energy, count) in ((e0, 1), (nothing, 2))
        seen = Int[]
        TMS.search(Dmrg(nsweeps = 10, tol = 1e-8), mpo, gs, Limits(maxdim = 8); done = 3, energy,
                   sweep_done = (k, _, _, _) -> push!(seen, k))
        @test seen == 4:3+count
    end
end

@testset "Steady state" begin
    # driven decaying qubits, H = Ω/2 X and a decay of rate γ down to Dn by Sm: the optical
    # Bloch equations give ⟨Z⟩ = -γ²/(γ² + 2Ω²), ⟨Y⟩ = 2γΩ/(γ² + 2Ω²) and ⟨X⟩ = 0, on two
    # qubits, dmrg needing two sites at least
    Ω, γ = 0.8, 0.5
    L = sum(-im * Ω / 2 * X(i) + Dissipator(sqrt(γ) * Sm)(i) for i in 1:2)
    v, ρ = steady_state(L, State{Mixed}(System(2, Qubit()), "FullyMixed"); algo = Dmrg(nsweeps = 10))
    @test abs(v) < 1e-12
    @test real(expect1(ρ, Z)) ≈ fill(-γ^2 / (γ^2 + 2Ω^2), 2) atol = 1e-10
    @test real(expect1(ρ, Y)) ≈ fill(2γ * Ω / (γ^2 + 2Ω^2), 2) atol = 1e-10
    @test real(expect1(ρ, X)) ≈ zeros(2) atol = 1e-10
    @test trace(ρ) ≈ 1
    @test hermiticity(ρ) ≈ 1

    # pumped and damped
    v, ρ = steady_state(Dissipator(Sp)(1) + Dissipator(Sm)(2), State{Mixed}(System(2, Qubit()), "+");
                        algo = Dmrg(nsweeps = 20), limits = Limits(cutoff = 1e-10, maxdim = 10))
    @test real(expect1(ρ, X)) ≈ zeros(2) atol = 1e-10
    @test real(expect1(ρ, Y)) ≈ zeros(2) atol = 1e-10
    @test real(expect1(ρ, Z)) ≈ [1, -1] atol = 1e-10
    # the state found is a density matrix, of trace one, where the eigenvector of L†L has norm
    # one and a sign of its own
    @test trace(ρ) ≈ 1

    # from a random state, a chain coupled along Z, every spin pumped up. The tolerances of these
    # searches hid in 1.x an error of the default Krylov search of ITensorMPS, 1e-5 here, which
    # the default of steady_state resolves, see "The Krylov parameters reach the searches"
    _, ρ = steady_state(-im * (-sum(Z(i) * Z(i+1) for i in 1:4)) + sum(Dissipator(Sp)(i) for i in 1:5),
                        random_state(Mixed, System(5, Qubit()), 10); algo = Dmrg(nsweeps = 40),
                        limits = Limits(maxdim = 10, cutoff = 1e-10))
    @test real(expect1(ρ, X)) ≈ zeros(5) atol = 1e-8
    @test real(expect1(ρ, Y)) ≈ zeros(5) atol = 1e-8
    @test real(expect1(ρ, Z)) ≈ ones(5) atol = 1e-8

    # an XX chain pumped at one end and damped at the other, which carries a current
    _, ρ = steady_state(-im * sum(X(i) * X(i+1) + Y(i) * Y(i+1) for i in 1:3) +
                        Dissipator(Sp)(1) + Dissipator(Sm)(4), State{Mixed}(System(4, Qubit()), "FullyMixed");
                        algo = Dmrg(nsweeps = 200), limits = Limits(cutoff = 1e-10, maxdim = 10))
    @test real(expect1(ρ, Z)) ≈ [1 / 17, 0., 0., -1 / 17] atol = 1e-7
    @test real(expect(ρ, [2 * (X(i) * Y(i+1) - Y(i) * X(i+1)) for i in 1:3])) ≈ fill(16 / 17, 3) atol = 1e-7
end

@testset "A steady state that is not unique" begin
    # dephasing keeps the number of fermions, so that the Lindbladian has a steady state in each
    # sector of it, of which dmrg returns the combination its starting state leads to, in
    # general not a density matrix: in 1.x, rounding took a state of two fermions to one that
    # was not hermitian, its trace complex and ⟨N⟩ far from 0.5, and nothing said so. A
    # starting state of complex weights on two sectors leads to such a combination, which the
    # warning shows. Conserved strongly, the number keeps the search in the sector of the
    # starting state, where the steady state is unique
    h = -sum(dag(C)(i) * C(i+1) + dag(C)(i+1) * C(i) for i in 1:3)
    lind = -im * h + sum(Dissipator(N)(i) for i in 1:4)
    lim = Limits(maxdim = 32)
    algo = Dmrg(nsweeps = 10)
    start(site) = State{Mixed}(System(4, site), ["Occ", "Emp", "Occ", "Emp"])
    one = State{Mixed}(System(4, Fermion()), ["Occ", "Emp", "Emp", "Emp"])
    # other warnings may come, those of KrylovKit in its oldest versions for instance
    @test_logs (:warn, r"not unique") match_mode = :any steady_state(lind, start(Fermion()) + im * one;
                                                                     algo, limits = lim)
    logs, (_, ρ) = Test.collect_test_logs() do
        steady_state(lind, start(Fermion(conserve = strong(N))); algo, limits = lim)
    end
    @test !any(l -> occursin("not unique", string(l.message)), logs)
    @test real(expect1(ρ, N)) ≈ fill(0.5, 4)
    @test trace(ρ) ≈ 1
end

@testset "The noise of the sweeps" begin
    # a number, or one value per sweep continued with its last, stored as the Float64 dmrg takes
    @test Dmrg(nsweeps = 2, noise = 0).noise === 0.
    @test Dmrg(nsweeps = 2, noise = [1, 0]).noise == [1., 0.]
    @test Dmrg(nsweeps = 2, noise = [1, 0]).noise isa Vector{Float64}
    L = Dissipator(Sp)(1) + Dissipator(Sm)(2)
    ρ0 = State{Mixed}(System(2, Qubit()), "FullyMixed")
    lim = Limits(maxdim = 16)
    for noise in (1e-6, [1e-6, 0.])
        _, ρ = steady_state(L, ρ0; algo = Dmrg(; nsweeps = 2, noise), limits = lim)
        @test real(expect1(ρ, Z)) ≈ [1, -1]
    end
    h = -sum(Z(i) for i in 1:4) - 0.3 * sum(X(i) * X(i + 1) for i in 1:3)
    r = random_state(Pure, System(4, Qubit()), 4)
    e, _ = ground_state(h, r; algo = Dmrg(nsweeps = 6), limits = Limits(maxdim = 16))
    @test first(ground_state(h, r; algo = Dmrg(nsweeps = 6, noise = [1e-5, 1e-7, 0.]),
                             limits = Limits(maxdim = 16))) ≈ e
end

@testset "The Krylov parameters reach the searches" begin
    # a Krylov space of a single vector holds nothing but the state it starts from, so that
    # the search does not leave it: dmrg ends where it started
    n = 4
    h = sum(X(i) * X(i + 1) for i in 1:n - 1)
    up = State{Pure}(System(n, Qubit()), "Up")
    lim = Limits(maxdim = 16)
    @test first(ground_state(h, up; algo = Dmrg(nsweeps = 3), limits = lim)) ≈ -3
    e, st = ground_state(h, up; algo = Dmrg(nsweeps = 3, krylov = Krylov(dim = 1)), limits = lim)
    @test e ≈ 0 atol = 1e-12
    @test real(expect1(st, Z)) ≈ ones(n)
    # and steady_state, from the fully mixed state; a single vector with its default of three
    # Krylov spaces fails in KrylovKit, which cannot restart from it
    L = Dissipator(Sp)(1) + Dissipator(Sm)(2)
    ρ0 = State{Mixed}(System(2, Qubit()), "FullyMixed")
    _, ρ = steady_state(L, ρ0; algo = Dmrg(nsweeps = 2, krylov = Krylov(dim = 1, maxiter = 1)),
                        limits = lim)
    @test real(expect1(ρ, Z)) ≈ [0, 0] atol = 1e-12
    # which the GroundState and SteadyState phases hand on, the latter with the same defaults
    sim = runTMS(SimData(phases = [
        CreateState{Pure}(n, Qubit(), "Up"),
        GroundState(hamiltonian = h, algo = Dmrg(nsweeps = 3, krylov = Krylov(dim = 1)), limits = lim)]);
        output = devnull)
    @test sim.energy ≈ 0 atol = 1e-12
    sim = runTMS(SimData(phases = [
        CreateState{Mixed}(2, Qubit(), "FullyMixed"),
        SteadyState(lindbladian = L, algo = Dmrg(nsweeps = 2, krylov = Krylov(dim = 1, maxiter = 1)),
                    limits = lim)]); output = devnull)
    @test real(expect1(sim.state, Z)) ≈ [0, 0] atol = 1e-12
    L2 = -im * (X(1) + 0.7 * X(2)) + Dissipator(Sm)(1)
    sim = runTMS(SimData(phases = [CreateState{Mixed}(2, Qubit(), "FullyMixed"),
                                   SteadyState(lindbladian = L2, algo = Dmrg(nsweeps = 3), limits = lim)]);
                 output = devnull)
    @test norm(sim.state - last(steady_state(L2, ρ0; algo = Dmrg(nsweeps = 3), limits = lim))) < 1e-12

    # a field left nothing takes the default of the method, field by field: KrylovKit's prefixed
    # for dmrg, and dim = 8 and maxiter = 3 for steady_state, which Krylov(tol = …) keeps where
    # it replaced them in 1.x
    kw = TMS.krylov_kwargs
    @test kw(Krylov()) == (;)
    @test kw(Krylov(dim = 8, maxiter = 3, tol = 1e-10)) == (krylovdim = 8, maxiter = 3, tol = 1e-10)
    @test kw(Krylov(dim = 8), Krylov(), "eigsolve_") == (eigsolve_krylovdim = 8,)
    @test kw(Krylov(tol = 1e-10), TMS.steady_krylov, "eigsolve_") ==
          (eigsolve_krylovdim = 8, eigsolve_maxiter = 3, eigsolve_tol = 1e-10)
    @test kw(Krylov(dim = 4), TMS.steady_krylov, "eigsolve_") ==
          (eigsolve_krylovdim = 4, eigsolve_maxiter = 3)
    # which the chain pumped up of "Steady state" shows in ten sweeps: the parameters of
    # ITensorMPS leave an error of a few percent, the default of steady_state 1e-9, and a
    # tolerance alone, the default of ITensorMPS, keeps the other two
    L = -im * (-sum(Z(i) * Z(i+1) for i in 1:4)) + sum(Dissipator(Sp)(i) for i in 1:5)
    r = random_state(Mixed, System(5, Qubit()), 10)
    z(krylov) = real(expect1(last(steady_state(L, r; algo = Dmrg(nsweeps = 10; krylov),
                                               limits = Limits(maxdim = 10, cutoff = 1e-10))), Z))
    @test z(Krylov()) ≈ ones(5) atol = 1e-7
    @test z(Krylov(tol = 1e-14)) ≈ z(Krylov()) atol = 1e-12
    @test maximum(abs.(z(Krylov(dim = 3, maxiter = 1)) .- 1)) > 1e-3
end

@testset "Ground and steady states on a charged system" begin
    h = -sum(dag(C)(i) * C(i+1) + dag(C)(i+1) * C(i) for i in 1:3)
    lind = -im * h + sum(Dissipator(sqrt(0.3) * N)(i) for i in 1:4)
    lim = Limits(cutoff = 1e-12, maxdim = 32)
    function start(site)
        sys = System(4, site)
        p(v) = State{Pure}(sys, v)
        return (p(["Occ", "Emp", "Occ", "Emp"]) + 0.5 * p(["Emp", "Occ", "Occ", "Emp"])) / sqrt(1.25)
    end
    # dmrg searches inside the sector its starting state lives in, and steady_state builds the
    # MPO of L†L, whose flux is zero whenever that of L is. Both must land where the
    # computation on the system conserving nothing does
    s0 = start(Fermion())
    e = first(ground_state(h, s0; algo = Dmrg(nsweeps = 3), limits = lim))
    # the two modes of negative energy filled, -2cos(π/5) - 2cos(2π/5)
    @test e ≈ -sqrt(5)
    z = real(first(steady_state(lind, mix(s0); algo = Dmrg(nsweeps = 2), limits = lim)))
    # an operator prepared once, for the system it was prepared on
    @test first(ground_state(PreMPO(s0, h), s0; algo = Dmrg(nsweeps = 3), limits = lim)) ≈ e
    ρ0 = mix(s0)
    @test real(first(steady_state(PreMPO(ρ0, lind), ρ0; algo = Dmrg(nsweeps = 2), limits = lim))) ≈ z atol = 1e-10
    for site in (Fermion(conserve = N), Fermion(conserve = strong(N)))
        @test first(ground_state(h, start(site); algo = Dmrg(nsweeps = 3), limits = lim)) ≈ e
        @test real(first(steady_state(lind, mix(start(site)); algo = Dmrg(nsweeps = 2), limits = lim))) ≈
              z atol = 1e-10
    end
end

@testset "The MPOs kept in the PreMPO" begin
    # ground_state keeps the MPO of the Hamiltonian, which a second search takes again, and
    # steady_state that of L†L for its limits and its algorithm
    n = 4
    h = sum(X(i) * X(i + 1) for i in 1:n - 1) + 0.5 * sum(Z(i) for i in 1:n)
    st = State{Pure}(System(n, Qubit()), "Up")
    pre = PreMPO(st, h)
    e1, _ = ground_state(pre, st; algo = Dmrg(nsweeps = 4), limits = Limits(maxdim = 16))
    mpo = last(pre.cache[:mpo])
    e2, _ = ground_state(pre, st; algo = Dmrg(nsweeps = 4), limits = Limits(maxdim = 16))
    @test last(pre.cache[:mpo]) === mpo
    @test e1 ≈ e2
    L = Dissipator(Sp)(1) + Dissipator(Sm)(2) - im * X(1) * X(2)
    ρ0 = State{Mixed}(System(2, Qubit()), "FullyMixed")
    pl = PreMPO(ρ0, L)
    _, ρ1 = steady_state(pl, ρ0; algo = Dmrg(nsweeps = 4))
    key, l2 = pl.cache[:squared_lindbladian]
    @test key == (Limits(), "naive")
    _, ρ2 = steady_state(pl, ρ0; algo = Dmrg(nsweeps = 4))
    @test last(pl.cache[:squared_lindbladian]) === l2
    @test norm(ρ1 - ρ2) < 1e-10
    # another truncation of L†L, or another algorithm of the product, builds it again, the MPO
    # of L being kept
    l = last(pl.cache[:mpo])
    for (mpo_limits, mpo_algo) in [(Limits(cutoff = 1e-14), "naive"), (Limits(), "zipup")]
        _, ρ3 = steady_state(pl, ρ0; algo = Dmrg(nsweeps = 4), mpo_limits, mpo_algo)
        @test first(pl.cache[:squared_lindbladian]) == (mpo_limits, mpo_algo)
        @test last(pl.cache[:squared_lindbladian]) !== l2
        @test norm(ρ3 - ρ1) < 1e-8
    end
    @test last(pl.cache[:mpo]) === l
end

@testset "The methods and their parameters" begin
    @test Dmrg(nsweeps = 3) isa SearchAlgorithm
    @test repr(Dmrg(nsweeps = 3)) == "Dmrg(nsweeps = 3, noise = 0.0, krylov = Krylov(), tol = 0.0)"
    @test repr(Dmrg(nsweeps = 3, noise = [1e-6, 0], krylov = Krylov(dim = 8), tol = 1e-10)) ==
          "Dmrg(nsweeps = 3, noise = [1.0e-6, 0.0], krylov = Krylov(dim = 8), tol = 1.0e-10)"
    @test repr(Krylov()) == "Krylov()"
    @test repr(Krylov(dim = 8, maxiter = 3, tol = 1e-10)) == "Krylov(dim = 8, maxiter = 3, tol = 1.0e-10)"
    # a tolerance is a real number, an integer one included
    @test Krylov(tol = 1).tol === 1.
    @test Dmrg(nsweeps = 2, tol = 1).tol === 1.
    # the number of sweeps has no default, and a search takes one at least
    @test_throws UndefKeywordError Dmrg()
    for nsweeps in (0, -1)
        @test_throws "a search takes at least one sweep, and nsweeps is $nsweeps" Dmrg(; nsweeps)
    end
end

@testset "What the searches refuse" begin
    h = -Z(1) * Z(2) - 0.5 * (X(1) + X(2))
    st = State{Pure}(System(2, Qubit()), "Up")
    ρ = mix(random_state(Pure, System(2, Qubit()), 2))
    algo = Dmrg(nsweeps = 2)
    # the method has no default, which the refusal writes
    @test_throws "ground_state needs its search method: algo = Dmrg(nsweeps = 10)" ground_state(h, st)
    @test_throws "steady_state needs its search method: algo = Dmrg(nsweeps = 10)" steady_state(
        Dissipator(Sm)(1), ρ)
    # a Hamiltonian on a mixed state, of which dmrg would minimise ρ ↦ Hρ + ρH with no
    # guarantee of ending on a density matrix, and a Lindbladian on a pure state
    @test_throws "ground_state searches a pure state: search it pure and mix it" ground_state(
        h, ρ; algo)
    @test_throws "steady_state searches a density matrix: mix the state" steady_state(
        Dissipator(Sm)(1), st; algo)
    # and so do the phases
    @test_throws "GroundState searches a pure state" runTMS(SimData(phases = [
        CreateState{Mixed}(2, Qubit(), "Up"), GroundState(; hamiltonian = h, algo)]); output = devnull)
    @test_throws "SteadyState searches a density matrix" runTMS(SimData(phases = [
        CreateState{Pure}(2, Qubit(), "Up"), SteadyState(; lindbladian = Dissipator(Sm)(1), algo)]);
        output = devnull)
    # an MPO, a PreMPO of another system or of the other representation, and what is no operator
    msg = "the solvers take an operator or its PreMPO, which keeps the MPOs they build, not an MPO"
    @test_throws msg ground_state(make_mpo(st, h), st; algo)
    @test_throws msg steady_state(make_mpo(ρ, -im * h), ρ; algo)
    other = System(3, Qubit())
    @test_throws "the PreMPO was prepared on another system than that of the state" ground_state(
        PreMPO(State{Pure}(other, "Up"), h), st; algo)
    @test_throws "the PreMPO was prepared on another system than that of the state" steady_state(
        PreMPO(State{Mixed}(other, "Up"), -im * h), ρ; algo)
    @test_throws "the PreMPO was prepared on a state of the other representation" ground_state(
        PreMPO(ρ, h), st; algo)
    @test_throws "the PreMPO was prepared on a state of the other representation" steady_state(
        PreMPO(st, h), ρ; algo)
    @test_throws "3 is not an operator or its PreMPO" ground_state(3, st; algo)
    @test_throws "mpo_algo is \"naive\" or \"zipup\", not \"fit\"" steady_state(
        Dissipator(Sm)(1), ρ; algo, mpo_algo = "fit")
    # the options of ITensorMPS are not passed on, `outputlevel` included
    @test_throws MethodError ground_state(X(1) * X(2), st; algo, outputlevel = 1)
    @test_throws MethodError steady_state(Dissipator(Sm)(1), ρ; algo, outputlevel = 1)
end

@testset "Excited states" begin
    # the bound of the width of the spectrum, by the inequalities of Weyl, term by term: exact on
    # a term of two sites, a constant counting for nothing, a positive term for its largest
    # eigenvalue, a hopping and its adjoint for the width of their sum, a term of more sites by
    # the norms of its blocks
    n = 6
    sys = System(n, Qubit())
    up = State{Pure}(sys, "Up")
    heis = sum(X(i) * X(i + 1) + Y(i) * Y(i + 1) + Z(i) * Z(i + 1) for i in 1:n-1)
    @test PreMPO(up, heis).width[] ≈ 6(n - 1)
    @test PreMPO(up, heis + 2.0 * Id(1)).width[] ≈ 6(n - 1)
    @test PreMPO(up, 2 * N(1) * N(2)).width[] ≈ 2
    @test PreMPO(up, Z(1) * Z(2) * Z(3)).width[] ≈ 2
    @test PreMPO(State{Pure}(System(4, Fermion()), "Emp"),
                 dag(C)(1) * C(3) + dag(C)(3) * C(1)).width[] ≈ 2
    # the string on the site of the first block counts, which a Majorana operator, odd and of no
    # definite form, shows: i a(1) a(3) is Y ⊗ Z ⊗ X, of width 2, where the bound was 0
    a = C + dag(C)
    @test PreMPO(State{Pure}(System(4, Fermion()), "Emp"), im * a(1) * a(3)).width[] ≈ 2
    # the lowest eigenstates in turn, each orthogonal to those before, the energy that of the
    # Hamiltonian alone, against its spectrum
    H = heis + 0.3 * sum(Z(i) for i in 1:n)
    LA = TMS.LinearAlgebra
    E = LA.eigvals(LA.Hermitian(Matrix(TMS.system_matrix(Pure, sys, TMS.lower(sys, H)))))
    @test PreMPO(up, H).width[] ≥ E[end] - E[1]
    r = random_state(Pure, sys, 4; rng = Xoshiro(1))
    lim = Limits(maxdim = 32)
    algo = Dmrg(nsweeps = 40, tol = 1e-13)
    e0, g = ground_state(H, r; algo, limits = lim)
    e1, x = ground_state(H, r; algo, limits = lim, orthogonal_to = [g])
    @test e0 ≈ E[1]
    @test e1 ≈ E[2] atol = 1e-9
    @test e1 ≈ real(expect(x, H))
    @test abs2(inner(g, x)) < 1e-12
    # the next one, with a weight above its gap, and one below, which falls back on the ground
    # state and warns
    e2, _ = ground_state(H, r; algo, limits = lim, orthogonal_to = [g, x], weight = 2(E[3] - E[1]))
    @test e2 ≈ E[3] atol = 1e-9
    @test_logs (:warn, r"keeps an overlap of .* the weight .* is below the gap") match_mode = :any ground_state(
        H, r; algo, limits = lim, orthogonal_to = [g], weight = 0.1 * (E[2] - E[1]))
    # a state of another sector of charge is orthogonal already, and left out
    qs = System(n, Qubit(conserve = N))
    a = State{Pure}(qs, [ isodd(i) ? "Up" : "Dn" for i in 1:n ])
    b = State{Pure}(qs, "Up")
    pre = PreMPO(a, heis)
    @test isempty(TMS.orthogonal_args(a, pre, [b], nothing))
    @test TMS.orthogonal_args(a, pre, [a], nothing).weight == pre.width[]
    @test TMS.orthogonal_args(a, pre, [a], 0.5).weight == 0.5
    @test_throws "on another system" ground_state(heis, a; algo, orthogonal_to = [up])
    @test_throws "orthogonal to pure states, not to a" ground_state(heis, a; algo, orthogonal_to = [mix(a)])
end
