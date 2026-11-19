# The representation Trajectory, a state held as one quantum trajectory: its types, the algebra it
# shares with Pure, its files, and what it refuses

const Trajectory = TensorMixedStates.Trajectory
const Ket = TensorMixedStates.Ket

@testset "Trajectory is a ket, as Pure, and holds its states as Pure does" begin
    @test Trajectory <: Ket && Pure <: Ket && !(Mixed <: Ket) && Ket <: TensorMixedStates.PM
    system = System(4, Fermion(conserve = strong(N)))
    @test all(siteind(Trajectory, system, i) == siteind(Pure, system, i) for i in 1:4)
    st = State{Trajectory}(system, [ isodd(i) ? "Occ" : "Emp" for i in 1:4 ])
    pure = State{Pure}(system, [ isodd(i) ? "Occ" : "Emp" for i in 1:4 ])
    @test st isa State{Trajectory} && TensorMixedStates.representation(st) === Trajectory
    overlap(a, b) = TensorMixedStates.ITensorMPS.inner(a.state, b.state)
    @test overlap(st, pure) ≈ 1
    mktempdir() do dir
        file = joinpath(dir, "t.h5")
        save_state(file, "t", st)
        back = load_state(file, "t")
        @test back isa State{Trajectory} && overlap(back, st) ≈ 1
    end
    @test_throws "ExactState holds Pure or Mixed states" ExactState{Trajectory}(System(2, Qubit()),
                                                                                "Up")
end

as_trajectory(p) = State{Trajectory}(p.system, p.state)

@testset "Trajectory measures what is linear in the density matrix exactly as Pure does" begin
    rng = Xoshiro(5)
    qubits = System(4, Qubit())
    fermions = System(4, Fermion(conserve = strong(N)))
    electrons = System(3, Electron(conserve = (Ntot, 2Sz)))
    cases = [
        (random_state(Pure, qubits, 4; rng), [X(1), Z], [(Z, Z)],
         sum(Z(i) * Z(i + 1) for i in 1:3) + sum(X(i) for i in 1:4), X(1) * X(2), Z(1)),
        (random_state(Pure, fermions, ["Occ", "Emp", "Occ", "Emp"], 4; rng), [N(2), N],
         [(dag(C), C)], sum(dag(C)(i) * C(i + 1) + dag(C)(i + 1) * C(i) for i in 1:3),
         dag(C)(1) * C(2) + dag(C)(2) * C(1), N(1)),
        (random_state(Pure, electrons, ["Up", "Dn", "Emp"], 4; rng), [Ntot(1), Nup],
         [(dag(Cup), Cup)],
         sum(dag(Cup)(i) * Cup(i + 1) + dag(Cup)(i + 1) * Cup(i) + Ntot(i) * Ntot(i + 1)
             for i in 1:2), Nupdn(2), Ntot(1)),
    ]
    for (p, ops, pairs, h, op, measured) in cases
        t = as_trajectory(p)
        @test expect(t, first(ops)) == expect(p, first(ops))
        @test expect1(t, last(ops)) == expect1(p, last(ops))
        @test expect2(t, only(pairs)) == expect2(p, only(pairs))
        @test probabilities(t, measured) == probabilities(p, measured)
        @test trace(t) == trace(p)
        @test reduced_density_matrix(t, [1, 2]) == reduced_density_matrix(p, [1, 2])
        @test sample(t; rng = Xoshiro(3)) == sample(p; rng = Xoshiro(3))
        ct, cp = collapse(t, measured; rng = Xoshiro(4)), collapse(p, measured; rng = Xoshiro(4))
        @test ct.result == cp.result && ct.probability == cp.probability
        @test ct.state isa State{Trajectory} && expect1(ct.state, last(ops)) ==
                                                 expect1(cp.state, last(ops))
        at, ap = apply(op, t), apply(op, p)
        @test at isa State{Trajectory} && expect1(at, last(ops)) == expect1(ap, last(ops))
        lim = Limits(maxdim = 8)
        et = evolve(-im * h, 0.1, t; nsteps = 2, limits = lim)
        ep = evolve(-im * h, 0.1, p; nsteps = 2, limits = lim)
        @test et isa State{Trajectory} && expect1(et, last(ops)) == expect1(ep, last(ops))
        @test measure(t, Measure(Trace, last(ops)); time = 0.) ==
              measure(p, Measure(Trace, last(ops)); time = 0.)
    end
    # a fidelity to a pure state is ⟨φ|ρ|φ⟩, linear
    p = random_state(Pure, qubits, 4; rng)
    q = random_state(Pure, qubits, 4; rng)
    t = as_trajectory(p)
    @test fidelity(t, q) == fidelity(p, q) && fidelity(q, t) == fidelity(q, p)
    # a pure state is taken as a trajectory
    c = TensorMixedStates.as_representation(Trajectory, p)
    @test c isa State{Trajectory} && expect1(c, Z) == expect1(p, Z)
end

@testset "Trajectory refuses what is not linear in the density matrix, and what needs draws" begin
    rng = Xoshiro(6)
    qubits = System(4, Qubit())
    p = random_state(Pure, qubits, 4; rng)
    t = as_trajectory(p)
    h = sum(Z(i) * Z(i + 1) for i in 1:3) + sum(X(i) for i in 1:4)
    TMS = TensorMixedStates
    for (f, message) in (
            (() -> purity(t), "purity is not linear"),
            (() -> renyi2(t), "purity is not linear"),
            (() -> renyi2(t, [1, 2]), "Renyi 2 entropy is not linear"),
            (() -> mutual_info_renyi2(t, [1]), "Renyi 2 entropy is not linear"),
            (() -> hermiticity(t), "hermiticity is not linear"),
            (() -> entanglement_entropy(t, 2), "entanglement entropy or spectrum"),
            (() -> entanglement_spectrum(t, 2), "entanglement entropy or spectrum"),
            (() -> entanglement_by_sector(t, 2), "entanglement entropy is not linear"),
            (() -> von_neumann_entropy(t, [1]), "von Neumann entropy is not linear"),
            (() -> log_negativity(t, [1], [2]), "logarithmic negativity is not linear"),
            (() -> inner(t, t), "inner product is not linear"),
            (() -> inner(t, p), "inner product is not linear"),
            (() -> inner(p, t), "inner product is not linear"),
            (() -> inner(t, X(1), t), "matrix element is not linear"),
            (() -> hs_fidelity(t, t), "hs_fidelity is not linear"),
            (() -> fidelity(t, t), "fidelity but to a pure state"),
            (() -> fidelity(t, mix(p)), "fidelity but to a pure state"),
            (() -> variance(t, X(1)), "measure A^2 and A"),
            (() -> measure(t, Measure(Norm); time = 0.), "norm is not linear"),
            (() -> measure(t, Measure(Overlap(p)); time = 0.), "not linear"),
            (() -> mix(t), "a trajectory is not a mixed state"),
            (() -> partial_trace(t; keep = [1]), "a trajectory is not a mixed state"),
            (() -> Vector(t), "the vector of one trajectory"),
            (() -> t ⊗ t, "the product of two trajectories"),
            (() -> t ⊗ p, "the product of a trajectory"),
            (() -> random_state(t, 2), "a random state comes with the draws"),
            (() -> random_state(Trajectory, qubits, 2), "a random state comes with the draws"),
            (() -> State{Trajectory}(qubits, Mixture([0.5 => "Up", 0.5 => "Dn"])),
             "a mixed initial state comes with the draws"),
            (() -> apply(Reset("Up")(1), t), "more than one Kraus operator on a trajectory"),
            (() -> ground_state(h, t), "a trajectory is not searched"),
            (() -> thermalize(h, 1., t), "a thermal state comes with the draws"),
            (() -> steady_state(-im * h, t), "a steady state comes with the draws"),
            (() -> ExactState(t), "ExactState holds Pure or Mixed states"),
            (() -> TMS.as_representation(Mixed, t), "a trajectory is not a pure or a mixed"),
            (() -> TMS.as_representation(Pure, t), "a trajectory is not a pure or a mixed"),
            (() -> TMS.as_representation(Trajectory, mix(p)), "a trajectory of a mixed state"))
        @test_throws message f()
    end
end

@testset "Trajectory weakens, compares and evolves by a PreMPO as Pure does" begin
    rng = Xoshiro(7)
    fermions = System(4, Fermion(conserve = strong(N)))
    p = random_state(Pure, fermions, ["Occ", "Emp", "Occ", "Emp"], 4; rng)
    t = as_trajectory(p)
    for target in ((), (N,))
        w = weaken(t, target)
        @test w isa State{Trajectory} && expect1(w, N) == expect1(weaken(p, target), N)
    end
    q = random_state(Pure, fermions, ["Emp", "Occ", "Occ", "Emp"], 4; rng)
    @test fidelity(t, weaken(q, ())) == fidelity(p, weaken(q, ()))
    # an operator measured on the state conserving nothing
    spins = System(6, Spin(1 / 2; conserve = 2Sz))
    ps = random_state(Pure, spins, ["Up", "Dn", "Up", "Dn", "Up", "Dn"], 4; rng)
    op = (Sx(1) + Sx(2)) * (Sx(3) + Sx(4)) * (Sx(5) + Sx(6))
    @test expect(as_trajectory(ps), op) == expect(ps, op)
    # a PreMPO prepared on a trajectory is for its wave function
    h = sum(dag(C)(i) * C(i + 1) + dag(C)(i + 1) * C(i) for i in 1:3)
    et = evolve(PreMPO(t, -im * h), 0.2, t; nsteps = 2)
    @test et isa State{Trajectory} && expect1(et, N) == expect1(evolve(PreMPO(p, -im * h), 0.2,
                                                                       p; nsteps = 2), N)
    # in a run, weakened then evolved
    run(R) = runTMS(SimData(; name = "w", phases = [
                 CreateState{R}(fermions, ["Occ", "Emp", "Occ", "Emp"]), Weaken(target = ()),
                 Evolve(evolver = -im * h, duration = 0.2, time_step = 0.1)],
                 final_measurements = Data("n") => N); output = devnull)
    @test run(Trajectory).data["n"] == run(Pure).data["n"]
end
