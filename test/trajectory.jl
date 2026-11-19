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

@testset "Trajectory leaves the dissipators, the channels and the exact states for later" begin
    rng = Xoshiro(8)
    qubits = System(3, Qubit())
    p = random_state(Pure, qubits, 4; rng)
    t = as_trajectory(p)
    h = sum(Z(i) * Z(i + 1) for i in 1:2) + sum(X(i) for i in 1:3)
    @test_throws "a dissipator on a trajectory comes with the draws" evolve(
        -im * h + Dissipator(X)(1), 0.1, t; nsteps = 1)
    @test_throws "a channel on a trajectory comes with the draws" evolve(
        depolarizing_channel(0.2)(1), 0.1, t; nsteps = 1)
    steady(R) = runTMS(SimData(; name = "s", phases = [CreateState{R}(qubits, "Up"),
                       SteadyState(lindbladian = -im * h, algo = Dmrg(nsweeps = 1))]);
                       output = devnull)
    @test_throws "a steady state comes with the draws" steady(Trajectory)
    # the exact states hold none, but give the fidelity to a pure state
    e = ExactState(p)
    @test_throws "ExactState holds Pure or Mixed states" TensorMixedStates.as_representation(
        Trajectory, e)
    @test_throws "ExactState holds Pure or Mixed states" random_state(ExactState{Trajectory},
                                                                      qubits, 2)
    convert(R) = runTMS(SimData(; name = "c", phases = [
                        CreateState{Pure}(qubits, "Up"; backend = ExactState),
                        ConvertState{R}()]); output = devnull)
    @test_throws "ExactState holds Pure or Mixed states" convert(Trajectory)
    q = random_state(Pure, qubits, 4; rng)
    @test fidelity(t, ExactState(q)) ≈ fidelity(t, q) rtol = 1e-12
    @test fidelity(ExactState(q), t) ≈ fidelity(q, t) rtol = 1e-12
    @test_throws "inner product is not linear" inner(t, e)
    @test_throws "fidelity but to a pure state" fidelity(t, mix(e))
end

@testset "Trajectory is built from a pure description, and State from a concrete one" begin
    qubits = System(2, Qubit())
    t = State{Trajectory}(qubits, Mixture([1.0 => "Up", 0 => "Dn"]))
    @test t isa State{Trajectory} && expect1(t, Z) == expect1(State{Pure}(qubits, "Up"), Z)
    @test occursin("State{Trajectory}", repr(t))
    @test_throws "a mixed initial state comes with the draws" State{Trajectory}(
        qubits, Mixture([0.5 => "Up", 0.5 => "Dn"]))
    @test_throws "a mixed initial state comes with the draws" State{Trajectory}(
        qubits, FullyMixed())
    @test_throws "a mixed initial state comes with the draws" State{Trajectory}(
        qubits, [[1 0; 0 0], "Up"])
    @test_throws "State{Ket} has no representation of its own" State{Ket}(qubits, "Up")
end

# whether the mean of the vectors `values` is `ref`, within `k` of its standard errors
function within_error(values::Vector{<:Vector}, ref::Vector; k = 4)
    n = length(values)
    m = sum(values) / n
    err = [ sqrt(sum(abs2(v[i] - m[i]) for v in values) / (n - 1) / n) for i in eachindex(ref) ]
    return all(abs.(m .- ref) .≤ k .* err .+ 1e-12)
end

@testset "Trajectory draws the Kraus operators of a channel, whose average is the channel" begin
    rng = Xoshiro(9)
    qubits = System(3, Qubit())
    p = random_state(Pure, qubits, 4; rng)
    t = as_trajectory(p)
    for kraus in (:reset, :pauli)
        ch = depolarizing_channel(0.4; kraus)(1) * amplitude_damping_channel(0.5, Sp)(2) *
             reset_channel(0.3, "Up")(3) * controlled(X)(1, 2)
        ref = real(expect1(apply(ch, mix(p)), Z))
        values = [ real(expect1(apply(ch, t; rng), Z)) for _ in 1:400 ]
        @test within_error(values, ref)
        same(seed) = expect1(apply(ch, t; rng = Xoshiro(seed)), Z)
        @test same(1) == same(1)
    end
    # on fermions conserving their number weakly, each trajectory keeps a definite number
    fermions = System(4, Fermion(conserve = N))
    pf = random_state(Pure, fermions, ["Occ", "Emp", "Occ", "Emp"], 4; rng)
    ch = amplitude_damping_channel(0.5, C)(2) * dephasing_channel(0.3)(3) *
         reset_channel(0.4, "FullyMixed")(1)
    ref = real(expect1(apply(ch, mix(pf)), N))
    values = [ real(expect1(apply(ch, as_trajectory(pf); rng), N)) for _ in 1:400 ]
    @test within_error(values, ref)
    # the Kraus operators as written: the identity and the Pauli strings, or the identity and the
    # resets to each basis state
    ks(c) = TensorMixedStates.kraus_operators(qubits, c(1))
    @test length(ks(depolarizing_channel(0.2; kraus = :pauli))) == 4
    @test length(ks(depolarizing_channel(0.2))) == 5
    @test length(TensorMixedStates.kraus_operators(qubits,
                 depolarizing_channel(0.2, 2; kraus = :pauli)(1, 2))) == 16
    # refused: Kraus operators of no definite charge, a channel not keeping the trace
    @test_throws "as written have no definite charge" apply(
        0.5 * (Gate(C + dag(C)) + Gate(im * (C - dag(C))))(2), as_trajectory(pf); rng)
    @test_throws "does not keep the trace" apply(0.5 * depolarizing_channel(0.2)(1), t; rng)
    @test_throws "does not keep the trace" apply(Gate(Proj("Up"))(1) + Gate(Proj("Up"))(1), t;
                                                 rng)
    @test_throws "kraus is :reset or :pauli" depolarizing_channel(0.2; kraus = :other)
    # drawn as written, never developed: a power, and a channel of many sites
    ch = (dephasing_channel(0.2, X)^10)(1)
    @test within_error([ real(expect1(apply(ch, t; rng), X)) for _ in 1:200 ],
                       real(expect1(apply(ch, mix(p)), X)))
    seven = System(7, Qubit())
    p7 = random_state(Pure, seven, 4; rng)
    ch = depolarizing_channel(0.3, 7)(1:7...)
    @test within_error([ real(expect1(apply(ch, as_trajectory(p7); rng), Z)) for _ in 1:200 ],
                       real(expect1(apply(ch, mix(p7)), Z)))
    # in a run, the generator of the simulation draws them
    gates(name) = runTMS(SimData(; name, seed = 3, phases = [
                    CreateState{Trajectory}(qubits, "Up"),
                    Gates(gates = depolarizing_channel(0.5)(1) * depolarizing_channel(0.5)(2))],
                    final_measurements = Data("z") => Z); output = devnull).data["z"]
    @test gates("a") == gates("b")
end
