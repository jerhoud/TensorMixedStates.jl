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
            (() -> random_state(t, 2), "randomising a trajectory is not implemented"),
            (() -> ground_state(h, t), "a trajectory is not searched"),
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
    run(R, d; kwargs...) = runTMS(SimData(; name = "w", kwargs..., phases = [
                 CreateState{R}(fermions, ["Occ", "Emp", "Occ", "Emp"]), Weaken(target = ()),
                 Evolve(evolver = -im * h, duration = 0.2, time_step = 0.1)],
                 final_measurements = d => N); output = devnull).data["n"][1]["N"]
    @test run(Trajectory, EachSample(Data("n")); samples = 1) == run(Pure, Data("n"))
end

@testset "Trajectory leaves the dissipators, the channels and the exact states for later" begin
    rng = Xoshiro(8)
    qubits = System(3, Qubit())
    p = random_state(Pure, qubits, 4; rng)
    t = as_trajectory(p)
    h = sum(Z(i) * Z(i + 1) for i in 1:2) + sum(X(i) for i in 1:3)
    @test_throws "given the generator rather than its PreMPO" PreMPO(t, -im * h +
                                                                      Dissipator(X)(1))
    @test_throws "a channel on a trajectory comes with the draws" evolve(
        depolarizing_channel(0.2)(1), 0.1, t; nsteps = 1)
    steady(R) = runTMS(SimData(; name = "s", samples = 1, phases = [
                       CreateState{R}(qubits, "Up"),
                       SteadyState(lindbladian = -im * h, algo = Dmrg(nsweeps = 1))]);
                       output = devnull)
    @test_throws "a steady state comes with the draws" steady(Trajectory)
    # the exact states hold none, but give the fidelity to a pure state
    e = ExactState(p)
    @test_throws "ExactState holds Pure or Mixed states" TensorMixedStates.as_representation(
        Trajectory, e)
    @test_throws "ExactState holds Pure or Mixed states" random_state(ExactState{Trajectory},
                                                                      qubits, 2)
    convert(R) = runTMS(SimData(; name = "c", samples = 1, phases = [
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
    # a pure state as a trajectory, its MPS copied
    ψ = random_state(Pure, qubits, 2; rng = Xoshiro(1))
    tψ = State{Trajectory}(ψ)
    @test tψ isa State{Trajectory} && expect1(tψ, X) == expect1(ψ, X)
    @test tψ.state !== ψ.state
    @test expect1(State{Trajectory}(tψ; rng = Xoshiro(2)), X) == expect1(ψ, X)
    @test_throws "a trajectory of a mixed state" State{Trajectory}(mix(ψ))
    # a projector of rank one is a pure state, whatever the draw
    @test expect1(State{Trajectory}(qubits, [[1 0; 0 0], "Up"]), Z) == [1, 1]
    @test_throws "State{Ket} has no representation of its own" State{Ket}(qubits, "Up")
    # the generator taken by every representation, which only a trajectory draws with
    for R in (Pure, Mixed), x in ("Up", Mixture([1.0 => "Up"]))
        @test expect1(State{R}(qubits, x; rng = Xoshiro(1)), Z) == expect1(State{R}(qubits, x), Z)
    end
    @test expect1(State{Mixed}(qubits, FullyMixed(); rng = Xoshiro(1)), Z) == [0, 0]
    ψ = State{Pure}(qubits, "Up")
    @test State{Pure}(qubits, ψ.state; rng = Xoshiro(1)).state === ψ.state
    run(name, state) = runTMS(SimData(; name, phases = [CreateState{Pure}(qubits, state)],
                                      final_measurements = Data("z") => Z); output = devnull)
    @test run("mps", ψ.state).data["z"] == run("name", "Up").data["z"]
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
        values = [ real(expect1(apply(ch, t; rng), Z)) for _ in 1:200 ]
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
    values = [ real(expect1(apply(ch, as_trajectory(pf); rng), N)) for _ in 1:200 ]
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
    @test within_error([ real(expect1(apply(ch, t; rng), X)) for _ in 1:100 ],
                       real(expect1(apply(ch, mix(p)), X)))
    seven = System(7, Qubit())
    p7 = random_state(Pure, seven, 4; rng)
    ch = depolarizing_channel(0.3, 7)(1:7...)
    @test within_error([ real(expect1(apply(ch, as_trajectory(p7); rng), Z)) for _ in 1:100 ],
                       real(expect1(apply(ch, mix(p7)), Z)))
    # in a run, the generator of the simulation draws them
    gates(name) = runTMS(SimData(; name, seed = 3, samples = 1, phases = [
                    CreateState{Trajectory}(qubits, "Up"),
                    Gates(gates = depolarizing_channel(0.5)(1) * depolarizing_channel(0.5)(2))],
                    final_measurements = EachSample(Data("z")) => Z); output = devnull).data["z"]
    @test gates("a") == gates("b")
end

# whether the density matrix averaged over the trajectories `ts` is that of `ρ`, element by
# element within `k` of its standard errors, at least `floor`: 1 / length(ts) for an element
# that the trajectories may never reach
function same_density(ts, ρ; k = 5, floor = 1e-10)
    vs = [ (v = dense_vector(State{Pure}(t.system, t.state)); v * v' / real(v' * v)) for t in ts ]
    n = length(vs)
    m = sum(vs) / n
    err = sqrt.(sum(v -> abs2.(v .- m), vs) / (n - 1) / n)
    ref = ket_bra_matrix(ρ)
    return all(abs.(m .- ref) .≤ k .* max.(err, floor))
end

@testset "Trajectory draws products in sums, coherences and charges as Mixed does" begin
    rng = Xoshiro(11)
    qubits = System(2, Qubit())
    p = random_state(Pure, qubits, 4; rng)
    t = as_trajectory(p)
    # a weak measurement with feedback, its parts products of gates of unequal weights
    ε = 0.8
    Mp, Mm = ((Id + ε * Z) / 2)^0.5, ((Id - ε * Z) / 2)^0.5
    for ch in (Gate(Mp)(1) * Gate(X)(2) + Gate(Mm)(1) * Gate(Z)(2),
               Gate(Proj("Up"))(1) * Gate(X)(2) + Gate(Proj("Dn"))(1),
               depolarizing_channel(0.4; kraus = :pauli)(1) * reset_channel(0.3, "Up")(2))
        @test same_density([ apply(ch, t; rng) for _ in 1:300 ], apply(ch, mix(p)))
    end
    # fermions: coherences and Jordan-Wigner signs, a reset to a degenerate target with rounding
    fermions = System(3, Fermion(conserve = N))
    pf = random_state(Pure, fermions, ["Occ", "Emp", "Occ"], 4; rng)
    ch = amplitude_damping_channel(0.5, C)(3) * Reset([0.5 1e-17 ; 1e-17 0.5])(2)
    @test same_density([ apply(ch, as_trajectory(pf); rng) for _ in 1:300 ], apply(ch, mix(pf)))
    # refused as Mixed refuses: Kraus operators changing a strong charge by different amounts
    strong_sys = System(3, Fermion(conserve = strong(N)))
    ps = random_state(Pure, strong_sys, ["Occ", "Emp", "Occ"], 4; rng)
    @test_throws "by different amounts" apply(amplitude_damping_channel(0.5, C)(2),
                                              as_trajectory(ps); rng)
    @test_throws "by different amounts" apply(Gate(C)(1) + Gate(dag(C))(2), as_trajectory(ps);
                                              rng)
    # refused: off the system, a reset to a state of no definite charge of a quantity conserved,
    # no sites
    @test_throws "which the system does not have" apply(depolarizing_channel(0.3)(5), t; rng)
    paired = as_trajectory(State{Pure}(System(2, Fermion(conserve = parity(N))), "Occ"))
    @test_throws "which a trajectory cannot hold" apply(Reset([1, 1])(1), paired; rng)
    @test_throws "1 site or more" depolarizing_channel(0.1, 0)
    # a run draws: a full depolarizing of |↑⟩ gives ↓ for some seed
    run(seed) = only(runTMS(SimData(; name = "d$seed", seed, samples = 1, phases = [
                    CreateState{Trajectory}(System(1, Qubit()), "Up"),
                    Gates(gates = depolarizing_channel(1.)(1))],
                    final_measurements = Data("z") => Z(1)); output = devnull).data["z"])["Z(1)"]
    @test any(seed -> run(seed) ≈ -1, 1:20)
end

# the trajectory of the scheme of `evolve` from the vector `ψ`, drawn with `rng` as it draws: `ks`
# the effective operators `(f, K)` of the terms, `decays` the jumps `(f, [L])` of those with jumps,
# `f` their function of time or `nothing`
function dense_trajectory(ψ, ks, decays, duration, nsteps, rng)
    dt = duration / nsteps
    c(f, t) = isnothing(f) ? 1.0 : f(t)
    pick(ps) = something(findfirst(≥(rand(rng) * sum(ps)), cumsum(ps)), length(ps))
    for k in 0:nsteps - 1
        t = k * dt
        rates = [ c(f, t) * real(ψ' * sum(L' * L for L in Ls) * ψ) / real(ψ' * ψ)
                  for (f, Ls) in decays ]
        n = LA.norm(ψ)
        if rand(rng) ≥ dt * sum(rates)
            ψ = exp(dt * sum(c(f, t + dt / 2) * K for (f, K) in ks)) * ψ
        else
            Ls = last(decays[pick(rates)])
            ψ = Ls[pick([ real(ψ' * L' * L * ψ) for L in Ls ])] * ψ
        end
        ψ *= n / LA.norm(ψ)
    end
    return ψ
end

# whether the trajectories evolved by `evolve` from `t` with the seeds 1:n are those of
# `dense_trajectory`
function same_trajectories(evolver, t, ks, decays; duration = 0.5, nsteps = 10, n = 10)
    ψ = raw_vector(State{Pure}(t.system, t.state))
    return all(1:n) do seed
        a = raw_vector(State{Pure}(t.system, evolve(evolver, duration, t; nsteps,
                                                    rng = Xoshiro(seed)).state))
        b = dense_trajectory(ψ, ks, decays, duration, nsteps, Xoshiro(seed))
        isapprox(a, b; rtol = 1e-8)
    end
end

struct HalfJumps <: TMS.EvolutionAlgorithm end

TMS.evolve_step(::HalfJumps, pre, state, time, dt, step, limits) =
    TMS.evolve_step(Tdvp(), pre, TMS.evolve_step(Tdvp(), pre, state, time, dt / 2, step, limits),
                    time + dt / 2, dt / 2, step, limits)

@testset "Trajectory draws the jumps of the dissipators, whose average is the evolution" begin
    rng = Xoshiro(12)
    qubits = System(3, Qubit())
    p = random_state(Pure, qubits, 4; rng)
    t = as_trajectory(p)
    h = sum(X(i) for i in 1:3) + sum(Z(i) * Z(i + 1) for i in 1:2)
    m(o) = lowered(qubits, o)
    # jumps of operators, of a channel but its identity, and of a collective operator; a rate
    # depending on time, of a term of its own
    dep = depolarizing_channel(0.5)(1)
    kdep = [ sqrt(0.4) * m(k) for k in TMS.kraus_operators(qubits, dep)
             if !(m(k) ≈ m(k)[1, 1] * LA.I) ]
    @test length(kdep) == 4
    f = s -> 1 + s
    g = [-im * h + 0.5 * Dissipator(Sm)(1) + 0.4 * Dissipator(dep) +
         0.3 * Dissipator(Sm(2) + Sm(3)), Dissipator(Z)(2) => f]
    jumps = [ [[sqrt(0.5) * m(Sm(1))] ; kdep ; [sqrt(0.3) * m(Sm(2) + Sm(3))]], [m(Z(2))] ]
    effective(Ls) = -0.5 * sum(L' * L for L in Ls)
    @test same_trajectories(g, t, [(nothing, -im * m(h) + effective(jumps[1])),
                                   (f, effective(jumps[2]))],
                            [(nothing, jumps[1]), (f, jumps[2])])
    # fermions: Jordan-Wigner strings in a jump across a site, a reset
    fermions = System(3, Fermion(conserve = N))
    pf = random_state(Pure, fermions, ["Occ", "Emp", "Occ"], 4; rng)
    tf = as_trajectory(pf)
    c(k) = jw_placed(fermions, k, matrix(C, fermions[k]))
    hf = sum(dag(C)(i) * C(i + 1) + dag(C)(i + 1) * C(i) for i in 1:2)
    hop = c(1)' * c(3)
    kres = [ sqrt(0.5) * lowered(fermions, k)
             for k in TMS.kraus_operators(fermions, Reset("Occ")(1)) ]
    gf = -im * hf + 0.5 * Dissipator(Reset("Occ"))(1) + 0.4 * Dissipator(dag(C)(1) * C(3)) +
         Dissipator(C)(2)
    jf = [ kres ; [sqrt(0.4) * hop, c(2)] ]
    @test same_trajectories(gf, tf, [(nothing, -im * lowered(fermions, hf) + effective(jf))],
                            [(nothing, jf)])
    # the average, at first order in the step
    ts = [ evolve(g, 0.5, t; nsteps = 10, rng) for _ in 1:100 ]
    @test all(x -> real(trace(x)) ≈ 1, ts)
    @test same_density(ts, evolve(g, 0.5, mix(p); nsteps = 10); floor = 1 / 100)
    # the identity of a channel is no jump: equal dissipators give equal trajectories
    same(d) = evolve(-im * h + d, 0.1, t; nsteps = 5, rng = Xoshiro(2))
    vector(x) = raw_vector(State{Pure}(x.system, x.state))
    @test vector(same(Dissipator(100 * depolarizing_channel(0.01))(1))) ≈
          vector(same(depolarizing_dissipator(1)(1)))
    # an algorithm of one's own drawing through the steps of Tdvp
    @test vector(evolve(g, 0.5, t; nsteps = 5, algo = HalfJumps(), rng = Xoshiro(3))) ≈
          vector(evolve(g, 0.5, t; nsteps = 10, rng = Xoshiro(3)))
    # the trace kept by the operators of the terms of one function together
    @test real(trace(evolve([-im * h + X(1), -X(1), Dissipator(Sm)(1)], 0.1, t; nsteps = 2,
                            rng))) ≈ 1
    # refused: a jump changing a strong charge, as Mixed refuses it, negative rates, a trace the
    # operators do not keep, steps backward or too long, jumps of no definite charge
    strong_sys = System(3, Fermion(conserve = strong(N)))
    ps = as_trajectory(random_state(Pure, strong_sys, ["Occ", "Emp", "Occ"], 4; rng))
    @test_throws "by different amounts" evolve(-im * hf + Dissipator(C)(2), 0.1, ps; nsteps = 1)
    @test real(trace(evolve(-im * hf + Dissipator(N)(2), 0.1, ps; nsteps = 2, rng))) ≈ 1
    @test_throws "negative rate, -0.1" evolve(-im * h - 0.1 * Dissipator(Sm)(1), 0.1, t;
                                             nsteps = 1)
    @test_throws "gives their jumps a negative rate" evolve(
        [-im * h, Dissipator(Sm)(1) => s -> -1], 0.1, t; nsteps = 1)
    @test_throws "does not keep the trace" evolve(Evolver(X)(1) + Dissipator(Sm)(1), 0.1, t;
                                                  nsteps = 1)
    @test_throws "does not keep the trace" evolve([-im * h + X(1), -X(1) => f,
                                                   Dissipator(Sm)(1)], 0.1, t; nsteps = 1)
    @test_throws "forward in real time" evolve(g, -0.1, t; nsteps = 1)
    @test_throws "above one: take smaller steps" evolve(-im * h + 100 * Dissipator(Sm)(1), 0.1,
                                                        t; nsteps = 1)
    @test_throws "no definite flux" evolve(-im * hf + Dissipator(C + dag(C))(2), 0.1, tf;
                                           nsteps = 1)
    # in a run, the generator of the simulation draws them, for each trajectory in turn
    run(name; on = [:state]) = runTMS(SimData(; name, seed = 5, samples = 1, phases = [
                    CreateState{Trajectory}(qubits, "Up"), CreateState{Trajectory}(qubits, "Up";
                                                                                   on = :b),
                    Evolve(evolver = g, duration = 0.5, nsteps = 10, on = on)],
                    final_measurements = EachSample(Data("z")) => Z); output = devnull).data["z"]
    @test run("a") == run("b")
    @test !isempty(run("c"; on = [:state, :b]))
end

@testset "Trajectory draws its initial state from a mixed description, whose average it is" begin
    rng = Xoshiro(13)
    draws(f) = [ f() for _ in 1:2000 ]
    same(f, ρ) = same_density(draws(f), ρ; floor = 1 / 2000)
    # local mixed states, by their eigenvectors, a mixture weighted by the traces of its states
    qubits = System(3, Qubit())
    locals = ["FullyMixed", [0.7 0.2; 0.2 0.3], "Up"]
    @test same(() -> State{Trajectory}(qubits, locals; rng), State{Mixed}(qubits, locals))
    m = Mixture([0.3 => ["Up", [1, 1], "Dn"], 0.7 => [[1 0; 0 1], "Dn", [2, 0]]])
    @test same(() -> State{Trajectory}(qubits, m; rng), State{Mixed}(qubits, m))
    @test same(() -> State{Trajectory}(qubits, FullyMixed(); rng),
               State{Mixed}(qubits, FullyMixed()))
    # fermions conserving their number weakly: each draw of one number, a sector drawn uniformly
    fermions = System(4, Fermion(conserve = N))
    fl = ["FullyMixed", [0.6 0; 0 0.4], "Occ", [0.5 0; 0 0.5]]
    ts = draws(() -> State{Trajectory}(fermions, fl; rng))
    @test same_density(ts, State{Mixed}(fermions, fl); floor = 1 / 2000)
    @test all(t -> isinteger(round(real(expect(t, sum(N(i) for i in 1:4))); digits = 10)), ts)
    @test same(() -> State{Trajectory}(fermions, FullyMixed(N => 2); rng),
               State{Mixed}(fermions, FullyMixed(N => 2)))
    # random states: the ancillas of the purification of random_state(Mixed) measured
    pq = TMS.purification(Trajectory, qubits, 4; rng = Xoshiro(3), eltype = ComplexF64)
    @test same(() -> TMS.projected(pq, qubits, rng), random_state(Mixed, qubits, 4;
                                                                  rng = Xoshiro(3)))
    states = ["Occ", "Emp", "Occ", "Emp"]
    pf = TMS.purification(Trajectory, fermions, states, 4; rng = Xoshiro(4), eltype = ComplexF64)
    @test same(() -> TMS.projected(pf, fermions, rng),
               random_state(Mixed, fermions, states, 4; rng = Xoshiro(4)))
    @test random_state(Trajectory, fermions, states, 4; rng) isa State{Trajectory}
    # refused: values conserved strongly spread over, coherences between parities, matrices that
    # are no density matrices
    strong_sys = System(2, Fermion(conserve = strong(N)))
    @test_throws "conserving it strongly forbids" State{Trajectory}(strong_sys, "FullyMixed")
    @test_throws "the terms have different values of N" State{Trajectory}(
        strong_sys, Mixture([0.5 => ["Occ", "Emp"], 0.5 => "Emp"]))
    @test State{Trajectory}(strong_sys, FullyMixed(N => 1); rng) isa State{Trajectory}
    # the rounding left in a local state does not decide its strongly conserved values
    ρ = [1e-18 0; 0 1]
    @test expect(State{Trajectory}(strong_sys, Mixture([0.5 => [ρ, ρ], 0.5 => "Occ"]); rng),
                 N(1)) ≈ 1
    @test_throws "cannot draw a random mixed state" random_state(Trajectory, strong_sys,
                                                                 ["Occ", "Emp"], 4)
    @test_throws "no definite charge" State{Trajectory}(System(2, Fermion(conserve = parity(N))),
                                                        [0.5 0.5; 0.5 0.5])
    @test_throws "is not Hermitian" State{Trajectory}(qubits, [0.5 0.5; 0 0.5])
    @test_throws "is not positive" State{Trajectory}(qubits, [1.5 0; 0 -0.5])
    @test_throws "drawn from the purification" random_state(Trajectory, qubits, FullyMixed(), 4)
    # a system conserving nothing imposes no parity on the fermions: a local state, a reset target
    # and Kraus operators superposing parities, the latter on the second site, under a string
    free = System(2, Fermion())
    coherent = [0.5 0.5; 0.5 0.5]
    @test same(() -> State{Trajectory}(free, coherent; rng), State{Mixed}(free, coherent))
    pfree = random_state(Pure, free, 4; rng = Xoshiro(5))
    for site in 1:2
        reset = Reset([1, 1])(site)
        @test same(() -> apply(reset, as_trajectory(pfree); rng), apply(reset, mix(pfree)))
    end
    # jumps of mixed parity, under a string
    lindblad = -im * (dag(C)(1) * C(2) + dag(C)(2) * C(1)) + Dissipator(Id + C)(2)
    jumped = [ evolve(lindblad, 0.2, as_trajectory(pfree); nsteps = 20, rng) for _ in 1:300 ]
    @test same_density(jumped, evolve(lindblad, 0.2, mix(pfree); nsteps = 20); floor = 1 / 300)
    u = (Id + im * (C + dag(C))) / sqrt(2)
    mixing = (0.5 * Gate(u) + 0.5 * Gate(Id))(2)
    @test same(() -> apply(mixing, as_trajectory(pfree); rng), apply(mixing, mix(pfree)))
    # the rounding between parities, or in a diagonal block, chooses neither the blocks nor the
    # eigenvectors drawn, nor the Kraus operators of a reset
    σe = [0.25 0 0 0.05; 0 0.3 0 0; 0 0 0.2 0; 0.05 0 0 0.25]
    σe[1, 2] = σe[2, 1] = 1e-14
    electron = System(1, Electron())
    @test all(_ -> abs(real(expect(State{Trajectory}(electron, σe; rng), TMS.F(1)))) ≈ 1, 1:50)
    @test all(_ -> abs(real(expect(State{Trajectory}(System(1, Qubit()),
                                                     [0.5 1e-14; 1e-14 0.5]; rng), Z(1)))) ≈ 1,
              1:20)
    parities = [0, 1, 1, 0]
    definite(m) = allequal(mod(parities[a] + parities[b], 2) for a in 1:4, b in 1:4
                           if abs(m[a, b]) > 1e-12)
    @test all(k -> definite(matrix(k, Electron())), TMS.reset_kraus(Reset(σe), Electron()))
    # in a run, the generator of the simulation draws them
    run(name; kwargs...) = runTMS(SimData(; name, seed = 5, samples = 1, phases = [
                    CreateState{Trajectory}(fermions, kwargs[:state];
                                            randomize = get(kwargs, :randomize, 0))],
                    final_measurements = EachSample(Data("n")) => N); output = devnull).data["n"]
    @test run("a"; state = "FullyMixed") == run("b"; state = "FullyMixed")
    @test run("c"; state = states, randomize = 4) == run("d"; state = states, randomize = 4)
    @test_throws "drawn from the purification" run("e"; state = FullyMixed(N => 2),
                                                   randomize = 4)
end

@testset "Trajectories run as the samples of a simulation, from the phase that makes them" begin
    qubits = System(3, Qubit())
    h = sum(Z(i) * Z(i + 1) for i in 1:2) + sum(X(i) for i in 1:3)
    g = -im * h + sum(0.5 * Dissipator(Sm)(i) for i in 1:3)
    run(name, phases; samples = 3) = runTMS(SimData(; name, seed = 9, samples, phases,
                    final_measurements = EachSample(Data("z")) => Z); output = devnull).data["z"]
    # as Realizations returning the phases from the first making a trajectory
    body = [CreateState{Trajectory}(qubits, ["Up", "FullyMixed", "Dn"]),
            Evolve(evolver = g, duration = 0.2, nsteps = 2)]
    @test run("a", body) == run("b", [Realizations(_ -> body)])
    prefix = CreateState{Pure}(qubits, "Up")
    tail = [ConvertState{Trajectory}(), Evolve(evolver = g, duration = 0.2, nsteps = 2)]
    @test run("c", [prefix; tail]) == run("d", [prefix, Realizations(_ -> tail)])
    rows = run("e", body)
    @test rows[1]["Z"] ≠ rows[2]["Z"]
    # their average that of the mixed state, for the same disorder
    disordered(R) = c -> begin
        J = 1 + 0.3 * randn(c.rng)
        [prefix, ConvertState{R}(), Evolve(evolver = -im * J * h + sum(0.5 * Dissipator(Sm)(i)
                                                                         for i in 1:3),
                                           duration = 0.5, nsteps = 20)]
    end
    average(R) = only(runTMS(SimData(; name = "m", seed = 4, samples = 80,
                                     phases = [Realizations(disordered(R))],
                                     final_measurements = Data("z") => Z); output = devnull).data["z"])
    t, m = average(Trajectory), average(Mixed)
    @test all(abs(t["Z($i)"] - m["Z($i)"]) < 4 * t["err(Z($i))"] for i in 1:3)
    # refused: trajectories without samples, or before Realizations, samples without either
    @test_throws "making trajectories needs samples" SimData(phases = body)
    @test_throws "comes before Realizations" SimData(samples = 2, phases = [
                     body[1], Realizations(_ -> body[2:2])])
    @test_throws "or makes trajectories" SimData(samples = 2, phases = [prefix])
    # a trajectory outside the samples refused when it appears, whatever makes it
    made = Update(c -> c.state = State{Trajectory}(c.state))
    evolved = Evolve(evolver = g, duration = 0.2, nsteps = 2)
    @test_throws "a trajectory outside the samples" runTMS(SimData(phases = [prefix, made]);
                                                           output = devnull)
    @test_throws "a trajectory outside the samples" runTMS(SimData(samples = 2, phases = [
                     prefix, made, Realizations(_ -> [evolved])]); output = devnull)
    @test length(run("f", [Realizations(_ -> [prefix, made, evolved])])) == 3
end

@testset "A trajectory weighs the averages by its trace, carried at each step" begin
    q = System(1, Qubit())
    K = Operator{1}("K", [1.0 0; 0 0.5])
    run(lists; samples = 3, more = []) = runTMS(SimData(; seed = 3, samples, phases = [
                    CreateState{Trajectory}(q, "FullyMixed"), Gates(gates = Gate(K)(1)), more...],
                    final_measurements = lists); output = devnull)
    rows = run(EachSample(Data("z")) => Z(1)).data["z"]
    @test all(r -> r["log_weight"] ≈ (r["Z(1)"] > 0 ? 0 : log(0.25)), rows)
    w = [ exp(r["log_weight"]) for r in rows ]
    z = [ r["Z(1)"] for r in rows ]
    avg = only(run(Data("z") => Z(1)).data["z"])
    @test avg["Z(1)"] ≈ sum(w .* z) / sum(w)
    @test avg["N_eff"] ≈ sum(w)^2 / sum(w .^ 2)
    @test avg["log_weight"] ≈ log(sum(w) / 3)
    # against the mixed state, its expectation value and its trace
    many = only(run(Data("z") => Z(1); samples = 400).data["z"])
    ρ = apply(Gate(K)(1), State{Mixed}(q, "FullyMixed"))
    @test abs(many["Z(1)"] - real(expect(ρ, Z(1)))) < 4 * many["err(Z(1))"]
    @test abs(many["log_weight"] - log(real(trace(ρ)))) < 4 * many["err(log_weight)"]
    # carried at each step: a thousand quarterings of the trace in one phase
    quarter = Update(c -> c.state = 0.5 * c.state; nsteps = 1000)
    rows = run(EachSample(Data("z")) => Z(1); more = [quarter]).data["z"]
    @test all(r -> r["log_weight"] ≈ 1000 * log(0.25) + (r["Z(1)"] > 0 ? 0 : log(0.25)), rows)
    # the column reserved, the value read
    @test_throws "one the averages or the samples write" run(Data("z") => ["log_weight" => Z(1)])
    @test only(run(Data("z") => ["lw" => :log_weight]).data["z"])["lw"] isa Real
    # outside a simulation, a trace near the bounds of the numbers refused
    @test_throws "near the bounds of the numbers" apply(Gate(1e-130 * Id)(1),
                                                         State{Trajectory}(q, "Up"))
    # weights hundreds of orders of magnitude apart, against sums in BigFloat
    lws = [0., 800., -800., 300., 800.5]
    xs = [0.3, -1.2, 5.0, 2.2, 0.7]
    cells(lw, x) = Pair{String, Any}["x" => x, "log_weight" => lw]
    sums = TMS.RowSum(0., cells(lws[1], xs[1]))
    for (lw, x) in zip(lws, xs)
        sums = TMS.added(sums, cells(lw, x), 0., (1, 1, 1), "d")
    end
    out = Dict(TMS.averaged_cells(sums))
    bw = exp.(big.(lws))
    @test out["x"] ≈ Float64(sum(bw .* xs) / sum(bw)) rtol = 1e-12
    @test out["log_weight"] ≈ Float64(log(sum(bw) / 5)) rtol = 1e-12
    @test out["N_eff"] ≈ Float64(sum(bw)^2 / sum(bw .^ 2)) rtol = 1e-12
end

@testset "A result imposed weighs a trajectory by its probability, a weight of zero stopping it" begin
    q = System(1, Qubit())
    # outside a simulation: the trace kept by a result drawn, times the probability of one imposed
    t = State{Trajectory}(q, [1, 1])
    @test real(trace(collapse(t, Z(1); rng = Xoshiro(1)).state)) ≈ 2
    r = collapse(t, Z(1) => 1)
    @test r.probability ≈ 0.5 && real(trace(r.state)) ≈ 1
    r = collapse(State{Trajectory}(q, "Up"), Z(1) => -1)
    @test r.probability == 0 && real(trace(r.state)) == 0
    @test_throws "of probability" collapse(State{Pure}(q, "Up"), Z(1) => -1)
    # in a simulation: the trajectories of probability zero stopped, counted in the mean weight
    ρ0 = [0.7 0; 0 0.3]
    post = Update(c -> c.state = collapse(c.state, Z(1) => 1).state)
    evolved = Evolve(evolver = -im * X(1), duration = 0.3, nsteps = 3,
                     measurements = Data("steps") => Z(1))
    sim = runTMS(SimData(seed = 4, samples = 400,
                         phases = [CreateState{Trajectory}(q, ρ0), post, evolved],
                         final_measurements = Data("z") => Z(1)); output = devnull)
    z = only(sim.data["z"])
    @test abs(z["log_weight"] - log(0.7)) < 4 * z["err(log_weight)"]
    @test abs(z["N"] - 280) < 4 * sqrt(400 * 0.7 * 0.3)
    @test all(r -> r["N"] == z["N"], sim.data["steps"])
    ρ = evolve(-im * X(1), 0.3, collapse(State{Mixed}(q, ρ0), Z(1) => 1).state; nsteps = 3)
    @test z["Z(1)"] ≈ real(expect(ρ, Z(1))) atol = 1e-8
end

@testset "Trajectories thermalize by their weights, whose average is the thermal state" begin
    n = 6
    sys = System(n, Qubit())
    h = -sum(Z(i) * Z(i + 1) for i in 1:n-1) - 0.7 * sum(X(i) for i in 1:n)
    # outside a simulation: normalized, the logarithm of its squared norm returned
    lt, t = thermalize(h, 0.5, State{Trajectory}(sys, "FullyMixed"; rng = Xoshiro(1)); nsteps = 4)
    @test real(trace(t)) ≈ 1 && isfinite(lt)
    # in a simulation, against the mixed state, the weights spreading as beta grows
    run(beta) = only(runTMS(SimData(seed = 2, samples = 100, phases = [
                    CreateState{Trajectory}(sys, "FullyMixed"),
                    Thermalize(; hamiltonian = h, beta, nsteps = 10)],
                    final_measurements = Data("e") => ["E" => h]); output = devnull).data["e"])
    neff = map((0.5, 1.5)) do beta
        r = run(beta)
        l, ρ = thermalize(h, beta, State{Mixed}(sys, "FullyMixed"); nsteps = 10)
        @test abs(r["E"] - real(expect(ρ, h))) < 4 * r["err(E)"]
        @test abs(r["log_weight"] - l) < 4 * r["err(log_weight)"]
        r["N_eff"]
    end
    @test neff[2] < neff[1]
    # refused: a trajectory other than the state, whose trace weighs the sample
    @test_throws "acts on :state" runTMS(SimData(samples = 1, phases = [
                     CreateState{Trajectory}(sys, "FullyMixed"),
                     CreateState{Trajectory}(sys, "FullyMixed"; on = :b),
                     Thermalize(hamiltonian = h, beta = 0.5, nsteps = 2, on = :b)]); output = devnull)
end
