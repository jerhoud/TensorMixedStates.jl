# The representation Trajectory, a state held as one quantum trajectory: its types, the algebra it
# shares with Pure, its files, and what it refuses

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

@testset "Trajectory measures the rest on its normalised ket, and refuses what needs draws" begin
    rng = Xoshiro(6)
    qubits = System(4, Qubit())
    p = random_state(Pure, qubits, 4; rng)
    q = random_state(Pure, qubits, 4; rng)
    # a weight of 9, which these measures do not see: those of the pure state p
    t = 3 * as_trajectory(p)
    for f in (purity, renyi2, s -> renyi2(s, [1, 2]), s -> mutual_info_renyi2(s, [1]), hermiticity,
              s -> entanglement_entropy(s, 2), s -> entanglement_spectrum(s, 2),
              s -> von_neumann_entropy(s, [1]), s -> log_negativity(s, [1], [2]),
              s -> inner(s, q), s -> inner(q, s), s -> inner(s, X(1), q),
              s -> fidelity(s, mix(q)), s -> variance(s, X(1)), Vector,
              s -> last.(measure(s, Measure([Norm, Overlap(q), Variance(X(2))]); time = 0.)))
        @test f(t) ≈ f(p) rtol = 1e-10
    end
    @test inner(t, t) ≈ 1
    for ((ct, a), (cp, b)) in zip(entanglement_by_sector(t, 2), entanglement_by_sector(p, 2))
        @test ct == cp && a.weight ≈ b.weight && a.entropy ≈ b.entropy && a.spectrum ≈ b.spectrum
    end
    h = sum(Z(i) * Z(i + 1) for i in 1:3) + sum(X(i) for i in 1:4)
    TMS = TensorMixedStates
    for (f, message) in (
            (() -> mix(t), "a trajectory is not a mixed state"),
            (() -> partial_trace(t; keep = [1]), "a trajectory is not a mixed state"),
            (() -> t ⊗ t, "the product of two trajectories"),
            (() -> t ⊗ p, "the product of a trajectory"),
            (() -> random_state(t, 2), "randomising a trajectory is not implemented"),
            (() -> ground_state(h, t), "a trajectory is not searched"),
            (() -> steady_state(-im * h, t), "a steady state comes with the draws"),
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

@testset "Trajectory leaves the dissipators and the channels to evolve and apply" begin
    rng = Xoshiro(8)
    qubits = System(3, Qubit())
    p = random_state(Pure, qubits, 4; rng)
    t = as_trajectory(p)
    h = sum(Z(i) * Z(i + 1) for i in 1:2) + sum(X(i) for i in 1:3)
    @test_throws "given the generator rather than its PreMPO" PreMPO(t, -im * h +
                                                                      Dissipator(X)(1))
    @test_throws "which a trajectory draws through apply" evolve(
        depolarizing_channel(0.2)(1), 0.1, t; nsteps = 1)
    steady(R) = runTMS(SimData(; name = "s", samples = 1, phases = [
                       CreateState{R}(qubits, "Up"),
                       SteadyState(lindbladian = -im * h, algo = Dmrg(nsweeps = 1))]);
                       output = devnull)
    @test_throws "a steady state comes with the draws" steady(Trajectory)
    # beside an exact state
    e = ExactState(p)
    q = random_state(Pure, qubits, 4; rng)
    @test fidelity(t, ExactState(q)) ≈ fidelity(t, q) rtol = 1e-12
    @test fidelity(ExactState(q), t) ≈ fidelity(q, t) rtol = 1e-12
    @test inner(t, e) ≈ inner(TensorMixedStates.as_pure(t), e) rtol = 1e-12
    @test fidelity(t, mix(e)) ≈ fidelity(TensorMixedStates.as_pure(t), mix(e)) rtol = 1e-12
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
    # refused: Kraus operators of no definite charge; a channel not keeping the trace weighs
    @test_throws "as written have no definite charge" apply(
        0.5 * (Gate(C + dag(C)) + Gate(im * (C - dag(C))))(2), as_trajectory(pf); rng)
    @test real(trace(apply(0.5 * depolarizing_channel(0.2)(1), t; rng))) ≈ 0.5
    @test real(trace(apply(Gate(Proj("Up"))(1) + Gate(Proj("Up"))(1), t; rng))) ≈
          2 * real(expect(p, Proj("Up")(1)))
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
    pick(ps) = something(findfirst(>(rand(rng) * sum(ps)), cumsum(ps)), findlast(>(0), ps))
    value(A, φ) = real(φ' * A * φ) / real(φ' * φ)
    rates(φ, t) = [ c(f, t) * value(sum(L' * L for L in Ls), φ) for (f, Ls) in decays ]
    stepped(φ, t, s) = exp(s * sum(c(f, t + s / 2) * K for (f, K) in ks)) * φ
    # the rest `τ` of a step from `t`: no jump up to s with the probability exp(-a0 s - b s²)
    function part(φ, t, τ)
        a0 = max(sum(rates(φ, t)), 0.0)
        whole = stepped(φ, t, τ)
        b = (max(sum(rates(whole, t + τ)), 0.0) - a0) / (2τ)
        lu = log(rand(rng))
        if lu < -τ * (a0 + b * τ)
            return whole
        end
        s = clamp(-2lu / (a0 + sqrt(max(a0^2 - 4b * lu, 0.0))), 0.0, τ)
        φ = stepped(φ, t, s)
        r = rates(φ, t + s)
        if sum(r) > 0
            Ls = last(decays[pick(r)])
            φ = Ls[pick([ value(L' * L, φ) for L in Ls ])] * φ
        end
        return s < τ ? part(φ, t + s, τ - s) : φ
    end
    n = LA.norm(ψ)
    for k in 0:nsteps - 1
        ψ = part(ψ, k * dt, dt)
        ψ *= n / LA.norm(ψ)
    end
    return ψ
end

# whether the trajectories evolved by `evolve` from `t` with the seeds 1:n are those of
# `dense_trajectory`
function same_trajectories(evolver, t, ks, decays; duration = 0.5, nsteps = 10, n = 10,
                           algo = Tdvp(), rtol = 1e-8)
    ψ = raw_vector(State{Pure}(t.system, t.state))
    return all(1:n) do seed
        a = raw_vector(State{Pure}(t.system, evolve(evolver, duration, t; nsteps, algo,
                                                    rng = Xoshiro(seed)).state))
        b = dense_trajectory(ψ, ks, decays, duration, nsteps, Xoshiro(seed))
        isapprox(a, b; rtol)
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
    # the average
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
    # ApproxW, whose steps before and after a jump take the place of the step in its cache
    @test same_trajectories(g, t, [(nothing, -im * m(h) + effective(jumps[1])),
                                   (f, effective(jumps[2]))],
                            [(nothing, jumps[1]), (f, jumps[2])]; algo = ApproxW(order = 4),
                            rtol = 1e-4)
    # steps holding several jumps, five on average
    xs = [sqrt(50) * m(X(1))]
    @test same_trajectories(-im * h + 50 * Dissipator(X)(1), t,
                            [(nothing, -im * m(h) + effective(xs))], [(nothing, xs)];
                            duration = 0.2, nsteps = 2)
    # a rate vanishing at the end of the evolution, which the end of the last step passes by
    # rounding
    @test real(trace(evolve([-im * h, Dissipator(Sm)(1) => s -> 0.3 - s], 0.3, t; nsteps = 10,
                            rng))) ≈ 1
    # a step so long that the norm without jump underflows: its jumps drawn, or zero when the
    # operators lose it
    q1 = System(1, Qubit())
    @test real(trace(evolve(-im * X(1) + 760 * Dissipator(Z)(1), 1.0,
                            State{Trajectory}(q1, "Up"); nsteps = 1, rng))) ≈ 1
    @test real(trace(evolve(-1000 * Proj("Dn")(1) + Dissipator(Sm)(1), 1.0,
                            State{Trajectory}(q1, "Dn"); nsteps = 1, rng))) == 0
    # refused: a jump changing a strong charge, as Mixed refuses it, negative rates, steps
    # backward, jumps of no definite charge
    strong_sys = System(3, Fermion(conserve = strong(N)))
    ps = as_trajectory(random_state(Pure, strong_sys, ["Occ", "Emp", "Occ"], 4; rng))
    @test_throws "by different amounts" evolve(-im * hf + Dissipator(C)(2), 0.1, ps; nsteps = 1)
    @test real(trace(evolve(-im * hf + Dissipator(N)(2), 0.1, ps; nsteps = 2, rng))) ≈ 1
    @test_throws "negative rate, -0.1" evolve(-im * h - 0.1 * Dissipator(Sm)(1), 0.1, t;
                                             nsteps = 1)
    @test_throws "gives their jumps a negative rate" evolve(
        [-im * h, Dissipator(Sm)(1) => s -> -1], 0.1, t; nsteps = 1)
    # a rate function negative within a step or at its end, refused whatever the steps, and one
    # zero at the end of the evolution, which rounding puts past it, 0.2 + 0.1 > 0.3
    for n in 1:3
        @test_throws "gives their jumps a negative rate" evolve(
            [-im * h, Dissipator(Sm)(1) => s -> 1 - 2s], 0.9, t; nsteps = n, rng)
    end
    @test real(trace(evolve([-im * h, Dissipator(Sm)(1) => s -> 0.3 - s], 0.3, t; nsteps = 3,
                            rng))) ≈ real(trace(t))
    @test_throws "forward in real time" evolve(g, -0.1, t; nsteps = 1)
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

@testset "Dephasing draws the projectors on the eigenspaces, or the phases between them" begin
    rng = Xoshiro(43)
    # the Kraus operators of the phases on a qutrit: unitaries over √3, the dephasing on average
    q3 = System(1, Qudit(3))
    ks = [ lowered(q3, k) for k in TMS.kraus_operators(q3, Dephase(Level; kraus = :phases)(1)) ]
    @test length(ks) == 3
    @test all(k -> 3 * k' * k ≈ LA.I, ks)
    ρ = (x = randn(rng, ComplexF64, 3, 3); x * x')
    @test sum(k * ρ * k' for k in ks) ≈ LA.Diagonal(LA.diag(ρ))
    # on qubits, the jump is Z, the average of the trajectories the evolution
    qubits = System(2, Qubit())
    p = random_state(Pure, qubits, 2; rng)
    t = as_trajectory(p)
    g = -im * (X(1) * X(2) + Z(1)) +
        sum(dephasing_dissipator(0.8; jumps = :phases)(i) for i in 1:2)
    ts = [ evolve(g, 0.5, t; nsteps = 5, rng) for _ in 1:200 ]
    @test same_density(ts, evolve(g, 0.5, mix(p); nsteps = 40); floor = 1 / 200)
    # the density matrix does not see the choice
    @test ket_bra_matrix(apply(dephasing_channel(0.3; kraus = :phases)(1), mix(p))) ≈
          ket_bra_matrix(apply(dephasing_channel(0.3)(1), mix(p)))
    # in the eigenbasis of a sum of operators of one site, products of phases of each site, those
    # of the dense matrix up to a phase each, which allow many sites
    tensor_sum(n) = sum(reduce(⊗, [ j == i ? N : Id for j in 1:n ]) for i in 1:n)
    q3s = System(3, Qubit())
    sum3 = Dephase(tensor_sum(3); kraus = :phases)
    products = [ lowered(q3s, k) for k in TMS.kraus_operators(q3s, sum3(1, 2, 3)) ]
    dense = [ lowered(q3s, k(1, 2, 3)) for k in TMS.generic_kraus(sum3, [ q3s[i] for i in 1:3 ]) ]
    @test length(products) == length(dense) == 4
    @test all(abs(LA.tr(a' * b)) ≈ LA.norm(a)^2 ≈ LA.norm(b)^2 for (a, b) in zip(products, dense))
    q14 = System(14, Qubit())
    g14 = dephasing_dissipator(0.5, tensor_sum(14); jumps = :phases)(1:14...)
    @test length(TMS.kraus_operators(q14, Dephase(tensor_sum(14); kraus = :phases)(1:14...))) == 15
    t14 = evolve(g14, 0.2, State{Trajectory}(q14, "+"); nsteps = 2, rng)
    @test real(trace(t14)) ≈ 1
    @test real(expect(t14, sum(N(i) for i in 1:14))) ≈ 7
    # the relaxation of a qubit passes its jumps to its dephasing: Z at the rate γ / 2
    q1 = System(1, Qubit())
    γ = 1 / 3. - 1 / 10.
    for (jumps, dephasing) in ((:projectors, [ sqrt(γ) * [1 0; 0 0], sqrt(γ) * [0 0; 0 1] ]),
                               (:phases, [ sqrt(γ / 2) * [1 0; 0 -1] ]))
        e = TMS.prepared_evolver(as_trajectory(State{Pure}(q1, "Up")),
                                 thermal_relaxation_dissipator(5., 3.; jumps)(1))
        js = [ lowered(q1, j) for d in e.decays for j in d.jumps ]
        @test length(js) == 1 + length(dephasing)
        @test any(j -> j ≈ sqrt(1 / 5.) * [0 1; 0 0], js)
        @test all(k -> any(j -> j ≈ k || j ≈ -k, js), dephasing)
    end
    @test occursin("kraus = :phases", repr(Dephase(Z; kraus = :phases)))
    @test_throws "the jumps of dephasing_dissipator are :projectors or :phases" dephasing_dissipator(
        0.1; jumps = :pauli)
    @test_throws "the Kraus operators of Dephase are :projectors or :phases" Dephase(Z; kraus = :x)
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
    # the jumps counted, a Poisson number of mean 2 for jumps X at the rate 2 for a time 1
    q1 = System(1, Qubit())
    counts = [ r["jumps"] for r in runTMS(SimData(; name = "j", seed = 3, samples = 200, phases = [
                   CreateState{Trajectory}(q1, "Up"),
                   Evolve(evolver = -im * Z(1) + 2 * Dissipator(X)(1), duration = 1.0, nsteps = 4)],
                   final_measurements = EachSample(Data("j")) => ["jumps" => :jumps]);
                   output = devnull).data["j"] ]
    @test all(c -> isinteger(c) && c ≥ 0, counts)
    @test abs(sum(counts) / 200 - 2) < 5 * sqrt(2 / 200)
    @test_throws "a name the run reserves" SimData(samples = 2, phases = [
                     CreateState{Trajectory}(q1, "Up"), CreateState{Pure}(q1, "Up"; on = :jumps)])
    # a new trajectory counts its own jumps, and carries its own weight
    each(phases, column) = [ r[column] for r in runTMS(SimData(; name = "r", seed = 5, samples = 20,
                                 phases, final_measurements = EachSample(Data("r")) =>
                                 ["jumps" => :jumps, "z" => Z(1)]); output = devnull).data["r"] ]
    renewed = [CreateState{Trajectory}(q1, "Up"),
               Evolve(evolver = 2 * Dissipator(X)(1), duration = 1.0, nsteps = 4),
               CreateState{Trajectory}(q1, "Up")]
    @test all(==(0), each(renewed, "jumps"))
    thermal = [CreateState{Trajectory}(q1, "FullyMixed"),
               Thermalize(hamiltonian = Z(1), beta = 1.0, nsteps = 2),
               CreateState{Trajectory}(q1, "Dn")]
    @test all(==(0), each(thermal, "log_weight"))
    # the run alone sets them
    @test_throws "jumps cannot be changed" each([CreateState{Trajectory}(q1, "Up"),
                                                 Update(c -> c.jumps = 7)], "jumps")
    @test_throws "log_weight cannot be changed" each([CreateState{Trajectory}(q1, "Up"),
                                                      Update(c -> c.log_weight = 3.)], "jumps")
    # names of their own in a simulation with samples making no trajectory
    pure = runTMS(SimData(; name = "p", samples = 2, phases = [Realizations(_ -> [
                      CreateState{Pure}(q1, "Up"), Update(c -> c.jumps = 1; provides = [:jumps])])],
                      final_measurements = Data("p") => ["jumps" => :jumps]); output = devnull)
    @test only(pure.data["p"])["jumps"] == 1
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

@testset "Channels and evolvers not keeping the trace weigh a trajectory by what they keep" begin
    rng = Xoshiro(21)
    qubits = System(2, Qubit())
    p = random_state(Pure, qubits, 2; rng)
    t = as_trajectory(p)
    # the average of trace × ⟨Z(1)⟩ and of the trace, against the channel applied to the mixed state
    function weighed(f, ρ; n = 400)
        values = [ (s = f(); [real(trace(s)), real(trace(s)) * real(expect(s, Z(1)))])
                   for _ in 1:n ]
        return within_error(values, [real(trace(ρ)), real(trace(ρ)) * real(expect(ρ, Z(1)))])
    end
    K = Operator{1}("K", [1.0 0; 0 0.5])
    for ch in ((0.5 * Gate(X) + 0.3 * Gate(Id))(1),
               Gate(K)(1) * depolarizing_channel(0.3)(2),
               0.7 * amplitude_damping_channel(0.4)(2))
        @test weighed(() -> apply(ch, t; rng), apply(ch, mix(p)))
    end
    # an evolver whose operators lose trace, with jumps, in long steps
    g = -im * (X(1) * X(2) + Z(1)) - 0.4 * Proj("Dn")(1) + Dissipator(sqrt(0.5) * Sm)(2)
    @test weighed(() -> evolve(g, 0.5, t; nsteps = 2, rng), evolve(g, 0.5, mix(p); nsteps = 40);
                  n = 300)
    # with a rate depending on time
    gt = [-im * (X(1) * X(2) + Z(1)) - 0.4 * Proj("Dn")(1),
          Dissipator(sqrt(0.5) * Sm)(2) => s -> 1 + 2s]
    @test weighed(() -> evolve(gt, 0.5, t; nsteps = 2, rng), evolve(gt, 0.5, mix(p); nsteps = 40);
                  n = 300)
end

@testset "Weights stay right at zero, under rounding, across rows and for one carrier" begin
    rng = Xoshiro(31)
    q = System(1, Qubit())
    q2 = System(2, Qubit())
    K = Operator{1}("K", [1.0 0; 0 0.1])
    # a part or a factor of probability zero leaves a null trajectory, not NaN nor an error
    t = State{Trajectory}(q2, ["Dn", "Up"])
    @test real(trace(apply(Gate(Proj("Up"))(1) + Gate(Proj("Up"))(1), t; rng))) == 0
    @test real(trace(apply(depolarizing_channel(0.1)(2) * Gate(Proj("Up"))(1), t; rng))) == 0
    # a probability zero but for rounding is no negative weight
    θ = 0.7
    Q = Operator{1}("Q", [sin(θ)^2 -sin(θ)*cos(θ); -sin(θ)*cos(θ) cos(θ)^2])
    @test real(trace(apply(Gate(Q)(1) + Gate(Q)(1), State{Trajectory}(q, [cos(θ), sin(θ)]);
                           rng))) < 1e-25
    # in a simulation, those trajectories stop, in a product or after two collapses in one step
    run(phases, lists; samples = 40) = runTMS(SimData(; seed = 7, samples, phases,
                                                      final_measurements = lists); output = devnull)
    r = only(run([CreateState{Trajectory}(q, "FullyMixed"),
                  Gates(gates = depolarizing_channel(0.1)(1) * Gate(Proj("Up"))(1))],
                 Data("z") => Z(1)).data["z"])
    @test abs(r["log_weight"] - log(0.5)) < 4 * r["err(log_weight)"]
    twice = Update(c -> (c.state = collapse(c.state, Z(1) => 1).state;
                         c.state = collapse(c.state, Z(2) => 1).state))
    r = only(run([CreateState{Trajectory}(q2, "FullyMixed"), twice], Data("z") => Z(1)).data["z"])
    @test abs(r["log_weight"] - log(0.25)) < 4 * r["err(log_weight)"]
    # the error of a weighted mean, that of a ratio of means
    rows = run([CreateState{Trajectory}(q, "FullyMixed"), Gates(gates = Gate(K)(1))],
               EachSample(Data("z")) => Z(1)).data["z"]
    avg = only(run([CreateState{Trajectory}(q, "FullyMixed"), Gates(gates = Gate(K)(1))],
                   Data("z") => Z(1)).data["z"])
    w = [ exp(r["log_weight"]) for r in rows ]
    x = [ r["Z(1)"] for r in rows ]
    m = sum(w .* x) / sum(w)
    n = length(w)
    @test avg["err(Z(1))"] ≈ sqrt(sum(w .^ 2 .* (x .- m) .^ 2) * n / (n - 1)) / sum(w)
    # a tolerance waits for the error of the mean weight too
    sim = runTMS(SimData(seed = 7, samples = 30, min_samples = 10, tol = Data("v") => 1e-6,
                         phases = [CreateState{Trajectory}(q, "FullyMixed"),
                                   Gates(gates = Gate(K)(1))],
                         final_measurements = Data("v") => ["one" => Proj("Up")(1) +
                                                                      Proj("Dn")(1)]);
                 output = devnull)
    @test only(sim.data["v"])["N"] == 30
    # every row of a simulation of trajectories weighs, 0 before the trajectory, in one text file
    in_tmp() do
        runTMS(SimData(name = "e", seed = 7, samples = 2, phases = [
                    CreateState{Pure}(q, "Up"; final_measurements = EachSample("e.txt") => Z(1)),
                    ConvertState{Trajectory}(final_measurements = EachSample("e.txt") => Z(1))]))
        lines = readlines("e/e.txt")
        @test occursin("log_weight", lines[1]) && allequal(length.(split.(lines, '\t')))
    end
    # :log_weight only in a simulation making trajectories; :log_trace refused on trajectories
    @test_throws "only a simulation making trajectories has" run(
        [CreateState{Pure}(q, "Up"), Realizations(c -> [Gates(gates = Gate(X)(1))])],
        Data("v") => ["lw" => :log_weight]; samples = 2)
    @test_throws "reads log_trace" SimData(samples = 2, phases = [
                     CreateState{Trajectory}(q, "FullyMixed"),
                     Thermalize(hamiltonian = Z(1), beta = 0.5, nsteps = 2,
                                measurements = Data("t") => :log_trace)])
    # only :state carries a weight
    @test_throws "which only :state carries" run([CreateState{Trajectory}(q, "FullyMixed"),
                                                   CreateState{Trajectory}(q, "FullyMixed";
                                                                           on = :b),
                                                   Gates(gates = Gate(K)(1), on = :b)],
                                                  Data("z") => Z(1); samples = 2)
    # truncation is no weight: an evolution or gates keeping the trace bring it back
    six = System(6, Qubit())
    h6 = sum(X(i) * X(i + 1) + Y(i) * Y(i + 1) for i in 1:5) + 0.5 * sum(Z(i) for i in 1:6)
    t6 = State{Trajectory}(six, "FullyMixed"; rng)
    @test real(trace(evolve(-im * h6, 1., t6; nsteps = 10, limits = Limits(maxdim = 2)))) ≈ 1
    circuit = prod(Gate(exp(-0.3im * (X ⊗ X + Y ⊗ Y)))(i, i + 1) for i in 1:5)
    @test real(trace(apply(circuit, evolve(-im * h6, 1., t6; nsteps = 4);
                           limits = Limits(maxdim = 2), rng))) ≈ 1
    r = run([CreateState{Trajectory}(six, "FullyMixed"),
             CreateState{Trajectory}(six, "FullyMixed"; on = :b),
             Evolve(evolver = -im * h6, duration = 1., nsteps = 10, limits = Limits(maxdim = 2),
                    on = [:state, :b])], Data("z") => Z(3); samples = 2).data["z"]
    @test only(r)["log_weight"] ≈ 0 atol = 1e-12
    # a local density matrix weighs by its trace, as a vector by its norm; negative ones refused
    @test real(trace(State{Trajectory}(q, [2 0.3; 0.3 1]; rng))) ≈ 3
    @test real(trace(State{Trajectory}(q, [0 0; 0 0]))) == 0
    @test_throws "so it is not a density matrix" State{Trajectory}(q, [-1 0; 0 0])
    # an operator keeping the trace, alone or in a product of channels, brings it back after a
    # truncation, its phase kept
    four = System(4, Qubit())
    u = exp(-0.7im * (X ⊗ X + Y ⊗ Y))
    t4 = State{Trajectory}(random_state(Pure, four, 4; rng))
    one_bond = Limits(maxdim = 1)
    @test real(trace(apply(u(2, 3), t4; limits = one_bond))) ≈ 1
    @test real(trace(apply(depolarizing_channel(0.1)(1) * Gate(u)(2, 3), t4; limits = one_bond,
                           rng))) ≈ 1
    # one that does not keep it multiplies it by ⟨A†A⟩, computed before, truncation or not; a
    # truncation keeps it
    A = X(2) * X(3) + 0.5 * Z(2)
    @test real(trace(apply(A, t4; limits = one_bond))) ≈ real(expect(t4, dag(A) * A))
    @test real(trace(truncate(t4; limits = one_bond))) ≈ 1
    pure(s) = raw_vector(State{Pure}(s.system, s.state))
    @test pure(apply(im * X(1), t4)) ≈ im * pure(apply(X(1), t4))
    # an Update too, though it provides the state
    halved = Update(c -> c.b = collapse(c.b, Z(1) => 1).state; provides = [:b])
    @test_throws "which only :state carries" run([CreateState{Trajectory}(q, "FullyMixed"),
                                                   CreateState{Trajectory}(q, [1, 1] / sqrt(2);
                                                                           on = :b),
                                                   halved], Data("z") => Z(1); samples = 2)
end

@testset "A run of trajectories notes the measurements it takes on each normalised ket" begin
    in_tmp() do
        q = System(3, Qubit())
        list = "m" => [Z, EntanglementEntropy(2)]
        runTMS(SimData(name = "t", seed = 1, samples = 1, final_measurements = list,
                       phases = [CreateState{Trajectory}(q, "Up")]))
        line = only(filter(l -> occursin("normalised ket", l), readlines("t/log")))
        @test endswith(line, "density matrix: EntanglementEntropy(2)")
        runTMS(SimData(name = "m", final_measurements = list,
                       phases = [CreateState{Mixed}(q, "Up")]))
        @test !occursin("normalised ket", read("m/log", String))
    end
end

@testset "ExactState holds a trajectory, drawn and measured as State draws and measures it" begin
    rng = Xoshiro(21)
    TMS = TensorMixedStates
    qubits = System(3, Qubit())
    p = random_state(Pure, qubits, 4; rng)
    t = 3 * as_trajectory(p)
    e = ExactState(t)
    pure(s) = raw_vector(State{Pure}(s.system, s.state))
    # the trace, its weight, kept both ways and in files
    @test e isa ExactState{Trajectory} && e.state ≈ pure(t) && trace(e) ≈ trace(t)
    back = State(e)
    @test back isa State{Trajectory} && pure(back) ≈ pure(t)
    mktempdir() do dir
        file = joinpath(dir, "e.h5")
        save_state(file, "e", e)
        r = load_state(file, "e")
        @test r isa ExactState{Trajectory} && r.state == e.state
    end
    @test TMS.as_representation(Trajectory, ExactState(p)) isa ExactState{Trajectory}
    @test ExactState{Trajectory}(qubits, "Up") isa ExactState{Trajectory}
    # linear in the density matrix, and on the normalised ket
    q = random_state(Pure, qubits, 4; rng)
    op = X(1) * Z(2) + 0.3 * Y(3)
    @test expect(e, op) ≈ expect(t, op) rtol = 1e-12
    @test expect1(e, Z) ≈ expect1(t, Z) rtol = 1e-12
    @test reduced_density_matrix(e, [1, 3]) ≈ reduced_density_matrix(t, [1, 3]) rtol = 1e-12
    @test von_neumann_entropy(e, [2]) ≈ von_neumann_entropy(t, [2]) rtol = 1e-10
    @test entanglement_entropy(e, 1) ≈ entanglement_entropy(t, 1) rtol = 1e-10
    @test log_negativity(e, [1], [2]) ≈ log_negativity(t, [1], [2]) rtol = 1e-10
    @test purity(e) == 1 && hermiticity(e) == 1
    @test variance(e, op) ≈ variance(t, op) rtol = 1e-10
    @test Vector(e) ≈ Vector(t) rtol = 1e-12
    @test fidelity(e, ExactState(q)) ≈ fidelity(t, q) rtol = 1e-12
    @test fidelity(e, mix(ExactState(q))) ≈ fidelity(t, mix(q)) rtol = 1e-12
    @test fidelity(e, ExactState(3 * as_trajectory(q))) ≈ fidelity(t, q) rtol = 1e-12
    @test inner(e, ExactState(q)) ≈ inner(t, q) rtol = 1e-12
    @test TMS.measured_norm(e) ≈ 1
    # applied and collapsed as State applies and collapses it, with the same draws
    dep = depolarizing_channel(0.6)(1) * amplitude_damping_channel(0.4)(2)
    @test all(1:10) do seed
        apply(dep, e; rng = Xoshiro(seed)).state ≈ pure(apply(dep, t; rng = Xoshiro(seed)))
    end
    A = X(2) * X(3) + 0.5 * Z(2)
    @test apply(A, e).state ≈ pure(apply(A, t))
    c = collapse(e, Z(1) => 1)
    @test c.probability ≈ collapse(t, Z(1) => 1).probability && trace(c.state) ≈ 9c.probability
    # refused as on a State
    h = sum(X(i) for i in 1:3) + sum(Z(i) * Z(i + 1) for i in 1:2)
    for (f, message) in (
            (() -> mix(e), "a trajectory is not a mixed state"),
            (() -> partial_trace(e; keep = [1]), "a trajectory is not a mixed state"),
            (() -> e ⊗ e, "the product of two trajectories"),
            (() -> e ⊗ ExactState(p), "the product of a trajectory"),
            (() -> random_state(e, 2), "randomising a trajectory is not implemented"),
            (() -> ground_state(h, e), "a trajectory is not searched"),
            (() -> steady_state(-im * h, e), "a steady state comes with the draws"),
            (() -> apply(Dissipator(X)(1), e), "is a generator"),
            (() -> TMS.as_representation(Pure, e), "a trajectory is not a pure or a mixed"),
            (() -> TMS.as_representation(Trajectory, mix(ExactState(p))),
             "a trajectory of a mixed state"))
        @test_throws message f()
    end
    # drawn at random from the purification, of a definite charge, that of the purification only
    # fixed by its start
    r = random_state(ExactState{Trajectory}, qubits, 2; rng)
    @test r isa ExactState{Trajectory} && real(trace(r)) ≈ 1
    fermions = System(3, Fermion(conserve = N))
    rf = random_state(ExactState{Trajectory}, fermions, ["Occ", "Emp", "Occ"], 2; rng)
    n = sum(expect1(rf, N))
    @test abs(n - round(n)) < 1e-12 && 1 ≤ round(n) ≤ 3
    @test_throws "without a sector" random_state(ExactState{Trajectory}, fermions, 2; rng)
end

@testset "ExactState evolves a trajectory exactly, its jumps drawn as State draws them" begin
    rng = Xoshiro(22)
    TMS = TensorMixedStates
    qubits = System(3, Qubit())
    p = random_state(Pure, qubits, 4; rng)
    e = ExactState(as_trajectory(p))
    h = sum(X(i) for i in 1:3) + sum(Z(i) * Z(i + 1) for i in 1:2)
    m(o) = lowered(qubits, o)
    # without dissipators, as a pure state, the time functions over the whole duration
    for g in (-im * h, [-im * h, -im * Z(2) => s -> cos(3s)])
        @test evolve(g, 1., e; nsteps = 20).state ≈
              evolve(g, 1., ExactState(p); nsteps = 20).state rtol = 1e-11
    end
    # the jumps replayed with the same seeds by an independent evaluation
    dep = depolarizing_channel(0.5)(1)
    kdep = [ sqrt(0.4) * m(k) for k in TMS.kraus_operators(qubits, dep)
             if !(m(k) ≈ m(k)[1, 1] * LA.I) ]
    g = -im * h + 0.5 * Dissipator(Sm)(1) + 0.4 * Dissipator(dep) +
        0.3 * Dissipator(Sm(2) + Sm(3))
    jumps = [[sqrt(0.5) * m(Sm(1))] ; kdep ; [sqrt(0.3) * m(Sm(2) + Sm(3))]]
    K = -im * m(h) - 0.5 * sum(L' * L for L in jumps)
    # weighed by their exact norm, see exact_norm
    @test all(1:10) do seed
        a = evolve(g, 0.5, e; nsteps = 10, rng = Xoshiro(seed)).state
        b = dense_trajectory(e.state, [(nothing, K)], [(nothing, jumps)], 0.5, 10, Xoshiro(seed))
        isapprox(a / LA.norm(a), b / LA.norm(b); rtol = 1e-8)
    end
    # their weighted average is the evolution of the density matrix
    ts = [ evolve(g, 0.5, e; nsteps = 10, rng).state for _ in 1:400 ]
    vs = [ v * v' / real(v' * v) for v in ts ]
    mean = sum(v * v' for v in ts) / sum(v -> real(v' * v), ts)
    err = sqrt.(sum(v -> abs2.(v .- mean), vs) / (length(vs) - 1) / length(vs))
    ρ = evolve(g, 0.5, mix(ExactState(p)); nsteps = 10).state
    @test all(abs.(mean .- ρ) .≤ 5 .* max.(err, 1 / 400))
    # the weights take the error of the rate taken linear out, even in a single step where it
    # biased the unweighted average by 7 error bars (measured)
    qubit = System(1, Qubit())
    g1 = -im * X(1) + 2 * Dissipator(Sm)(1)
    e1 = ExactState{Trajectory}(qubit, "Up")
    ts1 = [ evolve(g1, 1., e1; nsteps = 1, rng) for _ in 1:4000 ]
    w = [ real(trace(t)) for t in ts1 ]
    z = [ real(expect(t, Z(1))) for t in ts1 ]
    m = sum(w .* z) / sum(w)
    ref = real(expect(evolve(g1, 1., mix(ExactState{Pure}(qubit, "Up")); nsteps = 1), Z(1)))
    @test abs(m - ref) ≤ 5 * sqrt(sum(abs2, w .* (z .- m))) / sum(w)
    # in a run of samples, weighed by an operator not keeping the trace as a State is
    A = X(1) + 0.5 * Z(2)
    x = runTMS(SimData(; name = "x", samples = 200, seed = 3, backend = ExactState, phases = [
                       CreateState{Trajectory}(qubits, "Up"),
                       Evolve(evolver = g, duration = 0.5, time_step = 0.05), Gates(gates = A)],
                       final_measurements = Data("z") => [Z]); output = devnull).data["z"][1]
    ρ = apply(A, evolve(g, 0.5, mix(ExactState{Pure}(qubits, "Up")); nsteps = 10))
    @test all(abs.(x["Z"] .- expect1(ρ, Z)) .≤ 5 .* x["err(Z)"])
    @test x["err(log_weight)"] > 0
    @test abs(x["log_weight"] - log(real(trace(ρ)))) ≤ 5 * x["err(log_weight)"] + 1e-9
    # a second trajectory carries no weight: brought back to its trace, as on a State
    y = runTMS(SimData(; name = "y", samples = 3, seed = 3, backend = ExactState, phases = [
                       CreateState{Trajectory}(qubits, "Up"), CopyState(:state => :b),
                       Evolve(evolver = g, duration = 0.5, time_step = 0.1, on = [:state, :b])],
                       final_measurements = EachSample(Data("b")) => ["tb" => c -> real(trace(c.b))]);
               output = devnull).data["b"]
    @test length(y) == 3 && all(r -> r["tb"] ≈ 1, y)
end

@testset "ExactState thermalizes a trajectory exactly, its weight carried by a run" begin
    rng = Xoshiro(23)
    qubits = System(3, Qubit())
    h = sum(X(i) for i in 1:3) + sum(Z(i) * Z(i + 1) for i in 1:2) + 0.5 * sum(Z(i) for i in 1:3)
    p = random_state(Pure, qubits, 4; rng)
    # against the exponential of the dense matrix, from a trajectory of trace 3
    e = 3 * ExactState(as_trajectory(p))
    r = thermalize(h, 0.6, e; nsteps = 3)
    v = Vector(p)
    w = LA.exp(-0.3 * Matrix(matrix(h, qubits))) * (v / LA.norm(v))
    @test r.state isa ExactState{Trajectory} && real(trace(r.state)) ≈ 1
    @test r.state.state ≈ w / LA.norm(w) rtol = 1e-9
    @test r.log_trace ≈ log(LA.norm(w)^2) rtol = 1e-9
    # step by step, as the phase does
    l, st = 0., e
    for k in 1:3
        lk, st = thermalize(h, 0.2, st; nsteps = 1, step = k, total_beta = 0.6)
        l += lk
    end
    @test st.state ≈ r.state.state rtol = 1e-9
    @test l ≈ r.log_trace rtol = 1e-9
    # the evolver of the Hamiltonian as the Hamiltonian, a dissipator refused as on a State
    @test thermalize(Evolver(h), 0.6, e; nsteps = 3).state.state ≈ r.state.state rtol = 1e-12
    @test_throws "whose jumps a trajectory draws" thermalize(Dissipator(X)(1), 0.6, e; nsteps = 3)
    # in a run, a single trajectory weighed by 3X(1), of trace 9, before the phase
    one = runTMS(SimData(; name = "one", samples = 1, seed = 4, backend = ExactState, phases = [
                         CreateState{Trajectory}(qubits, "Up"), Gates(gates = 3 * X(1)),
                         Thermalize(hamiltonian = h, beta = 1., nsteps = 4)],
                         final_measurements = Data("w") => [Z]); output = devnull).data["w"][1]
    dn = thermalize(h, 1., ExactState{Trajectory}(qubits, ["Dn", "Up", "Up"]); nsteps = 4)
    @test one["log_weight"] ≈ log(9) + dn.log_trace rtol = 1e-9
    # the weighted average of trajectories drawn from the fully mixed state is the thermal state
    x = runTMS(SimData(; name = "th", samples = 300, seed = 4, backend = ExactState, phases = [
                       CreateState{Trajectory}(qubits, "FullyMixed"),
                       Thermalize(hamiltonian = h, beta = 1., nsteps = 4)],
                       final_measurements = Data("t") => [Z]); output = devnull).data["t"][1]
    ρ = thermalize(h, 1., ExactState{Mixed}(qubits, "FullyMixed"); nsteps = 4)
    @test all(abs.(x["Z"] .- expect1(ρ.state, Z)) .≤ 5 .* x["err(Z)"])
    @test x["err(log_weight)"] > 0
    @test abs(x["log_weight"] - ρ.log_trace) ≤ 5 * x["err(log_weight)"]
end

@testset "A trajectory lists its jumps: their time, their dissipator and its Kraus index" begin
    TMS = TensorMixedStates
    # X at the constant rate 2: the first jump at -log(u) / 2, u the first draw, unless u < e^-2
    t = State{Trajectory}(System(1, Qubit()), "Up")
    @test all(1:20) do seed
        js = []
        evolve(2 * Dissipator(X)(1), 1., t; nsteps = 1, rng = Xoshiro(seed), jumps = js)
        u = rand(Xoshiro(seed))
        u < exp(-2) ? isempty(js) :
            first(js).time ≈ -log(u) / 2 && TMS.obs_name(first(js).dissipator) == "X(1)" &&
            first(js).kraus == 1
    end
    # each dissipator as written, the Kraus operators of a channel by their index, on an MPS and
    # held whole, which draw the same jumps
    q = System(3, Qubit())
    p = random_state(Pure, q, 4; rng = Xoshiro(3))
    g = -im * sum(X(i) for i in 1:3) + 0.5 * Dissipator(Sm)(1) + depolarizing_dissipator(0.4)(2) +
        dephasing_dissipator(0.3)(3)
    names = Dict("Sm(1)" => 1:1, "Reset(\"FullyMixed\")(2)" => 1:4, "Dephase(Level)(3)" => 1:2)
    for st in (State{Trajectory}(p), ExactState(State{Trajectory}(p)))
        js = []
        evolve(g, 2., st; nsteps = 20, rng = Xoshiro(5), jumps = js)
        @test !isempty(js) && issorted([ j.time for j in js ]) && 0 < js[1].time && js[end].time < 2
        @test all(j -> j.kraus in names[TMS.obs_name(j.dissipator)], js)
    end
    # in a run, written for each trajectory, as many as :jumps counts
    in_tmp() do
        sim = runTMS(SimData(name = "j", seed = 1, samples = 3, phases = [
                             CreateState{Trajectory}(q, "Up"),
                             Evolve(evolver = g, duration = 1., time_step = 0.1)],
                             final_measurements = [EachSample("jumps.jsonl") => [:jump_list],
                                                   EachSample(Data("n")) => [:jumps]]))
        rows = [ TMS.JSON.parse(l) for l in readlines("j/jumps.jsonl") ]
        counts = Dict(r["sample"] => r["jumps"] for r in sim.data["n"])
        @test length(rows) == 3 && sum(values(counts)) > 0
        @test all(r -> length(r["jump_list"]) == counts[r["sample"]], rows)
        listed = reduce(vcat, [ r["jump_list"] for r in rows ])
        @test all(j -> j["kraus"] in names[j["dissipator"]], listed)
    end
end
