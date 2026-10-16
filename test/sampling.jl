# The projective measurements: probabilities, sample and collapse, of a site and of an operator
# placed on sites, on a result drawn or given, their fast paths for a product of involutions and
# for a sum of operators of one site checked against the dense one, and the measurement whose
# result is not read, the average of the states its results leave.

"""
the results of measuring `op` on `st` with their probabilities, as pairs, through the reduced
density matrix on the sites of `op` and the projectors on the eigenspaces of its matrix there:
the reference the fast paths of probabilities are checked against
"""
function dense_probabilities(st::State, op)
    sites, values, ps = TMS.measured_spectrum(st.system, op, "dense")
    ρ = reduced_density_matrix(st, sites)
    return [ λ => real(LA.tr(ρ * p)) for (λ, p) in zip(values, ps) ]
end

@testset "Qubit sampling" begin
    # `sample` takes the generator to draw from, and these tests pass one of their own
    # rather than leaning on the global one: the frequencies below then come out the same
    # whether the whole suite runs or only this group, so the tolerances are met or missed
    # once and for all instead of now and then
    rng = Xoshiro(20260920)
    sys = System(3, Qubit())
    st = State{Pure}(sys, ["Up", "Dn", "+"])  # site 3 is a 50/50 superposition
    # deterministic sites always give the same outcome
    for _ in 1:20
        @test sample(st, 1; rng) == 0
        @test sample(st, 2; rng) == 1
    end
    # superposed site: frequency should be close to 1/2
    n = 2000
    s3 = [sample(st, 3; rng) for _ in 1:n]
    @test all(x -> x in (0, 1), s3)
    @test isapprox(sum(s3) / n, 0.5; atol = 0.05)
    # sampling the whole state is consistent with per-site sampling
    for _ in 1:20
        r = sample(st; rng)
        @test r[1] == 0 && r[2] == 1 && r[3] in (0, 1)
    end
    # classical diagonal mixture with known probabilities
    p0 = 0.2
    stm = State{Mixed}(System(1, Qubit()), [p0 0. ; 0. 1 - p0])
    nm = 4000
    samples = [sample(stm, 1; rng) for _ in 1:nm]
    @test all(x -> x in (0, 1), samples)
    @test isapprox(count(==(0), samples) / nm, p0; atol = 0.03)
    samples_full = [sample(stm; rng)[1] for _ in 1:nm]
    @test isapprox(count(==(0), samples_full) / nm, p0; atol = 0.03)
    # the vector carried along a mixed chain held the probability of the outcomes drawn so far,
    # which fell below the smallest float after some 1074 qubits: every later site came out in
    # its last state. A fully mixed chain draws ones and zeros to its end
    far = sample(State{Mixed}(System(1200, Qubit()), "FullyMixed"); rng)
    @test 0 < sum(far[1101:1200]) < 100
end

@testset "Sampling a site of a charged state" begin
    # the site legs of a charged state have a direction, which the projector on a basis state
    # has to take
    sq = State{Pure}(System(3, Fermion(conserve = N)), ["Occ", "Emp", "Occ"])
    @test sample(sq, 1) == 1
    @test sample(sq, 2) == 0
end

@testset "Measuring a site" begin
    # the result is drawn with its probability, and the state left is the state projected onto
    # it and normalized, which the outcome on one half of a Bell pair fixes on the other.
    # collapse returns the triple (result, probability, state)
    rng = Xoshiro(20261007)
    plus = State{Pure}(System(2, Qubit()), "+")
    ps = probabilities(plus, 1)
    @test first.(ps) == [0, 1]
    @test last.(ps) ≈ [0.5, 0.5]
    ps = probabilities(plus, X(1))
    @test first.(ps) ≈ [-1, 1]
    @test last.(ps) ≈ [0, 1] atol = 1e-12
    @test sample(plus, X(1); rng) ≈ 1
    @test collapse(plus, 1; rng).result isa Int
    bell = apply(controlled(X)(1, 2) * H(1), State{Pure}(System(2, Qubit()), "Up"))
    for st in (bell, mix(bell)), _ in 1:10
        x, p, s = collapse(st, 1; rng)
        @test p ≈ 0.5
        @test real(expect(s, Z(2))) ≈ (x == 0 ? 1 : -1)
        @test real(trace(s)) ≈ 1
    end
    # the projector on a level is Proj(Level => x), a number denoting no state
    r = random_state(Pure, System(4, Qubit()), 4)
    for st in (r, mix(r)), _ in 1:4
        x, _, s = collapse(st, 2; rng)
        @test norm(s - normalize(apply(Proj(Level => x)(2), st))) < 1e-12
        λ, _, s = collapse(st, X(3); rng)
        @test norm(s - normalize(apply(Proj(X => λ)(3), st))) < 1e-12
    end
    up = State{Mixed}(System(1, Qubit()), "Up")
    @test count(≈(1), [ sample(up, X(1); rng) for _ in 1:2000 ]) / 2000 ≈ 0.5 atol = 0.05
    # Ntot does not tell Up from Dn, whose coherence the measurement keeps
    e = State{Pure}(System(1, Electron()), [[0., 1., 1., 0.] / sqrt(2)])
    for st in (e, mix(e))
        λ, p, s = collapse(st, Ntot(1); rng)
        @test λ ≈ 1
        @test p ≈ 1
        @test abs(expect(s, (dag(Cup) * Cdn)(1))) ≈ 0.5
    end
    # on a strong symmetry
    fs = mix(apply(sqrt(Swap)(1, 2),
                   State{Pure}(System(2, Qubit(conserve = strong(N))), ["Dn", "Up"])))
    x, _, s = collapse(fs, 1; rng)
    @test real(expect(s, N(1))) ≈ x
    @test real(expect(s, N(1) + N(2))) ≈ 1
end

@testset "Measuring an operator of several sites" begin
    # the results and their probabilities from the eigenspaces of the matrix of the operator on
    # its sites, and the projection by its projectors placed there; collapse returns the triple
    # (result, probability, state)
    rng = Xoshiro(20261008)
    # a parity measured on |++> leaves a Bell pair, and Swap tells the singlet from the triplets
    plus = State{Pure}(System(2, Qubit()), "+")
    ps = probabilities(plus, Z(1) * Z(2))
    @test first.(ps) ≈ [-1, 1]
    @test last.(ps) ≈ [0.5, 0.5]
    for st in (plus, mix(plus))
        s, _, after = collapse(st, Z(1) * Z(2); rng)
        @test real(expect(after, Z(1) * Z(2))) ≈ s
        @test real(expect(after, X(1) * X(2))) ≈ 1
    end
    ud = State{Pure}(System(2, Qubit()), ["Up", "Dn"])
    @test last.(probabilities(ud, Swap(1, 2))) ≈ [0.5, 0.5]
    # the number of excitations of three sites, from the diagonal of their density matrix
    r = random_state(Pure, System(5, Qubit()), 4)
    ρ = reduced_density_matrix(r, [1, 3, 5])
    counts = zeros(4)
    for (i, b) in enumerate(Iterators.product(0:1, 0:1, 0:1))
        counts[sum(b) + 1] += real(ρ[i, i])
    end
    ps = probabilities(r, N(1) + N(3) + N(5))
    @test first.(ps) ≈ 0:3
    @test last.(ps) ≈ counts
    # a hop over site 2, of spectrum -1, 0, 1: its projectors are polynomials of it, which
    # make_mpo places with their Jordan-Wigner strings, and so must the projection
    f = random_state(State{Pure}(System(4, Fermion(conserve = N)), ["Occ", "Occ", "Emp", "Emp"]), 4)
    h = dag(C)(1) * C(3) + dag(C)(3) * C(1)
    projs = Dict(-1 => (h * h - h) / 2, 0 => Id(1) - h * h, 1 => (h * h + h) / 2)
    for (λ, p) in probabilities(f, h)
        @test p ≈ real(expect(f, projs[round(Int, λ)])) atol = 1e-12
    end
    for st in (f, mix(f)), _ in 1:6
        λ, _, s = collapse(st, h; rng)
        ref = normalize(apply(make_mpo(f, projs[round(Int, λ)]), f))
        @test norm(s - (st isa State{Mixed} ? mix(ref) : ref)) < 1e-10
    end
    # an operator moving a charge is measured on the state conserving nothing, but the state it
    # leaves has no definite charges, which a conserving system cannot hold
    q = State{Pure}(System(3, Qubit(conserve = N)), ["Up", "Dn", "Up"])
    @test collapse(q, Z(1) * Z(2); rng).result ≈ -1
    @test last.(probabilities(q, X(1))) ≈ [0.5, 0.5]
    @test last.(probabilities(mix(q), X(1))) ≈ [0.5, 0.5]
    @test_throws "collapse weaken(state, ()) instead" collapse(q, X(1); rng)
    @test_throws "collapse weaken(state, ()) instead" collapse(mix(q), X(1) => 1)
    @test collapse(weaken(q, ()), X(1) => 1).probability ≈ 0.5
    # an operator that is not Hermitian is refused, one not commuting with the fermionic parity
    # measured, see "Measuring an operator of no definite fermionic parity"
    @test_throws "needs a Hermitian operator" probabilities(plus, Sp(1) * Z(2))
    # the dimension 2^64, which an Int wrapped around to 0, letting the matrix be built
    up64, xx = State{Pure}(System(64, Qubit()), "Up"), sum(X(i) * X(i + 1) for i in 1:63)
    @test_throws "of dimension 18446744073709551616, beyond the 4096" probabilities(up64, xx)
    @test_throws "needs an operator placed on sites" probabilities(plus, Id(1))
    @test_throws "which the system does not have" sample(plus, X(7))
    @test_throws "which the state does not have" collapse(plus, 3)
    # sample of a site checks it, which 1.x did not
    @test_throws "sample was given site 3, which the state does not have" sample(plus, 3)
    @test_throws "sample was given site 0" sample(mix(plus), 0)
end

@testset "Measuring an operator of no definite fermionic parity" begin
    # the projectors of the dense path are placed as operators, whose strings the descent lays:
    # on a state superposing the two parities the results of a Majorana operator, ±1, have the
    # probabilities its mean value gives, and the state projected on one of them has it as mean
    # value, in the pure and in the mixed representation
    φ = random_state(Pure, System(4, Fermion()), 4)
    γ3 = (C + dag(C))(3)
    @test abs(expect(φ, γ3)) > 1e-3
    for st in (φ, mix(φ))
        ps = probabilities(st, γ3)
        @test first.(ps) ≈ [-1, 1]
        @test sum(first.(ps) .* last.(ps)) ≈ real(expect(φ, γ3))
        λ, p, after = collapse(st, γ3 => 1)
        @test λ ≈ 1
        @test p ≈ last(last(ps))
        @test real(expect(after, γ3)) ≈ 1
        @test real(trace(after)) ≈ 1
    end
    # the product of two Majorana operators of two sites is anti-Hermitian
    @test_throws "needs a Hermitian operator" probabilities(φ, (C + dag(C))(1) * (C + dag(C))(2))
end

@testset "Collapse on a given result" begin
    # the probability of the result and the state projected onto it, on each way of measuring:
    # the result post-selected. collapse returns the triple (result, probability, state)
    plus = State{Pure}(System(2, Qubit()), "+")
    for st in (plus, mix(plus))
        s, p, after = collapse(st, Z(1) * Z(2) => -1)
        @test s ≈ -1
        @test p ≈ 0.5
        @test real(expect(after, Z(1) * Z(2))) ≈ -1
        @test real(expect(after, X(1) * X(2))) ≈ 1
    end
    # the projector of a count keeps a real state real, where it made it complex
    plus4 = State{Pure}(System(4, Qubit()), "+")
    for st in (plus4, mix(plus4))
        r = collapse(st, sum(Z(i) for i in 1:4) => 2)
        @test r.probability ≈ 1/4
        @test eltype(r.state.state[1]) == Float64
    end
    ψ = random_state(Pure, System(6, Qubit()), 6)
    P = Z(1) * X(3) * Y(5)
    op = N(1) + N(3) + N(5)
    _, values, projectors = TMS.measured_spectrum(ψ.system, op, "dense")
    for st in (ψ, mix(ψ))
        for (s, q) in probabilities(st, P)
            _, p, after = collapse(st, P => s)
            ref = normalize(apply(make_mpo(ψ, (Id(1) + s * P) / 2), ψ))
            @test p ≈ q
            @test norm(after - (st isa State{Mixed} ? mix(ref) : ref)) < 1e-12
        end
        for (v, q) in probabilities(st, op)
            _, p, after = collapse(st, op => v)
            proj = Operator{3}("P", projectors[findfirst(≈(v), values)])(1, 3, 5)
            @test p ≈ q
            @test norm(after - normalize(apply(proj, st))) < 1e-12
        end
    end
    _, p, singlet = collapse(State{Pure}(System(2, Qubit()), ["Up", "Dn"]), Swap(1, 2) => -1)
    @test p ≈ 0.5
    @test real(expect(singlet, Swap(1, 2))) ≈ -1
    bell = apply(controlled(X)(1, 2) * H(1), State{Pure}(System(2, Qubit()), "Up"))
    for st in (bell, mix(bell))
        x, p, after = collapse(st, 1 => 1)
        @test x === 1
        @test p ≈ 0.5
        @test real(expect(after, Z(2))) ≈ -1
    end
    # the level given as a real
    x, p, after = collapse(bell, 1 => 1.0)
    @test x === 1
    @test p ≈ 0.5
    @test real(expect(after, Z(2))) ≈ -1
    up = State{Pure}(System(2, Qubit()), "Up")
    @test_throws "0.5 is not a result of measuring X(1)" collapse(plus, X(1) => 0.5)
    @test_throws "of probability" collapse(up, Z(1) => -1)
    @test_throws "is not a result" collapse(up, 1 => 2)
    @test_throws "is not a result" collapse(up, 1 => 0.5)
    @test_throws "which the state does not have" collapse(up, 3 => 0)
end

@testset "Dephasing as the average of the collapses" begin
    # the measurement whose result is not read is the average of the states its results leave,
    # with the strings of the fermionic sites in between and the order of the sites; collapse
    # returns the triple (result, probability, state)
    averaged(st, op) = sum(c.probability * c.state for c in (collapse(st, op => λ)
                                                            for (λ, q) in probabilities(st, op)
                                                            if q > 1e-14))
    ρ = mix(random_state(Pure, System(4, Qubit()), 4))
    @test norm(apply(Dephase(Z ⊗ Z)(1, 3), ρ) - averaged(ρ, Z(1) * Z(3))) < 1e-12
    @test norm(apply(Dephase(Swap)(2, 4), ρ) - averaged(ρ, Swap(2, 4))) < 1e-12
    f = mix(random_state(State{Pure}(System(4, Fermion(conserve = N)),
                                     ["Occ", "Occ", "Emp", "Emp"]), 4))
    h = dag(C)(1) * C(3) + dag(C)(3) * C(1)
    j = im * (dag(C)(1) * C(3) - dag(C)(3) * C(1))
    hop = dag(C) ⊗ C - C ⊗ dag(C)
    current = im * (dag(C) ⊗ C + C ⊗ dag(C))
    for (D, op) in [(Dephase(hop)(1, 3), h), (Dephase(hop)(3, 1), h),
                    (Dephase(current)(1, 3), j), (Dephase(current)(3, 1), j)]
        @test norm(apply(D, f) - averaged(f, op)) < 1e-12
    end
end

@testset "Measuring a product of involutions" begin
    # measured by its mean value on any number of sites, it agrees with the matrix of the
    # operator on a few, and its projection, the state plus its image, with (1 + sP)/2;
    # collapse returns the triple (result, probability, state)
    rng = Xoshiro(20261009)
    ψ = random_state(Pure, System(6, Qubit()), 6)
    for op in (Z(1) * X(3) * Y(5), -2 * Z(1) * Z(2), X(2) * X(3) * X(4) * X(5)), st in (ψ, mix(ψ))
        @test !isnothing(TMS.involution(st.system, op))
        ps, ds = probabilities(st, op), dense_probabilities(st, op)
        @test first.(ps) ≈ first.(ds)
        @test last.(ps) ≈ last.(ds)
    end
    P = Z(1) * X(3) * Y(5)
    projected(λ) = normalize(apply(make_mpo(ψ, (Id(1) + λ * P) / 2), ψ))
    for _ in 1:4
        λ, _, after = collapse(ψ, P; rng)
        @test norm(after - projected(λ)) < 1e-12
        λ, _, after = collapse(mix(ψ), P; rng)
        @test norm(after - mix(projected(λ))) < 1e-12
    end
    p = last(last(probabilities(ψ, P)))
    @test count(_ -> collapse(ψ, P; rng).result > 0, 1:1000) / 1000 ≈ p atol = 0.05
    # on twenty sites, the bond dimension at most doubles
    big = random_state(Pure, System(20, Qubit()), 8)
    string20 = prod(X(i) for i in 1:2:20) * prod(Z(i) for i in 2:2:20)
    s, _, after = collapse(big, string20; rng)
    @test real(expect(after, string20)) ≈ s
    @test maxlinkdim(after) ≤ 2 * maxlinkdim(big)
    # a diagonal one on a conserving system, a parity of fermions, and a strong symmetry
    q = random_state(State{Pure}(System(6, Qubit(conserve = N)),
                                 ["Up", "Dn", "Up", "Dn", "Up", "Dn"]), 4)
    Pz = Z(1) * Z(4) * Z(6)
    @test last.(probabilities(q, Pz)) ≈ last.(dense_probabilities(q, Pz))
    @test isnothing(TMS.involution(q.system, X(1) * X(2)))
    s, _, after = collapse(q, Pz; rng)
    @test real(expect(after, Pz)) ≈ s
    f = random_state(State{Pure}(System(5, Fermion(conserve = N)),
                                 ["Occ", "Emp", "Occ", "Emp", "Occ"]), 4)
    Pf = F(1) * F(3) * F(4)
    @test last.(probabilities(f, Pf)) ≈ last.(dense_probabilities(f, Pf))
    for st in (f, mix(f))
        s, _, after = collapse(st, Pf; rng)
        @test real(expect(after, Pf)) ≈ s
        @test real(trace(after)) ≈ 1
    end
    fs = mix(apply(sqrt(Swap)(1, 2),
                   State{Pure}(System(3, Qubit(conserve = strong(N))), ["Dn", "Up", "Up"])))
    s, _, after = collapse(fs, Z(1) * Z(2); rng)
    @test real(trace(after)) ≈ 1
    @test real(expect(after, N(1) + N(2) + N(3))) ≈ 1
    # the factors must be even, so that they commute: i times two Majorana operators, an
    # involution whose factors anticommute, goes to the dense path
    φ = random_state(Pure, System(4, Fermion()), 4)
    m = im * (C + dag(C))(1) * (C + dag(C))(3)
    @test isnothing(TMS.involution(φ.system, (C + dag(C))(1) * (C + dag(C))(3)))
    for st in (φ, mix(φ))
        ps = probabilities(st, m)
        @test first.(ps) ≈ [-1, 1]
        @test last.(ps) ≈ last.(dense_probabilities(st, m))
    end
end

@testset "Measuring a sum of operators of one site" begin
    # its totals, integers or half integers, are counted on any number of sites, in agreement
    # with the matrix of the operator on a few, and its projection by an MPO carrying the
    # partial sums with the projector on its eigenspace; collapse returns the triple
    # (result, probability, state)
    rng = Xoshiro(20261010)
    ψ = random_state(Pure, System(6, Qubit()), 6)
    s = random_state(Pure, System([Qubit(), Spin(1), Qubit(), Spin(3/2)]), 6)
    for (st, op) in [(ψ, N(1) + N(3) + N(5)), (ψ, Sz(1) + Sz(2) + Sz(4)), (ψ, X(1) + X(2) + X(3)),
                     (ψ, 2 * (N(2) - N(6))), (s, Sz(1) + Sz(2) + Sz(3) + Sz(4))], x in (st, mix(st))
        @test !isnothing(TMS.counting(x.system, op, "test"))
        ps, ds = probabilities(x, op), dense_probabilities(x, op)
        @test first.(ps) ≈ first.(ds)
        @test last.(ps) ≈ last.(ds) atol = 1e-12
    end
    op = N(1) + N(3) + N(5)
    _, values, projectors = TMS.measured_spectrum(ψ.system, op, "dense")
    projected(st, v) =
        normalize(apply(Operator{3}("P", projectors[findfirst(≈(v), values)])(1, 3, 5), st))
    for st in (ψ, mix(ψ)), _ in 1:3
        v, _, after = collapse(st, op; rng)
        @test norm(after - projected(st, v)) < 1e-12
    end
    # on twenty sites, the number of particles of the region becomes sharp
    big = random_state(Pure, System(30, Qubit()), 8)
    NA = sum(N(i) for i in 6:25)
    ps = probabilities(big, NA)
    @test sum(last.(ps)) ≈ 1
    @test sum(first.(ps) .* last.(ps)) ≈ real(expect(big, NA))
    v, _, after = collapse(big, NA; rng)
    @test real(expect(after, NA)) ≈ v
    @test real(expect(after, NA * NA)) ≈ v^2
    # on charged fermions, the total stays, and on a strong symmetry
    f = random_state(State{Pure}(System(8, Fermion(conserve = N)),
                                 ["Occ", "Emp", "Occ", "Emp", "Occ", "Emp", "Occ", "Emp"]), 6)
    Nf = sum(N(i) for i in 2:6)
    @test last(last(probabilities(f, Nf))) ≈ 0 atol = 1e-12
    for st in (f, mix(f))
        v, _, after = collapse(st, Nf; rng)
        @test real(expect(after, Nf)) ≈ v
        @test real(expect(after, sum(N(i) for i in 1:8))) ≈ 4
        @test real(trace(after)) ≈ 1
    end
    fs = mix(apply(sqrt(Swap)(1, 2),
                   State{Pure}(System(3, Qubit(conserve = strong(N))), ["Dn", "Up", "Dn"])))
    v, _, after = collapse(fs, N(1) + N(3); rng)
    @test real(expect(after, N(1) + N(3))) ≈ v
    @test real(trace(after)) ≈ 1
    # the parts must be even, so that they commute: the sum of two Majorana operators, which
    # anticommute, has the results ±√2 and not -2, 0, 2, and goes to the dense path
    φ = random_state(Pure, System(4, Fermion()), 4)
    γγ = (C + dag(C))(1) + (C + dag(C))(3)
    @test isnothing(TMS.counting(φ.system, γγ, "test"))
    for st in (φ, mix(φ))
        ps = probabilities(st, γγ)
        @test first.(ps) ≈ [-√2, √2]
        @test last.(ps) ≈ last.(dense_probabilities(st, γγ))
    end
    # an impossible total has a probability of zero, which the sum of Fourier left at ±1e-14,
    # and collapse refuses it, where it projected on nothing
    f = random_state(State{Pure}(System(4, Fermion(conserve = N)), ["Occ", "Occ", "Emp", "Emp"]), 4)
    total = N(1) + N(2) + N(3) + N(4)
    for st in (f, mix(f))
        @test all(((t, p),) -> t == 2 ? p ≈ 1 : p == 0, probabilities(st, total))
        @test_throws "of probability 0.0" collapse(st, total => 3)
    end
end
