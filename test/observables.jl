# Expectation values: expect, through the descent and the compaction, expect1 and expect2 with
# their fast paths, inner(a, op, b), and apply on pure states, checked against dense matrices
# written with explicit Jordan-Wigner strings. And the measurements of a state that are not
# expectation values: purity, hermiticity, entropies, reduced density matrices, negativity,
# fidelities and variance. And the projective measurements: probabilities, sample and collapse,
# their fast paths checked against the dense one. And measure: the lists of measurements, their
# names, the kinds of their values, what they refuse, the measurements of the library, the
# functions of a context with the Check they return, and the preparation a Measure keeps. The
# channels are tested in channels.jl.

"""
the matrix, on `n` sites `site`, of the product of operators of one site `factors`, pairs
`op => i` in the order of the product, written with explicit Jordan-Wigner strings: an odd
operator on site `i` carries `F` on every site before it. Independent of the descent, so it can
be used as a reference for it.
"""
function jw_matrix(site, n, factors)
    mf, id = matrix(F, site), matrix(Id, site)
    odd(a) = TMS.parity_of(matrix(a, site), mf) == 1
    placed(a, j) = foldl(kron, [ k == j ? matrix(a, site) : (k < j && odd(a) ? mf : id) for k in 1:n ])
    return prod(placed(a, j) for (a, j) in factors)
end

"""
the exact correlation matrix ``\\langle c^\\dagger_i c_j \\rangle`` of a pure state of fermions,
see `jw_matrix`
"""
function exact_fermionic_correlations(state::State{Pure}, site)
    n = length(state)
    psi = dense_vector(state)
    return [ psi' * jw_matrix(site, n, [dag(C) => i, C => j]) * psi for i in 1:n, j in 1:n ]
end

"""
the vector of the pure state `st` in the full basis, site 1 most significant, not normalised
"""
raw_vector(st::State{Pure}) = vec(Array(IT.dense(prod(st.state)), reverse(siteinds(Pure, st.system))...))

"""
the dense matrix of the MPO of `op` on the sites of the pure state `st`
"""
function mpo_matrix(st::State{Pure}, op)
    idx = siteinds(Pure, st.system)
    a = Array(IT.dense(prod(make_mpo(st, op))), [ i' for i in reverse(idx) ]...,
              [ IT.dag(i) for i in reverse(idx) ]...)
    d = prod(IT.dim.(idx))
    return reshape(a, d, d)
end

"""
the results of measuring `op` on `st` with their probabilities, as pairs, through the reduced
density matrix on the sites of `op` and the projectors on the eigenspaces of its matrix there:
the reference the fast paths of probabilities are checked against
"""
function dense_probabilities(st::State, op)
    sites, values, ps = TMS.measured_spectrum(st.system, op, "dense")
    ρ = reduced_density_matrix(st, sites)
    return [ λ => real(TMS.LinearAlgebra.tr(ρ * p)) for (λ, p) in zip(values, ps) ]
end

@testset "Fermionic correlations" begin
    # expect2 lays the Jordan-Wigner strings itself: the whole correlation matrix, at every
    # distance, is checked against an exact computation, on a random state, whose correlations
    # do not vanish at any distance
    n = 6
    fe = Fermion()
    st = random_state(Pure, System(n, fe), 4)
    ref = exact_fermionic_correlations(st, fe)
    @test ref ≈ ref'
    @test [ ref[i, i] for i in 1:n ] ≈ expect1(st, N)
    @test expect2(st, (dag(C), C)) ≈ ref atol = 1e-10
    @test expect2(st, (C, dag(C))) ≈ [ i == j ? 1 - ref[i, i] : -ref[j, i] for i in 1:n, j in 1:n ] atol = 1e-10
    # an odd operator, and a pair of mixed parities, which 1.x refused, are measured through
    # expect with their strings; on a state superposing the two parities they do not vanish
    @test expect1(st, C) ≈ [ expect(st, C(i)) for i in 1:n ]
    @test maximum(abs, expect1(st, C)) > 1e-3
    # the odd part with the string the state keeps, the even part on the site: an operator of no
    # definite parity, several operators at once, a second call on the same state, a mixed state
    @test expect1(st, C + N) ≈ [ expect(st, (C + N)(i)) for i in 1:n ]
    e = expect1(st, [C, N, dag(C)])
    @test e[3] ≈ [ expect(st, dag(C)(i)) for i in 1:n ] && e[2] ≈ expect1(st, N)
    @test expect1(st, dag(C) + 0.5N) ≈ [ expect(st, (dag(C) + 0.5N)(i)) for i in 1:n ]
    ρ = random_state(Mixed, System(n, fe), 3)
    @test expect1(ρ, C + N) ≈ [ expect(ρ, (C + N)(i)) for i in 1:n ]
    @test expect2(st, (dag(C), N)) ≈
          [ i == j ? expect(st, (dag(C) * N)(i)) : expect(st, dag(C)(i) * N(j)) for i in 1:n, j in 1:n ]
    # any other pair row by row, each part of definite parity with its string: no definite parity,
    # a mixed state, and operators defined on some of the sites only
    pairwise(x, o1, o2, on1, on2) =
        [ (i in on1 && j in on2) ? expect(x, i == j ? (o1 * o2)(i) : o1(i) * o2(j)) : missing
          for i in 1:length(x.system), j in 1:length(x.system) ]
    same(a, b) = all(ismissing(x) ? ismissing(y) : isapprox(x, y; atol = 1e-12) for (x, y) in zip(a, b))
    @test same(expect2(st, (C + N, dag(C) + N)), pairwise(st, C + N, dag(C) + N, 1:n, 1:n))
    ρ = random_state(Mixed, System(n, fe), 3)
    @test same(expect2(ρ, (N, C)), pairwise(ρ, N, C, 1:n, 1:n))
    fb = random_state(Pure, System([fe, Boson(3), fe, Boson(3), fe]), 4)
    @test same(expect2(fb, (C, N)), pairwise(fb, C, N, [1, 3, 5], 1:5))
    # pairs of both parities in one call do not disturb each other
    m = expect2(st, [(N, N), (dag(C), C), (dag(C), N)])
    @test m[1] ≈ expect2(st, (N, N)) atol = 1e-10
    @test m[2] ≈ ref atol = 1e-10
    @test m[3] ≈ expect2(st, (dag(C), N)) atol = 1e-10
    # the mixed representation, computed by another algorithm, agrees
    @test expect2(mix(st), (dag(C), C)) ≈ ref atol = 1e-8
    b = random_state(Pure, System(4, Boson(3)), 4)
    # the entries (j, i) of equal operators, and of adjoint ones on a pure state, follow from (i, j)
    pairs = [(A, dag(A)), (N, A), (dag(A), N), (A, A), (N, N)]
    @test all(expect2(b, pairs) .≈ expect2(mix(b), pairs))
    ψp = random_state(State{Pure}(System(n, Fermion(conserve = parity(N))), "Emp"), 6)
    cc = [ expect(ψp, C(i) * C(j)) for i in 1:n, j in 1:n ]
    @test maximum(abs, cc) > 0.01
    @test expect2(ψp, (C, C)) ≈ cc atol = 1e-10
    # products written in any order, a sum on a site crossed by a string, renamed operators
    psi = dense_vector(st)
    exact(factors) = psi' * jw_matrix(fe, n, factors) * psi
    stm = mix(st)
    for d in 1:n-1
        @test expect(st, dag(C)(1) * C(1 + d)) ≈ ref[1, 1 + d] atol = 1e-10
        @test expect(stm, dag(C)(1) * C(1 + d)) ≈ ref[1, 1 + d] atol = 1e-8
        @test expect(st, C(1 + d) * dag(C)(1)) ≈ -ref[1, 1 + d] atol = 1e-10
        @test expect(st, C(1 + d) * (C + dag(C))(1)) ≈ exact([C => 1 + d, (C + dag(C)) => 1]) atol = 1e-10
        @test expect(stm, C(1 + d) * (C + dag(C))(1)) ≈ exact([C => 1 + d, (C + dag(C)) => 1]) atol = 1e-8
        @test expect(st, C(1 + d) * named(dag(C), "Cd")(1)) ≈ -ref[1, 1 + d] atol = 1e-10
        @test expect(st, dag(C)(1) * named(2C, "C2")(1 + d)) ≈ 2ref[1, 1 + d] atol = 1e-10
    end
end

@testset "Electron and Tj against Jordan-Wigner strings written by hand" begin
    # two species on a site, whose order on the site the matrices of Cup and Cdn hold, and
    # strings through the whole sites on the left
    n = 4
    for site in (Electron(), Tj())
        st = random_state(Pure, System(n, site), 4)
        psi = dense_vector(st)
        exact(factors) = psi' * jw_matrix(site, n, factors) * psi
        for a in (Cup, Cdn), b in (Cup, Cdn)
            ref = [ exact([dag(a) => i, b => j]) for i in 1:n, j in 1:n ]
            @test [ expect(st, dag(a)(i) * b(j)) for i in 1:n, j in 1:n ] ≈ ref atol = 1e-10
            @test expect2(st, (dag(a), b)) ≈ ref atol = 1e-10
        end
        # a pair hopping, four fermionic factors on two sites apart
        @test expect(st, dag(Cup)(1) * dag(Cdn)(1) * Cdn(3) * Cup(3)) ≈
              exact([dag(Cup) => 1, dag(Cdn) => 1, Cdn => 3, Cup => 3]) atol = 1e-10
    end
end

@testset "Superoperators of placed operators, unnormalized expectations and matrix elements" begin
    # decision 1: Left and Right are gone, a placed operator a acts on a density matrix as
    # Gate(a), ρ ↦ aρa†, or Evolver(a), ρ ↦ aρ + ρa†. A sum with its strings and its
    # coefficients, against the pure state: tr(b aρa†) = <ψ|a† b a|ψ> and
    # tr(b (aρ + ρa†)) = <ψ|b a|ψ> + <ψ|a† b|ψ>, expect times the trace being tr(b σ)
    sys = System(4, Fermion())
    ψ = random_state(Pure, sys, 4)
    ρ = mix(ψ)
    a = dag(C)(1) * C(3) + 0.5im * N(2)
    adjoint_a = dag(C)(3) * C(1) - 0.5im * N(2)
    b = dag(C)(3) * C(1) + N(4)
    σ = apply(make_mpo(ρ, Gate(a)), ρ)
    @test trace(σ) * expect(σ, b) ≈ inner(ψ, adjoint_a * b * a, ψ)
    σ = apply(make_mpo(ρ, Evolver(a)), ρ)
    @test trace(σ) * expect(σ, b) ≈ inner(ψ, b * a, ψ) + inner(ψ, adjoint_a * b, ψ)
    # the product of the gates of one site is the gate of the product, strings included
    @test norm(apply(Gate(dag(C))(1) * Gate(C)(3), ρ) - mix(apply(dag(C)(1) * C(3), ψ))) < 1e-12

    # `normalize = false` is gone: the value not divided by the trace is expect times the trace
    @test expect(2 * ρ, 3 * Id(1)) * trace(2 * ρ) ≈ 6
    @test expect(2 * ψ, N(1)) * trace(2 * ψ) ≈ 4 * expect(ψ, N(1))
    # 1.x measured c₁ρ, of trace zero, unnormalized. Without Left, its substance is the gate of
    # c₁, tr(X c₁ρc₁†) = tr(c₁† X c₁ρ), on sites charged weakly and strongly too, where the
    # evolver c₁ρ + ρc₁† has no definite flux
    for site in (Fermion(), Fermion(conserve = N), Fermion(conserve = strong(N)))
        s = System(4, site)
        φ = mix(normalize(State{Pure}(s, ["Occ", "Emp", "Occ", "Emp"]) +
                          0.5 * State{Pure}(s, ["Emp", "Occ", "Occ", "Emp"])))
        σ = apply(Gate(C(1)), φ)
        @test trace(σ) ≈ expect(φ, N(1))
        @test trace(σ) * expect(σ, [Id(1), N(3)]) ≈ [0.8, 0.8]
        if TMS.is_charged(s)
            @test_throws "no definite flux" make_mpo(φ, Evolver(C(1)))
        end
    end

    # the matrix element of an operator between two states, and of a channel between two mixed
    # ones: <<mix(u)| Gate(X) |mix(v)>> = |<u|X|v>|^2
    q = System(4, Qubit())
    u, v = random_state(Pure, q, 3), random_state(Pure, q, 3)
    m = X(1) * Z(3) + 0.3 * Y(2)
    @test inner(u, m, v) ≈ inner(u, apply(make_mpo(v, m), v))
    @test inner(mix(u), Gate(X)(1), mix(v)) ≈ abs2(inner(u, X(1), v))
end

@testset "The adjoint of a tensor product of fermions" begin
    # (A ⊗ B)(i, j) is A(i) * B(j), whose adjoint reverses the factors: two fermions
    # anticommute, so dag(C ⊗ C) is the opposite of dag(C) ⊗ dag(C)
    st = random_state(Pure, System(3, Fermion()), 4)
    for (i, j) in [(1, 2), (1, 3), (3, 1)]
        @test expect(st, dag(C ⊗ C)(i, j)) ≈ conj(expect(st, C(i) * C(j)))
        @test expect(st, dag(dag(C) ⊗ C)(i, j)) ≈ conj(expect(st, dag(C)(i) * C(j)))
    end
    # a hopping term written with its adjoint is self adjoint, and the dissipator of a
    # fermionic jump preserves the trace
    h = dag(C) ⊗ C
    @test abs(imag(expect(st, sum((h + dag(h))(i, i + 1) for i in 1:2)))) < 1e-12
    ρ = mix(st)
    for l in (C ⊗ C, dag(C) ⊗ C, C ⊗ dag(C))
        @test abs(trace(apply(make_mpo(ρ, Dissipator(l)(1, 2)), ρ))) < 1e-12
    end
end

@testset "Functions of a fermionic operator" begin
    # exp and mod of an operator that is not even, placed after a fermionic site, take the
    # string on their odd part, which is checked against the dense matrices of the explicit
    # c₁ = c ⊗ 1 ⊗ 1, c₂ = F ⊗ c ⊗ 1 and c₃ = F ⊗ F ⊗ c
    fe = Fermion()
    st = random_state(Pure, System(3, fe), 4)
    psi = dense_vector(st)
    ev(m) = psi' * m * psi
    c, f, id = matrix(C, fe), matrix(F, fe), matrix(Id, fe)
    cs = [kron(c, id, id), kron(f, c, id), kron(f, f, c)]
    g = 0.7 * (C + dag(C))
    for i in 1:3
        b = 0.7 * (cs[i] + cs[i]')
        @test expect(st, exp(g)(i)) ≈ ev(exp(b))
        @test expect(st, mod(g, 3)(i)) ≈ ev(exp(2im * π * b / 3))
        @test expect(st, exp(C)(i)) ≈ ev(kron(id, id, id) + cs[i])
        @test expect(mix(st), exp(g)(i)) ≈ ev(exp(b))
        # a non integer power is a function of the operator too
        @test expect(st, (g^0.5)(i)) ≈ ev(b^0.5)
    end
    # a projector of no definite parity, which 1.x refused, has its odd part strung
    pr = [1. 1. ; 1. 1.] / 2
    @test expect(st, C(3) * Proj([1., 1.] / sqrt(2))(1)) ≈ ev(cs[3] * kron(pr, id, id))
    pe, po = (pr + f * pr * f) / 2, (pr - f * pr * f) / 2
    @test expect(st, Proj([1., 1.] / sqrt(2))(2) * C(3)) ≈ ev((kron(id, pe, id) + kron(f, po, id)) * cs[3])
    p = Proj([0., 1.])
    @test expect(st, C(3) * p(1)) ≈ ev(cs[3] * kron([0. 0. ; 0. 1.], id, id))
    @test expect(st, p(1) * C(3)) ≈ ev(kron([0. 0. ; 0. 1.], id, id) * cs[3])
end

@testset "Entanglement and entropies" begin
    L2 = log(2)
    s2 = System(2, Qubit())
    s4 = System(4, Qubit())
    # a bond dimension above the Schmidt rank keeps singular values of exactly zero, whose
    # 0 * log(0) made the entropy NaN; a product state on links of dimension 3, where 1.x took
    # one from dmrg with a mindim of 3
    up = zeros(3, 2, 3)
    up[1, 1, 1] = 1
    g = mps_state(s4, [up[1:1, :, :], up, up, up[:, :, 1:1]])
    @test entanglement_entropy(g, 2)[1] ≈ 0 atol=1e-12
    bell = (State{Pure}(s2, "Up") + State{Pure}(s2, "Dn")) / sqrt(2)
    ghz = (State{Pure}(s4, "Up") + State{Pure}(s4, "Dn")) / sqrt(2)
    product = State{Pure}(s4, "Up")
    fullymixed = State{Mixed}(s4, "FullyMixed")

    # a maximally entangled pair carries one bit, with a flat spectrum
    e, spectrum = entanglement_entropy(bell, 1)
    @test e ≈ L2
    @test spectrum ≈ [0.5, 0.5]
    # a GHZ state carries that same one bit wherever it is cut
    for cut in 1:3
        @test entanglement_entropy(ghz, cut)[1] ≈ L2
    end
    # the cut is on the right of the site, so the last site separates the whole state from
    # nothing, and a position outside the chain is no cut at all
    @test entanglement_entropy(ghz, 4)[1] ≈ 0 atol=1e-12
    @test_throws "entanglement entropy at site 0 of a 4 site state" entanglement_entropy(ghz, 0)
    @test_throws "entanglement entropy at site 5 of a 4 site state" entanglement_entropy(ghz, 5)
    # a product state carries none
    e, spectrum = entanglement_entropy(product, 2)
    @test e ≈ 0 atol=1e-12
    @test spectrum ≈ [1.]

    # a pure state stays pure, whichever representation holds it (trace2 in 1.x)
    for st in [bell, ghz, mix(ghz)]
        @test purity(st) ≈ 1
        @test renyi2(st) ≈ 0 atol=1e-12
    end
    # but one of its sites, seen alone, is maximally mixed
    @test renyi2(mix(ghz), [1]) ≈ L2
    @test renyi2(mix(ghz), [1, 2]) ≈ L2
    # and the pure representation must answer the same: a subsystem of a pure state is not
    # pure, so this is an entanglement measure and not the 0. the whole state gives
    @test renyi2(ghz, [1]) ≈ L2
    @test renyi2(ghz, [1, 2]) ≈ L2
    @test renyi2(ghz, [1, 3]) ≈ L2                  # the sites need not be contiguous
    @test renyi2(product, [1, 2]) ≈ 0 atol=1e-12
    # across a cut the mutual information is twice it, the cut 1 of 1.x written as the
    # positions 1:1 (decision 9.2)
    @test mutual_info_renyi2(bell, 1:1) ≈ 2 * renyi2(bell, [1])

    # the infinite temperature state of n qubits has purity 2^-n and entropy n log 2
    @test purity(fullymixed) ≈ 1 / 16
    @test renyi2(fullymixed) ≈ 4L2
    @test renyi2(fullymixed, [1]) ≈ L2
    @test hermiticity(fullymixed) ≈ 1

    # mutual information: 2 log 2 for GHZ, none between independent halves, the cut 2 of 1.x
    # written as the positions 1:2 (decision 9.2)
    @test mutual_info_renyi2(ghz, 1:2) ≈ 2L2
    @test mutual_info_renyi2(mix(ghz), 1:2) ≈ 2L2
    @test mutual_info_renyi2(fullymixed, 1:2) ≈ 0 atol=1e-12

    # on a pure state the mutual information of a leading block is read from the entanglement
    # spectrum rather than from partial traces: both routes must agree, on a state whose
    # spectrum is not degenerate, a random one where 1.x took one evolved by tdvp
    k = 5
    spread = random_state(Pure, System(k, Qubit()), 4)
    for cut in 1:k - 1
        @test mutual_info_renyi2(spread, 1:cut) ≈ mutual_info_renyi2(mix(spread), 1:cut)
        # the leading block given as a range, the cut of 1.x (decision 9.2), and as a list
        @test mutual_info_renyi2(spread, 1:cut) ≈ mutual_info_renyi2(spread, collect(1:cut))
    end
    # a non contiguous part of a pure state is twice its entropy as well, one partial trace
    # standing for the two parts and the whole of the mixed state
    @test mutual_info_renyi2(spread, [1, 3]) ≈ mutual_info_renyi2(mix(spread), [1, 3])
    # and so on fermions, for a state of definite parity
    ψf = random_state(State{Pure}(System(6, Fermion(conserve = parity(N))), "Emp"), 8)
    @test mutual_info_renyi2(ψf, [2, 5]) ≈ mutual_info_renyi2(mix(ψf), [2, 5])

    # the same quantities as measurements must agree with the functions, the cuts written as
    # positions (decision 9.2)
    m(state, f) = only(last.(measure(state, Measure(f))))
    @test m(ghz, Trace) ≈ 1
    @test m(ghz, Norm) ≈ 1
    # Trace2, the exact duplicate of Purity, is gone (decision 7b), and TraceError and
    # HermiticityError, named like exceptions, are TraceDeviation and HermiticityDeviation
    @test m(ghz, Purity) ≈ 1
    @test m(ghz, Renyi2) ≈ 0 atol=1e-12
    @test m(ghz, Hermiticity) ≈ 1
    @test m(ghz, TraceDeviation) ≈ 0 atol=1e-12
    @test m(ghz, HermiticityDeviation) ≈ 0 atol=1e-12
    @test m(ghz, MaxLinkdim) == 2
    @test m(ghz, EntanglementEntropy(2)) ≈ L2
    @test m(fullymixed, Purity) ≈ 1 / 16
    @test m(fullymixed, Renyi2) ≈ 4L2
    @test m(fullymixed, HermiticityDeviation) ≈ 0 atol=1e-12
    @test m(mix(ghz), SubRenyi2([1])) ≈ L2
    @test m(ghz, MutualInfoRenyi2(1:2)) ≈ 2L2
    @test m(mix(ghz), MutualInfoRenyi2(1:2)) ≈ 2L2
    @test MutualInfoRenyi2(1:2).name == "MutualInfoRenyi2(1,2)"
    # MutualInfoRenyi2(3), of the cut 3, was named as the one site part [3], so that the two
    # could not be measured together
    st4 = random_state(Pure, System(4, Qubit()), 3)
    r = Dict(measure(st4, [MutualInfoRenyi2(1:3), MutualInfoRenyi2([3])]))
    @test r["MutualInfoRenyi2(1:3)"] ≈ mutual_info_renyi2(st4, [1, 2, 3])
    @test r["MutualInfoRenyi2(3)"] ≈ mutual_info_renyi2(st4, [3])
    # the spectrum asked for is written whole, the eigenvalues beyond the bond dimension at the
    # cut being zeros, so that the rows of a file keep their width
    v = only(measure(State{Pure}(System(3, Qubit()), "Up"), EntanglementEntropy(1, 4)))
    @test first(v) == "EntanglementEntropy(1,4)"
    @test length(last(v)) == 5
    @test last(v)[1] ≈ 0 atol=1e-12
    @test last(v)[2] ≈ 1
    @test last(v)[3:5] == zeros(3)
    # and the entropy is that of entanglement_entropy, followed by its spectrum
    @test last(only(measure(bell, EntanglementEntropy(1, 2)))) ≈ [L2, 0.5, 0.5]
end

@testset "Entanglement resolved by sector" begin
    # the sectors add up to the entanglement entropy, the charge fluctuating between the two
    # sides bringing the number entropy on top of what each sector holds
    total(d) = sum(v.weight * v.entropy for v in values(d)) -
               sum(v.weight * log(v.weight) for v in values(d))
    f = System(4, Fermion(conserve = N))
    p(v) = State{Pure}(f, v)
    st = (p(["Occ", "Emp", "Emp", "Occ"]) + 0.5 * p(["Emp", "Occ", "Occ", "Emp"])) / sqrt(1.25)
    for pos in 1:4
        @test total(entanglement_by_sector(st, pos)) ≈ first(entanglement_entropy(st, pos))
    end

    # cut after one site, the charge on the left is either one particle or none
    d = entanglement_by_sector(st, 1)
    @test sort(collect(keys(d)); by = string) == [IT.QN("N", 0), IT.QN("N", 1)]
    @test d[IT.QN("N", 1)].weight ≈ 0.8
    @test d[IT.QN("N", 0)].weight ≈ 0.2
    # cut after two, it is one particle in both terms, and the entanglement lives inside it
    d = entanglement_by_sector(st, 2)
    @test collect(keys(d)) == [IT.QN("N", 1)]
    @test d[IT.QN("N", 1)].weight ≈ 1
    @test d[IT.QN("N", 1)].spectrum ≈ [0.8, 0.2]
    @test d[IT.QN("N", 1)].entropy ≈ first(entanglement_entropy(st, 2))

    # several quantities at once, and a charge modulo 2 over a basis whose charges repeat
    e = System(3, Electron(conserve = (Ntot, 2Sz)))
    q(v) = State{Pure}(e, v)
    se = (q(["Up", "Dn", "Emp"]) + 0.5 * q(["UpDn", "Emp", "Emp"])) / sqrt(1.25)
    d = entanglement_by_sector(se, 1)
    @test d[IT.QN(("Ntot", 1), ("2Sz", 1))].weight ≈ 0.8
    @test d[IT.QN(("Ntot", 2), ("2Sz", 0))].weight ≈ 0.2
    b = System(3, Boson(4, conserve = parity(N)))
    r(v) = State{Pure}(b, v)
    d = entanglement_by_sector((r(["1", "2", "0"]) + 0.5 * r(["0", "3", "0"])) / sqrt(1.25), 1)
    @test d[IT.QN("parity(N)", 1, 2)].weight ≈ 0.8
    @test d[IT.QN("parity(N)", 0, 2)].weight ≈ 0.2

    # without charges there is a single sector holding the whole spectrum
    fd = System(4, Fermion())
    pd(v) = State{Pure}(fd, v)
    sd = (pd(["Occ", "Emp", "Emp", "Occ"]) + 0.5 * pd(["Emp", "Occ", "Occ", "Emp"])) / sqrt(1.25)
    d = entanglement_by_sector(sd, 2)
    @test collect(keys(d)) == [IT.QN()]
    @test d[IT.QN()].spectrum ≈ last(entanglement_entropy(sd, 2))

    # a mixed representation is refused, its links pairing the charge of a ket with a bra
    @test_throws "takes a pure state" entanglement_by_sector(mix(st), 1)
    @test_throws "between 1 and 4" entanglement_by_sector(st, 5)
end

@testset "Composed operators as measurements" begin
    # a measurement may be any operator of one site and not just a named one, named by how it
    # prints
    sys = System(3, Qubit())
    st = State{Pure}(sys, ["X+", "Y+", "Z+"])
    colnames(o) = first.(measure(st, [o]))
    # a channel is not an observable, and expect says so (decision 1: Left(X)(1) in 1.x)
    @test_throws "is a channel" expect(mix(st), Gate(X)(1))
    @test_throws "measure takes observables, and Gate(X)(1) is not one" measure(mix(st), Gate(X)(1))
    @test colnames(X) == ["X"]
    @test colnames(X * Y) == ["X*Y"]
    @test colnames(2X) == ["2X"]
    @test colnames(X + Y) == ["X+Y"]
    @test colnames(dag(S)) == ["dag(S)"]
    @test colnames(Id) == ["Id"]
    # a pair is named (A, B), where 1.x glued the two names, IdX and XY (layer 3, question 7)
    @test colnames((Id, X)) == ["(Id, X)"]
    @test colnames((X, Y)) == ["(X, Y)"]
    # an operator placed on sites by its compact printed form
    @test colnames(X(1) * Y(2) + 0.5 * Z(3)) == ["X(1)*Y(2)+0.5Z(3)"]
    # a coefficient rides outside the operator that `(c * A)(i)` builds, so it has to come off
    # where the tensor of one site is made and not only where the name is
    @test expect1(st, 2X) ≈ 2 * expect1(st, X)
    @test expect1(mix(st), 2X) ≈ 2 * expect1(mix(st), X)
    @test expect1(st, 0.5 * (X + Y)) ≈ 0.5 * (expect1(st, X) + expect1(st, Y))
    @test expect2(st, (2X, 3Y)) ≈ 6 * expect2(st, (X, Y))
end

@testset "Expectation values on an unnormalised state" begin
    # the environments do not carry the normalisation, which expect applies once at the end:
    # a state of norm 1/2 whose orthogonality centre is on site 4, so that the shortcut of the
    # left orthogonal sites is taken, gives the values of the normalised one
    sys = System(6, Qubit())
    base = State{Pure}(sys, ["Up", "+", "Up", "+", "+", "+"])
    st = State(base, 0.5 * TMS.ITensorMPS.orthogonalize(base.state, 4))
    @test TMS.ITensorMPS.leftlim(st.state) ≥ 1
    @test trace(st) ≈ 0.25
    @test expect1(st, Z) ≈ [1., 0., 1., 0., 0., 0.] atol = 1e-12
    @test expect1(st, X) ≈ [0., 1., 0., 1., 1., 1.] atol = 1e-12
    @test expect1(st, X) ≈ expect1(normalize(st), X)
    @test expect(st, Z(1) * Z(3)) ≈ expect(normalize(st), Z(1) * Z(3))
    @test expect2(st, (X, X)) ≈ expect2(normalize(st), (X, X))
    # the same built by two gates that are not unitary, as 1.x did: Sp|+> = |0>/sqrt(2)
    sp = apply(Sp(3) * Sp(1), State{Pure}(sys, "+"))
    @test trace(sp) ≈ 0.25
    @test expect1(sp, Z) ≈ [1., 0., 1., 0., 0., 0.] atol = 1e-12
    @test expect2(sp, (X, X)) ≈ expect2(normalize(sp), (X, X))
    # on a charged state the shortcut has arrows to respect
    q = State{Pure}(System(4, Qubit(conserve = N)), ["Up", "Up", "Dn", "Dn"])
    a = State(q, TMS.ITensorMPS.orthogonalize(q.state, 3))
    @test TMS.ITensorMPS.leftlim(a.state) ≥ 1
    @test expect1(a, N) ≈ [0., 0., 1., 1.] atol = 1e-12
    @test expect(a, N(3) * N(4)) ≈ 1
end

@testset "Positions given as a range" begin
    # the position arguments took a `Vector{Int}` and nothing else, so a range was refused
    # by a `MethodError`. The name of the measurement was built happily by
    # `compact_positions`, so the failure came only when the measurement was taken. The
    # partial trace of a range, and its refusal of a pure state, are in "A partial trace keeps
    # the trace, and what it refuses"
    sys = System(6, Qubit())
    stp = random_state(Pure, sys, 4)
    stm = mix(stp)
    @test renyi2(stm, 1:3) ≈ renyi2(stm, [1, 2, 3])
    @test renyi2(stp, 1:3) ≈ renyi2(stp, [1, 2, 3])
    @test mutual_info_renyi2(stm, 1:3) ≈ mutual_info_renyi2(stm, [1, 2, 3])
    # `measure` hands back a vector of name => value pairs, one per measurement
    @test last(only(measure(stm, SubRenyi2(1:3)))) ≈ renyi2(stm, [1, 2, 3])
    @test last(only(measure(stm, MutualInfoRenyi2(1:3)))) ≈ mutual_info_renyi2(stm, [1, 2, 3])
    @test first(only(measure(stm, SubRenyi2(1:3)))) == "SubRenyi2(1:3)"
end

@testset "A cut at either end" begin
    # a cut at 0 gave 0 on a mixed state and an error about the entanglement entropy on a
    # pure one, and a negative cut was taken for an empty part on a mixed state. The cuts of
    # 1.x are the positions 1:cut (decision 9.2): the cuts 0 and 3 an empty part and the
    # whole state, the cuts -1 and 4 sites the state does not have
    ψ = random_state(Pure, System(3, Qubit()), 2)
    for st in (ψ, mix(ψ))
        @test mutual_info_renyi2(st, 1:0) == 0
        @test mutual_info_renyi2(st, 1:3) == 0
        @test_throws "mutual_info_renyi2 was given site -1" mutual_info_renyi2(st, [-1])
        @test_throws "mutual_info_renyi2 was given site 4" mutual_info_renyi2(st, 1:4)
    end

    # SubRenyi2(3) was accepted and failed when measured: a leading block, on a pure state, is
    # read off the spectrum
    for cut in 1:2
        @test renyi2(ψ, 1:cut) ≈ renyi2(mix(ψ), collect(1:cut))
        @test renyi2(mix(ψ), 1:cut) ≈ renyi2(mix(ψ), collect(1:cut))
    end
    @test renyi2(ψ, 1:0) == 0
    # no site has no entropy, as the cut 0, where partial_trace refused to keep nothing
    @test renyi2(ψ, Int[]) == 0
    @test renyi2(mix(ψ), Int[]) == 0
    @test last(only(measure(ψ, SubRenyi2(Int[])))) == 0
    # positions in a vector of another element type failed at the first measurement
    @test last(only(measure(ψ, SubRenyi2([])))) == 0
    @test last(only(measure(mix(ψ), SubRenyi2([])))) == 0
    @test renyi2(ψ, Any[1, 2]) ≈ renyi2(ψ, [1, 2])
    @test mutual_info_renyi2(ψ, Any[1]) ≈ mutual_info_renyi2(ψ, [1])
    @test mutual_info_renyi2(ψ, []) == 0
    @test renyi2(mix(ψ), 1:3) ≈ renyi2(mix(ψ))
    @test_throws "renyi2 was given site 4" renyi2(ψ, 1:4)
    # the cut 2 of 1.x, SubRenyi2(2), written as the positions 1:2 (decision 9.2)
    @test SubRenyi2(1:2).name == "SubRenyi2(1,2)"
    @test last(only(measure(ψ, SubRenyi2(1:2)))) ≈ renyi2(mix(ψ), [1, 2])
end

@testset "Measuring a multiple of the identity" begin
    st = random_state(Pure, System(3, Qubit()), 4)
    for s in (st, mix(st))
        @test expect1(s, 2Id) ≈ fill(2, 3)
        @test expect1(s, (2X)^2) ≈ fill(4, 3)
        @test expect2(s, (2Id, X)) ≈ 2 * expect2(s, (Id, X))
        @test expect(s, 3 * Id(1)) ≈ 3
        @test last(only(measure(s, (0.5Z)^2))) ≈ fill(0.25, 3)
        @test last(only(measure(s, 3 * Id(2)))) ≈ 3
    end
end

@testset "Terms of one site next to a term over many sites" begin
    # the terms of one site are not taken for rounding next to a term whose norm of Frobenius
    # grows with its number of sites
    q = State{Pure}(System(100, Qubit()), "Up")
    h = sum([ Z(i) for i in 1:100 ]) + prod([ X(i) for i in 1:100 ])
    @test expect(q, h) ≈ 100
    @test inner(q, h, q) ≈ 100
    # a periodic Hubbard chain, whose hopping across the edge has a string over every site
    n, u = 50, 4.
    e = State{Pure}(System(n, Electron()), "UpDn")
    bonds = [ (i, mod1(i + 1, n)) for i in 1:n ]
    h = -sum([ dag(c)(i) * c(j) + dag(c)(j) * c(i) for c in (Cup, Cdn) for (i, j) in bonds ]) +
        u * sum([ Nupdn(i) for i in 1:n ])
    @test expect(e, h) ≈ u * n
    @test inner(e, h, e) ≈ u * n
end

@testset "Projector of several sites" begin
    # applied as a gate, with the Jordan-Wigner strings of the sites in between, the order of
    # its sites included; renamed by named, which takes no site (decision 2)
    ud = State{Pure}(System(2, Qubit()), ["Up", "Dn"])
    @test real(expect(ud, named(Proj(Swap => -1), "Singlet")(1, 2))) ≈ 0.5
    for st in (ud, mix(ud))
        @test real(expect(normalize(apply(Proj(Swap => -1)(1, 2), st)), Swap(1, 2))) ≈ -1
    end
    # 1.x refused to measure it as acting on several sites at once; the descent splits its
    # matrix on its sites (decision 2)
    @test expect(ud, Proj(Swap => -1)(1, 2)) ≈ 0.5
    fs = System(4, Fermion(conserve = N))
    f = random_state(State{Pure}(fs, ["Occ", "Occ", "Emp", "Emp"]), 4)
    projected(a, λ) = normalize(apply(make_mpo(f, (a * a + λ * a) / 2), f))
    h = dag(C)(1) * C(3) + dag(C)(3) * C(1)
    j = im * (dag(C)(1) * C(3) - dag(C)(3) * C(1))
    hop = dag(C) ⊗ C - C ⊗ dag(C)
    current = im * (dag(C) ⊗ C + C ⊗ dag(C))
    for (g, ref) in [(Proj(hop => 1)(1, 3), projected(h, 1)),
                     (Proj(hop => 1)(3, 1), projected(h, 1)),
                     (Proj(current => 1)(1, 3), projected(j, 1)),
                     (Proj(current => 1)(3, 1), projected(j, -1))]
        @test norm(normalize(apply(g, f)) - ref) < 1e-12
        @test norm(normalize(apply(g, mix(f))) - mix(ref)) < 1e-12
    end
end

@testset "A random fermionic oracle" begin
    # the Jordan-Wigner strings checked on random products against the dense matrices of the
    # strings written out, through expect on a pure and on a mixed state, the MPO, apply on a
    # pure and on a mixed state, the gates of functions of several sites, expect2, the matrix
    # of a tensor product and partial_trace. In the dense basis C(j) is
    # F ⊗ … ⊗ F ⊗ c ⊗ 1 ⊗ … ⊗ 1, and any other operator is written from those, so that the
    # reference knows nothing of the package
    fe = Fermion()
    mc, mf, mi = matrix(C, fe), matrix(F, fe), matrix(Id, fe)
    cd(j, n) = foldl(kron, [ k < j ? mf : k == j ? mc : mi for k in 1:n ])
    singles = [ (C, (j, n) -> cd(j, n)), (dag(C), (j, n) -> cd(j, n)'),
                (N, (j, n) -> cd(j, n)' * cd(j, n)),
                (C + 2dag(C), (j, n) -> cd(j, n) + 2cd(j, n)') ]
    count = Ref(0)
    # a factor on the positions 1:m, as an operator and as its dense matrix on the n sites of
    # the state, `site` taking a position to its site: an operator of one site, or a tensor
    # product of two, renamed or not, on two positions in either order
    function factor(m, site, n)
        a, da = rand(singles)
        i = rand(1:m)
        if m == 1 || rand(Bool)
            return a(i), da(site(i), n)
        end
        b, db = rand(singles)
        j = rand(setdiff(1:m, i))
        t = rand(Bool) ? a ⊗ b : named(a ⊗ b, "T$(count[] += 1)")
        return t(i, j), da(site(i), n) * db(site(j), n)
    end
    function product(m, site, n)
        fs = [ factor(m, site, n) for _ in 1:rand(1:3) ]
        return prod(first.(fs)), prod(last.(fs))
    end
    n = 4
    ψ = random_state(Pure, System(n, fe), 4)
    ρ = mix(ψ)
    v = dense_vector(ψ)
    value(d) = v' * d * v
    w = raw_vector(ψ)
    for _ in 1:30
        o, d = product(n, identity, n)
        @test expect(ψ, o) ≈ value(d) atol = 1e-10
        @test expect(ρ, o) ≈ value(d) atol = 1e-10
        @test mpo_matrix(ψ, o) ≈ d atol = 1e-10
        @test raw_vector(apply(o, ψ)) ≈ d * w atol = 1e-10
        @test norm(apply(o, ρ) - mix(apply(o, ψ))) < 1e-10
        keep = sort(TMS.Random.randperm(n)[1:rand(1:3)])
        o, d = product(length(keep), p -> keep[p], n)
        @test expect(partial_trace(ρ, keep; keep = true), o) ≈ value(d) atol = 1e-10
    end
    # a function of an even product applied as a gate, on sites in any order and apart
    odd = singles[[1, 2, 4]]
    for _ in 1:10
        fs = Any[ rand(odd), rand(odd) ]
        if rand(Bool)
            insert!(fs, rand(1:3), singles[3])
        end
        sites = TMS.Random.randperm(n)[1:length(fs)]
        a = reduce(⊗, first.(fs))
        da = prod(f[2](s, n) for (f, s) in zip(fs, sites))
        g = exp(0.3 * a - 0.3 * dag(a))(sites...)
        @test raw_vector(apply(g, ψ)) ≈ exp(0.3 * da - 0.3 * da') * w atol = 1e-10
        @test norm(apply(g, ρ) - mix(apply(g, ψ))) < 1e-10
    end
    # and on electrons, whose parity counts both spins
    el = Electron()
    mu, md, mfe, mie = matrix(Cup, el), matrix(Cdn, el), matrix(F, el), matrix(Id, el)
    ce(j, m) = foldl(kron, [ k < j ? mfe : k == j ? m : mie for k in 1:3 ])
    ψe = random_state(Pure, System(3, el), 4)
    ge = exp(-0.4im * (dag(Cup) ⊗ Cdn + dag(dag(Cup) ⊗ Cdn)))
    for (i, j) in [(1, 3), (3, 1), (1, 2)]
        h = ce(i, mu)' * ce(j, md) + ce(j, md)' * ce(i, mu)
        @test raw_vector(apply(ge(i, j), ψe)) ≈ exp(-0.4im * h) * raw_vector(ψe) atol = 1e-10
    end
    for (a, da) in singles, (b, db) in singles
        ref = [ value(i == j ? da(i, n) * db(i, n) : da(i, n) * db(j, n)) for i in 1:n, j in 1:n ]
        @test expect2(ψ, (a, b)) ≈ ref atol = 1e-10
        @test expect2(ρ, (a, b)) ≈ ref atol = 1e-10
    end
    # a tensor product on consecutive sites, of operators of one site and of renamed tensor
    # products of two, each with its width and its dense matrix from its first site
    for _ in 1:20
        fs = map(1:rand(2:3)) do _
            a, da = rand(singles)
            if rand(Bool)
                return (a, 1, da)
            end
            b, db = rand(singles)
            return (named(a ⊗ b, "T$(count[] += 1)"), 2, (j, n) -> da(j, n) * db(j + 1, n))
        end
        starts = cumsum([ 1; [ w for (_, w, _) in fs ] ])
        m = starts[end] - 1
        ref = prod(d(s, m) for ((_, _, d), s) in zip(fs, starts))
        @test matrix(reduce(⊗, first.(fs)), fe) ≈ ref atol = 1e-12
    end
end

@testset "apply cuts a product into pieces" begin
    # a product is applied factor by factor, never developed: a hundred Swaps, whose product
    # developed would have about 4^100 terms
    q = System(5, Qubit())
    st = State{Pure}(q, ["Dn", "Up", "Up", "Up", "Up"])
    swaps = [ Swap(mod1(i, 4), mod1(i, 4) + 1) for i in 1:100 ]
    op = prod(swaps)
    @test length(last(TMS.pieces(q, op))) == 100
    one_by_one = foldl((s, g) -> apply(g, s), reverse(swaps); init = st)
    @test norm(apply(op, st) - one_by_one) < 1e-10
    # gates far apart move sites, which come back to their places
    ψ = random_state(Pure, System(8, Fermion()), 4)
    far = (C ⊗ dag(C))(2, 7) * Swap(1, 5) * exp(-0.3im * (dag(C) ⊗ C + dag(dag(C) ⊗ C)))(3, 8)
    @test TMS.ITensorMPS.siteinds(apply(far, ψ).state) == TMS.ITensorMPS.siteinds(ψ.state)
    # factors of one site of definite parity gather into one piece, across factors they
    # commute with, so that their strings compensate
    fe = System(5, Fermion())
    @test length(last(TMS.pieces(fe, dag(C)(3) * N(1) * C(5)))) == 1
    @test length(last(TMS.pieces(fe, C(1) * Swap(2, 3) * dag(C)(4)))) == 2
    # but not across a factor of no definite parity, and across an odd one with a sign
    @test length(last(TMS.pieces(fe, C(1) * (C + N)(2) * dag(C)(4) * N(5)))) == 3
    ψ = random_state(Pure, fe, 4)
    for op in [C(1) * (C ⊗ N)(2, 3) * dag(C)(4), C(1) * (C + N)(2) * dag(C)(4) * N(5)]
        @test norm(apply(op, ψ) - apply(make_mpo(ψ, op), ψ)) < 1e-10
    end
    # a dense block on its sites when it lays no string elsewhere, its MPO otherwise, and
    # gates of one site for a single term
    h = dag(C) ⊗ C + dag(dag(C) ⊗ C)
    steps(op) = last(TMS.pure_steps(fe, op))
    @test length(only(steps(exp(-0.4im * h)(1, 2))).gates) == 1
    @test only(steps(exp(-0.4im * h)(1, 3))) isa TMS.WindowMPO
    @test only(steps((C + N)(3))) isa TMS.WindowMPO
    @test length(only(steps(C(4))).gates) == 4
    @test length(only(steps((C + dag(C))(3) * named(N + 2Id, "G")(5))).gates) == 4
    # a product that vanishes gives the null state, and an MPO is applied as it is
    ψ = random_state(Pure, fe, 4)
    @test norm(apply(C(2) * C(2), ψ)) == 0
    mpo = make_mpo(ψ, dag(C)(1) * C(3) + N(2))
    @test norm(apply(mpo, ψ) - apply(dag(C)(1) * C(3) + N(2), ψ)) < 1e-10
end

@testset "The algorithms of the product of an MPO by a state" begin
    # without truncation, "naive", the default, and "densitymatrix" give the exact product; on
    # a random state, where 1.x took one evolved by tdvp
    n = 4
    h = sum(X(i) * X(i + 1) for i in 1:n - 1)
    st = random_state(Pure, System(n, Qubit()), 4)
    m = make_mpo(st, h)
    @test norm(apply(m, st; apply_algo = "naive") - apply(m, st)) < 1e-12
    @test norm(apply(m, st; apply_algo = "densitymatrix") - apply(m, st)) < 1e-12
    up = State{Pure}(System(n, Qubit()), "Up")
    w(apply_algo) = evolve(-im * h, 0.5, up; algo = ApproxW(; order = 2, apply_algo), nsteps = 5)
    @test norm(w("naive") - w("densitymatrix")) < 1e-12
    # "fit" needs a number of sweeps of its own, "zipup" is missing for a state from the lowest
    # ITensorMPS accepted, and a name that is none of them is refused too
    for alg in ("fit", "zipup", "exact")
        @test_throws "apply_algo is" apply(m, st; apply_algo = alg)
        @test_throws "apply_algo is" ApproxW(order = 2, apply_algo = alg)
    end
    # and an Evolve phase with ApproxW hands it on, as evolve does
    sim = runTMS(SimData(phases = [CreateState{Pure}(n, Qubit(), "Up"),
                                   Evolve(evolver = -im * h, duration = 0.5, time_step = 0.1,
                                          algo = ApproxW(order = 2, apply_algo = "densitymatrix"))]);
                 output = devnull)
    @test norm(sim.state - w("densitymatrix")) < 1e-12
end

@testset "Gates on sites apart, and MPOs on the sites they act on" begin
    # a gate on sites apart goes through its MPO, decomposed by singular value decompositions,
    # the identity on the sites between, rather than bringing its sites together; a piece
    # developed goes through its MPO on the sites it acts on alone, from site 1 when its strings
    # start there; a single term, its strings laid, through gates of one site
    exact(st, op) = TMS.system_matrix(Pure, st.system, TMS.lower(st.system, op)) * raw_vector(st)
    function same(st, op)
        e = exact(st, op)
        return norm(raw_vector(apply(op, st)) - e) ≤ 1e-12 * max(1, norm(e))
    end
    q = random_state(Pure, System(8, Qubit()), 4)
    g4 = exp(-0.3im * (X ⊗ X + Y ⊗ Y + Z ⊗ Z))
    for op in [controlled(X)(2, 7), g4(1, 8), controlled(controlled(X))(1, 4, 8),
               controlled(X)(1, 6) * H(3) * controlled(Z)(2, 8) * X(5)]
        @test same(q, op)
    end
    @test only(last(TMS.pure_steps(q.system, controlled(X)(2, 7)))) isa TMS.WindowMPO
    @test only(last(TMS.pure_steps(q.system, controlled(X)(2, 3)))) isa TMS.GateBatch
    hop = exp(-0.7im * (dag(C) ⊗ C + dag(dag(C) ⊗ C)))
    f = random_state(Pure, System(7, Fermion()), 4)
    fq = apply(prod(hop(k, k + 1) for k in (1, 3, 5, 2, 4, 6)),
               State{Pure}(System(7, Fermion(conserve = N)), [ isodd(k) ? "Occ" : "Emp" for k in 1:7 ]))
    for st in (f, fq), op in [hop(2, 6), C(2) + 0.5 * C(5), (C ⊗ N)(1, 4), (N ⊗ C)(2, 5),
                              C(3) * dag(C)(6) * C(7)]
        @test same(st, op)
    end
    @test same(f, (C + dag(C))(4))
    # a result that vanishes through the MPO on its sites, on charged sites
    z = State{Pure}(System(5, Fermion(conserve = N)), ["Occ", "Emp", "Emp", "Emp", "Occ"])
    @test norm(apply(C(2) + C(3), z)) == 0
end

@testset "A result that vanishes on charged sites" begin
    # ITensors cannot factorise a tensor of charged indices that is exactly zero, which has no
    # block: whatever takes a state there gives the null state, through gates of one site or of
    # several, near or far apart, or through an MPO
    sys = System(5, Fermion(conserve = N))
    ψ = State{Pure}(sys, ["Emp", "Occ", "Emp", "Occ", "Emp"])
    h = dag(C) ⊗ C + dag(dag(C) ⊗ C)
    for a in ((C ⊗ C)(1, 2), (C ⊗ C)(1, 4), C(1) * Swap(1, 3), C(1) + C(3),
              exp(-0.3im * h)(1, 3) * (C ⊗ C)(1, 3))
        @test norm(apply(a, ψ)) == 0
    end
    @test norm(apply(make_mpo(ψ, C(1) + C(3)), ψ)) == 0
    @test norm(apply((Sp ⊗ Sp)(1, 2), State{Pure}(System(3, Qubit(conserve = N)), "Up"))) == 0
    # and what does not vanish is left as it was
    @test norm(apply((dag(C) ⊗ C)(1, 2), ψ)) ≈ 1
end

@testset "Inner products with an operator" begin
    q = System(4, Qubit())
    u, v = random_state(Pure, q, 3), random_state(Pure, q, 3)
    m = X(1) * Z(3) + 0.3 * Y(2) + Swap(2, 4)
    du, dv = dense_vector(u) * norm(u), dense_vector(v) * norm(v)
    @test inner(u, m, v) ≈ du' * TMS.system_matrix(Pure, q, TMS.lower(q, m)) * dv
    @test dot(u, m, v) ≈ inner(u, m, v)
    # between two density matrices a pure operator is either its gate or its generator
    @test_throws "is ambiguous" inner(mix(u), X(1), mix(v))
end

@testset "Inner products and fidelities" begin
    sys = System(3, Qubit())
    up    = State{Pure}(sys, "Up")
    dn    = State{Pure}(sys, "Dn")
    plus  = State{Pure}(sys, "+")
    mixup = State{Pure}(sys, ["+", "+", "Up"])
    icplx = State{Pure}(sys, ["i", "+", "+"])

    # <+++|++Up> = 1 * 1 * <+|Up> = 1/sqrt(2)
    @test inner(plus, mixup) ≈ 1/√2
    @test dot(plus, mixup) == inner(plus, mixup)          # dot is an alias
    @test inner(plus, plus) ≈ 1
    @test inner(up, dn) ≈ 0 atol = 1e-14
    # the first argument is the one conjugated
    @test inner(icplx, plus) ≈ conj(inner(plus, icplx))
    @test imag(inner(icplx, plus)) ≉ 0                    # the test would be empty otherwise

    # fidelity is normalised, so neither the norm nor the trace of its arguments matters
    @test fidelity(plus, mixup) ≈ 0.5
    @test fidelity(plus, plus) ≈ 1
    @test fidelity(up, dn) ≈ 0 atol = 1e-14
    @test fidelity(3 * plus, 2 * mixup) ≈ fidelity(plus, mixup)
    # against a mixed representation, in either order
    fm = State{Mixed}(sys, "FullyMixed")
    @test fidelity(plus, mix(mixup)) ≈ fidelity(plus, mixup)
    @test fidelity(mix(mixup), plus) ≈ fidelity(plus, mixup)
    @test fidelity(up, fm) ≈ 1/8
    @test fidelity(2 * up, fm) ≈ 1/8
    # the normalised Hilbert-Schmidt overlap of two mixed states
    @test hs_fidelity(mix(plus), mix(plus)) ≈ 1
    @test hs_fidelity(mix(plus), mix(mixup)) ≈ fidelity(plus, mixup)
    @test hs_fidelity(mix(up), fm) ≈ √2/4
    @test hs_fidelity(2 * mix(up), 3 * fm) ≈ hs_fidelity(mix(up), fm)

    # a pure and a mixed representation are vectors of different spaces
    @test_throws "pure and a mixed representation" inner(plus, fm)
    @test_throws "pure and a mixed representation" inner(fm, plus)
    # and two mixed ones have no fidelity, the refusal naming hs_fidelity
    @test_throws "no fidelity between two mixed" fidelity(fm, mix(plus))
    @test_throws "hs_fidelity" fidelity(fm, mix(plus))
end

@testset "Fidelity and Overlap as measurements" begin
    # the reference is written with the measurement, on a system of its own
    ref = State{Pure}(System(3, Qubit()), "+")
    st = State{Pure}(System(3, Qubit()), ["+", "+", "Up"])
    @test first(only(measure(st, Fidelity(ref)))) == "Fidelity"
    @test last(only(measure(st, Fidelity(ref)))) ≈ 0.5
    @test first(only(measure(st, Overlap(ref)))) == "Overlap"
    @test last(only(measure(st, Overlap(ref)))) ≈ 1/√2
    @test last(only(measure(mix(st), Fidelity(ref)))) ≈ 0.5
    # and through a run, where the system does not exist until the phase creates it, or with a
    # reference the run names, missing until there is one
    sim = runTMS(SimData(phases = [CreateState{Pure}(3, Qubit(), ["+", "+", "Up"];
                                                     final_measurements = Data("d") => [
                                                         Fidelity(ref), Overlap(ref), "named" => Fidelity(:ref)]),
                                   CopyState(:state => :ref; final_measurements = Data("d") => [
                                                 Fidelity(ref), Overlap(ref), "named" => Fidelity(:ref)])]);
                 output = devnull)
    first_row, second_row = sim.data["d"]
    @test last.(first_row[2:3]) ≈ [0.5, 1/√2]
    @test ismissing(last(first_row[4]))
    @test last(second_row[4]) ≈ 1

    # a reference built before the state was weakened is weakened in its turn to what the state
    # conserves, on a system mixing sites that conserve and sites that do not as well
    for sites in ([Fermion(conserve = strong(N)) for _ in 1:3],
                  [Fermion(conserve = strong(N)), Qubit(), Fermion(conserve = strong(N))])
        p = State{Pure}(System(sites), ["Occ", sites[2] isa Qubit ? "Up" : "Emp", "Occ"])
        @test last(only(measure(weaken(mix(p)), Fidelity(p)))) ≈ 1
        @test last(only(measure(weaken(mix(p), ()), Fidelity(p)))) ≈ 1
        @test last(only(measure(weaken(p), Overlap(p)))) ≈ 1
    end
    # and a pure reference on a state dropping part of what it conserves, a site left conserving
    # nothing being cut into a single block on a system still charged, as 1.x tested it
    sys = System([Qubit(conserve = N), Boson(3, conserve = named(N, "Nb")),
                  Boson(3, conserve = named(N, "Nb")), Qubit(conserve = N)])
    p(v) = State{Pure}(sys, v)
    ψ = normalize(p(["1", "2", "0", "0"]) + 0.5 * p(["0", "1", "1", "1"]) -
                  0.7im * p(["1", "0", "2", "0"]))
    for target in (N, named(N, "Nb"))
        @test last(only(measure(weaken(mix(ψ), target), Fidelity(ψ)))) ≈ 1
    end
    # but a reference conserving less than the state is refused, weakening the state at every
    # measurement costing a whole conversion of it
    strong_state = State{Pure}(System(2, Fermion(conserve = strong(N))), "Occ")
    weak_ref = State{Pure}(System(2, Fermion(conserve = N)), "Occ")
    @test_throws "the reference of Fidelity conserves less than the measured state" measure(
        mix(strong_state), Fidelity(weak_ref))
    @test_throws "which conserves strong(N)" measure(mix(strong_state), Fidelity(weak_ref))
    sq = State{Pure}(System(2, Qubit(conserve = N)), "Up")
    @test_throws "the reference of Overlap conserves less than the measured state" measure(
        sq, Overlap(State{Pure}(System(2, Qubit()), "Up")))
end

@testset "A named state measured, and the matrix elements of Overlap" begin
    # a list :name => [...] measures the state of that name with everything measure takes, its
    # names prefixed by the name of the state, nested lists included
    sys = System(4, Qubit())
    a = random_state(Pure, sys, 4; rng = Xoshiro(1))
    b = random_state(Pure, sys, 4; rng = Xoshiro(2))
    r = measure(a, [Z(1), :b => [Z, (X, X), Purity, "e" => c -> c.time, :a => [Z(2)]]];
                time = 0.5, a, b)
    @test first.(r) == ["Z(1)", "b.Z", "b.(X, X)", "b.Purity", "b.e", "b.a.Z(2)"]
    @test last(r[2]) ≈ real(expect1(b, Z))
    @test last(r[3]) ≈ real(expect2(b, (X, X)))
    @test last(r[5]) == 0.5
    @test last(r[6]) ≈ real(expect(a, Z(2)))
    # a run refuses a named state that no phase provides
    @test_throws "reads nowhere" SimData(phases = [CreateState{Pure}(2, Qubit(), "Up";
        final_measurements = "d" => [:nowhere => [Z]])])
    # the matrix elements <ref|op|state> of two pure states, a vector for an operator of one site
    @test last(only(measure(a, Overlap(:b, X(1) * Z(3)); b))) ≈ inner(b, apply(X(1) * Z(3), a))
    @test last(only(measure(a, Overlap(:a, X(1) * Z(3)); a))) ≈ expect(a, X(1) * Z(3))
    @test last(only(measure(a, Overlap(b, Z)))) ≈ [ inner(b, apply(Z(i), a)) for i in 1:4 ]
    @test_throws "matrix element of two pure states" measure(mix(a), Overlap(:b, Z); b = mix(b))
end

@testset "Energy variance" begin
    sys = System(6, Qubit())
    ising = -sum(Z(i) * Z(i+1) for i in 1:5)
    field = -sum(X(i) for i in 1:6)
    h = ising + field

    # zero exactly on an eigenstate, whichever one
    @test variance(State{Pure}(sys, "Up"), ising) ≈ 0 atol = 1e-12
    @test variance(State{Pure}(sys, ["Up", "Dn", "Up", "Dn", "Up", "Dn"]), ising) ≈ 0 atol = 1e-12
    @test variance(State{Pure}(sys, "+"), field) ≈ 0 atol = 1e-12
    # |+...+> is an eigenstate of the field but not of the couplings. <ZZ> is zero there and
    # so are the cross terms, so the variance is the number of couplings
    @test variance(State{Pure}(sys, "+"), h) ≈ 5

    # against the naive route, the one that forms H^2 and that the implementation avoids
    r = random_state(Pure, sys, 8)
    @test variance(r, h) ≈ real(expect(r, h * h)) - real(expect(r, h))^2
    # neither the norm nor a phase of the state changes it
    @test variance(3r, h) ≈ variance(r, h)
    @test variance(im * r, h) ≈ variance(r, h)
    # an mpo already built gives the same answer
    @test variance(r, make_mpo(r, h)) ≈ variance(r, h)

    # what it is for: a converged ground state, whose energy is the lowest eigenvalue of the
    # dense hamiltonian, has a variance the energy alone cannot show
    e, gs = ground_state(h, random_state(Pure, sys, 16); algo = Dmrg(nsweeps = 12),
                         limits = Limits(cutoff = 1e-14, maxdim = 32))
    on(m, i) = kron([ k == i ? m : [1. 0. ; 0. 1.] for k in 1:6 ]...)
    dense = -sum(on([1. 0. ; 0. -1.], i) * on([1. 0. ; 0. -1.], i + 1) for i in 1:5) -
            sum(on([0. 1. ; 1. 0.], i) for i in 1:6)
    @test e ≈ minimum(TMS.LinearAlgebra.eigvals(TMS.LinearAlgebra.Symmetric(dense))) atol = 1e-10
    @test variance(gs, h) < 1e-8
    @test variance(gs, h) < variance(r, h)

    # as a measurement, on the random state where 1.x took the ground state
    @test first(only(measure(r, Variance(h)))) == "Variance"
    @test last(only(measure(r, Variance(h)))) ≈ variance(r, h)

    # a hamiltonian on a density matrix has no variance in this sense
    @test_throws "needs a pure representation" variance(mix(r), h)
    # and Variance, which builds the MPO when it is prepared, refuses it there, before anything
    # is measured, the Measure keeping no preparation for it
    called = Ref(false)
    @test_throws "needs a pure representation" measure(mix(r), ["f" => c -> (called[] = true),
                                                                Variance(h)])
    @test !called[]
    m = Measure(Variance(h))
    @test_throws "needs a pure representation" measure(mix(r), m)
    @test isempty(m.prepared)
    @test last(only(measure(r, m))) ≈ variance(r, h)
end

@testset "Observables do not depend on the charges" begin
    # the same physics on sites that declare a conserved quantity and on sites that do not
    function chain(site)
        sys = System(4, site)
        p(v) = State{Pure}(sys, v)
        s = p(["Occ", "Emp", "Occ", "Emp"]) + 0.5 * p(["Emp", "Occ", "Occ", "Emp"]) -
            0.3 * p(["Occ", "Occ", "Emp", "Emp"])
        return s / norm(s)
    end
    d, q = chain(Fermion()), chain(Fermion(conserve = N))
    # an absolute tolerance, several of these being exactly zero, as the trace of a dissipator
    both(f) = @test isapprox(f(d), f(q); atol = 1e-12)
    both(s -> expect(s, N(2)))
    both(s -> expect1(s, N))
    both(s -> expect2(s, (N, N)))
    both(s -> expect2(s, (dag(C), C)))
    both(s -> [ expect(s, dag(C)(i) * C(j)) for i in 1:4, j in 1:4 ])
    # an operator changing the charges by no definite amount is measured on the state conserving
    # nothing
    γ, η = C + dag(C), C - dag(C)
    @test abs(expect(q, η(1) * γ(3))) > 0.1
    both(s -> expect(s, η(1) * γ(3)))
    both(s -> expect1(s, γ))
    both(s -> expect2(s, (η, γ)))
    both(s -> last.(probabilities(s, γ(2))))
    # and measure, whose values have the same kinds, an operator of no definite flux being found
    # Hermitian or not on the system conserving nothing
    ms = [N(2), N, (dag(C), C), dag(C)(1) * C(3), η(1) * γ(3), γ, (η, γ), Purity]
    flat(s) = [ x for (_, v) in measure(s, ms) for x in v ]
    both(flat)
    both(s -> flat(mix(s)))
    @test typeof.(last.(measure(d, ms))) == typeof.(last.(measure(q, ms)))
    @test last(measure(q, ms)[5]) isa Float64
    both(s -> norm(s))
    both(s -> entanglement_entropy(s, 2)[1])
    both(s -> collect(entanglement_entropy(s, 2)[2]))

    both(s -> trace(mix(s)))
    both(s -> purity(mix(s)))                                      # trace2 in 1.x
    both(s -> renyi2(mix(s)))
    both(s -> hermiticity(mix(s)))
    both(s -> expect(mix(s), N(2)))
    both(s -> expect2(mix(s), (dag(C), C)))
    both(s -> expect(mix(s), η(1) * γ(3)))
    both(s -> expect2(mix(s), (γ, N)))
    both(s -> trace(partial_trace(mix(s), [1, 3])))
    both(s -> mutual_info_renyi2(mix(s), 1:2))                     # the cut 2 (decision 9.2)
    both(s -> trace(apply(Gate(F)(1), mix(s))))
    both(s -> trace(apply(Dissipator(N)(2), mix(s))))
    both(s -> expect(apply(Reset("Occ")(2), mix(s)), N(2)))
end

@testset "Observables do not depend on the kind of symmetry" begin
    # a state conserving something strongly is measured through its weak form
    strong = TMS.strong
    function chain(site)
        sys = System(4, site)
        p(v) = State{Pure}(sys, v)
        s = p(["Occ", "Emp", "Occ", "Emp"]) + 0.5 * p(["Emp", "Occ", "Occ", "Emp"])
        return mix(s / norm(s))
    end
    w, s = chain(Fermion(conserve = N)), chain(Fermion(conserve = strong(N)))
    both(f) = @test isapprox(f(w), f(s); atol = 1e-12)
    @test trace(s) ≈ 1
    both(trace)
    both(purity)                                                   # trace2 in 1.x
    both(x -> expect(x, N(2)))
    both(x -> expect1(x, N))
    both(x -> expect2(x, (dag(C), C)))
    both(hermiticity)
    # a correlation whose two ends do not conserve the charge is measured all the same, the
    # two shifts it brings cancelling
    both(x -> [ expect(x, dag(C)(i) * C(j)) for i in 1:4, j in 1:4 ])
    both(x -> expect(x, (C - dag(C))(1) * (C + dag(C))(3)))
    # and measure, which prepares its list on the weak system
    ms = [N(2), N, (dag(C), C), dag(C)(1) * C(3), (C - dag(C))(1) * (C + dag(C))(3), Purity, Trace]
    both(x -> [ v for (_, value) in measure(x, ms) for v in value ])
    both(x -> trace(apply(Dissipator(N)(2), x)))
    # the adjoint exchanges ket and bra, which swaps the two charges a strong symmetry keeps
    # apart. It relabels the charges so that the exchange becomes a permutation of zero flux
    @test hermiticity(s) ≈ 1
    @test trace(hermitianize(s)) ≈ 1
    @test flux(dag(s).state) == flux(s.state)
    # on a Hermitian state the adjoint is invisible, so it is checked on one that is not,
    # holding a coherence between sites 1 and 3 that an off diagonal observable reads
    function skew(site)
        sys = System(4, site)
        p(v) = State{Pure}(sys, v)
        c1, c2 = ["Occ", "Emp", "Emp", "Occ"], ["Emp", "Emp", "Occ", "Occ"]
        a, b = p(c1) + 0.7im * p(c2), p(c2) - 0.4 * p(c1)
        return State(mix(a), mix(a).state + 0.6im * mix(b).state)
    end
    ws, ss = skew(Fermion(conserve = N)), skew(Fermion(conserve = strong(N)))
    adjoint_numbers(x) = [ expect(dag(x), dag(C)(1) * C(3)), expect(dag(x), N(3)),
                           inner(dag(x), x) ]
    @test isapprox(adjoint_numbers(ws), adjoint_numbers(ss); atol = 1e-12)
    @test abs(expect(dag(ss), dag(C)(1) * C(3))) > 0.1
    @test expect(dag(ss), dag(C)(1) * C(3)) ≈ conj(expect(ss, dag(C)(3) * C(1)))
    @test hermiticity(ss) < 0.99
    # tracing out part of the sites is refused, what is left being a mixture over several
    # sectors, which a state keeping the two charges apart cannot hold (decision 16)
    @test_throws "weaken it first" partial_trace(s, [1, 3])
    @test expect(partial_trace(weaken(s), [1, 3]), N(1)) ≈ expect(partial_trace(w, [1, 3]), N(1))
    # a number read off a partial trace is measured all the same, through the weak form
    both(x -> renyi2(x, [1, 2]))
    # which is also the cut 2 of 1.x, mutual_info_renyi2(x, 2) (decision 9.2)
    both(x -> mutual_info_renyi2(x, [1, 2]))
end

@testset "Sites where an operator is not defined, and chosen sites" begin
    # expect1 and expect2 act on the sites where their operators are defined, missing elsewhere
    ψ = random_state(Pure, System([Fermion(), Qubit(), Fermion()]), 4)
    x = expect1(ψ, X)
    @test ismissing(x[1]) && ismissing(x[3]) && x[2] ≈ expect(ψ, X(2))
    @test expect1(ψ, N) ≈ [ expect(ψ, N(i)) for i in 1:3 ]
    xx = expect2(ψ, (N, X))
    @test count(ismissing, xx) == 6 && xx[1, 2] ≈ expect(ψ, N(1) * X(2))
    # and so does measure, X being Hermitian where it is defined, and N and X not commuting on
    # the qubit
    v, vv = last.(measure(ψ, [X, (N, X)]))
    @test v isa Vector{Union{Missing, Float64}}
    @test ismissing(v[1]) && ismissing(v[3]) && v[2] ≈ real(x[2])
    @test vv isa Matrix{Union{Missing, ComplexF64}}
    @test count(ismissing, vv) == 6 && vv[1, 2] ≈ xx[1, 2]
    @test_throws "is defined on no site" expect1(State{Pure}(System(2, Fermion()), "Emp"), X)
    @test_throws "is defined on no site" measure(State{Pure}(System(2, Fermion()), "Emp"), X)
    # expect2 on chosen sites, missing elsewhere
    ψf = random_state(Pure, System(5, Fermion()), 4)
    m = expect2(ψf, (dag(C), C); sites = [1, 3, 4])
    @test count(ismissing, m) == 25 - 9
    @test m[3, 1] ≈ expect(ψf, dag(C)(3) * C(1)) && m[4, 4] ≈ expect(ψf, N(4))
    @test_throws "acts on site 9, which the system does not have" expect2(ψf, (N, N); sites = [1, 9])
end

@testset "A partial trace keeps the fermionic signs" begin
    # a fermion traced out is moved past the fermions kept on its right: on
    # (|110⟩ + |011⟩)/√2 the hopping from 3 to 1 crosses the fermion of site 2, and gives -1/2
    # on the reduced state as on the state
    fe = Fermion()
    s3 = System(3, fe)
    ρ = mix(normalize(State{Pure}(s3, ["Occ", "Occ", "Emp"]) + State{Pure}(s3, ["Emp", "Occ", "Occ"])))
    @test expect(partial_trace(ρ, [2]), dag(C)(1) * C(2)) ≈ expect(ρ, dag(C)(1) * C(3)) ≈ -0.5
    # a state of no definite parity takes the sign on a single fermionic operator as well
    ρ = mix(normalize(State{Pure}(s3, ["Emp", "Occ", "Emp"]) + State{Pure}(s3, ["Emp", "Occ", "Occ"])))
    @test expect(partial_trace(ρ, [1, 2]), C(1) + dag(C)(1)) ≈ expect(ρ, C(3) + dag(C)(3)) ≈ -1
    # every product of two operators of the sites kept, sites traced out on both sides of a
    # charged system
    ρ = mix(random_state(State{Pure}(System(5, Fermion(conserve = N)), ["Occ", "Occ", "Emp", "Emp", "Occ"]), 4))
    red = partial_trace(ρ, [1, 4]; keep = true)
    @test red.system == System(2, Fermion(conserve = N))
    for a in (Id, N, C, dag(C)), b in (Id, N, C, dag(C))
        @test expect(red, a(1) * b(2)) ≈ expect(ρ, a(1) * b(4)) atol = 1e-12
    end
    # two species on a site
    ρ = mix(random_state(Pure, System(3, Electron()), 4))
    @test expect(partial_trace(ρ, [2]), dag(Cdn)(1) * Cup(2)) ≈ expect(ρ, dag(Cdn)(1) * Cup(3))
    # no fermion kept on the right of one traced out: no sign, and the bond dimension stays
    ρ = mix(random_state(Pure, System(4, fe), 4))
    @test maxlinkdim(partial_trace(ρ, [3, 4])) ≤ maxlinkdim(ρ)
end

@testset "A partial trace keeps the trace, and what it refuses" begin
    # nothing is divided by the trace, which a state of trace zero keeps
    ρ = mix(random_state(Pure, System(4, Qubit()), 4))
    for k in ([1], [2], [4], [1, 3])
        @test trace(partial_trace(2 * ρ, k)) ≈ 2 * trace(ρ)
    end
    @test abs(trace(partial_trace(ρ - ρ, [1]))) < 1e-12
    @test length(partial_trace(ρ, 1:3; keep = true)) == 3
    # a position the state does not have, every site traced out, and a pure state, which is not
    # mixed behind the caller's back
    @test_throws "given site 5, which the state does not have" partial_trace(ρ, [5])
    @test_throws "given site 0" partial_trace(ρ, [0, 1]; keep = true)
    @test_throws "every site" partial_trace(ρ, 1:4)
    @test_throws "use mix(state) first" partial_trace(random_state(Pure, System(4, Qubit()), 4), [1, 2])
end

@testset "Reduced density matrix and von Neumann entropy" begin
    # tr(ρ_A M) is the expectation value of the operator whose matrix on the kept sites is M,
    # in the order of the sites and the Jordan-Wigner basis of the kept sites alone, which
    # also checks the order of kets and bras with operators that are not hermitian
    LA = TensorMixedStates.LinearAlgebra
    q = random_state(Pure, System(4, Qubit()), 4)
    ρ = reduced_density_matrix(q, [1, 3])
    @test LA.tr(ρ) ≈ 1
    @test LA.ishermitian(round.(ρ; digits = 12))
    for (a, b) in [(X, Z), (Sp, Sm), (Y, Id), (Sm, Sp)]
        m = matrix(a ⊗ b, Qubit(), Qubit())
        @test LA.tr(ρ * m) ≈ expect(q, a(1) * b(3)) atol = 1e-12
    end
    @test reduced_density_matrix(mix(q), [1, 3]) ≈ ρ
    @test reduced_density_matrix(q, Any[3, 1]) ≈ ρ
    s = random_state(Pure, System([Qubit(), Spin(1), Qubit()]), 4)
    ρs = reduced_density_matrix(s, [2, 3])
    @test size(ρs) == (6, 6)
    m = matrix(Sp ⊗ X, Spin(1), Qubit())
    @test LA.tr(ρs * m) ≈ expect(s, Sp(2) * X(3)) atol = 1e-12
    f = random_state(Pure, System(4, Fermion()), 4)
    ρf = reduced_density_matrix(f, [1, 3])
    for (a, b) in [(dag(C), C), (C, C), (N, dag(C))]
        m = matrix(a ⊗ b, Fermion(), Fermion())
        @test LA.tr(ρf * m) ≈ expect(f, a(1) * b(3)) atol = 1e-12
    end
    strong_ρ = State{Mixed}(System(3, Fermion(conserve = strong(N))), ["Occ", "Emp", "Occ"])
    @test real(LA.diag(reduced_density_matrix(strong_ρ, [1, 2]))) ≈ [0, 0, 1, 0]
    @test reduced_density_matrix(q, []) == ones(1, 1)
    @test_throws "given site 5, which the state does not have" reduced_density_matrix(q, [5])
    # the entropy of half a Bell pair is log 2, that of the pair 0, and on a cut it is the
    # entanglement entropy
    bell = apply(controlled(X)(1, 2) * H(1), State{Pure}(System(3, Qubit()), "Up"))
    @test vonneumann_entropy(bell, [1]) ≈ log(2)
    @test abs(vonneumann_entropy(bell, [1, 2])) < 1e-12
    @test vonneumann_entropy(q, 1:2) ≈ first(entanglement_entropy(q, 2))
    @test vonneumann_entropy(q, []) === 0.0
    measured = Dict(measure(q, [ReducedDensityMatrix([1, 3]), VonNeumannEntropy([1])]))
    @test measured["ReducedDensityMatrix(1,3)"] ≈ ρ
    @test measured["VonNeumannEntropy(1)"] ≈ vonneumann_entropy(q, [1])
end

@testset "Reduced density matrix of a pure state" begin
    # contracted on the pure state, it is the one of its mixed representation, fermionic signs
    # included: sites kept on both sides of sites traced out, and a site traced out first
    for (st, positions) in [
            (random_state(Pure, System(6, Qubit()), 8), [[6], [5, 2], [1, 3, 6]]),
            (random_state(State{Pure}(System(6, Fermion(conserve = N)),
                                      ["Occ", "Emp", "Occ", "Emp", "Occ", "Emp"]), 8),
             [[2, 4], [2, 3, 5]]),
            (random_state(State{Pure}(System(4, Electron(conserve = (Ntot, 2Sz))),
                                      ["Up", "Dn", "Emp", "UpDn"]), 6),
             [[1, 3], [2, 4]]),
            (random_state(Pure, System([Qubit(), Fermion(), Boson(3), Fermion(), Qubit()]), 6),
             [[2, 4], [1, 3, 5]])]
        for pos in positions
            @test reduced_density_matrix(st, pos) ≈ reduced_density_matrix(mix(st), pos)
        end
    end
end

@testset "Logarithmic negativity" begin
    # log 2 for a Bell pair, and for the Werner state p Bell + (1 - p) I/4, the depolarization
    # of the pair, log((1 + 3p)/2) above p = 1/3, separable below; the partial transpose of
    # either part gives the same, which lets a part without fermionic site take it
    bell = apply(controlled(X)(1, 2) * H(1), State{Pure}(System(3, Qubit()), "Up"))
    @test log_negativity(bell, [1], [2]) ≈ log(2)
    @test abs(log_negativity(bell, [1], [3])) < 1e-12
    for p in [0.6, 0.2]
        werner = apply(depolarizing_gate(1 - p, 2)(1, 2), mix(bell))
        @test log_negativity(werner, [1], [2]) ≈ max(0, log((1 + 3p) / 2)) atol = 1e-12
    end
    q = random_state(Pure, System(4, Qubit()), 4)
    @test log_negativity(q, [1], [3, 4]) ≈ log_negativity(q, [3, 4], [1])
    s = random_state(Pure, System([Fermion(), Qubit(), Fermion()]), 4)
    @test log_negativity(s, [1], [2]) ≈ log_negativity(s, [2], [1])
    @test_throws "a part without fermionic sites" log_negativity(s, [1], [3])
    @test_throws "sites in both parts" log_negativity(q, [1, 2], [2])
    @test log_negativity(q, [], [2]) == 0
    measured = only(measure(bell, LogNegativity([1], [2])))
    @test first(measured) == "LogNegativity(1;2)"
    @test last(measured) ≈ log(2)
    @test first(only(measure(q, LogNegativity(1:2, [4])))) == "LogNegativity(1,2;4)"
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
    # collapse returns the triple (result, probability, state), where 1.x returned the pair
    # (result, state) (decision 9.1)
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
    # the projector on a level is Proj(Level => x), a number denoting no state (decision 20)
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
    # (result, probability, state), where 1.x returned the pair (result, state) (decision 9.1)
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
    # 1.x refused an operator not commuting with the fermionic parity, now measured, see
    # "Measuring an operator of no definite fermionic parity"
    @test_throws "needs a Hermitian operator" probabilities(plus, Sp(1) * Z(2))
    @test_throws "needs an operator placed on sites" probabilities(plus, Id(1))
    @test_throws "which the system does not have" sample(plus, X(7))
    @test_throws "which the state does not have" collapse(plus, 3)
    # sample of a site checks it, which 1.x did not
    @test_throws "sample was given site 3, which the state does not have" sample(plus, 3)
    @test_throws "sample was given site 0" sample(mix(plus), 0)
end

@testset "Measuring an operator of no definite fermionic parity" begin
    # 1.x refused it; the projectors of the dense path are placed as operators, whose strings the
    # descent lays: on a state superposing the two parities the results of a Majorana operator,
    # ±1, have the probabilities its mean value gives, and the state projected on one of them
    # has it as mean value, in the pure and in the mixed representation
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
    # the result post-selected. collapse returns the triple (result, probability, state), where
    # 1.x returned the pair (probability, state) (decision 9.1)
    plus = State{Pure}(System(2, Qubit()), "+")
    for st in (plus, mix(plus))
        s, p, after = collapse(st, Z(1) * Z(2) => -1)
        @test s ≈ -1
        @test p ≈ 0.5
        @test real(expect(after, Z(1) * Z(2))) ≈ -1
        @test real(expect(after, X(1) * X(2))) ≈ 1
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
    # the level given as a real, where 1.x took a Pair{Int, Int}
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
    # returns the triple (result, probability, state) (decision 9.1)
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
    # collapse returns the triple (result, probability, state) (decision 9.1)
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
    # (result, probability, state) (decision 9.1)
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
end

@testset "What expect refuses" begin
    ψ = State{Pure}(System(2, Qubit()), "Up")
    @test_throws "expect takes an operator placed on sites, such as X(1)" expect(ψ, X)
    @test_throws "measure X(1) with expect" expect1(ψ, X(1))
    @test_throws "measure X(1) with expect" expect1(ψ, [Z, X(1)])
    @test_throws "(X, Y) rather than (X(1), Y(2))" expect2(ψ, (X(1), X(2)))
    @test_throws "(X, Y) rather than (X(1), Y(2))" expect2(ψ, [(X, X), (X(1), X(2))])
    @test_throws "acts on site 3, which the system does not have" expect(ψ, X(3))
    # an operator of several sites, which 1.x refused in expect, is measured
    @test expect(ψ, Swap(1, 2)) ≈ 1
end

@testset "An operator changing the charges by no definite amount" begin
    # measured on the state conserving nothing, where a product of such operators stays a single
    # term, rather than a sum over the charges each factor moves
    sys = System(2, Qubit(conserve = N))
    b = normalize(State{Pure}(sys, ["Up", "Dn"]) + State{Pure}(sys, ["Dn", "Up"]))
    for st in (b, mix(b))
        @test expect(st, X(1) * X(2)) ≈ 1
        @test abs(expect(st, X(1))) < 1e-12
        @test expect2(st, (X, X)) ≈ ones(2, 2)
        # and so does measure, the operator found Hermitian on the system conserving nothing
        xx, x, x1, xs = last.(measure(st, [X(1) * X(2), X(1), X, (X, X)]))
        @test xx isa Float64 && xx ≈ 1
        @test x isa Float64 && abs(x) < 1e-12
        @test x1 isa Vector{Float64} && maximum(abs, x1) < 1e-12
        @test xs isa Matrix{Float64} && xs ≈ ones(2, 2)
    end
    up = State{Pure}(System(20, Qubit(conserve = N)), "Up")
    @test abs(expect(up, prod(X(i) for i in 1:20))) < 1e-12
    # apply keeps refusing it, its result having no definite charges
    @test_throws "no definite flux" apply(X(1), b)
end

@testset "measure and its lists of measurements" begin
    # measure gives the pairs name => value of a list in its order, always a vector, even for one
    # measurement, rather than a bare value whose type would depend on the argument, which
    # decision 7 rejected (c). A list holds operators placed on sites, measured by expect,
    # operators of one site, by expect1, pairs of them, by expect2, the measurements of the
    # library, symbols, the values of their names given to measure, and named functions of the
    # context
    st = State{Pure}(System(3, Qubit()), ["X+", "Y+", "Z+"])
    r = measure(st, [X, (X, Y), Z(1) * Z(2), Purity, "E" => c -> c.time, :energy];
                time = 1.5, energy = -3.2)
    @test r isa Vector{Pair{String, Any}}
    @test first.(r) == ["X", "(X, Y)", "Z(1)*Z(2)", "Purity", "E", "energy"]
    @test last(r[1]) ≈ expect1(st, X)
    @test last(r[2]) ≈ expect2(st, (X, Y))
    @test last(r[3]) ≈ expect(st, Z(1) * Z(2))
    @test last(r[4]) ≈ 1
    @test last(r[5]) == 1.5
    @test last(r[6]) == -3.2
    single = measure(st, X(1))
    @test single isa Vector{Pair{String, Any}}
    @test first(only(single)) == "X(1)" && last(only(single)) ≈ 1
    # a Measure holds such a list, printed by the names of its measurements
    m = Measure(X, (X, Y), Z(1) * Z(2), Purity, "E" => :energy)
    @test repr(m) == "Measure(X, (X, Y), Z(1)*Z(2), Purity, E)"
    @test first.(measure(st, m; energy = 2)) == ["X", "(X, Y)", "Z(1)*Z(2)", "Purity", "E"]
    @test last(last(measure(st, m; energy = 2))) == 2
    # vectors and Measures in a list stand for their measurements, a matrix row by row, so that a
    # list can be built from parts: two Measures make one list, where 1.x gave one vector of pairs
    # for each, and a view is a vector, which 1.x took for a single measurement
    @test first.(measure(st, [[X(1), Z(2)], Measure(Y(1), (Z, Z)), [X(2) Y(2); Z(3) X(3)]])) ==
          ["X(1)", "Z(2)", "Y(1)", "(Z, Z)", "X(2)", "Y(2)", "Z(3)", "X(3)"]
    @test first.(measure(st, [Measure(X(1)), Measure(Z(1), Y(1))])) == ["X(1)", "Z(1)", "Y(1)"]
    @test first.(measure(st, view([X(1), X(2), Z(1)], 1:2))) == ["X(1)", "X(2)"]
    # a name given to any measurement, and to a function of the context, which needs one
    named = measure(st, ["a" => X, "b" => (X, Y), "c" => Purity, "d" => :e, "f" => X(1)]; e = 2)
    @test first.(named) == ["a", "b", "c", "d", "f"]
    @test last(named[1]) ≈ expect1(st, X)
    @test last(named[4]) == 2
    @test first.(measure(st, Measure("a" => X(1), "f" => c -> 1))) == ["a", "f"]
    @test first.(measure(st, [X(1), "f" => c -> 1])) == ["X(1)", "f"]
    # in a vector of named measurements, which Julia types Vector{Pair{String, Any}}, a named
    # function is told by its value, not by the type of its pair
    @test first.(measure(st, ["a" => X(1), "f" => c -> 1])) == ["a", "f"]
    # two measurements of a list cannot share a name (decision 6a)
    @test_throws "several measurements are named \"X(1)\": name them apart" Measure(X(1), X(1))
    @test_throws "named \"X(1)\"" Measure([X(1), [X(1)]])
    @test_throws "named \"X\"" Measure(X, "X" => Z(1))
    @test_throws "named \"X\"" Measure(X, Measure(X))
    @test_ok Measure(X, "X1" => X(1))
end

@testset "What a list of measurements refuses" begin
    st = State{Pure}(System(3, Qubit()), "Up")
    # when it is built (layer 3, question 2): a function needs a name, a constant and a label are
    # no measurements, a Check is returned by a function, a pair is made of two operators of one
    # site, an operator of several sites is placed, and a channel or a generator is no observable
    @test_throws "a measurement written as a function needs a name: \"name\" => c -> …" (
        Measure(c -> 1))
    @test_throws "a constant is not a measurement: write \"name\" => c -> 1.5" Measure(1.5)
    @test_throws "a constant is not a measurement" Measure("a" => 1)
    # a range of numbers as well, which 1.x measured as the vector of its constants
    @test_throws "a constant is not a measurement" measure(st, [1:2])
    @test_throws "\"label\" is not a measurement: a name goes before one" Measure("label")
    @test_throws "\"b\" is not a measurement" Measure("a" => "b")
    @test_throws "a Check is returned by a measurement written as a function" Measure(Check(1, 2))
    @test_throws "a pair measured on every pair of sites is (X, Y)" Measure((X, Y, Z))
    @test_throws "a pair measured on every pair of sites is (X, Y)" Measure((X(1), Y(2)))
    @test_throws "Swap acts on several sites, and is measured placed on them" Measure(Swap)
    @test_throws "measure takes observables, and Gate(X)(1) is not one" Measure(Gate(X)(1))
    @test_throws "measure takes observables, and Gate(X) is not one" Measure(Gate(X))
    @test_throws "measure takes observables, and Dissipator(X)(1) is not one" Measure(Dissipator(X)(1))
    @test_throws "a name goes before one measurement, not before a list" Measure("a" => [X, Y])
    @test_throws "nothing is not a measurement" Measure(nothing)
    @test_throws "is not a measurement" Measure(st)
    # when it is prepared on the system of the state, before anything is measured: a site the
    # system does not have, for an operator and for a measurement of the library, and an operator
    # of one site defined on none
    @test_throws "X(4) acts on site 4, which the system does not have: it has 3 sites" measure(st, X(4))
    @test_throws "SubRenyi2 was given site 5, which the state does not have: it has 3 sites" measure(
        st, SubRenyi2([5]))
    @test_throws "EntanglementEntropy was given site 4" measure(st, EntanglementEntropy(4))
    @test_throws "EntanglementEntropy was given site 0" measure(st, EntanglementEntropy(0, 2))
    @test_throws "LogNegativity was given site 7" measure(st, LogNegativity([1], [7]))
    @test_throws "ReducedDensityMatrix was given site 4" measure(st, ReducedDensityMatrix([1, 4]))
    @test_throws "is defined on no site" measure(st, C)
    # a symbol of the list needs its value, which 1.x left empty (layer 3, question 3), and is
    # refused before anything is measured
    called = Ref(false)
    @test_throws "measure was given no value of energy: measure(state, measurements; energy = …)" (
        measure(st, ["f" => c -> (called[] = true), :energy]))
    @test !called[]
    # and no value can be named state, which c.state reads
    @test_throws "a value given to measure cannot be named state" measure(st, X(1); state = 2)
end

@testset "Functions of the context and Check" begin
    # a measurement of one's own is a named function of the context c, which reads c.state,
    # c.time and the values given to measure by their names (decision 6b), where 1.x had
    # StateFunc and TimeFunc
    st = State{Pure}(System(3, Qubit()), "+")
    r = measure(st, ["x" => c -> real(expect(c.state, X(1))), "t" => c -> 2c.time,
                     "e" => c -> c.energy, "n" => c -> sort!(collect(propertynames(c)))];
                time = 0.5, energy = -1)
    @test first.(r) == ["x", "t", "e", "n"]
    @test last(r[1]) ≈ 1
    @test last(r[2]) == 1
    @test last(r[3]) == -1
    @test last(r[4]) == [:energy, :state, :time]
    # a complex time, as 1.x took
    @test last(only(measure(st, "t" => c -> c.time; time = 0.5im))) == 0.5im
    # what it reads and was not given is refused when it is read: the time, which has no default
    # (layer 3, question 3), and a value
    @test_throws "this measurement reads the time, which measure was not given" measure(
        st, "t" => c -> c.time)
    @test_throws "measure(state, measurements; time = t)" measure(st, "t" => c -> c.time)
    @test_throws "this measurement reads foo, which measure was not given" measure(st, "f" => c -> c.foo)
    # its value is given as it is computed, a complex one with no warning
    ph = apply(Phase(0.7)(1), State{Pure}(System(2, Qubit()), "+"))
    sp = @test_logs last(only(measure(ph, "sp" => c -> expect(c.state, Sp(1)))))
    @test sp ≈ exp(0.7im) / 2

    # a Check is what such a function returns to compare a value with a reference (decision 6b),
    # where 1.x had a measurement Check(name, a, b, tol): measure gives the value, the reference
    # and the norm of their difference
    ch = last(only(measure(st, "chk" => c -> Check(expect(c.state, X(1)), 1; tol = 1e-10))))
    @test keys(ch) == (:value, :reference, :difference)
    @test ch.value ≈ 1 && ch.reference == 1 && ch.difference < 1e-10
    ch = last(only(measure(st, "chk" => c -> Check(expect1(c.state, X), [1, 1, 1]; tol = 1e-10))))
    @test ch.value ≈ ones(3) && ch.reference == [1, 1, 1] && ch.difference < 1e-10
    # complex values are compared whole, each kept as it is computed
    ch = last(only(measure(ph, "chk" => c -> Check(expect(c.state, Sp(1)), exp(0.7im) / 2;
                                                     tol = 1e-10))))
    @test ch.value ≈ exp(0.7im) / 2 && ch.difference < 1e-10
    # beyond tol it is refused, the message giving the name, the time when measure is given one,
    # the values, their difference and tol, or only warned about with warn = true
    far(; kwargs...) = "chk" => c -> Check(1.0, 1.001; tol = 1e-6, kwargs...)
    @test_throws "Check chk failed: 1.0 and 1.001 differ by" measure(st, far())
    @test_throws "Check chk failed at time 2.5: 1.0 and 1.001 differ by" measure(st, far(); time = 2.5)
    @test_throws "beyond tol = 1.0e-6" measure(st, far())
    failed = r"^Check chk failed: 1\.0 and 1\.001 differ by .*, beyond tol = 1\.0e-6$"
    ch = @test_logs (:warn, failed) last(only(measure(st, far(warn = true))))
    @test ch.difference ≈ 0.001
    # without tol nothing is refused, the three being recorded
    ch = @test_logs last(only(measure(st, "chk" => c -> Check(1.0, 2.0))))
    @test ch == (value = 1.0, reference = 2.0, difference = 1.0)
    # values of different shapes are refused, tol or not
    @test_throws "Check chk compares values of different shapes, (2,) and (1, 2)" measure(
        st, "chk" => c -> Check([1.0, 2.0], [1.0 2.0]))
    @test_throws "different shapes, () and (1,)" measure(st, "chk" => c -> Check(1.0, [1.0]; tol = 1))
    # a value not given, which 1.x compared as an empty value, failing on a MethodError: the
    # function reading it is refused, and given the value, the Check compares it
    @test_throws "this measurement reads energy" measure(st, "e" => c -> Check(c.energy, -1.0))
    @test last(only(measure(st, "e" => c -> Check(c.energy, -1.0); energy = -1.0))) ==
          (value = -1.0, reference = -1.0, difference = 0.0)
end

@testset "Real and complex values" begin
    # the value of a Hermitian operator is real, decided numerically on its lowered terms when the
    # list is prepared, that of any other complex, whatever the element type of the state
    # (decision 7): a real value is a Float64 on a complex state, a complex one a ComplexF64 on a
    # real state. 1.x decided by a symbolic test, which took Sp(1)Sm(2) + Sm(1)Sp(2) for complex,
    # wrote an imaginary value under Im(name), and had RealValue, ImaginaryValue and ComplexValue
    # declare a kind, which are gone
    value(st, m) = last(only(measure(st, m)))
    up = State{Pure}(System(2, Qubit()), "Up")
    ph = apply(Phase(0.7)(1), State{Pure}(System(2, Qubit()), "+"))
    plus_i = State{Pure}(System(2, Qubit()), ["+", "i"])
    @test value(ph, X(1)) isa Float64
    @test value(ph, X(1)) ≈ cos(0.7)
    @test value(up, Sp(1)) isa ComplexF64
    @test value(ph, Sp(1)) ≈ exp(0.7im) / 2
    @test value(ph, Sp(1) * Sm(2) + Sm(1) * Sp(2)) isa Float64
    @test value(up, Id(2)) === 1.0
    iv = only(measure(plus_i, im * X(1) * Y(2)))
    @test first(iv) == "im*X(1)*Y(2)"
    @test last(iv) isa ComplexF64
    @test last(iv) ≈ im
    # a non-Hermitian operator stays complex on a real state, whose values are real
    fr = State{Pure}(System(3, Fermion()), ["Occ", "Emp", "Occ"])
    @test value(fr, dag(C)(1) * C(2)) isa ComplexF64
    @test value(fr, C) isa Vector{ComplexF64}
    @test value(fr, (dag(C), C)) isa Matrix{ComplexF64}
    @test value(fr, dag(C)(1) * C(1)) isa Float64
    # what is found Hermitian is: on states with complex amplitudes, the part dropped is rounding,
    # and nothing is warned about; what is not keeps its imaginary part, the whole value of an
    # anti-Hermitian operator
    hop = dag(C)(1) * C(2)
    γγ = (C + dag(C))(1) * (C + dag(C))(2)
    for (st, hermitian, anti, other) in (
            (random_state(Pure, System(3, Qubit()), 4),
             [X(1), Z(1) * Z(2) + X(3), Proj("Dn")(2), exp(0.3X)(1)],
             [im * X(1) * Y(2), (im * Z)(3), im * (X(1) * Y(2) - Y(1) * X(2))],
             [Sp(1), Sp(1) * Z(2)]),
            (random_state(Pure, System(3, Fermion()), 4),
             [hop + dag(C)(2) * C(1), im * (dag(C)(1) * C(3) - dag(C)(3) * C(1)), im * γγ],
             [hop - dag(C)(2) * C(1), γγ],
             [hop]))
        for x in (st, mix(st))
            for o in hermitian
                v = @test_logs value(x, o)
                @test v isa Float64
                @test v ≈ real(expect(x, o))
            end
            for o in [anti; other]
                @test value(x, o) isa ComplexF64
                @test value(x, o) ≈ expect(x, o)
            end
            @test all(o -> abs(real(value(x, o))) < 1e-12, anti)
            @test maximum(o -> abs(imag(value(x, o))), anti) > 1e-3
        end
    end

    # an operator of one site is real when its matrix is Hermitian on every site where it is
    # defined, its vector having missing elsewhere
    @test value(up, X) isa Vector{Float64}
    @test value(up, Sp) isa Vector{ComplexF64}
    @test value(up, Sp + Sm) isa Vector{Float64}
    # a pair is real when its operators are Hermitian on each site, commute there, and are not both
    # odd: the diagonal of (X, Y) is <XY> = i<Z>, which 1.x wrote as zero
    @test value(up, (X, X)) isa Matrix{Float64}
    @test value(up, (Z, N)) isa Matrix{Float64}
    @test value(up, (X, Y)) isa Matrix{ComplexF64}
    @test value(up, (X, Y)) ≈ [im 0; 0 im]
    @test value(up, (X, Z)) isa Matrix{ComplexF64}
    @test value(up, (Sp, Sm)) isa Matrix{ComplexF64}
    # two Majorana operators anticommute on two sites, which makes their product anti-Hermitian,
    # while one of them with an even operator of the other mode of the site gives a Hermitian one
    f = random_state(Pure, System(4, Fermion()), 4)
    γ = C + dag(C)
    g = value(f, (γ, γ))
    @test g isa Matrix{ComplexF64}
    @test maximum(abs(real(g[i, j])) for i in 1:4, j in 1:4 if i ≠ j) < 1e-12
    @test maximum(abs(imag(g[i, j])) for i in 1:4, j in 1:4 if i ≠ j) > 1e-3
    @test value(f, (N, N)) isa Matrix{Float64}
    @test value(f, (dag(C), C)) isa Matrix{ComplexF64}
    e = random_state(Pure, System(3, Electron()), 4)
    @test value(e, (Cup + dag(Cup), Ndn)) isa Matrix{Float64}
    @test value(e, (Cup + dag(Cup), Ndn)) ≈ real(expect2(e, (Cup + dag(Cup), Ndn)))
    @test value(e, (Cup + dag(Cup), Cdn + dag(Cdn))) isa Matrix{ComplexF64}

    # the measurements of the library and the functions give their values as they are computed
    # (decision 7), Trace keeping its imaginary part, which is what it is measured for; but Trace,
    # TraceDeviation, Overlap and ReducedDensityMatrix are complex on any state, a real one
    # included, where a text file would otherwise refuse its columns once the state turns complex
    ρ2 = (1 + 1im) * mix(State{Pure}(System(2, Qubit()), "+"))
    @test value(ρ2, Trace) ≈ 1 + 1im
    @test value(ρ2, TraceDeviation) ≈ -1im
    @test value(ph, Overlap(ph)) isa ComplexF64
    @test value(up, Overlap(up)) isa ComplexF64
    ρup = mix(up)
    @test value(ρup, Trace) isa ComplexF64
    @test value(up, Trace) isa ComplexF64
    @test value(ρup, TraceDeviation) isa ComplexF64
    @test value(ρup, ReducedDensityMatrix([1])) isa Matrix{ComplexF64}
    @test value(ρup, ReducedDensityMatrix([1])) ≈ [1 0; 0 0]
    sp = @test_logs value(ph, "sp" => c -> expect(c.state, Sp(1)))
    @test sp isa ComplexF64
end

@testset "The imaginary part dropped" begin
    # the imaginary part of the value of a Hermitian operator is dropped with a warning when it is
    # beyond rounding next to the scale of the operator, the sum over its terms of |coef| times
    # the operator norms of its blocks, and not next to the value (layer 3, question 4): a value
    # zero up to rounding is not warned about
    up = State{Pure}(System(2, Qubit()), "Up")
    @test (@test_logs measure(up, X(1))) == ["X(1)" => 0.0]
    r = random_state(Pure, System(5, Qubit()), 4)
    for s in (r, mix(r))
        o = X(2) * Z(3) - real(expect(s, X(2) * Z(3))) * Id(1)
        @test abs(expect(s, o)) < 1e-12
        @test_logs measure(s, [o, X, (X, Z), Z(1) * Z(5)])
    end
    # a density matrix that is not Hermitian gives the Hermitian X an imaginary part of 0.6
    ρ = State{Mixed}(System(1, Qubit()), [0.5 0.3im ; 0.3im 0.5])
    @test expect(ρ, X(1)) ≈ 0.6im
    msg = "the imaginary part of a Hermitian measurement is beyond rounding, dropped"
    x = @test_logs (:warn, msg) last(only(measure(ρ, X(1))))
    @test x isa Float64 && x == 0
    # the warning carries the measurement, the imaginary part and the time, in a group of its own,
    # which a run gathers into one warning
    logs, _ = Test.collect_test_logs(() -> measure(ρ, X(1); time = 2.))
    @test only(logs).group === :imaginary_part
    @test only(logs).kwargs[:measurement] == "X(1)"
    @test only(logs).kwargs[:imaginary] ≈ 0.6
    @test only(logs).kwargs[:time] == 2.
    # the measurements a call warns about, one warning for each value
    function warned(f)
        logs, _ = Test.collect_test_logs(f)
        @test all(l -> l.message == msg && l.group === :imaginary_part, logs)
        return [ l.kwargs[:measurement] for l in logs ]
    end
    @test warned(() -> measure(ρ, X)) == ["X"]
    # the same imaginary part next to an operator a hundred million times larger is rounding,
    # the value being nearly 0.6i all the same, and next to one a million times larger it is not
    @test_logs measure(ρ, X(1) + 1e8 * Z(1))
    @test warned(() -> measure(ρ, X(1) + 1e6 * Z(1))) == ["X(1)+1000000Z(1)"]
    # a state that is not normalised, whose trace is not real
    ρ2 = (1 + 1im) * mix(State{Pure}(System(2, Qubit()), "+"))
    @test warned(() -> measure(ρ2, X(1))) == ["X(1)"]
    # each value of an operator of one site and of a pair is checked
    @test warned(() -> measure(ρ2, X)) == ["X", "X"]
    @test warned(() -> measure(ρ2, (X, X))) == fill("(X, X)", 4)
end

@testset "A Measure keeps its preparation" begin
    # a Measure is prepared once for each representation and system it is measured on, and
    # measuring again reuses that preparation, where a list given as it is is prepared at each
    # call (decision 5)
    ψ = random_state(Pure, System(4, Qubit()), 4)
    m = Measure(X, (X, Y), Z(1) * Z(2), Purity)
    @test isempty(m.prepared)
    a = measure(ψ, m)
    p = only(values(m.prepared))
    @test measure(ψ, m) == a
    @test only(values(m.prepared)) === p
    # another state on an equal system shares it, a system being its content (decision 15)
    measure(random_state(Pure, System(4, Qubit()), 2), m)
    @test length(m.prepared) == 1
    # the mixed representation and another system have their own
    measure(mix(ψ), m)
    @test length(m.prepared) == 2
    measure(random_state(Pure, System(5, Qubit()), 2), m)
    @test length(m.prepared) == 3
    # a preparation refused leaves nothing, see "Energy variance"
end

@testset "Site indices out of the system" begin
    # nothing between writing X(10) and contracting its tensor compared that number with the
    # size of the system, and each path an indexed operator can take used to report a
    # BoundsError on an internal vector
    sys = System(4, Qubit())
    st = State{Pure}(sys, "Up")
    stm = mix(st)
    # the three entries, and both bounds
    for op in (X(5), X(0), X(-1))
        @test_throws "does not have" expect(st, op)
        @test_throws "does not have" measure(st, op)
        @test_throws "does not have" make_mpo(st, op)
        @test_throws "does not have" apply(op, st)
    end
    # the message names the factor at fault and the size of the system
    @test_throws "X(5) acts on site 5" expect(st, X(5))
    @test_throws "it has 4 sites" expect(st, X(5))
    # inside a product and inside a sum, the offending factor is the one named
    @test_throws "Z(9) acts on site 9" expect(st, Z(1) * Z(9))
    @test_throws "Z(7) acts on site 7" make_mpo(st, Z(1) * Z(2) + Z(3) * Z(7))
    # an operator of several sites is checked on each of its sites
    @test_throws "Swap(1,9) acts on site 9" apply(Swap(1, 9), st)
    # a time dependent evolver is a vector of terms, each one checked
    @test_throws "acts on site 6" PreMPO(st, [-im * X(1), -im * X(6)])
    # and a mixed representation goes through the same entries
    @test_throws "acts on site 8" expect(stm, X(8))
    @test_throws "acts on site 8" make_mpo(stm, Dissipator(Sm)(8))
    # an operator acting on several sites at once, which 1.x refused, is split on its sites by
    # the descent (decision 2), see operators.jl
    # a matrix of the wrong size for its site is refused by a message naming the operator and
    # the site, on the three paths. It failed on a DimensionMismatch from reshape
    b3 = Operator{1}("B3", [1. 0. 0. ; 0. 1. 0. ; 0. 0. 1.])
    @test_throws "B3 is given by a 3×3 matrix and cannot act on Qubit()" expect(st, b3(1))
    @test_throws "whose dimension is 2" make_mpo(st, b3(1))
    @test_throws "whose dimension is 2" apply(b3(1), st)
    @test_throws "whose dimension is 4" apply(Operator{2}("B9", ones(9, 9))(1, 2), st)
    @test_throws "does not have" renyi2(stm, [2, 7])
    # an operator placed on no site, or a channel, is refused by expect by name rather than by a
    # MethodError about iterating it
    @test_throws "such as Swap(1,2) rather than Swap" expect(st, Swap)
    @test_throws "expect takes an observable, and Gate(X) is a channel" expect(stm, Gate(X))

    # what is inside the system is untouched
    @test_ok expect(st, X(4) * Z(1))
    @test_ok make_mpo(st, sum(Z(i) * Z(i+1) for i in 1:3))
    @test_ok apply(Swap(1, 4), st)

    # a failed measurement used to leave the preobs cache longer than the state, with
    # undefined entries, so a second call on the same state met an UndefRefError rather than
    # the error it deserved
    s2 = State{Pure}(sys, "+")
    @test_throws "does not have" expect(s2, X(9))
    @test length(s2.preobs.left) ≤ length(s2)
    @test_throws "does not have" expect(s2, X(9))
    @test expect(s2, X(3)) ≈ 1
end

@testset "Measuring one state from several threads" begin
    # its caches were filled on demand with no lock: a thread found one sized and not yet
    # filled, and failed on an UndefRefError. Each state below is fresh, its caches empty,
    # pure and mixed, and on a single thread of Julia the tasks still run one after the other
    ψ = random_state(Pure, System(12, Qubit()), 8)
    ops = [Z(i) * Z(j) for i in 1:12 for j in i+1:12]
    ref = [expect(State(ψ, ψ.state), o) for o in ops]
    for s in (State(ψ, ψ.state), mix(State(ψ, ψ.state)))
        got = fetch.([Threads.@spawn expect(s, o) for o in ops])
        @test got ≈ ref
    end
    # and with one Measure, which the first task to need it prepares, under its lock, the
    # others finding that preparation: one for each representation
    m = Measure(ops..., Z, (X, X), Purity)
    ref = measure(State(ψ, ψ.state), Measure(ops..., Z, (X, X), Purity))
    same(a, b) = first.(a) == first.(b) && all(isapprox.(last.(a), last.(b); atol = 1e-10))
    for s in (State(ψ, ψ.state), mix(State(ψ, ψ.state)))
        got = fetch.([Threads.@spawn measure(s, m) for _ in 1:8])
        @test all(g -> same(g, ref), got)
    end
    @test length(m.prepared) == 2
end
