# Expectation values: expect, through the descent and the compaction, expect1 and expect2 with
# their fast paths, inner(a, op, b), and apply on pure states, checked against dense matrices
# written with explicit Jordan-Wigner strings.
#
# The testsets of 1.x on the other measurements, sampling, partial traces, entropies and the
# channels are still in v1/observables.jl; their assertions on expectation values and on apply
# to pure states are here.

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

@testset "The adjoint of a tensor product of fermions" begin
    # (A ⊗ B)(i, j) is A(i) * B(j), whose adjoint reverses the factors: two fermions
    # anticommute, so dag(C ⊗ C) is the opposite of dag(C) ⊗ dag(C)
    st = random_state(Pure, System(3, Fermion()), 4)
    for (i, j) in [(1, 2), (1, 3), (3, 1)]
        @test expect(st, dag(C ⊗ C)(i, j)) ≈ conj(expect(st, C(i) * C(j)))
        @test expect(st, dag(dag(C) ⊗ C)(i, j)) ≈ conj(expect(st, dag(C)(i) * C(j)))
    end
    # a hopping term written with its adjoint is self adjoint
    h = dag(C) ⊗ C
    @test abs(imag(expect(st, sum((h + dag(h))(i, i + 1) for i in 1:2)))) < 1e-12
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

@testset "Measuring a multiple of the identity" begin
    st = random_state(Pure, System(3, Qubit()), 4)
    for s in (st, mix(st))
        @test expect1(s, 2Id) ≈ fill(2, 3)
        @test expect1(s, (2X)^2) ≈ fill(4, 3)
        @test expect2(s, (2Id, X)) ≈ 2 * expect2(s, (Id, X))
        @test expect(s, 3 * Id(1)) ≈ 3
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

@testset "A random fermionic oracle" begin
    # the Jordan-Wigner strings checked on random products against the dense matrices of the
    # strings written out, through expect on a pure and on a mixed state, the MPO, expect2 and
    # the matrix of a tensor product. In the dense basis C(j) is F ⊗ … ⊗ F ⊗ c ⊗ 1 ⊗ … ⊗ 1, and
    # any other operator is written from those, so that the reference knows nothing of the
    # package
    fe = Fermion()
    mc, mf, mi = matrix(C, fe), matrix(F, fe), matrix(Id, fe)
    cd(j, n) = foldl(kron, [ k < j ? mf : k == j ? mc : mi for k in 1:n ])
    singles = [ (C, (j, n) -> cd(j, n)), (dag(C), (j, n) -> cd(j, n)'),
                (N, (j, n) -> cd(j, n)' * cd(j, n)),
                (C + 2dag(C), (j, n) -> cd(j, n) + 2cd(j, n)') ]
    count = Ref(0)
    # a factor on the sites 1:n, as an operator and as its dense matrix: an operator of one
    # site, or a tensor product of two, renamed or not, on two sites in either order
    function factor(n)
        a, da = rand(singles)
        i = rand(1:n)
        if rand(Bool)
            return a(i), da(i, n)
        end
        b, db = rand(singles)
        j = rand(setdiff(1:n, i))
        t = rand(Bool) ? a ⊗ b : named(a ⊗ b, "T$(count[] += 1)")
        return t(i, j), da(i, n) * db(j, n)
    end
    function product(n)
        fs = [ factor(n) for _ in 1:rand(1:3) ]
        return prod(first.(fs)), prod(last.(fs))
    end
    n = 4
    ψ = random_state(Pure, System(n, fe), 4)
    ρ = mix(ψ)
    v = dense_vector(ψ)
    value(d) = v' * d * v
    w = raw_vector(ψ)
    for _ in 1:30
        o, d = product(n)
        @test expect(ψ, o) ≈ value(d) atol = 1e-10
        @test expect(ρ, o) ≈ value(d) atol = 1e-10
        @test mpo_matrix(ψ, o) ≈ d atol = 1e-10
        @test raw_vector(apply(o, ψ)) ≈ d * w atol = 1e-10
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
    # a product that vanishes gives the null state, and an MPO is applied as it is
    ψ = random_state(Pure, fe, 4)
    @test norm(apply(C(2) * C(2), ψ)) == 0
    mpo = make_mpo(ψ, dag(C)(1) * C(3) + N(2))
    @test norm(apply(mpo, ψ) - apply(dag(C)(1) * C(3) + N(2), ψ)) < 1e-10
    @test_throws "apply_algo is" apply(mpo, ψ; apply_algo = "fit")
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
    both(f) = @test isapprox(f(d), f(q); atol = 1e-12)
    both(s -> expect(s, N(2)))
    both(s -> expect1(s, N))
    both(s -> expect2(s, (N, N)))
    both(s -> expect2(s, (dag(C), C)))
    both(s -> [ expect(s, dag(C)(i) * C(j)) for i in 1:4, j in 1:4 ])
    both(s -> expect(mix(s), N(2)))
    both(s -> expect2(mix(s), (dag(C), C)))
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
    both(x -> expect(x, N(2)))
    both(x -> expect1(x, N))
    both(x -> expect2(x, (dag(C), C)))
    # a correlation whose two ends do not conserve the charge is measured all the same, the
    # two shifts it brings cancelling
    both(x -> [ expect(x, dag(C)(i) * C(j)) for i in 1:4, j in 1:4 ])
end

@testset "Sites where an operator is not defined, and chosen sites" begin
    # expect1 and expect2 act on the sites where their operators are defined, missing elsewhere
    ψ = random_state(Pure, System([Fermion(), Qubit(), Fermion()]), 4)
    x = expect1(ψ, X)
    @test ismissing(x[1]) && ismissing(x[3]) && x[2] ≈ expect(ψ, X(2))
    @test expect1(ψ, N) ≈ [ expect(ψ, N(i)) for i in 1:3 ]
    xx = expect2(ψ, (N, X))
    @test count(ismissing, xx) == 6 && xx[1, 2] ≈ expect(ψ, N(1) * X(2))
    @test_throws "is defined on no site" expect1(State{Pure}(System(2, Fermion()), "Emp"), X)
    # expect2 on chosen sites, missing elsewhere
    ψf = random_state(Pure, System(5, Fermion()), 4)
    m = expect2(ψf, (dag(C), C); sites = [1, 3, 4])
    @test count(ismissing, m) == 25 - 9
    @test m[3, 1] ≈ expect(ψf, dag(C)(3) * C(1)) && m[4, 4] ≈ expect(ψf, N(4))
    @test_throws "acts on site 9, which the system does not have" expect2(ψf, (N, N); sites = [1, 9])
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
