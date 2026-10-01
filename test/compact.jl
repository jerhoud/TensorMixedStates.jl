# compact, and the MPO a com is laid as.
#
# Goes here: compact and the coms it builds, their place in the operator algebra, and the MPO
# PreMPO lays an operator as, compacted and reduced on the sites of the system with the
# relations its operators have there. The invariant to prefer is that the MPO stands for the
# operator it was made from: it gives the dense operator an MPO laid term by term gives, a
# channel for each term as PreMPO did before compacting. Bond dimensions are checked against the
# least any triangular MPO of the operator can have, computed here from the dense operator
# alone, so that the tests do not depend on how compact gets there.

"""
the operator `op` laid on the sites of `state` term by term, a channel for each term of several
sites, as PreMPO laid it before compacting: the reference the MPOs are checked against
"""
function naive_pre(state::State{R}, op) where R
    TMS = TensorMixedStates
    n = op isa Vector ? length(op) : 1
    s = TMS.removeMulti(simplify(TMS.adapt_representation(R, op)))
    return TMS.PreMPO!(TMS.PreMPO{R}(state.system, n), s)
end

"""
the dense matrix of `op` on the sites of `state`, laid term by term, or as PreMPO lays it,
`coefs` being the values of the time functions of a vector of operators
"""
exact(state::State, op, coefs = [1.]) = dense(state, make_mpo(naive_pre(state, op), coefs))
laid(state::State, op, coefs = [1.]) = dense(state, make_mpo(PreMPO(state, op), coefs))

"""
the least number of memory channels, those besides the term not yet begun and the term
finished, that a triangular MPO of the operator of dense matrix `h` can have on each link, the
sites having the dimensions `ds`: the rank of `h` across the link once its parts that are the
identity on the left or on the right are taken out
"""
function least_channels(h::AbstractMatrix, ds::Vector{Int})
    LA = TensorMixedStates.LinearAlgebra
    n = length(ds)
    return map(1:n-1) do k
        dl, dr = prod(ds[1:k]), prod(ds[k+1:n])
        m = reshape(permutedims(reshape(h, dr, dl, dr, dl), (2, 4, 1, 3)), dl^2, dr^2)
        vl = vec(Matrix{eltype(m)}(LA.I, dl, dl)) / sqrt(dl)
        vr = vec(Matrix{eltype(m)}(LA.I, dr, dr)) / sqrt(dr)
        m = m - vl * (vl' * m)
        m = m - (m * vr) * vr'
        count(>(1e-10 * norm(h)), LA.svdvals(m))
    end
end

"""
the number of memory channels of each link of the MPO PreMPO lays `op` as on `state`
"""
memory_channels(state::State, op) = PreMPO(state, op).linkdims .- 1

"""
check that PreMPO lays `op`, and `compact(op)`, as MPOs that stand for `op` on the sites of
`state` and have the least channels a triangular MPO of `op` can have there
"""
function check_compact(state::State{R}, op) where R
    h = exact(state, op)
    least = least_channels(h, [ dim(SysIndex{R}(state.system, k)) for k in 1:length(state) ])
    for a in (op, compact(op))
        @test norm(laid(state, a) - h) < 1e-12 * norm(h)
        @test memory_channels(state, a) == least
    end
end

@testset "compact keeps the operator" begin
    rng = Xoshiro(20261001)
    n = 5
    st = State{Pure}(System(n, Qubit()), "Up")
    c = randn(rng, n, n)
    for h in [Z(1) * X(4) + Z(1) * Y(4) + Z(2) * X(4) + Z(3) * Y(4),
              sum(c[i, j] * (X(i) * X(j) + Y(i) * Y(j) + Z(i) * Z(j)) for i in 1:n for j in i+1:n),
              sum(c[i, j] * c[j, k] * X(i) * Y(j) * Z(k) for i in 1:n for j in i+1:n for k in j+1:n),
              sum(0.6^(j - i) * Z(i) * Z(j) for i in 1:n for j in i+1:n) + sum(X(i) for i in 1:n) + 2 * Id(1)]
        hc = compact(h)
        # symbolically too, the relations of the sites being used only where a com is laid
        @test hc ≈ h
        @test norm(laid(st, hc) - exact(st, h)) < 1e-12 * norm(exact(st, h))
        # the same coefficients to the bit, which the fingerprint of a phase hashes
        @test compact(h) == hc
    end
end

@testset "An MPO has the least channels a triangular MPO can have" begin
    rng = Xoshiro(20261001)
    n = 5
    st = State{Pure}(System(n, Qubit()), "Up")
    c = randn(rng, n, n)
    for h in [sum(c[i, j] * Z(i) * Z(j) for i in 1:n for j in i+1:n),
              sum(c[i, j] * (X(i) * X(j) + Y(i) * Y(j) + Z(i) * Z(j)) for i in 1:n for j in i+1:n),
              sum(0.6^(j - i) * Z(i) * Z(j) for i in 1:n for j in i+1:n) + sum(X(i) for i in 1:n),
              sum(Z(i) * Z(j) for i in 1:n for j in i+1:min(n, i + 2)),
              sum(c[i, j] * c[j, k] * X(i) * Y(j) * Z(k) for i in 1:n for j in i+1:n for k in j+1:n),
              # a term the triangular form cuts where a smaller Schmidt rank would not
              X(1) * Y(2) + X(1) * Y(2) * Z(5)]
        check_compact(st, h)
    end
    # the atoms of a site are compared through their matrices where the com is laid: X*Y is
    # im times Z, Sx is half Sp + Sm, N is (1 - Z) / 2, S2 the identity times 3/4 and Sp*Sp
    # zero on a qubit, which compact alone does not know
    st3 = State{Pure}(System(3, Qubit()), "Up")
    for h in [(X * Y)(1) * Z(3) + Z(1) * X(3),
              Sx(1) * Sx(3) + Sp(1) * Sm(3) + Sm(1) * Sp(3),
              N(1) * N(3) + Z(1) * Z(3),
              S2(1) * Z(3) + Z(1) * Z(3),
              (Sp * Sp)(1) * X(3) + Z(1) * Z(3) + X(1) * Z(2) * (Sp * Sp)(3)]
        check_compact(st3, h)
    end
    @test memory_channels(st3, (X * Y)(1) * Z(3) + Z(1) * X(3)) == [1, 1]
    @test memory_channels(st3, N(1) * N(3) + Z(1) * Z(3)) == [1, 1]
    @test memory_channels(st, sum(0.7^(j - i) * (N(i) * N(j) + Z(i) * Z(j)) for i in 1:n for j in i+1:n)) ==
          fill(1, n - 1)
end

@testset "An MPO of fermions, spins, electrons and bosons" begin
    n = 5
    hf = sum(0.5^(j - i) * (dag(C)(i) * C(j) + dag(C)(j) * C(i)) + 0.3^(j - i) * N(i) * N(j)
             for i in 1:n for j in i+1:n)
    check_compact(State{Pure}(System(n, Fermion()), "Emp"), hf)
    check_compact(State{Pure}(System(n, Fermion(conserve = N)), [ isodd(i) ? "Occ" : "Emp" for i in 1:n ]), hf)
    m = 4
    hs = sum(0.7^(j - i) * (Sp(i) * Sm(j) + Sm(i) * Sp(j) + 0.5 * Sz(i) * Sz(j)) for i in 1:m for j in i+1:m)
    check_compact(State{Pure}(System(m, Spin(1, conserve = 2Sz)), [ isodd(i) ? "1" : "-1" for i in 1:m ]), hs)
    k = 3
    he = sum(0.5^(j - i) * (dag(Cup)(i) * Cup(j) + dag(Cup)(j) * Cup(i) + dag(Cdn)(i) * Cdn(j) + dag(Cdn)(j) * Cdn(i)) +
             0.2^(j - i) * Ntot(i) * Ntot(j) for i in 1:k for j in i+1:k) + sum(0.3 * Nup(i) * Ndn(i) for i in 1:k)
    check_compact(State{Pure}(System(k, Electron()), "Emp"), he)
    hb = sum(0.5^(j - i) * (dag(A)(i) * A(j) + dag(A)(j) * A(i)) + N(i) * N(j) for i in 1:k for j in i+1:k)
    check_compact(State{Pure}(System(k, Boson(3)), "0"), hb)
end

@testset "A com is laid on the sites it meets" begin
    # compact knows no system, and the same com is reduced with the matrices of each one it is
    # laid on: Sp*Sp vanishes on a qubit and not on a spin 1
    h = (Sp * Sp)(1) * Sx(3) + Sz(1) * Sz(3) + Sp(1) * Sm(3)
    hc = compact(h)
    q, s = State{Pure}(System(3, Qubit()), "Up"), State{Pure}(System(3, Spin(1)), "1")
    for st in (q, s)
        @test norm(laid(st, hc) - exact(st, h)) < 1e-12 * norm(exact(st, h))
    end
    @test memory_channels(q, hc) == least_channels(exact(q, h), [2, 2, 2]) == [2, 2]
    @test memory_channels(s, hc) == least_channels(exact(s, h), [3, 3, 3]) == [3, 3]
end

@testset "An MPO on a system of several kinds of sites" begin
    # each site reduces the com with its own matrices: Sp*Sp vanishes on the qubits and not on
    # the spins 1, and the string of C on the boson is its identity
    sys = System([Qubit(), Spin(1), Qubit(), Spin(1)])
    h = sum(0.6^(j - i) * (Sp(i) * Sm(j) + Sm(i) * Sp(j) + (Sp * Sp)(i) * Sz(j) + Sz(i) * Sz(j))
            for i in 1:4 for j in i+1:4)
    check_compact(State{Pure}(sys, ["Up", "1", "Up", "1"]), h)
    sys = System([Fermion(), Boson(3), Fermion()])
    h = dag(C)(1) * C(3) + dag(C)(3) * C(1) + 0.5 * N(1) * N(2) + 0.3 * (A + dag(A))(2) * N(3)
    check_compact(State{Pure}(sys, ["Occ", "1", "Emp"]), h)
end

@testset "An MPO in the mixed representation" begin
    m = 3
    ρ = State{Mixed}(System(m, Qubit()), "FullyMixed")
    h = sum(0.6^(j - i) * (X(i) * X(j) + N(i) * N(j)) for i in 1:m for j in i+1:m)
    dissipators = sum(Dissipator(sqrt(0.1) * Sm)(i) for i in 1:m)
    ev = exact(ρ, -im * h + dissipators)
    @test norm(laid(ρ, -im * compact(h) + dissipators) - ev) < 1e-12 * norm(ev)
    check_compact(ρ, -im * Evolver(h) + dissipators + Dissipator(X ⊗ X)(1, 3))
end

@testset "A compacted time dependent evolver" begin
    n = 5
    st = State{Pure}(System(n, Qubit()), "Up")
    h1 = sum(0.6^(j - i) * (X(i) * X(j) + N(i) * N(j)) for i in 1:n for j in i+1:n)
    h2 = sum(Z(i) * Z(i + 1) for i in 1:n-1) + sum(X(i) for i in 1:n)
    h = exact(st, 0.3 * h1 - 1.2 * h2)
    for a in ([h1, h2], compact([h1, h2]))
        @test norm(laid(st, a, [0.3, -1.2]) - h) < 1e-12 * norm(h)
    end
end

@testset "Measuring a com" begin
    n = 5
    h = sum(0.6^(j - i) * (X(i) * X(j) + N(i) * N(j)) + 0.3^(j - i) * Z(i) * Z(j) for i in 1:n for j in i+1:n) +
        0.7 * sum(Z(i) for i in 1:n)
    hc = compact(h)
    ψ = RandomState{Pure}(System(n, Qubit()), 4)
    # expect measures term by term what it is given, measure the compacted operator
    @test isapprox(expect(ψ, hc), expect(ψ, h); rtol = 1e-12)
    @test isapprox(last(only(measure(ψ, h))), expect(ψ, h); rtol = 1e-12)
    @test isapprox(variance(h, ψ), variance(make_mpo(naive_pre(ψ, h)), ψ); rtol = 1e-10, atol = 1e-12)
    ρ = RandomState{Mixed}(System(n, Qubit()), 4)
    @test isapprox(expect(ρ, hc), expect(ρ, h); rtol = 1e-12)
    hf = sum(0.5^(j - i) * (dag(C)(i) * C(j) + dag(C)(j) * C(i)) + 0.2 * N(i) * N(j) for i in 1:n for j in i+1:n)
    ψf = RandomState{Pure}(System(n, Fermion()), 4)
    @test isapprox(expect(ψf, compact(hf)), expect(ψf, hf); rtol = 1e-12)
    # the kind of the values is found on the operator, not on its coms
    @test last(only(measure(ψf, hf))) isa Real
    @test isapprox(last(only(measure(ψf, hf))), real(expect(ψf, hf)); rtol = 1e-12)
end

@testset "WI and WII of a com" begin
    n = 4
    st = State{Pure}(System(n, Qubit()), "Up")
    # with no relation through the identity, the terms keep their sites, and WI and WII are the
    # operators laid term by term gives
    h = sum(0.6^(j - i) * X(i) * X(j) + 0.3^(j - i) * Z(i) * Z(j) for i in 1:n for j in i+1:n) + 0.7 * sum(Z(i) for i in 1:n)
    for w in (make_approx_W1, make_approx_W2)
        a, b = dense(st, w(PreMPO(st, -im * h), 0.1)), dense(st, w(naive_pre(st, -im * h), 0.1))
        @test norm(a - b) < 1e-12 * norm(b)
    end
    # N is (1 - Z) / 2 where Z also opens terms, and the identity taken out moves part of the
    # operator to terms of one site: the approximations change, and remain of the first order,
    # their error over a step shrinking four times when the step is halved
    h = sum(0.7^(j - i) * (N(i) * N(j) + Z(i) * Z(j)) for i in 1:n for j in i+1:n) + 0.4 * sum(X(i) for i in 1:n)
    hd = exact(st, h)
    for w in (make_approx_W1, make_approx_W2)
        err(τ) = norm(dense(st, w(PreMPO(st, -im * h), τ)) - exp(-im * τ * hd))
        @test 3.5 < err(1e-3) / err(5e-4) < 4.5
    end
end

@testset "Merging coms" begin
    n = 5
    st = State{Pure}(System(n, Qubit()), "Up")
    h1 = sum(0.6^(j - i) * X(i) * X(j) for i in 1:n for j in i+1:n)
    h2 = sum(Z(i) * Z(i + 1) for i in 1:n-1) + sum(X(i) for i in 1:n)
    @test memory_channels(st, compact(compact(h1))) == memory_channels(st, h1)
    @test memory_channels(st, compact(h1) + h2) == memory_channels(st, h1 + h2)
    @test memory_channels(st, compact(h1) + compact(h1)) == memory_channels(st, h1)
    h = exact(st, h1 + h2)
    @test norm(laid(st, compact(h1) + h2) - h) < 1e-12 * norm(h)
end

@testset "What a com refuses" begin
    st = State{Pure}(System(4, Qubit()), "Up")
    hc = compact(sum(0.6^(j - i) * Z(i) * Z(j) for i in 1:4 for j in i+1:4) + sum(X(i) for i in 1:4))
    @test_throws "is compacted and cannot be multiplied" make_mpo(st, hc * X(1))
    @test_throws "is compacted and cannot be multiplied" make_mpo(st, hc^2)
    @test_throws "cannot apply sums as gates" apply(hc, st)
    @test_throws "which is placed on sites" Gate(compact(Z(1) * Z(2)))
    # a com takes its gate where a product by an operator on mixed states asks for it, which
    # refuses it as a product
    ρ = State{Mixed}(System(4, Qubit()), "FullyMixed")
    @test_throws "is compacted and cannot be multiplied" make_mpo(ρ, Gate(X)(1) * hc)
    @test_throws "compact takes an operator placed on sites" compact(X)
    # PreMPO compacts, and keeps its own refusal of a factor of several sites
    m2 = Operator{2}("M2", [1. 0 0 0 ; 0 0 1 0 ; 0 1 0 0 ; 0 0 0 1], plain_op)
    @test_throws "which an MPO cannot place" make_mpo(st, m2(1, 2))
end
