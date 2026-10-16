# The MPO of an operator: its terms compacted into as few channels as a triangular MPO can have,
# and laid on the sites of a system.
#
# The invariant to prefer is that the MPO stands for the operator it was made from: it gives the
# dense matrix of its terms. Bond dimensions are checked against the least any triangular MPO of
# the operator can have, computed here from the dense operator alone, so that the tests do not
# depend on how the compaction gets there. The MPO of channels and generators on a density matrix is
# tested in channels.jl.

"""
the number of memory channels of each link of the MPO PreMPO lays `op` as on `st`
"""
memory_channels(st::State, op) = PreMPO(st, op).linkdims .- 1

"""
check that the MPO of `op` on the sites of `st` stands for `op` and has the least channels a
triangular MPO of `op` can have there
"""
function check_mpo(st::State, op)
    h = lowered(st.system, op)
    least = least_channels(h, [ IT.dim(i) for i in siteinds(Pure, st.system) ])
    @test norm(mpo_matrix(st, op) - h) < 1e-12 * norm(h)
    @test memory_channels(st, op) == least
end

@testset "An MPO has the least channels a triangular MPO can have" begin
    rng = Xoshiro(20261001)
    n = 5
    st = State{Pure}(System(n, Qubit()), "Up")
    c = randn(rng, n, n)
    for h in [sum(c[i, j] * Z(i) * Z(j) for i in 1:n for j in i+1:n),
              sum(c[i, j] * (X(i) * X(j) + Y(i) * Y(j) + Z(i) * Z(j)) for i in 1:n for j in i+1:n),
              sum(0.6^(j - i) * Z(i) * Z(j) for i in 1:n for j in i+1:n) + sum(X(i) for i in 1:n),
              sum(0.6^(j - i) * Z(i) * Z(j) for i in 1:n for j in i+1:n) + 2 * Id(1),
              sum(Z(i) * Z(j) for i in 1:n for j in i+1:min(n, i + 2)),
              sum(c[i, j] * c[j, k] * X(i) * Y(j) * Z(k) for i in 1:n for j in i+1:n for k in j+1:n),
              # a term the triangular form cuts where a smaller Schmidt rank would not
              X(1) * Y(2) + X(1) * Y(2) * Z(5),
              sum(Swap(i, i + 1) for i in 1:n-1)]
        check_mpo(st, h)
    end
    # the blocks of a site are compared through their matrices: X*Y is im times Z, Sx is half
    # Sp + Sm, N is (1 - Z) / 2, S2 the identity times 3/4 and Sp*Sp zero on a qubit, and the
    # identity is taken out of each block by its trace
    st3 = State{Pure}(System(3, Qubit()), "Up")
    for h in [(X * Y)(1) * Z(3) + Z(1) * X(3),
              Sx(1) * Sx(3) + Sp(1) * Sm(3) + Sm(1) * Sp(3),
              N(1) * N(3) + Z(1) * Z(3),
              S2(1) * Z(3) + Z(1) * Z(3),
              (Sp * Sp)(1) * X(3) + Z(1) * Z(3) + X(1) * Z(2) * (Sp * Sp)(3)]
        check_mpo(st3, h)
    end
    @test memory_channels(st3, (X * Y)(1) * Z(3) + Z(1) * X(3)) == [1, 1]
    @test memory_channels(st3, N(1) * N(3) + Z(1) * Z(3)) == [1, 1]
    @test memory_channels(st, sum(0.7^(j - i) * (N(i) * N(j) + Z(i) * Z(j)) for i in 1:n for j in i+1:n)) ==
          fill(1, n - 1)
end

@testset "A factor holding the identity" begin
    # (Id + Z)(1) * Z(2) stands for Z(2) + Z(1) * Z(2), its part of one site laid on that site
    st = State{Pure}(System(3, Qubit()), "Up")
    a = ((Id + Z)^0.5)(1) * ((Id + Z)^0.5)(1) * Z(2)
    @test mpo_matrix(st, a) ≈ lowered(st.system, Z(2) + Z(1) * Z(2))
    check_mpo(st, a)
    sf = State{Pure}(System(4, Fermion()), "Emp")
    r = ((Id + 0.25N)^0.5)(3)
    @test mpo_matrix(sf, r * r * dag(C)(1) * C(4)) ≈ lowered(sf.system, (Id(3) + 0.25N(3)) * dag(C)(1) * C(4))
end

@testset "An MPO of fermions, spins, electrons and bosons" begin
    n = 5
    hf = sum(0.5^(j - i) * (dag(C)(i) * C(j) + dag(C)(j) * C(i)) + 0.3^(j - i) * N(i) * N(j)
             for i in 1:n for j in i+1:n)
    check_mpo(State{Pure}(System(n, Fermion()), "Emp"), hf)
    check_mpo(State{Pure}(System(n, Fermion(conserve = N)), [ isodd(i) ? "Occ" : "Emp" for i in 1:n ]), hf)
    # an odd operator, its strings in a channel of their own
    check_mpo(State{Pure}(System(n, Fermion()), "Emp"), C(3) + 0.5 * dag(C)(5))
    m = 4
    hs = sum(0.7^(j - i) * (Sp(i) * Sm(j) + Sm(i) * Sp(j) + 0.5 * Sz(i) * Sz(j)) for i in 1:m for j in i+1:m)
    check_mpo(State{Pure}(System(m, Spin(1, conserve = 2Sz)), [ isodd(i) ? "1" : "-1" for i in 1:m ]), hs)
    k = 3
    he = sum(0.5^(j - i) * (dag(Cup)(i) * Cup(j) + dag(Cup)(j) * Cup(i) + dag(Cdn)(i) * Cdn(j) + dag(Cdn)(j) * Cdn(i)) +
             0.2^(j - i) * Ntot(i) * Ntot(j) for i in 1:k for j in i+1:k) + sum(0.3 * Nup(i) * Ndn(i) for i in 1:k)
    check_mpo(State{Pure}(System(k, Electron()), "Emp"), he)
    hb = sum(0.5^(j - i) * (dag(A)(i) * A(j) + dag(A)(j) * A(i)) + N(i) * N(j) for i in 1:k for j in i+1:k)
    check_mpo(State{Pure}(System(k, Boson(3)), "0"), hb)
end

@testset "An MPO on a system of several kinds of sites" begin
    # each site has its own matrices: Sp*Sp vanishes on the qubits and not on the spins 1, and
    # the string of C on the boson is its identity
    sys = System([Qubit(), Spin(1), Qubit(), Spin(1)])
    h = sum(0.6^(j - i) * (Sp(i) * Sm(j) + Sm(i) * Sp(j) + (Sp * Sp)(i) * Sz(j) + Sz(i) * Sz(j))
            for i in 1:4 for j in i+1:4)
    check_mpo(State{Pure}(sys, ["Up", "1", "Up", "1"]), h)
    sys = System([Fermion(), Boson(3), Fermion()])
    h = dag(C)(1) * C(3) + dag(C)(3) * C(1) + 0.5 * N(1) * N(2) + 0.3 * (A + dag(A))(2) * N(3)
    check_mpo(State{Pure}(sys, ["Occ", "1", "Emp"]), h)
    # one operator on sites of two types: Sp*Sp vanishes on a qubit and not on a spin 1
    h = (Sp * Sp)(1) * Sx(3) + Sz(1) * Sz(3) + Sp(1) * Sm(3)
    q, s = State{Pure}(System(3, Qubit()), "Up"), State{Pure}(System(3, Spin(1)), "1")
    @test memory_channels(q, h) == least_channels(lowered(q.system, h), [2, 2, 2]) == [2, 2]
    @test memory_channels(s, h) == least_channels(lowered(s.system, h), [3, 3, 3]) == [3, 3]
end

@testset "Strings carried by the channels" begin
    # a term keeps its strings apart down to the compaction, which lays them on its blocks only:
    # a channel crossing a site without a block carries them there, F on a fermion and the
    # identity on a qubit
    n = 6
    sf = State{Pure}(System(n, Fermion()), "Emp")
    @test all(t -> length(t.blocks) == 2, TMS.lower(sf.system, dag(C)(1) * C(n) + dag(C)(n) * C(1)))
    check_mpo(sf, sum(dag(C)(i) * C(j) + dag(C)(j) * C(i) for i in 1:n for j in i+1:n))
    check_mpo(sf, C(5) + dag(C)(2) * N(4) * C(6) + C(1) * N(3))
    sys = System([Qubit(), Fermion(), Qubit(), Qubit(), Fermion(), Fermion()])
    sm = State{Pure}(sys, ["Up", "Emp", "Up", "Up", "Emp", "Emp"])
    check_mpo(sm, dag(C)(2) * C(6) + dag(C)(6) * C(2) + X(1) * C(5) + C(6) + X(3) * C(2) * X(4))
end

@testset "A time dependent operator" begin
    n = 5
    st = State{Pure}(System(n, Qubit()), "Up")
    h1 = sum(0.6^(j - i) * (X(i) * X(j) + N(i) * N(j)) for i in 1:n for j in i+1:n)
    h2 = sum(Z(i) * Z(i + 1) for i in 1:n-1) + sum(X(i) for i in 1:n)
    h = lowered(st.system, 0.3 * h1 - 1.2 * h2)
    @test norm(mpo_matrix(PreMPO(st, [h1, h2]), [0.3, -1.2]) - h) < 1e-12 * norm(h)
    pre = PreMPO(st, [h1, h2])
    @test_throws "takes as many time functions" make_mpo(pre, [1.])
    @test_throws "time functions take real values" make_mpo(pre, [1., im])
end

@testset "Measuring a compacted operator" begin
    # expect measures an operator through its compaction, which measure shares, and gives the
    # value of the dense operator on a pure state, a mixed one and fermions
    n = 5
    h = sum(0.6^(j - i) * (X(i) * X(j) + N(i) * N(j)) + 0.3^(j - i) * Z(i) * Z(j)
            for i in 1:n for j in i+1:n) + 0.7 * sum(Z(i) for i in 1:n)
    ψ = random_state(Pure, System(n, Qubit()), 4)
    v = dense_vector(ψ)
    @test isapprox(expect(ψ, h), v' * lowered(ψ.system, h) * v; rtol = 1e-12)
    @test isapprox(last(only(measure(ψ, h))), expect(ψ, h); rtol = 1e-12)
    @test isapprox(variance(ψ, h), real(v' * lowered(ψ.system, h)^2 * v) - real(expect(ψ, h))^2;
                   rtol = 1e-10, atol = 1e-12)
    LA = TMS.LinearAlgebra
    ρ = random_state(Mixed, System(n, Qubit()), 4)
    r = ket_bra_matrix(ρ)
    @test isapprox(expect(ρ, h), LA.tr(r * lowered(ρ.system, h)) / LA.tr(r); rtol = 1e-12)
    @test isapprox(last(only(measure(ρ, h))), real(expect(ρ, h)); rtol = 1e-12)
    hf = sum(0.5^(j - i) * (dag(C)(i) * C(j) + dag(C)(j) * C(i)) + 0.2 * N(i) * N(j)
             for i in 1:n for j in i+1:n)
    ψf = random_state(Pure, System(n, Fermion()), 4)
    vf = dense_vector(ψf)
    @test isapprox(expect(ψf, hf), vf' * lowered(ψf.system, hf) * vf; rtol = 1e-12)
    # the kind of the value is found on the operator, Hermitian and so real, its difference with
    # its adjoint compacting to nothing
    @test last(only(measure(ψf, hf))) isa Real
    @test last(only(measure(ψf, hf))) ≈ real(expect(ψf, hf)) rtol = 1e-12
end

@testset "The constant, the zero operator, the charges and the element type" begin
    st = State{Pure}(System(3, Qubit()), "Up")
    @test iszero(mpo_matrix(st, X(1) - X(1)))
    @test mpo_matrix(st, 2 * Id(1) + X(2)) ≈ lowered(st.system, 2 * Id(1) + X(2))
    # an operator whose terms change the charges differently is refused
    sn = State{Pure}(System(4, Fermion(conserve = N)), ["Occ", "Emp", "Occ", "Emp"])
    @test_throws "change the charges by" make_mpo(sn, C(1) + dag(C)(3))
    # a real MPO for a real operator, Y * Y included, a complex one otherwise
    @test eltype(make_mpo(st, Y(1) * Y(2))[1]) == Float64
    @test eltype(make_mpo(st, X(1) + Y(2))[1]) == ComplexF64
end

@testset "The algebra of PreMPOs" begin
    # the MPO of a product, built from those of its factors, stands for the product and has the
    # least channels, where developing it would give the product of the numbers of their terms,
    # and so does that of a linear combination; on qubits, on fermions with their strings and on
    # charged fermions, and the adjoint
    q = State{Pure}(System(5, Qubit()), "Up")
    f = State{Pure}(System(4, Fermion()), "Emp")
    cf = State{Pure}(System(4, Fermion(conserve = N)), ["Occ", "Emp", "Occ", "Emp"])
    for (st, ops) in [
            (q, [sum(X(i) * X(i + 1) + Y(i) * Y(i + 1) for i in 1:4) + 0.3 * Z(2),
                 sum(Z(i) * Z(i + 2) for i in 1:3) + X(1), 0.5 * X(3) * Y(5) + Z(1)]),
            (f, [sum(dag(C)(i) * C(i + 1) + dag(C)(i + 1) * C(i) for i in 1:3) + 0.4 * (C + dag(C))(2),
                 C(1) * C(4) + 0.3 * N(3), dag(C)(2) + 0.2 * C(3) * N(1)]),
            (cf, [sum(dag(C)(i) * C(i + 1) + dag(C)(i + 1) * C(i) for i in 1:3), dag(C)(2) * C(4),
                  N(1) + dag(C)(4) * C(1)])]
        pres = [ PreMPO(st, op) for op in ops ]
        ds = mpo_matrix.(pres)
        dims = [ IT.dim(i) for i in siteinds(Pure, st.system) ]
        p = prod(pres)
        @test norm(mpo_matrix(p) - prod(ds)) < 1e-12 * norm(prod(ds))
        @test p.linkdims .- 1 == least_channels(prod(ds), dims)
        @test norm(mpo_matrix(pres[1]') - ds[1]') < 1e-12 * norm(ds[1])
        s = pres[1] + 0.5 * pres[2] - pres[3] / 2
        S = ds[1] + 0.5 * ds[2] - ds[3] / 2
        @test norm(mpo_matrix(s) - S) < 1e-12 * norm(S)
        @test s.linkdims .- 1 == least_channels(S, dims)
        @test norm(mpo_matrix(sum(pres)) - sum(ds)) < 1e-12 * norm(sum(ds))
        @test norm(mpo_matrix(-2im * pres[1]) + 2im * ds[1]) < 1e-12 * norm(ds[1])
        z = pres[1] - pres[1]
        @test all(==(1), z.linkdims) && norm(mpo_matrix(z)) < 1e-14
        # the bounds add, a complex coefficient giving twice its norm for the width
        @test s.norm ≈ pres[1].norm + 0.5 * pres[2].norm + 0.5 * pres[3].norm
        @test s.width ≈ pres[1].width + 0.5 * pres[2].width + 0.5 * pres[3].width
        @test (2im * pres[1]).width ≈ 4 * pres[1].norm
    end
    # what cancels but for rounding opens no channel: an identity pushed onto an opening, in a
    # product that is 0.7 (0.1 + 0.2 - 0.3) (1 - X1) Z2 Z3, a product of pieces, M² = 0 on a
    # spin 1, and the merges of the sweeps, in random products of fermions
    q3 = State{Pure}(System(3, Qubit()), "Up")
    p = PreMPO(q3, (0.1 + 0.2) * X(1) * Z(2) * Z(3) + 0.3 * Z(2) * Z(3)) * PreMPO(q3, 0.7 * X(1) - 0.7 * Id(1))
    @test p.linkdims == [1, 1] && norm(mpo_matrix(p)) < 1e-15
    rot = [cos(1) -sin(1) 0; sin(1) cos(1) 0; 0 0 1] * [cos(2) 0 sin(2); 0 1 0; -sin(2) 0 cos(2)]
    M = Operator{1}("M", rot * [0. 1. 0.; 0. 0. 0.; 0. 0. 0.] * rot')
    s2 = State{Pure}(System(2, Spin(1)), "0")
    p = PreMPO(s2, M(1) * Sz(2)) * PreMPO(s2, M(1) * Sx(2))
    @test p.linkdims == [1] && norm(mpo_matrix(p)) < 1e-15
    rng = Xoshiro(20261016)
    for _ in 1:10
        ops = [ sum(round(2rand(rng) - 1; digits = 3) *
                    prod(rand(rng, [C, dag(C), N])(i) for i in sort(randperm(rng, 4)[1:rand(rng, 1:3)]))
                    for _ in 1:rand(rng, 2:6)) for _ in 1:3 ]
        pres = [ PreMPO(f, op) for op in ops ]
        ds = mpo_matrix.(pres)
        p = prod(pres)
        scale = prod(norm.(ds))
        @test norm(mpo_matrix(p) - prod(ds)) < 1e-12 * scale
        least = norm(prod(ds)) < 1e-12 * scale ? [0, 0, 0] : least_channels(prod(ds), [2, 2, 2, 2])
        @test p.linkdims .- 1 == least
    end
    # superoperators compose
    ρ = mix(q)
    l1, l2 = PreMPO(ρ, -im * X(1) * X(2) + Dissipator(Sm)(3)), PreMPO(ρ, Dissipator(Z)(2) + Dissipator(Sp)(4))
    @test norm(mpo_matrix(l1 * l2) - mpo_matrix(l1) * mpo_matrix(l2)) < 1e-12 * norm(mpo_matrix(l1)) * norm(mpo_matrix(l2))
    @test norm(mpo_matrix(l1 - 2 * l2) - (mpo_matrix(l1) - 2 * mpo_matrix(l2))) < 1e-12 * norm(mpo_matrix(l1))
    # what is refused: two systems, two representations, several time functions
    @test_throws "the PreMPOs of a product are on two systems" PreMPO(q, X(1)) * PreMPO(f, N(1))
    @test_throws "the PreMPOs of a product are of two representations" PreMPO(q, X(1)) * l1
    @test_throws "a PreMPO of several time functions is not multiplied" PreMPO(q, [X(1), Z(2)]) *
                                                                         PreMPO(q, X(1))
    @test_throws "the PreMPOs of a sum are on two systems" PreMPO(q, X(1)) + PreMPO(f, N(1))
    @test_throws "the PreMPOs of a sum are of two representations" PreMPO(q, X(1)) - l1
    @test_throws "a PreMPO of several time functions is not added or scaled" 2 * PreMPO(q, [X(1), Z(2)])
    # the bound of the norm, the sum of the sizes of the terms, multiplies, and that of the width
    # of the spectrum of a product is twice it, which gives a search of excited states its weight:
    # H shifted, so that H² has no level of two
    two = State{Pure}(System(2, Qubit()), "Up")
    h = PreMPO(two, Z(1) + 0.4 * X(1) * X(2) + 0.2 * Z(2) + 0.3 * Id(1))
    @test h.norm ≈ 1.9
    h2 = h * h
    @test h2.norm ≈ 1.9^2 && h2.width ≈ 2 * 1.9^2
    # a power by squaring, the power 0 the identity, a negative one refused
    @test norm(mpo_matrix(h^5) - mpo_matrix(h)^5) < 1e-12 * norm(mpo_matrix(h))^5
    @test mpo_matrix(h^0) ≈ TMS.LinearAlgebra.I(4)
    @test_throws "a PreMPO has no power -1: it has no inverse" h^-1
    E = TMS.LinearAlgebra.eigvals(TMS.LinearAlgebra.Hermitian(mpo_matrix(h2)))
    @test h2.width ≥ E[end] - E[1]
    # H conserves Z1Z2, and the lowest level of H² is not in the sector of "Up": the searches
    # start from "+"
    plus = State{Pure}(two.system, "+")
    algo = Dmrg(nsweeps = 10, tol = 1e-13)
    e0, g = ground_state(h2, plus; algo)
    @test e0 ≈ E[1] atol = 1e-10
    @test first(ground_state(h2, plus; algo, orthogonal_to = [g])) ≈ E[2] atol = 1e-8
end

@testset "Products of sums made of the PreMPOs of their factors" begin
    # PreMPO(state, op) builds a product that would develop into more terms than its factors have
    # from their PreMPOs multiplied, equal factors by squaring: the MPO of the operator, with the
    # least channels, as from its terms; on qubits, on fermions with their strings and on
    # charged fermions
    q = State{Pure}(System(5, Qubit()), "Up")
    mag = sum(Z(i) for i in 1:5)
    ham = sum(X(i) * X(i + 1) for i in 1:4) + sum(0.3 * Z(i) for i in 1:5)
    f = State{Pure}(System(4, Fermion()), "Emp")
    hop = sum(dag(C)(i) * C(i + 1) + dag(C)(i + 1) * C(i) for i in 1:3)
    γ = sum(c * (C + dag(C))(i) for (i, c) in enumerate([0.3, 0.7, 1.1, 0.13]))
    cf = State{Pure}(System(4, Fermion(conserve = N)), ["Occ", "Emp", "Occ", "Emp"])
    ntot = sum(N(i) for i in 1:4)
    for (st, op) in [(q, mag^4), (q, ham^3), (q, (ham^2 + ham)^2), (q, ham * mag * ham),
                     (q, 2im * ham^2 + mag), (f, hop^2 + γ^3), (f, (hop + γ)^3),
                     (cf, ntot^3 + hop^2)]
        @test any(TMS.by_blocks, TMS.addends(TMS.gathered(st.system, op)))
        check_mpo(st, op)
    end
    # a product of factors of one term each, or of two factors of two terms, keeps its terms
    @test !TMS.by_blocks(TMS.gathered(q.system, (X(1) + Z(2)) * (X(3) + Z(4))))
    @test !TMS.by_blocks(TMS.gathered(f.system, dag(C)(1) * C(3) * N(4)))
    # the estimates in Float64: a projector on a code space of 64 stabilizers has 2^64 terms,
    # which an Int wrapped around to 0, the product being developed
    p64 = prod((Id(2i - 1) + Z(2i - 1) * Z(2i)) / 2 for i in 1:64)
    q128 = State{Pure}(System(128, Qubit()), "Up")
    @test TMS.by_blocks(TMS.gathered(q128.system, p64))
    @test expect(q128, p64) ≈ 1
    # on charged fermions: a product that vanishes, (Σ C†)² being 0, does not count its flux
    # beside N(1), and a sum of two fluxes is refused
    ad = sum(dag(C)(i) for i in 1:4)
    @test norm(mpo_matrix(cf, ad^2 * hop + N(1)) - lowered(cf.system, N(1))) < 1e-12
    @test_throws "so it has no definite flux" PreMPO(cf, (ad * hop)^2 + N(1))
end

@testset "An MPO truncated in its triangular form" begin
    # limits cut, across each link, the singular values of the part of the operator acting on both
    # sides, for the Hilbert-Schmidt norm normalized on each site, as those of a bond of a state:
    # those of a power law fall off, cut by the cutoff or by maxdim, the bond dimension of the MPO
    # being two more than its channels; the cutoff of rounding leaves the operator itself, on
    # fermions with their strings and on a superoperator, and each time function is cut apart
    nnorm(m) = norm(m) / sqrt(size(m, 1))
    n = 8
    q = State{Pure}(System(n, Qubit()), "Up")
    h = sum(Z(i) * Z(j) / (j - i)^2 for i in 1:n for j in i+1:n) + sum(0.5 * X(i) for i in 1:n)
    full = mpo_matrix(PreMPO(q, h))
    @test PreMPO(q, h).linkdims .- 1 == [1, 2, 3, 4, 3, 2, 1]
    for (limits, links, error) in [(Limits(cutoff = 1e-14), [1, 2, 3, 4, 3, 2, 1], 1e-14),
                                   (Limits(cutoff = 1e-6), [1, 2, 3, 3, 3, 2, 1], 1e-4),
                                   (Limits(maxdim = 5), [1, 2, 3, 3, 3, 2, 1], 1e-4)]
        p = PreMPO(q, h; limits)
        @test p.linkdims .- 1 == links
        @test nnorm(mpo_matrix(p) - full) < error * nnorm(full)
    end
    f = State{Pure}(System(5, Fermion()), "Emp")
    hf = sum((dag(C)(i) * C(j) + dag(C)(j) * C(i)) / (j - i) for i in 1:5 for j in i+1:5) +
         0.3 * (C + dag(C))(3)
    @test nnorm(mpo_matrix(PreMPO(f, hf; limits = Limits(cutoff = 1e-24))) - mpo_matrix(PreMPO(f, hf))) < 1e-13
    ρ = State{Mixed}(System(4, Qubit()), "Up")
    L = -im * sum(Z(i) * Z(j) / (j - i)^2 for i in 1:4 for j in i+1:4) + Dissipator(sum(Sm(i) for i in 1:4))
    @test nnorm(mpo_matrix(PreMPO(ρ, L; limits = Limits(cutoff = 1e-24))) - mpo_matrix(PreMPO(ρ, L))) < 1e-13
    pv = PreMPO(q, [h, sum(X(i) * X(j) / (j - i)^3 for i in 1:n for j in i+1:n)])
    tv = PreMPO(q, [h, sum(X(i) * X(j) / (j - i)^3 for i in 1:n for j in i+1:n)];
                limits = Limits(cutoff = 1e-6))
    @test sum(tv.linkdims) < sum(pv.linkdims)
    for c in ([1., 0.], [0., 1.], [0.3, 0.7])
        @test nnorm(mpo_matrix(tv, c) - mpo_matrix(pv, c)) < 1e-4 * nnorm(mpo_matrix(pv, c))
        @test nnorm(mpo_matrix(TMS.truncated(pv, Limits(cutoff = 1e-6)), c) - mpo_matrix(tv, c)) < 1e-12 * nnorm(mpo_matrix(pv, c))
    end
end
