# The MPO of an operator: its terms compacted into as few channels as a triangular MPO can have,
# and laid on the sites of a system.
#
# The invariant to prefer is that the MPO stands for the operator it was made from: it gives the
# dense matrix of its terms. Bond dimensions are checked against the least any triangular MPO of
# the operator can have, computed here from the dense operator alone, so that the tests do not
# depend on how the compaction gets there. The testsets of 1.x on compact, its coms and the
# mixed representation are still in v1/compact.jl.

"""
the dense matrix of the operator `op` on the sites of `st`, from its terms
"""
exact(st::State, op) = TMS.system_matrix(Pure, st.system, TMS.lower(st.system, op))

"""
the dense matrix of the MPO `make_mpo` builds for `op` on the sites of `st`, `coefs` being the
values of the time functions of a vector of operators
"""
function laid(st::State, op, coefs = [1.])
    idx = siteinds(Pure, st.system)
    t = prod(make_mpo(PreMPO(st, op), coefs))
    a = Array(IT.dense(t), [ i' for i in reverse(idx) ]..., [ IT.dag(i) for i in reverse(idx) ]...)
    d = prod(IT.dim.(idx))
    return reshape(a, d, d)
end

"""
the least number of memory channels, those besides the term not yet begun and the term
finished, that a triangular MPO of the operator of dense matrix `h` can have on each link, the
sites having the dimensions `ds`: the rank of `h` across the link once its parts that are the
identity on the left or on the right are taken out
"""
function least_channels(h::AbstractMatrix, ds::Vector{Int})
    LA = TMS.LinearAlgebra
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
the number of memory channels of each link of the MPO PreMPO lays `op` as on `st`
"""
memory_channels(st::State, op) = PreMPO(st, op).linkdims .- 1

"""
check that the MPO of `op` on the sites of `st` stands for `op` and has the least channels a
triangular MPO of `op` can have there
"""
function check_mpo(st::State, op)
    h = exact(st, op)
    least = least_channels(h, [ IT.dim(i) for i in siteinds(Pure, st.system) ])
    @test norm(laid(st, op) - h) < 1e-12 * norm(h)
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
    @test laid(st, a) ≈ exact(st, Z(2) + Z(1) * Z(2))
    check_mpo(st, a)
    sf = State{Pure}(System(4, Fermion()), "Emp")
    r = ((Id + 0.25N)^0.5)(3)
    @test laid(sf, r * r * dag(C)(1) * C(4)) ≈ exact(sf, (Id(3) + 0.25N(3)) * dag(C)(1) * C(4))
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
    @test memory_channels(q, h) == least_channels(exact(q, h), [2, 2, 2]) == [2, 2]
    @test memory_channels(s, h) == least_channels(exact(s, h), [3, 3, 3]) == [3, 3]
end

@testset "A time dependent operator" begin
    n = 5
    st = State{Pure}(System(n, Qubit()), "Up")
    h1 = sum(0.6^(j - i) * (X(i) * X(j) + N(i) * N(j)) for i in 1:n for j in i+1:n)
    h2 = sum(Z(i) * Z(i + 1) for i in 1:n-1) + sum(X(i) for i in 1:n)
    h = exact(st, 0.3 * h1 - 1.2 * h2)
    @test norm(laid(st, [h1, h2], [0.3, -1.2]) - h) < 1e-12 * norm(h)
    pre = PreMPO(st, [h1, h2])
    @test_throws "takes as many time functions" make_mpo(pre, [1.])
    @test_throws "time functions take real values" make_mpo(pre, [1., im])
end

@testset "The constant, the zero operator, the charges and the element type" begin
    st = State{Pure}(System(3, Qubit()), "Up")
    @test iszero(laid(st, X(1) - X(1)))
    @test laid(st, 2 * Id(1) + X(2)) ≈ exact(st, 2 * Id(1) + X(2))
    # an operator whose terms change the charges differently is refused
    sn = State{Pure}(System(4, Fermion(conserve = N)), ["Occ", "Emp", "Occ", "Emp"])
    @test_throws "change the charges by" make_mpo(sn, C(1) + dag(C)(3))
    # a real MPO for a real operator, Y * Y included, a complex one otherwise
    @test eltype(make_mpo(st, Y(1) * Y(2))[1]) == Float64
    @test eltype(make_mpo(st, X(1) + Y(2))[1]) == ComplexF64
end
