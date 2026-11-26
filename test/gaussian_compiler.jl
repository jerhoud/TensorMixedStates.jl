# The translation of operators into polynomials in Majoranas for the Gaussian backend (compiler of
# src/gaussian/), against dense matrices built here from the matrices of the sites, with the
# Jordan-Wigner strings laid by hand, independently of the lowering of TMS, qubits and spins 1/2
# along their Pauli axes, compared on the even states, the only ones of the backend.

if !isdefined(Main, :dense_majorana)
    include("gaussian_dense.jl")
end
using .LA: I, Diagonal, eigen, Hermitian

const gc_rng = Random.Xoshiro(20261127)

# the site matrix m on site k of sys, after F on the sites before it when m is odd
function gc_embed(sys, k, m)
    p = fermion_parities(sys.sites[k])
    odd = any(m[a, b] != 0 && isodd(p[a] + p[b]) for a in axes(m, 1), b in axes(m, 2))
    mats = [j < k && odd ? TensorMixedStates.matrix(TensorMixedStates.F, s) :
            j == k ? m : Matrix{ComplexF64}(I, dim(s), dim(s)) for (j, s) in enumerate(sys.sites)]
    reduce(kron, mats)
end
gc_op(sys, O, k) = gc_embed(sys, k, TensorMixedStates.matrix(O, sys.sites[k]))

# the kind of modes of a site, as in the compiler
gc_kind(s) = !isnothing(fermion_species(s)) ? :fermion : dim(s) == 2 ? :qubit : nothing

# the unitary from the occupation basis of a qubit along a to the basis of the site: |0⟩ the +1
# eigenvector of σ_a, |1⟩ = σ_b|0⟩, so that σ_b, σ_c, σ_a read X, Y, Z there
function gc_axis_unitary(a)
    a = a / LA.norm(a)
    σb, _ = TensorMixedStates.qubit_majoranas(a)
    E = eigen(Hermitian(TensorMixedStates.pauli(a)))
    v0 = E.vectors[:, argmax(E.values)]
    hcat(v0, σb * v0)
end

# perm such that a vector in the order of the modes (mode 1 the slowest) is x[perm] for x in the
# order of TMS; the modes are taken block by block, a block per kind in the order the kinds appear;
# an Electron has its modes ↑ then ↓, and its index n↑ + 2n↓ in TMS
function gc_to_modes(sys)
    n = length(sys.sites)
    kinds = [gc_kind(s) for s in sys.sites]
    order = [k for b in unique(filter(!isnothing, kinds)) for k in 1:n if kinds[k] == b]
    nm = [kinds[k] == :fermion ? length(fermion_species(sys.sites[k])) : 1 for k in 1:n]
    M = sum(nm[order])
    perm = zeros(Int, 2^M)
    for idx in 0:2^M-1
        bits = dense_bits(idx, M)
        loc = zeros(Int, n)
        pos = 0
        for k in order
            b = bits[pos+1:pos+nm[k]]
            pos += nm[k]
            loc[k] = nm[k] == 1 ? Int(b[1]) : Int(b[1]) + 2 * Int(b[2])
        end
        tms = 0
        for k in 1:n
            tms = tms * dim(sys.sites[k]) + loc[k]
        end
        perm[idx+1] = tms + 1
    end
    perm
end

# the polynomial of op against the matrix O of TMS: I ⊗ O_even + Y ⊗ O_odd with the ancilla, on the
# even states, the only ones of the backend
function gc_matches(sys, op, O; axes = Dict{Int, Symbol}())
    W = reduce(kron, [gc_kind(s) == :qubit ? gc_axis_unitary(TensorMixedStates.axis_vector(get(axes, k, :Z))) :
                      Matrix{ComplexF64}(I, dim(s), dim(s)) for (k, s) in enumerate(sys.sites)])
    perm = gc_to_modes(sys)
    O = (W' * O * W)[perm, perm]
    M = trailing_zeros(size(O, 1))
    P = Diagonal([(-1)^count(dense_bits(k - 1, M)) for k in 1:2^M])
    Oe = (O + P * O * P) / 2
    expected = kron(Matrix{ComplexF64}(I, 2, 2), Oe) + kron(dense_Y, O - Oe)
    even = Diagonal([iseven(count(dense_bits(k - 1, M + 1))) ? 1.0 : 0.0 for k in 1:2^(M + 1)])
    got = dense_poly(TensorMixedStates.majorana_poly(sys, op; axes), M + 1)
    isapprox(got * even, expected * even; atol = 1e-12)
end

@testset "fermionic operators" begin
    sys = System([Fermion(), Fermion(), Fermion()])
    for k in 1:3, O in (C, dag(C), N, TensorMixedStates.F)
        @test gc_matches(sys, O(k), gc_op(sys, O, k))
    end
    c(k) = gc_op(sys, C, k)
    cd(k) = gc_op(sys, dag(C), k)
    @test gc_matches(sys, dag(C)(1) * C(3) + dag(C)(3) * C(1), cd(1) * c(3) + cd(3) * c(1))
    @test gc_matches(sys, C(1) * C(2) + 0.3im * dag(C)(2) * dag(C)(3), c(1) * c(2) + 0.3im * cd(2) * cd(3))
    @test gc_matches(sys, N(1) * N(3) - 2 * C(2) + 0.5 * N(2), gc_op(sys, N, 1) * gc_op(sys, N, 3) - 2 * c(2) + 0.5 * gc_op(sys, N, 2))
    ee = [Electron(), Electron()]
    sys = System(ee)
    for k in 1:2, O in (Cup, Cdn, dag(Cup), dag(Cdn), Nup, Ndn, Ntot, Nupdn, Sz, Sp, Sm, Sx, Sy, TensorMixedStates.F)
        @test gc_matches(sys, O(k), gc_op(sys, O, k))
    end
    u(k) = gc_op(sys, Cup, k)
    d(k) = gc_op(sys, Cdn, k)
    @test gc_matches(sys, dag(Cup)(1) * Cup(2) + dag(Cdn)(2) * Cdn(1) + 4 * Nupdn(1) + Cup(1) * Cdn(1),
                     u(1)' * u(2) + d(2)' * d(1) + 4 * gc_op(sys, Nupdn, 1) + u(1) * d(1))
    @test gc_matches(sys, Sp(1) * Sm(2) + Cdn(2), gc_op(sys, Sp, 1) * gc_op(sys, Sm, 2) + d(2))
    sys = System([Fermion(), Electron(), Fermion()])
    @test gc_matches(sys, dag(C)(1) * Cdn(2) + C(3) * Cup(2) + Nupdn(2) * N(3),
                     gc_op(sys, dag(C), 1) * gc_op(sys, Cdn, 2) + gc_op(sys, C, 3) * gc_op(sys, Cup, 2) +
                     gc_op(sys, Nupdn, 2) * gc_op(sys, N, 3))
    # random sums of products of operators on random sites
    ops = (C, dag(C), N)
    sys = System([Fermion() for _ in 1:4])
    for _ in 1:20
        op, O = 0.0 * N(1), zeros(ComplexF64, 16, 16)
        for _ in 1:3
            ks = rand(gc_rng, 1:4, rand(gc_rng, 1:3))
            Os = rand(gc_rng, ops, length(ks))
            a = randn(gc_rng, ComplexF64)
            op += a * prod(o(k) for (o, k) in zip(Os, ks))
            O += a * prod(gc_op(sys, o, k) for (o, k) in zip(Os, ks))
        end
        @test gc_matches(sys, op, O)
    end
end

# the degree of a polynomial, without the ancilla's Majorana 1 or 2 counted apart
gc_degree(sys, op; axes = Dict{Int, Symbol}()) =
    TensorMixedStates.GaussianKernels.degree(TensorMixedStates.majorana_poly(sys, op; axes))

@testset "qubits and spins 1/2" begin
    sys = System([Qubit() for _ in 1:4])
    q(O, k) = gc_embed(sys, k, TensorMixedStates.matrix(O, sys.sites[k]))
    for k in 1:4, O in (X, Y, Z, N, Sp, Sm)
        @test gc_matches(sys, O(k), q(O, k))
    end
    # quadratic: Z, neighbouring XX and YX, X on the first qubit (ancilla) and on the last (even
    # states); X on the second is not
    @test gc_degree(sys, Z(2)) == 2 && gc_degree(sys, X(2) * X(3)) == 2 && gc_degree(sys, Y(3) * X(4)) == 2
    @test gc_degree(sys, X(1)) == 2 && gc_degree(sys, Y(4)) == 2 && gc_degree(sys, X(2)) == 4
    @test gc_degree(sys, X(1) * X(3)) > 2
    # axes X, Y and Z, random sums of products of Paulis
    axes = Dict(1 => :X, 2 => :Y, 4 => :X)
    for _ in 1:10
        op, O = 0.0 * Z(1), zeros(ComplexF64, 16, 16)
        for _ in 1:3
            ks = rand(gc_rng, 1:4, rand(gc_rng, 1:3))
            Os = rand(gc_rng, (X, Y, Z), length(ks))
            a = randn(gc_rng, ComplexF64)
            op += a * prod(o(k) for (o, k) in zip(Os, ks))
            O += a * prod(q(o, k) for (o, k) in zip(Os, ks))
        end
        @test gc_matches(sys, op, O; axes)
    end
    # Ising along X: the field is a local parity, ZZ between neighbours a quadratic, XX quartic
    xaxes = Dict(k => :X for k in 1:4)
    @test gc_degree(sys, X(1) + Z(2) * Z(3); axes = xaxes) == 2 && gc_degree(sys, X(2) * X(3); axes = xaxes) == 4
    @test_throws ArgumentError TensorMixedStates.majorana_poly(sys, X(1); axes = Dict(1 => :W))
    sys = System([Spin(1/2), Qubit(), Spin(1/2)])
    for k in (1, 3), O in (Sx, Sy, Sz, Sp, Sm)
        @test gc_matches(sys, O(k), gc_op(sys, O, k))
    end
end

@testset "fermions and qubits mixed" begin
    for sites in ([Fermion(), Qubit(), Fermion(), Qubit()], [Qubit(), Fermion(), Qubit(), Fermion()],
                  [Qubit(), Electron()])
        sys = System(sites)
        fs = [k for (k, s) in enumerate(sites) if gc_kind(s) == :fermion]
        qs = [k for (k, s) in enumerate(sites) if gc_kind(s) == :qubit]
        fops = sites[fs[1]] isa Electron ? (Cup, Cdn, Nup) : (C, dag(C), N)
        D = prod(dim.(sites))
        for _ in 1:10
            op, O = 0.0 * Z(qs[1]), zeros(ComplexF64, D, D)
            for _ in 1:3
                ks = vcat(rand(gc_rng, fs, rand(gc_rng, 0:2)), rand(gc_rng, qs, rand(gc_rng, 0:2)))
                if isempty(ks)
                    continue
                end
                Os = [k in fs ? rand(gc_rng, fops) : rand(gc_rng, (X, Y, Z)) for k in ks]
                a = randn(gc_rng, ComplexF64)
                op += a * prod(o(k) for (o, k) in zip(Os, ks))
                O += a * prod(gc_op(sys, o, k) for (o, k) in zip(Os, ks))
            end
            @test gc_matches(sys, op, O)
        end
    end
end

@testset "sites that have no mode" begin
    # a string over a qubit is the identity there
    p = TensorMixedStates.majorana_poly(System([Fermion(), Qubit(), Fermion()]), dag(C)(1) * C(3))
    @test p == TensorMixedStates.majorana_poly(System([Fermion(), Fermion(), Qubit()]), dag(C)(1) * C(2))
    @test_throws ArgumentError TensorMixedStates.majorana_poly(System([Fermion(), Boson(3)]), C(1) * A(2))
    @test TensorMixedStates.majorana_poly(System([Fermion(), Boson(3)]), N(1)) ==
          TensorMixedStates.majorana_poly(System([Fermion()]), N(1))
    @test_throws ArgumentError TensorMixedStates.MajoranaLayout(System([Tj(), Fermion()]))
end
