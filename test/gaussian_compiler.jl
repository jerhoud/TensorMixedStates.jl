# The translation of operators into polynomials in Majoranas for the Gaussian backend (compiler of
# src/gaussian/), against dense matrices built here from the matrices of the sites, with the
# Jordan-Wigner strings laid by hand, independently of the lowering of TMS.

if !isdefined(Main, :dense_majorana)
    include("gaussian_dense.jl")
end
using .LA: I, Diagonal

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

# perm such that a vector in the order of the modes (mode 1 the slowest) is x[perm] for x in the
# order of TMS; an Electron has its modes ↑ then ↓, and its index n↑ + 2n↓ in TMS
function gc_to_modes(sys)
    nm = [length(fermion_species(s)) for s in sys.sites]
    M = sum(nm)
    perm = zeros(Int, 2^M)
    for idx in 0:2^M-1
        bits = dense_bits(idx, M)
        tms, pos = 0, 0
        for (k, s) in enumerate(sys.sites)
            b = bits[pos+1:pos+nm[k]]
            pos += nm[k]
            tms = tms * dim(s) + (nm[k] == 1 ? Int(b[1]) : Int(b[1]) + 2 * Int(b[2]))
        end
        perm[idx+1] = tms + 1
    end
    perm
end

# the polynomial of op against the matrix O of TMS: I ⊗ O_even + Y ⊗ O_odd with the ancilla
function gc_matches(sys, op, O)
    perm = gc_to_modes(sys)
    O = O[perm, perm]
    M = trailing_zeros(size(O, 1))
    P = Diagonal([(-1)^count(dense_bits(k - 1, M)) for k in 1:2^M])
    Oe = (O + P * O * P) / 2
    expected = kron(Matrix{ComplexF64}(I, 2, 2), Oe) + kron(dense_Y, O - Oe)
    isapprox(dense_poly(TensorMixedStates.majorana_poly(sys, op), M + 1), expected; atol = 1e-12)
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

@testset "sites that are not fermionic" begin
    # a string over a qubit is the identity there; an operator on it is not accepted yet
    p = TensorMixedStates.majorana_poly(System([Fermion(), Qubit(), Fermion()]), dag(C)(1) * C(3))
    @test p == TensorMixedStates.majorana_poly(System([Fermion(), Fermion()]), dag(C)(1) * C(2))
    @test_throws ArgumentError TensorMixedStates.majorana_poly(System([Fermion(), Qubit()]), C(1) * X(2))
    @test_throws ArgumentError TensorMixedStates.MajoranaLayout(System([Tj(), Fermion()]))
end
