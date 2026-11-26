# The translation of operators into polynomials in Majoranas for the Gaussian backend (compiler of
# src/gaussian/), against dense matrices built here from the matrices of the sites, with the
# Jordan-Wigner strings laid by hand, independently of the lowering of TMS, qubits and spins 1/2
# along their Pauli axes, compared on the even states, the only ones of the backend.

if !isdefined(Main, :dense_majorana)
    include("gaussian_dense.jl")
end
using .LA: I, Diagonal, eigen, Hermitian
using TensorMixedStates.GaussianKernels: MajoranaPoly, quadratic_ground_state

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

# the ground energy of H from its quadratic form, H = c + (1/4) Σ k γγ = c + (i/4) Σ h γγ, h = -ik,
# in the even sector of the backend, which holds all the states of the system
function gc_ground_energy(sys, H; axes = Dict{Int, Symbol}())
    c, k, _ = TensorMixedStates.gaussian_generator(sys, H; axes)
    h = -im * k
    @assert all(iszero, imag(h))
    real(c) + quadratic_ground_state(real(h), 1).energy
end

gc_lowest(O) = eigen(Hermitian((O + O') / 2)).values[1]

@testset "generators" begin
    # quadratic_form against the dense matrices
    G = TensorMixedStates.GaussianKernels
    p = MajoranaPoly(Int[] => 0.3, [1, 2] => 1.0 - 2im, [2, 5] => 0.5, [3, 6] => -1im)
    c, k = G.quadratic_form(p, 6)
    γγ(a, b) = dense_matrix(x -> dense_majorana(a, dense_majorana(b, x)), 3)
    @test dense_poly(p, 3) ≈ c * I + sum(k[a, b] * γγ(a, b) for a in 1:6, b in 1:6 if a != b) / 4
    @test_throws ArgumentError G.quadratic_form(MajoranaPoly([1] => 1), 2)
    @test_throws ArgumentError G.quadratic_form(MajoranaPoly([1, 2, 3, 4] => 1), 4)
    # ground energies: fermions with hopping, pairing and a chemical potential
    sys = System([Fermion() for _ in 1:4])
    H = sum(-(dag(C)(j) * C(j + 1) + dag(C)(j + 1) * C(j)) + 0.4 * (C(j) * C(j + 1) + dag(C)(j + 1) * dag(C)(j)) for j in 1:3) +
        0.3 * sum(N(j) for j in 1:4)
    cd(k) = gc_op(sys, dag(C), k)
    c(k) = gc_op(sys, C, k)
    Hd = sum(-(cd(j) * c(j + 1) + cd(j + 1) * c(j)) + 0.4 * (c(j) * c(j + 1) + cd(j + 1) * cd(j)) for j in 1:3) +
         0.3 * sum(gc_op(sys, N, j) for j in 1:4)
    @test gc_ground_energy(sys, H) ≈ gc_lowest(Hd) atol = 1e-12
    # Ising along X with longitudinal fields on both ends (ancilla, even states)
    sys = System([Qubit() for _ in 1:5])
    q(O, k) = gc_embed(sys, k, TensorMixedStates.matrix(O, sys.sites[k]))
    H = -sum(Z(j) * Z(j + 1) for j in 1:4) - 0.7 * sum(X(j) for j in 1:5) + 0.3 * Z(1) - 0.2 * Z(5)
    Hd = -sum(q(Z, j) * q(Z, j + 1) for j in 1:4) - 0.7 * sum(q(X, j) for j in 1:5) + 0.3 * q(Z, 1) - 0.2 * q(Z, 5)
    xaxes = Dict(k => :X for k in 1:5)
    @test gc_ground_energy(sys, H; axes = xaxes) ≈ gc_lowest(Hd) atol = 1e-12
    _, k, _ = TensorMixedStates.gaussian_generator(sys, -im * H; axes = xaxes)
    @test all(iszero, imag(k))
    # electrons: hopping of both spins and a field on the spin
    sys = System([Electron(), Electron()])
    H = -(dag(Cup)(1) * Cup(2) + dag(Cup)(2) * Cup(1) + dag(Cdn)(1) * Cdn(2) + dag(Cdn)(2) * Cdn(1)) + 0.5 * Sz(1)
    e(O, k) = gc_op(sys, O, k)
    Hd = -(e(Cup, 1)' * e(Cup, 2) + e(Cup, 2)' * e(Cup, 1) + e(Cdn, 1)' * e(Cdn, 2) + e(Cdn, 2)' * e(Cdn, 1)) + 0.5 * e(Sz, 1)
    @test gc_ground_energy(sys, H) ≈ gc_lowest(Hd) atol = 1e-12
end

@testset "generators refused" begin
    message(sys, op; axes = Dict{Int, Symbol}()) =
        try
            TensorMixedStates.gaussian_generator(sys, op; axes)
            ""
        catch e
            e isa ArgumentError ? e.msg : rethrow()
        end
    # Ising along Z: ZZ is quartic, the error proposes X
    sys = System([Qubit() for _ in 1:4])
    m = message(sys, -sum(Z(j) * Z(j + 1) for j in 1:3) - 0.7 * sum(X(j) for j in 1:4))
    @test occursin("degree 4", m) && occursin("2 => :X", m) && occursin("3 => :X", m)
    # qubits coupled across another: no axis, reorder the sites
    m = message(System([Qubit() for _ in 1:3]), X(1) * X(3) + Y(1) * Y(3))
    @test occursin("No Pauli axis", m) && occursin("Qubits 1 and 3", m)
    # axes that only the ends allow: XY along X, and X₁Y₂, a chain from the ancilla when the axis of
    # qubit 1 is X (or Z₁Z₂ one when it is Z); a rule per site (couplings exclude their Paulis) finds
    # none on qubit 1 (on 2 qubits every generator is quadratic, so the chain is longer)
    sys = System([Qubit() for _ in 1:4])
    H = sum(Y(j) * Y(j + 1) + Z(j) * Z(j + 1) for j in 1:3) + X(1) * Y(2) + 0.5 * X(3)
    m = message(sys, H)
    proposed = Dict(parse(Int, r[1]) => Symbol(r[2]) for r in eachmatch(r"(\d+) => :([XYZ])", m))
    @test length(proposed) == 4 && gc_degree(sys, H; axes = proposed) == 2
    @test gc_degree(sys, H; axes = Dict(k => :X for k in 1:4)) == 2
    @test message(System([Qubit(), Qubit()]), Y(1) * X(2) + Z(1) * Z(2) + X(1) * Z(2)) == ""
    # fermions and qubits: the qubits come second, the last one is the last mode
    sys = System([Fermion(), Fermion(), Qubit(), Qubit()])
    m = message(sys, dag(C)(1) * C(2) + dag(C)(2) * C(1) + Z(3) * Z(4) + X(3) + X(4))
    @test occursin("make it quadratic", m)
    # Hubbard: fermions only, no advice on axes
    m = message(System([Electron(), Electron()]), Nupdn(1) + dag(Cup)(1) * Cup(2))
    @test occursin("sites 1", m) && occursin("degree 4", m) && !occursin("axes", m)
end

@testset "checks of the compiler" begin
    G = TensorMixedStates.GaussianKernels
    # an axis given for a site that has no qubit is refused, not ignored
    sys = System([Fermion(), Qubit()])
    @test_throws ArgumentError TensorMixedStates.majorana_poly(sys, Z(2); axes = Dict(1 => :X))
    @test_throws ArgumentError TensorMixedStates.majorana_poly(sys, Z(2); axes = Dict(3 => :X))
    # a residue of cancelling terms is kept, and shown in the message
    sys = System([Fermion() for _ in 1:3])
    m = try
        TensorMixedStates.gaussian_generator(sys, 0.1 * N(1) * N(2) + 0.2 * N(1) * N(2) - 0.3 * N(1) * N(2) + N(3))
        ""
    catch e
        e.msg
    end
    @test occursin("degree 4", m) && occursin("e-17", m)
    # hash and quadratic_form
    p = MajoranaPoly([1, 2] => 1.0, Int[] => 2.0)
    @test hash(p) == hash(MajoranaPoly(Int[] => 2.0, [1, 2] => 1.0))
    @test_throws ArgumentError G.quadratic_form(MajoranaPoly([1, 7] => 1), 6)
end
