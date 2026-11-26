# The states of the Gaussian backend (GaussianState) against dense states: product states made by
# TMS, carried to the basis of the modes and lifted with the ancilla, Ψ = |0⟩ψ_even + i|1⟩ψ_odd.

if !isdefined(Main, :dense_majorana)
    include("gaussian_dense.jl")
end
if !isdefined(Main, :gc_to_modes)
    include("gaussian_tms.jl")
end

const gs_rng = Random.Xoshiro(20261128)

# the dense vector of a product state of TMS, in the basis of the modes, lifted with the ancilla
function gs_lifted(sys, ψ; axes = Dict{Int, Symbol}())
    W = reduce(kron, [gc_kind(s) == :qubit ? gc_axis_unitary(TensorMixedStates.axis_vector(get(axes, k, :Z))) :
                      Matrix{ComplexF64}(I, dim(s), dim(s)) for (k, s) in enumerate(sys.sites)])
    x = (W' * ψ)[gc_to_modes(sys)]
    M = trailing_zeros(length(x))
    even = [iseven(count(dense_bits(k - 1, M))) ? x[k] : 0 for k in eachindex(x)]
    vcat(even, im .* (x - even))
end

using .LA: dot, eigvals

gs_random_vector(d) = (v = randn(gs_rng, ComplexF64, d); v / LA.norm(v))

@testset "creation from product states" begin
    for (sites, axes) in (([Fermion(), Electron(), Qubit(), Fermion()], Dict(3 => :X)),
                          ([Qubit(), Qubit(), Qubit(), Qubit()], Dict(1 => :X, 2 => :Y)),
                          ([Qubit(), Fermion(), Qubit()], Dict{Int, Symbol}()),
                          ([Electron()], Dict{Int, Symbol}()))
        sys = System(sites)
        for _ in 1:5
            vs = [gs_random_vector(dim(s)) for s in sites]
            g = TensorMixedStates.GaussianState(State{Pure}(sys, vs); axes)
            Ψ = gs_lifted(sys, reduce(kron, vs); axes)
            @test g.cov ≈ dense_covariance(Ψ) atol = 1e-12
            @test g.cov * transpose(g.cov) ≈ I atol = 1e-12
            @test !g.definite
        end
    end
    # states of definite parity, Fock states and eigenstates of the axes
    sys = System([Fermion(), Electron(), Qubit()])
    g = TensorMixedStates.GaussianState(State{Pure}(sys, ["Occ", "Up", "Dn"]))
    @test g.definite
    @test g.cov ≈ dense_covariance(gs_lifted(sys, reduce(kron, [[0, 1.0], [0, 1.0, 0, 0], [0, 1.0]]))) atol = 1e-12
    # refused: a state that is not a product, a mixed state
    sys = System([Qubit() for _ in 1:4])
    @test_throws ArgumentError TensorMixedStates.GaussianState(random_state(Pure, sys, 2))
    @test_throws ArgumentError TensorMixedStates.GaussianState(State{Mixed}(sys, "Up"))
end

# a random Gaussian state of sys (axes Z): with a definite parity, an even Gaussian of the modes; else
# a Gaussian Ψ of the modes and the ancilla, whose physical state is ψ_even - iψ_odd; the GaussianState
# holding it and its physical vector in the basis of TMS
function gs_random(sys; definite::Bool)
    M = sum(gc_kind(s) == :fermion ? length(fermion_species(s)) : 1 for s in sys.sites)
    if definite
        x = dense_gaussian(zeros(Bool, M), gs_rng)
        Ψ = vcat(x, zeros(ComplexF64, 2^M))
    else
        Ψ = dense_gaussian(zeros(Bool, M + 1), gs_rng)
        x = Ψ[1:2^M] - im * Ψ[2^M+1:end]
    end
    ψ = similar(x)
    ψ[gc_to_modes(sys)] = x
    layout = TensorMixedStates.MajoranaLayout(sys)
    g = TensorMixedStates.GaussianState{Pure}(sys, Dict{Int, Symbol}(), layout, dense_covariance(Ψ), definite)
    g, ψ
end

# the reduced density matrix of the sites kept of the vector ψ (the first site the slowest)
function gs_partial(ψ, dims, kept)
    n = length(dims)
    T = reshape(ψ, reverse(dims)...)
    rows = reverse([n + 1 - k for k in kept])
    A = reshape(permutedims(T, vcat(rows, setdiff(1:n, rows))), prod(dims[kept]), :)
    A * A'
end

gs_entropy(ρ) = (λ = max.(eigvals(Hermitian(ρ)), 0); -sum(x * log(x) for x in λ if x > 0; init = 0.0))

@testset "values and entropies" begin
    for (sites, ops) in (([Fermion() for _ in 1:5], (C, dag(C), N)),
                         ([Qubit() for _ in 1:5], (X, Y, Z)),
                         ([Fermion(), Qubit(), Fermion(), Qubit()], nothing),
                         ([Electron(), Fermion(), Electron()], nothing)),
        definite in (true, false)
        sys = System(sites)
        g, ψ = gs_random(sys; definite)
        n = length(sites)
        dims = [dim(s) for s in sites]
        opsof(k) = gc_kind(sites[k]) == :qubit ? (X, Y, Z) : sites[k] isa Electron ? (Cup, Cdn, Nup) : (C, dag(C), N)
        for _ in 1:5
            ks = rand(gs_rng, 1:n, rand(gs_rng, 1:3))
            Os = [rand(gs_rng, opsof(k)) for k in ks]
            op = prod(o(k) for (o, k) in zip(Os, ks))
            O = prod(gc_op(sys, o, k) for (o, k) in zip(Os, ks))
            @test TensorMixedStates.expect_value(g, op) ≈ dot(ψ, O * ψ) atol = 1e-12
        end
        for kept in ([1], [1, 2], [2, 3], [1, 3], [n - 1, n], [n])
            ρ = gs_partial(ψ, dims, kept)
            @test reduced_density_matrix(g, kept) ≈ ρ atol = 1e-12
            @test von_neumann_entropy(g, kept) ≈ gs_entropy(ρ) atol = 1e-10
        end
        for cut in 0:n
            @test entanglement_entropy(g, cut) ≈ gs_entropy(gs_partial(ψ, dims, 1:cut)) atol = 1e-10
        end
    end
end
