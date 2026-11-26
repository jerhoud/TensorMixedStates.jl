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
