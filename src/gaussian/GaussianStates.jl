# The Gaussian backend: states of free fermions and qubits held exactly as fermionic Gaussian states,
# by the covariance of the Majoranas of the backend (the ancilla first, see compiler/Majoranas.jl).
# A state of one component, whose phase is not observable, keeps its covariance only.

"""
    GaussianState(state::State{Pure}; axes = Dict(), limits = Limits())

the product state `state` held exactly as a fermionic Gaussian state; `axes` maps a qubit (a site of
dimension 2 that is not fermionic) to the Pauli `:X`, `:Y` or `:Z` along which it is occupied, `:Z`
by default; `limits` is without use. A state that is not a product, an MPS of link dimension above
1, is refused.

    CreateState(…; backend = GaussianState)
"""
struct GaussianState{R <: PM} <: AbstractState{R}
    system::System
    axes::Dict{Int, Symbol}
    layout::MajoranaLayout
    cov::Matrix{Float64}
    # whether the physical state has a definite parity, see the entropies
    definite::Bool
end

function GaussianState(state::State{Pure}; axes = Dict{Int, Symbol}(), limits::Limits = Limits())
    d = maxlinkdim(state.state)
    if d > 1
        throw(ArgumentError("the Gaussian backend starts from a product state, not from an MPS of " *
                            "link dimension $d"))
    end
    system = state.system
    axes = Dict{Int, Symbol}(axes)
    layout = MajoranaLayout(system; axes)
    vectors = local_vectors(state)
    cov, definite = product_covariance(system, layout, vectors)
    GaussianState{Pure}(system, axes, layout, cov, definite)
end

GaussianState(state::State{R}; kwargs...) where R =
    throw(ArgumentError("the Gaussian backend holds pure states only for now, not a State{$R}"))

# the normalised vector of each site of a product state
function local_vectors(state::State)
    ψ = state.state
    n = length(ψ)
    vectors = Vector{Vector{ComplexF64}}(undef, n)
    for j in 1:n
        t = ψ[j]
        for k in (j - 1, j + 1)
            if 1 <= k <= n
                l = commonind(ψ[j], ψ[k])
                t *= onehot(dag(l) => 1)
            end
        end
        v = Vector{ComplexF64}(Array(t, siteind(ψ, j)))
        vectors[j] = v / norm(v)
    end
    vectors
end

# the covariance of Ψ = |0⟩ψ_even + i|1⟩ψ_odd, ψ the product of the vectors, on the Majoranas of the
# backend, and whether ψ has a definite parity. Ψ = √2 Π χ, χ = |+y⟩ ⊗ ψ the product with the
# ancilla in (|0⟩ + i|1⟩)/√2 and Π = (1 + P)/2 the projector on the even states, so that
# ⟨Ψ|γ_a γ_b|Ψ⟩ = ⟨χ|γ_a γ_b|χ⟩ + ⟨χ|γ_a γ_b P|χ⟩, each a product of values on the units (the
# ancilla, then the sites with modes in the order of the modes): γ_a of unit i is the parities of
# the units before it times its local Majorana.
function product_covariance(system::System, layout::MajoranaLayout, vectors)
    order = sort([k for k in eachindex(system.sites) if !isempty(layout.majoranas[k])]; by = k -> layout.first[k])
    units = [(ComplexF64[1, im] / sqrt(2), [pauli_x, pauli_y], pauli_z)]
    for k in order
        push!(units, (vectors[k], layout.majoranas[k], layout.parity[k]))
    end
    ev(v, m) = dot(v, m * v)
    parities = [ev(v, P) for (v, _, P) in units]
    definite = all(abs(p) == 1 for p in parities[2:end])
    nu = length(units)
    # products of the parities before and after each unit
    before = ones(ComplexF64, nu)
    after = ones(ComplexF64, nu)
    for u in 2:nu
        before[u] = before[u - 1] * parities[u - 1]
    end
    for u in nu-1:-1:1
        after[u] = after[u + 1] * parities[u + 1]
    end
    first = cumsum([0; [length(γs) for (_, γs, _) in units]]) .+ 1
    m = first[end] - 1
    M = zeros(m, m)
    for i in 1:nu
        vi, γi, Pi = units[i]
        for (r, ga) in enumerate(γi), (s, gb) in enumerate(γi)
            if r < s
                v = ev(vi, ga * gb) + before[i] * after[i] * ev(vi, ga * gb * Pi)
                M[first[i] + r - 1, first[i] + s - 1] = real(-im * v)
            end
        end
        between = one(ComplexF64)
        for j in i+1:nu
            vj, γj, Pj = units[j]
            for (r, ga) in enumerate(γi), (s, gb) in enumerate(γj)
                v = ev(vi, ga * Pi) * between * ev(vj, gb) + before[i] * ev(vi, ga) * ev(vj, gb * Pj) * after[j]
                M[first[i] + r - 1, first[j] + s - 1] = real(-im * v)
            end
            between *= parities[j]
        end
    end
    M - transpose(M), definite
end

representation(::GaussianState{R}) where R = R
show(io::IO, s::GaussianState{R}) where R = print(io, "GaussianState{$R}($(s.system))")
ignored_options(::GaussianState) = (:algo, :limits, :mpo_limits, :randomize)
trace(::GaussianState{Pure}) = 1.0
norm(::GaussianState{Pure}) = 1.0

as_representation(::Type{Pure}, state::GaussianState{Pure}; limits::Limits = Limits()) = state
as_representation(::Type{R}, ::GaussianState{Pure}; limits::Limits = Limits()) where R =
    throw(ArgumentError("the Gaussian backend holds pure states only for now, not $R"))
