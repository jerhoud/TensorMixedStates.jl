# The Gaussian backend: states of free fermions and qubits held exactly as fermionic Gaussian states,
# by the covariance of the Majoranas of the backend (the ancilla first, see compiler/Majoranas.jl).
# A state of one component, whose phase is not observable, keeps its covariance only.

import .GaussianKernels: GaussianSum, PureGaussian, LogComplex, evolve!, fix!, relink!, norm2, apply_monomial,
                         cross_term, describe, overlap, fock_state, logpfaffian, expect_monomial, fock_density

"""
    GaussianState(state::State{R}; axes = :Z, max_terms = 64, limits = Limits())

the product state `state`, pure, a trajectory or mixed (`R = Pure`, `Trajectory` or `Mixed`), held exactly as a
fermionic Gaussian state, on fermionic sites and sites of dimension 2 (qubits); `limits` is without
use. A state that is not a product, an MPS of link dimension above 1, is refused. A trajectory evolves
under dissipators whose effective operator -iH - (1/2) Σ L†L is quadratic; when it is not unitary, the
trajectory is of one component and its jumps keep it Gaussian: gates exp(A) of quadratic A, projectors on
occupations, linear operators (`C`, `dag(C)`, `Sp`, `Sm`). A mixed state is on each site a Gaussian state
of its modes, which commutes with their parity (`"FullyMixed"`, mixtures of occupations; on a site of
one mode any state that commutes with its parity), of one component; it takes Gaussian unitaries, products
of Majoranas, projectors on occupations, linear operators and channels exp(L), see `apply`, converts to a
`State` through its density matrix, and evolves under a constant Lindbladian, -iH with H quadratic and
dissipators of operators linear in its Majoranas (`C`, `dag(C)` of fermions placed before any qubit, `Sm` and
`Sp` of the first or the last qubit of a system of qubits along `:Z`), exactly, with its `steady_state`, and
under -iH(t) without dissipators as a pure state does; `thermalize` takes it to e^{-βH/2} ρ e^{-βH/2}, exactly.

- `axes`: the Pauli `:X`, `:Y` or `:Z` along which the qubits are occupied, one for all, a vector of
  one per site (`nothing` on the others) or a `Dict` from sites to Paulis
- `max_terms`: the largest number of Gaussian components of the state, which operators that are not
  Gaussian multiply (default 64); a measurement of K components keeps K(K-1)/2 matrices of the size
  of the covariance

    CreateState(…; backend = GaussianState)
    CreateState(GaussianState(State{Pure}(System(4, Qubit()), "X+"); axes = :X))
"""
struct GaussianState{R <: PM} <: AbstractState{R}
    system::System
    axes::Dict{Int, Symbol}
    layout::MajoranaLayout
    # one component: its covariance, its phase free; several: their sum, with their phases
    cov::Union{Nothing, Matrix{Float64}}
    sum::Union{Nothing, GaussianSum}
    # the physical parity ±1 of each component (the occupation of the ancilla), nothing when not known:
    # the sectors of the evolutions that need them, and the entropies of a state of one component
    parities::Union{Nothing, Vector{Int}}
    # the squared norm of the state
    trace::Float64
    max_terms::Int
end

# the system, the axes by site and the layout of the product state `state`, refused when it is not a
# product or has a site without modes
function product_layout(state::State, axes)
    d = maxlinkdim(state.state)
    if d > 1
        throw(ArgumentError("the Gaussian backend starts from a product state, not from an MPS of " *
                            "link dimension $d"))
    end
    system = state.system
    axes = axes_by_site(system, axes)
    layout = MajoranaLayout(system; axes)
    for k in eachindex(system.sites)
        if isempty(layout.majoranas[k])
            throw(ArgumentError("the Gaussian backend holds fermionic sites and sites of dimension 2, " *
                                "not site $k, a $(typeof(system[k]))"))
        end
    end
    system, axes, layout
end

function GaussianState(state::State{R}; axes = Dict{Int, Symbol}(), max_terms::Int = 64,
                       limits::Limits = Limits()) where {R <: Ket}
    system, axes, layout = product_layout(state, axes)
    vectors = local_vectors(state)
    cov, definite = product_covariance(system, layout, vectors)
    parities = definite ? [Int(sign(cov[1, 2]))] : nothing
    GaussianState{R}(system, axes, layout, cov, nothing, parities, norm(state.state)^2, max_terms)
end

# the generator of the references of the sums: fixed, so that a run can be reproduced (the reference
# changes the results by rounding only)
const reference_seed = 20261126

# a state of one component like `state`, of covariance `cov` and physical parity `parities`
function single(state::GaussianState{R}, cov::Matrix{Float64}; parities = state.parities,
                trace::Real = state.trace) where R
    if !(isnothing(parities) || length(parities) == 1)
        throw(ArgumentError("a state of one component has one parity, not $(length(parities))"))
    end
    GaussianState{R}(state.system, state.axes, state.layout, cov, nothing, parities, trace, state.max_terms)
end

# whether the state is of one component and of a known parity
definite(state::GaussianState) = isnothing(state.sum) && !isnothing(state.parities)
# a mixed state of the backend commutes with the parity
definite(::GaussianState{Mixed}) = true

# whether the states are even, the ancilla carrying the parity of the system, so that the compiler may take
# the total parity as 1 (`polynomial`); not a mixed state, its ancilla empty
even_states(::GaussianState) = true
even_states(::GaussianState{Mixed}) = false

# the terms ts as a polynomial on the Majoranas of the state
state_polynomial(state::GaussianState, ts::Vector{Term}) =
    polynomial(state.system, state.layout, ts; even = even_states(state))

# the parity of a component after the monomial α: changed by a single Majorana of the ancilla (1 or 2),
# which changes its occupation
flipped(p::Int, α) = isodd(count(<=(2), α)) ? -p : p

# a state like `state` held by the sum `s`, refused beyond max_terms components, of one component
# when it has one; `trace` its squared norm when known, `parities` those of the components when known
function superposed(state::GaussianState{R}, s::GaussianSum; trace::Union{Nothing, Real} = nothing,
                    parities::Union{Nothing, Vector{Int}} = nothing) where {R <: Ket}
    K = length(s.comps)
    if !(isnothing(parities) || length(parities) == K)
        throw(ArgumentError("$(length(parities)) parities for $K components"))
    end
    if K > state.max_terms
        throw(ArgumentError("the state would have $K Gaussian components, above max_terms = $(state.max_terms): " *
                            "raise it with GaussianState(…; max_terms)"))
    end
    # computed only when not given, K² overlaps; zero when the components cancel to rounding beside the
    # scale of their coefficients, as the null state an operator annihilating the state leaves
    t = trace
    if isnothing(t)
        t = real(ComplexF64(norm2(s)))
        scale = sum(abs, s.coefs; init = 0.0)^2
        if abs(t) <= rounding_tol * scale
            t = 0.0
        end
    end
    if K == 1
        fix!(s)
        return single(state, only(s.comps).cov; parities, trace = t)
    end
    GaussianState{R}(state.system, state.axes, state.layout, nothing, s, parities, t, state.max_terms)
end

# the state as a sum that can be changed, its phases not brought up to date
function as_sum(state::GaussianState{<:Ket})
    if isnothing(state.sum)
        g = describe(state.cov)
        return GaussianSum([g], [sqrt(state.trace)], Random.Xoshiro(reference_seed))
    end
    deepcopy(state.sum)
end


# a mixed product state whose sites commute with their parity, held as |0⟩⟨0| ⊗ ρ: the ancilla empty, the
# covariance of each site from its density matrix, none between sites (their odd Majoranas average to zero);
# the ancilla of the pure states (|0⟩ψ_even + i|1⟩ψ_odd) would not leave a mixed state Gaussian
function GaussianState(state::State{Mixed}; axes = Dict{Int, Symbol}(), max_terms::Int = 64,
                       limits::Limits = Limits())
    system, axes, layout = product_layout(state, axes)
    t = nonzero_trace(state, "GaussianState")
    M = zeros(2 + layout.count, 2 + layout.count)
    M[1, 2] = 1
    for k in eachindex(system.sites)
        ρ = local_density(state, k)
        # the ancilla first, then the Majoranas of site k from layout.first[k] + 2
        js = layout.first[k] .+ (1:length(layout.majoranas[k])) .+ 1
        M[js, js] = site_covariance(layout, k, ρ / tr(ρ))
    end
    GaussianState{Mixed}(system, axes, layout, M - transpose(M), nothing, nothing, t, max_terms)
end

# the density matrix of site k of the mixed product `state`, not normalised, from its tensor: the other sites
# contribute their traces, and the signs of a fermionic partial trace only the elements odd in the parity,
# which `site_covariance` refuses
function local_density(state::State{Mixed}, k::Int)
    a = mps_array(state.state, k, siteind(Mixed, state.system, k))
    pairs_matrix(vec(a[:, invperm(mixed_elements(state.system, k)), :]), [dim(state.system[k])])
end

# the upper triangle of the covariance of the modes of site k in its density matrix ρ, normalised; refused
# when ρ does not commute with their parity, or is not the Gaussian state of that covariance, its values on
# the monomials of four Majoranas and more off Wick's theorem by more than 1e-12 (on a site of several modes:
# an Electron in (n↑ + n↓ - 2n↑n↓)/2 has the covariance of a Gaussian state, not its ⟨n↑n↓⟩)
function site_covariance(layout::MajoranaLayout, k::Int, ρ::AbstractMatrix)
    γs, P = layout.majoranas[k], layout.parity[k]
    if norm(ρ * P - P * ρ) > 1e-13 * norm(ρ)
        throw(ArgumentError("the Gaussian backend holds mixed states that commute with the parity: site $k has " *
                            "coherences between even and odd occupations"))
    end
    q = length(γs)
    m = zeros(q, q)
    for r in 1:q, t in r+1:q
        m[r, t] = real(-im * tr(ρ * γs[r] * γs[t]))
    end
    a = m - transpose(m)
    for (β, c) in local_monomials(γs, size(ρ, 1))
        if iseven(length(β)) && length(β) >= 4
            gap = abs(tr(ρ * c) - ComplexF64(expect_monomial(a, β)))
            if gap > 1e-12
                throw(ArgumentError("the Gaussian backend holds Gaussian mixed states, and the state of site $k " *
                                    "is not one: its correlations of $(length(β)) Majoranas miss Wick's theorem " *
                                    "by $(round(gap; sigdigits = 2))"))
            end
        end
    end
    m
end

# the axes given to GaussianState as a Dict from the qubits to their Paulis
function axes_by_site(system::System, axes)
    qubits = [k for k in eachindex(system.sites) if isnothing(fermion_species(system[k])) && dim(system[k]) == 2]
    if axes isa Symbol
        return Dict{Int, Symbol}(k => axes for k in qubits)
    elseif axes isa AbstractVector
        if length(axes) != length(system.sites)
            throw(ArgumentError("a vector of axes has one entry per site, $(length(system.sites)) here, not $(length(axes))"))
        end
        return Dict{Int, Symbol}(k => a for (k, a) in enumerate(axes) if !isnothing(a))
    end
    Dict{Int, Symbol}(axes)
end

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
    units = [(ComplexF64[1, im] / sqrt(2), [pauli_x, pauli_y], pauli_z)]
    for k in layout.order
        push!(units, (vectors[k], layout.majoranas[k], layout.parity[k]))
    end
    ev(v, m) = dot(v, m * v)
    parities = [ev(v, P) for (v, _, P) in units]
    # a site of definite parity when its amplitude in the other sector, ‖(1 - sP)v‖/2 with s the sign
    # of ⟨P⟩, is within 1e-13 of its norm, as in `in_sectors`: 1 - |⟨P⟩| sees it only squared
    definite = all(norm(v - sign(real(p)) * (P * v)) / 2 <= 1e-13 * norm(v)
                   for ((v, _, P), p) in zip(units[2:end], parities[2:end]))
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
trace(state::GaussianState) = state.trace
ncomponents(state::GaussianState) = isnothing(state.sum) ? 1 : length(state.sum.comps)
norm(state::GaussianState{<:Ket}) = sqrt(state.trace)

# the state as the representation R, its data unchanged
with_representation(::Type{R}, s::GaussianState) where {R <: Ket} =
    GaussianState{R}(s.system, s.axes, s.layout, s.cov, s.sum, s.parities, s.trace, s.max_terms)

# whether the backend keeps the norm of a trajectory exactly, see `exact_norm`
exact_norm(::GaussianState) = true

as_representation(::Type{R}, state::GaussianState{R}; limits::Limits = Limits()) where {R <: PM} = state
as_representation(::Type{Trajectory}, state::GaussianState{Pure}; limits::Limits = Limits()) =
    with_representation(Trajectory, state)
as_representation(::Type{Pure}, ::GaussianState{Trajectory}; limits::Limits = Limits()) =
    error("a trajectory is not a pure or a mixed state: run the simulation in Pure or Mixed")
as_representation(::Type{Pure}, ::GaussianState{Mixed}; limits::Limits = Limits()) =
    error("a mixed state cannot be made pure: ask for Mixed rather than Pure")
as_representation(::Type{Trajectory}, ::GaussianState{Mixed}; limits::Limits = Limits()) =
    not_yet("a trajectory of a mixed state")

# |ψ⟩⟨ψ| for ψ of one component and a definite parity, which commutes with the parity: its covariance, the
# ancilla made empty, its terms with the system zero by the definite parity
function as_representation(::Type{Mixed}, state::GaussianState{Pure}; limits::Limits = Limits())
    if !isnothing(state.sum)
        throw(ArgumentError("the Gaussian backend holds a mixed state of one Gaussian component, not the " *
                            "$(ncomponents(state)) of this state"))
    elseif isnothing(state.parities)
        throw(ArgumentError("the Gaussian backend holds mixed states that commute with the parity, and the parity " *
                            "of this pure state is not known: project it on a sector of the parity first"))
    end
    M = copy(state.cov)
    M[1:2, :] .= 0
    M[:, 1:2] .= 0
    M[1, 2], M[2, 1] = 1, -1
    GaussianState{Mixed}(state.system, state.axes, state.layout, M, nothing, nothing, state.trace, state.max_terms)
end
as_representation(::Type{Mixed}, ::GaussianState{Trajectory}; limits::Limits = Limits()) =
    error("a trajectory is not a pure or a mixed state: run the simulation in Pure or Mixed")

mix(state::GaussianState{Pure}; limits::Limits = Limits(), max_dense::Union{Nothing, Int} = nothing) =
    as_representation(Mixed, state)
mix(state::GaussianState{Mixed}; limits::Limits = Limits(), max_dense::Union{Nothing, Int} = nothing) = state
mix(::GaussianState{Trajectory}; kwargs...) = not_mixed()

# the value on the state, normalised, of a polynomial in its Majoranas: Σ c_α Pf(iM[α, α]) for one
# component, Σ c_α Σ_kl conj(a_k) a_l ⟨G_k|c(α)|G_l⟩ / ⟨ψ|ψ⟩ for several (see `sum_monomial`); the
# values of the monomials already met kept in cache
function poly_value(state::GaussianState, p::MajoranaPoly,
                    cache::Dict{Vector{Int}, ComplexF64} = Dict{Vector{Int}, ComplexF64}())
    nonzero_trace(state, "expect")
    value(α) = isnothing(state.sum) ? ComplexF64(expect_monomial(state.cov, α)) : sum_monomial(state, α)
    sum((c * get!(() -> value(α), cache, α) for (α, c) in p.terms); init = zero(ComplexF64))
end

# ⟨G_1|c(α)|G_2⟩ = ω Pf([W[α, α] Z[α, :]; -Z[α, :]ᵀ C]), Z and C empty for a transition matrix
bordered_pfaffian(W, Z, C, α) = logpfaffian([W[α, α] Z[α, :]; -transpose(Z[α, :]) C])

# (ω, W, Z, C) of `bordered_pfaffian` for the pair (G_1, G_2), nothing where the precise way is needed. The
# transition matrix -4i P₊Σ⁻¹ of Bravyi and Gosset (eq. (28)), Σ = M_1 + M_2 and P₊ = (I + iM_1)/2, splits on
# the singular vectors of Σ: W on the regular ones V_R, and Z(-4iG⁻¹)Zᵀ on those of the singular values at
# rounding level V_S, which components orthogonal by structure have (G and c(β)G: gates T, jumps, the branches
# of a ring), Z an orthonormal basis of P₊V_S and G = ZᵀΣZ; C = (i/4)G keeps these small singular values without
# dividing by them, after the robust Wick theorem of Chen and Scuseria (J. Chem. Phys. 158, 231102 (2023),
# eqs. (28), (32)-(34)).
function pair_transition(g1::PureGaussian, g2::PureGaussian)
    M1 = g1.cov
    Σ = M1 + g2.cov
    m = size(Σ, 1)
    F = svd(Σ)
    # the split is exact for zero singular values only, in groups of four (a pair of modes each); the others at
    # least 1e-3 of the largest, which bounds the error of W (2.3e-13 measured at a reciprocal condition 2.3e-3)
    k = count(<=(1e-13 * F.S[1]), F.S)
    if k % 4 != 0 || (k < m && F.S[m-k] < 1e-3 * F.S[1])
        return nothing
    end
    VS, VR = F.V[:, m-k+1:m], F.V[:, 1:m-k]
    W = -2im * (I + im * M1) * (VR * ((VR' * Σ * VR) \ VR'))
    if k == 0
        o = overlap(g1, g2)
        return iszero(o) ? nothing : (o, W, zeros(ComplexF64, m, 0), zeros(ComplexF64, 0, 0))
    end
    E = eigen(Hermitian(im * (VS' * M1 * VS)))
    if count(>(0), E.values) != k ÷ 2
        return nothing
    end
    Z = VS * E.vectors[:, E.values .> 0]
    C = (im / 4) * (transpose(Z) * Σ * Z)
    # ω by the precise way on the rows of Z of largest volume; |Pf| below 1e-8 left to the precise way as a
    # precaution (errors measured within 1.6e-14 down to |Pf| = 1.5e-6)
    α0 = sort(qr(transpose(Z), ColumnNorm()).p[1:k÷2])
    b = bordered_pfaffian(W, Z, C, α0)
    if b.logabs < log(1e-8)
        return nothing
    end
    cross_term(g1, α0, g2) / b, W, Z, C
end

# `pair_transition` for each pair k < l of a sum, brought up to date first. Computed once for a sum, which no
# operation changes once brought up to date.
const pair_cache = WeakKeyDict{GaussianSum, Matrix{Any}}()

function pair_transitions(s::GaussianSum)
    fix!(s)
    get!(pair_cache, s) do
        K = length(s.comps)
        P = Matrix{Any}(nothing, K, K)
        for k in 1:K, l in k+1:K
            P[k, l] = pair_transition(s.comps[k], s.comps[l])
        end
        P
    end
end

# Σ_kl conj(a_k) a_l ⟨G_k|c(α)|G_l⟩ / ⟨ψ|ψ⟩: the diagonal by Pf(iM_k[α, α]), the pairs l < k from k < l
# by ⟨G_k|c(α)|G_l⟩ = conj(⟨G_l|c(α)†|G_k⟩), c(α)† = ± c(α), α strictly increasing as in a MajoranaPoly.
# Where the parities of the components are known, a term where c(α)G_l lies in the other sector of the
# ancilla from G_k is zero exactly.
function sum_monomial(state::GaussianState{<:Ket}, α::Vector{Int})
    s = state.sum
    fix!(s)
    p = state.parities
    P = pair_transitions(s)
    w = length(α)
    adjoint_sign = iseven(w * (w - 1) ÷ 2) ? 1 : -1
    total = zero(LogComplex)
    for (k, (gk, ak)) in enumerate(zip(s.comps, s.coefs))
        if isnothing(p) || p[k] == flipped(p[k], α)
            total += conj(ak) * ak * expect_monomial(gk.cov, α)
        end
        for l in k+1:length(s.comps)
            if !isnothing(p) && p[k] != flipped(p[l], α)
                continue
            end
            gl, al = s.comps[l], s.coefs[l]
            if isnothing(P[k, l])
                v = cross_term(gk, α, gl)
            else
                ω, W, Z, C = P[k, l]
                v = ω * bordered_pfaffian(W, Z, C, α)
            end
            total += conj(ak) * al * v + conj(al) * ak * conj(v) * adjoint_sign
        end
    end
    ComplexF64(total) / state.trace
end

# measured on what conserves nothing, as on a State, an operator of no definite flux too
expect_value(state::GaussianState, a::IndexedOp{IsOperator}) = poly_value(state, cached_polynomial(state, a))

# the polynomial of the operator a on the layout of the state, kept by the identity of both: the rates of the jumps
# of a trajectory, measured at each step, are lowered once (72 % of a monitored trajectory went to lowering them
# again, profile P11); emptied beyond 1000 operators, which a program building new ones in a loop would reach;
# behind a lock, trajectories running on several threads
const polynomial_cache = IdDict{Any, MajoranaPoly}()
const polynomial_lock = ReentrantLock()

function cached_polynomial(state::GaussianState, a::IndexedOp{IsOperator})
    key = (a, state.layout, state.system, even_states(state))
    p = lock(() -> get(polynomial_cache, key, nothing), polynomial_lock)
    if isnothing(p)
        p = state_polynomial(state, lower(weaken(state.system, ()), a))
        lock(polynomial_lock) do
            if length(polynomial_cache) > 1000
                empty!(polynomial_cache)
            end
            polynomial_cache[key] = p
        end
    end
    p
end

# the matrix of size d with a 1 at (r, c)
function matrix_unit(d::Int, r::Int, c::Int)
    m = zeros(ComplexF64, d, d)
    m[r, c] = 1
    m
end

# the local indices of the basis states of sites of dimensions dims, in the order of kron (the first
# site the slowest)
function kron_indices(dims::Vector{Int})
    out = Vector{Int}[]
    for idx in 0:prod(dims)-1
        local_ = zeros(Int, length(dims))
        r = idx
        for i in length(dims):-1:1
            local_[i] = r % dims[i] + 1
            r ÷= dims[i]
        end
        push!(out, local_)
    end
    out
end

# ρ_A[a, b] = ⟨|b⟩⟨a|⟩, |b⟩⟨a| an operator on the sites given by its matrix and placed, its
# Jordan-Wigner strings laid by `lower`, so that tr(ρ_A matrix(op)) = ⟨op⟩ as for a State
function reduced_density_matrix(state::GaussianState, positions::AbstractVector{Int};
                                max_dense::Union{Nothing, Int} = dense_bound(state))
    check_positions(state, positions, "reduced_density_matrix")
    check_increasing(positions, "reduced_density_matrix")
    if isempty(positions)
        return ones(1, 1)
    end
    kept = collect(positions)
    check_reduced(state, kept, max_dense)
    system = weaken(state.system, ())
    dims = [dim(system[k]) for k in kept]
    basis = kron_indices(dims)
    ρ = zeros(ComplexF64, length(basis), length(basis))
    cache = Dict{Vector{Int}, ComplexF64}()
    for (ia, a) in enumerate(basis), (ib, b) in enumerate(basis)
        unit = reduce(kron, [matrix_unit(d, b[i], a[i]) for (i, d) in enumerate(dims)])
        op = Operator{length(kept)}("E($b,$a)", unit)(kept...)
        ρ[ia, ib] = poly_value(state, state_polynomial(state, lower(system, op)), cache)
    end
    ρ
end

"""
    von_neumann_entropy(state::GaussianState, positions; max_dense)

by the eigenvalues ±iν of the covariance of the modes of the sites, Σ h((1 + ν)/2), when the state has
a definite parity and the sites are fermionic, or qubits only and the first or the last ones (the
reduced state is then the Gaussian one of their modes); otherwise from `reduced_density_matrix`.
"""
function von_neumann_entropy(state::GaussianState, positions::AbstractVector{Int};
                             max_dense::Union{Nothing, Int} = dense_bound(state))
    check_positions(state, positions, "von_neumann_entropy")
    layout, order = state.layout, state.layout.order
    kept = sort(unique(positions); by = k -> layout.first[k])
    fermions = all(!isnothing(fermion_species(state.system[k])) for k in order)
    qubits = all(isnothing(fermion_species(state.system[k])) for k in order)
    ends = kept == order[1:length(kept)] || kept == order[end-length(kept)+1:end]
    if !(definite(state) && (fermions || (qubits && ends)))
        return invoke(von_neumann_entropy, Tuple{AbstractState, AbstractVector{Int}}, state, positions; max_dense)
    end
    majoranas = reduce(vcat, [layout.first[k] .+ (0:length(layout.majoranas[k])-1) .+ 2 for k in kept]; init = Int[])
    ν = eigvals(Hermitian(im * state.cov[majoranas, majoranas]))
    ps = (1 .+ ν[length(ν)÷2+1:end]) ./ 2
    shannon_entropy(vcat(ps, 1 .- ps))
end

function entanglement_entropy(state::GaussianState{<:Ket}, cut::Int)
    check_cut(length(state.system.sites), cut)
    von_neumann_entropy(state, 1:cut)
end

entanglement_entropy(::GaussianState{Mixed}, cut::Int) =
    throw(ArgumentError("the Gaussian backend does not compute the operator space entanglement entropy of a " *
                        "mixed state; von_neumann_entropy gives the entropy of a part"))

truncate(state::GaussianState; limits::Limits) = state

# the basis of the occupations of the modes of site k as the columns of a matrix, column 1 + Σ_r n_r 2^(r-1)
# being (c_1†)^n_1 ⋯ (c_q†)^n_q |0⟩ with c_r = (γ_{2r-1} + iγ_{2r})/2 from its local Majoranas, the
# vacuum the vector all c_r annihilate
function occupation_basis(layout::MajoranaLayout, k::Int)
    γs = layout.majoranas[k]
    q = length(γs) ÷ 2
    cs = [(γs[2r-1] + im * γs[2r]) / 2 for r in 1:q]
    vacuum = nullspace(reduce(vcat, cs))[:, 1]
    B = zeros(ComplexF64, length(vacuum), 2^q)
    for o in 0:2^q-1
        B[:, o + 1] = foldr((r, v) -> (o >> (r - 1)) & 1 == 1 ? cs[r]' * v : v, 1:q; init = vacuum)
    end
    B
end

# the array a with the matrix B applied to its axis k
function on_axis(a::Array, k::Int, B::AbstractMatrix)
    n = ndims(a)
    p = vcat(k, setdiff(1:n, k))
    b = reshape(B * reshape(permutedims(a, p), size(a, k), :), size(B, 1), size(a)[setdiff(1:n, k)]...)
    permutedims(b, invperm(p))
end

# the vector of the physical state in the basis of TMS, its phase free: Ψ = Σ_x ⟨x|Ψ⟩ |x⟩ over the
# even Fock states of the modes, the ancilla first (⟨x|Ψ⟩ = Σ_k a_k ⟨x|G_k⟩ for several components),
# and ψ = Ψ_even - iΨ_odd, an array of the occupations of the sites carried to their bases
function physical_vector(state::GaussianState{<:Ket})
    system, layout = state.system, state.layout
    order = layout.order
    dims = [dim(s) for s in system.sites]
    nmodes = layout.count ÷ 2 + 1
    if isnothing(state.sum)
        comps, coefs = [describe(state.cov)], [one(LogComplex)]
    else
        fix!(state.sum)
        comps, coefs = state.sum.comps, state.sum.coefs
    end
    a = zeros(ComplexF64, dims...)
    for idx in 0:2^nmodes-1
        x = [((idx >> (nmodes - j)) & 1) == 1 for j in 1:nmodes]
        if isodd(count(x))
            continue
        end
        amp = ComplexF64(sum((c * overlap(fock_state(x), g) for (g, c) in zip(comps, coefs)); init = zero(LogComplex)))
        occupations = zeros(Int, length(dims))
        pos = 2
        for k in order
            q = length(layout.majoranas[k]) ÷ 2
            occupations[k] = sum(x[pos + r - 1] << (r - 1) for r in 1:q) + 1
            pos += q
        end
        a[occupations...] = x[1] ? -im * amp : amp
    end
    for k in eachindex(system.sites)
        a = on_axis(a, k, occupation_basis(layout, k))
    end
    in_kron(a)
end

# the density matrix of a mixed state in the basis of TMS, in the order of kron: that of its modes
# (`fock_density`, the empty ancilla left out), its indices carried to the occupations of the sites as in
# `physical_vector`, then to the bases of the sites on both sides
function physical_matrix(state::GaussianState{Mixed})
    system, layout = state.system, state.layout
    dims = [dim(s) for s in system.sites]
    n, m = length(dims), layout.count ÷ 2
    # the ancilla held empty by every operation, up to rounding (4.5e-16 measured after a purification)
    M = state.cov
    if abs(M[1, 2] - 1) > 1e-10 || maximum(abs, M[1:2, 3:end]; init = 0.0) > 1e-10
        error("the ancilla of this mixed Gaussian state is not empty, which no operation should leave: " *
              "please report it")
    end
    ρ = fock_density(M[3:end, 3:end])
    place = Vector{CartesianIndex{n}}(undef, 2^m)
    for i in 1:2^m
        x = [((i - 1) >> (m - j)) & 1 == 1 for j in 1:m]
        occupations = zeros(Int, n)
        pos = 1
        for k in layout.order
            q = length(layout.majoranas[k]) ÷ 2
            occupations[k] = sum(x[pos + r - 1] << (r - 1) for r in 1:q) + 1
            pos += q
        end
        place[i] = CartesianIndex(occupations...)
    end
    a = zeros(ComplexF64, dims..., dims...)
    for i in 1:2^m, j in 1:2^m
        a[place[i], place[j]] = ρ[i, j]
    end
    for k in 1:n
        B = occupation_basis(layout, k)
        a = on_axis(on_axis(a, k, B), n + k, conj(B))
    end
    reshape(permutedims(a, [n:-1:1; 2n:-1:n+1]), prod(dims), prod(dims))
end

"""
    State(state::GaussianState; limits)

the state as an MPS, through its dense vector, from an overlap for each Fock state of its modes:
refused beyond 16 modes, where it would take minutes; a mixed state through its dense density
matrix, refused beyond `TensorMixedStates.dense_elements[]` elements, 4^m for m modes (2^24 at
first, 12 modes).
"""
function State(state::GaussianState{R}; limits::Limits = Limits()) where {R <: Ket}
    m = state.layout.count ÷ 2
    if m > 16
        throw(ArgumentError("the Gaussian backend converts to a State through 2^m overlaps for m modes, " *
                            "up to 16 modes, not $m"))
    end
    ψ = physical_vector(state)
    State(ExactState{R}(state.system, sqrt(state.trace) * ψ / norm(ψ), Val(:raw)); limits)
end

function State(state::GaussianState{Mixed}; limits::Limits = Limits())
    d = prod(float(dim(s)) for s in state.system.sites)
    check_elements(d^2, "the density matrix of the state", dense_bound(state); settable = false)
    State(ExactState{Mixed}(state.system, state.trace * physical_matrix(state), Val(:raw)); limits)
end

# one component: its covariance; several: for each, its covariance, its support, its amplitude and
# its coefficient, the numbers as their logarithm and phase
function write_state(g, state::GaussianState{R}) where R
    attributes(g)["representation"] = string(nameof(R))
    attributes(g)["axes"] = join(("$k:$P" for (k, P) in sort(collect(state.axes))), ",")
    if !isnothing(state.parities)
        g["parities"] = state.parities
    end
    attributes(g)["trace"] = state.trace
    attributes(g)["max_terms"] = state.max_terms
    if isnothing(state.sum)
        g["covariance"] = state.cov
        return nothing
    end
    s = state.sum
    fix!(s)
    for (k, (c, a)) in enumerate(zip(s.comps, s.coefs))
        h = create_group(g, "component$k")
        h["covariance"] = c.cov
        h["support"] = Int.(c.x)
        h["amplitude"] = [c.amp.logabs, real(c.amp.phase), imag(c.amp.phase)]
        h["coefficient"] = [a.logabs, real(a.phase), imag(a.phase)]
    end
    attributes(g)["components"] = length(s.comps)
    return nothing
end

function read_state(::Type{<:GaussianState}, g, sites)
    name = read(attributes(g)["representation"])
    representations = Dict("Pure" => Pure, "Trajectory" => Trajectory, "Mixed" => Mixed)
    if !haskey(representations, name)
        error("a GaussianState of unknown representation \"$name\"")
    end
    R = representations[name]
    system = System(sites)
    text = read(attributes(g)["axes"])
    axes = Dict{Int, Symbol}(parse(Int, first(split(e, ":"))) => Symbol(last(split(e, ":"))) for e in split(text, ",") if !isempty(e))
    layout = MajoranaLayout(system; axes)
    max_terms = haskey(attributes(g), "max_terms") ? Int(read(attributes(g)["max_terms"])) : 64
    trace = read(attributes(g)["trace"])
    parities = haskey(g, "parities") ? Vector{Int}(read(g["parities"])) : nothing
    if haskey(g, "covariance")
        cov = read(g["covariance"])
        # files of before the parities of the components: a flag for a state of one component
        if haskey(attributes(g), "definite") && read(attributes(g)["definite"]) == 1
            parities = [Int(sign(cov[1, 2]))]
        end
        return GaussianState{R}(system, axes, layout, cov, nothing, parities, trace, max_terms)
    end
    number(v) = LogComplex(v[1], complex(v[2], v[3]))
    hs = [g["component$k"] for k in 1:read(attributes(g)["components"])]
    comps = [PureGaussian(read(h["covariance"]), Bool.(read(h["support"])), number(read(h["amplitude"]))) for h in hs]
    coefs = [number(read(h["coefficient"])) for h in hs]
    s = GaussianSum(comps, coefs, Random.Xoshiro(reference_seed))
    GaussianState{R}(system, axes, layout, nothing, s, parities, trace, max_terms)
end
