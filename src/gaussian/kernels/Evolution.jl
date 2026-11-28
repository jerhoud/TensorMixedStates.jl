# Gaussian unitaries, and superpositions of Gaussian states that keep their relative phases under
# them. U = exp(-iHt) with H = (i/4) Σ h_ab γ_a γ_b fixes U γ_a U† only up to a phase on each
# state, so the phases are carried by links to a generic reference F evolved with the state:
# ⟨UF|UG_k⟩ = ⟨F|G_k⟩. This is the random reference of Guaita, Hackl and Quella (Ann. Phys. 495,
# 170757 (2026), arXiv:2409.11628, §III E), with descriptions of Dias and König; F_o = γ₁F serves
# the odd components.

"""
    propagator(h, t)

the rotation R = exp(h t), under which the covariance of exp(-iHt)|ψ⟩, H = (i/4) Σ h_ab γ_a γ_b
with `h` real antisymmetric, is R M Rᵀ.
"""
propagator(h::AbstractMatrix{<:Real}, t::Real) = exp(h * t)

antisymmetric(M) = (M - transpose(M)) / 2

# the normal modes of a real antisymmetric matrix of real Schur form Q T Qᵀ, as pairs (a, b) of columns
# of Q: a 2 × 2 block on (i, i + 1), then the real eigenvalues (zero to rounding) paired in their order;
# the pairing follows the blocks LAPACK returns, not the size of the eigenvalues
function normal_pairs(T::AbstractMatrix)
    m = size(T, 1)
    pairs = Tuple{Int, Int}[]
    real_eigenvalues = Int[]
    i = 1
    while i <= m
        if i < m && T[i+1, i] != 0
            push!(pairs, (i, i + 1))
            i += 2
        else
            push!(real_eigenvalues, i)
            i += 1
        end
    end
    for k in 1:2:length(real_eigenvalues)
        push!(pairs, (real_eigenvalues[k], real_eigenvalues[k+1]))
    end
    pairs
end

"""
    unitary_expectation(M, h, t, c = 0)

⟨ψ|e^{t(c + K)}|ψ⟩ with its phase, as a `LogComplex`, for the state ψ of covariance `M` and
K = (1/4) Σ h_ab γ_a γ_b, `h` real antisymmetric, so that e^{tK} = exp(-iHt) for H = (i/4) Σ h_ab γ_a γ_b.
In the normal modes of h, e^{tK} = Π_j (cos φ_j + sin φ_j γ'_a γ'_b), φ_j = ε_j t/2, whose expectation is
by Wick's theorem the pfaffian Pf(D iM' D + C), M' the covariance in the modes, D the sines on the first
Majorana of each mode, C the cosines on the pairs: no square root, no choice of sign.
"""
function unitary_expectation(M::AbstractMatrix{<:Real}, h::AbstractMatrix{<:Real}, t::Real, c::Number = 0)
    m = size(M, 1)
    if size(M) != (m, m) || size(h) != (m, m) || isodd(m)
        throw(DimensionMismatch("M and h are of the same size 2n, not $(size(M)) and $(size(h))"))
    end
    F = schur(Matrix{Float64}(antisymmetric(h)))
    T = F.T
    pairs = normal_pairs(T)
    O = F.Z[:, reduce(vcat, [[a, b] for (a, b) in pairs]; init = Int[])]
    d = ones(m)
    C = zeros(m, m)
    for (j, (a, b)) in enumerate(pairs)
        φ = (T[a, b] - T[b, a]) / 4 * t
        d[2j-1] = sin(φ)
        C[2j-1, 2j], C[2j, 2j-1] = cos(φ), -cos(φ)
    end
    A = Diagonal(d) * (im .* (transpose(O) * M * O)) * Diagonal(d) + C
    LogComplex(real(t * c), cis(imag(t * c))) * logpfaffian(A)
end

"""
    rotate(g, R)

the state U|ψ⟩ for the Gaussian unitary U of rotation `R`, up to a phase.
"""
rotate(g::PureGaussian, R::AbstractMatrix{<:Real}) = describe(antisymmetric(R * g.cov * transpose(R)))

"""
    GaussianSum(comps, coefs, rng)

the state Σ_k coefs[k] |comps[k]⟩, with a reference drawn from `rng`, which `evolve!` evolves and
`fix!` brings up to date.
"""
mutable struct GaussianSum
    comps::Vector{PureGaussian}
    coefs::Vector{LogComplex}
    # F, F_o (γ₁F when created), u_c = ⟨F|γ_c|F_o⟩ and the links ⟨F|G_k⟩ or ⟨F_o|G_k⟩, as at
    # the anchor
    ref::PureGaussian
    refodd::PureGaussian
    u::Vector{ComplexF64}
    links::Vector{LogComplex}
    # the rotation since the anchor, applied by fix!, and whether there is one
    R::Matrix{Float64}
    rotated::Bool
end

function GaussianSum(comps::AbstractVector{PureGaussian}, coefs::AbstractVector, rng)
    if length(comps) != length(coefs)
        throw(DimensionMismatch("$(length(comps)) components and $(length(coefs)) coefficients"))
    end
    if isempty(comps)
        throw(ArgumentError("a GaussianSum needs at least one component"))
    end
    n = nmodes(first(comps))
    for g in comps
        if nmodes(g) != n
            throw(DimensionMismatch("components of $n and $(nmodes(g)) modes"))
        end
    end
    O = Matrix(qr(randn(rng, 2n, 2n)).Q)
    if det(O) < 0
        O[:, 1] .*= -1
    end
    F = describe(antisymmetric(O * fock_covariance(zeros(Bool, n)) * transpose(O)))
    u = im .* F.cov[:, 1]
    u[1] = 1
    s = GaussianSum(collect(comps), LogComplex.(coefs), F, apply_monomial([1], F), u, LogComplex[],
                    Matrix{Float64}(I, 2n, 2n), false)
    relink!(s)
    s
end

nmodes(s::GaussianSum) = nmodes(s.ref)

reference(s::GaussianSum, g::PureGaussian) = iseven(count(g.x)) ? s.ref : s.refodd

# the links of the components to the references, after the components changed
function relink!(s::GaussianSum)
    if s.rotated
        throw(ArgumentError("relink! needs a GaussianSum brought up to date by fix!"))
    end
    s.links = [overlap(reference(s, g), g) for g in s.comps]
    if any(iszero, s.links)
        throw(ErrorException("a component is orthogonal to the reference of its GaussianSum"))
    end
    s
end

"""
    evolve!(s, R)

applies the Gaussian unitary of rotation `R`, of determinant 1 as those of `propagator`, to `s`;
`fix!` computes its effect.
"""
function evolve!(s::GaussianSum, R::AbstractMatrix{<:Real})
    s.R = R * s.R
    s.rotated = true
    s
end

# the phase of g changed so that ⟨ref|c(α)|g⟩ has the phase of target
function align(g::PureGaussian, ref::PureGaussian, α, target::Number)
    o = cross_term(ref, α, g)
    if iszero(o)
        throw(ErrorException("a component is orthogonal to the reference of its GaussianSum"))
    end
    PureGaussian(g.cov, g.x, g.amp * LogComplex(0.0, target) / LogComplex(0.0, o.phase))
end

"""
    fix!(s)

the components of `s` brought up to date with their phases, which become the new anchor.
"""
function fix!(s::GaussianSum)
    if !s.rotated
        return s
    end
    R = s.R
    if det(R) < 0
        throw(ArgumentError("a GaussianSum evolves by rotations of determinant 1, as those of propagator"))
    end
    F = rotate(s.ref, R)
    # ⟨UF|γ_b|UF_o⟩ = Σ_c R_bc ⟨F|γ_c|F_o⟩, taken at its largest component
    w = R * s.u
    b = argmax(abs.(w))
    Fo = align(rotate(s.refodd, R), F, [b], w[b])
    s.ref, s.refodd, s.u = F, Fo, w
    s.comps = [align(rotate(g, R), reference(s, g), Int[], l.phase) for (g, l) in zip(s.comps, s.links)]
    s.R = Matrix{Float64}(I, size(R))
    s.rotated = false
    s
end

"""
    evolve_branches!(s, branch, steps)

applies to each component k of `s` the operator e^{t(c + K)} of its branch, `(h, t, c) = steps[branch[k]]`,
K = (1/4) Σ h_ab γ_a γ_b with `h` real antisymmetric: a Gaussian unitary times e^{tc}, different from one
branch to another. The phase each branch takes is that of ⟨F|e^{t(c + K)}|F⟩ for the reference F
(`unitary_expectation`); the links to F are computed again, in O(K n³).
"""
function evolve_branches!(s::GaussianSum, branch::AbstractVector{<:Integer}, steps::AbstractDict)
    if length(branch) != length(s.comps)
        throw(DimensionMismatch("$(length(branch)) branches for $(length(s.comps)) components"))
    end
    fix!(s)
    F, Fo = s.ref, s.refodd
    moved = Dict{Int, Tuple{PureGaussian, PureGaussian, Matrix{Float64}, Float64}}()
    for p in unique(branch)
        h, t, c = steps[p]
        R = propagator(h, t)
        # U_p F with its phase, then U_p F_o from ⟨U_p F|γ_b|U_p F_o⟩ = Σ_c R_bc ⟨F|γ_c|F_o⟩, as in fix!
        Fp = align(rotate(F, R), F, Int[], unitary_expectation(F.cov, h, t, c).phase)
        w = R * s.u
        b = argmax(abs.(w))
        moved[p] = (Fp, align(rotate(Fo, R), Fp, [b], w[b]), R, real(t * c))
    end
    for k in eachindex(s.comps)
        Fp, Fop, R, loss = moved[branch[k]]
        g = s.comps[k]
        # ⟨U_p F|U_p G⟩ = ⟨F|G⟩, the link
        s.comps[k] = align(rotate(g, R), iseven(count(g.x)) ? Fp : Fop, Int[], s.links[k].phase)
        s.coefs[k] *= LogComplex(loss, 1)
    end
    relink!(s)
end
