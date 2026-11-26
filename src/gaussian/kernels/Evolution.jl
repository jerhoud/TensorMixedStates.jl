# Gaussian unitaries, and superpositions of Gaussian states that keep their relative phases under
# them. U = exp(-iHt) with H = (i/4) Σ h_ab γ_a γ_b fixes U γ_a U† only up to a phase on each
# state, so the phases are carried by links to a generic reference F evolved with the state:
# ⟨UF|UG_k⟩ = ⟨F|G_k⟩. This is the random reference of Guaita, Hackl and Quella (arXiv:2409.11628,
# §III E), with descriptions of Dias and König; F_o = γ₁F serves the odd components.

"""
    propagator(h, t)

the rotation R = exp(h t), under which the covariance of exp(-iHt)|ψ⟩, H = (i/4) Σ h_ab γ_a γ_b
with `h` real antisymmetric, is R M Rᵀ.
"""
propagator(h::AbstractMatrix{<:Real}, t::Real) = exp(h * t)

antisymmetric(M) = (M - transpose(M)) / 2

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
    # the rotation since the anchor, applied by fix!
    R::Matrix{Float64}
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
                    Matrix{Float64}(I, 2n, 2n))
    relink!(s)
    s
end

nmodes(s::GaussianSum) = nmodes(s.ref)

reference(s::GaussianSum, g::PureGaussian) = iseven(count(g.x)) ? s.ref : s.refodd

# the links of the components to the references, after the components changed
function relink!(s::GaussianSum)
    if s.R != I
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
    if s.R == I
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
    s
end
