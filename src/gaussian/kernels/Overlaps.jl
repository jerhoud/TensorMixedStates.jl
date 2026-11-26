# Overlaps with their phase, monomials c(α) = γ_{α_1} ⋯ γ_{α_w} on a state, and Wick's theorem.

function check_modes(g1::PureGaussian, g2::PureGaussian)
    if nmodes(g1) != nmodes(g2)
        throw(DimensionMismatch("states of $(nmodes(g1)) and $(nmodes(g2)) modes"))
    end
end

"""
    apply_monomial(α, g)

the state c(α)|ψ⟩ = γ_{α_1} ⋯ γ_{α_w} |ψ⟩, exactly.
"""
function apply_monomial(α::AbstractVector{<:Integer}, g::PureGaussian)
    D = ones(2nmodes(g))
    x = g.x
    amp = g.amp
    for b in Iterators.reverse(α)
        # γ_b M γ_b flips the row and the column b of the covariance, and γ_b|y⟩ = coef|x⟩
        D[b] = -D[b]
        y = copy(x)
        y[(b + 1) ÷ 2] = !y[(b + 1) ÷ 2]
        coef, _ = majorana_on_fock(b, y)
        amp *= conj(coef)
        x = y
    end
    PureGaussian(D .* g.cov .* transpose(D), x, amp)
end

"""
    overlap(g1, g2)

⟨ψ₁|ψ₂⟩ as a `LogComplex`, zero for opposite parities.
"""
function overlap(g1::PureGaussian, g2::PureGaussian)
    check_modes(g1, g2)
    n = nmodes(g1)
    if isodd(count(g1.x)) != isodd(count(g2.x))
        return zero(LogComplex)
    end
    # c(α)†|x₁⟩ = coef |x₂⟩, then ⟨ψ₁|x₁⟩⟨x₁|c(α)|ψ₂⟩ = uv
    α = [2j - 1 for j in 1:n if g1.x[j] != g2.x[j]]
    coef = one(ComplexF64)
    y = g1.x
    for b in α
        c, y = majorana_on_fock(b, y)
        coef *= c
    end
    uv = conj(g1.amp) * conj(coef) * g2.amp
    # Bravyi and Gosset, eq. (27), with (φ₀, φ₁, φ₂) = (ψ₁, |x₁⟩, ψ₂) as Dias and König do:
    # ⟨ψ₂|ψ₁⟩ uv = σ iⁿ 4⁻ⁿ Pf(R), and Pf(R) = Pf(iM₀) Pf(S) = iⁿ σ Pf(S)
    conj(LogComplex(-n * log(4), (-1)^n) * logpfaffian(overlap_matrix(g1.cov, g1.x, g2.cov, α)) / uv)
end

# S, the Schur complement of the block iM₀ in the matrix R_α of Bravyi and Gosset, eq. (27):
# R_α = [iM₀ -I I 0; I iM₁ -I 0; -I I iDM₂D Jᵀ+iDM₂Jᵀ; 0 0 -J+iJM₂D iJM₂Jᵀ], M₁ the covariance
# of |x₁⟩, J the rows α of the identity, D = I - JᵀJ. The block is eliminated exactly, with
# (iM₀)⁻¹ = iM₀ for a pure state, rather than the next one i(M₀ + M₁) too: a good support bounds
# its pfaffian, not its condition number.
function overlap_matrix(M0, x1, M2, α)
    m = size(M0, 1)
    w = length(α)
    M1 = fock_covariance(x1)
    keep = ones(m)
    keep[α] .= 0
    DM2D = keep .* M2 .* transpose(keep)
    S = zeros(ComplexF64, 2m + w, 2m + w)
    A, B, C = 1:m, m+1:2m, 2m+1:2m+w
    S[A, A] .= im .* (M1 .+ M0)
    S[A, B] .= -im .* M0
    S[B, A] .= -im .* M0
    S[B, B] .= im .* (DM2D .+ M0)
    for k in 1:m
        S[A[k], B[k]] -= 1
        S[B[k], A[k]] += 1
    end
    for (i, a) in enumerate(α)
        S[B, C[i]] .= im .* keep .* M2[:, a]
        S[B[a], C[i]] += 1
        S[C[i], B] .= im .* M2[a, :] .* keep
        S[C[i], B[a]] -= 1
        S[C[i], C] .= im .* M2[a, α]
    end
    S
end

"""
    cross_term(g1, α, g2)

⟨ψ₁|c(α)|ψ₂⟩ = ⟨ψ₁|γ_{α_1} ⋯ γ_{α_w}|ψ₂⟩ as a `LogComplex`.
"""
cross_term(g1::PureGaussian, α::AbstractVector{<:Integer}, g2::PureGaussian) =
    overlap(g1, apply_monomial(α, g2))

"""
    expect_monomial(M, α)

⟨ψ|c(α)|ψ⟩ = Pf(iM[α, α]) as a `LogComplex` for the state of covariance `M`, the indices of `α`
distinct.
"""
function expect_monomial(M::AbstractMatrix{<:Real}, α::AbstractVector{<:Integer})
    if !allunique(α)
        throw(ArgumentError("expect_monomial needs distinct Majoranas, reduce γ_b γ_b = 1 first"))
    end
    logpfaffian(im .* M[α, α])
end

"""
    transition_matrix(M1, M2)

the matrix W with ⟨ψ₁|γ_a γ_b|ψ₂⟩ = ⟨ψ₁|ψ₂⟩ W_ab (Bravyi and Gosset, eq. (28)), then
⟨ψ₁|c(α)|ψ₂⟩ = ⟨ψ₁|ψ₂⟩ Pf(W[α, α]); its error grows with the condition number of M1 + M2, which
is singular when ψ₁ and ψ₂ are orthogonal.
"""
function transition_matrix(M1::AbstractMatrix{<:Real}, M2::AbstractMatrix{<:Real})
    W = im .* conj.((-2I + im * M1 - im * M2) / (M1 + M2))
    W[diagind(W)] .= 1
    W
end
