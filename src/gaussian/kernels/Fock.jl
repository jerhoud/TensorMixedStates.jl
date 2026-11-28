# Fock states |x⟩, x[j] the occupation of mode j.

"""
    fock_covariance(x)

the covariance of the Fock state |x⟩.
"""
function fock_covariance(x::AbstractVector{Bool})
    n = length(x)
    M = zeros(2n, 2n)
    for j in 1:n
        s = x[j] ? -1.0 : 1.0
        M[2j-1, 2j] = s
        M[2j, 2j-1] = -s
    end
    M
end

"""
    majorana_on_fock(b, x)

`(coef, y)` such that γ_b |x⟩ = coef |y⟩.
"""
function majorana_on_fock(b::Integer, x::AbstractVector{Bool})
    j = (b + 1) ÷ 2
    s = iseven(count(view(x, 1:j-1))) ? 1.0 : -1.0
    coef = isodd(b) ? complex(s) : s * im * (x[j] ? -1 : 1)
    y = Vector{Bool}(x)
    y[j] = !y[j]
    coef, y
end

"""
    fock_density(M)

the density matrix of the mixed Gaussian state of covariance `M` (real antisymmetric, of eigenvalues
in [-1, 1]) on the Fock states |x⟩ of its modes, the index of |x⟩ having x as its bits, the first mode the
most significant.
"""
function fock_density(M::AbstractMatrix{<:Real})
    m = size(M, 1) ÷ 2
    N = 2^m
    index(x) = foldl((i, b) -> 2i + b, x; init = 0) + 1
    # γ_c as the signed permutation |x⟩ → coefs[c][i] |targets[c][i]⟩, i the index of x
    targets = [zeros(Int, N) for _ in 1:2m]
    coefs = [zeros(ComplexF64, N) for _ in 1:2m]
    for i in 1:N
        x = [((i - 1) >> (m - j)) & 1 == 1 for j in 1:m]
        for c in 1:2m
            coef, y = majorana_on_fock(c, x)
            targets[c][i] = index(y)
            coefs[c][i] = coef
        end
    end
    # (Σ_c z_c γ_c ρ)ᵀ from ρᵀ, column by column
    function combined(z, W)
        V = zeros(ComplexF64, N, N)
        for c in eachindex(z)
            if z[c] != 0
                for i in 1:N
                    view(V, :, targets[c][i]) .+= (z[c] * coefs[c][i]) .* view(W, :, i)
                end
            end
        end
        V
    end
    # ρ = 2^{-m} Π_k (1 - i μ_k η_a η_b) over the blocks [0 μ_k; -μ_k 0] of the real Schur form M = Z T Zᵀ,
    # η = Zᵀγ, the factors commuting; a block of one is a zero eigenvalue, a factor 1
    F = schur(Matrix{Float64}(M))
    T, Z = F.T, F.Z
    W = Matrix{ComplexF64}(I, N, N) / N
    a = 1
    while a <= 2m
        if a < 2m && T[a + 1, a] != 0
            μ = (T[a, a + 1] - T[a + 1, a]) / 2
            W = W - im * μ * combined(view(Z, :, a), combined(view(Z, :, a + 1), W))
            a += 2
        else
            a += 1
        end
    end
    transpose(W)
end
