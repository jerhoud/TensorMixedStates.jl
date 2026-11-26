# Pure Gaussian states with their phase, after Dias and König (Quantum 8, 1350 (2024),
# arXiv:2307.12912): the covariance fixes the state up to a phase, which the amplitude on a Fock
# state carries.

"""
    PureGaussian(cov, x, amp)

the pure Gaussian state |ψ⟩ of covariance `cov` with `amp = ⟨x|ψ⟩`, where |x⟩ is chosen with
|⟨x|ψ⟩|² ≥ 2⁻ⁿ so that the phase is known to rounding.
"""
struct PureGaussian
    cov::Matrix{Float64}
    x::Vector{Bool}
    amp::LogComplex
end

nmodes(g::PureGaussian) = length(g.x)

"""
    find_support(M)

`(x, logp)` with |x⟩ a Fock state of probability `exp(logp) ≥ 2⁻ⁿ` in the pure Gaussian state of
covariance `M`, found by taking in turn the more likely occupation of each mode given the ones
before (Dias and König, findsupport).
"""
function find_support(M0::AbstractMatrix{<:Real})
    M = Matrix{Float64}(M0)
    n = size(M, 1) ÷ 2
    x = zeros(Bool, n)
    logp = 0.0
    for j in 1:n
        c, d = 2j - 1, 2j
        q = (1 + M[c, d]) / 2
        # either outcome is valid when q is near 1/2, so where the comparison falls is harmless
        s = q >= 1/2 ? 1.0 : -1.0
        p = s > 0 ? q : 1 - q
        x[j] = s < 0
        logp += log(p)
        # the covariance after finding mode j in that occupation, by Wick's theorem
        wc = M[:, c]
        wd = M[:, d]
        f = s / (2p)
        for b in axes(M, 2), a in axes(M, 1)
            M[a, b] += f * (wd[a] * wc[b] - wc[a] * wd[b])
        end
        M[c, :] .= 0
        M[:, c] .= 0
        M[d, :] .= 0
        M[:, d] .= 0
        M[c, d] = s
        M[d, c] = -s
    end
    x, logp
end

"""
    describe(M)

the pure Gaussian state of covariance `M`, its phase chosen so that its amplitude on its
support is positive.
"""
function describe(M::AbstractMatrix{<:Real})
    x, logp = find_support(M)
    PureGaussian(Matrix{Float64}(M), x, LogComplex(logp / 2, 1))
end
