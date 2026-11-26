# Measurements of the occupation of a mode, with the phase of the state after them.

fock_state(y) = PureGaussian(fock_covariance(y), Vector{Bool}(y), one(LogComplex))

# the probability of finding mode j empty (s = 1) or occupied (s = -1) in the pure state of
# covariance M; the less likely outcome uses 1 - M_cd² = Σ_{a ∉ {c,d}} (M_ca² + M_da²)/2 rather
# than 1 - |M_cd|, which would turn the rounding of M_cd into a probability and its root into a
# spurious amplitude
function occupation_probability(M, j, s)
    c, d = 2j - 1, 2j
    m = M[c, d]
    if s * m >= 0
        return (1 + abs(m)) / 2
    end
    r = 0.0
    for a in axes(M, 1)
        if a != c && a != d
            r += M[c, a]^2 + M[d, a]^2
        end
    end
    r / (4 * (1 + abs(m)))
end

"""
    project(g, j, occupied)

`(p, g′)`: the probability `p` of finding mode `j` occupied (or empty) and the normalised state after
it, with its phase; an error if `p = 0`. The error on g′ grows as `p` becomes small, not its
weight √p.
"""
function project(g::PureGaussian, j::Integer, occupied::Bool)
    s = occupied ? -1.0 : 1.0
    p = occupation_probability(g.cov, j, s)
    if p == 0
        throw(ArgumentError("mode $j cannot be found $(occupied ? "occupied" : "empty")"))
    end
    M = project_covariance!(copy(g.cov), j, s, p)
    sqrtp = LogComplex(log(p) / 2, 1)
    if g.x[j] == occupied
        # ⟨x|Πψ⟩ = ⟨x|ψ⟩, and |⟨x|ψ⟩|²/p ≥ 2⁻ⁿ: x stays a support
        return p, PureGaussian(M, copy(g.x), g.amp / sqrtp)
    end
    y, _ = find_support(M)
    p, PureGaussian(M, y, overlap(fock_state(y), g) / sqrtp)
end

"""
    project!(s, j, occupied)

replaces `s` by Π s, not normalised, Π the projector on mode `j` occupied (or empty): each component is
projected and its coefficient multiplied by the square root of its probability; the components whose
weight √p |a| is below 1e-13 times the largest are dropped, the projection of the others exact (a
component that cannot be found in that occupation has p zero to rounding, and its projection would be
noise); an error if all are. Returns the indices of the components kept.
"""
function project!(s::GaussianSum, j::Integer, occupied::Bool)
    fix!(s)
    sign = occupied ? -1.0 : 1.0
    ps = [occupation_probability(g.cov, j, sign) for g in s.comps]
    logw = [p > 0 ? log(p) / 2 + a.logabs : -Inf for (p, a) in zip(ps, s.coefs)]
    if all(==(-Inf), logw)
        throw(ArgumentError("mode $j cannot be found $(occupied ? "occupied" : "empty")"))
    end
    largest = maximum(logw)
    comps = PureGaussian[]
    coefs = LogComplex[]
    kept = Int[]
    for (k, (g, a, w)) in enumerate(zip(s.comps, s.coefs, logw))
        if w > -Inf && w >= largest + log(1e-13)
            p, h = project(g, j, occupied)
            push!(comps, h)
            push!(coefs, a * LogComplex(log(p) / 2, 1))
            push!(kept, k)
        end
    end
    s.comps, s.coefs = comps, coefs
    relink!(s)
    kept
end

"""
    norm2(s)

⟨ψ|ψ⟩ for the state `s`, as a `LogComplex`, by the K(K-1)/2 overlaps between its components;
brings `s` up to date.
"""
function norm2(s::GaussianSum)
    fix!(s)
    total = zero(LogComplex)
    for k in eachindex(s.comps)
        total += conj(s.coefs[k]) * s.coefs[k]
        for l in k+1:length(s.comps)
            t = conj(s.coefs[k]) * s.coefs[l] * overlap(s.comps[k], s.comps[l])
            total += t + conj(t)
        end
    end
    total
end
