# The WI and WII approximations of the exponential of an evolver, as MPOs built from its PreMPO,
# their combinations of orders 1 to 4, and the values of time functions they are taken at.

"""
    w_charges(pre, coefs)

the charges of the links of the approximations WI and WII, where one channel serves both the
term not yet begun and the term finished, so that the operator must have zero flux
"""
function w_charges(pre::PreMPO, coefs::Vector{Float64})
    q = mpo_charges(pre, coefs)
    total = q[1][end]
    if total ≠ QN()
        error("the approximations WI and WII need an operator of zero flux, and this one " *
              "carries $(-total), so it would move the charge its system conserves")
    end
    return q
end

"""
    make_approx_W1(pre, tau[, coefs])

the MPO of the approximation WI of the exponential of `tau` times the operator of `pre`,
`coefs` being the values of its time functions, see `make_mpo`
"""
function make_approx_W1(pre::PreMPO{R}, tau::Number, coefs::AbstractVector = [1.]) where R
    coefs = check_coefs(pre, coefs)
    sys = pre.system
    elt = promote_type(mpo_eltype(pre), typeof(tau))
    links = mpo_links(pre, coefs, w_charges, 0)
    ts = map(1:length(sys)) do i
        idx = siteind(R, sys, i)
        llink, rlink = links[i], links[i+1]
        a = site_array(elt, idx, llink, rlink)
        add_block!(a, 1, 1, Matrix{elt}(I, dim(idx), dim(idx)))
        for (l, r, _, ref, m) in pre.terms[i]
            c = coefs[ref]
            if c ≠ 0
                # the coefficient and the time step on the closing piece alone, as in make_mpo
                add_block!(a, l, r, m, r == 1 ? c * tau : one(c))
            end
        end
        return site_tensor(a, idx, llink, rlink)
    end
    return close_ends(ts, links, 1)
end

"""
    duhamel(m, x)
    duhamel(m, x, y)

``\\int_0^1 e^{sm} x e^{(1-s)m} ds``, and ``\\int_{0<s<t<1} e^{(1-t)m} x e^{(t-s)m} y e^{sm}
ds\\,dt``, read in the exponential of a block triangular matrix, as in Van Loan, Computing
integrals involving the matrix exponential, IEEE Trans. Autom. Control 23, 395 (1978)
"""
function duhamel(m::Matrix, x::Matrix)
    k = size(m, 1)
    z = zero(m)
    return exp([m x; z m])[1:k, k+1:2k]
end

function duhamel(m::Matrix, x::Matrix, y::Matrix)
    k = size(m, 1)
    z = zero(m)
    return exp([m x z; z m y; z z m])[1:k, 2k+1:3k]
end

"""
    make_approx_W2(pre, tau[, coefs])

the MPO of the approximation WII of Zaletel et al., Phys. Rev. B 91, 165112 (2015), of the
exponential of `tau` times the operator of `pre`, `coefs` being the values of its time
functions, see `make_mpo`. Its error, of order ``\\tau^2``, comes from the products of terms
that cross a same link.
"""
function make_approx_W2(pre::PreMPO{R}, tau::Number, coefs::AbstractVector = [1.]) where R
    coefs = check_coefs(pre, coefs)
    sys = pre.system
    elt = promote_type(mpo_eltype(pre), typeof(tau))
    dims = channel_counts(pre)
    links = mpo_links(pre, coefs, w_charges, 0)
    ts = map(1:length(sys)) do i
        idx = siteind(R, sys, i)
        llink, rlink = links[i], links[i+1]
        ldim, rdim = dims[i], dims[i+1]
        k = dim(idx)
        v = Matrix{Union{Nothing, Matrix{elt}}}(nothing, ldim, rdim)
        for (l, r, _, ref, u) in pre.terms[i]
            c = coefs[ref]
            if c ≠ 0
                # the coefficient and the time step on the closing piece alone, as in make_mpo
                m = (r == 1 ? c * tau : one(c)) * u
                v[l, r] = isnothing(v[l, r]) ? m : v[l, r] + m
            end
        end
        # their equation 11: each block is read in the exponential of the terms of one site D,
        # the transport A from channel l to channel r, the closing B of l and the opening C of
        # r, each taken once at most. On dense matrices, whose zeros outside the flux of a
        # block stay exact zeros
        d = something(v[1, 1], zeros(elt, k, k))
        closing = v[:, 1]
        opening = v[1, :]
        a = site_array(elt, idx, llink, rlink)
        add_block!(a, 1, 1, exp(d))
        for l in 2:ldim
            if !isnothing(closing[l])
                add_block!(a, l, 1, duhamel(d, closing[l]))
            end
        end
        for r in 2:rdim
            if !isnothing(opening[r])
                add_block!(a, 1, r, duhamel(d, opening[r]))
            end
        end
        for l in 2:ldim, r in 2:rdim
            block = isnothing(v[l, r]) ? nothing : duhamel(d, v[l, r])
            if !isnothing(closing[l]) && !isnothing(opening[r])
                both = duhamel(d, closing[l], opening[r]) + duhamel(d, opening[r], closing[l])
                block = isnothing(block) ? both : block + both
            end
            if !isnothing(block)
                add_block!(a, l, r, block)
            end
        end
        return site_tensor(a, idx, llink, rlink)
    end
    return close_ends(ts, links, 1)
end

"""
    approx_W_mpos(pre, tau, algo[, coefs])

the MPOs of the approximation of `algo`, an `ApproxW`, of the exponential of `tau` times the
operator of `pre`, to apply in turn
"""
function approx_W_mpos(pre::PreMPO, tau::Number, algo::ApproxW, coefs::AbstractVector = [1.])
    make = algo.w == 1 ? make_approx_W1 : make_approx_W2
    return [ make(pre, tau * c, coefs) for c in w_approx_coefs[algo.order] ]
end

"""
    step_coefs(coefs, t, dt, order)

the values of the time functions `coefs` for each exponential of a step of `dt` from `t`, in
the order they are applied: at the middle of the step up to order 2, and from order 3 the
commutator-free Magnus integrator of order 4 of Blanes and Moan, two exponentials combining
the values at the two Gauss points of the step
"""
function step_coefs(coefs, t::Number, dt::Number, order::Int)
    if order ≤ 2
        return [ map(f -> f(t + dt / 2), coefs) ]
    end
    v1, v2 = (map(f -> f(t + (1/2 + s * √3/6) * dt), coefs) for s in (-1, 1))
    w1, w2 = (3 + 2√3) / 12, (3 - 2√3) / 12
    return [ w1 * v1 + w2 * v2, w2 * v1 + w1 * v2 ]
end
