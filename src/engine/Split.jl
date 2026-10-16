# The split of an operator of several sites into a sum of tensor products of operators of one
# site, each of a definite charge and a definite fermionic parity, which the lowering places on a
# system with their Jordan-Wigner strings. A superoperator, a channel or a generator, is split
# the same way on the density matrices of its sites, with a parity for the kets and one for the
# bras.

"""
    block_dim(R, site)

the dimension of the matrix of an operator on `site` for the representation `R`: that of the
site for `Pure`, its square for `Mixed`, a superoperator acting on its density matrices
"""
block_dim(::Type{Pure}, site::AbstractSite) = dim(site)
block_dim(::Type{Mixed}, site::AbstractSite) = dim(site)^2

"""
    Grade

the charge and the fermionic parities of an element of an operator on a site: for an operator,
its parity, 0 or 1; for a superoperator, that of the kets plus twice that of the bras, two bits
added separately
"""
struct Grade
    q::QN
    p::Int
end

(a::Grade + b::Grade) = Grade(a.q + b.q, a.p ⊻ b.p)
(a::Grade - b::Grade) = Grade(a.q - b.q, a.p ⊻ b.p)
(a::Grade == b::Grade) = a.q == b.q && a.p == b.p
hash(a::Grade, h::UInt) = hash(a.p, hash(a.q, h))

"""
    element_grades(R, site)

the grade of each element ``|a\\rangle\\langle b|`` of `site`, in the order of a vectorised
matrix: the charge and the parity of an operator of one site made of that element alone. For
`Mixed`, `a` and `b` are elements ``(k, b)`` of the density matrix, the ket faster, of the
charge of the mixed index, that of the ket minus that of the bra, starred where it is strong.
"""
function element_grades(::Type{Pure}, site::AbstractSite)
    c = basis_charges(site)
    p = fermion_parities(site)
    return vec([ Grade(c[a] - c[b], mod(p[a] + p[b], 2)) for a in eachindex(c), b in eachindex(c) ])
end

function element_grades(::Type{Mixed}, site::AbstractSite)
    c = basis_charges(site)
    p = fermion_parities(site)
    st = strong_names(site)
    q = vec([ c[k] - star(c[b], st) for k in eachindex(c), b in eachindex(c) ])
    π = vec([ mod(p[k], 2) + 2mod(p[b], 2) for k in eachindex(c), b in eachindex(c) ])
    return vec([ Grade(q[e] - q[f], π[e] ⊻ π[f]) for e in eachindex(q), f in eachindex(q) ])
end

"""
    string_odd(outs, ins)

whether an element of an operator on several sites, of parities `outs` and `ins` on each site,
changes sign when a tensor product of matrices becomes a tensor product of operators:
`kron(A, B)` is the matrix of `(A * F^p) ⊗ B`, `p` the parity of `B`, so an element
``|a\\rangle\\langle b|`` of a site changes sign when `b` is odd and so are the elements of
the sites after it
"""
function string_odd(outs, ins)
    odd, after = false, 0
    for j in length(outs):-1:1
        odd = odd ⊻ (isodd(ins[j]) && isodd(after))
        after += outs[j] + ins[j]
    end
    return odd
end

"""
    string_signs(R, sites)

the sign of each element of an operator on `sites`, one axis per site, see `string_odd`. A
superoperator on `Mixed` sites takes the sign of its kets times that of its bras: it is a sum of
`ρ ↦ A ρ B†`, whose elements are those of `A` times those of `conj(B)`.
"""
function string_signs(::Type{Pure}, sites::Vector{<:AbstractSite})
    d = [ dim(s) for s in sites ]
    p = [ fermion_parities(s) for s in sites ]
    js = eachindex(d)
    return map(CartesianIndices(Tuple(d .^ 2))) do c
        outs = [ p[j][mod1(c[j], d[j])] for j in js ]
        ins = [ p[j][cld(c[j], d[j])] for j in js ]
        return string_odd(outs, ins) ? -1 : 1
    end
end

function string_signs(::Type{Mixed}, sites::Vector{<:AbstractSite})
    d = [ dim(s) for s in sites ]
    p = [ fermion_parities(s) for s in sites ]
    js = eachindex(d)
    ket(e, j) = p[j][mod1(e, d[j])]
    bra(e, j) = p[j][cld(e, d[j])]
    return map(CartesianIndices(Tuple(d .^ 4))) do c
        # the elements (k, b) of the density matrix going out of and into each site
        eo = [ mod1(c[j], d[j]^2) for j in js ]
        ei = [ cld(c[j], d[j]^2) for j in js ]
        odd = string_odd([ ket(eo[j], j) for j in js ], [ ket(ei[j], j) for j in js ]) ⊻
              string_odd([ bra(eo[j], j) for j in js ], [ bra(ei[j], j) for j in js ])
        return odd ? -1 : 1
    end
end

"""
    along(f, y, j)

the linear map `f` applied along axis `j` of the array `y`
"""
function along(f::AbstractMatrix, y::AbstractArray, j::Int)
    p = [j; setdiff(1:ndims(y), j)]
    z = permutedims(y, p)
    w = reshape(f * reshape(z, size(z, 1), :), size(f, 1), size(z)[2:end]...)
    return permutedims(w, invperm(p))
end

"""
    part_on(x, d, on)

the part of an operator acting on the sites where `on` is true and as the identity on the
others, with an axis for each of the former: the component along the identity is projected
out on the former, the latter are traced out and divided by their dimension
"""
function part_on(x::AbstractArray, d, on)
    y = x
    # from the last site, so that dropping an axis leaves those still to come in place
    for j in length(d):-1:1
        e = vec(Matrix{Float64}(I, d[j], d[j]))
        if on[j]
            y = along(I - e * transpose(e) / d[j], y, j)
        else
            y = dropdims(along(transpose(e) / d[j], y, j); dims = j)
        end
    end
    return y
end

"""
    svd_terms(y, grades, g, tol)

the terms of an operator acting on each of its sites, as pairs of a coefficient and of one
vectorised factor per site: `y` has an axis per site, `grades` gives the grade of each value of
each axis and `g` is the grade of the operator. The sites are split off one by one by singular
value decompositions, made grade by grade so that every factor has a definite one.
"""
function svd_terms(y::AbstractArray, grades, g::Grade, tol)
    if ndims(y) == 1
        # what rounding left outside the grade of the operator is dropped, so that the
        # factor carries exactly that grade
        v = [ grades[1][k] == g ? y[k] : zero(eltype(y)) for k in eachindex(y) ]
        c = norm(v)
        return c ≤ tol ? [] : [ (c, [ v / c ]) ]
    end
    rest = size(y)[2:end]
    rg = vec([ sum(grades[l + 1][i[l]] for l in eachindex(rest)) for i in CartesianIndices(rest) ])
    ym = reshape(y, size(y, 1), :)
    terms = []
    for g1 in unique(grades[1])
        rows = findall(==(g1), grades[1])
        cols = findall(==(g - g1), rg)
        if isempty(cols)
            continue
        end
        f = svd(ym[rows, cols])
        for k in eachindex(f.S)
            if f.S[k] ≤ tol
                break
            end
            v = zeros(eltype(f.Vt), length(rg))
            v[cols] = f.Vt[k, :]
            sub = svd_terms(reshape(v, rest), grades[2:end], g - g1, tol / f.S[k])
            if isempty(sub)
                continue
            end
            u = zeros(eltype(f.U), size(y, 1))
            u[rows] = f.U[:, k]
            for (c, fs) in sub
                push!(terms, (f.S[k] * c, [ [u] ; fs ]))
            end
        end
    end
    return terms
end

"""
    split_matrix(R, name, m, sites)

the operator `name` of matrix `m` on `sites`, as `matrix` gives it, written as a sum of tensor
products of operators of one site, each of a definite charge and a definite fermionic parity:
pairs `(coef, factors)`, `factors[j]` the matrix of the term on site `j`, `nothing` for the
identity. Each part acting on some of the sites and as the identity on the others is split on
its own, so that no term spans more sites than it acts on. An operator of no definite charge is
refused; one of no definite parity has its parts split apart. For `Mixed`, `m` is the matrix of
a superoperator, each factor that of a superoperator of one site.
"""
function split_matrix(::Type{R}, name, m::AbstractMatrix, sites::Vector{<:AbstractSite}) where R
    n = length(sites)
    d = [ block_dim(R, s) for s in sites ]
    # real factors for a real operator, which give a real MPO, cheaper to contract
    if eltype(m) <: Complex && all(x -> iszero(imag(x)), m)
        m = real(m)
    end
    m = float(m)
    tol = rounding_tol * norm(m)
    # one axis per site, holding the vectorised matrix of an operator of one site: a matrix
    # reshaped puts the last site first and every output before every input, so the two axes
    # of each site are brought together, the output varying faster
    x = reshape(permutedims(reshape(m, (reverse(d)..., reverse(d)...)),
                            [ k for j in 1:n for k in (n - j + 1, 2n - j + 1) ]),
                Tuple(d .^ 2)) .* string_signs(R, sites)
    grades = [ element_grades(R, s) for s in sites ]
    totals = unique(sum(grades[k][c[k]] for k in 1:n)
                    for c in CartesianIndices(x) if abs(x[c]) > tol)
    # a FluxError, which expect answers by measuring on the state conserving nothing
    qs = unique(g.q for g in totals)
    if length(qs) > 1
        names = unique(q.name for s in sites for q in conserved(s).quantities)
        charge_refusal(name, unique(String[ x for s in sites for x in strong_names(s) ]), qs,
                       "$name carries no definite charge of $(join(names, ", ")), so it cannot " *
                       "act on $(join(sites, " ⊗ "))")
    end
    terms = Tuple{Number, Vector{Union{Nothing, Matrix}}}[]
    for on in Iterators.product(fill((false, true), n)...)
        y = part_on(x, d, on)
        js = findall(collect(on))
        if isempty(js)
            if abs(y[]) > tol
                push!(terms, (y[], fill(nothing, n)))
            end
            continue
        end
        for g in totals, (c, vs) in svd_terms(y, grades[js], g, tol)
            fs = Vector{Union{Nothing, Matrix}}(nothing, n)
            for (j, v) in zip(js, vs)
                fs[j] = reshape(v, d[j], d[j])
            end
            push!(terms, (c, fs))
        end
    end
    return terms
end
