# compact, which writes the terms of several sites of an operator as coms, the blocks of
# channels of Operators.jl, with no more channels than its terms need, and ≈, which compares
# two operators through the products of atoms they expand into, those of a com included.

export compact

################### Atoms ###################

"""
    add!(d, k, x)

add `x` to the value of `d` at `k`, zero when absent, and return `d`
"""
function add!(d::AbstractDict, k, x)
    d[k] = get(d, k, zero(x)) + x
    return d
end

"""
    side_terms(a)

the terms of `a`: those of its argument taken to the side of `a` for a `Left` or a `Right`,
which are linear, and `a` itself for any other operator
"""
side_terms(a::Left) = [ Left(x) for x in sumsubs(a.arg) ]
side_terms(a::Right) = [ Right(x) for x in sumsubs(a.arg) ]
side_terms(a::Op) = [a]

"""
    linearize(o)

the operator of one site `o` as a combination of its atoms, the operators it sums that are not
sums themselves, their coefficients taken apart: `2X + Y` gives `[X => 2, Y => 1]`, and
`Left(X + Y)` gives `[Left(X) => 1, Left(Y) => 1]`, gathered and sorted by `simplify_sum`.
"""
function linearize(o::Op)
    s = simplify_sum([ scalarcoef(t) * u for t in sumsubs(o) for u in side_terms(scalararg(t)) ])
    return Pair{Op, Number}[ scalararg(t) => scalarcoef(t) for t in sumsubs(s) if scalarcoef(t) ≠ 0 ]
end

"""
    coef_type(xs)

the floating point type the coefficients `xs` are computed in, `Float64` at least
"""
coef_type(xs) = float(mapreduce(typeof, promote_type, xs; init = Float64))

"""
    com_coefs(a)

the coefficients of the pieces of the com `a`, as `linearize` gives them
"""
com_coefs(a::ComOp) = (x for p in a.pieces for (_, _, o) in p for (_, x) in linearize(o))


################### Rank decisions ###################

"""
    interpolative(rows, scales, tol)

the positions of the rows kept among `rows`, in increasing order, and the coefficients `t`
expressing every row in terms of them, `rows[i] ≈ sum(t[i, j] * rows[kept[j]])`.

The rows are taken one at a time, the one with the largest residual on those already taken
first, as long as that residual exceeds `tol` times the scale of its row. Every other row is a
combination of the rows taken, computed from its projections, and a row itself under that
bound is dropped whole. A coefficient whose part in its row is under the bound is dropped too,
so that rounding leaves no piece of the order of 1e-17 in the MPO.

The pivoting is what keeps the coefficients bounded: taking the rows in their order, a row
taken while small makes the following ones its large multiples, and `λ^(j-i)` gets coefficients
growing as `λ^(-i)`.
"""
function interpolative(rows::Vector{Vector{T}}, scales::Vector{<:Real}, tol::Real) where T
    n = length(rows)
    res = copy.(rows)
    # the rows that are neither negligible nor taken yet
    cand = [ i for i in 1:n if norm(rows[i]) > tol * scales[i] ]
    basis = Vector{T}[]
    kept = Int[]
    while true
        best, bestnorm = 0, zero(real(T))
        for i in cand
            ν = norm(res[i])
            if ν > tol * scales[i] && ν > bestnorm
                best, bestnorm = i, ν
            end
        end
        if best == 0
            break
        end
        # orthogonalized once more, the residuals having been updated one vector at a time
        v = copy(res[best])
        for b in basis
            v .-= dot(b, v) .* b
        end
        q = v ./ norm(v)
        push!(basis, q)
        push!(kept, best)
        filter!(≠(best), cand)
        for i in cand
            res[i] .-= dot(q, res[i]) .* q
        end
    end
    k = length(kept)
    # the rows kept on the basis, triangular in the order they were taken: the coefficients of
    # a row are its projections solved against it, not the projections themselves
    r = [ dot(basis[a], rows[kept[b]]) for a in 1:k, b in 1:k ]
    knorms = [ norm(rows[i]) for i in kept ]
    t = zeros(T, n, k)
    for (a, i) in enumerate(kept)
        t[i, a] = one(T)
    end
    for i in cand
        c = [ dot(basis[a], rows[i]) for a in 1:k ]
        for a in k:-1:1
            for b in a+1:k
                c[a] -= r[a, b] * c[b]
            end
            c[a] /= r[a, a]
        end
        for a in 1:k
            if abs(c[a]) * knorms[a] > tol * scales[i]
                t[i, a] = c[a]
            end
        end
    end
    order = sortperm(kept)
    return kept[order], t[:, order]
end

"""
    dense_rows(entries)

the sorted keys the dictionaries `entries` hold, and the dictionaries as dense rows over them
"""
function dense_rows(entries::Vector{Dict{K, T}}) where {K, T}
    ks = sort!(unique!(K[ k for e in entries for k in keys(e) ]))
    pos = Dict(k => i for (i, k) in enumerate(ks))
    rows = Vector{T}[ zeros(T, length(ks)) for _ in entries ]
    for (g, e) in zip(rows, entries), (k, x) in e
        g[pos[k]] = x
    end
    return ks, rows
end


################### Blocks of channels ###################

"""
    SitePieces{T}

the pieces of a site of a block of channels: each pair of channels `(l, r)`, numbered as in
`ComOp`, mapped to the combination of atoms laid there, with coefficients of type `T`
"""
const SitePieces{T} = Dict{Tuple{Int, Int}, Dict{Op, T}}

"""
    struct ComBlock{T}

a com being built or reduced, its pieces kept as combinations of atoms with coefficients of
type `T`: `pieces[k]` holds those of the `k`-th site from `start`, and `dims[k]` is the number
of channels of the link on the right of that site.
"""
struct ComBlock{T}
    start::Int
    dims::Vector{Int}
    pieces::Vector{SitePieces{T}}
end

"""
    direct_pass(R, terms, T, tol)

the block of channels of the terms of several sites `terms`, each given by its coefficient and
its factors in the order of their sites, as `(site, combination of atoms)`. `R` is the
representation of the operator and `T` the type of the coefficients.

It is built from left to right in a single pass. On each site, the rows are the channels of the
link on the left continued by the identity or by an atom of the site, and the terms beginning
there; the columns are what remains of the terms on the right of the site, their suffixes. The
rows `interpolative` keeps are the channels of the link on the right, and the coefficients of
all the rows are the pieces of the site. A term thus opens on its first site and closes on its
last, and the channels are independent on their left: their number is the rank of the coupling
of the terms crossing the link whenever the suffixes are independent operators, which
`right_sweep` sees to otherwise.
"""
function direct_pass(::Type{R}, terms, ::Type{T}, tol::Real) where {R, T}
    id = IdentityOp{R, Generic, 1}()
    # the suffixes of the terms, each a node (site, factor, rest) of a tree the terms share, the
    # factors divided by their first coefficient so that proportional ones are one; and the
    # terms by their first site
    F = Vector{Pair{Op, T}}
    nodes = Tuple{Int, F, Int}[]
    node_ids = Dict{Tuple{Int, F, Int}, Int}()
    starts = Dict{Int, Vector{Tuple{T, F, Int}}}()
    for (c, fs) in terms
        c = T(c)
        normed = Tuple{Int, F}[]
        for (s, comb) in fs
            if any(x -> first(x) isa IdentityOp, comb)
                error("bug: the factor $comb of site $s holds the identity")
            end
            β = T(last(comb[1]))
            c *= β
            push!(normed, (s, [ a => T(x) / β for (a, x) in comb ]))
        end
        rest = 0
        for p in length(normed):-1:2
            key = (normed[p]..., rest)
            rest = get!(node_ids, key) do
                push!(nodes, key)
                length(nodes)
            end
        end
        s, f = normed[1]
        push!(get!(starts, s, Tuple{T, F, Int}[]), (c, f, rest))
    end
    first_site = minimum(keys(starts))
    last_site = maximum(first(last(fs)) for (_, fs) in terms)
    n = last_site - first_site + 1
    pieces = [ SitePieces{T}() for _ in 1:n ]
    dims = zeros(Int, n)
    # the channels of the link on the left, over the suffixes `cols`
    q = Vector{T}[]
    cols = Int[]
    for (j, k) in enumerate(first_site:last_site)
        m = length(q)
        rows = Dict{Tuple{Int, Op}, Dict{Int, T}}()
        for r in 1:m, (x, v) in zip(q[r], cols)
            if iszero(x)
                continue
            end
            s, f, rest = nodes[v]
            if s == k
                for (a, β) in f
                    if rest == 0
                        add!(get!(Dict{Op, T}, pieces[j], (r, 0)), a, x * β)
                    else
                        add!(get!(Dict{Int, T}, rows, (r, a)), rest, x * β)
                    end
                end
            else
                add!(get!(Dict{Int, T}, rows, (r, id)), v, x)
            end
        end
        for (c, f, rest) in get(starts, k, Tuple{T, F, Int}[]), (a, β) in f
            add!(get!(Dict{Int, T}, rows, (0, a)), rest, c * β)
        end
        # the channels continued first, in order, then the terms beginning here
        rowkeys = sort!(collect(keys(rows)); by = x -> (x[1] == 0 ? m + 1 : x[1], x[2]))
        newcols, grows = dense_rows([ rows[key] for key in rowkeys ])
        # a row continuing a channel is measured against all that the channel holds, so that
        # a remainder of rounding is not taken for a channel of its own. The pieces of the site
        # are only the closings yet, one per channel
        content = zeros(real(T), m)
        for ((r, _), comb) in pieces[j], (_, x) in sort!(collect(comb); by = first)
            content[r] += abs2(x)
        end
        for (key, g) in zip(rowkeys, grows)
            if key[1] > 0
                content[key[1]] += real(dot(g, g))
            end
        end
        scales = [ key[1] > 0 ? sqrt(content[key[1]]) : norm(g) for (key, g) in zip(rowkeys, grows) ]
        kept, t = interpolative(grows, scales, tol)
        for (i, (r, a)) in enumerate(rowkeys), p in eachindex(kept)
            if !iszero(t[i, p])
                add!(get!(Dict{Op, T}, pieces[j], (r, p)), a, t[i, p])
            end
        end
        q = grows[kept]
        cols = newcols
        dims[j] = length(kept)
    end
    if pop!(dims) ≠ 0
        error("bug: channels left open on the last site of the terms")
    end
    return ComBlock{T}(first_site, dims, pieces)
end

"""
    block_of(a, c, T)

the com `a` times `c` as a block of channels with coefficients of type `T`, the coefficient
laid on its openings, which every term takes once
"""
function block_of(a::ComOp, c::Number, ::Type{T}) where T
    pieces = map(a.pieces) do p
        d = SitePieces{T}()
        for (l, r, o) in p, (atom, β) in linearize(o)
            add!(get!(Dict{Op, T}, d, (l, r)), atom, T(l == 0 ? c * β : β))
        end
        d
    end
    return ComBlock{T}(a.start, copy(a.linkdims), pieces)
end

"""
    direct_sum(blocks)

the blocks of channels `blocks` side by side in one, the channels of each link numbered block
after block
"""
function direct_sum(blocks::Vector{ComBlock{T}}) where T
    first_site = minimum(b -> b.start, blocks)
    last_site = maximum(b -> b.start + length(b.pieces) - 1, blocks)
    n = last_site - first_site + 1
    dims = zeros(Int, n - 1)
    pieces = [ SitePieces{T}() for _ in 1:n ]
    for b in blocks
        o = b.start - first_site
        offs = [ dims[o + j] for j in eachindex(b.dims) ]
        for j in eachindex(b.dims)
            dims[o + j] += b.dims[j]
        end
        for (j, p) in enumerate(b.pieces), ((l, r), comb) in p
            pieces[o + j][(l == 0 ? 0 : offs[j - 1] + l, r == 0 ? 0 : offs[j] + r)] = comb
        end
    end
    return ComBlock{T}(first_site, dims, pieces)
end

"""
    left_sweep!(b, tol)

reduce the block `b`, from left to right, to channels whose left operators are independent,
and return it: on each link, a channel is a combination of the pieces entering it, and one
that is a combination of the others is merged into them, its pieces on the next site carried
over with the coefficients of the combination
"""
function left_sweep!(b::ComBlock{T}, tol::Real) where T
    for j in 1:length(b.pieces) - 1
        m = b.dims[j]
        entries = [ Dict{Tuple{Int, Op}, T}() for _ in 1:m ]
        for ((l, r), comb) in b.pieces[j], (a, x) in comb
            if r > 0
                entries[r][(l, a)] = x
            end
        end
        _, rows = dense_rows(entries)
        kept, t = interpolative(rows, norm.(rows), tol)
        if length(kept) == m
            continue
        end
        here = SitePieces{T}()
        for ((l, r), comb) in b.pieces[j]
            p = r == 0 ? 0 : findfirst(==(r), kept)
            if !isnothing(p)
                here[(l, p)] = comb
            end
        end
        next = SitePieces{T}()
        for ((l, r), comb) in sort!(collect(b.pieces[j + 1]); by = first)
            if l == 0
                next[(l, r)] = comb
                continue
            end
            for p in eachindex(kept)
                if !iszero(t[l, p])
                    d = get!(Dict{Op, T}, next, (p, r))
                    for (a, x) in comb
                        add!(d, a, t[l, p] * x)
                    end
                end
            end
        end
        b.pieces[j] = here
        b.pieces[j + 1] = next
        b.dims[j] = length(kept)
    end
    return b
end

"""
    mirror(b)

the block `b` read from right to left: its sites in the reverse order, the piece `(l, r)` of a
site becoming `(r, l)`. Mirroring twice gives `b` back, `start` being kept, which `left_sweep!`
does not read.
"""
mirror(b::ComBlock{T}) where T =
    ComBlock{T}(b.start, reverse(b.dims), [ SitePieces{T}((r, l) => c for ((l, r), c) in p) for p in reverse(b.pieces) ])

"""
    right_sweep(b, tol)

the block `b` reduced, from right to left, to channels whose right operators are independent:
`left_sweep!` on its mirror, where a channel is a combination of the pieces leaving it. After
`left_sweep!`, or after `direct_pass`, whose channels are already independent on their left,
the channels of each link are as few as the rank of the coupling across it allows.
"""
right_sweep(b::ComBlock, tol::Real) = mirror(left_sweep!(mirror(b), tol))

"""
    coms_of(R, b)

the block `b` as coms of the representation `R`, one for each run of links with channels, a
piece whose coefficients all vanish being left out
"""
function coms_of(::Type{R}, b::ComBlock) where R
    coms = ComOp{R}[]
    # the links without channels cut the block into its coms
    cuts = [ 0; findall(iszero, b.dims); length(b.pieces) ]
    for (c, d) in zip(cuts, cuts[2:end])
        if d - c < 2
            continue
        end
        pieces = map(c+1:d) do s
            p = Tuple{Int, Int, GenericOp{R, 1}}[]
            for (lr, comb) in sort!(collect(b.pieces[s]); by = first)
                o = simplify_sum(GenericOp{R, 1}[ x * a for (a, x) in comb ])
                if scalarcoef(o) ≠ 0
                    push!(p, (lr..., o))
                end
            end
            p
        end
        push!(coms, ComOp{R}(b.start + c, b.dims[c+1:d-1], pieces))
    end
    return coms
end


################### Reduction on the sites of a system ###################

"""
    atom_basis(atoms, site, tol)

the atoms of one site `atoms` written on a basis of them, compared through their matrices on
`site`: a dictionary from each atom to its coefficient on the identity, its combination of the
atoms kept, and the norm of its matrix.

The part of each matrix along the identity, its trace over that of the identity, is taken out
first, and `interpolative` keeps the atoms whose remainders are independent, each measured
against its whole matrix. The atoms kept are then independent together with the identity, and
an atom that is a number times the identity, as `S2` on a spin, or zero, as `Sp*Sp` on a spin
1/2, keeps none of them.
"""
function atom_basis(atoms::Vector{Op}, site::AbstractSite, tol::Real)
    ms = [ matrix(a, site) for a in atoms ]
    S = float(mapreduce(eltype, promote_type, ms; init = Float64))
    d = size(ms[1], 1)
    τ = [ S(sum(m[i, i] for i in 1:d) / d) for m in ms ]
    norms = [ norm(S.(vec(m))) for m in ms ]
    rows = [ S.(vec(m)) for m in ms ]
    for (g, x) in zip(rows, τ), i in 1:d
        g[(i - 1) * d + i] -= x
    end
    kept, t = interpolative(rows, norms, tol)
    basis = Dict{Op, Tuple{S, Vector{Pair{Op, S}}, real(S)}}()
    for (j, a) in enumerate(atoms)
        c = τ[j] - sum(t[j, p] * τ[kept[p]] for p in eachindex(kept); init = zero(S))
        basis[a] = (c, Pair{Op, S}[ atoms[kept[p]] => t[j, p] for p in eachindex(kept) if !iszero(t[j, p]) ], norms[j])
    end
    return basis
end

"""
    on_basis(comb, basis, id, tol)

the combination of atoms `comb` written on the atoms `basis` keeps and on the identity `id`. A
coefficient whose part is under `tol` times that of the whole combination is dropped, as
`interpolative` drops one, so that what cancels up to rounding leaves nothing.
"""
function on_basis(comb::Dict{Op, T}, basis, id::Op, tol::Real) where T
    out = Dict{Op, T}()
    for (a, x) in sort!(collect(comb); by = first)
        c, lin, _ = basis[a]
        if !iszero(c)
            add!(out, id, x * c)
        end
        for (b, y) in lin
            add!(out, b, x * y)
        end
    end
    scale = sqrt(sum(abs2(x) * basis[a][3]^2 for (a, x) in comb; init = 0.0))
    return filter!(p -> abs(last(p)) * basis[first(p)][3] > tol * scale, out)
end

"""
    push_identities!(b, id)

take the identity `id` out of the openings of the block `b`, from left to right, and return the
terms of one site this leaves, a combination for each site of `b`. A term opened by the
identity on a site begins in fact on the next one: `α` times the identity opening channel `r`
becomes `α` times each piece leaving `r` on the next site, an opening there, which may hold the
identity in turn, or a term of one site when the piece closes.
"""
function push_identities!(b::ComBlock{T}, id::Op) where T
    singles = [ Dict{Op, T}() for _ in b.pieces ]
    for j in 1:length(b.pieces) - 1
        opened = sort!([ (r, pop!(comb, id)) for ((l, r), comb) in b.pieces[j] if l == 0 && haskey(comb, id) ];
                       by = first)
        filter!(p -> !isempty(last(p)), b.pieces[j])
        for (r, α) in opened, ((l, s), comb) in sort!(collect(b.pieces[j + 1]); by = first)
            if l == r
                to = s == 0 ? singles[j + 1] : get!(Dict{Op, T}, b.pieces[j + 1], (0, s))
                for (a, x) in comb
                    add!(to, a, α * x)
                end
            end
        end
    end
    return singles
end

"""
    reduce_on_sites(R, a, c, system)

the com `a` times `c` as a block of channels on `system`, as few as its operators allow there,
and the terms of one site this leaves, a combination of atoms for each site of the block. It is
how `PreMPO` lays a com.

`compact` takes the atoms of a site as independent, `X*Y` and `Z` being two operators to it.
Here they are compared through their matrices on the site, see `atom_basis`, every piece is
written on a basis of them that the identity completes, and `push_identities!` takes the
identity out of the openings and of the closings. On a qubit, `N(1)*N(3) + Z(1)*Z(3)` has `N`
written `(1 - Z) / 2` and keeps `5/4 * Z(1)*Z(3)` as its part of several sites, the rest going
to the terms of one site and to the constant. The sweeps then reduce what the relations made
dependent. On each link the channels are as few as the rank of the operator across it, once
its parts that are the identity on either side are taken out, which no triangular MPO goes
below.
"""
function reduce_on_sites(::Type{R}, a::ComOp{R}, c::Number, sys::System) where R
    id = IdentityOp{R, Generic, 1}()
    bases = map(zip(com_sites(a), a.pieces)) do (k, p)
        atoms = unique!(sort!(Op[ id; [ atom for (_, _, o) in p for (atom, _) in linearize(o) ] ]))
        atom_basis(atoms, sys[k], rounding_tol)
    end
    T = coef_type([ c; collect(com_coefs(a)); [ first(v) for basis in bases for v in values(basis) ] ])
    b = block_of(a, c, T)
    for (j, basis) in enumerate(bases)
        b.pieces[j] = filter!(p -> !isempty(last(p)),
                              SitePieces{T}(lr => on_basis(comb, basis, id, rounding_tol) for (lr, comb) in b.pieces[j]))
    end
    singles = push_identities!(b, id)
    m = mirror(b)
    closed = push_identities!(m, id)
    b = right_sweep(left_sweep!(mirror(m), rounding_tol), rounding_tol)
    return b, [ mergewith(+, s, t) for (s, t) in zip(singles, reverse(closed)) ]
end


################### compact ###################

"""
    holds_identity(f)

whether the factor `f` of a product has the identity among its atoms, see `linearize`
"""
holds_identity(f) = f isa AtIndex && any(x -> first(x) isa IdentityOp, linearize(f.op))

"""
    expand_identities(s)

the simplified operator `s` with each term of several sites that has a factor holding the
identity written out, that factor replaced by the sum of its atoms on its site: `(Id + Z)(1) *
Z(2)`, which a merge of two factors of one site gives, becomes `Z(2) + Z(1) * Z(2)`. A term of
several sites is laid on the sites it acts on, which a part of it acting on fewer would not be.
"""
function expand_identities(s::IndexedOp)
    return simplify_sum(map(sumsubs(s)) do t
        fs = prodsubs(scalararg(t))
        if count(f -> !(f isa IdentityOp), fs) < 2 || !any(holds_identity, fs)
            return t
        end
        written = map(fs) do f
            holds_identity(f) ? simplify_sum([ x * a(f.index...) for (a, x) in linearize(f.op) ]) : f
        end
        return scalarcoef(t) * removeMulti(simplify_prod(written))
    end)
end

"""
    compact_simplified(s, tol, what)

`compact` of the operator `s`, simplified and with its Jordan-Wigner strings spelled out by
`removeMulti` already, `what` naming the caller in the refusal of a factor of several sites.
`PreMPO` and `make_obs` simplify the operator themselves, the one to keep its own messages,
the other to find the kind of the operator on that same simplification.
"""
function compact_simplified(s::IndexedOp{R}, tol::Real, what::String) where R
    s = expand_identities(s)
    kept = IndexedOp{R}[]
    coms = Tuple{Number, ComOp{R}}[]
    terms = Tuple{Number, Vector{Tuple{Int, Vector{Pair{Op, Number}}}}}[]
    for t in sumsubs(s)
        c, a = scalarcoef(t), scalararg(t)
        if a isa ComOp
            push!(coms, (c, a))
            continue
        end
        fs = filter(x -> !(x isa IdentityOp), prodsubs(a))
        foreach(o -> check_one_site(o, what), fs)
        if length(fs) < 2
            push!(kept, t)
        else
            push!(terms, (c, [ (only(f.index), linearize(f.op)) for f in fs ]))
        end
    end
    if isempty(terms) && isempty(coms)
        return s
    end
    T = coef_type([ [ c for (c, _) in terms ]; [ x for (_, fs) in terms for (_, comb) in fs for (_, x) in comb ];
                    [ c for (c, _) in coms ]; [ x for (_, a) in coms for x in com_coefs(a) ] ])
    blocks = ComBlock{T}[]
    if !isempty(terms)
        push!(blocks, direct_pass(R, terms, T, tol))
    end
    for (c, a) in coms
        push!(blocks, block_of(a, c, T))
    end
    b = right_sweep(left_sweep!(direct_sum(blocks), tol), tol)
    return simplify_sum([ kept; coms_of(R, b) ])
end

"""
    compact(op; tol = rounding_tol)

`op` written so that its MPO has the smallest bond dimension: its terms of one site and its
constant as they are, and its terms of several sites gathered into coms, blocks of channels
in which the terms share what they have in common, printed `com(sites,linkdims)`. `make_mpo`,
`PreMPO`, `tdvp`, `dmrg`, `approx_W`, `measure` and `expect` take the result as they take
`op`, and all but `expect` compact an operator themselves: calling `compact` saves doing it
again at each call, and speeds up `expect`, which measures an operator as it is given. An array
of operators is compacted element by element, which is how a time dependent evolver is
compacted, each term keeping its time function.

`op` is simplified first, and the atoms of its factors, the operators of one site they sum,
are taken as independent: `compact` does not know that `X*Y` is `im * Z` on a qubit, nor that
`N` is `(1 - Z) / 2`. `PreMPO` does: when it lays a com, it compares them through their
matrices on the sites of the system. The bond dimension of the MPO on each link is then 2 plus
the rank of the operator across it, once its parts that are the identity on either side are
taken out. In `N(1)*N(3) + Z(1)*Z(3)`, `N` is then written `(1 - Z) / 2`, which moves part of
the operator to terms of one site and to the constant: the operator is the same, its
approximations WI and WII are not.

`tol` decides what counts as zero: a channel whose part is below `tol` times all it holds is
dropped, and so is a term below `tol` times itself. The default, `rounding_tol`, drops only
what the package takes as rounding, and the operator stays exact; a larger value gives an
approximation of it, which is not the best one of its size.

A com can be added to other operators, multiplied by a number, measured and lifted to a
mixed representation, but neither multiplied by another operator nor made a gate: take the
product or the gate first, and compact it, as in `compact(Gate(A)(1, 2))`. `compact(op) ≈ op`
compares the two, term by term.

# Examples

    H = compact(sum(0.6^(j - i) * Z(i) * Z(j) for i in 1:20 for j in i+1:20))
    maxlinkdim(make_mpo(state, H))      # 3, as for the sum itself, which make_mpo compacts:
                                        # compact spares doing it again at each call
"""
compact(op::IndexedOp; tol::Real = rounding_tol) =
    compact_simplified(removeMulti(simplify(op)), tol, "compact")

compact(a::GenericOp; kwargs...) =
    error("compact takes an operator placed on sites, such as $(a((1:nsites(a))...)) rather than $a")

compact(a::AbstractArray; kwargs...) = map(x -> compact(x; kwargs...), a)


################### Comparison ###################

"""
    AtomProduct

a product of atoms, as `com_expansion` and `monomials` give them: its factors `(sites, atom)`
in the order of the sites
"""
const AtomProduct = Vector{Tuple{Tuple, Op}}

"""
    com_expansion(a)

the terms of the com `a` as products of atoms, a dictionary from each product to its
coefficient. Coefficients meant to cancel may leave products with a coefficient of the order of
rounding: the expansion is `a` to that precision, not to the bit.
"""
function com_expansion(a::ComOp)
    T = coef_type(com_coefs(a))
    result = Dict{AtomProduct, T}()
    open = Dict{AtomProduct, T}[]
    for (j, (k, ps)) in enumerate(zip(com_sites(a), a.pieces))
        next = [ Dict{AtomProduct, T}() for _ in 1:get(a.linkdims, j, 0) ]
        for (l, r, o) in ps
            from = l == 0 ? Dict(AtomProduct() => one(T)) : open[l]
            to = r == 0 ? result : next[r]
            lin = linearize(o)
            for (p, x) in from, (atom, β) in lin
                if !iszero(x)
                    add!(to, atom isa IdentityOp ? p : [ p; ((k,), atom) ], x * T(β))
                end
            end
        end
        open = next
    end
    return result
end

"""
    monomials(op)

the operator placed on sites `op` as the dictionary of its products of atoms to their
coefficients: its terms once simplified, with a product of combinations of atoms expanded into
products of atoms and a com into its terms, see `com_expansion`.
"""
function monomials(op::IndexedOp)
    d = Dict{AtomProduct, Number}()
    for t in sumsubs(removeMulti(simplify(op)))
        c, a = scalarcoef(t), scalararg(t)
        if a isa ComOp
            for (w, x) in com_expansion(a)
                add!(d, w, c * x)
            end
            continue
        end
        fs = filter(x -> !(x isa IdentityOp), prodsubs(a))
        combs = [ [ (f.index, atom) => β for (atom, β) in linearize(f.op) ] for f in fs ]
        for choice in distribute(combs)
            # the identity left out, as com_expansion does
            add!(d, AtomProduct([ first(x) for x in choice if !(last(first(x)) isa IdentityOp) ]),
                 c * prod(last, choice; init = 1))
        end
    end
    return d
end

"""
    a ≈ b

for two operators placed on sites, whether their terms are equal up to the tolerances `atol`
and `rtol`, as `isapprox` compares two vectors, `rtol` defaulting to zero when `atol` is
given: both are simplified and written as sums of products of atoms, see `monomials`, those of
a com included. The comparison is symbolic, as `==` is: `X(1) ≈ (Sp + Sm)(1)` is false, on a
qubit as well.
"""
function isapprox(a::IndexedOp{R}, b::IndexedOp{R};
                  atol::Real = 0, rtol::Real = atol > 0 ? 0 : sqrt(eps(Float64))) where R
    da, db = monomials(a), monomials(b)
    nrm(d) = sqrt(sum(abs2, values(d); init = 0.0))
    # a product of `db` alone comes with its sign, which the norm does not see
    return nrm(mergewith(-, da, db)) ≤ max(atol, rtol * max(nrm(da), nrm(db)))
end
