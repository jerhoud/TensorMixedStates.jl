# The compaction of the terms of an operator on a system into a block of channels, as few on each
# link as the rank of the operator across it once its parts acting as the identity on either
# side are taken out. The atoms of a site are the matrices the terms have there, their strings
# laid on them, and those of the strings of the terms crossing it without a factor, compared
# through their matrices.

"""
    add!(d, k, x)

add `x` to the value of `d` at `k`, zero when absent, and return `d`
"""
function add!(d::AbstractDict, k, x)
    d[k] = get(d, k, zero(x)) + x
    return d
end

"""
    interpolative(rows, scales, tol)

the positions of the rows kept among `rows`, in increasing order, and the coefficients `t`
expressing every row in terms of them, `rows[i] ≈ sum(t[i, j] * rows[kept[j]])`.

The row of largest residual is taken first, as long as that residual exceeds `tol` times the
scale of its row; a row, or a coefficient, under that bound is dropped. This pivoting keeps
the coefficients bounded.
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
        # orthogonalized again: the residuals were updated one vector at a time
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
    # triangular in the order the rows were taken: the coefficients of a row are its
    # projections solved against it, not the projections themselves
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

"""
    Combination{T}

an operator of one site as a combination of the matrices of its site in a block of channels,
see `ComBlock`, numbered from 1, the identity being 0, with coefficients of type `T`
"""
const Combination{T} = Dict{Int, T}

"""
    SitePieces{T}

the pieces of a site of a block of channels: each pair of channels `(l, r)` mapped to the
combination laid there, `l == 0` for a term opening on the site, `r == 0` for one closing
"""
const SitePieces{T} = Dict{Tuple{Int, Int}, Combination{T}}

"""
    struct ComBlock{T}

a block of channels: `pieces[k]` holds those of the `k`-th site from `start`, `atoms[k]` the
matrices their combinations number on that site, its atoms, and their basis once `on_basis!` has
written the pieces on it, and `dims[k]` is the number of channels of the link on the right of
that site, for every site but the last.
"""
struct ComBlock{T}
    start::Int
    dims::Vector{Int}
    pieces::Vector{SitePieces{T}}
    atoms::Vector{Vector{Matrix{T}}}
end

"""
    direct_pass(terms, strings, tol)

the first site, the numbers of channels and the pieces of the block of channels of the terms of
several sites `terms`, see `ComBlock`, each given by its coefficient and its factors in the
order of their sites, as `(site, atom, g)`, `g` being the strings of the term on the sites
between the factor before and this one, `strings(k, g)` their atom on site `k`, 0 for the
identity.

Built from left to right: on each site, the channels on the left continued by an atom, or by
the strings of their terms where these have no factor, and the terms beginning there, are rows
over the suffixes of the terms; the rows `interpolative` keeps are the channels on the right,
and the coefficients of all the rows the pieces of the site. A term opens on its first site and
closes on its last. The atoms are taken as independent: the relations `on_basis!` finds between
them, and the identity `push_identities!` moves out, leave channels dependent on either side,
which the sweeps of `compacted` reduce.
"""
function direct_pass(terms::Vector{Tuple{T, Vector{Tuple{Int, Int, Int}}}}, strings,
                     tol::Real) where T
    # the suffixes of the terms, nodes (site, atom, strings before, rest) of a tree the terms share
    nodes = Tuple{Int, Int, Int, Int}[]
    node_ids = Dict{Tuple{Int, Int, Int, Int}, Int}()
    starts = Dict{Int, Vector{Tuple{T, Int, Int}}}()
    for (c, fs) in terms
        rest = 0
        for p in length(fs):-1:2
            key = (fs[p]..., rest)
            rest = get!(node_ids, key) do
                push!(nodes, key)
                length(nodes)
            end
        end
        s, a, _ = fs[1]
        push!(get!(starts, s, Tuple{T, Int, Int}[]), (c, a, rest))
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
        rows = Dict{Tuple{Int, Int}, Dict{Int, T}}()
        for r in 1:m, (x, v) in zip(q[r], cols)
            if iszero(x)
                continue
            end
            s, a, g, rest = nodes[v]
            if s == k
                if rest == 0
                    add!(get!(Combination{T}, pieces[j], (r, 0)), a, x)
                else
                    add!(get!(Dict{Int, T}, rows, (r, a)), rest, x)
                end
            else
                add!(get!(Dict{Int, T}, rows, (r, strings(k, g))), v, x)
            end
        end
        for (c, a, rest) in get(starts, k, Tuple{T, Int, Int}[])
            add!(get!(Dict{Int, T}, rows, (0, a)), rest, c)
        end
        # the channels continued first, in order, then the terms beginning here
        rowkeys = sort!(collect(keys(rows)); by = x -> (x[1] == 0 ? m + 1 : x[1], x[2]))
        newcols, grows = dense_rows([ rows[key] for key in rowkeys ])
        # a row continuing a channel is measured against all the channel holds, so that a
        # remainder of rounding is not taken for a channel of its own; the pieces of the site
        # hold only the closings yet
        content = zeros(real(T), m)
        for ((r, _), comb) in pieces[j]
            content[r] += sum(abs2, values(comb))
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
                add!(get!(Combination{T}, pieces[j], (r, p)), a, t[i, p])
            end
        end
        q = grows[kept]
        cols = newcols
        dims[j] = length(kept)
    end
    if pop!(dims) ≠ 0
        error("bug: channels left open on the last site of the terms")
    end
    return first_site, dims, pieces
end

"""
    on_basis!(b, tol)

write the pieces of the block `b`, and return it, on a basis of each site: the identity and the
parts of the atoms of the site orthogonal to it, as few as they span. Each atom being the
identity times its trace over the dimension plus that part, the identity is taken out of every
piece as far as it goes, which `push_identities!` then moves out of the openings and the
closings. A coefficient under `tol` times its whole combination is dropped.
"""
function on_basis!(b::ComBlock{T}, tol::Real) where T
    for (j, atoms) in enumerate(b.atoms)
        if isempty(atoms)
            continue
        end
        d = size(atoms[1], 1)
        τ = [ tr(m) / d for m in atoms ]
        parts = [ m - x * I for (m, x) in zip(atoms, τ) ]
        anorms = norm.(atoms)
        kept, t = interpolative(vec.(parts), anorms, tol)
        basis = parts[kept]
        norms = [ sqrt(d); norm.(basis) ]
        written = SitePieces{T}()
        for (lr, comb) in b.pieces[j]
            out = Combination{T}()
            for (a, x) in sort!(collect(comb); by = first)
                if a == 0
                    add!(out, 0, x)
                    continue
                end
                if !iszero(τ[a])
                    add!(out, 0, x * τ[a])
                end
                for p in eachindex(kept)
                    if !iszero(t[a, p])
                        add!(out, p, x * t[a, p])
                    end
                end
            end
            scale = sqrt(sum(abs2(x) * (a == 0 ? d : anorms[a]^2) for (a, x) in comb; init = 0.0))
            filter!(p -> abs(last(p)) * norms[first(p) + 1] > tol * scale, out)
            if !isempty(out)
                written[lr] = out
            end
        end
        b.pieces[j] = written
        b.atoms[j] = basis
    end
    return b
end

"""
    push_identities!(b)

take the identity out of the openings of the block `b`, from left to right, and return the
terms of one site this leaves, a combination for each site of `b`: `α` times the identity
opening channel `r` becomes `α` times each piece leaving `r` on the next site, an opening
there, which may hold the identity in turn, or a term of one site when the piece closes.
"""
function push_identities!(b::ComBlock{T}) where T
    singles = [ Combination{T}() for _ in b.pieces ]
    for j in 1:length(b.pieces) - 1
        # the identity of each opening, taken out
        opened = Dict{Int, T}()
        for ((l, r), comb) in b.pieces[j]
            if l == 0 && haskey(comb, 0)
                opened[r] = pop!(comb, 0)
            end
        end
        filter!(p -> !isempty(last(p)), b.pieces[j])
        for ((l, s), comb) in sort!(collect(b.pieces[j + 1]); by = first)
            α = get(opened, l, nothing)
            if isnothing(α)
                continue
            end
            to = s == 0 ? singles[j + 1] : get!(Combination{T}, b.pieces[j + 1], (0, s))
            for (a, x) in comb
                add!(to, a, α * x)
            end
        end
    end
    return singles
end

"""
    left_sweep!(b, tol)

reduce the block `b`, from left to right, to channels whose left operators are independent,
and return it, a channel that is a combination of the others being merged into them on the
next site
"""
function left_sweep!(b::ComBlock{T}, tol::Real) where T
    for j in 1:length(b.pieces) - 1
        m = b.dims[j]
        entries = [ Dict{Tuple{Int, Int}, T}() for _ in 1:m ]
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
                    d = get!(Combination{T}, next, (p, r))
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

the block `b` read from right to left, the piece `(l, r)` of a site becoming `(r, l)`, its
combinations shared. `start` is kept, which `push_identities!` and `left_sweep!` do not read, so
that mirroring twice gives `b` back.
"""
mirror(b::ComBlock{T}) where T =
    ComBlock{T}(b.start, reverse(b.dims),
                [ SitePieces{T}((r, l) => c for ((l, r), c) in p) for p in reverse(b.pieces) ],
                reverse(b.atoms))

"""
    combination_matrix(b, j, comb, d)

the matrix of the combination `comb` of the atoms of the `j`-th site of the block `b`, a site
of dimension `d`
"""
function combination_matrix(b::ComBlock{T}, j::Int, comb::Combination{T}, d::Int) where T
    m = zeros(T, d, d)
    for (a, x) in comb
        if a == 0
            m += x * I
        else
            m += x * b.atoms[j][a]
        end
    end
    return m
end

"""
    normalized_norm(m)

the Hilbert-Schmidt norm of the matrix `m` over that of the identity, 1 for any unitary, so
that the norm of a tensor product, the product of those of its factors, does not grow with
their number
"""
normalized_norm(m::AbstractMatrix) = norm(m) / sqrt(size(m, 1))

"""
    string_factors(R, system, t)

the factors of the term `t`, as `(site, matrix, g)`: its blocks with its strings laid on them,
see `laid_strings`, and `g` the strings it has on the sites between the factor before and this
one, see `strings_on`, the same all the way since a string changes only on a block or on site
1, and 0 when none of those sites is fermionic. A string from site 1 takes a factor of its own
on the first fermionic site.
"""
function string_factors(::Type{R}, system::System, t::Term) where R
    fs = Tuple{Int, Matrix, Int}[]
    k = next_fermion(system, 1)
    if !isempty(t.blocks) && k < first(first(t.blocks)) && strings_on(t, k) ≠ 0
        push!(fs, (k, string_matrix(R, matrix(F, system[k]), strings_on(t, k)), 0))
    end
    for (k, m) in t.blocks
        g = 0
        if !isempty(fs)
            j = first(last(fs))
            g = next_fermion(system, j + 1) < k ? strings_on(t, j) : 0
        end
        f, s = matrix(F, system[k]), strings_on(t, k)
        push!(fs, (k, f == I || s == 0 ? m : m * string_matrix(R, f, s), g))
    end
    return fs
end

"""
    compacted(R, system, terms)

the terms of an operator on `system` as its constant, its terms of one site, a matrix for each
site, and the block of channels of its terms of several sites, `nothing` when it has none, the
matrices being those of superoperators for `R = Mixed`. The strings of the terms are laid on
their factors, see `string_factors`, and carried by the channels between them. The block, built
by `direct_pass`, is written on a basis of each site, see `on_basis!`, the identity moved out of
its openings and its closings into terms of one site, see `push_identities!`, and the sweeps,
`left_sweep!` from the left then on its mirror from the right, reduce what these relations made
dependent: each link is left with as few channels as the rank of the coupling across it. A term
of one site that rounding left in place of zero is dropped, the terms measured by
`normalized_norm`.
"""
function compacted(::Type{R}, system::System, terms::Vector{Term}) where R
    scale = maximum(t -> abs(t.coef) * prod(normalized_norm(m) for (_, m) in t.blocks; init = 1.0),
                    terms; init = 0.0)
    T = promote_type(Float64, (typeof(t.coef) for t in terms)...,
                     (eltype(m) for t in terms for (_, m) in t.blocks)...)
    constant = zero(T)
    singles = Dict{Int, Matrix{T}}()
    atoms = Dict{Int, Vector{Matrix{T}}}()
    atom_ids = Dict{Tuple{Int, Matrix}, Int}()
    # the number of the atom of the matrix `m` on site `k`, from 1
    atom(k, m) = get!(atom_ids, (k, m)) do
        push!(get!(atoms, k, Matrix{T}[]), T.(m))
        length(atoms[k])
    end
    # the atom of the strings `g` on site `k`, where a term crosses it without a factor
    string_ids = Dict{Tuple{Int, Int}, Int}()
    string_atom(k, g) = g == 0 || next_fermion(system, k) ≠ k ? 0 :
        get!(() -> atom(k, string_matrix(R, matrix(F, system[k]), g)), string_ids, (k, g))
    multi = Tuple{T, Vector{Tuple{Int, Int, Int}}}[]
    for t in terms
        fs = string_factors(R, system, t)
        if isempty(fs)
            constant += t.coef
        elseif length(fs) == 1
            k, m, _ = only(fs)
            add!(singles, k, T.(t.coef * m))
        else
            push!(multi, (T(t.coef), [ (k, atom(k, m), g) for (k, m, g) in fs ]))
        end
    end
    b = nothing
    if !isempty(multi)
        start, dims, pieces = direct_pass(multi, string_atom, rounding_tol)
        # the atoms read after the pass, which adds those of the strings
        b = ComBlock{T}(start, dims, pieces,
                        [ get(atoms, k, Matrix{T}[]) for k in start:start + length(pieces) - 1 ])
        b = on_basis!(b, rounding_tol)
        opened = push_identities!(b)
        b = mirror(b)
        closed = push_identities!(b)
        b = left_sweep!(mirror(b), rounding_tol)
        b = mirror(left_sweep!(mirror(b), rounding_tol))
        for (j, comb) in enumerate(mergewith(+, s, c) for (s, c) in zip(opened, reverse(closed)))
            if !isempty(comb)
                k = b.start + j - 1
                add!(singles, k, combination_matrix(b, j, comb, block_dim(R, system[k])))
            end
        end
    end
    return constant, filter(p -> normalized_norm(last(p)) > rounding_tol * scale, singles), b
end
