# The MPO of an operator: PreMPO, which gathers the terms of a simplified operator, and
# make_mpo, which builds the MPO from it, with make_approx_W1 and make_approx_W2 for the WI and
# WII approximations of an exponential.

export PreMPO, make_mpo, make_approx_W1, make_approx_W2

struct PreMPO{R <: PM}
    system::System
    linkdims::Vector{Int}
    # the pieces on each site, see `add_term!`
    terms::Vector{Vector{Tuple{Int, Int, ITensor, Int, Matrix}}}
    # the number of time functions the terms take, one per element of a time dependent
    # evolver: it cannot be counted from `terms`, a whole element vanishing on its sites
    nterms::Int
    function PreMPO{R}(system::System, nterms::Int = 1) where R
        n = length(system)
        return new{R}(system, fill(1, n - 1),
                      [ Tuple{Int, Int, ITensor, Int, Matrix}[] for _ in 1:n ], nterms)
    end
end

"""
    site_matrix(u, idx)

the matrix of the one site operator `u` on the site of index `idx`, rows on `idx'`. It is read
through `dense`, which takes a delta or a block sparse tensor as any other: with the NDTensors
of ITensors 0.7, `Array` fails on the delta of a charged site and `denseblocks` on a dense
tensor.
"""
site_matrix(u::ITensor, idx::Index) = Array(dense(u), idx', dag(idx))

"""
    add_term!(pre, k, l, r, u, ref)

add to site `k` of `pre` the piece `u`, going from channel `l` to channel `r` and taking the
time function `ref`, with its matrix on the index of the site, read once here rather than at
every MPO built from `pre`
"""
add_term!(pre::PreMPO{R}, k::Int, l::Int, r::Int, u::ITensor, ref::Int) where R =
    push!(pre.terms[k], (l, r, u, ref, site_matrix(u, SysIndex{R}(pre.system, k))))

"""
    check_coefs(pre, coefs)

refuse `coefs` unless it holds one real value per time function of `pre`. A pure term lifted
for a mixed state, ``A \\rho + \\rho A^\\dagger``, would need a complex value on one side and
its conjugate on the other, so a complex function is written as its real and imaginary parts,
each with its own term, on pure states as well, which keeps one interface for both.
"""
function check_coefs(pre::PreMPO, coefs)
    if length(coefs) ≠ pre.nterms
        error("an evolver of $(pre.nterms) terms takes as many time functions, got $(length(coefs))")
    end
    for (k, c) in enumerate(coefs)
        if !(c isa Real)
            error("time functions take real values, and term $k got $c: " *
                  "write f * A as real(f) * A + imag(f) * (im * A)")
        end
    end
end

"""
    com_tensor(system, o, k)

the tensor on site `k` of `system` of the piece `o` of a com, an operator of one site times a
coefficient, built from its matrix as `tensor` builds that of a placed operator: placing it is
not possible for the identity, which has no site once placed.
"""
com_tensor(sys::System, o::GenericOp{R, 1}, k::Int) where R =
    scalarcoef(o) * legs_on(scalararg(o), [sys[k]], [SysIndex{Pure}(sys, k)], [SysIndex{R}(sys, k)])

"""
    PreMPO!(pre, coef, term[, ref])
    PreMPO!(pre, coef, com[, ref])
    PreMPO!(pre, op[, ref])
    PreMPO!(pre, ops)

add to `pre`, and return it, `coef` times the constant or the term of one site `term`, the
com `com` times `coef`, or the terms of the simplified and compacted operator `op`, `ref`
numbering their time function: the terms of several sites come gathered in coms, see
`compact_simplified`, and a com takes as few channels as the operators of its sites allow, see
`reduce_on_sites`. Each operator of the vector `ops` gets the time function of its position.
"""
function PreMPO!(pre::PreMPO{R}, coef::Number, a::IndexedOp{R}, ref::Int = 1) where R
    sys = pre.system
    # the identity has no site, and a term made of it alone is laid on the first one
    subs = filter(o -> !(o isa IdentityOp), prodsubs(a))
    if isempty(subs)
        kdx = SysIndex{R}(sys, 1)
        add_term!(pre, 1, 1, 1, coef * delta(kdx', dag(kdx)), ref)
        return pre
    end
    if length(subs) > 1
        error("bug: the term $a of several sites reaches PreMPO! without being compacted")
    end
    o = only(subs)
    u = tensor(sys, o)
    # a term whose factor vanishes on its site, as C(1)*C(1) or Sp(1)*Sp(1) on a spin 1/2, is
    # dropped here, where the sites are known: simplify cannot tell, one name standing for
    # operators of different algebras on different sites. Kept, on a charged system its
    # tensor would have no block, hence no flux
    if !iszero(u)
        add_term!(pre, only(o.index), 1, 1, coef * u, ref)
    end
    return pre
end

function PreMPO!(pre::PreMPO{R}, coef::Number, a::ComOp{R}, ref::Int = 1) where R
    sys = pre.system
    b, singles = reduce_on_sites(R, a, coef, sys)
    ks = b.start .+ (0:length(b.pieces) - 1)
    ld = pre.linkdims
    # the channels of the com numbered after those its links already have
    offs = [ ld[k] for k in ks[1:end-1] ]
    for (j, k) in enumerate(ks[1:end-1])
        ld[k] += b.dims[j]
    end
    for (j, k) in enumerate(ks)
        laid = sort!(collect(b.pieces[j]); by = first)
        if !isempty(singles[j])
            push!(laid, (0, 0) => singles[j])
        end
        for ((l, r), comb) in laid
            o = simplify_sum(GenericOp{R, 1}[ x * atom for (atom, x) in comb ])
            if scalarcoef(o) ≠ 0
                add_term!(pre, k, l == 0 ? 1 : offs[j-1] + l, r == 0 ? 1 : offs[j] + r,
                          com_tensor(sys, o, k), ref)
            end
        end
    end
    return pre
end

PreMPO!(pre::PreMPO{R}, a::IndexedOp{R}, ref::Int = 1) where {R <: PM} =
    PreMPO!(pre, scalarcoef(a), scalararg(a), ref)

function PreMPO!(pre::PreMPO{R}, s::SumOp{R, Indexed}, ref::Int=1) where {R <: PM}
    for p in s.subs
        PreMPO!(pre, p, ref)
    end
    return pre
end

function PreMPO!(pre::PreMPO, as)
    for (i, a) in enumerate(as)
        PreMPO!(pre, a, i)
    end
    return pre
end

"""
    adapt_representation(::Type{R}, op)

the operator `op` adapted to the representation `R` of the MPO. A pure operator for a mixed
state is an evolver, `-im * H`, and is lifted with `Evolver`; a mixed one for a pure state is
refused. The terms of a time dependent evolver, a vector, are adapted one by one.
"""
adapt_representation(::Type{Pure}, a::IndexedOp{Mixed}) =
    error("cannot build a pure MPO from the mixed operator $a, " *
          "the state must be in mixed representation (see ToMixed)")
adapt_representation(::Type{Mixed}, a::IndexedOp{Pure}) = Evolver(a)
adapt_representation(::Type{R}, a::Vector) where R = map(x -> adapt_representation(R, x), a)
adapt_representation(::Type{R}, a) where R = a

"""
    PreMPO(::State, op)

the operator `op` preprocessed for the representation of the state, to be turned into an MPO
by `make_mpo`, `make_approx_W1` or `make_approx_W2`, or passed to `tdvp`, `approx_W`, `dmrg`
or `steady_state` in place of the operator, which saves preprocessing it again: a phase of
one's own evolving one step at a time prepares its evolver once. `op` may also be a vector of
operators, the terms of a time dependent evolver, each multiplied by its own real time function.
The terms of several sites are compacted, see `compact`, so that the bond dimension of the MPO
is the least any triangular MPO of the operator can have.

A pure operator ``A`` given for a mixed state is lifted with `Evolver` to
``\\rho \\mapsto A \\rho + \\rho A^\\dagger``: the hamiltonian part of an evolver must already be
written `-im * H`.

# Examples

    pre = PreMPO(state, [sum(X(i) for i in 1:10), sum(Z(i) for i in 1:10)])
    mpo = make_mpo(pre, [1., 0.5])
"""
function PreMPO(state::State{R}, a) where R
    # on the operator as it was written, so that the message names what the caller wrote
    # and not what `simplify` made of it
    check_indices(state.system, a)
    n = a isa Vector ? length(a) : 1
    s = removeMulti(simplify(adapt_representation(R, a)))
    gather(x) = compact_simplified(x, rounding_tol, "an MPO")
    return PreMPO!(PreMPO{R}(state.system, n), s isa Vector ? map(gather, s) : gather(s))
end

"""
    mpo_eltype(::PreMPO, coefs)

the element type of the tensors of the MPO built from `pre` and `coefs`, `Float64` at least.
Real matrices with real coefficients give a real MPO, whose contractions are cheaper. The
approximations WI and WII promote it with the type of their time step.
"""
mpo_eltype(pre::PreMPO, coefs) =
    promote_type(
        Float64,
        mapreduce(typeof, promote_type, coefs; init = Bool),
        mapreduce(t -> eltype(t[3]), promote_type, Iterators.flatten(pre.terms); init = Bool))

"""
    charge_mismatch(system, a, b, shown)

raise the error for an operator on `system` whose terms do not all carry the same charge, `a`
and `b` being the charges of two of them, printed when `shown`, or those of two partial terms
entering the same channel of a com, whose values would mean nothing to the reader.
"""
function charge_mismatch(sys::System, a::QN, b::QN, shown::Bool)
    st = strong_names(sys)
    only_strong = !isempty(st) && weak_qn(a, st, String[]) == weak_qn(b, st, String[])
    error("the terms of this operator do not all carry the same charge, " *
          (shown ? "$(-a) and $(-b), " : "") *
          (only_strong ?
              "which only a strong symmetry tells apart: drop `strong` or weaken the state" :
              "so it has no definite flux and cannot be put on a system that conserves it"))
end

"""
    channel_counts(pre)

the number of channels the terms of `pre` take on each link, from the one on the left of the
first site to the one on the right of the last, the term not yet begun included.
"""
channel_counts(pre::PreMPO) = [ 1; pre.linkdims; 1 ]

"""
    mpo_charges(pre, coefs)

the charge of every channel of every link of the MPO, `q[i + 1][k]` for channel `k` of the
link on the right of site `i`.

A channel stands for a term partly placed, the sites on its left having contributed their
factors, and its charge is minus the sum of their fluxes. The first channel is the term not
yet begun and has no charge; the last is the term finished and has minus the flux of the whole
operator, which every term must share or the operator is refused. A channel of a com is
entered by several pieces, which must agree. Terms whose time function is zero are left out.
"""
function mpo_charges(pre::PreMPO{R}, coefs) where R
    n = length(pre.system)
    tm = pre.terms
    q = [ Union{Nothing, QN}[ QN(); fill(nothing, d) ] for d in channel_counts(pre) ]
    total = nothing
    for i in 1:n
        for (l, r, u, ref, _) in tm[i]
            if coefs[ref] == 0
                continue
            end
            left = q[i][l]
            if isnothing(left)
                error("bug: channel $l of link $(i-1) has no charge")
            end
            # the link loses what the operator carries, which is the sign that makes the
            # tensor of the site a flux of zero once its two ends are oriented
            c = left - flux(u)
            if r == 1
                if !isnothing(total) && total ≠ c
                    charge_mismatch(pre.system, total, c, true)
                end
                total = c
            else
                if !isnothing(q[i+1][r]) && q[i+1][r] ≠ c
                    charge_mismatch(pre.system, q[i+1][r], c, false)
                end
                q[i+1][r] = c
            end
        end
    end
    total = something(total, QN())
    for qi in q
        qi[end] = total
        # a channel no term goes through, kept so that the numbering does not move
        replace!(qi, nothing => total)
    end
    return q
end

"""
    w_charges(pre, coefs)

the charges of the links of the approximations WI and WII, where one channel serves both the
term not yet begun and the term finished. The two must then have the same charge, so the
operator must have zero flux, as does any generator of an evolution keeping the state in its
sector.
"""
function w_charges(pre::PreMPO, coefs)
    q = mpo_charges(pre, coefs)
    total = q[1][end]
    if total ≠ QN()
        error("the approximations WI and WII need an operator of zero flux, and this one " *
              "carries $(-total), so it would move the charge its system conserves")
    end
    return q
end

"""
    mpo_links(pre, coefs, charges, extra)

the links of the MPO of `pre`, from the one on the left of the first site to the one on the
right of the last, each with `extra` channels besides those of `channel_counts`: one, for the
term finished, in `make_mpo`, and none in WI and WII, where the term not yet begun serves as
finished. On a charged system the channels carry the charges `charges(pre, coefs)` gives,
`charges` being `mpo_charges` or `w_charges`; otherwise the links are plain indices.
"""
function mpo_links(pre::PreMPO, coefs, charges, extra::Int)
    ds = channel_counts(pre) .+ extra
    if !is_charged(pre.system)
        return [ Index(d, "Link,l=$(i-1)") for (i, d) in enumerate(ds) ]
    end
    q = charges(pre, coefs)
    return [ Index([ q[i][k] => 1 for k in 1:d ]...; tags = "Link,l=$(i-1)")
             for (i, d) in enumerate(ds) ]
end

"""
    close_ends(ts, links, k)

the MPO of the tensors `ts` of the sites, its outer links `links[1]` and `links[end]` fixed on
channel 1, the term not yet begun, and on channel `k`, the term finished.
"""
function close_ends(ts::Vector{ITensor}, links::Vector{<:Index}, k::Int)
    ts[1] *= onehot(links[1] => 1)
    ts[end] *= onehot(dag(links[end]) => k)
    return MPO(ts)
end

"""
    site_array(elt, idx, llink, rlink)

the array, zero, of the tensor of the site of index `idx` in an MPO of links `llink` and
`rlink`, its axes those of `idx'`, `dag(idx)`, `llink` and `rlink`: filled block by block with
`add_block!`, it is made a tensor at once by `site_tensor`
"""
site_array(elt::Type, idx::Index, llink::Index, rlink::Index) =
    zeros(elt, dim(idx), dim(idx), dim(llink), dim(rlink))

"""
    add_block!(a, l, r, m[, c])

add `c` times the matrix `m` of a one site operator, rows on `idx'`, to the block of channels
`l` and `r` of the array `a` of `site_array`, and return `a`
"""
function add_block!(a::Array{<:Number, 4}, l::Int, r::Int, m::AbstractMatrix, c::Number = 1)
    @views a[:, :, l, r] .+= c .* m
    return a
end

"""
    site_tensor(a, idx, llink, rlink)

the tensor of the array `a` of `site_array`, made at once: writing the elements of a block
sparse tensor one by one made most of the cost of an MPO. On a charged system only the blocks
that are not zero are kept, and they must all have the flux the charges of the links give.
"""
site_tensor(a::Array{<:Number, 4}, idx::Index, llink::Index, rlink::Index) =
    ITensor(a, idx', dag(idx), dag(llink), rlink)

"""
    make_mpo(::PreMPO[, coefs])
    make_mpo(::State, op)

the MPO of an operator, in the representation of the state. For a time dependent evolver,
`coefs` holds the value of each time function, one real number per term; it defaults to
`[1.]`, a single operator. The form taking a `State` builds the MPO of a single operator: that of a vector of
terms is built from its `PreMPO`, with its `coefs`.

A pure operator ``A`` given for a mixed state becomes its `Evolver`,
``\\rho \\mapsto A \\rho + \\rho A^\\dagger``, see `PreMPO`, and not the gate
``\\rho \\mapsto A \\rho A^\\dagger`` that `apply` makes of it.

# Examples

    mpo = make_mpo(state, sum(Z(i) * Z(i + 1) for i in 1:9))
"""
function make_mpo(pre::PreMPO{R}, coefs=[1.]) where R
    check_coefs(pre, coefs)
    sys = pre.system
    elt = mpo_eltype(pre, coefs)
    dims = channel_counts(pre)
    # the links of a charged MPO carry the charge each channel has accumulated, without which
    # the tensor of a site would hold several fluxes
    links = mpo_links(pre, coefs, mpo_charges, 1)
    ts = map(1:length(sys)) do i
        idx = SysIndex{R}(sys, i)
        llink, rlink = links[i], links[i+1]
        a = site_array(elt, idx, llink, rlink)
        id = Matrix{elt}(I, dim(idx), dim(idx))
        add_block!(a, 1, 1, id)
        add_block!(a, 1 + dims[i], 1 + dims[i+1], id)
        for (l, r, _, ref, m) in pre.terms[i]
            c = coefs[ref]
            if c ≠ 0
                # the coefficient of a term goes on its closing piece alone, the only one
                # with r == 1: laid on every piece, a term of k sites took it to the power k
                if r == 1
                    add_block!(a, l, 1 + dims[i+1], m, c)
                else
                    add_block!(a, l, r, m)
                end
            end
        end
        return site_tensor(a, idx, llink, rlink)
    end
    return close_ends(ts, links, 2)
end

make_mpo(state::State, a) = make_mpo(PreMPO(state, a))

"""
    make_approx_W1(::PreMPO, tau[, coefs])
    make_approx_W1(::State, op, tau)

the MPO of the approximation WI of the exponential of `tau` times the operator, `coefs` being
as for `make_mpo`, and the form taking a `State` being for a single operator as well. On a
charged system the operator must have zero flux.
"""
function make_approx_W1(pre::PreMPO{R}, tau::Number, coefs=[1.]) where R
    check_coefs(pre, coefs)
    sys = pre.system
    elt = promote_type(mpo_eltype(pre, coefs), typeof(tau))
    links = mpo_links(pre, coefs, w_charges, 0)
    ts = map(1:length(sys)) do i
        idx = SysIndex{R}(sys, i)
        llink, rlink = links[i], links[i+1]
        a = site_array(elt, idx, llink, rlink)
        add_block!(a, 1, 1, Matrix{elt}(I, dim(idx), dim(idx)))
        for (l, r, _, ref, m) in pre.terms[i]
            c = coefs[ref]
            if c ≠ 0
                # the coefficient and the time step go on the closing piece of a term alone,
                # as in make_mpo
                add_block!(a, l, r, m, r == 1 ? c * tau : one(c))
            end
        end
        return site_tensor(a, idx, llink, rlink)
    end
    return close_ends(ts, links, 1)
end

make_approx_W1(state::State, a, tau::Number) = make_approx_W1(PreMPO(state, a), tau)

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
    make_approx_W2(::PreMPO, tau[, coefs])
    make_approx_W2(::State, op, tau)

the MPO of the approximation WII of the exponential of `tau` times the operator, `coefs` being
as for `make_mpo`, and the form taking a `State` being for a single operator as well. On a
charged system the operator must have zero flux.

It is the WII of Zaletel et al., Phys. Rev. B 91, 165112 (2015). It keeps every product of
terms of which no two cross the same link, the terms of one site included, wherever a term of
several sites goes through their site; its error, of order ``\\tau^2``, comes from the terms
that cross a same link.
"""
function make_approx_W2(pre::PreMPO{R}, tau::Number, coefs=[1.]) where R
    check_coefs(pre, coefs)
    sys = pre.system
    elt = promote_type(mpo_eltype(pre, coefs), typeof(tau))
    dims = channel_counts(pre)
    links = mpo_links(pre, coefs, w_charges, 0)
    ts = map(1:length(sys)) do i
        idx = SysIndex{R}(sys, i)
        llink, rlink = links[i], links[i+1]
        ldim, rdim = dims[i], dims[i+1]
        k = dim(idx)
        v = Matrix{Union{Nothing, Matrix{elt}}}(nothing, ldim, rdim)
        for (l, r, _, ref, u) in pre.terms[i]
            c = coefs[ref]
            if c ≠ 0
                # the coefficient and the time step go on the closing piece of a term alone,
                # as in make_mpo
                m = (r == 1 ? c * tau : one(c)) * u
                v[l, r] = isnothing(v[l, r]) ? m : v[l, r] + m
            end
        end
        # their equation 11: each block is read in the exponential of the terms of one site D,
        # the transport A from channel l to channel r, the closing B of l and the opening C of
        # r, each taken once at most, which puts exp(D) around every piece and takes a closing
        # and an opening on the same site in both orders. Computed on dense matrices, whose
        # zeros outside the flux of a block stay exact zeros
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

make_approx_W2(state::State, a, tau::Number) = make_approx_W2(PreMPO(state, a), tau)
    