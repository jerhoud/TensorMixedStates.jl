# The MPO of an operator: PreMPO, which gathers the terms of a simplified operator, and
# make_mpo, which builds the MPO from it, with make_approx_W1 and make_approx_W2 for the WI and
# WII approximations of an exponential.

export PreMPO, make_mpo, make_approx_W1, make_approx_W2

struct PreMPO{R <: PM}
    system::System
    linkdims::Vector{Int}
    terms::Vector{Vector{Tuple{Int, Int, ITensor, Int}}}
    # the number of time functions the terms take, one per element of a time dependent
    # evolver: it cannot be counted from `terms`, a whole element vanishing on its sites
    nterms::Int
    function PreMPO{R}(system::System, nterms::Int = 1) where R
        n = length(system)
        return new{R}(system, fill(1, n - 1), [ Tuple{Int, Int, ITensor, Int}[] for _ in 1:n ], nterms)
    end
end

"""
    check_coefs(pre, coefs)

refuse `coefs` unless it holds one value per time function of `pre`.
"""
function check_coefs(pre::PreMPO, coefs)
    if length(coefs) ≠ pre.nterms
        error("an evolver of $(pre.nterms) terms takes as many time functions, got $(length(coefs))")
    end
end

"""
    PreMPO!(pre, coef, factors[, ref])
    PreMPO!(pre, op[, ref])
    PreMPO!(pre, ops)

add to `pre`, and return it, the term `coef` times the product of the one site `factors`, or
the terms of the simplified operator `op`, `ref` numbering their time function. A term of
several sites takes a channel of its own on every link it spans. Each operator of the vector
`ops` gets the time function of its position.
"""
function PreMPO!(pre::PreMPO{R}, coef::Number, subs::Vector{<:IndexedOp{R}}, ref::Int=1) where R
    foreach(o -> check_one_site(o, "an MPO"), subs)
    sys = pre.system
    # the identity has no site, and a term made of it alone is laid on the first one
    subs = filter(o -> !(o isa IdentityOp), subs)
    if isempty(subs)
        kdx = SysIndex{R}(sys, 1)
        push!(pre.terms[1], (1, 1, coef * delta(kdx', dag(kdx)), ref))
        return pre
    end
    # a term whose factor vanishes on its site, as C(1)*C(1) or Sp(1)*Sp(1) on a spin 1/2, is
    # dropped here, where the sites are known: simplify cannot tell, one name standing for
    # operators of different algebras on different sites. Kept, it took a channel on every
    # link it spans, and on a charged system its tensor has no block, hence no flux
    if any(o -> iszero(tensor(sys, o)), subs)
        return pre
    end
    ld = pre.linkdims
    tm = pre.terms
    fst = subs[1].index[1]
    lst = subs[end].index[1]
    if fst == lst
        push!(tm[fst], (1, 1, coef * tensor(sys, subs[1]), ref))
    else
        for k in fst:lst-1
            ld[k] += 1
        end    
        push!(tm[fst], (1, ld[fst], coef * tensor(sys, subs[1]), ref))
        i = fst
        for ind in subs[2:end-1]
            j = ind.index[1]
            for k in i+1:j-1
                kdx = SysIndex{R}(sys, k)
                push!(tm[k],(ld[k-1], ld[k], delta(kdx', dag(kdx)), ref))
            end
            push!(tm[j], (ld[j-1], ld[j], tensor(sys, ind), ref))
            i = j
        end
        for k in i+1:lst-1
            kdx = SysIndex{R}(sys, k)
            push!(tm[k],(ld[k-1], ld[k], delta(kdx', dag(kdx)), ref))
        end
        push!(tm[lst], (ld[lst-1], 1, tensor(sys, subs[end]), ref))
    end
    return pre
end

PreMPO!(pre::PreMPO{R}, a::IndexedOp{R}, ref::Int = 1) where {R <: PM} =
    PreMPO!(pre, scalarcoef(a), prodsubs(a), ref)

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
by `make_mpo`, `make_approx_W1` or `make_approx_W2`, or passed to `tdvp` or `approx_W` in
place of the operator, which saves preprocessing it again. `op` may also be a vector of
operators, the terms of a time dependent evolver, each multiplied by its own time function.

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
    return PreMPO!(PreMPO{R}(state.system, n), removeMulti(simplify(adapt_representation(R, a))))
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
    mpo_charges(pre, coefs)

the charge of every channel of every link of the MPO, `q[i + 1][k]` for channel `k` of the
link on the right of site `i`.

A channel stands for a term partly placed, the sites on its left having contributed their
factors, and its charge is minus the sum of their fluxes. The first channel is the term not
yet begun and has no charge; the last is the term finished and has minus the flux of the whole
operator, which every term must share or the operator is refused. Terms whose time function
is zero are left out.
"""
function mpo_charges(pre::PreMPO{R}, coefs) where R
    n = length(pre.system)
    ld = pre.linkdims
    tm = pre.terms
    rdims = [ i == n ? 1 : ld[i] for i in 1:n ]
    q = [ Vector{Union{Nothing, QN}}(nothing, 1 + (i == 0 ? 1 : rdims[i])) for i in 0:n ]
    q[1][1] = QN()
    total = nothing
    for i in 1:n
        q[i+1][1] = QN()
        for (l, r, u, ref) in tm[i]
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
                    st = strong_names(pre.system)
                    only_strong = !isempty(st) &&
                        weak_qn(total, st, String[]) == weak_qn(c, st, String[])
                    error("the terms of this operator do not all carry the same charge, " *
                          "$(-total) and $(-c), " *
                          (only_strong ?
                              "which only a strong symmetry tells apart: drop `strong` or " *
                              "weaken the state" :
                              "so it has no definite flux and cannot be put on a system " *
                              "that conserves it"))
                end
                total = c
            else
                q[i+1][r] = c
            end
        end
    end
    if isnothing(total)
        total = QN()
    end
    for i in 0:n
        q[i+1][end] = total
        for k in eachindex(q[i+1])
            if isnothing(q[i+1][k])
                # a channel no term goes through, kept so that the numbering does not move
                q[i+1][k] = total
            end
        end
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
    link_maker(pre, coefs, charges)

a function of `(i, d)` giving the link on the right of site `i` with `d` channels. On a
charged system the channels carry the charges `charges(pre, coefs)` gives, `charges` being
`mpo_charges` or `w_charges`; otherwise the link is a plain index.
"""
function link_maker(pre::PreMPO, coefs, charges)
    if !is_charged(pre.system)
        return (i, d) -> Index(d, "Link,l=$i")
    end
    q = charges(pre, coefs)
    return (i, d) -> Index([ q[i+1][k] => 1 for k in 1:d ]...; tags = "Link,l=$i")
end

"""
    close_end(w, link, k)

the tensor `w` of the first or the last site, with its outer link fixed on channel `k`.
"""
close_end(w::ITensor, link::Index, k::Int) = w * onehot(link => k)

"""
    add_block!(w, llink, l, rlink, r, u, idx[, c])

add `c` times the one site operator `u` to the block of channels `l` and `r` of `w`, the
tensor of the site of index `idx`, and return `w`. Zeros are skipped: a block sparse tensor
refuses an element outside its flux, even a zero.
"""
function add_block!(w::ITensor, llink::Index, l::Int, rlink::Index, r::Int, u::ITensor,
                    idx::Index, c::Number = 1)
    for j in eachindval(idx, idx')
        v = c * u[j...]
        if !iszero(v)
            w[llink => l, rlink => r, j...] += v
        end
    end
    return w
end

"""
    make_mpo(::PreMPO[, coefs])
    make_mpo(::State, op)

the MPO of an operator, in the representation of the state. For a time dependent evolver,
`coefs` holds the value of each time function, one per term; it defaults to `[1.]`, a single
operator.

A pure operator ``A`` given for a mixed state becomes its `Evolver`,
``\\rho \\mapsto A \\rho + \\rho A^\\dagger``, see `PreMPO`, and not the gate
``\\rho \\mapsto A \\rho A^\\dagger`` that `apply` makes of it.

# Examples

    mpo = make_mpo(state, sum(Z(i) * Z(i + 1) for i in 1:9))
"""
function make_mpo(pre::PreMPO{R}, coefs=[1.]) where R
    check_coefs(pre, coefs)
    sys = pre.system
    ld = pre.linkdims
    tm = pre.terms
    n = length(sys)
    ts = Vector{ITensor}(undef, n)
    elt = mpo_eltype(pre, coefs)
    # the links of a charged MPO carry the charge each channel has accumulated, without which
    # the tensor of a site would hold several fluxes
    mklink = link_maker(pre, coefs, mpo_charges)
    rdim = 1
    rlink = mklink(0, 2)
    for i in 1:n
        idx = SysIndex{R}(sys, i)
        ldim = rdim
        llink = rlink
        if i == n
            rdim = 1
        else
            rdim = ld[i]
        end
        rlink = mklink(i, 1 + rdim)
        w = ITensor(elt, idx', dag(idx), dag(llink), rlink)
        id = delta(dag(idx), idx')
        add_block!(w, llink, 1, rlink, 1, id, idx)
        add_block!(w, llink, 1 + ldim, rlink, 1 + rdim, id, idx)
        for (l, r, u, ref) in tm[i]
            c = coefs[ref]
            if c ≠ 0
                # the coefficient of a term goes on its closing piece alone, the only one
                # with r == 1: laid on every piece, a term of k sites took it to the power k
                if r == 1
                    add_block!(w, llink, l, rlink, r + rdim, u, idx, c)
                else
                    add_block!(w, llink, l, rlink, r, u, idx)
                end
            end
        end
        if i == 1
            w = close_end(w, llink, 1)
        end
        if i == n
            w = close_end(w, dag(rlink), 2)
        end
        ts[i] = w
    end
    return MPO(ts)
end

make_mpo(state::State, a) = make_mpo(PreMPO(state, a))

"""
    make_approx_W1(::PreMPO, tau[, coefs])
    make_approx_W1(::State, op, tau)

the MPO of the approximation WI of the exponential of `tau` times the operator, `coefs` being
as for `make_mpo`. On a charged system the operator must have zero flux.
"""
function make_approx_W1(pre::PreMPO{R}, tau::Number, coefs=[1.]) where R
    check_coefs(pre, coefs)
    sys = pre.system
    ld = pre.linkdims
    tm = pre.terms
    n = length(sys)
    ts = Vector{ITensor}(undef, n)
    elt = promote_type(mpo_eltype(pre, coefs), typeof(tau))
    mklink = link_maker(pre, coefs, w_charges)
    rdim = 1
    rlink = mklink(0, 1)
    for i in 1:n
        idx = SysIndex{R}(sys, i)
        llink = rlink
        if i == n
            rdim = 1
        else
            rdim = ld[i]
        end
        rlink = mklink(i, rdim)
        w = ITensor(elt, idx', dag(idx), dag(llink), rlink)
        add_block!(w, llink, 1, rlink, 1, delta(dag(idx), idx'), idx)
        for (l, r, u, ref) in tm[i]
            c = coefs[ref]
            if c ≠ 0
                # the coefficient and the time step go on the closing piece of a term alone,
                # as in make_mpo
                add_block!(w, llink, l, rlink, r, u, idx, r == 1 ? c * tau : one(c))
            end
        end
        if i == 1
            w = close_end(w, llink, 1)
        end
        if i == n
            w = close_end(w, dag(rlink), 1)
        end
        ts[i] = w
    end
    return MPO(ts)
end

make_approx_W1(state::State, a, tau::Number) = make_approx_W1(PreMPO(state, a), tau)

"""
    make_approx_W2(::PreMPO, tau[, coefs])
    make_approx_W2(::State, op, tau)

the MPO of the approximation WII of the exponential of `tau` times the operator, `coefs` being
as for `make_mpo`. On a charged system the operator must have zero flux.
"""
function make_approx_W2(pre::PreMPO{R}, tau::Number, coefs=[1.]) where R
    check_coefs(pre, coefs)
    sys = pre.system
    ld = pre.linkdims
    tm = pre.terms
    n = length(sys)
    ts = Vector{ITensor}(undef, n)
    elt = promote_type(mpo_eltype(pre, coefs), typeof(tau))
    mklink = link_maker(pre, coefs, w_charges)
    rdim = 1
    rlink = mklink(0, 1)
    for i in 1:n 
        idx = SysIndex{R}(sys, i)
        ldim = rdim
        if i == n
            rdim = 1
        else
            rdim = ld[i]
        end
        llink = rlink
        rlink = mklink(i, rdim)
        v = fill(ITensor(), (ldim, rdim))
        for (l, r, u, ref) in tm[i]
            c = coefs[ref]
            if c ≠ 0
                # the coefficient and the time step go on the closing piece of a term alone,
                # as in make_mpo
                v[l, r] += (r == 1 ? c * tau : one(c)) * u
            end
        end
        d = v[1, 1]
        if isempty(d)
            e = delta(dag(idx), idx')
        else
            e = exp(d)
        end
        v[1, 1] = e
        for l in 2:ldim, r in 2:rdim
            vl = v[l, 1]
            vr = v[1, r]
            if !isempty(vl) && !isempty(vr)
                v[l, r] += replaceprime(vl'' * e' * vr, 3=>1)
            end
        end
        for r in 2:rdim
            vr = v[1, r]
            if !isempty(vr)
                v[1, r] = replaceprime(e' * vr, 2=>1)
            end
        end
        for l in 2:ldim
            vl = v[l, 1]
            if !isempty(vl)
                v[l, 1] = replaceprime(vl' * e, 2=>1)
            end
        end
        
        w = ITensor(elt, idx', dag(idx), dag(llink), rlink)
        for l in 1:ldim, r in 1:rdim
            if !isempty(v[l, r])
                add_block!(w, llink, l, rlink, r, v[l, r], idx)
            end
        end
        if i == 1
            w = close_end(w, llink, 1)
        end
        if i == n
            w = close_end(w, dag(rlink), 1)
        end
        ts[i] = w
    end
    return MPO(ts)
end

make_approx_W2(state::State, a, tau::Number) = make_approx_W2(PreMPO(state, a), tau)
    