export PreMPO, make_mpo, make_approx_W1, make_approx_W2

struct PreMPO{R <: PM}
    system::System
    linkdims::Vector{Int}
    terms::Vector{Vector{Tuple{Int, Int, ITensor, Int}}}
    function PreMPO{R}(system::System) where R
        n = length(system)
        return new{R}(system, fill(1, n - 1), [ Tuple{Int, Int, ITensor, Int}[] for _ in 1:n ])
    end
end

function PreMPO!(pre::PreMPO{R}, coef::Number, subs::Vector{<:IndexedOp{R}}, ref::Int=1) where R
    foreach(o -> check_one_site(o, "an MPO"), subs)
    sys = pre.system
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

adapt an operator to the representation the MPO is built in. A pure operator given for a
mixed state is an evolver (`-im * hamiltonian`) and is lifted with `Evolver`, which is what
makes `-im * H` work on a mixed state. A time dependent evolver is a vector of terms and
each term is lifted on its own.
"""
adapt_representation(::Type{Pure}, a::IndexedOp{Mixed}) =
    error("cannot build a pure MPO from the mixed operator $a, " *
          "the state must be in mixed representation (see ToMixed)")
adapt_representation(::Type{Mixed}, a::IndexedOp{Pure}) = Evolver(a)
adapt_representation(::Type{R}, a::Vector) where R = map(x -> adapt_representation(R, x), a)
adapt_representation(::Type{R}, a) where R = a

"""
    PreMPO(::State, op)

preprocess an operator, or a vector of operators for a time dependent evolver, in which
case each one is a term whose coefficient is given by the matching time function.
The result can be passed wherever an operator that must be turned into an MPO is expected.
The operator is first adapted to the representation of the state, see `adapt_representation`.
"""
function PreMPO(state::State{R}, a) where R
    # on the operator as it was written, so that the message names what the caller wrote
    # and not what `simplify` made of it
    check_indices(state.system, a)
    return PreMPO!(PreMPO{R}(state.system), removeMulti(simplify(adapt_representation(R, a))))
end

"""
    mpo_eltype(::PreMPO, coefs)

return the element type needed for the tensors of the MPO built from `pre` and `coefs`.
The approximations WI and WII promote this further with the type of their time step.
Real operators (built only from real matrices with real coefficients) give a real MPO,
which makes all subsequent ITensor contractions about twice as fast.
`Float64` is used as a floor so that integer or boolean data never reaches the tensors.
"""
mpo_eltype(pre::PreMPO, coefs) =
    promote_type(
        Float64,
        mapreduce(typeof, promote_type, coefs; init = Bool),
        mapreduce(t -> eltype(t[3]), promote_type, Iterators.flatten(pre.terms); init = Bool))

"""
    mpo_charges(pre, coefs)

the charge of every channel of every link of the MPO, or `nothing` on a system without
charges.

A channel stands for a term of the operator partly placed: the sites on its left have
contributed their factors and the ones on its right have not. Its charge is therefore what
those factors carry, accumulated from the left, and that is what the link has to record for
the MPO to be a tensor of definite flux. Channel one is the term not yet begun and carries
nothing; the last one is the term finished and carries the flux of the whole operator, which
every term must agree on.
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
                          "$(total) and $(c), " *
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

the charges of the links of the approximations WI and WII, which share one channel between
the term not yet begun and the term finished. That only makes sense when the two carry the
same charge, so the operator has to have a flux of zero, which an evolution generator has
anyway: one that moved the charge would not keep the state in its sector.
"""
function w_charges(pre::PreMPO, coefs)
    q = mpo_charges(pre, coefs)
    total = q[1][end]
    if total ≠ QN()
        error("the approximations WI and WII need an operator of zero flux, and this one " *
              "carries $total, so it would move the charge its system conserves")
    end
    return q
end

"""
    link_maker(pre, coefs, charges)

a function giving the link on the right of site `i` with `d` channels.

The three builders draw their links the same way and differed only in which charges they
ask for, `mpo_charges` or `w_charges`, so they share this. Without charges it is the plain
index it always was.
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

the tensor of the first or last site, with its dangling link fixed on channel `k`
"""
close_end(w::ITensor, link::Index, k::Int) = w * onehot(link => k)

"""
    make_mpo(::PreMPO[, coefs])
    make_mpo(::State, operator)

build an mpo representing an operator
"""
function make_mpo(pre::PreMPO{R}, coefs=[1.]) where R
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
        # zeros are skipped rather than written: a block sparse tensor refuses an element
        # outside its flux even when what is written there is nothing
        for j in eachindval(idx, idx')
            v = id[j...]
            if !iszero(v)
                w[llink => 1, rlink => 1, j...] = v
                w[llink => 1 + ldim, rlink => 1 + rdim, j...] = v
            end
        end
        for (l, r, u, ref) in tm[i]
            c = coefs[ref]
            if c ≠ 0
                if r == 1
                    r += rdim
                end
                for j in eachindval(idx, idx')
                    v = c * u[j...]
                    if !iszero(v)
                        w[llink=>l, rlink=>r, j...] += v
                    end
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
    make_approx_W1(::State, operator, tau)

build MPO representing approximation WI of a given operator and time step
"""
function make_approx_W1(pre::PreMPO{R}, tau::Number, coefs=[1.]) where R
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
        id = delta(dag(idx), idx')
        for j in eachindval(idx, idx')
            v = id[j...]
            if !iszero(v)
                w[llink=>1, rlink=>1, j...] = v
            end
        end
        for (l, r, u, ref) in tm[i]
            c = coefs[ref]
            if c ≠ 0
                if r == 1
                    c *= tau
                end
                for j in eachindval(idx, idx')
                    v = c * u[j...]
                    if !iszero(v)
                        w[llink=>l, rlink=>r, j...] += v
                    end
                end
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
    make_approx_W2(::State, operator, tau)

build MPO representing approximation WII of a given operator and time step
"""
function make_approx_W2(pre::PreMPO{R}, tau::Number, coefs=[1.]) where R
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
                if r == 1
                    c *= tau
                end
                v[l, r] += c * u
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
            u = v[l, r]
            if !isempty(u)
                for j in eachindval(idx, idx')
                    x = u[j...]
                    if !iszero(x)
                        w[llink=>l, rlink=>r, j...] += x
                    end
                end
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
    