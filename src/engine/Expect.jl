# The expectation values of operators on a state: expect, through the descent and the
# compaction, expect1 and expect2, for an operator of one site on every site or a pair of them
# on every pair of sites, with their fast paths, and inner(a, op, b). The environments on the
# left and on the right of each site, and those on the left carrying the Jordan-Wigner string,
# are kept in the cache `preobs` of the state.

export expect, expect1, expect2

"""
    on_trace(state, t, i)

the tensor `t` of a pure operator on site `i` carried onto the mixed index there, so that
contracting it with the density matrix applies the operator and traces the site out, on a
state conserving nothing strongly, see `weak_form`
"""
function on_trace(state::State{Mixed}, t::ITensor, i::Int)
    s = state.system
    if !isempty(strong_names(s[i]))
        error("bug: site $i of a strongly conserving state measured without going through weak_form")
    end
    j = siteind(Pure, s, i)
    c = last(mixer(j, siteind(Mixed, s, i), s[i]))
    # daggered to meet the mixed index of the state
    return t * dag(c)
end

"""
    identity_at(state, i)

the tensor measuring the identity on site `i`: the identity of the site on a pure state, the
trace over the site on a mixed one
"""
function identity_at(state::State{Pure}, i::Int)
    j = siteind(Pure, state.system, i)
    return delta(dag(j), j')
end

function identity_at(state::State{Mixed}, i::Int)
    j = siteind(Pure, state.system, i)
    return on_trace(state, denseblocks(delta(dag(j), j')), i)
end

"""
    obs_tensor(state, m, i)

the tensor measuring the operator of matrix `m` on site `i`: its tensor on a pure state,
carried onto the trace on a mixed one, see `on_trace`
"""
function obs_tensor(state::State{Pure}, m::AbstractMatrix, i::Int)
    j = siteind(Pure, state.system, i)
    return charged_itensor(m, [j', dag(j)])
end

function obs_tensor(state::State{Mixed}, m::AbstractMatrix, i::Int)
    j = siteind(Pure, state.system, i)
    return on_trace(state, charged_itensor(m, [j', dag(j)]), i)
end

"""
    create_loc!(l, state)

fill `l` with what `get_loc` returns for every site, mixed states only
"""
create_loc!(_, ::State{Pure}) = error("bug: get_loc on pure states")

function create_loc!(l, state::State{Mixed})
    n = length(state)
    resize!(l, n)
    for i in 1:n
        l[i] = identity_at(state, i) * state.state[i]
    end
    return l
end

"""
    get_loc(state, i)

the tensor of the density matrix on site `i` with that site traced out, cached for all sites
on first use, mixed states only
"""
get_loc(state::State, i::Int) = cached!(create_loc!, state.preobs.loc, state)[i]

"""
    create_right!(r, state)

fill `r` with what `get_right` returns for every site; on a pure state, the right orthogonal
sites from `rightlim` on give the identity, uncontracted
"""
function create_right!(r, state::State{Pure})
    st = state.state
    n = length(state)
    rl = ITensorMPS.rightlim(st) - 1
    resize!(r, n)
    r[n] = dag(st[n]')
    for i in n-1:-1:1
        rlink = commonind(st[i], st[i+1])
        v = if i >= rl
            delta(dag(rlink), rlink')
        else
            r[i+1] * identity_at(state, i+1) * st[i+1]
        end
        r[i] = v * dag(st[i]')
    end
    return r
end

function create_right!(r, state::State{Mixed})
    n = length(state)
    resize!(r, n)
    t = r[n] = ITensor(1.)
    for i in n-1:-1:1
        t = r[i] = get_loc(state, i+1) * t
    end
    return r
end

"""
    get_right(state, i)

the environment on the right of site `i`: the sites after `i` traced out on a mixed state,
ket with bra on a pure one, the bra of site `i` included; cached for all sites on first use
"""
get_right(state::State, i::Int) = cached!(create_right!, state.preobs.right, state)[i]

"""
    extend_left!(l, state, i)

extend `l`, which holds what `get_left` returns for its first sites, the first one at least,
up to site `i`; on a pure state, the left orthogonal sites up to `leftlim` give the identity,
uncontracted
"""
function extend_left!(l, state::State{Pure}, i::Int)
    st = state.state
    ll = ITensorMPS.leftlim(st) + 1
    # pushed one by one: an interruption leaves what is there complete
    for k in length(l)+1:i
        llink = commonind(st[k-1], st[k])
        v = if k <= ll
            # the directions the branch below leaves, the ket link as st[k-1] holds it and the
            # bra one daggered and primed: only charges tell the other way round apart
            delta(llink, dag(llink)')
        else
            l[k-1] * identity_at(state, k-1) * dag(st[k-1]')
        end
        push!(l, v * st[k])
    end
    return l
end

function extend_left!(l, state::State{Mixed}, i::Int)
    st = state.state
    for k in length(l)+1:i
        push!(l, l[k-1] * identity_at(state, k-1) * st[k])
    end
    return l
end

"""
    get_left(state, i)

the environment on the left of site `i`: the sites before `i` traced out on a mixed state, ket
with bra on a pure one, times the tensor of the state on `i`, not divided by the trace;
computed up to `i` on demand and cached
"""
function get_left(state::State, i::Int)
    l = state.preobs.left
    # read under the lock too, the vector growing as it is extended
    return lock(state.preobs.lock) do
        if isempty(l)
            push!(l, state.state[1])
        end
        if length(l) < i
            extend_left!(l, state, i)
        end
        l[i]
    end
end

"""
    string_at(state, i)

the tensor measuring the Jordan-Wigner string `F` on site `i`, the identity on a site that is
not fermionic
"""
function string_at(state::State, i::Int)
    f = matrix(F, state.system[i])
    return f == I ? identity_at(state, i) : obs_tensor(state, f, i)
end

"""
    extend_strung!(l, state, i)

extend `l`, which holds what `get_strung_left` returns for its first sites, the first one at
least, up to site `i`
"""
function extend_strung!(l, state::State{Pure}, i::Int)
    st = state.state
    for k in length(l)+1:i
        push!(l, l[k-1] * string_at(state, k-1) * dag(st[k-1]') * st[k])
    end
    return l
end

function extend_strung!(l, state::State{Mixed}, i::Int)
    st = state.state
    for k in length(l)+1:i
        push!(l, l[k-1] * string_at(state, k-1) * st[k])
    end
    return l
end

"""
    get_strung_left(state, i)

the environment on the left of site `i` as `get_left` gives it, with the Jordan-Wigner string
`F` on every site before `i`, which an odd operator on `i` takes; computed up to `i` on demand
and cached
"""
function get_strung_left(state::State, i::Int)
    l = state.preobs.strung
    return lock(state.preobs.lock) do
        if isempty(l)
            push!(l, state.state[1])
        end
        if length(l) < i
            extend_strung!(l, state, i)
        end
        l[i]
    end
end

"""
    Expector(pos, t)
    Expector()

an expectation value being contracted from left to right: `t` holds the sites up to `pos`,
the tensor of the state on `pos` included, its site left open for an operator.
`Expector()` is the empty one, at `pos = 0`.
"""
struct Expector
    pos::Int
    t::ITensor
end

Expector() = Expector(0, ITensor())

"""
    zip_between(state, t, from, to)

the contraction `t` of an expector at `from`, its operator placed there, carried to `to`, the
sites in between traced out and the site of `to` left open
"""
function zip_between(state::State{Pure}, t::ITensor, from::Int, to::Int)
    st = state.state
    t *= dag(st[from]')
    for k in from+1:to-1
        t *= st[k] * identity_at(state, k)
        t *= dag(st[k]')
    end
    return t * st[to]
end

function zip_between(state::State{Mixed}, t::ITensor, from::Int, to::Int)
    for k in from+1:to-1
        t *= get_loc(state, k)
    end
    return t * state.state[to]
end

"""
    zipto(state, a, i)

the expector `a` carried to site `i`, at or after `a.pos`, see `zip_between`; from the empty
expector, `get_left(state, i)`
"""
function zipto(state::State, a::Expector, i::Int)
    if a.pos == 0
        return Expector(i, get_left(state, i))
    elseif a.pos == i
        return a
    end
    return Expector(i, zip_between(state, a.t, a.pos, i))
end

"""
    zipend(state, a)

the expector `a` completed with the environment on the right of `a.pos`, see `get_right`
"""
zipend(state::State, a::Expector) = Expector(a.pos, a.t * get_right(state, a.pos))

"""
    expectfactor(state, a, t, i)

the expector `a` carried to site `i` and multiplied by the tensor `t` measuring an operator
there
"""
function expectfactor(state::State, a::Expector, t::ITensor, i::Int)
    a = zipto(state, a, i)
    return Expector(a.pos, a.t * t)
end

"""
    times_piece(state, t, b, j, comb)

the expector tensor `t`, carried to the `j`-th site of the block `b`, times the piece `comb`
measured there; on a pure state the identity only renames the index of the site
"""
function times_piece(state::State, t::ITensor, b::ComBlock, j::Int, comb::Combination)
    k = b.start + j - 1
    if state isa State{Pure} && collect(keys(comb)) == [0]
        return comb[0] * prime(t, siteind(Pure, state.system, k))
    end
    return t * obs_tensor(state, combination_matrix(b, j, comb, dim(state.system[k])), k)
end

"""
    expect_block(state, b)

the sum of the terms of the block of channels `b` measured on `state`, not divided by its
trace, with one expector per channel
"""
function expect_block(state::State, b::ComBlock)
    open = Union{Nothing, ITensor}[]
    total = 0.
    for (j, ps) in enumerate(b.pieces)
        k = b.start + j - 1
        carried = [ isnothing(t) ? nothing : zip_between(state, t, k - 1, k) for t in open ]
        next = Vector{Union{Nothing, ITensor}}(nothing, get(b.dims, j, 0))
        for ((l, r), comb) in sort!(collect(ps); by = first)
            from = l == 0 ? get_left(state, k) : carried[l]
            if isnothing(from)
                continue
            end
            t = times_piece(state, from, b, j, comb)
            if r == 0
                total += scalar(t * get_right(state, k))
            else
                # out of place: a piece may share its storage with the cache of `get_left`
                next[r] = isnothing(next[r]) ? t : next[r] + t
            end
        end
        open = next
    end
    return total
end

"""
    measured_compacted(system, op)

the operator `op` compacted on `system`, see `compacted_op`, and the system it is on: `system`
itself or, when `op` changes the charges by no definite amount, see `FluxError`, the same sites
conserving nothing, where it is measured on the state conserving nothing, see `plain_form`
"""
measured_compacted(system::System, a::IndexedOp{IsOperator}) =
    or_on_plain(s -> (compacted_op(s, a), s), system)

"""
    PreparedOp

an operator placed on sites made ready to be measured on the states of a system, see
`prepare_op`: its constant, its terms of one site and its block of several, see `compacted_op`,
on the system or, when `plain`, on the same sites conserving nothing, and `scale`, the bound of
its norm, which bounds its expectation values
"""
struct PreparedOp
    plain::Bool
    constant::Number
    singles::Dict
    block::Union{Nothing, ComBlock}
    scale::Float64
end

PreparedOp(c::CompactedOp, plain::Bool) = PreparedOp(plain, c.constant, c.singles, c.block, c.norm)

"""
    prepare_op(system, op)

the operator `op` made ready to be measured on the states of `system`, see `PreparedOp`
"""
function prepare_op(system::System, a::IndexedOp{IsOperator})
    c, on = measured_compacted(system, a)
    return PreparedOp(c, on !== system)
end

"""
    prepared_value(state, p)

the expectation value of the operator `p` prepared by `prepare_op` on the system of `state`,
which `weak_form` gives, divided by the real part of its trace
"""
function prepared_value(state::State, p::PreparedOp)
    if p.plain
        state = plain_form(state)
    end
    v = p.constant * trace(state)
    for k in sort!(collect(keys(p.singles)))
        v += scalar(get_left(state, k) * obs_tensor(state, p.singles[k], k) * get_right(state, k))
    end
    if !isnothing(p.block)
        v += expect_block(state, p.block)
    end
    return v / nonzero_trace(state, "expect")
end

"""
    expect(state, op)

the expectation value of the operator `op`, placed on sites, or the array or the tuple of
those of an array or a tuple of them, on the state of an extension as well, divided by the real
part of the trace of the state: ``\\langle\\psi|op|\\psi\\rangle
/ \\langle\\psi|\\psi\\rangle`` on a pure state, ``\\mathrm{tr}(op\\,\\rho) / \\mathrm{tr}\\,\\rho`` on a
mixed one. Its factors may come in any order, act on several sites, and the Jordan-Wigner
strings of fermionic operators are laid for you. On a system conserving charges, an operator
changing them by no definite amount is measured on the state conserving nothing, see `weaken`.
A channel, a generator and an operator not placed on sites, `X` rather than `X(1)`, see
`expect1`, are refused.

# Examples

    expect(state, X(1) * Y(2) + Y(1) * Z(3))
    expect(state, [X(1) * Y(2), X(3), Swap(1, 3)])
    expect(state, C(3) * dag(C)(1))
"""
function expect(state::State, a::IndexedOp{IsOperator})
    state = weak_form(state)
    return prepared_value(state, prepare_op(state.system, a))
end

expect(state::AbstractState, ops::Union{AbstractArray, Tuple}) = map(o -> expect(state, o), ops)

expect(::State, a::Op{IsChannel}) = error("expect takes an observable, and $a is a channel")

expect(::State, a::Op{IsGenerator}) = error("expect takes an observable, and $a is a generator")

expect(::State, a::GenericOp{IsOperator, N}) where N =
    error("expect takes an operator placed on sites, such as $(a((1:N)...)) rather than $a, " *
          "and expect1 measures an operator of one site on every site")

"""
    defined_on(op, site)

whether the generic operator `op` is defined on `site`: every operator it is made of, and every
operator their definitions are made of, is declared there, every state it projects on, and a
matrix it is given by is of the dimension of the site
"""
function defined_on(a::Operator{N}, site::AbstractSite) where N
    if a.expr isa GenericOp
        return defined_on(a.expr, site)
    elseif a.expr isa Matrix
        return size(a.expr, 1) == dim(site)^N
    elseif !isnothing(a.expr)
        return true
    end
    def = operator_definition(site, Val(Symbol(a.name)))
    return def isa GenericOp ? defined_on(def, site) : !isnothing(def)
end

defined_on(a::Union{SumOp, ProdOp, TensorOp}, site::AbstractSite) =
    all(x -> defined_on(x, site), a.subs)
defined_on(a::Union{ScalarOp, DagOp, ExpOp, IntPowOp, GenPowOp, ModOp}, site::AbstractSite) =
    defined_on(a.arg, site)
defined_on(::IdentityOp, ::AbstractSite) = true

function defined_on(a::Proj, site::AbstractSite)
    if a.state isa Pair
        return defined_on(first(a.state), site)
    elseif a.state isa String
        return !isnothing(state_definition(site, Val(Symbol(a.state)))) ||
               !isnothing(string_state(site, a.state)) || haskey(common_states, a.state)
    end
    return true
end

"""
    defined_sites(system, op)

the sites of `system` where the generic operator `op` is defined, see `defined_on`, refused
when there is none
"""
function defined_sites(system::System, a::GenericOp)
    sites = [ i for i in 1:length(system) if defined_on(a, system[i]) ]
    if isempty(sites)
        error("$a is defined on no site of $system")
    end
    return sites
end

"""
    site_value(state, m, i)

the expectation value, not divided by the trace, of the operator of matrix `m` on site `i`, its
odd part with the Jordan-Wigner string on the sites before
"""
function site_value(state::State, m::AbstractMatrix, i::Int)
    left(p) = p == 0 ? get_left(state, i) : get_strung_left(state, i)
    return sum(scalar(left(p) * obs_tensor(state, part, i) * get_right(state, i))
               for (p, part) in graded_parts(m, matrix(F, state.system[i])))
end

"""
    definite_on(system, m, sites)

whether the operator of matrix `m` on `sites` changes their charges by a definite amount, which
a tensor on charged indices needs; always so on a system conserving nothing
"""
function definite_on(system::System, m::AbstractMatrix, sites)
    idx = [ siteind(Pure, system, k) for k in sites ]
    return has_definite_flux(on_legs(Matrix(m), [ i' for i in idx ], [ dag(i) for i in idx ])...)
end

"""
    definite_everywhere(system, ops, sites)

whether each operator of one site of `ops` changes the charges by a definite amount on each site
of `sites` where it is defined, see `definite_on`
"""
definite_everywhere(system::System, ops, sites) =
    !is_charged(system) || all(definite_on(system, matrix(o, system[i]), [i])
                               for o in ops for i in intersect(sites, defined_sites(system, o)))

"""
    expect1(state, op)

the expectation values of the operator of one site `op` on every site, normalised as by
`expect`: a vector indexed by site, `missing` where `op` is not defined, or, for an array of
operators, an array of such vectors arranged as the operators. The even part of the operator on
a site is measured there on the density matrix of the site, its odd part with the Jordan-Wigner
string on the sites before, both from environments the state keeps. On the state of an
extension, each value is measured by its `expect`.

# Examples

    expect1(state, X)
    expect1(state, [X, Y, Z])
"""
function expect1(state::State, ops::AbstractArray{<:SimpleOp})
    state = weak_form(state)
    if !definite_everywhere(state.system, ops, 1:length(state))
        state = plain_form(state)
    end
    sys = state.system
    n = length(sys)
    scale = nonzero_trace(state, "expect1")
    # the density matrix of each site, shared by the operators measured there
    locals = Vector{Union{Nothing, ITensor}}(nothing, n)
    return map(ops) do o
        on = defined_sites(sys, o)
        map(1:n) do i
            if !(i in on)
                return missing
            end
            v = 0.
            for (p, part) in graded_parts(matrix(o, sys[i]), matrix(F, sys[i]))
                if p == 1
                    v += scalar(get_strung_left(state, i) * obs_tensor(state, part, i) *
                                get_right(state, i))
                    continue
                end
                if isnothing(locals[i])
                    locals[i] = get_left(state, i) * get_right(state, i)
                end
                v += scalar(locals[i] * obs_tensor(state, part, i))
            end
            return v / scale
        end
    end
end

# on the state of an extension, through its expect
function expect1(state::AbstractState, ops::AbstractArray{<:SimpleOp})
    n = length(state)
    return map(ops) do o
        on = defined_sites(state.system, o)
        return [ i in on ? expect(state, o(i)) : missing for i in 1:n ]
    end
end

expect1(state::AbstractState, op::SimpleOp) = only(expect1(state, [op]))

expect1(::AbstractState, a::IndexedOp) =
    error("expect1 measures an operator on every site, X rather than X(1): measure $a with expect")

expect1(::AbstractState, a::Op) =
    error("expect1 measures an operator of one site, and $a is not one")

function expect1(state::AbstractState, ops::AbstractArray)
    for o in ops
        if !(o isa SimpleOp)
            expect1(state, o)
        end
    end
    return expect1(state, convert(Array{SimpleOp}, ops))
end

"""
    pair_parity(system, (o1, o2), sites)

the fermionic parity both operators of the pair have on every site of `sites`, or `nothing`
when they do not have one and the same there
"""
function pair_parity(sys::System, (o1, o2), sites)
    ps = unique(parity_of(matrix(o, sys[i]), matrix(F, sys[i])) for o in (o1, o2) for i in sites)
    return length(ps) == 1 ? only(ps) : nothing
end

"""
    swap_rule(state, (o1, o2), sites, parity)

the function giving the entry `(j, i)`, `j > i`, of the pair `(o1, o2)` from its entry `(i, j)`,
or `nothing`: for `o1 == o2`, the same value, of opposite sign when `parity` is odd, the two
sites being different; on a pure state, for `o2` the adjoint of `o1` on every site of `sites`,
its conjugate, `(o1(i) o2(j))†` being `o1(j) o2(i)`. On a mixed state, which truncation leaves
only nearly hermitian, the second is computed.
"""
function swap_rule(state::State, (o1, o2), sites, parity::Int)
    sys = state.system
    if o1 == o2
        return parity == 1 ? (-) : identity
    elseif state isa State{Pure} && all(i -> nearly(matrix(o2, sys[i]), collect(matrix(o1, sys[i])')), sites)
        return conj
    end
    return nothing
end

"""
    closed_row(state, e, m, i, parity, upto)

the operator of matrix `m` placed on site `i` of the expector `e`, with its string when it is
odd, and the result completed on each later site `j` up to `upto`, that site left open for a
second operator
"""
function closed_row(state::State, e::Expector, m::AbstractMatrix, i::Int, parity::Int, upto::Int)
    sys = state.system
    f(k) = matrix(F, sys[k])
    env = expectfactor(state, e, obs_tensor(state, parity == 1 ? m * f(i) : m, i), i)
    row = Vector{ITensor}(undef, upto)
    for j in i+1:upto
        row[j] = zipend(state, zipto(state, env, j)).t
        if j < upto
            next = parity == 1 ? obs_tensor(state, f(j), j) : identity_at(state, j)
            env = zipto(state, expectfactor(state, env, next, j), j + 1)
        end
    end
    return row
end

"""
    fast_pairs!(r, state, pairs, parities, sites)

fill the entries `(i, j)` of the matrices `r` of the pairs of operators `pairs`, for `i` and `j`
in `sites`, each pair being defined on all of them with one parity there, its entry of
`parities`, see `pair_parity`, and return `r`
"""
function fast_pairs!(r, state::State{Pure}, pairs, parities, sites)
    sys = state.system
    scale = nonzero_trace(state, "expect2")
    rules = [ swap_rule(state, p, sites, q) for (p, q) in zip(pairs, parities) ]
    upto = maximum(sites)
    for i in sites
        e = zipto(state, Expector(), i)
        e = Expector(e.pos, e.t / scale)
        t = zipend(state, e).t
        later = filter(>(i), sites)
        # the rows of an operator placed on i, shared by the pairs it opens
        rows = Dict{Tuple{SimpleOp, Int}, Vector{ITensor}}()
        row(o, q) = get!(() -> closed_row(state, e, matrix(o, sys[i]), i, q, upto), rows, (o, q))
        for (k, ((o1, o2), q, m)) in enumerate(zip(pairs, parities, rules))
            r[k][i, i] = scalar(t * obs_tensor(state, matrix(o1 * o2, sys[i]), i))
            if isempty(later)
                continue
            end
            row1 = row(o1, q)
            row2 = isnothing(m) ? row(o2, q) : nothing
            for j in later
                v = scalar(row1[j] * obs_tensor(state, matrix(o2, sys[j]), j))
                r[k][i, j] = v
                if !isnothing(m)
                    r[k][j, i] = m(v)
                else
                    # swapping two odd operators costs a sign
                    w = scalar(row2[j] * obs_tensor(state, matrix(o1, sys[j]), j))
                    r[k][j, i] = q == 1 ? -w : w
                end
            end
        end
    end
    return r
end

function fast_pairs!(r, state::State{Mixed}, pairs, parities, sites)
    sys = state.system
    scale = nonzero_trace(state, "expect2")
    rules = [ swap_rule(state, p, sites, q) for (p, q) in zip(pairs, parities) ]
    need = unique(parities)
    f(k) = matrix(F, sys[k])
    upto = maximum(sites)
    for i in sites
        e = zipto(state, Expector(), i)
        start = Expector(e.pos, e.t / scale)
        t = zipend(state, start).t
        for (k, (o1, o2)) in enumerate(pairs)
            r[k][i, i] = scalar(t * obs_tensor(state, matrix(o1 * o2, sys[i]), i))
        end
        # the environments carried from i, through the identity for even pairs and through F
        # for odd ones, site i left open
        env = Dict(q => start for q in need)
        for j in i+1:upto
            if j in sites
                ends = Dict(q => zipend(state, zipto(state, env[q], j)) for q in need)
                correlation(a, b, q) =
                    scalar(ends[q].t * obs_tensor(state, q == 1 ? matrix(a, sys[i]) * f(i) : matrix(a, sys[i]), i) *
                           obs_tensor(state, matrix(b, sys[j]), j))
                for (k, ((o1, o2), q, m)) in enumerate(zip(pairs, parities, rules))
                    v = correlation(o1, o2, q)
                    r[k][i, j] = v
                    if !isnothing(m)
                        r[k][j, i] = m(v)
                    else
                        w = correlation(o2, o1, q)
                        r[k][j, i] = q == 1 ? -w : w
                    end
                end
            end
            if j < upto
                for q in need
                    next = q == 1 ? obs_tensor(state, f(j), j) : identity_at(state, j)
                    env[q] = zipto(state, expectfactor(state, env[q], next, j), j + 1)
                end
            end
        end
    end
    return r
end

"""
    general_pairs!(r, state, (o1, o2), on1, on2)

fill the entries `(i, j)` of the matrix `r` of the pair `(o1, o2)`, for `i` in `on1` and `j` in
`on2`, whatever the parities of the two operators, and return `r`. The part of parity `p1` of
`o1` on `i` and that of parity `p2` of `o2` on a later `j` take the string `F` before `i` when
`p1 + p2` is odd, the part on `i` followed by `F^p2`, and `F^p2` between `i` and `j`: one row
of `closed_row` for each. An entry with `j` before `i` is that of the pair written the other way
round, with the sign of swapping two odd parts.
"""
function general_pairs!(r, state::State, (o1, o2), on1, on2)
    sys = state.system
    scale = nonzero_trace(state, "expect2")
    parts(o, k) = graded_parts(matrix(o, sys[k]), matrix(F, sys[k]))
    for i in on1, j in on2
        r[i, j] = i == j ? site_value(state, matrix(o1 * o2, sys[i]), i) / scale : 0
    end
    # the rows of `a` on an earlier site, closed by `b` on each later one; written the other way
    # round, `b(j) a(i)` is `(-1)^(p q) a(i) b(j)` for parts of parities p and q
    function rows!(a, b, ona, onb, swapped)
        upto = maximum(onb; init = 0)
        for i in ona
            later = filter(>(i), onb)
            if isempty(later)
                continue
            end
            closing = Dict(j => parts(b, j) for j in later)
            for (p, x) in parts(a, i), q in unique(q for j in later for (q, _) in closing[j])
                left = isodd(p + q) ? get_strung_left(state, i) : get_left(state, i)
                row = closed_row(state, Expector(i, left / scale), x, i, q, upto)
                for j in later, (q2, y) in closing[j]
                    if q2 ≠ q
                        continue
                    end
                    v = scalar(row[j] * obs_tensor(state, y, j))
                    if swapped
                        r[j, i] += isodd(p * q) ? -v : v
                    else
                        r[i, j] += v
                    end
                end
            end
        end
    end
    rows!(o1, o2, on1, on2, false)
    rows!(o2, o1, on2, on1, true)
    return r
end

"""
    expect2(state, (o1, o2); sites)
    expect2(state, [(o1, o2), ...]; sites)

the matrix whose entry `(i, j)` is the expectation value of `o1(i) * o2(j)`, normalised as by
`expect`, the diagonal being that of `(o1 * o2)(i)`, for `i` and `j` in `sites` (all sites by
default), `missing` elsewhere and where `o1` is not defined on `i` or `o2` on `j`; for a
vector of pairs, a vector of such matrices. The Jordan-Wigner strings are laid for you. A pair
whose operators have one and the same parity on all those sites is measured by sharing the
environments of each site, any other row by row, each part of definite parity with its string.
On the state of an extension, each entry is measured by its `expect`.

# Examples

    expect2(state, (X, X))
    expect2(state, [(X, Y), (C, dag(C))]; sites = 1:10)
"""
function expect2(state::State, pairs::Vector{<:Tuple{SimpleOp, SimpleOp}};
                 sites::AbstractVector{Int} = 1:length(state))
    state = weak_form(state)
    sys = state.system
    n = length(sys)
    for i in sites
        check_index(sys, i, "expect2")
    end
    sites = sort!(unique(collect(sites)))
    if !definite_everywhere(sys, [ o for p in pairs for o in p ], sites)
        state = plain_form(state)
        sys = state.system
    end
    r = [ Matrix{Union{Missing, Number}}(missing, n, n) for _ in pairs ]
    fast, parities = Int[], Int[]
    for (k, (o1, o2)) in enumerate(pairs)
        on1, on2 = intersect(sites, defined_sites(sys, o1)), intersect(sites, defined_sites(sys, o2))
        q = on1 == on2 == sites ? pair_parity(sys, (o1, o2), sites) : nothing
        if !isnothing(q)
            push!(fast, k)
            push!(parities, q)
            continue
        end
        general_pairs!(r[k], state, (o1, o2), on1, on2)
    end
    if !isempty(fast)
        rs = fast_pairs!([ r[k] for k in fast ], state, pairs[fast], parities, sites)
        for (k, m) in zip(fast, rs)
            r[k] = m
        end
    end
    return [ map(identity, m) for m in r ]
end

# on the state of an extension, through its expect
function expect2(state::AbstractState, pairs::Vector{<:Tuple{SimpleOp, SimpleOp}};
                 sites::AbstractVector{Int} = 1:length(state))
    sys = state.system
    n = length(sys)
    for i in sites
        check_index(sys, i, "expect2")
    end
    return map(pairs) do (o1, o2)
        r = Matrix{Union{Missing, Number}}(missing, n, n)
        on1, on2 = intersect(sites, defined_sites(sys, o1)), intersect(sites, defined_sites(sys, o2))
        for i in on1, j in on2
            r[i, j] = expect(state, i == j ? (o1 * o2)(i) : o1(i) * o2(j))
        end
        return map(identity, r)
    end
end

expect2(state::AbstractState, pair::Tuple{SimpleOp, SimpleOp}; kwargs...) =
    only(expect2(state, [pair]; kwargs...))

# ::Any rather than _, which Julia 1.10 refuses before keywords
expect2(::AbstractState, ::Any; kwargs...) =
    error("expect2 correlates two operators on every pair of sites, (X, Y) rather than " *
          "(X(1), Y(2)): measure X(1) * Y(2) with expect")

"""
    inner(a::State, op, b::State)
    dot(a::State, op, b::State)

``\\langle a | op | b \\rangle`` for the operator `op`, placed on sites, between two pure states of
one system, not normalised, through the MPO of `op`, an operator changing the charges by no
definite amount being taken between the states conserving nothing, as by `expect`. Between two
density matrices, for a channel or a generator, ``\\mathrm{tr}(a^\\dagger\\, op(b))``; a pure
operator is ambiguous there, and refused.
"""
function inner(a::State{Pure}, op::IndexedOp{IsOperator}, b::State{Pure})
    same_system(a, b)
    return or_on_plain((x, y) -> inner(x.state', make_mpo(y, op), y.state), a, b)
end

function inner(a::State{Mixed}, op::IndexedOp{<:Union{IsChannel, IsGenerator}}, b::State{Mixed})
    same_system(a, b)
    return inner(a.state', make_mpo(b, op), b.state)
end

inner(::State{Mixed}, op::IndexedOp{IsOperator}, ::State{Mixed}) =
    error("between density matrices the operator $op is ambiguous: write Gate($op) for the map " *
          "ρ ↦ KρK† or Evolver($op) for ρ ↦ Kρ + ρK†")

dot(a::State, op::IndexedOp, b::State) = inner(a, op, b)
