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
        error("bug: site $i of a strongly conserving state measured without going through " *
              "weak_form")
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
sites, those after `ortho_lims`, give the identity, uncontracted
"""
function create_right!(r, state::State{Pure})
    st = state.state
    n = length(state)
    rl = last(ortho_lims(st))
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
    string_at(state, i)

the tensor measuring the Jordan-Wigner string `F` on site `i`, the identity on a site that is
not fermionic
"""
function string_at(state::State, i::Int)
    f = matrix(F, state.system[i])
    return f == I ? identity_at(state, i) : obs_tensor(state, f, i)
end

"""
    extend_left!(l, state, i, at)

extend `l`, which holds the environments on the left of the first sites, the first one at
least, up to site `i`, `at(state, k)` measuring each site before: the identity, see
`identity_at`, or the string, see `string_at`; on a pure state, the left orthogonal sites, those
before `ortho_lims`, give the identity, uncontracted, for the identity
"""
function extend_left!(l, state::State{Pure}, i::Int, at)
    st = state.state
    ll = at === identity_at ? first(ortho_lims(st)) : 0
    # pushed one by one: an interruption leaves what is there complete
    for k in length(l)+1:i
        v = if k <= ll
            llink = commonind(st[k-1], st[k])
            # the directions the branch below leaves, the ket link as st[k-1] holds it and the
            # bra one daggered and primed: only charges tell the other way round apart
            delta(llink, dag(llink)')
        else
            l[k-1] * at(state, k-1) * dag(st[k-1]')
        end
        push!(l, v * st[k])
    end
    return l
end

function extend_left!(l, state::State{Mixed}, i::Int, at)
    st = state.state
    for k in length(l)+1:i
        push!(l, l[k-1] * at(state, k-1) * st[k])
    end
    return l
end

"""
    left_environment(state, l, at, i)

the environment on the left of site `i` kept in the cache `l` of `state`, computed up to `i` on
demand, see `extend_left!`
"""
function left_environment(state::State, l::Vector, at, i::Int)
    # read under the lock too, the vector growing as it is extended
    return lock(state.preobs.lock) do
        if isempty(l)
            push!(l, state.state[1])
        end
        if length(l) < i
            extend_left!(l, state, i, at)
        end
        l[i]
    end
end

"""
    get_left(state, i)

the environment on the left of site `i`: the sites before `i` traced out on a mixed state, ket
with bra on a pure one, times the tensor of the state on `i`, not divided by the trace;
computed up to `i` on demand and cached
"""
get_left(state::State, i::Int) = left_environment(state, state.preobs.left, identity_at, i)

"""
    get_strung_left(state, i)

the environment on the left of site `i` as `get_left` gives it, with the Jordan-Wigner string
`F` on every site before `i`, which an odd operator on `i` takes; computed up to `i` on demand
and cached
"""
get_strung_left(state::State, i::Int) = left_environment(state, state.preobs.strung, string_at, i)

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
    is_hermitian(system, terms, scale)

whether the operator of the `terms` on `system` is Hermitian: their difference with their
adjoints, see `dag(::Term)`, compacts to nothing once its terms of rounding size next to
`scale`, see `term_size`, are dropped
"""
function is_hermitian(system::System, terms::Vector{Term}, scale::Real)
    diff = reduced([ terms ; scaled(-1, map(dag, terms)) ])
    kept = [ t for t in diff if term_size(t) > rounding_tol * scale ]
    if isempty(kept)
        return true
    end
    constant, singles, b = compacted(Pure, system, kept)
    return abs(constant) ≤ rounding_tol * scale && isempty(singles) && isnothing(b)
end

"""
    hermitian_blocks(c)

whether the operator compacted `c`, see `compacted_op`, is Hermitian: zero, or of no flux, an
adjoint changing the charges the other way, and minus its adjoint vanishing, see `vanishes`
"""
hermitian_blocks(c::CompactedOp) =
    isnothing(c.flux) || (c.flux == QN() && vanishes(combined(c.system, [(c, 1), (c', -1)])))

"""
    PreparedOp(op, plain)

an operator placed on sites made ready to be measured on the states of a system, see
`prepare_op`: compacted, see `CompactedOp`, the bound of its norm bounding its expectation
values, on the system or, when `plain`, on the same sites conserving nothing, and whether it is
Hermitian, see `is_hermitian` and `hermitian_blocks`, which makes its values real
"""
struct PreparedOp
    op::CompactedOp{Pure}
    plain::Bool
    hermitian::Bool
end

PreparedOp(c::CompactedOp{Pure}, plain::Bool) =
    PreparedOp(c, plain,
               isnothing(c.terms) ? hermitian_blocks(c) : is_hermitian(c.system, c.terms, c.norm))

"""
    worth_keeping(p)

whether the operator prepared `p` is worth keeping for the next call, see `kept_for!`: of
several terms, or held as the MPO of its factors, where a product of local operators prepares
in about the time finding it would take
"""
worth_keeping(p::PreparedOp) = isnothing(p.op.terms) || length(p.op.terms) > 1

"""
    prepare_op(system, op)

the operator `op` made ready to be measured on the states of `system`, see `PreparedOp`: compacted
on `system`, see `compacted_op`, or, when it changes the charges by no definite amount, see
`FluxError`, its part changing none, see `flux_zero_terms`, the only one with a value on a state
of `system`, Hermitian or not as `op` itself. When that part has too many terms, or `op` holds a
product compacted by blocks, see `by_blocks`, it is compacted on the same sites conserving
nothing, where it is measured on the state conserving nothing, see `plain_form`, a copy that
weighed 452 MiB, against 95 MiB for `Z(20)`, measuring `X(20)` on 40 qubits of link dimension
200. `system` keeps it for the next calls once it is asked for twice, if it is `worth_keeping`,
see `kept_for!`
"""
function prepare_op(system::System, a::IndexedOp{IsOperator})
    return kept_for!(system, system.prepared, a; worth = worth_keeping) do
        try
            return PreparedOp(compacted_op(system, a), false)
        catch e
            if !(e isa FluxError)
                rethrow()
            end
        end
        plain = conserving_nothing(system)
        g = gathered(plain, a)
        if !any(by_blocks, addends(g))
            terms = lower(plain, a, g)
            zero_part = flux_zero_terms(system, terms)
            if !isnothing(zero_part)
                hermitian = is_hermitian(plain, terms, sum(term_size, terms; init = 0.0))
                return PreparedOp(from_terms(Pure, system, a, zero_part), false, hermitian)
            end
        end
        return PreparedOp(compacted_op(plain, a), true)
    end::PreparedOp
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
    c = p.op
    v = c.constant * trace(state)
    for k in sort!(collect(keys(c.singles)))
        v += scalar(get_left(state, k) * obs_tensor(state, c.singles[k], k) * get_right(state, k))
    end
    if !isnothing(c.block)
        v += expect_block(state, c.block)
    end
    return v / nonzero_trace(state, "expect")
end

"""
    site_hermitian(system, o)

whether the operator of one site `o` is Hermitian on every site of `system` where it is
defined, and the largest of its operator norms there, the scale of its values
"""
function site_hermitian(system::System, o::SimpleOp)
    ms = [ matrix(o, system[i]) for i in defined_sites(system, o) ]
    return (all(m -> nearly(m, m'), ms), maximum(opnorm, ms; init = 0.0))
end

"""
    pair_hermitian(system, (a, b))

whether every product `a(i) * b(j)` of the operators of one site `a` and `b` on `system` is
Hermitian, and the product of their largest operator norms, the scale of its values
"""
function pair_hermitian(system::System, (a, b)::Tuple{SimpleOp, SimpleOp})
    on(o) = [ (i, matrix(o, system[i])) for i in defined_sites(system, o) ]
    as, bs = on(a), on(b)
    hermitian(ms) = all(nearly(m, m') for (_, m) in ms)
    even(ms) = all(parity_of(m, matrix(F, system[i])) == 0 for (i, m) in ms)
    commuting = all(nearly(x * y, y * x) for (i, x) in as for (j, y) in bs if i == j)
    # two odd operators anticommute on two sites, which makes their product anti-Hermitian
    h = hermitian(as) && hermitian(bs) && commuting && (even(as) || even(bs))
    norms(ms) = maximum(opnorm(m) for (_, m) in ms; init = 0.0)
    return (h, norms(as) * norms(bs))
end

"""
    imaginary_tol

how large the imaginary part of the expectation value of a Hermitian operator may be, relative
to the scale of the operator, see `term_size`, before `expect` or `measure` warns as it drops it
"""
const imaginary_tol = 1e-8

"""
    real_value(x, scale, name, time)

the real part of the value `x` of the Hermitian operator or measurement `name`, element by
element, with a warning when the imaginary part dropped is beyond rounding next to `scale`, of
the group `:imaginary_part` and carrying `measurement`, `imaginary` and `time`, so that a logger
may gather such warnings
"""
function real_value(x::Number, scale, name, time)
    if abs(imag(x)) > imaginary_tol * scale
        @warn("the imaginary part of a Hermitian measurement is beyond rounding, dropped",
              measurement = name, imaginary = imag(x), time, _group = :imaginary_part)
    end
    return real(x)
end

real_value(::Missing, _, _, _) = missing
real_value(x::AbstractArray, args...) = map(y -> real_value(y, args...), x)

"""
    complex_value(x)

the value `x` of an operator or a measurement that is not Hermitian as a complex number,
element by element, so that its type does not depend on the state measured
"""
complex_value(x::Number) = complex(x)
complex_value(::Missing) = missing
complex_value(x::AbstractArray) = map(complex_value, x)

"""
    ruled_value(x, hermitian, scale, name)

the value `x` of the operator `name` as `expect`, `expect1` and `expect2` give it: real when the
operator is Hermitian, see `real_value`, and complex otherwise, see `complex_value`
"""
ruled_value(x, hermitian::Bool, scale, name) =
    hermitian ? real_value(x, scale, name, nothing) : complex_value(x)

"""
    expect(state, op; complex = false)

the expectation value of the operator `op`, placed on sites, or the array or the tuple of those
of an array or a tuple of them, on the state of an extension as well, through its
`expect_value`, divided by the real part of the trace of the state:
``\\langle\\psi|op|\\psi\\rangle / \\langle\\psi|\\psi\\rangle`` on a pure state,
``\\mathrm{tr}(op\\,\\rho) / \\mathrm{tr}\\,\\rho`` on a mixed one. Its factors may come in any
order, act on several sites, and the Jordan-Wigner strings of fermionic operators are laid
for you. On a system conserving charges, an operator changing them by no definite
amount is measured by its part changing none, the only one with a value on a state of definite
charges, or on the state conserving nothing, see `weaken`, when that part has many more terms
than the operator, as for a long product of `Sx`. A channel, a generator and an operator not
placed on sites, `X` rather than `X(1)`, see `expect1`, are refused.

The value of a Hermitian operator is real, as `measure` gives it, its imaginary part dropped,
with a warning when it is more than rounding, as a density matrix that lost its hermiticity
gives; that of any other operator is complex, and an array of values of both kinds complex.
With `complex = true`, the value is the complex number computed, imaginary part included. The
system of the state keeps an operator of several terms prepared once it is asked for twice.

# Examples

    expect(state, X(1) * Y(2) + Y(1) * Z(3))
    expect(state, [X(1) * Y(2), X(3), Swap(1, 3)])
    expect(state, C(3) * dag(C)(1))
    expect(ρ, X(1); complex = true)
"""
function expect(state::State, a::IndexedOp{IsOperator}; complex::Bool = false)
    state = weak_form(state)
    p = prepare_op(state.system, a)
    x = prepared_value(state, p)
    return complex ? complex_value(x) : ruled_value(x, p.hermitian, p.op.norm, obs_name(a))
end

"""
    expect_value(state, op)

the expectation value of the operator `op` placed on sites, as the state computes it, divided
by the real part of its trace: the number `expect` gives with `complex = true`, and makes real
for a Hermitian operator otherwise. The state of an extension adds a method, which `expect`,
`expect1`, `expect2` and `measure` go through, see `AbstractState`.
"""
function expect_value(state::State, a::IndexedOp{IsOperator})
    state = weak_form(state)
    return prepared_value(state, prepare_op(state.system, a))
end

# on the state of an extension, through its expect_value, Hermitian or not on its physical sites
function expect(state::AbstractState, a::IndexedOp{IsOperator}; complex::Bool = false)
    x = expect_value(state, a)
    if complex
        return complex_value(x)
    end
    p = prepare_op(state.system, a)
    return ruled_value(x, p.hermitian, p.op.norm, obs_name(a))
end

function expect(state::AbstractState, ops::AbstractArray; kwargs...)
    values = map(o -> expect(state, o; kwargs...), ops)
    # real and complex values make an array of complex numbers, as [1.0, 2im], where map would
    # give an array of Number
    T = mapreduce(typeof, promote_type, values; init = Union{})
    return T === Union{} ? values : convert(Array{T}, values)
end

expect(state::AbstractState, ops::Tuple; kwargs...) = map(o -> expect(state, o; kwargs...), ops)

expect(::State, a::Op{IsChannel}; kwargs...) =
    error("expect takes an observable, and $(obs_name(a)) is a channel")

expect(::State, a::Op{IsGenerator}; kwargs...) =
    error("expect takes an observable, and $(obs_name(a)) is a generator")

expect(::State, a::GenericOp{IsOperator, N}; kwargs...) where N =
    error("expect takes an operator placed on sites, such as $(obs_name(a((1:N)...))) rather " *
          "than $(obs_name(a)), and expect1 measures an operator of one site on every site")

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
        error("$(obs_name(a)) is defined on no site of $system")
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
    measured_part(system, m, i)

the matrix `m` of an operator on site `i` of `system`, or, when it changes the charges by no
definite amount, its part changing none, see `flux_zero_part`: the only one with a value on a
state of `system`, and one a tensor on its charged indices can hold
"""
measured_part(system::System, m::AbstractMatrix, i::Int) =
    definite_on(system, m, [i]) ? m : flux_zero_part(m, system[i])

"""
    expect1(state, op; complex = false)

the expectation values of the operator of one site `op` on every site, normalised as by
`expect`: a vector indexed by site, `missing` where `op` is not defined, or, for an array of
operators, an array of such vectors arranged as the operators. The values of an operator
Hermitian on every site where it is defined are real, the others complex, and all of them the
complex numbers computed with `complex = true`, as for `expect`. The even part of the operator
on a site is measured there on the density matrix of the site, its odd part with the
Jordan-Wigner string on the sites before, both from environments the state keeps. On the state
of an extension, each value is measured by its `expect_value`.

# Examples

    expect1(state, X)
    expect1(state, [X, Y, Z])
"""
function expect1(state::AbstractState, ops::AbstractArray{<:SimpleOp}; complex::Bool = false)
    values = expect1_computed(state, ops)
    return map(ops, values) do o, v
        complex ? complex_value(v) : ruled_value(v, site_hermitian(state.system, o)..., obs_name(o))
    end
end

"""
    expect1_computed(state, ops)

the values of `expect1` as they are computed, before they are made real or complex
"""
function expect1_computed(state::State, ops::AbstractArray{<:SimpleOp})
    state = weak_form(state)
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
            for (p, part) in graded_parts(measured_part(sys, matrix(o, sys[i]), i), matrix(F, sys[i]))
                if iszero(part)
                    continue
                end
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

# on the state of an extension, through its expect_value
function expect1_computed(state::AbstractState, ops::AbstractArray{<:SimpleOp})
    n = length(state)
    return map(ops) do o
        on = defined_sites(state.system, o)
        return [ i in on ? expect_value(state, o(i)) : missing for i in 1:n ]
    end
end

expect1(state::AbstractState, op::SimpleOp; kwargs...) = only(expect1(state, [op]; kwargs...))

expect1(::AbstractState, a::IndexedOp; kwargs...) =
    error("expect1 measures an operator on every site, X rather than X(1): measure " *
          "$(obs_name(a)) with expect")

expect1(::AbstractState, a::Op; kwargs...) =
    error("expect1 measures an operator of one site, and $(obs_name(a)) is not one")

function expect1(state::AbstractState, ops::AbstractArray; kwargs...)
    for o in ops
        if !(o isa SimpleOp)
            expect1(state, o)
        end
    end
    return expect1(state, convert(Array{SimpleOp}, ops); kwargs...)
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
    elseif state isa State{Pure} &&
           all(i -> nearly(matrix(o2, sys[i]), collect(matrix(o1, sys[i])')), sites)
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
                function correlation(a, b, q)
                    x = q == 1 ? matrix(a, sys[i]) * f(i) : matrix(a, sys[i])
                    return scalar(ends[q].t * obs_tensor(state, x, i) *
                                  obs_tensor(state, matrix(b, sys[j]), j))
                end
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
`on2`, whatever the parities and the charges of the two operators, and return `r`. The part of
parity `p1` of `o1` on `i` and that of parity `p2` of `o2` on a later `j` take the string `F`
before `i` when `p1 + p2` is odd, the part on `i` followed by `F^p2`, and `F^p2` between `i` and
`j`: one row of `closed_row` for each. Each part is split by charge as well, see `charge_parts`,
and only parts whose charges cancel are paired, the others having no value on a state of
definite charges; the diagonal is the part of `o1 * o2` changing no charge, see
`measured_part`. An entry with `j` before `i` is that of the pair written the other way round,
with the sign of swapping two odd parts.
"""
function general_pairs!(r, state::State, (o1, o2), on1, on2)
    sys = state.system
    scale = nonzero_trace(state, "expect2")
    # the parts of each parity and of each charge, see `charge_parts`, a single charge on a site
    # conserving nothing
    parts(o, k) = [ (p, c, x) for (p, y) in graded_parts(matrix(o, sys[k]), matrix(F, sys[k]))
                    for (c, x) in charge_parts(y, sys[k]) if !iszero(x) ]
    # a float zero, which the entries no pair of parts reaches keep: an integer one would come out
    # as a Complex{Int}
    for i in on1, j in on2
        m = i == j ? measured_part(sys, matrix(o1 * o2, sys[i]), i) : nothing
        r[i, j] = isnothing(m) || iszero(m) ? 0. : site_value(state, m, i) / scale
    end
    # the rows of `a` on an earlier site, closed by `b` on each later one; written the other way
    # round, `b(j) a(i)` is `(-1)^(p q) a(i) b(j)` for parts of parities p and q. Parts whose
    # charges do not cancel have no value on a state of definite charges, and are not paired
    function rows!(a, b, ona, onb, swapped)
        upto = maximum(onb; init = 0)
        for i in ona
            later = filter(>(i), onb)
            if isempty(later)
                continue
            end
            closing = Dict(j => parts(b, j) for j in later)
            for (p, ca, x) in parts(a, i),
                q in unique(q for j in later for (q, cb, _) in closing[j] if ca + cb == QN())
                left = isodd(p + q) ? get_strung_left(state, i) : get_left(state, i)
                row = closed_row(state, Expector(i, left / scale), x, i, q, upto)
                for j in later, (q2, cb, y) in closing[j]
                    if q2 ≠ q || ca + cb ≠ QN()
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
    expect2(state, (o1, o2); sites, complex = false)
    expect2(state, [(o1, o2), ...]; sites, complex = false)

the matrix whose entry `(i, j)` is the expectation value of `o1(i) * o2(j)`, normalised as by
`expect`, the diagonal being that of `(o1 * o2)(i)`, for `i` and `j` in `sites` (all sites by
default), `missing` elsewhere and where `o1` is not defined on `i` or `o2` on `j`; for a
vector of pairs, a vector of such matrices. The matrix of a pair whose products are all
Hermitian, `o1` and `o2` Hermitian, commuting on a site and not both fermionic, is real, the
others complex, and all of them the complex numbers computed with `complex = true`, as for
`expect`. The Jordan-Wigner strings are laid for you. A pair whose operators have one and the
same parity on all those sites, and change the charges by definite amounts there, is measured
by sharing the environments of each site, any other row by row, each part of definite parity
and charge with its string, the parts whose charges do not cancel left out. On the state of an
extension, each entry is measured by its `expect_value`.

# Examples

    expect2(state, (X, X))
    expect2(state, [(X, Y), (C, dag(C))]; sites = 1:10)
"""
function expect2(state::AbstractState, pairs::Vector{<:Tuple{SimpleOp, SimpleOp}};
                 sites::AbstractVector{Int} = 1:length(state), complex::Bool = false)
    values = expect2_computed(state, pairs; sites)
    return map(pairs, values) do (a, b), v
        complex ? complex_value(v) :
                  ruled_value(v, pair_hermitian(state.system, (a, b))...,
                              "($(obs_name(a)), $(obs_name(b)))")
    end
end

"""
    expect2_computed(state, pairs; sites)

the values of `expect2` as they are computed, before they are made real or complex
"""
function expect2_computed(state::State, pairs::Vector{<:Tuple{SimpleOp, SimpleOp}};
                          sites::AbstractVector{Int} = 1:length(state))
    state = weak_form(state)
    sys = state.system
    n = length(sys)
    for i in sites
        check_index(sys, i, "expect2")
    end
    sites = sort!(unique(collect(sites)))
    r = [ Matrix{Union{Missing, Number}}(missing, n, n) for _ in pairs ]
    fast, parities = Int[], Int[]
    for (k, (o1, o2)) in enumerate(pairs)
        on1 = intersect(sites, defined_sites(sys, o1))
        on2 = intersect(sites, defined_sites(sys, o2))
        # a pair changing the charges by no definite amount goes row by row, its parts of
        # charges cancelling alone, see `general_pairs!`
        q = on1 == on2 == sites && definite_everywhere(sys, [o1, o2], sites) ?
            pair_parity(sys, (o1, o2), sites) : nothing
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

# on the state of an extension, through its expect_value
function expect2_computed(state::AbstractState, pairs::Vector{<:Tuple{SimpleOp, SimpleOp}};
                          sites::AbstractVector{Int} = 1:length(state))
    sys = state.system
    n = length(sys)
    for i in sites
        check_index(sys, i, "expect2")
    end
    return map(pairs) do (o1, o2)
        r = Matrix{Union{Missing, Number}}(missing, n, n)
        on1 = intersect(sites, defined_sites(sys, o1))
        on2 = intersect(sites, defined_sites(sys, o2))
        for i in on1, j in on2
            r[i, j] = expect_value(state, i == j ? (o1 * o2)(i) : o1(i) * o2(j))
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
the same sites, taken on what both conserve, not normalised, through the MPO of `op`, an
operator changing the charges by no definite amount being taken between the states conserving
nothing, as by `expect`. Between two density matrices, for a channel or a generator,
``\\mathrm{tr}(a^\\dagger\\, op(b))``; a pure operator is ambiguous there, and refused.
"""
function inner(a::State{Pure}, op::IndexedOp{IsOperator}, b::State{Pure})
    a, b = on_common_system(a, b)
    return or_on_plain((x, y) -> inner(x.state', make_mpo(y, op), y.state), a, b)
end

function inner(a::State{Mixed}, op::IndexedOp{<:Union{IsChannel, IsGenerator}}, b::State{Mixed})
    a, b = on_common_system(a, b)
    return inner(a.state', make_mpo(b, op), b.state)
end

function inner(::State{Mixed}, op::IndexedOp{IsOperator}, ::State{Mixed})
    name = obs_name(op)
    error("between density matrices the operator $name is ambiguous: write Gate($name) for the " *
          "map ρ ↦ KρK† or Evolver($name) for ρ ↦ Kρ + ρK†")
end

dot(a::State, op::IndexedOp, b::State) = inner(a, op, b)
