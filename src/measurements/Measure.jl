# The measurement system: Measure, a list of measurements prepared once for each system it is
# measured on, the measurements of the library, the measurements written as functions of a
# context and the Check they may return, and measure, which computes a list on a state.

export Measure, measure, Check, Purity, Trace, TraceDeviation, Norm, Hermiticity,
       HermiticityDeviation, Renyi2, MaxLinkdim, MemoryUsage, SubRenyi2, MutualInfoRenyi2,
       EntanglementEntropy, ReducedDensityMatrix, VonNeumannEntropy, LogNegativity, Fidelity,
       Overlap, Variance

"""
    Context

what a measurement written as a function receives: `c.state`, `c.time`, and the values given
to `measure` by their names, as `c.energy`; reading one that was not given is refused
"""
struct Context
    state::AbstractState
    time::Union{Nothing, Number}
    values::Dict{Symbol, Any}
end

function getproperty(c::Context, name::Symbol)
    if name === :state
        return getfield(c, :state)
    elseif name === :time
        t = getfield(c, :time)
        if isnothing(t)
            error("this measurement reads the time, which measure was not given: " *
                  "measure(state, measurements; time = t)")
        end
        return t
    end
    values = getfield(c, :values)
    if !haskey(values, name)
        error("this measurement reads $name, which measure was not given: " *
              "measure(state, measurements; $name = …)")
    end
    return values[name]
end

propertynames(c::Context) = (:state, :time, keys(getfield(c, :values))...)

"""
    Check(value, reference; tol = nothing, warn = false)

what a measurement written as a function returns to compare `value` with `reference`, numbers
or arrays of one shape: `measure` gives `(value, reference, difference)`, the difference being
`norm(value - reference)` over the entries both have, `missing` when one misses an entry the
other has, and refuses a difference beyond `tol` or missing, or only warns with `warn = true`.
Without `tol`, the three are recorded and nothing is refused.

# Examples

    "check_cos" => c -> Check(expect(c.state, X(1)), cos(2c.time); tol = 1e-6)
"""
struct Check
    value
    reference
    tol::Union{Nothing, Real}
    warn::Bool
end

Check(value, reference; tol::Union{Nothing, Real} = nothing, warn::Bool = false) =
    Check(value, reference, tol, warn)

"""
    checked(name, check, time)

the value of the `Check` returned by the measurement `name` at `time`, `nothing` when measure
was given no time, refused or warned about when the check fails
"""
function checked(name::String, ch::Check, time)
    a, b = ch.value, ch.reference
    shape(x) = x isa AbstractArray ? size(x) : ()
    if shape(a) ≠ shape(b)
        error("Check $name compares values of different shapes, $(shape(a)) and $(shape(b))")
    end
    same = ismissing.(a) == ismissing.(b)
    d = same ? sqrt(sum(abs2, skipmissing(a isa AbstractArray ? a .- b : [a - b]); init = 0.0)) :
               missing
    if !isnothing(ch.tol) && (!same || d > ch.tol)
        when = isnothing(time) ? "" : " at time $time"
        message = same ? "Check $name failed$when: $a and $b differ by $d, beyond tol = $(ch.tol)" :
                         "Check $name failed$when: $a and $b miss different entries"
        if ch.warn
            @warn message
        else
            error(message)
        end
    end
    return (value = a, reference = b, difference = d)
end

"""
    StateMeasurement(name, prepare)

a measurement of the library, written under `name`: `prepare(state)` gives the function of a
state measuring it on the states of the representation and the system of `state`, see `Measure`
"""
struct StateMeasurement
    name::String
    prepare::Function
end

show(io::IO, m::StateMeasurement) = print(io, m.name)

"""
    ContextMeasurement(name, f, reads)

a measurement of the library, written under `name`, that the function `f` of the context gives,
reading the values `reads`, as a reference state named in the simulation
"""
struct ContextMeasurement
    name::String
    f::Function
    reads::Vector{Symbol}
end

show(io::IO, m::ContextMeasurement) = print(io, m.name)

"""
    OnNamed(state, entry)

the measurement `entry` made on the state named `state`, as a list `:φ => [Z, Purity]` gives
it, its value `missing` while there is no such state
"""
struct OnNamed
    state::Symbol
    entry
end

"""
    measurement_reads(m)

the values the measurement `m` reads: a symbol itself, a measurement of the library reading a
named state that name, a measurement of a named state that name and what it reads
"""
measurement_reads(s::Symbol) = [s]
measurement_reads(m::ContextMeasurement) = m.reads
measurement_reads(m::OnNamed) = [m.state; measurement_reads(m.entry)]
measurement_reads(_) = Symbol[]

"""
    on_state(name, f)

the measurement of the library `name` that `f`, a function of the state, gives, with nothing
to prepare
"""
on_state(name::String, f) = StateMeasurement(name, _ -> f)

"""
    compact_positions(positions)

a set of site positions written in short for a name, runs of three or more consecutive sites
becoming ranges, two different sets giving two different names

# Examples

    compact_positions(1:20)        # "1:20"
    compact_positions([1,2,3,7,8]) # "1:3,7,8"
"""
function compact_positions(p)
    v = sort(unique(collect(p)))
    if isempty(v)
        return "[]"
    end
    parts = String[]
    i = 1
    while i ≤ length(v)
        j = i
        while j < length(v) && v[j + 1] == v[j] + 1
            j += 1
        end
        push!(parts, j > i + 1 ? "$(v[i]):$(v[j])" : join(v[i:j], ","))
        i = j + 1
    end
    return join(parts, ",")
end

"""
    on_positions(name, f, positions...)

the measurement of the library `name` that `f(state, positions...)` gives, named after the
positions, which are checked when it is prepared
"""
function on_positions(name::String, f, positions...)
    full = "$name($(join(map(compact_positions, positions), ";")))"
    return StateMeasurement(full, function (state)
        for p in positions
            check_positions(state, p, name)
        end
        return st -> f(st, positions...)
    end)
end

"""
    Purity

the measurement of the purity of the state, see `purity`
"""
const Purity = on_state("Purity", purity)

"""
    Trace

the measurement of the trace of the state, see `trace`, with its imaginary part: a complex
number, whatever the state, so that its columns do not change during a run
"""
const Trace = on_state("Trace", st -> complex(trace(st)))

"""
    TraceDeviation

the measurement of `1 - trace(state)`, see `Trace`, which numerical errors tend to move, a
complex number as `Trace`
"""
const TraceDeviation = on_state("TraceDeviation", st -> complex(1 - trace(st)))

"""
    Norm

the measurement of the norm of the state, see `norm`
"""
const Norm = on_state("Norm", norm)

"""
    Hermiticity

the measurement of how Hermitian the density matrix is, see `hermiticity`
"""
const Hermiticity = on_state("Hermiticity", hermiticity)

"""
    HermiticityDeviation

the measurement of `1 - hermiticity(state)`, see `Hermiticity`: the squared norm of the
anti-Hermitian part of the density matrix relative to that of the whole
"""
const HermiticityDeviation = on_state("HermiticityDeviation", st -> 1 - hermiticity(st))

"""
    Renyi2

the measurement of the Rényi-2 entropy of the state, see `renyi2`
"""
const Renyi2 = on_state("Renyi2", renyi2)

"""
    MaxLinkdim

the measurement of the largest bond dimension of the state, see `maxlinkdim`
"""
const MaxLinkdim = on_state("MaxLinkdim", maxlinkdim)

"""
    MemoryUsage

the measurement of the memory the state occupies, in bytes, the caches the measurements fill
included
"""
const MemoryUsage = on_state("MemoryUsage", Base.summarysize)

"""
    SubRenyi2(positions)

the measurement of the Rényi-2 entropy of the sites at `positions`, see `renyi2`

# Examples

    SubRenyi2(1:3)
"""
SubRenyi2(positions::AbstractVector) = on_positions("SubRenyi2", renyi2, positions)

"""
    MutualInfoRenyi2(positions)

the measurement of the Rényi-2 mutual information between the sites at `positions` and the
others, see `mutual_info_renyi2`
"""
MutualInfoRenyi2(positions::AbstractVector) =
    on_positions("MutualInfoRenyi2", mutual_info_renyi2, positions)

"""
    ReducedDensityMatrix(positions)

the measurement of the density matrix of the sites at `positions`, see
`reduced_density_matrix`, a complex matrix whatever the state, as `Trace`
"""
ReducedDensityMatrix(positions::AbstractVector) =
    on_positions("ReducedDensityMatrix", (st, p) -> complex(reduced_density_matrix(st, p)),
                 positions)

"""
    VonNeumannEntropy(positions)

the measurement of the von Neumann entropy of the sites at `positions`, see
`von_neumann_entropy`
"""
VonNeumannEntropy(positions::AbstractVector) =
    on_positions("VonNeumannEntropy", von_neumann_entropy, positions)

"""
    LogNegativity(a, b)

the measurement of the logarithmic negativity between the sites at positions `a` and those at
`b`, see `log_negativity`
"""
LogNegativity(a::AbstractVector, b::AbstractVector) =
    on_positions("LogNegativity", log_negativity, a, b)

"""
    EntanglementEntropy(cut)
    EntanglementEntropy(cut, n)

the measurement of the entanglement entropy across the cut between sites `cut` and `cut + 1`,
see `entanglement_entropy`, followed, given `n`, by the first `n` values of the spectrum, see
`entanglement_spectrum`, padded with zeros beyond the bond dimension of the cut
"""
EntanglementEntropy(cut::Int) =
    StateMeasurement("EntanglementEntropy($cut)", function (state)
        check_positions(state, [cut], "EntanglementEntropy")
        return st -> entanglement_entropy(st, cut)
    end)

EntanglementEntropy(cut::Int, n::Int) =
    StateMeasurement("EntanglementEntropy($cut,$n)", function (state)
        check_positions(state, [cut], "EntanglementEntropy")
        return function (st)
            # the entropy of the spectrum, which one decomposition gives
            sp = entanglement_spectrum(st, cut)
            return [ shannon_entropy(sp) ; sp[1:min(length(sp), n)] ;
                     zeros(max(0, n - length(sp))) ]
        end
    end)

"""
    reference_on(state, ref, what)

the reference state `ref` of the measurement `what` weakened to what the measured `state`
conserves, a reference conserving less being refused
"""
function reference_on(state::AbstractState, ref::State, what::String)
    target = symmetries(state.system)
    source = symmetries(ref.system)
    for (name, strong) in target.names
        k = findfirst(q -> first(q) == name, source.names)
        if isnothing(k) || (strong && !last(source.names[k]))
            error("the reference of $what conserves less than the measured state, which " *
                  "conserves $target")
        end
    end
    return weaken(ref, target)
end

"""
    on_named(what, name, f)

the measurement `what` of the state against the reference named `name`, `f(state, ref)`, the
reference weakened to what the state conserves, `missing` while there is none
"""
on_named(what::String, name::Symbol, f) = ContextMeasurement(what, function (c)
    ref = getproperty(c, name)
    return ismissing(ref) ? missing : f(c.state, reference_on(c.state, ref, what))
end, [name])

"""
    Fidelity(ref)
    Fidelity(name)

the measurement of the fidelity with the reference state `ref`, see `fidelity`, `ref` being
weakened to what the measured state conserves; given a `name`, the reference is the value of
that name, a state named in the simulation, `missing` while there is none
"""
Fidelity(ref::State) = StateMeasurement("Fidelity", function (state)
    r = reference_on(state, ref, "Fidelity")
    return st -> fidelity(st, r)
end)

Fidelity(name::Symbol) = on_named("Fidelity", name, (st, r) -> fidelity(st, r))

"""
    applied_element(ref, op, state)

``\\langle ref | op | state \\rangle`` for the operator `op` placed on sites, applied to `state`,
or between the states conserving nothing when it changes the charges by no definite amount, as
`expect` takes it
"""
applied_element(ref::State{Pure}, a::IndexedOp{IsOperator}, state::State{Pure}) =
    or_on_plain((r, s) -> complex(inner(r, apply(a, s))), ref, state)

"""
    matrix_element(ref, op, state)

``\\langle ref | op | state \\rangle`` for two pure states, with an operator placed on sites, or
the vector of its values on each site for an operator of one site, `missing` where it is not
defined
"""
function matrix_element(ref::State, op, state::State)
    if !(ref isa State{Pure} && state isa State{Pure})
        error("Overlap with an operator is the matrix element of two pure states")
    end
    if op isa SimpleOp
        on = defined_sites(state.system, op)
        return [ i in on ? applied_element(ref, op(i), state) : missing for i in 1:length(state) ]
    end
    return applied_element(ref, op, state)
end

"""
    Overlap(ref[, op])
    Overlap(name[, op])

the measurement of the inner product with the reference state `ref`, see `inner`, `ref` being
weakened to what the measured state conserves, a complex number whatever the states, as
`Trace`; given a `name`, the reference is the value of that name, as for `Fidelity`. With an
operator, the matrix element ``\\langle ref | op | state \\rangle`` of two pure states, of an
operator placed on sites, or the vector of those of an operator of one site on every site: from
a reference `A|ψ⟩` evolved with the state `|ψ⟩`, the correlations ``\\langle A^\\dagger(0) op(t)
\\rangle`` at two times.

# Examples

    Overlap(:excited)
    Overlap(:kicked, Z)            # ⟨kicked|Z_i|state⟩ for every site i
"""
Overlap(ref::State) = StateMeasurement("Overlap", function (state)
    r = reference_on(state, ref, "Overlap")
    return st -> complex(inner(r, st))
end)

Overlap(name::Symbol) = on_named("Overlap", name, (st, r) -> complex(inner(r, st)))

Overlap(ref::State, op::Union{IndexedOp{IsOperator}, SimpleOp}) =
    StateMeasurement("Overlap", function (state)
        r = reference_on(state, ref, "Overlap")
        return st -> matrix_element(r, op, st)
    end)

Overlap(name::Symbol, op::Union{IndexedOp{IsOperator}, SimpleOp}) =
    on_named("Overlap", name, (st, r) -> matrix_element(r, op, st))


"""
    variance_of(state, op)

the function measuring the variance of `op` on the states of the system of `state`, its MPO
built once, or by the `variance` of the state of an extension; a mixed state is refused here,
before any measurement
"""
function variance_of(state::State{Pure}, op)
    on, mpo = or_on_plain(s -> (s, make_mpo(s, op)), state)
    # the states to come measured on the same form as `state`
    return on === state ? (st -> variance(st, mpo)) : (st -> variance(plain_form(st), mpo))
end

variance_of(state::State{Mixed}, op) = variance(state, op)

variance_of(::AbstractState, op) = st -> variance(st, op)

"""
    Variance(op)

the measurement of the variance of the operator `op`, see `variance`, its MPO built once for
each system
"""
Variance(op) = StateMeasurement("Variance", state -> variance_of(state, op))

"""
    checked_measurement(x)

`x` checked to be a measurement: an operator placed on sites, an operator of one site, a pair
of them, a measurement of the library, or a symbol; anything else is refused
"""
checked_measurement(a::IndexedOp{IsOperator}) = a
checked_measurement(a::SimpleOp) = a
checked_measurement(p::Tuple{SimpleOp, SimpleOp}) = p
checked_measurement(m::StateMeasurement) = m
checked_measurement(m::ContextMeasurement) = m
checked_measurement(m::OnNamed) = m
checked_measurement(s::Symbol) = s

checked_measurement(::Function) =
    error("a measurement written as a function needs a name: \"name\" => c -> …")

checked_measurement(x::Number) = error("a constant is not a measurement: write \"name\" => c -> $x")

checked_measurement(x::String) =
    error("\"$x\" is not a measurement: a name goes before one, \"$x\" => measurement")

checked_measurement(::Check) =
    error("a Check is returned by a measurement written as a function: \"name\" => c -> Check(…)")

checked_measurement(::Tuple) =
    error("a pair measured on every pair of sites is (X, Y), two operators of one site")

checked_measurement(::AbstractArray) =
    error("a name goes before one measurement, not before a list")

checked_measurement(a::GenericOp{IsOperator}) =
    error("$(obs_name(a)) acts on several sites, and is measured placed on them")

checked_measurement(a::Op) = error("measure takes observables, and $(obs_name(a)) is not one")

checked_measurement(x) = error("$(repr(x)) is not a measurement")

"""
    named_measurement(x)

`x` checked to be what may follow a name: a function of the context, or a measurement, see
`checked_measurement`; told apart by its value, a vector of named measurements having pairs of
type `Pair{String, Any}`
"""
named_measurement(f::Function) = f
named_measurement(x) = checked_measurement(x)

"""
    entry_name(m)

the name of the measurement `m` when none is given: the compact form of an operator, `(A, B)`
for a pair, the name of a measurement of the library, a symbol as itself
"""
entry_name(a::Op) = obs_name(a)
entry_name((a, b)::Tuple{SimpleOp, SimpleOp}) = "($(obs_name(a)), $(obs_name(b)))"
entry_name(m::Union{StateMeasurement, ContextMeasurement}) = m.name
entry_name(s::Symbol) = string(s)

"""
    entries(list)

the measurements of `list` as pairs `name => measurement`, the vectors and the `Measure` it
holds flattened, a matrix row by row, each named as given or by `entry_name`, and checked, see
`checked_measurement`
"""
entries(x::AbstractVector) = reduce(vcat, [ entries(y) for y in x ]; init = Pair{String, Any}[])
entries(x::AbstractMatrix) = entries(vec(permutedims(x)))
entries((name, x)::Pair{String}) = Pair{String, Any}[ name => named_measurement(x) ]
entries((state, x)::Pair{Symbol}) =
    Pair{String, Any}[ "$state.$n" => OnNamed(state, m) for (n, m) in entries(x) ]

function entries(x)
    m = checked_measurement(x)
    return Pair{String, Any}[ entry_name(m) => m ]
end

"""
    Measure(measurements...)
    Measure([measurements...])

a list of measurements, measured together by `measure`, which prepares it once for each system
it is measured on and keeps that preparation. Each measurement is one of

- an operator placed on sites, `X(1) * X(2)`, measured by `expect`
- an operator of one site, `X`, measured on every site where it is defined by `expect1`
- a pair of them, `(X, Y)`, measured on every pair of sites by `expect2`
- a measurement of the library, as `Purity` or `EntanglementEntropy(3)`
- a symbol, `:energy`, the value of that name given to `measure`
- `"name" => m`, `m` one of these or a function of the context `c`: `c.state`, `c.time` and the
  values given to `measure`, `c.energy`
- `:φ => list`, the measurements of the list made on the state named `φ`, a value given to
  `measure`, their names prefixed by `φ.`

A vector, or a `Measure`, inside the list stands for its measurements. Two measurements of a
list cannot share a name.

# Examples

    Measure(X, (X, Y), Z(1) * Z(2), Purity, "E" => H, "cos" => c -> cos(2c.time), :energy)
"""
struct Measure
    entries::Vector{Pair{String, Any}}
    prepared::Dict{Tuple{DataType, System}, Vector{Any}}
    lock::ReentrantLock
    function Measure(list::AbstractVector)
        es = entries(list)
        names = first.(es)
        twice = unique(n for n in names if count(==(n), names) > 1)
        if !isempty(twice)
            error("several measurements are named $(join(repr.(twice), ", ")): name them apart, " *
                  "\"name\" => measurement")
        end
        return new(es, Dict{Tuple{DataType, System}, Vector{Any}}(), ReentrantLock())
    end
end

Measure(args...) = Measure(collect(Any, args))

entries(m::Measure) = m.entries

show(io::IO, m::Measure) = print(io, "Measure(", join(first.(m.entries), ", "), ")")

"""
    OpValue(prepared, op)

an operator placed on sites as `measure` prepares it, see `prepare_op`, measured on the state
of an extension by its `expect_value`; its value is made real when it is Hermitian
"""
struct OpValue
    prepared::PreparedOp
    op::IndexedOp{IsOperator}
end

"""
    SiteValues(op, hermitian, scale)

an operator of one site as `measure` prepares it, measured by `expect1`: whether it is
Hermitian on every site where it is defined, and the largest of its operator norms there
"""
struct SiteValues
    op::SimpleOp
    hermitian::Bool
    scale::Float64
end

"""
    PairValues(pair, hermitian, scale)

a pair of operators of one site as `measure` prepares it, measured by `expect2`: whether each
product is Hermitian, and the product of the largest operator norms of the two
"""
struct PairValues
    pair::Tuple{SimpleOp, SimpleOp}
    hermitian::Bool
    scale::Float64
end

"""
    OnState(f)

a measurement of the library as `measure` prepares it, `f` a function of the state
"""
struct OnState
    f::Function
end

"""
    OnNamedItem(state, entry, measure)

a measurement of a named state as `measure` prepares it: the name of the state, the measurement,
and the `Measure` of it, made at its first value under the name of the measurement, which keeps
its preparation for the system of that state
"""
mutable struct OnNamedItem
    state::Symbol
    entry
    measure::Union{Nothing, Measure}
end

"""
    prepared_entry(state, m)

the measurement `m` made ready for the states of the representation and of the system of
`state`, see `OpValue`, `SiteValues`, `PairValues` and `OnState`, a function or a symbol
staying as it is
"""
function prepared_entry(state::AbstractState, a::IndexedOp{IsOperator})
    system = state isa State ? weak_form(state).system : state.system
    return OpValue(prepare_op(system, a), a)
end

prepared_entry(state::AbstractState, o::SimpleOp) =
    SiteValues(o, site_hermitian(state.system, o)...)

prepared_entry(state::AbstractState, pair::Tuple{SimpleOp, SimpleOp}) =
    PairValues(pair, pair_hermitian(state.system, pair)...)

prepared_entry(state::AbstractState, m::StateMeasurement) = OnState(m.prepare(state))
prepared_entry(::AbstractState, f::Function) = f
prepared_entry(::AbstractState, m::ContextMeasurement) = m.f
prepared_entry(::AbstractState, s::Symbol) = s
prepared_entry(::AbstractState, m::OnNamed) = OnNamedItem(m.state, m.entry, nothing)

"""
    prepared(m, state)

the measurements of the `Measure` `m` made ready for the states of the representation and of
the system of `state`, see `prepared_entry`, prepared once and kept in `m`
"""
prepared(m::Measure, state::AbstractState) =
    lock(m.lock) do
        get!(m.prepared, (typeof(state), state.system)) do
            Any[ prepared_entry(state, last(e)) for e in m.entries ]
        end
    end

"""
    entry_value(item, name, c, computed)

the value of the prepared measurement `item`, named `name`, in the context `c`, `computed`
holding the values of `expect1` and of `expect2` for the operators and the pairs of the list
"""
function entry_value(it::OpValue, name, c::Context, _)
    state = c.state
    x = state isa State ? prepared_value(weak_form(state), it.prepared) : expect_value(state, it.op)
    return it.prepared.hermitian ? real_value(x, it.prepared.op.norm, name, getfield(c, :time)) :
                                   complex_value(x)
end

function entry_value(it::Union{SiteValues, PairValues}, name, c::Context, computed)
    x = it isa SiteValues ? computed.sites[it.op] : computed.pairs[it.pair]
    return it.hermitian ? real_value(x, it.scale, name, getfield(c, :time)) : complex_value(x)
end

entry_value(it::OnState, _, c::Context, _) = it.f(c.state)

function entry_value(f::Function, name, c::Context, _)
    v = f(c)
    return v isa Check ? checked(name, v, getfield(c, :time)) : v
end

entry_value(s::Symbol, _, c::Context, _) = getproperty(c, s)

function entry_value(it::OnNamedItem, name, c::Context, _)
    state = getproperty(c, it.state)
    if ismissing(state)
        return missing
    end
    if isnothing(it.measure)
        it.measure = Measure(Pair{String, Any}[ name => it.entry ])
    end
    values = getfield(c, :values)
    return last(only(measure(state, it.measure; time = getfield(c, :time), values...)))
end

"""
    measure(state, measurement; time, values...)
    measure(state, measurements; time, values...)

the value of one measurement on `state`, or those of a list of them, a vector, a `Measure` or
`:φ => list`, as a vector of pairs `name => value` in their order, as `v[i]` and `v[[i]]` give
an element and a vector. An operator placed on sites gives its expectation value, an
operator of one site a vector indexed by site and a pair a matrix, as `expect1` and `expect2`
give them, made real when the operator, or each product of the pair, is Hermitian, with a
warning when the imaginary part dropped is beyond rounding, and complex otherwise; any other
value is given as it is computed, a `Check` as `(value, reference, difference)`. `time` and
the `values` are what a measurement written as a function reads as `c.time` and `c.energy`,
and what a symbol stands for, a symbol of the list without a value being refused before
anything is measured. A `Measure` keeps its preparation, which a list given as it is makes
anew at each call.

# Examples

    measure(state, Purity)
    measure(state, [X, (X, Y), Z(1) * Z(2), Purity])
    measure(state, ["E" => H, "cos" => c -> cos(2c.time), :energy]; time = 1.5, energy = -3.2)
"""
function measure(state::AbstractState, m::Measure; time::Union{Nothing, Number} = nothing,
                 values...)
    given = Dict{Symbol, Any}(values)
    if haskey(given, :state)
        error("a value given to measure cannot be named state, which c.state reads")
    end
    items = prepared(m, state)
    for it in items
        # the time is given apart from the values, and read as they are
        if it isa Symbol && !haskey(given, it) && (it ≠ :time || isnothing(time))
            error("measure was given no value of $it: measure(state, measurements; $it = …)")
        end
    end
    ops = unique(SimpleOp[ it.op for it in items if it isa SiteValues ])
    pairs = unique(Tuple{SimpleOp, SimpleOp}[ it.pair for it in items if it isa PairValues ])
    # as computed: the values are made real or complex below, under their names
    sites = isempty(ops) ? [] : expect1_computed(state, ops)
    correlations = isempty(pairs) ? [] : expect2_computed(state, pairs)
    computed = (sites = Dict(zip(ops, sites)), pairs = Dict(zip(pairs, correlations)))
    c = Context(state, time, given)
    return Pair{String, Any}[ name => entry_value(it, name, c, computed)
                              for ((name, _), it) in zip(m.entries, items) ]
end

measure(state::AbstractState,
        list::Union{AbstractArray, Pair{Symbol, <:Union{AbstractArray, Measure}}}; kwargs...) =
    measure(state, Measure(list); kwargs...)

measure(state::AbstractState, m; kwargs...) = last(only(measure(state, Measure(m); kwargs...)))
