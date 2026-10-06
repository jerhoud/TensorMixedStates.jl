# The measurements: their types (Trace, EntanglementEntropy, Check...), their names and the kind
# of their values, and measure, which computes a whole set of them at once.

export StateFunc, TimeFunc, Check, Measure, Trace, TraceError, Trace2, Purity, Norm, Hermiticity, HermiticityError, Renyi2, SubRenyi2
export EntanglementEntropy, MutualInfoRenyi2, measure
export MaxLinkdim, MemoryUsage
export Fidelity, Overlap, Variance
export RealValue, ImaginaryValue, ComplexValue

"""
    struct StateFunc
    StateFunc(name, func)

a measurement given by a function of the state, written under `name`. The function returns a
number, or a vector or a matrix of numbers, taken as real unless declared otherwise, see
`RealValue`.

# Examples

    halftrace = StateFunc("HalfTrace", st -> trace(st) / 2)
    measure(state, halftrace)
    measurements = "data.dat" => halftrace
"""
struct StateFunc
    name::String
    obs::Function
end

show(io::IO, s::StateFunc) =
    print(io, s.name)

"""
    struct TimeFunc
    TimeFunc(name, obs)

a measurement given by a function of the simulation time, or by a constant, written under
`name`. A number or a function given as a measurement becomes one, and a string a label with
no value. A function given as it is is named `"func"`: two of them in one set need a
`TimeFunc` each, two measurements of a set not sharing a name.

# Examples

    measurements = "data" => [TimeFunc("sinus", t -> sin(t)), TimeFunc("cosinus", t -> cos(t))]
"""
struct TimeFunc
    name::String
    obs::Union{Number, Vector, Function}
end

show(io::IO, s::TimeFunc) =
    print(io, s.name)

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
    struct ObsOp
    ObsOp(name, op, obs)

an operator measurement: `op` as written, for the `expect` of a representation of one's own,
and `obs` its terms, simplified and compacted, for the `expect_norm` of a `State`
"""
struct ObsOp
    name::String
    op::IndexedOp{Pure}
    obs::Vector{IndexedOp{Pure}}
end
ObsOp(name::String, op::IndexedOp{Pure}, s::IndexedOp{Pure}) = ObsOp(name, op, sumsubs(s))

"""
    struct ObsExp1
    ObsExp1(name, obs)

a one site operator measured on every site with `expect1`, for `measure`.
"""
struct ObsExp1
    name::String
    obs ::SimpleOp
end

"""
    struct ObsExp2
    ObsExp2(name, obs)

a pair of one site operators whose correlations are measured on every pair of sites with
`expect2`, for `measure`.
"""
struct ObsExp2
    name::String
    obs::Tuple{SimpleOp, SimpleOp}
end


"""
    struct Check
    Check(name, obs1, obs2[, tol])

a measurement comparing two measurements: its value is theirs followed by the norm of their
difference, and an error is raised when that difference exceeds `tol`, if given. They are
compared as computed, before their kind (see `RealValue`) is applied.

# Examples

    measure(state, Check("check", X(1)X(2), t -> sin(2t)), 0.8)
    measurements = "data" => Check("trace", Trace, 1., 1e-6)
"""
struct Check
    name::String
    obs1
    obs2
    tol::Union{Nothing, Number}
    Check(name, o1, o2, tol=nothing) = new(name, o1, o2, tol)
end

"""
    @enum ValueKind

the kind of the values of a measurement, which decides how `measure` gives them.

# Enumeration values

- `real_kind`: a real number, its imaginary part being rounding
- `imaginary_kind`: a purely imaginary number, given as its imaginary part under the name
  `Im(name)`
- `complex_kind`: a complex number
"""
@enum ValueKind real_kind imaginary_kind complex_kind

"""
    struct Declared
    Declared(obs, kind)

a measurement together with the kind of its values, built by `RealValue`, `ImaginaryValue`,
`ComplexValue`, and `make_obs` for every other measurement
"""
struct Declared
    obs
    kind::ValueKind
end

"""
    RealValue(measurement)
    ImaginaryValue(measurement)
    ComplexValue(measurement)

declare the values of a measurement real, purely imaginary or complex, when the kind
`measure` gives it by itself is not the one wanted. A real value is written in one column,
an imaginary one as its imaginary part in a column named `Im(name)`, and a complex one in
two, real then imaginary part. A warning is written when the part dropped is more than
rounding.

By itself, an operator is real when it is self adjoint, imaginary when its adjoint is its
opposite, and complex when the symbolic test cannot prove either, as for
`Sp(1)Sm(2) + Sm(1)Sp(2)`. A state function or a function of time is real, and a number takes
the kind of its type.

A declaration holds for each measurement of a vector or a `Check`, one inside it keeping its
own.

# Examples

    measurements = "data" => [RealValue(Sp(1)Sm(2) + Sm(1)Sp(2)), ComplexValue(my_state_function)]
"""
RealValue(obs) = Declared(obs, real_kind)

"""
    ImaginaryValue(measurement)

declare the values of a measurement purely imaginary, see `RealValue`.
"""
ImaginaryValue(obs) = Declared(obs, imaginary_kind)

"""
    ComplexValue(measurement)

declare the values of a measurement complex, see `RealValue`.
"""
ComplexValue(obs) = Declared(obs, complex_kind)

"""
    kind_constructors

the name of the function declaring each `ValueKind`, for printing a `Declared`.
"""
const kind_constructors = Dict(real_kind => "RealValue", imaginary_kind => "ImaginaryValue",
                               complex_kind => "ComplexValue")

show(io::IO, d::Declared) = print(io, kind_constructors[d.kind], "(", d.obs, ")")

"""
    operator_kind(s)

the kind of the values of an operator already simplified: real when `simplify_dag` leaves it
unchanged, imaginary when it gives its opposite, complex otherwise. Comparing canonical forms,
it can only err towards complex, never losing an imaginary part.
"""
function operator_kind(s)
    d = simplify_dag(s)
    if d == s
        return real_kind
    elseif d == simplify(-s)
        return imaginary_kind
    else
        return complex_kind
    end
end

"""
    pair_kind(a, b)

the kind of the correlation matrix of `a` and `b`: the one shared by `a * b` on one site and
the pair on two sites in either order, complex when they differ. The Jordan-Wigner strings of
fermionic entries, in between and self adjoint, do not change it.
"""
function pair_kind(a::SimpleOp, b::SimpleOp)
    ks = unique([ operator_kind(simplify(p)) for p in ((a * b)(1), a(1) * b(2), a(2) * b(1)) ])
    return length(ks) == 1 ? only(ks) : complex_kind
end

"""
    value_kind(leaf)

the kind of the values of a measurement that is not an operator: real for a function, and
for a constant that of its type, see `constant_kind`, so that how it is written does not
depend on a value
"""
value_kind(o::TimeFunc) = o.obs isa Function ? real_kind : constant_kind(o.obs)
value_kind(_) = real_kind

"""
    constant_kind(x)

the kind of a constant measurement: complex for a complex number or an array holding one,
real otherwise.
"""
constant_kind(::Complex) = complex_kind
constant_kind(x::AbstractArray) = any(y -> constant_kind(y) == complex_kind, x) ? complex_kind : real_kind
constant_kind(_) = real_kind

"""
    obs_op(o, s)

the `ObsOp` measuring the operator `o`, given simplified as `s`, its terms of several sites
compacted, see `compact`
"""
obs_op(o::IndexedOp{Pure}, s::IndexedOp{Pure}) =
    ObsOp(obs_name(o), o, compact_simplified(removeMulti(s), rounding_tol, "expect"))

"""
    make_leaf(o)

the measurement `o` stands for: an `ObsOp` for an operator placed on sites, an `ObsExp1` for
a one site operator, an `ObsExp2` for a pair of them, a `TimeFunc` for a number, a string or
a function, and `o` itself otherwise.
"""
make_leaf(o::IndexedOp{Pure}) =
    obs_op(o, simplify(o))
make_leaf(o::Tuple{SimpleOp, SimpleOp}) =
    ObsExp2(obs_name(first(o)) * obs_name(last(o)), o)
make_leaf(o::SimpleOp) =
    ObsExp1(obs_name(o), o)
make_leaf(o::Number) =
    TimeFunc(string(o), o)
make_leaf(o::String) =
    TimeFunc(o, [])
make_leaf(o::Function) =
    TimeFunc("func", o)
make_leaf(o) = o

"""
    make_obs(o)

the measurement `o` made ready for `measure`: its leaf, see `make_leaf`, in a `Declared` with
the kind declared or found by `operator_kind`, `pair_kind` or `value_kind`. An array, a range
included, becomes a `Vector` or a `Matrix` of them, and a `Check` a `Check` of them.
"""
make_obs(o::AbstractArray) = make_obs.(collect(o))
make_obs(o::Check) =
    Check(o.name, make_obs(o.obs1), make_obs(o.obs2), o.tol)
make_obs(o::Declared) = declare(o.obs, o.kind)
# the kind is taken before compacting: the adjoint of a compacted operator may swap its
# channels, which operator_kind, comparing structures, would take as complex
function make_obs(o::IndexedOp{Pure})
    s = simplify(o)
    return Declared(obs_op(o, s), operator_kind(s))
end
make_obs(o::SimpleOp) =
    Declared(make_leaf(o), operator_kind(simplify(o(1))))
make_obs(o::Tuple{SimpleOp, SimpleOp}) =
    Declared(make_leaf(o), pair_kind(o...))
function make_obs(o)
    leaf = make_leaf(o)
    return Declared(leaf, value_kind(leaf))
end

"""
    declare(o, kind)

the measurement `o` made ready as by `make_obs` with the given kind, for each measurement of
an array or a `Check`, a declaration inside `o` keeping its own
"""
declare(o::Declared, ::ValueKind) = make_obs(o)
declare(o::AbstractArray, kind::ValueKind) = map(x -> declare(x, kind), collect(o))
declare(o::Check, kind::ValueKind) =
    Check(o.name, declare(o.obs1, kind), declare(o.obs2, kind), o.tol)
declare(o, kind::ValueKind) = Declared(make_leaf(o), kind)

"""
    row_major(x)

the elements of the vector or matrix `x` in a vector, row by row, as a matrix is written
"""
row_major(x::AbstractVector) = x
row_major(x::AbstractMatrix) = vec(permutedims(x))

"""
    flat_measurements(x)

the measurements of `x` in a flat vector, nested arrays being unrolled, a matrix row by row.
"""
flat_measurements(x::AbstractArray) =
    reduce(vcat, [ flat_measurements(y) for y in row_major(x) ]; init = [])
flat_measurements(x) = [x]

"""
    written_name(kind, name)

the name a value of the given kind is written under: `Im(name)` for an imaginary one, `name`
otherwise.
"""
written_name(kind::ValueKind, name) = kind == imaginary_kind ? "Im($name)" : name

"""
    measure_names(o)

the names the measurements of `o` are written under, a `Symbol` being named after itself.
"""
measure_names(o::Declared) = [ written_name(o.kind, n) for n in measure_names(o.obs) ]
measure_names(o::Union{ObsOp, ObsExp1, ObsExp2, StateFunc, TimeFunc, Check}) = [o.name]
measure_names(o::Symbol) = [string(o)]
measure_names(_) = String[]

"""
    struct Measure
    Measure(args...)

a set of measurements, what a destination is given, a vector inside it standing for its
measurements. Two measurements of a set may not share a name, names being column headers.

Building one is only needed to measure several sets at once,
`measure(state, [Measure(...), Measure(...)])`, which returns one group of results per set and
computes a product shared by them only once; `output` builds them, one per destination.
"""
struct Measure
    measurements::Vector
    function Measure(obs::Vector)
        # flattened, a name standing for one value; a vector inside a `Check` stays one
        measurements = flat_measurements(make_obs.(obs))
        # two long operators may be abbreviated to the same name
        ns = reduce(vcat, measure_names.(measurements); init = String[])
        dup = unique([n for n in ns if count(==(n), ns) > 1])
        if !isempty(dup)
            error("several measurements of the same set are named $(join(repr.(dup), ", ")). " *
                  "Names are used as column headers, so they must differ: split them between " *
                  "destinations, or name them explicitly.")
        end
        return new(measurements)
    end
end

Measure(args...) = Measure([args...])

# a set inside a set stands for its measurements, which are already made
make_obs(o::Measure) = o.measurements


"""
    leaves(T, o)

the measurements of type `T` in `o`, found through arrays, sets, declarations and checks.
"""
leaves(T, o::Union{Vector, Matrix}) = reduce(vcat, [ leaves(T, x) for x in o ]; init = T[])
leaves(T, o::Measure) = leaves(T, o.measurements)
leaves(T, o::Declared) = leaves(T, o.obs)
leaves(T, o::Check) = [ leaves(T, o.obs1); leaves(T, o.obs2) ]
leaves(T, o) = o isa T ? T[o] : T[]

"""
    get_prods(state, o)
    get_exp1(o)
    get_exp2(o)

what the measurements of `o` ask for, for `measure` to compute each kind at once: the terms
of the operator measurements on a `State`, the operators as written on another representation,
see `prod_values`, the operators of the `ObsExp1` for `expect1`, and the pairs of the `ObsExp2`
for `expect2`
"""
get_prods(::State, o) = reduce(vcat, [ l.obs for l in leaves(ObsOp, o) ]; init = IndexedOp{Pure}[])
get_prods(::AbstractState, o) = IndexedOp{Pure}[ l.op for l in leaves(ObsOp, o) ]
get_exp1(o) = SimpleOp[ l.obs for l in leaves(ObsExp1, o) ]
get_exp2(o) = Tuple{SimpleOp, SimpleOp}[ l.obs for l in leaves(ObsExp2, o) ]

"""
    prod_values(state, prods)

the values of what `get_prods` gives: on a `State`, by `expect_norm` on the terms already
simplified, divided by the trace as `expect` divides, and on another representation by its
`expect`
"""
prod_values(state::State, terms) = expect_norm(state, terms) ./ real(trace(state))
prod_values(state::AbstractState, ops) = expect(state, ops)

"""
    Trace

a state function measuring the trace of the density matrix, see `trace` and `StateFunc`. Only
its real part is given: `ComplexValue(Trace)` keeps the imaginary part numerical errors may
add, see `TraceError`.
"""
const Trace = StateFunc("Trace", trace)

"""
    TraceError

a state function measuring `1 - trace(state)`, see `Trace`, a good check of the accuracy of a
simulation, numerical inaccuracies tending to move the trace. Only its real part is given:
`ComplexValue(TraceError)` keeps the imaginary part, which comes from the anti-Hermitian part
of the density matrix, see `HermiticityError`.
"""
const TraceError = StateFunc("TraceError", st -> 1. - trace(st))

"""
    Trace2
    Purity

state functions measuring the purity ``\\mathrm{tr}(\\rho^2)``, see `trace2` and `StateFunc`.
"""
const Trace2 = StateFunc("Trace2", trace2)

"""
    Purity

the state function `Trace2` under the name `Purity`, see `Trace2`.
"""
const Purity = StateFunc("Purity", trace2)

"""
    Norm

a state function measuring the norm of the state, see `norm` and `StateFunc`.
"""
const Norm = StateFunc("Norm", norm)

"""
    Hermiticity

a state function measuring how Hermitian the density matrix is, from 0 when it is
anti-Hermitian to 1 when it is Hermitian, see `hermiticity` and `StateFunc`.
"""
const Hermiticity = StateFunc("Hermiticity", hermiticity)

"""
    HermiticityError

a state function measuring `1 - hermiticity(state)`, see `Hermiticity`: the squared norm of the
anti-Hermitian part of the density matrix, relative to that of the whole. When the exact state
is Hermitian, as under a Lindbladian, gates and noisy gates, and for a thermal or a steady
state, its square root is a lower bound of the error relative to the norm of the state, a
criterion of convergence, see [Checking the accuracy](@ref).
"""
const HermiticityError = StateFunc("HermiticityError", st -> 1. - hermiticity(st))

"""
    Renyi2

a state function measuring the Rényi-2 entropy of the state, see `renyi2` and `StateFunc`.
"""
const Renyi2 = StateFunc("Renyi2", renyi2)

"""
    SubRenyi2(cut)
    SubRenyi2([positions...])

a state function measuring the Rényi-2 entropy of the sites at `positions`, see `renyi2`. An
integer `cut` stands for the sites `1:cut` and is named after them, the site `k` alone being
`SubRenyi2([k])`. On a pure representation, a cut is as cheap as `EntanglementEntropy`, while
positions mix the state first, which is much more expensive.

# Examples

    measurements = "data" => [SubRenyi2(3), SubRenyi2([1, 4])]
"""
SubRenyi2(pos) = StateFunc("SubRenyi2($(compact_positions(pos)))", st -> renyi2(st, pos))

SubRenyi2(cut::Int) =
    StateFunc("SubRenyi2($(compact_positions(1:cut)))", st -> renyi2(st, cut))

"""
    EntanglementEntropy(cut)
    EntanglementEntropy(cut, spectrum)

a state function measuring the entanglement entropy, or the OSEE on a mixed representation,
across the cut between sites `cut` and `cut + 1`, see `entanglement_entropy`. Given
`spectrum`, the entropy is followed by the first `spectrum` values of the spectrum, padded
with zeros beyond the bond dimension of the cut, so that every row has the same width.

# Examples

    measurements = "data" => [EntanglementEntropy(3), EntanglementEntropy(5, 4)]
"""
EntanglementEntropy(cut) = StateFunc("EntanglementEntropy($cut)",
    st-> begin
        ee, _ = entanglement_entropy(st, cut)
        return ee
    end)
EntanglementEntropy(cut, spectrum) = StateFunc("EntanglementEntropy($cut,$spectrum)",
    st-> begin
        ee, sp = entanglement_entropy(st, cut)
        return [[ee]; sp[1:min(length(sp), spectrum)]; zeros(max(0, spectrum - length(sp)))]
    end)

"""
    MutualInfoRenyi2(cut)
    MutualInfoRenyi2([positions...])

a state function measuring the Rényi-2 mutual information between the sites at `positions`
and the rest, see `mutual_info_renyi2`. An integer `cut` stands for the sites `1:cut` and is
named after them: `MutualInfoRenyi2(1:3)` for `cut = 3`, `MutualInfoRenyi2(1,2)` for
`cut = 2`.
"""
MutualInfoRenyi2(part) = StateFunc("MutualInfoRenyi2($(compact_positions(part)))", st -> mutual_info_renyi2(st, part))

MutualInfoRenyi2(cut::Int) =
    StateFunc("MutualInfoRenyi2($(compact_positions(1:cut)))", st -> mutual_info_renyi2(st, cut))

"""
    reference_on(st, ref)

the reference state `ref` put on the system of the measured state `st`, weakened first to
what `st` conserves, a simulation possibly weakening its state after `ref` was built; a
reference conserving less than `st` is refused
"""
function reference_on(st::AbstractState, ref::State)
    target = symmetries(st.system)
    source = symmetries(ref.system)
    # the reference rather than the measured state, which would be converted at every
    # measurement
    for (name, strong) in target.names
        k = findfirst(q -> q[1] == name, source.names)
        if isnothing(k) || (strong && !source.names[k][2])
            error("the reference of Fidelity or Overlap conserves less than the measured state, " *
                  "which conserves $target: give a reference that conserves as much, or weaken " *
                  "the state with a Weaken phase")
        end
    end
    return State(st.system, weaken(ref, target))
end

"""
    Fidelity(ref)

a state function measuring the fidelity with the reference state `ref`, see `fidelity`. `ref`
is put on the system of the measured state, weakened to what it conserves, so that it can still
be measured after a `Weaken` phase. Between two mixed representations it is refused, see
`fidelity`, at the first measurement.

# Examples

    measurements = "data" => [Fidelity(ground_state), Purity]
"""
Fidelity(ref::State) = StateFunc("Fidelity", st -> fidelity(st, reference_on(st, ref)))

"""
    Overlap(ref)

a state function measuring the inner product with the reference state `ref`,
``\\langle ref | \\psi \\rangle`` on a pure representation and
``\\mathrm{tr}(ref^\\dagger \\rho)`` on a mixed one, see `inner`. Unlike `Fidelity` it is not
normalised, and it is a `ComplexValue`, written in two columns. `ref` is put on the system of
the measured state as for `Fidelity`, and a reference of the other representation is refused
at the first measurement.

# Examples

    measurements = "data" => Overlap(initial_state)
"""
Overlap(ref::State) = ComplexValue(StateFunc("Overlap", st -> inner(reference_on(st, ref), st)))

"""
    Variance(hamiltonian)

a state function measuring the variance of the energy of `hamiltonian`, zero exactly when
the state is one of its eigenstates, see `variance`. A mixed representation is refused at the
first measurement. It costs a `dmrg` sweep: ask for it in `final_measurements` or under
a large `measurements_period`.

# Examples

    final_measurements = "data" => Variance(hamiltonian)
"""
Variance(h) = StateFunc("Variance", st -> variance(st, h))

"""
    MaxLinkdim

a state function measuring the maximum bond dimension of the state, see `maxlinkdim`.
"""
const MaxLinkdim = StateFunc("MaxLinkdim", maxlinkdim)

"""
    MemoryUsage

a state function measuring the memory the state occupies, in bytes, including the caches
filled by the measurements already made on it, so that it depends on them.
"""
const MemoryUsage = StateFunc("MemoryUsage", Base.summarysize)

"""
    thresh_warn_dropped

the size above which the part `measure` drops from a value is reported, below it being
rounding.
"""
const thresh_warn_dropped = 1e-6

"""
    dropped_hint(o)

the end of the warning on a dropped part: for a state function or a function of time, whose
kind is only a default, the declaration that keeps it.
"""
dropped_hint(::Union{StateFunc, TimeFunc}) = ", ComplexValue keeps it"
dropped_hint(_) = ""

"""
    kept_value(declared, name, x, t)

the value `x` of a measurement given its kind: its real part, its imaginary part, or a complex
number even when `x` is real, so that the values of a measurement have one type whatever the
state. A dropped part larger than rounding is reported with `@warn`.
"""
kept_value(o::Declared, name, x, t) = kept_value(o.kind, name, x, t, dropped_hint(o.obs))

function kept_value(kind::ValueKind, name, x::Number, t, hint)
    if kind == complex_kind
        return complex(x)
    end
    kept, dropped, part = kind == real_kind ? (real(x), imag(x), "imaginary") : (imag(x), real(x), "real")
    if abs(dropped) > thresh_warn_dropped
        @warn("large $part part: time $t, $name " *
              @sprintf("%8.1e (rel %8.1e)", dropped, dropped / abs(x)) * hint)
    end
    return kept
end

function kept_value(kind::ValueKind, name, x::AbstractArray, t, hint)
    # the empty value of a symbol not given keeps its type
    if isempty(x)
        return x
    end
    return map(keys(x)) do ij
        kept_value(kind, "$name $(Tuple(ij))", x[ij], t, hint)
    end
end

kept_value(::ValueKind, _, x, _, _) = x

"""
    get_val(o, v, state, t; kwargs...)

the pairs `name => value` of the measurements of `o`, grouped as in `o`: `v` holds the
expectation values `measure` has computed, `t` is the time, and `kwargs` give the values of
the `Symbol` measurements, empty for one it does not give.
"""
get_val(o::Vector{Measure}, v::Dict, st::AbstractState, t::Number; kwargs...) =
    [get_val(x, v, st, t; kwargs...) for x in o]
get_val(o::Measure, v::Dict, st::AbstractState, t::Number; kwargs...) =
    [get_val(x, v, st, t; kwargs...) for x in o.measurements]
function get_val(o::Declared, v::Dict, st::AbstractState, t::Number; kwargs...)
    name, x = get_val(o.obs, v, st, t; kwargs...)
    return written_name(o.kind, name) => kept_value(o, name, x, t)
end
get_val(o::Union{ObsExp1, ObsExp2}, v::Dict, ::AbstractState, ::Number; kwargs...) =
    o.name => v[o.obs]
get_val(o::ObsOp, v::Dict, ::State, ::Number; kwargs...) =
    o.name => sum(v[p] for p in o.obs)
get_val(o::ObsOp, v::Dict, ::AbstractState, ::Number; kwargs...) =
    o.name => v[o.op]
get_val(o::TimeFunc, ::Dict, ::AbstractState, t::Number; kwargs...) =
    if o.obs isa Function
        o.name => o.obs(t)
    else
        o.name => o.obs
    end
get_val(o::StateFunc, ::Dict, st::AbstractState, ::Number; kwargs...) = o.name => o.obs(st)
get_val(o::Symbol, ::Dict, st::AbstractState, ::Number; kwargs...) =
    if haskey(kwargs, o)
        string(o) => kwargs[o]
    else
        string(o) => []
    end

"""
    part_value(o, v, state, t; kwargs...)

the value of a part of a `Check`, as computed: without its name and before its kind applies.
"""
part_value(o::Declared, args...; kwargs...) = last(get_val(o.obs, args...; kwargs...))
part_value(o::Union{Vector, Matrix}, args...; kwargs...) = map(x -> part_value(x, args...; kwargs...), o)
part_value(o, args...; kwargs...) = last(get_val(o, args...; kwargs...))

"""
    written_part(o, x, t)

the value `x` of a part of a `Check` as written given its kind, except that an imaginary part
stays a complex number, a part having no `Im(name)` of its own on the line of its check
"""
function written_part(o::Declared, x, t)
    y = kept_value(o, only(measure_names(o.obs)), x, t)
    return o.kind == imaginary_kind ? on_numbers(z -> im * z, y) : y
end
written_part(o::Union{Vector, Matrix}, x, t) = map((p, y) -> written_part(p, y, t), o, x)
written_part(_, x, _) = x

"""
    on_numbers(f, x)

`f` applied to every number of `x`, through nested arrays, anything else being left as it is.
"""
on_numbers(f, x::Number) = f(x)
on_numbers(f, x::AbstractArray) = map(y -> on_numbers(f, y), x)
on_numbers(_, x) = x

function get_val(o::Check, v::Dict, st::AbstractState, t::Number; kwargs...)
    # compared as computed: what the kind of a part drops must not decide whether the check
    # passes
    v1 = part_value(o.obs1, v, st, t; kwargs...)
    v2 = part_value(o.obs2, v, st, t; kwargs...)
    # a symbol not given here has an empty value, which nothing compares to
    if isempty(v1) || isempty(v2)
        if !isnothing(o.tol)
            error("Check $(o.name) has nothing to compare, a symbol it measures not being given here")
        end
        return o.name => Any[v1, v2, []]
    end
    d = norm(v1 - v2)
    if !isnothing(o.tol) && d > o.tol
        error("Check $(o.name) failed with values $v1, $v2 and difference $d")
    end
    # not a vector literal, which would convert the three to a common type and change their
    # columns: `map` narrows the container as reading it back from a checkpoint does
    return o.name => map(identity, Any[written_part(o.obs1, v1, t), written_part(o.obs2, v2, t), d])
end

"""
    measure(state, measurements[, t]; kwargs...)
    measure(state, ::Measure[, t]; kwargs...)
    measure(state, ::Vector{Measure}[, t]; kwargs...)

the measurements asked for on `state` at simulation time `t`, 0 by default: a measurement or
a vector of them gives a vector of pairs `name => value`, and a vector of `Measure` one such
vector per set. Asking for all measurements in one call is more efficient.

Each value has the kind of its measurement, see `RealValue`, while `expect`, `expect1` and
`expect2` give values as computed. A `Symbol` takes its value from the keyword argument of
that name, and is empty without one: the phases pass `sweep`, `energy` for a dmrg search, and
`beta` and `log_trace` for `thermal_state`.

# Examples

    measure(state, X(1))     # compute observable X(1)
    measure(state, X)        # compute observable X on all sites
    measure(state, (X, Y))   # compute correlations XY on all pairs of sites
    measure(state, Check("check", X(1)X(2), t->sin(2t)), 0.8) # compute and check the given observable against a computed value
    measure(state, [X(2), Y, (X, Y)]) # several measurements together
"""
measure(state::AbstractState, args, t::Number = 0.; kwargs...) =
    measure(state, Measure(args), t; kwargs...)

measure(state::AbstractState, m::Measure, t::Number = 0.; kwargs...) =
    measure(state, [m], t; kwargs...)[1]

function measure(state::AbstractState, m::Vector{Measure}, t::Number = 0.; kwargs...)
    vals = Dict()
    for (items, compute) in ((get_prods(state, m), prod_values), (get_exp1(m), expect1),
                             (get_exp2(m), expect2))
        u = unique(items)
        if !isempty(u)
            push!(vals, (u .=> compute(state, u))...)
        end
    end
    return get_val(m, vals, state, t; kwargs...)
end