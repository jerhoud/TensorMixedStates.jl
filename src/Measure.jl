export StateFunc, TimeFunc, Check, Measure, Trace, TraceError, Trace2, Purity, Norm, Hermiticity, HermiticityError, Renyi2, SubRenyi2
export EE, EntanglementEntropy, MutualInfoRenyi2, Mutual_Info_Renyi2, measure
export Linkdim, MaxLinkdim, MemoryUsage
export Fidelity, Overlap, Variance
export RealValue, ImaginaryValue, ComplexValue

"""
    struct StateFunc
    StateFunc(name, func)

a data type to represent a function of `State`. This is used by `measure`.
This how `Trace`, `Purity` ... are implemented.

# Examples

    Trace = StateFunc("Trace", trace)
    measure(state, Trace)
    measures = "data.dat" => Trace
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

a data type to represent a function of simulation time, the sibling of `StateFunc`. This
is what `measure` builds for a number, a string or a function given as a measurement.

Giving the name matters as soon as there are two of them: an anonymous function is named
`"func"`, and two measurements of one set may not share a name, so `t -> sin(t)` and
`t -> cos(t)` together are refused. Naming them here is the way out.

# Examples

    measures = "data" => [TimeFunc("sinus", t -> sin(t)), TimeFunc("cosinus", t -> cos(t))]
"""
struct TimeFunc
    name::String
    obs::Union{Number, Vector, Function}
end

show(io::IO, s::TimeFunc) =
    print(io, s.name)

"""
    compact_positions(positions)

write a set of site positions in a short form, contiguous runs becoming ranges, so that a
measurement made on many sites still has a readable name. Names are used as column headers
and as keys, so two different sets must keep two different names.

# Examples

    compact_positions(1:20)        # "1:20"
    compact_positions([1,2,3,7,8]) # "1:3,7:8"
"""
compact_positions(p::Int) = string(p)

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
    ObsOp(name, obs)

a data type to represent an observable defined by quantum operators. This is used by `measure`.
"""
struct ObsOp
    name::String
    obs::Vector{IndexedOp{Pure}}
end
ObsOp(name::String, o::IndexedOp{Pure}) = ObsOp(name, sumsubs(o))

"""
    struct ObsExp1
    ObsExp1(name, obs)

a data type to represent an observable applied on all sites. This is used by `measure`.
"""
struct ObsExp1
    name::String
    obs ::SimpleOp
end

"""
    struct ObsExp2
    ObsExp2(name, obs)

a data type to represent a correlation applied on all pairs of site. This is used by `measure`.
"""
struct ObsExp2
    name::String
    obs::Tuple{SimpleOp, SimpleOp}
end


"""
    struct Check
    Check(name, obs1, obs2[, tol])

a measurement that checks the equality between two measurements. It throws an error if the difference is larger than tol.

The two are compared as they are computed, before `measure` gives them their kind (see
`RealValue`), which only decides how they are written.
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

the kind of the values of a measurement, which decides how `measure` gives them

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

a measurement together with the kind of its values. `RealValue`, `ImaginaryValue` and
`ComplexValue` build one for the user, and `make_obs` one for every other measurement, with
the kind `operator_kind` finds for an operator or `value_kind` gives anything else.
"""
struct Declared
    obs
    kind::ValueKind
end

"""
    RealValue(measurement)
    ImaginaryValue(measurement)
    ComplexValue(measurement)

declare the values of a measurement real, purely imaginary or complex, for when the kind
`measure` gives it by itself is not the one wanted. A real value is written in one column,
an imaginary one as its imaginary part, in one column named `Im(name)`, and a complex one in
two, its real part then its imaginary part. The part a declaration drops is checked, and a
warning is written when it is more than rounding.

`measure` finds the kind of an operator by itself: real when it is self adjoint, imaginary
when its adjoint is its opposite, complex otherwise. The test is symbolic and only says real
or imaginary when it can prove it, so what it misses comes out complex, with a column of
rounding: `Sp(1)Sm(2) + Sm(1)Sp(2)` for instance, `Sm` being a name whose relation to `Sp`
the test does not know. A state function or a function of time is real unless declared
otherwise, and a number takes the kind of its type.

A declaration given a vector or a `Check` holds for each of its measurements, and one given
inside another is the one that holds for what it contains.

# Examples

    measures = "data" => [RealValue(Sp(1)Sm(2) + Sm(1)Sp(2)), ComplexValue(my_state_function)]
"""
RealValue(obs) = Declared(obs, real_kind)

@doc (@doc RealValue)
ImaginaryValue(obs) = Declared(obs, imaginary_kind)

@doc (@doc RealValue)
ComplexValue(obs) = Declared(obs, complex_kind)

const kind_constructors = Dict(real_kind => "RealValue", imaginary_kind => "ImaginaryValue",
                               complex_kind => "ComplexValue")

show(io::IO, d::Declared) = print(io, kind_constructors[d.kind], "(", d.obs, ")")

"""
    operator_kind(s)

the kind of the values of an operator already simplified: real when `simplify_dag` leaves it
unchanged, imaginary when it gives its opposite, complex otherwise.

Canonical forms are compared, so equal means equal, and an operator whose adjoint the rules of
`simplify_dag` cannot bring back to it comes out complex. That is the one direction a mistake
can go: a column of rounding too many, never an imaginary part lost. The rules rest on the
`OpType` of the named operators, as they do for the dissipators.
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

# an entry of a correlation matrix is the expectation value of one of three operators, the
# product on one site or the pair on two sites in either order: the Jordan-Wigner strings of
# the fermionic entries lie on the sites in between and are self adjoint. The matrix takes
# the kind all three share, complex when they differ, since a kind per entry would give its
# rows columns of their own
function pair_kind(a::SimpleOp, b::SimpleOp)
    ks = unique([ operator_kind(simplify(p)) for p in ((a * b)(1), a(1) * b(2), a(2) * b(1)) ])
    return length(ks) == 1 ? only(ks) : complex_kind
end

# the symbolic test does not see into a function, which is real unless declared otherwise.
# A constant, often the reference of a `Check`, takes the kind of its type, which is part of
# how the measurement is written and not a property of a value computed along the way
value_kind(o::TimeFunc) = o.obs isa Function ? real_kind : constant_kind(o.obs)
value_kind(_) = real_kind

constant_kind(::Complex) = complex_kind
constant_kind(x::AbstractArray) = any(y -> constant_kind(y) == complex_kind, x) ? complex_kind : real_kind
constant_kind(_) = real_kind

make_leaf(o::IndexedOp{Pure}) =
    ObsOp(obs_name(o), simplify(o))
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

make_obs(o::Union{Vector, Matrix}) = make_obs.(o)
make_obs(o::Check) =
    Check(o.name, make_obs(o.obs1), make_obs(o.obs2), o.tol)
make_obs(o::Declared) = declare(o.obs, o.kind)
# simplified once, for the measurement and for the test of its kind
function make_obs(o::IndexedOp{Pure})
    s = simplify(o)
    return Declared(ObsOp(obs_name(o), s), operator_kind(s))
end
make_obs(o::SimpleOp) =
    Declared(make_leaf(o), operator_kind(simplify(o(1))))
make_obs(o::Tuple{SimpleOp, SimpleOp}) =
    Declared(make_leaf(o), pair_kind(o...))
function make_obs(o)
    leaf = make_leaf(o)
    return Declared(leaf, value_kind(leaf))
end

declare(o::Declared, ::ValueKind) = make_obs(o)
declare(o::Union{Vector, Matrix}, kind::ValueKind) = map(x -> declare(x, kind), o)
declare(o::Check, kind::ValueKind) =
    Check(o.name, declare(o.obs1, kind), declare(o.obs2, kind), o.tol)
declare(o, kind::ValueKind) = Declared(make_leaf(o), kind)

# row by row, the order in which the lines of a matrix are written
row_major(x::AbstractVector) = x
row_major(x::AbstractMatrix) = vec(permutedims(x))

flat_measures(x::AbstractArray) = reduce(vcat, [ flat_measures(y) for y in row_major(x) ]; init = [])
flat_measures(x) = [x]

"""
    struct Measure
    Measure(args...)

a data type to hold a set of measurements, which is what a destination is given.

Building one by hand is only needed to measure several sets at once, `measure(state,
[Measure(...), Measure(...)])`, which returns one group of results per set and computes a
product shared by two of them only once. A single set is written as a plain vector, and
`output` builds these for you, one per destination.

A vector inside a set stands for its measurements, each one given on its own.
"""
struct Measure
    measures::Vector
    function Measure(obs::Vector)
        # a single name for several values would have to be a vector, which is neither a
        # column header nor a key. Inside a `Check` a vector stays one, compared element
        # by element
        m = new(flat_measures(make_obs.(obs)))
        # names become column headers and keys, so two measurements sharing one would be
        # written on top of each other. It takes a long operator, abbreviated to the same
        # text as another, to get there, so the way out is left to the caller: name the
        # measurements apart or put them in different destinations.
        ns = measure_names(m)
        dup = unique([n for n in ns if count(==(n), ns) > 1])
        if !isempty(dup)
            error("several measurements of the same set are named $(join(repr.(dup), ", ")). " *
                  "Names are used as column headers, so they must differ: split them between " *
                  "destinations, or name them explicitly.")
        end
        return m
    end
end

written_name(kind::ValueKind, name) = kind == imaginary_kind ? "Im($name)" : name

measure_names(o::Measure) = reduce(vcat, measure_names.(o.measures); init = String[])
measure_names(o::Declared) = [ written_name(o.kind, n) for n in measure_names(o.obs) ]
measure_names(o::Union{ObsOp, ObsExp1, ObsExp2, StateFunc, TimeFunc, Check}) = [o.name]
measure_names(o::Symbol) = [string(o)]
measure_names(_) = String[]

Measure(args...) = Measure([args...])

get_prods(o::Union{Vector, Matrix}) = vcat(get_prods.(o)...)
get_prods(o::Measure) = vcat(get_prods.(o.measures)...)
get_prods(o::Declared) = get_prods(o.obs)
get_prods(o::ObsOp) = o.obs
get_prods(o::Check) = [get_prods(o.obs1); get_prods(o.obs2)]
get_prods(_) = IndexedOp{Pure}[]

get_exp1(o::Union{Vector, Matrix}) = vcat(get_exp1.(o)...)
get_exp1(o::Measure) = vcat(get_exp1.(o.measures)...)
get_exp1(o::Declared) = get_exp1(o.obs)
get_exp1(o::ObsExp1) = [o.obs]
get_exp1(o::Check) = [get_exp1(o.obs1); get_exp1(o.obs2)]
get_exp1(_) = SimpleOp[]

get_exp2(o::Union{Vector, Matrix}) = vcat(get_exp2.(o)...)
get_exp2(o::Measure) = vcat(get_exp2.(o.measures)...)
get_exp2(o::Declared) = get_exp2(o.obs)
get_exp2(o::ObsExp2) = [o.obs]
get_exp2(o::Check) = [get_exp2(o.obs1); get_exp2(o.obs2)]
get_exp2(_) = Tuple{SimpleOp, SimpleOp}[]

"""
    Trace

a state function to measure the trace of the system (density matrix). See also `StateFunc` and `trace`.
"""
const Trace = StateFunc("Trace", trace)

"""
    TraceError

a state function to measure the deviation to trace 1 of the system (density matrix).
See also `StateFunc` and `Trace`.
This is good way to measure the coherence of a simulation as numerical
inaccuracies tend to change the trace of the density matrix (which should stay 1)
"""
const TraceError = StateFunc("TraceError", st -> 1. - trace(st))

"""
    Trace2
    Purity

state functions to measure the trace of the square of the density matrix. See also `StateFunc` and `trace2`.
"""
const Trace2 = StateFunc("Trace2", trace2)

@doc (@doc Trace2)
const Purity = StateFunc("Purity", trace2)

"""
    Norm

a state function to measure the norm of the state. See also `StateFunc` and `norm`.
"""
const Norm = StateFunc("Norm", norm)

"""
    Hermiticity

a state function to measure the degree of hermiticity of the density matrix.
Return 1 if density matrix is Hermitian, 0 for anti Hermitian, in between otherwise.
Sea also `StateFunc` and `hermiticity`
"""
const Hermiticity = StateFunc("Hermiticity", hermiticity)

"""
    HermiticityError

a state function to measure the deviation of the Hermiticity from 1.
See `StateFunc`, `Hermiticity` and `hermiticity`.
"""
const HermiticityError = StateFunc("HermiticityError", st -> 1. - hermiticity(st))

"""
    Renyi2

a state function to measure the Renyi-2 entropy of the system.
See also `StateFunc` and `renyi2`.
"""
const Renyi2 = StateFunc("Renyi2", renyi2)

"""
    SubRenyi2([positions...])

a state function to measure the Renyi-2 entropy of a subsystem described by the positions
given. On a pure representation this measures how much that subsystem is entangled with
the rest, and the state is first turned into its mixed representation to do it, which is
much more expensive than the other state functions.
See also `StateFunc`, `Renyi2` and `renyi2`.
"""
SubRenyi2(pos) = StateFunc("SubRenyi2($(compact_positions(pos)))", st -> renyi2(st, pos))

"""
    EntanglementEntropy(pos)
    EntanglementEntropy(pos, spectrum)

a state function to measure entanglement entropy / OSEE and associated spectrum. The cut
is on the right of `pos`, between sites `pos` and `pos + 1`, and the spectrum is that of
the reduced density matrix of the sites up to `pos`, so it is made of squared singular
values summing to one. `spectrum` is how many of them to write out.
See also `StateFunc` and `entanglement_entropy`.
"""
EntanglementEntropy(pos) = StateFunc("EntanglementEntropy($pos)",
    st-> begin
        ee, _ = entanglement_entropy(st, pos)
        return ee
    end)
EntanglementEntropy(pos, spectrum) = StateFunc("EntanglementEntropy($pos,$spectrum)",
    st-> begin
        ee, sp = entanglement_entropy(st, pos)
        return [[ee]; sp[1:min(length(sp), spectrum)]]
    end)

# the docstring goes through `@doc` rather than sitting above the call, because the macro
# expands to a toplevel block and a docstring cannot be attached to one
Base.@deprecate EE(pos) EntanglementEntropy(pos) false
Base.@deprecate EE(pos, spectrum) EntanglementEntropy(pos, spectrum) false

@doc """
    EE(pos)
    EE(pos, spectrum)

deprecated, use [`EntanglementEntropy`](@ref) instead. The label written to the output
files follows the new name, so a column that read `EE(3)` now reads
`EntanglementEntropy(3)`.
""" EE

"""
    MutualInfoRenyi2(link)
    MutualInfoRenyi2([positions...])

a state function to measure the Renyi-2 mutual information of the given subsystems.
See also `StateFunc` and `mutual_info_renyi2`.
"""
MutualInfoRenyi2(part) = StateFunc("MutualInfoRenyi2($(compact_positions(part)))", st -> mutual_info_renyi2(st, part))

# the docstring goes through `@doc` rather than sitting above the call, because the macro
# expands to a toplevel block and a docstring cannot be attached to one
Base.@deprecate Mutual_Info_Renyi2(part) MutualInfoRenyi2(part) false

@doc """
    Mutual_Info_Renyi2(part)

deprecated, use [`MutualInfoRenyi2`](@ref) instead.

It forwards to `MutualInfoRenyi2`, so the measurement is the same one, but the label written
to the output files is the new name. The deprecation warning only shows with
`--depwarn=yes`, which is what running the tests does; an ordinary run stays silent, so this
line is the notice.
""" Mutual_Info_Renyi2

"""
    reference_on(st, ref)

the reference state `ref` put on the system of the measured state `st`, weakened first to
what `st` conserves: a simulation may weaken its state after the reference was built, and
the two would no longer have the same sites. Weakening is exact, so nothing measured
changes, and it costs nothing when the two already conserve the same.
"""
reference_on(st::State, ref::State) = State(st.system, weaken(ref, symmetries(st.system)))

"""
    Fidelity(ref)

a state function to measure the fidelity with the reference state `ref`.
See also `StateFunc` and `fidelity`.

`ref` is put on the system of the state being measured, which the strictness of `fidelity`
would otherwise refuse: a measurement is written when the simulation is described, before
the system it will run on exists. It is weakened first to what that state conserves, so that
it goes on being measured after a `Weaken` phase.

# Examples

    measures = "data" => [Fidelity(ground_state), Purity]
"""
Fidelity(ref::State) = StateFunc("Fidelity", st -> fidelity(st, reference_on(st, ref)))

"""
    Overlap(ref)

a state function to measure the inner product with the reference state `ref`, that is
``\\langle ref | \\psi \\rangle`` on a pure representation. Unlike `Fidelity` it is not
normalised and it is complex, so it is written as two columns: it is declared
`ComplexValue`.
See also `StateFunc` and `inner`.

`ref` is put on the system of the state being measured, as for `Fidelity`.

# Examples

    measures = "data" => Overlap(initial_state)
"""
Overlap(ref::State) = ComplexValue(StateFunc("Overlap", st -> inner(reference_on(st, ref), st)))

"""
    Variance(hamiltonian)

a state function to measure the variance of the energy of `hamiltonian`, which is zero
exactly when the state is one of its eigenstates.
See also `StateFunc` and `variance`.

The MPO is built at every measurement, since the system the simulation runs on does not
exist when the measurement is written. That is cheap next to the two contractions the
variance itself costs, which are those of a `dmrg` sweep: ask for this in `final_measures`
or under a large `measures_period`, not at every sweep.

# Examples

    final_measures = "data" => Variance(hamiltonian)
"""
Variance(h) = StateFunc("Variance", st -> variance(h, st))

"""
    MaxLinkdim

a state function to measure the maximum bond dimension, named after the `maxlinkdim` it
measures, as the other state functions are after theirs.
See also `StateFunc` and `maxlinkdim`.
"""
const MaxLinkdim = StateFunc("MaxLinkdim", maxlinkdim)

# the docstring goes through `@doc` rather than sitting above the call, because the macro
# expands to a toplevel block and a docstring cannot be attached to one
Base.@deprecate_binding Linkdim MaxLinkdim false ", use MaxLinkdim instead."

@doc """
    Linkdim

deprecated, use [`MaxLinkdim`](@ref) instead. The label written to the output files follows
the new name.

Note that a deprecated *binding* only warns on a qualified access,
`TensorMixedStates.Linkdim`; after `using TensorMixedStates` it is silent, so this line and
the changelog are the notice.
""" Linkdim

"""
    MemoryUsage

a state function to measure the memory used by the state. See also `StateFunc`.
"""
const MemoryUsage = StateFunc("MemoryUsage", Base.summarysize)

# the part `measure` drops from a value is reported above this size, below it is rounding
const thresh_warn_dropped = 1e-6

# for a function, a kind is only a default, which a declaration overrides
dropped_hint(::Union{StateFunc, TimeFunc}) = ", ComplexValue keeps it"
dropped_hint(_) = ""

"""
    kept_value(declared, name, x, t)

the value `x` of a measurement given the kind it is declared of: its real part, its imaginary
part, or a complex number even when `x` is real, so that the values of a measurement all have
one type whatever the state. The part dropped is reported with `@warn` when it is more than
rounding, which `output` sends to the log of the simulation.
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
    # the empty value of a symbol the running algorithm does not provide keeps its type
    if isempty(x)
        return x
    end
    return map(keys(x)) do ij
        kept_value(kind, "$name $(Tuple(ij))", x[ij], t, hint)
    end
end

kept_value(::ValueKind, _, x, _, _) = x

get_val(o::Vector{Measure}, v::Dict, st::State, t::Number; kwargs...) =
    [get_val(x, v, st, t; kwargs...) for x in o]
get_val(o::Measure, v::Dict, st::State, t::Number; kwargs...) =
    [get_val(x, v, st, t; kwargs...) for x in o.measures]
function get_val(o::Declared, v::Dict, st::State, t::Number; kwargs...)
    name, x = get_val(o.obs, v, st, t; kwargs...)
    return written_name(o.kind, name) => kept_value(o, name, x, t)
end
get_val(o::Union{ObsExp1, ObsExp2}, v::Dict, ::State, ::Number; kwargs...) = o.name => v[o.obs]
get_val(o::ObsOp, v::Dict, ::State, ::Number; kwargs...) = o.name => sum(v[p] for p in o.obs)
get_val(o::TimeFunc, ::Dict, ::State, t::Number; kwargs...) =
    if o.obs isa Function
        o.name => o.obs(t)
    else
        o.name => o.obs
    end
get_val(o::StateFunc, ::Dict, st::State, ::Number; kwargs...) = o.name => o.obs(st)
get_val(o::Symbol, ::Dict, st::State, ::Number; kwargs...) =
    if haskey(kwargs, o)
        string(o) => kwargs[o]
    else
        string(o) => []
    end

part_value(o::Declared, args...; kwargs...) = last(get_val(o.obs, args...; kwargs...))
part_value(o::Union{Vector, Matrix}, args...; kwargs...) = map(x -> part_value(x, args...; kwargs...), o)
part_value(o, args...; kwargs...) = last(get_val(o, args...; kwargs...))

# a part has no name of its own on the line of its check, where `Im(name)` could not say
# that a number is an imaginary part: an imaginary part stays a complex number
function written_part(o::Declared, x, t)
    y = kept_value(o, only(measure_names(o.obs)), x, t)
    return o.kind == imaginary_kind ? on_numbers(z -> im * z, y) : y
end
written_part(o::Union{Vector, Matrix}, x, t) = map((p, y) -> written_part(p, y, t), o, x)
written_part(_, x, _) = x

on_numbers(f, x::Number) = f(x)
on_numbers(f, x::AbstractArray) = map(y -> on_numbers(f, y), x)
on_numbers(_, x) = x

function get_val(o::Check, v::Dict, st::State, t::Number; kwargs...)
    # compared as computed: the kind of a part decides how it is written, and what it drops
    # must not decide whether the check passes, a complex reference given as a function of
    # time, real unless declared otherwise, for instance
    v1 = part_value(o.obs1, v, st, t; kwargs...)
    v2 = part_value(o.obs2, v, st, t; kwargs...)
    d = norm(v1 - v2)
    if !isnothing(o.tol) && d > o.tol
        error("Check $(o.name) failed with values $v1, $v2 and difference $d")
    end
    # a vector literal would convert the three to a common type, a real reference or the
    # distance to a complex number, and change the columns they take. `map` leaves each as
    # it is and narrows the container the way reading it back from a checkpoint does
    return o.name => map(identity, Any[written_part(o.obs1, v1, t), written_part(o.obs2, v2, t), d])
end

"""
    measure(state, args[, t])
    measure(state, measure[, t])
    measure(state, [measures...][, t])

compute the requested measurements on the given state and simulation time

It is more efficient to ask all required measurements in one call

Each value is given the kind of its measurement, see `RealValue`: a real number, a complex
number, or the imaginary part of a purely imaginary one under the name `Im(name)`. `expect`,
`expect1` and `expect2` give the values as they are computed.

# Examples
    measure(state, X(1))     # compute observable X(1)
    measure(state, X)        # compute observable X on all sites
    measure(state, (X, Y))   # compute correlations XY on all pairs of sites
    measure(state, Check("check", X(1)X(2), t->sin(2t)), 0.8) # compute and check the given observable against a computed value
    measure(state, [X(2), Y, (X, Y)]) # several measures together
"""
measure(state::State, args, t::Number = 0.; kwargs...) =
    measure(state, Measure(args), t; kwargs...)

measure(state::State, m::Measure, t::Number = 0.; kwargs...) =
    measure(state, [m], t; kwargs...)[1]

function measure(state::State, m::Vector{Measure}, t::Number = 0.; kwargs...)
    vals = Dict()
    prods = collect(Set(get_prods(m)))
    if !isempty(prods)
        # make_obs simplified them once, there is nothing left for expect to normalise
        push!(vals, (prods .=> expect_norm(state, prods))...)
    end
    exp1s = collect(Set(get_exp1(m)))
    if !isempty(exp1s)
        push!(vals, (exp1s .=> expect1(state, exp1s))...)
    end
    exp2s = collect(Set(get_exp2(m)))
    if !isempty(exp2s)
        push!(vals, (exp2s .=> expect2(state, exp2s))...)
    end
    return get_val(m, vals, state, t; kwargs...)
end