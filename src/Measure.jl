export StateFunc, TimeFunc, Check, Measure, Trace, TraceError, Trace2, Purity, Norm, Hermiticity, HermiticityError, Renyi2, SubRenyi2
export EE, EntanglementEntropy, MutualInfoRenyi2, Mutual_Info_Renyi2, measure
export Linkdim, MaxLinkdim, MemoryUsage
export Fidelity, Overlap, Variance

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
"""
struct Check
    name::String
    obs1
    obs2
    tol::Union{Nothing, Number}
    Check(name, o1, o2, tol=nothing) = new(name, o1, o2, tol)
end

make_obs(o::IndexedOp{Pure}) =
    ObsOp(obs_name(o), simplify(o))
make_obs(o::Union{Vector, Matrix}) = make_obs.(o)
make_obs(o::Tuple{SimpleOp, SimpleOp}) =
    ObsExp2(obs_name(first(o)) * obs_name(last(o)), o)
make_obs(o::SimpleOp) =
    ObsExp1(obs_name(o), o)
make_obs(o::Number) =
    TimeFunc(string(o), o)
make_obs(o::String) =
    TimeFunc(o, [])
make_obs(o::Function) =
    TimeFunc("func", o)
make_obs(o::Check) =
    Check(o.name, make_obs(o.obs1), make_obs(o.obs2), o.tol)
make_obs(o) = o

"""
    struct Measure
    Measure(args...)

a data type to hold a set of measurements, which is what a destination is given.

Building one by hand is only needed to measure several sets at once, `measure(state,
[Measure(...), Measure(...)])`, which returns one group of results per set and computes a
product shared by two of them only once. A single set is written as a plain vector, and
`output` builds these for you, one per destination.
"""
struct Measure
    measures::Vector
    function Measure(obs::Vector)
        m = new(make_obs.(obs))
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

measure_names(o::Measure) = reduce(vcat, measure_names.(o.measures); init = String[])
measure_names(o::Union{Vector, Matrix}) = reduce(vcat, measure_names.(o); init = String[])
measure_names(o::Union{ObsOp, ObsExp1, ObsExp2, StateFunc, TimeFunc, Check}) = [o.name]
measure_names(o::Symbol) = [string(o)]
measure_names(_) = String[]

Measure(args...) = Measure([args...])

get_prods(o::Union{Vector, Matrix}) = vcat(get_prods.(o)...)
get_prods(o::Measure) = vcat(get_prods.(o.measures)...)
get_prods(o::ObsOp) = o.obs
get_prods(o::Check) = [get_prods(o.obs1); get_prods(o.obs2)]
get_prods(_) = IndexedOp{Pure}[]

get_exp1(o::Union{Vector, Matrix}) = vcat(get_exp1.(o)...)
get_exp1(o::Measure) = vcat(get_exp1.(o.measures)...)
get_exp1(o::ObsExp1) = [o.obs]
get_exp1(o::Check) = [get_exp1(o.obs1); get_exp1(o.obs2)]
get_exp1(_) = SimpleOp[]

get_exp2(o::Union{Vector, Matrix}) = vcat(get_exp2.(o)...)
get_exp2(o::Measure) = vcat(get_exp2.(o.measures)...)
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
normalised and it is complex, so it is written as two columns.
See also `StateFunc` and `inner`.

`ref` is put on the system of the state being measured, as for `Fidelity`.

# Examples

    measures = "data" => Overlap(initial_state)
"""
Overlap(ref::State) = StateFunc("Overlap", st -> inner(reference_on(st, ref), st))

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

get_val(o::Vector{Measure}, v::Dict, st::State, t::Number; kwargs...) =
    [get_val(x, v, st, t; kwargs...) for x in o]
get_val(o::Measure, v::Dict, st::State, t::Number; kwargs...) =
    [get_val(x, v, st, t; kwargs...) for x in o.measures]
function get_val(o::Union{Vector, Matrix}, v::Dict, st::State, t::Number; kwargs...)
    gvs = [get_val(x, v, st, t; kwargs...) for x in o]
    return first.(gvs) => last.(gvs)
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

function get_val(o::Check, v::Dict, st::State, t::Number; kwargs...)
    v1 = last(get_val(o.obs1, v, st, t; kwargs...))
    v2 = last(get_val(o.obs2, v, st, t; kwargs...))
    d = norm(v1 - v2)
    if !isnothing(o.tol) && d > o.tol
        error("Check $(o.name) failed with values $v1, $v2 and difference $d")
    end
    return o.name => [v1, v2, d]
end

"""
    measure(state, args[, t])
    measure(state, measure[, t])
    measure(state, [measures...][, t])

compute the requested measurements on the given state and simulation time

It is more efficient to ask all required measurements in one call

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