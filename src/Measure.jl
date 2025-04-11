export StateFunc, TimeFunc, Check, Measure, Trace, TraceError, Trace2, Purity, Norm, EE, Linkdim, MemoryUsage, measure

"""
    struct StateFunc
    StateFunc(name, obs)

a data type to represent a function of `State`. This is used by `measure`.
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

a data type to represent a function of simulation time. This is used by `measure`.
"""
struct TimeFunc
    name::String
    obs
end

show(io::IO, s::TimeFunc) =
    print(io, s.name)

"""
    struct ObsLit
    ObsLit(name, obs)

a data type to represent an observable defined by quantum operators. This is used by `measure`.
"""
struct ObsOp
    name::String
    obs::Vector{ExprIndexed{Pure}}
end
ObsOp(name::String, o::ExprIndexed{Pure}) = ObsOp(name, sumsubs(o))

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
    tol
    Check(name, o1, o2, tol=nothing) = new(name, o1, o2, tol)
end

function make_obs(o::ExprIndexed{Pure})
    io = IOBuffer()
    print(IOContext(io, :compact => true), o)
    seek(io, 0)
    name = read(io, String)
    ObsOp(name, process(o))
end
make_obs(o::Tuple{SimpleOp, SimpleOp}) =
    ObsExp2(first(o).name * last(o).name, o)
make_obs(o::SimpleOp) =
    ObsExp1(o.name, o)
make_obs(o::Union{Number, <:Vector}) =
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

a data type to hold a set of measurements
"""
struct Measure
    measures::Vector
    Measure(obs::Vector) = new(make_obs.(obs))
end

Measure(args...) = Measure([args...])

get_prods(o::Vector{Measure}) = vcat(get_prods.(o)...)
get_prods(o::Measure) = vcat(get_prods.(o.measures)...)
get_prods(o::ObsOp) = o.obs
get_prods(o::Check) = [get_prods(o.obs1); get_prods(o.obs2)]
get_prods(_) = ExprIndexed{Pure}[]

get_exp1(o::Vector{Measure}) = vcat(get_exp1.(o)...)
get_exp1(o::Measure) = vcat(get_exp1.(o.measures)...)
get_exp1(o::ObsExp1) = [o.obs]
get_exp1(o::Check) = [get_exp1(o.obs1); get_exp1(o.obs2)]
get_exp1(_) = SimpleOp[]

get_exp2(o::Vector{Measure}) = vcat(get_exp2.(o)...)
get_exp2(o::Measure) = vcat(get_exp2.(o.measures)...)
get_exp2(o::ObsExp2) = [o.obs]
get_exp2(o::Check) = [get_exp2(o.obs1); get_exp2(o.obs2)]
get_exp2(_) = Tuple{SimpleOp, SimpleOp}[]

Trace = StateFunc("Trace", trace)
TraceError = StateFunc("TraceError", st -> 1. - trace(st))
Trace2 = StateFunc("Trace2", trace2)
Purity = StateFunc("Purity", trace2)
Norm = StateFunc("Norm", norm)
EE(pos) = StateFunc("EE($pos)",
    st-> begin
        ee, _ = entanglement_entropy(st, pos)
        return ee
    end)
EE(pos, spectre) = StateFunc("EE($pos,$spectre)",
    st-> begin
        ee, sp = entanglement_entropy(st, pos)
        return [[ee]; sp[1:min(length(sp), spectre)]]
    end)
Linkdim = StateFunc("Linkdim", maxlinkdim)
MemoryUsage = StateFunc("MemoryUsage", Base.summarysize)

get_val(o::Vector{Measure}, v::Dict, st::State, t::Number; kwargs...) =
    [get_val(x, v, st, t; kwargs...) for x in o]
get_val(o::Measure, v::Dict, st::State, t::Number; kwargs...) =
    [get_val(x, v, st, t; kwargs...) for x in o.measures]
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

# Examples
    measure(state, X(1))     # compute observable X(1)
    measure(state, X)        # compute observable X on all sites
    measure(state, (X, Y))   # compute correlations XY on all pairs of sites
    measure(state, Check("check", X(1)X(2), sin(2t)), 0.8) # compute and check the given observable against a computed value
"""
measure(state::State, args, t::Number = 0.; kwargs...) =
    measure(state, Measure(args), t; kwargs...)

measure(state::State, m::Measure, t::Number = 0.; kwargs...) =
    measure(state, [m], t; kwargs...)[1]

function measure(state::State, m::Vector{Measure}, t::Number = 0.; kwargs...)
    vals = Dict()
    prods = collect(Set(get_prods(m)))
    if !isempty(prods)
        push!(vals, (prods .=> expect(state, prods))...)
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