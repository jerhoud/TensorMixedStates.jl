export StateFunc, TimeFunc, Check, Measure, Trace, Trace2, Purity, Norm, EE, Linkdim, MemoryUsage, measure

struct StateFunc
    name
    obs
end

struct TimeFunc
    name
    obs
end

struct ObsLit
    name
    obs
end

struct ObsExp1
    name
    obs 
end

struct ObsExp2
    name
    obs
end

make_obs(o::SumLit) =
    ObsLit(string(o), insertFfactors(o).ps)
make_obs(o::ProdLit) =
    ObsLit(string(o), [insertFfactors(o)])
make_obs(o::Tuple{<:Function, <:Function}) =
    ObsExp2(first(o)()*last(o)(), o)
make_obs(o::Function) =
    ObsExp1(o(), o)
make_obs(o::Union{Number, String}) =
    TimeFunc(string(o), [])
make_obs(o) = o

wrap_check(o) = make_obs(o)
wrap_check(o::Number) = TimeFunc("$o", o)
wrap_check(o::Vector) = TimeFunc("$o", o)
wrap_check(o::Function) =
    try
        ObsExp1(o(), o)
    catch
        TimeFunc("func", o)
    end

struct Check
    name
    obs1
    obs2
    tol
    Check(name, o1, o2, tol=nothing) =
        new(name, wrap_check(o1), wrap_check(o2), tol)
end

struct Measure
    measures::Vector
    Measure(obs::Vector) = new(make_obs.(obs))
end

Measure(args...) = Measure([args...])

get_prods(o::Vector{Measure}) = vcat(get_prods.(o)...)
get_prods(o::Measure) = vcat(get_prods.(o.measures)...)
get_prods(o::ObsLit) = o.obs
get_prods(o::Check) = [get_prods(o.obs1); get_prods(o.obs2)]
get_prods(_) = []

get_exp1(o::Vector{Measure}) = vcat(get_exp1.(o)...)
get_exp1(o::Measure) = vcat(get_exp1.(o.measures)...)
get_exp1(o::ObsExp1) = [o.obs]
get_exp1(o::Check) = [get_exp1(o.obs1); get_exp1(o.obs2)]
get_exp1(_) = []

get_exp2(o::Vector{Measure}) = vcat(get_exp2.(o)...)
get_exp2(o::Measure) = vcat(get_exp2.(o.measures)...)
get_exp2(o::ObsExp2) = [o.obs]
get_exp2(o::Check) = [get_exp2(o.obs1); get_exp2(o.obs2)]
get_exp2(_) = []

Trace = StateFunc("Trace", trace)
Trace2 = StateFunc("Trace2", trace2)
Purity = StateFunc("Purity", trace2)
Norm = StateFunc("Norm", norm)
EE(pos) = StateFunc("EE($pos)",
    st-> begin
        ee, _ = entanglement_entropy!(st, pos)
        return ee
    end)
EE(pos, spectre) = StateFunc("EE($pos,$spectre)",
    st-> begin
        ee, sp = entanglement_entropy!(st, pos)
        return [[ee]; sp[1:min(length(sp), spectre)]]
    end)
Linkdim = StateFunc("Linkdim", maxlinkdim)
MemoryUsage = StateFunc("MemoryUsage", Base.summarysize)

get_val(o::Vector{Measure}, v::Dict, st::State, t::Number; kwargs...) =
    [get_val(x, v, st, t; kwargs...) for x in o]
get_val(o::Measure, v::Dict, st::State, t::Number; kwargs...) =
    [get_val(x, v, st, t; kwargs...) for x in o.measures]
get_val(o::Union{ObsExp1, ObsExp2}, v::Dict, ::State, ::Number; kwargs...) = o.name => v[o.obs]
get_val(o::ObsLit, v::Dict, ::State, ::Number; kwargs...) = o.name => sum(v[p] for p in o.obs)
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
        string(0) => []
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

measure(state::State, args, t::Number = 0.; kwargs...) =
    measure(state, Measure(args), t; kwargs...)

measure(state::State, m::Measure, t::Number = 0.; kwargs...) =
    measure(state, [m], t; kwargs...)[1]

function measure(state::State, m::Vector{Measure}, t::Number = 0.; kwargs...)
    st = copy(state)
    prep = PreObs(st)
    vals = Dict()
    prods = collect(Set(get_prods(m)))
    if !isempty(prods)
        push!(vals, (prods .=> expect!(st, prods, prep))...)
    end
    exp1s = collect(Set(get_exp1(m)))
    if !isempty(exp1s)
        push!(vals, (exp1s .=> expect1!(st, exp1s, prep))...)
    end
    exp2s = collect(Set(get_exp2(m)))
    if !isempty(exp2s)
        push!(vals, (exp2s .=> expect2!(st, exp2s, prep))...)
    end
    return get_val(m, vals, st, t; kwargs...)
end