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
    measures::Dict
    Measure() = new(Dict())
end

Measure(args...) = Measure!(Measure(), args...)

function Measure!(output::Measure, outs::Vector{<:Pair})
    foreach(outs) do out
        Measure!(output, out)
    end
    return output
end

function Measure!(output::Measure, out::Pair{<:Any, <:Vector})
    file = first(out)
    foreach(last(out)) do o
        Measure!(output, file => o)
    end
    return output
end
    
function Measure!(output::Measure, out::Pair)
    file_ops = get!(output.measures, first(out), [])
    push!(file_ops, make_obs(last(out)))
    return output
end

get_prods(o::Measure) = vcat(get_prods.(values(o.measures))...)
get_prods(o::Vector) = vcat(map(get_prods, o)...)
get_prods(o::ObsLit) = o.obs
get_prods(o::Check) = [get_prods(o.obs1); get_prods(o.obs2)]
get_prods(_) = []

get_exp1(o::Measure) = vcat(get_exp1.(values(o.measures))...)
get_exp1(o::Vector) = vcat(map(get_exp1, o)...)
get_exp1(o::ObsExp1) = [o.obs]
get_exp1(o::Check) = [get_exp1(o.obs1); get_exp1(o.obs2)]
get_exp1(_) = []

get_exp2(o::Measure) = vcat(get_exp2.(values(o.measures))...)
get_exp2(o::Vector) = vcat(map(get_exp2, o)...)
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
        return [[ee]; sp[1:spectre]]
    end)
Linkdim = StateFunc("Linkdim", maxlinkdim)
MemoryUsage = StateFunc("MemoryUsage", Base.summarysize)

get_val(o::Measure, v::Dict, st::State, t::Number) =
    Dict(k => get_val(x, v, st, t) for (k, x) in o.measures)
get_val(o, v::Dict, st::State,  t::Number) =
    map(o) do out
        get_val(out, v, st, t)
    end
get_val(o::Union{ObsExp1, ObsExp2}, v::Dict, ::State, ::Number) = o.name => v[o.obs]
get_val(o::ObsLit, v::Dict, ::State, ::Number) = o.name => sum(v[p] for p in o.obs)
get_val(o::TimeFunc, ::Dict, ::State, t::Number) =
    if o.obs isa Function
        o.name => o.obs(t)
    else
        o.name => o.obs
    end
get_val(o::StateFunc, ::Dict, st::State, ::Number) = o.name => o.obs(st)

function get_val(o::Check, v::Dict, st::State, t::Number)
    v1 = last(get_val(o.obs1, v, st, t))
    v2 = last(get_val(o.obs2, v, st, t))
    d = norm(v1 - v2)
    if !isnothing(o.tol) && d > o.tol
        error("Check $(o.name) failed with values $v1, $v2 and difference $d")
    end
    return o.name => [v1, v2, d]
end

function measure(state::State, args::Vector, t::Number = 0.)
    r = measure(state, Measure("" => args), t)
    return last.(r[""])
end

function measure(state::State, arg, t::Number = 0)
    r = measure(state, Measure("" => arg), t)
    return last(r[""][1])
end

function measure(state::State, m::Measure, t::Number = 0.)
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
    return get_val(m, vals, st, t)
end