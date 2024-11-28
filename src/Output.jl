export Func, Check, Output, Trace, Trace2, Norm, EE, measure

struct Func
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

wrap_check(o::Number) = Func("$o", o)
wrap_check(o::Function) = Func("func", o)
wrap_check(o) = o

struct Check
    name
    obs1
    obs2
    tol
    Check(name, o1, o2, tol=nothing) =
        new(name, wrap_check(o1), wrap_check(o2), tol)
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


struct Output
    outputs::Dict
    Output() = new(Dict())
end

Output(args...) = Output!(Output(), args...)

function Output!(output::Output, outs::Vector{<:Pair})
    foreach(outs) do out
        Output!(output, out)
    end
    return output
end

function Output!(output::Output, out::Pair{<:Any, <:Vector})
    file = first(out)
    foreach(last(out)) do o
        Output!(output, file => o)
    end
    return output
end
    
function Output!(output::Output, out::Pair)
    file_ops = get!(output.outputs, first(out), [])
    push!(file_ops, make_obs(last(out)))
    return output
end

get_prods(o::Output) = vcat(get_prods.(values(o.outputs))...)
get_prods(o::Vector) = vcat(map(get_prods, o)...)
get_prods(o::ObsLit) = o.obs
get_prods(o::Check) = [get_prods(o.obs1); get_prods(o.obs2)]
get_prods(_) = []

get_exp1(o::Output) = vcat(get_exp1.(values(o.outputs))...)
get_exp1(o::Vector) = vcat(map(get_exp1, o)...)
get_exp1(o::ObsExp1) = [o.obs]
get_exp1(o::Check) = [get_exp1(o.obs1); get_exp1(o.obs2)]
get_exp1(_) = []

get_exp2(o::Output) = vcat(get_exp2.(values(o.outputs))...)
get_exp2(o::Vector) = vcat(map(get_exp2, o)...)
get_exp2(o::ObsExp2) = [o.obs]
get_exp2(o::Check) = [get_exp2(o.obs1); get_exp2(o.obs2)]
get_exp2(_) = []

Trace = Func("Trace", trace)
Trace2 = Func("Trace2", trace2)
Norm = Func("Norm", norm)
EE(pos, spectre=0) = Func("EE($pos,$spectre)",
    st-> begin
        ee, sp = entanglement_entropy!(st, pos)
        return [[ee]; sp[1:spectre]]
    end)

get_val(o::Output, v::Dict, t::Number) =
    Dict(k => get_val(x, v, t) for (k, x) in o.outputs)
get_val(o, v::Dict,  t::Number) =
    map(o) do out
        get_val(out, v, t)
    end
get_val(o::Union{ObsExp1, ObsExp2}, v::Dict, ::Number) = o.name => v[o.obs]
get_val(o::ObsLit, v::Dict, ::Number) = o.name => sum(v[p] for p in o.obs)
get_val(o::Func, ::Dict, t::Number) =
    if o.obs isa Function
        o.name => o.obs(t)
    else
        o.name => o.obs
    end

function get_val(o::Check, v::Dict, t::Number)
    v1 = last(get_val(o.obs1, v, t))
    v2 = last(get_val(o.obs2, v, t))
    d = abs(v1 - v2)
    if !isnothing(o.tol) && d > o.tol
        error("Check $(o.name) failed with values $v1, $v2 and difference $d")
    end
    return o.name => [v1, v2, d]
end

function measure(state::State, args::Vector)
    r = measure(state, Output("" => args))
    return last.(r[""])
end

function measure(state::State, arg)
    r = measure(state, Output("" => arg))
    return last(r[""][1])
end

function measure(state::State, output::Output)
    st = copy(state)
    prep = PreObs(st)
    vals = Dict()
    prods = collect(Set(get_prods(output)))
    if !isempty(prods)
        push!(vals, (prods .=> expect!(st, prods, prep))...)
    end
    exp1s = collect(Set(get_exp1(output)))
    if !isempty(exp1s)
        push!(vals, (exp1s .=> expect1!(st, exp1s, prep))...)
    end
    exp2s = collect(Set(get_exp2(output)))
    if !isempty(exp2s)
        push!(vals, (exp2s .=> expect2!(st, exp2s, prep))...)
    end
    return get_val(output, vals, state.time)
end