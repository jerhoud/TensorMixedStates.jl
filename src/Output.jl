export Func, Check, Output, Trace, Trace2, Norm, EE

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

struct Check
    obs1
    obs2
    tol
end

Check(o, v::Number; tol=nothing) =
    Check(o, Func("$v", v), tol)

Check(o, f::Function; tol=nothing) =
    Check(o, Func("func", f), tol)

Check(o1, o2; tol=nothing) =
    Check(o1, o2, tol)

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

Output!(output::Output, out::Vector) = 
    Output!(output, "" => out)

Output!(output::Output, outs...) = Output!(output::Output, [outs...])

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

get_prods(o::Output) = vcat(get_prods.(last.(o))...)
get_prods(o::Vector) = vcat(getmap(get_prods, o)...)
get_prods(o::ObsLit) = o.obs
get_prods(o::Check) = [get_prods(o.obs1); get_prods(o.obs2)]
get_prods(_) = []

get_exp1(o::Output) = vcat(get_exp1.(last.(o))...)
get_exp1(o::Vector) = vcat(map(get_exp1, o)...)
get_exp1(o::ObsExp1) = [o.obs]
get_exp1(o::Check) = [get_exp1(o.obs1); get_exp1(o.obs2)]
get_exp1(o) = []

get_exp2(o::Output) = vcat(get_exp2.(last.(o))...)
get_exp2(o::Vector) = vcat(map(get_exp2, o)...)
get_exp2(o::ObsExp2) = [o.obs]
get_exp2(o::Check) = [get_exp2(o.obs1); get_exp2(o.obs2)]
get_exp2(o) = []

Trace = Func("Trace", trace)
Trace2 = Func("Trace2", trace2)
Norm = Func("Norm", norm)
EE(pos, spectre=0) = Func("EE($pos,$spectre)",
    st-> begin
        ee, sp = entanglement_entropy(st, pos)
        return [[ee]; sp[1:spectre]]
    end)

function get_observables(output::Output, state::State)
    t = state.time
    st = copy(state)
    prep = PreObs(st)
    prods = collect(Set(get_prods(output)))
    exp1s = collect(Set(get_exp1(output)))
    exp2s = collect(Set(get_exp2(output)))
    prodvals = expect(st, prods, prep)
    exp1vals = expect1(st, exp1s, prep)
    exp2vals = expect2(st, exp2s, prep)
    vals = Dict()
    push!(vals, (prods .=> prodvals)...)
    push!(vals, (exp1s .=> exp1vals)...)
    push!(vals, (exp2s .=> exp2vals)...)
end