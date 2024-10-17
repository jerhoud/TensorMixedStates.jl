export site

function combinerto(i1::Index, i2::Index, i3::Index)
    c = combiner(i1, i2)
    i = combinedind(c)
    replaceind(c, i, i3)
end

recombine(r, i, j, k) =
    *(r, combinerto.(k, j, i)..., combinerto.(k', j', i')...)

function make_indices(withdim::Bool, tp::String, i::Index...)
    if withdim
        d = dim(i[1]) / 2
        j = [siteind(tp; dim = d) for _ in i]
        k = [siteind(tp; dim = d) for _ in i]
    else
        j = [siteind(tp) for _ in i]
        k = [siteind(tp) for _ in i]
    end
    return (j, k)
end

function make_state(withdim::Bool, tp::String, st::String, i::Index; kwargs...)
    if withdim
        j = siteind(tp; dim = dim(i) / 2)
    else
        j = siteind(tp)
    end
    s = state(j, st; kwargs...)
    return s * dag(s') * combinerto(j', j, i)
end

function make_observable(withdim::Bool, tp::String, name::String, i::Index...; kwargs...)
    j, k = make_indices(withdim, tp, i...)
    o = op(name, j...; kwargs...)
    return recombine(*(o, delta.(k, k')...), i, j, k)
end

function make_operator(withdim::Bool, tp::String, name::String, i::Index...; kwargs...)
    j, k = make_indices(withdim, tp, i...)
    oj = op(name, j...; kwargs...)
    ok = op(name, k...; kwargs...)
    r = *(oj, delta.(k, k')...) - *(dag(ok), delta.(j, j')...)
    return recombine(r, i, j, k)
end

function make_gate(withdim::Bool, tp::String, name::String, i::Index...; kwargs...)
    j, k = make_indices(withdim, tp, i...)
    oj = op(name, j...; kwargs...)
    ok = op(name, k...; kwargs...)
    return recombine(oj * dag(ok), i, j, k)
end

function make_dissipator(withdim::Bool, tp::String, name::String, i::Index...; kwargs...)
    j, k = make_indices(withdim, tp, i...)
    oj = op(name, j...; kwargs...)
    aoj = swapprime(dag(oj'), 1=>2)
    ok = op(name, k...; kwargs...)
    aok = swapprime(dag(ok'), 1=>2)
    r = oj * dag(ok) -
        0.5 * *(replaceprime(aoj * oj, 2 => 1), delta.(k, k')...) -
        0.5 * *(replaceprime(aok * ok, 2 => 1), delta.(j, j')...)
    return recombine(r, i, j, k)
end

site(type::String; kwargs...) =
    tp ->
    if tp == Pure
        siteind(type; kwargs...)
    else
        siteind("Mixed" * type; kwargs...)
    end

macro mixer(
    type_exp,
    withdim_exp,
    states_exp,
    obs_exp,
    dissip_exp
)
    type = eval(type_exp)
    withdim = eval(withdim_exp)
    states = eval(states_exp)
    obs = eval(obs_exp)
    dissip = eval(dissip_exp)

    mixed = "Mixed" * type
    stype = Symbol(type)

    e = Expr(:block)

    if withdim ≠ 0
        push!(e.args,
        quote
            ITensors.space(::(@SiteType_str($mixed)); dim = $withdim) = 2 * space((@SiteType_str($type))(); dim)
        end)
    else
        push!(e.args,
        quote
            ITensors.space(::(@SiteType_str($mixed))) = 2 * space((@SiteType_str($type))())
        end)
    end

    parse

    withdim = (withdim ≠ false)

    if withdim
        quote
            ITensors.state(::StateName{N}, ::(@SiteType_str($mixed)), i::Index; kwargs...) where {N} =
            make_state($withdim, $type, string(N), i; kwargs...)
        end
    end

    foreach(states) do st
        push!(e.args,
        quote
            ITensors.state(::(@StateName_str($st)), ::(@SiteType_str($mixed)), i::Index; kwargs...) =
            make_state($withdim, $type, $st, i; kwargs...)
        end)
    end

    push!(e.args,
    quote
        function ITensors.op(::OpName"obs", ::(@SiteType_str($mixed)), i::Index)
            if $withdim
                j = siteind($type; dim = dim(i) / 2)
            else
                j = siteind($type)
            end
            return combinerto(j', j, i) * dense(delta(j, j'))
        end
    end)

    foreach(obs) do o
        if o isa String
            so = bso = o
        else
            so = o[1]
            bso = o[2]
            push!(e.args,
            quote
                ITensors.op(::(@StateName_str($so)), ::(@SiteType_str($type)), i::Index; kwargs...) =
                op((@StateName_str($bso))(), (@SiteType_str($type))(), i; kwargs...)
            end)
        end
        sobs = "obs" * so
        soper = "oper" * so
        sgate = "gate" * so
        push!(e.args,
        quote
            ITensors.op(::(@OpName_str($sobs)), ::(@SiteType_str($mixed)), i::Index...; kwargs...) =
            make_observable($withdim, $type, $bso, i...; kwargs...)
            ITensors.op(::(@OpName_str($soper)), ::(@SiteType_str($mixed)), i::Index...; kwargs...) =
            make_operator($withdim, $type, $bso, i...; kwargs...)
            ITensors.op(::(@OpName_str($sgate)), ::(@SiteType_str($mixed)), i::Index...; kwargs...) =
            make_gate($withdim, $type, $bso, i...; kwargs...)
            @opLit($so)
        end)
    end

    foreach(dissip) do o
        so = o[1]
        bso = o[2]
        push!(e.args,
        quote
            ITensors.op(::(@OpName_str($so)), ::(@SiteType_str($mixed)), i::Index...; kwargs...) =
            make_dissipator($withdim, $type, $bso, i...; kwargs...)
            @opLit($so)
        end)
    end

    return e
end