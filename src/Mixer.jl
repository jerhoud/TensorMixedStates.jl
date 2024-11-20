export site

struct MixObservable end
struct MixEvolve end
struct MixEvolve2 end
struct MixGate end
struct MixDissipator end

mixed_index(i::Index) =
    addtags(Index(2dim(i), tags(i)), "Mixed")

pure_index(i::Index) =
    Index(dim(i)÷2, tags(i))

function combinerto(i1::Index, i2::Index, i3::Index)
    c = combiner(i1, i2)
    i = combinedind(c)
    replaceind(c, i, i3)
end

function make_state(::Type{Pure}, st::String, i::Index)
    return state(i, st)
end

function make_state(::Type{Mixed}, st::String, i::Index)
    j = pure_index(i)
    s = state(j, st)
    return s * dag(s') * combinerto(j', j, i)
end

function make_operator(::Type{Pure}, tj::ITensor, ::Index)
    return tj
end

function make_operator(::Type{MixObservable}, tj::ITensor, i::Index)
    j = noprime(ind(tj, 1))
    return tj * combinerto(j', j, i)
end

function obs(i::Index)
    j = pure_index(i)
    return dense(delta(j, j')) * combinerto(j', j, i)
end

function make_operator(::Type{MixEvolve}, tj::ITensor, i::Index)
    j = noprime(ind(tj, 1))
    k = sim(j)
    tk = replaceinds(tj, (j, j'), (k, k'))
    return *(tj, delta(k, k'), combinerto(k, j, i), combinerto(k', j', i ))
end

function make_operator(::Type{MixEvolve2}, tj::ITensor, i::Index)
    j = noprime(ind(tj, 1))
    k = sim(j)
    tk = replaceinds(tj, (j, j'), (k, k'))
    return *(tk, delta(j, j'), combinerto(k, j, i), combinerto(k', j', i ))
end

function make_operator(::Type{MixGate}, tj::ITensor, i::Index...)
    j = filter(t->hasplev(t, 0), inds(tj))
    k = sim(j)
    tk = replaceinds(tj, (j..., j'...), (k..., k'...))
    return *(tj, dag(tk), combinerto.(k, j, i)..., combinerto.(k', j', i')...)
end

function make_operator(::Type{MixDissipator}, tj::ITensor, i::Index)
    j = noprime(ind(tj, 1))
    k = sim(j)
    tk = replaceinds(tj, (j, j'), (k, k'))
    atj = swapprime(dag(tj'), 1=>2)
    atk = swapprime(dag(tk'), 1=>2)
    r = tj * dag(tk) -
        0.5 * *(replaceprime(atj * tj, 2 => 1), delta.(k, k')) -
        0.5 * *(replaceprime(atk * tk, 2 => 1), delta.(j, j'))
    return *(r, combinerto(k, j, i), combinerto(k', j', i'))
end

site(::Type{Pure}, type::String; kwargs...) =
    siteind(type; kwargs...)

site(::Type{Mixed}, type::String; kwargs...) =
    mixed_index(siteind(type; kwargs...))

site(type::String; kwargs...) =
    tp -> site(tp, type; kwargs...)

macro add_operators(names...)
    quote
        @opLits(false, false, $names...)
    end
end

macro add_fermionic_operators(names...)
    quote
        @opLits(true, false, $name...)
    end
end

macro add_dissipators(names...)
    quote
        @opLits(false, true, $names...)
    end
end

@add_operators(
    "X", "Y", "Z", "Sp", "Sm",
    "H", "P", "S", "T", "SqrtNOT" => "√NOT", "Rx", "Ry", "Rz",
    "CNOT", "CX", "CY", "CZ", "Swap", "CPhase", "CRx", "CRy", "CRz", "SqrtSwap" => "√SWAP",
    "Rxx", "Rxy", "Ryy", "Rzz",
    "CCNOT", "CSwap",
    "CCCNOT",
    "F", "A", "Adag", "N",
    "Aup", "Adagup", "Adn", "Adagdn",
    "Fup", "Fdn", "Nup", "Ndn", "Nupdn", "Ntot",
    "Sx", "Sy", "Sz", "Sx2", "Sy2", "Sz2", "S2",
)

@add_fermionic_operators(
    "C", "Cdag",
    "Cup", "Cdagup", "Cdn", "Cdagdn",
)

@add_dissipators(
    "DUp" => "S+",
    "DDn" => "S-",
    "DX" => "X",
    "DY" => "Y",
    "DZ" => "Z",
    "DPhase" => "Z",
)