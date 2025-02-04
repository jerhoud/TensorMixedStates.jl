export @add_operators, @add_fermionic_operators, @add_dissipators, @add_fermionic_dissipators
export make_operator, MixObservable, MixGate, MixEvolve, MixEvolve2, MixDissipator, MixDissipatorF, create_mixed_gate

struct MixObservable end
struct MixGate end
struct MixEvolve end
struct MixEvolve2 end
struct MixDissipator end
struct MixDissipatorF end


function make_operator(::TPure, ::System, t::ITensor, ::Int...)
    return t
end
    
function make_operator(::Type{MixObservable}, system::System, t::ITensor, i::Int)
    idx = system.pure_sites[i]
    jdx = sim(idx)
    kdx = system.mixed_sites[i]
    return t * delta(jdx, jdx') * combinerto(jdx, idx, kdx) * combinerto(jdx', idx', kdx')
end

function make_operator(::Type{MixGate}, system::System, t::ITensor, i::Int...)
    idx = map(k->system.pure_sites[k], i)
    jdx = sim(idx)
    kdx = map(k->system.mixed_sites[k], i)
    ti = t
    tj = replaceinds(ti, (idx..., idx'...), (jdx..., jdx'...))
    return *(ti, dag(tj), combinerto.(jdx, idx, kdx)..., combinerto.(jdx', idx', kdx')...)
end

function make_operator(::Type{MixEvolve}, system::System, t::ITensor, i::Int)
    idx = system.pure_sites[i]
    jdx = sim(idx)
    kdx = system.mixed_sites[i]
    return t * delta(jdx, jdx') * combinerto(jdx, idx, kdx) * combinerto(jdx', idx', kdx')
end

function make_operator(::Type{MixEvolve2}, system::System, t::ITensor, i::Int)
    idx = system.pure_sites[i]
    jdx = sim(idx)
    kdx = system.mixed_sites[i]
    return dag(t) * delta(jdx, jdx') * combinerto(idx, jdx, kdx) * combinerto(idx', jdx', kdx')
end

function make_operator(::Type{MixDissipator}, system::System, t::ITensor, i::Int)
    idx = system.pure_sites[i]
    jdx = sim(idx)
    kdx = system.mixed_sites[i]
    ti = t
    tj = replaceinds(ti, (idx, idx'), (jdx, jdx'))
    ati = swapprime(dag(ti'), 1=>2)
    atj = swapprime(dag(tj'), 1=>2)
    r = ti * dag(tj) -
        0.5 * replaceprime(ati * ti, 2 => 1) * delta.(jdx, jdx') -
        0.5 * replaceprime(atj * tj, 2 => 1) * delta.(idx, idx')
    return r * combinerto(jdx, idx, kdx) * combinerto(jdx', idx', kdx')
end

function make_operator(::Type{MixDissipatorF}, system::System, t::ITensor, i::Int)
    idx = system.pure_sites[i]
    jdx = sim(idx)
    kdx = system.mixed_sites[i]
    ti = t
    tj = replaceinds(ti, (idx, idx'), (jdx, jdx'))
    ati = swapprime(dag(ti'), 1=>2)
    atj = swapprime(dag(tj'), 1=>2)
    r = - 0.5 * replaceprime(ati * ti, 2 => 1) * delta.(jdx, jdx') -
        0.5 * replaceprime(atj * tj, 2 => 1) * delta.(idx, idx')
    return r * combinerto(jdx, idx, kdx) * combinerto(jdx', idx', kdx')
end

(oper::Operator)(sitename::String; tp = Pure, kwargs...) =
    oper(Site(sitename); tp, kwargs...)

function (oper::Operator)(site::Site; tp = Pure, kwargs...)
    if oper.dissipator
        tp = MixDissipator
    end
    d = oper.dim
    s = System(d, site)
    o = make_operator(tp, s, oper(s.pure_sites...; kwargs...), (1:d)...)
    idx = filter(i-> hasplev(i, 0), inds(o))
    if d == 1
        return Array(o, idx[1]', idx[1])
    else
        c = combiner(idx...)
        j = combinedind(c)
        return Array(c' * o * c, j', j)
    end
end

create_mixed_gate(name::String, i::Int, ops::Vector{Operator}, weights::Vector{Float64}; kwargs...) =
    Lit(
        Operator(name, name, 1, false, false),
        system -> begin
            t = ITensor()
            for (op, w) in zip(ops, weights)
                t += w * make_operator(MixGate, system, op(system.pure_sites[i]), i)
            end
            return t
        end,
        (i,),
        NamedTuple(kwargs)
    )

"""
    @add_operators(names::Vector)
    @add_fermionic_operators(names::Vector)

Defines the given names as operators

# Examples
    @add_operators(["Foo"])    #defines operator Foo
    Foo.name          # return "Foo"
    Foo.dim           # return number of indices
    Foo(1, 3, 7)      # represent operator Foo applied on sites 1, 3, 7
    Foo(i1, i2, i3)   # return ITensor of operator Foo applied on ITensor indices i1, i2
    Foo("Qubit")      # return matrix of operator Foo applied on Qubit
"""
macro add_operators(names, dim::Int = 1)
    quote
        @make_operators($names, $dim, false, false)
    end
end,

macro add_fermionic_operators(names)
    quote
        @make_operators($names, 1, true, false)
    end
end

"""
    @add_dissipators(names::Vector{Pairs})
    @add_fermionic_dissipators(names::Vector{Pairs})

Defines the given names as dissipators

# Examples
    @add_operators(["DFoo" => "Foo"])    #defines dissipator DFoo based on ITensor operator "Foo"
    DFoo.opname  # return base operator name (here "Foo")
    DFoo(1)      # represent operator Foo applied on sites 1
    DFoo(idx)    # return ITensor of operator DFoo applied on ITensor indices idx
"""
macro add_dissipators(names)
    quote
        @make_operators($names, 1, false, true)
    end
end

macro add_fermionic_dissipators(names)
    quote
        @make_operators($names, 1, true, true)
    end
end

@add_operators(
    [
    "F", "I", "Id",
    "X", "Y", "Z", "Sp" => "S+", "Sm" => "S-", "ProjUp", "ProjDn",
    "H", "P", "S", "T", "SqrtNOT" => "√NOT", "Rx", "Ry", "Rz",
    "A", "Adag", "N",
    "Aup", "Adagup", "Adn", "Adagdn",
    "Fup", "Fdn", "Nup", "Ndn", "Nupdn", "Ntot",
    "Sx", "Sy", "Sz", "Sx2", "Sy2", "Sz2", "S2",
    ]
)

@add_operators(
    [
    "CNOT", "CX", "CY", "CZ", "Swap", "CPhase", "CRx", "CRy", "CRz", "SqrtSwap" => "√SWAP",
    "Rxx", "Rxy", "Ryy", "Rzz",
    ],
    2
)

@add_operators([ "CCNOT", "CSwap"], 3)
@add_operators([ "CCCNOT" ], 4)

@add_fermionic_operators(
    [
    "C", "Cdag",
    "Cup", "Cdagup", "Cdn", "Cdagdn",
    ]
)

@add_dissipators(
    [
    "DUp" => "S+",
    "DDn" => "S-",
    "DX" => "X",
    "DY" => "Y",
    "DZ" => "Z",
    "DPhase" => "Z",
    "DN" => "N",
    "DA" => "A",
    "DAdag" => "Adag",
    "DAup" => "Aup",
    "DAdagup" => "Adagup",
    "DAdn" => "Adn",
    "DAdagdn" => "Adagdn",
    "DNup" => "Nup",
    "DNdn" => "Ndn",
    "DNtot" => "Ntot",
    ]
)

@add_fermionic_dissipators(
    [
    "DC" => "C",
    "DCdag" => "Cdag",
    "DCup" => "Cup",
    "DCdagup" => "Cdagup",
    "DCdn" => "Cdn",
    "DCdagdn" => "Cdagdn",
    ]
)
