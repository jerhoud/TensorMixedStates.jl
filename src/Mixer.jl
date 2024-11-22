import ITensors: state
export @add_operators, @add_fermionic_operators, @add_dissipators

macro add_operators(names)
    quote
        @opLit($names, false, false)
    end
end

macro add_fermionic_operators(names)
    quote
        @opLit($names, true, false)
    end
end

macro add_dissipators(names)
    quote
        @opLit($names, false, true)
    end
end

@add_operators(
    [
    "F", "I", "Id",
    "X", "Y", "Z", "Sp", "Sm",
    "H", "P", "S", "T", "SqrtNOT" => "√NOT", "Rx", "Ry", "Rz",
    "CNOT", "CX", "CY", "CZ", "Swap", "CPhase", "CRx", "CRy", "CRz", "SqrtSwap" => "√SWAP",
    "Rxx", "Rxy", "Ryy", "Rzz",
    "CCNOT", "CSwap",
    "CCCNOT",
    "A", "Adag", "N",
    "Aup", "Adagup", "Adn", "Adagdn",
    "Fup", "Fdn", "Nup", "Ndn", "Nupdn", "Ntot",
    "Sx", "Sy", "Sz", "Sx2", "Sy2", "Sz2", "S2",
    ]
)

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
    ]
)

state(::StateName"Half", ::SiteType"MixedQubit") = [0.5 0 0 0.5]