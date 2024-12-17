export @add_operators, @add_fermionic_operators, @add_dissipators

"""
    @add_operators(names::Vector)
    @add_fermionic_operators(names::Vector)

Defines the given names as operator functions

# Examples
    @add_operators(["Foo"])    #defines function Foo
    Foo()          # return "Foo"
    Foo(1, 3, 7)   # represent operator Foo applied on sites 1, 3, 7
    Foo(i1, i2)    # return ITensor of operator Foo applied on ITensor indices i1, i2
"""
macro add_operators(names)
    quote
        @opLit($names, false, false)
    end
end,

macro add_fermionic_operators(names)
    quote
        @opLit($names, true, false)
    end
end

"""
    @add_dissipators(names::Vector{Pairs})

Defines the given names as dissipator functions

# Examples
    @add_operators(["DFoo" => "Foo"])    #defines dissipator DFoo based on ITensor operator "Foo"
    DFoo()       # return "DFoo"
    DFoo(1)      # represent operator Foo applied on sites 1
    DFoo(idx)    # return ITensor of operator DFoo applied on ITensor indices idx
"""
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
