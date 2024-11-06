@mixer(
    # name of base type
    "Qubit",
    # whether it takes an optional dim parameter (if yes its default value)
    false,
    # states to be converted (vector or tuple)
    [
        "0", "1", "+", "-", "i", "-i", "Up", "Dn", "↑", "↓", 
        "X+", "X-", "Y+", "Y-", "Z+", "Z-",
       "Tetra1", "Tetra2", "Tetra3", "Tetra4",
    ],
    # observables / operators to be converted (vector or tuple), names must be valid julia identifiers
    # change a name using "newname" => "basename" (newname must not be an existing name of the base type)
    [
        "I", "Id",
        "X", "Y", "Z", "Sp", "Sm",
        "H", "P", "S", "T", "SqrtNOT" => "√NOT", "Rx", "Ry", "Rz",
        "CNOT", "CX", "CY", "CZ", "Swap", "CPhase", "CRx", "CRy", "CRz", "SqrtSwap" => "√SWAP",
        "Rxx", "Rxy", "Ryy", "Rzz",
        "CCNOT", "CSwap",
        "CCCNOT" 
    ],
    # dissipators to be created (vector or tuple) "dissipator_name" => "base operator"
    # dissipator_name must be a valid julia identifier
    [
        "DUp" => "S+", 
        "DDn" => "S-",
        "DPhase" => "Z",
        "DX" => "X",
        "DY" => "Y",
        "DZ" => "Z",
    ],
    # list of operators having has_fermion_string set to true
    []
)

ITensors.state(::StateName"Half", ::SiteType"MixedQubit") = [ 0.5 0 0 0.5 ]

export graph_state

function graph_state(tp::Union{Type{Pure}, Type{Mixed}}, g::Vector{Tuple{Int, Int}}; output::Output = no_output, limits::Limits = no_limits)
    n = graph_base_size(g)
    (
        CreateState(
            name = "creating initial state + with $n qubits for graph state",
            type = tp,
            system = fill(site("Qubit"), n),
            state = "+",
        ),
        Gates(
            name = "building graph state",
            gates = prod(CZ(i, j) for (i, j) in g),
            output = output,
            limits = limits,
        )
    )
end
