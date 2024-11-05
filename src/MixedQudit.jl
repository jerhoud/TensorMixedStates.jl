@mixer(
    # name of base type
    "Qudit",
    # whether it takes an optional dim parameter (if yes its default value)
    2,
    # states to be converted (vector or tuple)
    [],
    # observables / operators to be converted (vector or tuple), names must be valid julia identifiers
    # change a name using "newname" => "basename" (newname must not be an existing name of the base type)
    [
        "A", "Adag", "N",
    ],
    # dissipators to be created (vector or tuple) "dissipator_name" => "base operator"
    # dissipator_name must be a valid julia identifier
    [
        "DA" => "A",
        "DAdag" => "Adag",
    ]
)
