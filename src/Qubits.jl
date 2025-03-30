export Qubits

struct Qubit <: AbstractSite end

dim(::Qubit) = 2

controlled(op::ExprOp{Pure, N}) where N =
    TensorOp{Pure, N+1}(ExprOp[[Proj("Up")] ; [Id for _ in 1:N]]) + (Proj("Dn") ⊗ op)

@def_states(Qubit(),
[
    ["Up", "Z+"] => [1., 0.],
    ["Dn", "Z-"] => [0., 1.],
    ["+", "X+"] => [1., 1.] / √2,
    ["-", "X-"] => [1., -1.] / √2,
    ["i", "Y+"] => [1., im] / √2,
    ["-i", "Y-"] => [1., -im] / √2,
])


@def_operators(Qubit(),
[
    F = Id,
    X = [0. 1. ; 1. 0.],
    Y = [0. -im; im 0.],
    Z = [1. 0. ; 0. -1.],
    Sp = [0. 1. ; 0. 0.],
    Sm = [0. 0. ; 1. 0.],
    Sx = X / 2,
    Sy = Y / 2,
    Sz = Z / 2,
    S2 = 0.75 * Id,
    H = [1. 1. ; 1. -1] / √2,
    S = [1. 0. ; 0. im],
    Swap = (Id ⊗ Id + X ⊗ X + Y ⊗ Y + Z ⊗ Z) / 2
])

module Qubits

import ..Qubit, ..controlled, ..X, ..Y, ..Z, ..Sx, ..Sy, ..Sz, ..S2, ..Sp, ..Sm, ..H, ..S, ..Swap
export Qubit, controlled, X, Y, Z, Sx, Sy, Sz, S2, Sp, Sm, H, S, Swap

end