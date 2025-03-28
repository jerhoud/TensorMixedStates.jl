export Qubits

struct Qubit <: AbstractSite end

dim(::Qubit) = 2

controlled(op::ExprOp{Pure, N}) where N =
    TensorOp{Pure, N+1}(ExprOp[[PUp] ; [Id for _ in 1:N]]) + (PDn ⊗ op)

@def_states(Qubit(),
[
    ["Up", "0", "Z+"] => [1., 0.],
    ["Dn", "1", "Z-"] => [0., 1.],
    ["+", "X+"] => [1., 1.] / √2,
    ["-", "X-"] => [1., -1.] / √2,
    ["i", "Y+"] => [1., im] / √2,
    ["-i", "Y-"] => [1., -im] / √2,
])


@def_operators(Qubit(),
[
    Id = [1. 0. ; 0. 1.],
    F = [1. 0. ; 0. 1.],
    X = [0. 1. ; 1. 0.],
    Y = [0. -im; im 0.],
    Z = [1. 0. ; 0. -1.],
    PUp = [1. 0. ; 0. 0.],
    PDn = [0. 0. ; 0. 1.],
    Sp = [0. 1. ; 0. 0.],
    Sm = [0. 0. ; 1. 0.],
    NOT = X,
    CX = controlled(X),
    CNOT = CX,
    CY = controlled(Y),
    CZ = controlled(Z),
])

module Qubits

import ..Qubit, ..controlled, ..Id, ..F, ..X, ..Y, ..Z, ..PUp, ..PDn, ..Sp, ..Sm, ..NOT, ..CX, ..CNOT, ..CY, ..CZ
export Qubit, controlled, Id, F, X, Y, Z, PUp, PDn, Sp, Sm, NOT, CX, CNOT, CY, CZ

end