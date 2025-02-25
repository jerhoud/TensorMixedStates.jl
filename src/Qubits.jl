export Qubits

struct Qubit <: AbstractSite end

site_dim(::Qubit) = 2
site_name(::Qubit) = "Qubit"

controlled(op::ExprOp{N}) where N =
    TensorOp{N+1}(ExprOp[ [PUp] ; [I for _ in 1:N]]) + (PDn ⊗ op)

@def_operators(Qubit(),
[
    (I = [1. 0. ; 0. 1.], "Identity operator"),
    (F = [1. 0. ; 0. 1.], "Jordan Wigner F, identity in this case"),
    (X = [0. 1. ; 1. 0.], "Pauli operator X"),
    (Y = [0. -im; im 0.], "Pauli operator Y"),
    (Z = [1. 0. ; 0. -1.], "Pauli operator Z"),
    (PUp = [1. 0. ; 0. 0.], "Projector on up"),
    (PDn = [0. 0. ; 0. 1.], "Projector on down"),
    (Sp = [0. 1. ; 0. 0.], "S+"),
    (Sm = [0. 0. ; 1. 0.], "S-"),
    (NOT = X, "NOT gate, same as X"),
    (CX = controlled(X), "Controlled X gate, same as CNOT"),
    (CNOT = CX, "Controlled NOT gate, same as CX"),
    (CY = controlled(Y), "Controlled Y gate"),
    (CZ = controlled(Z), "Controlled Z gate"),
])
    
module Qubits

import ..Qubit, ..controlled, ..I, ..F, ..X, ..Y, ..Z, ..PUp, ..PDn, ..Sp, ..Sm, ..CX, ..CNOT, ..CY, ..CZ
export Qubit, controlled, I, F, X, Y, Z, PUp, PDn, Sp, Sm, CX, CNOT, CY, CZ

end