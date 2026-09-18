# Operators

## Usage

There are two kinds of operators: generic (like `X`) and indexed (like `X(3)`). Indexed operators are applied to specific site numbers.

- Operators can be used to define Hamiltonians, for example

    hamiltonian = - j * sum(X(i)X(i+1) + Y(i)Y(i+1) for i in 1:n-1) - h * sum(Z(i) for i in 1:n)

or Lindbladian dissipators like

    dissipators = sum(Dissipator(Sp)(i) for i in 1:n)

to build Lindbladian

    lindbladian = -im * hamiltonian + dissipators

Note the factor `-im` for the Hamiltonian.

- Operators can be used to define quantum gates like

    gates = H(1)Swap(1, 2)H(1)

Noisy gates can be defined using the `Gate` constructor, for example

    noisygate = 0.7Gate(Id) + 0.1Gate(X) + 0.1Gate(Y) + 0.1Gate(Z)

- Operators can be used to define observables

    obs = X(1)X(2)Z(3)

## Reference

Complex operators can be built from a rich set of functions, for example

    Rxy(t) = exp(-im * t * (X⊗X + Y⊗Y) / 4)

Operators can be added and multiplied using usual operators (`+`, `-`, `*`, `/`, `^`).

```@docs
Operator
AtIndex
⊗
Proj
Dissipator
Gate
SetState
Left
Right
Evolver
Identity
JW
JW_F
Multi_F
dag(::GenericOp{Pure})
isfermionic
has_fermionic
matrix
tensor
simplify
removeMulti
```

## Operator types

These types describe the operators themselves, they are mostly useful when writing functions
operating on operators.

```@docs
PM
Pure
Mixed
GI
Generic
Indexed
Op
GenericOp
IndexedOp
SimpleOp
OpType
plain_op
fermionic_op
selfadjoint_op
involution_op
```
