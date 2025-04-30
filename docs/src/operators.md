# Operators

## Usage

There are two kinds of operators generic (like `X`) and indexed (like `X(3)`). Indexed operators are applied to specific site numbers.

- Operators can be used to defined Hamiltonians, for example

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

Complex operators can be build from a rich set of functions, for example

    Rxy(t) = exp(-im * t * (X⊗X + Y⊗Y) / 4)

Operators can be added and multiplied using usual operators (`+`, `-`, `*`, `/`, `^`).

```@docs
Pure
Mixed
ExprOp
Operator
Indexed
⊗
Proj
Dissipator
Gate
dag
fermionic
matrix
tensor
simplify
insertFfactors
process
```
