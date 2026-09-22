# Operators

```@contents
Pages = ["operators.md"]
Depth = 2
```

## Usage

The examples of this page use qubits on a six site system, that is

```@setup operators
using TensorMixedStates
using .Qubits
n = 6
j = 1.0
h = 0.5
```

```julia
using TensorMixedStates, .Qubits

n = 6
j = 1.0
h = 0.5
```

There are two kinds of operators: generic (like `X`) and indexed (like `X(3)`). Indexed operators are applied to specific site numbers.

- Operators can be used to define Hamiltonians, for example

```@example operators
hamiltonian = - j * sum(X(i)X(i+1) + Y(i)Y(i+1) for i in 1:n-1) - h * sum(Z(i) for i in 1:n)
```

or Lindbladian dissipators like

```@example operators
dissipators = sum(Dissipator(Sp)(i) for i in 1:n)
```

to build Lindbladian

```@example operators
lindbladian = -im * hamiltonian + dissipators
```

Note the factor `-im` for the Hamiltonian.

- Operators can be used to define quantum gates like

```@example operators
gates = H(1)Swap(1, 2)H(1)
```

Noisy gates can be defined using the `Gate` constructor, for example

```@example operators
noisygate = 0.7Gate(Id) + 0.1Gate(X) + 0.1Gate(Y) + 0.1Gate(Z)
```

- Operators can be used to define observables

```@example operators
obs = X(1)X(2)Z(3)
```

## Reference

Complex operators can be built from a rich set of functions, for example

```@example operators
Rxy(t) = exp(-im * t * (X⊗X + Y⊗Y) / 4)
```

```@example operators
Rxy(0.2)(2, 5)
```

Operators can be added and multiplied using usual operators (`+`, `-`, `*`, `/`, `^`).

```@docs
Operator
⊗
Proj
Dissipator
Gate
SetState
Left
Right
Evolver
named
parity
mod(::GenericOp{Pure}, ::Integer)
dag(::GenericOp{Pure})
isfermionic
has_fermionic
flux(::SimpleOp, ::AbstractSite)
matrix
tensor
simplify
```

## Operator types

These types describe the operators themselves, they are mostly useful when writing functions
operating on operators.

```@docs
Pure
Mixed
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
