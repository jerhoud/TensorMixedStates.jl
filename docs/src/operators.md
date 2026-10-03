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

Operators define Hamiltonians, for example

```@example operators
hamiltonian = - j * sum(X(i)X(i+1) + Y(i)Y(i+1) for i in 1:n-1) - h * sum(Z(i) for i in 1:n)
```

or Lindbladian dissipators like

```@example operators
dissipators = sum(Dissipator(Sp)(i) for i in 1:n)
```

to build a Lindbladian

```@example operators
lindbladian = -im * hamiltonian + dissipators
```

Note the factor `-im` for the Hamiltonian.

They define quantum gates, like

```@example operators
gates = H(1)Swap(1, 2)H(1)
```

Noisy gates can be defined using the `Gate` constructor, for example

```@example operators
noisygate = 0.7Gate(Id) + 0.1Gate(X) + 0.1Gate(Y) + 0.1Gate(Z)
```

And they define observables

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
Operator{N}(::String, ::Union{Matrix, Function, GenericOp{Pure, N}}, ::OpType, ::AbstractSite, ::AbstractSite...) where N
⊗(::GenericOp{Pure, N}, ::GenericOp{Pure, M}) where {N, M}
Proj
Dissipator
Gate
SetState
Left
Right
Evolver
named
parity
mod(::GenericOp{Pure}, ::Int)
dag(::GenericOp{Pure})
exp(::GenericOp{Pure})
sqrt(::GenericOp)
isfermionic
has_fermionic
matrix
tensor
simplify
```

## Lindblad and Kraus forms

An evolver and gates written by the user can be read back as the terms a representation of
one's own needs to unravel them, quantum trajectories for instance, see
[Representations of one's own](@ref): the hamiltonian and the jump operators of an evolver, and
the Kraus operators of each channel of a product of gates. `map_sites` places these operators,
and those it measures, on the system its tensors lie on.

```@docs
lindblad_terms
kraus_operators
map_sites
```

## Operator types

These types describe the operators themselves, they are mostly useful when writing functions
operating on operators.

```@docs
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
