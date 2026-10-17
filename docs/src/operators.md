# Operators

```@contents
Pages = ["operators.md"]
Depth = 2
```

## Usage

The examples of this page use qubits on a six site system, that is

```@setup operators
using TensorMixedStates
using .Qubits, .Circuits
n = 6
j = 1.0
h = 0.5
```

```julia
using TensorMixedStates, .Qubits, .Circuits

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

They define quantum gates, the module `Circuits` giving `H`, the rotations and the controlled
gates, like

```@example operators
gates = H(1)Swap(1, 2)H(1)
```

Noisy gates can be defined using the `Gate` constructor, for example

```@example operators
noisygate = 0.7Gate(Id) + 0.1Gate(X) + 0.1Gate(Y) + 0.1Gate(Z)
```

The reset of sites of any type to a state is the channel [`reset_gate`](@ref), and
its Lindblad generator, for an evolver, [`reset_dissipator`](@ref); to the fully mixed
state, it is the depolarization, [`depolarizing_gate`](@ref) and
[`depolarizing_dissipator`](@ref). The dephasing channel and its
generator are [`dephasing_gate`](@ref) and [`dephasing_dissipator`](@ref), built on
[`Dephase`](@ref), and the decay of the excitations an operator destroys,
[`amplitude_damping_gate`](@ref) and [`amplitude_damping_dissipator`](@ref).

Two conventions differ from one text to another. `depolarizing_gate(p)` replaces the state of a
site by the fully mixed one with probability `p`, ``(1 - p)\rho + p\,I/d``: on a qubit, the
channel ``(1 - p)\rho + \frac{p}{3}(X \rho X + Y \rho Y + Z \rho Z)``, which applies each
Pauli operator with probability ``p/3``, is `depolarizing_gate(4p/3)`. On a qubit, `"Up"` is
``|0\rangle`` and `"Dn"` is ``|1\rangle``, the excitation that `N` counts and that `Sp`
destroys: `amplitude_damping_gate` and `amplitude_damping_dissipator` take `"Dn"` to `"Up"`.

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
⊗(::TensorMixedStates.GenericOp{K, N}, ::TensorMixedStates.GenericOp{K, M}) where {K, N, M}
Proj
Level
Swap
controlled
Dissipator
Gate
Reset
Dephase
reset_gate
reset_dissipator
depolarizing_gate
depolarizing_dissipator
dephasing_gate
dephasing_dissipator
amplitude_damping_gate
amplitude_damping_dissipator
Evolver
named
parity
mod(::TensorMixedStates.GenericOp{TensorMixedStates.IsOperator}, ::Int)
dag(::TensorMixedStates.GenericOp{TensorMixedStates.IsOperator})
exp(::TensorMixedStates.GenericOp{TensorMixedStates.IsOperator})
sqrt(::TensorMixedStates.GenericOp)
matrix
tensor
```

## Gates of circuits

The module `Circuits`, `using .Circuits`, gathers the gates of quantum circuits: `H`, `S` and
`T`, which `H` and `S` extend to qudits, `Phase`, `Sum`, the rotations and the controlled gates.
A rotation by `θ` has `θ/2` in its exponent, `Rx(θ)` being ``e^{-i\theta X/2}`` on a qubit and
`Rxx(θ)` ``e^{-i\theta X \otimes X/2}`` on two, as in OpenQASM, Qiskit or PennyLane.

```@docs
TensorMixedStates.Circuits
Phase
Sum
Rx
Ry
Rz
Rxx
Ryy
Rzz
CX
CNOT
CY
CZ
CCX
Toffoli
```

## Operators on another system

A representation of one's own that keeps its tensors on another system, see
[States of one's own](@ref), places the operators it measures or applies there with
`map_sites`.

```@docs
map_sites
```

## Operator types

These types describe the operators themselves, they are mostly useful when writing functions
operating on operators. They are not exported: `TensorMixedStates.GenericOp`.

```@docs
TensorMixedStates.Op
TensorMixedStates.OpKind
TensorMixedStates.GenericOp
TensorMixedStates.IndexedOp
TensorMixedStates.SimpleOp
```
