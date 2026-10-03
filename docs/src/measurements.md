# Measurements

## Kinds of measurements

`measure` takes operators, state functions, symbols, checks, numbers or functions of the
time, and strings, which are written as a label with no value. The examples below use a four
qubit state

```@setup measurements
using TensorMixedStates
using .Qubits
mystate = State{Mixed}(System(4, Qubit()), ["+", "Up", "+", "Dn"])
```

```julia
using TensorMixedStates, .Qubits

mystate = State{Mixed}(System(4, Qubit()), ["+", "Up", "+", "Dn"])
```

With an indexed operator, we get the expectation of the corresponding observable

```@example measurements
measure(mystate, X(1)X(3)Z(4))
```

With a generic operator (acting on one site only), we get the expectation of the operator on each site

```@example measurements
measure(mystate, X)
```

With a couple of generic operators, we get the correlation matrix of those observables

```@example measurements
measure(mystate, (X, Y))
```

There are some state functions predefined:

| state function | value |
|:---|:---|
| `Trace` | the trace of the density matrix |
| `TraceError` | `1 - trace`, to monitor how far the trace drifts |
| `Trace2`, `Purity` | the trace of the square of the density matrix |
| `Norm` | the norm of the state |
| `Hermiticity` | 1 for a Hermitian density matrix, down to 0 for an anti-Hermitian one |
| `HermiticityError` | `1 - hermiticity` |
| `Renyi2` | the Rényi entropy of order 2 of the state |
| `SubRenyi2(sites)`, `SubRenyi2(link)` | the Rényi entropy of order 2 of the sites given, or of those up to the link |
| `MutualInfoRenyi2(sites)`, `MutualInfoRenyi2(link)` | the Rényi-2 mutual information between the sites given, or those up to the link, and the rest |
| `EntanglementEntropy(l)` | the entanglement entropy across the cut between sites `l` and `l+1`, the operator space entanglement entropy (OSEE) on a mixed representation |
| `EntanglementEntropy(l, n)` | the same, followed by the first `n` values of the spectrum it is computed from, see `entanglement_entropy` |
| `Fidelity(ref)` | the fidelity with the reference state `ref`, refused between two mixed representations |
| `Overlap(ref)` | the inner product with the reference state `ref` |
| `Variance(hamiltonian)` | the variance of the energy of `hamiltonian`, on a pure representation |
| `MaxLinkdim` | the maximum bond dimension of the representation |
| `MemoryUsage` | the memory the state occupies in bytes, the caches of the measurements included |

They are used like this

```@example measurements
measure(mystate, TraceError)
```

New state functions may be defined by

```@example measurements
stfunc = StateFunc("half_trace", st -> trace(st) / 2)
```

Some quantities produced by the algorithm itself, rather than computed from the state, are
asked for with a symbol

```julia
measures = "sweeps.dat" => :sweep
```

`:sweep` is the sweep number and is available in every phase that sweeps, that is `Evolve`,
`GroundState` and `SteadyState`, and in a phase of your own that passes it to `output`, as
`output(sim, measures; sweep = k)`. `:energy` is the current energy and is available in
`GroundState` and `SteadyState`, where it is the value dmrg minimises, zero at the steady
state. Asking for a symbol that the running algorithm does not provide is not an error: the
measurement produces an empty value, so `:energy` in an `Evolve` phase writes its name and
the time with no value in a file, and empty values in a json file or a `Data` object. The
symbols are given to the `measures` of a phase, not to its `final_measures`, which are taken
once it is over: the energy a `GroundState` ends with is measured as the last line of its
`measures`, or as the expectation value of its hamiltonian.

Checks can be performed (useful for coherence tests)

```julia
measure(mystate, Check("purity", Purity, 1))            # the two values and their distance
measure(mystate, Check("purity", Purity, 1, 1e-8))      # an error if the distance is above 1e-8
measure(mystate, Check("cos", X(1), t -> cos(2t)), 0.3) # against a function of time, at t = 0.3
```

The two measurements compared may be operators, state functions, constants, vectors or
functions of time, like `t -> cos(2t)`, in which case the simulation time must be given to
`measure` as its third argument. They may be symbols as well: a check on a symbol that is not
given has empty values, and is an error when it has a tolerance, having nothing to pass.

## Real, imaginary and complex values

Each value is given the kind its measurement calls for. An operator is real when it is self
adjoint, like `X(1)` or `Z(1)Z(2) + X(1)`, purely imaginary when its adjoint is its opposite,
like `dag(C)(1)C(2) - dag(C)(2)C(1)`, and complex otherwise, like `Sp(1)`. A correlation
matrix takes one kind for all its entries, so `(X, Y)`, whose diagonal is
``\langle XY \rangle = i \langle Z \rangle``, is complex. A state function or a function of
time is real, except `Overlap`, and a number takes the kind of its type. `Trace` and
`TraceError` thus give the real part of the trace, the imaginary part being a numerical
error, see [`TraceError`](@ref).

A real value is written in one column, its imaginary part, which is rounding, being dropped.
An imaginary one is written in one column too, as its imaginary part, under the name
`Im(name)`, and a complex one in two, its real part then its imaginary part. A part dropped
that is more than rounding is reported in the log.

The test giving the kind of an operator is symbolic and only says real or imaginary when it
can prove it, so what it misses comes out complex: `Sp(1)Sm(2) + Sm(1)Sp(2)`, for instance,
which is real. `RealValue`, `ImaginaryValue` and `ComplexValue` declare the kind by hand, as
they do for a state function whose values are complex

```julia
measures = "data" => [RealValue(Sp(1)Sm(2) + Sm(1)Sp(2)), ComplexValue(stfunc)]
```

## Reference

```@docs
measure
StateFunc
TimeFunc
Check
Measure
RealValue
ImaginaryValue
ComplexValue
expect
expect1
expect2
variance
sample(::State{Pure})
entanglement_entropy
entanglement_by_sector
renyi2
mutual_info_renyi2
Trace
TraceError
Trace2
Purity
Norm
Hermiticity
HermiticityError
Renyi2
SubRenyi2
Fidelity
Overlap
Variance
MutualInfoRenyi2
Mutual_Info_Renyi2
EntanglementEntropy
EE
MaxLinkdim
Linkdim
MemoryUsage
```

## Output

```@docs
output
log_msg
```
