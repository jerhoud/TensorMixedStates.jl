# Measurements

## measure

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

## Types of measurements

There are five types of measurements we can ask to measure. The examples below use a four
qubit state

```@setup measurements
using TensorMixedStates
using .Qubits
mystate = State{Mixed}(System(4, Qubit()), "Up")
```

```julia
using TensorMixedStates, .Qubits

mystate = State{Mixed}(System(4, Qubit()), "Up")
```

With an indexed operator, we get the expectation of the corresponding observable

```@example measurements
measure(mystate, X(1)X(2)Z(4))
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

```julia
Trace                   # returns the trace of the state
Trace2                  # returns the trace of the square of the state (alternate name: Purity)
TraceError              # returns 1 - trace, useful for monitoring trace deviations
Renyi2                  # returns the Renyi entropy of order 2 of the state
SubRenyi2(sub)          # returns the Renyi entropy of order 2 of the subsystem (sub is a vector containing the indices of the sites of the subsystem)
MutualInfoRenyi2(sub)   # returns the Renyi2 mutual information of the two subsystem (you either give one subsystem as a vector of indices or a splitting link)
EntanglementEntropy(l)     # entanglement entropy for the cut on the right of site l, between l and l+1
EntanglementEntropy(l, n)  # the same, followed by the first n eigenvalues of the reduced density matrix of sites 1 to l
Hermiticity             # returns 1 if density matrix is really Hermitian and down to 0 for anti Hermitian density matrix
HermiticityError        # returns 1 - Hermiticity for monitoring hermiticity deviation
MaxLinkdim              # returns the maximum bond dimension of the representation
```

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
`GroundState` and `SteadyState`. `:energy` is the current energy and is available in
`GroundState` only. Asking for a symbol that the running algorithm does not provide is not
an error: the measurement silently produces an empty value, so `:energy` in an `Evolve`
phase writes nothing.

Checks can be performed (useful for coherence tests)

```julia
measure(mystate, Check(name, obs1, obs2))      # returns 3 values obs1, obs2 and |obs2 - obs1|
measure(mystate, Check(name, obs1, obs2, tol)) # if |obs2 - obs1|>tol throw an error
measure(mystate, Check(name, obs1, obs2), t)
```

In `Check`, obs may also be constants, vectors and function of time (like `t -> sin(t)`), in this case the simulation
time must be fed to `measure` as 3rd argument.

## Real, imaginary and complex values

Each value is given the kind its measurement calls for. An operator is real when it is self
adjoint, like `X(1)` or `Z(1)Z(2) + X(1)`, purely imaginary when its adjoint is its opposite,
like `dag(C)(1)C(2) - dag(C)(2)C(1)`, and complex otherwise, like `Sp(1)`. A correlation
matrix takes one kind for all its entries, so `(X, Y)`, whose diagonal is
``\langle XY \rangle = i \langle Z \rangle``, is complex. A state function or a function of
time is real, except `Overlap`, and a number takes the kind of its type.

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

## Output

```@docs
output
log_msg
```