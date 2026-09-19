# Measurements

## measure

```@docs
measure
StateFunc
TimeFunc
Check
Measure
expect
expect1
expect2
sample(::State{Pure})
entanglement_entropy
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
Mutual_Info_Renyi2
EE
Linkdim
MemoryUsage
```

## Types of measurements

There are five types of measurements we can ask to measure. The examples below use a four
qubit state

```@setup measurements
using TensorMixedStates
using .Qubits
state = State{Mixed}(System(4, Qubit()), "Up")
```

```julia
using TensorMixedStates, .Qubits

state = State{Mixed}(System(4, Qubit()), "Up")
```

With an indexed operator, we get the expectation of the corresponding observable

```@example measurements
measure(state, X(1)X(2)Z(4))
```

With a generic operator (acting on one site only), we get the expectation of the operator on each site

```@example measurements
measure(state, X)
```

With a couple of generic operators, we get the correlation matrix of those observables

```@example measurements
measure(state, (X, Y))
```

There are some state functions predefined:

```julia
Trace                   # returns the trace of the state
Trace2                  # returns the trace of the square of the state (alternate name: Purity)
TraceError              # returns 1 - trace, useful for monitoring trace deviations
Renyi2                  # returns the Renyi entropy of order 2 of the state
SubRenyi2(sub)          # returns the Renyi entropy of order 2 of the subsystem (sub is a vector containing the indices of the sites of the subsystem)
Mutual_Info_Renyi2(sub) # returns the Renyi2 mutual information of the two subsystem (you either give one subsystem as a vector of indices or a splitting link)
EE(l)                   # returns entanglement entropy for the cut on the right of site l, between l and l+1
EE(l, n)                # the same, followed by the first n eigenvalues of the reduced density matrix of sites 1 to l
Hermiticity             # returns 1 if density matrix is really Hermitian and down to 0 for anti Hermitian density matrix
HermiticityError        # returns 1 - Hermiticity for monitoring hermiticity deviation
Linkdim                 # returns the maximum bond dimension of the representation
```

They are used like this

```@example measurements
measure(state, TraceError)
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
measure(state, Check(name, obs1, obs2))      # returns 3 values obs1, obs2 and |obs2 - obs1|
measure(state, Check(name, obs1, obs2, tol)) # if |obs2 - obs1|>tol throw an error
measure(state, Check(name, obs1, obs2), t)
```

In `Check`, obs may also be constants, vectors and function of time (like `t -> sin(t)`), in this case the simulation
time must be fed to `measure` as 3rd argument.

## Output

```@docs
output
log_msg
```