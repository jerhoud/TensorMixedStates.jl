# Measurements

## measure

```@docs
measure
StateFunc
Check
Measure
expect
expect1
expect2
entanglement_entropy
renyi2
mutual_info_renyi2
Trace
TraceError
Trace2
Purity
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

There are five types of measurements we can ask to measure

With an indexed operator, we get the expectation of the corresponding observable

    mesure(state, X(1)X(2)Z(4))

With a generic operator (acting on one site only), we get the expectation of the operator on each site

    measure(state, X)

With a couple of generic operators, we get the correlation matrix of those observators

    measure(state, (X, Y))

There are a some state functions predefined:

    Trace              # returns the trace of the state
    Trace2                  # returns the trace of the square of the state (alternate name: Purity)
    TraceError              # returns 1 - trace, usefull for monitoring trace deviations
    Renyi2                  # returns the Renyi entropy of order 2 of the state
    SubRenyi2(sub)          # returns the Renyi entropy of order 2 of the subsystem (sub is a vector containing the indices of the sites of the subsystem)
    Mutual_Info_Renyi2(sub) # returns the Renyi2 mutual information of the two subsystem (you either give one subsystem as a vector of indices or a splitting link)
    EE(l, n)                # returns entanglement entropy at site l and first n singular values
    Hermiticity             # returns 1 if density matrix is really Hermitian and down to 0 for anti Hermitian density matrix
    HermiticityError        # returns 1 - Hermiticity for monitoring hermiticity deviation
    Linkdim                 # returns the maximum bond dimension of the representation

They are used like this

    measure(state, TraceError)

New state functions may be defined by

    stfunc = StateFunc(name, state->... )

Checks can be performed (useful for coherence tests)

    measure(state, Check(name, obs1, obs2))      # returns 3 values ob1, obs2 and |obs2 - obs1|
    messure(state, Check(name, obs1, obs2, tol)) # if |obs2 - obs1|>tol throw an error
    measure(state, Check(name, obs1, obs2), t)

In `Check`, obs may also be constants, vectors and function of time (like `t -> sin(t)`), in this case the simulation
time must be fed to `measure` as 3rd argument.

## Output

```@docs
output
log_msg
```