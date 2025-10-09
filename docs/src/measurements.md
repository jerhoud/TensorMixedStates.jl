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

    measure(state, Trace)              # returns the trace of the state
    measure(state, Trace2)             # returns the trace of the square of the state (alternate name: Purity)
    measure(state, TraceError)         # returns 1 - trace, usefull for monitoring trace deviations
    measure(state, Renyi2)             # returns the Renyi entropy of order 2 of the state
    measure(state, SubRenyi2([1, 3]))  # returns the Renyi entropy of order 2 of the subsystem (here sites 1 and 3)
    measure(state, EE(l, n))           # returns entanglement entropy at site l and first n singular values
    measure(state, Hermitianity)       # returns 1 if density matrix is really Hermitian and down to 0 for anti Hermitian density matrix
    measure(state, HermitianityError)  # return 1 - Hermitianity for monitoring hermitianity deviation

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