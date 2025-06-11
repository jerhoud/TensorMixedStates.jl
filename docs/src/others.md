# Others

## Graphs

graphs are useful to describe interactions or gates to apply 

```@docs
line_graph
circle_graph
complete_graph
graph_base_size
```

## MPO

Matrix Product Operators (MPO) are used under the hood by TMS to operate on MPS (inner state representation).
Except for `apply` and `measure` all operators are converted to MPO internally.

```@docs
PreMPO
make_mpo
make_approx_W1
make_approx_W2
```

## Time dependent operators

For time evolution, it may be useful to have time dependent operators. In TMS, a time dependent operator
is described in the following way: a vector of indexed operators and a vector of time functions. For example,
to describe

```math
    h(t) = - e^{-t} \sum_{i=1}^{n-1} \sigma_x^i \sigma_x^{i+1} - \sin(t) \sum_{i=1}^n \sigma_z^i  
```
we use

    hs = -im * [ -sum(X(i)X(i+1) for i in 1:n-1), -sum(Z(i) for i in 1:n)]

and

    coefs = [ t -> exp(-t), t -> sin(t) ]

`h` is passed as usual to `tdvp` or `approx_W` and `coefs` is passed as a keyword argument called `coefs`.
When using `Simulation` the simulation time is used for `t`, for `State` the initial simulation time is passed as a keyword argument called `time_start` (which default to 0)

    tdvp(hs, duration, initial_state; coefs, time_start)

With the high level interface, one can use time dependent evolver for the Evolve phase with the following syntax

    evolver = hs => coefs
