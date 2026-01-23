# Algorithms

## Ground state computation

Using DMRG, we can compute ground states for States and Simulations

```@docs
dmrg
DmrgObserver
```

## Steady state computation

Using DMRG on the square Lindbladian we can compute steady states for mixed States and Simulations

```@docs
steady_state
```

## Time evolution

Using TDVP or ApproxW, we can do time evolution (Hamiltonian or Lindbladian)

```@docs
tdvp
approx_W
TdvpObserver
ApproxWObserver
```

## Gate application

```@docs
apply
```