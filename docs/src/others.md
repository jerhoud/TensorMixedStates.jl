# Others

## Graphs

graphs are useful to describe interactions or gates to apply 

```@docs
line_graph
circle_graph
complete_graph
square_lattice
graph_base_size
```

## MPO

Matrix Product Operators (MPO) are used under the hood by TMS to operate on MPS (inner state representation).
Except for `apply` and `measure` all operators are converted to MPO internally.

### How the MPO is built, and what it costs

TMS gives every term of a sum a channel of its own. A term acting on sites `i` to `j` opens
a channel at `i`, carries it across the links in between and closes it at `j`. Two more
channels run along the whole chain, one for the terms not yet started and one for those
already finished. The bond dimension of the MPO on a given link is therefore

```
2 + the number of terms crossing that link
```

For an operator whose terms are short ranged this is the best one can do. The Ising chain
`sum(Z(i)Z(i+1) for i in 1:n-1)` has exactly one term crossing each link, so its MPO has
bond dimension 3 whatever the number of sites.

For a long ranged operator it is another matter, because every term is kept separate rather
than merged with its neighbours. Measured on `sum(Z(i)Z(j) for i in 1:n-1 for j in i+1:n)`:

| sites | terms | MPO bond dimension in TMS | with an SVD compressed construction |
|------:|------:|--------------------------:|------------------------------------:|
| 10    | 45    | 27                        | 3                                   |
| 20    | 190   | 102                       | 3                                   |
| 40    | 780   | 402                       | 3                                   |

The last column is what `ITensorMPS`' `OpSum` gives on the same operator: it compresses the
construction with an SVD and finds the three channels that suffice. TMS does not compress,
so its bond dimension grows like the number of terms crossing the middle of the chain, here
`n²/4`. If you come from `ITensorMPS` and build a long ranged Hamiltonian, this is the
difference you will see, and the contraction cost follows it.

This is a deliberate trade. The construction leaves the MPO in the triangular form that the
WI and WII approximations of [`ApproxW`](@ref) need, which a compressed MPO no longer has.
And because each term keeps its own identity, its coefficient can be changed from one sweep
to the next without rebuilding anything, which is what makes time dependent evolvers
possible at no extra cost.

What follows from it in practice: the cost is paid per term, so it is worth writing an
operator with as few terms as possible. `simplify` is applied on the way to the MPO and will
merge what it can, but it cannot know that two terms you wrote separately describe the same
physics.

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

```@setup others
using TensorMixedStates
using .Qubits
n = 6
```

```@example others
hs = -im * [ -sum(X(i)X(i+1) for i in 1:n-1), -sum(Z(i) for i in 1:n)]
```

and

```@example others
coefs = [ t -> exp(-t), t -> sin(t) ]
```

`hs` is passed as usual to `tdvp` or `approx_W` and `coefs` is passed as a keyword argument called `coefs`.
When using `Simulation` the simulation time is used for `t`, for `State` the initial simulation time is passed as a keyword argument called `time_start` (which default to 0)

```julia
tdvp(hs, duration, initial_state; coefs, time_start)
```

With the high level interface, one can use time dependent evolver for the Evolve phase with the following syntax

```julia
evolver = hs => coefs
```
