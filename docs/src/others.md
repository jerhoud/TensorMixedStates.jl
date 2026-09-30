# Graphs and MPOs

## Graphs

Graphs are useful to describe interactions or gates to apply.

```@docs
line_graph
circle_graph
complete_graph
square_lattice
graph_base_size
```

## MPO

Matrix Product Operators (MPO) are used under the hood by TMS to act on the states, which it
represents as matrix product states. Every operator is converted to an MPO internally, except in `apply`, which applies gates one
by one, and in `expect`, `expect1`, `expect2` and `measure`, which contract the state
directly.

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
to the next without simplifying the operator again: only the tensors of the MPO are rebuilt,
from the terms kept, which is what makes time dependent evolvers cheap.

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
