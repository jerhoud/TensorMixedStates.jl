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

When it builds an MPO, TMS compacts the operator: its terms of several sites are gathered into blocks of channels in which they share what they have in common, and
the operators of each site are compared through their matrices on the sites of the system, so
that `X*Y` and `im * Z` on a qubit, or `N` and `(1 - Z) / 2`, are known to be related. The MPO
keeps the triangular form: one channel runs along the chain for the terms not yet started,
another for those already finished, and in between each link has as many channels as the rank
of the operator across it, once its parts that are the identity on either side are taken out.
No MPO of this form has fewer.

For an operator whose terms are short ranged, this is what one channel per term gives already:
the Ising chain `sum(Z(i)Z(i+1) for i in 1:n-1)` has bond dimension 3 whatever the number of
sites. For a long ranged one, the difference is large. Measured on
`sum(Z(i)Z(j) for i in 1:n-1 for j in i+1:n)`:

| sites | terms | one channel per term | compacted |
|------:|------:|---------------------:|----------:|
| 10    | 45    | 27                   | 3         |
| 20    | 190   | 102                  | 3         |
| 40    | 780   | 402                  | 3         |

Three channels is also what `ITensorMPS`' `OpSum` finds there, with an SVD. Each term of a time
dependent evolver is compacted on its own, so that it keeps its time function, and changing
the coefficients from one step to the next rebuilds the tensors of the MPO alone.

A product of sums would develop into as many terms as the product of theirs: the moment `M^4` of
`M = sum(Z(i) for i in 1:n)` has of order n⁴. When that outnumbers the terms of its factors
together, TMS builds it from the MPOs of its factors instead, multiplied site by site and
compacted the same way, equal factors in a row by squaring, `H^4` being `(H^2)^2`; on a density
matrix, a gate or a dissipator of an operator of many terms is made likewise of its sides `Kρ`
and `ρK†`. The cost is then linear in the number of sites: `expect(ψ, M^4)` on 24 sites takes
0.01 s, where its 10 903 terms took 2.2 s. The same algebra is at hand on the `PreMPO`s
themselves: `pre'`, `a * b`, `a ^ p`, `a + b`, `a - b`, `c * a` and `sum(v)`.

Two things follow from it in practice:

- The WI and WII approximations of [`ApproxW`](@ref) need the triangular form, and their error,
  of order τ², depends on how the operator is split into terms, not on the operator alone.
  Compacting keeps the first and the last site of every term, unless an operator of a site is a
  combination of other ones and of the identity there, as `N` when `Z` is used too, or `N` with
  the Jordan-Wigner strings `F` of fermions: part of a term then goes to terms of fewer sites,
  and the approximations change, while remaining of the same order.
- `expect` prepares the operator it measures, compacted, at each call, and a `Measure` keeps
  that preparation for the calls of `measure` that follow: on a long ranged operator measured
  many times, a `Measure` built once saves most of the cost, as a `PreMPO` does for the MPO of
  an evolver.

```@docs
PreMPO
make_mpo
```
