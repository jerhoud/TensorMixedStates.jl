# Choosing and checking the parameters

## Choosing them

Three parameters set the error of a result, and the user chooses all three: the truncation of
the states, `limits`, the time step of an evolution, `nsteps` or `time_step`, and the number of
sweeps of a search, `Dmrg(nsweeps)`. No value is right in general; the way to choose them is
the same everywhere:

1. Start small: a `maxdim` the machine handles in seconds, the default `cutoff`, which only
   drops what rounding leaves, a time step a tenth of the shortest time scale of the evolver,
   and for a search a few sweeps with `maxdim` growing along them, as
   `Limits(maxdim = [10, 20, 50, 100])`.
2. Measure along the run what tells the error, `MaxLinkdim`, `EntanglementEntropy`, and on a
   density matrix `TraceDeviation`, `HermiticityDeviation` and `Purity`, see below.
3. Improve one parameter at a time, doubling `maxdim`, halving the time step, adding sweeps,
   and run again: the digits that do not move are those the result has. When `MaxLinkdim`
   reaches `maxdim`, the truncation is what limits the result.
4. Keep the cheapest parameters giving the digits needed, and say which checks they passed.

A density matrix needs a larger bond dimension than a wave function, that of a pure state
already the square of the bond dimension of its wave function, and a conserved quantity declared on the sites,
see [Conserving a quantity](@ref), makes a given bond dimension cheaper.

## Checking the accuracy

Three main approximations set the error of a result, and the user chooses all three: the
truncation of the states, given by `limits`, the time step of an evolution or of
`thermalize`, and the convergence of a search, the sweeps of its `Dmrg`. A
result can be trusted once it no longer moves when each of them is improved: larger limits, a
smaller time step, more sweeps. The measurements taken along a run tell which of them limits
the accuracy. Two more have defaults: the tolerance of the Krylov method of each local step,
see `Krylov`, and the truncation `mpo_limits` of the ``L^\dagger L`` of `steady_state`, none by
default; truncated, it makes the value `steady_state` returns the residual of the truncated
operator.

### Truncation

`MaxLinkdim` tells when the bond dimension reaches `maxdim`, from which point the results
depend on it, and `EntanglementEntropy(cut)`, the operator space entanglement entropy on a
mixed state, shows the entanglement that makes the bond dimension grow. A truncation by
`cutoff` shows in neither, and only running again tells it. `Limits()` discards the singular
values below about 1e-13 of the norm, which keeps rounding from making the bond dimension
grow and may leave an error of that order: the truncation that makes a computation affordable
is a choice, checked by running again with a larger `maxdim` and a smaller `cutoff`.

### Time step

An evolution is checked by halving its time step, and by comparing `Tdvp` with `ApproxW`,
whose errors have different origins. `ApproxW` of order `k` makes an error of order
``\tau^k`` in the time step ``\tau`` when the evolver does not depend on time. When it does,
`Tdvp` and `ApproxW` of order 1 or 2 take its functions at the middle of each step, which caps
them at order 2, while `ApproxW` of order 3 or 4 keeps the order of its approximation, see
[Time dependent evolvers](@ref). `Tdvp` projects the evolution on the states of the bond
dimension the state has, its steps on two sites letting that dimension grow, and it enlarges
the bonds of a state whose bond dimension is below 5 by a global Krylov expansion before a
step: from a state of small bond dimension, as a product state, or for `thermalize`, compare
with `expand_period = 1`, which expands before every step. A thermal state has checks of its
own, see [Thermal states](@ref).

### Indicators of a density matrix

The exact density matrix is hermitian under a Lindbladian, under gates and noisy gates, which
act on both of its sides, and for a thermal or a steady state, and its trace stays one when the
evolution preserves the trace, as a Lindbladian does. What deviates from that is error, which
three measurements of the state follow along a run:

- `TraceDeviation`, ``1 - \mathrm{tr}\,\rho``, the drift of the trace. Expectation values are
  divided by the trace, which hides the drift but not the error it comes from.
- `HermiticityDeviation`, the squared norm of the anti-hermitian part of ``\rho``, relative to
  that of ``\rho``. That part being error, the square root of `HermiticityDeviation` is a lower
  bound of the error relative to the norm of ``\rho``, in the Hilbert-Schmidt norm: a value of
  `1e-6` proves an error of `1e-3` at least. It bounds the error from below only, an error
  without an anti-hermitian part going unseen. `hermitianize_period` removes the
  anti-hermitian part, and with it what `HermiticityDeviation` shows: measure it before deciding
  to hermitianize.
- `Purity`, ``\mathrm{tr}(\rho^\dagger \rho) / (\mathrm{tr}\,\rho)^2``, lies between ``1/D``
  and 1 for a density matrix of dimension ``D``. Above 1, the state is no longer a density
  matrix: it has lost its positivity, which nothing in a matrix product state guarantees,
  unless `HermiticityDeviation` shows that it has lost its hermiticity.

### Conserved quantities

A quantity declared conserved on the sites, see [Conserving a quantity](@ref), keeps the
errors out of the sectors the exact state never reaches: the tensors have no blocks there, so
that neither truncation nor rounding can put anything in them. A pure state keeps its charge
exactly, a mixed state under a weak symmetry keeps no coherence between two sectors, and
under a strong one stays in its sector. A quantity the exact evolution conserves without its
being declared, as the energy under a Hamiltonian, which no site can declare, is still a
check: its drift along a run is error, though its constancy proves nothing.

### Searches

A search runs the sweeps of its `Dmrg`, and stops earlier once the energy changes by less than
its `tol` between two sweeps, which a search stuck in a
metastable state does as well. `Variance(hamiltonian)`, zero exactly for an eigenstate, checks
the state found, see `variance`, but does not tell the ground state from another eigenstate,
which searches from other states, or in other sectors, do. Extrapolating the energy to
zero variance as `maxdim` grows estimates the exact energy, and the error of the last one.

The value `steady_state` returns, the residual ``\|L\rho\|^2`` for ``\rho`` of unit
Hilbert-Schmidt norm, is zero for a steady state. The state it finds is hermitian when the
steady state is unique and the search has converged, and it warns when its `HermiticityDeviation`
exceeds `1e-6`: a state that is not hermitian comes from a search that has not converged, or
from a Lindbladian with several steady states, of which dmrg returns any combination: one in
each sector of a quantity that its Hamiltonian and its jump operators all commute with, for
instance. The one a system reaches depends on the state it starts from: it is the limit of
its evolution in time, when that converges, which `evolve` gives. When the quantity is one a
site can conserve, declaring it `strong` keeps `steady_state` in the sector of the starting
state, see [Conserving a quantity](@ref). A hermitian state does not prove the steady state
unique.
