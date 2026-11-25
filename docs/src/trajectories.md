# Quantum trajectories

A density matrix is an MPS of the squares of the site dimensions, and its bond dimension grows
with the entanglement of the mixed state, which may be far larger than that of the pure states
it is made of. A quantum trajectory holds one of these pure states, a wave function drawn at
random: the average of ``|\psi\rangle\langle\psi|`` over the trajectories is the density matrix
that `State{Mixed}` would hold, the dissipators and the channels drawing jumps or Kraus
operators where the density matrix takes their sum. The representation is `Trajectory`, as in
`State{Trajectory}`: a trajectory is a ket, and what holds for a pure state holds for it, but
for the draws and the weights, and the measurements that are not linear in the density matrix
are those of each pure state, see [What a trajectory measures](@ref).

## Making a trajectory

`State{Trajectory}` takes the descriptions of `State{Mixed}` and draws a pure state from them,
with the random number generator `rng`, the global one by default:

```julia
using TensorMixedStates, .Qubits, .Fermions, Random

sys = System(4, Qubit())
t = State{Trajectory}(sys, ["Up", "FullyMixed", [0.7 0.2; 0.2 0.3], "Dn"]; rng = Xoshiro(1))
t = State{Trajectory}(sys, Mixture([0.5 => "Up", 0.5 => "Dn"]))
t = State{Trajectory}(System(4, Fermion(conserve = N)), FullyMixed(N => 2))
t = random_state(Trajectory, sys, 16)
t = State{Trajectory}(ψ)            # the pure state ψ, a ground state for instance
```

- A pure local state is taken as it is; a density matrix gives one of its eigenvectors, with its
  eigenvalue for probability.
- A `Mixture` gives one of its terms, with its weight times the traces of its local states, then
  its local states as above.
- `FullyMixed(N => 2)` gives a basis state of the sector, all with the same probability.
- `random_state(Trajectory, …)` takes the arguments of `random_state(Mixed, …)`, draws the same
  purification, and measures its ancillas: the average of the trajectory is the random mixed
  state drawn from that purification.

On a system conserving a quantity, a trajectory has a definite charge of it, like a pure state:
a local density matrix with coherences between its charges is refused, and so is, when it is
conserved strongly, a description spreading over several of its values. A system conserving
nothing imposes nothing, not even the parity of fermions: to keep real fermions to states of
definite parity, declare it, `Fermion(conserve = parity(N))`, see
[Conserving a quantity](@ref).

## Evolving and applying

`evolve` takes the evolver of a mixed state, `-im * H` plus dissipators. Between jumps, the
trajectory evolves under the effective operator ``-iH - \frac12 \sum_L L^\dagger L``, and a jump
happens at the time ``s`` where the probability of none, ``\exp(-\int_0^s r)``, reaches a number
drawn at random, ``r = \sum_L \langle L^\dagger L\rangle`` being the rate of the jumps, taken
linear over each step between its values at both ends; the operator ``L`` is drawn among those
of the dissipators with its rate ``\langle L^\dagger L\rangle`` at that time, and the rest of the
step drawn again, so that a step may hold several jumps. The trajectory is then brought back to
its trace or, when the operators do not keep it, weighed by what they keep.
Besides the error of the method of the steps, the average of the trajectories has one from the
rate taken linear, which halving the step shows: a rate that rises and falls back within a
step, nearly zero at both ends, draws no jump in it. An `ExactState{Trajectory}`, whose norm no
truncation lowers, weighs its trajectory `:state` by the norm the evolution leaves over the
probability drawn with, which corrects the rate taken linear where it is not zero, at the cost
of weights different from one. A rate zero at both ends of a step but not within it still
draws no jump there, and one nearly zero at both ends gives a few jumps of very large weights:
take steps short enough for the rate. Its other trajectories carry no weight and are brought
back to their trace, as on a `State`.

```julia
L = -im * sum(X(i) for i in 1:4) + sum(Dissipator(sqrt(0.2) * Sp)(i) for i in 1:4)
t = evolve(L, 1.0, State{Trajectory}(sys, "Dn"); nsteps = 50, rng = Xoshiro(2))
```

`apply` of a channel draws one of its Kraus operators with its probability, `rng` giving the
draws. The jumps and the Kraus operators are those the expression is written with: `L` for
`Dissipator(L)`, the gates of a sum of `Gate`, the projectors of a `Dephase`, or its phases
with `kraus = :phases`, the operators of
a `Reset`; a channel given by its matrix alone is decomposed by its Choi matrix, block by block
of the charges conserved. Two ways of writing one channel give the same average with different
variances: `depolarizing_channel(p; kraus = :reset)`, the default, resets the sites to the fully
mixed state with the probability `p`, `kraus = :pauli` applies the identity or a Pauli gate, on
qubits; `dephasing_dissipator(γ)` jumps to the eigenstates of the basis at the rate `γ`,
`dephasing_dissipator(γ; jumps = :phases)` applies `Z` at the rate `γ / 2` on a qubit, which
leaves the populations as they are. On an Ising ring of 8 qubits with disordered couplings in a
transverse field, with ``\gamma = 0.1``, the second gave ``\langle Z_1 Z_5\rangle`` at
``t = 3.2`` a spread 4 times smaller from one trajectory to another.

The algebra of states acts on the kets, as for pure states: `t1 + t2` is a superposition, not a
mixture, which prepares a state before any draw.

## In a simulation

A simulation making trajectories needs `samples`: each sample is one trajectory, running the
phases from the first `CreateState{Trajectory}` or `ConvertState{Trajectory}` on, the phases
before them running once, as the phases before `Realizations` do, see
[Samples and disorder](@ref). `samples = 1` runs a single trajectory.

```julia
using TensorMixedStates, .Qubits

n = 6
sys = System(n, Qubit())
ham = -sum(Z(i) * Z(i + 1) for i in 1:n-1) - 0.8 * sum(X(i) for i in 1:n)
damping = sum(Dissipator(sqrt(0.2) * Sp)(i) for i in 1:n)    # "Dn" to "Up", at the rate 0.2
runTMS(SimData(
    name = "damped",
    seed = 1,
    samples = 500,
    phases = [
        CreateState{Trajectory}(sys, "Dn"),
        Evolve(evolver = -im * ham + damping, duration = 2.0, time_step = 0.02,
               limits = Limits(maxdim = 32), measurements = "z" => [Z, :step]),
    ],
    final_measurements = EachSample("final") => [Z(1)],
))
```

The files hold the averages over the trajectories, with their errors, columns `c`, `err(c)`,
`std(c)`, `N` and `N_eff`, and `EachSample("final")` a row for each trajectory. Everything the
samples of a disorder do, the trajectories do: `tol` stops them once the errors are small
enough, `TMS_SAMPLES` prolongs a finished run, and the program launched several times at once
shares them between the processes. A pure state prepared once, a ground state for instance,
becomes the starting point of every trajectory by `ConvertState{Trajectory}()` after the phase
that computes it. The value `:jumps` counts the jumps the dissipators of the `Evolve` phases drew
in the trajectory `:state` since a phase put it there, 0 before the first: `"jumps" => :jumps`
gives their mean number, and `EachSample` that of each trajectory. The value `:jump_list` lists
them, each as `(time, dissipator, kraus)`: its time, the dissipator as written and placed, as
`"Sm(1)"`, and the index of the operator drawn among the Kraus operators of a channel,
`TensorMixedStates.kraus_operators(system, dissipator)`, 1 for an operator;
`EachSample("jumps.jsonl") => [:jump_list]` writes it for each trajectory, in JSON, and
`evolve(…; jumps = list)` gives the same list outside a simulation. The values `:log_weight`,
`:jumps` and `:jump_list` belong to the run: no phase of a simulation making trajectories may set
them.

An operator that does not keep the trace, a `Gate(K)` of a `K` that is not unitary or a
projection, a channel or an evolver that loses trace, changes the trace of the trajectory, which
is its weight. The averages weigh each trajectory by it, their errors being those of a ratio of
weighted means, `N_eff` saying how many trajectories count, and the column `log_weight` gives the
logarithm of the mean weight, an estimate of the logarithm of the trace of the mixed state, with
its error. Every row of a simulation of trajectories has that column, 0 before the trajectory,
`EachSample` giving the weight of each, and the measurements may read it as the value
`:log_weight`. The run carries the trace into that value at each step where it is far from one,
so that it neither underflows nor overflows. Only the trajectory `:state` carries the weight: a
phase changing the trace of another one is refused. `SaveState` saves the trajectory `:state`
with its weight, which `LoadState` gives back to a trajectory loaded into `:state`. A local
density matrix weighs the trajectory
drawn from it by its trace, as the norm of a vector does.

Every operation multiplies the trace of a trajectory by what it keeps of the trace of the mixed
state, computed before it and not read on the ket it leaves: one for what keeps the trace, the
sum of ``\langle K^\dagger K\rangle`` for a channel, ``\langle A^\dagger A\rangle`` for an
operator ``A``, the probability of a result imposed, so that truncation and rounding are never a
weight; a truncation keeps the trace. One exception, an evolution whose operators lose trace,
whose weight has no closed form, takes the norm its steps leave, truncation included, as does a
`PreMPO` given to `evolve` or `Evolve`, which is taken as it is: give the evolver itself.

`thermalize` and `Thermalize` take a trajectory drawn from a mixed state, as
`CreateState{Trajectory}(system, "FullyMixed")`, to ``e^{-\beta H/2} \psi`` normalized, its
weight multiplied by the squared norm: the averages are those of the thermal state, and the mean
weight is ``Z / \prod_i d_i``. The weights spread as ``\beta`` grows, which `N_eff` shows: at low
temperature, few trajectories count. The value `:log_trace`, the logarithm of the weight of one
trajectory, is refused in a measurement, its average not being that of the thermal state: the
column `log_weight` gives it.

A result imposed by `collapse`, as `collapse(t, Z(1) => 1)` in an `Update`, multiplies the
weight by its probability, a postselection whose probability `log_weight` then estimates. A
trajectory of probability zero is left null and stops its sample, which the log says: it counts
in the mean weight, but in no row after, whose column `N` counts the trajectories left; when
every trajectory stops, the rows after are written by none.

Outside a simulation, the trace of a trajectory is its weight: average `trace(t) * expect(t, A)`
and divide by the average of `trace(t)`. `thermalize` normalizes the trajectory it is given, as
it does a mixed state: the weight is then the exponential of its `log_trace` times the trace the
trajectory had.

With disorder, the function of `Realizations` returns the phases of the trajectory, from its
`CreateState{Trajectory}` or `ConvertState{Trajectory}`: each sample is then one realization of
the disorder and one trajectory. A phase making a trajectory before `Realizations` is refused,
since every sample would go on from the same trajectory, and so is a trajectory loaded by
`LoadState` or set by `Update` outside the samples: a trajectory is loaded among the phases
`Realizations` returns.

## What a trajectory measures

The average over the trajectories gives what is linear in the density matrix: the expectation
values, the probabilities of `sample` and `collapse`, the trace, the fidelity to a pure state.
Anything else is measured on the pure state each trajectory holds, its ket normalised, and
averaged with the weights: the entropies, the entanglement spectrum and the logarithmic
negativity of each trajectory, the average entanglement of a monitored system, not that of the
density matrix; `Purity` and `Norm` are 1; `Variance(A)` is the mean of the variances of `A` in
each trajectory, smaller than that of the density matrix, which the averages of `A^2` and `A`
give; `inner`, `Vector` and `fidelity` take the normalised ket, and `hs_fidelity`, between
mixed representations, refuses it. Such averages depend on the
way the channels and the dissipators are unravelled, see `dephasing_dissipator(γ; jumps)`. A
trajectory is not mixed by `mix`. An `ExactState{Trajectory}` holds its ket whole, without
truncation, as `backend = ExactState` makes it in a simulation. A ground state is searched in
`Pure`, then made a trajectory by `State{Trajectory}(ψ)` or `ConvertState{Trajectory}()`.

Not available yet: `SteadyState`, `ConvertState{Trajectory}` of a mixed
state, the product `⊗` of trajectories and `random_state(t, n)` of a trajectory `t`.
