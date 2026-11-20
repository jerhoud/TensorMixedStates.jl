# Quantum trajectories

A density matrix is an MPS of the squares of the site dimensions, and its bond dimension grows
with the entanglement of the mixed state, which may be far larger than that of the pure states
it is made of. A quantum trajectory holds one of these pure states, a wave function drawn at
random: the average of ``|\psi\rangle\langle\psi|`` over the trajectories is the density matrix
that `State{Mixed}` would hold, the dissipators and the channels drawing jumps or Kraus
operators where the density matrix takes their sum. The representation is `Trajectory`, as in
`State{Trajectory}`: a trajectory is a ket, and what holds for a pure state holds for it, but
for the draws and for the measurements that are not linear in the density matrix, see
[What a trajectory refuses](@ref).

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

`evolve` takes the evolver of a mixed state, `-im * H` plus dissipators. At each step of length
``dt``, a jump happens with the probability ``\delta p = dt \sum_L \langle L^\dagger L\rangle``,
the operator ``L`` drawn among those of the dissipators with its rate ``\langle L^\dagger
L\rangle``, and replaces the step; otherwise the step evolves the trajectory under the effective
operator ``-iH - \frac12 \sum_L L^\dagger L``. The trajectory is then brought back to its trace.
This scheme is of first order in ``dt``: the average of the trajectories has an error of order
``dt``, which halving the step shows, and a step whose ``\delta p`` exceeds one is refused.

```julia
L = -im * sum(X(i) for i in 1:4) + sum(Dissipator(sqrt(0.2) * Sp)(i) for i in 1:4)
t = evolve(L, 1.0, State{Trajectory}(sys, "Dn"); nsteps = 50, rng = Xoshiro(2))
```

`apply` of a channel draws one of its Kraus operators with its probability, `rng` giving the
draws. The jumps and the Kraus operators are those the expression is written with: `L` for
`Dissipator(L)`, the gates of a sum of `Gate`, the projectors of a `Dephase`, the operators of
a `Reset`; a channel given by its matrix alone is decomposed by its Choi matrix, block by block
of the charges conserved. Two ways of writing one channel give the same average with different
variances: `depolarizing_channel(p; kraus = :reset)`, the default, resets the sites to the fully
mixed state with the probability `p`, `kraus = :pauli` applies the identity or a Pauli gate, on
qubits.

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
that computes it.

With disorder, the function of `Realizations` returns the phases of the trajectory, from its
`CreateState{Trajectory}` or `ConvertState{Trajectory}`: each sample is then one realization of
the disorder and one trajectory. A phase making a trajectory before `Realizations` is refused,
since every sample would go on from the same trajectory, and so is a trajectory loaded by
`LoadState` or set by `Update` outside the samples: a trajectory is loaded among the phases
`Realizations` returns.

## What a trajectory refuses

The average over the trajectories gives what is linear in the density matrix: the expectation
values, the probabilities of `sample` and `collapse`, the trace, the fidelity to a pure state.
The rest is refused with a message saying so: the purity, the entropies and the entanglement
spectrum, the logarithmic negativity, the fidelity to a mixed state, the measurement `Norm`,
inner products and matrix elements between states. `Variance(A)` is refused too: measure `A^2`
and `A`, whose averages give it. A trajectory is not mixed by `mix`, has no `Vector`, and is not
held by an `ExactState`. A ground state is searched in `Pure`, then made a trajectory by
`State{Trajectory}(ψ)` or `ConvertState{Trajectory}()`.

Not available yet: `SteadyState` and `Thermalize`, `ConvertState{Trajectory}` of a mixed
state, the product `⊗` of trajectories, `random_state(t, n)` of a trajectory `t`, and a channel
in one product with a gate that does not keep the trace.
