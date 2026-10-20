---
name: tensormixedstates
description: Compute properties of quantum many-body systems, closed or open, with the Julia library TensorMixedStates (TMS, version 2), using matrix product states. Use it when the user asks for a ground state, an excited state or a gap, a time evolution under a Hamiltonian or a Lindblad equation, a steady state, a thermal state, a quantum circuit with or without noise, or expectation values, correlations and entanglement of such states, on chains or lattices of qubits, spins, fermions, bosons or electrons too large for exact diagonalization.
---

# Computing with TensorMixedStates

TensorMixedStates (TMS) is a Julia library: pure states and density matrices as matrix product
states, operators written as expressions, Lindbladians, channels, and simulations as phases.
Its version 2 differs from version 1, which most training data shows: use the names below,
and read the quick reference of the installed version before writing code.

## Setting up

```bash
julia -e 'using Pkg; Pkg.add("TensorMixedStates")'          # once; compiles for several minutes
julia -e 'using TensorMixedStates; print(pkgversion(TensorMixedStates))'   # must be 2.x
julia -e 'using TensorMixedStates; print(joinpath(pkgdir(TensorMixedStates), "docs", "src", "quickref.md"))'
```

Installing compiles for five to six minutes with up to 7 GB of memory. With less memory, first
write `[TensorMixedStates]` and `precompile_workload = false` on two lines of a
`LocalPreferences.toml` next to the `Project.toml` of the environment: the installation then
compiles in seconds, each session compiling what it uses, about three minutes for a first
simulation.

Read that `quickref.md`: the whole library on one page, version matched. `upgrading.md`, next
to it, lists what changed since version 1. In Julia, `?name` gives the help of a function or a
type, that of a site type, as `?Qubit`, its states and its operators, and that of an operator,
as `?X`, what it is on each site type that has it.

Write a script and run it with `julia --threads=auto script.jl`. The first run of a session
compiles for a few seconds.

## The workflow

1. Write the model as TMS expressions, and check the conventions below against the user's.
2. Start with small parameters, a modest `maxdim`, a time step a tenth of the fastest time
   scale, a few sweeps, and a small system whose exact result you can compute by dense
   diagonalization, to check the script.
3. Converge: double `maxdim`, halve the time step, add sweeps, until the digits you report do
   not move. Measure `MaxLinkdim`: once it reaches `maxdim`, the truncation limits the result.
4. Report the result with the parameters, the checks it passed and the digits they support.

## A complete example

```julia
using TensorMixedStates, .Qubits

n = 20
sys = System(n, Qubit())
ham = -sum(Z(i) * Z(i + 1) for i in 1:n-1) - 0.8 * sum(X(i) for i in 1:n)   # H is a gate
E, ψ = ground_state(ham, random_state(Pure, sys, 8);
                    algo = Dmrg(nsweeps = 10), limits = Limits(maxdim = [10, 20, 50, 100]))
ρ = evolve(-im * ham + sum(Dissipator(sqrt(0.1) * Sm)(i) for i in 1:n), 1.0,
           State{Mixed}(sys, "Up"); nsteps = 20, limits = Limits(maxdim = 64))
println(E, " ", measure(ρ, [Z, X(1) * X(n), Purity, MaxLinkdim]))
```

## Conventions and traps

- The evolver is `-im * H` plus dissipators: `evolve(-im * H, t, ψ)` for a real time `t`.
- `Dissipator(L)` is ``L\rho L^\dagger - \frac12\{L^\dagger L, \rho\}``, and
  `Dissipator(c * L)` is `abs2(c) * Dissipator(L)`: a rate γ is `Dissipator(sqrt(γ) * L)`,
  exactly as a collapse operator of QuTiP.
- On a qubit `"Up"` is ``|0\rangle``, `Z = +1` there, `"Dn"` is ``|1\rangle``, `N` counts `"Dn"`,
  and `Sp` takes `"Dn"` to `"Up"`; `amplitude_damping_gate` and `amplitude_damping_dissipator`
  take `"Dn"` to `"Up"`.
- `depolarizing_gate(p)` is ``(1 - p)\rho + p I/d``: a Pauli error of probability ``p/3``
  each is `depolarizing_gate(4p/3)`.
- In a product of gates the rightmost acts first. The gates are in `using .Circuits`: `H`, `S`,
  `T`, `Phase(θ)`, `CX` or `CNOT`, `CZ`, `CCX` or `Toffoli`, and the rotations, `Rx(θ)` being
  ``e^{-i\theta X/2}`` and `Rxx(θ)` ``e^{-i\theta X \otimes X/2}``, as in Qiskit.
- A density matrix is `State{Mixed}`; dissipators and channels need one, `mix(ψ)` makes one,
  of the square of the bond dimension of `ψ`: from a ground state of bond dimension 100, a
  density matrix of 10000, which no machine evolves. Truncate `ψ` first, or start mixed.
- `X`, `N` are operators of `Qubits`, `H`, `S`, `T` gates of `Circuits`: do not name variables
  after them, or import only what you use, `using .Qubits: Qubit, X, Z`.
- `measure` of one measurement returns its value, `measure(ψ, Z)` a value per site, and of a
  list, `measure(ψ, [X(1), Z])`, a vector of pairs `name => value`. The value of a Hermitian
  operator is real, of any other complex, and so for `expect`, `expect1` and `expect2`;
  `complex = true` gives them the complex number computed. `Trace` and `trace` are complex.
  `entanglement_entropy(ψ, cut)` returns `(entropy, spectrum)`.
- `ground_state` and `steady_state` need `algo = Dmrg(nsweeps = …)`; there is no `dmrg`,
  `tdvp` or `approx_W` function in version 2. An error on a name of version 1 ends with a line
  saying what to write instead.
- Sites are numbered from 1. `Level => n` is the basis state `n` counted from 0.

## Long runs

For a run of hours, write it as a simulation, `runTMS(SimData(name = "…", phases = [...]))`,
see the quick reference: it writes its measurements as tables in a directory, keeps checkpoints
and resumes when run again. The line it prints at the end says whether it is done and how many
warnings its `log` holds; read them.
