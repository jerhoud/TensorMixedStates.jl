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

Read that `quickref.md`: the whole library on one page, version matched. `upgrading.md`, next
to it, lists what changed since version 1. In Julia, `?name` gives the help of any name.

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
ham = -sum(Z(i) * Z(i + 1) for i in 1:n-1) - 0.8 * sum(X(i) for i in 1:n)   # not H: H is Hadamard
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
- In a product of gates the rightmost acts first.
- A density matrix is `State{Mixed}`; dissipators and channels need one, `mix(ψ)` makes one,
  of the square of the bond dimension of `ψ`: from a ground state of bond dimension 100, a
  density matrix of 10000, which no machine evolves. Truncate `ψ` first, or start mixed.
- `H`, `X`, `S`, `T`, `N` are operators of `Qubits`: do not name variables after them.
- `measure` returns a vector of pairs `name => value`, even for one measurement, and an
  operator of one site, `Z`, gives one pair `"Z" => [a value per site]`. The value of a
  Hermitian operator is real, of any other complex; `Trace` is complex. For one number call
  the function: `expect(ψ, X(1) * X(2))`, `expect1(ψ, Z)`, `purity(ρ)`,
  `entanglement_entropy(ψ, cut)`, which returns `(entropy, spectrum)`. `expect`, `expect1` and
  `trace` are complex as soon as the state is, even for a Hermitian operator: take `real`.
- `ground_state` and `steady_state` need `algo = Dmrg(nsweeps = …)`; there is no `dmrg`,
  `tdvp` or `approx_W` function in version 2. An error on a name of version 1 ends with a line
  saying what to write instead.
- Sites are numbered from 1. `Level => n` is the basis state `n` counted from 0.

## Long runs

For a run of hours, write it as a simulation, `runTMS(SimData(name = "…", phases = [...]))`,
see the quick reference: it writes its measurements as tables in a directory, keeps checkpoints
and resumes when run again. The line it prints at the end says whether it is done and how many
warnings its `log` holds; read them.
