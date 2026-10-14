# Quick reference

The whole library on one page, for a reader who knows quantum mechanics and wants a correct
computation quickly, a person or an AI. `?name` at the Julia prompt gives the details of any
name below.

## Loading

```julia
using TensorMixedStates, .Qubits, .Fermions, .Bosons   # the module of each site type: also
                                                      # .Spins, .Electrons, .Tjs, .Qudits, .Qbosons
```

Run scripts with `julia --threads=auto script.jl`. The first call of a session compiles, a few
seconds; the precompiled workload covers the common paths.

## Systems and states

```julia
sys = System(10, Qubit())                          # 10 qubits
mixed_sys = System([Qubit(), Boson(4), Fermion()]) # sites of several types
ψ = State{Pure}(sys, "Up")                         # a product state, the same on every site
ψ = State{Pure}(sys, [isodd(i) ? "Up" : "Dn" for i in 1:10])
ρ = State{Mixed}(sys, "FullyMixed")                # density matrices are Mixed
ρ = mix(ψ)                                         # |ψ⟩⟨ψ|
ψ = random_state(Pure, sys, 8)                     # random, bond dimension 8
ψ = ghz_state(sys, "Up", "Dn")                     # also dicke_state, w_state, dimer_state,
                                                   # fermi_sea, slater_state, fully_mixed,
                                                   # graph_state (Qubits), aklt_state (Spins)
```

A local state is a name of the site (`"Up"`, `"Dn"`, `"+"`, `"0"`, `"1"` on a qubit), a vector,
a density matrix, or `Level => n`, the basis state `n` counted from 0. Conservation: a site
`Qubit(conserve = N)`, `Fermion(conserve = N)`, on a mixed state also `conserve = strong(N)`;
`weaken(state)` drops a level of it.

## Operators

```julia
ham = -sum(Z(i) * Z(i + 1) for i in 1:9) - 0.5 * sum(X(i) for i in 1:10) # X(i): X on site i;
                                          # not H, the Hadamard gate of Qubits
hop = -sum(dag(C)(i) * C(i + 1) + dag(C)(i + 1) * C(i) for i in 1:9)      # fermions: strings
                                                                         # laid automatically
L = -im * ham + sum(Dissipator(sqrt(0.1) * Sm)(i) for i in 1:10)          # a Lindbladian
gates = controlled(X)(1, 2) * H(1)        # a product acts right to left: H first
noisy = depolarizing_gate(0.01)(3)        # a channel, for a mixed state
myop = named([1 1 ; 1 -1] / √2, "MyOp")   # an operator of one's own
```

`X`, `Y`, `Z`, `Sp`, `Sm`, `N` on qubits; `C`, `N` on fermions; `A`, `N` on bosons; `Sx`, `Sy`,
`Sz`, `Sp`, `Sm` on spins; `Id`, `F`, `Level`, `Swap` and `controlled` on every site. `⊗` makes an
operator of several sites, `X ⊗ X`, placed by `(X ⊗ X)(2, 5)`; `exp`, `sqrt`, `dag`,
`Proj("Up")` and `Proj(Z => 1)` work on operators.

## Conventions

| what | convention |
|:---|:---|
| the evolver | `-im * H` plus dissipators: ``d\psi/dt = -iH\psi``, ``\hbar = 1`` |
| a jump operator | `Dissipator(L)`, ``L\rho L^\dagger - \frac12\{L^\dagger L, \rho\}``; `Dissipator(c * L)` is `abs2(c) * Dissipator(L)`, so a rate ``\gamma`` is `γ * Dissipator(L)` |
| qubit states | `"Up"` = ``\vert 0\rangle``, `"Dn"` = ``\vert 1\rangle``; `Z` is +1 on `"Up"`; `N` counts `"Dn"`; `Sp` takes `"Dn"` to `"Up"` |
| amplitude damping | `amplitude_damping_gate(p)`, `amplitude_damping_dissipator(γ)`: `"Dn"` decays to `"Up"` |
| depolarization | `depolarizing_gate(p)` is ``(1-p)\rho + p\,I/d``: each Pauli error with probability ``p/3`` is `depolarizing_gate(4p/3)` |
| a product of gates | the rightmost factor acts first |
| values measured | real for a Hermitian operator, complex otherwise; `Trace` complex |

## Solvers

```julia
lim = Limits(maxdim = 100)          # truncation: cutoff (default 1e-26, rounding), maxdim, mindim
E, gs = ground_state(ham, random_state(Pure, sys, 8); algo = Dmrg(nsweeps = 10), limits = lim)
E1, ex = ground_state(ham, random_state(Pure, sys, 8); algo = Dmrg(nsweeps = 20), limits = lim,
                      orthogonal_to = [gs])                     # the first excited state
ψt = evolve(-im * ham, 1.0, gs; nsteps = 20, limits = lim)        # algo = Tdvp() by default
ρt = evolve(L, 1.0, mix(gs); algo = ApproxW(order = 4), nsteps = 20, limits = lim)
v, ρs = steady_state(L, mix(gs); algo = Dmrg(nsweeps = 10), limits = lim)
logZ, ρβ = thermalize(ham, 0.5, State{Mixed}(sys, "FullyMixed"); nsteps = 10, limits = lim)
ψg = apply(gates, State{Pure}(sys, "Up"); limits = lim)
```

`evolve` takes a time dependent evolver as `[A, B] => [f, g]`, real functions of the time, read
from `time` (keyword, default 0). The error is set by `limits`, the time step (`nsteps`) and the
sweeps: a result is trusted once it does not move when each is improved, see
[Checking the accuracy](@ref).

## Measurements

```julia
measure(ψt, [X(1) * X(10), Z, (X, X), Purity, EntanglementEntropy(5)])
# pairs: an operator on sites, an operator of one site on every site, a correlation matrix,
# measurements of the state
measure(ψt, ["E" => c -> real(expect(c.state, ham)), "t" => c -> c.time]; time = 1.0)
expect(ψt, X(1) * X(10)); expect1(ψt, Z); expect2(ψt, (X, X))
probabilities(ψt, Z(1) * Z(2)); sample(ψt, Z(1)); collapse(ψt, Z(1))  # (result, probability, state)
```

## Simulations

```julia
runTMS(SimData(
    name = "quench",                      # the directory of its files
    seed = 1,                             # the random draws, reproducible
    phases = [
        CreateState{Pure}(10, Qubit(), "Up"),
        GroundState(hamiltonian = ham, algo = Dmrg(nsweeps = 10), limits = lim),
        CopyState(:state => :ground),     # a named state, for Fidelity(:ground)
        Evolve(evolver = -im * (ham + 0.5 * sum(Z(i) for i in 1:10)), duration = 1.0,
               time_step = 0.05, limits = lim,
               measurements = "data" => [Z(1), Fidelity(:ground), :step]),
    ],
))
```

Other phases: `ToMixed`, `Weaken`, `PartialTrace`, `Gates`, `Thermalize`, `SteadyState`,
`SaveState`, `LoadState`, `Update`. A phase acts on `:state`, or on a named state by `on = :φ`,
and `Evolve` and `Gates` on several, `on = [:state, :φ]`; a list `:φ => [Z, Purity]` measures
the named state, and `Overlap(:φ, X)` gives ``\langle φ \vert X_i \vert state \rangle`` on every
site, the correlations at two times of `φ = A|ψ⟩` evolved with `|ψ⟩`. A destination is a file name (a table), `"name.jsonl"` (JSON
Lines), `"stdout"`, or `Data("name")` (kept in `sim.data`, `data_to_frame`). Running the program
again resumes a stopped run from its checkpoint; `runTMS(sd; max_time, checkpoint_interval)`
or the variables `TMS_MAX_TIME`, `TMS_CHECKPOINT_INTERVAL` set them; a file `restart` in the
directory, or `cleanTMS(sd)`, starts over. `runTMS(sd; output = devnull)` writes no files.

## When something fails

An error on a name or a keyword of version 1.x ends with a line saying what 2.0 writes instead,
see [Upgrading from 1.x](@ref). The log of a run, in its directory, holds its warnings, which
the line `runTMS` prints at the end counts.
