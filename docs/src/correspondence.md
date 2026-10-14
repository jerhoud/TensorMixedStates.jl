# Coming from QuTiP or ITensors

## From QuTiP

QuTiP holds states and operators as dense matrices, TMS as matrix product states and operators,
which reach much larger systems at the price of a truncation, see
[Checking the accuracy](@ref). The conventions are the same: on a qubit, QuTiP's `basis(2, 0)`
is `"Up"`, with ``\sigma_z = +1``, `sigmap()` is `Sp` and `sigmam()` is `Sm`, and a collapse
operator `c` of `mesolve` is `Dissipator(c)` in an evolver.

| QuTiP | TMS |
|:---|:---|
| `basis(2, 0)`, `basis(2, 1)` on each of `n` qubits | `State{Pure}(System(n, Qubit()), "Up")`, `"Dn"` |
| `tensor(basis(2, 0), basis(2, 1))` | `State{Pure}(System(2, Qubit()), ["Up", "Dn"])` |
| `ket2dm(ψ)` | `mix(ψ)` |
| `sigmax()`, `sigmay()`, `sigmaz()` on qubit `i` | `X(i)`, `Y(i)`, `Z(i)` |
| `sigmap()`, `sigmam()` | `Sp`, `Sm` |
| `destroy(d)` on boson `i` | `A(i)` on `Boson(d)` |
| `qeye(d)` | `Id` |
| `tensor(sigmaz(), sigmaz())` acting on qubits `i` and `j` | `Z(i) * Z(j)` |
| `mesolve(H, ρ0, tlist, c_ops = [c1, c2], e_ops = [o])` | `evolve(-im * H + Dissipator(c1) + Dissipator(c2), t, ρ0; nsteps)`, then `measure(ρ, o)` |
| `sesolve(H, ψ0, tlist)` | `evolve(-im * H, t, ψ0; nsteps)` |
| `steadystate(H, c_ops)` | `steady_state(-im * H + sum(Dissipator.(c_ops)), ρ; algo = Dmrg(nsweeps))` |
| `H.groundstate()` | `ground_state(H, ψ; algo = Dmrg(nsweeps))` |
| `expect(o, ρ)` | `expect(ρ, o)`, or `measure(ρ, o)` |
| `ρ.ptrace([0, 2])` | `partial_trace(ρ, [1, 3]; keep = true)`, sites counted from 1 |
| `entropy_vn(ρ.ptrace([0]))` | `vonneumann_entropy(ρ, [1])` |

A collapse operator carries its rate, `sqrt(γ) * c` in QuTiP as in TMS: `Dissipator(sqrt(γ) * c)`
is `γ * Dissipator(c)`. A QuTiP operator is a matrix on the whole space; a TMS operator is placed
on sites, `Z(3)`, and an operator of one site given alone, `Z`, is measured on every site. Where
`mesolve` returns the values at each time of `tlist`, a TMS program evolves step by step and
measures between the steps, or measures in an `Evolve` phase, see
[High Level Interface](@ref).

## From ITensors

TMS is built on ITensors and ITensorMPS: a pure state holds an `MPS`, `state.state`, on the
indices `siteinds(Pure, system)`, and its algorithms call those of ITensorMPS. What it adds is
the density matrix, held as an MPS too, the operators written as expressions rather than
`OpSum`, the Lindbladians, the channels, the measurements and the simulations.

| ITensors and ITensorMPS | TMS |
|:---|:---|
| `s = siteinds("Qubit", n)` | `System(n, Qubit())` |
| `MPS(s, "0")` | `State{Pure}(System(n, Qubit()), "Up")` |
| `random_mps(s; linkdims = 8)` | `random_state(Pure, system, 8)` |
| `os = OpSum(); os += "Z", 1, "Z", 2; MPO(os, s)` | `Z(1) * Z(2)`, made an MPO when it is used, or `make_mpo(state, op)` |
| `dmrg(H, ψ0; nsweeps, maxdim, cutoff)` | `ground_state(H, ψ0; algo = Dmrg(nsweeps), limits = Limits(maxdim, cutoff))` |
| `tdvp(H, -im * t, ψ0; nsteps)` | `evolve(-im * H, t, ψ0; nsteps)` |
| `expect(ψ, "Z")`, `correlation_matrix(ψ, "Z", "Z")` | `expect1(ψ, Z)`, `expect2(ψ, (Z, Z))`, or `measure(ψ, [Z, (Z, Z)])` |
| `apply(gates, ψ; cutoff)` | `apply(gates, ψ; limits)`, the gates written as operators, `controlled(X)(1, 2) * H(1)` |
| `inner(ψ, ϕ)` | `inner(ψ, ϕ)` |

`tdvp` of ITensorMPS computes ``e^{t A}`` for its operator ``A`` and time ``t``, so that a real
time evolution under ``H`` takes `-im * t`; `evolve` takes the evolver `-im * H` and the real
time `t`, so that the same `-im` is written on the operator. A Lindbladian is an evolver as
well, on a mixed state, which ITensorMPS has no notion of.
