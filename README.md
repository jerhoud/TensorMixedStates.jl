<div align="center">
  <img src="docs/src/assets/logo.png" alt="TensorMixedStates" width="200">
</div>

# TensorMixedStates.jl

[![CI](https://github.com/jerhoud/TensorMixedStates.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/jerhoud/TensorMixedStates.jl/actions/workflows/CI.yml)
[![Documentation](https://github.com/jerhoud/TensorMixedStates.jl/actions/workflows/documentation.yml/badge.svg)](https://github.com/jerhoud/TensorMixedStates.jl/actions/workflows/documentation.yml)
[![pages-build-deployment](https://github.com/jerhoud/TensorMixedStates.jl/actions/workflows/pages/pages-build-deployment/badge.svg)](https://github.com/jerhoud/TensorMixedStates.jl/actions/workflows/pages/pages-build-deployment)
[![codecov](https://codecov.io/gh/jerhoud/TensorMixedStates.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/jerhoud/TensorMixedStates.jl)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![License: GPL v3+](https://img.shields.io/badge/license-GPL--3.0--or--later-blue.svg)](LICENSE)

A Julia library to simulate closed and open quantum systems, pure and mixed states, using
matrix product states.

TMS represents the density matrix of the system as a matrix product state. It offers
manipulations of systems and states, a rich set of sites and operators that the user can
extend, and the algorithms that go with them: ground states by DMRG, Hamiltonian and
Lindbladian evolution by TDVP and others, and the application of gates, noisy ones included.
Being built on ITensor, it computes fast and runs in parallel naturally. Its syntax for
operators is meant to make observables, gates, Hamiltonians and Lindbladians easy to write,
and an optional high level interface expresses a whole simulation in a few lines.

## Installation

```julia
]add TensorMixedStates
```

TMS requires Julia 1.10.5 or later. Note the `]`, which enters the package manager.

## A first simulation

What TMS is for is open systems: the state it carries is the density matrix, so dissipation
is written directly into the evolution rather than approximated away. Here six qubits evolve
under a transverse field Ising hamiltonian while each one decays, and the high level
interface describes the run as a list of phases, executes them, takes the measurements and
writes everything to disk.

```julia
using TensorMixedStates, .Qubits

runTMS(SimData(
    name = "dissipative_ising",
    description = "six qubits under a transverse field Ising hamiltonian, each decaying at rate 0.2",
    phases = [
        CreateState(type = Mixed(), system = System(6, Qubit()), state = "Up"),
        Evolve(
            duration = 2.0,
            time_step = 0.1,
            algo = Tdvp(),
            limits = Limits(maxdim = 64),
            evolver = -im * (-sum(Z(i)Z(i + 1) for i in 1:5) - sum(X(i) for i in 1:6))
                      + sum(Dissipator(sqrt(0.2) * Sm)(i) for i in 1:6),
            measures = "data" => [Z, Purity],
        ),
    ],
))
```

`type = Mixed()` is what makes the state a density matrix, and the `Dissipator` terms added
to the hamiltonian are what turn the evolution into a Lindblad equation. The run writes a
`dissipative_ising` directory: `log` for what happened, `data` for the measurements,
`description` and `stamp` for what was asked and when, and `prog.jl`, a copy of the script
that produced it. The magnetization on the six sites and the purity, at the start and at the
end of the evolution:

```
Z        0.1    0.94098941    0.94117682   0.94117682   0.94117682   0.94117682    0.94098941
Purity   0.1    0.78980701
...
Z        2     -0.095182096  -0.19011196  -0.14498043  -0.14498043  -0.19011196   -0.095182096
Purity   2      0.1086794
```

The qubits start pure and pointing up; by the end the magnetization has reversed and the
purity has fallen to 0.11, which is a state no pure state code could have represented.

Everything is available directly too, without the phases: build a `State`, call `tdvp` or
`dmrg` on it, and read `expect1`, `expect2` or `entanglement_entropy` off the result. The
manual covers both styles, and [`examples/`](examples) holds eleven working scripts.

## Documentation

<https://jerhoud.github.io/TensorMixedStates.jl>

## Article of reference

TMS is described in *TensorMixedStates: A Julia library for simulating pure and mixed
quantum states using matrix product states*,
[SciPost Phys. Codebases **72** (2026)](https://doi.org/10.21468/SciPostPhysCodeb.72).
SciPost asks that it be cited together with the
[codebase release](https://doi.org/10.21468/SciPostPhysCodeb.72-r1.28) it documents. BibTeX
entries for both are in the [documentation](https://jerhoud.github.io/TensorMixedStates.jl/stable/#References),
and GitHub's "Cite this repository" button reads them from [CITATION.cff](CITATION.cff).

If TMS helped produce the data of a publication, a word in your acknowledgements would be
warmly appreciated. And please do not hesitate to get in touch, at the address given in the
article, to say what you are doing with TMS, how you find it and what you would need: such
messages are what shapes what comes next.

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md), and [CHANGELOG.md](CHANGELOG.md) for what changed.

## Licence

TensorMixedStates is distributed under the [GNU General Public License, version 3 or
later](LICENSE). It builds on ITensor, which is distributed under the Apache License 2.0.
