# Measurements

## What is measured

`measure(state, measurement)` returns the value of one measurement, and
`measure(state, measurements; time, values...)` those of a list, in their order, as pairs
`name => value`, as `v[i]` and `v[[i]]` give an element and a vector. A measurement is an
operator, a measurement of the library, a symbol, or a function of the context, each possibly
named, and a list of them is a vector or a `Measure`, which keeps what it prepares for the
next call. The examples below use a four qubit state

```@setup measurements
using TensorMixedStates
using .Qubits
mystate = State{Mixed}(System(4, Qubit()), ["+", "Up", "+", "Dn"])
```

```julia
using TensorMixedStates, .Qubits

mystate = State{Mixed}(System(4, Qubit()), ["+", "Up", "+", "Dn"])
```

With an operator placed on sites, we get its expectation value

```@example measurements
measure(mystate, X(1)X(3)Z(4))
```

With an operator of one site, its expectation value on each site

```@example measurements
measure(mystate, X)
```

With a pair of operators of one site, the matrix of their correlations

```@example measurements
measure(mystate, (X, Y))
```

The library measures the state as a whole:

| measurement | value |
|:---|:---|
| `Trace` | the trace of the state, a complex number |
| `TraceDeviation` | `1 - trace`, which numerical errors tend to move, a complex number |
| `Purity` | ``\mathrm{tr}(\rho^\dagger\rho)/(\mathrm{tr}\,\rho)^2``, ``\mathrm{tr}(\rho^2)`` for a density matrix of trace one, 1 on a pure state |
| `Norm` | the norm of the state |
| `Hermiticity` | 1 for a Hermitian density matrix, down to 0 for an anti-Hermitian one |
| `HermiticityDeviation` | `1 - hermiticity`, the squared norm of the anti-Hermitian part relative to that of the whole |
| `Renyi2` | the Rényi entropy of order 2 of the state |
| `SubRenyi2(positions)` | the Rényi entropy of order 2 of the sites at `positions`, `1:k` for the first `k` |
| `MutualInfoRenyi2(positions)` | the Rényi-2 mutual information between the sites at `positions` and the others |
| `VonNeumannEntropy(positions)` | the von Neumann entropy of the sites at `positions`, from their reduced density matrix |
| `ReducedDensityMatrix(positions)` | the density matrix of the sites at `positions`, a complex matrix, see `reduced_density_matrix` |
| `LogNegativity(positions_a, positions_b)` | the logarithmic negativity between the sites of `positions_a` and those of `positions_b` |
| `EntanglementEntropy(cut)` | the entanglement entropy across the cut between sites `cut` and `cut+1`, the operator space entanglement entropy (OSEE) on a mixed representation |
| `EntanglementEntropy(cut, n)` | the same, followed by the first `n` values of the spectrum it is computed from, see `entanglement_entropy` |
| `Fidelity(ref)` | the fidelity with the reference state `ref`, or with the state of that name in a simulation, `Fidelity(:ground)` |
| `Overlap(ref)` | the inner product with the reference state `ref`, a complex number, or with a state named in a simulation |
| `Overlap(ref, op)` | the matrix element ``\langle ref \vert op \vert state \rangle`` of two pure states, or its vector for an operator of one site on every site, the correlations at two times of a reference evolved with the state |
| `Variance(op)` | the variance of the operator `op` |
| `MaxLinkdim` | the maximum bond dimension of the state |
| `MemoryUsage` | the memory the state occupies in bytes, the caches of the measurements included |

```@example measurements
measure(mystate, [Trace, Purity, EntanglementEntropy(2)])
```

Anything else is measured by a function of the context `c`, named: `c.state` is the state,
`c.time` the time given to `measure`, and the other values given to it are read as `c.energy`,
for instance. A symbol stands for the value of that name:

```@example measurements
measure(mystate, ["half_trace" => c -> trace(c.state) / 2, "cos" => c -> cos(2c.time), :energy];
        time = 0.3, energy = -1.5)
```

A list given after the name of a state, `:φ => [Z, (X, X), Purity]`, measures the state of that
name with all of the above, its measurements named after it, `φ.Z`, `φ.(X, X)` and `φ.Purity`,
as their columns are:

```@example measurements
measure(mystate, [Z(1), :reference => [Z(1), Purity]];
        reference = State{Pure}(System(4, Qubit()), "Dn"))
```

In a simulation, the values are those its phases provide: `:step` the number of steps the
phase has done, the sweeps of a search, `:energy` the energy of each sweep of `GroundState`, and
for `SteadyState` the eigenvalue of ``L^\dagger L``, zero for a steady state, and `:beta` and `:log_trace` the inverse temperature
and the logarithm of the trace of `Thermalize`, see `thermalize`. A state named in a
simulation, by `CopyState` for instance, is a value as well, and `Update` provides the values it
declares. A measurement reading a value that no phase provides is refused when the simulation is
written, see `SimData`.

A function of the context may return a `Check`, which compares a value with a reference:

```@example measurements
measure(mystate, "check" => c -> Check(real(expect(c.state, X(1))), 1.0; tol = 1e-8))
```

The value is then the value, the reference and their distance, and a distance beyond `tol` is an
error, or only a warning with `warn = true`. Without `tol`, the three are recorded and nothing
is refused.

## Real and complex values

The value of an operator is real when the operator is Hermitian, which is tested numerically,
and complex otherwise, whatever the state, so that the columns of a file do not change during a
run. A pair `(A, B)` is Hermitian when `A` and `B` are, commute on a site and are not both
fermionic: `(X, Y)`, whose diagonal is ``\langle XY \rangle = i \langle Z \rangle``, is
complex. The imaginary part of the value of a Hermitian operator is dropped, with a warning when
it is more than rounding, which a density matrix that lost its hermiticity gives; a simulation
gathers these warnings into one at the end of each phase. `Trace`, `TraceDeviation`, `Overlap`
and `ReducedDensityMatrix` are complex, the other measurements of the library real, and a
function of the context gives what it returns.

## Reference

```@docs
measure
Check
Measure
expect
expect1
expect2
variance
sample(::State{Pure})
probabilities
collapse
entanglement_entropy
entanglement_by_sector
renyi2
mutual_info_renyi2
reduced_density_matrix
vonneumann_entropy
log_negativity
Trace
TraceDeviation
Purity
Norm
Hermiticity
HermiticityDeviation
Renyi2
SubRenyi2
Fidelity
Overlap
Variance
MutualInfoRenyi2
VonNeumannEntropy
ReducedDensityMatrix
LogNegativity
EntanglementEntropy
MaxLinkdim
MemoryUsage
```

## [Output](@id measure-output)

`output` takes measurements and writes them to destinations, as the phases of `runTMS` do with
their `measurements`, given as pairs `destination => measurements`. The destinations are

- a file name: a text table, a line naming the columns, then a line per time measured, the
  numbers written with the format of the simulation, `SimData(format = "%.8g")` by default. A
  file holds one list of measurements, and another list needs another file. Its columns are
  those of its first line: a value that changes its cells later, as that of a state named after
  it, missing at first, needs a JSON Lines file.

  ```julia
  "magnetization" => [X, Z(1)]
  ```

- `"stdout"` or `"-"`, and `"stderr"`: a text table on the stream, written again with its new
  columns when the list changes; `""` writes nothing.

- a file name ending in `.json` or `.jsonl`: JSON Lines, an object for each time measured, the
  time and each measurement under its name, a vector as an array, a matrix as the array of its
  rows and a complex number as an object, as `{"time":0.1,"X":[0.5,0.5],"Trace":{"re":1.0,"im":0.0}}`.

  ```julia
  "data.jsonl" => [Purity, X(2)Z(3), (X, Y)]
  ```

- `Data(name)`: the rows are kept in memory, in the `data` of the `Simulation` that `runTMS`
  returns, which `data_to_frame` turns into a `DataFrame`.

  ```julia
  Data("mydata") => [TraceDeviation, X(1), Y]
  ```

- a `TensorMixedStates.Destination` of one's own.

The columns of a table are `time`, then a column per value: the name of the measurement, as
`X(1)*X(3)` for an operator placed on sites, `X(i)` for each site of an operator of one site,
`X(i)Y(j)` for each pair of sites of a pair, `m(i)` and `m(i,j)` when they are named `m`,
`m[i]` and `m[i,j]` for the entries of a vector or of a matrix of another measurement, and
`m.field` for those of a named tuple. A complex value takes the two columns `Re(m)` and
`Im(m)`. The files `runTMS` writes itself in the directory of the simulation cannot be
destinations: `log`, `stop`, `restart`, `error`, `running`, `stamp`, `description`, `prog.jl`,
`prog_args.json` and the checkpoint files.

```julia
mysim = runTMS(simdata)
df = data_to_frame(mysim.data["mydata"])
```

```@docs
output
log_message
TensorMixedStates.Destination
```
