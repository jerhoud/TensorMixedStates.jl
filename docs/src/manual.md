# Manual

```@contents
Pages = ["manual.md"]
Depth = 3
```

## Import

To use TMS, you must first import it with

```@example manual
using TensorMixedStates
```

## Sites and Systems

The first step in using TMS is the definition of your quantum system. In TMS, a system is composed of a finite
number of sites numbered from 1 (1, 2, ..., N). These sites may be all identical or not.

There are eight different predefined types of site: `Qubit`, `Qudit`, `Fermion`, `Boson`, `Spin`, `Electron`, `Tj` and `Qboson`.

To use each of these sites and the corresponding predefined operators you need first to import the corresponding module.
For example to use qubits, you need to write

```@example manual
using .Qubits
```

Note the "." before the name and the "s" at the end.

To define a site just call the corresponding creator for example

```@example manual
s = Qubit()
```

Four site creators need an argument: `Qudit(dim)` and `Boson(dim)` for the dimension of
the local Hilbert space, `Spin(s)` for the spin, and `Qboson(q, dim)` for the deformation
parameter and the dimension. For example

```@example manual
using .Bosons, .Spins

s = Boson(4)
```

```@example manual
s = Spin(3/2)
```

You can now define a quantum system by declaring the sites it contains:

```@example manual
system1 = System(10, Qubit())
nothing # hide
```

gives you a system with 10 qubits. Systems may have different types of site, in this case you must feed `System` with an array of sites

```@example manual
using .Fermions

system2 = System([Qubit(), Boson(4), Fermion()])
nothing # hide
```

gives you a three site system.

## States

States may be in pure or mixed representation, these two possibilities are represented in TMS, by `Pure` or `Mixed`.

To create a state, we call the State creator

```@example manual
state1 = State{Pure}(system1, "Up")
nothing # hide
```

returns a pure up state in a 10 qubit system. Predefined local states are designated by their name. Here `"Up"` is a predefined state of site `Qubit`.

Sites need not be all in the same local state, in which case we give State an array of local states

```@example manual
state2 = State{Mixed}(system2, ["+", "2", "Occ"])
nothing # hide
```

Here we choose a mixed representation.

States may be added or multiplied by a number (they need to be based on the same system). For example

```@example manual
ghz = (State{Pure}(system1, "Up") + State{Pure}(system1, "Dn")) / sqrt(2)
nothing # hide
```

We can transform a pure representation into a mixed representation by

```@example manual
mixedstate = mix(state1)
nothing # hide
```

For mixed states there is a local mixed state `"FullyMixed"` which correspond to a density matrix proportional to the identity matrix (that is the infinite temperature state).

If you need a local state which is not predefined, it is possible to pass its vector (or matrix for mixed states) directly, For example, we could also define `state1` by

```@example manual
state1 = State{Pure}(system1, [1., 0.])
nothing # hide
```

## Limits

TMS uses Matrix Product State to internally represent quantum states. It is important to control the parameters of this approximation, in particular the maximum bond dimension and the cutoff on singular values. To achieve this, many functions accept a `Limits` object as keyword argument containing those parameters. It is built thus

```@example manual
lim = Limits(cutoff = 1e-10, maxdim = 50)
```

each (or both) of the arguments may be omitted in which case it corresponds to an absence of constraint for this parameter. In particular, `Limits()` represents no constraint.

To apply the constraints on a state, one uses

```@example manual
newstate = truncate(ghz; limits = lim)
nothing # hide
```

Many functions accept such an argument. For example, when adding states instead of

```@example manual
state = (state1 + ghz) / 2
nothing # hide
```

One can write

```@example manual
state = +(state1, ghz; limits = lim) / 2
nothing # hide
```

## Operators

In TMS, there are two kinds of operators: generic operators and indexed operators. For example,

```@example manual
X
```

represents the ``\sigma_x`` Pauli operator for qubits. This is a *generic operator*, it is not applied to a specific site.

```@example manual
X(3)
```

represents the ``\sigma_x`` Pauli operator applied to the system site number 3. This is an *indexed operator*.
 
Note that all predefined operator names start with a capital letter, so it is better to keep your own identifiers lowercase to prevent name collisions.

The operator system is very rich and flexible. For example, if you want to use this Hamiltonian

```math
H = \sum_{i=1}^{n-1} \sigma_x(i) \sigma_x(i+1)
```

you will simply write

```@example manual
n = 10
h = sum(X(i)X(i+1) for i in 1:n-1)
```

Many operations are defined on generic operators:

- addition, multiplication and power by a number
- tensor product: `X⊗X` is a two site operator (`⊗` is usually obtained by typing \otimes in your editor, just in case, one can also write `tensor(X, X)`)
- `dag` represents the adjoint operator, for example `C` is the `c` operator for fermions and `dag(C)` is ``c^\dagger``.
- `Dissipator` represents a Lindblad dissipator, for example `Dissipator(Sp)` is the jump operator that may flip a qubit toward up (`Sp` is the ``S^+`` operator)
- `Gate` represents an operator to be applied as a gate on a mixed state. It is useful to define noisy gate operators, for example `0.9Gate(Id) + 0.1Gate(X)` is a noisy gate operator that will apply an ``\sigma_x`` gate 10 percent of the time.
- `Proj` represents an operator that projects on the given state, for example `Proj("Up")` projects qubits on the up state.
- the functions `exp` and `sqrt`: for example `sqrt(Swap)`
- `controlled` for qubits makes controlled gates: `CX = controlled(X)`

For example one can define the Rxy 2-site operator by

```@example manual
Rxy(t) = exp(-im * t * (X⊗X + Y⊗Y) / 4)
```

If this is not enough to define your favorite operator you can create new ones by specifying their matrix.
The number in braces is the number of sites on which the operator must be applied.

```@example manual
myop = Operator{1}("MyOp", [1 1 ; 1 -1] / √2, involution_op)
```

```@example manual
myswap = Operator{2}("MySwap", [1 0 0 0 ; 0 0 1 0 ; 0 1 0 0 ; 0 0 0 1], involution_op)
```

Finally from generic operators, we define indexed operators by simply applying them to the corresponding sites

```@example manual
Rxy(0.2)(2, 5)
```

```@example manual
myop(3)
```

```@example manual
myswap(4, 7)
```

In the case of Hamiltonian or Lindbladian evolution the Hamiltonian part is to be multiplied by -im:

```julia
evolver = -im * hamiltonian + dissipators
```

## Algorithms

We can now work with states and operators.

We can apply gates with `apply`

```julia
newstate = apply(gates, oldstate; limits)
```

the `gates` argument is an indexed operator representing the gates to apply

the keyword argument `limits` fixes the constraints to apply

We can compute ground states with `dmrg`

```julia
energy, groundstate = dmrg(hamiltonian, startstate; options...)
```

the options are `limits` to set constraints and `nsweeps` to fix the number of sweeps among others.

We can do time evolution with `tdvp` and `approx_W`

```julia
newstate = tdvp(evolver, time, oldstate; options...)
newstate = approx_W(evolver, time, oldstate; options...)
```

the options are `limits` for the constraints, `nsweeps` for the number of steps to do and for `approx_W`, `order` and `w` for the parameters of the algorithm (`order = 4, w = 2` are usually good)

For more details, see the reference or the inline help.

## Measurements

Once we have created a state, we may want to measure it. Take for example a three qubit state

```@example manual
state = State{Pure}(System(3, Qubit()), "Up")
nothing # hide
```

```@example manual
result = measure(state, X(1)X(3))
```

will give ``\langle \psi | \sigma_x^1 \sigma_x^3 | \psi \rangle``

```@example manual
result = measure(state, X)
```

will give the array of the ``\langle \psi | \sigma_x^i | \psi \rangle``

```@example manual
result = measure(state, (X, Y))
```

will give the matrix of the ``\langle \psi | \sigma_x^i \sigma_y^j | \psi \rangle``

We can also measure other properties with
- `Trace` : the trace of the density matrix, this should be one, so it is a good indicator for accumulated error
- `TraceError`: measure the deviation from trace 1
- `Trace2`, `Purity`: measure the trace of the square of the density matrix
- `Hermiticity`: measure how well the density matrix is Hermitian, return 1 if Hermitian, 0 if anti-Hermitian
or any value in between
- `HermiticityError` measure the deviation from Hermiticity 1
- `Renyi2`: measure the Renyi entropy of order 2 of the system
- `SubRenyi2`: measure the Renyi entropy of order 2 of a subsystem
- `EE`: entanglement entropy for pure representation, OSEE for mixed
- `Linkdim`: the maximum bond dimension of the representation
- `MemoryUsage`: the memory used to store the representation

We can also ask for several measurements at the same time

```@example manual
results = measure(state, [X, X(2)Z(3), (X, Y), Trace, MemoryUsage])
```

For more details see the reference or the inline help.

## High Level Interface

### Framework

Most simulations follow the same pattern: start from some simple state, make some evolution and make measurements during or after the evolution and save the results to file. For these simple cases, TMS presents a simpler interface.

A simple simulation follows a single state through a certain number of phases which act in a simple way on the state and make measurements during and/or after the evolution and save the results to file.

The following phases are available:

- `CreateState` : create a simple state
- `GroundState` : compute the ground state using dmrg (requires a pure state)
- `ToMixed` : go from pure representation to mixed representation
- `Evolve` : do Hamiltonian or Lindbladian evolution
- `Gates` : apply some gates
- `PartialTrace` : trace the system over some sites (requires a mixed state)
- `SteadyState` : compute the steady state of a Lindblad equation (still experimental, requires a mixed state)
- `SaveState` : write the state to disk in a hdf5 file
- `LoadState` : read back a state written by `SaveState`

with these phases we define a `SimData` object that describes the simulation and finally, we call

```julia
runTMS(simdata)
```

which executes the simulation.

Phases that sweep take their measurements at every step by default. The
`measures_period` field of `Evolve`, `GroundState` and `SteadyState` raises that interval:
`measures_period = 10` measures one step out of ten, which is what long runs usually want.

The `phases` field is a list, but that list may contain lists, to any depth, and is
flattened before the simulation starts. This is meant for programs that build their phases
in pieces, a helper returning the several phases it needs rather than a single one, as
`create_graph_state` does.

### Example

What TMS is for is open systems, so here is one: six qubits evolving under a transverse
field Ising hamiltonian while each of them decays. `CreateState{Mixed}` is what makes the
state a density matrix, and the `Dissipator` terms added to the hamiltonian are what turn
the evolution into a Lindblad equation.

```julia
using TensorMixedStates, .Qubits

runTMS(SimData(
    name = "dissipative_ising",
    description = "six qubits under a transverse field Ising hamiltonian, each decaying at rate 0.2",
    phases = [
        CreateState{Mixed}(6, Qubit(), "Up"),
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

The magnetization on the six sites and the purity, at the start and at the end of the run:

```
Z        0.1    0.94098941    0.94117682   0.94117682   0.94117682   0.94117682    0.94098941
Purity   0.1    0.78980701
...
Z        2     -0.095182096  -0.19011196  -0.14498043  -0.14498043  -0.19011196   -0.095182096
Purity   2      0.1086794
```

The qubits start pure and pointing up; by the end the magnetization has reversed and the
purity has fallen to 0.11, which is a state no pure state code could have represented.

A longer one, the tight binding chain of fermions with dephasing noise:

```julia
using TensorMixedStates, .Fermions

hamiltonian(n) = -sum(dag(C)(i)C(i+1)+dag(C)(i+1)C(i) for i in 1:n-1)
dissipators(n, gamma) = sum(Dissipator(sqrt(4gamma) * N)(i) for i in 1:n)

sim_data(n, gamma, step) = SimData(
    name = "fermion_chain_with_dephasing",
    phases = [
        CreateState{Mixed}(n, Fermion(), [ iseven(i) ? "Occ" : "Emp" for i in 1:n ]),
        Evolve(
            duration = 4,
            time_step = step,
            algo = Tdvp(),
            evolver = -im*hamiltonian(n) + dissipators(n, gamma),
            limits = Limits(cutoff = 1e-30, maxdim = 100),
            measures = [
                "density.dat" => N,
                "OSEE.dat" => EE(div(n, 2))
            ]
        )
    ]
)

runTMS(sim_data(40, 1., 0.05))
```

### Output

`runTMS` creates a directory named after the `SimData` object `name` field and puts the output files there. In particular, it produces a `log` file showing the progression of the computation, a `prog.jl` file containing a copy of the script, a `description` file containing the content of the `SimData` `description` field, a `stamp` file containing version and date info, a `running` empty file is present during the computation, in case of error an empty `error` file is created.

Three keyword arguments may be given `restart` (default `false`) erases the directory before starting, `clean` (default `false`) erases the directory and does not run the simulation, `output` (default `nothing`) if set, does not create the directory nor any output files and redirect all output to the given io channel (useful values are stdout and devnull). 

### Long runs, checkpoints and stopping

A simulation meant to run for hours or days can save its progress, so that a crash, a
batch system killing the job, or a deliberate stop does not throw the computation away.
Two fields of `SimData` control it.

```julia
SimData(
    name = "my_simulation",
    checkpoint_interval = 600,       # seconds between two checkpoints
    max_time = 3.5 * 3600,           # stop cleanly after this long
    phases = [...],
)
```

`checkpoint_interval` is the time between two saves, `0` (the default) disables
checkpointing entirely. `max_time` is a wall clock budget: once it is past, the simulation
writes a checkpoint and returns instead of carrying on. Set it comfortably below the limit
of your batch job, since a checkpoint is only taken between two sweeps: a sweep that lasts
ten minutes delays the stop by up to ten minutes.

A checkpoint is written to `checkpoint.h5` and `checkpoint.json` in the simulation
directory. Both are written to temporary files and moved into place, so an interruption
during the save leaves the previous checkpoint intact.

#### Resuming

`runTMS` resumes on its own: run the same program again and it picks up where it left off,
skipping the phases that were finished and restarting the interrupted one at the sweep it
had reached. There is nothing to pass and nothing to change in the program. Running it
once more after the simulation completed does nothing.

Output files are cut back to the length they had at the checkpoint before the simulation
continues, so the measurements written between the last checkpoint and the interruption
are not duplicated. The result is the same file as an uninterrupted run would have
produced.

Use `restart = true` to ignore an existing checkpoint and start over, as it erases the
directory.

A checkpoint records which phases it belongs to, and `runTMS` refuses to resume one that
was written by a different simulation rather than mixing the two. So editing the phases of
a program and running it again under the same name reports an error instead of quietly
continuing something else, which matters when trying things out interactively. Give the
simulation another name, or pass `restart = true`.

Functions are the blind spot of that check: changing the coefficients of a time dependent
evolver, or the body of a `StateFunc`, leaves the phases looking the same, and the
simulation resumes from a checkpoint computed with the old ones. Restart such a run rather
than resume it.

#### Stopping on purpose

Three things ask a running simulation to stop, and all three write a checkpoint first:

- `max_time` running out,
- the file `stop` appearing in the simulation directory, typically with `touch
  my_simulation/stop`, from the shell or from the epilogue of a batch job,
- an interrupt, that is `Ctrl-C`.

The `stop` file is the one to reach for in batch, since it does not depend on how the
queueing system signals its jobs. It is removed when the simulation next starts, so it
never blocks a later run.

Note that a simulation writing to a directory turns `Ctrl-C` into a clean stop rather than
an immediate exit, for the whole program.

#### What can be resumed inside a phase

`Evolve`, `GroundState` and `SteadyState` are resumed at the sweep they reached. The other
phases are short enough to be replayed, and a checkpoint is taken between phases whenever
one is due.

The exception is `Gates`, which hands its whole list of gates to the tensor network library
in one go and therefore cannot be cut in the middle. A deep circuit is better written as
several `Gates` phases, which gives resume points at no cost.

### Measurements

Measurements are specified in the `measures` or `final_measures` fields. They take the form of a pair or list of pairs.

```julia
measures = destination => measurements
measures = [ dest1 => meas1, dest2 => meas2, ...]
```

The possible measurements are described in the measurements section of this manual. There are three types of destinations:

- filenames: writes the specified measurements to the given file as they are made. Special filenames are "stdout" (or "-"), "stderr", "" (for devnull)

  ```julia
  "file.dat" => X
  ```

- json filenames: filenames ending by ".json" are treated differently: data is accumulated during the simulation and written at the end in the JSON format.

  ```julia
  "file.json" => [Purity, X(2)Z(3), (X, Y)]
  ```

- Data object: data is accumulated during the simulation and stored in the `data` field of the `Simulation` object returned by `runTMS`. This is useful for analyzing the data inside the program.

  ```julia
  Data("mydata") => [TraceError, X(1), Y]
  ```

The `DataToFrame` function can be used on the result to get a `DataFrame` object (the `DataFrames` package must be imported first)

```julia
sim = runTMS(simdata)
df = DataToFrame(sim.data["mydata"])
```


For more information, see the reference or inline help for each phase, `SimData` and `runTMS`.
