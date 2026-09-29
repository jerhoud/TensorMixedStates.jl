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

To define a site, call the corresponding constructor, for example

```@example manual
s = Qubit()
```

Four site constructors need an argument: `Qudit(dim)` and `Boson(dim)` for the dimension of
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

gives you a system with 10 qubits. Systems may have sites of different types, in which case `System` is given an array of sites

```@example manual
using .Fermions

system2 = System([Qubit(), Boson(4), Fermion()])
nothing # hide
```

gives you a three site system.

## States

States may be in pure or mixed representation, these two possibilities are represented in TMS, by `Pure` or `Mixed`.

To create a state, we call the `State` constructor

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

For mixed states there is a local mixed state `"FullyMixed"` which corresponds to a density matrix proportional to the identity matrix (that is the infinite temperature state).

If you need a local state which is not predefined, it is possible to pass its vector (or matrix for mixed states) directly. For example, we could also define `state1` by

```@example manual
state1 = State{Pure}(system1, [1., 0.])
nothing # hide
```

## Conserved quantities

A site can be told that a quantity is conserved, for instance the number of fermions:

```@example manual
mysystem = System(10, Fermion(conserve = N))
mystate = State{Pure}(mysystem, [isodd(i) ? "Occ" : "Emp" for i in 1:10])
nothing # hide
```

The tensors are then block sparse, which makes a large computation smaller and faster, and in
exchange the state stays in the sector it was built in: a local state such as `"+"`, which
superposes two numbers of particles, is refused, and so is an operator of no definite charge.
For a mixed state, a quantity may also be conserved strongly, `Fermion(conserve = strong(N))`.
All of this is described in [Conserving a quantity](@ref).

## Limits

TMS uses Matrix Product States to represent quantum states internally. It is important to control the parameters of this approximation, in particular the maximum bond dimension and the cutoff on singular values. To achieve this, many functions accept a `Limits` object as keyword argument containing those parameters. It is built thus

```@example manual
lim = Limits(cutoff = 1e-10, maxdim = 50)
```

Any of the arguments may be omitted in which case it corresponds to an absence of constraint for this parameter. In particular, `Limits()` represents no constraint.

A third parameter, `mindim`, sets the bond dimension the truncation is not allowed to go below, as in `Limits(cutoff = 1e-10, maxdim = 50, mindim = 4)`. `maxdim` keeps the last word when the two ask for opposite things.

To apply the constraints on a state, one uses

```@example manual
newstate = truncate(ghz; limits = lim)
nothing # hide
```

Many functions accept such an argument. For example, when adding states instead of

```@example manual
state3 = (state1 + ghz) / 2
nothing # hide
```

one can write

```@example manual
state3 = +(state1, ghz; limits = lim) / 2
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
 
Note that all predefined operator names start with a capital letter, so it is better to keep your own identifiers lowercase to prevent name collisions with them.

Lowercase is not a safe harbour by itself, though. `using TensorMixedStates` also brings in
some sixty lowercase names, among them `state`, `output`, `measure`, `trace`, `dim`,
`matrix`, `tensor`, `apply`, `norm` and `sample`. Assigning to one of them at the top level
of your program shadows the function for the rest of the file, and if you happen to have
used it before assigning to it, Julia 1.10 and 1.11 refuse the assignment outright with
`cannot assign a value to imported variable`. Inside a function there is no such issue,
where `state = ...` is an ordinary local variable. This manual prefixes its own variables
with `my`, as in `mystate` and `myop`, which is one way of staying clear.

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
- `Dissipator` represents a Lindblad dissipator, for example `Dissipator(Sp)` is the dissipator whose jump operator, `Sp`, the ``S^+`` operator, flips a qubit toward up
- `Gate` represents an operator to be applied as a gate on a mixed state. It is useful to define noisy gate operators, for example `0.9Gate(Id) + 0.1Gate(X)` is a noisy gate operator that will apply a ``\sigma_x`` gate 10 percent of the time.
- `Proj` represents an operator that projects on the given state, for example `Proj("Up")` projects qubits on the up state.
- the functions `exp` and `sqrt`: for example `sqrt(Swap)`
- `controlled` for qubits makes controlled gates: `CX = controlled(X)`

For example one can define the Rxy 2-site operator by

```@example manual
Rxy(t) = exp(-im * t * (X⊗X + Y⊗Y) / 4)
```

If this is not enough to define your favorite operator you can create new ones with `named`,
from an expression, a matrix or a function of the sites. The type of the operator, which
`simplify` reasons with (see `OpType`), is read off its matrix: here an involution.

```@example manual
myop = named([1 1 ; 1 -1] / √2, "MyOp")
```

```@example manual
myop.type
```

A matrix of several sites given alone does not say how many sites it acts on, which is then
given in braces, with the type:

```@example manual
myswap = Operator{2}("MySwap", [1 0 0 0 ; 0 0 1 0 ; 0 1 0 0 ; 0 0 0 1], involution_op)
```

Whichever way it is given, the type is checked against the matrix each time the operator is
placed on a site.

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

An operator of several sites defined by a matrix, as `myswap` is, or by a function of such an
operator, as the exponential in `Rxy`, can only be applied as a gate. To measure it or to put it in a hamiltonian, give
the sites it acts on when creating it, one per index or a single one for identical sites,
whose number is then read off the size of the matrix:

```@example manual
myswap2 = named([1 0 0 0 ; 0 0 1 0 ; 0 1 0 0 ; 0 0 0 1], "MySwap2", Qubit())
```

It is then split into a sum of products of one site operators, which becomes its
definition, the way `Swap` is defined by an expression:

```@example manual
myswap2.expr
```

The factors carry a definite charge of what their sites conserve. On a fermionic site the
matrix has to commute with `F`, since it is taken as it is, with no Jordan-Wigner string:
an operator moving fermions between sites is written with `C` and `dag(C)` instead.

In the case of Hamiltonian or Lindbladian evolution the Hamiltonian part is to be multiplied by -im:

```julia
evolver = -im * hamiltonian + dissipators
```

### Fermions

A fermionic operator is written as it is: `C(i)` destroys a fermion on site `i` and
`dag(C)(j)` creates one on site `j`, the Jordan-Wigner strings being inserted for you wherever
the operator is used, in expectation values, in MPOs and in gates. A hopping term is thus
simply

```@example manual
hopping = sum(dag(C)(i)C(i+1) + dag(C)(i+1)C(i) for i in 1:n-1)
nothing # hide
```

The expectation value of a single fermionic operator vanishes on any state of definite parity,
so `measure` refuses to measure a fermionic operator on every site, as `C`, and so does
`expect1`; placed on a site, `C(3)`, it is accepted, for a state that superposes parities. The
two operators of a correlation, as `(dag(C), C)`, must both be fermionic or both not.

## Algorithms

We can now work with states and operators.

### Gates

Gates are applied with `apply`. In a product of gates the rightmost acts first: here the
Hadamard gate on the first qubit, then the controlled not, which turns two up qubits into a
Bell pair

```@example manual
mybell = apply(controlled(X)(1, 2) * H(1), State{Pure}(System(2, Qubit()), "Up"))
measure(mybell, [Z(1)Z(2), X(1)X(2)])
```

Applying all the gates in a single call is much more efficient than one by one, and the
keyword argument `limits` constrains the truncations made along the way.

### Ground states

The ground state of a hamiltonian is computed by `dmrg`, from a starting state, here a random
state of bond dimension 8. It returns the energy and the ground state, which can then be
measured. For the Ising chain in a transverse field

```@example manual
mysystem = System(10, Qubit())
hamiltonian = -sum(Z(i)Z(i + 1) for i in 1:9) - sum(X(i) for i in 1:10)
energy, ground = dmrg(hamiltonian, RandomState{Pure}(mysystem, 8); nsweeps = 10,
                      limits = Limits(maxdim = 20))
energy
```

```@example manual
measure(ground, X)
```

`nsweeps` is the number of sweeps, and `limits` may give one value per sweep, as in
`Limits(maxdim = [10, 20, 50])`.

### Time evolution

Time evolution is done by `tdvp` or `approx_W`, which take an evolver, the time to evolve for
and the state. The evolver follows the convention

```julia
evolver = -im * hamiltonian + dissipators
```

dissipators being accepted on a mixed state only. Under ``H = \sum_i \sigma_z^i``, qubits
starting in `"+"` precess, with ``\langle \sigma_x \rangle = \cos 2t``, which is
``\cos 1 \approx 0.5403`` at ``t = 0.5``:

```@example manual
myevolved = tdvp(-im * sum(Z(i) for i in 1:4), 0.5, State{Pure}(System(4, Qubit()), "+");
                 nsweeps = 10)
measure(myevolved, X)
```

`nsweeps` is the number of steps, each of length `t / nsweeps`, and `limits` constrains the
state as before. With no hamiltonian and `Dissipator(Sm)` on each site, the qubits of a mixed
state decay from up to down, with ``\langle \sigma_z \rangle = 2e^{-t} - 1``, which is
``-0.2642`` at ``t = 1``:

```@example manual
myrho = State{Mixed}(System(4, Qubit()), "Up")
mydecayed = tdvp(sum(Dissipator(Sm)(i) for i in 1:4), 1.0, myrho; nsweeps = 10)
measure(mydecayed, Z)
```

`approx_W` is called the same way, with in addition the order of its approximation, from 1
to 4, and `w`, 1 or 2, for which `order = 4, w = 2` is usually a good choice:

```@example manual
mydecayed = approx_W(sum(Dissipator(Sm)(i) for i in 1:4), 1.0, myrho; order = 4, w = 2,
                     nsweeps = 10)
measure(mydecayed, Z)
```

An evolver may also depend on time, see [Time dependent evolvers](@ref).

### Steady states

The steady state of an open system is computed by `steady_state`, from a mixed state to start
from. It returns a value, zero for a steady state, and the state. Qubits decaying toward down
at rate 1 and pumped toward up at rate 0.5 settle at ``\langle \sigma_z \rangle = -1/3``:

```@example manual
value, mysteady = steady_state(sum(Dissipator(Sm)(i) + 0.5Dissipator(Sp)(i) for i in 1:4),
                               myrho; nsweeps = 10)
measure(mysteady, Z)
```

For more details, see the [Algorithms](algorithms.md) page of the reference or the inline
help.

## Measurements

Once we have created a state, we may want to measure it. Take for example a three qubit state

```@example manual
mystate = State{Pure}(System(3, Qubit()), ["+", "Up", "+"])
nothing # hide
```

```@example manual
result = measure(mystate, X(1)X(3))
```

will give ``\langle \psi | \sigma_x^1 \sigma_x^3 | \psi \rangle``

```@example manual
result = measure(mystate, X)
```

will give the array of the ``\langle \psi | \sigma_x^i | \psi \rangle``

```@example manual
result = measure(mystate, (X, Y))
```

will give the matrix of the ``\langle \psi | \sigma_x^i \sigma_y^j | \psi \rangle``, complex
since its diagonal is ``\langle \sigma_x \sigma_y \rangle = i \langle \sigma_z \rangle``

We can also measure other properties with
- `Trace`: the trace of the density matrix, which should be one, so it is a good indicator of
  accumulated error
- `TraceError`: the deviation of the trace from 1
- `Trace2`, `Purity`: the trace of the square of the density matrix
- `Norm`: the norm of the state
- `Hermiticity`: how Hermitian the density matrix is, 1 if Hermitian, 0 if anti-Hermitian, or
  any value in between
- `HermiticityError`: the deviation of the hermiticity from 1
- `Renyi2`: the Rényi entropy of order 2 of the system
- `SubRenyi2`: the Rényi entropy of order 2 of a subsystem
- `MutualInfoRenyi2`: the Rényi-2 mutual information between a subsystem and the rest
- `EntanglementEntropy`: the entanglement entropy for a pure representation, the OSEE for a
  mixed one
- `Fidelity`, `Overlap`: the fidelity with a reference state, and the inner product with it
- `Variance`: the variance of the energy of a hamiltonian
- `MaxLinkdim`: the maximum bond dimension of the representation
- `MemoryUsage`: the memory the state occupies, including the caches of the measurements
  already made on it

We can also ask for several measurements at the same time

```@example manual
results = measure(mystate, [X, X(2)Z(3), (X, Y), Trace, MemoryUsage])
```

For more details see the reference or the inline help.

## High Level Interface

### Framework

Most simulations follow the same pattern: start from a simple state, evolve it, measure it
during or after the evolution and save the results to files. For these cases TMS offers a
higher level interface, in which a simulation follows a single state through a sequence of
phases, each acting on the state and making its measurements.

The following phases are available:

- `CreateState` : create a simple state
- `GroundState` : compute the ground state using dmrg (requires a pure state)
- `ToMixed` : go from pure representation to mixed representation
- `Evolve` : do Hamiltonian or Lindbladian evolution
- `Gates` : apply some gates
- `PartialTrace` : trace the system over some sites (requires a mixed state)
- `Weaken` : conserve less, for instance a strong symmetry asked for weakly (see `weaken`)
- `SteadyState` : compute the steady state of a Lindblad equation (requires a mixed state)
- `SaveState` : write the state to disk in an HDF5 file
- `LoadState` : read back a state written by `SaveState`

With these phases we define a `SimData` object that describes the simulation and finally, we call

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

### Examples

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

A longer one, the tight binding chain of fermions with dephasing noise. The hopping and the
dephasing both conserve the number of fermions, so the sites are told to conserve it,
`Fermion(conserve = N)`, and the tensors become block sparse, see [What it saves](@ref).
Since the jump operators, `N`, commute with that number, it could also be conserved
strongly, see [Weak and strong symmetries](@ref).

```julia
using TensorMixedStates, .Fermions

hamiltonian(n) = -sum(dag(C)(i)C(i+1)+dag(C)(i+1)C(i) for i in 1:n-1)
dissipators(n, gamma) = sum(Dissipator(sqrt(4gamma) * N)(i) for i in 1:n)

sim_data(n, gamma, step) = SimData(
    name = "fermion_chain_with_dephasing",
    phases = [
        CreateState{Mixed}(n, Fermion(conserve = N), [ iseven(i) ? "Occ" : "Emp" for i in 1:n ]),
        Evolve(
            duration = 4,
            time_step = step,
            algo = Tdvp(),
            evolver = -im*hamiltonian(n) + dissipators(n, gamma),
            limits = Limits(cutoff = 1e-30, maxdim = 100),
            measures = [
                "density.dat" => N,
                "OSEE.dat" => EntanglementEntropy(div(n, 2))
            ]
        )
    ]
)

runTMS(sim_data(40, 1., 0.05))
```

### Output

`runTMS` creates a directory named after the `name` field of the `SimData` object and puts
the output files there, in particular:

- `log`: the progression of the computation;
- `prog.jl`: a copy of the script;
- `description`: the content of the `description` field of the `SimData` object;
- `stamp`: the versions, the date, the BLAS library and the thread settings of the run;
- `running`: an empty file present during the computation;
- `error`: an empty file created in case of error.

Three keyword arguments may be given:

- `restart` (default `false`): erase the directory before starting;
- `clean` (default `false`): erase the directory and do not run the simulation;
- `output` (default `nothing`): if set, create neither the directory nor any output file, and
  redirect all output to the given stream, `stdout` or `devnull` for instance.

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

`checkpoint_interval` is the time between two saves, `0` (the default) disables the periodic
checkpoints: a stop, from `max_time` or the `stop` file, or an interrupt, still writes one, so
that the simulation can be resumed. `max_time` is a wall clock budget: once it is past, the simulation
writes a checkpoint and returns instead of carrying on. Set it comfortably below the limit
of your batch job, since a checkpoint is only taken between two sweeps: a sweep that lasts
ten minutes delays the stop by up to ten minutes.

A checkpoint is written in the simulation directory, the state to `checkpoint-1.h5` or
`checkpoint-2.h5` and the rest to `checkpoint.json`, which names the state file. The state
goes to the file the previous checkpoint does not use, and `checkpoint.json` is moved into
place last, so an interruption at any point of the save leaves the previous checkpoint
intact.

Without a directory, when `runTMS` is given `output`, nothing is saved: an interrupt goes on
to the caller, and `max_time` stops the simulation with a message saying that it cannot be
resumed.

#### Resuming

`runTMS` resumes on its own: run the same program again and it picks up where it left off,
skipping the phases that were finished and restarting the interrupted one at the sweep it
had reached. There is nothing to pass and nothing to change in the program. With periodic
checkpoints on, one is written after the last phase, so that running it once more after the
simulation completed does nothing.

Output files are cut back to the length they had at the checkpoint before the simulation
continues, so the measurements written between the last checkpoint and the interruption
are not duplicated. The result is the same file as an uninterrupted run would have
produced.

Use `restart = true` to ignore an existing checkpoint and start over, as it erases the
directory.

The random number generator is not part of a checkpoint. What draws random numbers after the
point a run resumes from, a `CreateState` with `randomize` or a measurement calling `sample`,
draws other numbers than the uninterrupted run, and the `seed` of a `CreateState` finished
before that point is not applied again. The results are as valid, but they are not the same
numbers. A `CreateState` that draws its state after its own `seed` is reproduced, since it is
replayed whole; samples measured during an evolution are not.

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

- filenames: writes the specified measurements to the given file as they are made. Special filenames are "stdout" (or "-"), "stderr", "" (for devnull). The files `runTMS` writes itself in the simulation directory cannot be destinations: `log`, `stop`, `error`, `running`, `stamp`, `description`, `prog.jl` and the checkpoint files

  ```julia
  "file.dat" => X
  ```

- json filenames: filenames ending in ".json" are treated differently: data is accumulated during the simulation and written at the end in the JSON format.

  ```julia
  "file.json" => [Purity, X(2)Z(3), (X, Y)]
  ```

- Data object: data is accumulated during the simulation and stored in the `data` field of the `Simulation` object returned by `runTMS`. This is useful for analyzing the data inside the program.

  ```julia
  Data("mydata") => [TraceError, X(1), Y]
  ```

A complex value takes two columns in a file, its real part then its imaginary part, a json
file writes it as `{"re": …, "im": …}`, and a `Data` object holds it as a complex number. Which
values are complex is described in the measurements section.

A json file and a `Data` object hold, for each measurement, the lists `"times"`, `"data"` and
`"events"`: the time of each value, the value, and the number of the measurement set it
belongs to, counted for each destination, one per time the destination was written. Values
measured together share their event, where their time alone repeats over the sweeps of a
ground state search, over a circuit, or once the time is set back. A matrix is written in a
json file as the list of its rows, as a file writes it row by row.

The `data_to_frame` function can be used on the result to get a `DataFrame` object, with one row per event (the `DataFrames` package must be imported first)

```julia
mysim = runTMS(simdata)
df = data_to_frame(mysim.data["mydata"])
```


For more information, see the reference or inline help for each phase, `SimData` and `runTMS`.
