# Manual

```@contents
Pages = ["manual.md"]
Depth = 3
```

## How to read this manual

The sections up to [Measurements](@ref manual-measurements) present the objects of TMS one after the other, sites,
states, operators, then the algorithms and the measurements, with small examples run directly
from Julia. [High Level Interface](@ref) shows how to write a whole simulation as a list of
phases, which writes its results to files and can be stopped and resumed; most programs are
written that way. [Threads and performance](@ref) is for when a computation gets large. The
reference pages that follow give every function in detail.

The examples assume TMS has been imported:

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

A state is either a wave function, its pure representation, written `Pure`, or a density
matrix, its mixed representation, written `Mixed`, which open systems need. Both are stored as
a matrix product state (MPS): a chain of tensors, one per site, joined by links whose size,
the bond dimension, measures how much entanglement, or for a density matrix how many
correlations, the state can hold. A density matrix is stored as the vector of its entries, so
each site has dimension ``d^2`` instead of ``d``, its norm is the Hilbert-Schmidt norm
``\sqrt{\mathrm{tr}\,\rho^\dagger\rho}`` rather than its trace, and it costs more than a wave
function of the same system.

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
exchange the state stays in the sector it was built in: a local state such as `[1, 1] / √2`,
which superposes two numbers of particles, is refused, and so is an operator of no definite charge.
For a mixed state, a quantity may also be conserved strongly, `Fermion(conserve = strong(N))`.
All of this is described in [Conserving a quantity](@ref).

## Limits

TMS uses Matrix Product States to represent quantum states internally. It is important to control the parameters of this approximation, in particular the maximum bond dimension and the cutoff, the total weight of the singular values a
truncation of the tensors may drop. To achieve this, many functions accept a `Limits` object as keyword argument containing those parameters. It is built thus

```@example manual
lim = Limits(cutoff = 1e-10, maxdim = 50)
```

Any of the arguments may be omitted in which case it corresponds to an absence of constraint for this parameter. In particular, `Limits()` represents no constraint.

A third parameter, `mindim`, sets the bond dimension the truncation is not allowed to go below, as in `Limits(cutoff = 1e-10, maxdim = 50, mindim = 4)`. `maxdim` takes precedence when the two conflict.

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
seventy lowercase names, among them `state`, `output`, `measure`, `trace`, `dim`,
`matrix`, `tensor`, `apply`, `norm` and `sample`. Assigning to one of them at the top level
of your program shadows the function for the rest of the file, and if you happen to have
used it before assigning to it, Julia 1.10 and 1.11 refuse the assignment outright with
`cannot assign a value to imported variable`. There is no such issue inside a function,
where `state = ...` is an ordinary local variable. This manual prefixes its own variables with
`my`, as in `mystate` and `myop`, which is one way of staying clear.

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
from an expression, a matrix or a function of the sites, see
[Operators and measurements of one's own](@ref):

```@example manual
myop = named([1 1 ; 1 -1] / √2, "MyOp")
```

Finally from generic operators, we define indexed operators by simply applying them to the corresponding sites

```@example manual
Rxy(0.2)(2, 5)
```

```@example manual
myop(3)
```

An operator of several sites that is a function of an operator, as the exponential in `Rxy`,
can only be applied as a gate, unless it is created with the sites it acts on, see
[Operators and measurements of one's own](@ref).

### Fermions

A fermionic operator is written as it is: `C(i)` destroys a fermion on site `i` and
`dag(C)(j)` creates one on site `j`, the Jordan-Wigner strings being inserted for you wherever
the operator is used, in expectation values, in MPOs and in gates. A hopping term is thus
simply

```@example manual
hopping = sum(dag(C)(i)C(i+1) + dag(C)(i+1)C(i) for i in 1:n-1)
nothing # hide
```

The strings follow the order of the sites: `C(j)` is the matrix of `C` on site `j` with `F`
on every site before it. A product of two fermionic operators placed in increasing order of
the sites, `A(i) * B(j)` with `i < j`, is thus `A * F` on site `i`, `F` on every site in
between and `B` on site `j`, and in the other order it takes the sign of the swap. A matrix
given for an operator of several sites, see [Operators and measurements of one's own](@ref),
is read in the basis of the sites as it is, with no string.

A function of an even fermionic operator of several sites, as the exponential of a hopping
term, `exp(-0.1im * (dag(C) ⊗ C + dag(dag(C) ⊗ C)))(2, 5)`, is applied as a gate with its
strings, on sites in any order and apart. A function of an odd operator mixes the two
parities and is refused.

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

`nsweeps` is the number of sweeps, passes of dmrg along the chain and back, and `limits` may
give one value per sweep, as in `Limits(maxdim = [10, 20, 50])`.

### Time evolution

Time evolution is done by `tdvp` or `approx_W`, which take an evolver, the time to evolve for
and the state. The evolver follows the convention

```julia
evolver = -im * hamiltonian + dissipators
```

dissipators being accepted on a mixed state only. On a pure state, it integrates
``\frac{d}{dt}|\psi\rangle = -iH|\psi\rangle``, and on a mixed state the Lindblad equation

```math
\frac{d\rho}{dt} = -i[H, \rho]
    + \sum_k \left(L_k \rho L_k^\dagger - \tfrac{1}{2}\{L_k^\dagger L_k, \rho\}\right)
```

each term `Dissipator(L)` giving one jump operator ``L_k``, with ``\hbar = 1``. A rate
``\gamma`` is written `γ * Dissipator(L)`, or equally `Dissipator(sqrt(γ) * L)`, since
`Dissipator(c * L)` is `abs2(c) * Dissipator(L)`. Under ``H = \sum_i \sigma_z^i``, qubits
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

`approx_W` approximates the exponential of the evolver over a step by MPOs, built from the
approximations WI or WII of Zaletel et al. It is called the same way, with in addition the
order of its approximation, from 1 to 4, and `w`, 1 or 2 for WI or WII, for which
`order = 4, w = 2` is usually a good choice:

```@example manual
mydecayed = approx_W(sum(Dissipator(Sm)(i) for i in 1:4), 1.0, myrho; order = 4, w = 2,
                     nsweeps = 10)
measure(mydecayed, Z)
```

An evolver may also depend on time, see [Time dependent evolvers](@ref).

### Steady states

`steady_state` looks for the state ``\rho`` with ``L\rho = 0``, ``L`` being the evolver,
starting from a mixed state, by running dmrg on ``L^\dagger L``. It returns two things: the
residual ``\|L\rho\|^2``, for ``\rho`` of unit Hilbert-Schmidt norm, which is zero at a true
steady state and tells how well the search converged, and the steady state itself, normalized
to trace one. Qubits decaying toward down at rate 1 and pumped toward up at rate
0.5 settle at ``\langle \sigma_z \rangle = -1/3``:

```@example manual
value, mysteady = steady_state(sum(Dissipator(Sm)(i) + 0.5Dissipator(Sp)(i) for i in 1:4),
                               myrho; nsweeps = 10)
measure(mysteady, Z)
```

For more details, see the [Algorithms](algorithms.md) page of the reference or the inline
help.

## [Measurements](@id manual-measurements)

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

Fermionic operators follow two rules. A single one given for every site, as `C`, is refused by
`measure` and `expect1`, since its value vanishes on every state of definite parity; placed on
a site, as `C(3)`, it is accepted, for a state that superposes parities. In a correlation, as
`(dag(C), C)`, the two operators must be both fermionic or both not.

We can also measure properties of the state as a whole, with state functions such as
`Trace`, `Purity`, `EntanglementEntropy(l)` or `Fidelity(ref)`: the
[Measurements](measurements.md) page has the table of them all.

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

- `CreateState`: create a simple state
- `GroundState`: compute the ground state using dmrg (requires a pure state)
- `ToMixed`: go from pure representation to mixed representation
- `Evolve`: do Hamiltonian or Lindbladian evolution
- `Gates`: apply some gates
- `PartialTrace`: trace the system over some sites (requires a mixed state)
- `Weaken`: conserve less, for instance a strong symmetry asked for weakly (see `weaken`)
- `SteadyState`: compute the steady state of a Lindblad equation (requires a mixed state)
- `SaveState`: write the state to disk in an HDF5 file
- `LoadState`: read back a state written by `SaveState`

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

TMS is made for open systems, so here is one: six qubits evolving under a transverse
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

The magnetization on the six sites and the purity after the first time step and at the end
of the run, `Evolve` measuring after each step:

```
Z        0.1    0.94098941    0.94117682   0.94117682   0.94117682   0.94117682    0.94098941
Purity   0.1    0.78980701
...
Z        2     -0.095182096  -0.19011196  -0.14498043  -0.14498043  -0.19011196   -0.095182096
Purity   2      0.1086794
```

The qubits start pure and pointing up; by the end the magnetization has reversed and the
purity has fallen to 0.11: a mixed state, which no single wave function can represent.

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

`runTMS` creates a directory named after the `name` field of the `SimData` object and runs
the phases in it: a relative file name, of a destination or of `SaveState` and `LoadState`, is
taken there. The directory holds in particular:

- `log`: the progression of the computation;
- `prog.jl`: a copy of the script;
- `description`: the content of the `description` field of the `SimData` object;
- `stamp`: the versions, the date, the BLAS library and the thread settings of the run;
- `running`: an empty file present during the computation;
- `error`: an empty file created in case of error.

`runTMS` takes three keyword arguments:

- `restart` (default `false`): erase the directory before starting;
- `clean` (default `false`): erase the directory and do not run the simulation;
- `output` (default `nothing`): if set, create neither the directory nor the files `runTMS`
  writes in it, and redirect all output to the given stream, `stdout` or `devnull` for
  instance. A `SaveState` still writes its file, in the current directory.

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

`checkpoint_interval` is the time between two saves; `0`, the default, means no periodic save.
Even then, a stop, from `max_time`, the `stop` file or an interrupt, writes a checkpoint, so
that the simulation can always be resumed. `max_time` is a wall clock budget: once it is past,
the simulation writes a checkpoint and returns instead of carrying on. Set it comfortably below the limit
of your batch job, since a checkpoint is only taken between two sweeps, two steps or two
phases: a sweep that lasts ten minutes delays the stop by up to ten minutes.

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
had reached. There is nothing to pass and nothing to change in the program. A checkpoint is
also written after the last phase, when periodic checkpoints are on or when one is already on
the disk, left by a stop, so that running the program once more after the simulation
completed does nothing. A simulation that never wrote one runs again from the start.

Output files are cut back to the length they had at the checkpoint before the simulation
continues, so the measurements written between the last checkpoint and the interruption
are not duplicated. The result is the same file as an uninterrupted run would have
produced. The log is the exception: it keeps the history of every run, so what an
interrupted run wrote after its last checkpoint stays, with the line saying why it stopped,
followed by the line marking the resume and by the steps done again from the checkpoint.

Use `restart = true` to ignore an existing checkpoint and start over, as it erases the
directory.

The random number generator is not part of a checkpoint, so after a resume the random draws,
of a `CreateState` with `randomize` or of a measurement calling `sample`, differ from those of
the uninterrupted run: the results are as valid, but they are not the same numbers. A
`CreateState` that the resume runs again applies its `seed` again and draws the same state,
while one finished before the resume point does not, and samples measured during an evolution
are not reproduced.

A checkpoint records which phases it belongs to, and `runTMS` refuses to resume one that
was written by a different simulation rather than mixing the two. So editing the phases of
a program and running it again under the same name reports an error instead of quietly
continuing something else, which matters when trying things out interactively. Give the
simulation another name, or pass `restart = true`.

Functions are the blind spot of that check: changing the coefficients of a time dependent
evolver, or the body of a `StateFunc`, leaves the phases looking the same, and the
simulation resumes from a checkpoint computed with the old ones. Restart such a run rather
than resume it. The other way round, a function that captures another function, as
`t -> g(t)` built inside a function taking `g`, looks different when the program is included
again in the same Julia session, and its checkpoint is then refused as another simulation's:
resume it from a new session.

#### Stopping on purpose

Three things ask a running simulation to stop, and all three write a checkpoint first:

- `max_time` running out,
- the file `stop` appearing in the simulation directory, typically with `touch
  my_simulation/stop`, from the shell or from the epilogue of a batch job,
- an interrupt, that is `Ctrl-C`.

The `stop` file is the one to reach for in batch, since it does not depend on how the
queueing system signals its jobs. It is removed when the simulation next starts, so it
never blocks a later run.

Note that while a simulation writing to a directory runs, `Ctrl-C` stops it cleanly rather
than ending the program. When `runTMS` returns, `Ctrl-C` gets back the behaviour Julia gives
it by default.

#### What can be resumed inside a phase

`Evolve`, `GroundState` and `SteadyState` are resumed at the sweep they reached, and a phase
of your own written with `run_steps` at the step it reached, see [Phases of one's own](@ref own-phases).
The other phases are short enough to be replayed, and a checkpoint is taken between phases
whenever one is due.

The exception is `Gates`, which hands its whole list of gates to the tensor network library
in one go and therefore cannot be cut in the middle. A deep circuit is better written as
several `Gates` phases, which gives resume points at no cost.

### Measurements

Measurements are specified in the `measures` field of `Evolve`, `GroundState` and
`SteadyState`, taken as the phase sweeps, and in the `final_measures` field of every phase and
of `SimData`, taken at its end. They take the form of a pair or list of pairs.

```julia
measures = destination => measurements
measures = [ dest1 => meas1, dest2 => meas2, ...]
```

The possible measurements are described on the [Measurements](measurements.md) page. A
destination is a file, a json file or a `Data` object, which keeps the values in the program;
the destinations and the format of what they hold are described in
[Output](@ref measure-output).

For more information, see the reference or inline help for each phase, `SimData` and `runTMS`.

### Phases, algorithms and representations of one's own

Phases, time evolution algorithms and representations of a state of your own are described
in [Extending TMS](@ref).

## Threads and performance

Most of the running time goes into the tensor contractions of ITensors. How they use the cores
of the machine comes down to three choices.

**Starting Julia.** Start it with as many threads as the machine has:

```sh
julia --threads=auto my_script.jl
```

or `julia -t auto my_script.jl` for short; `--threads=4`, or `-t 4`, gives exactly four. The
garbage collector then runs on as many threads from Julia 1.12 on, and on half as many
before, which speeds up the runs, and they are what the `:blocks` mode below runs on. Started
without them, Julia still runs the products of dense matrices on several cores, BLAS having
threads of its own, but everything else on a single one. The threads of Julia are fixed when
it starts: they cannot be added from within a program.

**The BLAS library.** On `x86_64` machines running Linux or Windows, MKL, which replaces
OpenBLAS for the products of dense matrices, is often faster, see [BLAS backend](@ref). This
choice is independent of the two others.

**The mode.** The `threading` field of `SimData`, or `set_threading`, chooses how the
contractions use the threads of Julia:

- `:dense`: BLAS runs each product of matrices on several threads, and Strided, which ITensors
  uses for the permutations of dense tensors, runs on a single one, as ITensors recommends.
  Julia starts Strided on as many threads as it has itself, where they compete with those of
  BLAS. `SimData` applies this mode by default; a program calling the functions of TMS directly
  should start with `set_threading(:dense)`;
- `:blocks` asks ITensors to run the products of the blocks of block sparse tensors in parallel
  on the threads of Julia, BLAS running each of them on a single thread. The tensors of a
  system are block sparse when it conserves something, see [Conserving a quantity](@ref), so
  this mode is meant for such systems: it is worth trying on them, with OpenBLAS above all;
- `:auto`, for `SimData` only, chooses before each phase: `:blocks` when the system of the
  state conserves something and Julia has several threads, `:dense` otherwise and until there
  is a state.
  `set_threading(mysystem)` makes the same choice once.

```julia
using MKL                      # if MKL is installed, on x86_64 Linux or Windows: before TMS
using TensorMixedStates
set_threading(:dense)          # first thing in a program calling the functions of TMS directly

old = set_threading(:blocks)   # set_threading returns the threading it replaces,
set_threading(old)             # which it takes back

SimData(name = "my_simulation", threading = :blocks, phases = [...])
```

The threading `set_threading` returns is a `ThreadingState`, which gives the threads of BLAS and
of Strided and whether block sparse multithreading is on. One can be built by hand to try other
settings: `set_threading(ThreadingState(blas = 2))` puts BLAS on two threads and leaves the rest
as it is.

The `stamp` file of a simulation records all these settings, so that the running times of two
runs can be compared.

On a laptop with four cores and an Intel processor, for ground state searches and time
evolutions of chains at bond dimensions from 64 to 768, the running times changed as follows:

| | dense tensors | block sparse tensors |
|:---|:---|:---|
| `--threads=auto` with `:dense`, rather than `julia` alone | −2 to −11 % | −8 to −14 % |
| MKL rather than OpenBLAS | −6 to −21 % | −5 to −14 % |
| Strided left on the threads of Julia, rather than on one as `:dense` puts it | +1 to +61 % | no effect |
| `:blocks` rather than `:dense`, with OpenBLAS | +20 to +65 % | −19 to +1 % |
| `:blocks` rather than `:dense`, with MKL | +35 to +67 % | −1 to +31 % |

Strided weighs most at moderate bond dimensions, and hardly at all from 512 on.

The documentations of Julia and of ITensors agree that the way to find the best settings is to
try them on a few sweeps of your own calculation.

### BLAS backend

Most of the running time of TMS is spent in the tensor contractions of ITensors, which
themselves call BLAS. Julia ships with OpenBLAS and TMS uses it as it comes: switching the
BLAS backend affects the whole Julia session, so that choice is left to you rather than
made by a library you load.

On `x86_64` machines running Linux or Windows, Intel's MKL is often noticeably faster on
these contractions. To use it, add `MKL` to your project and load it before TMS:

```julia
using MKL
using TensorMixedStates
```

MKL is not distributed for macOS nor for ARM machines, where OpenBLAS is the only option.
