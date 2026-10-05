# High Level Interface

## Framework

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
- `Thermalize`: take a mixed state towards the thermal state of a hamiltonian, by evolution in
  imaginary time
- `SaveState`: write the state to disk in an HDF5 file
- `LoadState`: read back a state written by `SaveState`

With these phases we define a `SimData` object that describes the simulation and finally, we call

```julia
runTMS(simdata)
```

which executes the simulation.

Phases that sweep take their measurements at every step by default. The
`measurements_period` field of `Evolve`, `GroundState`, `SteadyState` and `Thermalize` raises that
interval:
`measurements_period = 10` measures one step out of ten, which is what long runs usually want.

The `phases` field is a list, but that list may contain lists, to any depth, and is
flattened before the simulation starts. This is meant for programs that build their phases
in pieces, a helper returning the several phases it needs rather than a single one, as
`create_graph_state` does.

## Examples

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
            measurements = "data" => [Z, Purity],
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
            measurements = [
                "density.dat" => N,
                "OSEE.dat" => EntanglementEntropy(div(n, 2))
            ]
        )
    ]
)

runTMS(sim_data(40, 1., 0.05))
```

## Output

`runTMS` creates a directory named after the `name` field of the `SimData` object and runs
the phases in it: a relative file name, of a destination or of `SaveState` and `LoadState`, is
taken there. The directory holds in particular:

- `log`: the progression of the computation;
- `prog.jl`: a copy of the script;
- `prog_args.json`: the command line arguments it was given, if any;
- `description`: the content of the `description` field of the `SimData` object;
- `stamp`: the versions, the date, the BLAS library and the thread settings of the run;
- `running`: present during the computation, naming the machine and the process running it.
  `runTMS` refuses to start in a directory where another run may be going on. A file left by
  a process of the same machine that no longer exists, killed for instance, is replaced with
  a warning in the log; any other is removed by hand;
- `error`: an empty file created in case of error.

`runTMS` takes three keyword arguments:

- `restart` (default `false`): erase the directory before starting;
- `clean` (default `false`): erase the directory and do not run the simulation;
- `output` (default `nothing`): if set, create neither the directory nor the files `runTMS`
  writes in it, and redirect all output to the given stream, `stdout` or `devnull` for
  instance. A `SaveState` still writes its file, in the current directory.

## Long runs, checkpoints and stopping

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

### Resuming

`runTMS` resumes on its own: run the same program again and it picks up where it left off,
skipping the phases that were finished and restarting the interrupted one at the sweep it
had reached. There is nothing to pass and nothing to change in the program. A checkpoint is
also written after the last phase, when periodic checkpoints are on or when one is already on
the disk, left by a stop, so that running the program once more after the simulation
completed does nothing. A simulation that never wrote one runs again from the start.

Output files are cut back to the length they had at the checkpoint before the simulation
continues, so the measurements written between the last checkpoint and the interruption
are not duplicated. The result is the same file as an uninterrupted run would have
produced. A json file and a `Data` object are put back as they were at the checkpoint, the
values measured after it dropped. The log is the exception: it keeps the history of every run, so what an
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

A checkpoint is resumed only by the program that wrote it: the program run must be the same,
byte for byte, as the copy `prog.jl` in the directory, and be given the same command line
arguments, which `prog_args.json` holds when there are some. Editing the program and running it again under the
same name thus reports an error instead of quietly continuing something else. To resume with
the edited program all the same, more time or other threads for instance, copy it onto
`prog.jl`, and its arguments into `prog_args.json`; to start over, pass `restart = true` or
give the simulation another name. Only these two files are compared: a change in a file the
program includes, or in a parameter it reads elsewhere, from an environment variable or a data
file, goes unseen. A simulation run from the REPL, with no program file, is always resumed.

### Stopping on purpose

Three things ask a running simulation to stop, and all three write a checkpoint first:

- `max_time` running out,
- the file `stop` appearing in the simulation directory, typically with `touch
  my_simulation/stop` from the shell; in a batch job, `max_time` set a little below the time
  limit does the same without a file,
- an interrupt, that is `Ctrl-C`.

`max_time` is the one to reach for in batch, since it does not depend on how the queueing
system signals its jobs. The `stop` file is removed when the simulation next starts, so it
never blocks a later run.

A stopped run returns the simulation it resumes from, and [`stopped`](@ref) tells it from one
that completed:

```julia
sim = runTMS(sim_data)
if stopped(sim)
    println("stopped at time ", sim.time, ", run the program again to resume")
end
```

Note that while a simulation writing to a directory runs, `Ctrl-C` stops it cleanly rather
than ending the program. When `runTMS` returns, `Ctrl-C` gets back the behaviour Julia gives
it by default.

### What can be resumed inside a phase

`Evolve`, `GroundState`, `SteadyState` and `Thermalize` are resumed at the sweep they reached, and a phase
of your own written with `run_steps` at the step it reached, see [Phases of one's own](@ref own-phases).
The other phases are short enough to be replayed, and a checkpoint is taken between phases
whenever one is due.

The exception is `Gates`, which hands its whole list of gates to the tensor network library
in one go and therefore cannot be cut in the middle. A deep circuit is better written as
several `Gates` phases, which gives resume points at no cost.

## Measurements

Measurements are specified in the `measurements` field of `Evolve`, `GroundState`,
`SteadyState` and `Thermalize`, taken as the phase sweeps, and in the `final_measurements` field of every phase and
of `SimData`, taken at its end. They take the form of a pair or list of pairs.

```julia
measurements = destination => measurements
measurements = [ dest1 => meas1, dest2 => meas2, ...]
```

The possible measurements are described on the [Measurements](measurements.md) page. A
destination is a file, a json file or a `Data` object, which keeps the values in the program;
the destinations and the format of what they hold are described in
[Output](@ref measure-output).

For more information, see the reference or inline help for each phase, `SimData` and `runTMS`.

## Phases, algorithms and representations of one's own

Phases, time evolution algorithms and representations of a state of your own are described
in [Extending TMS](@ref).
