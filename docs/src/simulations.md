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
- `CopyState`: keep a copy of the state under a name, or put a named one back, for the
  measurements and the phases after it, a search of excited states for instance
- `Update`: change the simulation by a function of your own, in a few steps

A phase acts on the state of the simulation, `:state`, or on the state named by its field `on`,
`GroundState(…; on = :ground)`, which leaves the others as they are; `Evolve` and `Gates` take
several names, `on = [:state, :φ]`, evolving each state alike, the time advancing once, as the
correlations at two times need, see the example `examples/high_level/two_time_correlations.jl`.
The measurements read `:state`, and a named state by a list `:φ => [...]`, see
[Measurements](measurements.md).

With these phases we define a `SimData` object that describes the simulation and finally, we call

```julia
runTMS(simdata)
```

which executes the simulation.

Phases made of steps, the steps of an evolution or the sweeps of a search, take their
measurements after every step by default. The `measurements_period` field of `Evolve`,
`GroundState`, `SteadyState`, `Thermalize` and `Update` raises that interval:
`measurements_period = 10` measures one step out of ten, which is what long runs usually want.

The `phases` field is a list, but that list may contain lists, to any depth, and is
flattened before the simulation starts. This is meant for programs that build their phases
in pieces, a helper returning the several phases it needs rather than a single one.

## Examples

TMS is made for open systems, so here is one: six qubits evolving under a transverse
field Ising hamiltonian while each of them jumps from up to down. `CreateState{Mixed}` is what
makes the state a density matrix, and the `Dissipator` terms added to the hamiltonian are what
turn the evolution into a Lindblad equation.

```julia
using TensorMixedStates, .Qubits

runTMS(SimData(
    name = "dissipative_ising",
    description = "six qubits under a transverse field Ising hamiltonian, each jumping from up to down",
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

The file `data` is a table, a line naming its columns and then a line per time step, `Evolve`
measuring after each step: the magnetization on the six sites and the purity, here after the
first step and at the end of the run,

```
time  Z(1)               Z(2)               Z(3)               Z(4)               Z(5)               Z(6)               Purity
0.1   0.9409894108307    0.941176821604     0.9411768236482    0.9411768236508    0.9411768216059    0.9409894108332    0.7898070092902
...
2    -0.09518209547961  -0.1901119643657   -0.1449804250078   -0.1449804249806   -0.190111964564    -0.0951820957316    0.1086794010911
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

`runTMS` creates a directory named after the `name` field of the `SimData` object and writes
the files of the simulation there: a relative file name, of a destination or of `SaveState`
and `LoadState`, is taken there. The working directory of the process does not change, so a
file opened by a phase of your own with a relative name is taken from it. The directory holds
in particular:

- `log`: the progression of the computation;
- `prog.jl`: a copy of the script;
- `prog_args.json`: the command line arguments it was given, if any;
- `description`: the content of the `description` field of the `SimData` object;
- `stamp`: the versions, the date, the seed of the random draws, the BLAS library and the
  thread settings of the run;
- `running`: present during the computation, naming the machine and the process running it.
  `runTMS` refuses to start in a directory where another run may be going on. A file left by
  a process of the same machine that no longer exists, killed for instance, is replaced with
  a warning in the log; any other is removed by hand;
- `error`: the error that ended the run, if one did.

When `runTMS` returns, a line on the console says whether the simulation is done or stopped,
where its results are and how many warnings its log holds, which a program, or an AI running it,
reads without opening the log.

A file `restart` put in the directory makes the next run start over, erasing the directory
first, and `cleanTMS(sim_data)` erases it at once. `runTMS` takes four keyword arguments:

- `max_time` and `checkpoint_interval`, see [Long runs, checkpoints and stopping](@ref);
- `threading`, see [Threads and performance](@ref);
- `output` (default `nothing`): if set, create neither the directory nor the files `runTMS`
  writes in it, and redirect all output to the given stream, `stdout` or `devnull` for
  instance. A `SaveState` still writes its file, in the current directory.

The random draws of the phases use the generator of the run, seeded by `SimData(seed = …)`, or
by a seed drawn when there is none and written in `stamp`.

## Long runs, checkpoints and stopping

A simulation meant to run for hours or days can save its progress, so that a crash, a
batch system killing the job, or a deliberate stop does not throw the computation away.
Two keyword arguments of `runTMS` control it, or the environment variables
`TMS_CHECKPOINT_INTERVAL` and `TMS_MAX_TIME`, which a batch script sets without editing the
program.

```julia
runTMS(sim_data;
       checkpoint_interval = 600,    # seconds between two checkpoints
       max_time = 3.5 * 3600)        # stop cleanly after this long
```

`checkpoint_interval` is the time between two saves; `0`, the default, means no periodic save.
Even then, a stop, from `max_time`, the `stop` file or an interrupt, writes a checkpoint, so
that the simulation can always be resumed. `max_time` is a wall clock budget: once it is past,
the simulation writes a checkpoint and returns instead of carrying on. Set it comfortably below the limit
of your batch job, since a checkpoint is only taken between two sweeps, two steps or two
phases: a sweep that lasts ten minutes delays the stop by up to ten minutes.

A checkpoint is written in the simulation directory, the state and the rest of the
simulation, its time, its values and where its files were, to `checkpoint-1.h5` or
`checkpoint-2.h5`, and `checkpoint.json` names that file. It goes to the file the previous
checkpoint does not use, and `checkpoint.json` is moved into place last, so an interruption at
any point of the save leaves the previous checkpoint intact. A checkpoint is refused when it
was written in the layout of another version of TensorMixedStates, or by another minor version of
Julia.

Without a directory, when `runTMS` is given `output`, nothing is saved: an interrupt goes on
to the caller, and `max_time` stops the simulation with a message saying that it cannot be
resumed.

### Resuming

`runTMS` resumes on its own: run the same program again and it picks up where it left off,
skipping the phases that were finished and restarting the interrupted one at the step it
had reached. There is nothing to pass and nothing to change in the program. A checkpoint is
also written after the last phase, when periodic checkpoints are on or when one is already on
the disk, left by a stop, so that running the program once more after the simulation
completed does nothing. A simulation that never wrote one runs again from the start.

Output files are cut back to the length they had at the checkpoint before the simulation
continues, so the measurements written between the last checkpoint and the interruption
are not duplicated. The result is the file an uninterrupted run would have produced, with the
same values up to rounding: a search resumed between two sweeps rebuilds what dmrg kept from
one sweep to the next, which may change the last digits of a small value, as a residual close
to zero. A `Data` object is put back as it was at the checkpoint, the values measured after
it dropped. The log is the exception: it keeps the history of every run, so what an
interrupted run wrote after its last checkpoint stays, with the line saying why it stopped,
followed by the line marking the resume and by the steps done again from the checkpoint.

A file `restart` in the directory makes the run ignore an existing checkpoint and start over,
erasing the directory.

The generator of the random draws of the run is part of a checkpoint, so that a resumed run
draws the same numbers as the uninterrupted one, in a `CreateState` with `randomize`, or in a
measurement or a phase drawing with `c.rng` or `sim.rng`. A draw made with another generator,
the global one of Julia for instance, is not reproduced.

A checkpoint is resumed only by the program that wrote it: the program run must be the same,
byte for byte, as the copy `prog.jl` in the directory, and be given the same command line
arguments, which `prog_args.json` holds when there are some. Editing the program and running it again under the
same name thus reports an error instead of quietly continuing something else. To resume with
the edited program all the same, more time or other threads for instance, copy it onto
`prog.jl`, and its arguments into `prog_args.json`, or set the environment variable
`TMS_ACCEPT_PROGRAM=1`; to start over, put a file `restart` in the directory or give the
simulation another name. Only these two files are compared: a change in a file the
program includes, or in a parameter it reads elsewhere, from an environment variable or a data
file, goes unseen. A simulation run from the REPL, with no program file, is always resumed.

A program running several simulations copies itself into the directory of each, so that once
it is edited none of them resumes, those it did not change included, until it is copied onto
the `prog.jl` of each or they start over. Two simulations of one program must have different
names: under the same name, they share a directory, and the second resumes the checkpoint of
the first, of the same program, giving its results back as its own.

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

`Evolve`, `GroundState`, `SteadyState`, `Thermalize` and `Update` are resumed at the step or
the sweep they reached, `Gates` at the step of `apply` it reached, a batch of gates or an MPO
on some sites, and a phase of your own made of steps at the step it reached, see
[Phases of one's own](@ref own-phases).
The other phases are short enough to be replayed, and a checkpoint is taken between phases
whenever one is due.

## Measurements

Measurements are specified in the `measurements` field of `Evolve`, `GroundState`,
`SteadyState`, `Thermalize` and `Update`, taken after its steps, and in the
`final_measurements` field of every phase and of `SimData`, taken at its end. They take the form of a pair or list of pairs.

```julia
measurements = destination => measurements
measurements = [ dest1 => meas1, dest2 => meas2, ...]
```

The possible measurements are described on the [Measurements](measurements.md) page, the
values the phases provide, `:step` or `:energy`, included. A destination is a file, a JSON Lines
file or a `Data` object, which keeps the values in the program;
the destinations and the format of what they hold are described in
[Output](@ref measure-output).

For more information, see the reference or inline help for each phase, `SimData` and `runTMS`.

## Phases, methods, destinations and states of one's own

Phases, methods of evolution and of search, destinations and states of your own are described
in [Extending TMS](@ref).
