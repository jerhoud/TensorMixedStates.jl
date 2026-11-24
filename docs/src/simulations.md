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
- `ConvertState{R}`: go to representation `R`, of which `ToMixed` is `ConvertState{Mixed}`
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

- `log`: the progression of the computation: the machine and the versions of each launch, the
  bond dimension of the MPO of each solver, and for each step its values, the largest bond
  dimension of the state and its duration;
- `prog.jl`: a copy of the script;
- `prog`: a copy of the files the script includes, if any, at their paths relative to the
  script, a `..` becoming `_parent`;
- `prog_args.json`: the command line arguments it was given, if any;
- `description`: the content of the `description` field of the `SimData` object;
- `stamp`: the versions of Julia, of TensorMixedStates, of the packages doing the computations
  and of those the script loads, the date, the machine, the seed of the random draws, the BLAS
  library and the thread settings of the run;
- `running`: present during the computation, naming the machine and the process running it.
  `runTMS` refuses to start in a directory where another run may be going on. A file left by
  a process of the same machine that no longer exists, killed for instance, is replaced with
  a warning in the log; any other is removed by hand;
- `error`: the error that ended the run, if one did, that of a sample too, whose part holds the
  whole error;
- `lock`: guards the files several processes share;
- with samples, `branch.h5`, `queue.json`, `merged.jls`, `parts` and `pending`: the branch point,
  the queue of the samples, their merge, the directories of the parts and of the rows not merged
  yet. `merged.jls` is needed to prolong a run or to read its `Data` again.

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
are not duplicated. The result is the same file as an uninterrupted run would have
produced. A `Data` object is put back as it was at the checkpoint, the values measured after
it dropped. The log is the exception: it keeps the history of every run, so what an
interrupted run wrote after its last checkpoint stays, with the line saying why it stopped,
followed by the line marking the resume and by the steps done again from the checkpoint.

A file `restart` in the directory makes the run ignore an existing checkpoint and start over,
erasing the directory.

The generator of the random draws of the run is part of a checkpoint, so that a resumed run
draws the same numbers as the uninterrupted one, in a `CreateState` with `randomize`, or in a
measurement or a phase drawing with `c.rng` or `sim.rng`. A draw made with another generator,
the global one of Julia for instance, is not reproduced. A seed gives the same draws for one minor
version of Julia, which `stamp` records, not necessarily for another.

A checkpoint is resumed only by the program that wrote it: the program run must be the same,
byte for byte, as the copy `prog.jl` in the directory, include the same files as the copies in
`prog`, and be given the same command line arguments, which `prog_args.json` holds when there
are some. Editing the program or a file it includes and running it again under the same name
thus reports an error instead of quietly continuing something else. To resume with the edited
program all the same, more time or other threads for instance, copy it onto `prog.jl`, the
files it includes into `prog`, and its arguments into `prog_args.json`, or set the environment
variable `TMS_ACCEPT_PROGRAM=1`, which the log notes; to start over, put a file `restart` in
the directory or give the simulation another name. The files included before
`using TensorMixedStates` are not seen, nor the packages the program loads, whose versions
`stamp` gives, nor a parameter it reads elsewhere, from an environment variable or a data file.
A simulation run from the REPL, with no program file, is always resumed.

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

## Samples and disorder

A simulation with `samples` runs its last phases many times, once for each sample, each sample
with its own draw of the disorder. The phases before the last phase, `Realizations`, run once
and lead to the branch point; each sample starts again from the state they reached and runs the
phases the function of `Realizations` returns. A simulation of quantum trajectories branches in
the same way at its first `CreateState{Trajectory}` or `ConvertState{Trajectory}`, without
`Realizations`, see [Quantum trajectories](@ref).

```julia
runTMS(SimData(
    name = "disordered_ground",
    seed = 1,
    samples = 1000,
    phases = [
        CreateState{Pure}(10, Qubit(), "Up"),
        Realizations() do c
            J = [1 + 0.1 * randn(c.rng) for i in 1:9]
            H = -sum(J[i] * Z(i) * Z(i + 1) for i in 1:9) - sum(X(i) for i in 1:10)
            GroundState(hamiltonian = H, algo = Dmrg(nsweeps = 20, tol = 1e-10),
                        limits = Limits(maxdim = 64))
        end,
    ],
    final_measurements = ["energy" => [:energy, Z], EachSample("energies") => :energy],
))
```

The function is given `c.sample`, the index of the sample, and `c.rng`, the generator it draws
the disorder with, and returns the phases of that sample. It must draw with `c.rng` alone, so
that a sample can be made again, and return the same phases, by type and name, for every
sample: only their operators and values change. It is called once when the `SimData` is
written, so that its phases are checked before the run, and the phases of each sample are
checked again before it runs.

Each sample draws with generators of its own, made from the seed of the run and the index of the
sample: `c.rng` for its disorder, and `sim.rng` for its phases and measurements. A sample is
thus the same whichever run computes it, and two simulations of the same `seed` whose functions
draw alike get the same disorder in each sample, whatever their phases: put the draw in a
function both programs include to compare them on the same realizations.

The final measurements of the `SimData` are made at the end of each sample. A file or a `Data`
then holds averages over the samples: a row for each row of a sample, matched by its phase and
its step, the final measurements of a phase whatever its number of steps, with the columns `N`,
the samples in it, `N_eff`, their effective number, equal to `N` here, and for each column `c`
of the values, `c`, its mean, `err(c)`, the error on the mean, and `std(c)`, the standard
deviation of the samples, `NaN` for one sample. The file `energy` above holds one row, the
average ground state energy and magnetization. A row of a `Data` holds these cells, `row["Z(1)"]`
and `row["err(Z(1))"]`, and gives the vector of the cells of a value of every site or of a
vector, `row["Z"]` and `row["err(Z)"]`, as a row without samples gives the value; `row["N"]` is
the number of samples, and a value named `N`, as the operator `N` of every site, is read by its
cells, `row["N(1)"]`, or named otherwise, `"n" => N`. A row that some samples lack, a step that
only the slow searches reach for instance, averages those that have it. The samples of a row
must have the same time, but for rounding, and the same columns, and a value that is not a
number cannot be averaged. A row written before the samples is shared by all, and written with
no error and no deviation. A measurement written averaged or for each sample cannot be named
`sample`, `weight` or `trajectory`, which the samples write or will. The means and deviations
are updated sample after sample by the algorithm of West (1979), each value shifted by its value
in the first sample (Chan, Golub and LeVeque 1983), accurate even when the deviation is small
against the mean.

`EachSample(d)` writes to `d` a row for each sample instead, after the time the column `sample`
giving its index, 0 before the samples: the file `energies` above holds the energy of each
sample, for distributions or analyses of one's own after the run. The value `:sample` gives the
index of the sample to any measurement. A stream, `"stdout"` for instance, and a destination of
one's own are written for each sample, without that column.

The phases before the branch point run once, in the directory of the simulation, whose log
tells them, and end with the file `branch.h5`. The samples are then run by parts, each in a
directory of its own, `parts/1`, `parts/2`, …, with its log, which tells the phases and steps of
its first sample only and a line for each sample: a part takes the next sample no part has
taken, runs it, and leaves its rows in the directory `pending`. The samples ended are merged into
the files and the `Data` of the simulation in their order, up to the first one not ended, so
that they hold samples 1 to m, the same to the last bit whichever process ran each, unless a
sample reads `sim.data`, which holds the samples merged when it runs, and an average never
favours the samples that end first. An averaged file is written whole at each
merge.

Running the program several times at once, `julia prog.jl & julia prog.jl & …`, runs the samples
on as many cores, each process with its own memory: a process waits while another runs the
phases before the branch point, takes a part, and ends when no sample is left, its `Simulation`
holding the state of the last sample it ran and the `Data` merged so far, which `sim.data` also
gives in a sample; the console says how many samples are merged and how many processes still
run. A sample runs on a single thread, BLAS and Strided included, so that it is the same
whatever the threads of the process, which change the rounding of the products: launch a
process for each core instead. The `threading` of `runTMS` applies to the phases before the
branch point; `threading = nothing` leaves the threading of the samples as it is, and their
results may then differ by rounding from one process to another. A file a sample writes, by
`SaveState` for instance, needs a name of its own for each sample, `"state$(c.sample).h5"`: the
processes write at once.

A part stopped, or left by a process that died, is taken over by the next process, which resumes
it where its checkpoint is. The file `stop` stops all the parts in the directory of the
simulation, one part in its own; a run finding no other process running removes it. A sample
that fails, or that cannot be merged, makes the other processes fail with its error once they
end their sample, and the next run tries it again. The environment variable `TMS_SAMPLES`
overrides `samples`: running the program again with `TMS_SAMPLES=2000` after 1000 samples goes
on with sample 1001 from the branch point, without running the phases before it again, and
writes the same files as 2000 samples in one run. With fewer samples than those merged, a run
changes nothing. A run is prolonged with the minor version of Julia that started it; another
version of TensorMixedStates goes on when it reads its files, which the log notes, the samples
then coming from two versions.

A tolerance stops the samples once their averages are precise enough: with
`tol = "energy" => 1e-3`, or `"energy" => (atol = 1e-3, rtol = 1e-2)`, or a list of them,
the samples stop at the first whose errors in `energy`, every `err(c)` at every row, are within
`atol + rtol * |c|`, every row holding at least `min_samples` samples, 10 by default, so that an
error underestimated by chance does not stop them; `samples` is then the most that runs. A
quantity constant over the first samples, as a rare event, has an error of zero there: raise
`min_samples` above the samples it takes to show. The tolerance names an averaged file or
`Data`, a `Data` alone in a run without a directory. The samples merged are the same whichever
process runs them, the samples running beyond the first that reaches the tolerance ending
without being merged. A program run again with a tighter tolerance goes on from the samples
merged; with a looser one, it keeps them.

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
