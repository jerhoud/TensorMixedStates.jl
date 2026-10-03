# Checkpointing: stopping a simulation and starting it again from where it stopped.
#
# Goes here: `Checkpointer`, the phase fingerprint `phases_id`, `flatten_phases` and the
# resume path of `runTMS`. The tests below create output directories, so they run inside a
# temporary directory and are the only ones of the suite to touch the disk. They also
# change the working directory of the process while they run, which is safe as long as the
# test groups stay sequential.

"""
phases exercising a resume, along with the three counters arming them. A counter fires
when it counts down to zero and then stays disarmed, `0` never fires and a negative value
fires every time.

`stop_in` asks the run to stop: the `stop` file is what a user creates to end a simulation
cleanly, but `runTMS` erases a leftover one when it starts, so it has to appear while the
run is going. A measurement is the one piece of user code called at every sweep, and it
runs with the simulation directory as working directory, which is exactly where the file is
looked for. This makes the stop fall on a known sweep instead of depending on a clock,
which a loaded machine would make unreliable.

`fail_in` raises an `InterruptException`, taking the same path as a real interrupt without
involving a signal, and `crash_in` raises an ordinary error, which stands for the simulation
being killed outright: no checkpoint is written, so the run resumes from an older one.

The measurements return the same constant whether armed or not, so that every run writes the
same columns and the outputs can be compared as they are.
"""
function resume_phases(stop_in::Ref{Int}, fail_in::Ref{Int}, crash_in::Ref{Int})
    function fire!(r)
        if r[] < 0
            return true
        elseif r[] == 0
            return false
        end
        r[] -= 1
        return r[] == 0
    end
    stopper = StateFunc("Stopper", _ -> begin
        if fire!(stop_in)
            touch("stop")
        end
        0.
    end)
    breaker = StateFunc("Breaker", _ -> begin
        if fire!(fail_in)
            throw(InterruptException())
        end
        if fire!(crash_in)
            error("simulated kill")
        end
        0.
    end)
    # the sweep of each measurement is written too: a resumed dmrg phase numbered its own from
    # 1 again, ITensorMPS counting afresh the sweeps it is asked for
    measures = ["data" => [X(1), Y(1), Z(2), :sweep, stopper, breaker]]
    evolve(op) = Evolve(; duration = 0.3, time_step = 0.1, algo = Tdvp(), evolver = -im * op,
                        limits = Limits(maxdim = 10, cutoff = 1e-15), measures)
    # a dmrg phase resumes on its sweep count rather than on a simulation time, which is a
    # path of its own through `DmrgObserver`
    ground = GroundState(; hamiltonian = sum(-Z(i) for i in 1:3), nsweeps = 3,
                         limits = Limits(maxdim = 10, cutoff = 1e-15), measures)
    return [CreateState{Pure}(3, Qubit(), "X+"), evolve(Z(1)), evolve(Z(2)), ground]
end

# phases of one's own driving tdvp with the observer of the package, as the docstring of
# `run_phase` describes: the first does not read its resume point, the second does
Base.@kwdef struct PlainEvolve
    name::String = "plain evolution"
    time_start = nothing
    final_measures = []
    measures = []
end

TensorMixedStates.run_phase(sim::Simulation, p::PlainEvolve) =
    tdvp(-im * X(1), 0.4, sim; nsweeps = 4, limits = Limits(maxdim = 4, cutoff = 1e-15),
         observer! = TdvpObserver(sim, p.measures, 1))

Base.@kwdef struct ResumingEvolve
    name::String = "resuming evolution"
    time_start = nothing
    final_measures = []
    measures = []
end

function TensorMixedStates.run_phase(sim::Simulation, p::ResumingEvolve)
    done, _ = resume_step(sim)
    return tdvp(-im * X(1), 0.4, sim; nsweeps = 4, first_sweep = done + 1,
                limits = Limits(maxdim = 4, cutoff = 1e-15),
                observer! = TdvpObserver(sim, p.measures, 1))
end

# a phase of one's own written as a loop of steps, each a kick on the first qubit and the time
# it takes, measured after it
Base.@kwdef struct Kicks
    name::String = "kicks"
    time_start = nothing
    final_measures = []
    nkicks::Int = 4
    measures = []
end

TensorMixedStates.run_phase(sim::Simulation, p::Kicks) =
    run_steps(sim, p.nkicks) do sim, k
        sim = apply(exp(-0.3im * X)(1), sim)
        sim = Simulation(sim, sim.state, sim.time + 0.1)
        output(sim, p.measures)
        return sim
    end

@testset "A SimData is not a phase" begin
    # A SimData inside `phases` used to be accepted and to silently skip phases: the loop it
    # opened shared the phase counter of the loop around it. It cost the first inner phase
    # on the very first run, with no checkpoint involved and no message. Grouping is done
    # with vectors.
    p = CreateState{Pure}(2, Qubit(), "Up")
    inner = SimData(name = "inner", phases = [p, p])
    @test_throws "cannot be a phase" SimData(name = "outer", phases = [p, inner])
    # what grouping is for, and it still flattens to any depth
    @test length(SimData(name = "flat", phases = [p, [p, [p, p]]]).phases) == 4
end

@testset "A simulation starts with its state" begin
    # it has none before its first phase, which used to fail on `nothing` inside its solver
    @test_throws "first phase must be" SimData(phases = [Gates(gates = X(1))])
    @test_throws "first phase must be" SimData(phases = [[], [Gates(gates = X(1))]])
    @test_throws "first phase must be" SimData(phases = [])
end

@testset "Resuming reproduces an uninterrupted run" begin
    mktempdir() do dir
        cd(dir) do
            stop_in, fail_in, crash_in = Ref(0), Ref(0), Ref(0)
            phases = resume_phases(stop_in, fail_in, crash_in)
            runTMS(SimData(; name = "ref", phases))
            reference = read("ref/data", String)

            # one sweep per run: the measurement of every sweep asks for a stop, so each
            # run does a single sweep and checkpoints. Runs made after the last phase is
            # over find a checkpoint pointing past the end and do nothing, so the count
            # only has to be large enough.
            stop_in[] = -1
            sim_data = SimData(; name = "chk", phases, checkpoint_interval = 1e-9)
            for _ in 1:14
                runTMS(sim_data)
            end
            @test read("chk/data", String) == reference

            # an interrupt in the middle of the last phase, then a plain resume
            stop_in[] = 0
            fail_in[] = 5
            sim_data = SimData(; name = "int", phases, checkpoint_interval = 1e-9)
            runTMS(sim_data)
            @test fail_in[] == 0                          # the interrupt did happen
            @test read("int/data", String) ≠ reference
            runTMS(sim_data)
            @test read("int/data", String) == reference

            # a simulation killed outright writes no checkpoint, so it resumes from an
            # older one and the measurements made in between are on the file twice unless
            # they are cut back. The interval here is long enough that the only checkpoint
            # is the one the stop writes.
            stop_in[] = 2
            sim_data = SimData(; name = "kill", phases, checkpoint_interval = 1e9)
            runTMS(sim_data)
            checkpointed = length(readlines("kill/data"))
            crash_in[] = 3
            @test_throws ErrorException runTMS(sim_data)
            @test crash_in[] == 0                         # the kill did happen
            @test isfile("kill/error")
            # measurements were written past the last checkpoint, which is what the resume
            # has to cut back before running them again
            @test length(readlines("kill/data")) > checkpointed
            runTMS(sim_data)
            @test read("kill/data", String) == reference
            # the marker describes the last run, which succeeded
            @test !isfile("kill/error")
        end
    end
end

@testset "Interrupting a phase without a solver" begin
    mktempdir() do dir
        cd(dir) do
            fail_in = Ref(0)
            breaker = StateFunc("Breaker", _ -> begin
                if fail_in[] > 0
                    fail_in[] -= 1
                    if fail_in[] == 0
                        throw(InterruptException())
                    end
                end
                0.
            end)
            measures = ["data" => [X(1), Y(1)]]
            evolve(op) = Evolve(; duration = 0.3, time_step = 0.1, algo = Tdvp(),
                                evolver = -im * op,
                                limits = Limits(maxdim = 10, cutoff = 1e-15), measures)
            # a `Gates` phase runs no sweep, so nothing records its progress while it runs.
            # An interrupt in the second one has to resume from the state the first one
            # produced and from sweep 0, not from the state and the sweep count the
            # evolution before them left behind. `measure` computes every measurement
            # before writing any of them, so the breaker throws without a partial line.
            phases = [CreateState{Pure}(3, Qubit(), "X+"), evolve(Z(1)),
                      Gates(; gates = Z(1)),
                      Gates(; gates = X(1), final_measures = ["data" => [breaker]]),
                      evolve(Z(2))]
            runTMS(SimData(; name = "ref", phases))
            reference = read("ref/data", String)

            # no checkpoint is due at the phase boundaries, so the only one written is the
            # one the interrupt asks for
            fail_in[] = 1
            sim_data = SimData(; name = "chk", phases, checkpoint_interval = 1e9)
            runTMS(sim_data)
            @test fail_in[] == 0                  # the interrupt did happen
            runTMS(sim_data)
            @test read("chk/data", String) == reference
        end
    end
end

@testset "Accumulating destinations survive a resume" begin
    mktempdir() do dir
        cd(dir) do
            stop_in = Ref(0)
            stopper = StateFunc("Stopper", _ -> begin
                if stop_in[] > 0
                    stop_in[] -= 1
                    if stop_in[] == 0
                        touch("stop")
                    end
                end
                0.
            end)
            # a text destination is continued from the position the checkpoint recorded,
            # but a json one and a `Data` one accumulate in memory and are only handed over
            # at the end, so the checkpoint has to carry what they hold
            measures = ["data" => [X(1), stopper], "out.json" => [X(1)], Data("d") => [X(1)]]
            phases = [CreateState{Pure}(3, Qubit(), "X+"),
                      Evolve(duration = 0.6, time_step = 0.1, algo = Tdvp(),
                             evolver = -im * Z(1),
                             limits = Limits(maxdim = 10, cutoff = 1e-15); measures)]
            ref = runTMS(SimData(; name = "ref", phases))
            reference = read("ref/data", String)
            ref_json = TensorMixedStates.JSON.parsefile("ref/out.json")

            stop_in[] = 3                                  # stop half way through the six
            sim_data = SimData(; name = "chk", phases, checkpoint_interval = 1e-9)
            runTMS(sim_data)
            @test stop_in[] == 0                           # the stop did happen
            sim = runTMS(sim_data)
            @test read("chk/data", String) == reference
            @test TensorMixedStates.JSON.parsefile("chk/out.json") == ref_json
            @test sim.data["d"]["X(1)"]["times"] == ref.data["d"]["X(1)"]["times"]
        end
    end
end

@testset "Complex values and matrices survive a resume" begin
    # json holds neither complex numbers nor matrices, so the checkpoint marks them in a
    # `Data` destination and rebuilds them: the resumed run hands back the values, element
    # types included, that an uninterrupted one gives
    mktempdir() do dir
        cd(dir) do
            stop_in = Ref(0)
            stopper = StateFunc("Stopper", _ -> begin
                if stop_in[] > 0
                    stop_in[] -= 1
                    if stop_in[] == 0
                        touch("stop")
                    end
                end
                0.
            end)
            ms = [Sp(1), (X, Z), (Sp, Sm)]
            measures = ["data" => [ms; stopper], "out.json" => ms, Data("d") => ms]
            phases = [CreateState{Pure}(2, Qubit(), "X+"),
                      Evolve(duration = 0.6, time_step = 0.1, algo = Tdvp(),
                             evolver = -im * (Z(1) + X(1)X(2)),
                             limits = Limits(maxdim = 10, cutoff = 1e-15); measures)]
            ref = runTMS(SimData(; name = "ref", phases))
            stop_in[] = 3
            sim_data = SimData(; name = "chk", phases, checkpoint_interval = 1e-9)
            runTMS(sim_data)
            @test stop_in[] == 0
            sim = runTMS(sim_data)
            @test read("chk/data", String) == read("ref/data", String)
            @test TensorMixedStates.JSON.parsefile("chk/out.json") ==
                  TensorMixedStates.JSON.parsefile("ref/out.json")
            for name in ["Sp(1)", "XZ", "SpSm"]
                resumed, whole = sim.data["d"][name]["data"], ref.data["d"][name]["data"]
                @test map(typeof, resumed) == map(typeof, whole)
                @test all(resumed .≈ whole)
            end
        end
    end
end

@testset "A complex time and an older checkpoint" begin
    mktempdir() do dir
        cd(dir) do
            # a complex simulation time is marked in the checkpoint as well
            phases = [CreateState{Pure}(2, Qubit(), "Up"),
                      Gates(gates = X(1), time_start = 0.5im, final_measures = Data("d") => Z(1)),
                      Gates(gates = X(1))]
            sim_data = SimData(; name = "ctime", phases, checkpoint_interval = 1e-9)
            runTMS(sim_data)
            k = TensorMixedStates.load_checkpoint("ctime")
            data = TensorMixedStates.restore_series(k.outputs["data"]["d"])
            @test only(data["Z(1)"]["times"]) == 0.5im
            @test k.phase_time == 0.5im
            # and a checkpoint of an earlier version is refused: its values do not say which
            # call of output they came from
            meta = "ctime/checkpoint.json"
            write(meta, replace(read(meta, String), "\"version\":3" => "\"version\":2"))
            @test_throws "has version 2" runTMS(sim_data)
        end
    end
end

# a measurement that creates the stop file when it has been taken `stop_in[]` times
function stopper_at(stop_in::Ref{Int})
    return StateFunc("Stopper", _ -> begin
        if stop_in[] > 0
            stop_in[] -= 1
            if stop_in[] == 0
                touch("stop")
            end
        end
        0.
    end)
end

@testset "A resumed ground state search stops where it would have" begin
    # the tolerance is checked against the energy of the sweep before, which the checkpoint
    # carries, a stop on it records the phase as done, and the state checkpointed is the one
    # of the sweep it is labelled with: a resume, wherever it falls, runs the same sweeps
    # and writes the same lines as the uninterrupted search
    mktempdir() do dir
        cd(dir) do
            stop_in = Ref(0)
            h = -sum(Z(i) * Z(i + 1) for i in 1:3) - 0.8 * sum(X(i) for i in 1:4)
            phases = [CreateState{Pure}(4, Qubit(), "Up"),
                      GroundState(hamiltonian = h, nsweeps = 30, tolerance = 1e-10,
                                  limits = Limits(maxdim = 8, cutoff = 1e-14),
                                  measures = "data" => [stopper_at(stop_in), Z(1)],
                                  final_measures = "final" => [Z(1)])]
            runTMS(SimData(; name = "ref", phases))
            reference = read("ref/data", String)
            sweeps = count(l -> startswith(l, "Z(1)"), split(reference, '\n'))
            @test 1 < sweeps < 30                  # the tolerance stops the search
            for k in 1:sweeps
                stop_in[] = k
                sim_data = SimData(; name = "chk$k", phases, checkpoint_interval = 1e-9)
                runTMS(sim_data)
                @test stop_in[] == 0
                runTMS(sim_data)
                @test read("chk$k/data", String) == reference
                @test read("chk$k/final", String) == read("ref/final", String)
            end
        end
    end
end

@testset "A resume after the last phase" begin
    # a stop asked for in the final measurements of the last phase checkpoints past it: the
    # resume runs no phase, and the final measurements of the simulation, which the stopped
    # run leaves to it, are taken at the time the simulation reached
    mktempdir() do dir
        cd(dir) do
            stop_in = Ref(0)
            phases = [CreateState{Pure}(2, Qubit(), "Up"),
                      Evolve(duration = 0.3, time_step = 0.1, algo = Tdvp(), evolver = -im * X(1),
                             limits = Limits(maxdim = 4, cutoff = 1e-15),
                             final_measures = "data" => [stopper_at(stop_in)])]
            ref = runTMS(SimData(; name = "ref", phases, final_measures = "fin" => Z(1)))
            stop_in[] = 1
            sim_data = SimData(; name = "chk", phases, final_measures = "fin" => Z(1),
                               checkpoint_interval = 1e9)
            runTMS(sim_data)
            @test stop_in[] == 0
            sim = runTMS(sim_data)
            @test sim.time ≈ ref.time
            @test read("chk/fin", String) == read("ref/fin", String)
            @test read("chk/data", String) == read("ref/data", String)
        end
    end
end

@testset "An output the checkpoint does not know is created anew" begin
    # a file opened after the checkpoint a run resumes from holds what a killed attempt left:
    # the uninterrupted run creates it, so the resumed one does too, where a file the
    # checkpoint knows is continued from where it was cut back
    mktempdir() do dir
        cd(dir) do
            write("late", "left by a killed attempt\n")
            write("known", "kept\ncut back\n")
            sim = Simulation(nothing)
            TensorMixedStates.restore_outputs!(sim.outputs,
                Dict("files" => Dict("known" => Dict("text" => 5)), "data" => Dict()))
            println(get_sim_file(sim, "late"), "new")
            println(get_sim_file(sim, "known"), "more")
            TensorMixedStates.close_sim_files(sim)
            @test read("late", String) == "new\n"
            @test read("known", String) == "kept\nmore\n"
        end
    end
end

@testset "Values that are not finite" begin
    # json has no number for them: the checkpoint marks them and gives them back, and a json
    # file writes them as Julia prints them
    mktempdir() do dir
        cd(dir) do
            phases = [CreateState{Pure}(2, Qubit(), "Up"),
                      Gates(gates = X(1), final_measures = [Data("d") => [NaN, -Inf], "out.json" => [Inf]]),
                      Gates(gates = X(1))]
            runTMS(SimData(; name = "nonfinite", phases, checkpoint_interval = 1e-9))
            @test only(TensorMixedStates.JSON.parsefile("nonfinite/out.json")["Inf"]["data"]) == "Inf"
            k = TensorMixedStates.load_checkpoint("nonfinite")
            data = TensorMixedStates.restore_series(k.outputs["data"]["d"])
            @test isnan(only(data["NaN"]["data"]))
            @test only(data["-Inf"]["data"]) == -Inf
        end
    end
end

@testset "Running a completed simulation again" begin
    # with periodic checkpoints on, one is written after the last phase, so that the run
    # resumes past every phase: nothing is computed again and the files stay as they were
    mktempdir() do dir
        cd(dir) do
            phases = [CreateState{Pure}(2, Qubit(), "Up"),
                      Evolve(duration = 0.3, time_step = 0.1, algo = Tdvp(), evolver = -im * X(1),
                             limits = Limits(maxdim = 4, cutoff = 1e-15), measures = "data" => [Z(1)])]
            sim_data = SimData(; name = "done", phases, checkpoint_interval = 1e9,
                               final_measures = "fin" => Z(1))
            first = runTMS(sim_data)
            data, fin = read("done/data", String), read("done/fin", String)
            again = runTMS(sim_data)
            @test again.time ≈ first.time
            @test read("done/data", String) == data
            @test read("done/fin", String) == fin
            @test occursin("Resuming from checkpoint: phase 3", read("done/log", String))
        end
    end
end

@testset "Running again a simulation that was stopped and completed" begin
    # without periodic checkpoints, the checkpoint a stop wrote stayed on the disk once the
    # resumed run had completed, and the next run resumed from it: the files were cut back to
    # the stop and their end computed again
    mktempdir() do dir
        cd(dir) do
            stop_in = Ref(2)
            phases = [CreateState{Pure}(2, Qubit(), "Up"),
                      Evolve(duration = 0.3, time_step = 0.1, algo = Tdvp(), evolver = -im * X(1),
                             limits = Limits(maxdim = 4, cutoff = 1e-15),
                             measures = "data" => [Z(1), stopper_at(stop_in)])]
            sim_data = SimData(; name = "done", phases)
            runTMS(sim_data)
            @test stop_in[] == 0                          # the stop did happen
            first = runTMS(sim_data)                      # resumed and completed
            data = read("done/data", String)
            again = runTMS(sim_data)
            @test again.time ≈ first.time
            @test read("done/data", String) == data
            @test occursin("Resuming from checkpoint: phase 3", read("done/log", String))
        end
    end
end

@testset "A resume logs the time it resumes from" begin
    # it logged the time the interrupted phase had started from, 0 here, rather than the time
    # its last committed sweep had reached
    mktempdir() do dir
        cd(dir) do
            stop_in = Ref(2)
            phases = [CreateState{Pure}(2, Qubit(), "Up"),
                      Evolve(duration = 0.3, time_step = 0.1, algo = Tdvp(), evolver = -im * X(1),
                             limits = Limits(maxdim = 4, cutoff = 1e-15),
                             measures = "data" => [Z(1), stopper_at(stop_in)])]
            sim_data = SimData(; name = "chk", phases)
            runTMS(sim_data)
            @test stop_in[] == 0
            runTMS(sim_data)
            m = match(r"Resuming from checkpoint: phase 2, sweep 2, simulation time (\S+)",
                      read("chk/log", String))
            @test !isnothing(m)
            @test parse(Float64, m[1]) ≈ 0.2
        end
    end
end

@testset "Per sweep schedules" begin
    rs = TensorMixedStates.resume_schedule
    @test rs(1e-8, 3) == 1e-8                       # one value covers every sweep
    @test rs([1, 2, 3, 4], 2) == [3, 4]
    @test rs([1, 2, 3], 5) == [3]                   # a schedule that ran out keeps its last
    resumed = rs(Limits(cutoff = 1e-14, maxdim = [2, 4, 8], mindim = [1, 2, 3]), 1)
    @test resumed.maxdim == [4, 8]
    @test resumed.mindim == [2, 3]

    # the evolution solvers pick their value out sweep by sweep instead, so that the
    # sweep numbers of the phase, which a resume keeps, land on the right one
    sv, sl = TensorMixedStates.sweep_value, TensorMixedStates.sweep_limits
    @test sv(1e-8, 3) == 1e-8                       # one value covers every sweep
    @test sv([2, 4, 8], 2) == 4
    @test sv([2, 4, 8], 7) == 8                     # a schedule that ran out keeps its last
    @test sl(Limits(cutoff = 1e-14, maxdim = [2, 4, 8], mindim = [1, 2, 3]), 3) ==
          Limits(cutoff = 1e-14, maxdim = 8, mindim = 3)

    # `first_sweep` is what hides that difference from the phases: dmrg is asked for the
    # sweeps that are left and handed the tail of its schedules, so resuming at sweep 3 of
    # 4 has to be the same run as asking for 2 sweeps with that tail written out by hand
    sys = System(6, Qubit())
    ham = sum(-Z(i) * Z(i + 1) for i in 1:5) - sum(1. * X(i) for i in 1:6)
    st = State{Pure}(sys, "Z+")
    e1, _ = dmrg(ham, st; nsweeps = 4, first_sweep = 3,
                 limits = Limits(cutoff = 1e-14, maxdim = [2, 2, 8, 8]),
                 noise = [1e-2, 1e-3, 0., 0.])
    e2, _ = dmrg(ham, st; nsweeps = 2,
                 limits = Limits(cutoff = 1e-14, maxdim = [8, 8]), noise = [0., 0.])
    @test e1 == e2

    # a ground state resuming half way has to go on with the maxdim its sweep was due, not
    # start the schedule over, which would truncate a state the uninterrupted run kept
    mktempdir() do dir
        cd(dir) do
            stop_in = Ref(0)
            stopper = StateFunc("Stopper", _ -> begin
                if stop_in[] > 0
                    stop_in[] -= 1
                    if stop_in[] == 0
                        touch("stop")
                    end
                end
                0.
            end)
            phases = [
                CreateState{Pure}(6, Qubit(), "Z+"),
                GroundState(; hamiltonian = sum(-Z(i) * Z(i + 1) for i in 1:5) -
                                            sum(1. * X(i) for i in 1:6),
                            nsweeps = 4, limits = Limits(maxdim = [2, 2, 8, 8], cutoff = 1e-14),
                            noise = [1e-2, 1e-3, 0., 0.],
                            measures = ["data" => [MaxLinkdim, X(1), Z(1)Z(2), stopper]]),
            ]
            runTMS(SimData(; name = "ref", phases))
            stop_in[] = 1
            sim_data = SimData(; name = "chk", phases, checkpoint_interval = 1e-9)
            runTMS(sim_data)
            runTMS(sim_data)
            @test read("chk/data", String) == read("ref/data", String)
        end
    end

    # an evolution counts its sweeps from the start of the phase and resumes on that
    # count, so its schedule has to be read at the sweep being run, not restarted
    mktempdir() do dir
        cd(dir) do
            stop_in = Ref(0)
            stopper = StateFunc("Stopper", _ -> begin
                if stop_in[] > 0
                    stop_in[] -= 1
                    if stop_in[] == 0
                        touch("stop")
                    end
                end
                0.
            end)
            phases = [
                CreateState{Pure}(6, Qubit(), "X+"),
                Evolve(duration = 0.4, time_step = 0.1, algo = Tdvp(),
                       evolver = -im * (sum(-Z(i) * Z(i + 1) for i in 1:5) -
                                        sum(1. * X(i) for i in 1:6)),
                       limits = Limits(maxdim = [2, 2, 8, 8], cutoff = 1e-14),
                       measures = ["data" => [MaxLinkdim, X(1), Z(1)Z(2), stopper]]),
            ]
            runTMS(SimData(; name = "ref", phases))
            # the link dimension really does follow the schedule, otherwise the run below
            # would agree with the reference for want of anything to disagree about
            @test length(unique(l -> split(l)[1] == "MaxLinkdim" ? split(l)[3] : "",
                                readlines("ref/data"))) > 2
            stop_in[] = 2
            sim_data = SimData(; name = "chk", phases, checkpoint_interval = 1e-9)
            runTMS(sim_data)
            runTMS(sim_data)
            @test read("chk/data", String) == read("ref/data", String)
        end
    end
end

@testset "A resumed steady state has trace one" begin
    # the eigenvector dmrg gives has norm one and a sign of its own, and the solver normalizes
    # it only when it returns. A steady state checkpointed on its last sweep, whose resume
    # does not run the solver, was handed to the next phases and saved with a trace of -1.33
    mktempdir() do dir
        cd(dir) do
            stop_in = Ref(0)
            lind = -im * (X(1) + 0.7 * X(2) + Z(1) * Z(2)) + Dissipator(Sm)(1) + Dissipator(Sp)(2)
            phases = [CreateState{Mixed}(3, Qubit(), "Up"),
                      SteadyState(lindbladian = lind, nsweeps = 4,
                                  limits = Limits(cutoff = 1e-12, maxdim = 32),
                                  measures = "data" => [stopper_at(stop_in)]),
                      Evolve(duration = 0.2, time_step = 0.1, algo = Tdvp(), evolver = lind,
                             limits = Limits(cutoff = 1e-12, maxdim = 32)),
                      SaveState(file = "saved.h5")]
            ref = runTMS(SimData(; name = "ref", phases))
            for k in (4, 2)
                stop_in[] = k
                sim_data = SimData(; name = "chk$k", phases, checkpoint_interval = 1e-9)
                stopped = runTMS(sim_data)
                @test stop_in[] == 0
                @test trace(stopped.state) ≈ 1
                sim = runTMS(sim_data)
                @test trace(sim.state) ≈ 1
                @test trace(load_state("chk$k/saved.h5", "state")) ≈ 1
                @test expect1(sim.state, Z) ≈ expect1(ref.state, Z)
            end
        end
    end
end

@testset "A resume keeps the system of the phases" begin
    # the state of a checkpoint came back on a system of its own, so that a measurement
    # comparing with a state built on the system of the phases failed on every resume
    mktempdir() do dir
        cd(dir) do
            sys = System(2, Qubit())
            ref = State{Pure}(sys, "Up")
            fid = StateFunc("Fid", st -> fidelity(ref, st))
            for mixed in (false, true)
                stop_in = Ref(0)
                evolver = mixed ? -im * X(1) + Dissipator(Sm)(2) : -im * X(1)
                phases = [CreateState(type = Pure(), system = sys, state = "Up");
                          mixed ? [ToMixed()] : [];
                          Evolve(; duration = 0.3, time_step = 0.1, algo = Tdvp(), evolver,
                                 measures = "data" => [fid, stopper_at(stop_in)])]
                runTMS(SimData(; name = "ref$mixed", phases))
                stop_in[] = 2
                sim_data = SimData(; name = "chk$mixed", phases, checkpoint_interval = 1e-9)
                runTMS(sim_data)
                @test stop_in[] == 0
                sim = runTMS(sim_data)
                @test sim.state.system === sys
                @test read("chk$mixed/data", String) == read("ref$mixed/data", String)
            end
        end
    end
end

@testset "stopped tells a stopped run from a completed one" begin
    # the number of lines written, or the time reached, was the only way to tell
    mktempdir() do dir
        cd(dir) do
            stop_in = Ref(2)
            phases = [CreateState{Pure}(2, Qubit(), "Up"),
                      Evolve(duration = 0.3, time_step = 0.1, algo = Tdvp(), evolver = -im * X(1),
                             measures = "data" => [stopper_at(stop_in)])]
            sim_data = SimData(; name = "sim", phases, checkpoint_interval = 1e-9)
            @test stopped(runTMS(sim_data))
            @test !stopped(runTMS(sim_data))
            @test !stopped(Simulation(State{Pure}(System(2, Qubit()), "Up")))
        end
    end
end

@testset "The log keeps the history of every run" begin
    # cut back to the checkpoint on a resume, as the measurements are, the log lost the line
    # saying why a run had stopped
    mktempdir() do dir
        cd(dir) do
            stop_in = Ref(2)
            phases = [CreateState{Pure}(2, Qubit(), "Up"),
                      Evolve(duration = 0.3, time_step = 0.1, algo = Tdvp(), evolver = -im * X(1),
                             measures = "data" => [stopper_at(stop_in)])]
            sim_data = SimData(; name = "sim", phases, checkpoint_interval = 1e-9)
            runTMS(sim_data)
            runTMS(sim_data)
            log = read("sim/log", String)
            stop = findfirst("the simulation can be resumed", log)
            @test !isnothing(stop)
            @test first(stop) < first(findfirst("Resuming from checkpoint", log))
            # stopped in the course of its evolution, which it resumes from where it was
            @test occursin("Stopping in phase 2", log)
            @test occursin("***** Stopping phase \"Time evolution\"", log)
            @test count("Evolving state from simulation time 0.0 to", log) == 1
        end
    end
end

@testset "SaveState refuses a file of the simulation" begin
    # a state saved as a checkpoint file was destroyed by the next checkpoint
    mktempdir() do dir
        cd(dir) do
            phases(file) = [CreateState{Pure}(2, Qubit(), "Up"), SaveState(; file)]
            @test_throws "a file of the simulation directory" runTMS(SimData(name = "a",
                phases = phases("checkpoint-1.h5")))
            runTMS(SimData(name = "b", phases = phases("mine.h5")))
            @test isfile("b/mine.h5")
        end
    end
end

@testset "Two names of one file" begin
    # each name opened the file, emptying what the other had written
    mktempdir() do dir
        cd(dir) do
            runTMS(SimData(name = "sim", phases = [
                CreateState{Pure}(2, Qubit(), "Up"),
                Evolve(duration = 0.3, time_step = 0.1, algo = Tdvp(), evolver = -im * X(1),
                       measures = ["data" => X(1), "./data" => Z(1), "sub/../data" => Y(1)])]))
            @test length(readlines("sim/data")) == 9
        end
    end
end

@testset "The copy of the program runs again" begin
    # run again from the directory of the simulation, prog.jl was copied onto itself, which cp
    # refused, and the run failed
    mktempdir() do dir
        cd(dir) do
            write("prog.jl", """
                using TensorMixedStates, .Qubits
                runTMS(SimData(name = "sim", phases = [CreateState{Pure}(2, Qubit(), "Up")]))
                """)
            include(joinpath(dir, "prog.jl"))
            @test isfile("sim/prog.jl")
            include(joinpath(dir, "sim", "prog.jl"))
            @test !isfile("sim/error")
            @test read("sim/prog.jl", String) == read("prog.jl", String)
        end
    end
end

@testset "A restart does not remove the current directory" begin
    # with ".", an ancestor or the absolute path of the current directory, rm emptied it, the
    # program included, before failing on the directory itself
    mktempdir() do dir
        cd(dir) do
            phases = [CreateState{Pure}(2, Qubit(), "Up")]
            mkpath("sub")
            touch("sub/precious")
            cd("sub") do
                refusal = "which contains the current directory"
                for name in (".", "..", pwd())
                    @test_throws refusal runTMS(SimData(; name, phases); restart = true)
                    @test_throws refusal runTMS(SimData(; name, phases); clean = true)
                end
                @test isfile("precious")
                # an ordinary name is removed and run again
                runTMS(SimData(; name = "sim", phases))
                touch("sim/old")
                runTMS(SimData(; name = "sim", phases); restart = true)
                @test !isfile("sim/old")
                @test isfile("sim/log")
            end
        end
    end
end

@testset "A resumed dmrg with a measurement period" begin
    # a deadline already past stops every run after one sweep or one phase. A stop measured a
    # sweep its period skips, the line of a checkpointed sweep was cut from the log, and a
    # search checkpointed on its last sweep lost the line it ends with
    mktempdir() do dir
        cd(dir) do
            h = -sum(Z(i) * Z(i + 1) for i in 1:3) - 0.8 * sum(X(i) for i in 1:4)
            phases = [CreateState{Pure}(4, Qubit(), "Up"),
                      GroundState(hamiltonian = h, nsweeps = 5,
                                  limits = Limits(maxdim = 8, cutoff = 1e-14),
                                  measures = "data" => [Z(1), :sweep], measures_period = 2)]
            runTMS(SimData(; name = "ref", phases))
            sim_data = SimData(; name = "chk", phases, max_time = -1)
            runTMS(sim_data)
            runTMS(sim_data)
            # stopped after its first sweep, the search is not done
            @test !occursin("Done, dmrg", read("chk/log", String))
            for _ in 1:8
                runTMS(sim_data)
            end
            @test read("chk/data", String) == read("ref/data", String)
            lines(f) = filter(l -> startswith(l, "sweep") || startswith(l, "Done"), readlines(f))
            @test lines("chk/log") == lines("ref/log")
        end
    end
end

@testset "A phase of one's own resumes correctly" begin
    # its sweeps were checkpointed, and the resume handed it the state it had reached, on
    # which it ran all its sweeps again: it ended at <Z(1)> = 0.362 instead of 0.697
    mktempdir() do dir
        cd(dir) do
            for (name, P) in (("plain", PlainEvolve), ("resuming", ResumingEvolve))
                stop_in = Ref(0)
                phases = [CreateState{Pure}(2, Qubit(), "Up"),
                          P(measures = "data" => [Z(1), stopper_at(stop_in)])]
                ref = runTMS(SimData(; name = "ref$name", phases))
                stop_in[] = 2
                sim_data = SimData(; name = "chk$name", phases, checkpoint_interval = 1e-9)
                stopped = runTMS(sim_data)
                @test stop_in[] == 0
                # a stopped run hands back what it resumes from
                @test stopped.time ≈ (name == "plain" ? 0 : 0.2)
                sim = runTMS(sim_data)
                @test read("chk$name/data", String) == read("ref$name/data", String)
                @test real(expect(sim.state, Z(1))) ≈ real(expect(ref.state, Z(1)))
                @test sim.time ≈ ref.time
            end
        end
    end
end

# a loop of one's own whose steps each run tdvp with the observer of the package
Base.@kwdef struct Evolutions
    name::String = "evolutions"
    time_start = nothing
    final_measures = []
    measures = []
end

TensorMixedStates.run_phase(sim::Simulation, p::Evolutions) =
    run_steps(sim, 2) do sim, k
        tdvp(-im * X(1), 0.2, sim; nsweeps = 2, limits = Limits(maxdim = 4, cutoff = 1e-15),
             observer! = TdvpObserver(sim, p.measures, 1))
    end

@testset "A loop of one's own resumes after its last step" begin
    # run_steps commits each step once it is measured: a stop between two steps resumes after
    # the last one, from the state and the simulation time it had reached
    mktempdir() do dir
        cd(dir) do
            stop_in = Ref(0)
            phases = [CreateState{Pure}(2, Qubit(), "Up"),
                      Kicks(measures = "data" => [Z(1), stopper_at(stop_in)],
                            final_measures = "final" => Z(1))]
            ref = runTMS(SimData(; name = "ref", phases))
            stop_in[] = 2
            sim_data = SimData(; name = "chk", phases)
            stopped = runTMS(sim_data)
            @test stop_in[] == 0
            # a stopped run hands back what it resumes from, and takes no final measurement
            @test stopped.time ≈ 0.2
            @test !isfile("chk/final")
            sim = runTMS(sim_data)
            @test read("chk/data", String) == read("ref/data", String)
            @test read("chk/final", String) == read("ref/final", String)
            @test sim.time ≈ ref.time ≈ 0.4
            @test real(expect(sim.state, Z(1))) ≈ real(expect(ref.state, Z(1)))
            # outside runTMS the steps simply run
            @test run_steps((s, k) -> Simulation(s, s.state, s.time + 1), Simulation(ref.state), 3).time ≈ 3
            # a step has to hand back the simulation
            @test_throws "has to return the simulation" run_steps((s, k) -> s.state, Simulation(ref.state), 1)
        end
    end
end

@testset "A solver within a step of one's own" begin
    # the sweeps of a solver run within a step were committed as steps of the phase, and a stop
    # falling in the middle of the solver committed the unfinished step: the resume went on
    # after it, from a state half evolved
    mktempdir() do dir
        cd(dir) do
            stop_in = Ref(0)
            phases = [CreateState{Pure}(2, Qubit(), "Up"),
                      Evolutions(measures = "data" => [Z(1), stopper_at(stop_in)])]
            ref = runTMS(SimData(; name = "ref", phases))
            # the third measurement is the first sweep of the second step
            stop_in[] = 3
            sim_data = SimData(; name = "chk", phases)
            stopped = runTMS(sim_data)
            @test stop_in[] == 0
            @test stopped.time ≈ 0.2                      # the last step done, the first
            sim = runTMS(sim_data)
            @test read("chk/data", String) == read("ref/data", String)
            @test sim.time ≈ ref.time ≈ 0.4
            @test real(expect(sim.state, Z(1))) ≈ real(expect(ref.state, Z(1)))
        end
    end
end

@testset "A simulation built by hand writes its json files once closed" begin
    mktempdir() do dir
        cd(dir) do
            sim = Simulation(State{Pure}(System(2, Qubit()), "Up"))
            output(sim, "d.json" => Z(1))
            @test !isfile("d.json")
            close_sim_files(sim)
            @test only(TensorMixedStates.JSON.parsefile("d.json")["Z(1)"]["data"]) ≈ 1
        end
    end
end

@testset "A destination is not a file of the simulation" begin
    # a destination called stop stopped the simulation at its first sweep and was erased by
    # the next run, one called checkpoint.json overwrote the checkpoint, and so on
    mktempdir() do dir
        cd(dir) do
            p = CreateState{Pure}(2, Qubit(), "Up")
            for name in ("stop", "log", "checkpoint.json", "./running")
                @test_throws "a file of the simulation directory" runTMS(SimData(name = "s",
                    phases = [p, Gates(gates = X(1), final_measures = name => [Z(1)])]))
            end
            # every file a run and its checkpoints leave in the directory is one of those
            runTMS(SimData(name = "all", description = "d", checkpoint_interval = 1e-9,
                           phases = [p, Evolve(duration = 0.2, time_step = 0.1, algo = Tdvp(),
                                               evolver = -im * X(1), measures = "data" => Z(1))]))
            @test issubset(setdiff(readdir("all"), ["data"]), TensorMixedStates.simulation_files)
            # without a directory nothing is written there, and any name goes
            @test_ok runTMS(SimData(phases = [p, Gates(gates = X(1), final_measures = "stop" => Z(1))]);
                            output = devnull)
        end
    end
end

@testset "An interrupt with no directory reaches the caller" begin
    # nothing can be saved, so nothing can be resumed: the interrupt was taken for a
    # checkpoint that was never written, and the run returned as if it had completed
    breaker = StateFunc("Breaker", _ -> throw(InterruptException()))
    phases = [CreateState{Pure}(2, Qubit(), "Up"),
              Gates(gates = X(1), final_measures = "data" => [breaker])]
    @test_throws InterruptException runTMS(SimData(; phases); output = devnull)
end

@testset "A checkpoint is replaced whole" begin
    # the state file the metadata names is the one read, whatever a crash left beside it: the
    # state and the metadata were renamed one after the other, and a kill between the two paired
    # the new state with the previous counts
    mktempdir() do dir
        cd(dir) do
            stop_in = Ref(0)
            phases = [CreateState{Pure}(3, Qubit(), "X+"),
                      Evolve(duration = 0.4, time_step = 0.1, algo = Tdvp(), evolver = -im * Z(1),
                             limits = Limits(maxdim = 10, cutoff = 1e-15),
                             measures = "data" => [X(1), stopper_at(stop_in)])]
            runTMS(SimData(; name = "ref", phases))
            stop_in[] = 2
            sim_data = SimData(; name = "chk", phases, checkpoint_interval = 1e-9)
            runTMS(sim_data)
            h5 = filter(endswith(".h5"), readdir("chk"))
            @test length(h5) == 1
            other = only(h5) == "checkpoint-1.h5" ? "checkpoint-2.h5" : "checkpoint-1.h5"
            write(joinpath("chk", other), "left half written")
            runTMS(sim_data)
            @test read("chk/data", String) == read("ref/data", String)
        end
    end
end

@testset "A complex time with no imaginary part" begin
    # it came back real, and the time took one column instead of two after the resume
    mktempdir() do dir
        cd(dir) do
            stop_in = Ref(0)
            phases = [CreateState{Pure}(2, Qubit(), "Up"),
                      Gates(gates = X(1), time_start = complex(0.5),
                            final_measures = "data" => [Z(1), stopper_at(stop_in)]),
                      Gates(gates = X(1), final_measures = "data" => [Z(1)])]
            runTMS(SimData(; name = "ref", phases))
            stop_in[] = 1
            sim_data = SimData(; name = "chk", phases, checkpoint_interval = 1e9)
            runTMS(sim_data)
            @test stop_in[] == 0
            runTMS(sim_data)
            @test read("chk/data", String) == read("ref/data", String)
        end
    end
end

@testset "Phase fingerprint" begin
    # the phases of a simulation are what `SimData` made of them, flattened
    id(phases) = TensorMixedStates.phases_id(SimData(; phases).phases)
    base = [CreateState{Pure}(3, Qubit(), "X+"),
            Evolve(duration = 1., time_step = 0.1, algo = Tdvp(), evolver = -im * Z(1))]

    # how the phases are grouped says nothing about what is computed
    @test id(base) == id([first(base), [last(base)]])
    @test id(base) == id([[[first(base)]], last(base)])
    # and it does not depend on the indices a session happens to draw
    @test id(base) == id(deepcopy(base))
    # a dictionary or a set is what it holds, in no order, and no longer raises
    ph = TensorMixedStates.phase_hash
    @test ph(UInt(0), Dict("a" => 1, "b" => 2)) == ph(UInt(0), Dict("b" => 2, "a" => 1))
    @test ph(UInt(0), Dict("a" => 1)) ≠ ph(UInt(0), Dict("a" => 2))
    @test ph(UInt(0), Set([1, 2])) == ph(UInt(0), Set([2, 1]))
    # FNV-1a, whose definition is fixed, where `Base.hash` changes between versions of Julia
    # and refused a checkpoint after an upgrade as belonging to another simulation
    o = TensorMixedStates.fnv_offset
    @test TensorMixedStates.fnv(o, codeunits("a")) == 0xaf63dc4c8601ec8c   # its published value
    @test ph(o, 1.5) == 0xaa95e93229a27c80
    @test ph(o, "abc") == 0xc11ab6d2519bc2b2
    @test ph(o, 3) == 0xc7c2bf3b330983e6
    @test ph(o, 1 + 2im) == 0x7717980363c8e066
    # a type is written by its full path, not as the program happens to see it: `Qubit` or
    # `TensorMixedStates.Qubit`, depending on the `using`, gave the same simulation two
    # fingerprints
    tk = TensorMixedStates.type_key
    @test tk(Qubit) == "TensorMixedStates.Qubit"
    @test tk(State{Pure}) == "TensorMixedStates.State{TensorMixedStates.Pure}"
    @test tk(Vector{Float64}) == "Core.Array{Core.Float64,1}"
    @test tk(Union{Int, Nothing}) == "Union{Core.Int64,Core.Nothing}"

    # a checkpoint of another simulation must be refused, so anything a phase says has to
    # count. Resuming into the wrong simulation is silent, which is what makes it serious.
    @test id(base) ≠ id([CreateState{Pure}(3, Qubit(), "X-"), last(base)])
    @test id(base) ≠ id([CreateState{Mixed}(3, Qubit(), "X+"), last(base)])
    @test id(base) ≠ id([CreateState{Pure}(4, Qubit(), "X+"), last(base)])
    @test id(base) ≠ id([CreateState{Pure}(3, Boson(2), "X+"), last(base)])
    @test id(base) ≠ id([first(base), Evolve(duration = 2., time_step = 0.1,
                                             algo = Tdvp(), evolver = -im * Z(1))])
    @test id(base) ≠ id([first(base), Evolve(duration = 1., time_step = 0.1,
                                             algo = Tdvp(), evolver = -im * Z(2))])
    @test id(base) ≠ id([first(base), Evolve(duration = 1., time_step = 0.1,
                                             algo = ApproxW(order = 2), evolver = -im * Z(1))])
    # the parameters of the solvers count too
    @test id(base) ≠ id([first(base), Evolve(duration = 1., time_step = 0.1,
                                             algo = Tdvp(krylov = Krylov(tol = 1e-10)),
                                             evolver = -im * Z(1))])
    @test id([first(base), Evolve(duration = 1., time_step = 0.1, algo = ApproxW(order = 2),
                                  evolver = -im * Z(1))]) ≠
          id([first(base), Evolve(duration = 1., time_step = 0.1,
                                  algo = ApproxW(order = 2, apply_algo = "naive"),
                                  evolver = -im * Z(1))])
    @test id(base) ≠ id([first(base), Evolve(duration = 1., time_step = 0.1, algo = Tdvp(),
                                             evolver = -im * Z(1), measures = ["f" => X])])
    @test id([base; Gates(gates = X(1)); Gates(gates = Z(1))]) ≠
          id([base; Gates(gates = Z(1)); Gates(gates = X(1))])
    @test id(base) ≠ id(base[1:1])

    # an anonymous function counts by the variables it captures, not by the name of its type,
    # which a counter gives: the same program included again in a session named its
    # functions anew, and its checkpoint was refused as another simulation's. Two copies of
    # one function are two such names
    timed(c) = [first(base), Evolve(duration = 1., time_step = 0.1, algo = Tdvp(),
                                    evolver = [-im * Z(1)] => [c])]
    @test id(timed(t -> exp(-t))) == id(timed(t -> exp(-t)))
    decay1(g) = t -> exp(-g * t)
    decay2(g) = t -> exp(-g * t)
    @test id(timed(decay1(1.))) == id(timed(decay2(1.)))
    # a function with a name counts by it
    @test id(timed(sin)) ≠ id(timed(cos))

    # a State given as is, rather than described, is part of what the simulation computes
    sys = System(2, Qubit())
    with(st) = [CreateState(type = Pure(), system = sys, state = st)]
    @test id(with(State{Pure}(sys, "Z+"))) == id(with(State{Pure}(sys, "Z+")))
    @test id(with(State{Pure}(sys, "Z+"))) ≠ id(with(State{Pure}(sys, "Z-")))
end

@testset "Nested phases" begin
    flatten = TensorMixedStates.flatten_phases
    @test flatten([1, [2, [3, 4]], 5]) == [1, 2, 3, 4, 5]
    @test flatten([[], [1], []]) == [1]
    @test flatten(1) == [1]

    # phases built in pieces, as `create_graph_state` does, must run like a flat list. A
    # checkpoint numbers the phases, and numbering them without flattening first would
    # quietly skip the ones inside a sublist.
    mktempdir() do dir
        cd(dir) do
            phases = resume_phases(Ref(0), Ref(0), Ref(0))
            runTMS(SimData(; name = "flat", phases))
            runTMS(SimData(; name = "nested", phases = [[phases[1]], [phases[2:end]]],
                           checkpoint_interval = 1e-9))
            @test read("nested/data", String) == read("flat/data", String)
        end
    end
end

@testset "Checkpointer" begin
    C = TensorMixedStates.Checkpointer
    @test !TensorMixedStates.checkpoint_due(C())                     # 0 disables it
    @test !TensorMixedStates.checkpoint_due(C(interval = -1))        # and so does any value below
    @test TensorMixedStates.checkpoint_due(C(interval = 1e-9))
    @test !TensorMixedStates.stop_requested(C())
    @test TensorMixedStates.stop_requested(C(max_time = -1))         # deadline already past

    # the resume point belongs to the phase it was written in, and is read once
    c = C()
    sim = Simulation(nothing; checkpoint = c)
    TensorMixedStates.commit!(c, sim.outputs, 2, 0, 0., 0., nothing)
    c.resume = TensorMixedStates.Commit(2, 4, 0., 0., nothing, -1.5, (files = Dict(), data = Dict()))
    @test resume_step(sim) == (4, -1.5)
    @test resume_step(sim) == (0, nothing)
    c.resume = TensorMixedStates.Commit(3, 4, 0., 0., nothing, nothing, (files = Dict(), data = Dict()))
    @test resume_step(sim) == (0, nothing)                         # another phase
end
