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
    measures = ["data" => [X(1), Y(1), Z(2), stopper, breaker]]
    evolve(op) = Evolve(; duration = 0.3, time_step = 0.1, algo = Tdvp(), evolver = -im * op,
                        limits = Limits(maxdim = 10, cutoff = 1e-15), measures)
    # a dmrg phase resumes on its sweep count rather than on a simulation time, which is a
    # path of its own through `DmrgObserver`
    ground = GroundState(; hamiltonian = sum(-Z(i) for i in 1:3), nsweeps = 3,
                         limits = Limits(maxdim = 10, cutoff = 1e-15), measures)
    return [CreateState{Pure}(3, Qubit(), "X+"), evolve(Z(1)), evolve(Z(2)), ground]
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
            # measurements were written past the last checkpoint, which is what the resume
            # has to cut back before running them again
            @test length(readlines("kill/data")) > checkpointed
            runTMS(sim_data)
            @test read("kill/data", String) == reference
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

@testset "Per sweep schedules" begin
    rs = TensorMixedStates.resume_schedule
    @test rs(1e-8, 3) == 1e-8                       # one value covers every sweep
    @test rs([1, 2, 3, 4], 2) == [3, 4]
    @test rs([1, 2, 3], 5) == [3]                   # a schedule that ran out keeps its last
    @test rs(Limits(cutoff = 1e-14, maxdim = [2, 4, 8]), 1).maxdim == [4, 8]

    # the evolution solvers pick their value out sweep by sweep instead, so that the
    # sweep numbers of the phase, which a resume keeps, land on the right one
    sv, sl = TensorMixedStates.sweep_value, TensorMixedStates.sweep_limits
    @test sv(1e-8, 3) == 1e-8                       # one value covers every sweep
    @test sv([2, 4, 8], 2) == 4
    @test sv([2, 4, 8], 7) == 8                     # a schedule that ran out keeps its last
    @test sl(Limits(cutoff = 1e-14, maxdim = [2, 4, 8]), 3) ==
          Limits(cutoff = 1e-14, maxdim = 8)

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
    @test id(base) ≠ id([first(base), Evolve(duration = 1., time_step = 0.1, algo = Tdvp(),
                                             evolver = -im * Z(1), measures = ["f" => X])])
    @test id(base) ≠ id(reverse(base))
    @test id(base) ≠ id(base[1:1])

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
    @test SimData(phases = [1, [2, [3, 4]], 5]).phases == [1, 2, 3, 4, 5]

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

    c = C()
    c.skip = 4
    @test TensorMixedStates.first_sweep!(c) == 5
    @test TensorMixedStates.first_sweep!(c) == 1                     # consumed by the first phase
end
