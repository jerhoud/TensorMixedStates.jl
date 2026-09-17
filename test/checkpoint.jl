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
    fire!(r) = r[] < 0 || (r[] > 0 && (r[] -= 1) == 0)
    stopper = StateFunc("Stopper", _ -> (fire!(stop_in) && touch("stop"); 0.))
    breaker = StateFunc("Breaker", _ -> begin
        fire!(fail_in) && throw(InterruptException())
        fire!(crash_in) && error("simulated kill")
        0.
    end)
    measures = ["data" => [X(1), Y(1), Z(2), stopper, breaker]]
    evolve(op) = Evolve(; duration = 0.3, time_step = 0.1, algo = Tdvp(), evolver = -im * op,
                        limits = Limits(maxdim = 10, cutoff = 1e-15), measures)
    return [CreateState{Pure}(3, Qubit(), "X+"), evolve(Z(1)), evolve(Z(2))]
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
            for _ in 1:10
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
    Base.exit_on_sigint(true)   # runTMS turned it off, leave the process as it was found
end

@testset "Phase fingerprint" begin
    id(phases) = TensorMixedStates.phases_id(TensorMixedStates.flatten_phases(phases))
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

    # phases built in pieces, as `create_graph_state` does, must run like a flat list. A
    # checkpoint numbers the phases, and numbering them without flattening first would
    # quietly skip the ones inside a sublist.
    mktempdir() do dir
        cd(dir) do
            phases = resume_phases(Ref(0), Ref(0), Ref(0))
            runTMS(SimData(; name = "flat", phases))
            runTMS(SimData(; name = "nested", phases = [[phases[1]], [phases[2:3]]],
                           checkpoint_interval = 1e-9))
            @test read("nested/data", String) == read("flat/data", String)
        end
    end
end

@testset "Checkpointer" begin
    C = TensorMixedStates.Checkpointer
    @test !TensorMixedStates.checkpoint_due(C())                     # 0 disables it
    @test TensorMixedStates.checkpoint_due(C(interval = 1e-9))
    @test !TensorMixedStates.stop_requested(C())
    @test TensorMixedStates.stop_requested(C(max_time = -1))         # deadline already past

    c = C()
    c.skip = 4
    @test TensorMixedStates.first_sweep!(c) == 5
    @test TensorMixedStates.first_sweep!(c) == 1                     # consumed by the first phase
end
