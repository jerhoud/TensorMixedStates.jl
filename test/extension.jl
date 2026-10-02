# Extending the package from outside.
#
# Goes here: what a package built on TensorMixedStates relies on. A site type it declares at
# its top level is only checked by a package written for the occasion, precompiled and loaded
# in a process of its own: that code runs while the package is precompiled, and only what
# Julia keeps of that run is there once the package is loaded, while in the process of the
# tests the same code would run when included and show nothing. A representation of its own,
# with its state type, needs no such thing and is checked here directly.

"""
    run_with_packages(script, dir)

the output of the Julia code `script`, run in a process of its own which loads the packages of
the directory `dir` besides those of the environment of the tests, and precompiles them in a
depot of its own, which leaves nothing in the depot of the user.

A package of `dir` is written without a project file, which lets it load any package of that
environment: with one, Julia would look for its dependencies in `dir` alone.
"""
function run_with_packages(script::String, dir::String)
    sep = Sys.iswindows() ? ";" : ":"
    mktempdir() do depot
        cmd = `$(Base.julia_cmd()) --project=$(Base.active_project()) -e $script`
        return read(addenv(cmd, "JULIA_LOAD_PATH" => join([dir; LOAD_PATH], sep),
                                "JULIA_DEPOT_PATH" => join([depot; DEPOT_PATH], sep)), String)
    end
end

@testset "A site type declared in a package" begin
    # @def_states and @def_operators used to fill dictionaries of TensorMixedStates, which a
    # package writes into while it is precompiled and which Julia does not keep: once loaded,
    # the package had its site type and its operator names, but neither states nor matrices
    mktempdir() do dir
        src = joinpath(dir, "SiteProbe", "src")
        mkpath(src)
        write(joinpath(src, "SiteProbe.jl"), """
            module SiteProbe

            using TensorMixedStates, TensorMixedStates.Qubits
            import TensorMixedStates: dim

            export Probe, Pz

            struct Probe <: AbstractSite end

            dim(::Probe) = 2

            @def_states(Probe(), [ "a" => [1., 0.], ["b", "β"] => [0., 1.] ])

            @def_operators(Probe(), [ involution_op => [ Pz = [1. 0. ; 0. -1.] ] ])

            # an operator declared from outside for a site type of TensorMixedStates
            @def_operators(Qubit(), [ plain_op => [ Probed = [0. 0. ; 1. 0.] ] ])

            end
            """)
        # each check prints true, or the error it raised, which the comparison below then shows
        out = run_with_packages("""
            using TensorMixedStates, TensorMixedStates.Qubits, SiteProbe
            check(f) = try f() catch e sprint(showerror, e) end
            println(check(() -> state(Probe(), "β") == [0., 1.]))
            println(check(() -> matrix(Pz, Probe()) == [1. 0. ; 0. -1.]))
            println(check(() -> real(expect(State{Pure}(System(2, Probe()), ["a", "b"]),
                                            Pz(2))) ≈ -1))
            println(check(() -> matrix(SiteProbe.Probed, Qubit()) == [0. 0. ; 1. 0.]))
            """, dir)
        @test out == "true\ntrue\ntrue\ntrue\n"
    end
end

# A representation of one's own: a state wrapping a pure `State`, which goes through the
# phases of the package, its measurements, its state files and its checkpoints like a `State`.
struct Wrapped <: Representation end

struct WrappedState <: AbstractState
    system::System
    inner::State{Pure}
end

WrappedState(inner::State{Pure}) = WrappedState(inner.system, inner)

TensorMixedStates.run_phase(sim::Simulation, phase::CreateState{Wrapped}) =
    Simulation(sim, WrappedState(State{Pure}(phase.system, phase.state)))

TensorMixedStates.apply(op::IndexedOp, st::WrappedState; kwargs...) =
    WrappedState(apply(op, st.inner; kwargs...))

TensorMixedStates.expect_norm(st::WrappedState, terms::Vector) =
    TensorMixedStates.expect_norm(st.inner, terms)
TensorMixedStates.expect1(st::WrappedState, ops) = expect1(st.inner, ops)
TensorMixedStates.expect2(st::WrappedState, pairs) = expect2(st.inner, pairs)

TensorMixedStates.write_state(g, st::WrappedState) = TensorMixedStates.write_state(g, st.inner)
TensorMixedStates.read_state(::Type{WrappedState}, g, sites, system) =
    WrappedState(TensorMixedStates.read_state(State{Pure}, g, sites, system))

# a state with no way of being written
struct Unwritable <: AbstractState
    system::System
end

@testset "A representation of one's own" begin
    sys = System(3, Qubit())
    phases = [CreateState(type = Wrapped(), system = sys, state = "Up"),
              Gates(gates = X(1) * H(2))]
    measures = Data("m") => [Z, Z(1) * Z(3), (Z, Z)]
    # the threading :auto reads the system of the state before each phase
    sim = runTMS(SimData(phases = phases, final_measures = measures, threading = :auto);
                 output = devnull)
    @test sim.state isa WrappedState
    @test last(sim.data["m"]["Z"]["data"]) ≈ [-1, 0, 1]
    vals = measure(sim.state, [Z(1) * Z(3), (Z, Z)])
    @test last(vals[1]) ≈ -1
    @test last(vals[2]) ≈ [1 0 -1 ; 0 1 0 ; -1 0 1]
    mktempdir() do dir
        file = joinpath(dir, "wrapped.h5")
        save_state(file, "w", sim.state)
        back = load_state(file, "w"; system = sim.state.system)
        @test back isa WrappedState
        @test abs(inner(back.inner, sim.state.inner)) ≈ 1
        # refused before the file is touched, which leaves the state saved under that name
        @test_throws "has no method of TensorMixedStates.write_state" save_state(file, "w",
                                                                               Unwritable(sys))
        @test load_state(file, "w") isa WrappedState
        # a run stopped after its first phase writes the state in its checkpoint, and the run
        # resumed from it reads it back and finishes as the uninterrupted one
        cd(dir) do
            stopped = runTMS(SimData(name = "wrapped", phases = phases,
                                     final_measures = measures, max_time = 0))
            @test stopped.state isa WrappedState
            @test !haskey(stopped.data, "m")
            resumed = runTMS(SimData(name = "wrapped", phases = phases, final_measures = measures))
            @test resumed.state isa WrappedState
            @test last(resumed.data["m"]["Z"]["data"]) ≈ [-1, 0, 1]
        end
    end
end
