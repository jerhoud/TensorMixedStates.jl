# Extending the package from outside.
#
# Goes here: what a package built on TensorMixedStates relies on. A site type it declares at
# its top level is only checked by a package written for the occasion, precompiled and loaded
# in a process of its own: that code runs while the package is precompiled, and only what
# Julia keeps of that run is there once the package is loaded, while in the process of the
# tests the same code would run when included and show nothing. A representation or an
# algorithm of its own needs no such thing and is checked here directly.

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

@testset "data_to_frame without DataFrames" begin
    # a MethodError naming no remedy: a hint says what to load, and only while it is not loaded
    script = """
        using TensorMixedStates
        try
            data_to_frame(Dict())
        catch e
            print(occursin("needs the DataFrames package", sprint(showerror, e)))
        end
        """
    @test read(`$(Base.julia_cmd()) --project=$(Base.active_project()) -e $script`, String) == "true"
    message = try data_to_frame(1) catch e sprint(showerror, e) end
    @test !occursin("needs the DataFrames package", message)
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

TensorMixedStates.expect(st::WrappedState, op::IndexedOp) = expect(st.inner, op)
TensorMixedStates.expect1(st::WrappedState, ops) = expect1(st.inner, ops)
TensorMixedStates.expect2(st::WrappedState, pairs) = expect2(st.inner, pairs)

TensorMixedStates.write_state(g, st::WrappedState) = TensorMixedStates.write_state(g, st.inner)
TensorMixedStates.read_state(::Type{WrappedState}, g, sites, system) =
    WrappedState(TensorMixedStates.read_state(State{Pure}, g, sites, system))

# a state whose type has a parameter, which a state file does not record: write_state writes it
# in the group and a method of read_state for the type without parameters reads it back
struct Tagged{T} <: AbstractState
    system::System
    inner::State{Pure}
end

function TensorMixedStates.write_state(g, st::Tagged{T}) where T
    g["tag"] = T
    return TensorMixedStates.write_state(g, st.inner)
end

function TensorMixedStates.read_state(::Type{<:Tagged}, g, sites, system)
    inner = TensorMixedStates.read_state(State{Pure}, g, sites, system)
    return Tagged{read(g, "tag")}(inner.system, inner)
end

# and one read by a method for its types with parameters only, which a file cannot give
struct Untagged{T} <: AbstractState
    system::System
    inner::State{Pure}
end

TensorMixedStates.write_state(g, st::Untagged) = TensorMixedStates.write_state(g, st.inner)
function TensorMixedStates.read_state(::Type{Untagged{T}}, g, sites, system) where T
    inner = TensorMixedStates.read_state(State{Pure}, g, sites, system)
    return Untagged{T}(inner.system, inner)
end

# a state with no way of being written
struct Unwritable <: AbstractState
    system::System
end

@testset "What a state of one's own has no method of" begin
    # LinearAlgebra takes any argument for these, and they failed on iterate
    st = Unwritable(System(2, Qubit()))
    @test_throws "no method matching norm(::Unwritable)" norm(st)
    @test_throws "no method matching normalize(::Unwritable)" normalize(st)
    @test_throws "no method matching dot(::Unwritable, ::Unwritable)" dot(st, st)
    @test_throws "no method matching norm(::Unwritable)" measure(st, Norm)
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
        # a type with parameters is recorded without them
        save_state(file, "t", Tagged{3}(sys, State{Pure}(sys, "Up")))
        @test load_state(file, "t") isa Tagged{3}
        save_state(file, "u", Untagged{3}(sys, State{Pure}(sys, "Up")))
        @test_throws "is read by a method for the type without them" load_state(file, "u")
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
        # LoadState truncates only when it is given limits, and WrappedState has no truncate
        loaded = runTMS(SimData(phases = [LoadState(file = file, statename = "w")]);
                        output = devnull)
        @test loaded.state isa WrappedState
        @test_throws MethodError runTMS(SimData(phases = [
                LoadState(file = file, statename = "w", limits = Limits(maxdim = 2))]);
            output = devnull)
    end
end

# A representation on another system: a pure state whose physical sites are the odd sites of a
# doubled system, each followed by an ancilla left empty, as a purification interleaves them.
# It receives the operators as written and places them on its system with map_sites.
struct InterleavedState <: AbstractState
    system::System
    inner::State{Pure}
end

interleaved(i) = 2i - 1

TensorMixedStates.apply(op::IndexedOp, st::InterleavedState; kwargs...) =
    InterleavedState(st.system, apply(map_sites(interleaved, op), st.inner; kwargs...))

TensorMixedStates.expect(st::InterleavedState, op::IndexedOp) =
    expect(st.inner, map_sites(interleaved, op))

# the state on the doubled system goes in a subgroup, read back on the sites it lies on
TensorMixedStates.write_state(g, st::InterleavedState) =
    TensorMixedStates.write_state(TensorMixedStates.HDF5.create_group(g, "wide"), st.inner)

function TensorMixedStates.read_state(::Type{InterleavedState}, g, sites, system)
    wide = TensorMixedStates.read_state(State{Pure}, g["wide"],
                                        reduce(vcat, [ [s, s] for s in sites ]), nothing)
    return InterleavedState(isnothing(system) ? System(sites) : system, wide)
end

@testset "A representation on another system" begin
    sys = System(3, Fermion())
    hop(θ) = exp(-im * θ * (dag(C) ⊗ C + dag(dag(C) ⊗ C)))
    # the second hop goes past site 2, occupied or not, which the strings must account for
    gates = hop(0.4)(1, 3) * hop(0.7)(2, 3)
    st = apply(gates, State{Pure}(sys, ["Occ", "Occ", "Emp"]))
    wide = InterleavedState(sys, State{Pure}(System(6, Fermion()),
                                            ["Occ", "Emp", "Occ", "Emp", "Emp", "Emp"]))
    wide = apply(gates, wide)
    ops = [dag(C)(3) * C(1), C(1) * dag(C)(3), N(1) * N(3) + 0.5 * dag(C)(2) * C(3), 2.]
    # the values are real or complex as declared: in a vector of Number, ≈ compares them exactly
    @test [last.(measure(wide, ops))...] ≈ [last.(measure(st, ops))...]
    @test expect(wide, [N(2), dag(C)(1) * C(3)]) ≈ expect(st, [N(2), dag(C)(1) * C(3)])
    # a tuple of operators gives a tuple of values, on a State as on a state of one's own
    pair = (N(2), C(1) * dag(C)(3))
    @test expect(st, pair) isa Tuple
    @test collect(expect(wide, pair)) ≈ collect(expect(st, pair))
    mktempdir() do dir
        file = joinpath(dir, "wide.h5")
        save_state(file, "wide", wide)
        back = load_state(file, "wide"; system = sys)
        @test back isa InterleavedState
        @test [last.(measure(back, ops))...] ≈ [last.(measure(st, ops))...]
    end
end

# An algorithm of one's own: tdvp one step at a time, `run_steps` doing the bookkeeping. It
# makes the same steps as `Tdvp`, and so writes the same lines, interrupted or not.
struct Stepwise <: Algo end

function TensorMixedStates.evolve(::Stepwise, ::State, sim::Simulation, phase::Evolve;
                                  evolver, coefs, nsteps, kwargs...)
    dt = phase.duration / nsteps
    return run_steps(sim, nsteps) do sim, k
        sim = tdvp(evolver, dt, sim; phase.limits)
        if mod(k, phase.measures_period) == 0
            output(sim, phase.measures; sweep = k)
        end
        return sim
    end
end

# Tdvp on the state of the representation above, through the tdvp of the state it wraps
function TensorMixedStates.evolve(algo::Tdvp, st::WrappedState, sim::Simulation, phase::Evolve;
                                  kwargs...)
    s = TensorMixedStates.evolve(algo, st.inner, Simulation(sim, st.inner), phase; kwargs...)
    return Simulation(s, WrappedState(s.state))
end

@testset "An algorithm of one's own" begin
    mktempdir() do dir
        cd(dir) do
            # a measurement that creates the stop file when it has been taken `stop_in[]` times
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
            start = CreateState{Pure}(2, Qubit(), "Up")
            stepped(algo) = Evolve(duration = 0.3, time_step = 0.1, algo = algo,
                                   evolver = -im * (X(1) + Z(1) * Z(2)),
                                   limits = Limits(maxdim = 4, cutoff = 1e-15),
                                   measures = "data" => [X(1), Z(2), :sweep, stopper])
            runTMS(SimData(name = "tdvp", phases = [start, stepped(Tdvp())]))
            reference = read("tdvp/data", String)
            runTMS(SimData(name = "own", phases = [start, stepped(Stepwise())]))
            @test read("own/data", String) == reference
            # stopped after each of its steps and resumed, it goes on from where it stopped
            for k in 1:3
                stop_in[] = k
                data = SimData(name = "chk$k", phases = [start, stepped(Stepwise())])
                runTMS(data)
                @test stop_in[] == 0
                runTMS(data)
                @test read("chk$k/data", String) == reference
            end
        end
    end
    # Tdvp on the state of an extension, and an algorithm with no method for it
    sys = System(2, Qubit())
    brief(algo) = Evolve(duration = 0.2, time_step = 0.1, algo = algo, evolver = -im * X(1),
                         final_measures = Data("z") => Z(1))
    plain = runTMS(SimData(phases = [CreateState(type = Pure(), system = sys, state = "Up"),
                                     brief(Tdvp())]); output = devnull)
    wrapped = runTMS(SimData(phases = [CreateState(type = Wrapped(), system = sys, state = "Up"),
                                       brief(Tdvp())]); output = devnull)
    @test wrapped.state isa WrappedState
    @test last(wrapped.data["z"]["Z(1)"]["data"]) ≈ last(plain.data["z"]["Z(1)"]["data"])
    @test_throws "has no method for ApproxW on a WrappedState" runTMS(SimData(phases = [
            CreateState(type = Wrapped(), system = sys, state = "Up"),
            brief(ApproxW(order = 2))]); output = devnull)
end
