# Extending the package from outside.
#
# Goes here: what a package built on TensorMixedStates relies on. A site type it declares at
# its top level is only checked by a package written for the occasion, precompiled and loaded
# in a process of its own: that code runs while the package is precompiled, and only what
# Julia keeps of that run is there once the package is loaded, while in the process of the
# tests the same code would run when included and show nothing. A representation of its own
# needs no such thing and is checked here directly.
#
# The algorithms and the phases of one's own run by runTMS are here too, their resume with the
# others in runner.jl.

using TensorMixedStates: AbstractState

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
    # the package had its site type and its operator names, but neither states nor matrices.
    # They now declare on the site type (decision 19)
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

            @def_states Probe [ "a" => [1., 0.], ["b", "β"] => [0., 1.] ]

            @def_operators Probe [ Pz = [1. 0. ; 0. -1.] ]

            # an operator declared from outside for a site type of TensorMixedStates
            @def_operators Qubit [ Probed = [0. 0. ; 1. 0.] ]

            end
            """)
        # each check prints true, or the error it raised, which the comparison below then shows
        out = run_with_packages("""
            using TensorMixedStates, TensorMixedStates.Qubits, SiteProbe
            check(f) = try f() catch e sprint(showerror, e) end
            println(check(() -> TensorMixedStates.local_state(Probe(), "β") == [0., 1.]))
            println(check(() -> matrix(Pz, Probe()) == [1. 0. ; 0. -1.]))
            println(check(() -> real(expect(State{Pure}(System(2, Probe()), ["a", "b"]),
                                            Pz(2))) ≈ -1))
            println(check(() -> matrix(SiteProbe.Probed, Qubit()) == [0. 0. ; 1. 0.]))
            """, dir)
        @test out == "true\ntrue\ntrue\ntrue\n"
    end
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

# A representation of one's own: a state wrapping a pure `State`, which goes through the state
# files like a `State`
struct WrappedState <: AbstractState
    system::System
    inner::State{Pure}
end

TensorMixedStates.write_state(g, st::WrappedState) = TensorMixedStates.write_state(g, st.inner)

function TensorMixedStates.read_state(::Type{WrappedState}, g, sites)
    inner = TensorMixedStates.read_state(State{Pure}, g, sites)
    return WrappedState(inner.system, inner)
end

# phases of one's own creating a WrappedState and turning its first qubit, step by step
Base.@kwdef struct CreateWrapped <: TensorMixedStates.AbstractPhase
    name::String = "wrapped"
    time_start = nothing
    final_measurements = []
end

TensorMixedStates.creates_state(::CreateWrapped) = true

function TensorMixedStates.run_phase(::CreateWrapped, sim::Simulation, ::Int)
    st = State{Pure}(System(3, Qubit()), "Up")
    return Simulation(sim, WrappedState(st.system, st))
end

Base.@kwdef struct TurnWrapped <: TensorMixedStates.AbstractPhase
    name::String = "turn"
    time_start = nothing
    final_measurements = []
    measurements = []
end

TensorMixedStates.prepare_phase(p::TurnWrapped, _) = p
TensorMixedStates.step_count(::TurnWrapped) = 3
TensorMixedStates.run_step(::TurnWrapped, sim, _) =
    Simulation(sim, WrappedState(sim.state.system, apply(exp(-0.4im * X)(1), sim.state.inner)),
               sim.time + 1)

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

function TensorMixedStates.read_state(::Type{<:Tagged}, g, sites)
    inner = TensorMixedStates.read_state(State{Pure}, g, sites)
    return Tagged{read(g, "tag")}(inner.system, inner)
end

# and one read by a method for its types with parameters only, which a file cannot give
struct Untagged{T} <: AbstractState
    system::System
    inner::State{Pure}
end

TensorMixedStates.write_state(g, st::Untagged) = TensorMixedStates.write_state(g, st.inner)

function TensorMixedStates.read_state(::Type{Untagged{T}}, g, sites) where T
    inner = TensorMixedStates.read_state(State{Pure}, g, sites)
    return Untagged{T}(inner.system, inner)
end

@testset "A representation of one's own" begin
    sys = System(3, Qubit())
    st = WrappedState(sys, apply(X(1) * H(2), State{Pure}(sys, "Up")))
    # run by phases of one's own, measured, and resumed through its state files
    mktempdir() do dir
        cd(dir) do
            z = "z" => c -> real(expect(c.state.inner, Z(1)))
            phases = [CreateWrapped(), TurnWrapped(measurements = "z" => z)]
            ref = runTMS(SimData(; name = "ref", phases))
            sim = runTMS(SimData(; name = "chk", phases); max_time = -1)
            while stopped(sim)
                sim = runTMS(SimData(; name = "chk", phases); max_time = -1)
            end
            @test sim.state isa WrappedState
            @test read("chk/z", String) == read("ref/z", String)
            @test sim.time == 3
            # LoadState reads it back, not truncated at the default limits, which it has no
            # method for, and :auto chooses the threading from its system
            save_state("w.h5", "w", st)
            sim = runTMS(SimData(name = "load", phases = [LoadState(file = "../w.h5", statename = "w")]);
                         threading = :auto)
            @test sim.state isa WrappedState
            @test_throws MethodError runTMS(SimData(name = "cut", phases = [
                LoadState(file = "../w.h5", statename = "w", limits = Limits(maxdim = 2))]))
        end
    end
    mktempdir() do dir
        file = joinpath(dir, "wrapped.h5")
        save_state(file, "w", st)
        back = load_state(file, "w")
        @test back isa WrappedState
        @test abs(inner(back.inner, st.inner)) ≈ 1
        # refused before the file is touched, which leaves the state saved under that name
        @test_throws "has no method of TensorMixedStates.write_state" save_state(file, "w",
                                                                               Unwritable(sys))
        @test load_state(file, "w") isa WrappedState
        # a type with parameters is recorded without them, and a type read only with them
        # has, for the file, no method of read_state
        save_state(file, "t", Tagged{3}(sys, State{Pure}(sys, "Up")))
        @test load_state(file, "t") isa Tagged{3}
        save_state(file, "u", Untagged{3}(sys, State{Pure}(sys, "Up")))
        @test_throws "Untagged, its type has no method of TensorMixedStates.read_state" load_state(
            file, "u")
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

function TensorMixedStates.read_state(::Type{InterleavedState}, g, sites)
    wide = TensorMixedStates.read_state(State{Pure}, g["wide"],
                                        reduce(vcat, [ [s, s] for s in sites ]))
    return InterleavedState(System(sites), wide)
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
    ops = [dag(C)(3) * C(1), C(1) * dag(C)(3), N(1) * N(3) + 0.5 * dag(C)(2) * C(3)]
    @test [ expect(wide, o) for o in ops ] ≈ [ expect(st, o) for o in ops ]
    # measure goes through its expect, the values of a Hermitian operator real there too
    @test [last.(measure(wide, ops))...] ≈ [last.(measure(st, ops))...]
    @test last(only(measure(wide, dag(C)(1) * C(3) + dag(C)(3) * C(1)))) isa Float64
    @test expect(wide, [N(2), dag(C)(1) * C(3)]) ≈ expect(st, [N(2), dag(C)(1) * C(3)])
    # expect1 and expect2 measure each value through its expect, on all sites or chosen ones, where
    # in 1.x a representation of one's own defined its own, and so does measure for an operator of
    # one site and a pair
    @test expect1(wide, N) ≈ expect1(st, N)
    @test expect1(wide, [N, C]) ≈ expect1(st, [N, C])
    @test expect2(wide, (dag(C), C)) ≈ expect2(st, (dag(C), C))
    @test expect2(wide, [(N, N), (C, C)]) ≈ expect2(st, [(N, N), (C, C)])
    m = expect2(wide, (dag(C), C); sites = [1, 3])
    @test count(ismissing, m) == 5 && m[3, 1] ≈ expect(st, dag(C)(3) * C(1))
    @test_throws "expect2 acts on site 9, which the system does not have" expect2(wide, (N, N);
                                                                                sites = [1, 9])
    vals = measure(wide, [N, (dag(C), C), (N, N)])
    @test first.(vals) == ["N", "(dag(C), C)", "(N, N)"]
    @test last(vals[1]) isa Vector{Float64} && last(vals[1]) ≈ expect1(st, N)
    @test last(vals[2]) isa Matrix{ComplexF64} && last(vals[2]) ≈ expect2(st, (dag(C), C))
    @test last(vals[3]) isa Matrix{Float64} && last(vals[3]) ≈ expect2(st, (N, N))
    # a tuple of operators gives a tuple of values, on a State as on a state of one's own
    pair = (N(2), C(1) * dag(C)(3))
    @test expect(st, pair) isa Tuple
    @test collect(expect(wide, pair)) ≈ collect(expect(st, pair))
    mktempdir() do dir
        file = joinpath(dir, "wide.h5")
        save_state(file, "wide", wide)
        # read back on a system of the same sites, which it compares with (decision 15)
        back = load_state(file, "wide")
        @test back isa InterleavedState && back.system == sys
        @test [ expect(back, o) for o in ops ] ≈ [ expect(st, o) for o in ops ]
    end
end

# An evolution of one's own: each step by the exponential of the dense matrix of the MPO of the
# evolver, its time functions read at the middle of the step, exact on a small system
struct DenseStep <: EvolutionAlgorithm end

function TensorMixedStates.evolve_step(::DenseStep, pre::PreMPO, coefs, state::State{R}, time, dt,
                                       step, limits) where R
    mpo = isnothing(coefs) ? make_mpo(pre) : make_mpo(pre, [ f(time + dt / 2) for f in coefs ])
    s = siteinds(R, state.system)
    u = exp(dt * prod(mpo), [ i' for i in s ], s)
    return State(state, TensorMixedStates.ITensorMPS.MPS(IT.noprime(u * prod(state.state)), s))
end

# one that records what each step is given and multiplies the state by `factor`
struct StepRecorder <: EvolutionAlgorithm
    calls::Vector{Any}
    factor::Float64
end

StepRecorder(factor = 1.) = StepRecorder([], factor)

function TensorMixedStates.evolve_step(a::StepRecorder, pre, coefs, state, time, dt, step, limits)
    push!(a.calls, (time, dt, step, limits.maxdim, isnothing(coefs) ? nothing : [ f(time) for f in coefs ]))
    return a.factor * state
end

# one that has no step
struct Stepless <: EvolutionAlgorithm end

# A search of one's own: the lowest eigenvector of the dense matrix of the MPO
struct DenseSearch <: SearchAlgorithm end

function TensorMixedStates.search(::DenseSearch, mpo, state::State{R}, limits; sweep_done = nothing,
                                  krylov = Krylov(), done = 0, energy = nothing) where R
    LA = TensorMixedStates.LinearAlgebra
    s = siteinds(R, state.system)
    d = prod(IT.dim.(s))
    f = LA.eigen(LA.Hermitian(reshape(Array(prod(mpo), [ i' for i in s ]..., s...), d, d)))
    st = TensorMixedStates.ITensorMPS.MPS(IT.ITensor(f.vectors[:, 1], s...), s)
    return (f.values[1], State(state, st))
end

@testset "An algorithm of one's own" begin
    # an evolution of one's own goes through evolve and thermalize, which give each step the
    # PreMPO of the evolver, its time functions, the time, the step, its number and its limits.
    # Tdvp, at full bond dimension, is exact on three sites, see evolve.jl
    sys = System(3, Qubit())
    h = X(1) * X(2) + Z(2) * Z(3) + 0.5 * sum(Z(i) for i in 1:3)
    ψ = State{Pure}(sys, ["Up", "+", "Dn"])
    lim = Limits(cutoff = 0)
    for (st, ev) in [(ψ, -im * h), (mix(ψ), -im * h + Dissipator(sqrt(0.3) * Sm)(1))]
        @test norm(evolve(ev, 0.5, st; algo = DenseStep()) - evolve(ev, 0.5, st; nsteps = 4, limits = lim)) <
              1e-10
    end
    # and an Evolve phase, resumed after each step as the uninterrupted run
    mktempdir() do dir
        cd(dir) do
            phases = [CreateState{Pure}(3, Qubit(), ["Up", "+", "Dn"]),
                      Evolve(evolver = -im * h, duration = 0.4, time_step = 0.1, algo = DenseStep(),
                             measurements = "data" => [X(1), Z])]
            runTMS(SimData(; name = "ref", phases))
            while stopped(runTMS(SimData(; name = "chk", phases); max_time = -1))
            end
            @test read("chk/data", String) == read("ref/data", String)
        end
    end
    hs = [-im * X(1) * X(2), -im * sum(Z(i) for i in 1:3)]
    fs = [t -> 1., t -> cos(2t)]
    @test norm(evolve(hs => fs, 0.5, ψ; algo = DenseStep(), nsteps = 5) -
               evolve(hs => fs, 0.5, ψ; nsteps = 5, limits = lim)) < 1e-10
    # thermalize, against e^{-βH}/Z
    x, z, o = [0. 1. ; 1. 0.], [1. 0. ; 0. -1.], [1. 0. ; 0. 1.]
    hd = kron(x, x, o) + kron(o, z, z) + 0.5 * (kron(z, o, o) + kron(o, z, o) + kron(o, o, z))
    g = exp(-1.2 * hd)
    trg = sum(g[k, k] for k in 1:8)
    ρ0 = State{Mixed}(sys, "FullyMixed")
    l, ρ = thermalize(h, 1.2, ρ0; algo = DenseStep(), nsteps = 2)
    @test l ≈ log(trg / 8)
    @test norm(ρ - dense_state(sys, g / trg)) < 1e-10
    # what each step is given
    rec = StepRecorder()
    st = evolve(-im * h, 0.75, ψ; algo = rec, nsteps = 3, time = 1., step = 5,
                limits = Limits(maxdim = [1, 2, 3, 4, 5, 6, 7]))
    @test rec.calls == [(1., 0.25, 5, 5, nothing), (1.25, 0.25, 6, 6, nothing), (1.5, 0.25, 7, 7, nothing)]
    @test norm(st - ψ) < 1e-14
    rec = StepRecorder()
    evolve([-im * h] => [t -> 2t], 0.5, ψ; algo = rec, nsteps = 2)
    @test rec.calls == [(0., 0.25, 1, nothing, [0.]), (0.25, 0.25, 2, nothing, [0.5])]
    # by thermalize, steps of -dβ/2 from a time 0, the state brought back to trace one and its
    # trace added to the logarithm
    rec = StepRecorder(2.)
    l, ρ = thermalize(h, 0.75, ρ0; algo = rec, nsteps = 3, step = 2)
    @test rec.calls == [ (0., -0.125, k, nothing, nothing) for k in 2:4 ]
    @test l ≈ 3 * log(2)
    @test norm(ρ - ρ0) < 1e-14
    @test_throws "the trace of the state reached -1.0 at step 2 of thermalize: take smaller steps" thermalize(
        h, 0.75, ρ0; algo = StepRecorder(-1.), step = 2)
    @test_throws MethodError evolve(-im * h, 0.5, ψ; algo = Stepless())

    # a search of one's own goes through ground_state and steady_state
    LA = TensorMixedStates.LinearAlgebra
    e, gs = ground_state(h, ψ; algo = DenseSearch())
    @test e ≈ minimum(LA.eigvals(LA.Symmetric(hd)))
    @test real(expect(gs, h)) ≈ e
    @test norm(gs) ≈ 1
    Ω, γ = 0.8, 0.5
    L = sum(-im * Ω / 2 * X(i) + Dissipator(sqrt(γ) * Sm)(i) for i in 1:2)
    v, ρ = steady_state(L, State{Mixed}(System(2, Qubit()), "FullyMixed"); algo = DenseSearch())
    @test abs(v) < 1e-12
    @test trace(ρ) ≈ 1
    @test real(expect1(ρ, Z)) ≈ fill(-γ^2 / (γ^2 + 2Ω^2), 2)
end
