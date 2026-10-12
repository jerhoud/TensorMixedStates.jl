# Extending the package from outside.
#
# Goes here: what a package built on TensorMixedStates relies on. A site type it declares at
# its top level is only checked by a package written for the occasion, precompiled and loaded
# in a process of its own: that code runs while the package is precompiled, and only what
# Julia keeps of that run is there once the package is loaded, while in the process of the
# tests the same code would run when included and show nothing. A representation of its own
# needs no such thing and is checked here directly.
#
# An assertion that needs a later layer is kept as `@test_skip`, marked with that layer. The
# algorithms of one's own of 1.x, which run in phases, are in v1/extension.jl.

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
    @test_skip "CreateState and Gates of a representation of one's own, measured and resumed"   # couche 5
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

@testset "An algorithm of one's own" begin
    @test_skip "tdvp step by step through run_steps, resumed as the uninterrupted run"   # couche 5
end
