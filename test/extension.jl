# Extending the package from a package of one's own.
#
# Goes here: what a package built on TensorMixedStates declares at its top level. That code
# runs while the package is precompiled, and only what Julia keeps of that run is there once
# the package is loaded, so these tests write a package, have it precompiled and load it in a
# process of its own: in the process of the tests, the same code would run when included and
# show nothing.

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
            println(check(() -> real(expect(State{Pure}(System(2, Probe()), ["a", "b"]), Pz(2))) ≈ -1))
            println(check(() -> matrix(SiteProbe.Probed, Qubit()) == [0. 0. ; 1. 0.]))
            """, dir)
        @test out == "true\ntrue\ntrue\ntrue\n"
    end
end
