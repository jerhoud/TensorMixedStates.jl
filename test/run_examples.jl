# Runs the examples of `examples/high_level`, to catch the day one of them stops working.
#
# This is not a test group: it is absent from `GROUPS` in `runtests.jl` and the suite does
# not run it. A step of its own in CI does, on one job of the matrix, because it costs
# minutes rather than seconds. Run it by hand with
#
#     julia --project=. test/run_examples.jl
#
# The examples of `examples/article` are deliberately left out. They are the ones of the
# companion article and run at the sizes it published, hours rather than minutes.
#
# Each example is included in a module of its own, so that two of them binding the same
# name do not collide, and all of them share one session, which pays the compilation once
# instead of five times. They write output directories and `runTMS` changes the working
# directory while it runs, so the whole thing happens in a temporary directory.

using TensorMixedStates
using Test

const EXAMPLES = ["dmrg", "gates", "precession", "ising_quench", "complete_graph_tdvp"]
const EXAMPLE_DIR = normpath(joinpath(@__DIR__, "..", "examples", "high_level"))

cd(mktempdir()) do
    @testset verbose = true "examples" begin
        for name in EXAMPLES
            @testset "$name" begin
                t = @elapsed @eval module $(Symbol("Example_", name))
                    include($(joinpath(EXAMPLE_DIR, name * ".jl")))
                end
                # reached only if the example ran to the end without throwing
                @test true
                println("    $name ran in $(round(t; digits = 1)) s")
            end
        end
    end
end
