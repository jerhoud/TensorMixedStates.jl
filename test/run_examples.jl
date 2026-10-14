# Runs the examples of `examples/high_level` and `examples/low_level`, to catch the day one of
# them stops working, or, the low level ones checking their results against exact values, starts
# giving a wrong result.
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
# instead of once per example. They write their output directories in the working directory, so
# the whole thing happens in a temporary directory.

using TensorMixedStates
using Test

const EXAMPLES = [
    "high_level" => ["dmrg", "gates", "precession", "ising_quench", "complete_graph_tdvp",
                     "fermion_chain_conserved", "excited_states", "measured_trajectory",
                     "thermal_phase"],
    "low_level" => ["spin_one_chain", "free_fermions", "domain_wall_quench",
                    "driven_damped_qubits", "thermal_state", "noisy_ghz", "entanglement",
                    "strong_then_weak"],
]
const EXAMPLE_DIR = normpath(joinpath(@__DIR__, "..", "examples"))

cd(mktempdir()) do
    @testset verbose = true "examples" begin
        for (dir, names) in EXAMPLES, name in names
            @testset "$dir/$name" begin
                t = @elapsed @eval module $(Symbol("Example_", name))
                    include($(joinpath(EXAMPLE_DIR, dir, name * ".jl")))
                end
                # reached only if the example ran to the end without throwing
                @test true
                println("    $dir/$name ran in $(round(t; digits = 1)) s")
            end
        end
    end
end
