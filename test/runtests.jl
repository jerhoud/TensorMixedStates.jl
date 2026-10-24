using TensorMixedStates, .Qubits, .Qudits, .Fermions, .Bosons, .Spins, .Electrons, .Tjs, .Qbosons,
      .Circuits
using Test, Random, Aqua

# Every source of randomness in the suite is pinned, so that a run can be reproduced and a
# failure cannot come and go. This covers `random_state`, which draws from the global
# generator by default. The statistical tests of `sampling.jl` use a generator of their own on top
# of this, so that they do not depend on how much randomness the groups before them drew.
# Note that a fixed seed pins the stream for one Julia version, not across versions.
Random.seed!(20260920)

# runTMS ends with a line of information, which hundreds of runs would make of the output of the
# tests: only the warnings are shown, the tests of that line catching it themselves
TensorMixedStates.Logging.global_logger(
    TensorMixedStates.Logging.ConsoleLogger(stderr, TensorMixedStates.Logging.Warn))

include("utils.jl")

"""
the groups of tests, each one lives in the file of the same name, whose header comment
says what belongs there. To add a group, create the file and add its name here.
"""
const GROUPS = ["building", "channels", "compact", "docs", "evolve", "extension", "fermions",
                "hints", "observables", "operators", "outputs", "prepared", "runner", "sampling",
                "searches", "sites", "states_io"]

# without arguments everything runs, otherwise only the groups given on the command
# line, as in `Pkg.test(test_args = ["observables", "sites"])`
selected = isempty(ARGS) ? GROUPS : ARGS
for g in selected
    if !(g in GROUPS)
        error("unknown test group \"$g\", expected one of $(join(GROUPS, ", "))")
    end
end

# one testset per group, which the summary shows as one line, unfolded by Test when one of
# its tests fails; a line as each group starts and ends shows the progress of a long run, a
# failure printed between those of its group
@testset verbose=true "TensorMixedStates.jl" begin
    for (k, g) in enumerate(selected)
        label = "[$k/$(length(selected))] $g"
        println(label)
        flush(stdout)
        t = time()
        @testset "$g" begin
            include("$g.jl")
        end
        println("$label: $(round(Int, time() - t)) s")
        flush(stdout)
    end
end
