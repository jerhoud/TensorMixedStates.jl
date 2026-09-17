using TensorMixedStates, .Qubits, .Qudits, .Fermions, .Bosons, .Spins, .Electrons, .Tjs, .Qbosons
using Test, Aqua, DataFrames

include("utils.jl")

"""
the groups of tests, each one lives in the file of the same name, whose header comment
says what belongs there. To add a group, create the file and add its name here.
"""
const GROUPS = ["building", "operators", "sites", "observables", "states_io", "evolve", "algorithms"]

# without arguments everything runs, otherwise only the groups given on the command
# line, as in `Pkg.test(test_args = ["observables", "sites"])`
selected = isempty(ARGS) ? GROUPS : ARGS
for g in selected
    g in GROUPS || error("unknown test group \"$g\", expected one of $(join(GROUPS, ", "))")
end

@testset verbose=true "TensorMixedStates.jl" begin
    for g in selected
        include("$g.jl")
    end
end
