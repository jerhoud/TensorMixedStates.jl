# Produces `state_v1.h5`, a state file in version 1 of the format, the one written by
# TensorMixedStates up to and including 1.3.0. `load_state` has to go on reading it, and a
# checkpoint left by such a version holds a state in that format, so a resume depends on it
# too.
#
# As for the other scripts of this directory, the suite does not run this one: it is here so
# that the file can be regenerated and argued with. It has to run against the released 1.3.0,
# not against the working tree, whose format is no longer version 1:
#
#     git worktree add /tmp/tms-v130 v1.3.0
#     cp Manifest.toml /tmp/tms-v130/
#     julia --project=/tmp/tms-v130 test/reference/make_state_v1.jl
#     git worktree remove /tmp/tms-v130
#
# The sites are chosen to exercise the parameters as version 1 wrote them, that is as
# `Float64`: a site with none, one with an integer, one with a float, and one with both.

using TensorMixedStates, .Qubits, .Bosons, .Spins, .Qbosons

file = get(ARGS, 1, joinpath(@__DIR__, "state_v1.h5"))
rm(file; force = true)

sites = [Qubit(), Boson(4), Spin(3/2), Qboson(0.1, 3)]
pure = State{Pure}(System(sites), ["Up", "2", "1/2", "1"])
save_state(file, "pure", pure)
save_state(file, "mixed", mix(pure))

println("wrote $file")
println("  sites        : ", sites)
println("  expect Z(1)  : ", expect(pure, Z(1)))
println("  expect N(2)  : ", expect(pure, N(2)))
println("  expect Sz(3) : ", expect(pure, Sz(3)))
