# Produces `state_v2.h5`, a state file in version 2 of the format, written by the released
# 1.6.0. `save_state` promises that a file it writes stays readable by every later version,
# and this file holds the reading of the current format to that promise.
#
# As for the other scripts of this directory, the suite does not run this one: it is here so
# that the file can be regenerated and argued with. It has to run against the released 1.6.0,
# not against the working tree:
#
#     git worktree add /tmp/tms-v160 v1.6.0
#     cp Manifest.toml /tmp/tms-v160/
#     julia --project=/tmp/tms-v160 test/reference/make_state_v2.jl
#     git worktree remove /tmp/tms-v160
#
# The states exercise what version 2 added to version 1, and nothing more, each one adding 20
# to 40 kB to the file, most of it HDF5 structure: a pure state conserving two quantities on a
# site, a mixed one conserving a quantity strongly, whose correlation checks that its numbers
# are read back intact, and site fields of every kind a file accepts, on `Kindly`, the site of
# `test/states_io.jl`.

using TensorMixedStates, .Fermions, .Electrons, .Qubits
strong = TensorMixedStates.strong

struct Kindly <: AbstractSite
    conserve::Symbol
    on::Bool
    label::String
    n::Int
    x::Float64
    void::Nothing
end

TensorMixedStates.dim(::Kindly) = 2

@def_states(Kindly(:none, true, "a", 3, 0.5, nothing), [ "1" => [0., 1.] ])

file = get(ARGS, 1, joinpath(@__DIR__, "state_v2.h5"))
rm(file; force = true)

sys = System(2, Fermion(conserve = strong(N)))
ρ = mix((State{Pure}(sys, ["Occ", "Emp"]) + 0.5 * State{Pure}(sys, ["Emp", "Occ"])) / sqrt(1.25))
save_state(file, "fermions_strong", ρ)
save_state(file, "electrons",
           State{Pure}(System(2, Electron(conserve = (Ntot, 2Sz))), ["Up", "Dn"]))
save_state(file, "kindly",
           State{Pure}(System([Kindly(:none, true, "a", 3, 0.5, nothing), Qubit()]), ["1", "Up"]))

println("wrote $file")
println("  N             : ", real.(expect1(ρ, N)))
println("  c1†c2         : ", expect(ρ, dag(C)(1) * C(2)))
