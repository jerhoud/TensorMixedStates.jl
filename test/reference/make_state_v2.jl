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
# The states exercise what version 2 added to version 1: conserved quantities, weak and
# strong, on one species and on two, a partial trace of a charged state, and site fields of
# every kind a file accepts, on `Kindly`, the site of `test/states_io.jl`.

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

spread(sys) = (State{Pure}(sys, ["Occ", "Emp", "Occ", "Emp"]) +
               0.5 * State{Pure}(sys, ["Emp", "Occ", "Occ", "Emp"])) / sqrt(1.25)
ψ = spread(System(4, Fermion(conserve = N)))
save_state(file, "fermions", ψ)
save_state(file, "fermions_mixed", mix(ψ))
save_state(file, "fermions_strong", mix(spread(System(4, Fermion(conserve = strong(N))))))
save_state(file, "electrons", State{Pure}(System(2, Electron(conserve = (Ntot, 2Sz))), ["Up", "Dn"]))
save_state(file, "partial_trace", partial_trace(mix(ψ), [1, 3]))
save_state(file, "kindly",
           State{Pure}(System([Kindly(:none, true, "a", 3, 0.5, nothing), Qubit()]), ["1", "Up"]))

println("wrote $file")
println("  N             : ", real.(expect1(ψ, N)))
println("  c1†c2         : ", expect(ψ, dag(C)(1) * C(2)))
