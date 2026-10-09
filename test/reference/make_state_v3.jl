# Produces `state_v3.h5`, a state file in version 3 of the format, the one of version 2.
# `save_state` promises that a file it writes stays readable by every later version, and this
# file holds the reading of the current format to that promise.
#
# As for the other scripts of this directory, the suite does not run this one: it is here so
# that the file can be regenerated and argued with. Until 2.0.0 is released, it runs against the
# working tree, and the file is written again when the format changes:
#
#     julia --project=. test/reference/make_state_v3.jl
#
# The states exercise what version 3 added to version 2, what the sites conserve written typed
# apart from their other fields, and nothing more: a pure state on a site conserving a quantity
# next to another field, one conserving two, and one conserving nothing, and a mixed one
# conserving a parity strongly, whose numbers are read back intact.

using TensorMixedStates, .Bosons, .Electrons, .Fermions, .Qubits

file = get(ARGS, 1, joinpath(@__DIR__, "state_v3.h5"))
rm(file; force = true)

save_state(file, "charged",
           State{Pure}(System([Boson(4, conserve = N), Qubit(), Electron(conserve = (Ntot, 2Sz))]),
                       ["2", "Up", "Dn"]))
sys = System(2, Fermion(conserve = strong(parity(N))))
ψ = (State{Pure}(sys, ["Emp", "Emp"]) + 0.5 * State{Pure}(sys, ["Occ", "Occ"])) / sqrt(1.25)
save_state(file, "parity_strong", mix(ψ))

println("wrote $file")
