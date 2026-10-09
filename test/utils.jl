# Helpers shared by every test group. Nothing here runs a test by itself: this file is
# included before the groups, outside of the global testset, so that the macros it defines
# exist by the time the group files are parsed. The helpers of 1.x the later layers need are
# in v1/utils.jl until those layers come.

# the package and ITensors under short names, for what a test reaches that is not exported
const TMS = TensorMixedStates
const IT = TensorMixedStates.ITensors

"""
test whether a statement executes without throwing an exception
"""
macro test_ok(a)
    :(@test ($(esc(a)); true))
end

"""
test a statement (with `@test_ok`) with `type = Pure` and `type = Mixed`
"""
macro test_pm(a)
    tp = :type
    quote
        @test_ok (($(esc(tp)) -> $(esc(a)))(Pure))
        @test_ok (($(esc(tp)) -> $(esc(a)))(Mixed))
    end
end

"""
the vector of a pure state in the full basis, site 1 most significant, normalised. Only usable
on small systems.
"""
function dense_vector(state::State{Pure})
    idx = siteinds(Pure, state.system)
    psi = vec(Array(prod(state.state), reverse(idx)...))
    return psi / norm(psi)
end

"""
the density matrix of a mixed state of qubits in the full basis, site 1 most significant. Only
usable on small systems.
"""
function dense_matrix(ρ::State{Mixed})
    n = length(ρ.system)
    v = Array(prod(ρ.state), reverse(siteinds(Mixed, ρ.system))...)
    perm = vcat([2k - 1 for k in 1:n], [2k for k in 1:n])
    return reshape(permutedims(reshape(v, ntuple(_ -> 2, 2n)...), perm), 2^n, 2^n)
end
