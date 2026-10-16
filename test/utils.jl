# Helpers shared by every test group. Nothing here runs a test by itself: this file is
# included before the groups, outside of the global testset, so that the macros it defines
# exist by the time the group files are parsed.

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
the vector of the pure state `st` in the full basis, site 1 most significant, not normalised
"""
raw_vector(st::State{Pure}) = vec(Array(IT.dense(prod(st.state)), reverse(siteinds(Pure, st.system))...))

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

"""
the density matrix of the mixed state `ρ`, its rows the kets and its columns the bras, site 1
most significant, read through the combiners its mixed indices are made with
"""
function ket_bra_matrix(ρ::State{Mixed})
    sys = ρ.system
    t = prod(ρ.state)
    js = siteinds(Pure, sys)
    bs = IT.Index[]
    for i in 1:length(sys)
        b, c = TMS.mixer(js[i], siteind(Mixed, sys, i), sys[i])
        t *= IT.dag(c)
        push!(bs, b')
    end
    d = prod(IT.dim.(js))
    return reshape(Array(IT.dense(t), reverse(js)..., reverse(bs)...), d, d)
end

"""
the dense matrix of the placed operator `op` on the whole of the system `sys`, from the terms
of its descent
"""
lowered(sys::System, op) = TMS.system_matrix(Pure, sys, TMS.lower(sys, op))

"""
the dense superoperator of the placed channel or generator `a` on the whole of the system
`sys`, from the terms of its descent
"""
superop(sys::System, a) = TMS.system_matrix(Mixed, sys, TMS.lower(sys, a))

"""
the dense matrix of the MPO of `op` on the sites of the pure state `st`
"""
function mpo_matrix(st::State{Pure}, op)
    idx = siteinds(Pure, st.system)
    a = Array(IT.dense(prod(make_mpo(st, op))), [ i' for i in reverse(idx) ]...,
              [ IT.dag(i) for i in reverse(idx) ]...)
    d = prod(IT.dim.(idx))
    return reshape(a, d, d)
end

"""
the least number of memory channels, those besides the term not yet begun and the term
finished, that a triangular MPO of the operator of dense matrix `h` can have on each link, the
sites having the dimensions `ds`: the rank of `h` across the link once its parts that are the
identity on the left or on the right are taken out
"""
function least_channels(h::AbstractMatrix, ds::Vector{Int})
    LA = TMS.LinearAlgebra
    n = length(ds)
    return map(1:n-1) do k
        dl, dr = prod(ds[1:k]), prod(ds[k+1:n])
        m = reshape(permutedims(reshape(h, dr, dl, dr, dl), (2, 4, 1, 3)), dl^2, dr^2)
        vl = vec(Matrix{eltype(m)}(LA.I, dl, dl)) / sqrt(dl)
        vr = vec(Matrix{eltype(m)}(LA.I, dr, dr)) / sqrt(dr)
        m = m - vl * (vl' * m)
        m = m - (m * vr) * vr'
        count(>(1e-10 * norm(h)), LA.svdvals(m))
    end
end
