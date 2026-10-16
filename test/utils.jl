# Helpers shared by every test group. Nothing here runs a test by itself: this file is
# included before the groups, outside of the global testset, so that the macros it defines
# exist by the time the group files are parsed.

# the package, ITensors and LinearAlgebra under short names, for what a test reaches that is
# not exported
const TMS = TensorMixedStates
const IT = TensorMixedStates.ITensors
const LA = TensorMixedStates.LinearAlgebra

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
the vector of the pure state `st` in the full basis, site 1 most significant, not normalised
"""
raw_vector(st::State{Pure}) = vec(Array(IT.dense(prod(st.state)), reverse(siteinds(Pure, st.system))...))

"""
the vector of the pure state `st` in the full basis, site 1 most significant, normalised. Only
usable on small systems.
"""
function dense_vector(st::State{Pure})
    v = raw_vector(st)
    return v / norm(v)
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
the matrix on the whole of `system`, site 1 most significant, of the matrix `m` of site `k`, its
odd part carrying `F` on every site before `k`
"""
function jw_placed(system, k, m)
    sites = system.sites
    f = matrix(F, sites[k])
    id(i) = matrix(Id, sites[i])
    laid(a, strings) = foldl(kron, [ i < k ? (strings ? matrix(F, sites[i]) : id(i)) :
                                     i == k ? a : id(i) for i in eachindex(sites) ])
    return laid((m + f * m * f) / 2, false) + laid((m - f * m * f) / 2, true)
end

"""
the matrix `m` of one qubit on site `i` of `n` qubits, site 1 most significant
"""
on_qubit(m, i, n) = kron([ k == i ? m : [1. 0. ; 0. 1.] for k in 1:n ]...)

"""
    mpo_matrix(m::MPO, idx)
    mpo_matrix(pre::PreMPO, coefs = [1.])
    mpo_matrix(st::State, m::MPO)
    mpo_matrix(st::State, op)

the dense matrix, site 1 most significant, of the MPO `m` on the indices `idx` of its sites, of
the MPO of `pre` for the values `coefs` of its time functions, or of the MPO `m` or the MPO of
`op` on the sites of `st`
"""
function mpo_matrix(m::TMS.ITensorMPS.MPO, idx)
    # the 16 indices of 8 sites, which ITensors warns about from 14
    t = IT.@disable_warn_order prod(m)
    a = Array(IT.dense(t), [ i' for i in reverse(idx) ]..., [ IT.dag(i) for i in reverse(idx) ]...)
    d = prod(IT.dim.(idx))
    return reshape(a, d, d)
end
mpo_matrix(pre::PreMPO{R}, coefs = [1.]) where R =
    mpo_matrix(make_mpo(pre, coefs), siteinds(R, pre.system))
mpo_matrix(st::State{R}, m::TMS.ITensorMPS.MPO) where R <: TMS.PM =
    mpo_matrix(m, siteinds(R, st.system))
mpo_matrix(st::State, op) = mpo_matrix(PreMPO(st, op))

"""
    in_tmp(f)

`f()` with a temporary directory as working directory
"""
in_tmp(f) = mktempdir(dir -> cd(f, dir))

"""
    run_through(sd; kwargs...)

the simulation `sd` run again and again until it is not stopped
"""
function run_through(sd; kwargs...)
    for _ in 1:500
        sim = runTMS(sd; kwargs...)
        if !stopped(sim)
            return sim
        end
    end
    error("the simulation did not end")
end

"""
the least number of memory channels, those besides the term not yet begun and the term
finished, that a triangular MPO of the operator of dense matrix `h` can have on each link, the
sites having the dimensions `ds`: the rank of `h` across the link once its parts that are the
identity on the left or on the right are taken out
"""
function least_channels(h::AbstractMatrix, ds::Vector{Int})
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
