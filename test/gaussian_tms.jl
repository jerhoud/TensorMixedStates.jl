# Dense matrices of the operators of TMS for the tests of the Gaussian backend: the matrices of the
# sites with the Jordan-Wigner strings laid by hand, and the change from the basis of TMS to that of
# the modes of the backend (axes of the qubits, order of the modes).

using .LA: I, eigen, Hermitian

# the site matrix m on site k of sys, after F on the sites before it when m is odd
function gc_embed(sys, k, m)
    p = fermion_parities(sys.sites[k])
    odd = any(m[a, b] != 0 && isodd(p[a] + p[b]) for a in axes(m, 1), b in axes(m, 2))
    mats = [j < k && odd ? TensorMixedStates.matrix(TensorMixedStates.F, s) :
            j == k ? m : Matrix{ComplexF64}(I, dim(s), dim(s)) for (j, s) in enumerate(sys.sites)]
    reduce(kron, mats)
end
gc_op(sys, O, k) = gc_embed(sys, k, TensorMixedStates.matrix(O, sys.sites[k]))

# the kind of modes of a site, as in the compiler
gc_kind(s) = !isnothing(fermion_species(s)) ? :fermion : dim(s) == 2 ? :qubit : nothing

# the unitary from the occupation basis of a qubit along a to the basis of the site: |0⟩ the +1
# eigenvector of σ_a, |1⟩ = σ_b|0⟩, so that σ_b, σ_c, σ_a read X, Y, Z there
function gc_axis_unitary(a)
    a = a / LA.norm(a)
    σb, _ = TensorMixedStates.qubit_majoranas(a)
    E = eigen(Hermitian(TensorMixedStates.pauli(a)))
    v0 = E.vectors[:, argmax(E.values)]
    hcat(v0, σb * v0)
end

# perm such that a vector in the order of the modes (mode 1 the slowest) is x[perm] for x in the
# order of TMS; the modes are taken block by block, a block per kind in the order the kinds appear;
# an Electron has its modes ↑ then ↓, and its index n↑ + 2n↓ in TMS
function gc_to_modes(sys)
    n = length(sys.sites)
    kinds = [gc_kind(s) for s in sys.sites]
    order = [k for b in unique(filter(!isnothing, kinds)) for k in 1:n if kinds[k] == b]
    nm = [kinds[k] == :fermion ? length(fermion_species(sys.sites[k])) : 1 for k in 1:n]
    M = sum(nm[order])
    perm = zeros(Int, 2^M)
    for idx in 0:2^M-1
        bits = dense_bits(idx, M)
        loc = zeros(Int, n)
        pos = 0
        for k in order
            b = bits[pos+1:pos+nm[k]]
            pos += nm[k]
            loc[k] = nm[k] == 1 ? Int(b[1]) : Int(b[1]) + 2 * Int(b[2])
        end
        tms = 0
        for k in 1:n
            tms = tms * dim(sys.sites[k]) + loc[k]
        end
        perm[idx+1] = tms + 1
    end
    perm
end

