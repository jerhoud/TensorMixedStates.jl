# Ground states of quadratic Hamiltonians H = (i/4) Σ h_ab γ_a γ_b, from the normal form of h:
# in the real Schur form h = Q T Qᵀ, a block [0 ε; -ε 0] on columns (i, j) of Q is a mode of energy
# ε (H = Σ ε (n - 1/2)), empty in the ground state; the real eigenvalues, zero to rounding, are
# paired into modes of zero energy, which any occupation fills.

"""
    quadratic_ground_state(h, parity = nothing)

`(cov, energy, energies, occupied)`: the covariance and energy of a ground state of H = (i/4) Σ h_ab γ_a γ_b,
of the given parity (±1) or of any, the energies of its normal modes in increasing order, and whether
its lowest mode is occupied (to reach the parity).
"""
function quadratic_ground_state(h::AbstractMatrix{<:Real}, parity::Union{Nothing, Integer} = nothing)
    if !(parity === nothing || parity == 1 || parity == -1)
        throw(ArgumentError("a parity is 1 or -1, not $parity"))
    end
    if size(h, 1) != size(h, 2) || isodd(size(h, 1))
        throw(ArgumentError("h is a square matrix of size 2n for n modes, not $(size(h))"))
    end
    if isempty(h) && parity == -1
        throw(ArgumentError("a state of no mode has parity 1"))
    end
    F = schur(Matrix{Float64}(antisymmetric(h)))
    T, Q = F.T, F.Z
    m = size(T, 1)
    modes = Tuple{Int, Int, Float64}[]
    for (a, b) in normal_pairs(T)
        if b == a + 1 && T[b, a] != 0
            ε = (T[a, b] - T[b, a]) / 2
            push!(modes, ε >= 0 ? (a, b, ε) : (b, a, -ε))
        else
            push!(modes, (a, b, max(abs(T[a, a]), abs(T[b, b]))))
        end
    end
    sort!(modes; by = last)
    M = zeros(m, m)
    for (a, b, _) in modes
        M .+= Q[:, a] .* transpose(Q[:, b]) .- Q[:, b] .* transpose(Q[:, a])
    end
    energies = [ε for (_, _, ε) in modes]
    energy = -sum(energies; init = 0.0) / 2
    occupied = parity !== nothing && real(logpfaffian(M).phase) * parity < 0
    if occupied
        # the other sector: the lowest mode occupied
        a, b, ε = first(modes)
        M .-= 2 .* (Q[:, a] .* transpose(Q[:, b]) .- Q[:, b] .* transpose(Q[:, a]))
        energy += ε
    end
    (cov = M, energy = energy, energies = energies, occupied = occupied)
end
