# Non-unitary Gaussian operations on a pure state of one component, held by its covariance: steps under
# a complex quadratic generator and linear operators, with the factor they put on the squared norm. The
# state is handled through its annihilators b = z·γ, z in the range of P₋ = (I - iM)/2, as the monitored
# free fermions of Cao, Tilloy and De Luca (SciPost Phys. 7, 024 (2019)) and Le Gal, Turkeshi and Schirò
# (PRX Quantum 5, 030329 (2024), appendix F) are held by their occupied orbitals, here in Majoranas.

# an orthonormal basis Z (2n × n) of the coefficients of the annihilators of the pure state of covariance M,
# the eigenvectors of iM of eigenvalue -1
function annihilators(M::AbstractMatrix{<:Real})
    E = eigen(Hermitian(im .* M))
    E.vectors[:, E.values .< 0]
end

# the covariance of the pure state whose annihilators have the coefficients spanned by the orthonormal Q,
# real when Q spans an isotropic subspace (QᵀQ = 0); an error when it does not to 1e-10, the precision of
# the steps lost
function covariance_of(Q::AbstractMatrix)
    iM = im .* (2 .* (Q * Q') - I)
    defect = maximum(abs, imag(iM); init = 0.0)
    if defect > 1e-10
        throw(ErrorException("a non-unitary Gaussian step lost its precision (defect $(round(defect; sigdigits = 2)) " *
                             "of the covariance): take smaller steps"))
    end
    antisymmetric(real(iM))
end

"""
    nonunitary_pieces(k, t)

the number of pieces of a step e^{tK}, K = (1/4) Σ k_ab γ_a γ_b, in `nonunitary_step`: each grows a vector
by at most e, ‖e^{(t/m)k}‖ ≤ e^{(t/m)‖Im k‖}, so that the QR of each loses little.
"""
nonunitary_pieces(k::AbstractMatrix, t::Real) = max(1, ceil(Int, abs(t) * opnorm(imag(k))))

"""
    nonunitary_step(M, R, m = 1)

the covariance of e^{mK}|ψ⟩, normalised, and log ‖e^{mK}|ψ⟩‖² for the normalised pure state ψ of covariance
`M`, K = (1/4) Σ k_ab γ_a γ_b with k complex antisymmetric and R = exp(k) complex orthogonal (`propagator`
for k real), applied m times: the annihilators Z of ψ become Q of RZ = QR each time, and the squared norm is
multiplied by Π_j |R_jj|. That factor, eq. (F12) of Le Gal, Turkeshi and Schirò for occupied orbitals of a
state of fixed number (Π |R_jj|²), is here the one of the annihilators of any pure Gaussian state (measured
against the dense exponential); R should grow vectors by a bounded factor, see `nonunitary_pieces`.
"""
function nonunitary_step(M::AbstractMatrix{<:Real}, R::AbstractMatrix, m::Int = 1)
    n = size(M, 1) ÷ 2
    Q = annihilators(M)
    lognorm = 0.0
    for _ in 1:m
        F = qr(R * Q)
        Q = Matrix(F.Q)[:, 1:n]
        lognorm += sum(log ∘ abs, diag(F.R))
    end
    covariance_of(Q), lognorm
end

# the part P₊u of u that creates on the pure state of covariance M, its annihilators Z: the rest annihilates
creating_part(Z::AbstractMatrix, u::AbstractVector) = u - Z * (Z' * u)

"""
    linear_annihilates(M, u)

whether the linear operator L = Σ_a u_a γ_a annihilates the pure state of covariance `M`, its creating
part P₊u within 1e-13 of u in norm, where rounding would make the direction of L|ψ⟩.
"""
linear_annihilates(M::AbstractMatrix{<:Real}, u::AbstractVector) =
    norm(creating_part(annihilators(M), u)) <= 1e-13 * norm(u)

"""
    linear_applied(M, u)

the covariance of L|ψ⟩, normalised, and ‖L|ψ⟩‖² = u†(I + iM)u = 2‖P₊u‖² for the normalised pure state ψ of
covariance `M` and the linear operator L = Σ_a u_a γ_a; `(nothing, 0.0)` when L annihilates ψ, see
`linear_annihilates`. Only the part w = P₊u of u creates: L|ψ⟩ = (w·γ)|ψ⟩, annihilated by w·γ and by the
annihilators z·γ of ψ with w̄†z = 0. Its precision is that of w, about eps ‖u‖ / ‖w‖.
"""
function linear_applied(M::AbstractMatrix{<:Real}, u::AbstractVector)
    Z = annihilators(M)
    w = creating_part(Z, u)
    if norm(w) <= 1e-13 * norm(u)
        return nothing, 0.0
    end
    p = 2 * norm(w)^2
    w ./= norm(w)
    # the annihilators of ψ orthogonal to w̄, the last ones of a basis of Z starting along Z†w̄
    a = Z' * conj(w)
    rest = (qr(reshape(a, :, 1)).Q * Matrix{ComplexF64}(I, length(a), length(a)))[:, 2:end]
    covariance_of(hcat(Z * rest, w)), p
end

"""
    lyapunov_propagator(X, Y, t)

`(E, G)` with M(t) = E M(0) Eᵀ + G the solution of dM/dt = X M + M Xᵀ + Y, the covariance of a mixed
Gaussian state under a Lindbladian of quadratic H and linear jumps (Bravyi and König 2012, eq. (32)):
E = e^{Xt} and G = ∫₀ᵗ e^{Xs} Y e^{Xᵀs} ds = F₁₂ Eᵀ, F₁₂ the upper right block of the exponential of
[X Y; 0 -Xᵀ] t (Van Loan 1978); the block e^{-Xᵀt} grows as e^{‖X‖t}, which `lyapunov_pieces` bounds, the
instability Bravyi and König note for long times.
"""
function lyapunov_propagator(X::AbstractMatrix{<:Real}, Y::AbstractMatrix{<:Real}, t::Real)
    n = size(X, 1)
    F = exp([X Y; zeros(n, n) -transpose(X)] * t)
    E = F[1:n, 1:n]
    E, F[1:n, n+1:end] * transpose(E)
end

"""
    lyapunov_pieces(X, t)

the number of pieces of a step of `lyapunov_propagator` over `t`, so that each grows its blocks by at most e.
"""
lyapunov_pieces(X::AbstractMatrix{<:Real}, t::Real) = max(1, ceil(Int, abs(t) * opnorm(X)))

"""
    purified_covariance(M)

the covariance [M S; -S -M] of a pure Gaussian state on twice the Majoranas whose first half has the
covariance `M` of a mixed state, S = √(I + M²): symmetric and a function of M², so that it commutes with M
and the result squares to -I; the eigenvalues of I + M² below zero by rounding are taken as zero. A step
on the first half, `nonunitary_step` of R ⊕ I, is then one of the mixed state, ρ → e^K ρ e^{K†}, the factor
on the norm its trace (the trace of the operator on the purifying half being one).
"""
function purified_covariance(M::AbstractMatrix{<:Real})
    E = eigen(Symmetric(I + M * M))
    S = E.vectors * Diagonal(sqrt.(max.(E.values, 0.0))) * transpose(E.vectors)
    [M S; -S -M]
end
