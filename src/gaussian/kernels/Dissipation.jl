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

# the covariance of the pure state whose annihilators have the coefficients spanned by the orthonormal Q
covariance_of(Q::AbstractMatrix) = antisymmetric(real(im .* (2 .* (Q * Q') - I)))

"""
    nonunitary_step(M, R)

the covariance of e^K|ψ⟩, normalised, and log ‖e^K|ψ⟩‖² for the normalised pure state ψ of covariance `M`
and K = (1/4) Σ k_ab γ_a γ_b with k complex antisymmetric, R = exp(k) complex orthogonal (`propagator` for
k real): with Z the annihilators of ψ and RZ = QR, those of e^K|ψ⟩ are Q and ‖e^K|ψ⟩‖² = Π_j |R_jj|,
the Majorana form of eq. (F12) of Le Gal, Turkeshi and Schirò.
"""
function nonunitary_step(M::AbstractMatrix{<:Real}, R::AbstractMatrix)
    n = size(M, 1) ÷ 2
    F = qr(R * annihilators(M))
    covariance_of(Matrix(F.Q)[:, 1:n]), sum(log ∘ abs, diag(F.R))
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
covariance `M` and the linear operator L = Σ_a u_a γ_a; an error when L annihilates ψ, see
`linear_annihilates`. Only the part w = P₊u of u creates: L|ψ⟩ = (w·γ)|ψ⟩, annihilated by w·γ and by the
annihilators z·γ of ψ with w̄†z = 0.
"""
function linear_applied(M::AbstractMatrix{<:Real}, u::AbstractVector)
    Z = annihilators(M)
    w = creating_part(Z, u)
    if norm(w) <= 1e-13 * norm(u)
        throw(ArgumentError("the linear operator annihilates the state"))
    end
    p = 2 * norm(w)^2
    w ./= norm(w)
    # the annihilators of ψ orthogonal to w̄, the last ones of a basis of Z starting along Z†w̄
    a = Z' * conj(w)
    rest = (qr(reshape(a, :, 1)).Q * Matrix{ComplexF64}(I, length(a), length(a)))[:, 2:end]
    covariance_of(hcat(Z * rest, w)), p
end
