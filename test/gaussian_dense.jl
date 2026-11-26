# An exact evaluator of fermionic states as vectors of 2ⁿ amplitudes, written from the
# definitions and apart from GaussianKernels, to check it: the Majoranas are Pauli strings by
# Jordan-Wigner, |1⟩ is an occupied mode, and mode 1 is the most significant bit of the index.

const dense_X = ComplexF64[0 1; 1 0]
const dense_Y = ComplexF64[0 -im; im 0]
const dense_Z = ComplexF64[1 0; 0 -1]

dense_nmodes(ψ) = trailing_zeros(length(ψ))
dense_index(x) = 1 + sum(x[j] << (length(x) - j) for j in eachindex(x); init = 0)
dense_bits(k, n) = [((k >> (n - j)) & 1) == 1 for j in 1:n]

function dense_fock(x)
    ψ = zeros(ComplexF64, 2^length(x))
    ψ[dense_index(x)] = 1
    ψ
end

# the 2×2 matrix O on mode k
function dense_site(O, k, ψ)
    n = dense_nmodes(ψ)
    T = reshape(ψ, 2^(n - k), 2, 2^(k - 1))
    R = zeros(ComplexF64, size(T))
    for i in 1:2, l in 1:2
        R[:, i, :] .+= O[i, l] .* view(T, :, l, :)
    end
    vec(R)
end

function dense_majorana(a, ψ)
    j = (a + 1) ÷ 2
    φ = dense_site(isodd(a) ? dense_X : dense_Y, j, ψ)
    for k in 1:j-1
        φ = dense_site(dense_Z, k, φ)
    end
    φ
end

# M_ab = -i ⟨γ_a γ_b⟩ for a ≠ b, of a normalised ψ
function dense_covariance(ψ)
    m = 2 * dense_nmodes(ψ)
    φ = [dense_majorana(b, ψ) for b in 1:m]
    M = zeros(m, m)
    for a in 1:m, b in 1:m
        if a != b
            v = -im * dot(φ[a], φ[b])
            @assert abs(imag(v)) < 1e-12
            M[a, b] = real(v)
        end
    end
    M
end

# random rotations exp(θ γ_a γ_b / 2) = cos(θ/2) + sin(θ/2) γ_a γ_b of n modes, as (a, b, θ)
function dense_rotations(n, rng; nrot = 4 * n^2 + 4)
    rots = Tuple{Int, Int, Float64}[]
    for _ in 1:nrot
        a, b = rand(rng, 1:2n), rand(rng, 1:2n)
        if a != b
            push!(rots, (a, b, 2π * rand(rng)))
        end
    end
    rots
end

function dense_rotate(ψ, rots)
    for (a, b, θ) in rots
        ψ = cos(θ / 2) .* ψ .+ sin(θ / 2) .* dense_majorana(a, dense_majorana(b, ψ))
    end
    ψ
end

# a random pure Gaussian state with a random phase
dense_gaussian(x, rng) =
    dense_rotate(dense_fock(x) .* cis(2π * rand(rng)), dense_rotations(length(x), rng))

# c(α)ψ = γ_{α_1} ⋯ γ_{α_w} ψ
function dense_monomial(α, ψ)
    for b in Iterators.reverse(α)
        ψ = dense_majorana(b, ψ)
    end
    ψ
end

# the pfaffian by expansion along the first row
function dense_pfaffian(A)
    m = size(A, 1)
    if m == 0
        return one(eltype(A))
    elseif isodd(m)
        return zero(eltype(A))
    end
    s = zero(eltype(A))
    for j in 2:m
        rest = [k for k in 2:m if k != j]
        s += (-1)^j * A[1, j] * dense_pfaffian(A[rest, rest])
    end
    s
end

# the matrix of a linear map f on vectors of n modes
function dense_matrix(f, n)
    A = zeros(ComplexF64, 2^n, 2^n)
    for k in 1:2^n
        e = zeros(ComplexF64, 2^n)
        e[k] = 1
        A[:, k] = f(e)
    end
    A
end

dense_pair(M, a, b, ψ) = M[a, b] .* dense_majorana(a, dense_majorana(b, ψ))

# H = (i/4) Σ h_ab γ_a γ_b
dense_hamiltonian(h) =
    dense_matrix(ψ -> (im / 4) .* sum(dense_pair(h, a, b, ψ) for a in axes(h, 1), b in axes(h, 2) if a != b),
                 size(h, 1) ÷ 2)

# the state of a description: the top eigenvector of Σ_{a<b} M_ab (-i γ_a γ_b), of eigenvalue n,
# with the phase of its amplitude on g.x
function dense_vector(g)
    M = g.cov
    m = size(M, 1)
    O = dense_matrix(ψ -> -im .* sum(dense_pair(M, a, b, ψ) for a in 1:m, b in 1:m if a < b), m ÷ 2)
    v = eigen(Hermitian(O)).vectors[:, end]
    k = dense_index(g.x)
    v .* cis(angle(ComplexF64(g.amp)) - angle(v[k]))
end
