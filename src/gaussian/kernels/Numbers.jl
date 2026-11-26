# Complex numbers held by their logarithm, and pfaffians computed as such: the overlaps of
# states of hundreds of modes are products of hundreds of factors, which leave the range of
# Float64 (lot 0, E6: |Pf| ~ 1e-800 at 2n = 400).

"""
    LogComplex(logabs, phase)
    LogComplex(z)

the complex number `exp(logabs) * phase`, `phase` normalised to modulus 1, zero when
`logabs == -Inf` or `phase == 0`.
"""
struct LogComplex
    logabs::Float64
    phase::ComplexF64
    function LogComplex(logabs::Real, phase::Number)
        if logabs == -Inf || phase == 0
            return new(-Inf, zero(ComplexF64))
        end
        if !isfinite(logabs)
            throw(DomainError(logabs, "the logarithm of a LogComplex is finite or -Inf"))
        end
        new(logabs, phase / abs(phase))
    end
end

LogComplex(z::Number) = z == 0 ? zero(LogComplex) : LogComplex(log(abs(z)), z)
LogComplex(z::LogComplex) = z

Base.zero(::Type{LogComplex}) = LogComplex(-Inf, 0)
Base.one(::Type{LogComplex}) = LogComplex(0.0, 1)
Base.iszero(a::LogComplex) = a.logabs == -Inf
Base.ComplexF64(a::LogComplex) = exp(a.logabs) * a.phase
Base.abs(a::LogComplex) = exp(a.logabs)

Base.:*(a::LogComplex, b::LogComplex) = LogComplex(a.logabs + b.logabs, a.phase * b.phase)
Base.:*(a::LogComplex, z::Number) = a * LogComplex(z)
Base.:*(z::Number, a::LogComplex) = LogComplex(z) * a
Base.conj(a::LogComplex) = LogComplex(a.logabs, conj(a.phase))
Base.:-(a::LogComplex) = LogComplex(a.logabs, -a.phase)

function Base.inv(a::LogComplex)
    if iszero(a)
        throw(DivideError())
    end
    LogComplex(-a.logabs, conj(a.phase))
end
Base.:/(a::LogComplex, b::LogComplex) = a * inv(b)
Base.:/(a::LogComplex, z::Number) = a / LogComplex(z)

function Base.:+(a::LogComplex, b::LogComplex)
    if iszero(a)
        return b
    elseif iszero(b)
        return a
    end
    m = max(a.logabs, b.logabs)
    v = exp(a.logabs - m) * a.phase + exp(b.logabs - m) * b.phase
    LogComplex(m + log(abs(v)), v)
end
Base.:-(a::LogComplex, b::LogComplex) = a + (-b)

"""
    logpfaffian(A)

the pfaffian of the antisymmetric part of `A` as a `LogComplex`, `A` being antisymmetric up to
rounding.
"""
function logpfaffian(A::AbstractMatrix)
    if isempty(A)
        return one(LogComplex)
    end
    l, s = parlett_reid((A - transpose(A)) / 2)
    LogComplex(l, s)
end

# (log|Pf(A)|, phase of Pf(A)) for A antisymmetric, by the Parlett-Reid elimination with partial pivoting
# (Wimmer, ACM TOMS 38, 30 (2012), PFAPACK). The updates keep A exactly antisymmetric, so that a column that
# vanishes gives Pf = 0 (SkewLinearAlgebra 1.1, with fast math and the pivot tested on the other side of the
# diagonal, gave NaN on the overlap of two orthogonal states)
function parlett_reid(A0::AbstractMatrix{<:Number})
    n = size(A0, 1)
    A = Matrix{float(eltype(A0))}(A0)
    phase = one(complex(eltype(A)))
    if isodd(n)
        return -Inf, zero(phase)
    end
    logabs = 0.0
    τ = similar(A, n)
    c = similar(A, n)
    @inbounds for k in 1:2:n-1
        # pivot: the largest element of row k beyond the diagonal
        kp = k + 1
        largest = abs(A[k, k+1])
        for j in k+2:n
            if abs(A[k, j]) > largest
                kp, largest = j, abs(A[k, j])
            end
        end
        if kp != k + 1
            for j in k:n
                A[k+1, j], A[kp, j] = A[kp, j], A[k+1, j]
            end
            for i in k:n
                A[i, k+1], A[i, kp] = A[i, kp], A[i, k+1]
            end
            phase = -phase
        end
        p = A[k, k+1]
        if iszero(p)
            return -Inf, zero(phase)
        end
        logabs += log(abs(p))
        phase *= p / abs(p)
        # A[i, j] += τ_i c_j - c_i τ_j: the same rounding on (i, j) and (j, i), which stay opposite
        for i in k+2:n
            τ[i] = A[k, i] / p
            c[i] = A[i, k+1]
        end
        for j in k+2:n, i in k+2:n
            A[i, j] += τ[i] * c[j] - c[i] * τ[j]
        end
    end
    logabs, phase
end

"""
    logpfaffian(A, blocks)

the pfaffian of `A` by eliminating in turn its leading diagonal blocks of the sizes `blocks`,
with Pf([P B; -Bᵀ C]) = Pf(P) Pf(C + Bᵀ P⁻¹ B); an error if one of these blocks is singular.
"""
function logpfaffian(A::AbstractMatrix, blocks)
    pf = one(LogComplex)
    for k in blocks
        P = A[1:k, 1:k]
        B = A[1:k, k+1:end]
        p = logpfaffian(P)
        if iszero(p)
            throw(ArgumentError("a block eliminated by logpfaffian is singular"))
        end
        pf *= p
        A = A[k+1:end, k+1:end] + transpose(B) * (P \ B)
    end
    pf * logpfaffian(A)
end
