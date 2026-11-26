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
    l, s = logabspfaffian((A - transpose(A)) / 2)
    LogComplex(l, s)
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
