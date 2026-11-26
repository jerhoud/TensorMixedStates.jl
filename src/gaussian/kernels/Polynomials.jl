# Polynomials in Majoranas, Σ_α c_α c(α) with c(α) = γ_{α_1} ⋯ γ_{α_w} and α strictly increasing.

"""
    MajoranaPoly(α => c, ...)
    MajoranaPoly(terms::AbstractDict)

the polynomial Σ c c(α) in Majoranas, each α a strictly increasing list of positive indices (the
empty list is the identity); pairs with the same α add up, and exactly zero coefficients are dropped.
"""
struct MajoranaPoly
    terms::Dict{Vector{Int}, ComplexF64}
    function MajoranaPoly(terms::AbstractDict)
        d = Dict{Vector{Int}, ComplexF64}()
        for (α, c) in terms
            if any(<(1), α) || any(α[k] >= α[k+1] for k in 1:length(α)-1)
                throw(ArgumentError("a monomial is a strictly increasing list of positive indices, not $α"))
            end
            if c != 0
                d[Vector{Int}(α)] = c
            end
        end
        new(d)
    end
end

function MajoranaPoly(ps::Pair...)
    d = Dict{Vector{Int}, ComplexF64}()
    for (α, c) in ps
        d[α] = get(d, α, 0) + c
    end
    MajoranaPoly(d)
end

Base.zero(::Type{MajoranaPoly}) = MajoranaPoly()
Base.iszero(p::MajoranaPoly) = isempty(p.terms)
Base.length(p::MajoranaPoly) = length(p.terms)
Base.:(==)(p::MajoranaPoly, q::MajoranaPoly) = p.terms == q.terms
Base.hash(p::MajoranaPoly, h::UInt) = hash(p.terms, h)

"""
    degree(p)

the largest number of Majoranas in a monomial of `p`, 0 for the zero polynomial.
"""
degree(p::MajoranaPoly) = maximum(length, keys(p.terms); init = 0)

"""
    parity(p)

0 if all monomials of `p` are even, 1 if all are odd, `nothing` if both occur; 0 for zero.
"""
function parity(p::MajoranaPoly)
    ps = unique(length(α) % 2 for α in keys(p.terms))
    isempty(ps) ? 0 : length(ps) == 1 ? only(ps) : nothing
end

# c(α) c(β) = sign c(α Δ β): each γ_b of β moves left past the γ_a of α with a > b, and γ_b γ_b = 1
function monomial_product(α::Vector{Int}, β::Vector{Int})
    crossings = sum(count(>(b), α) for b in β; init = 0)
    (iseven(crossings) ? 1 : -1), sort!(symdiff(α, β))
end

function Base.:+(p::MajoranaPoly, q::MajoranaPoly)
    d = copy(p.terms)
    for (β, c) in q.terms
        d[β] = get(d, β, 0) + c
    end
    MajoranaPoly(d)
end

Base.:*(z::Number, p::MajoranaPoly) = MajoranaPoly(Dict(α => z * c for (α, c) in p.terms))
Base.:*(p::MajoranaPoly, z::Number) = z * p
Base.:-(p::MajoranaPoly) = -1 * p
Base.:-(p::MajoranaPoly, q::MajoranaPoly) = p + (-q)

function Base.:*(p::MajoranaPoly, q::MajoranaPoly)
    d = Dict{Vector{Int}, ComplexF64}()
    for (α, a) in p.terms, (β, b) in q.terms
        s, γ = monomial_product(α, β)
        d[γ] = get(d, γ, 0) + s * a * b
    end
    MajoranaPoly(d)
end

# c(α)† = γ_{α_w} ⋯ γ_{α_1} = (-1)^{w(w-1)/2} c(α)
Base.adjoint(p::MajoranaPoly) =
    MajoranaPoly(Dict(α => (iseven(length(α) * (length(α) - 1) ÷ 2) ? 1 : -1) * conj(c) for (α, c) in p.terms))
