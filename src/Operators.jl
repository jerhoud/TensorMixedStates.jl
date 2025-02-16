export Exps, ProdOf, SumOf, TensorOf, ExpOf, Operator, ⊗

abstract type Exps{T} end

struct ProdOf{T} <: Exps{T}
    coef::Number
    factors::Vector{<:Exps{T}}
end

(a::Exps{T} * b::Exps{T}) where T = ProdOf(1, [a, b])
(a::Number * b::Exps{T}) where T = ProdOf(a, [b])
(a::Exps{T} * b::Number) where T = b * a
(a::Exps{T} / b::Number) where T = inv(b) * a
-(a::Exps{T}) = -1 * a

struct SumOf{T} <: Exps{T}
    terms::Vector{<:Exps{T}}
end

(a::Exps{T} + b::Exps{T}) where T = SumOf([a, b])

struct TensorOf{T} <: Exps{T}
    factors::Vector{<:Exps{T}}
end

(a::Exps{T} ⊗ b::Exps{T}) where T = TensorOf([a, b])

struct ExpOf{T} <: Exps{T}
    arg::Exps{T}
end

exp(a::Exps{T}) where T = ExpOf(a)

struct BaseOp end
struct Operator <: Exps{BaseOp}
    name::String
end
