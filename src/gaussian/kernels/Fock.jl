# Fock states |x⟩, x[j] the occupation of mode j.

"""
    fock_covariance(x)

the covariance of the Fock state |x⟩.
"""
function fock_covariance(x::AbstractVector{Bool})
    n = length(x)
    M = zeros(2n, 2n)
    for j in 1:n
        s = x[j] ? -1.0 : 1.0
        M[2j-1, 2j] = s
        M[2j, 2j-1] = -s
    end
    M
end

"""
    majorana_on_fock(b, x)

`(coef, y)` such that γ_b |x⟩ = coef |y⟩.
"""
function majorana_on_fock(b::Integer, x::AbstractVector{Bool})
    j = (b + 1) ÷ 2
    s = iseven(count(view(x, 1:j-1))) ? 1.0 : -1.0
    coef = isodd(b) ? complex(s) : s * im * (x[j] ? -1 : 1)
    y = Vector{Bool}(x)
    y[j] = !y[j]
    coef, y
end
