struct MixEvolve end
struct MixEvolve2 end
struct MixDissipator end

function make_operator(::Type{MixEvolve}, tj::ITensor, i::Index)
    j = noprime(ind(tj, 1))
    k = sim(j)
    tk = replaceinds(tj, (j, j'), (k, k'))
    return *(tj, delta(k, k'), combinerto(k, j, i), combinerto(k', j', i ))
end

function make_operator(::Type{MixEvolve2}, tj::ITensor, i::Index)
    j = noprime(ind(tj, 1))
    k = sim(j)
    tk = replaceinds(tj, (j, j'), (k, k'))
    return *(tk, delta(j, j'), combinerto(k, j, i), combinerto(k', j', i ))
end

function make_operator(::Type{MixDissipator}, tj::ITensor, i::Index)
    j = noprime(ind(tj, 1))
    k = sim(j)
    tk = replaceinds(tj, (j, j'), (k, k'))
    atj = swapprime(dag(tj'), 1=>2)
    atk = swapprime(dag(tk'), 1=>2)
    r = tj * dag(tk) -
        0.5 * *(replaceprime(atj * tj, 2 => 1), delta.(k, k')) -
        0.5 * *(replaceprime(atk * tk, 2 => 1), delta.(j, j'))
    return *(r, combinerto(k, j, i), combinerto(k', j', i'))
end





