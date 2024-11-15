struct PreMPO
    linkdims::Vector{Int}
    links::Vector{Vector{Tuple{Int, Int, Int}}}
    tensors::Dict{Tuple{Int, Int}, ITensor}
end


function prepare_mpo(tp, a::SumLit, sites)
    n = lenght(sites)
    linkdims = fill(0, n + 1)
    links = [ Tuple{Int, Int, Int}[] for _ in 1:n ]
    for (i, p) in enumerate(a.ps)
        fst = p.ls[1].index[1]
        lst = p.ls[end].index[1]
        for j in fst + 1:lst
            linkdims[j] += 1
        end
        for j in fst:lst
            push!(links[j], (i, linkdims[j], linkdims[j+1]))
        end
    end
    tensors = Dict{Tuple{Int, Int}, ITensor}()
    for (nt, p) in enumerate(a.ps)
        i = 0
        dissip = false
        t = ITensor()
        idx = Index()
        c = p.coef
        for l in p.ls
            j = l.index[1]
            if i ≠ j
                if i ≠ 0
                    tensors[(nt, i)] = t
                end
                dissip = l.dissipator
                idx = sites[j]
                if tp == Mixed && !dissip 
                    idx = Index(dim(idx)÷2, tags(idx))
                end
                t = c * op(l.name, idx; l.param...)
                c = 1
                i = j
            elseif dissip || l.dissipator
                die("It is not allowed to multiply a dissipator by another operator")
            else
                t = replaceprime(t' * op(l.name, idx; l.param...), 2 => 1)
            end
        end
        if i ≠ 0
            tensors[(nt, i)] = t
        end
    end
    return PreMPO(linkdims, links, tensors)
end

struct MixEvolve end
struct MixGate end

function make_operator(::Type{MixEvolve}, tj::ITensor, i::Index, j::Index)
    k = sim(j)
    tk = replaceinds(tj, (j, j'), (k, k'))
    return (*(tj, delta(k, k'), combinerto(k, j, i) * combinerto(k', j', i )),
            *(tk, delta(j, j'), combinerto(k, j, i) * combinerto(k', j', i )))
end

function make_operator(::Type{GateEvolve}, tj::ITensor, i::Index, j::Index)
    k = sim(j)
    tk = replaceinds(tj, (j, j'), (k, k'))
    return *(tj, dag(tk), combinerto(k, j, i), combinerto(k', j', i'))
end

function mpo_abcd(tp, sl::SumLit, i::Int, sites, pre::PreMPO)
    idx = sites[i]
    n = length(sites)
    ldim = pre.linkdims[i]
    rdim = pre.linkdims[i + 1]
    ts = pre.tensors
    if tp == MixEvolve
        twice = true
        ldim *= 2
        rdim *= 2
    else
        twice = false
    end
    a = fill(ITensor(), (ldim, rdim))
    b = fill(ITensor(), ldim)
    c = fill(ITensor(), rdim)
    d = ITensor()

    for (term, rlink, llink) in pre.links[i]
        t = ts[(term, i)]
        tidx = noprime(t.inds[1])
        samedim = dim(tidx) == dim(idx)
        if samedim # for pure systems or dissipators
            t = replaceinds(t, (tidx, tidx'), (idx, idx'))
        else
            t = make_operator(tp, t, idx, tidx)
        end
        if rlink == 0
            if llink == 0
                d += t
            else
                if twice
                    b[2 * llink - 1] += t[1]
                    b[2 * llink] += t[2]
                else
                    b[llink] += t
                end
            end
        else
            if llink == 0
                if twice
                    c[2 * rlink - 1] += t[1]
                    c[2 * rlink] += t[2]      
                else
                    c[rlink] += t
                end
            else
                if twice
                    a[2 * llink - 1, 2 * rlink - 1] += t[1]
                    a[2 * llink, 2 * rlink] += t[2]
                else
                    a[llink, rlink] += t
                end
            end
        end
    end 

    return (a, b, c, d)
end