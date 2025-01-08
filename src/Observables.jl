import ITensors: norm, dag
import ITensorMPS: expect, normalize

export trace, trace2, norm, normalize, dag, symmetrize, normsym
export PreObs, expect, expect1, expect2
export entanglement_entropy

function mixed_obs(state::State, t::ITensor, i::Int)
    j = state.system.pure_sites[i]
    k = state.system.mixed_sites[i]
    return t * combinerto(j', j, k)
end

function mixed_obs(state::State, i::Int)
    j = state.system.pure_sites[i]
    k = state.system.mixed_sites[i]
    state.state[i] * (dense(delta(j', j)) * combinerto(j', j, k))
end

function mixed_dag(state::State, i::Int)
    j = state.system.pure_sites[i]
    k = state.system.mixed_sites[i]
    return dag(state.state[i]) * combinerto(j', j, k) * combinerto(j, j', k)
end

"""
    trace(::State)

Return the trace of the system, mostly usefull for mixed representations.
This should be one.
"""
function trace(state::State)
    if state.type == Pure
        norm(state.state)^2
    else
        scalar(prod(mixed_obs(state, i) for i in 1:length(state)))
    end
end

"""
    trace2(::State)

Return the trace of the square density matrix, mostly usefull for mixed representations.
Should be one for pure representation.
"""
trace2(state::State) =
    if state.type == Pure
        norm(state.state)^4
    else
        norm(state.state)^2
    end

"""
    norm(::State)

Return the norm of the state, mostly usefull for pure representation.
This should be one for pure representation.
"""
norm(state::State) =
    norm(state.state)


"""
    normalize(::State)

normalize state so that norm = 1 for pure state and trace = 1 for mixed state
"""
normalize(state::State) =
    if state.type == Pure
        return State(state, normalize(state.state))
    else
        return State(state, state.state * (1. / trace(state)))
    end



dag(state::State) =
    if state.type == Pure
        state
    else
        st = copy(state.state)
        for i in 1:length(st)
            st[i] = mixed_dag(state, i)
        end
        return State(state, st)
    end


"""
    symmetrize(::State)

symmetrize state so that it is hermitian (only for mixed state)
"""
symmetrize(state::State) =
    if state.type == Pure
        state
    else
        State(state, 0.5*(state.state + dag(state).state))
    end

"""
    normsym(::State)

equivalent to (normalyze(symmetrize(state)))
"""
normsym(state::State) = normalize(symmetrize(state))

"""
    struct PreObs
    PreObs(::State)

A data structure to hold preprocessing data for observable expectation computations.
"""
struct PreObs
    loc::Vector{ITensor}
    left::Vector{ITensor}
    right::Vector{ITensor}
end

PreObs(::TPure, ::State) = nothing
    
function PreObs(::TMixed, state::State)
    n = length(state)
    st = state.state
    vloc = [ mixed_obs(state, i) for i in 1:n]
    v = ITensor(1)
    vleft = [[st[1]]; [(v *= vloc[i]; v * st[i+1]) for i in 1:n-1]]
    v = ITensor(1)
    vright = reverse([[ITensor(1)]; [v *= vloc[i] for i in n:-1:2]])
    return PreObs(vloc, vleft, vright)
end

PreObs(state::State) = PreObs(state.type, state)





unroll(x) =
    if x[1] isa Number
        x
    else
        v = [ map(t->t[i], x) for i in 1:length(x[1]) ]
        if x[1] isa Matrix
            return reshape(v, size(x[1]))
        else
            return v
        end
    end

function expect!(::TPure, state::State, p::ProdLit, ::Nothing)
    if p.coef == 0
        return 0
    end
    sites = state.system.pure_sites
    st = state.state
    imin = p.ls[1].index[1]
    imax = p.ls[end].index[1]
    if isortho(st)
        c = orthocenter(st)
        if c < imin
            orthogonalize!(st, imin)
        elseif c > imax
            orthogonalize!(st, imax)
        end
    else
        orthogonalize!(st, imin)
    end
    r = ITensor(1)
    o = ITensor()
    j = 0
    foreach(p.ls) do l
        i = l.index[1]
        if j == i
            idx = sites[i]
            o = replaceprime(o' * op(l.opname, idx; l.param...), 2=>1)
        else
            if j == 0
                if i ≠ 1
                    idx = commonind(st[i-1], st[i])
                    r = delta(idx, idx')
                end
            else
                r *= st[j]
                r *= o
                r *= dag(st[j]')
                for k in j+1:i-1
                    kdx = sites[k]
                    r *= st[k] * delta(kdx, kdx')
                    r *= dag(st[k]')
                end
            end
            j = i
            o = op(l.opname, sites[i]; l.param...)
        end
    end
    r *= st[j]
    if j ≠ length(st)
        idx = commonind(st[j], st[j+1])
        r *= delta(idx, idx')
    end
    r *= o
    r *= dag(st[j]')
    return p.coef * scalar(r)
end

function expect!(::TMixed, state::State, p::ProdLit, pre::PreObs)
    psites = state.system.pure_sites
    st = state.state
    r = ITensor(1)
    o = ITensor()
    j = 0
    foreach(p.ls) do l
        i = l.index[1]
        if j == i
            o = replaceprime(o' * op(l.opname, psites[i]; l.param...), 2=>1)
        else
            if j == 0
                r = pre.left[i]
            else
                r *= mixed_obs(state, o, j)
                for k in j+1:i-1
                    r *= pre.loc[k]
                end
                r *= st[i]
            end
            j = i
            o = op(l.opname, psites[i]; l.param...)
        end
    end
    if j ≠ 0
        r *= mixed_obs(state, o, j)
        r *= pre.right[j]
    end
    return p.coef * scalar(r)
end

"""
    expect(::State, obs[, ::PreObs])

Compute expectation values of `obs` on the given state.

# Examples
    expect(state, X(1)*Y(2) + Y(1)*Z(3))
    expect(state, [X(1)*Y(2), X(3), Z(1)*X(2)])

    pre = PreObs(state)
    expect(state, X(1)*Z(3), pre)
    expect...
"""
expect(state::State, args...) = expect!(copy(state), args...)

expect!(state::State, p::ProdLit, pre=PreObs(state)) =
    expect!(state.type, state, p, pre)
    
expect!(state::State, op::SumLit, pre=PreObs(state)) =
    sum(op.ps; init=0) do p
        expect!(state, p, pre)
    end

expect!(state::State, op, pre=PreObs(state)) =
    map(op) do o
        expect!(state, o, pre)
    end





function expect1_one(state::State, o, i::Int, t::ITensor)
    idx = state.system.pure_sites[i]
    op = o(idx)
    if state.type == Pure
        return scalar(op * t)
    else
        return scalar(mixed_obs(state, op, i) * t)
    end
end
    
expect1_one(state::State, op::Vector, i::Int, t::ITensor) =
    map(op) do o
        expect1_one(state, o, i, t)
    end

function expect1!(::TPure, state::State, op, ::Nothing)
    n = length(state)
    st = state.state
    r = Vector(undef, n)
    for i in 1:n
        orthogonalize!(st, i)
        t = st[i]
        if i > 1
            idx = commonind(st[i-1], st[i])
            t *= delta(idx, idx')
        end
        if i < n
            idx = commonind(st[i], st[i+1])
            t *= delta(idx, idx')
        end
        t *= dag(st[i]')
        r[i] = expect1_one(state, op, i, t)
    end
    return unroll(r)
end
    
function expect1!(::TMixed, state::State, op, pre::PreObs)
    n = length(state)
    r = [ expect1_one(state, op, i, pre.left[i] * pre.right[i]) for i in 1:n ]
    return unroll(r)
end

"""
    expect1(::State, op[, ::PreObs])

Compute the expectation values of the given operators on all sites.

# Examples
    expect1(state, X)
    expect1(state, [X, Y, Z])

    pre = PreObs(state)
    expect1(state, X, pre)
    ...
"""
expect1(state::State, args...) = expect1!(copy(state), args...)

expect1!(state::State, op, pre=PreObs(state)) =
    expect1!(state.type, state, op, pre)





function expect2_one(state::State, ops::Tuple{Any, Any}, i1::Int, i2::Int, t::ITensor, rev::Bool)
    if rev
        op2, op1 = ops
    else
        op1, op2 = ops
    end
    o1 = op1(state.system.pure_sites[i1])
    o2 = op2(state.system.pure_sites[i2])
    if i1 == i2
        o = replaceprime(o1' * o2, 2=>1)
        if state.type == Pure
            return scalar(o * t)
        else
            return scalar(mixed_obs(state, o, i1) * t)
        end
    else
        if state.type == Pure
             return scalar(o1 * t * o2)
        else
            return scalar(mixed_obs(state, o1, i1) * t * mixed_obs(state, o2, i2))
        end
    end
end

expect2_one(state::State, ops, i1::Int, i2::Int, t::ITensor, rev::Bool) =
    map(ops) do op
        expect2_one(state, op, i1, i2, t, rev)
    end

function expect2!(::TPure, state::State, ops, ::Nothing)
    n = length(state)
    st = state.state
    sites = state.system.pure_sites
    c = Matrix(undef, n, n)
    for i in 1:n
        orthogonalize!(st, i)
        l = st[i]
        if i > 1
            lidx = commonind(st[i-1], st[i])
            l *= delta(lidx, lidx')
        end
        l *= dag(st[i]')
        t = l
        if i < n
            ridx = commonind(st[i], st[i+1])
            t *= delta(ridx, ridx')
        end
        c[i, i] = expect2_one(state, ops, i, i, t, false)
        for j in i+1:n
            l *= st[j]
            t = l
            if j < n
                rjdx = commonind(st[j], st[j+1])
                t *= delta(rjdx, rjdx')
            end
            t *= dag(st[j]')
            c[i, j] = expect2_one(state, ops, i, j, t, false)
            c[j, i] = expect2_one(state, ops, i, j, t, true)
            jdx = sites[j]
            l *= delta(jdx, jdx') * dag(st[j]')
        end
    end
    return unroll(c)
end
    
function expect2!(::TMixed, state::State, ops, pre::PreObs)
    n = length(state)
    st = state.state
    r = Matrix(undef, n, n)
    for i in 1:n
        l = pre.left[i]
        r[i, i] = expect2_one(state, ops, i, i, l * pre.right[i], false)
        for j in i+1:n
            t = l * (st[j] * pre.right[j])
            r[i, j] = expect2_one(state, ops, i, j, t, false)
            r[j, i] = expect2_one(state, ops, i, j, t, true)
            l *= pre.loc[j]
        end
    end
    return unroll(r)
end

"""
    expect2(::State, op_pairs[, ::PreObs])

Compute the expectation values of the given pairs of operators on all sites.

# Examples
    expect2(state, (X, X))
    expect2(state, [(X, Y), (X, Z), (Y, Z)])

    pre = PreObs(state)
    expect2(state, (X, Y), pre)
    ...
"""
expect2(state::State, args...) = expect2!(copy(state), args...)

expect2!(state::State, ops, pre=PreObs(state)) =
    expect2!(state.type, state, ops, pre)



"""
    entanglement_entropy(::State, ::Int)

Return the entanglement entropy of the given state at the given site.
Also return the associated spectrum.
"""
entanglement_entropy(state::State, args...) = entanglement_entropy!(copy(state), args...)

function entanglement_entropy!(state::State, pos::Int)
    s = orthogonalize!(state.state, pos)
    _, S = svd(s[pos], (linkinds(s, pos-1)..., siteinds(s, pos)...))
    sp = [ S[i,i]^2 for i in 1:dim(S, 1) ]
    ee = -sum(p * log(p) for p in sp)
    return (ee, sp)
end
