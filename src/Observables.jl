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
function trace(state::State; prep = nothing)
    if state.type == Pure
        norm(state.state)^2
    elseif isnothing(prep)
        scalar(prod(mixed_obs(state, i) for i in 1:length(state)))
    else
        prep.norm
    end
end

"""
    trace2(::State)

Return the trace of the square density matrix, mostly usefull for mixed representations.
Should be one for pure representation.
"""
trace2(state::State; prep = nothing) =
    if state.type == Pure
        1.
    else
        (norm(state.state) / real(trace(state; prep))) ^ 2
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
        return State(state, state.state / real(trace(state)))
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
symmetrize(state::State; cutoff = 0., maxdim = typemax(Int)) =
    if state.type == Pure
        state
    else
        State(state, 0.5*(+(state.state, dag(state).state; cutoff, maxdim)))
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
    norm::Number
end

function PreObs(::TPure, state::State)
    n = length(state)
    st = state.state
    sites = state.system.pure_sites
    vleft = Vector{ITensor}(undef, n)
    vleft[1] = st[1]
    ll = ITensorMPS.leftlim(st) + 1
    for i in 2:n
        llink = commonind(st[i-1], st[i])
        v = if i <= ll
            delta(llink, llink')
        else
            k = sites[i-1]
            vleft[i-1] * delta(k, k') * dag(st[i-1]')
        end
        vleft[i] = v * st[i]
    end
    vnorm = scalar(vleft[n] * delta(sites[n], sites[n]') * dag(st[n]'))
    vright = Vector{ITensor}(undef, n)
    vright[n] = dag(st[n]') / real(vnorm)
    rl = ITensorMPS.rightlim(st) - 1
    for i in n-1:-1:1
        rlink = commonind(st[i], st[i+1])
        v = if i >= rl
            delta(rlink, rlink') / real(vnorm)
        else
            k = sites[i+1]
            vright[i+1] * delta(k, k') * st[i+1]
        end
        vright[i] = v * dag(st[i]')
    end
    return PreObs([], vleft, vright, vnorm)
end
    
function PreObs(::TMixed, state::State)
    n = length(state)
    st = state.state
    vloc = [ mixed_obs(state, i) for i in 1:n]
    v = ITensor(1)
    vleft = [[st[1]]; [(v *= vloc[i]; v * st[i+1]) for i in 1:n-1]]
    vnorm = scalar(v * vloc[n])
    v = ITensor(1. / real(vnorm))
    vright = reverse([[v]; [v *= vloc[i] for i in n:-1:2]])
    return PreObs(vloc, vleft, vright, vnorm)
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

function expect(::TPure, state::State, p::ProdLit, pre::PreObs)
    if p.coef == 0.
        return 0.
    end
    sys = state.system
    sites = sys.pure_sites
    st = state.state
    r = ITensor(1)
    o = ITensor()
    j = 0
    foreach(p.ls) do l
        i = l.index[1]
        if j == i
            o = replaceprime(o' * l(Pure, sys), 2=>1)
        else
            if j == 0
                r = pre.left[i]
            else
                r *= o * dag(st[j]')
                for k in j+1:i-1
                    idk = sites[k]
                    r *= st[k] * delta(idk, idk')
                    r *= dag(st[k]')
                end
                r *= st[i]
            end
            j = i
            o = l(Pure, sys)
        end
    end
    r *= o * pre.right[j]
    return p.coef * scalar(r)
end

function expect(::TMixed, state::State, p::ProdLit, pre::PreObs)
    if p.coef == 0.
        return 0.
    end
    sys = state.system
    st = state.state
    r = ITensor(1)
    o = ITensor()
    j = 0
    foreach(p.ls) do l
        i = l.index[1]
        if j == i
            o = replaceprime(o' * l(Pure, sys), 2=>1)
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
            o = l(Pure, sys)
        end
    end
    r *= mixed_obs(state, o, j)
    r *= pre.right[j]
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
expect(state::State, p::ProdLit, pre=PreObs(state)) =
    expect(state.type, state, p, pre)
    
expect(state::State, op::SumLit, pre=PreObs(state)) =
    sum(op.ps; init=0) do p
        expect(state, p, pre)
    end

expect(state::State, op, pre=PreObs(state)) =
    map(op) do o
        expect(state, o, pre)
    end


function expect1_one(::TPure, state::State, o::Operator, i::Int, t::ITensor)
    idx = state.system.pure_sites[i]
    return scalar(o(idx) * t)
end

function expect1_one(::TMixed, state::State, o::Operator, i::Int, t::ITensor)
    idx = state.system.pure_sites[i]
    return scalar(mixed_obs(state, o(idx), i) * t)
end

expect1_one(tp, state::State, op, i::Int, t::ITensor) =
    map(op) do o
        expect1_one(tp, state, o, i, t)
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
function expect1(state::State, op, pre::PreObs=PreObs(state))
    n = length(state)
    r = [ expect1_one(state.type, state, op, i, pre.left[i] * pre.right[i]) for i in 1:n ]
    return unroll(r)
end


function expect2(::TPure, state::State, ops::Vector{Tuple{Operator, Operator}}, pre::PreObs)
    oplist = collect(Set([first.(ops) ; last.(ops)]))
    n = length(state)
    sys = state.system
    st = state.state
    r = Matrix(undef, n, n)
    for i in 1:n
        idx = sys.pure_sites[i]
        ldict = Dict{Operator, ITensor}()
        ti = pre.left[i] * pre.right[i]
        r[i, i] = map(ops) do (o1, o2)
            t = replaceprime(o1(idx) * o2(idx)' + o1(idx)' * o2(idx), 2 => 1)
            scalar(ti * t) / 2
        end
        for op in oplist
            t = op(idx)
            if op.fermionic
                t = replaceprime(t' * F(idx), 2 => 1)
            end
            ldict[op] = pre.left[i] * t * dag(st[i]')
        end
        rdict = Dict{Operator, ITensor}()
        for j in i+1:n
            jdx = sys.pure_sites[j]
            for op in oplist
                rdict[op] = pre.right[j] * op(jdx) * st[j]
            end
            r[i, j] = map(ops) do (o1, o2)
                scalar(ldict[o1] * rdict[o2])
            end
            r[j, i] = map(ops) do (o1, o2)
                v = scalar(ldict[o2] * rdict[o1])
                if o1.fermionic
                    v *= -1
                end
                v
            end
            if j < n
                for op in oplist
                    p = if op.fermionic
                        F(jdx)
                    else
                        delta(jdx, jdx')
                    end
                    ldict[op] *= st[j] * p
                    ldict[op] *= dag(st[j]')
                end
            end    
        end
    end
    return unroll(r)
end

function expect2(::TMixed, state::State, ops::Vector{Tuple{Operator, Operator}}, pre::PreObs)
    oplist = collect(Set([first.(ops) ; last.(ops)]))
    n = length(state)
    st = state.state
    sys = state.system
    r = Matrix(undef, n, n)
    for i in 1:n
        idx = sys.pure_sites[i]
        ldict = Dict{Operator, ITensor}()
        ti = pre.left[i] * pre.right[i]
        r[i, i] = map(ops) do (o1, o2)
            t = replaceprime(o1(idx) * o2(idx)' + o1(idx)' * o2(idx), 2 => 1)
            scalar(ti * mixed_obs(state, t, i)) / 2
        end
        for op in oplist
            t = op(idx)
            if op.fermionic
                t = replaceprime(t' * F(idx), 2 => 1)
            end
            ldict[op] = pre.left[i] * mixed_obs(state, t, i)
        end
        rdict = Dict{Operator, ITensor}()
        for j in i+1:n
            jdx = sys.pure_sites[j]
            q = st[j] * pre.right[j]
            for op in oplist
                rdict[op] = q * mixed_obs(state, op(jdx), j)
            end
            r[i, j] = map(ops) do (o1, o2)
                scalar(ldict[o1] * rdict[o2])
            end
            r[j, i] = map(ops) do (o1, o2)
                v = scalar(ldict[o2] * rdict[o1])
                if o1.fermionic
                    v *= -1
                end
                v
            end
            if j < n
                for op in oplist
                    p = if op.fermionic
                        st[j] * mixed_obs(state, F(jdx), j)
                    else
                        pre.loc[j]
                    end
                    ldict[op] *= p
                end    
            end
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
expect2(state::State, ops::Tuple{Operator, Operator}, pre=PreObs(state)) =
    expect2(state, [ops], pre)

expect2(state::State, ops::Vector{Tuple{Operator, Operator}}, pre=PreObs(state)) =
    expect2(state.type, state, ops, pre)



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
    sump = sum(sp)
    ee = -sum(p * log(p) for p in sp) / sump + log(sump)
    return (ee, sp)
end
