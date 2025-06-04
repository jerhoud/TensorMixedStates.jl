export trace, trace2, norm, normalize, dag, symmetrize, normsym, hermitianity
export expect, expect1, expect2
export entanglement_entropy, partial_trace

function tensor_trace(state::State, i::Int)
    s = state.system
    j = s[Pure(), i]
    k = s[Mixed(), i]
    return dense(delta(j', j)) * combinerto(k, j', j)
end

function tensor_obs(state::State, ind::Indexed{Pure, 1})
    s = state.system
    t = tensor(s, ind)
    j = s[Pure(), ind.index...]
    k = s[Mixed(), ind.index...]
    return t * combinerto(k, j, j')
end

tensor_obs(state::State, i::Int) =
    tensor_trace(state, i) * state.state[i]

function tensor_dag(state::State, i::Int)
    s = state.system
    j = s[Pure(), i]
    k = s[Mixed(), i]
    return dag(state.state[i]) * combinerto(k, j', j) * combinerto(k, j, j')
end

function get_loc(state::State, i::Int)
    l = state.preobs.loc
    if isempty(l)
        create_loc!(l, state.type, state)
    end
    return l[i]
end

function get_right(state::State, i::Int)
    r = state.preobs.right
    if isempty(r)
        create_right!(r, state.type, state)
    end
    return r[i]
end

function get_left(state::State, i::Int)
    l = state.preobs.left
    if length(l) < i
        create_left!(l, state.type, state, i)
    end
    return l[i]
end

"""
    trace(::State)

Return the trace of the system, mostly usefull for mixed representations.
This should be one.
"""
function trace(state::State)
    t = state.preobs.trace
    if isempty(t)
        create_trace!(t, state.type, state)
    end
    return t[1]
end

create_loc!(_, ::Pure, ::State) = error("bug: get_loc on pure states")
function create_loc!(l, ::Mixed, state::State)
    n = length(state)
    resize!(l, n)
    for i in 1:n
        l[i] = tensor_obs(state, i)
    end
    return l
end

function create_right!(r, ::Pure, state::State)
    st = state.state
    s = state.system
    n = length(state)
    rl = ITensorMPS.rightlim(st) - 1
    resize!(r, n)
    r[n] = dag(st[n]')
    for i in n-1:-1:1
        rlink = commonind(st[i], st[i+1])
        v = if i >= rl
            delta(rlink, rlink')
        else
            k = s[Pure(), i+1]
            r[i+1] * delta(k, k') * st[i+1]
        end
        r[i] = v * dag(st[i]')
    end
    return r
end

function create_right!(r, ::Mixed, state::State)
    n = length(state)
    resize!(r, n)
    t = r[n] = ITensor(1.)
    for i in n-1:-1:1
        t = r[i] = get_loc(state, i+1) * t
    end
    return r
end

function create_trace!(t, ::Pure, state::State)
    resize!(t, 1)
    k = state.system[Pure(), 1]
    t[1] = scalar(get_right(state, 1) * delta(k, k') * state.state[1])
    return t
end

function create_trace!(t, ::Mixed, state::State)
    resize!(t, 1)
    t[1] = scalar(get_loc(state, 1) * get_right(state, 1))
    return t
end

function create_left!(l, ::Pure, state::State, i::Int)
    st = state.state
    s = state.system
    ll = ITensorMPS.leftlim(st) + 1
    j = length(l)
    resize!(l, i)
    if j == 0
        l[1] = st[1] / real(trace(state))
        j = 1
    end
    for k in j+1:i
        llink = commonind(st[k-1], st[k])
        v = if k <= ll
            delta(llink, llink')
        else
            idx = s[Pure(), k-1]
            l[k-1] * delta(idx, idx') * dag(st[k-1]')
        end
        l[k] = v * st[k]
    end
    return l
end

function create_left!(l, ::Mixed, state::State, i::Int)
    st = state.state
    j = length(l)
    resize!(l, i)
    if j == 0
        l[1] = st[1] / real(trace(state))
        j = 1
    end
    for k in j+1:i
        l[k] = l[k-1] * tensor_trace(state, k-1) * st[k]
    end
    return l
end

"""
    trace_error(::State)

Return the deviaton to 1 of the trace of the system, mostly usefull for mixed representations.
This should be zero.
"""
trace_error(state::State; preobs = nothing) = 1 - trace(state; preobs) 

"""
    trace2(::State)

Return the trace of the square density matrix, mostly usefull for mixed representations.
This is one for pure representations.
"""
trace2(state::State) =
    if state.type isa Pure
        1.
    else
        (norm(state.state) / real(trace(state))) ^ 2
    end

"""
    norm(::State)

Return the norm of the state, mostly usefull for pure representations.
This should be one for pure representation.
"""
norm(state::State) =
    norm(state.state)


"""
    normalize(::State)

normalize state so that norm = 1 for pure state and trace = 1 for mixed state
"""
normalize(state::State) =
    if state.type isa Pure
        return State(state, normalize(state.state))
    else
        return State(state, state.state / real(trace(state)))
    end


"""
    dag(::State)

adjoint of density matrix for mixed representation
"""
dag(state::State) =
    if state.type isa Pure
        state
    else
        n = length(state)
        st = MPS(n)
        for i in 1:n
            st[i] = tensor_dag(state, i)
        end
        return State(state, st)
    end


"""
    symmetrize(::State)

symmetrize state so that it is hermitian (only for mixed state)
"""
symmetrize(state::State; limits::Limits=Limits()) =
    if state.type isa Pure
        state
    else
        State(state, 0.5*(+(state.state, dag(state).state; limits.cutoff, limits.maxdim)))
    end



"""
    normsym(::State)

equivalent to normalize(symmetrize(state))
"""
normsym(state::State) = normalize(symmetrize(state))

"""
    hermitianity(::State)

hermitianity measure whether density matrix for mixed state is Hermitian as it should.

return a value from 0 (anti Hermitian) to 1 (Hermitian)

return 1 for pure state
"""
hermitianity(state::State) =
    if state.type isa Pure
        1.
    else
        0.5 + 0.5 * real(dot(state.state, dag(state).state)) / norm(state.state)^2
    end


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

function expect(::Pure, state::State, coef::Number, subs::Vector{<:ExprIndexed{Pure}})
    if coef == 0.
        return 0.
    end
    sys = state.system
    st = state.state
    local r::ITensor
    j = 0
    foreach(subs) do ind
        i = ind.index[1]
        if j == 0
            r = get_left(state, i)
        else
            r *= dag(st[j]')
            for k in j+1:i-1
                idk = sys[Pure(), k]
                r *= st[k] * delta(idk, idk')
                r *= dag(st[k]')
            end
            r *= st[i]
        end
        r *= tensor(sys, ind)
        j = i
    end
    r *= get_right(state, j)
    return coef * scalar(r)
end

function expect(::Mixed, state::State, coef::Number, subs::Vector{<:ExprIndexed{Pure}})
    if coef == 0.
        return 0.
    end
    st = state.state
    local r::ITensor
    j = 0
    foreach(subs) do ind
        i = ind.index[1]
        if j == 0
            r = get_left(state, i)
        else
            for k in j+1:i-1
                r *= get_loc(state, k)
            end
            r *= st[i]
        end
        r *= tensor_obs(state, ind)
        j = i
    end
    r *= get_right(state, j)
    return coef * scalar(r)
end

"""
    expect(::State, obs)

Compute expectation values of `obs` on the given state.

# Examples
    expect(state, X(1)*Y(2) + Y(1)*Z(3))
    expect(state, [X(1)*Y(2), X(3), Z(1)*X(2)])

"""
expect(state::State, p::ExprIndexed{Pure}) =
    expect(state.type, state, prodcoef(p), prodsubs(p))
    
expect(state::State, op::SumOp{Pure, IndexOp}) =
    sum(op.ps; init=0) do p
        expect(state, p)
    end

expect(state::State, op) =
    map(op) do o
        expect(state, o)
    end


expect1_one(::Pure, state::State, op::SimpleOp, i::Int, t::ITensor) =
    scalar(t * tensor(state.system, op(i)))

expect1_one(::Mixed, state::State, op::SimpleOp, i::Int, t::ITensor) =
    scalar(t * tensor_obs(state, op(i)))

expect1_one(tp, state::State, ops, i::Int, t::ITensor) =
    map(ops) do o
        expect1_one(tp, state, o, i, t)
    end


"""
    expect1(::State, op)

Compute the expectation values of the given operators on all sites.

# Examples
    expect1(state, X)
    expect1(state, [X, Y, Z])

"""
function expect1(state::State, op)
    n = length(state)
    r = [ expect1_one(state.type, state, op, i, get_left(state, i) * get_right(state, i)) for i in 1:n ]
    return unroll(r)
end


function expect2(::Pure, state::State, ops::Vector{<:Tuple{SimpleOp, SimpleOp}})
    oplist = collect(Set([first.(ops) ; last.(ops)]))
    n = length(state)
    sys = state.system
    st = state.state
    r = Matrix(undef, n, n)
    for i in 1:n
        ldict = Dict{SimpleOp, ITensor}()
        ti = get_left(state, i) * get_right(state, i)
        r[i, i] = map(ops) do (o1, o2)
            scalar(ti * tensor(sys, (o1 * o2)(i)))
        end
        for op in oplist
            opF = fermionic(op) ? op * F : op
            ldict[op] = get_left(state, i) * tensor(sys, opF(i)) * dag(st[i]')
        end
        rdict = Dict{SimpleOp, ITensor}()
        for j in i+1:n
            for op in oplist
                rdict[op] = get_right(state, j) * tensor(sys, op(j)) * st[j]
            end
            r[i, j] = map(ops) do (o1, o2)
                scalar(ldict[o1] * rdict[o2])
            end
            r[j, i] = map(ops) do (o1, o2)
                v = scalar(ldict[o2] * rdict[o1])
                if fermionic(o1)
                    v *= -1
                end
                v
            end
            if j < n
                for op in oplist
                    p = if fermionic(op)
                        tensor(sys, F(j))
                    else
                        jdx = sys[Pure(), j]
                        delta(jdx', jdx)
                    end
                    ldict[op] *= st[j] * p
                    ldict[op] *= dag(st[j]')
                end
            end    
        end
    end
    return unroll(r)
end

function expect2(::Mixed, state::State, ops::Vector{<:Tuple{SimpleOp, SimpleOp}})
    oplist = collect(Set([first.(ops) ; last.(ops)]))
    n = length(state)
    st = state.state
    r = Matrix(undef, n, n)
    for i in 1:n
        ldict = Dict{SimpleOp, ITensor}()
        ti = get_left(state, i) * get_right(state, i)
        r[i, i] = map(ops) do (o1, o2)
            scalar(ti * tensor_obs(state, (o1 * o2)(i)))
        end
        for op in oplist
            opF = fermionic(op) ? op * F : op
            ldict[op] = get_left(state, i) * tensor_obs(state, opF(i))
        end
        rdict = Dict{SimpleOp, ITensor}()
        for j in i+1:n
            q = st[j] * get_right(state, j)
            for op in oplist
                rdict[op] = q * tensor_obs(state, op(j))
            end
            r[i, j] = map(ops) do (o1, o2)
                scalar(ldict[o1] * rdict[o2])
            end
            r[j, i] = map(ops) do (o1, o2)
                v = scalar(ldict[o2] * rdict[o1])
                if fermionic(o1)
                    v *= -1
                end
                v
            end
            if j < n
                for op in oplist
                    p = if fermionic(op)
                        st[j] * tensor_obs(state, F(j))
                    else
                        get_loc(state, j)
                    end
                    ldict[op] *= p
                end    
            end
        end
    end
    return unroll(r)
end

"""
    expect2(::State, op_pairs)

Compute the expectation values of the given pairs of operators on all sites.

# Examples
    expect2(state, (X, X))
    expect2(state, [(X, Y), (X, Z), (Y, Z)])
    ...
"""
expect2(state::State, ops::Tuple{SimpleOp, SimpleOp}) =
    expect2(state, [ops])[1]

expect2(state::State, ops::Vector{<:Tuple{SimpleOp, SimpleOp}}) =
    expect2(state.type, state, ops)



"""
    entanglement_entropy(::State, ::Int)

Return the entanglement entropy of the given state at the given site.
Also return the associated spectrum.
"""
function entanglement_entropy(state::State, pos::Int)
    s = orthogonalize(state.state, pos)
    _, S = svd(s[pos], (linkinds(s, pos-1)..., siteinds(s, pos)...))
    sp = [ S[i,i]^2 for i in 1:dim(S, 1) ]
    sp /= sum(sp)
    ee = -sum(p * log(p) for p in sp)
    return (ee, sp)
end

"""
    partial_trace(::State, ::Vector{Int} [; keepers = true])

return the state partially traced at the given positions,
alternatively one can give the positions to keep by setting `keepers = true`
"""
function partial_trace(state::State, pos::Vector{Int}; keepers::Bool = false)
    if state.type isa Pure
        error("cannot partial trace a pure representation")
    end
    n = length(state)
    if keepers
        keep = pos
    else
        keep = filter(e->e ∉ pos, 1:n)
    end
    kn = length(keep)
    if kn == 0
        error("partial_trace cannot trace all sites of a state")
    end
    mps = state.state
    sys = state.system
    j = 0
    t = Vector{ITensor}(undef, kn)
    s = Vector{AbstractSite}(undef, kn)
    sp = Vector{Index}(undef, kn)
    sm = Vector{Index}(undef, kn)
    for (i, k) in enumerate(keep)
        s[i] = sys[k] 
        sp[i] = sys[Pure(), k]
        sm[i] = sys[Mixed(), k]
        if i == 1
            t[1] = get_left(state, k)
        else
            for l in j+1:k-1
                t[i-1] *= get_loc(state, l)
            end 
            t[i] = mps[k]
        end
        j = k
    end
    if j == 0
        error("partial_trace cannot trace all sites of a state")
    else
        t[kn] *= get_right(state, keep[kn])
    end
    for i in 1:kn-1
        idx = commonind(t[i], t[i+1])
        jdx = settags(idx, "Link, l=$i")
        replaceind!(t[i], idx, jdx)
        replaceind!(t[i+1], idx, jdx)
    end
    return State(Mixed(), System(s, sp, sm), MPS(t))
end