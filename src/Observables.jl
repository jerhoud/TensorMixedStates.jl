import ITensors: norm
import ITensorMPS: expect
export random_state, trace, trace2, norm
export Preprocess, preprocess, expect, expect1, expect2

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


function random_state(::TPure, system::System, linkdims::Int; start_time::Float64=0.0)
    st = random_mps(ComplexF64, system.pure_sites; linkdims)
    return State(Pure, system, st, start_time)
end

function random_state(::TPure, state::State, linkdims::Int; start_time::Float64=state.time)
    if state.type ≠ Pure
        error("cannot produce a random state from a mixed state")
    end
    st = random_mps(ComplexF64, state.system.pure_sites, state.state; linkdims)
    return State(Pure, system, st, start_time)
end

function random_state(::TMixed, system::System, linkdims::Int; start_time::Float64=0.0)
    n = length(system)
    psites = system.pure_sites
    msites = system.mixed_sites
    super = System(vcat(psites, sim.(psites)), vcat(msites, sim.(msites)))
    super_pure = random_state(Pure, super, 1 + linkdims ÷ 2)
    super_mixed = truncate(mix_state(super_pure), maxdim = linkdims, cutoff = 0)
    t = ITensor(1)
    for i in 2n:-1:n+1
        t *= mixed_obs(super_mixed, i)
    end
    super_mixed.state[n] *= t
    return State(Mixed, system, MPS(super_mixed.state[1:n]), start_time)
end

function trace(state::State)
    if state.type == Pure
        return 1.0
    else
        scalar(prod(mixed_obs(state, i) for i in 1:length(state)))
    end
end

trace2(state::State) =
    if state.type == Pure
        return 1.0
    else
        return norm(state.state)^2
    end

norm(state::State) =
    norm(state.state)

struct Preprocess
    loc::Vector{ITensor}
    left::Vector{ITensor}
    right::Vector{ITensor}
end


preprocess(::TPure, ::State) = nothing
    
function preprocess(::TMixed, state::State)
    n = length(state)
    st = state.state
    vloc = [ mixed_obs(state, i) for i in 1:n]
    v = ITensor(1)
    vleft = vcat([st[1]], [(v *= vloc[i]; v * st[i+1]) for i in 1:n-1])
    v = ITensor(1)
    vright = reverse(vcat([ITensor(1)], [v *= vloc[i] for i in n:-1:2]))
    return Preprocess(vloc, vleft, vright)
end

preprocess(state::State) = preprocess(state.type, state)

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

function expect(::TPure, state::State, p::ProdLit, ::Nothing)
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
                    r *= st[k]
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

function expect(::TMixed, state::State, p::ProdLit, prep::Preprocess)
    st = state.state
    system = state.system
    psites = system.pure_sites
    r = ITensor(1)
    o = ITensor()
    j = 0
    foreach(p.ls) do l
        i = l.index[1]
        if j == i
            o = replaceprime(o' * op(l.opname, psites[i]; l.param...), 2=>1)
        else
            if j == 0
                r = prep.left[i]
            else
                r *= mixed_obs(state, o, j)
                for k in j+1:i-1
                    r *= prep.loc[k]
                end
                r *= st[i]
            end
            j = i
            o = op(l.opname, psites[i]; l.param...)
        end
    end
    if j ≠ 0
        r *= mixed_obs(state, o, j)
        r *= prep.right[j]
    end
    return p.coef * scalar(r)
end

expect(state::State, p::ProdLit, prep=preprocess(state)) =
    expect(state.type, state, p, prep)
    
expect(state::State, op::SumLit, prep=preprocess(state)) =
    sum(op.ps; init=0) do p
        expect(state, p, prep)
    end

expect(state::State, op, prep=preprocess(state)) =
    map(op) do o
        expect(state, o, prep)
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

function expect1(::TPure, state::State, op, ::Nothing)
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
    
function expect1(::TMixed, state::State, op, prep::Preprocess)
    n = length(state)
    r = [ expect1_one(state, op, i, prep.left[i] * prep.right[i]) for i in 1:n ]
    return unroll(r)
end

expect1(state::State, op, prep=preprocess(state)) =
    expect1(state.type, state, op, prep)



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
        expect2_one(tp, op, i1, i2, t, rev)
    end

function expect2(::TPure, state::State, ops, prep::Nothing)
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
            if i < n
                ridx = commonind(st[i], st[i+1])
                t *= delta(ridx, ridx')
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
    
function expect2(::TMixed, state::State, ops, prep::Preprocess)
    n = length(state)
    st = state.state
    r = Matrix(undef, n, n)
    for i in 1:n
        l = prep.left[i]
        r[i, i] = expect2_one(state, ops, i, i, l * prep.right[i], false)
        for j in i+1:n
            t = l * (st[j] * prep.right[j])
            r[i, j] = expect2_one(state, ops, i, j, t, false)
            r[j, i] = expect2_one(state, ops, i, j, t, true)
            l *= prep.loc[j]
        end
    end
    return unroll(r)
end

expect2(state::State, ops, prep=preprocess(state)) =
    expect2(state.type, state, ops, prep)