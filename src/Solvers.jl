export tdvp, dmrg
import ITensorMPS: tdvp, dmrg


tdvp(mpo::MPO, t::Number, state::State; observer = NoObserver(), observer! = observer, kwargs...) =
    State(state, tdvp(mpo, t, state.state; observer!, kwargs...))

tdvp(op::Lits, t::Number, state::State; kwargs...) =
    tdvp(make_mpo(state, op), t, state; kwargs...)

function tdvp(op::Vector{<:Lits}, t::Number, state::State; coefs::Vector{<:Function}, nsweeps = 1, kwargs...)
    pre = PreMPO(state, op)
    st = state.state
    dt = t / nsweeps
    for sweep in 1:nsweeps
        tsweep = t + (sweep - 1) * dt
        tf = tsweep + dt/2
        cs = map(f->f(tf), coefs)
        st = tdvp(make_mpo(pre, cs), dt, st; nsweeps = 1, time_start = tsweep, kwargs...)
    end
    return State(state, st)
end


function dmrg(mpo::MPO, state::State; kwargs...)
    e, st = dmrg(mpo, state.state; kwargs...)
    return (e, State(state, st)) 
end

dmrg(op::Lits, state::State; kwargs...) =
    dmrg(make_mpo(state, op), state; kwargs...)

const w_approx_coefs = Vector{ComplexF64}[
    [
        1.
    ],
    [
        0.5 + 0.5im,
        0.5 - 0.5im
    ],
    [
        0.10566243270259355 - 0.39433756729740643im,
        0.39433756729740643 + 0.10566243270259355im,
        0.39433756729740643 - 0.10566243270259355im,
        0.10566243270259355 + 0.39433756729740643im
    ],
    [
        0.2588533986109182 + 0.0447561340111419im,
        -0.03154685814880379 + 0.24911905427556322im,
        0.1908290521106672 - 0.23185374923210605im,
        0.16372881485443674,
        0.1908290521106672 + 0.23185374923210605im,        
        -0.03154685814880379 - 0.24911905427556322im,
        0.2588533986109182 - 0.0447561340111419im,
    ]
]

function make_approx_W(pre::PreMPO, t::Number; order::Int, w::Int, coefs = [1.])
    if order < 1 || order > length(w_approx_coefs)
        error("W approximation of order $order is not implemented")
    end
    if w == 1
        return map(c->make_approx_W1(pre, t * c, coefs), w_approx_coefs[order])
    elseif w == 2
        return map(c->make_approx_W2(pre, t * c, coefs), w_approx_coefs[order])
    else 
        error("W approximation is only defined for w=1 or 2 (not $w)")
    end
end

make_approx_W(op::Lits, t::Number, state::State; order::Int, w::Int) =
    make_approx_W(PreMPO(state, op), t; order, w)

function evolve(op::Lits, t::Number, state::State;
    nsweeps::Int = 1, order::Int = 1, w::Int = 1, observer = NoObserver(), time_start = zero(t),  kwargs...)
    dt = t / nsweeps
    mpos = make_approx_W(op, dt, state; order, w)
    for sweep in 1:nsweeps
        current_time = time_start + sweep * dt
        for mpo in mpos
            state = apply(mpo, state; kwargs...)
        end
        measure!(obsevrer; sweep, current_time, state)
    end
    return state
end

function evolve(op::Vector{Lits}, t::Number, state::State; coefs::Vector{<:Function},
    nsweeps::Int = 1, order::Int = 1, w::Int = 1, observer = NoObserver(), time_start = zero(t),  kwargs...)
    dt = t / nsweeps
    pre = PreMPO(state, op)
    for sweep in 1:nsweeps
        current_time = time_start + sweep * dt
        tf = current_time - dt / 2
        mpos = make_approx_W(pre, dt; order, w, coefs = map(f->f(tf), coefs))
        for mpo in mpos
            state = apply(mpo, state; kwargs...)
        end
        measure!(obsevrer; sweep, current_time, state)
    end
    return state
end
