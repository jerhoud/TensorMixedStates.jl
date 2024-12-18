export tdvp, dmrg, approx_W
import ITensorMPS: tdvp, dmrg

function tdvp(pre::PreMPO, t::Number, state::State;
    observer! = NoObserver(), coefs=nothing, n_expand = 0, nsweeps = 1, time_start = zero(t), kwargs...)
    time_dep = !isnothing(coefs)
    if n_expand == 0 && !time_dep
        return State(state, tdvp(make_mpo(pre), t, state.state; observer!, nsweeps, time_start, kwargs...))
    else
        st = state.state
        dt = t / nsweeps
        if !time_dep
            mpo = make_mpo(pre)
        end
        for sweep in 1:nsweeps
            current_time = time_start + sweep * dt
            if time_dep
                tf = current_time - dt / 2
                mpo = make_mpo(pre, map(f->f(tf), coefs))
            end
            st = tdvp(mpo, dt, st; nsweeps = 1, kwargs...)
            if n_expand ≠ 0 && mod(sweep, n_expand) == 0
                st = expand(st, mpo; alg="global_krylov")
            end
            measure!(observer!; sweep, st, current_time)
        end
        return State(state, st)    
    end
end

tdvp(op, t::Number, state::State; kwargs...) =
    tdvp(PreMPO(state, op), t, state; kwargs...)


function dmrg(mpo::MPO, state::State; nsweeps = 1, observer! = NoObserver(), kwargs...)
    e, st = dmrg(mpo, state.state; outputlevel = 0, nsweeps, observer = observer!, kwargs...)
    return (e, State(state, st))
end

dmrg(op, state::State; kwargs...) =
    dmrg(make_mpo(state, op, MixObservable), state; kwargs...)

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

make_approx_W(op, t::Number, state::State; order::Int, w::Int) =
    make_approx_W(PreMPO(state, op), t; order, w)

correct_approx_w(state::State) = state

function approx_W(pre::PreMPO, t::Number, state::State; coefs = nothing, n_correct::Int = 0,
    nsweeps::Int = 1, order::Int = 1, w::Int = 1, observer! = NoObserver(), time_start = zero(t),  kwargs...)
    st = state.state
    dt = t / nsweeps
    time_dep = !isnothing(coefs)
    if !time_dep
        mpos = make_approx_W(pre, dt; order, w)
    end
    for sweep in 1:nsweeps
        current_time = time_start + sweep * dt
        if time_dep
            tf = current_time - dt / 2
            mpos = make_approx_W(pre, dt; order, w, coefs = map(f->f(tf), coefs))
        end
        for mpo in mpos
            st = apply(mpo, st; kwargs...)
        end
        if n_correct ≠ 0 && mod(sweep, n_correct) == 0
            st = correct_approx_w(st)
        end
        measure!(observer!; sweep, st, current_time)
    end
    return State(state, st)    
end

approx_W(op, t::Number, state::State; kwargs...) =
    approx_W(PreMPO(state, op), t, state; kwargs...)
