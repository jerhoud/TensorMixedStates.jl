export tdvp, dmrg, approx_W, steady_state

"""
    tdvp(evolver, t, ::State; options...)
    tdvp(evolver, t, ::Simulation; options...)

do time evolution with tdvp algorithm on a state / sim for the given time t. Also see `TdvpObserver`

# Options

- `nsweeps`: number sweeps to do (time step = t / nsweeps) 
- `coefs`: coefficients for time dependent evolver
- `n_expand`: do expansion steps every n_expand steps (default 0 means no expansion)
- `n_hermitianize`: make hermitian (for mixed states) every n_hermitianize steps (default 0 for no corrections)
- others are identical to ITensorMPS.tdvp
"""
function tdvp(pre::PreMPO{R}, t::Number, state::State{R};
    observer! = NoObserver(), coefs=nothing, n_expand = 0, n_hermitianize = 0,
    nsweeps = 1, time_start = zero(t), limits::Limits=Limits(), kwargs...) where {R <: PM}
    time_dep = !isnothing(coefs)
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
        st = tdvp(mpo, dt, st; nsweeps = 1, limits.cutoff, limits.maxdim, kwargs...)
        if n_hermitianize ≠ 0 && mod(sweep, n_hermitianize) == 0
            st = hermitianize(State(state, st); limits).state
        end    
        measure!(observer!; sweep, state = st, current_time, mpo)
        if n_expand ≠ 0 && mod(sweep, n_expand) == 0
            st = expand(st, mpo; alg="global_krylov")
        end
    end
    return State(state, st)
end

tdvp(op, t::Number, state::State; kwargs...) =
    tdvp(PreMPO(state, op), t, state; kwargs...)

"""
    dmrg(hamiltonian, ::State; options...)
    dmrg(hamiltonian, ::Simulation; options...)

optimize for ground state of the given Hamiltonian starting with state / simulation using dmrg.

Note that Dmrg does not work for mixed representations.

# Options

- `nsweeps`: number of sweeps
- `observer!`: observer (see `DmrgObserver`)
- `limits`: constraints on the mps (`cutoff` and `maxdim` may be vectors with different values for each sweep)
- others identical to ITensorMPS.dmrg
"""
function dmrg(mpo::MPO, state::State; nsweeps = 1, observer! = NoObserver(), limits::Limits=Limits(), kwargs...)
    e, st = dmrg(mpo, state.state; outputlevel = 0, nsweeps, observer = observer!, limits.cutoff, limits.maxdim, kwargs...)
    return (e, State(state, st))
end

dmrg(op, state::State; kwargs...) =
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

make_approx_W(op, t::Number, state::State; order::Int, w::Int) =
    make_approx_W(PreMPO(state, op), t; order, w)

"""
    approx_W(evolver, t, ::State; options...)
    approx_W(evolver, t, ::Simulation; options...)

time evolution using approximation WI or WII at a given order. Also see `ApproxWObserver`

# Options

- `coefs`: coefficients for time dependent evolution
- `n_hermitianize`: make hermitian (for mixed states) every n_hermitianize steps (default 0 for no corrections)
- `nsweeps`: number of steps (time step is t / nsweeps)
- `order`: order of approximation
- `w`: 1 or 2 for WI or WII
- `observer!`: observer (see ApproxWObserver)
- `time_start`: the simulation time at the beginning of evolution
- `limits`: MPS constraints
"""
function approx_W(pre::PreMPO{R}, t::Number, state::State{R}; coefs = nothing, n_hermitianize::Int = 0,
    nsweeps::Int = 1, order::Int = 1, w::Int = 1, observer! = NoObserver(), time_start = zero(t),
    limits::Limits=Limits(), kwargs...) where {R <: PM}
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
            st = apply(mpo, st; limits.cutoff, limits.maxdim, kwargs...)
        end
        if n_hermitianize ≠ 0 && mod(sweep, n_hermitianize) == 0
            st = hermitianize(State(state, st); limits).state;
        end
        measure!(observer!; sweep, state = st, current_time, mpos)
    end
    return State(state, st)    
end

approx_W(op, t::Number, state::State; kwargs...) =
    approx_W(PreMPO(state, op), t, state; kwargs...)

"""
    steady_state(lindbladian, state; kwargs...)

compute the steady state of the given Lindbladian starting on the given mixed state using DMRG on (L+)L

return achieved "energy" (which should be zero) and computed steady state

# Options
- `nsweeps`: number of sweeps
- `observer!`: observer (see `DmrgObserver`)
- `limits`: constraints on the mps (`cutoff` and `maxdim` may be vectors with different values for each sweep)
- `mpo_limits`: sets the limit on the MPO of (L+)L (default is no truncation)
- `alg`: is "naive"(default) or "zipup": alorithm to compute (L+)L 
- others identical to ITensorMPS.dmrg

"""
function steady_state(op::IndexedOp{Mixed}, state::State{Mixed};
    limits::Limits = Limits(), nsweeps::Int = 1,
    observer! = NoObserver(), mpo_limits::Limits = Limits(), alg::String = "naive", kwargs...)
    truncate = (mpo_limits != Limits())
    l = make_mpo(state, op)
    l2 = apply(replaceprime(dag(l)', 2=>0), l; mpo_limits.cutoff, mpo_limits.maxdim, alg, truncate)
    return dmrg(l2, state; nsweeps, limits, observer!, kwargs...)
end
