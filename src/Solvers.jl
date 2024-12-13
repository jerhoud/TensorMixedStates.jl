export tdvp, dmrg, TdvpObserver, DmrgObserver

import ITensorMPS: tdvp, dmrg, measure!, checkdone!

struct TdvpObserver <: AbstractObserver
    sim::Simulation
    measurements
    periodicity
end

mutable struct DmrgObserver <: AbstractObserver
    sim::Simulation
    measurements
    periodicity
    tol
    energy
    DmrgObserver(s, m, p, t) = new(s, m, p, t, 0.)
end

function measure!(o::TdvpObserver; sweep, half_sweep_is_done, half_sweep, current_time, state, kwargs...)
    if half_sweep_is_done && half_sweep == 2 && mod(sweep, o.periodicity) == 0
        st = State(o.sim.state, state)
        sim = Simulation(o.sim, st, o.sim.time + current_time)
        output(sim, o.measurements; sweep)
    end
    return nothing
end

function checkdone!(o::DmrgObserver; energy, sweep, state, kwargs...)
    if sweep ≠ 1 && abs(o.energy - energy) < o.tol
        return true
    elseif mod(sweep, o.periodicity) == 0
        st = State(o.sim.state, state)
        sim = Simulation(o.sim, st)
        output(sim, o.measurements; energy, sweep)
    end
    o.energy = energy
    return false
end

tdvp(state::State, mpo::MPO, t::Number, n::Int; kwargs...) =
    State(state, tdvp(mpo, n * t, state.state; time_step = t, nsteps = n, kwargs...))

function dmrg(state::State, mpo::MPO, n::Int; kwargs...)
    e, st = dmrg(mpo, state.state; nsweeps = n, kwargs...)
    return (e, State(state, st)) 
end

function evolve(state::State, mpo::MPO, t::Number, n::Int; order, kwargs...)

end