export TdvpObserver, DmrgObserver, EvolveObserver
import ITensorMPS: measure!, checkdone!

struct TdvpObserver <: AbstractObserver
    sim::Simulation
    measurements
    period
end

struct ApproxWObserver <: AbstractObserver
    sim::Simulation
    measurements
    period
end

mutable struct DmrgObserver <: AbstractObserver
    sim::Simulation
    measurements
    period
    tol
    energy
    DmrgObserver(sim, measurements, period, tol) = new(sim, measurements, period, tol, 0.)
end

function measure!(o::TdvpObserver; sweep, half_sweep_is_done, half_sweep, current_time, state, kwargs...)
    if half_sweep_is_done && half_sweep == 2 && mod(sweep, o.period) == 0
        st = State(o.sim.state, state)
        sim = Simulation(o.sim, st, current_time)
        output(sim, o.measurements; sweep)
    end
    return nothing
end

function measure!(o::ApproxWObserver; sweep, current_time, state, kwargs...)
    if mod(sweep, o.period) == 0
        st = State(o.sim.state, state)
        sim = Simulation(o.sim, st, current_time)
        output(sim, o.measurements; sweep)
    end
    return nothing
end

function checkdone!(o::DmrgObserver; energy, sweep, psi, kwargs...)
    if o.sim.state.type == Mixed
        energy /= 2.
    end
    stop = false
    if sweep ≠ 1 && abs(o.energy - energy) < o.tol
        stop = true
    end
    if stop || mod(sweep, o.period) == 0
        st = normalize(State(o.sim.state, psi))
        sim = Simulation(o.sim, st)
        output(sim, o.measurements; energy, sweep)
    end
    o.energy = energy
    return stop
end
