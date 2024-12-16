export TdvpObserver, DmrgObserver, EvolveObserver
import ITensorMPS: measure!, checkdone!

struct TdvpObserver <: AbstractObserver
    sim::Simulation
    measurements
    periodicity
end

struct ApproxWObserver <: AbstractObserver
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
    DmrgObserver(sim, measurements, periodicity, tol) = new(sim, measurements, periodicity, tol, 0.)
end

function measure!(o::TdvpObserver; sweep, half_sweep_is_done, half_sweep, current_time, state, kwargs...)
    if half_sweep_is_done && half_sweep == 2 && mod(sweep, o.periodicity) == 0
        st = State(o.sim.state, state)
        sim = Simulation(o.sim, st, current_time)
        output(sim, o.measurements; sweep)
    end
    return nothing
end

function measure!(o::ApproxWObserver; sweep, current_time, state, kwargs...)
    if mod(sweep, o.periodicity) == 0
        st = State(o.sim.state, state)
        sim = Simulation(o.sim, st, current_time)
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
