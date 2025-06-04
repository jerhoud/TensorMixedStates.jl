export TdvpObserver, DmrgObserver, ApproxWObserver

"""
    struct TdvpObserver
    TdvpObserver(sim, measurements, period)

an observer for tdvp which make measurements every period steps
"""
struct TdvpObserver <: AbstractObserver
    sim::Simulation
    measurements
    period
end

"""
    struct ApproxWObserver
    ApproxWObserver(sim, measurements, period)

an observer for approx_W which make measurements every period steps
"""
struct ApproxWObserver <: AbstractObserver
    sim::Simulation
    measurements
    period
end

"""
    struct DmrgObserver
    DmrgObserver(sim, measurements, period, tol)

an observer for dmrg which makes and outputs measurements every period steps and stops it when energy improvements are smaller than tol
"""
mutable struct DmrgObserver <: AbstractObserver
    sim::Simulation
    measurements
    period
    tol
    energy
    DmrgObserver(sim, measurements, period, tol) = new(sim, measurements, period, tol, 0.)
end

function measure!(o::TdvpObserver; sweep, current_time, state, mpo, kwargs...)
    if mod(sweep, o.period) == 0
        st = State(o.sim.state, state)
        sim = Simulation(o.sim, st, current_time)
        output(sim, o.measurements; sweep)
    end
    if sweep == 1
        log_msg(o.sim, "Tdvp MPO: maxlinkdim=$(maxlinkdim(mpo)), memory=$(Base.summarysize(mpo))")
    end
    log_msg(o.sim, "sim_time $(round(current_time; digits=8))")
    return nothing
end

function measure!(o::ApproxWObserver; sweep, current_time, state, mpos, kwargs...)
    if mod(sweep, o.period) == 0
        st = State(o.sim.state, state)
        sim = Simulation(o.sim, st, current_time)
        output(sim, o.measurements; sweep)
    end
    if sweep == 1
        log_msg(o.sim, "Approx_W MPOS: maxlinkdim=$(maxlinkdim(mpos[1])), memory=$(Base.summarysize(mpos))")
    end
   log_msg(o.sim, "sim_time $(round(current_time; digits=8))")
    return nothing
end

function checkdone!(o::DmrgObserver; energy, sweep, psi, kwargs...)
    if o.sim.state.type isa Mixed 
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
    log_msg(o.sim, "sweep $sweep")
    return stop
end
