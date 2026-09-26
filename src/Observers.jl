export TdvpObserver, DmrgObserver, ApproxWObserver

"""
    struct TdvpObserver
    TdvpObserver(sim, measurements, period)

an observer for tdvp which make measurements every period steps
"""
struct TdvpObserver <: AbstractObserver
    sim::Simulation
    measurements::Union{Vector, Pair}
    period::Int
end

"""
    struct ApproxWObserver
    ApproxWObserver(sim, measurements, period)

an observer for approx_W which make measurements every period steps
"""
struct ApproxWObserver <: AbstractObserver
    sim::Simulation
    measurements::Union{Vector, Pair}
    period::Int
end

"""
    struct DmrgObserver
    DmrgObserver(sim, measurements, period, tol)

an observer for dmrg which makes and outputs measurements every period steps and stops it
when energy improvements are smaller than tol. `done` is the number of sweeps already done
before this run, non zero when the phase resumes from a checkpoint.
"""
mutable struct DmrgObserver <: AbstractObserver
    sim::Simulation
    measurements::Union{Vector, Pair}
    period::Int
    tol::Number
    energy::Float64
    done::Int
    DmrgObserver(sim, measurements, period, tol, done = 0) = new(sim, measurements, period, tol, 0., done)
end

function measure!(o::TdvpObserver; sweep, current_time, state, mpo, kwargs...)
    if sweep_due(o.period, sweep)
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
    if sweep_due(o.period, sweep)
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
    c = o.sim.checkpoint
    # ITensorMPS counts the sweeps of a resumed run from 1 again, and what is measured, logged
    # and checkpointed goes by those of the phase, as it does for tdvp and approx_W. The
    # energy is compared with the sweep before in this run, the only one known
    s = sweep + o.done
    stop = sweep ≠ 1 && abs(o.energy - energy) < o.tol
    # the measurements of a sweep are written before the checkpoint records how far the
    # output files go, so that resuming at the next sweep does not cut them away. This is
    # the order tdvp and approx_W have by construction, their observer being asked to
    # measure before it is asked whether to stop.
    if stop || stop_requested(c) || sweep_due(o.period, s)
        st = normalize(State(o.sim.state, psi))
        sim = Simulation(o.sim, st)
        output(sim, o.measurements; energy, sweep = s)
    end
    # a dmrg sweep does not change the simulation time, so the state is checkpointed as it
    # is and the sweep count is what a resume needs
    if checkpoint_step!(c, o.sim, State(o.sim.state, psi), o.sim.time, s)
        stop = true
    end
    o.energy = energy
    log_msg(o.sim, "sweep $s")
    return stop
end

# tdvp and approx_W ask their observer whether to stop, the same way dmrg does. This is
# where checkpointing lives: the observer holds the simulation, the solvers do not have to
# know anything about it.
checkdone!(o::Union{TdvpObserver, ApproxWObserver}; sweep, state, current_time, kwargs...) =
    checkpoint_step!(o.sim.checkpoint, o.sim, State(o.sim.state, state), current_time, sweep)
