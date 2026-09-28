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
    DmrgObserver(sim, measurements, period, tol[, done; nsweeps, energy])

an observer for dmrg which makes and outputs measurements every period steps and stops it
when energy improvements are smaller than tol. `done` is the number of sweeps already done
before this run, non zero when the phase resumes from a checkpoint, `energy` the energy of
the last of them, which the first sweep is compared with, and `nsweeps` the sweeps of the
phase, which a stop on the tolerance records as done.
"""
mutable struct DmrgObserver <: AbstractObserver
    sim::Simulation
    measurements::Union{Vector, Pair}
    period::Int
    tol::Real
    energy::Union{Nothing, Float64}
    done::Int
    nsweeps::Int
    DmrgObserver(sim, measurements, period, tol, done = 0; nsweeps = typemax(Int), energy = nothing) =
        new(sim, measurements, period, tol, energy, done, nsweeps)
end

"""
    sweep_commit!(sim, state, time, sweep; energy)

the end of a sweep of the phase being run, once its measurements and its log are written:
commit it, write the commit if a checkpoint is due or a stop is asked for, and return whether
the solver has to stop.

The order is what keeps a checkpoint and the outputs in step, so it is here, once, rather than
in each observer: what is written after the commit is written again by the resumed run, and
what is written before it is kept. The sweep is committed only in a phase that has read its
resume point, see `resume_sweeps!`; in any other, the commit stays the start of the phase,
while a stop and an interrupt are honoured all the same.
"""
function sweep_commit!(sim::Simulation, state::State, t::Number, sweep::Int; energy = nothing)
    c = sim.checkpoint
    if c.sweeps
        k = c.last
        # a solver run on a simulation of one's own, outside `runTMS`, has no phase around it
        phase, phase_time = isnothing(k) ? (1, t) : (k.phase, k.phase_time)
        commit!(c, sim.outputs, phase, sweep, phase_time, t, state; energy)
    end
    stop = stop_requested(c)
    if stop || checkpoint_due(c)
        write_checkpoint(c, sim.outputs)
    end
    c.stopping = stop
    return stop
end

"""
    sweep_done!(observer; sweep, state, current_time, kwargs...)

the end of a sweep of `tdvp` or `approx_W`, once all its work is done, expansion and
hermitianization included, and whether the solver has to stop. The observers of the package
write the measurements and the log of the sweep and then commit it, see `sweep_commit!`, so
that a checkpoint holds exactly what is written up to it. Any other observer is handed the
sweep as ITensorMPS hands it one, `measure!` then `checkdone!`.
"""
function sweep_done!(o; kwargs...)
    measure!(o; kwargs...)
    return checkdone!(o; kwargs...)
end

# the two evolutions differ only in the operators they apply, which the first sweep logs
function evolution_sweep_done!(o::Union{TdvpObserver, ApproxWObserver}, sweep, current_time,
                               state, label, mpo, operators)
    st = State(o.sim.state, state)
    if sweep_due(o.period, sweep)
        output(Simulation(o.sim, st, current_time), o.measurements; sweep)
    end
    if sweep == 1
        log_msg(o.sim, "$label: maxlinkdim=$(maxlinkdim(mpo)), memory=$(Base.summarysize(operators))")
    end
    log_msg(o.sim, "sim_time $(round(current_time; digits=8))")
    return sweep_commit!(o.sim, st, current_time, sweep)
end

sweep_done!(o::TdvpObserver; sweep, current_time, state, mpo, kwargs...) =
    evolution_sweep_done!(o, sweep, current_time, state, "Tdvp MPO", mpo, mpo)

sweep_done!(o::ApproxWObserver; sweep, current_time, state, mpos, kwargs...) =
    evolution_sweep_done!(o, sweep, current_time, state, "Approx_W MPOS", mpos[1], mpos)

function checkdone!(o::DmrgObserver; energy, sweep, psi, kwargs...)
    # ITensorMPS counts the sweeps of a resumed run from 1 again, and what is measured, logged
    # and checkpointed goes by those of the phase, as it does for tdvp and approx_W. The
    # energy is compared with the sweep before, which a resumed phase reads from the
    # checkpoint, so that it stops where the uninterrupted run stops
    s = sweep + o.done
    stop = !isnothing(o.energy) && abs(o.energy - energy) < o.tol
    # a stop asked for is not a reason to measure: the uninterrupted run did not measure a
    # sweep its period skips, and the resumed one continues after it
    if stop || sweep_due(o.period, s)
        st = normalize(State(o.sim.state, psi))
        output(Simulation(o.sim, st), o.measurements; energy, sweep = s)
    end
    log_msg(o.sim, "sweep $s")
    o.energy = energy
    # a dmrg sweep does not change the simulation time, so the state is committed as it is
    # and the sweep count is what a resume needs. The state is copied, dmrg going on with the
    # next sweep in the same MPS, and a stop on the tolerance records the phase as done, not
    # to be run again
    if sweep_commit!(o.sim, State(o.sim.state, copy(psi)), o.sim.time, stop ? o.nsweeps : s; energy)
        stop = true
    end
    return stop
end
