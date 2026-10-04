# The observers given to tdvp, dmrg and approx_W, which output measurements every given number of
# steps, log the progress, and stop the algorithm when the simulation is asked to stop.

export TdvpObserver, DmrgObserver, ApproxWObserver, ThermalObserver

"""
    TdvpObserver(sim, measurements, period)

an observer for `tdvp` that outputs `measurements`, given as `output` takes them, every
`period` steps, at the time each step reaches; a `period` below one means never. It also
logs the simulation time of every step and, within `runTMS`, stops the evolution when the
simulation is asked to stop.

# Examples

    tdvp(evolver, 1., sim; nsweeps = 10, observer! = TdvpObserver(sim, "data" => [X, Z(1)], 2))
"""
struct TdvpObserver <: AbstractObserver
    sim::Simulation
    measurements::Union{Vector, Pair}
    period::Int
end

"""
    ApproxWObserver(sim, measurements, period)

an observer for `approx_W` that outputs `measurements`, given as `output` takes them, every
`period` steps, at the time each step reaches; a `period` below one means never. It also
logs the simulation time of every step and, within `runTMS`, stops the evolution when the
simulation is asked to stop.

# Examples

    approx_W(evolver, 1., sim; order = 4, nsweeps = 10,
             observer! = ApproxWObserver(sim, "data" => [X, Z(1)], 2))
"""
struct ApproxWObserver <: AbstractObserver
    sim::Simulation
    measurements::Union{Vector, Pair}
    period::Int
end

"""
    ThermalObserver(sim, measurements, period)

an observer for `thermal_state` that outputs `measurements`, given as `output` takes them,
every `period` steps, the symbols `:beta` and `:log_trace` taking the inverse temperature and
the logarithm of the trace each step reaches; a `period` below one means never. It also logs
the inverse temperature of every step and, within `runTMS`, stops the computation when the
simulation is asked to stop.

# Examples

    thermal_state(H, 2.0, sim; nsteps = 20,
                  observer! = ThermalObserver(sim, "data" => [:beta, :log_trace, Z(1)], 1))
"""
struct ThermalObserver <: AbstractObserver
    sim::Simulation
    measurements::Union{Vector, Pair}
    period::Int
end

"""
    DmrgObserver(sim, measurements, period, tol[, done; nsweeps, energy])

an observer for `dmrg` and `steady_state` that stops the search when the energy changes by
less than `tol` from one sweep to the next. It outputs `measurements`, given as `output`
takes them, on the normalized state every `period` sweeps and on the sweep the tolerance
stops at; a `period` below one means only then. It also logs every sweep and, within
`runTMS`, stops the search when the simulation is asked to stop.

The other arguments serve a search resumed from a checkpoint: `done` is the number of sweeps
already done, `energy` the energy of the last of them, which the first sweep is compared with,
and `nsweeps` the sweeps of the whole phase, which a stop on the tolerance records as done.

# Examples

    dmrg(H, sim; nsweeps = 10, observer! = DmrgObserver(sim, "data" => [X, Z(1)], 1, 1e-8))
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
    sweep_commit!(sim, state, time, sweep; carried)

close a sweep of the phase being run, once its measurements and its log are written: commit
it, with the value `carried` the phase carries to its next sweep, see `Commit`, write the commit
if a checkpoint is due or a stop is asked for, and return whether the solver has to stop.

This order keeps a checkpoint and the outputs in step, so it lives here rather than in each
observer: what is written after the commit is written again by the resumed run, and what is
written before it is kept. The sweep is committed only in a phase that has read its resume
point, see `resume_step`; in any other the commit stays the start of the phase, while a
stop and an interrupt are honoured all the same.
"""
function sweep_commit!(sim::Simulation, state::AbstractState, t::Number, sweep::Int;
                       carried = nothing)
    c = sim.checkpoint
    if c.sweeps
        k = c.last
        # a solver run on a simulation of one's own, outside `runTMS`, has no phase around it
        phase, phase_time = isnothing(k) ? (1, t) : (k.phase, k.phase_time)
        commit!(c, sim.outputs, phase, sweep, phase_time, t, state; carried)
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
hermitianization included, returning whether the solver has to stop. The observers of the
package write the measurements and the log of the sweep and then commit it, see
`sweep_commit!`, so that a checkpoint holds exactly what is written up to it. Any other
observer is handed the sweep as ITensorMPS hands it one, `measure!` then `checkdone!`.
"""
function sweep_done!(o; kwargs...)
    measure!(o; kwargs...)
    return checkdone!(o; kwargs...)
end

"""
    evolution_sweep_done!(observer, sweep, current_time, state, label, mpo, operators)

`sweep_done!` for a `TdvpObserver` or an `ApproxWObserver`. The two differ only in the
operators they apply, which the first sweep logs under `label`.
"""
function evolution_sweep_done!(o::Union{TdvpObserver, ApproxWObserver}, sweep, current_time,
                               state, label, mpo, operators)
    st = State(o.sim.state, state)
    if sweep_due(o.period, sweep)
        output(Simulation(o.sim, st, current_time), o.measurements; sweep)
    end
    if sweep == 1
        log_message(o.sim, "$label: maxlinkdim=$(maxlinkdim(mpo)), memory=$(Base.summarysize(operators))")
    end
    log_message(o.sim, "sim_time $(round(current_time; digits=8))")
    return sweep_commit!(o.sim, st, current_time, sweep)
end

sweep_done!(o::TdvpObserver; sweep, current_time, state, mpo, kwargs...) =
    evolution_sweep_done!(o, sweep, current_time, state, "Tdvp MPO", mpo, mpo)

sweep_done!(o::ApproxWObserver; sweep, current_time, state, mpos, kwargs...) =
    evolution_sweep_done!(o, sweep, current_time, state, "Approx_W MPOS", mpos[1], mpos)

function sweep_done!(o::ThermalObserver; sweep, state, beta, log_trace, kwargs...)
    st = State(o.sim.state, state)
    if sweep_due(o.period, sweep)
        output(Simulation(o.sim, st), o.measurements; sweep, beta, log_trace)
    end
    log_message(o.sim, "beta $(round(beta; digits=8))")
    # the simulation time does not move, and the logarithm of the trace is carried in the
    # commit, as the energy of a dmrg sweep is, for a resumed computation to go on from it
    return sweep_commit!(o.sim, st, o.sim.time, sweep; carried = log_trace)
end

function checkdone!(o::DmrgObserver; energy, sweep, psi, kwargs...)
    # ITensorMPS counts the sweeps of a resumed run from 1 again, and what is measured, logged
    # and checkpointed goes by those of the phase, as it does for tdvp and approx_W. The
    # energy is compared with the sweep before, which a resumed phase reads from the
    # checkpoint, so that it stops where the uninterrupted run stops
    s = sweep + o.done
    stop = !isnothing(o.energy) && abs(o.energy - energy) < o.tol
    # normalized as the solver normalizes what it returns: an eigenvector of L†L has norm one
    # and a sign of its own, and a steady state checkpointed on its last sweep, whose resume
    # does not run the solver, was handed on and saved as it was, of trace -1.33. Normalizing
    # makes a new MPS, dmrg going on with the next sweep in its own
    st = normalize(State(o.sim.state, psi))
    # a stop asked for is not a reason to measure: the uninterrupted run did not measure a
    # sweep its period skips, and the resumed one continues after it
    if stop || sweep_due(o.period, s)
        output(Simulation(o.sim, st), o.measurements; energy, sweep = s)
    end
    log_message(o.sim, "sweep $s")
    o.energy = energy
    # a dmrg sweep does not change the simulation time, so the sweep count is what a resume
    # needs, and a stop on the tolerance records the phase as done, not to be run again
    if sweep_commit!(o.sim, st, o.sim.time, stop ? o.nsweeps : s; carried = energy)
        stop = true
    end
    return stop
end
