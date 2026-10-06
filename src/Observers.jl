# The observers given to the solvers, which output measurements every given number of steps,
# log the progress, and stop the solver when the simulation is asked to stop.

export TdvpObserver, DmrgObserver, ApproxWObserver, ThermalObserver

"""
    TdvpObserver(sim, measurements, period)

an observer for `tdvp` that outputs `measurements`, given as `output` takes them, every
`period` steps, at the time each step reaches; a `period` below one means never. It also
logs the simulation time of every step and, within `runTMS`, stops the evolution when the
simulation is asked to stop.

# Examples

    tdvp(evolver, 1., sim; nsteps = 10, observer! = TdvpObserver(sim, "data" => [X, Z(1)], 2))
"""
struct TdvpObserver <: AbstractObserver
    sim::Simulation
    measurements::Vector
    period::Int
    TdvpObserver(sim, measurements, period) = new(sim, measure_sets(measurements), period)
end

"""
    ApproxWObserver(sim, measurements, period)

an observer for `approx_W` that outputs `measurements`, given as `output` takes them, every
`period` steps, at the time each step reaches; a `period` below one means never. It also
logs the simulation time of every step and, within `runTMS`, stops the evolution when the
simulation is asked to stop.

# Examples

    approx_W(evolver, 1., sim; order = 4, nsteps = 10,
             observer! = ApproxWObserver(sim, "data" => [X, Z(1)], 2))
"""
struct ApproxWObserver <: AbstractObserver
    sim::Simulation
    measurements::Vector
    period::Int
    ApproxWObserver(sim, measurements, period) = new(sim, measure_sets(measurements), period)
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
    measurements::Vector
    period::Int
    ThermalObserver(sim, measurements, period) = new(sim, measure_sets(measurements), period)
end

"""
    DmrgObserver(sim, measurements, period, tol[, done; nsweeps, energy])

an observer for `dmrg` and `steady_state` that stops the search when the energy changes by
less than `tol` from one sweep to the next. It outputs `measurements`, given as `output`
takes them, on the normalized state every `period` sweeps and on the sweep the tolerance
stops at; a `period` below one means only then. It also logs every sweep and, within
`runTMS`, stops the search when the simulation is asked to stop.

The other arguments serve a search resumed from a checkpoint:

- `done`: the number of sweeps already done
- `energy`: the energy of the last of them, which the first sweep is compared with
- `nsweeps`: the sweeps of the whole phase, which a stop on the tolerance records as done

# Examples

    dmrg(H, sim; nsweeps = 10, observer! = DmrgObserver(sim, "data" => [X, Z(1)], 1, 1e-8))
"""
mutable struct DmrgObserver <: AbstractObserver
    sim::Simulation
    measurements::Vector
    period::Int
    tol::Real
    energy::Union{Nothing, Float64}
    done::Int
    nsweeps::Int
    DmrgObserver(sim, measurements, period, tol, done = 0; nsweeps = typemax(Int), energy = nothing) =
        new(sim, measure_sets(measurements), period, tol, energy, done, nsweeps)
end

"""
    sweep_commit!(sim, state, time, sweep; carried)

close a sweep once its measurements and its log are written: commit it with the value
`carried` to its next sweep, see `Commit`, write the checkpoint if one is due or a stop is
asked for, and return whether the solver has to stop. This order keeps a checkpoint and the
outputs in step: what is written before the commit is kept on resume, what is written after
it is written again. The sweep is committed only in a phase that has read its resume point, see
`resume_step`, a stop being honoured in any phase.
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

the end of a sweep of `tdvp`, `approx_W` or `thermal_state`, all its work done, returning
whether the solver has to stop. The observers of the package write the measurements and the
log of the sweep, then commit it, see `sweep_commit!`; any other observer gets `measure!`
then `checkdone!`, as from ITensorMPS.
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
    # log_trace is carried in the commit for a resumed computation to go on from it
    return sweep_commit!(o.sim, st, o.sim.time, sweep; carried = log_trace)
end

function checkdone!(o::DmrgObserver; energy, sweep, psi, kwargs...)
    # ITensorMPS counts the sweeps of a resumed search from 1: those of the phase are measured,
    # logged and checkpointed, and the energy of the sweep before, read from the checkpoint,
    # makes a resumed search stop where an uninterrupted one does
    s = sweep + o.done
    stop = !isnothing(o.energy) && abs(o.energy - energy) < o.tol
    # normalized as the solver normalizes what it returns, an eigenvector of L†L having a sign
    # of its own: a search resumed after its last sweep hands this state on without running
    # the solver. normalize makes a new MPS, leaving that of dmrg untouched
    st = normalize(State(o.sim.state, psi))
    # a stop asked for is no reason to measure: the resumed run continues after this sweep,
    # which the uninterrupted run measures only when due
    if stop || sweep_due(o.period, s)
        output(Simulation(o.sim, st), o.measurements; energy, sweep = s)
    end
    log_message(o.sim, "sweep $s")
    o.energy = energy
    # a stop on the tolerance records the phase as done, not to be run again on resume
    if sweep_commit!(o.sim, st, o.sim.time, stop ? o.nsweeps : s; carried = energy)
        stop = true
    end
    return stop
end
