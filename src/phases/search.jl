# What the phases searching a state by dmrg, GroundState and SteadyState, share: running the
# search from the sweep a resumed run had reached.

"""
    run_search(solve, sim, phase, what, final_line)

run a phase that searches a state by dmrg, `GroundState` or `SteadyState`, `solve(sim;
options...)` calling its solver from the sweep a resumed run had reached, logging `what` and,
unless the run stops for a checkpoint, `final_line(e)` with the value reached, that recorded by
the checkpoint when no sweep is left
"""
function run_search(solve, sim::Simulation, phase, what::String, final_line)
    done, e = resume_step(sim)
    if done < phase.nsweeps
        log_message(sim, "$what with $(phase.nsweeps - done) sweeps of Dmrg")
        e, sim = solve(sim; phase.nsweeps, first_sweep = done + 1, phase.limits,
            observer! = DmrgObserver(sim, phase.measurements, phase.measurements_period,
                                     phase.tol, done; phase.nsweeps, energy = e))
    end
    # a search stopped for a checkpoint is not done, and its resume writes the line
    if !stopped(sim)
        log_message(sim, final_line(e))
    end
    return sim
end
