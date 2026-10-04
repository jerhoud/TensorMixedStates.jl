# The SteadyState phase, which searches the steady state of a Lindbladian by dmrg.

export SteadyState

"""
    SteadyState(; lindbladian, limits, nsweeps, noise, tolerance, measures, options...)

a phase that searches the steady state of a Lindbladian, see `steady_state`, on a mixed state.

# Fields

- `name`, `time_start`, `final_measures`: the fields every phase has, see `AbstractPhase`
- `lindbladian`: the Lindbladian ``L`` whose steady state is searched, of the form
  `-im * hamiltonian + dissipators`
- `mpo_limits`: the truncation of the MPO of ``L^\\dagger L`` (default `Limits()`)
- `mpo_algo`: the algorithm computing ``L^\\dagger L``, `"naive"` (default) or `"zipup"`
- `limits`: constraints on the state, see `Limits`, required
- `nsweeps`: the maximum number of sweeps, required
- `noise`: the noise to apply, a number or one value per sweep (default 0)
- `krylov`: the parameters of the Krylov search of each local step, see `Krylov` (default
  `Krylov(dim = 8, maxiter = 3)`, see `steady_state`)
- `measures`: the measurements to make during the search, see `output` (default `[]`)
- `measures_period`: the number of sweeps between two measurements (default 1)
- `tolerance`: the search stops when the dmrg energy changes by less than this from one sweep
  to the next (default 0, never)

# Examples

    SteadyState(
        lindbladian = -im * hamiltonian + dissipators,
        limits = Limits(cutoff = 1e-20, maxdim = [10, 20, 50]),
        nsweeps = 10,
        tolerance = 1e-5,
    )
"""
@kwdef struct SteadyState <: AbstractPhase
    name::String = "Steady state optimization"
    time_start::Union{Nothing, Number} = nothing
    final_measures = []
    lindbladian::IndexedOp{Mixed}
    mpo_limits::Limits = Limits()
    mpo_algo::String = "naive"
    limits::Limits
    nsweeps::Int
    noise::Union{Float64, Vector{Float64}} = 0.
    krylov::Krylov = Krylov(dim = 8, maxiter = 3)
    measures = []
    measures_period::Int = 1
    tolerance::Real = 0.
    # a noise is stored as the Float64 dmrg takes, as in GroundState, and mpo_algo is checked
    # when the phase is written, as the fields of ApproxW
    function SteadyState(name, time_start, final_measures, lindbladian, mpo_limits, mpo_algo,
                         limits, nsweeps, noise::Union{Real, AbstractVector{<:Real}}, krylov,
                         measures, measures_period, tolerance)
        check_mpo_algo(mpo_algo)
        return new(name, time_start, final_measures, lindbladian, mpo_limits, mpo_algo, limits,
                   nsweeps, noise isa AbstractVector ? Vector{Float64}(noise) : Float64(noise),
                   krylov, measures, measures_period, tolerance)
    end
end

function run_phase(sim::Simulation, phase::SteadyState)
    if sim.state isa State{Pure}
        error("state must be in mixed representation for computing steady state")
    end
    return run_search(
        (sim; kwargs...) -> steady_state(phase.lindbladian, sim; phase.mpo_limits,
                                         phase.mpo_algo, phase.noise, phase.krylov, kwargs...),
        sim, phase, "Searching for steady state",
        e -> "Done, dmrg final value is $e (0 for steady state)")
end
