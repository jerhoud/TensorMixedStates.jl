# The GroundState phase, which searches the ground state of a hamiltonian by dmrg.

export GroundState

"""
    GroundState(; hamiltonian, limits, nsweeps, noise, tol, measurements, options...)

a phase that searches the ground state of a hamiltonian by dmrg, see `dmrg`, on a pure state.

# Fields

- `name`, `time_start`, `final_measurements`: the fields every phase has, see `AbstractPhase`
- `hamiltonian`: the hamiltonian whose ground state is searched
- `limits`: constraints on the state, see `Limits`, required
- `nsweeps`: the maximum number of sweeps, required
- `noise`: the noise to apply, a number or one value per sweep (default 0)
- `krylov`: the parameters of the Krylov search of each local step, see `Krylov` (default
  `Krylov()`)
- `measurements`: the measurements to make during the search, see `output` (default `[]`)
- `measurements_period`: the number of sweeps between two measurements (default 1)
- `tol`: the search stops when the energy changes by less than this from one sweep to
  the next (default 0, never)

# Examples

    GroundState(hamiltonian = X(1)X(2), nsweeps = 10,
                limits = Limits(cutoff = 1e-10, maxdim = [10, 20, 30]), tol = 1e-6)
"""
@kwdef struct GroundState <: AbstractPhase
    name::String = "Ground state computation using Dmrg"
    time_start::Union{Nothing, Number} = nothing
    final_measurements = []
    hamiltonian::IndexedOp{Pure}
    limits::Limits
    nsweeps::Int
    noise::Union{Float64, Vector{Float64}} = 0.
    krylov::Krylov = Krylov()
    measurements = []
    measurements_period::Int = 1
    tol::Real = 0.
    # a noise is a real number, an integer one included, stored as the Float64 dmrg takes
    GroundState(name, time_start, final_measurements, hamiltonian, limits, nsweeps,
                noise::Union{Real, AbstractVector{<:Real}}, krylov, measurements,
                measurements_period, tol) =
        new(name, time_start, final_measurements, hamiltonian, limits, nsweeps,
            noise isa AbstractVector ? Vector{Float64}(noise) : Float64(noise),
            krylov, measurements, measurements_period, tol)
end

run_phase(sim::Simulation, phase::GroundState) =
    run_search((sim; kwargs...) -> dmrg(phase.hamiltonian, sim; phase.noise, phase.krylov,
                                        kwargs...),
               sim, phase, "Optimizing state", e -> "Done, dmrg final energy is $e")
