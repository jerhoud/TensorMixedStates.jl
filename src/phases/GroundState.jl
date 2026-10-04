# The GroundState phase, which searches the ground state of a hamiltonian by dmrg, and its
# deprecated name Dmrg.

export GroundState, Dmrg

"""
    GroundState(; hamiltonian, limits, nsweeps, noise, tolerance, measures, options...)

a phase that searches the ground state of a hamiltonian by dmrg, see `dmrg`, on a pure state.

# Fields

- `name`, `time_start`, `final_measures`: the fields every phase has, see `AbstractPhase`
- `hamiltonian`: the hamiltonian whose ground state is searched
- `limits`: constraints on the state, see `Limits`, required
- `nsweeps`: the maximum number of sweeps, required
- `noise`: the noise to apply, a number or one value per sweep (default 0)
- `krylov`: the parameters of the Krylov search of each local step, see `Krylov` (default
  `Krylov()`)
- `measures`: the measurements to make during the search, see `output` (default `[]`)
- `measures_period`: the number of sweeps between two measurements (default 1)
- `tolerance`: the search stops when the energy changes by less than this from one sweep to
  the next (default 0, never)

# Examples

    GroundState(hamiltonian = X(1)X(2), nsweeps = 10,
                limits = Limits(cutoff = 1e-10, maxdim = [10, 20, 30]), tolerance = 1e-6)
"""
@kwdef struct GroundState <: AbstractPhase
    name::String = "Ground state computation using Dmrg"
    time_start::Union{Nothing, Number} = nothing
    final_measures = []
    hamiltonian::IndexedOp{Pure}
    limits::Limits
    nsweeps::Int
    noise::Union{Float64, Vector{Float64}} = 0.
    krylov::Krylov = Krylov()
    measures = []
    measures_period::Int = 1
    tolerance::Real = 0.
    # a noise is a real number, an integer one included, stored as the Float64 dmrg takes
    GroundState(name, time_start, final_measures, hamiltonian, limits, nsweeps,
                noise::Union{Real, AbstractVector{<:Real}}, krylov, measures, measures_period,
                tolerance) =
        new(name, time_start, final_measures, hamiltonian, limits, nsweeps,
            noise isa AbstractVector ? Vector{Float64}(noise) : Float64(noise),
            krylov, measures, measures_period, tolerance)
end

# the docstring goes through `@doc` rather than sitting above the call, because the macro
# expands to a toplevel block and a docstring cannot be attached to one. It comes before the
# deprecation: after it, attaching the docstring used the deprecated name and warned
@doc """
    Dmrg

deprecated, use [`GroundState`](@ref) instead.

`Dmrg` is an alias of `GroundState`, which works until it is removed in a future version. Its
deprecation warning only shows with `--depwarn=yes`, as when running the tests, so this line
is the notice.
""" Dmrg

Base.@deprecate_binding Dmrg GroundState false ", use GroundState instead."

run_phase(sim::Simulation, phase::GroundState) =
    run_search((sim; kwargs...) -> dmrg(phase.hamiltonian, sim; phase.noise, phase.krylov,
                                        kwargs...),
               sim, phase, "Optimizing state", e -> "Done, dmrg final energy is $e")
