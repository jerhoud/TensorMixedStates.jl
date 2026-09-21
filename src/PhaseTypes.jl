export Phases, Algo, CreateState, LoadState, SaveState, ToMixed, Tdvp, ApproxW, Evolve, Gates, GroundState, Dmrg, PartialTrace, SteadyState


"""
A phase type to create the simulation state

# Fields

- `name`: the name of the phase
- `time_start`: the initial simulation time (default 0., use `nothing` to keep the current simulation time)
- `final_measures`: the measurements to make at the end of the phase see `measure` and `output`
- `type`: the type of state to create `Pure()` or `Mixed()`
- `system`: a System object to describe the system (see `System`) (unused if a State object is given)
- `state`: a description of the state (or a State object)
- `randomize`: the link dimension for the random state to create (default 0 for no randomizing)
- `seed`: set the random generator seed for randomize (default nothing)

# Examples
    CreateState(type = Pure(), system = System(10, Qubit()), state = "Up")
    CreateState(type = Mixed(), system = System(3, Qubit()), state = ["Up", "Dn", "Up"])
    CreateState(type = Pure(), system = System(10, Qubit()), randomize = 50)
    CreateState(type = Pure(), system = System(10, Qubit()), state = "Up", randomize = 50)
    CreateState{Pure}(10, Qubit(), "Up")                                      # simple form
    CreateState{Mixed}([Qubit(), Boson(4), Fermion()], ["Up", "2", "Occ"])    # other simple form
"""
@kwdef struct CreateState{R <: PM}
    name::String = "Creating state"
    time_start::Union{Nothing, Number} = 0.
    final_measures = []
    type::R
    system::Union{Nothing, System} = nothing
    state = nothing
    randomize::Int = 0
    seed::Union{Nothing, Int} = nothing
end

CreateState{R}(n, site, state; kwargs...) where R =
    CreateState(;type = R(), system = System(n, site), state, kwargs...)
CreateState{R}(sites, state; kwargs...) where R =
    CreateState(;type = R(), system = System(sites), state, kwargs...)


"""
A phase type to save the state to disk in a hdf5 file (see `save_state`)

SaveState(file = "myfile.h5")
SaveState(file = "myfile.h5", statename = "after_evolution")

several states can be saved in the same file under different `statename`,
saving under a name already present in the file replaces it

# Fields

- `name`: the name of the phase
- `time_start`: the simulation time at the start of the phase (`nothing` keeps the current time)
- `final_measures`: the measurements to make at the end of the phase
- `file`: the name of the hdf5 file to write to
- `statename`: the name under which the state is stored in the file
"""
@kwdef struct SaveState
    name::String = "Saving state"
    time_start::Union{Nothing, Number} = nothing
    final_measures = []
    file::String
    statename::String = "state"
end



"""
A phase type to load the state from a hdf5 file written by `SaveState` (see `load_state`)

LoadState(file = "myfile.h5")
LoadState(file = "myfile.h5", statename = "after_evolution")

# Fields

- `name`: the name of the phase
- `time_start`: the simulation time at the start of the phase (`nothing` keeps the current time)
- `final_measures`: the measurements to make at the end of the phase
- `file`: the name of the hdf5 file to read from
- `statename`: the name under which the state is stored in the file
- `limits`: the truncation applied to the state after loading
"""
@kwdef struct LoadState
    name::String = "Loading state"
    time_start::Union{Nothing, Number} = nothing
    final_measures = []
    file::String
    statename::String = "state"
    limits::Limits = Limits()
end


"""
A phase type to switch to mixed representation

# Fields

- `name`: the name of the phase
- `time_start`: the simulation time to use (no much use here)
- `final_measures`: the measurements to make at the end of the phase see `measure` and `output`
- `limits` : constraints on the final state

# Examples

    ToMixed()
    ToMixed(limits = Limits(cutoff = 1e-10, maxdim = 10))
"""
@kwdef struct ToMixed
    name::String = "Switching to mixed state representation"
    time_start::Union{Nothing, Number} = nothing
    final_measures = []
    limits::Limits = Limits()
end




"""
An algorithm type for `Evolve`

# Examples
    Tdvp()
    Tdvp(n_expand = 5)        # tdvp with expansion steps every 5 steps
    Tdvp(n_hermitianize = 3)  # tdvp, make hermitian every 3 steps
"""
@kwdef struct Tdvp
    n_expand::Int = 0
    n_hermitianize::Int = 0
end

show(io::IO, s::Tdvp) =
    print(io, "Tdvp(n_expand = $(s.n_expand), n_hermitianize = $(s.n_hermitianize))")

"""
An algorithm type for `Evolve`

This corresponds to time evolution with exponential approximation WI or WII combined to obtained approximation of the given order

# Examples
    ApproxW(order = 2)                     # order 2, WII
    ApproxW(order = 4, w = 1)              # order 4, WI
    ApproxW(order = 4, n_hermitianize = 3) # order 4, make hermitian every 3 steps
"""
@kwdef struct ApproxW
    order::Int
    w::Int = 2
    n_hermitianize::Int = 0
end

show(io::IO, s::ApproxW) =
    print(io, "ApproxW(order = $(s.order), w = $(s.w), n_hermitianize = $(s.n_hermitianize))")


"""
    Algo = Union{Tdvp, ApproxW}

the type of the time evolution algorithms accepted by the `algo` field of the `Evolve` phase
"""
const Algo = Union{Tdvp, ApproxW}


"""
A phase type for time evolution

# Examples

    Evolve(duration = 2., time_step = 0.1, algo = Tdvp(), evolver = -im*(Z(1)Z(2)+(Z(2)Z(3))), measures = [X, Y, Z])

# Fields

- `name`: the name of the phase
- `time_start`: the initial simulation time
- `final_measures`: the measurements to make at the end of the phase see `measure` and `output`
- `limits`: a Limits object to set cutoff, maxdim and mindim (see `Limits`)
- `duration`: the duration of the time evolution
- `time_step`: the time step
- `algo`: the algorithm used (one of `Tdvp()` or `ApproxW(...)`)
- `evolver`: the hamiltonian (evolver = -im * H) with a possible dissipator (evolver = -im * H + D)
- `measures`: the measurement to make (default [])
- `measures_period`: number of time steps between measurements (default 1)
"""
@kwdef struct Evolve
    name::String = "Time evolution"
    time_start::Union{Nothing, Number} = nothing
    final_measures = []
    limits::Limits = Limits()
    duration::Number
    time_step::Number
    algo::Algo
    evolver::Union{IndexedOp, Pair}
    measures_period::Int = 1
    measures = []
end


"""
A phase type for applying gates

# Fields
- `name`: the name of the phase
- `time_start`: the simulation time to use at the start of the phase
- `final_measures`: the measurements to make at the end of the phase see `measure` and `output`
- `limits`: constraints to enforce at each step of the computation
- `gates`: the gates to apply

# Examples

    Gates(gates = controlled(X)(1, 3)*controlled(Z)(2, 4), limits = Limits(cutoff=1e-10, maxdim = 20))
"""
@kwdef struct Gates
    name::String = "Applying gates"
    time_start::Union{Nothing, Number} = nothing
    final_measures = []
    gates::IndexedOp
    limits::Limits = Limits()
end


"""
A phase type for computing the ground state using Dmrg

# Examples
    GroundState(hamiltonian = X(1)X(2), nsweeps = 10, limits = Limits(cutoff = 1e-10, maxdim = [10, 20, 30]), tolerance = 1e-6)

# Fields
- `name`: the name of the phase
- `time_start`: the initial simulation time
- `final_measures`: the measurements to make at the end of the phase see `measure` and `output`
- `hamiltonian`: the Hamiltonian whose ground state is requested
- `limits`: the limits on the state
- `nsweeps`: the maximum number of sweeps
- `noise`: the noise to apply (either a number or a vector of numbers, ITensor dmrg documentation)
- `measures`: measurements to make during the computation
- `measures_period`: the interval at which the measurements are made
- `tolerance`: computation is stopped if the progression in the energy between sweeps is lower than this number
"""
@kwdef struct GroundState
    name::String = "Ground state computation using Dmrg"
    time_start::Union{Nothing, Number} = nothing
    final_measures = []
    hamiltonian::IndexedOp{Pure}
    limits::Limits
    nsweeps::Int
    noise::Union{Float64, Vector{Float64}} = 0.
    measures = []
    measures_period::Int = 1
    tolerance::Number = 0.
end

# the docstring goes through `@doc` rather than sitting above the call, because the macro
# expands to a toplevel block and a docstring cannot be attached to one
Base.@deprecate_binding Dmrg GroundState false ", use GroundState instead."

@doc """
    Dmrg

deprecated, use [`GroundState`](@ref) instead.

`Dmrg` is an alias of `GroundState` and goes on working, but it is marked deprecated in the
runtime and will be removed in a future version. The deprecation warning only shows with
`--depwarn=yes`, which is what running the tests does; an ordinary run stays silent, so this
line is the notice.
""" Dmrg



"""
a phase type for applying a partial trace

# Fields
- `name`: the name of the phase
- `time_start`: the initial simulation time
- `final_measures`: the measurements to make at the end of the phase see `measure` and `output`
- `trace_positions`: an array of site numbers on which to trace
- `keep_positions`: an array of site numbers which are not traced (all the others are)

# Examples

    PartialTrace(trace_positions = [2, 3, 6])
    PartialTrace(keep_positions = [1, 4, 5])
"""
@kwdef struct PartialTrace
    name::String = "Computing partial trace"
    time_start::Union{Nothing, Number} = nothing
    final_measures = []
    trace_positions::Union{Nothing, Vector{Int}} = nothing
    keep_positions::Union{Nothing, Vector{Int}} = nothing
end


"""
a phase to compute the steady state of a Lindbladian

# Fields

- `name`: the name of the phase
- `time_start`: the initial simulation time
- `final_measures`: the measurements to make at the end of the phase see `measure` and `output`
- `lindbladian`: the Lindbladian whose steady state is requested (should be of the form -im * hamiltonian + dissipators)
- `mpo_limits`: limits on the resulting MPO (nothing for no truncation)
- `mpo_algo`: algorithm for computing (L+)L: "naive" (default) or "zipup"
- `limits`: limits on the state MPS
- `nsweeps`: maximum number of sweeps
- `measures`: measurements to be made during the computation
- `measures_period`: the interval at which the measurements are made
- `tolerance`: computation is stopped if the progression in the energy between sweeps is lower than this number

# Examples

    SteadyState(
        lindbladian = -im * hamiltonian + dissipators,
        limits = Limits(cutoff = 1e-20, maxdim = [10, 20, 50]),
        nsweeps = 10,
        tolerance = 1e-5,
    )
"""
@kwdef struct SteadyState
    name::String = "Steady state optimization"
    time_start::Union{Nothing, Number} = nothing
    final_measures = []
    lindbladian::IndexedOp{Mixed}
    mpo_limits::Limits = Limits()
    mpo_algo::String = "naive"
    limits::Limits
    nsweeps::Int
    measures = []
    measures_period::Int = 1
    tolerance::Number = 0.
end


"""   
    Phases = Union{CreateState, SaveState, LoadState, ToMixed, Evolve, Gates, GroundState, PartialTrace, SteadyState}

A type that contains all possible phase types for SimData and runTMS.
Each of the types contains at least the three following fields (like SimData).

- `name`: the name of the phase
- `time_start`: the simulation time to use at the start of the phase
- `final_measures`: the measurements to make at the end of the phase see `measure` and `output`
"""
const Phases = Union{CreateState, SaveState, LoadState, ToMixed, Evolve, Gates, GroundState, PartialTrace, SteadyState}

# The phases are printed field by field, read from the type rather than written out one by
# one, so that a field added to a phase shows up in the log and in `prog.jl` without
# anything else to change.
function show(io::IO, s::Phases)
    t = typeof(s)
    print(io, "\n", nameof(t))
    if !isempty(t.parameters)
        print(io, "{", join(nameof.(t.parameters), ", "), "}")
    end
    print(io, "(")
    fs = fieldnames(t)
    for (i, f) in enumerate(fs)
        print(io, "\n    ", f, " = ", repr(getfield(s, f)), i < length(fs) ? "," : ")")
    end
end