export Limits, no_limits, Phases, CreateState, LoadState, SaveState, ToMixed, Tdvp, ApproxW, Evolve, Gates, Dmrg

"""
A type to hold MPS limits

# Fields
- `cutoff`: the cutoff under which singular values are neglected
- `maxdim`: the maximum bond dimension
"""
@kwdef struct Limits
    cutoff::Float64
    maxdim::Int
end
  
const no_limits = Limits(cutoff = 0, maxdim = typemax(Int))


"""
A phase type to create the simulation state

# Fields

- `type`: the type of state to create Pure or Mixed
- `system`: a System object to describe the system (see `System`)
- `state`: a description of the state (or a State object to randomize)
- `randomize`: the link dimension for the random state to create (default 0 for no randomizing)

# Examples
    CreateState(type = Pure, system = System(10, "Qubit"), state = "Up")
    CreateState(type = Mixed, system = System(3, "Qubit"), state = ["Up", "Dn", "Up])
    CreateState(type = Pure, system = System(10, "Qubit"), randomize = 50)
    CreateState(type = Pure, system = System(10, "Qubit"), state = "Up", randomize = 50)
"""
@kwdef struct CreateState
    name::String = "Creating state"
    time_start = 0.
    final_measures = []
    type::Union{Nothing, TPure, TMixed} = nothing
    system::Union{Nothing, System} = nothing
    state::Union{Nothing, String, Vector{String}, State} = nothing
    randomize::Int = 0
end

show(io::IO, s::CreateState) = 
    print(io,
        """
        CreateState(
            name = $(repr(s.name)),
            time_start = $(s.time_start),
            final_measures = $(s.final_measures),
            type = $(s.type),
            system = $(s.system),
            state = $(repr(s.state)),
            randomize = $(s.randomize))""")

"""
A phase type to save the state to disk

SaveState(file = "myfile")
"""
@kwdef struct SaveState
    name::String = "Saving state"
    time_start = nothing
    final_measures = []
    file::String
    group_name::String = "state"
end

show(io::IO, s::SaveState) = 
    print(io,
    """
    SaveState(
        name = $(repr(s.name)),
        time_start = $(s.time_start),
        final_measures = $(s.final_measures),
        file = $(repr(s.file)),
        group_name = $(repr(s.group_name)))"""
    )


"""
A phase type to load the state from disk

LoadState(file = "myfile")
"""
@kwdef struct LoadState
    name::String = "Loading state"
    time_start = nothing
    final_measures = []
    file::String
    group_name::String = "state"
    limits::Limits = no_limits
end

show(io::IO, s::LoadState) = 
    print(io,
    """
    LoadState(
        name = $(repr(s.name)),
        time_start = $(s.time_start),
        final_measures = $(s.final_measures),
        file = $(repr(s.file)),
        group_name = $(repr(s.group_name)),
        limits = $(s.limits))"""
    )

"""
A phase type to switch to mixed representation

# Examples
    ToMixed()
    ToMixed(limits = Limits(cutoff = 1e-10, maxdim = 10))
"""
@kwdef struct ToMixed
    name::String = "Switching to mixed state representation"
    time_start = nothing
    final_measures = []
    limits::Limits = no_limits
end

show(io::IO, s::ToMixed) = 
    print(io,
    """
    ToMixed(
        name = $(repr(s.name)),
        time_start = $(s.time_start),
        final_measures = $(s.final_measures),
        limits = $(s.limits))"""
    )



"""
An algorithm type for `Evolve`
"""
struct Tdvp end

"""
An algorithm type for `Evolve`

This corresponds to time evolution with exponential approximation WI or WII combined to obtained approximation of the given order

# Examples
    ApproxW(order = 2, w = 2)     order 2, WII
    ApproxW(order = 4, w = 1)     order 4, WI
"""
@kwdef struct ApproxW
    order
    w = 2
end

show(io::IO, s::ApproxW) =
    print(io, "ApproxW(order = $(s.order), w = $(s.w))")

"""
A phase type for time evolution

# Examples

    Evolve(duration = 2., time_step = 0.1, algo = Tdvp, evolver = -im*(Z(1)Z(2)+(Z(2)Z(3))), measures = [X, Y, Z])

# Fields
- `limits`: a Limits object to set cutoff and maxdim (see `Limits`)
- `duration`: the duration of the time evolution
- `time_step`: the time step
- `algo`: the algorithm used (one of `Tdvp` or `ApproxW`)
- `evolver`: the hamiltonian (evolver = -im * H) with a possible dissipator (evolver = -im * H + D)
- `measures`: the measurement to make (default [])
- `measures_period`: number of time steps between measurments (default 1)
"""
@kwdef struct Evolve
    name::String = "Time evolution"
    time_start = nothing
    final_measures = []
    limits::Limits = no_limits
    duration::Number
    time_step::Number
    algo
    evolver
    measures_period::Int = 1
    measures = []
    corrections::Int = 0
end

show(io::IO, s::Evolve) = 
    print(io,
    """
    Evolve(
        name = $(repr(s.name)),
        time_start = $(s.time_start),
        final_measures = $(s.final_measures),
        limits = $(s.limits),
        duration = $(s.duration),
        time_step = $(s.time_step),
        algo = $(s.algo),
        evolver = $(s.evolver),
        measures_period = $(s.measures_period),
        measures = $(s.measures),
        corrections = $(s.corrections))"""
    )

"""
A phase type for applying gates

# Examples
    Gates(gates = CNOT(1, 3)*CZ(2,4), limits = Limits(cutoff=1e-10, maxdim = 20))
"""
@kwdef struct Gates
    name::String = "Applying gates"
    time_start = nothing
    final_measures = []
    gates::ProdLit
    limits::Limits
  end

show(io::IO, s::Gates) = 
    print(io,
    """
    Gates(
        name = $(repr(s.name)),
        time_start = $(s.time_start),
        final_measures = $(s.final_measures),
        limits = $(s.limits),
        gates = $(s.gates))"""
    )

"""
A phase type for optimizing with Dmrg

# Examples
    Dmrg(hamiltonian = X(1)X(2), nsweeps = 10, cutoff = 1e-10, maxim = [10, 20, 30], tolerance =1e-6)
"""
@kwdef struct Dmrg
    name::String = "Dmrg optimization"
    time_start = nothing
    final_measures = []
    hamiltonian::Lits
    cutoff::Float64
    maxdim::Vector{Int}
    nsweeps::Int
    measures = []
    measures_period::Int = 1
    tolerance::Number = 0.
end

show(io::IO, s::Dmrg) = 
    print(io,
    """
    Dmrg(
        name = $(repr(s.name)),
        time_start = $(s.time_start),
        final_measures = $(s.final_measures),
        nsweeps = $(s.nsweeps),
        hamiltonian = $(s.hamiltonian),
        cutoff = $(s.cutoff),
        maxdim = $(s.maxdim),
        measures_period = $(s.measures_period),
        measures = $(s.measures),
        tolerance = $(s.tolerance))"""
    )


"""
    Phases = Union{CreateState, SaveState, LoadState, ToMixed, Evolve, Gates, Dmrg}

A type that contains all possible phase types for SimData and runTMS.
Each of the types contains at least the three following fields (like SimData).

- `name`: the name of the phase
- `time_start`: the simulation time to use at the start of the phase
- `final_measures`: the measurements to make at the end of the phase see `Measure` and `output`
"""
Phases = Union{CreateState, SaveState, LoadState, ToMixed, Evolve, Gates, Dmrg}