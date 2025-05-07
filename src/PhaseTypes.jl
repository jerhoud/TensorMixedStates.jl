export Limits, Phases, Algo, CreateState, LoadState, SaveState, ToMixed, Tdvp, ApproxW, Evolve, Gates, Dmrg, Partial_Trace


"""
A phase type to create the simulation state

# Fields

- `name`: the name of the phase
- `time_start`: the initial simulation time
- `final_measures`: the measurements to make at the end of the phase see `measure` and `output`
- `type`: the type of state to create `Pure()` or `Mixed()`
- `system`: a System object to describe the system (see `System`) (unused if a State object is given)
- `state`: a description of the state (or a State object)
- `randomize`: the link dimension for the random state to create (default 0 for no randomizing)

# Examples
    CreateState(type = Pure(), system = System(10, Qubit()), state = "Up")
    CreateState(type = Mixed(), system = System(3, Qubit()), state = ["Up", "Dn", "Up])
    CreateState(type = Pure(), system = System(10, Qubit()), randomize = 50)
    CreateState(type = Pure(), system = System(10, Qubit()), state = "Up", randomize = 50)
    CreateState(Pure(), 10, Qubit(), "Up")                                      # simple form
    CreateState(Mixed(), [Qubit(), Boson(4), Fermion()], ["Up", "2", "Occ"])    # other simple form
"""
@kwdef struct CreateState
    name::String = "Creating state"
    time_start::Number = 0.
    final_measures = []
    type::Union{Nothing, Pure, Mixed} = nothing
    system::Union{Nothing, System} = nothing
    state = nothing
    randomize::Int = 0
end

CreateState(type, n, site, state; kwargs...) = CreateState(;type, system = System(n, site), state, kwargs...)
CreateState(type, sites, state; kwargs...) = CreateState(;type, system = System(sites), state, kwargs...)

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

not implemented
"""
@kwdef struct SaveState
    name::String = "Saving state"
    time_start::Union{Nothing, Number} = nothing
    final_measures = []
    file::String
    statename::String = "state"
end

show(io::IO, s::SaveState) = 
    print(io,
    """

    SaveState(
        name = $(repr(s.name)),
        time_start = $(s.time_start),
        final_measures = $(s.final_measures),
        file = $(repr(s.file)),
        statename = $(repr(s.statename)))"""
    )


"""
A phase type to load the state from disk

LoadState(file = "myfile")

not implemented
"""
@kwdef struct LoadState
    name::String = "Loading state"
    time_start::Union{Nothing, Number} = nothing
    final_measures = []
    file::String
    statename::String = "state"
    limits::Limits = Limits()
end

show(io::IO, s::LoadState) = 
    print(io,
    """

    LoadState(
        name = $(repr(s.name)),
        time_start = $(s.time_start),
        final_measures = $(s.final_measures),
        file = $(repr(s.file)),
        statename = $(repr(s.statename)),
        limits = $(s.limits))"""
    )

"""
A phase type to switch to mixed representation

# Fields

- `name`: the name of the phase
- `time_start`: the simulation time to use (no much use here)
- `final_measures`: the measurements to make at the end of the phase see `measure` and `output`
' `limits` : constraints on the final state

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

# Examples
    Tdvp()
    Tdvp(n_expand = 5)     tdvp with expansion steps every 5 steps
    Tdvp(n_symmetrize = 3) tdvp, make hermitian every 3 steps
"""
@kwdef struct Tdvp
    n_expand::Int = 0
    n_symmetrize::Int = 0
end

show(io::IO, s::Tdvp) =
    print(io, "Tdvp(n_expand = $(s.n_expand), n_symmetrize = $(s.n_symmetrize))")

"""
An algorithm type for `Evolve`

This corresponds to time evolution with exponential approximation WI or WII combined to obtained approximation of the given order

# Examples
    ApproxW(order = 2)                   order 2, WII
    ApproxW(order = 4, w = 1)            order 4, WI
    ApproxW(order = 4, n_symmetrize = 3) order 4, make hermitian every 3 steps
"""
@kwdef struct ApproxW
    order::Int
    w::Int = 2
    n_symmetrize::Int = 0
end

show(io::IO, s::ApproxW) =
    print(io, "ApproxW(order = $(s.order), w = $(s.w), n_symmetrize = $(s.n_symmetrize))")


Algo = Union{Tdvp, ApproxW}


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
    time_start::Union{Nothing, Number} = nothing
    final_measures = []
    limits::Limits = Limits()
    duration::Number
    time_step::Number
    algo::Algo
    evolver::Union{ExprIndexed, Pair{ExprIndexed, Vector}}
    measures_period::Int = 1
    measures = []
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
        measures = $(s.measures))"""
    )

"""
A phase type for applying gates

- `name`: the name of the phase
- `time_start`: the simulation time to use at the start of the phase
- `final_measures`: the measurements to make at the end of the phase see `measure` and `output`
- `limits`: constraints to enforce at each step of the computation
- `gates`: the gates to apply

# Examples

    Gates(gates = CNOT(1, 3)*CZ(2,4), limits = Limits(cutoff=1e-10, maxdim = 20))
"""
@kwdef struct Gates
    name::String = "Applying gates"
    time_start::Union{Nothing, Number} = nothing
    final_measures = []
    gates::ExprIndexed
    limits::Limits = Limits()
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
    time_start::Union{Nothing, Number} = nothing
    final_measures = []
    hamiltonian::ExprIndexed{Pure}
    limits::Limits
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
a phase type for applying a partial trace

# Examples

    Partial_Trace(trace_positions = [2, 3, 6])
"""
@kwdef struct Partial_Trace
    name::String = "Dmrg optimization"
    time_start::Union{Nothing, Number} = nothing
    final_measures = []
    trace_positions::Union{Nothing, Vector{Int}} = nothing
    keep_positions::Union{Nothing, Vector{Int}} = nothing
end

show(io::IO, s::Partial_Trace) = 
    print(io,
    """
    
    Partial_Trace(
        name = $(repr(s.name)),
        time_start = $(s.time_start),
        final_measures = $(s.final_measures),
        trace_positions = $(s.trace_positions),
        keep_positions = $(s.keep_positions))""")

"""   
    Phases = Union{CreateState, SaveState, LoadState, ToMixed, Evolve, Gates, Dmrg, Partial_Trace}

A type that contains all possible phase types for SimData and runTMS.
Each of the types contains at least the three following fields (like SimData).

- `name`: the name of the phase
- `time_start`: the simulation time to use at the start of the phase
- `final_measures`: the measurements to make at the end of the phase see `measure` and `output`
"""
Phases = Union{CreateState, SaveState, LoadState, ToMixed, Evolve, Gates, Dmrg, Partial_Trace}