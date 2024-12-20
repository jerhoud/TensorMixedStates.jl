export Limits, no_limits, CreateState, LoadState, SaveState, ToMixed, Tdvp, ApproW, Evolve, Gates, Dmrg

@kwdef struct Limits
    cutoff::Float64
    maxdim::Int
end
  
const no_limits = Limits(cutoff = 0, maxdim = typemax(Int))

@kwdef struct CreateState
    name::String = "Creating state"
    time_start = 0.
    final_measures = []
    type::Union{Nothing, TPure, TMixed} = nothing
    system::Union{Nothing, System} = nothing
    state::Union{Nothing, String, Vector{String}, State} = nothing
    randomize::Int = 0
end

@kwdef struct SaveState
    name::String = "Saving state"
    time_start = nothing
    final_measures = []
    file::String
    group_name::String = "state"
end
  
@kwdef struct LoadState
    name::String = "Loading state"
    time_start = nothing
    final_measures = []
    file::String
    group_name::String = "state"
    limits::Limits = no_limits
end
  
@kwdef struct ToMixed
    name::String = "Switching to mixed state representation"
    time_start = nothing
    final_measures = []
    limits::Limits = no_limits
end

struct Tdvp end

@kwdef struct ApproxW
    order
    w = 2
end

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

@kwdef struct Gates
    name::String = "Applying gates"
    time_start = nothing
    final_measures = []
    gates::ProdLit
    limits::Limits
  end

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
