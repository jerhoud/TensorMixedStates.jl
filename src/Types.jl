export Limits, no_limits, CreateState

@kwdef struct Limits
    cutoff::Float64
    maxdim::Int
end
  
const no_limits = Limits(cutoff = 0, maxdim = typemax(Int))

@kwdef struct CreateState
    name::String = "Creating state"
    end_measures = []
    type::Union{Type{Pure}, Type{Mixed}}
    system
    state::Union{Nothing, String, Vector{String}} = nothing
    randomize::Int = 0
end

@kwdef struct ToMixed
    name::String = "Switching to mixed state representation"
    end_measures = []
    limits::Limits = no_limits
end

struct Tdvp end

@kwdef struct ApproxW
    order
    w = 2
end

@kwdef struct Evolve
    name::String = "Time evolution"
    end_measures = []
    Limits::Limits = no_limits
    time_start = nothing
    time::Number
    time_step::Number
    algo
    evolver
    measures_periodicity::Int = 1
    measures = []
    corrections::Int = 0
end