export RandomState

"""
    RandomState{Pure|Mixed}([eltype, ]::System, linkdims::Int)
    RandomState{Pure|Mixed}([eltype, ]::Int, ::AbstractSite, linkdims::Int)
    RandomState{Pure|Mixed}([eltype, ]::Vector{<:AbstractSite}, linkdims::Int)
    RandomState([eltype, ]::State, linkdims::Int)

Return a random state, with the specified link dimension.
If a State is given, randomize the given state.

A mixed one is built as a purification on twice the system, whose mixed representation
squares the link dimension, so the pure state it starts from is built one size up and the
result truncated back to what was asked for.

`eltype` is the element type of the tensors and defaults to `ComplexF64`.
Passing `Float64` gives a real state, which makes every later contraction
significantly cheaper, but is only correct when the whole computation stays real.
"""
struct RandomState{R <: PM} end

function RandomState{Pure}(elt::Type{<:Number}, system::System, linkdims::Int)
    st = random_mps(elt, system.pure_indices; linkdims)
    return State{Pure}(system, st)
end

function RandomState{Mixed}(elt::Type{<:Number}, system::System, linkdims::Int)
    n = length(system)
    super = system ⊗ system
    # rounding up rather than down: `floor` would land on the square below, so asking for
    # 50 gave 49 and asking for 10 gave 9
    super_rand = mix(RandomState{Pure}(elt, super, ceil(Int, sqrt(linkdims))))
    return truncate(partial_trace(super_rand, collect(1:n); keepers = true);
                    limits = Limits(maxdim = linkdims))
end

RandomState{R}(system::System, linkdims::Int) where R =
    RandomState{R}(ComplexF64, system, linkdims)

RandomState{R}(elt::Type{<:Number}, size::Int, site::AbstractSite, linkdims::Int) where R =
    RandomState{R}(elt, System(size, site), linkdims)

RandomState{R}(size::Int, site::AbstractSite, linkdims::Int) where R =
    RandomState{R}(ComplexF64, System(size, site), linkdims)

RandomState{R}(elt::Type{<:Number}, sites::Vector{<:AbstractSite}, linkdims::Int) where R =
    RandomState{R}(elt, System(sites), linkdims)

RandomState{R}(sites::Vector{<:AbstractSite}, linkdims::Int) where R =
    RandomState{R}(ComplexF64, System(sites), linkdims)

function RandomState(elt::Type{<:Number}, state::State{Pure}, linkdims::Int)
    st = copy(state.state)
    for i in eachindex(st)
        st[i] = ITensors.convert_eltype(elt, st[i])
    end
    ITensorMPS.randomizeMPS!(elt, st, state.system.pure_indices, linkdims)
    return State{Pure}(state.system, st)
end

RandomState(state::State{Pure}, linkdims::Int) =
    RandomState(ComplexF64, state, linkdims)

RandomState(::State{Mixed}, ::Int) =
    error("randomizing a mixed state is not implemented")

RandomState(::Type{<:Number}, ::State{Mixed}, ::Int) =
    error("randomizing a mixed state is not implemented")
