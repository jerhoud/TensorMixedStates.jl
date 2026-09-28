export RandomState

"""
    RandomState{Pure|Mixed}([eltype, ]::System, linkdims::Int)
    RandomState{Pure|Mixed}([eltype, ]::Int, ::AbstractSite, linkdims::Int)
    RandomState{Pure|Mixed}([eltype, ]::Vector{<:AbstractSite}, linkdims::Int)
    RandomState([eltype, ]::State, linkdims::Int)
    RandomState{Mixed}([eltype, ]::System, states, linkdims::Int)

Return a random state, with the specified link dimension.
If a State is given, randomize the given state.

A mixed one is built as a purification on twice the system, whose mixed representation
squares the link dimension: its link dimension is the largest square not above the one asked
for, `isqrt(linkdims)^2`, 9 for 10.

`eltype` is the element type of the tensors and defaults to `ComplexF64`.
Passing `Float64` gives a real state, which makes every later contraction
significantly cheaper, but is only correct when the whole computation stays real.

On a system whose sites conserve something there is no sector to draw a state in, so the
forms taking a system alone do not apply. A pure one is drawn by randomising a state you
already have, which leaves it in the sector it was in; a mixed one is drawn by naming the
`states` its purification starts from, written as for `State`, since what is left after
tracing half of that purification is a mixture over the sectors around the one named.

A quantity conserved strongly leaves no room for a random mixed state: the purification
would have to be traced in half, and what that leaves spreads over several sectors. A pure one
is still drawn by randomising a state you have.
"""
struct RandomState{R <: PM} end

function RandomState{Pure}(elt::Type{<:Number}, system::System, linkdims::Int)
    if is_charged(system)
        error("cannot draw a random pure state on a system that conserves something, there " *
              "being no sector to draw it in: randomise a state you have instead, with " *
              "RandomState(state, linkdims), which keeps the sector that one lives in")
    end
    st = random_mps(elt, system.pure_indices; linkdims)
    return State{Pure}(system, st)
end

# the purification lives on twice the system, whose mixed representation squares the link
# dimension, so the pure state it starts from is drawn at the square root of what was asked,
# rounded down, and nothing is truncated. Truncating the density matrix to a size in between
# cut into a spectrum with no small tail, which left a trace off one and negative eigenvalues,
# and a state a little below the size asked for is what the first step of an evolution or of
# dmrg, truncating to that size, keeps whole, where one above it gets truncated at once.
# The traced state is put back on `system`: the first half of `system ⊗ system` keeps its
# indices, and those are what the partial trace leaves
function purify(elt::Type{<:Number}, start::State{Pure}, system::System, linkdims::Int)
    super_rand = mix(RandomState(elt, start, isqrt(linkdims)))
    ρ = partial_trace(super_rand, collect(1:length(system)); keepers = true)
    return State{Mixed}(system, ρ.state)
end

# refused first, so that a system conserving something strongly is not sent to the form taking
# states, which refuses it as well
function refuse_strong(system::System)
    if !isempty(strong_names(system))
        error("cannot draw a random mixed state on a system conserving something strongly: " *
              "tracing half of its purification leaves a mixture over several sectors, " *
              "which such a system cannot hold")
    end
end

function RandomState{Mixed}(elt::Type{<:Number}, system::System, linkdims::Int)
    refuse_strong(system)
    if is_charged(system)
        error("cannot draw a random mixed state on a system that conserves something " *
              "without a sector to start from: name the states its purification starts " *
              "from, with RandomState{Mixed}(system, states, linkdims)")
    end
    super = system ⊗ system
    return purify(elt, RandomState{Pure}(elt, super, isqrt(linkdims)), system, linkdims)
end

"""
    double(states, n)

the local states of a purification, the same on each half of the doubled system
"""
double(states::Vector, ::Int) = [ states ; states ]
# one local state for every site, in a list of type Any as `State` builds it, so that an index or
# an amplitude vector repeated is not itself taken for the amplitudes of a single site
double(states, n::Int) = Any[ states for _ in 1:2n ]
double(states::Vector{<:Number}, n::Int) = Any[ states for _ in 1:2n ]

function RandomState{Mixed}(elt::Type{<:Number}, system::System, states, linkdims::Int)
    refuse_strong(system)
    n = length(system)
    super = system ⊗ system
    return purify(elt, State{Pure}(super, double(states, n)), system, linkdims)
end

RandomState{Mixed}(system::System, states, linkdims::Int) =
    RandomState{Mixed}(ComplexF64, system, states, linkdims)

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
        # a tensor of its own: convert_eltype hands back the same one when the type already
        # matches, and randomizeMPS! writes into a single site in place
        st[i] = ITensors.convert_eltype(elt, copy(st[i]))
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
