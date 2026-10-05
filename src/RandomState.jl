# RandomState, random pure or mixed states on a system with a given link dimension.

export RandomState

"""
    RandomState{Pure|Mixed}([eltype, ]::System, linkdims::Int)
    RandomState{Pure|Mixed}([eltype, ]::Int, ::AbstractSite, linkdims::Int)
    RandomState{Pure|Mixed}([eltype, ]::Vector{<:AbstractSite}, linkdims::Int)
    RandomState([eltype, ]::State, linkdims::Int)
    RandomState{Mixed}([eltype, ]::System, states, linkdims::Int)

a random state of link dimension `linkdims` or, given a pure `State`, that state randomised.

A mixed state is built as a purification on twice the system, which squares the link
dimension: its link dimension is the largest square not above `linkdims`,
`isqrt(linkdims)^2`, 9 for 10.

`eltype`, `ComplexF64` by default, is the element type of the tensors. `Float64` gives a real
state, which makes every later contraction significantly cheaper, but is correct only if the
whole computation stays real.

On a system conserving something there is no sector to draw a state in, so the forms without
`states` or a `State` are refused. A pure state is drawn by randomising a state you have, which
stays in its sector; a mixed one by naming the `states` its purification starts from, written
as for `State`, the result being a mixture over the sectors around the one named. The
purification draws its correlations through the sites it can change: a site whose charge the
named states pin, on its two copies, stays pure, and no correlation crosses it, so that the
link dimension may be lower than asked. No mixed
state can be drawn when a quantity is conserved strongly, since tracing half of the
purification leaves a mixture over several sectors. Randomising a mixed state is not
implemented.

# Examples

    RandomState{Pure}(system, 20)
    RandomState{Mixed}(Float64, 10, Qubit(), 16)
    RandomState(state, 20)
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

"""
    purify(elt, start, system, linkdims)

a random mixed state on `system`: the pure state `start`, on `system ⊗ system`, randomised at
link dimension `isqrt(linkdims)`, made mixed, which squares its link dimension, and traced over
its second half.

Nothing is truncated: truncating the density matrix to a size in between would cut into a
spectrum with no small tail, leaving a trace off one and negative eigenvalues. A state a
little below the size asked for is also kept whole by the first truncation of an evolution or
of dmrg, where one above it would be truncated at once. The traced state is put back on
`system`, whose indices are those of the first half, which the partial trace keeps.
"""
function purify(elt::Type{<:Number}, start::State{Pure}, system::System, linkdims::Int)
    super_rand = mix(RandomState(elt, start, isqrt(linkdims)))
    ρ = partial_trace(super_rand, collect(1:length(system)); keep = true)
    return State{Mixed}(system, ρ.state)
end

"""
    refuse_strong(system)

refuse a system conserving something strongly, on which tracing half of a purification would
leave a mixture over several sectors. It is checked first, so that such a system is not told
to name the states of its purification, which would be refused as well.
"""
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

the local states of a purification on the doubled system, `states` on each half, for a
system of `n` sites.
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
