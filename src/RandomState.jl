# RandomState, random pure or mixed states on a system with a given link dimension.

export RandomState

"""
    RandomState{Pure|Mixed}([eltype, ]::System, linkdims::Int)
    RandomState{Pure|Mixed}([eltype, ]::Int, ::AbstractSite, linkdims::Int)
    RandomState{Pure|Mixed}([eltype, ]::Vector{<:AbstractSite}, linkdims::Int)
    RandomState([eltype, ]::State, linkdims::Int)
    RandomState{Mixed}([eltype, ]::System, states, linkdims::Int)

a random state of link dimension `linkdims` or, given a pure `State`, that state randomised.

A mixed state is drawn as a purification on twice the system: its link dimension is the
largest square not above `linkdims`, `isqrt(linkdims)^2`, 9 for 10.

`eltype`, `ComplexF64` by default, is the element type of the tensors: `Float64` gives a real
state, significantly cheaper to contract, but correct only if the whole computation stays real.

On a system conserving something, the forms without `states` or a `State` are refused, there
being no sector to draw in. A pure state is drawn by randomising a state you have, which stays
in its sector; a mixed one by naming the `states` its purification starts from, written as for
`State`, the result being a mixture over the sectors around the one named. A site whose charge
the named states pin on its two copies stays pure, no correlation crossing it, so that the link
dimension may be lower than asked. A mixed state cannot be drawn on a system conserving
something strongly, and randomising a mixed state is not implemented.

A system conserving nothing knows nothing of the parity of its fermionic sites, which every
physical state keeps: a pure state drawn on it has no definite parity, and a mixed one holds
coherences between the two. Conserving the parity, as `Fermion(conserve = parity(N))`, or the
number of fermions keeps it: a pure state is then drawn by randomising one of the parity
wanted, a mixed one by naming the states its purification starts from.

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
link dimension `isqrt(linkdims)`, made mixed and traced over its second half, then put on
`system`, whose indices are those of the first half. Nothing is truncated: the spectrum has no
small tail, and a cut would leave a trace off one and negative eigenvalues.
"""
function purify(elt::Type{<:Number}, start::State{Pure}, system::System, linkdims::Int)
    super_rand = mix(RandomState(elt, start, isqrt(linkdims)))
    ρ = partial_trace(super_rand, collect(1:length(system)); keep = true)
    return State{Mixed}(system, ρ.state)
end

"""
    refuse_strong(system)

refuse a system conserving something strongly, on which tracing half of a purification would
leave a mixture over several sectors; checked first, so that such a system is not told to name
the states of its purification.
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
# a list of type Any, as in `State`: a `Vector{Int}` would be taken for the amplitudes of a
# single site
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
        # copied: convert_eltype may hand back the same tensor, which randomizeMPS! may write
        # in place
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
