# Random states on a system: a random pure state of a given link dimension, a state randomised
# in its sector, and a random mixed state drawn as a purification.

export random_state

"""
    random_state(Pure, system, linkdims; rng, eltype = ComplexF64)
    random_state(state, linkdims; rng, eltype = ComplexF64)
    random_state(Mixed, system, linkdims; rng, eltype = ComplexF64)
    random_state(Mixed, system, states, linkdims; rng, eltype = ComplexF64)

a random state of link dimension `linkdims`: a pure state, the pure `state` randomised in its
sector, or a mixed state drawn as a purification on twice the system, of link dimension
`isqrt(linkdims)^2`. On a system conserving something, a pure state is drawn by randomising
one of the sector wanted, and a mixed one from the local `states` its purification starts
from, written as for `State`; a mixed state cannot be drawn on a system conserving something
strongly. `rng` is the random number generator, the global one by default; `eltype = Float64`
gives a real state, cheaper to contract, correct only if the whole computation stays real.

# Examples

    random_state(Pure, system, 20)
    random_state(Mixed, system, 16; eltype = Float64)
    random_state(State{Pure}(System(6, Fermion(conserve = N)), "Occ"), 8)
"""
function random_state(::Type{Pure}, system::System, linkdims::Int;
                      rng::AbstractRNG = Random.default_rng(), eltype::Type{<:Number} = ComplexF64)
    if is_charged(system)
        error("cannot draw a random pure state on a system that conserves something, there " *
              "being no sector to draw it in: randomise a state of the sector wanted, " *
              "random_state(state, linkdims)")
    end
    return State{Pure}(system, random_mps(rng, eltype, siteinds(Pure, system); linkdims))
end

function random_state(state::State{Pure}, linkdims::Int;
                      rng::AbstractRNG = Random.default_rng(), eltype::Type{<:Number} = ComplexF64)
    st = copy(state.state)
    for i in eachindex(st)
        # copied: convert_eltype may hand back the same tensor, which randomizeMPS! may write
        # in place; only the conversion tells whether complex amplitudes are complex indeed
        st[i] = try
            ITensors.convert_eltype(eltype, copy(st[i]))
        catch e
            e isa InexactError || rethrow()
            error("a state of complex amplitudes is not randomised into one of eltype $eltype")
        end
    end
    ITensorMPS.randomizeMPS!(rng, eltype, st, siteinds(Pure, state.system), linkdims)
    return State{Pure}(state.system, st)
end

random_state(::State{Mixed}, ::Int; kwargs...) = error("randomising a mixed state is not implemented")

"""
    refuse_strong(system)

refuse a system conserving something strongly, on which tracing half of a purification would
leave a mixture over several sectors
"""
function refuse_strong(system::System)
    if !isempty(strong_names(system))
        error("cannot draw a random mixed state on a system conserving something strongly: " *
              "tracing half of its purification leaves a mixture over several sectors, " *
              "which such a system cannot hold")
    end
end

function random_state(::Type{Mixed}, system::System, linkdims::Int;
                      rng::AbstractRNG = Random.default_rng(), eltype::Type{<:Number} = ComplexF64)
    refuse_strong(system)
    if is_charged(system)
        error("cannot draw a random mixed state on a system that conserves something without a " *
              "sector to start from: name the local states its purification starts from, " *
              "random_state(Mixed, system, states, linkdims)")
    end
    purification = random_state(Pure, system ⊗ system, isqrt(linkdims); rng, eltype)
    # the ancillas lie on the right, so that tracing them out takes no fermionic sign
    return traced(mix(purification), 1:length(system))
end

"""
    double(states, n)

the local states of a purification on the doubled system, `states` on each half, for a system
of `n` sites
"""
function double(states::Vector, n::Int)
    if length(states) ≠ n
        error("incompatible sizes between system ($n) and states ($(length(states)))")
    end
    return [ states ; states ]
end

# a list of type Any, as in `State`: a vector of numbers is the amplitudes of one local state
double(states, n::Int) = Any[ states for _ in 1:2n ]
double(states::Vector{<:Number}, n::Int) = Any[ states for _ in 1:2n ]

function random_state(::Type{Mixed}, system::System, states, linkdims::Int;
                      rng::AbstractRNG = Random.default_rng(), eltype::Type{<:Number} = ComplexF64)
    refuse_strong(system)
    start = State{Pure}(system ⊗ system, double(states, length(system)))
    return traced(mix(random_state(start, isqrt(linkdims); rng, eltype)), 1:length(system))
end
