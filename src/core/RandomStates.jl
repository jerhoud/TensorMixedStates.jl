# Random states on a system: a random pure state of a given link dimension, a state randomised
# in its sector, and a random mixed state drawn as a purification.

export random_state

"""
    random_state(Pure, system, linkdims; rng, eltype = ComplexF64)
    random_state(Pure, system, states, linkdims; rng, eltype = ComplexF64)
    random_state(state, linkdims; rng, eltype = ComplexF64)
    random_state(Mixed, system, linkdims; rng, eltype = ComplexF64)
    random_state(Mixed, system, states, linkdims; rng, eltype = ComplexF64)
    random_state(Trajectory, system, linkdims; rng, eltype = ComplexF64)
    random_state(Trajectory, system, states, linkdims; rng, eltype = ComplexF64)

a random state of link dimension `linkdims`: a pure state, the pure `state` randomised in its
sector, or a mixed state drawn as a purification on twice the system, of link dimension
`isqrt(linkdims)^2`. The ancillas lie on the right of the system, so that tracing them out
takes no fermionic sign: a mixed state has a rank of at most `isqrt(linkdims)`, and a purity of
at least `1 / isqrt(linkdims)`, whatever the size of the system. On a system conserving
something, a pure state is drawn in the sector of the local `states`, written as for `State`,
or by randomising a state of the sector wanted, and a mixed one from the local `states` its
purification starts from; a mixed state cannot be drawn on a system conserving something
strongly. A trajectory is drawn as the mixed state is, its ancillas then measured, so that its
average over these measurements is the mixed state of that purification. `rng` is the random
number generator, the global one by default; `eltype = Float64` gives a real state, cheaper to
contract, correct only if the whole computation stays real.

# Examples

    random_state(Pure, system, 20)
    random_state(Mixed, system, 16; eltype = Float64)
    random_state(Pure, System(4, Fermion(conserve = N)), ["Occ", "Emp", "Occ", "Emp"], 8)  # N = 2
    random_state(ψ, 8)     # in the sector of ψ
"""
function random_state(::Type{Pure}, system::System, linkdims::Int;
                      rng::AbstractRNG = Random.default_rng(), eltype::Type{<:Number} = ComplexF64)
    if is_charged(system)
        error("cannot draw a random pure state on a system that conserves something, there " *
              "being no sector to draw it in: give the local states of the sector wanted, " *
              "random_state(Pure, system, states, linkdims), or randomise a state of it, " *
              "random_state(state, linkdims)")
    end
    return State{Pure}(system, random_mps(rng, eltype, siteinds(Pure, system); linkdims))
end

random_state(::Type{Pure}, system::System, states, linkdims::Int;
             rng::AbstractRNG = Random.default_rng(), eltype::Type{<:Number} = ComplexF64) =
    random_state(State{Pure}(system, states), linkdims; rng, eltype)

function random_state(state::State{Pure}, linkdims::Int;
                      rng::AbstractRNG = Random.default_rng(), eltype::Type{<:Number} = ComplexF64)
    st = copy(state.state)
    for i in eachindex(st)
        # copied: convert_eltype may hand back the same tensor, which randomizeMPS! may write
        # in place; only the conversion tells whether complex amplitudes are complex indeed
        st[i] = try
            ITensors.convert_eltype(eltype, copy(st[i]))
        catch e
            if !(e isa InexactError)
                rethrow()
            end
            error("a state of complex amplitudes is not randomised into one of eltype $eltype")
        end
    end
    ITensorMPS.randomizeMPS!(rng, eltype, st, siteinds(Pure, state.system), linkdims)
    return State{Pure}(state.system, st)
end

random_state(::State{Mixed}, ::Int; kwargs...) =
    error("randomising a mixed state is not implemented")

random_state(::State{Trajectory}, ::Int; kwargs...) =
    error("randomising a trajectory is not implemented: draw one with " *
          "random_state(Trajectory, system, states, linkdims)")

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

"""
    purification(R, system, linkdims; rng, eltype)
    purification(R, system, states, linkdims; rng, eltype)

the random pure state on twice `system` whose ancillas, on the right, purify the random mixed
state of `random_state(R, …)`, drawn from the local `states` on each half when they are given
"""
function purification(R, system::System, linkdims::Int; rng::AbstractRNG, eltype::Type{<:Number})
    refuse_strong(system)
    if is_charged(system)
        error("cannot draw a random mixed state on a system that conserves something without a " *
              "sector to start from: name the local states its purification starts from, " *
              "random_state($(nameof(R)), system, states, linkdims)")
    end
    return random_state(Pure, system ⊗ system, isqrt(linkdims); rng, eltype)
end

function purification(R, system::System, states, linkdims::Int; rng::AbstractRNG,
                      eltype::Type{<:Number})
    refuse_strong(system)
    start = State{Pure}(system ⊗ system, double(states, length(system)))
    return random_state(start, isqrt(linkdims); rng, eltype)
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

# the ancillas lie on the right, so that tracing them out takes no fermionic sign
random_state(::Type{Mixed}, system::System, linkdims::Int;
             rng::AbstractRNG = Random.default_rng(), eltype::Type{<:Number} = ComplexF64) =
    traced(mix(purification(Mixed, system, linkdims; rng, eltype)), 1:length(system))

random_state(::Type{Mixed}, system::System, states, linkdims::Int;
             rng::AbstractRNG = Random.default_rng(), eltype::Type{<:Number} = ComplexF64) =
    traced(mix(purification(Mixed, system, states, linkdims; rng, eltype)), 1:length(system))

random_state(::Type{Trajectory}, system::System, linkdims::Int;
             rng::AbstractRNG = Random.default_rng(), eltype::Type{<:Number} = ComplexF64) =
    projected(purification(Trajectory, system, linkdims; rng, eltype), system, rng)

random_state(::Type{Trajectory}, system::System, states, linkdims::Int;
             rng::AbstractRNG = Random.default_rng(), eltype::Type{<:Number} = ComplexF64) =
    projected(purification(Trajectory, system, states, linkdims; rng, eltype), system, rng)

no_random_mixture() =
    error("a random mixed state is drawn from the purification of a product state: give its " *
          "local states")

random_state(::Type{Mixed}, ::System, ::Union{Mixture, FullyMixed}, ::Int; kwargs...) =
    no_random_mixture()
random_state(::Type{Trajectory}, ::System, ::Union{Mixture, FullyMixed}, ::Int; kwargs...) =
    no_random_mixture()

"""
    projected(p, system, rng)

the state of `system` left by the purification `p` of a mixed state, its ancillas on the right
measured one after the other in their basis, each outcome drawn with its probability given
those before, see `purification`: the average of its density matrix over the draws is the
mixed state
"""
function projected(p::State{Pure}, system::System, rng::AbstractRNG)
    n = length(system)
    # the system on the left of the centre is orthonormal: the weights of the ancillas are those
    # of the tensors on its right
    st = orthogonalize(p.state, n + 1)
    a = st[n + 1]
    for k in n+1:2n
        s = siteind(Pure, p.system, k)
        ps = [ norm(a * onehot(dag(s) => j))^2 for j in 1:dim(s) ]
        j = draw_index(ps, rng)
        a = a * onehot(dag(s) => j) / sqrt(ps[j] / sum(ps))
        if k < 2n
            a *= st[k + 1]
        end
    end
    ts = [ st[i] for i in 1:n ]
    ts[n] *= a
    mps = MPS(ts)
    if hasqns(ts[1]) && !is_charged(system)
        mps = dense(mps)
    end
    return State{Trajectory}(system, replace_siteinds(mps, siteinds(Pure, system)))
end
