# The partial trace of a density matrix, its fermionic signs included: the state of some of its
# sites, on the system they make.

export partial_trace

"""
    check_positions(state, positions, what)

refuse a position of `positions` the state does not have, `what` naming the function given it
"""
function check_positions(state::AbstractState, positions, what)
    n = length(state)
    for p in positions
        if !(1 ≤ p ≤ n)
            error("$what was given site $p, which the state does not have: it has $n sites")
        end
    end
    return nothing
end

"""
    sign_blocks(R, system, i, f)

the matrix, on the index of site `i` of `system` in the representation `R`, of the parity `f`
of the site taken on both sides, ``\\rho \\mapsto f\\rho f``, and the function giving that of a
projector `p` taken on the kets alone, ``\\rho \\mapsto p\\rho``: `f` and `p` themselves on a
pure state
"""
sign_blocks(::Type{Pure}, ::System, ::Int, f::AbstractMatrix) = (f, p -> p)

function sign_blocks(::Type{Mixed}, system::System, i::Int, f::AbstractMatrix)
    o = mixed_elements(system, i)
    return (kron(f, f)[o, o], p -> kron(one(p), p)[o, o])
end

"""
    trace_signs(state, kept)

the state whose partial trace keeping the increasing positions `kept` is the reduced state of
`state`, its fermionic signs included: ``U\\rho U^\\dagger``, or ``U|\\psi\\rangle`` on a pure
state, ``U = (-1)^{\\sum n_l n_k}``, `l` running over the fermionic sites traced out, `k` over
the fermionic sites kept on their right, and `n` the parity of a site. It is the product by an
MPO whose links, of dimension 2, carry the parity traced out on the left, read on the ket side
alone of a site traced out, which the trace makes that of the bra. A state needing no sign is
given back as it is.
"""
function trace_signs(state::State{R}, kept) where R
    system = state.system
    n = length(system)
    fermionic(i) = next_fermion(system, i) == i
    last_kept = findlast(fermionic, kept)
    if isnothing(last_kept) ||
       !any(l -> fermionic(l) && !insorted(l, kept), 1:kept[last_kept]-1)
        return state
    end
    links = [ is_charged(system) ? Index(QN() => 2; tags = "Parity,l=$i") :
                                   Index(2; tags = "Parity,l=$i") for i in 0:n ]
    st = state.state
    ts = map(1:n) do i
        x = siteind(R, system, i)
        d = dim(system[i])
        one = Matrix{Float64}(I, d, d)
        f = matrix(F, system[i])
        both, on_kets = sign_blocks(R, system, i, f)
        id = Matrix{Float64}(I, dim(x), dim(x))
        # the maps taking the parity p carried in to the parity carried out
        a = zeros(dim(x), dim(x), 2, 2)
        for p in 1:2
            if !fermionic(i)
                a[:, :, p, p] = id
            elseif insorted(i, kept)
                a[:, :, p, p] = p == 1 ? id : both
            else
                for q in 0:1
                    a[:, :, p, xor(p - 1, q) + 1] = on_kets((one + (-1)^q * f) / 2)
                end
            end
        end
        u = charged_itensor(a, [x', dag(x), dag(links[i]), links[i + 1]])
        return noprime(u * st[i])
    end
    # no parity carried in on the left, either one carried out on the right
    ts[1] *= onehot(links[1] => 1)
    ts[n] *= onehot(dag(links[n + 1]) => 1) + onehot(dag(links[n + 1]) => 2)
    for i in 1:n-1
        l = commonind(st[i], st[i + 1])
        c = combiner(l, links[i + 1]; tags = tags(l))
        ts[i] *= c
        ts[i + 1] *= dag(c)
    end
    return State{R}(system, MPS(ts))
end

"""
    partial_trace(state, positions [; keep = false])

the mixed `state` with the sites at `positions` traced out, or, with `keep = true`, all the
others: a mixed state on the system of the sites kept, in their order, with the trace of
`state`. An operator of the sites kept has on it the value it has on `state`, fermionic signs
included; tracing out a fermionic site with fermionic sites kept on its right may double the
bond dimension.

A pure state is refused (`mix` it first), so is a state conserving something strongly
(`weaken` it first), and so is tracing out every site.

# Examples

    partial_trace(ρ, [1, 2])
    partial_trace(ρ, 3:5; keep = true)
"""
function partial_trace(state::State{Mixed}, positions::AbstractVector{Int}; keep::Bool = false)
    if !isempty(strong_names(state.system))
        error("partial_trace cannot trace out part of a state conserving something strongly, " *
              "which leaves a mixture over several sectors: weaken it first")
    end
    check_positions(state, positions, "partial_trace")
    kept = keep ? sort(unique(positions)) : setdiff(1:length(state), positions)
    if isempty(kept)
        error("partial_trace cannot trace out every site of a state")
    end
    return traced(trace_signs(state, kept), kept)
end

# not mixed behind the caller's back, the result being a state
partial_trace(::State{Pure}, ::AbstractVector{Int}; keep::Bool = false) =
    error("partial_trace needs a mixed representation: use mix(state) first")
