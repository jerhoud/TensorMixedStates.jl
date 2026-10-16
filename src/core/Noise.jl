# The channels and generators of the common noises on sites of any type: resetting,
# depolarizing, dephasing and amplitude damping, each as a gate and as a dissipator.

export reset_gate, reset_dissipator, depolarizing_gate, depolarizing_dissipator,
    dephasing_gate, dephasing_dissipator, amplitude_damping_gate, amplitude_damping_dissipator

"""
    check_probability(what, p)

refuse a probability `p` of the channel `what` outside [0, 1]
"""
function check_probability(what::String, p::Real)
    if !(0 ≤ p ≤ 1)
        error("$what takes a probability between 0 and 1, not $p")
    end
end

"""
    tensor_power(a, n)

the tensor product of `n` copies of the generic operator `a`
"""
tensor_power(a::GenericOp, n::Int) = reduce(⊗, fill(a, n))

"""
    resets(n, state)
    resets(states)

the channel resetting `n` sites to `state` each, or as many sites as states to a state each, a
vector of numbers being the amplitudes of a single state
"""
resets(n::Int, state) = tensor_power(Reset(state), n)
resets(state) = Reset(state)
resets(states::Vector) = reduce(⊗, [ Reset(s) for s in states ])
resets(states::Vector{<:Number}) = Reset(states)

"""
    reset_dissipator(γ, state)
    reset_dissipator(γ, n, state)
    reset_dissipator(γ, [state1, state2, ...])

the generator resetting a site, `n` sites, or as many sites as states, at rate `γ` to `state`,
or to a state per site, the sites being reset together: `γ * Dissipator(Reset(state))`,
``\\gamma\\,(\\sigma \\otimes \\mathrm{tr}_S\\rho - \\rho)``. A state is given as for `Reset`,
and a vector of numbers is the amplitudes of a single state. Evolving under it for a time `t`
is `reset_gate(1 - exp(-γt), ...)`.

# Examples

    -im * H + sum(reset_dissipator(0.1, "Up")(i) for i in 1:10)
"""
reset_dissipator(γ::Real, s...) = γ * Dissipator(resets(s...))

"""
    reset_gate(p, state)
    reset_gate(p, n, state)
    reset_gate(p, [state1, state2, ...])

the channel resetting a site, `n` sites, or as many sites as states, with probability `p` to
`state`, or to a state per site, the sites being reset together:
``(1 - p)\\,\\rho + p\\,\\sigma \\otimes \\mathrm{tr}_S\\rho``.

# Examples

    reset_gate(0.1, "Up")(3)
    reset_gate(0.1, ["Up", "Dn"])(1, 2)
"""
function reset_gate(p::Real, s...)
    check_probability("reset_gate", p)
    r = resets(s...)
    return (1 - p) * Gate(tensor_power(Id, nsites(r))) + p * r
end

"""
    depolarizing_dissipator(γ, n = 1)

the generator depolarizing `n` sites together at rate `γ`, `reset_dissipator(γ, n,
"FullyMixed")`, ``\\gamma\\,(\\mathrm{tr}_S\\rho \\otimes I/d^n - \\rho)``. Evolving under it
for a time `t` is `depolarizing_gate(1 - exp(-γt), n)`.

# Examples

    -im * H + sum(depolarizing_dissipator(0.1)(i) for i in 1:10)
"""
depolarizing_dissipator(γ::Real, n::Int = 1) = reset_dissipator(γ, n, "FullyMixed")

"""
    depolarizing_gate(p, n = 1)

the channel depolarizing `n` sites together with probability `p`, `reset_gate(p, n,
"FullyMixed")`, ``(1 - p)\\,\\rho + p\\,\\mathrm{tr}_S\\rho \\otimes I/d^n``: on a qubit,
`(1 - 3p/4) * Gate(Id) + p/4 * (Gate(X) + Gate(Y) + Gate(Z))`. Two sites depolarized together
are not two sites depolarized each on its own, `depolarizing_gate(p) ⊗ depolarizing_gate(p)`.

# Examples

    prod(depolarizing_gate(0.01)(i) for i in 1:10)
    depolarizing_gate(0.02, 2)(1, 2) * controlled(Z)(1, 2)
"""
function depolarizing_gate(p::Real, n::Int = 1)
    check_probability("depolarizing_gate", p)
    return reset_gate(p, n, "FullyMixed")
end

"""
    dephasing_dissipator(γ, A = Level)

the generator dephasing at rate `γ` in the eigenbasis of `A`, by default the basis of the site,
`γ * Dissipator(Dephase(A))`, ``\\gamma\\,(\\Delta(\\rho) - \\rho)``, `Δ` being `Dephase(A)`.
Every coherence between two eigenspaces decays at the same rate `γ`, where `Dissipator(A)`
damps it at ``(\\lambda_a - \\lambda_b)^2 / 2``, faster for distant eigenvalues: the two agree
when `A` has two eigenvalues, `Dissipator(Z)` being `dephasing_dissipator(2, Z)`. Evolving
under it for a time `t` is `dephasing_gate(1 - exp(-γt), A)`.

# Examples

    -im * H + sum(dephasing_dissipator(0.1)(i) for i in 1:10)
    dephasing_dissipator(0.1, Z ⊗ Z)(1, 2)      # the parity of two qubits
"""
dephasing_dissipator(γ::Real, a::GenericOp{IsOperator} = Level) = γ * Dissipator(Dephase(a))

"""
    dephasing_gate(p, A = Level)

the channel dephasing with probability `p` in the eigenbasis of `A`, by default the basis of
the site, ``(1 - p)\\,\\rho + p\\,\\Delta(\\rho)``, `Δ` being `Dephase(A)`.

# Examples

    prod(dephasing_gate(0.01)(i) for i in 1:10)
    dephasing_gate(0.1, Z ⊗ Z)(1, 2)      # the parity of two qubits
"""
function dephasing_gate(p::Real, a::GenericOp{IsOperator, N} = Level) where N
    check_probability("dephasing_gate", p)
    return (1 - p) * Gate(tensor_power(Id, N)) + p * Dephase(a)
end

"""
    amplitude_damping_dissipator(γ, L = Sp)

the generator of the decay at rate `γ` of the excitations the operator `L` of one site
destroys, `γ * Dissipator(L)`: `Sp` on a `Qubit`, from `"Dn"` to `"Up"`, `C` on a `Fermion`,
`A` on a `Boson`. Evolving under it for a time `t` is `amplitude_damping_gate(1 - exp(-γt), L)`.

# Examples

    -im * H + sum(amplitude_damping_dissipator(0.1)(i) for i in 1:10)
    amplitude_damping_dissipator(0.1, A)(3)     # the loss of photons
"""
amplitude_damping_dissipator(γ::Real, a::SimpleOp = Sp) = γ * Dissipator(a)

"""
    DecayEnd(a)

the function of a site giving the density matrix every state of the site ends in under the
dissipator of `a`, refused when there is none: its generator must have a single steady state,
every other mode decaying
"""
struct DecayEnd <: Function
    jump::SimpleOp
end

function (d::DecayEnd)(site::AbstractSite)
    e = eigen(matrix(Dissipator(d.jump), site))
    tol = spectral_tol(e.values)
    steady = findall(λ -> abs(λ) ≤ tol, e.values)
    if length(steady) ≠ 1 || any(λ -> abs(λ) > tol && real(λ) ≥ -tol, e.values)
        error("$(obs_name(Dissipator(d.jump))) does not take every state of $site to a single " *
              "one, so its amplitude damping has no probability 1")
    end
    σ = reshape(e.vectors[:, only(steady)], dim(site), dim(site))
    return σ / tr(σ)
end

show(io::IO, d::DecayEnd) = show_func(io, "DecayEnd", d.jump)

"""
    amplitude_damping_gate(p, L = Sp)

the channel of `amplitude_damping_dissipator(γ, L)` over a time `t`, `p = 1 - exp(-γt)` being the
probability that an excitation of a `Qubit`, a `Fermion` or a `Boson` decays. On a qubit, its
Kraus operators are `Id - (1 - sqrt(1 - p)) * N` and `sqrt(p) * Sp`. For `p = 1`, an infinite
time, it is the reset to the state the decay ends in, refused when there is none.

# Examples

    prod(amplitude_damping_gate(0.01)(i) for i in 1:10)
    amplitude_damping_gate(0.1, A)(3)           # each photon lost with probability 0.1
"""
function amplitude_damping_gate(p::Real, a::SimpleOp = Sp)
    check_probability("amplitude_damping_gate", p)
    return p < 1 ? exp(-log1p(-p) * Dissipator(a)) : Reset(DecayEnd(a))
end
