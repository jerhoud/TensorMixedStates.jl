# The methods of the solvers, one structure each, shared by the functions and the phases: Krylov,
# the parameters of the local steps, Tdvp and ApproxW for the evolutions, Dmrg for the searches,
# and the abstract types an algorithm of one's own extends.

export Krylov, Tdvp, ApproxW, Dmrg

"""
    Krylov(; dim = nothing, maxiter = nothing, tol = nothing)

the parameters of the Krylov method of each local step: the exponentiation of tdvp, the
eigenvector search of dmrg. A field left `nothing` takes the default of the method, field by
field, that of KrylovKit for tdvp (30, 100, 1e-12) and of ITensorMPS for dmrg (3, 1, 1e-14),
`steady_state` taking `dim = 8` and `maxiter = 3`.

- `dim`: the largest dimension of a Krylov space
- `maxiter`: the number of Krylov spaces built in turn
- `tol`: the tolerance; for tdvp, it sets the number of vectors, `dim` only bounding it

# Examples

    Krylov(tol = 1e-10)            # tdvp: fewer vectors, a lower precision
    Krylov(dim = 8, maxiter = 3)   # dmrg: a more accurate local step
"""
@kwdef struct Krylov
    dim::Union{Nothing, Int} = nothing
    maxiter::Union{Nothing, Int} = nothing
    tol::Union{Nothing, Float64} = nothing
end

show(io::IO, k::Krylov) =
    print(io, "Krylov(", join((string(f, " = ", repr(getfield(k, f))) for f in fieldnames(Krylov)
                               if !isnothing(getfield(k, f))), ", "), ")")

"""
    krylov_kwargs(k, defaults = Krylov(), prefix = "")

the fields of `k` that are given, those of `defaults` in place of the others, as keyword
arguments of KrylovKit prefixed by `prefix`
"""
function krylov_kwargs(k::Krylov, defaults::Krylov = Krylov(), prefix = "")
    value(f) = isnothing(getfield(k, f)) ? getfield(defaults, f) : getfield(k, f)
    return NamedTuple(Symbol(prefix, f == :dim ? :krylovdim : f) => value(f)
                      for f in fieldnames(Krylov) if !isnothing(value(f)))
end

"""
    EvolutionAlgorithm

the supertype of the methods of `evolve` and `thermalize`, `Tdvp` and `ApproxW`. An algorithm of
one's own is a subtype with a method of `TensorMixedStates.evolve_step`, which makes one step,
see `evolve_step`.
"""
abstract type EvolutionAlgorithm end

"""
    SearchAlgorithm

the supertype of the methods of `ground_state` and `steady_state`, `Dmrg`. An algorithm of one's
own is a subtype with a method of `TensorMixedStates.search`, see `search`.
"""
abstract type SearchAlgorithm end

# the fields are read from the type, as for a phase
show(io::IO, a::Union{EvolutionAlgorithm, SearchAlgorithm}) =
    print(io, nameof(typeof(a)), "(",
          join(("$f = $(repr(getfield(a, f)))" for f in fieldnames(typeof(a))), ", "), ")")

"""
    Tdvp(; krylov = Krylov(), expand_period = 0, hermitianize_period = 0)

time evolution by tdvp, for `evolve` and `thermalize`. The bond dimension is enlarged by a global
Krylov expansion before any step that starts from a state whose largest bond dimension is below
`min(5, maxdim)`, as a product state: tdvp is held to the bonds of the state it starts from, and
would leave the error of the first steps, 5e-5 instead of 1e-7 on a Heisenberg chain. A single
site, which tdvp does not take, is exponentiated whole.

- `krylov`: the Krylov exponentiation of each local step, see `Krylov`
- `expand_period`: also enlarge it before the first step and then every `expand_period` steps
  (default 0, never)
- `hermitianize_period`: make a mixed state hermitian every `hermitianize_period` steps
  (default 0, never)

# Examples

    Tdvp()
    Tdvp(expand_period = 5)
    Tdvp(krylov = Krylov(tol = 1e-10))
"""
@kwdef struct Tdvp <: EvolutionAlgorithm
    krylov::Krylov = Krylov()
    expand_period::Int = 0
    hermitianize_period::Int = 0
end

"""
    w_approx_coefs

for each order from 1 to 4, the fractions of the time step whose W approximations, applied in
turn, make up the approximation of that order
"""
const w_approx_coefs = Vector{ComplexF64}[
    [
        1.
    ],
    [
        0.5 + 0.5im,
        0.5 - 0.5im
    ],
    [
        0.10566243270259355 - 0.39433756729740643im,
        0.39433756729740643 + 0.10566243270259355im,
        0.39433756729740643 - 0.10566243270259355im,
        0.10566243270259355 + 0.39433756729740643im
    ],
    [
        0.2588533986109182 + 0.0447561340111419im,
        -0.03154685814880379 + 0.24911905427556322im,
        0.1908290521106672 - 0.23185374923210605im,
        0.16372881485443674,
        0.1908290521106672 + 0.23185374923210605im,
        -0.03154685814880379 - 0.24911905427556322im,
        0.2588533986109182 - 0.0447561340111419im,
    ]
]

"""
    ApproxW(; order, w = 2, apply_algo = "naive", hermitianize_period = 0)

time evolution by the WI or WII approximations of the exponential of the evolver, combined into
an approximation of the given order, for `evolve` and `thermalize`.

- `order`: the order of the approximation, from 1 to 4. In imaginary time, as for
  `thermalize`, two of the seven sub steps of order 4 go backwards, which amplifies the high
  energies by about ``e^{0.03\\,d\\beta\\,\\Delta E}``, `ΔE` the width of the spectrum: unstable
  unless `dβ ΔE` is small, where orders 1 to 3 only go forwards
- `w`: 1 or 2 for WI or WII (default 2)
- `apply_algo`: the algorithm of the product of the state by each MPO, see `apply`: `"naive"`
  (default), exact then truncated, or `"densitymatrix"`, faster at large bond dimensions without
  conservation (2.4x at 128 on a Heisenberg chain) but with the cutoff `eps()` it needs, which
  stops the error near 1e-7 and lets it grow with the number of steps
- `hermitianize_period`: make a mixed state hermitian every `hermitianize_period` steps
  (default 0, never): from order 2, the complex fractions of the step keep a density matrix
  hermitian only up to the order of the approximation

# Examples

    ApproxW(order = 2)
    ApproxW(order = 4, w = 1)
"""
struct ApproxW <: EvolutionAlgorithm
    order::Int
    w::Int
    apply_algo::String
    hermitianize_period::Int
    function ApproxW(order, w, apply_algo, hermitianize_period)
        if !(1 ≤ order ≤ length(w_approx_coefs))
            error("ApproxW has orders 1 to $(length(w_approx_coefs)), not $order")
        end
        if w ∉ (1, 2)
            error("ApproxW has the approximations WI and WII, w = 1 or 2, not $w")
        end
        check_apply_algo(apply_algo)
        return new(order, w, apply_algo, hermitianize_period)
    end
end

ApproxW(; order::Int, w::Int = 2, apply_algo::String = "naive",
        hermitianize_period::Int = 0) = ApproxW(order, w, apply_algo, hermitianize_period)

"""
    Dmrg(; nsweeps, noise = 0, krylov = Krylov(), tol = 0)

the search of an eigenvector by dmrg, for `ground_state` and `steady_state`.

- `nsweeps`: the number of sweeps, at most
- `noise`: the noise of each sweep, a number or one value per sweep, continued with its last
- `krylov`: the Krylov search of each local step, see `Krylov`
- `tol`: the search stops when the energy changes by less than `tol` from one sweep to the next
  (default 0, never)

# Examples

    Dmrg(nsweeps = 10)
    Dmrg(nsweeps = 20, noise = [1e-6, 1e-7, 0], tol = 1e-10)
"""
struct Dmrg <: SearchAlgorithm
    nsweeps::Int
    noise::Union{Float64, Vector{Float64}}
    krylov::Krylov
    tol::Float64
    function Dmrg(nsweeps, noise, krylov, tol)
        if nsweeps < 1
            error("a search takes at least one sweep, and nsweeps is $nsweeps")
        end
        return new(nsweeps, noise isa AbstractVector ? Vector{Float64}(noise) : Float64(noise),
                   krylov, tol)
    end
end

Dmrg(; nsweeps::Int, noise::Union{Real, AbstractVector{<:Real}} = 0., krylov::Krylov = Krylov(),
     tol::Real = 0.) = Dmrg(nsweeps, noise, krylov, tol)

"""
    due(period, k)

whether something done every `period` steps is due at step `k`, never for a period below one
"""
due(period::Int, k::Int) = period > 0 && k % period == 0

"""
    step_value(x, k)
    step_limits(limits, k)

the value of step `k`, from 1, of a parameter given as a number or one value per step, continued
with its last, and the limits of step `k`
"""
step_value(x::AbstractVector, k::Int) = x[min(k, end)]
step_value(x, _) = x

step_limits(l::Limits, k::Int) =
    Limits(step_value(l.cutoff, k), step_value(l.maxdim, k), step_value(l.mindim, k))

"""
    later_steps(x, done)

the values of a parameter given as one value per step or per sweep, for the steps after the
first `done`, a number staying as it is
"""
later_steps(x::AbstractVector, done::Int) = x[min(done + 1, end):end]
later_steps(x, _) = x

later_steps(l::Limits, done::Int) =
    Limits(later_steps(l.cutoff, done), later_steps(l.maxdim, done), later_steps(l.mindim, done))
