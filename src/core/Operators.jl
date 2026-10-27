# The algebra of operators, channels and generators: their types, generic or placed on sites,
# and how they are built, combined, compared and printed, before any site is involved. The later
# layers walk an expression with `scalarcoef`, `scalararg`, `sumsubs`, `prodsubs`,
# `factor_sites`, `scale`, `power`, `is_pure` and `obs_name`.

export Pure, Mixed
export Operator, Id, Proj, Gate, Reset, Dephase, Dissipator, Evolver, named, parity, dag, ⊗,
    tensor, map_sites

############# Types ################

"""
    abstract type Representation

the supertype of the representations of a state: `Pure`, `Mixed`, and those an extension
defines.
"""
abstract type Representation end

"""
    abstract type PM <: Representation

the supertype of `Pure` and `Mixed`, the representations of a state this package defines
"""
abstract type PM <: Representation end

"""
    Pure

the representation of a state as a wave function, as in `State{Pure}`.
"""
struct Pure <: PM end

"""
    Mixed

the representation of a state as a density matrix, as in `State{Mixed}`.
"""
struct Mixed <: PM end

"""
    abstract type GI

the supertype of `Generic`, for an operator not yet placed on sites, and `Indexed`, for one
placed on sites
"""
abstract type GI end
struct Generic <: GI end
struct Indexed <: GI end

"""
    abstract type OpKind
    IsOperator
    IsChannel
    IsGenerator

the kinds of operators: `IsOperator` for an operator, a linear map of the states of its sites;
`IsChannel` for a channel, a completely positive map of density matrices; `IsGenerator` for a
generator, the Lindblad generator of an evolution
"""
abstract type OpKind end
struct IsOperator <: OpKind end
struct IsChannel <: OpKind end
struct IsGenerator <: OpKind end

"""
    Op{K <: OpKind, T <: GI, N}

the supertype of all operators: `K` their kind, `IsOperator`, `IsChannel` or `IsGenerator`, `T`
`Generic` or `Indexed`, `N` the number of sites of a generic operator, 1 for a placed one.
"""
abstract type Op{K <: OpKind, T <: GI, N} end

"""
    GenericOp{K, N}

the operators not yet placed on sites, `Op{K, Generic, N}`.
"""
const GenericOp{K, N} = Op{K, Generic, N}

"""
    IndexedOp{K}

the operators placed on sites, `Op{K, Indexed, 1}`.
"""
const IndexedOp{K} = Op{K, Indexed, 1}

"""
    SimpleOp

the generic operators of one site, `GenericOp{IsOperator, 1}`.
"""
const SimpleOp = GenericOp{IsOperator, 1}

"""
    no_signed_zero(x)

`x` with its signed zeros made positive: `-0.0 == 0.0`, but `isequal` and `hash` tell them
apart
"""
no_signed_zero(x::Union{AbstractFloat, Complex{<:AbstractFloat}}) = x + zero(x)
no_signed_zero(x::AbstractArray) = map(no_signed_zero, x)
no_signed_zero(x) = x

"""
    obs_name(op)

the compact printed form of the operator, a long sum, product, tensor product or list of sites
cut down to its first and its last elements: the default name of a measurement or of a
conserved quantity, and how an error message names an operator
"""
obs_name(op) = sprint(print, op; context = :compact => true)

############### Operator ###############

"""
    Operator{N}(name, def = nothing)

a named operator of `N` sites. `def` defines it: a matrix, a function of the sites, an
expression of other operators, or `nothing` for an operator each site type declares with
`@def_operators`. To define an operator of one's own, see `named`.

# Examples

    Operator{1}("Z", [1 0 ; 0 -1])
    Operator{2}("Swap", [1 0 0 0 ; 0 0 1 0 ; 0 1 0 0 ; 0 0 0 1])
"""
struct Operator{N} <: GenericOp{IsOperator, N}
    name::String
    expr::Union{Nothing, Matrix, Function, GenericOp{IsOperator, N}}
    Operator{N}(name::String, expr::Union{Nothing, AbstractMatrix, Function,
                                          GenericOp{IsOperator, N}} = nothing) where N =
        new{N}(name, no_signed_zero(expr isa AbstractMatrix ? Matrix(expr) : expr))
end

############ Identity ############

"""
    IdentityOp{K, T, N}

the identity, one value for each kind of operator, so that every way of writing it, as
`Id ⊗ Id` or `Id(3)`, is recognized by its type. Placed, it has no site.
"""
struct IdentityOp{K, T, N} <: Op{K, T, N} end

IdentityOp(::Op{K, T, N}) where {K, T, N} = IdentityOp{K, T, N}()

"""
    Id

the identity operator, defined on every site type.
"""
const Id = IdentityOp{IsOperator, Generic, 1}()

################ Numbers and sums #############

"""
    float_coef(x)

the number `x` as the coefficient of an operator: an integer of fixed width, or a complex
number of them, as a float, whose products cannot silently wrap around
"""
float_coef(x::Union{Bool, Base.BitInteger, Complex{<:Union{Bool, Base.BitInteger}}}) = float(x)
float_coef(x::Number) = x

"""
    check_coef(K, c)

the number `c` as a coefficient of an operator of kind `K`: any number for an operator, a real
one for a generator, a real one of zero or above for a channel
"""
check_coef(::Type{IsOperator}, c::Number) = c

check_coef(::Type{IsGenerator}, c::Number) =
    if isreal(c)
        real(c)
    else
        error("a generator takes real coefficients, not $c: write -im * H + Dissipator(L), H an " *
              "operator")
    end

check_coef(::Type{IsChannel}, c::Number) =
    if isreal(c) && real(c) ≥ 0
        real(c)
    else
        error("a channel takes weights of zero or above, not $c")
    end

"""
    ScalarOp{K, T, N}

a number times an operator, never itself a `ScalarOp`: `2 * (3X)` is `6X`. It is built by
`scale`.
"""
struct ScalarOp{K, T, N} <: Op{K, T, N}
    coef::Number
    arg::Op{K, T, N}
end

"""
    scalarcoef(a)

the coefficient of an operator, 1 unless it is a `ScalarOp`
"""
scalarcoef(a::ScalarOp) = a.coef
scalarcoef(::Op) = 1

"""
    scalararg(a)

the operator without its coefficient
"""
scalararg(a::ScalarOp) = a.arg
scalararg(a::Op) = a

"""
    SumOp{K, T, N}

a sum of operators, nested sums flattened and terms of coefficient zero left out. An empty
sum is `0 * Id`, a sum of one term that term.
"""
struct SumOp{K, T, N} <: Op{K, T, N}
    subs::Vector{<:Op{K, T, N}}
    function SumOp(subs::Vector{<:Op{K, T, N}}) where {K, T, N}
        # in one pass: a sum built term by term, as `sum` of a generator does, copies its terms
        # at each one
        s = Op{K, T, N}[]
        for x in subs
            if x isa SumOp
                append!(s, x.subs)
            elseif scalarcoef(x) ≠ 0
                push!(s, x)
            end
        end
        if isempty(s)
            return scale(0., IdentityOp{K, T, N}())
        elseif length(s) == 1
            return s[1]
        else
            return new{K, T, N}(s)
        end
    end
end

"""
    scale(c, a)

the operator `a` times the number `c`: `a` itself for 1, `0 * Id` for 0, each term multiplied
for a sum, `c` checked by `check_coef`
"""
function scale(c::Number, a::Op{K, T, N}) where {K, T, N}
    c = check_coef(K, c)
    return if c == 0
        ScalarOp{K, T, N}(0., IdentityOp{K, T, N}())
    elseif c == 1
        a
    elseif a isa SumOp
        SumOp(map(x -> scale(c, x), a.subs))
    else
        ScalarOp{K, T, N}(no_signed_zero(float_coef(c * scalarcoef(a))), scalararg(a))
    end
end

"""
    sumsubs(a)

the terms of a sum, an operator that is not a sum being its only term
"""
sumsubs(a::SumOp) = a.subs
sumsubs(a::Op) = [a]

############### Products ###############

"""
    ProdOp{K, T, N}

a product of operators, nested products flattened and the coefficients of the factors
gathered in front of it. An empty product is `Id`, a product of one factor that factor.
"""
struct ProdOp{K, T, N} <: Op{K, T, N}
    subs::Vector{<:Op{K, T, N}}
    # the factors as they are given, which ProdOp(subs) below flattens first
    ProdOp{K, T, N}(subs::Vector{<:Op{K, T, N}}) where {K, T, N} = new{K, T, N}(subs)
end

"""
    prodsubs(a)

the factors of a product, without their coefficients, an operator that is not a product being
its only factor
"""
prodsubs(a::ProdOp) = a.subs
prodsubs(a::ScalarOp) = prodsubs(a.arg)
prodsubs(a::Op) = [a]

function ProdOp(subs::Vector{<:Op{K, T, N}}) where {K, T, N}
    c = prod(scalarcoef.(subs); init = 1)
    if isempty(subs)
        return scale(c, IdentityOp{K, T, N}())
    end
    s = reduce(vcat, prodsubs.(subs))
    return scale(c, length(s) == 1 ? s[1] : ProdOp{K, T, N}(s))
end

############### Tensor products ############

"""
    TensorOp{K, N}

a tensor product of generic operators acting on `N` sites, the sum of theirs, the
coefficients of the factors gathered in front of it
"""
struct TensorOp{K, N} <: GenericOp{K, N}
    subs::Vector{<:GenericOp{K}}
    TensorOp{K, N}(subs::Vector{<:GenericOp{K}}) where {K, N} =
        if length(subs) == 1
            subs[1]
        else
            c = prod(scalarcoef.(subs))
            s = scalararg.(subs)
            if all(x -> x isa IdentityOp, s)
                return scale(c, IdentityOp{K, Generic, N}())
            end
            scale(c, new{K, N}(s))
        end
end

"""
    tensorsubs(a)

the factors of a tensor product, an operator that is not one being its only factor
"""
tensorsubs(a::TensorOp) = a.subs
tensorsubs(a::GenericOp) = [a]

"""
    nsites(a)

the number of sites the generic operator `a` acts on
"""
nsites(::GenericOp{K, N}) where {K, N} = N

"""
    factor_sites(a::TensorOp)

the ranges of the sites each factor of the tensor product `a` acts on, in order
"""
function factor_sites(a::TensorOp)
    stop = 0
    return map(a.subs) do o
        start = stop + 1
        stop += nsites(o)
        start:stop
    end
end

############ Proj ################

"""
    Proj(state)
    Proj(A => λ)

the projector ``|s\\rangle\\langle s|`` on a local pure state of one site, given as `State`
takes it, or the projector on the eigenspace of eigenvalue `λ` of the Hermitian operator `A`,
of one site or several. An operator placed on sites places the projector there:
`Proj(N(3) => 1)` is `Proj(N => 1)(3)`.

# Examples

    Proj("Up")
    Proj(Level => 1)       # the second basis state
    Proj(Sz => 0)          # on a Spin(1), the state of Sz = 0
    Proj(Swap => -1)       # the antisymmetric part of two sites, the singlet of two qubits
    Proj(Swap(1, 3) => -1) # that of sites 1 and 3, Proj(Swap => -1)(1, 3)
"""
struct Proj{N} <: GenericOp{IsOperator, N}
    state::Union{String, Vector, Function, Pair{<:GenericOp{IsOperator}, <:Real}}
    Proj(state::Union{String, AbstractVector, Function}) =
        new{1}(no_signed_zero(state isa AbstractVector ? Vector(state) : state))
    Proj(p::Pair{<:GenericOp{IsOperator, N}, <:Real}) where N =
        new{N}(first(p) => no_signed_zero(last(p)))
end

Proj(n::Integer) =
    error("a number does not denote a state: the projector on the basis state $n, counted from " *
          "0, is Proj(Level => $n)")

Proj(p::Pair{<:IndexedOp{IsOperator}, <:Real}) =
    if first(p) isa AtIndex
        Proj(first(p).op => last(p))(first(p).index...)
    else
        error("cannot project on an eigenspace of $(obs_name(first(p))) once placed: take the " *
              "projector of the operator before placing it, Proj((N ⊗ Id + Id ⊗ N) => 1)(1, 2) " *
              "rather than Proj(N(1) + N(2) => 1)")
    end

############ AtIndex ################

"""
    AtIndex{K, N}

an operator placed on sites, as `X(1)` or `Swap(2, 4)`, on distinct sites
"""
struct AtIndex{K, N} <: IndexedOp{K}
    op::GenericOp{K, N}
    index::NTuple{N, Int}
    function AtIndex(op::GenericOp{K, N}, index::NTuple{N, Int}) where {K, N}
        if !allunique(index)
            error("$(obs_name(op)) acts on $N sites and cannot be placed on $index, which " *
                  "repeats a site")
        end
        if scalararg(op) isa IdentityOp
            return scale(scalarcoef(op), IdentityOp{K, Indexed, 1}())
        end
        return scale(scalarcoef(op), new{K, N}(scalararg(op), index))
    end
end

(op::GenericOp{K, N})(index::Vararg{Int, N}) where {K, N} = AtIndex(op, index)

############## Channels and generators ###########

"""
    Gate(A)

the channel ``\\rho \\mapsto A\\rho A^\\dagger`` of the operator `A`, generic or placed on
sites: `Gate(X(1))` is `Gate(X)(1)`, `Gate(c * A)` is `abs2(c) * Gate(A)`, and the gate of a
product is the product of the gates. Channels combine with weights of zero or above, which is
how a noisy gate is written. An operator multiplied by a channel, or in a tensor product with
one, becomes its gate.

# Examples

    0.9 * Gate(Id) + 0.1 * Gate(X)      # X with probability 0.1
    Gate(X(1) + X(2))
"""
struct Gate{T, N} <: Op{IsChannel, T, N}
    arg::Op{IsOperator, T, N}
    function Gate(a::Op{IsOperator, T, N}) where {T, N}
        c, x = abs2(scalarcoef(a)), scalararg(a)
        if x isa IdentityOp
            return scale(c, IdentityOp{IsChannel, T, N}())
        elseif x isa AtIndex
            return scale(c, AtIndex(Gate(x.op), x.index))
        elseif x isa ProdOp
            return scale(c, ProdOp(Op{IsChannel, T, N}[ Gate(y) for y in x.subs ]))
        end
        return scale(c, new{T, N}(x))
    end
end

"""
    Reset(state)

the channel resetting a site to the local state `state`, given as `State` takes it,
``\\rho \\mapsto \\sigma \\otimes \\mathrm{tr}_i \\rho``, ``\\sigma`` the density matrix
of `state` and ``\\mathrm{tr}_i`` the trace over the site. A vector or a matrix is taken as it
is, without normalization. On a fermionic site this holds on states of definite parity: on the
others, an odd operator `O` of the other sites takes the value of
``\\langle F\\rangle_\\sigma F_i O``, ``F_i`` the parity of the site. An operator of one site
placed on it places the reset there: `Reset(N(3) => 1)` is `Reset(N => 1)(3)`.

# Examples

    Reset("Up")(3)
    Reset(Level => 0)(3)
    Reset([0.3 0. ; 0. 0.7])(3)
"""
struct Reset <: GenericOp{IsChannel, 1}
    state::Union{String, Vector, Matrix, Function, Pair{<:SimpleOp, <:Real}}
    Reset(state::Union{String, AbstractVector, AbstractMatrix, Function}) =
        new(no_signed_zero(state isa AbstractArray ? Array(state) : state))
    Reset(p::Pair{<:SimpleOp, <:Real}) = new(first(p) => no_signed_zero(last(p)))
end

Reset(n::Integer) =
    error("a number does not denote a state: the reset to the basis state $n, counted from 0, " *
          "is Reset(Level => $n)")

Reset(p::Pair{<:IndexedOp{IsOperator}, <:Real}) =
    if first(p) isa AtIndex{IsOperator, 1}
        Reset(first(p).op => last(p))(only(first(p).index))
    else
        error("cannot reset to an eigenstate of $(obs_name(first(p))), which is not an operator " *
              "of one site placed on it: write Reset(A => λ)(i), A an operator of one site")
    end

"""
    Dephase(A)

the channel dephasing in the eigenbasis of the Hermitian operator `A`, of one site or several,
``\\rho \\mapsto \\sum_\\lambda P_\\lambda \\rho P_\\lambda``, ``P_\\lambda`` the projector
on the eigenspace of `A` of eigenvalue ``\\lambda``: the measurement of `A` whose result is not
read. It erases the coherences between eigenspaces and keeps those within one. An operator
placed on sites places the dephasing there: `Dephase(N(1))` is `Dephase(N)(1)`.

# Examples

    Dephase(Level)(3)       # in the basis of the site
    Dephase(Ntot)(2)        # on an Electron, keeps the coherence of Up and Dn
    Dephase(Z ⊗ Z)(1, 2)    # the parity of two qubits, not read
"""
struct Dephase{N} <: GenericOp{IsChannel, N}
    arg::GenericOp{IsOperator, N}
    Dephase(arg::GenericOp{IsOperator, N}) where N = new{N}(arg)
end

Dephase(a::IndexedOp{IsOperator}) =
    if a isa AtIndex
        Dephase(a.op)(a.index...)
    else
        error("cannot dephase in the eigenbasis of $(obs_name(a)) once placed: dephase in " *
              "that of the operator before placing it, Dephase(N ⊗ Id + Id ⊗ N)(1, 2) rather " *
              "than Dephase(N(1) + N(2))")
    end

"""
    Dissipator(L)
    Dissipator(Φ)

the Lindblad dissipator of the jump operator `L`, generic or placed on sites,
``\\rho \\mapsto L\\rho L^\\dagger - \\frac{1}{2}\\{L^\\dagger L, \\rho\\}``, or that of
the channel `Φ`, the sum of those of its Kraus operators, ``\\Phi - \\mathrm{id}`` for `Φ`
preserving the trace. `Dissipator(c * L)` is `abs2(c) * Dissipator(L)`, and
`Dissipator(Sm(1))` is `Dissipator(Sm)(1)`.

# Examples

    Dissipator(Sm)(1)
    Dissipator(Sm(1) + Sm(2))           # a collective decay
    0.1 * Dissipator(Reset("Up"))(3)
"""
struct Dissipator{T, N} <: Op{IsGenerator, T, N}
    arg::Op{<:Union{IsOperator, IsChannel}, T, N}
    function Dissipator(a::Op{<:Union{IsOperator, IsChannel}, T, N}) where {T, N}
        x = scalararg(a)
        c = a isa Op{IsOperator} ? abs2(scalarcoef(a)) : scalarcoef(a)
        if x isa IdentityOp
            return scale(0, IdentityOp{IsGenerator, T, N}())
        elseif x isa Gate
            return scale(c, Dissipator(x.arg))
        elseif x isa AtIndex
            return scale(c, AtIndex(Dissipator(x.op), x.index))
        end
        return scale(c, new{T, N}(x))
    end
end

"""
    Evolver(K)

the generator ``\\rho \\mapsto K\\rho + \\rho K^\\dagger`` of the operator `K`, generic or
placed on sites, ``d\\psi/dt = K\\psi`` on a pure state: for `K = -im * H`, the Hamiltonian
part ``-i[H, \\rho]`` of a generator. An operator added to a generator becomes its
`Evolver`. A real coefficient comes out of it, a complex one stays inside.

# Examples

    -im * H + Dissipator(Sm)(1)         # Evolver(-im * H) + Dissipator(Sm)(1)
"""
struct Evolver{T, N} <: Op{IsGenerator, T, N}
    arg::Op{IsOperator, T, N}
    function Evolver(a::Op{IsOperator, T, N}) where {T, N}
        c, x = scalarcoef(a), scalararg(a)
        if x isa AtIndex
            return AtIndex(Evolver(scale(c, x.op)), x.index)
        elseif isreal(c)
            return scale(real(c), new{T, N}(x))
        end
        return new{T, N}(a)
    end
end

"""
    map_sites(f, op)

the operator placed on sites `op` with each factor moved from its sites to their images by
`f`: `map_sites(i -> 2i - 1, X(1) * Y(2))` is `X(1) * Y(3)`.
"""
map_sites(f, a::AtIndex) = AtIndex(a.op, map(f, a.index))
map_sites(f, a::SumOp) = SumOp(map(x -> map_sites(f, x), a.subs))
map_sites(f, a::ProdOp) = ProdOp(map(x -> map_sites(f, x), a.subs))
map_sites(f, a::ScalarOp) = scale(a.coef, map_sites(f, a.arg))
map_sites(_, a::IdentityOp) = a
map_sites(f, a::Gate{Indexed}) = Gate(map_sites(f, a.arg))
map_sites(f, a::Dissipator{Indexed}) = Dissipator(map_sites(f, a.arg))
map_sites(f, a::Evolver{Indexed}) = Evolver(map_sites(f, a.arg))

############## Arithmetic ###########

(a::Number * b::Op) = scale(a, b)
(a::Op * b::Number) = scale(b, a)
(a::Op / b::Number) = scale(inv(b), a)
-(a::Op) = scale(-1, a)
(a::Op{K, T, N} + b::Op{K, T, N}) where {K, T, N} = SumOp([a, b])
(a::Op - b::Op) = a + (-b)
(a::Op{K, T, N} * b::Op{K, T, N}) where {K, T, N} = ProdOp([a, b])

(a::Op{IsGenerator, T, N} * b::Op{IsGenerator, T, N}) where {T, N} =
    error("generators are added, not multiplied: $(obs_name(a)) * $(obs_name(b)) is not a " *
          "generator")

# an operator beside a generator is its Hamiltonian part, beside a channel its gate
(a::Op{IsOperator, T, N} + b::Op{IsGenerator, T, N}) where {T, N} = Evolver(a) + b
(a::Op{IsGenerator, T, N} + b::Op{IsOperator, T, N}) where {T, N} = a + Evolver(b)
(a::Op{IsOperator, T, N} * b::Op{IsChannel, T, N}) where {T, N} = Gate(a) * b
(a::Op{IsChannel, T, N} * b::Op{IsOperator, T, N}) where {T, N} = a * Gate(b)

"""
    op1 ⊗ op2
    tensor(op1, op2, ...)

the tensor product of generic operators, or of channels, acting on the sites of `op1`
followed by those of `op2`: `(A ⊗ B)(i, j)` is `A(i) * B(j)`. `⊗` is typed `\\otimes`, and
`tensor` is the same function under a name in ASCII.

# Examples

    Rxy(t) = exp(-im * t * (X ⊗ X + Y ⊗ Y) / 4)
"""
(a::GenericOp{K, N} ⊗ b::GenericOp{K, M}) where {K, N, M} =
    TensorOp{K, N + M}([tensorsubs(a) ; tensorsubs(b)])
(a::GenericOp{IsGenerator} ⊗ b::GenericOp{IsGenerator}) =
    error("generators are added, not multiplied: $(obs_name(a)) ⊗ $(obs_name(b)), placed as a " *
          "product, is not a generator")
(a::GenericOp{IsOperator} ⊗ b::GenericOp{IsChannel}) = Gate(a) ⊗ b
(a::GenericOp{IsChannel} ⊗ b::GenericOp{IsOperator}) = a ⊗ Gate(b)
⊗(a::GenericOp, b::GenericOp, c::GenericOp, d::GenericOp...) = foldl(⊗, (a, b, c, d...))
const tensor = ⊗

############## Functions of an operator ###########

"""
    is_natural(p)

whether the exponent `p` is an integer of zero or above, whatever its type
"""
is_natural(p::Number) = isreal(p) && isinteger(real(p)) && real(p) ≥ 0

"""
    IntPowOp{K, N}
    GenPowOp{K, N}

the powers of an operator: `IntPowOp` for an integer exponent of zero or above, the product
kept whole; `GenPowOp` for any other, the principal power, from which only the modulus of a
coefficient comes out
"""
struct IntPowOp{K, N} <: GenericOp{K, N}
    arg::GenericOp{K, N}
    expo::Int
    function IntPowOp(arg::GenericOp{K, N}, n::Int) where {K, N}
        c = scalarcoef(arg)
        a = scalararg(arg)
        if n == 0 || a isa IdentityOp
            return scale(c^n, IdentityOp{K, Generic, N}())
        elseif n == 1
            return arg
        elseif a isa IntPowOp
            return scale(c^n, IntPowOp(a.arg, a.expo * n))
        end
        return scale(c^n, new{K, N}(a, n))
    end
end

struct GenPowOp{K, N} <: GenericOp{K, N}
    arg::GenericOp{K, N}
    expo::Number
    function GenPowOp(arg::GenericOp{K, N}, p::Number) where {K, N}
        c = scalarcoef(arg)
        a = scalararg(arg)
        m = abs(c)
        if a isa IdentityOp
            if iszero(c) && real(p) ≤ 0
                error("$(obs_name(new{K, N}(arg, p))) does not exist: $(obs_name(arg)) is not " *
                      "invertible")
            end
            return scale(complex(c)^p, a)
        elseif !iszero(m) && !isone(m)
            # the modulus commutes with everything and leaves the logarithm exactly, the
            # phase would turn the eigenvalues across its cut and stays inside
            return scale(m^p, GenPowOp(scale(c / m, a), p))
        end
        return new{K, N}(arg, no_signed_zero(p))
    end
end

"""
    power(a, p)

`a^p`, an `IntPowOp` when `p` is natural, a `GenPowOp` otherwise
"""
power(a::GenericOp, p::Number) = is_natural(p) ? IntPowOp(a, Int(real(p))) : GenPowOp(a, p)

power(a::GenericOp{IsChannel}, p::Number) =
    if is_natural(p)
        IntPowOp(a, Int(real(p)))
    else
        error("$(obs_name(a)) has no power $p: only the integer powers of a channel are channels")
    end

power(a::GenericOp{IsGenerator}, ::Number) =
    error("$(obs_name(a)) has no power: generators are added, not multiplied")

(a::GenericOp ^ p::Number) = power(a, p)

# a literal negative exponent would go through `inv`, which an operator has no method of
Base.literal_pow(::typeof(^), a::Op, ::Val{p}) where p = a ^ p

# the power of a placed operator is that of its operator; anything else placed only takes an
# integer power, its product
(a::AtIndex ^ p::Number) = power(a.op, p)(a.index...)

(a::IdentityOp{K, Indexed} ^ ::Number) where K = a

(a::IndexedOp ^ p::Number) =
    if !is_natural(p)
        error("$(obs_name(a)) has no power $p once placed: take the power of the operator " *
              "before placing it")
    elseif iszero(p)
        IdentityOp(a)
    else
        prod(fill(a, Int(real(p))))
    end

"""
    sqrt(::GenericOp)

the principal square root of a generic operator, `op^0.5`
"""
sqrt(a::GenericOp) = a ^ 0.5

"""
    ExpOp{K, N}

the exponential of a generic operator, of kind `K`: an operator for an operator, a channel for
a generator, see `exp`
"""
struct ExpOp{K, N} <: GenericOp{K, N}
    arg::GenericOp{<:Union{IsOperator, IsGenerator}, N}
end

"""
    exp(::GenericOp)

the exponential of a generic operator, or of a generic generator, the channel of the evolution
it generates.

# Examples

    exp(-im * t * (X ⊗ X + Y ⊗ Y) / 4)
    exp(t * (-im * (X ⊗ X) + Dissipator(Sm ⊗ Sm)))
"""
exp(a::GenericOp{IsOperator, N}) where N = ExpOp{IsOperator, N}(a)
exp(a::GenericOp{IsGenerator, N}) where N = ExpOp{IsChannel, N}(a)

"""
    DagOp{N}

the adjoint of a generic operator on pure states, see `dag`
"""
struct DagOp{N} <: GenericOp{IsOperator, N}
    arg::GenericOp{IsOperator, N}
    DagOp(arg::DagOp) = arg.arg
    DagOp(arg::ScalarOp{IsOperator, Generic}) = scale(conj(arg.coef), DagOp(arg.arg))
    DagOp(arg::GenericOp{IsOperator, N}) where N = new{N}(arg)
end

"""
    dag(op)
    op'

the adjoint of an operator on pure states, generic or placed on sites, a placed one having its
coefficients conjugated and its products reversed.

# Examples

    dag(C) * C
    h = sum(C'(i) * C(i + 1) for i in 1:9); H = -(h + h')
"""
dag(a::GenericOp{IsOperator}) = DagOp(a)
dag(a::AtIndex{IsOperator}) = AtIndex(dag(a.op), a.index)
dag(a::ScalarOp{IsOperator, Indexed}) = scale(conj(a.coef), dag(a.arg))
dag(a::SumOp{IsOperator, Indexed}) = SumOp(IndexedOp{IsOperator}[ dag(x) for x in a.subs ])
dag(a::ProdOp{IsOperator, Indexed}) =
    ProdOp(IndexedOp{IsOperator}[ dag(x) for x in reverse(a.subs) ])
dag(a::IdentityOp{IsOperator, Indexed}) = a

Base.adjoint(a::Op{IsOperator}) = dag(a)

"""
    ModOp{N}

an operator `A` taken modulo an integer `m` of at least 2, ``e^{2i\\pi A/m}``, see `mod`
"""
struct ModOp{N} <: GenericOp{IsOperator, N}
    arg::GenericOp{IsOperator, N}
    modulus::Int
    ModOp(arg::GenericOp{IsOperator, N}, modulus::Int) where N =
        if modulus < 2
            error("a modulus is at least 2, got $modulus")
        else
            new{N}(arg, modulus)
        end
end

"""
    mod(A, m)

the operator ``e^{2i\\pi A/m}``. Conserved, `conserve = mod(N, 3)` gives the charges
``n \\bmod 3``, added modulo 3 over the sites.
"""
mod(a::GenericOp{IsOperator}, m::Int) = ModOp(a, m)

"""
    parity(A)

the operator ``(-1)^A``, `mod(A, 2)`: conserved, a charge of ``\\mathbb{Z}_2``.

# Examples

    Fermion(conserve = parity(N))
"""
parity(a::GenericOp{IsOperator}) = ModOp(a, 2)

############### Purity ###############

"""
    is_pure(a)

whether the channel or generator `a` acts on pure states: a channel of a single Kraus operator,
a generator with no dissipator. A sum of channels is taken as a mixture.
"""
is_pure(::Union{IdentityOp, Gate, Evolver}) = true
is_pure(::Union{Reset, Dephase, Dissipator, SumOp{IsChannel}}) = false
is_pure(a::Union{SumOp{IsGenerator}, ProdOp, TensorOp}) = all(is_pure, a.subs)
is_pure(a::Union{ScalarOp, IntPowOp, ExpOp}) = is_pure(a.arg)
is_pure(a::AtIndex) = is_pure(a.op)

"""
    pure_part(a)

the operator `K` of the generator `a` with no dissipator,
``\\rho \\mapsto K\\rho + \\rho K^\\dagger``
"""
pure_part(a::Evolver) = a.arg
pure_part(::IdentityOp{IsGenerator, T, N}) where {T, N} = scale(0.5, IdentityOp{IsOperator, T, N}())
pure_part(a::ScalarOp{IsGenerator}) = scale(a.coef, pure_part(a.arg))
pure_part(a::SumOp{IsGenerator, T, N}) where {T, N} =
    SumOp(Op{IsOperator, T, N}[ pure_part(x) for x in a.subs ])
pure_part(a::AtIndex{IsGenerator}) = AtIndex(pure_part(a.op), a.index)

"""
    kraus(a)

the Kraus operator of the channel `a` of a single one: `K` for `Gate(K)`, ``e^{K}`` for the
exponential of `Evolver(K)`
"""
kraus(a::Gate) = a.arg
kraus(::IdentityOp{IsChannel, T, N}) where {T, N} = IdentityOp{IsOperator, T, N}()
kraus(a::ScalarOp{IsChannel}) = scale(sqrt(a.coef), kraus(a.arg))
kraus(a::ProdOp{IsChannel, T, N}) where {T, N} =
    ProdOp(Op{IsOperator, T, N}[ kraus(x) for x in a.subs ])
kraus(a::TensorOp{IsChannel}) = reduce(⊗, map(kraus, a.subs))
kraus(a::IntPowOp{IsChannel}) = kraus(a.arg)^a.expo
kraus(a::ExpOp{IsChannel}) = exp(pure_part(a.arg))
kraus(a::AtIndex{IsChannel}) = AtIndex(kraus(a.op), a.index)

############### Naming ###############

"""
    named(def, name)

the operator `name` defined by `def`: an expression of operators, a matrix, or a function of
the sites returning a matrix. A matrix or a function acts on one site, an expression on as many
as it does. Sites conserving operators of the same name share one charge, so renaming keeps
two quantities apart.

# Examples

    named([1 1 ; 1 -1] / √2, "MyH")
    Fermion(conserve = named(N, "Nf"))
"""
named(op::GenericOp{IsOperator, N}, name::String) where N = Operator{N}(name, op)
named(def::Union{AbstractMatrix, Function}, name::String) = Operator{1}(name, def)

# an operator broadcasts as a single value, as in `Level .=> bits`
Base.broadcastable(a::Op) = Ref(a)

################## Equality #################

# An operator is compared by what it is made of: the default `==` of an immutable struct
# compares a `Vector` field by identity. `hash` goes with it.

(a::Op == b::Op) =
    typeof(a) == typeof(b) &&
    all(f -> getfield(a, f) == getfield(b, f), fieldnames(typeof(a)))

hash(a::Op, h::UInt) =
    foldl((h, f) -> hash(getfield(a, f), h), fieldnames(typeof(a)); init = hash(typeof(a), h))

############## Printing ###############

"""
    show_func(io, name, args)

print `name(args)`, `args` being a vector of arguments or a single one
"""
function show_func(io::IO, name, args)
    print(io, name, "(")
    if args isa Vector
        join(io, (repr(a) for a in args), ",")
    else
        show(io, args)
    end
    print(io, ")")
end

"""
    print_number(io, x)

print the real number `x`, a float of integer value as that integer below `maxintfloat`
"""
print_number(io::IO, x::AbstractFloat) =
    if isinteger(x) && abs(x) < maxintfloat(Float64)
        print(io, Int(x))
    else
        print(io, x)
    end
print_number(io::IO, x::Real) = print(io, x)

"""
    print_coef(io, a)

print the number `a` in front of an operator: nothing for 1, `-` for -1, `2im*` for an
imaginary number, a rational or a complex number in parentheses
"""
print_coef(io::IO, a::Number) =
    if a ≠ 1
        if a == -1
            print(io, "-")
        elseif isa(a, Complex)
            if imag(a) == 0
                print_coef(io, real(a))
            elseif real(a) == 0
                print_coef(io, imag(a))
                print(io, "im*")
            elseif isfinite(a)
                print(io, "(")
                print_number(io, real(a))
                print(io, imag(a) < 0 ? " - " : " + ")
                print_number(io, abs(imag(a)))
                print(io, "im)")
            else
                print(io, "(", a, ")")
            end
        elseif isa(a, Rational)
            # (1//2)X and not 1//2X, which reads 1//(2X)
            print(io, "(", a, ")")
        else
            print_number(io, a)
        end
    end

"""
    infix(io, subs, op)

print `subs` separated by `op`, the operands after the first a precedence higher, `*` and `⊗`
being left associative
"""
function infix(io::IO, subs, op::String)
    p = get(io, :precedence, 0)
    for (k, s) in enumerate(subs)
        if k > 1
            print(io, op)
        end
        print(k == 1 ? io : IOContext(io, :precedence => p + 1), s)
    end
end

"""
    paren(f, io, out_prec, in_prec = out_prec)

print through `f(io)` an expression of precedence `out_prec`, in parentheses when the
surrounding precedence is higher, its inside at precedence `in_prec`
"""
function paren(f, io::IO, out_prec::Int, in_prec::Int = out_prec)
    ext_prec = get(io, :precedence, 0)
    if out_prec < ext_prec
        print(io, "(")
    end
    f(IOContext(io, :precedence => in_prec))
    if out_prec < ext_prec
        print(io, ")")
    end
end

"""
    shortened(io, xs)

`xs`, cut down to its first three and its last two elements around `"..."` when `io` is
compact and `xs` has more than six
"""
function shortened(io::IO, xs::AbstractVector)
    n = length(xs)
    return n > 6 && get(io, :compact, false) ? [xs[1:3]; "..."; xs[n-1:n]] : xs
end

show(io::IO, op::Operator) = print(io, op.name)
show(io::IO, ::IdentityOp{IsOperator, Generic, N}) where N = print(io, join(fill("Id", N), "⊗"))
show(io::IO, ::IdentityOp{IsOperator, Indexed}) = print(io, "Id")
show(io::IO, ::IdentityOp{IsChannel, Generic, N}) where N =
    print(io, "Gate(", join(fill("Id", N), "⊗"), ")")
show(io::IO, ::IdentityOp{IsChannel, Indexed}) = print(io, "Gate(Id)")
# ρ ↦ ρ, as a generator
show(io::IO, ::IdentityOp{IsGenerator, Generic, N}) where N =
    print(io, "Evolver(0.5", join(fill("Id", N), "⊗"), ")")
show(io::IO, ::IdentityOp{IsGenerator, Indexed}) = print(io, "Evolver(0.5Id)")
show(io::IO, a::Proj) = show_func(io, "Proj", [a.state])
show(io::IO, a::Reset) = show_func(io, "Reset", [a.state])
function show(io::IO, ind::AtIndex)
    print(IOContext(io, :precedence => 500), ind.op)
    print(io, "(", join(shortened(io, collect(ind.index)), ","), ")")
end

show(io::IO, a::ScalarOp) =
    paren(io, Base.operator_precedence(:*)) do io
        print_coef(io, a.coef)
        print(io, a.arg)
    end

show(io::IO, a::SumOp) =
    paren(io, Base.operator_precedence(:+)) do io
        n = length(a.subs)
        subs = if n > 6 && get(io, :compact, false)
            dots = startswith(sprint(print, a.subs[4]; context = io), '-') ? "-..." : "..."
            [a.subs[1:3]; dots; a.subs[n-1:n]]
        else
            a.subs
        end
        for (k, x) in enumerate(subs)
            s = sprint(print, x; context = io)
            if k > 1 && s[1] ≠ '-'
                print(io, "+")
            end
            print(io, s)
        end
    end

show(io::IO, a::ProdOp) =
    paren(io, Base.operator_precedence(:*)) do io
        infix(io, shortened(io, a.subs), "*")
    end

show(io::IO, a::TensorOp) =
    paren(io, Base.operator_precedence(:⊗)) do io
        infix(io, shortened(io, a.subs), "⊗")
    end

show(io::IO, a::Union{IntPowOp, GenPowOp}) =
    paren(io, Base.operator_precedence(:^)) do io
        # the base a precedence higher, ^ being right associative: (X^0.5)^0.5
        print(IOContext(io, :precedence => Base.operator_precedence(:^) + 1), a.arg)
        print(io, "^", a.expo isa Union{Complex, Rational} ? "($(a.expo))" : a.expo)
    end

show(io::IO, a::ExpOp) =
    paren(io, 1000, 0) do io
        show_func(io, "exp", a.arg)
    end

show(io::IO, a::DagOp) =
    paren(io, 1000, 0) do io
        show_func(io, "dag", a.arg)
    end

show(io::IO, a::Union{Gate, Dephase, Dissipator, Evolver}) =
    paren(io, 1000, 0) do io
        show_func(io, nameof(typeof(a)), a.arg)
    end

show(io::IO, a::ModOp) =
    paren(io, 1000, 0) do io
        if a.modulus == 2
            show_func(io, "parity", a.arg)
        else
            show_func(io, "mod", [a.arg, a.modulus])
        end
    end
