# Site types: the declarations of a site type of one's own, the matrices of the operators of
# the built in types, their local states, what they conserve and the indices they give, and
# what operators measure and do on their states, whether the sites conserve something or not.

# the conserved quantities of a site, as tuples
charges(site) = [ (q.name, q.modulus, q.charges, q.strong) for q in TMS.conserved(site).quantities ]

# Defining a site type from outside the package, which is what a user does and what the built
# in types cannot check: they are declared inside TensorMixedStates.

struct Dummit <: AbstractSite end

TensorMixedStates.dim(::Dummit) = 2

@def_states Dummit [ "1" => [0., 1.] ]

# `N` reaches this module through the `using` of runtests.jl: `@def_operators` must register
# N for Dummit without binding the name again
@def_operators Dummit [ N = [0. 0. ; 0. 1.] ]

# a charge modulo 2 declared for the site, whose definition lives in the library and not in the
# operator: its eigenvalues ±1 were read as an integer charge
@def_operators Dummit [ Pd = parity(N) ]

@create_site_module(Dummits, [Dummit, N])

struct Dummit2 <: AbstractSite end

TensorMixedStates.dim(::Dummit2) = 2

# an operator with a definition of its own, whose name a site declaration must not take
const Renamed = named(parity(N), "Renamed")

# a site carrying nothing in a field that is not the last one: printing must keep it
struct Middling <: AbstractSite
    a::Union{Nothing, Int}
    b::Int
end

TensorMixedStates.dim(::Middling) = 2

# a site whose field happens to be named conserve without being a Conservation
struct Pretender <: AbstractSite
    conserve::String
end

TensorMixedStates.dim(::Pretender) = 2

# a site declaring a state spelt as a generic form, and one of its own for FullyMixed
struct Declarer <: AbstractSite end

TensorMixedStates.dim(::Declarer) = 2

@def_states Declarer [ "0" => [0., 1.], "FullyMixed" => [0.9 0. ; 0. 0.1] ]

# a site whose rule of names is interrupted
struct Interrupter <: AbstractSite end

TensorMixedStates.dim(::Interrupter) = 2

TensorMixedStates.string_state(::Interrupter, ::String) = throw(InterruptException())

# a fermionic site of one's own, its parities giving F
struct Oddling <: AbstractSite end

TensorMixedStates.dim(::Oddling) = 2

TensorMixedStates.fermion_parities(::Oddling) = (0, 1)

@def_operators Oddling [ Co = [0. 1. ; 0. 0.] ]

TensorMixedStates.fermion_species(::Oddling) = (Co,)

# sites whose declarations do not fit them, refused at their first use
struct Misdeclared <: AbstractSite end

TensorMixedStates.dim(::Misdeclared) = 2

@def_states Misdeclared [ "Big" => [1., 0., 0.] ]

struct Misoperated <: AbstractSite end

TensorMixedStates.dim(::Misoperated) = 2

@def_operators Misoperated [ Big = [1. 0. 0. ; 0. 1. 0. ; 0. 0. 1.] ]

struct Misparity <: AbstractSite end

TensorMixedStates.dim(::Misparity) = 2

TensorMixedStates.fermion_parities(::Misparity) = (0, 2)

# a site type declaring Sp and Sz only, and one declaring A only, which get the relations declared
# once for every site type
struct Spinner <: AbstractSite end

TensorMixedStates.dim(::Spinner) = 2

@def_operators Spinner [ Sp = [0. 1. ; 0. 0.], Sz = [0.5 0. ; 0. -0.5] ]

struct Laddered <: AbstractSite end

TensorMixedStates.dim(::Laddered) = 3

@def_operators Laddered [ A = [0. 1. 0. ; 0. 0. √2 ; 0. 0. 0.] ]

@testset "Qubit measuring" begin
    @test norm(matrix(Sx^2 + Sy^2 + Sz^2 - S2, Qubit())) ≈ 0 atol = 1e-12
    # X, Y and Z, alone and in pairs, on their eigenstates, pure and mixed
    for rep in (identity, mix)
        st(n, v) = rep(State{Pure}(System(n, Qubit()), v))
        @test expect(st(1, "Z+"), [X(1), Y(1), Z(1)]) ≈ [0, 0, 1]
        @test expect1(st(6, ["X+", "Y+", "Z+", "X-", "Y-", "Z-"]), [X, Y, Z]) ≈
            [[1, 0, 0, -1, 0, 0], [0, 1, 0, 0, -1, 0], [0, 0, 1, 0, 0, -1]]
        @test expect2(st(3, ["X+", "Y-", "Z-"]), [(X, Y), (Z, Z)]) ≈
            [[0 -1 0; 0 0 0; 0 0 -im], [1 0 0; 0 1 0; 0 0 1]]
        @test expect(st(4, ["Z+", "X-", "Z-", "Y-"]), [Z(1)Y(4), X(2)Z(3), Z(1)Z(3), X(2)Y(4)]) ≈
            [-1, 1, -1, 1]
    end
end

@testset "Qubit excitation number" begin
    q = Qubit()
    # N counts "Dn" as the occupied state, so that it agrees with Sz rather than against it
    @test matrix(N, q) ≈ [0. 0. ; 0. 1.]
    @test matrix(N, q) ≈ matrix(Proj(Level => 1), q)
    @test matrix(N, q) ≈ matrix((Id - Z) / 2, q)
    @test matrix(N, q) ≈ matrix(0.5 * Id - Sz, q)
    # the name is shared with the sites that already had it
    @test Qubits.N === Bosons.N === Fermions.N
    # conserving it is the same conservation as conserving Sz, with integer charges
    @test charges(Qubit(conserve = N)) == [("N", 1, [0, 1], false)]
    @test charges(Qubit(conserve = 2Sz)) == [("2Sz", 1, [1, -1], false)]
    # "Up" being empty, it is Sm that creates an excitation
    @test flux(Sm, Qubit(conserve = N)) == IT.QN("N", 1)
    @test flux(Sp, Qubit(conserve = N)) == IT.QN("N", -1)
    @test_throws "no definite flux" flux(X, Qubit(conserve = N))
end

@testset "Spin excitation number" begin
    # N counts the excitations above the state of maximal Sz
    for x in (1/2, 1, 3/2, 2)
        site = Spin(x)
        @test matrix(N, site) ≈ matrix(x * Id - Sz, site)
        @test all(k -> matrix(N, site)[k, k] ≈ k - 1, 1:dim(site))
    end
    @test matrix(N, Spin(1/2)) ≈ matrix(N, Qubit())
    @test Spins.N === Qubits.N === Bosons.N
    @test charges(Spin(1/2, conserve = N)) == [("N", 1, [0, 1], false)]
    @test charges(Spin(3/2, conserve = N)) == [("N", 1, [0, 1, 2, 3], false)]
    @test charges(Spin(3/2, conserve = 2Sz)) == [("2Sz", 1, [3, 1, -1, -3], false)]
    @test flux(Sm, Spin(1, conserve = N)) == IT.QN("N", 1)
    @test flux(Sp, Spin(1, conserve = N)) == IT.QN("N", -1)
end

@testset "Level numbers the basis states of every site" begin
    for site in (Qubit(), Fermion(), Boson(4), Qboson(0.1, 4), Qudit(5), Spin(3/2))
        @test matrix(Level, site) ≈ matrix(N, site)
    end
    numbers(d) = [ i == j ? i - 1. : 0. for i in 1:d, j in 1:d ]
    @test matrix(Level, Electron()) ≈ numbers(4)
    @test matrix(Level, Tj()) ≈ numbers(3)
    @test matrix(Level, Dummit()) ≈ numbers(2)
    @test flux(Level, Electron(conserve = (Ntot, 2Sz))) == IT.QN()
end

@testset "The parity of a site gives F" begin
    # decision 19: F is no longer declared, the parities of the basis states give it
    @test matrix(F, Fermion()) == [1. 0. ; 0. -1.]
    @test matrix(F, Electron()) == [1. 0 0 0 ; 0 -1 0 0 ; 0 0 -1 0 ; 0 0 0 1]
    @test matrix(F, Tj()) == [1. 0 0 ; 0 -1 0 ; 0 0 -1]
    @test matrix(F, Oddling()) == [1. 0. ; 0. -1.]
    for site in (Qubit(), Boson(3), Spin(1), Qudit(3), Qboson(0.1, 3), Dummit())
        @test matrix(F, site) == matrix(Id, site)
    end
    # F is the product of the parities of the two spins, parity(Nup) and parity(Ndn)
    for site in (Electron(), Tj())
        @test matrix(F, site) ≈ matrix(parity(Nup) * parity(Ndn), site)
    end
    @test fermion_species(Oddling()) == (Co,)
    @test fermion_species(Qubit()) === nothing
    @test_throws "F is not declared" @macroexpand @def_operators Dummit2 [ F = [1. 0. ; 0. -1.] ]
end

@testset "Swap exchanges the states of two sites of one type" begin
    # the permutation of the two states on every site type, a sign coming between two odd ones
    for site in (Qubit(), Qudit(3), Boson(4), Spin(1))
        d = dim(site)
        m = matrix(Swap, site)
        for i in 1:d, j in 1:d
            v = zeros(d^2)
            v[(i - 1) * d + j] = 1
            w = zeros(d^2)
            w[(j - 1) * d + i] = 1
            @test m * v == w
        end
    end
    q = Qubit()
    @test matrix(Swap, q) ≈ matrix((Id ⊗ Id + Z ⊗ Z) / 2 + Sp ⊗ Sm + Sm ⊗ Sp, q)
    # the fermionic swap exchanges the modes, a string included: P (C ⊗ Id) P is Id ⊗ C
    for (site, c) in ((Fermion(), C), (Electron(), Cup), (Electron(), Cdn), (Tj(), Cdn))
        p = matrix(Swap, site)
        @test p * matrix(c ⊗ Id, site) * p ≈ matrix(Id ⊗ c, site)
    end
    # its matrix knows no charge, and placed on a system it keeps the charges only when both
    # sites carry the same ones
    @test matrix(Swap, Fermion(conserve = N), Fermion()) == matrix(Swap, Fermion())
    mixed = State{Pure}(System([Fermion(conserve = N), Fermion()]), ["Occ", "Emp"])
    @test_throws "carries no definite charge of N" make_mpo(mixed, Swap(1, 2))
    @test_throws "of one type and one dimension" matrix(Swap, Qubit(), Fermion())
    @test_throws "of one type and one dimension" matrix(Swap, Boson(3), Boson(4))
end

@testset "Relations declared once for every site type" begin
    # a site type declaring Sp and Sz gets Sm, Sx, Sy and S2, one declaring A gets Q and P
    @test matrix(Sm, Spinner()) == [0. 0. ; 1. 0.]
    @test matrix(Sx, Spinner()) ≈ matrix(Sx, Qubit())
    @test matrix(Sy, Spinner()) ≈ matrix(Sy, Qubit())
    @test matrix(S2, Spinner()) ≈ 0.75 * matrix(Id, Spinner())
    @test matrix(Q, Laddered()) ≈ matrix(Q, Boson(3))
    @test matrix(P, Laddered()) ≈ matrix(P, Boson(3))
    # and so do the site types of the library
    for (site, ops) in ((Qubit(), (Sm, Sx, Sy, S2)), (Spin(3/2), (Sm, Sx, Sy, S2)),
                        (Electron(), (Sm, Sx, Sy, S2, Nup, Ndn, Ntot)),
                        (Tj(), (Sm, Sx, Sy, S2, Nup, Ndn, Ntot)),
                        (Boson(4), (Q, P)), (Qboson(0.5, 4), (Q, P)))
        for op in ops
            @test size(matrix(op, site)) == (dim(site), dim(site))
        end
    end
    @test matrix(Ntot, Electron()) ≈ matrix(Nup + Ndn, Electron())
    # S2 is written (Sp * Sm + Sm * Sp) / 2 + Sz^2, which is Sx^2 + Sy^2 + Sz^2 and real where Sp
    # and Sz are
    for site in (Qubit(), Electron(), Tj(), Spinner())
        @test matrix(S2, site) isa Matrix{Float64}
        @test matrix(S2, site) ≈ matrix(Sx^2 + Sy^2 + Sz^2, site)
    end
    @test matrix(S2, Qubit()) == 0.75 * matrix(Id, Qubit())
    # the declaration of a site type wins over the relation, and a relation holds only where
    # what it is made of is declared, the check at first use leaving it alone
    @test matrix(S2, Spin(1)) == 2 * matrix(Id, Spin(1))
    @test_throws "operator Sp is not defined for site Fermion" matrix(Sx, Fermion())
    @test_ok System(2, Spinner())
    @test_ok System(2, Fermion())
end

@testset "Fermion measuring" begin
    # N on the states named by their level, pure and mixed
    for rep in (identity, mix)
        @test expect(rep(State{Pure}(System(1, Fermion()), "1")), N(1)) ≈ 1
        @test expect1(rep(State{Pure}(System(2, Fermion()), ["0", "1"])), N) ≈ [0, 1]
    end
end

@testset "Boson measuring" begin
    for rep in (identity, mix)
        st = rep(State{Pure}(System(4, Boson(4)), ["0", "1", "2", "3"]))
        @test expect1(st, N) ≈ [0, 1, 2, 3]
    end
end

@testset "Spin measuring" begin
    @test norm(matrix(Sx^2 + Sy^2 + Sz^2 - S2, Spin(5/2))) ≈ 0 atol = 1e-12
    @test norm(matrix(Sx^2 + Sy^2 + Sz^2 - S2, Spin(4))) ≈ 0 atol = 1e-12
    # the eigenstates named by a value of Sz, Sx or Sy
    value(m) = isinteger(m) ? string(Int(m)) : "$(Int(2m))/2"
    for s in (1, 3/2), axis in ("X", "Y", "Z"), m in -s:s
        v = TMS.local_state(Spin(s), axis * value(m))
        a = matrix(axis == "X" ? Sx : axis == "Y" ? Sy : Sz, Spin(s))
        @test a * v ≈ m * v atol = 1e-12
    end
    @test_throws "state 5/2 is not defined for site Spin" TMS.local_state(Spin(3/2), "5/2")
    # and Sx, Sy and Sz measured on those states, pure and mixed
    for rep in (identity, mix)
        st(s, v) = rep(State{Pure}(System(length(v), Spin(s)), v))
        @test expect1(st(3/2, ["-3/2", "-1/2", "1/2", "3/2"]), [Sx, Sy, Sz]) ≈
            [[0, 0, 0, 0], [0, 0, 0, 0], [-3/2, -1/2, 1/2, 3/2]]
        on_axes(ms) = [ a * m for a in ("X", "Y", "Z") for m in ms ]
        @test expect1(st(1, on_axes(["-1", "0", "1"])), [Sx, Sy, Sz]) ≈
            [[-1, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, -1, 0, 1, 0, 0, 0],
             [0, 0, 0, 0, 0, 0, -1, 0, 1]]
        @test expect1(st(3/2, on_axes(["-3/2", "-1/2", "1/2", "3/2"])), [Sx, Sy, Sz]) ≈
            [[-3/2, -1/2, 1/2, 3/2, 0, 0, 0, 0, 0, 0, 0, 0],
             [0, 0, 0, 0, -3/2, -1/2, 1/2, 3/2, 0, 0, 0, 0],
             [0, 0, 0, 0, 0, 0, 0, 0, -3/2, -1/2, 1/2, 3/2]]
    end
end

@testset "Electron measuring" begin
    for rep in (identity, mix)
        st = rep(State{Pure}(System(4, Electron()), ["Emp", "Up", "Dn", "UpDn"]))
        @test expect1(st, [Nup, Ndn, Nupdn, Ntot]) ≈
            [[0, 1, 0, 1], [0, 0, 1, 1], [0, 0, 0, 1], [0, 1, 1, 2]]
    end
end

@testset "Tj measuring" begin
    for rep in (identity, mix)
        st = rep(State{Pure}(System(3, Tj()), ["Emp", "Up", "Dn"]))
        @test expect1(st, [Nup, Ndn, Ntot]) ≈ [[0, 1, 0], [0, 0, 1], [0, 1, 1]]
    end
end

@testset "One electron of fully mixed spin" begin
    for site in (Electron(conserve = strong(Ntot)), Tj(conserve = strong(Ntot)))
        st = State{Mixed}(System(2, site), ["MixedSpin", "↑|↓"])
        @test trace(st) ≈ 1
        @test real(expect(st, Ntot(1))) ≈ 1
        @test real(expect(st, Nup(2))) ≈ 0.5
        @test abs(expect(st, Sz(1))) < 1e-14
        @test real(purity(st)) ≈ 1 / 4                     # trace2 in 1.x
        @test_throws "spreads over several charges" State{Mixed}(System(2, site), "FullyMixed")
    end
end

@testset "Site parameters out of their range are refused" begin
    @test_throws "half integer of 0 or more" Spin(-1/2)
    @test_throws "1 level at least" Boson(0)
    @test_throws "dimension of 1 at least" Qudit(0)
    @test_throws "1 - q^n < 0" Qboson(1.5, 3)
    @test_throws "1 - q^n < 0" Qboson(-2.0, 3)
    @test dim(Spin(0)) == 1
    @test matrix(A, Qboson(-2.0, 2)) ≈ [0 sqrt(3); 0 0]
end

@testset "Qboson measuring" begin
    for rep in (identity, mix)
        st = rep(State{Pure}(System(4, Qboson(0.1, 4)), ["0", "1", "2", "3"]))
        @test expect1(st, N) ≈ [0, 1, 2, 3]
    end
end

@testset "State mixing" begin
    # a pure state and its density matrix give the same expectation values, on random states
    # of several bond dimensions
    for (site, ops) in ((Qubit(), [X, Y, Z]), (Fermion(), [N]), (Boson(3), [N]),
                        (Spin(3/2), [Sx, Sy, Sz]),
                        (Electron(), [Nup, Ndn, Nupdn, Ntot, Sx, Sy, Sz]),
                        (Tj(), [Nup, Ndn, Ntot, Sx, Sy, Sz]))
        for χ in (1, 5, 15)
            ψ = random_state(Pure, System(10, site), χ)
            @test expect1(ψ, ops) ≈ expect1(mix(ψ), ops)
        end
    end
end

@testset "Mix ordering" begin
    # a gate applied then the state mixed, or the state mixed then the gate applied
    ψ = State{Pure}(System(2, Qubit()), "0")
    g = (Id ⊗ X)(1, 2)
    for ρ in (mix(apply(g, ψ)), apply(g, mix(ψ)))
        @test expect(ρ, Z(1)) ≈ 1
        @test expect(ρ, Z(2)) ≈ -1
    end
end

@testset "Qudit operators" begin
    # the clock and shift operators are unitary, of order d, and obey the Weyl relation
    for d in (2, 3, 5)
        s = Qudit(d)
        ω = exp(2im * π / d)
        z = matrix(Z, s)
        x = matrix(X, s)
        id = matrix(Id, s)
        @test z^d ≈ id
        @test x^d ≈ id
        @test z * z' ≈ id
        @test x * x' ≈ id
        @test z * x ≈ ω * (x * z)
        @test matrix(N, s) ≈ [ i == j ? i - 1. : 0. for i in 1:d, j in 1:d ]
        h = matrix(H, s)
        @test h * h' ≈ id
        @test h^4 ≈ id
        @test h * x * inv(h) ≈ z
        @test h * z * inv(h) ≈ inv(x)
        p = matrix(S, s)
        @test p * p' ≈ id
        a = p * x * inv(p)
        b = x * z
        @test a ≈ (a[findfirst(!=(0), a)] / b[findfirst(!=(0), a)]) * b
        @test p^(isodd(d) ? d : 2d) ≈ id
    end
    @test matrix(Z, Qudit(2)) ≈ matrix(Z, Qubit())
    @test matrix(X, Qudit(2)) ≈ matrix(X, Qubit())
    @test matrix(H, Qudit(2)) ≈ matrix(H, Qubit())
    @test matrix(S, Qudit(2)) ≈ matrix(Qubits.S, Qubit())
    @test matrix(X, Qudit(1)) == ones(1, 1)
    @test matrix(Z, Qudit(1)) == ones(1, 1)
    @test matrix(Sum, Qudit(2), Qudit(2)) ≈ matrix(controlled(X), Qubit(), Qubit())
    # Sum adds the level of the first qudit to the second, modulo the dimension of the second,
    # which it reads on its sites
    m = matrix(Sum, Qudit(3))
    for i in 0:2, j in 0:2
        @test m[3i + mod(i + j, 3) + 1, 3i + j + 1] ≈ 1
    end
    m = matrix(Sum, Qudit(2), Qudit(3))
    for i in 0:1, j in 0:2
        @test m[3i + mod(i + j, 3) + 1, 3i + j + 1] ≈ 1
    end
    # on states: Sum adds the level of the first qudit to the second, and Z is the phase of the level
    st = apply(Sum(1, 2), State{Pure}(System(2, Qudit(3)), ["1", "2"]))
    @test expect1(st, N) ≈ [1, 0]
    @test expect(State{Pure}(System(1, Qudit(3)), "1"), Z(1)) ≈ exp(2im * π / 3)
    # on sites next to each other or apart, and, given by a function of its sites, in an MPO
    st = State{Pure}(System(3, Qudit(3)), ["1", "1", "0"])
    @test real(expect1(st, N)) ≈ [1, 1, 0]
    @test real(expect1(apply(Sum(1, 2), st), N)) ≈ [1, 2, 0]
    @test real(expect1(apply(Sum(1, 3), st), N)) ≈ [1, 1, 1]
    @test maxlinkdim(make_mpo(st, Sum(1, 2))) == 4
    for rep in (identity, mix)
        @test expect1(rep(State{Pure}(System(3, Qudit(3)), ["0", "1", "2"])), Z) ≈
            [exp(2im * π * n / 3) for n in 0:2]
    end
end

@testset "Unknown operator and state names" begin
    @test_throws "operator X is not defined for site Boson" matrix(X, Boson(3))
    @test_throws "operator A is not defined for site Fermion" matrix(A, Fermion())
    @test_throws "state Zorglub is not defined for site Qubit" State{Pure}(System(2, Qubit()), "Zorglub")
    @test_throws "state Zorglub is not defined for site Dummit" State{Pure}(System(2, Dummit()), "Zorglub")
    @test_throws "a state of 3 components cannot be one of Qubit()" State{Pure}(System(2, Qubit()), [1., 0., 0.])
    @test_throws "a 3×3 density matrix cannot be one of Qubit()" State{Mixed}(
        System(2, Qubit()), [1. 0. 0. ; 0. 0. 0. ; 0. 0. 0.])
    @test_throws "a state of 3 components cannot be one of Qubit()" matrix(Reset([1., 0., 0.]), Qubit())
end

@testset "A state a site declares is the one used" begin
    ls = TMS.local_state
    @test ls(Declarer(), "0") == [0., 1.]
    @test ls(Declarer(), "FullyMixed") == [0.9 0. ; 0. 0.1]
    # its rule of names is still there, and FullyMixed for a site that declares none
    @test ls(Declarer(), "1") == [0., 1.]
    @test ls(Qubit(), "FullyMixed") == [0.5 0. ; 0. 0.5]
    @test_throws "state Zorglub is not defined for site Declarer" ls(Declarer(), "Zorglub")
    # decision 20: an error of the rule of names is no longer taken for a name it does not read
    @test_throws InterruptException ls(Interrupter(), "Up")
    # a copy of what the site declares, which changing changed it for everyone
    v = ls(Qubit(), "Up")
    v[1], v[2] = 0, 1
    @test ls(Qubit(), "Up") == [1, 0]
    m = matrix(X, Qubit())
    m[1, 2] = 5
    @test matrix(X, Qubit()) == [0 1 ; 1 0]
end

@testset "Local states" begin
    ls = TMS.local_state
    # decision 20: A => λ is an eigenstate; a basis vector for a diagonal operator and a simple
    # eigenvalue, the uniform mixture of the eigenspace otherwise
    @test ls(Qudit(4), Level => 2) == [0., 0., 1., 0.]
    @test ls(Spin(1), Sz => -1) == [0., 0., 1.]
    @test ls(Electron(), Ntot => 1) ≈ [0 0 0 0 ; 0 0.5 0 0 ; 0 0 0.5 0 ; 0 0 0 0]
    @test ls(Qubit(), X => 1) ≈ [0.5 0.5 ; 0.5 0.5]
    @test ls(Qubit(), Id => 1) ≈ ls(Qubit(), "FullyMixed")
    @test_throws "not an eigenvalue" ls(Qubit(), Z => 2)
    @test_throws "needs a Hermitian operator" ls(Qubit(), Sp => 0)
    @test_throws "a number does not denote a state" ls(Qubit(), 1)
    @test_throws "is not a local state" ls(Qubit(), :Up)
    @test ls(Qubit(), s -> [0., 1.]) == [0., 1.]
    # in a pure state, a density matrix is refused, with what A => λ needs to give a vector
    @test_throws "A => λ gives a vector only" State{Pure}(System(2, Qubit()), X => 1)
    @test_ok State{Mixed}(System(2, Qubit()), X => 1)
    @test_ok State{Pure}(System(2, Spin(1)), Sz => 0)
end

@testset "Custom site type" begin
    @test dim(Dummit()) == 2
    # N arrives here by `using` from another site module: declaring it again for a new site
    # neither fails nor disturbs the site it came from, and the name goes on standing for one
    # and the same operator
    @test matrix(N, Dummit()) == [0. 0. ; 0. 1.]
    @test matrix(N, Boson(4)) ≈ TMS.LinearAlgebra.diagm(0:3)
    @test Dummits.Dummit === Dummit
    @test Dummits.N === N
    @test N === Bosons.N
    @test real(expect1(State{Pure}(System(2, Dummit()), "1"), N)) ≈ [1, 1]
    @test real(expect1(State{Pure}(System(2, Boson(4)), ["1", "3"]), N)) ≈ [1, 3]
    # a name standing for an operator with a definition of its own cannot be declared
    @test_throws "definition of its own" @def_operators Dummit2 [ Renamed = [0. 0. ; 0. 1.] ]
    @test_throws "operator Renamed is not defined for site Dummit2" matrix(Operator{1}("Renamed"), Dummit2())
    # declarations take the type, and a name once
    @test_throws "takes the site type, not a site" @macroexpand @def_states Dummit2() [ "1" => [0., 1.] ]
    @test_throws "declares the operator X twice" @macroexpand @def_operators Dummit2 [ X = [0. 1. ; 1. 0.], X = [0. 1. ; 1. 0.] ]
    @test_throws "declares the state 1 twice" @macroexpand @def_states Dummit2 [ "1" => [0., 1.], ["0", "1"] => [1., 0.] ]
    @test_throws "a literal string" @macroexpand @def_states Dummit2 [ name => [0., 1.] ]
    # the declarations of a site are checked at its first use, in a system
    @test_throws "a state of 3 components cannot be one of Misdeclared()" System(2, Misdeclared())
    @test_throws "operator Big is given by a 3×3 matrix and cannot act on Misoperated()" System(2, Misoperated())
    @test_throws "fermion_parities gives [0, 2]" System(2, Misparity())
    @test_ok System(2, Oddling())
end

@testset "Index tags" begin
    # the tag says which site type the index belongs to, whatever the modules in scope
    for site in [Qubit(), Fermion(), Boson(4), Spin(3/2), Electron(), Tj(), Qboson(0.1, 3),
                 Qudit(3), Dummit()]
        i = siteind(Pure, System(1, site), 1)
        @test IT.hastags(i, string(nameof(typeof(site))))
        @test IT.hastags(i, "Site")
    end
end

@testset "Declaring a conserved quantity" begin
    # a power or the adjoint of a charge modulo d is one modulo d
    @test charges(Qudit(2, conserve = dag(Z))) == [("dag(Z)", 2, [0, 1], false)]
    @test charges(Qudit(4, conserve = Z^2)) == [("Z^2", 4, [0, 2, 0, 2], false)]
    @test flux(X, Qudit(2, conserve = dag(Z))) == IT.QN("dag(Z)", 1, 2)
    # X changes Z^2 by 2 modulo 4 whatever the level, so X ⊗ X keeps it
    z2 = State{Pure}(System(3, Qudit(4, conserve = Z^2)), "0")
    @test maxlinkdim(make_mpo(z2, sum((X ⊗ X)(i, i + 1) for i in 1:2))) == 3
    @test charges(Boson(3, conserve = N^2)) == [("N^2", 1, [0, 1, 4], false)]
    # a quantity given twice is refused
    @test_throws "Ntot is given twice in conserve" Electron(conserve = (Ntot, Ntot))
    @test_throws "Ntot is given twice in conserve" Electron(conserve = (strong(Ntot), Ntot))
    # large charges are compared to a relative tolerance
    @test_ok Boson(80, conserve = dag(A) * A)
    @test_ok Spin(10, conserve = Sp * Sm - Sm * Sp)
    # the charges of every basis state, with the modulus
    @test charges(Fermion(conserve = N)) == [("N", 1, [0, 1], false)]
    @test charges(Fermion(conserve = parity(N))) == [("parity(N)", 2, [0, 1], false)]
    @test charges(Fermion(conserve = named(N, "Nf"))) == [("Nf", 1, [0, 1], false)]
    # a renamed quantity keeps the modulus of what it renames
    @test charges(Fermion(conserve = named(parity(N), "P"))) == [("P", 2, [0, 1], false)]
    paired = State{Pure}(System(4, Fermion(conserve = named(parity(N), "P"))), "Emp")
    @test_ok make_mpo(paired, sum(C(i) * C(i + 1) + dag(C)(i + 1) * dag(C)(i) for i in 1:3))
    @test charges(Boson(6, conserve = mod(N, 3))) == [("mod(N,3)", 3, [0, 1, 2, 0, 1, 2], false)]
    @test charges(Spin(1, conserve = Sz)) == [("Sz", 1, [1, 0, -1], false)]
    @test charges(Spin(1/2, conserve = 2Sz)) == [("2Sz", 1, [1, -1], false)]
    @test charges(Electron(conserve = (Ntot, 2Sz))) ==
        [("2Sz", 1, [0, 1, -1, 0], false), ("Ntot", 1, [0, 1, 1, 2], false)]
    @test isempty(charges(Qubit()))
    # the modulus of Z coming from its definition, mod(N, d)
    for d in (3, 5, 8)
        @test charges(Qudit(d, conserve = Z)) == [("Z", d, collect(0:d-1), false)]
    end
    @test charges(Qudit(2, conserve = Z)) == [("Z", 2, [0, 1], false)]
    @test TMS.site_charges(Z, Qudit(2)) == (2, [0, 1])
    # what cannot be a charge
    @test_throws "is not diagonal" Spin(1, conserve = Sx)
    @test_throws "2Sz rather than Sz" Spin(1/2, conserve = Sz)
    @test_throws "already carries a charge modulo" Qudit(3, conserve = mod(Z, 2))
    # a unitary given by its matrix is not read off its spectrum: its charge and the modulus
    # are written with mod
    @test_throws "a charge modulo m is written mod(Q, m)" Qubit(conserve = S)
    # decision 18: the charges are typed, so a name may hold any character but a final star,
    # which names the bra of a strong quantity
    for name in ("a%2", "a:b", "a;b", "a!")
        @test charges(Boson(3, conserve = named(N, name))) == [(name, 1, [0, 1, 2], false)]
    end
    @test_throws "a name ending in *" Boson(3, conserve = named(N, "a*"))
    # one name for two moduli on the same system
    @test_throws "modulo 2 on one site and modulo 3" System([Qudit(2, conserve = Z), Qudit(3, conserve = Z)])
    # a conserved quantity is carried by one site
    @test_throws "acting on one site" Fermion(conserve = Swap)
    @test_throws "acting on one site" Qubit(conserve = mod(Swap, 2))
    # the four cases, read directly
    @test TMS.site_charges(N, Fermion()) == (1, [0, 1])
    @test TMS.site_charges(parity(N), Boson(4)) == (2, [0, 1, 0, 1])
    @test TMS.site_charges(Z, Qudit(3)) == (3, [0, 1, 2])
    @test TMS.site_charges(Pd, Dummit()) == (2, [0, 1])
    # the string of versions 1 and 2 is read back as it was written
    @test TMS.decode_conserve("2Sz:0,1,-1,0;Ntot:0,1,1,2") == TMS.conserved(Electron(conserve = (Ntot, 2Sz)))
    @test TMS.decode_conserve("") == Conservation()
    # a site prints under the names of its charges
    @test repr(Qubit()) == "Qubit()"
    @test repr(Boson(4)) == "Boson(4)"
    @test repr(Spin(3/2)) == "Spin(1.5)"
    @test repr(Fermion(conserve = N)) == "Fermion(conserve = N)"
    @test repr(Boson(4, conserve = N)) == "Boson(4, conserve = N)"
    @test repr(Electron(conserve = (Ntot, 2Sz))) == "Electron(conserve = (2Sz, Ntot))"
    @test repr(Fermion(conserve = parity(N))) == "Fermion(conserve = parity(N))"
    # declaring a conservation makes a different site, and two declarations agree
    @test Fermion(conserve = N) == Fermion(conserve = N)
    @test Fermion(conserve = N) ≠ Fermion()
end

@testset "Printing a site" begin
    @test repr(Dummit()) == "Dummit()"
    @test TMS.conserved(Dummit()) == Conservation()
    # only the trailing fields carrying nothing are left out
    @test repr(Middling(nothing, 3)) == "Middling(nothing, 3)"
    @test repr(Middling(7, 3)) == "Middling(7, 3)"
    # a field named conserve that is not a Conservation is an ordinary field
    @test repr(Pretender("hello")) == "Pretender(\"hello\")"
    @test repr(Pretender("")) == "Pretender()"
    @test TMS.conserved(Pretender("N:0,1")) == Conservation()
end

@testset "Charged indices" begin
    hf = TMS.has_definite_flux
    i = siteind(Pure, System(1, Qubit(conserve = N)), 1)
    @test hf([1., 0.], [i])
    @test hf([0., 1.], [i])
    @test !hf([1., 1.] / √2, [i])
    @test hf([1e-15, 1.], [i])
    @test hf([1., 1.], [IT.Index(2)])
    sec = IT.space
    # a system whose sites declare nothing is dense
    s0 = System(3, Qubit())
    @test !IT.hasqns(siteind(Pure, s0, 1))
    @test !IT.hasqns(siteind(Mixed, s0, 1))
    @test dim(siteind(Pure, s0, 1)) == 2 && dim(siteind(Mixed, s0, 1)) == 4
    # the sectors are those the site recorded, one block per basis state
    s1 = System(2, Fermion(conserve = N))
    @test sec(siteind(Pure, s1, 1)) == [IT.QN("N", 0) => 1, IT.QN("N", 1) => 1]
    # the mixed index carries differences of charges: |m><n| has q(m) - q(n)
    @test sec(siteind(Mixed, s1, 1)) == [IT.QN("N", -1) => 1, IT.QN("N", 0) => 2, IT.QN("N", 1) => 1]
    # blocks are not merged, or a basis whose equal charges are not contiguous would be
    # reordered: parity on a boson gives 0, 1, 0, 1
    s2 = System(2, Boson(4, conserve = parity(N)))
    @test sec(siteind(Pure, s2, 1)) == [IT.QN("parity(N)", c, 2) => 1 for c in (0, 1, 0, 1)]
    s3 = System(2, Qudit(3, conserve = Z))
    @test sec(siteind(Pure, s3, 1)) == [IT.QN("Z", c, 3) => 1 for c in 0:2]
    s4 = System(2, Electron(conserve = (Ntot, 2Sz)))
    @test length(sec(siteind(Pure, s4, 1))) == 4
    @test dim(siteind(Pure, s4, 1)) == 4
    # one site declaring something makes every index charged, a site declaring nothing taking
    # a trivial charge
    s5 = System([Fermion(conserve = N), Qubit(), Fermion(conserve = N)])
    @test IT.hasqns(siteind(Pure, s5, 2))
    @test sec(siteind(Pure, s5, 2)) == [IT.QN() => 2]
    @test TMS.is_charged(s5)
    @test !TMS.is_charged(s0)
    # read once off the sites when the system is built, with the names conserved strongly
    @test TMS.strong_names(s5) == String[] && TMS.strong_names(s0) == String[]
    s6 = System([Electron(conserve = strong(Ntot)), Qubit(), Fermion(conserve = N)])
    @test TMS.is_charged(s6) && TMS.strong_names(s6) == ["Ntot"]
    # and its fermionic sites, the first from each site on
    @test [ TMS.next_fermion(s6, k) for k in 1:4 ] == [1, 3, 3, 4]
    @test TMS.next_fermion(System(3, Qubit()), 1) == 4
end

@testset "Operator tensors on a charged system" begin
    # the ITensor of an operator on a system is internal (decision 2): what 1.x read on its flux
    # is read on what the operator does. Rounding between charges, which an eigendecomposition
    # leaves behind, is not a charge, and the flux of an operator and its placement agree on it.
    # What counts as rounding is relative to the operator, so that a very small one is not taken
    # for the zero operator, which carries no charge
    b = Boson(5, conserve = parity(N))
    g = exp(0.2 * (A^2 + dag(A)^2) + 0.1 * N)
    @test flux(g, b) == IT.QN("parity(N)", 0, 2)
    function placed(site)
        s = State{Pure}(System(3, site), ["1", "2", "0"])
        return [expect(s, g(2)), norm(apply(g(2), s)), norm(apply(make_mpo(s, g(2)), s))]
    end
    @test placed(b) ≈ placed(Boson(5))
    @test_throws "no definite flux" flux(1e-20 * X, Qubit(conserve = N))
    # a site conserving nothing, between sites that do, takes a trivial charge, on which every
    # operator remains available
    ψ = State{Pure}(System([Fermion(conserve = N), Qubit(), Fermion(conserve = N)]),
                    ["Occ", "Up", "Emp"])
    for st in (ψ, mix(ψ))
        @test expect(apply(X(2), st), Z(2)) ≈ -1
        @test expect(apply(X(2), st), N(1)) ≈ 1
    end
    @test expect(apply(make_mpo(ψ, X(2)), ψ), Z(2)) ≈ -1
    # a channel has no flux: Gate(C) and Dissipator(C), which carry no charge, are applied to
    # charged sites in channels.jl
end

@testset "Weakening to the same symmetries" begin
    sys = System([Electron(conserve = (Ntot, 2Sz)), Electron(conserve = (2Sz, Ntot))])
    st = State{Pure}(sys, ["Up", "Dn"])
    t = symmetries(sys)
    @test weaken(sys, t) === sys
    @test weaken(st, t) === st
    @test_throws "Ntot is given twice" weaken(st, (Ntot, Ntot))
    @test_throws "Ntot is given twice" weaken(sys, (strong(Ntot), Ntot))
end

@testset "The order of the conserved quantities" begin
    @test Electron(conserve = (Ntot, 2Sz)) == Electron(conserve = (2Sz, Ntot))
    sys = System(3, Electron(conserve = (strong(Ntot), 2Sz)))
    @test repr(symmetries(sys)) == "(2Sz, strong(Ntot))"
    st = State{Mixed}(sys, "Up")
    for t in ((2Sz, strong(Ntot)), (strong(Ntot), 2Sz))
        @test weaken(sys, t) === sys
        @test weaken(st, t) === st
    end
end

@testset "Declaring a strong symmetry" begin
    mixed(s) = IT.space(siteind(Mixed, System(1, s), 1))
    # the strength is recorded on the quantity and not on the site
    @test charges(Fermion(conserve = strong(N))) == [("N", 1, [0, 1], true)]
    @test charges(Electron(conserve = (strong(Ntot), 2Sz))) ==
        [("2Sz", 1, [0, 1, -1, 0], false), ("Ntot", 1, [0, 1, 1, 2], true)]
    @test TMS.decode_conserve("N!:0,1") == TMS.conserved(Fermion(conserve = strong(N)))
    # a name ending in a character of several bytes, inside which the star or the ! was cut
    nu = Fermion(conserve = strong(named(N, "N↑")))
    @test TMS.decode_conserve("N↑!:0,1") == TMS.conserved(nu)
    @test TMS.decode_conserve("N↑%2:0,1").quantities[1].name == "N↑"
    ρu = mix(State{Pure}(System(2, nu), ["Occ", "Emp"]))
    @test expect(ρu, N(1)) ≈ 1
    @test trace(weaken(ρu)) ≈ 1
    @test repr(Fermion(conserve = strong(N))) == "Fermion(conserve = strong(N))"
    @test repr(Electron(conserve = (strong(Ntot), 2Sz))) == "Electron(conserve = (2Sz, strong(Ntot)))"
    # nothing changes on the pure side
    @test IT.space(siteind(Pure, System(1, Fermion(conserve = strong(N))), 1)) ==
        IT.space(siteind(Pure, System(1, Fermion(conserve = N)), 1))
    @test flux(N, Fermion(conserve = strong(N))) == IT.QN("N", 0)
    # keeping the two sides apart gives one block per pair of charges
    @test length(mixed(Fermion(conserve = N))) == 3
    @test length(mixed(Fermion(conserve = strong(N)))) == 4
    @test length(mixed(Boson(4, conserve = N))) == 7
    @test length(mixed(Boson(4, conserve = strong(N)))) == 16
    # a name ending in ! is an ordinary name now that the charges are typed
    @test_ok Fermion(conserve = strong(named(N, "N!")))
    # ITensors takes charge names of sixteen characters at most, the bra of a strong one
    # taking a star
    @test_ok Fermion(conserve = parity(dag(C) * C))
    @test_throws "at most 15 characters" Fermion(conserve = strong(parity(dag(C) * C)))
end

@testset "Charges that cannot live together" begin
    @test_throws "strongly on one site and weakly on another" System(
        [Fermion(conserve = N), Fermion(conserve = strong(N))])
    @test_ok System(2, Electron(conserve = (strong(Ntot), strong(2Sz))))
    @test_throws "of the four components" System(
        [Electron(conserve = (strong(Ntot), strong(2Sz))), Fermion(conserve = N)])
    @test_ok System(2, Fermion(conserve = strong(N)))
    @test_ok System([Fermion(conserve = strong(N)), Qubit(), Fermion(conserve = strong(N))])
end

@testset "Superoperators under a strong symmetry" begin
    ss = System(2, Fermion(conserve = strong(N)))
    ρ = mix(random_state(State{Pure}(ss, ["Occ", "Emp"]), 2))
    # a jump commuting with the charge is what a strong symmetry asks for: F on one site flips
    # the sign of a hopping that leaves it
    @test expect(apply(Gate(F)(1), ρ), dag(C)(1) * C(2)) ≈ -expect(ρ, dag(C)(1) * C(2))
    # one that moves the charge does not, and the refusal points at the way out
    @test_throws "conserving it strongly forbids: drop `strong`" apply(Dissipator(C)(1), ρ)
    # a jump of no definite charge at all is refused for what it lacks, under a weak symmetry as
    # under a strong one, where dropping `strong` would not help it
    for site in (Boson(3, conserve = N), Boson(3, conserve = strong(N)))
        ρb = mix(State{Pure}(System(2, site), "1"))
        @test_throws "no definite flux" apply(Gate(Q)(1), ρb)
    end
    # and so is an operator on several sites, which no one site check sees
    ψ = State{Pure}(System(2, Qubit(conserve = N)), "Up")
    @test_throws "carries no definite charge of N" make_mpo(ψ, exp(-0.3im * (X ⊗ X))(1, 2))
    # applied, its dense gate is refused in the same words, where the message of ITensors came
    # out, "Fluxes not all equal"
    for st in (ψ, mix(ψ))
        @test_throws "carries no definite charge of N" apply(exp(-0.3im * (X ⊗ X))(1, 2), st)
        @test_throws "no definite flux" apply((X ⊗ Id)(1, 2), st)
    end
    @test norm(apply(Swap(1, 2), ψ)) ≈ 1
end

@testset "Operators on several sites of a charged system" begin
    # an operator on several sites is laid on the index of each of them, and a channel on the
    # mixed index of each, which is where charged indices went wrong in 1.x. Each case is
    # compared with the same sites conserving nothing, on a quantity that reads a relative
    # phase or a population, which a tensor laid in a permuted basis would change
    same(f, charged, dense) = @test f(charged) ≈ f(dense)

    # a gate written as an expression, and one written as an exponential
    same(Qubit(conserve = N), Qubit()) do site
        sys = System(3, site)
        s = State{Pure}(sys, ["Up", "Dn", "Up"]) + 0.5 * State{Pure}(sys, ["Dn", "Up", "Up"])
        s = s / norm(s)
        a = apply(Swap(1, 2), s)
        b = apply(exp(-0.3im * (Z ⊗ Z))(2, 3), s)
        m = apply(Swap(2, 3), mix(s))
        [expect(a, N(1)), expect(b, Sp(1) * Sm(2) + Sm(1) * Sp(2)), expect(m, N(3))]
    end

    # Swap in expect and in MPOs, which place its factors one by one: written with X and Y, which
    # carry no charge of their own, it was refused on qubits that conserve something
    swaps = sum(Swap(i, i + 1) for i in 1:3)
    function swapped(site)
        sys = System(4, site)
        s = State{Pure}(sys, ["Up", "Dn", "Up", "Dn"]) +
            0.5 * State{Pure}(sys, ["Dn", "Up", "Up", "Dn"])
        return s / norm(s)
    end
    same(Qubit(conserve = N), Qubit()) do site
        s = swapped(site)
        [expect(s, Swap(2, 3)); inner(s.state', make_mpo(s, swaps), s.state);
         expect(s, controlled(Swap)(1, 2, 3)); trace(apply(make_mpo(mix(s), swaps), mix(s)))]
    end
    @test last(only(measure(swapped(Qubit(conserve = N)), Swap(1, 2)))) ≈
        last(only(measure(swapped(Qubit()), Swap(1, 2))))
    same(Qubit(conserve = strong(N)), Qubit()) do site
        ρ = mix(swapped(site))
        trace(apply(make_mpo(ρ, swaps), ρ))
    end
    # and its factors being real, so is its MPO, which Y made complex
    @test eltype(make_mpo(swapped(Qubit()), swaps)[2]) == Float64

    # a matrix obtained through an eigendecomposition, as the exponential of a hermitian one
    # is, keeps between charges the rounding the exact matrix does not have: the imaginary time
    # step of the Bose-Hubbard chain was refused on bosons conserving their number
    same(Boson(4, conserve = N), Boson(4)) do site
        s = State{Pure}(System(3, site), ["1", "2", "0"])
        g = exp(-0.1 * (A ⊗ dag(A) + dag(A) ⊗ A + 0.3 * N ⊗ N))(1, 2)
        a, m = apply(g, s), apply(g, mix(s))
        [norm(a); [expect(a, N(i)) for i in 1:3]; trace(m); [expect(m, N(i)) for i in 1:3]]
    end

    # charges that do not increase along the basis, with a channel and generators. The
    # dissipator is added to the state rather than applied alone, what it gives on its own
    # having no trace. Left and Right being gone (decision 1), Evolver, K ρ + ρ K†, acts on
    # both sides where they acted on one each
    same(Boson(4, conserve = parity(N)), Boson(4)) do site
        sys = System(3, site)
        p(v) = State{Pure}(sys, v)
        ρ = mix((p(["1", "2", "3"]) + 0.5 * p(["3", "0", "1"])) / sqrt(1.25))
        ρ = ρ + 0.1 * apply(Dissipator(0.3 * (A ⊗ A))(1, 2), ρ)
        ρ = apply(Gate(exp(-0.2im * (N ⊗ N)))(2, 3), ρ)
        ρ = apply(Evolver(exp(-0.1im * (N ⊗ N)))(1, 2), ρ)
        [trace(ρ); [expect(ρ, N(i)) for i in 1:3];
         expect(ρ, (dag(A) * dag(A))(1) * (A * A)(2) * (A * A)(3))]
    end
    same(Electron(conserve = (Ntot, 2Sz)), Electron()) do site
        sys = System(2, site)
        s = State{Pure}(sys, ["Up", "Dn"]) + 0.5 * State{Pure}(sys, ["Dn", "Up"])
        s = s / norm(s)
        g = exp(-0.3im * (Sz ⊗ Ntot + Ntot ⊗ Ntot))
        exchange = Sp(1) * Sm(2) + Sm(1) * Sp(2)
        [expect(apply(Gate(g)(1, 2), mix(s)), exchange), expect(apply(g(1, 2), s), exchange)]
    end

    # a site conserving nothing beside one that conserves, which used to put a plain index and
    # a charged one in the same tensor
    same([Qubit(), Fermion(conserve = N)], [Qubit(), Fermion()]) do sites
        s = State{Pure}(System(sites), ["Up", "Occ"])
        g = exp(-0.4im * (X ⊗ N))
        [expect(apply(Gate(g)(1, 2), mix(s)), Z(1)), expect(apply(g(1, 2), s), Y(1))]
    end

    same(Fermion(conserve = N), Fermion()) do site
        sys = System(3, site)
        p(v) = State{Pure}(sys, v)
        ρ = mix((p(["Occ", "Emp", "Occ"]) + 0.5 * p(["Emp", "Occ", "Occ"])) / sqrt(1.25))
        ρ = apply(Reset("Emp")(3), ρ)
        [expect(ρ, N(i)) for i in 1:3]
    end

    # under a strong symmetry, where the bra of each site lives apart from its ket
    same(Fermion(conserve = strong(N)), Fermion()) do site
        sys = System(3, site)
        p(v) = State{Pure}(sys, v)
        ρ = mix((p(["Occ", "Occ", "Emp"]) + 0.5 * p(["Occ", "Emp", "Occ"])) / sqrt(1.25))
        ρ = ρ + 0.1 * apply(Dissipator(0.3 * (N ⊗ N))(1, 2), ρ)
        ρ = apply(Gate(exp(-0.2im * (N ⊗ N)))(2, 3), ρ)
        [trace(ρ); [expect(ρ, N(i)) for i in 1:3]; expect(ρ, dag(C)(2) * C(3))]
    end

    # the refusal names the quantity even when the first site of the operator conserves nothing
    hyb = mix(State{Pure}(System([Qubit(), Fermion(conserve = strong(N))]), ["Up", "Occ"]))
    @test_throws "changes N by different amounts" apply(Dissipator(X ⊗ C)(1, 2), hyb)
    # `tensor` being ⊗ now, the ITensor of an operator given only its sites is internal, and
    # what 1.x read on it, one index for its sites and its flux, is not ported
end

@testset "A term vanishing on its site" begin
    # a product whose factor is zero on its site, as C(1)*C(1), is a term equal to zero, whether
    # a square vanishes depending on the site. On a charged system its tensor had no block,
    # hence no flux, and building the MPO failed on an ArgumentError from ITensors. Each
    # operator is compared with the same one written without that term, which must act alike
    # and take no more channels
    conf = ["Occ", "Emp", "Occ", "Emp"]
    for site in (Fermion(conserve = N), Fermion())
        s = random_state(State{Pure}(System(4, site), conf), 4)
        @test norm(make_mpo(s, C(1) * C(1))) == 0
        # and so does a sum whose terms all have a coefficient zero, which leaves none
        g = 0.0
        @test norm(make_mpo(s, g * C(1) + g * dag(C)(2))) == 0
        @test maxlinkdim(make_mpo(s, C(1) * C(1) * N(3) + N(1) * N(3))) ==
              maxlinkdim(make_mpo(s, N(1) * N(3)))
        # the approximations WI and WII are built from the same terms
        for w in (TMS.make_approx_W1, TMS.make_approx_W2)
            ev(a) = apply(w(PreMPO(s, -im * a), 0.1), s)
            @test norm(ev(C(1) * C(1) * N(3) + N(1) * N(3)) - ev(N(1) * N(3))) < 1e-12
        end
        # the square of a hopping hamiltonian, which the naive route to a variance forms, holds
        # such terms: each hopping term squared puts C twice on one site
        hop = sum(dag(C)(i) * C(i + 1) + dag(C)(i + 1) * C(i) for i in 1:3)
        @test inner(s.state', make_mpo(s, hop * hop), s.state) ≈ expect(s, hop * hop)
        # on a density matrix it is lifted with Evolver, H ρ + ρ H, whose trace is twice <H>
        @test trace(apply(make_mpo(mix(s), hop * hop), mix(s))) ≈ 2 * expect(s, hop * hop)
    end
    same(st, a, b) = @test norm(apply(make_mpo(st, a), st) - apply(make_mpo(st, b), st)) < 1e-12
    for site in (Fermion(conserve = N), Fermion())
        s = random_state(State{Pure}(System(4, site), conf), 4)
        same(s, C(1) * C(1) + N(2), N(2))
        same(s, C(1) * C(1) * N(3) + N(1) * N(3), N(1) * N(3))
        same(mix(s), C(1) * C(1) * N(3) + N(1) * N(3), N(1) * N(3))
    end
    s = random_state(State{Pure}(System(4, Fermion()), conf), 4)
    # what is dropped is what vanishes on the site at hand: the square of a fermionic operator
    # is not zero when it is a Majorana one, and the same Sp squares to zero on a spin 1/2 and
    # not on a spin 1
    G = named(C + dag(C), "G")
    same(s, G(1) * G(1) * N(3), N(3))
    for site in (Spin(1/2, conserve = 2Sz), Spin(1/2))
        @test norm(make_mpo(State{Pure}(System(3, site), "1/2"), Sp(1) * Sp(1) * Sz(2))) == 0
    end
    for site in (Spin(1, conserve = 2Sz), Spin(1))
        @test norm(make_mpo(State{Pure}(System(3, site), "0"), Sp(1) * Sp(1) * Sz(2))) > 0
    end
end

@testset "Measuring on every kind of charged site" begin
    # conserving a quantity must not change a single measured number, whatever the site and
    # the shape of the charge: an index keeps its blocks in the order of the basis, one per
    # state, while a combiner sorts and merges them, the two agreeing only when the charges of
    # the basis increase and never repeat, which they do not on electrons, t-J sites, or under
    # parity
    cases = ((Qubit(conserve = N),              Qubit(),        "Dn",  N),
             (Qubit(conserve = 2Sz),            Qubit(),        "Dn",  Z),
             (Qubit(conserve = parity(N)),      Qubit(),        "Dn",  N),
             (Spin(3/2, conserve = 2Sz),        Spin(3/2),      "1/2", Sz),
             (Boson(4, conserve = N),           Boson(4),       "2",   N),
             (Boson(4, conserve = parity(N)),   Boson(4),       "2",   N),
             (Fermion(conserve = N),            Fermion(),      "Occ", N),
             (Electron(conserve = Ntot),        Electron(),     "Up",  Ntot),
             (Electron(conserve = (Ntot, 2Sz)), Electron(),     "Up",  Ntot),
             (Tj(conserve = (Ntot, 2Sz)),       Tj(),           "Up",  Ntot),
             (Qboson(0.1, 4, conserve = N),     Qboson(0.1, 4), "2",   N),
             (Qudit(4, conserve = Z),           Qudit(4),       "1",   Z))    # Zd is Z now
    for (charged, dense, st, op) in cases, rep in (identity, mix)
        state(site) = rep(State{Pure}(System(3, site), fill(st, 3)))
        @test expect(state(charged), op(2)) ≈ expect(state(dense), op(2))
    end
    # an operator of two sites applied as a gate is laid on the indices of both at once, and
    # the norm of the result reads off the eigenvalue it picks, which a permuted basis would
    # change
    for (charged, dense, st, count) in ((Fermion(conserve = N), Fermion(), "Occ", N),
                                        (Boson(4, conserve = parity(N)), Boson(4), "2", N),
                                        (Electron(conserve = (Ntot, 2Sz)), Electron(), "Up", Ntot))
        for op in (count ⊗ count, Id ⊗ count)
            two(site) = State{Pure}(System(3, site), fill(st, 3))
            @test norm(apply(op(1, 2), two(charged))) ≈ norm(apply(op(1, 2), two(dense)))
        end
    end
end
