# Site types: the declarations of a site type of one's own, the matrices of the operators of
# the built in types, their local states, what they conserve and the indices they give.
#
# An assertion that needs a later layer is kept as `@test_skip`, marked with that layer; the
# testsets of later layers are still in v1/sites.jl.

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

@testset "The spin of a qubit" begin
    @test norm(matrix(Sx^2 + Sy^2 + Sz^2 - S2, Qubit())) ≈ 0 atol = 1e-12
    @test_skip "the measurements of Qubit measuring"        # couche 5
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
    @test_skip "the measurements of Spin measuring"         # couche 5
end

@testset "One electron of fully mixed spin" begin
    for site in (Electron(conserve = strong(Ntot)), Tj(conserve = strong(Ntot)))
        st = State{Mixed}(System(2, site), ["MixedSpin", "↑|↓"])
        @test trace(st) ≈ 1
        @test_skip real(expect(st, Ntot(1))) ≈ 1          # couche 2
        @test_skip real(trace2(st)) ≈ 1 / 4                # couche 3
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
    @test matrix(Sumd(2), Qudit(2), Qudit(2)) ≈ matrix(controlled(X), Qubit(), Qubit())
    # Sum adds the level of the first qudit to the second, modulo d
    m = matrix(Sumd(3), Qudit(3))
    for i in 0:2, j in 0:2
        @test m[3i + mod(i + j, 3) + 1, 3i + j + 1] ≈ 1
    end
    @test_skip "the states Sumd and Z act on, through apply and expect"     # couche 2
end

@testset "Unknown operator and state names" begin
    @test_throws "operator X is not defined for site Boson" matrix(X, Boson(3))
    @test_throws "operator A is not defined for site Fermion" matrix(A, Fermion())
    @test_throws "state Zorglub is not defined for site Qubit" State{Pure}(System(2, Qubit()), "Zorglub")
    @test_throws "state Zorglub is not defined for site Dummit" State{Pure}(System(2, Dummit()), "Zorglub")
    @test_throws "a state of 3 components cannot be one of Qubit()" State{Pure}(System(2, Qubit()), [1., 0., 0.])
    @test_throws "a 3×3 density matrix cannot be one of Qubit()" State{Mixed}(
        System(2, Qubit()), [1. 0. 0. ; 0. 0. 0. ; 0. 0. 0.])
    @test_skip "a Reset to a state of 3 components on a qubit"                  # couche 2
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
    @test_skip "the measurements of N on Dummit and Boson"                      # couche 2
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
    @test_skip "an MPO of X ⊗ X on Qudit(4, conserve = Z^2)"                # couche 2
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
    @test_skip "a pairing term in an MPO under named(parity(N))"              # couche 2
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
