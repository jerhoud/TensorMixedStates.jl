# Defining a site type from outside the package, which is what a user does and what the
# built in types cannot check: they are declared inside TensorMixedStates, where dim is
# local and the operator names are not imported. The site below is deliberately as
# trivial as possible; everything else about sites is covered by the seven real types.

struct Dummit <: AbstractSite end

TensorMixedStates.dim(::Dummit) = 2

@def_states(Dummit(), [ "1" => [0., 1.] ])

# `N` reaches this module through the `using` of runtests.jl, which is the interesting case:
# `@def_operators` must register N for Dummit without binding the name again
@def_operators(Dummit(), [ selfadjoint_op => [ N = [0. 0. ; 0. 1.] ] ])

@create_site_module(Dummits, [Dummit, N])

# a second site, to check that the same name declared with another OpType is refused
struct Dummit2 <: AbstractSite end

TensorMixedStates.dim(::Dummit2) = 2

# a site carrying nothing in a field that is not the last one, to check that printing does
# not drop it: the call would no longer line up with the fields
struct Middling <: AbstractSite
    a::Union{Nothing, Int}
    b::Int
end

TensorMixedStates.dim(::Middling) = 2

# a site whose field happens to be named conserve without being one. The field being
# optional, nothing declares this wrong, and printing must not raise on it
struct Pretender <: AbstractSite
    conserve::String
end

TensorMixedStates.dim(::Pretender) = 2

@testset "Qubit measuring" begin
    @test_pm test_phases(CreateState{type}(1, Qubit(), "Z+"; 
        final_measures = check([X(1), Y(1), Z(1)], [0, 0, 1])))
    @test_pm test_phases(CreateState{type}(6, Qubit(), ["X+", "Y+", "Z+", "X-", "Y-", "Z-"];
        final_measures = check([X, Y, Z], [[1, 0, 0, -1, 0, 0], [0, 1, 0, 0, -1, 0], [0, 0, 1, 0, 0, -1]])))
    @test_pm test_phases(CreateState{type}(3, Qubit(), ["X+", "Y-", "Z-"];
        final_measures = check([(X, Y), (Z, Z)], [[0 -1 0; 0 0 0; 0 0 -im], [1 0 0; 0 1 0; 0 0 1]])))
    @test_pm test_phases(CreateState{type}(4, Qubit(), ["Z+", "X-", "Z-", "Y-"];
        final_measures = check([Z(1)Y(4), X(2)Z(3), Z(1)Z(3), X(2)Y(4)], [-1, 1, -1, 1])))
    @test norm(matrix(Sx^2+Sy^2+Sz^2-S2, Qubit()))≈0 atol=1e-12
end

@testset "Qubit excitation number" begin
    q = Qubit()
    # N counts "Dn" as the occupied state, so that it agrees with Sz rather than against it
    @test matrix(N, q) ≈ [0. 0. ; 0. 1.]
    @test matrix(N, q) ≈ matrix(Proj(1), q)
    @test matrix(N, q) ≈ matrix((Id - Z) / 2, q)
    @test matrix(N, q) ≈ matrix(0.5 * Id - Sz, q)

    # the name is shared with the sites that already had it, as the whole scheme intends
    @test Qubits.N === Bosons.N === Fermions.N

    # conserving it is the same conservation as conserving Sz, said in the language of
    # excitations and with integer charges rather than the doubled ones of a half integer spin
    @test Qubit(conserve = N).conserve == "N:0,1"
    @test Qubit(conserve = 2Sz).conserve == "2Sz:1,-1"

    # "Up" being empty, it is Sm that creates an excitation
    @test flux(Sm, Qubit(conserve = N)) == TensorMixedStates.ITensors.QN("N", 1)
    @test flux(Sp, Qubit(conserve = N)) == TensorMixedStates.ITensors.QN("N", -1)
    @test_throws "no definite flux" flux(X, Qubit(conserve = N))
end

@testset "Spin excitation number" begin
    # N counts the excitations above the state of maximal Sz, which is the Holstein-Primakoff
    # counting and the same convention as for a qubit
    for x in (1/2, 1, 3/2, 2)
        site = Spin(x)
        @test matrix(N, site) ≈ matrix(x * Id - Sz, site)
        # integer for every spin, half integer ones included, which is what spares the
        # doubling that 2Sz needs
        @test all(k -> matrix(N, site)[k, k] ≈ k - 1, 1:dim(site))
    end

    # a qubit and a spin one half are the same system, and now say so
    @test matrix(N, Spin(1/2)) ≈ matrix(N, Qubit())
    @test Spins.N === Qubits.N === Bosons.N

    @test Spin(1/2, conserve = N).conserve == "N:0,1"
    @test Spin(3/2, conserve = N).conserve == "N:0,1,2,3"
    # the same conservation as 2Sz, counted from the other end
    @test Spin(3/2, conserve = 2Sz).conserve == "2Sz:3,1,-1,-3"

    @test flux(Sm, Spin(1, conserve = N)) == TensorMixedStates.ITensors.QN("N", 1)
    @test flux(Sp, Spin(1, conserve = N)) == TensorMixedStates.ITensors.QN("N", -1)
end

@testset "Fermion measuring" begin
    @test_pm test_phases(CreateState{type}(1, Fermion(), "1";
        final_measures = check(N(1), 1)))
    @test_pm test_phases(CreateState{type}(2, Fermion(), ["0", "1"];
        final_measures = check(N, [0, 1])))
end

@testset "Boson measuring" begin
    @test_pm test_phases(CreateState{type}(4, Boson(4), ["0", "1", "2", "3"];
        final_measures = check(N, [0, 1, 2, 3])))
end

@testset "Spin measuring" begin
    @test_pm test_phases(CreateState{type}(4, Spin(3/2), ["-3/2", "-1/2", "1/2", "3/2"];
        final_measures = check([Sx, Sy, Sz], [[0, 0, 0, 0], [0, 0, 0, 0], [-3/2, -1/2, 1/2, 3/2]])))
    @test norm(matrix(Sx^2+Sy^2+Sz^2-S2,Spin(5/2)))≈0 atol=1e-12
    @test norm(matrix(Sx^2+Sy^2+Sz^2-S2,Spin(4)))≈0 atol=1e-12
    @test_pm test_phases(CreateState{type}(9, Spin(1),
            ["X-1", "X0", "X1", "Y-1", "Y0", "Y1", "Z-1", "Z0", "Z1"];
        final_measures = check([Sx, Sy, Sz],
            [[-1, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, -1, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, -1, 0, 1]])))
    @test_pm test_phases(CreateState{type}(12, Spin(3/2),
            ["X-3/2", "X-1/2", "X1/2", "X3/2", "Y-3/2", "Y-1/2", "Y1/2", "Y3/2", "Z-3/2", "Z-1/2", "Z1/2", "Z3/2"];
        final_measures = check([Sx, Sy, Sz],
            [[-3/2, -1/2, 1/2, 3/2, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, -3/2, -1/2, 1/2, 3/2, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, -3/2, -1/2, 1/2, 3/2]])))
end

@testset "Electron measuring" begin
    @test_pm test_phases(CreateState{type}(4, Electron(), ["Emp", "Up", "Dn", "UpDn"];
        final_measures = check([Nup, Ndn, Nupdn, Ntot], [[0, 1, 0, 1], [0, 0, 1, 1], [0, 0, 0, 1], [0, 1, 1, 2]])))        
end

@testset "Tj measuring" begin
    @test_pm test_phases(CreateState{type}(3, Tj(), ["Emp", "Up", "Dn"];
        final_measures = check([Nup, Ndn, Ntot], [[0, 1, 0], [0, 0, 1], [0, 1, 1]])))        
end

@testset "Qboson measuring" begin
    @test_pm test_phases(CreateState{type}(4, Qboson(0.1, 4), ["0", "1", "2", "3"];
        final_measures = check(N, [0, 1, 2, 3])))
end

@testset "State mixing" begin
    @test_ok check_mix(
        [1, 5, 15],
        [Qubit(), Fermion(), Boson(3), Spin(3/2), Electron(), Tj()],
        [[X, Y, Z], [N], [N], [Sx, Sy, Sz], [Nup, Ndn, Nupdn, Ntot, Sx, Sy, Sz], [Nup, Ndn, Ntot, Sx, Sy, Sz]])
end

@testset "Mix ordering" begin
    sys = System(2, Qubit())
    ux2 = Id ⊗ X
    psi = State{Pure}(sys, "0")
    rho_from_pure = mix(apply(ux2(1, 2), psi))
    rho_direct = apply(ux2(1, 2), mix(psi))

    @test expect(rho_from_pure, Z(1)) ≈ 1
    @test expect(rho_direct, Z(1)) ≈ 1
    @test expect(rho_from_pure, Z(2)) ≈ -1
    @test expect(rho_direct, Z(2)) ≈ -1
end

@testset "Qudit measuring" begin
    # the clock and shift operators are unitary, of order d, and obey the Weyl relation
    for d in (2, 3, 5)
        s = Qudit(d)
        ω = exp(2im * π / d)
        z = matrix(Zd, s)
        x = matrix(Xd, s)
        id = matrix(Id, s)
        @test z^d ≈ id
        @test x^d ≈ id
        @test z * z' ≈ id
        @test x * x' ≈ id
        @test z * x ≈ ω * (x * z)
        # the level operator, and the Fourier operator that exchanges the two bases
        @test matrix(N, s) ≈ [ i == j ? i - 1. : 0. for i in 1:d, j in 1:d ]
        h = matrix(Hd, s)
        @test h * h' ≈ id
        @test h^4 ≈ id
        @test h * x * inv(h) ≈ z
        @test h * z * inv(h) ≈ inv(x)
        # the phase operator is Clifford: it maps Xd onto Xd Zd up to a phase
        p = matrix(S, s)
        @test p * p' ≈ id
        a = p * x * inv(p)
        b = x * z
        @test a ≈ (a[findfirst(!=(0), a)] / b[findfirst(!=(0), a)]) * b
        # its order is d for odd d and 2d for even d, which is intrinsic
        @test p^(isodd(d) ? d : 2d) ≈ id
    end
    # a qudit of dimension 2 is a qubit
    @test matrix(Zd, Qudit(2)) ≈ matrix(Z, Qubit())
    @test matrix(Xd, Qudit(2)) ≈ matrix(X, Qubit())
    @test matrix(Hd, Qudit(2)) ≈ matrix(H, Qubit())
    @test matrix(S, Qudit(2)) ≈ matrix(Qubits.S, Qubit())
    @test matrix(Sumd(2), Qudit(2), Qudit(2)) ≈ matrix(controlled(X), Qubit(), Qubit())
    # Sum adds the level of the first qudit to the second, modulo d
    st = State{Pure}(System(3, Qudit(3)), ["1", "1", "0"])
    @test real(expect1(st, N)) ≈ [1, 1, 0]
    @test real(expect1(apply(Sumd(3)(1, 2), st), N)) ≈ [1, 2, 0]
    @test real(expect1(apply(Sumd(3)(1, 3), st), N)) ≈ [1, 1, 1]
    # being an expression rather than a matrix, Sum also goes into an MPO
    @test maxlinkdim(make_mpo(st, Sumd(3)(1, 2))) == 4
    # the states are named by their level, and Zd reads that level back as a phase
    sys = System(3, Qudit(3))
    @test expect1(State{Pure}(sys, ["0", "1", "2"]), Zd) ≈ [exp(2im * π * n / 3) for n in 0:2]
    @test_pm test_phases(CreateState{type}(3, Qudit(3), ["0", "1", "2"];
        final_measures = check(Zd, [exp(2im * π * n / 3) for n in 0:2])))
end

@testset "Unknown operator and state names" begin
    # a mistyped name is the most ordinary mistake there is with these libraries, so it
    # has to name what was not found rather than surface as a dictionary error. X is a
    # perfectly good name, just not one a boson has
    @test_throws "operator X is not defined for site Boson" matrix(X, Boson(3))
    @test_throws "state Zorglub is not defined for site Qubit" State{Pure}(
        System(2, Qubit()), "Zorglub")
    # a site type declared outside the package goes through the same library
    @test_throws "state Zorglub is not defined for site Dummit" State{Pure}(
        System(2, Dummit()), "Zorglub")
end

@testset "Custom site type" begin
    @test dim(Dummit()) == 2
    # N arrives here by `using` from another site module: declaring it again for a new site
    # must neither fail nor disturb the site it came from, and the name must go on standing
    # for one and the same operator
    @test real(expect1(State{Pure}(System(2, Dummit()), "1"), N)) ≈ [1, 1]
    @test real(expect1(State{Pure}(System(2, Boson(4)), ["1", "3"]), N)) ≈ [1, 3]
    @test Dummits.Dummit === Dummit
    @test Dummits.N === N
    @test N === Bosons.N
    # the same name with another OpType would change the meaning of N for every site
    # already using it, so it is refused, and the refusal leaves the library untouched
    @test_throws "must agree on the OpType" @def_operators(Dummit2(),
        [ plain_op => [ N = [0. 0. ; 0. 1.] ] ])
    @test_throws "operator N is not defined for site Dummit2" matrix(N, Dummit2())
end

@testset "Index tags" begin
    # the tag says which site type the index belongs to, and must not depend on which site
    # modules the user has imported: `string(typeof(site))` printed the module prefix when
    # the module was not in scope, and ITensors cuts a tag at 16 characters, so every site
    # type came out tagged "TensorMixedState"
    for site in [Qubit(), Fermion(), Boson(4), Spin(3/2), Electron(), Tj(), Qboson(0.1, 3),
                 Qudit(3), Dummit()]
        # `hastags` belongs to ITensors, which the test environment does not import
        @test TensorMixedStates.ITensors.hastags(Index(site), string(nameof(typeof(site))))
        @test TensorMixedStates.ITensors.hastags(Index(site), "Site")
    end
end

@testset "Declaring a conserved quantity" begin
    # a site records the name of each conserved quantity, its modulus when that is not 1,
    # and the charge of every basis state: an operator cannot be written to a state file,
    # and the charges are all that is needed afterwards
    @test Fermion(conserve = N).conserve == "N:0,1"
    @test Fermion(conserve = parity(N)).conserve == "parity(N)%2:0,1"
    @test Fermion(conserve = named(N, "Nf")).conserve == "Nf:0,1"
    @test Boson(6, conserve = mod(N, 3)).conserve == "mod(N,3)%3:0,1,2,0,1,2"
    @test Spin(1, conserve = Sz).conserve == "Sz:1,0,-1"
    @test Spin(1/2, conserve = 2Sz).conserve == "2Sz:1,-1"
    @test Electron(conserve = (Ntot, 2Sz)).conserve == "Ntot:0,1,1,2;2Sz:0,1,-1,0"
    @test Qubit().conserve == ""

    # the modulus of Zd is read off its spectrum, the eigenvalues being genuine roots of
    # unity, so a clock operator needs no modulus written anywhere
    for d in (3, 5, 8)
        @test Qudit(d, conserve = Zd).conserve == "Zd%$d:" * join(0:d-1, ",")
    end

    # what cannot be a charge, and the message says how far it is from being one
    @test_throws "is not diagonal" Spin(1, conserve = Sx)
    @test_throws "2Sz rather than Sz" Spin(1/2, conserve = Sz)
    @test_throws "already carries a charge modulo" Qudit(3, conserve = mod(Zd, 2))

    # a conserved quantity is carried by one site. Without this the multi site operator went
    # into matrix, which expanded its definition and complained about an operator of another
    # site type, naming nothing the caller had written
    @test_throws "acts on one site" Fermion(conserve = Swap)
    @test_throws "acts on one site" Qubit(conserve = mod(Swap, 2))

    # the four cases, read directly
    @test TensorMixedStates.site_charges(N, Fermion()) == (1, [0, 1])
    @test TensorMixedStates.site_charges(parity(N), Boson(4)) == (2, [0, 1, 0, 1])
    @test TensorMixedStates.site_charges(Zd, Qudit(3)) == (3, [0, 1, 2])

    # the recorded form is read back as it was written
    c = Electron(conserve = (Ntot, 2Sz)).conserve
    @test TensorMixedStates.decode_conserve(c) ==
        [("Ntot", 1, [0, 1, 1, 2], false), ("2Sz", 1, [0, 1, -1, 0], false)]
    @test TensorMixedStates.decode_conserve("") == Tuple{String, Int, Vector{Int}, Bool}[]

    # a site conserving nothing prints as it always did, and one that conserves prints
    # under the name of its charges rather than under the charges themselves
    @test repr(Qubit()) == "Qubit()"
    @test repr(Boson(4)) == "Boson(4)"
    @test repr(Spin(3/2)) == "Spin(1.5)"
    @test repr(Fermion(conserve = N)) == "Fermion(conserve = N)"
    @test repr(Boson(4, conserve = N)) == "Boson(4, conserve = N)"
    @test repr(Electron(conserve = (Ntot, 2Sz))) == "Electron(conserve = (Ntot, 2Sz))"
    @test repr(Fermion(conserve = parity(N))) == "Fermion(conserve = parity(N))"

    # declaring a conservation makes a different site, and two declarations agree
    @test Fermion(conserve = N) == Fermion(conserve = N)
    @test Fermion(conserve = N) ≠ Fermion()
end


@testset "Printing a site" begin
    # a site declaring no conservation prints as it always did, and one declaring none
    # either because it has no such field at all
    @test repr(Dummit()) == "Dummit()"
    @test TensorMixedStates.conserved(Dummit()) == ""
    @test TensorMixedStates.conserved(Fermion(conserve = N)) == "N:0,1"

    # only the trailing fields carrying nothing are left out
    @test repr(Middling(nothing, 3)) == "Middling(nothing, 3)"
    @test repr(Middling(7, 3)) == "Middling(7, 3)"

    # show has to print something whatever a site put in a field named conserve, an error
    # raised while printing being far worse than an odd looking site
    @test repr(Pretender("hello")) == "Pretender(conserve = hello)"
    @test repr(Pretender("")) == "Pretender()"
end

@testset "Charged indices" begin
    Q = TensorMixedStates.ITensors
    sec(i) = Q.space(i)
    pure(sys, k) = TensorMixedStates.SysIndex{Pure}(sys, k)
    mixed(sys, k) = TensorMixedStates.SysIndex{Mixed}(sys, k)

    # a system whose sites declare nothing is dense, exactly as before
    s0 = System(3, Qubit())
    @test !Q.hasqns(pure(s0, 1))
    @test !Q.hasqns(mixed(s0, 1))
    @test dim(pure(s0, 1)) == 2 && dim(mixed(s0, 1)) == 4

    # the sectors are those the site recorded, one block per basis state
    s1 = System(2, Fermion(conserve = N))
    @test sec(pure(s1, 1)) == [Q.QN("N", 0) => 1, Q.QN("N", 1) => 1]

    # the mixed index carries differences of charges, not sums: |m><n| has q(m) - q(n)
    @test sec(mixed(s1, 1)) == [Q.QN("N", -1) => 1, Q.QN("N", 0) => 2, Q.QN("N", 1) => 1]

    # blocks are not merged, or a basis whose equal charges are not contiguous would be
    # reordered: parity on a boson gives 0, 1, 0, 1
    s2 = System(2, Boson(4, conserve = parity(N)))
    @test sec(pure(s2, 1)) == [Q.QN("parity(N)", c, 2) => 1 for c in (0, 1, 0, 1)]

    # a modulus read off the spectrum, and charges with two components
    s3 = System(2, Qudit(3, conserve = Zd))
    @test sec(pure(s3, 1)) == [Q.QN("Zd", c, 3) => 1 for c in 0:2]
    s4 = System(2, Electron(conserve = (Ntot, 2Sz)))
    @test length(sec(pure(s4, 1))) == 4
    @test dim(pure(s4, 1)) == 4

    # the mode is a property of the whole list: one site declaring something makes every
    # index charged, a site declaring nothing taking a trivial charge rather than staying
    # dense, since an MPS cannot mix the two kinds
    s5 = System([Fermion(conserve = N), Qubit(), Fermion(conserve = N)])
    @test Q.hasqns(pure(s5, 2))
    @test sec(pure(s5, 2)) == [Q.QN() => 2]
    @test TensorMixedStates.is_charged(s5.sites)
    @test !TensorMixedStates.is_charged(s0.sites)

    # a site on its own says what it declares and nothing more
    @test !Q.hasqns(Index(Qubit()))
    @test Q.hasqns(Index(Fermion(conserve = N)))
end

@testset "Operator tensors on a charged system" begin
    Q = TensorMixedStates.ITensors
    ten(sys, op) = TensorMixedStates.tensor(sys, op)

    # the flux of an operator placed on a charged system is the charge it carries
    sys = System(3, Fermion(conserve = N))
    @test flux(ten(sys, N(1))) == Q.QN("N", 0)
    @test flux(ten(sys, Left(N)(1))) == Q.QN("N", 0)
    @test flux(ten(sys, Left(C)(1))) == Q.QN("N", -1)
    @test flux(ten(sys, Right(C)(1))) == Q.QN("N", 1)
    @test flux(ten(sys, Gate(C)(1))) == Q.QN("N", 0)
    @test flux(ten(sys, Dissipator(C)(1))) == Q.QN("N", 0)

    # a site conserving nothing, inside a system where another one does, takes a trivial
    # index, on which the whole vocabulary of operators remains available
    mixed = System([Fermion(conserve = N), Qubit(), Fermion(conserve = N)])
    @test flux(ten(mixed, X(2))) == Q.QN()
    @test flux(ten(mixed, Left(X)(2))) == Q.QN()
    @test flux(ten(mixed, N(1))) == Q.QN("N", 0)

    # an operator carrying no flux is refused by a message naming it, rather than by the
    # `Fluxes not all equal` of ITensors, raised where neither operator nor site is in sight
    q = System(2, Qubit(conserve = 2Sz))
    @test_throws "X on site Qubit connects charges" ten(q, X(1))
    @test_throws "no definite flux" ten(q, Left(X)(1))
end

@testset "Declaring a strong symmetry" begin
    Q = TensorMixedStates.ITensors
    strong = TensorMixedStates.strong
    mixed(s) = Q.space(mix(Index(s), s))

    # the strength is recorded on the quantity and not on the site, so one site may hold
    # both kinds, and it travels in the string the site keeps
    @test Fermion(conserve = strong(N)).conserve == "N!:0,1"
    @test Electron(conserve = (strong(Ntot), 2Sz)).conserve ==
        "Ntot!:0,1,1,2;2Sz:0,1,-1,0"
    @test TensorMixedStates.decode_conserve("N!:0,1") == [("N", 1, [0, 1], true)]

    # and a site goes on printing as the call that built it
    @test repr(Fermion(conserve = strong(N))) == "Fermion(conserve = strong(N))"
    @test repr(Electron(conserve = (strong(Ntot), 2Sz))) ==
        "Electron(conserve = (strong(Ntot), 2Sz))"

    # nothing changes on the pure side: only the bra of the mixed index takes another name
    @test Q.space(Index(Fermion(conserve = strong(N)))) ==
        Q.space(Index(Fermion(conserve = N)))
    @test flux(N, Fermion(conserve = strong(N))) == Q.QN("N", 0)

    # keeping the two sides apart gives one block per pair of charges instead of one per
    # difference, which is the whole gain of a strong symmetry
    @test length(mixed(Fermion(conserve = N))) == 3
    @test length(mixed(Fermion(conserve = strong(N)))) == 4
    @test length(mixed(Boson(4, conserve = N))) == 7
    @test length(mixed(Boson(4, conserve = strong(N)))) == 16

    # a name ending in ! could not be told from the mark a site puts on a strong symmetry
    @test_throws "cannot be told from the mark" Fermion(conserve = named(N, "N!"))
end

@testset "Charges that cannot live together" begin
    strong = TensorMixedStates.strong

    # the same name would stand for the charge of the ket on one site and for a difference
    # on the other, and the flux of a state would add the two
    @test_throws "strongly on one site and weakly on another" System(
        [Fermion(conserve = N), Fermion(conserve = strong(N))])

    # a strong quantity takes two of the four components ITensors allows, so two of them
    # fit exactly and anything more does not
    @test_ok System(2, Electron(conserve = (strong(Ntot), strong(2Sz))))
    @test_throws "of the four components" System(
        [Electron(conserve = (strong(Ntot), strong(2Sz))), Fermion(conserve = N)])

    # while sites that agree, or say nothing, live together
    @test_ok System(2, Fermion(conserve = strong(N)))
    @test_ok System([Fermion(conserve = strong(N)), Qubit(), Fermion(conserve = strong(N))])
end

@testset "Superoperators under a strong symmetry" begin
    Q = TensorMixedStates.ITensors
    strong = TensorMixedStates.strong
    ten(sys, op) = TensorMixedStates.tensor(sys, op)
    sys = System(2, Fermion(conserve = strong(N)))

    # a jump commuting with the charge is what a strong symmetry asks for, and it passes
    @test flux(ten(sys, Dissipator(N)(1))) == Q.QN(("N", 0), ("N*", 0))
    @test flux(ten(sys, Gate(F)(1))) == Q.QN(("N", 0), ("N*", 0))

    # one that moves the charge does not, and the refusal names it and the way out rather
    # than leaving the `Fluxes not all equal` of ITensors through
    @test_throws "changes N between its two sides" ten(sys, Dissipator(C)(1))
    @test_throws "drop `strong`" ten(sys, Dissipator(C)(1))

    # the very same jump is fine when the quantity is conserved weakly, which is the whole
    # difference between the two: a weak symmetry lets the charge move, a strong one does not
    @test flux(ten(System(2, Fermion(conserve = N)), Dissipator(C)(1))) == Q.QN("N", 0)
end
