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
