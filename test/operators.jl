# The algebra of operators: building, comparing, ordering and printing them, their matrix on
# sites, the Jordan-Wigner strings of a tensor product included, the split of an operator of
# several sites into factors of one site, and the descent of an operator onto a system.
#
# An assertion that needs a later layer is kept as `@test_skip`, marked with that layer, those
# of `simplify`, internal from layer 2 on (decision 2), among them. The testsets of `simplify`,
# of the superoperators, of the lowering onto a system and of the classification by OpType,
# which decision 2 removes, are still in v1/operators.jl; their assertions of the first layer
# are here.

@testset "The adjoint is an involution" begin
    # which only shows on operators that are not self adjoint
    q = Qubit()
    for a in [S, T, S * T, 2im * S, X + im * Y, exp(im * X)]
        @test matrix(dag(dag(a)), q) ≈ matrix(a, q)
    end
    @test matrix(dag(dag(dag(S))), q) ≈ matrix(dag(S), q)
end

@testset "The zero operator" begin
    # a term of coefficient zero, and a sum left with no term at all, as an empty one, are the
    # zero operator
    q, f = Qubit(), Fermion()
    @test matrix(X - X, q) ≈ zeros(2, 2)
    @test matrix(0 * X, q) ≈ zeros(2, 2)
    g = 0.0
    @test matrix(g * C + g * dag(C), f) == zeros(2, 2)
    @test matrix(TMS.SumOp(TMS.Op{IsOperator, TMS.Generic, 1}[]), f) == zeros(2, 2)
end

@testset "A matrix is taken on one site at least" begin
    # which every method of matrix asks for: given none, the refusal came from the factor the
    # recursion reached, X for 2X, or as a MethodError on F
    for a in (X, 2X, exp(X), X * Y, X ⊗ Y, F, 2F, F * X, Id, Proj("Up"), Proj(X => 1))
        @test_throws MethodError matrix(a)
    end
    @test_throws "X acts on 1 site and was given 2" matrix(X, Qubit(), Qubit())
end

@testset "Global ordering of operators" begin
    # `isless(::Op, ::Op)` ranks the types first and, for two of the same rank, asks `isless`
    # again on the pair. A type with a ranking but no `isless` of its own therefore recurses on
    # itself instead of comparing anything
    for (a, b) in [(exp(X), exp(Y)), (dag(X), dag(Y)), (parity(X), parity(Y)), (X^2, Y^2),
                   (X^0.5, Y^0.5), (X ⊗ Y, Y ⊗ X), (X + Y, Y + Z), (X * Y, Y * Z), (2X, 2Y),
                   (X(1), X(2)), (Proj("Dn"), Proj("Up"))]
        @test isless(a, b) && !isless(b, a)
    end
    # a Proj holds a name, a vector or an eigenspace, which have no order between them
    # (decision 20: a number is not a state)
    @test length(sort([Proj("1"), Proj("Up"), Proj([1., 0.]), Proj("0")])) == 4
    @test_skip matrix(simplify(Proj("1") + Proj("Up")), Qubit()) ≈ matrix(Id, Qubit())   # couche 2
    for (a, b) in [(Gate(X), Gate(Y)), (Dephase(X), Dephase(Y)), (Dissipator(X), Dissipator(Y)),
                   (Evolver(X(1)), Evolver(Y(1))), (Reset("Dn"), Reset("Up"))]
        @test isless(a, b) && !isless(b, a)
    end
    # a Reset holds a name, a vector, a matrix or an eigenspace, which have no order between them
    @test length(sort([Reset("Up"), Reset([1., 0.]), Reset("Dn"), Reset([1. 0. ; 0. 0.]),
                       Reset(Level => 1)])) == 5
    # two of the same type but of different kinds compare by what they hold, without recursing
    @test isless(exp(X), exp(Dissipator(X))) || isless(exp(Dissipator(X)), exp(X))
    # the ranking puts the types in the order it declares
    @test isless(Id, X)                     # Identity before Operator
    @test isless(X(1), Gate(X)(1))          # AtIndex before Gate
    @test isless(Gate(X)(1), Dissipator(X)(1))
    @test isless(Dissipator(X)(1), Evolver(X(1)))
end

@testset "A projector is self adjoint" begin
    # `matrix` refuses a projector on a mixed state, so every one is self adjoint, and its
    # adjoint simplifies to it whatever it projects on (decision 20: Proj("1"), not Proj(1))
    for p in (Proj("Up"), Proj("1"), Proj([1, im] / √2), Proj(X => -1))
        @test_skip simplify(dag(p)) == p                    # couche 2
        @test matrix(dag(p), Qubit()) ≈ matrix(p, Qubit())
    end
end

@testset "Projector on an eigenspace" begin
    # the eigenvalue is matched to rounding: eigen gives those of Sx on a spin 1 as
    # -0.9999999999999984, about 1e-15 and 0.9999999999999999
    s = Spin(1)
    v = [1., 0., -1.] / √2
    @test matrix(Proj(Sx => 0), s) ≈ v * v'
    @test sum(matrix(Proj(Sx => λ), s) for λ in -1:1) ≈ matrix(Id, s)
    # an eigenspace of several states, which the projector on a state cannot give
    @test matrix(Proj(Ntot => 1), Electron()) ≈ matrix(Proj("Up") + Proj("Dn"), Electron())
    @test matrix(Proj(Level => 1), Qubit()) ≈ matrix(Proj("1"), Qubit())       # decision 20
    # equal, and ordered among the other projectors, by its operator and its eigenvalue
    @test Proj(X => -1) == Proj(X => -1.0)
    @test hash(Proj(X => -1)) == hash(Proj(X => -1.0))
    @test hash(Proj(X => -0.0)) == hash(Proj(X => 0.0))
    @test length(sort([Proj(X => 1), Proj("1"), Proj(X => -1), Proj("Up"), Proj(Z => 1)])) == 5
    # of a definite charge when the operator is
    @test flux(Proj(N => 1), Qubit(conserve = N)) == IT.QN("N", 0)
    @test_throws "no definite flux" flux(Proj(X => 1), Qubit(conserve = N))
    @test_throws "0.5 is not an eigenvalue of Sx" matrix(Proj(Sx => 0.5), s)
    @test_throws "Proj needs a Hermitian operator" matrix(Proj(Sp => 1), Qubit())
    # decision 2: an operator of no definite fermionic parity has its projectors, their strings
    # being laid by the lowering
    @test matrix(Proj(C + dag(C) => 1), Fermion()) ≈ [1. 1. ; 1. 1.] / 2
    @test_throws "Proj(N => 1)(3) rather than Proj(N(3) => 1)" Proj(X(3) => 1)
    # of several sites, the singlet and the spin 2 of two spins one, printed as before whatever
    # its number of sites
    @test Proj("Up") isa SimpleOp
    @test Proj(Swap => -1) isa GenericOp{IsOperator, 2}
    @test repr(Proj("Up")) == "Proj(\"Up\")"
    @test repr(Proj(Swap => -1)) == "Proj(Swap => -1)"
    @test matrix(Proj(Swap => -1), Qubit()) ≈ matrix((Id ⊗ Id - Swap) / 2, Qubit())
    SS = Sx ⊗ Sx + Sy ⊗ Sy + Sz ⊗ Sz
    @test matrix(named(Proj(SS => 1), "P2"), s) ≈ matrix(SS / 2 + SS * SS / 6 + (Id ⊗ Id) / 3, s)
    @test_throws "0.5 is not an eigenvalue of Swap on Qubit() ⊗ Qubit()" matrix(Proj(Swap => 0.5), Qubit())
    @test_throws "Proj needs a Hermitian operator" matrix(Proj(Sp ⊗ Sm => 0), Qubit())
end

@testset "Non integer powers of a scaled operator" begin
    # the power is the principal one. A positive coefficient comes out of it unchanged, any
    # other phase does not and stays inside, only the modulus coming out: a negative one is
    # the phase -1, and so is -1 + 0im. The power decides this when it is built and simplify
    # when a sum collapses into a scaled operator afterwards, and the two have to reach the
    # same form
    q = Qubit()
    for (c, a) in [(-1., -X - X), (im, im * X + im * X), (2., 2X + 2X),
                   (1. + im, (1 + im) * X + (1 + im) * X),
                   (-1. + 0im, (-1. + 0im) * X + (-1. + 0im) * X)]
        @test_skip simplify(a^0.5) == (2c * X)^0.5          # couche 2
        @test matrix((2c * X)^0.5, q) ≈ sqrt(2c * matrix(X, q))
    end
    @test (2X)^0.5 == sqrt(2.) * X^0.5
    @test (-2X)^0.5 == sqrt(2.) * (-X)^0.5
    # the sign really does leave the coefficient for the operator
    @test_skip string(simplify((-X - X)^0.5)) == string(sqrt(2.) * (-X)^0.5)   # couche 2
    # an integer exponent needs none of this care
    @test_skip simplify((-X - X)^2) == 4Id                  # couche 2
end

@testset "Operators as dictionary keys" begin
    # a type with an `==` of its own needs a matching `hash`: `Set` and `Dict` pick their
    # bucket by hash and only compare within it, so two equal operators would otherwise
    # land apart. The default hash follows the identity of the `subs` vector rather than
    # its contents, so it does not do. What matters as much is that every form is covered,
    # including the wrappers that merely hold one of those vectors: they used to fall back
    # to `===` and never compared equal to themselves
    for (a, b) in [(X(1) * Y(2), X(1) * Y(2)),                          # ProdOp
                   (X(1) + Y(2), X(1) + Y(2)),                          # SumOp
                   (X ⊗ Y, X ⊗ Y),                                      # TensorOp
                   ((X * Y)(3), (X * Y)(3)),                            # AtIndex of a product
                   (X(1) * Y(2) + (X * Y)(3), X(1) * Y(2) + (X * Y)(3)),
                   (2 * X(1) * Y(2), 2 * X(1) * Y(2)),                  # ScalarOp
                   (dag(X * Y), dag(X * Y)),                            # DagOp
                   ((X * Y)^2, (X * Y)^2),                              # IntPowOp
                   ((X * Y)^0.5, (X * Y)^0.5),                          # GenPowOp
                   (exp(X * Y), exp(X * Y)),                            # ExpOp
                   (parity(X * Y), parity(X * Y)),                      # ModOp
                   (Proj(X => 1), Proj(X => 1.0)),                      # Proj
                   (Gate(X + Y), Gate(X + Y)),                          # Gate
                   (Gate(X(1) + Y(2)), Gate(X(1) + Y(2))),
                   (Reset([1. 0. ; 0. 0.]), Reset([1. 0. ; 0. 0.])),    # Reset
                   (Dephase(X ⊗ Y), Dephase(X ⊗ Y)),                    # Dephase
                   (Dissipator(X * Y), Dissipator(X * Y)),              # Dissipator
                   (Evolver(X(1) * Y(2)), Evolver(X(1) * Y(2))),        # Evolver
                   (exp(Dissipator(X)), exp(Dissipator(X))),            # ExpOp of a generator
                   # an Operator built afresh each call, whose expr is a matrix for one
                   # and a whole expression for the other
                   (Phase(0.3), Phase(0.3)),
                   (controlled(Z), controlled(Z))]
        @test a == b
        @test hash(a) == hash(b)
        @test length(Set([a, b])) == 1
    end
    @test_skip "a product the descent gave strings"   # couche 2
    # and operators that differ must keep differing
    @test length(Set([X(1) * Y(2), X(1) * Z(2)])) == 2
    # what it is for: a product shared by two measurements is asked for only once
    @test_skip length(Set(TMS.get_prods(State{Pure}(System(3, Qubit()), "Up"),
                                        Measure([X(1) * Y(2) + Z(1), X(1) * Y(2) + Z(3)])))) == 3   # couche 3
end

@testset "Modulo and parity operators" begin
    f, b = Fermion(), Boson(6)

    # mod(A, m) is exp(2iπA/m), so its eigenvalues are the m-th roots of unity of those of A
    @test matrix(parity(N), f) ≈ [1 0 ; 0 -1]
    @test matrix(parity(N), b) ≈ [i == j ? (-1.)^(i - 1) : 0. for i in 1:6, j in 1:6]
    @test matrix(mod(N, 3), b) ≈
        [i == j ? exp(2im * π * mod(i - 1, 3) / 3) : 0im for i in 1:6, j in 1:6]
    @test matrix(parity(N), b) ≈ matrix(mod(N, 2), b)

    # simplify must never change the operator it stands for
    for a in [parity(N), mod(N, 3), dag(parity(N)), dag(mod(N, 3)), parity(2N), parity(N)^2]
        @test_skip matrix(simplify(a), b) ≈ matrix(a, b)    # couche 2
    end

    # the adjoint puts the sign into the operator instead of leaving a dag outside
    @test_skip matrix(simplify(dag(mod(N, 3))), b) ≈ adjoint(matrix(mod(N, 3), b))    # couche 2
    @test matrix(dag(mod(N, 3)), b) ≈ adjoint(matrix(mod(N, 3), b))
    @test matrix(dag(dag(mod(N, 3))), b) ≈ matrix(mod(N, 3), b)

    # (-1)^N is hermitian, whatever form simplify settles on
    @test matrix(parity(N), b) ≈ adjoint(matrix(parity(N), b))

    @test_throws "a modulus is at least 2, got 1" mod(N, 1)
    # decision 2: an operator of no definite fermionic parity has its parity operator,
    # exp(iπC) = 1 + iπC
    @test matrix(parity(C), f) ≈ [1 im * π ; 0 1]
end

@testset "Renaming an operator" begin
    f = Fermion()
    # the definition is kept and only the label changes: a plain relabelling would send the
    # lookup after a name no site defines
    @test matrix(named(N, "Nf"), f) ≈ matrix(N, f)
    @test repr(named(N, "Nf")) == "Nf"
    @test repr(named(N, "Nf")(3)) == "Nf(3)"
    @test_skip matrix(simplify(named(N, "Nf")), f) ≈ matrix(N, f)    # couche 2
    # an expression has no name of its own and can be renamed too
    @test matrix(named(2Sz, "SzA"), Spin(1/2)) ≈ matrix(2Sz, Spin(1/2))
    # decision 2: a renamed operator has no type, the lowering reading its fermionic parity off
    # its matrix, so one of no definite parity is named as any other
    @test matrix(named(C + N, "x"), f) ≈ matrix(C + N, f)
end

@testset "Flux of an operator" begin
    # the flux is the charge an operator carries, read off the charges of the states it
    # connects. A site conserving nothing puts no constraint on anything
    @test flux(X, Qubit()) == IT.QN()
    @test flux(N, Fermion(conserve = N)) == IT.QN("N", 0)
    @test flux(C, Fermion(conserve = N)) == IT.QN("N", -1)
    @test flux(dag(C), Fermion(conserve = N)) == IT.QN("N", 1)

    # in units of the declared charge, which is why a half integer one is written doubled
    @test flux(Sp, Qubit(conserve = 2Sz)) == IT.QN("2Sz", 2)
    @test flux(Sm, Qubit(conserve = 2Sz)) == IT.QN("2Sz", -2)

    # the difference is taken modulo the charge, without which the shift operator, which
    # generates the very symmetry Z records, would be refused on its wrap around
    @test flux(Z, Qudit(3, conserve = Z)) == IT.QN("Z", 0, 3)
    @test flux(X, Qudit(3, conserve = Z)) == IT.QN("Z", 1, 3)
    @test flux(X^2, Qudit(3, conserve = Z)) == IT.QN("Z", 2, 3)
    @test flux(A, Boson(4, conserve = parity(N))) == IT.QN("parity(N)", 1, 2)

    # several charges at once
    @test flux(Cup, Electron(conserve = (Ntot, 2Sz))) == IT.QN(("Ntot", -1), ("2Sz", -1))
    @test flux(Ntot, Electron(conserve = (Ntot, 2Sz))) == IT.QN(("Ntot", 0), ("2Sz", 0))

    # an operator raising and lowering at once carries no flux, and the message says which
    # two differences it found
    @test_throws "differing by 2Sz=2 and by 2Sz=-2" flux(X, Qubit(conserve = 2Sz))
    @test_throws "no definite flux" flux(H, Qudit(3, conserve = Z))
    # under parity the same X is fine, the two differences becoming one modulo 2
    @test flux(X, Qubit(conserve = parity(2Sz))) == IT.QN("parity(2Sz)", 0, 2)

    @test_throws "only defined for one site operators" flux(Swap, Qubit(conserve = 2Sz))
end

@testset "A matrix knows no charge" begin
    # the matrix of an operator is written in the basis of its sites whatever they conserve,
    # an operator on several charged sites included
    f, fq, fs = Fermion(), Fermion(conserve = N), Fermion(conserve = strong(N))
    @test matrix(C, fq) == matrix(C, f)
    @test matrix(C, fs) == matrix(C, f)
    for op in (Gate(C), Dissipator(C), Evolver(C), Reset("Occ"), Dephase(N))
        @test matrix(op, fq) == matrix(op, f)
        @test matrix(op, fs) == matrix(op, f)
    end
    e, eq = Electron(), Electron(conserve = (Ntot, 2Sz))
    @test matrix(Ntot ⊗ Sz, eq, eq) == matrix(Ntot ⊗ Sz, e, e)
    b, bq = Boson(4), Boson(4, conserve = parity(N))
    @test matrix(N ⊗ N, bq, bq) == matrix(N ⊗ N, b, b)

    # an operator carrying no definite flux on a site still has a matrix: what is refused is
    # placing it on a system that conserves
    @test matrix(X, Qubit(conserve = N)) == matrix(X, Qubit())

    # a single site stands for as many identical ones as the operator needs
    q = Qubit()
    @test matrix(Swap, q) == matrix(Swap, q, q)
    @test matrix(Gate(Swap), q) == matrix(Gate(Swap), q, q)
end

@testset "Operators of several sites given by a matrix" begin
    # what was left of splitting them on their sites once their sites go: the matrix of an
    # operator given by a matrix, a function of the sites or an expression is what it was given
    # (decision 2: an operator has no type, and its sites are those of its matrix)
    q, s1, fe = Qubit(), Spin(1), Fermion()
    msw = [1. 0 0 0 ; 0 0 1 0 ; 0 1 0 0 ; 0 0 0 1]
    @test matrix(Operator{2}("Sw", msw), q) ≈ msw
    ss = Sx ⊗ Sx + Sy ⊗ Sy + Sz ⊗ Sz
    p2e = Operator{2}("P2e", ss / 2 + ss * ss / 6 + (Id ⊗ Id) / 3)
    mp2 = matrix(p2e, s1)
    @test matrix(Operator{2}("P2", mp2), s1) ≈ mp2
    # sites that differ, given in the order of the indices, and more than two of them
    mk = matrix(Sz ⊗ X + Sp ⊗ Z + 0.3 * Sx ⊗ Id, s1, q)
    @test matrix(Operator{2}("K", mk), s1, q) ≈ mk
    m3 = matrix(X ⊗ Z ⊗ Y + Z ⊗ Id ⊗ Z + 0.5 * Id ⊗ X ⊗ Id + 0.2 * Id ⊗ Id ⊗ Id, q)
    @test matrix(Operator{3}("T3", m3), q) ≈ m3
    # a function of the sites, an expression kept whole, and a complex matrix
    zx = Operator{2}("ZX", (a, b) -> kron(matrix(Z, a), matrix(X, b)))
    @test matrix(zx, q) ≈ kron(matrix(Z, q), matrix(X, q))
    rxx = exp(-0.3im * (X ⊗ X))
    @test matrix(Operator{2}("R", rxx), q) ≈ matrix(rxx, q)
    @test matrix(Operator{2}("Y2", matrix(Y ⊗ Sp, q)), q) ≈ matrix(Y ⊗ Sp, q)
    # its MPO, expect and apply act as the matrix it was given by, on sites in any order, and on
    # a system mixing its sites
    chain(op) = sum(op(i, i + 1) for i in 1:3)
    st = random_state(Pure, System(4, q), 4)
    sw = Operator{2}("Sw", msw)
    @test maxlinkdim(make_mpo(st, chain(sw))) == maxlinkdim(make_mpo(st, chain(Swap)))
    @test norm(apply(make_mpo(st, chain(sw)), st) - apply(make_mpo(st, chain(Swap)), st)) < 1e-12
    st1 = random_state(Pure, System(4, s1), 4)
    p2 = Operator{2}("P2", mp2)
    # the eight terms of the projector of the AKLT chain across a link, and the two channels
    # of the term not yet begun and of the term finished
    @test maxlinkdim(make_mpo(st1, chain(p2))) == 10
    for (i, j) in [(1, 3), (4, 2), (2, 1)]
        @test expect(st1, p2(i, j)) ≈ expect(st1, p2e(i, j))
        @test norm(apply(p2(i, j), st1) - apply(p2e(i, j), st1)) < 1e-12
    end
    k = Operator{2}("K", mk)
    sk = random_state(Pure, System([s1, q, s1, q]), 4)
    for (i, j) in [(1, 2), (3, 2), (1, 4)]
        @test norm(apply(k(i, j), sk) - apply(make_mpo(sk, k(i, j)), sk)) < 1e-12
    end
    # a single site, a fermionic one included
    @test matrix(Operator{1}("Ex", exp(0.3im * X)), q) ≈ matrix(exp(0.3im * X), q)
    @test matrix(Operator{1}("Fs", s -> matrix(Sz, s)), s1) ≈ matrix(Sz, s1)
    @test matrix(Operator{1}("Cf", C), fe) ≈ matrix(C, fe)
    # what cannot be: sites of another number, or a matrix, or a function of the sites, of a size
    # other than theirs
    @test_throws "acts on 2 sites and was given 3" matrix(Operator{2}("A", msw), q, q, q)
    @test_throws "whose dimension is 6" matrix(Operator{2}("A", msw), s1, q)
    @test_throws "whose dimension is 3" matrix(Operator{1}("A", [1. 0 ; 0 1]), s1)
    @test_throws "whose dimension is 3" matrix(Operator{1}("Fn", s -> [1. 0 ; 0 1]), s1)
end

@testset "Operators of several sites split into factors of one site" begin
    # the lowering writes an operator of several sites as a sum of tensor products of operators
    # of one site, `(coef, factors)` with `nothing` for the identity
    q, s1, fe = Qubit(), Spin(1), Fermion()
    factor(f) = isnothing(f) ? Id : Operator{1}("f", f)
    # the terms give back the matrix, the strings of ⊗ included, and every factor has a definite
    # charge, which flux would refuse otherwise, and a definite parity
    exact(ts, m, sites) =
        sum(c * matrix(⊗(map(factor, fs)...), sites...) for (c, fs) in ts) ≈ m &&
        all(isnothing(f) || flux(factor(f), s) isa IT.QN &&
                            !isnothing(TMS.parity_of(f, matrix(F, s)))
            for (_, fs) in ts for (f, s) in zip(fs, sites))
    spans(ts) = sort([ findall(!isnothing, fs) for (_, fs) in ts ])
    crossing(ts) = count(t -> count(!isnothing, t[2]) == 2, ts)

    # the identity is taken out on each site before the split, so that the fewest terms cross
    # the link: three for Swap, as its expression has, and eight for the projector of the AKLT
    # chain, where its expression takes twelve, the same with a charge conserved
    msw = [1. 0 0 0 ; 0 0 1 0 ; 0 1 0 0 ; 0 0 0 1]
    sw = TMS.split_matrix(Pure, "Sw", msw, [q, q])
    @test exact(sw, msw, [q, q]) && crossing(sw) == 3
    ss = Sx ⊗ Sx + Sy ⊗ Sy + Sz ⊗ Sz
    mp2 = matrix(ss / 2 + ss * ss / 6 + (Id ⊗ Id) / 3, s1)
    p2 = TMS.split_matrix(Pure, "P2", mp2, [s1, s1])
    @test exact(p2, mp2, [s1, s1]) && crossing(p2) == 8
    s1q = Spin(1, conserve = Sz)
    p2q = TMS.split_matrix(Pure, "P2", mp2, [s1q, s1q])
    @test exact(p2q, mp2, [s1q, s1q]) && crossing(p2q) == 8
    # sites that differ, and more than two of them: no term spans more sites than it acts on
    mk = matrix(Sz ⊗ X + Sp ⊗ Z + 0.3 * Sx ⊗ Id, s1, q)
    @test exact(TMS.split_matrix(Pure, "K", mk, [s1, q]), mk, [s1, q])
    m3 = matrix(X ⊗ Z ⊗ Y + Z ⊗ Id ⊗ Z + 0.5 * Id ⊗ X ⊗ Id + 0.2 * Id ⊗ Id ⊗ Id, q)
    t3 = TMS.split_matrix(Pure, "T3", m3, [q, q, q])
    @test exact(t3, m3, [q, q, q]) && spans(t3) == [Int[], [1, 2, 3], [1, 3], [2]]
    # real factors for a real operator, even given as a complex matrix
    rsw = TMS.split_matrix(Pure, "Sw", complex(msw), [q, q])
    @test all(f -> isnothing(f) || eltype(f) <: Real, (f for (_, fs) in rsw for f in fs))
    rxx = matrix(exp(-0.3im * (X ⊗ X)), q)
    @test exact(TMS.split_matrix(Pure, "R", rxx, [q, q]), rxx, [q, q])
    # a charge taken modulo, and an operator of no definite charge
    bp = Boson(4, conserve = parity(N))
    mb = matrix(A ⊗ dag(A) + dag(A) ⊗ A + 0.5 * (A * A) ⊗ N, Boson(4))
    @test exact(TMS.split_matrix(Pure, "Bp", mb, [bp, bp]), mb, [bp, bp])
    qn = Qubit(conserve = N)
    @test_throws "no definite charge of N" TMS.split_matrix(Pure, "XX", matrix(X ⊗ X, q), [qn, qn])

    # on fermions the matrix holds the strings between its sites, which the split takes out:
    # C ⊗ dag(C), refused by 1.x, is a single term, and the string that C on the second site
    # lays on the first is no factor there
    hop = matrix(C ⊗ dag(C), fe)
    ts = TMS.split_matrix(Pure, "H", hop, [fe, fe])
    @test exact(ts, hop, [fe, fe]) && length(ts) == 1
    c2 = matrix(Id ⊗ C, fe)
    ts = TMS.split_matrix(Pure, "C2", c2, [fe, fe])
    @test exact(ts, c2, [fe, fe]) && spans(ts) == [[2]]
    # an operator of no definite parity has its two parts split apart
    ic = matrix(Id ⊗ Id, fe) + c2
    ts = TMS.split_matrix(Pure, "IC", ic, [fe, fe])
    @test exact(ts, ic, [fe, fe]) && spans(ts) == [Int[], [2]]
    # the fermionic Swap, a function of fermionic operators, three fermionic sites, a site of
    # another type between two, two species on a site, and a charge conserved
    h = dag(C) ⊗ C + dag(dag(C) ⊗ C)
    el, fn = Electron(), Fermion(conserve = N)
    for (m, sites) in [(matrix(Swap, fe), [fe, fe]),
                       (matrix(exp(-0.3 * (h * h)), fe), [fe, fe]),
                       (matrix(C ⊗ N ⊗ dag(C), fe), [fe, fe, fe]),
                       (matrix(C ⊗ X ⊗ dag(C), fe, q, fe), [fe, q, fe]),
                       (matrix(dag(Cup) ⊗ Cup + dag(Cdn) ⊗ Cdn + Ntot ⊗ Ntot, el), [el, el]),
                       (hop, [fn, fn])]
        @test exact(TMS.split_matrix(Pure, "M", m, sites), m, sites)
    end
end

@testset "The descent of an operator onto a system" begin
    # the descent writes a placed operator as terms, each a coefficient and one matrix per site,
    # the Jordan-Wigner strings laid: their sum is checked against the explicit matrix of each
    # factor of one site and of definite parity, F ⊗ … ⊗ F ⊗ a ⊗ 1 ⊗ … for an odd one
    sites(sys) = [ sys[k] for k in 1:length(sys) ]
    identity_on(sys) = foldl(kron, [ matrix(Id, s) for s in sites(sys) ])
    whole(sys, t) =
        foldl(kron, [ get(Dict(t.blocks), k, matrix(Id, s)) for (k, s) in enumerate(sites(sys)) ])
    lowered(sys, op) =
        sum((t.coef * whole(sys, t) for t in TMS.lower(sys, op)); init = 0 * identity_on(sys))
    function jw(sys, j, a)
        odd = TMS.parity_of(a, matrix(F, sys[j])) == 1
        return foldl(kron, [ k < j && odd ? matrix(F, s) : k == j ? a : matrix(Id, s)
                             for (k, s) in enumerate(sites(sys)) ])
    end
    q, fe, el = Qubit(), Fermion(), Electron()
    x, y, z = matrix(X, q), matrix(Y, q), matrix(Z, q)
    s3 = System(3, q)
    @test lowered(s3, X(1) + 2 * Z(2) * Y(3)) ≈ jw(s3, 1, x) + 2 * jw(s3, 2, z) * jw(s3, 3, y)

    # fermions: a hop over a fermion, a product written backwards, and Swap or a placed tensor
    # product on sites in decreasing order, which need no sign of reordering
    s4 = System(4, fe)
    c = matrix(C, fe)
    cd, n = collect(c'), c' * c
    hop = jw(s4, 1, cd) * jw(s4, 3, c) + jw(s4, 3, cd) * jw(s4, 1, c)
    @test lowered(s4, dag(C)(1) * C(3) + dag(C)(3) * C(1)) ≈ hop
    @test lowered(s4, C(3) * dag(C)(1)) ≈ jw(s4, 3, c) * jw(s4, 1, cd)
    sw = identity_on(s4) - jw(s4, 1, n) - jw(s4, 3, n) + hop
    @test lowered(s4, Swap(1, 3)) ≈ sw
    @test lowered(s4, Swap(3, 1)) ≈ sw
    @test lowered(s4, (C ⊗ dag(C))(3, 1)) ≈ jw(s4, 3, c) * jw(s4, 1, cd)
    @test lowered(s4, (C ⊗ N ⊗ dag(C))(4, 2, 1)) ≈ jw(s4, 4, c) * jw(s4, 2, n) * jw(s4, 1, cd)
    @test lowered(s4, dag(C ⊗ C)(1, 2)) ≈ (jw(s4, 1, c) * jw(s4, 2, c))'
    # a function of several fermionic sites, placed both ways, with a fermion between them
    h = dag(C) ⊗ C + dag(dag(C) ⊗ C)
    @test lowered(s4, exp(-0.3 * (h * h))(1, 3)) ≈ exp(-0.3 * hop^2)
    @test lowered(s4, exp(-0.3 * (h * h))(3, 1)) ≈ exp(-0.3 * hop^2)
    # a block of no definite parity is split into its even and odd parts, which 1.x refused
    @test lowered(s4, (C + N)(2)) ≈ jw(s4, 2, c) + jw(s4, 2, n)
    @test lowered(s4, exp(C)(2)) ≈ identity_on(s4) + jw(s4, 2, c)
    # strings on the fermionic sites only, and two species on a site
    sm = System([fe, q, fe])
    @test lowered(sm, C(1) * X(2) * dag(C)(3)) ≈ jw(sm, 1, c) * jw(sm, 2, x) * jw(sm, 3, cd)
    se = System(3, el)
    cu, cn, nt = matrix(Cup, el), matrix(Cdn, el), matrix(Ntot, el)
    @test lowered(se, dag(Cup)(1) * Cup(3) + dag(Cdn)(3) * Cdn(1) + (Ntot ⊗ Ntot)(2, 3)) ≈
          jw(se, 1, collect(cu')) * jw(se, 3, cu) + jw(se, 3, collect(cn')) * jw(se, 1, cn) +
          jw(se, 2, nt) * jw(se, 3, nt)
    @test lowered(se, (dag(Cup) ⊗ Cup + dag(Cdn) ⊗ Cdn)(3, 1)) ≈
          jw(se, 3, collect(cu')) * jw(se, 1, cu) + jw(se, 3, collect(cn')) * jw(se, 1, cn)
    # an operator of several sites given by a function of its sites
    q3 = Qudit(3)
    @test lowered(System(2, q3), Qudits.Sum(1, 2)) ≈ matrix(Qudits.Sum, q3, q3)

    # the reduction: an identity goes into the coefficient, a product that vanishes drops its
    # term, blocks equal up to a coefficient are gathered, and a real operator keeps real terms
    s2 = System(2, q)
    terms(op) = [ (t.coef, first.(t.blocks)) for t in TMS.lower(s2, op) ]
    @test terms(X(1) * X(1)) == [(1., Int[])]
    @test terms(X(1) * Y(1)) == [(im, [1])]
    @test only(TMS.lower(s2, X(1) * Y(1))).blocks == [1 => z]
    @test terms(X(1) + X(1)) == [(2., [1])]
    @test isempty(terms(X(1) - X(1))) && isempty(terms(Sp(1) * Sp(1)))
    yy = only(TMS.lower(s2, Y(1) * Y(2)))
    @test yy.coef === -1. && all(m -> eltype(m) == Float64, last.(yy.blocks))
    # the strings of a hop cancel on its left: each term stays on its two sites
    chain = TMS.lower(System(6, fe), sum(dag(C)(i) * C(i + 1) + dag(C)(i + 1) * C(i) for i in 1:5))
    @test sort([ first.(t.blocks) for t in chain ]) == sort([ [i, i + 1] for i in 1:5 for _ in 1:2 ])

    # on a charged system, every block has a definite flux and every term the same one
    fn = Fermion(conserve = N)
    hops = dag(C)(1) * C(3) + dag(C)(3) * C(1)
    @test lowered(System(3, fn), hops) ≈ lowered(System(3, fe), hops)
    @test_throws "change the charges by" TMS.lower(System(3, fn), C(1) + dag(C)(2))
    @test_throws "no definite flux" TMS.lower(System(2, Qubit(conserve = N)), X(1))
    # a site the system does not have, named as written
    @test_throws "X(5) acts on site 5, which the system does not have: it has 3 sites" TMS.lower(s3, X(5))
end

@testset "Strings kept apart while the terms are developed" begin
    # the strings of the odd factors are laid only at the end: products whose strings cross
    # blocks of both parities, in both orders, and an adjoint, give the matrix of the operator
    # built from the matrices of its factors placed
    fe = Fermion()
    fs = System(5, fe)
    f, c, n, id = matrix(F, fe), matrix(C, fe), matrix(N, fe), matrix(Id, fe)
    placed(m, k, odd) = foldl(kron, [ j < k ? (odd ? f : id) : j == k ? m : id for j in 1:5 ])
    sw = kron(id, matrix(Swap, fe, fe), id, id)
    for (op, m) in [(dag(C)(1) * N(3) * C(5), placed(c', 1, true) * placed(n, 3, false) * placed(c, 5, true)),
                    (C(5) * dag(C)(2) * C(3) * dag(C)(4),
                     placed(c, 5, true) * placed(c', 2, true) * placed(c, 3, true) * placed(c', 4, true)),
                    ((C + dag(C))(2) * C(4) * N(1),
                     placed(c + c', 2, true) * placed(c, 4, true) * placed(n, 1, false)),
                    (C(1) * Swap(2, 3) * dag(C)(4), placed(c, 1, true) * sw * placed(c', 4, true)),
                    (dag(C ⊗ Id ⊗ C)(1, 2, 4), (placed(c, 1, true) * placed(c, 4, true))')]
        @test TMS.system_matrix(Pure, fs, TMS.lower(fs, op)) ≈ m
    end
    # and so a hopping on a long chain takes no string outside the sites it hops between
    long = System(1000, fe)
    ts = TMS.lower(long, sum(dag(C)(i) * C(i + 1) + dag(C)(i + 1) * C(i) for i in 1:999))
    @test length(ts) == 1998
    @test all(t -> length(t.blocks) == 2 && last(first.(t.blocks)) == first(first.(t.blocks)) + 1, ts)
end

@testset "Products gathered before the descent" begin
    # the first pass of the descent develops nothing: it sorts the factors of a product by
    # their sites as far as they (anti)commute, and merges those sharing a site, the sites of
    # one among those of the other, when they can be brought together, which leaves fewer terms
    # to develop and the operator unchanged
    plain(sys, op) = TMS.system_matrix(Pure, sys, TMS.reduced(TMS.op_terms(sys, op, ())))
    same(sys, op) = TMS.system_matrix(Pure, sys, TMS.lower(sys, op)) ≈ plain(sys, op)
    nterms(sys, op) = length(TMS.lower(sys, op))
    q, fe = Qubit(), Fermion()
    s5, f4, f5 = System(5, q), System(4, fe), System(5, fe)
    @test TMS.gathered(s5, Swap(1, 2) * X(5) * Swap(1, 2)) == X(5)
    @test nterms(s5, Swap(1, 2) * Z(1) * Z(2) * Swap(1, 2)) == 1
    @test TMS.gathered(s5, Swap(1, 2) * Swap(2, 1)) isa TMS.IdentityOp
    @test TMS.gathered(f5, Swap(1, 3) * Swap(3, 1)) isa TMS.IdentityOp
    @test TMS.gathered(s5, Swap(1, 2) * X(5) * Swap(1, 2) + Swap(3, 4) * Z(1) * Swap(3, 4)) ==
          X(5) + Z(1)
    # two odd factors that change places take a sign
    g = TMS.gathered(f4, C(1) * C(2) * dag(C)(1))
    @test TMS.scalarcoef(g) == -1 && last(TMS.prodsubs(g)) == C(2)
    # a factor goes back past one that commutes with it, even placed between its sites
    g = TMS.gathered(s5, Swap(1, 3) * X(2) * Z(3))
    @test length(TMS.prodsubs(g)) == 2 && last(TMS.prodsubs(g)) == X(2)
    # a factor linked to both keeps them apart, and a merge that would give more terms is
    # refused: the two hops would merge into n₁(1 - n₄), cut into 4 terms by its identity parts
    @test TMS.gathered(s5, Swap(1, 2) * Swap(2, 3) * Swap(1, 2)) == Swap(1, 2) * Swap(2, 3) * Swap(1, 2)
    hops = (dag(C) ⊗ C)(1, 4) * Swap(2, 3) * (dag(C) ⊗ C)(4, 1)
    @test TMS.gathered(f5, hops) == (dag(C) ⊗ C)(1, 4) * (dag(C) ⊗ C)(4, 1) * Swap(2, 3)
    @test nterms(f5, hops) == 4

    # the result is a normal form, which the pass leaves as it is, and the operator stays the
    # same
    for (sys, op) in [(s5, Swap(2, 3) * Swap(1, 2) * Swap(3, 4)), (s5, X(1) * Swap(2, 3) * Swap(1, 2) * X(1)),
                      (s5, Z(1) * Swap(1, 2) * Y(1) * Swap(2, 3) * X(2)), (f4, C(1) * C(2) * dag(C)(1)),
                      (f5, N(2) * C(5) * N(2) * dag(C)(2)), (f5, hops)]
        g = TMS.gathered(sys, op)
        @test TMS.gathered(sys, g) == g
        @test same(sys, op)
    end
    @test TMS.gathered(s5, Swap(2, 3) * Swap(1, 2) * Swap(3, 4)) == Swap(2, 3) * Swap(1, 2) * Swap(3, 4)
end

@testset "tensor is ⊗ under a name in ASCII" begin
    @test tensor === ⊗
    @test tensor(X, Y) == X ⊗ Y
    @test tensor(X, Y, Z, X) == ((X ⊗ Y) ⊗ Z) ⊗ X
    @test tensor(System(1, Qubit()), System(2, Qubit())) == System(3, Qubit())
end

@testset "The matrix of a tensor product of fermions" begin
    # (A ⊗ B) on consecutive sites is A(1) * B(2), Jordan-Wigner strings included: its matrix
    # is checked against the explicit c₁ = c ⊗ 1, c₂ = F ⊗ c and c₃ = F ⊗ F ⊗ c
    fe = Fermion()
    c, f, id = matrix(C, fe), matrix(F, fe), matrix(Id, fe)
    c1, c2 = kron(c, id), kron(f, c)
    @test matrix(C ⊗ dag(C), fe) ≈ c1 * c2'
    @test matrix(dag(C) ⊗ C, fe) ≈ c1' * c2
    @test matrix(C ⊗ C, fe) ≈ c1 * c2
    @test matrix(N ⊗ C, fe) ≈ c1' * c1 * c2
    @test matrix((C + N) ⊗ C, fe) ≈ kron(c + c' * c, id) * c2
    @test matrix(C ⊗ Id ⊗ dag(C), fe) ≈ kron(c, id, id) * kron(f, f, c)'
    @test matrix(dag(C ⊗ C), fe) ≈ -matrix(dag(C) ⊗ dag(C), fe)
    # a factor of no definite parity after a fermionic site makes the product a sum
    @test_throws "no definite fermionic parity" matrix(C ⊗ (C + N), fe)

    # an operator of several sites given by such an expression is the operator the expression
    # stands for (decision 2: with no type)
    h = dag(C) ⊗ C + dag(dag(C) ⊗ C)
    hm = c1' * c2 + c2' * c1
    @test matrix(h, fe) ≈ hm
    @test matrix(Operator{2}("H2", h * h), fe) ≈ hm^2
    e = Operator{2}("E", exp(-0.3 * (h * h)))
    @test matrix(e, fe) ≈ exp(-0.3 * hm^2)
    @test expect(State{Pure}(System(3, fe), ["Occ", "Emp", "Emp"]), e(1, 2)) ≈ exp(-0.3)
    @test expect(State{Pure}(System(3, fe), ["Occ", "Emp", "Emp"]), e(1, 3)) ≈ exp(-0.3)

    # a factor of several sites has the parity of its pieces placed, a renamed expression as
    # well as a tensor product. Both were taken as even, which lost the sign of the adjoint and
    # the strings of the matrix, and made a dissipator change the trace
    cn = named(C ⊗ N, "CN")
    a2 = C ⊗ Id + Id ⊗ C
    @test_skip "the norm of a(1, 2, 3, 4) applied to a state, for cn ⊗ cn and a2 ⊗ a2"   # couche 2
    @test matrix(C ⊗ cn, fe) ≈ matrix(C ⊗ C ⊗ N, fe)
    @test matrix(C ⊗ a2, fe) ≈ matrix(C ⊗ C ⊗ Id + C ⊗ Id ⊗ C, fe)
    @test_throws "no definite fermionic parity" matrix(C ⊗ (C ⊗ Id + N ⊗ Id), fe)
    @test_skip "a dissipator of cn ⊗ cn keeps the trace"     # couche 2
end

@testset "Powers" begin
    # a literal exponent goes through literal_pow, which called inv for a negative one
    p = -1
    @test X^-1 == X^p
    @test (Id + 0.5X)^-1 == (Id + 0.5X)^p
    # an integer power is a product, with the adjoint and the strings of that product; any
    # other is the principal power, a function of the operator as exp is, of which only the
    # modulus of a coefficient comes out. Decision 2: an involution is not reduced, its matrix
    # saying what it is
    q, fe, s1 = Qubit(), Fermion(), Spin(1)
    @test matrix(X^10000, q) ≈ matrix(Id, q)
    @test matrix(X^10001, q) ≈ matrix(X, q)
    @test (2X)^3 == 8 * X^3
    # an integer coefficient is stored as a float, whose products do not wrap around as those
    # of integers silently do: (2X)^64 was 0Id, prod(2Sz(i) for i in 1:63) measured -1 on
    # sites up, and (2X)^-1 raised the DomainError of a negative power of an integer
    @test (2X)^64 == 2.0^64 * X^64
    @test (2X)^-1 == 0.5 * X^-1
    @test expect(State{Pure}(System(2, q), "Up"), (2^32 * Sz(1)) * (2^32 * Sz(2))) ≈ 2.0^62
    @test_skip simplify(Sz^2 * Sz) == Sz^3                  # couche 2
    for (a, s) in [(X^3, q), ((X + Y)^3, q), (Sz^2, s1), ((Sp + Sm)^3, s1),
                   ((C + dag(C))^2, fe), (C^3, fe), (dag((X + im * Y)^2), q)]
        @test_skip matrix(simplify(a), s) ≈ matrix(a, s)    # couche 2
    end
    @test matrix(dag((X + im * Y)^2), q) ≈ matrix((X + im * Y)^2, q)'

    # a phase stays inside: the principal square root of -Z has i on its eigenvalue -1
    @test matrix(sqrt(-Z), q) ≈ [im 0 ; 0 1]
    up = State{Pure}(System(2, q), "Up")
    @test expect(up, sqrt(-Z)(1)) ≈ im
    @test expect(up, ((im * X)^0.5)(1)) ≈ matrix((im * X)^0.5, q)[1, 1]
    @test matrix(X^(0.5im), q) ≈ matrix(X, q)^(0.5im)
    # a power that does not exist is refused, where Julia gives zero for C^0.5
    @test_throws "does not exist" matrix(sqrt(C), fe)
    @test_throws "does not exist" matrix(Proj("0")^(-0.5), q)     # decision 20

    # powers of the same operator merge, and what they give leaves no coefficient inside the
    # product nor an F out of its place
    @test_skip simplify((-X)^0.5 * (-X)^0.5) == -X          # couche 2
    @test_skip simplify(X^0.5 * X^0.5) == X                 # couche 2
    @test_skip simplify(F^0.5 * F^0.5 * C) == simplify(F * C)   # couche 2
    # and placed, a merge giving a multiple of the identity leaves only its coefficient: -Id
    # stayed in the product, which the sort by site then failed on
    a = ((im * X)^0.5)(1) * ((im * X)^1.5)(1)
    @test_skip simplify(a * Z(2)) == simplify(Z(2) * a) == simplify(-Z(2))   # couche 2

    # a placed operator takes powers, an integer one of anything placed and any one of an
    # operator on its sites
    @test X(1)^2 == (X^2)(1)
    @test X(1)^0.5 == (X^0.5)(1)
    @test (X(1) + Z(2))^2 == (X(1) + Z(2)) * (X(1) + Z(2))
    st = random_state(Pure, System(3, q), 4)
    @test expect(st, X(1)^2) ≈ 1
    @test expect(st, (X(1) + Z(2))^2) ≈ expect(st, (X(1) + Z(2)) * (X(1) + Z(2)))
    @test expect(st, X(1)^0.5) ≈ expect(st, (X^0.5)(1))
    @test_throws "has no power 0.5 once placed" (X(1) + Z(2))^0.5

    # a non integer power of several sites is kept whole, as exp is: it is refused where it
    # would have to be split, and applied as a gate, a fermion beside it included
    @test_skip "a non integer power of several sites refused by expect and applied as a gate"   # couche 2
    # a string is not moved past a factor of several sites kept whole, which it commutes with
    # only if the factor is even. A factor with no string in the way is still sorted
    xi = Operator{2}("XI", kron([0. 1. ; 1. 0.], [1. 0. ; 0. 1.]))
    @test_skip !(xi(1, 2) * dag(C)(4) ≈ dag(C)(4) * xi(1, 2))   # couche 2
    @test_skip xi(1, 2) * N(4) ≈ N(4) * xi(1, 2)            # couche 2
end

@testset "The identity is one value of each kind" begin
    # every construction of an identity gives the same value, told by its type alone: on
    # several sites, and placed, where it has no site, being the identity of the whole system
    q = Qubit()
    @test Id ⊗ Id isa TMS.IdentityOp
    @test Gate(Id) isa TMS.IdentityOp{IsChannel}
    @test Gate(Id ⊗ Id) == Gate(Id) ⊗ Gate(Id) == Id ⊗ Gate(Id)
    @test Gate(Id)^2 == Gate(Id)
    @test Gate(3Id) == 9 * Gate(Id)
    @test Gate(Id(2)) == Gate(Id)(1)
    # a dissipator of the identity is zero, an evolver of it twice the identity
    @test iszero(TMS.scalarcoef(Dissipator(Id)))
    @test iszero(TMS.scalarcoef(Dissipator(2 * Id(1))))
    @test iszero(TMS.scalarcoef(Dissipator(Gate(Id))))
    @test matrix(Evolver(Id), q) ≈ 2 * matrix(Gate(Id), q)
    @test Id(3) == Id(1)
    @test (Id ⊗ Id)(1, 3) == Id(2)
    @test (Id ⊗ Id)^3 == Id ⊗ Id
    @test matrix((-(Id ⊗ Id))^0.5, q, q) ≈ matrix(im * (Id ⊗ Id), q, q)
    @test X(3)^0 == Id(1)
    @test Id(2)^0.5 == Id(1)
    st = random_state(Pure, System(3, q), 2)
    @test expect(st, Id(2) - X(2) * X(2)) ≈ 0 atol = 1e-12
    # it is its own adjoint, and so measured as real, a number as any expectation value is
    @test_skip last(only(measure(st, Id(2)))) === 1.0       # couche 3
    # its tensor, where one is needed, goes on the site the caller names
    @test expect1(st, Id) ≈ ones(3)
    @test expect1(mix(st), Id) ≈ ones(3)
    @test expect2(st, (Id, X))[1, 3] ≈ expect(st, X(3))
    @test_skip "the MPO of X(1) * X(2) + 0.5 * Id(3), a product by 2Id(2) applied, a gate of Id"   # couche 2
    # having no site, the identity is the identity of any system, whatever site it was given
    @test expect(st, Id(5)) ≈ 1
end

@testset "Operators of several sites that were refused" begin
    # a renamed operator is its definition, developed as the definition itself is
    q = Qubit()
    st = State{Pure}(System(2, q), ["Up", "Dn"])
    @test expect(st, named(exp(-0.3im * Swap), "R")(1, 2)) ≈ expect(st, exp(-0.3im * Swap)(1, 2))
    # one site given for identical ones, a function of the sites included
    zx = Operator{2}("ZX", (a, b) -> kron(matrix(Z, a), matrix(X, b)))
    @test matrix(zx, q) ≈ kron(matrix(Z, q), matrix(X, q))
    # a fermionic operator controlled is its expression (decision 20: Proj("0"), not Proj(0))
    fe = Fermion()
    @test matrix(controlled(C), q, fe) ≈ matrix(Proj("0") ⊗ Id + Proj("1") ⊗ C, q, fe)
end

@testset "Every operator of the libraries has its matrix" begin
    # what is left of checking the type of an operator against its matrix, which decision 2
    # removes: F, which every parity is read against, squares to the identity
    for s in (Qubit(), Fermion(), Electron(), Tj(), Boson(3))
        @test matrix(F, s)^2 ≈ matrix(Id, s)
    end
    # and every operator a site type declares has a matrix of its size, on sites of every size
    for s in [Qubit(), Fermion(), Electron(), Tj(), Spin(0), Spin(1/2), Spin(1), Qudit(1),
              Qudit(2), Qudit(3), Boson(2), Boson(4), Qboson(1., 2), Qboson(0.5, 4)]
        for name in TMS.declared_names(TMS.operator_definition, s)
            @test size(matrix(Operator{1}(name), s)) == (dim(s), dim(s))
        end
    end
end

@testset "named defines an operator" begin
    # decision 2: an operator has no type, and named takes no site, the matrix being read on
    # the sites the operator is given
    @test matrix(named([1 1; 1 -1] / √2, "MyH"), Qubit()) ≈ [1 1; 1 -1] / √2
    @test matrix(named(s -> [0. 1.; 1. 0.], "Fx"), Qubit()) == [0. 1.; 1. 0.]
    @test matrix(named(Sp + Sm, "Sx2"), Spin(1)) ≈ 2 * matrix(Sx, Spin(1))
    @test_throws "given by a 3×3 matrix and cannot act on Qubit()" matrix(named(rand(3, 3), "R"), Qubit())
    # an operator of several sites is given its number of sites, a single site then standing for
    # as many as it acts on, and it can be measured
    sw = matrix(Swap, Qubit())
    s1 = Operator{2}("MySwap", sw)
    @test matrix(s1, Qubit()) ≈ sw
    @test expect(State{Pure}(System(2, Qubit()), ["Up", "X+"]), s1(1, 2)) ≈ 0.5
end

@testset "Signed zeros" begin
    # -0.0 and 0.0 are equal to `==` but not to `isless` or `hash`: simplify, which sorts the
    # terms before merging the equal ones, left a pair apart when a third term sorted between
    # them, and measure computed an expectation value twice
    c1 = -1 * (0.0 - 1.0im)       # -0.0 + 1.0im
    c3 = -1 * (0.0 - 2.0im)       # -0.0 + 2.0im
    f(c) = exp(c * X ⊗ X)(1, 2)
    @test_skip simplify(f(c1) - f(1.0im) + f(c3)) == simplify(f(c3))   # couche 2
    # the coefficient of a product is 0.25 - 0.0im
    a = (-0.5im * X(1)) * (0.5im * Y(2))
    b = 0.25 * X(1) * Y(2)
    @test hash(a) == hash(b)
    @test length(Set([a, b])) == 1
    @test isequal(((X + Z)^complex(-0.0, 0.5)).expo, 0.5im)
    # and so in every number an operator stores: the matrix of an Operator, the state of a
    # Proj or of a Reset (decision 2: an Operator has no type)
    m(x) = Operator{1}("M", [1.0 x; 0.0 1.0])
    @test hash(m(-0.0)) == hash(m(0.0))
    p, q = Proj([1.0, -0.0]), Proj([1.0, 0.0])
    @test hash(p) == hash(q)
    @test !(isless(p, q) || isless(q, p))
    @test hash(Reset([1.0 -0.0; 0.0 0.0])) == hash(Reset([1.0 0.0; 0.0 0.0]))
    # a projector is ordered by its values, which tie where `==` holds: its printed form told
    # [1, 0] and [1.0, 0.0] apart, and [1, 1] sorted between the two
    g(v) = exp(im * Proj(v) ⊗ X)(1, 2)
    @test_skip simplify(g([1, 0]) - g([1.0, 0.0]) + g([1, 1])) == simplify(g([1, 1]))   # couche 2
    @test isless(Proj("1"), Proj("Up")) && isless(Proj("Up"), Proj([1, 0]))      # decision 20
end

@testset "Printing reads back as the operator" begin
    # a name is a column header: it has to denote the operator measured. The right operand of
    # * and ⊗ and the base of ^ lost their parentheses, and so did rational numbers
    for op in (X ⊗ (Y * Z), (Y * Z) ⊗ X, (X ⊗ Y) * (Z ⊗ Z), (X^0.5)^0.5, (1//2) * X,
               X^(1//2), (-1.0 + 0im) * X, 0.5im * X, (X + Z)^0.5, X(1) * Y(2), 2.0X,
               (1 + 2im) * X, (1.5 - 2im) * X, 1e16 * X)
        @test eval(Meta.parse(TMS.obs_name(op))) == op
    end
    # and so do channels and generators
    for op in (0.9 * Gate(Id) + 0.1 * Gate(X), Gate(X(1) + Y(2)), Reset("Up")(1),
               Reset([0.5 0. ; 0. 0.5])(2), Dephase(Z ⊗ Z)(1, 2), Dissipator(Reset(Level => 0)),
               -im * X(1) * X(2) + Dissipator(Sm)(1), 0.5 * Dissipator(Sm(1) + Sm(2)),
               exp(0.1 * (-im * (X ⊗ X) + Dissipator(Sm ⊗ Sm))), Gate(Id ⊗ Id))
        @test eval(Meta.parse(TMS.obs_name(op))) == op
    end
    @test TMS.obs_name(X ⊗ (Y * Z)) == "X⊗(Y*Z)"
    @test TMS.obs_name((X^0.5)^0.5) == "(X^0.5)^0.5"
    # an integer coefficient, stored as a float, prints as the integer it was, and so does a
    # float of integer value, as long as it stands for a single integer
    @test repr(2X) == repr(2.0X) == "2X"
    @test repr((1 + 2im) * X) == "(1 + 2im)X"
    @test repr(1e16 * X) == "1.0e16X"
end

@testset "Long sums and products print compactly" begin
    # the terms left out take the sign of the first of them
    s = X(1) + X(2) + X(3) - X(4) + sum(X(i) for i in 5:9) - X(10)
    @test TMS.obs_name(s) == "X(1)+X(2)+X(3)-...+X(9)-X(10)"
    @test eval(Meta.parse(repr(s))) == s
    @test TMS.obs_name(sum(X(i) for i in 1:10)) == "X(1)+X(2)+X(3)+...+X(9)+X(10)"
    @test TMS.obs_name(sum(X(i) for i in 1:6)) == "X(1)+X(2)+X(3)+X(4)+X(5)+X(6)"

    p = prod(Z(i) for i in 1:10)
    @test TMS.obs_name(p) == "Z(1)*Z(2)*Z(3)*...*Z(9)*Z(10)"
    @test eval(Meta.parse(repr(p))) == p
    @test TMS.obs_name(prod(Z(i) for i in 1:6)) == "Z(1)*Z(2)*Z(3)*Z(4)*Z(5)*Z(6)"
    # the factors kept are parenthesized as in the full form
    @test TMS.obs_name(prod(Z(i) for i in 1:6) * (X(1) + X(2))) ==
        "Z(1)*Z(2)*Z(3)*...*Z(6)*(X(1)+X(2))"
    @test TMS.obs_name(X(1) + 2p) == "X(1)+2Z(1)*Z(2)*Z(3)*...*Z(9)*Z(10)"
end

@testset "map_sites moves the factors of an operator" begin
    # compared as they are built, ≈ on operators needing simplify
    twice(i) = 2i
    @test map_sites(twice, X(1) * Y(2)) == X(2) * Y(4)
    @test map_sites(twice, 0.5 * Swap(1, 3) + Z(2) - 3 * Id(1)) == 0.5 * Swap(2, 6) + Z(4) - 3 * Id(1)
    # the factors keep their order, and so the sign of a product of fermionic operators
    @test map_sites(twice, C(3) * dag(C)(1)) == C(6) * dag(C)(2)
    @test map_sites(twice, C(3) * dag(C)(1)) ≠ dag(C)(2) * C(6)
    # a channel and a generator, with what they hold placed
    @test map_sites(twice, Gate(X(1) + Z(2)) * Gate(Y)(1)) == Gate(X(2) + Z(4)) * Gate(Y)(2)
    @test map_sites(twice, -im * Z(2) + Dissipator(X(1) + Y(3))) == -im * Z(4) + Dissipator(X(2) + Y(6))
    @test_throws "repeats a site" map_sites(i -> 1, Swap(1, 2))
end
