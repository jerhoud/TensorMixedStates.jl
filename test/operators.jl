# Operator algebra.
#
# Goes here: everything that manipulates operators without touching a state, that is
# simplify, removeMulti, dag and the fermionic classification. The invariant to prefer
# is that rewriting an operator must not change the matrix it stands for: it holds
# whatever normal form simplify settles on, so these tests survive changes in the
# simplifier that assertions on the printed form would not.

@testset "Operator simplification" begin
    q = Qubit()
    # simplify must never change the operator it stands for
    for a in [X * X, X * Y, X + Y, 2X, 0 * X, X - X, X^2, X^3, sqrt(X * X), (X + Y)^2,
              dag(S), dag(X * Y), dag(X + im * Y), dag(dag(X)), exp(im * X),
              Id * X, X * Id, (X + Y) * Z, (2X) * (3Y), H * H, S * S, T * T,
              Proj(0) + Proj(1)]
        @test matrix(simplify(a), q) ≈ matrix(a, q)
    end
    for a in [X ⊗ Y, Swap, Swap * Swap, (X ⊗ Y) * (Y ⊗ X), controlled(X), X ⊗ Id + Id ⊗ X]
        @test matrix(simplify(a), q, q) ≈ matrix(a, q, q)
    end
    # and a power of a superoperator, for which simplify had no method: an integer one is the
    # composition repeated, a non integer one is kept whole
    for a in [Left(X)^2, Right(S)^3, Gate(X)^2, Gate(H)^0.5, Dissipator(Sp)^2,
              (0.9 * Gate(Id) + 0.1 * Gate(X))^3]
        @test matrix(simplify(a), q) ≈ matrix(a, q)
    end
    # an F crossing a factor takes the sign of its parity, a sum of odd operators included,
    # and stops before a factor that has none
    f = Fermion()
    for a in [F * C, F * dag(C), F * (C + dag(C)), (C + dag(C)) * F, F * (2C + im * dag(C)),
              F * C * F, F * (C * N), F * (C + N), F * exp(N), F * C^3]
        @test matrix(simplify(a), f) ≈ matrix(a, f)
    end
    # the adjoint of a non integer power, which is not the power of the adjoint as soon as
    # the operator has a negative eigenvalue
    for a in [dag(sqrt(X)), dag(Sz^0.5), dag((X + Y)^0.3), dag(sqrt(N)), dag(X^3)]
        @test matrix(simplify(a), q) ≈ matrix(a, q)
    end
    # a few normal forms that the simplifier is expected to reach
    @test simplify(X * X) == Id
    @test simplify(X^2) == Id
    @test simplify(X * Id) == X
    # a run collapsing to the identity leaves its two neighbours adjacent, and they have to
    # be merged in their turn. The simplifier used to make a single pass over the factors
    # and stopped at X*X here, which also made it depend on how many times it was called
    @test simplify(X * Y * Y * X) == Id
    @test simplify(X * Y * Z * Z * Y * X) == Id
    @test simplify(Z * X * X * Z * Y) == Y
    @test simplify(simplify(H * X * Y * Y * X * H)) == simplify(H * X * Y * Y * X * H)
    # what must not move: two different bases on one site do not commute
    @test simplify(X * Y * X) == X * Y * X
    @test matrix(simplify(X - X), q) ≈ zeros(2, 2)
    @test matrix(simplify(0 * X), q) ≈ zeros(2, 2)
    # dag is an involution, which only shows on non hermitian operators
    for a in [S, T, S * T, 2im * S, X + im * Y, exp(im * X)]
        @test matrix(dag(dag(a)), q) ≈ matrix(a, q)
    end
    @test matrix(dag(dag(dag(S))), q) ≈ matrix(dag(S), q)
end

@testset "Fermionic classification" begin
    # a product of two fermionic operators is not fermionic, and the parity of the
    # number of fermionic factors is what decides
    for a in [C, dag(C), 2C, -C, C + C, C * N, C^3, dag(C) * C * dag(C), dag(dag(C))]
        @test isfermionic(a)
    end
    for a in [N, Id, dag(C) * C, C * C, N * N, C^2]
        @test !isfermionic(a)
    end
    # these combinations have no defined fermionic nature and must be rejected
    @test_throws ErrorException isfermionic(C + N)
    @test_throws ErrorException isfermionic(exp(C))
    # a term of coefficient zero, the zero operator, has every parity: it is left out of a sum
    # rather than making `0C + dag(C)` a mix of the two. A sum left with no term at all, as an
    # empty one, is the zero operator
    g = 0.0
    @test isfermionic(g * C + dag(C))
    @test !isfermionic(g * C + g * dag(C))
    @test matrix(g * C + g * dag(C), Fermion()) == zeros(2, 2)
    @test matrix(TensorMixedStates.SumOp(TensorMixedStates.Op{Pure, TensorMixedStates.Generic, 1}[]),
                 Fermion()) == zeros(2, 2)
    # has_fermionic asks the other question: whether a factor still needs its
    # Jordan-Wigner string. A product of two fermionic operators is not fermionic, but
    # both of its factors are, so it does.
    @test has_fermionic(dag(C)(1))
    @test has_fermionic(dag(C)(1) * C(3))
    @test has_fermionic(2 * C(2) + C(4))
    @test !has_fermionic(N(1))
    @test !has_fermionic(N(1) * N(2))
    @test !has_fermionic(Id(1))
    # simplify is what inserts the strings, and so what makes the answer false
    @test !has_fermionic(simplify(dag(C)(3)))
    @test !has_fermionic(simplify(dag(C)(1) * C(3)))
end

@testset "Multi_F removal" begin
    # simplify represents a Jordan-Wigner string of two sites or more by a single
    # Multi_F, removeMulti spells it out as one F per site
    op = simplify(dag(C)(1) * C(4))
    rm = TensorMixedStates.removeMulti(op)
    @test occursin("Multi_F", string(op))
    @test !occursin("Multi_F", string(rm))
    @test TensorMixedStates.removeMulti(rm) == rm
    # both forms must measure the same thing
    n = 5
    sys = System(n, Fermion())
    st = tdvp(-im * sum(dag(C)(i)C(i + 1) + dag(C)(i + 1)C(i) for i in 1:n - 1), 0.7,
              State{Pure}(sys, ["1", "0", "1", "0", "1"]); limits = Limits(maxdim = 32))
    for j in 2:n
        a = simplify(dag(C)(1) * C(j))
        @test expect(st, a) ≈ expect(st, TensorMixedStates.removeMulti(a))
    end
end

@testset "Global ordering of operators" begin
    # `isless(::Op, ::Op)` ranks the types first and, for two of the same rank, asks
    # `isless` again on the pair. A type with a ranking but no `isless` of its own
    # therefore recurses on itself instead of comparing anything
    @test isless(Dissipator(X), Dissipator(Y))
    @test !isless(Dissipator(Y), Dissipator(X))
    # a Proj holds a number, a name or a vector, which have no order between them either
    @test length(sort([Proj(1), Proj("Up"), Proj([1., 0.]), Proj(0)])) == 4
    @test matrix(simplify(Proj(1) + Proj("Up")), Qubit()) ≈ identity_operator(2)
    @test isless(Evolver(X(1)), Evolver(Y(1)))
    @test !isless(Evolver(Y(1)), Evolver(X(1)))
    # a SetState holds a name, a vector or a matrix, and two of them have to be ordered
    # whatever they hold: there is no order between those types, nor between two matrices
    @test isless(SetState("Dn"), SetState("Up"))
    @test length(sort([SetState("Up"), SetState([1., 0.]), SetState("Dn"),
                       SetState([1. 0. ; 0. 0.])])) == 4
    # the ranking puts the types in the order it declares
    @test isless(Id, X)                     # Identity before Operator
    @test isless(X(1), Gate(X)(1))          # AtIndex before Gate
    @test isless(Gate(X)(1), Dissipator(X)(1))
    @test isless(Dissipator(X)(1), Evolver(X(1)))
end

@testset "A projector is self adjoint" begin
    # `matrix` refuses a projector on a mixed state, so every one is self adjoint, and its
    # adjoint simplifies to it whatever it projects on. `measure` relies on this to find
    # that a projector gives real values
    for p in (Proj("Up"), Proj(1), Proj([1, im] / √2))
        @test simplify(dag(p)) == p
        @test matrix(dag(p), Qubit()) ≈ matrix(p, Qubit())
    end
end

@testset "Non integer powers of a scaled operator" begin
    # the power is the principal one. A positive coefficient comes out of it unchanged, any
    # other phase does not and stays inside, only the modulus coming out: a negative one is
    # the phase -1, and so is -1 + 0im. `PowOp` decides this when the power is built and
    # `simplify_pow` when a sum collapses into a scaled operator afterwards, and the two have
    # to reach the same form
    q = Qubit()
    for (c, a) in [(-1., -X - X), (im, im * X + im * X), (2., 2X + 2X),
                   (1. + im, (1 + im) * X + (1 + im) * X),
                   (-1. + 0im, (-1. + 0im) * X + (-1. + 0im) * X)]
        @test simplify(a^0.5) == (2c * X)^0.5
        @test matrix((2c * X)^0.5, q) ≈ sqrt(2c * matrix(X, q))
    end
    # the sign really does leave the coefficient for the operator
    @test string(simplify((-X - X)^0.5)) == string(sqrt(2.) * (-X)^0.5)
    # an integer exponent needs none of this care
    @test simplify((-X - X)^2) == 4Id
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
                   ((X * Y)^2, (X * Y)^2),                              # PowOp
                   (exp(X * Y), exp(X * Y)),                            # ExpOp
                   (parity(X * Y), parity(X * Y)),                      # ModOp
                   (Left(X * Y), Left(X * Y)),                          # Left
                   (Right(X * Y), Right(X * Y)),                        # Right
                   (Gate(X * Y), Gate(X * Y)),                          # Gate
                   (Dissipator(X * Y), Dissipator(X * Y)),              # Dissipator
                   (Evolver(X(1) * Y(2)), Evolver(X(1) * Y(2))),        # Evolver
                   (TensorMixedStates.Multi_F{Pure}(2, 4, false, false),   # Multi_F
                    TensorMixedStates.Multi_F{Pure}(2, 4, false, false)),
                   # an Operator built afresh each call, whose expr is a matrix for one
                   # and a whole expression for the other
                   (Phase(0.3), Phase(0.3)),
                   (controlled(Z), controlled(Z)),
                   (simplify(C(3) * dag(C)(1)), simplify(C(3) * dag(C)(1)))]
        @test a == b
        @test hash(a) == hash(b)
        @test length(Set([a, b])) == 1
    end
    # and operators that differ must keep differing
    @test length(Set([X(1) * Y(2), X(1) * Z(2)])) == 2
    # what it is for: a product shared by two measurements is asked for only once
    m = Measure([X(1) * Y(2) + Z(1), X(1) * Y(2) + Z(3)])
    @test length(Set(TensorMixedStates.get_prods(m))) == 3
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
        @test matrix(simplify(a), b) ≈ matrix(a, b)
    end

    # the adjoint puts the sign into the operator instead of leaving a dag outside
    @test matrix(simplify(dag(mod(N, 3))), b) ≈ adjoint(matrix(mod(N, 3), b))
    @test matrix(dag(dag(mod(N, 3))), b) ≈ matrix(mod(N, 3), b)

    # (-1)^N is hermitian, whatever form simplify settles on
    @test matrix(parity(N), b) ≈ adjoint(matrix(parity(N), b))

    @test_throws "a modulus is at least 2, got 1" mod(N, 1)
    @test_throws "exponentiates the fermionic operator" isfermionic(parity(C))
end

@testset "Renaming an operator" begin
    f = Fermion()
    # the definition is kept and only the label changes: a plain relabelling would send the
    # lookup after a name no site defines
    @test matrix(named(N, "Nf"), f) ≈ matrix(N, f)
    @test repr(named(N, "Nf")) == "Nf"
    @test repr(named(N, "Nf")(3)) == "Nf(3)"
    @test matrix(simplify(named(N, "Nf")), f) ≈ matrix(N, f)
    # an expression has no name of its own and can be renamed too
    @test matrix(named(2Sz, "SzA"), Spin(1/2)) ≈ matrix(2Sz, Spin(1/2))
    # a fermionic expression stays fermionic, which is what gives it its Jordan-Wigner string,
    # and one of no definite fermionic nature cannot be named at all
    @test named(2C, "C2").type == fermionic_op
    @test named(dag(C), "Cd").type == fermionic_op
    @test named(dag(C) * C, "n").type == plain_op
    @test_throws "cannot sum fermionic and non fermionic operators" named(C + N, "x")
end

@testset "Flux of an operator" begin
    IT = TensorMixedStates.ITensors

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
    # generates the very symmetry Zd records, would be refused on its wrap around
    @test flux(Zd, Qudit(3, conserve = Zd)) == IT.QN("Zd", 0, 3)
    @test flux(Xd, Qudit(3, conserve = Zd)) == IT.QN("Zd", 1, 3)
    @test flux(Xd^2, Qudit(3, conserve = Zd)) == IT.QN("Zd", 2, 3)
    @test flux(A, Boson(4, conserve = parity(N))) == IT.QN("parity(N)", 1, 2)

    # several charges at once
    @test flux(Cup, Electron(conserve = (Ntot, 2Sz))) == IT.QN(("Ntot", -1), ("2Sz", -1))
    @test flux(Ntot, Electron(conserve = (Ntot, 2Sz))) == IT.QN(("Ntot", 0), ("2Sz", 0))

    # an operator raising and lowering at once carries no flux, and the message says which
    # two differences it found
    @test_throws "differing by 2Sz=2 and by 2Sz=-2" flux(X, Qubit(conserve = 2Sz))
    @test_throws "no definite flux" flux(Hd, Qudit(3, conserve = Zd))
    # under parity the same X is fine, the two differences becoming one modulo 2
    @test flux(X, Qubit(conserve = parity(2Sz))) == IT.QN("parity(2Sz)", 0, 2)

    @test_throws "only defined for one site operators" flux(Swap, Qubit(conserve = 2Sz))
end

@testset "Superoperators in the dense basis" begin
    # without charges the mixed basis orders the ket fastest, which gives a superoperator an
    # exact form owing nothing to the way TMS assembles it
    q, f = Qubit(), Fermion()
    i2 = identity_operator(2)
    for op in [Sp, Sm, S, T, X, Z, H]
        mat = matrix(op, q)
        @test matrix(Left(op), q) ≈ kron(i2, mat)
        @test matrix(Right(op), q) ≈ kron(conj(mat), i2)
        @test matrix(Gate(op), q) ≈ kron(conj(mat), mat)
    end

    # a dissipator is what its definition says, AρA† - (A†Aρ + ρA†A)/2
    mat = matrix(C, f)
    aa = adjoint(mat) * mat
    @test matrix(Dissipator(C), f) ≈
        kron(conj(mat), mat) - 0.5 * (kron(i2, aa) + kron(conj(aa), i2))

    # and on a site of another dimension
    b = Boson(3)
    i3 = identity_operator(3)
    mb = matrix(A, b)
    @test matrix(Left(A), b) ≈ kron(i3, mb)
    @test matrix(Right(A), b) ≈ kron(conj(mb), i3)
end

@testset "A matrix knows no charge" begin
    # the matrix of an operator is written in the basis of its sites whatever they conserve.
    # A superoperator read off the mixed index of a charged site came out in the order the
    # charges sort its sectors, and an operator on several charged sites not at all
    strong = TensorMixedStates.strong
    f, fq, fs = Fermion(), Fermion(conserve = N), Fermion(conserve = strong(N))
    for op in (Left(C), Right(C), Gate(C), Dissipator(C), SetState("Occ"))
        @test matrix(op, fq) == matrix(op, f)
        @test matrix(op, fs) == matrix(op, f)
    end
    e, eq = Electron(), Electron(conserve = (Ntot, 2Sz))
    @test matrix(Ntot ⊗ Sz, eq, eq) == matrix(Ntot ⊗ Sz, e, e)
    b, bq = Boson(4), Boson(4, conserve = parity(N))
    @test matrix(Left(N ⊗ N), bq, bq) == matrix(Left(N ⊗ N), b, b)

    # an operator carrying no definite flux on a site still has a matrix: what is refused is
    # placing it on a system that conserves
    @test matrix(Left(X), Qubit(conserve = N)) == matrix(Left(X), Qubit())

    # a single site stands for as many identical ones as the operator needs
    q = Qubit()
    @test matrix(Left(Swap), q) == matrix(Left(Swap), q, q)
    @test matrix(Gate(Swap), q) == matrix(Gate(Swap), q, q)
end

@testset "Operators of several sites split on their sites" begin
    # an operator of several sites that simplify cannot develop, given with its sites, is
    # split into a sum of products of one site operators, which becomes its definition.
    # Splitting must not change the matrix it stands for, whatever it was given by
    TMS = TensorMixedStates
    q, s1 = Qubit(), Spin(1)
    msw = [1. 0 0 0 ; 0 0 1 0 ; 0 1 0 0 ; 0 0 0 1]
    sw = Operator{2}("Sw", msw, involution_op, q)
    @test matrix(sw, q) ≈ msw
    ss = Sx ⊗ Sx + Sy ⊗ Sy + Sz ⊗ Sz
    p2e = Operator{2}("P2e", ss / 2 + ss * ss / 6 + (Id ⊗ Id) / 3, selfadjoint_op)
    mp2 = matrix(p2e, s1)
    p2m = Operator{2}("P2m", mp2, selfadjoint_op)
    p2 = Operator{2}("P2", mp2, selfadjoint_op, s1)
    @test matrix(p2, s1) ≈ mp2
    # sites that differ, given in the order of the indices, and more than two of them
    mk = matrix(Sz ⊗ X + Sp ⊗ Z + 0.3 * Sx ⊗ Id, s1, q)
    @test matrix(Operator{2}("K", mk, plain_op, s1, q), s1, q) ≈ mk
    m3 = matrix(X ⊗ Z ⊗ Y + Z ⊗ Id ⊗ Z + 0.5 * Id ⊗ X ⊗ Id + 0.2 * Id ⊗ Id ⊗ Id, q)
    @test matrix(Operator{3}("T3", m3, plain_op, q), q) ≈ m3
    # a function of the sites, an expression simplify keeps whole, and a complex matrix
    zx = Operator{2}("ZX", (a, b) -> kron(matrix(Z, a), matrix(X, b)), plain_op, q)
    @test matrix(zx, q) ≈ kron(matrix(Z, q), matrix(X, q))
    rxx = exp(-0.3im * (X ⊗ X))
    @test matrix(Operator{2}("R", rxx, plain_op, q), q) ≈ matrix(rxx, q)
    @test matrix(Operator{2}("Y2", matrix(Y ⊗ Sp, q), plain_op, q), q) ≈ matrix(Y ⊗ Sp, q)

    # simplify develops it into factors of one site. They are named after the operator, the
    # index 0 being the part acting on a site alone, the others pairing across the sites
    one_site(a) = all(t -> all(f -> f isa TMS.AtIndex{Pure, 1}, TMS.prodsubs(t)),
                      TMS.sumsubs(a))
    @test one_site(simplify(sw(1, 2)))
    @test one_site(simplify(p2(3, 1)))
    terms(a) = [ TMS.scalararg(t).subs for t in TMS.sumsubs(a.expr) ]
    factor_names(a) = Set(f.name for t in terms(a) for f in t if f isa Operator)
    fe = Fermion()
    nn = Operator{2}("NN", matrix(N ⊗ N, fe), plain_op, fe)
    @test factor_names(nn) == Set(["NN¹₀", "NN²₀", "NN¹₁", "NN²₁"])
    @test all(t -> t[1].name[end] == t[2].name[end],
              filter(t -> all(f -> f isa Operator, t), terms(sw)))
    # the identity is taken out on each site before the split, so that the fewest terms
    # cross the link: three for Swap, as its expression has, and eight for the projector
    # of the AKLT chain, where its expression takes twelve
    crossing(a) = count(t -> count(f -> f isa Operator, t) == 2, terms(a))
    @test crossing(sw) == 3
    @test crossing(p2) == 8

    # it then goes into an MPO or expect as an operator given by an expression does, and
    # the MPO acts as the matrix it was split from
    chain(op) = sum(op(i, i + 1) for i in 1:3)
    st = RandomState{Pure}(System(4, q), 4)
    @test maxlinkdim(make_mpo(st, chain(sw))) == maxlinkdim(make_mpo(st, chain(Swap)))
    @test norm(apply(make_mpo(st, chain(sw)), st) - apply(make_mpo(st, chain(Swap)), st)) < 1e-12
    st1 = RandomState{Pure}(System(4, s1), 4)
    @test maxlinkdim(make_mpo(st1, chain(p2))) == 10
    @test maxlinkdim(make_mpo(st1, chain(p2e))) == 14
    @test norm(apply(make_mpo(st1, chain(p2)), st1) - apply(make_mpo(st1, chain(p2e)), st1)) < 1e-12
    # on sites that are not neighbours, or given in the other order
    for (i, j) in [(1, 3), (4, 2), (2, 1)]
        @test expect(st1, p2(i, j)) ≈ expect(st1, p2e(i, j))
        @test norm(apply(p2m(i, j), st1) - apply(make_mpo(st1, p2(i, j)), st1)) < 1e-12
    end
    # sites that differ, placed on a system that mixes them
    k = Operator{2}("K", mk, plain_op, s1, q)
    km = Operator{2}("Km", mk, plain_op)
    sk = RandomState{Pure}(System([s1, q, s1, q]), 4)
    for (i, j) in [(1, 2), (3, 2), (1, 4)]
        @test norm(apply(km(i, j), sk) - apply(make_mpo(sk, k(i, j)), sk)) < 1e-12
    end
    # a density matrix, with the operator in a gate, a dissipator and a hamiltonian. The MPOs
    # are compared whole on three sites: their links are near a hundred channels wide, and
    # applying them to a random density matrix takes gigabytes
    ρ = mix(State{Pure}(System(3, s1), "0"))
    for p in (x -> Gate(x)(1, 3), x -> Dissipator(x)(2, 3), x -> -im * x(1, 2))
        @test norm(prod(make_mpo(ρ, p(p2))) - prod(make_mpo(ρ, p(p2e)))) < 1e-12
    end

    # on sites that conserve, every factor carries a definite charge, which the expression
    # written with Sx and Sy does not, and the operator acts on the same sites weakened
    s1q = Spin(1, conserve = Sz)
    p2q = Operator{2}("P2q", mp2, selfadjoint_op, s1q)
    conf(sys) = normalize(State{Pure}(sys, ["1", "0", "-1", "0"]) +
                          0.5 * State{Pure}(sys, ["0", "1", "0", "-1"]))
    stq, std = conf(System(4, s1q)), conf(System(4, s1))
    @test expect(stq, chain(p2q)) ≈ expect(std, chain(p2e))
    @test expect(weaken(stq, ()), chain(p2q)) ≈ expect(std, chain(p2e))
    @test maxlinkdim(make_mpo(stq, chain(p2q))) == 10
    # the factors split charge by charge recombine into the matrix they came from, a charge
    # taken modulo included
    @test matrix(p2q, s1q) ≈ mp2
    bp = Boson(4, conserve = parity(N))
    mb = matrix(A ⊗ dag(A) + dag(A) ⊗ A + 0.5 * (A * A) ⊗ N, Boson(4))
    @test matrix(Operator{2}("Bp", mb, plain_op, bp), bp) ≈ mb
    @test_throws "no definite flux" make_mpo(stq, chain(p2e))
    @test_throws "no definite charge of N" Operator{2}("XX", matrix(X ⊗ X, q), plain_op,
                                                       Qubit(conserve = N))

    # a matrix is taken as it is, with no Jordan-Wigner string, so on a fermionic site it has
    # to be even: its factors are then crossed by the strings of other operators with no sign
    nne = Operator{2}("NNe", N ⊗ N, plain_op)
    sf = RandomState{Pure}(System(4, fe), 4)
    for p in (x -> x(1, 3), x -> C(4) * x(1, 3) * dag(C)(2), x -> dag(C)(1) * x(4, 2) * C(3))
        @test expect(sf, p(nn)) ≈ expect(sf, p(nne))
    end
    hop = matrix(C ⊗ dag(C), fe)
    @test_throws "does not commute with F on its site 1" Operator{2}("H", hop, plain_op, fe)

    # a single site has nothing to split, its definition is only computed on it once, and
    # its type is kept, fermionic included
    ex = Operator{1}("Ex", exp(0.3im * X), plain_op, q)
    @test ex.expr isa Matrix && ex.expr ≈ matrix(exp(0.3im * X), q)
    @test Operator{1}("Fs", s -> matrix(Sz, s), selfadjoint_op, s1).expr ≈ matrix(Sz, s1)
    cf = Operator{1}("Cf", C, fermionic_op, fe)
    @test isfermionic(cf) && matrix(cf, fe) ≈ matrix(C, fe)
    @test_throws "has dimension 3" Operator{1}("A", [1. 0 ; 0 1], plain_op, s1)

    # what cannot be split
    @test_throws "acts on 2 sites and was given 3" Operator{2}("A", msw, plain_op, q, q, q)
    @test_throws "has dimension 6" Operator{2}("A", msw, plain_op, s1, q)
end
