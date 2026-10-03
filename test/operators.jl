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
    @test length(Set(TensorMixedStates.get_prods(State{Pure}(System(3, Qubit()), "Up"), m))) == 3
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
    @test_throws "parity(C) has no definite fermionic parity" isfermionic(parity(C))
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
    @test_throws "C+N has no definite fermionic parity" named(C + N, "x")
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
    # the identity, which has no site, counts as one site
    one_site(a) = all(t -> all(f -> f isa TMS.AtIndex{Pure, 1} || f isa TMS.IdentityOp,
                               TMS.prodsubs(t)), TMS.sumsubs(a))
    @test one_site(simplify(sw(1, 2)))
    @test one_site(simplify(p2(3, 1)))
    terms(a) = [ TMS.scalararg(t).subs for t in TMS.sumsubs(a.expr)
                 if !(TMS.scalararg(t) isa TMS.IdentityOp) ]
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
    # its eight terms give the MPO its eight channels, which the twelve of the expression come
    # down to once compacted
    @test maxlinkdim(make_mpo(st1, chain(p2))) == 10
    @test maxlinkdim(make_mpo(st1, chain(p2e))) == 10
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
    @test_throws "whose dimension is 3" Operator{1}("A", [1. 0 ; 0 1], plain_op, s1)

    # what cannot be split
    @test_throws "acts on 2 sites and was given 3" Operator{2}("A", msw, plain_op, q, q, q)
    @test_throws "whose dimension is 6" Operator{2}("A", msw, plain_op, s1, q)
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

    # an operator of several sites given by such an expression and its sites is the operator
    # the expression stands for
    h = dag(C) ⊗ C + dag(dag(C) ⊗ C)
    hm = c1' * c2 + c2' * c1
    @test matrix(h, fe) ≈ hm
    @test matrix(Operator{2}("H2", h * h, selfadjoint_op, fe), fe) ≈ hm^2
    e = Operator{2}("E", exp(-0.3 * (h * h)), plain_op, fe)
    @test matrix(e, fe) ≈ exp(-0.3 * hm^2)
    st = State{Pure}(System(3, fe), ["Occ", "Emp", "Emp"])
    @test expect(st, e(1, 2)) ≈ exp(-0.3)
    @test expect(st, e(1, 3)) ≈ exp(-0.3)

    # a factor of several sites has the parity of its pieces placed, a renamed expression as
    # well as a tensor product. Both were taken as even, which lost the sign of the adjoint and
    # the strings of the matrix, and made a dissipator change the trace
    cn = named(C ⊗ N, "CN")
    a2 = C ⊗ Id + Id ⊗ C
    ψ = RandomState{Pure}(System(4, fe), 4)
    for a in (cn ⊗ cn, a2 ⊗ a2)
        aψ = apply(make_mpo(ψ, a(1, 2, 3, 4)), ψ)
        @test expect(ψ, dag(a)(1, 2, 3, 4) * a(1, 2, 3, 4)) ≈ norm(aψ)^2
    end
    @test matrix(C ⊗ cn, fe) ≈ matrix(C ⊗ C ⊗ N, fe)
    @test matrix(C ⊗ a2, fe) ≈ matrix(C ⊗ C ⊗ Id + C ⊗ Id ⊗ C, fe)
    @test_throws "no definite fermionic parity" matrix(C ⊗ (C ⊗ Id + N ⊗ Id), fe)
    ρ = mix(RandomState{Pure}(System(4, fe), 4))
    @test abs(trace(apply(make_mpo(ρ, Dissipator(cn ⊗ cn)(1, 2, 3, 4)), ρ))) < 1e-12
end

@testset "Products differing by their last factor are gathered" begin
    # P*X(i) + P*Y(i) is P*(X+Y)(i): the terms of an odd sum of one site share their
    # Jordan-Wigner string, and an MPO laid term by term carries each product once. The metric
    # is the number of factors simplify leaves, the terms of such an MPO, measured on sums
    # written both ways: PreMPO compacts, which would gather the products whatever simplify
    # did. The operator must not change
    TMS = TensorMixedStates
    factors(a) = sum(t -> count(f -> !(f isa TMS.IdentityOp), TMS.prodsubs(TMS.scalararg(t))),
                     TMS.sumsubs(TMS.removeMulti(simplify(a))))
    fe = Fermion()
    st = RandomState{Pure}(System(5, fe), 2)
    @test factors(sum((C + dag(C))(i) for i in 1:5)) == 15
    @test factors(sum(C(i) + dag(C)(i) for i in 1:5)) == 15
    # compacted, the strings share a single channel
    @test maxlinkdim(make_mpo(st, sum((C + dag(C))(i) for i in 1:5))) == 3
    @test factors(sum(Dissipator(C + dag(C))(i) for i in 1:4)) == 13
    q = RandomState{Pure}(System(3, Qubit()), 2)
    a = X(1) * Z(2) + X(1) * Y(2)
    @test factors(a) == 2
    @test norm(apply(make_mpo(q, a), q) - apply(make_mpo(q, X(1) * Z(2)), q) -
               apply(make_mpo(q, X(1) * Y(2)), q)) < 1e-12
end

@testset "Powers" begin
    # a literal exponent goes through literal_pow, which called inv for a negative one
    p = -1
    @test X^-1 == X^p
    @test (Id + 0.5X)^-1 == (Id + 0.5X)^p
    # an integer power is a product, with the adjoint, the parity and the strings of that
    # product, an involution reduced at once; any other is the principal power, a function of
    # the operator as exp is, of which only the modulus of a coefficient comes out
    q, fe, s1 = Qubit(), Fermion(), Spin(1)
    @test X^10000 == Id
    @test X^10001 == X
    @test (2X)^3 == 8X
    @test simplify(Sz^2 * Sz) == Sz^3
    for (a, s) in [(X^3, q), ((X + Y)^3, q), (Sz^2, s1), ((Sp + Sm)^3, s1),
                   ((C + dag(C))^2, fe), (C^3, fe), (dag((X + im * Y)^2), q)]
        @test matrix(simplify(a), s) ≈ matrix(a, s)
    end
    @test matrix(dag((X + im * Y)^2), q) ≈ matrix((X + im * Y)^2, q)'
    @test isfermionic(C^3)
    @test !isfermionic(C^2)

    # a phase stays inside: the principal square root of -Z has i on its eigenvalue -1
    @test matrix(sqrt(-Z), q) ≈ [im 0 ; 0 1]
    up = State{Pure}(System(2, q), "Up")
    @test expect(up, sqrt(-Z)(1)) ≈ im
    @test expect(up, ((im * X)^0.5)(1)) ≈ matrix((im * X)^0.5, q)[1, 1]
    @test matrix(X^(0.5im), q) ≈ matrix(X, q)^(0.5im)
    # a power that does not exist is refused, where Julia gives zero for C^0.5
    @test_throws "does not exist" matrix(sqrt(C), fe)
    @test_throws "does not exist" matrix(Proj(0)^(-0.5), q)

    # powers of the same operator merge, and what they give leaves no coefficient inside the
    # product nor an F out of its place
    @test simplify((-X)^0.5 * (-X)^0.5) == -X
    @test simplify(X^0.5 * X^0.5) == X
    @test simplify(F^0.5 * F^0.5 * C) == simplify(F * C)
    # and placed, a merge giving a multiple of the identity leaves only its coefficient: -Id
    # stayed in the product, which the sort by site then failed on
    a = ((im * X)^0.5)(1) * ((im * X)^1.5)(1)
    @test simplify(a * Z(2)) == simplify(Z(2) * a) == simplify(-Z(2))

    # a placed operator takes powers, an integer one of anything placed and any one of an
    # operator on its sites
    st = RandomState{Pure}(System(3, q), 4)
    @test expect(st, X(1)^2) ≈ 1
    @test expect(st, (X(1) + Z(2))^2) ≈ expect(st, (X(1) + Z(2)) * (X(1) + Z(2)))
    @test expect(st, X(1)^0.5) ≈ expect(st, (X^0.5)(1))
    @test_throws "has no power 0.5 once placed" (X(1) + Z(2))^0.5

    # a non integer power of several sites is kept whole, as exp is: it is refused where it
    # would have to be split, and applied as a gate, a fermion beside it included
    @test_throws "acts on several sites" expect(st, sqrt(Swap)(1, 2))
    sf = State{Pure}(System([q, q, fe]), ["Up", "Dn", "Emp"])
    @test_ok apply(sqrt(Swap)(1, 2) * dag(C)(3), sf)
    # a string is not moved past a factor of several sites kept whole, which it commutes with
    # only if the factor is even: this one, odd on site 1, used to compare equal to the product
    # in the other order, its opposite. A factor with no string in the way is still sorted
    xi = Operator{2}("XI", kron([0. 1. ; 1. 0.], [1. 0. ; 0. 1.]), plain_op)
    @test !(xi(1, 2) * dag(C)(4) ≈ dag(C)(4) * xi(1, 2))
    @test xi(1, 2) * N(4) ≈ N(4) * xi(1, 2)
end

@testset "The identity is one value of each kind" begin
    # every construction of an identity gives the same value, told by its type alone: on
    # several sites, on a density matrix, and placed, where it has no site, being the identity
    # of the whole system
    q = Qubit()
    @test Id ⊗ Id isa TensorMixedStates.IdentityOp
    @test Left(Id) == Right(Id) == Gate(Id)
    @test Id(3) == Id(1)
    @test (Id ⊗ Id)(1, 3) == Id(2)
    @test (Id ⊗ Id)^3 == Id ⊗ Id
    @test Left(Id)^2 == Left(Id)
    @test matrix((-(Id ⊗ Id))^0.5, q, q) ≈ matrix(im * (Id ⊗ Id), q, q)
    @test X(3)^0 == Id(1)
    @test Id(2)^0.5 == Id(1)
    @test iszero(TensorMixedStates.scalarcoef(simplify(Dissipator(Id)(1))))
    @test iszero(TensorMixedStates.scalarcoef(simplify(Evolver(-im * Id(1)))))
    st = RandomState{Pure}(System(3, q), 2)
    ρ = mix(st)
    @test expect(st, Id(2) - X(2) * X(2)) ≈ 0 atol = 1e-12
    # it is its own adjoint, and so measured as real, a number as any expectation value is
    @test last(only(measure(st, Id(2)))) === 1.0
    # its tensor, where one is needed, goes on the site the caller names
    @test expect1(st, Id) ≈ ones(3)
    @test expect1(ρ, Id) ≈ ones(3)
    @test expect2(st, (Id, X))[1, 3] ≈ expect(st, X(3))
    h = X(1) * X(2) + 0.5 * Id(3)
    @test real(inner(st.state', make_mpo(st, h), st.state)) / real(inner(st.state, st.state)) ≈ real(expect(st, h))
    @test norm(apply(2Id(2) * X(1), st) - 2 * apply(X(1), st)) < 1e-12
    @test trace(apply(Gate(Id)(2), ρ)) ≈ 1
    # having no site, the identity is the identity of any system, whatever site it was given
    @test expect(st, Id(5)) ≈ 1
end


@testset "Operators of several sites that were refused" begin
    # a renamed operator is its definition, developed as the definition itself is
    q = Qubit()
    st = State{Pure}(System(2, q), ["Up", "Dn"])
    @test expect(st, named(exp(-0.3im * Swap), "R")(1, 2)) ≈ expect(st, exp(-0.3im * Swap)(1, 2))
    # one site given for identical ones, a function of the sites included
    zx = Operator{2}("ZX", (a, b) -> kron(matrix(Z, a), matrix(X, b)), plain_op)
    @test matrix(zx, q) ≈ kron(matrix(Z, q), matrix(X, q))
    # a fermionic operator controlled is plain, and is its expression
    fe = Fermion()
    @test matrix(controlled(C), q, fe) ≈ matrix(Proj(0) ⊗ Id + Proj(1) ⊗ C, q, fe)
end

@testset "The tensor of a matrix" begin
    IT = TensorMixedStates.ITensors
    # one site given for several identical ones is laid on as many indices of that site,
    # charges included: it was laid on a bare index, which a system that conserves refused
    q, qn = Qubit(), Qubit(conserve = N)
    m = matrix(Sp ⊗ Sm, q, q)
    @test IT.dims(tensor(m, q)) == (4, 4)
    @test IT.flux(tensor(m, qn)) == IT.QN("N", 0)
    @test IT.array(IT.dense(tensor(m, qn))) ≈ IT.array(IT.dense(tensor(Sp ⊗ Sm, qn, qn)))
    @test_throws "no definite charge" tensor(matrix(X, q), qn)
    # a size matching no number of sites was taken for as many basis states of a new index
    @test_throws "acts on no number of Qubit()" tensor(rand(3, 3), q)
    @test_throws "whose dimension is 4" tensor(rand(2, 2), q, q)
end

@testset "The type of an operator is checked against its matrix" begin
    # simplify reasons with the type, so a type the matrix belies gave a wrong result with no
    # message: this involution squared to Id where its square is zero
    q, f = Qubit(), Fermion()
    @test_throws "declared involution_op but is not self adjoint on Qubit()" matrix(
        Operator{1}("Bad", [0. 1.; 0. 0.], involution_op), q)
    @test_throws "its square is not the identity" matrix(
        Operator{1}("TwoX", [0. 2.; 2. 0.], involution_op), q)
    @test_throws "declared selfadjoint_op but is not self adjoint" matrix(
        Operator{1}("A", [0. 1.; 0. 0.], selfadjoint_op), q)
    # an odd operator declared plain crosses F with no sign, and a fermionic one on a site
    # with no F with one
    k = Operator{1}("K", [0. 1.; 0. 0.], plain_op)
    @test_throws "declared plain_op but does not commute with F on Fermion()" expect(
        State{Pure}(System(2, f), "Emp"), k(2))
    @test_throws "declared fermionic_op but does not anticommute with F on Qubit()" matrix(
        Operator{1}("Cq", [0. 1.; 0. 0.], fermionic_op), q)
    # F is what every parity is read against, and what simplify squares to Id
    for s in (q, f, Electron(), Tj(), Boson(3))
        @test matrix(F, s)^2 ≈ matrix(Id, s)
    end
    # named checks a matrix at once, simplify reading its type before any matrix is computed:
    # W(2)^2 was measured 1 and V(2) - dag(V)(2) simplified to zero
    @test_throws "W is declared involution_op but is not self adjoint" named(
        [0. 1.; 0. 0.], "W"; type = involution_op)
    @test_throws "V is declared selfadjoint_op but is not self adjoint" named(
        [0. 1.; 0. 0.], "V"; type = selfadjoint_op)
    @test named([0. 1.; 1. 0.], "Xn").type == involution_op
    # checked when the operator is built on its sites as well
    @test_throws "declared selfadjoint_op but is not self adjoint" Operator{2}("P",
        kron([0. 1.; 0. 0.], [1. 0.; 0. 1.]), selfadjoint_op, q)
    # every operator of the libraries has the type its matrices say, on sites of every size.
    # What a site type declares are the methods of operator_definition for that type
    declared(t) = [ string(m.sig.parameters[3].parameters[1])
                    for m in methods(TensorMixedStates.operator_definition)
                    if m.sig.parameters[2] == t ]
    for s in [Qubit(), Fermion(), Electron(), Tj(), Spin(0), Spin(1/2), Spin(1), Qudit(1),
              Qudit(2), Qudit(3), Boson(2), Boson(4), Qboson(1., 2), Qboson(0.5, 4)]
        for name in declared(typeof(s))
            if name ≠ "F" && isdefined(TensorMixedStates, Symbol(name))
                o = getfield(TensorMixedStates, Symbol(name))
                if o isa Operator
                    @test_ok matrix(o, s)
                end
            end
        end
    end
end

@testset "named defines an operator" begin
    # without sites, a matrix is read on its own
    @test named([1 1; 1 -1] / √2, "MyH").type == involution_op
    @test named([0. 2.; 2. 0.], "TwoX").type == selfadjoint_op
    @test named([0. 1.; 0. 0.], "Up").type == plain_op
    @test named(s -> [0. 1.; 1. 0.], "Fx").type == plain_op
    # an expression keeps the type of what it renames
    @test named(N, "Nf").type == selfadjoint_op
    @test named(C, "D").type == fermionic_op
    # given its sites, the type is read off the matrix there, F included
    @test named(Sp + Sm, "Sx2", Spin(1)).type == selfadjoint_op
    @test matrix(named(Sp + Sm, "Sx2", Spin(1)), Spin(1)) ≈ 2 * matrix(Sx, Spin(1))
    @test named(s -> [0. 1.; 1. 0.], "Fx", Qubit()).type == involution_op
    @test named([0. 1.; 1. 0.], "Xf", Fermion()).type == fermionic_op
    @test_throws "M has no definite fermionic parity on Fermion()" named([1. 1.; 1. 1.], "M",
                                                                           Fermion())
    # a single site stands for as many as the size of the matrix asks for, and the operator
    # can then be measured
    sw = matrix(Swap, Qubit())
    s1 = named(sw, "MySwap", Qubit())
    @test s1 isa Operator{2}
    @test s1.type == involution_op
    @test matrix(s1, Qubit()) ≈ sw
    st = State{Pure}(System(2, Qubit()), ["Up", "X+"])
    @test expect(st, s1(1, 2)) ≈ 0.5
    @test_throws "a 3×3 matrix acts on no number of Qubit()" named(rand(3, 3), "R", Qubit())
    # a type given is kept, and checked
    @test named(X, "X2", Qubit(); type = plain_op).type == plain_op
    @test_throws "declared selfadjoint_op but is not self adjoint" named(Sp, "Sp2", Qubit();
                                                                          type = selfadjoint_op)
end

@testset "Signed zeros" begin
    # -0.0 and 0.0 are equal to `==` but not to `isless` or `hash`: simplify, which sorts the
    # terms before merging the equal ones, left a pair apart when a third term sorted between
    # them, and measure computed an expectation value twice
    c1 = -1 * (0.0 - 1.0im)       # -0.0 + 1.0im
    c3 = -1 * (0.0 - 2.0im)       # -0.0 + 2.0im
    f(c) = exp(c * X ⊗ X)(1, 2)
    @test simplify(f(c1) - f(1.0im) + f(c3)) == simplify(f(c3))
    a = simplify((-0.5im * X(1)) * (0.5im * Y(2)))
    b = simplify(0.25 * X(1) * Y(2))
    @test hash(a) == hash(b)
    @test length(Set([a, b])) == 1
    @test isequal(((X + Z)^complex(-0.0, 0.5)).expo, 0.5im)
    # and so in every number an operator stores: the matrix of an Operator, the state of a
    # Proj or of a SetState
    m(x) = Operator{1}("M", [1.0 x; 0.0 1.0], plain_op)
    @test hash(m(-0.0)) == hash(m(0.0))
    p, q = Proj([1.0, -0.0]), Proj([1.0, 0.0])
    @test hash(p) == hash(q)
    @test !(isless(p, q) || isless(q, p))
    @test hash(SetState([1.0 -0.0; 0.0 0.0])) == hash(SetState([1.0 0.0; 0.0 0.0]))
    # a projector is ordered by its values, which tie where `==` holds: its printed form told
    # [1, 0] and [1.0, 0.0] apart, and [1, 1] sorted between the two
    g(v) = exp(im * Proj(v) ⊗ X)(1, 2)
    @test simplify(g([1, 0]) - g([1.0, 0.0]) + g([1, 1])) == simplify(g([1, 1]))
    @test isless(Proj(1), Proj("Up")) && isless(Proj("Up"), Proj([1, 0]))
end

@testset "Printing reads back as the operator" begin
    # a name is a column header: it has to denote the operator measured. The right operand of
    # * and ⊗ and the base of ^ lost their parentheses, and so did rational numbers
    for op in (X ⊗ (Y * Z), (Y * Z) ⊗ X, (X ⊗ Y) * (Z ⊗ Z), (X^0.5)^0.5, (1//2) * X,
               X^(1//2), (-1.0 + 0im) * X, 0.5im * X, (X + Z)^0.5, X(1) * Y(2))
        @test eval(Meta.parse(TensorMixedStates.obs_name(op))) == op
    end
    @test TensorMixedStates.obs_name(X ⊗ (Y * Z)) == "X⊗(Y*Z)"
    @test TensorMixedStates.obs_name((X^0.5)^0.5) == "(X^0.5)^0.5"
end

@testset "Long sums and products print compactly" begin
    # the terms left out take the sign of the first of them
    s = X(1) + X(2) + X(3) - X(4) + sum(X(i) for i in 5:9) - X(10)
    @test TensorMixedStates.obs_name(s) == "X(1)+X(2)+X(3)-...+X(9)-X(10)"
    @test eval(Meta.parse(repr(s))) == s
    @test TensorMixedStates.obs_name(sum(X(i) for i in 1:10)) == "X(1)+X(2)+X(3)+...+X(9)+X(10)"
    @test TensorMixedStates.obs_name(sum(X(i) for i in 1:6)) == "X(1)+X(2)+X(3)+X(4)+X(5)+X(6)"

    p = prod(Z(i) for i in 1:10)
    @test TensorMixedStates.obs_name(p) == "Z(1)*Z(2)*Z(3)*...*Z(9)*Z(10)"
    @test eval(Meta.parse(repr(p))) == p
    @test TensorMixedStates.obs_name(prod(Z(i) for i in 1:6)) == "Z(1)*Z(2)*Z(3)*Z(4)*Z(5)*Z(6)"
    # the factors kept are parenthesized as in the full form
    @test TensorMixedStates.obs_name(prod(Z(i) for i in 1:6) * (X(1) + X(2))) ==
        "Z(1)*Z(2)*Z(3)*...*Z(6)*(X(1)+X(2))"
    @test TensorMixedStates.obs_name(X(1) + 2p) == "X(1)+2Z(1)*Z(2)*Z(3)*...*Z(9)*Z(10)"
end


@testset "map_sites moves the factors of an operator" begin
    twice(i) = 2i
    @test map_sites(twice, X(1) * Y(2)) ≈ X(2) * Y(4)
    @test map_sites(twice, 0.5 * Swap(1, 3) + Z(2) - 3 * Id(1)) ≈ 0.5 * Swap(2, 6) + Z(4) - 3 * Id(1)
    # the factors keep their order, and so the sign of a product of fermionic operators
    @test map_sites(twice, C(3) * dag(C)(1)) ≈ C(6) * dag(C)(2)
    @test !(map_sites(twice, C(3) * dag(C)(1)) ≈ dag(C)(2) * C(6))
    @test map_sites(twice, Gate(X)(1) - im * Z(2)) ≈ Gate(X)(2) - im * Z(4)
    @test_throws "repeats a site" map_sites(i -> 1, Swap(1, 2))
end
