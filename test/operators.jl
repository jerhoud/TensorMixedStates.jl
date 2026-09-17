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
    # a few normal forms that the simplifier is expected to reach
    @test simplify(X * X) == Id
    @test simplify(X^2) == Id
    @test simplify(X * Id) == X
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
end

@testset "Multi_F removal" begin
    # simplify represents a Jordan-Wigner string of two sites or more by a single
    # Multi_F, removeMulti spells it out as one F per site
    op = simplify(dag(C)(1) * C(4))
    rm = removeMulti(op)
    @test occursin("Multi_F", string(op))
    @test !occursin("Multi_F", string(rm))
    @test removeMulti(rm) == rm
    # both forms must measure the same thing
    n = 5
    sys = System(n, Fermion())
    st = tdvp(-im * sum(dag(C)(i)C(i + 1) + dag(C)(i + 1)C(i) for i in 1:n - 1), 0.7,
              State{Pure}(sys, ["1", "0", "1", "0", "1"]); limits = Limits(maxdim = 32))
    for j in 2:n
        a = simplify(dag(C)(1) * C(j))
        @test expect(st, a) ≈ expect(st, removeMulti(a))
    end
end
