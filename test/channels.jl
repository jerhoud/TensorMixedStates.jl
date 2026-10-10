# Channels and generators: their normal forms, their algebra and what it refuses, their matrices
# on the density matrices of sites, their purity, and the library of noises.
#
# The invariant to prefer is that the matrix of a channel or a generator is the map it stands
# for, written on vec(ρ) in an order spelled out here element by element, so that the tests do
# not depend on how TMS lays it out. Their descent onto a system, their MPO and their
# application to a state come with the later steps of layer 2.

"""
the identity matrix of dimension `d`
"""
eye(d) = Matrix{Float64}(TMS.LinearAlgebra.I, d, d)

"""
vec(ρ) of a density matrix of two sites of dimensions `d1` and `d2` in the order of the matrices
of channels: the last site fastest and, on each site, the ket faster than the bra
"""
mixed_vec(ρ, d1, d2) =
    vec([ ρ[(k1 - 1) * d2 + k2, (b1 - 1) * d2 + b2] for k2 in 1:d2, b2 in 1:d2, k1 in 1:d1, b1 in 1:d1 ])

@testset "Channels and generators in their normal form" begin
    # a gate and a dissipator take the modulus squared of a coefficient, a gate distributes over
    # a product, and placing what they hold is placing them
    @test Gate(X(1)) == Gate(X)(1)
    @test Gate((1 + im) * X) == 2 * Gate(X)
    @test Gate(X * Y) == Gate(X) * Gate(Y)
    @test Gate(2X(1) * Y(2)) == 4 * Gate(X)(1) * Gate(Y)(2)
    @test Gate(X(1) + Y(2)) isa Gate{TMS.Indexed}
    @test Dissipator(Sm(1)) == Dissipator(Sm)(1)
    @test Dissipator(2im * Sm) == 4 * Dissipator(Sm)
    # the dissipator of a gate is that of its operator, that of a weighted channel weighted
    @test Dissipator(Gate(Sm)) == Dissipator(Sm)
    @test Dissipator(Gate(Sm)(1)) == Dissipator(Sm)(1)
    @test Dissipator(0.3 * Reset("Up")) == 0.3 * Dissipator(Reset("Up"))
    # a real coefficient comes out of an evolver, a complex one stays inside
    @test Evolver(X(1)) == Evolver(X)(1)
    @test Evolver(-im * X(1)) == Evolver(-im * X)(1)
    @test Evolver(-2X) == -2 * Evolver(X)
    @test TMS.scalarcoef(Evolver(-im * X)) == 1
    @test Evolver(-im * X).arg == -im * X
end

@testset "An operator beside a channel or a generator" begin
    # added to a generator, an operator is its Hamiltonian part; multiplied by a channel or in a
    # tensor product with one, its gate
    @test -im * X(1) + Dissipator(Sm)(1) == Evolver(-im * X(1)) + Dissipator(Sm)(1)
    @test Dissipator(Sm)(1) - X(1) == Dissipator(Sm)(1) - Evolver(X)(1)
    @test X(1) * Gate(Y)(2) == Gate(X)(1) * Gate(Y)(2)
    @test Gate(Y)(2) * X(1) == Gate(Y)(2) * Gate(X)(1)
    @test X ⊗ Reset("Up") == Gate(X) ⊗ Reset("Up")
    @test Reset("Up") ⊗ X == Reset("Up") ⊗ Gate(X)
    # nothing else mixes the kinds
    @test_throws MethodError X(1) + Gate(Y)(1)
    @test_throws MethodError Gate(X)(1) + Dissipator(Y)(1)
    @test_throws MethodError X(1) * Dissipator(Y)(1)
end

@testset "Weights of channels and coefficients of generators" begin
    # a channel is weighted by real numbers of zero or above, a generator by real numbers
    @test_throws "a channel takes weights of zero or above, not -1" -Gate(X)
    @test_throws "a channel takes weights of zero or above, not -1" Gate(Id) - Gate(X)
    @test_throws "a channel takes weights of zero or above, not im" im * Gate(X)
    @test_throws "a generator takes real coefficients, not im" im * Dissipator(X)
    @test_throws "a generator takes real coefficients" -im * (X(1) + Dissipator(Sm)(1))
    @test (1 + 0im) * Gate(X) == Gate(X)
    @test TMS.scalarcoef(-Dissipator(X)) == -1
    @test iszero(TMS.scalarcoef(0im * Dissipator(X)))
    @test iszero(TMS.scalarcoef(0 * Gate(X)))
end

@testset "Products and powers of channels and generators" begin
    # generators are added, never multiplied, and a channel takes integer powers only
    @test_throws "generators are added, not multiplied" Dissipator(X)(1) * Dissipator(Y)(1)
    @test_throws "generators are added, not multiplied" Dissipator(X) ⊗ Dissipator(Y)
    @test_throws "has no power" Dissipator(X)^2
    @test_throws "has no power" Dissipator(X)(1)^2
    @test_throws "only the integer powers of a channel are channels" Gate(X)^0.5
    @test_throws "only the integer powers of a channel are channels" sqrt(Gate(X))
    @test_throws "only the integer powers of a channel are channels" Gate(X)(1)^0.5
    q = Qubit()
    g = 0.9 * Gate(Id) + 0.1 * Gate(X)
    @test matrix(Gate(X)^2, q) ≈ matrix(Gate(X), q)^2
    @test matrix(g^3, q) ≈ matrix(g, q)^3
    @test matrix(Gate(X) * Gate(S), q) ≈ matrix(Gate(X * S), q)
end

@testset "Superoperators in the dense basis" begin
    # on one site the order of vec(ρ) is the ket fastest, which gives every superoperator an
    # exact form owing nothing to the way TMS assembles it
    q, f, b = Qubit(), Fermion(), Boson(3)
    for (op, s) in [(Sp, q), (Sm, q), (S, q), (X, q), (H, q), (C, f), (A, b), (A + 0.5im * N, b)]
        m = matrix(op, s)
        i = eye(size(m, 1))
        mm = m' * m
        @test matrix(Gate(op), s) ≈ kron(conj(m), m)
        @test matrix(Evolver(op), s) ≈ kron(i, m) + kron(conj(m), i)
        @test matrix(Dissipator(op), s) ≈ kron(conj(m), m) - (kron(i, mm) + kron(transpose(mm), i)) / 2
    end
    # a reset gives its state whatever the state it acts on, a dephasing erases the coherences
    # between the eigenspaces and keeps those within one
    σ = [0.3 0.2im ; -0.2im 0.7]
    @test matrix(Reset(σ), q) ≈ vec(σ) * transpose(vec(eye(2)))
    @test matrix(Reset("Dn"), q) ≈ matrix(Reset(Level => 1), q)
    @test matrix(Dephase(Level), q) ≈ [1 0 0 0 ; 0 0 0 0 ; 0 0 0 0 ; 0 0 0 1]
    e = Electron()
    ρ = randn(ComplexF64, 4, 4)
    out = reshape(matrix(Dephase(Ntot), e) * vec(ρ), 4, 4)
    n = real(TMS.LinearAlgebra.diag(matrix(Ntot, e)))
    @test out ≈ [ n[i] == n[j] ? ρ[i, j] : 0 for i in 1:4, j in 1:4 ]
    @test_throws "a state of 3 components cannot be one of Qubit()" matrix(Reset([1., 0., 0.]), q)
    @test_throws "a number does not denote a state" Reset(1)
    @test_throws "write Dephase(N)(1) rather than Dephase(N(1))" Dephase(N(1))
end

@testset "The matrix of a channel acts on vec(ρ) in the order of the sites" begin
    # on two sites of different dimensions, so that no order can pass for another
    q, b = Qubit(), Boson(3)
    k = randn(ComplexF64, 6, 6)
    kop = Operator{2}("K", k)
    h = k + k'
    ρ = randn(ComplexF64, 6, 6)
    v(x) = mixed_vec(x, 2, 3)
    @test matrix(Gate(kop), q, b) * v(ρ) ≈ v(k * ρ * k')
    @test matrix(Evolver(kop), q, b) * v(ρ) ≈ v(k * ρ + ρ * k')
    @test matrix(Dissipator(kop), q, b) * v(ρ) ≈ v(k * ρ * k' - (k' * k * ρ + ρ * k' * k) / 2)
    ps = TMS.eigenspaces(h, "", "")[2]
    @test matrix(Dephase(Operator{2}("H", h)), q, b) * v(ρ) ≈ v(sum(p * ρ * p for p in ps))
    # a tensor product of channels acts on each site
    x, a = matrix(X, q), matrix(A + dag(A), b)
    @test matrix(Gate(X) ⊗ Gate(A + dag(A)), q, b) * v(ρ) ≈ v(kron(x, a) * ρ * kron(x, a)')
    @test matrix(Gate(X ⊗ (A + dag(A))), q, b) ≈ matrix(Gate(X) ⊗ Gate(A + dag(A)), q, b)
end

@testset "The tensor product of channels on fermions" begin
    # (A ⊗ B)(i, j) is A(i) * B(j) for channels too: each factor takes the strings of the odd
    # factors after it, on the kets and on the bras, as the operators of a gate do
    f = Fermion()
    ρ = randn(ComplexF64, 4, 4)
    v(x) = mixed_vec(x, 2, 2)
    for k in (C ⊗ C, C ⊗ dag(C), (C + dag(C)) ⊗ C, N ⊗ C)
        m = matrix(k, f)
        @test matrix(Gate(k), f) * v(ρ) ≈ v(m * ρ * m')
        @test matrix(Gate(k), f) ≈ matrix(Gate(k.subs[1]) ⊗ Gate(k.subs[2]), f)
    end
    @test matrix(Gate(C ⊗ N ⊗ C), f) ≈ matrix(Gate(C) ⊗ Gate(N) ⊗ Gate(C), f)
    # a reset keeps the parity of nothing on its kets, but strings on both sides cancel in it
    @test matrix(Reset("Emp") ⊗ Gate(C), f) ≈ kron(matrix(Reset("Emp"), f), matrix(Gate(C), f))
    @test_throws "no definite fermionic parity" matrix(Gate(C) ⊗ Reset("Emp"), f)
end

@testset "The dissipator of a channel" begin
    # the sum of those of its Kraus operators: Φ - id when Φ keeps the trace, and otherwise the
    # anticommutator of the sum of the K†K its trace gives
    q, b = Qubit(), Boson(3)
    @test matrix(Dissipator(Reset("Up")), q) ≈ matrix(Reset("Up"), q) - eye(4)
    @test matrix(Dissipator(0.5 * Gate(Sm) + 0.3 * Gate(Sp)), q) ≈
        0.5 * matrix(Dissipator(Sm), q) + 0.3 * matrix(Dissipator(Sp), q)
    @test matrix(Dissipator(Gate(X) ⊗ Reset("Up")), q) ≈ matrix(Gate(X) ⊗ Reset("Up"), q) - eye(16)
    @test matrix(Dissipator(Gate(Sm) ⊗ Gate(A)), q, b) ≈ matrix(Dissipator(Sm ⊗ A), q, b)
end

@testset "A channel keeps the trace, a generator keeps it constant" begin
    # tr ρ is the row vec(I)ᵀ, which a channel of the library leaves as it is and a generator
    # sends to zero
    q, b = Qubit(), Boson(3)
    t1, t2 = vec(eye(2))', mixed_vec(eye(4), 2, 2)'
    for c in (reset_gate(0.3, "Dn"), depolarizing_gate(0.2), dephasing_gate(0.4), Dephase(X),
              amplitude_damping_gate(0.3), thermal_relaxation_gate(5., 3., 0.7))
        @test t1 * matrix(c, q) ≈ t1
    end
    @test t2 * matrix(reset_gate(0.3, ["Up", "Dn"]), q) ≈ t2
    @test transpose(vec(eye(3))) * matrix(amplitude_damping_gate(0.3, A), b) ≈ transpose(vec(eye(3)))
    for g in (-im * X + Dissipator(Sm), reset_dissipator(0.4, "Up"), dephasing_dissipator(0.2),
              depolarizing_dissipator(0.3), amplitude_damping_dissipator(0.5),
              thermal_relaxation_dissipator(5., 3.))
        @test norm(t1 * matrix(g, q)) < 1e-14
    end
end

@testset "The exponential of a generator is the channel of its evolution" begin
    q = Qubit()
    γ, t = 0.7, 0.4
    p = 1 - exp(-γ * t)
    ev(L) = matrix(exp(t * L), q)
    @test exp(Dissipator(X)) isa TMS.ExpOp{IsChannel}
    @test ev(Evolver(-im * (X + 0.3Z))) ≈ matrix(Gate(exp(-im * t * (X + 0.3Z))), q)
    @test ev(amplitude_damping_dissipator(γ)) ≈ matrix(amplitude_damping_gate(p), q)
    @test ev(reset_dissipator(γ, ["Up", "Dn"])) ≈ matrix(reset_gate(p, ["Up", "Dn"]), q)
    @test ev(depolarizing_dissipator(γ, 2)) ≈ matrix(depolarizing_gate(p, 2), q)
    @test ev(dephasing_dissipator(γ, Z ⊗ Z)) ≈ matrix(dephasing_gate(p, Z ⊗ Z), q)
    @test ev(thermal_relaxation_dissipator(5., 3.)) ≈ matrix(thermal_relaxation_gate(5., 3., t), q)
    b = Boson(5)
    @test matrix(exp(t * amplitude_damping_dissipator(γ, A)), b) ≈ matrix(amplitude_damping_gate(p, A), b)
end

@testset "The library of noises" begin
    q, f, b = Qubit(), Fermion(), Boson(6)
    p = 0.3
    # the qubit forms of 1.x
    @test matrix(depolarizing_gate(p), q) ≈
        matrix((1 - 3p/4) * Gate(Id) + p/4 * (Gate(X) + Gate(Y) + Gate(Z)), q)
    @test matrix(Dissipator(Z), q) ≈ matrix(dephasing_dissipator(2, Z), q)
    @test matrix(amplitude_damping_gate(p), q) ≈ matrix(Gate(Id - (1 - sqrt(1 - p)) * N) + p * Gate(Sp), q)
    @test matrix(thermal_relaxation_gate(5., 3., 0.4), q) ≈
        matrix(dephasing_gate(1 - exp(-0.4 * (1/3 - 1/10))) * amplitude_damping_gate(1 - exp(-0.4 / 5)), q)
    # sites reset together, or each to its state, a vector of numbers being one state
    @test reset_gate(p, ["Up", "Dn"]) == (1 - p) * Gate(Id ⊗ Id) + p * Reset("Up") ⊗ Reset("Dn")
    @test reset_gate(p, 2, "Up") == (1 - p) * Gate(Id ⊗ Id) + p * Reset("Up") ⊗ Reset("Up")
    @test reset_gate(p, [1., 0.]) == (1 - p) * Gate(Id) + p * Reset([1., 0.])
    # amplitude damping on other sites: the same Kraus operators on a fermion, and on a boson
    # each photon lost on its own, the number of those left binomial
    @test matrix(amplitude_damping_gate(p, C), f) ≈ matrix(Gate(Id - (1 - sqrt(1 - p)) * N) + p * Gate(C), f)
    n, d = 4, 6
    ρ = zeros(d, d)
    ρ[n + 1, n + 1] = 1
    out = reshape(matrix(amplitude_damping_gate(p, A), b) * vec(ρ), d, d)
    @test real(TMS.LinearAlgebra.diag(out))[1:n + 1] ≈ [ binomial(n, k) * (1 - p)^k * p^(n - k) for k in 0:n ]
    # with p = 1, an infinite time, the reset to the state the decay ends in, refused when there
    # is none
    @test matrix(amplitude_damping_gate(1.), q) ≈ matrix(Reset("Up"), q)
    @test repr(amplitude_damping_gate(1.)) == "Reset(DecayEnd(Sp))"
    @test amplitude_damping_gate(1., A) == amplitude_damping_gate(1., A)
    @test matrix(amplitude_damping_gate(1.), q) ≈ matrix(Gate(Id - N) + Gate(Sp), q)
    @test matrix(amplitude_damping_gate(1., A), b) ≈ matrix(Reset("0"), b)
    @test matrix(amplitude_damping_gate(1.), Spin(1)) ≈ matrix(Reset(Level => 0), Spin(1))
    @test_throws "does not take every state of Electron() to a single one" matrix(amplitude_damping_gate(1.), Electron())
    @test matrix(thermal_relaxation_gate(5., 3., 1e4), q) ≈ matrix(Reset("Up"), q)
    # parameters out of their domain
    for g in (() -> reset_gate(1.5, "Up"), () -> depolarizing_gate(-0.1), () -> dephasing_gate(2.),
              () -> amplitude_damping_gate(1.1))
        @test_throws "takes a probability between 0 and 1" g()
    end
    @test_throws "must be positive" thermal_relaxation_gate(0., 3., 1.)
    @test_throws "is above 2T1" thermal_relaxation_dissipator(1., 3.)
    @test_throws "takes a time of zero or above" thermal_relaxation_gate(5., 3., -1.)
end

@testset "Purity of channels and generators" begin
    # a channel of a single Kraus operator, a generator with no dissipator, act on pure states;
    # a sum of channels is taken as a mixture
    pure = TMS.is_pure
    @test pure(Gate(X)(1) * Gate(Y)(2))
    @test pure(Gate(X(1) + Y(2)))
    @test pure(Gate(Id))
    @test pure(0.5 * Gate(X) ⊗ Gate(Y))
    @test pure(Gate(X)^3)
    @test pure(-im * X(1) + -im * Z(2) + Evolver(Y)(1))
    @test pure(exp(0.1 * Evolver(-im * X)))
    @test !pure(0.9 * Gate(Id) + 0.1 * Gate(X))
    @test !pure(Reset("Up")(1))
    @test !pure(Dephase(Z)(1) * Gate(X)(1))
    @test !pure(-im * X(1) + Dissipator(Sm)(1))
    @test !pure(exp(0.1 * Dissipator(X)))
end

@testset "Dense matrices of functions are bounded" begin
    # a function of an operator or a superoperator has to be built dense: beyond dimension 4096
    # it is refused rather than exhausting memory
    q = Qubit()
    x13 = reduce(⊗, fill(X, 13))
    @test_throws "needs a dense matrix of dimension 8192" matrix(exp(-im * x13), q)
    @test_throws "needs a dense matrix of dimension 8192" matrix(x13^2, q)
    @test_throws "needs a dense matrix of dimension 8192" matrix(Proj(reduce(⊗, fill(Z, 13)) => 1), q)
    x7 = reduce(⊗, fill(X, 7))
    @test_throws "needs a dense matrix of dimension 16384" matrix(exp(Dissipator(x7)), q)
    @test_throws "needs a dense matrix of dimension 16384" matrix(Dephase(x7), q)
end
