# Channels and generators: their normal forms, their algebra and what it refuses, their matrices
# on the density matrices of sites, their purity, the library of noises, and the pure ones on
# pure states.
#
# The invariant to prefer is that the matrix of a channel or a generator is the map it stands
# for, written on vec(ρ) in an order spelled out here element by element, so that the tests do
# not depend on how TMS lays it out; on a pure state, that a pure channel acts as its Kraus
# operator; on a system, that the terms of a channel or a generator give the superoperator
# built from the matrices of the operators it is made of on the whole system, Jordan-Wigner
# strings included, and that its MPO applied to a density matrix, read back through the
# combiners of the mixed indices, gives that map written on dense matrices, whatever order the
# charges give those indices, and so does applying them to it.

"""
the identity matrix of dimension `d`
"""
eye(d) = Matrix{Float64}(LA.I, d, d)

"""
vec(ρ) of a density matrix of two sites of dimensions `d1` and `d2` in the order of the matrices
of channels: the last site fastest and, on each site, the ket faster than the bra
"""
mixed_vec(ρ, d1, d2) =
    vec([ ρ[(k1 - 1) * d2 + k2, (b1 - 1) * d2 + b2]
          for k2 in 1:d2, b2 in 1:d2, k1 in 1:d1, b1 in 1:d1 ])

"""
the dense matrix of the MPO of `a` on the mixed state `ρ`, `coefs` the values of the time
functions of a vector, its rows and columns in the order of the elements of the density
matrices, see `TMS.mixed_elements`
"""
function laid_mixed(ρ::State{Mixed}, a, coefs = [1.])
    sys = ρ.system
    idx = siteinds(Mixed, sys)
    p = foldl(kron, [ eye(IT.dim(idx[k]))[:, TMS.mixed_elements(sys, k)] for k in 1:length(sys) ])
    return p * mpo_matrix(PreMPO(ρ, a), coefs) * p'
end

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
        @test matrix(Dissipator(op), s) ≈
              kron(conj(m), m) - (kron(i, mm) + kron(transpose(mm), i)) / 2
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
    n = real(LA.diag(matrix(Ntot, e)))
    @test out ≈ [ n[i] == n[j] ? ρ[i, j] : 0 for i in 1:4, j in 1:4 ]
    # the eigenvalues are told apart relative to the largest, however small the operator
    @test matrix(Dephase(1e-14 * Z), q) ≈ matrix(Dephase(Z), q)
    @test matrix(Proj(1e-14 * Z => -1e-14), q) ≈ matrix(Proj(Z => -1), q)
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
    # and a channel of no definite parity after a fermionic site is the sum of its parts that
    # have one: a reset, the sum of the gates of its Kraus operators |Emp⟩⟨i|, one even, one odd
    ks = [ Operator{1}("K$i", [ j == 1 && k == i ? 1. : 0. for j in 1:2, k in 1:2 ]) for i in 1:2 ]
    @test matrix(Gate(C) ⊗ Reset("Emp"), f) ≈ sum(matrix(Gate(C ⊗ k), f) for k in ks)
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
    t3 = transpose(vec(eye(3)))
    @test t3 * matrix(amplitude_damping_gate(0.3, A), b) ≈ t3
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
    @test ev(dephasing_dissipator(γ)) ≈ matrix(dephasing_gate(p), q)
    @test ev(reset_dissipator(γ, ["Up", "Dn"])) ≈ matrix(reset_gate(p, ["Up", "Dn"]), q)
    @test ev(depolarizing_dissipator(γ, 2)) ≈ matrix(depolarizing_gate(p, 2), q)
    @test ev(dephasing_dissipator(γ, Z ⊗ Z)) ≈ matrix(dephasing_gate(p, Z ⊗ Z), q)
    @test ev(thermal_relaxation_dissipator(5., 3.)) ≈ matrix(thermal_relaxation_gate(5., 3., t), q)
    b = Boson(5)
    @test matrix(exp(t * amplitude_damping_dissipator(γ, A)), b) ≈
          matrix(amplitude_damping_gate(p, A), b)
end

@testset "The library of noises" begin
    q, f, b = Qubit(), Fermion(), Boson(6)
    p = 0.3
    # the forms of a qubit
    @test matrix(depolarizing_gate(p), q) ≈
        matrix((1 - 3p/4) * Gate(Id) + p/4 * (Gate(X) + Gate(Y) + Gate(Z)), q)
    @test matrix(Dissipator(Z), q) ≈ matrix(dephasing_dissipator(2, Z), q)
    @test matrix(amplitude_damping_gate(p), q) ≈
        matrix(Gate(Id - (1 - sqrt(1 - p)) * N) + p * Gate(Sp), q)
    @test matrix(thermal_relaxation_gate(5., 3., 0.4), q) ≈
        matrix(dephasing_gate(1 - exp(-0.4 * (1/3 - 1/10))) *
               amplitude_damping_gate(1 - exp(-0.4 / 5)), q)
    # sites reset together, or each to its state, a vector of numbers being one state
    @test reset_gate(p, ["Up", "Dn"]) == (1 - p) * Gate(Id ⊗ Id) + p * Reset("Up") ⊗ Reset("Dn")
    @test reset_gate(p, 2, "Up") == (1 - p) * Gate(Id ⊗ Id) + p * Reset("Up") ⊗ Reset("Up")
    @test reset_gate(p, [1., 0.]) == (1 - p) * Gate(Id) + p * Reset([1., 0.])
    # amplitude damping on other sites: the same Kraus operators on a fermion, and on a boson
    # each photon lost on its own, the number of those left binomial
    @test matrix(amplitude_damping_gate(p, C), f) ≈
        matrix(Gate(Id - (1 - sqrt(1 - p)) * N) + p * Gate(C), f)
    n, d = 4, 6
    ρ = zeros(d, d)
    ρ[n + 1, n + 1] = 1
    out = reshape(matrix(amplitude_damping_gate(p, A), b) * vec(ρ), d, d)
    @test real(LA.diag(out))[1:n + 1] ≈ [ binomial(n, k) * (1 - p)^k * p^(n - k) for k in 0:n ]
    # with p = 1, an infinite time, the reset to the state the decay ends in, refused when there
    # is none
    @test matrix(amplitude_damping_gate(1.), q) ≈ matrix(Reset("Up"), q)
    @test repr(amplitude_damping_gate(1.)) == "Reset(DecayEnd(Sp))"
    @test amplitude_damping_gate(1., A) == amplitude_damping_gate(1., A)
    @test matrix(amplitude_damping_gate(1.), q) ≈ matrix(Gate(Id - N) + Gate(Sp), q)
    @test matrix(amplitude_damping_gate(1., A), b) ≈ matrix(Reset("0"), b)
    @test matrix(amplitude_damping_gate(1.), Spin(1)) ≈ matrix(Reset(Level => 0), Spin(1))
    @test_throws "does not take every state of Electron() to a single one" matrix(
        amplitude_damping_gate(1.), Electron())
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
    z13 = reduce(⊗, fill(Z, 13))
    @test_throws "needs a dense matrix of dimension 8192" matrix(Proj(z13 => 1), q)
    x7 = reduce(⊗, fill(X, 7))
    @test_throws "needs a dense matrix of dimension 16384" matrix(exp(Dissipator(x7)), q)
    @test_throws "needs a dense matrix of dimension 16384" matrix(Dephase(x7), q)
    # the dimension 2^64, which an Int wrapped around to 0, letting the matrix be built
    x64, x32 = reduce(⊗, fill(X, 64)), reduce(⊗, fill(X, 32))
    @test_throws "needs a dense matrix of dimension 18446744073709551616" matrix(exp(-im * x64), q)
    @test_throws "needs a dense matrix of dimension 18446744073709551616" matrix(Dephase(x32), q)
end

@testset "Pure channels on pure states" begin
    # a channel of a single Kraus operator K takes |ψ⟩⟨ψ| to K|ψ⟩⟨ψ|K†, the pure state K|ψ⟩,
    # each form of a channel giving its K
    q = Qubit()
    ψ = random_state(Pure, System(4, q), 4)
    for (Φ, K) in [(Gate(X)(2), X(2)), (Gate(X(1) + Z(3)), X(1) + Z(3)),
                   (0.25 * Gate(S)(1) * Gate(H)(4), 0.5 * S(1) * H(4)), (Gate(Id)(2), Id(1)),
                   ((Gate(X) ⊗ Gate(Y))(1, 3), (X ⊗ Y)(1, 3)), ((Gate(H)^3)(2), (H^3)(2)),
                   (exp(0.3 * Evolver(-im * (X ⊗ X)))(2, 3), exp(0.3 * (-im * (X ⊗ X)))(2, 3)),
                   (exp(0.2 * (-im * X + Evolver(Id)))(1), exp(0.2 * (-im * X + Id))(1))]
        @test TMS.kraus(Φ) == K
        @test norm(apply(Φ, ψ) - apply(K, ψ)) < 1e-12
    end
    # and so against the dense vector
    x, h, i = matrix(X, q), matrix(H, q), matrix(Id, q)
    v = raw_vector(ψ)
    @test raw_vector(apply(0.25 * Gate(X)(1) * Gate(H)(3), ψ)) ≈ 0.5 * kron(x, i, h, i) * v
    @test raw_vector(apply(exp(0.3 * Evolver(-im * (X ⊗ X)))(2, 3), ψ)) ≈
        kron(i, exp(-0.3im * kron(x, x)), i) * v
    # on fermions, the strings of the Kraus operator
    φ = random_state(Pure, System(4, Fermion()), 4)
    for (Φ, K) in [(Gate(C)(1) * Gate(dag(C))(3), C(1) * dag(C)(3)),
                   (Gate(C(2) + dag(C)(4)), C(2) + dag(C)(4)),
                   ((Gate(C) ⊗ Gate(C))(1, 3), (C ⊗ C)(1, 3))]
        @test norm(apply(Φ, φ) - apply(K, φ)) < 1e-12
    end
    # a mixture, a reset, a generator have no single Kraus operator, and a channel no matrix
    # element between two pure states
    @test_throws "has more than one Kraus operator" apply(Reset("Up")(1), ψ)
    @test_throws "has more than one Kraus operator" apply(0.9 * Gate(Id)(1) + 0.1 * Gate(X)(1), ψ)
    @test_throws "and Dissipator(X)(1) is a generator" apply(Dissipator(X)(1), ψ)
    @test_throws "and Evolver(X)(1) is a generator" apply(Evolver(X)(1), ψ)
    @test_throws MethodError inner(ψ, Gate(X)(1), ψ)
end

@testset "The MPO of a pure generator on a pure state" begin
    # a pure state follows the operator K of a generator with no dissipator, dψ/dt = Kψ
    ψ = random_state(Pure, System(3, Qubit()), 4)
    L = -im * (X(1) * X(2) + Z(3)) + Evolver(0.1 * Z)(2)
    K = -im * (X(1) * X(2) + Z(3)) + 0.1 * Z(2)
    @test TMS.pure_part(L) == K
    @test norm(make_mpo(ψ, L) - make_mpo(ψ, K)) < 1e-12
    @test norm(make_mpo(PreMPO(ψ, [L, Evolver(X)(1)]), [1., 0.5]) -
               make_mpo(PreMPO(ψ, [K, X(1)]), [1., 0.5])) < 1e-12
    # the identity of generators, ρ ↦ ρ, is the evolver of half the identity
    idg = TMS.IdentityOp{IsGenerator, TMS.Generic, 1}()
    @test TMS.pure_part(idg) == 0.5 * Id
    @test matrix(Evolver(TMS.pure_part(idg)), Qubit()) ≈ matrix(idg, Qubit())
    @test_throws "holds a dissipator, which a pure state cannot follow" make_mpo(ψ,
        -im * X(1) + Dissipator(Sm)(1))
    # a channel is refused as such, and an operator beside a generator is taken, where both
    # failed on a MethodError
    @test_throws "is a channel, which acts on a density matrix" make_mpo(ψ, Gate(X(1)))
    @test_throws "is a channel, which acts on a density matrix" evolve(Gate(X(1)), 0.1, ψ)
    @test norm(make_mpo(PreMPO(ψ, [-im * X(1), Evolver(-im * Z(2))]), [1., 0.5]) -
               make_mpo(PreMPO(ψ, [-im * X(1), -im * Z(2)]), [1., 0.5])) < 1e-12
end

@testset "Gates, evolvers and dissipators on a density matrix" begin
    # each is made of the operators it holds on each side of the density matrix, their strings
    # included, on sites apart
    sided = TMS.sided
    fs, qs = System(3, Fermion()), System(3, Qubit())
    for (sys, k) in [(fs, C(1) + 0.5 * dag(C)(3)), (fs, dag(C)(1) * C(3) + N(2)),
                     (qs, X(1) * Y(3) + 0.3 * Z(2)), (fs, C(2)), (fs, (C ⊗ dag(C))(3, 1))]
        K, ss = lowered(sys, k), sys.sites
        i, m = one(K), K' * K
        @test superop(sys, Gate(k)) ≈ sided(K, K', ss)
        @test superop(sys, Evolver(k)) ≈ sided(K, i, ss) + sided(i, K', ss)
        @test superop(sys, Dissipator(k)) ≈
              sided(K, K', ss) - (sided(m, i, ss) + sided(i, m, ss)) / 2
    end
    # a Hamiltonian added to a dissipator
    h = dag(C)(1) * C(3) + dag(C)(3) * C(1)
    H, c1, ss = lowered(fs, h), lowered(fs, C(1)), fs.sites
    i, m = one(H), c1' * c1
    @test superop(fs, -im * h + Dissipator(C)(1)) ≈
        sided(-im * H, i, ss) + sided(i, im * H, ss) + sided(c1, c1', ss) -
        (sided(m, i, ss) + sided(i, m, ss)) / 2
end

@testset "Resets and dephasings on a density matrix" begin
    sided = TMS.sided
    fs, qs = System(3, Fermion()), System(3, Qubit())
    # on sites with no string, the reset of a site alone
    @test superop(qs, Reset("Up")(2)) ≈ kron(eye(4), matrix(Reset("Up"), Qubit()), eye(4))
    # a fermion between fermions reset to a state of eigenvalues λ and eigenvectors v: the sum of
    # the gates of its Kraus operators |v⟩⟨i|, each with the strings of its parity
    for (st, σ) in [("Occ", [0. 0. ; 0. 1.]), ("FullyMixed", [0.5 0. ; 0. 0.5])]
        e = LA.eigen(LA.Hermitian(σ))
        ks = [ λ => lowered(fs, Operator{1}("K", v * [ j == i ? 1. : 0. for j in 1:2 ]')(2))
               for (λ, v) in zip(e.values, eachcol(e.vectors)) for i in 1:2 ]
        @test superop(fs, Reset(st)(2)) ≈ sum(λ * sided(k, k', fs.sites) for (λ, k) in ks)
    end
    # a dephasing on sites apart, the sum of the gates of its projectors
    for (sys, a) in [(qs, Z ⊗ Z), (fs, im * (C + dag(C)) ⊗ (C + dag(C))), (fs, N ⊗ N)]
        vals = TMS.eigenspaces(matrix(a, sys[1], sys[3]), "", "")[1]
        ps = [ lowered(sys, Proj(a => λ)(1, 3)) for λ in vals ]
        @test superop(sys, Dephase(a)(1, 3)) ≈ sum(sided(p, p, sys.sites) for p in ps)
    end
    # on an Electron, whose parity is that of both spins, the reset of a site leaves the others
    # in their state, which their partial trace shows
    es = System(3, Electron())
    e = mix((State{Pure}(es, ["Up", "Up", "Emp"]) + State{Pure}(es, ["Emp", "Up", "Up"])) / sqrt(2))
    @test norm(partial_trace(apply(Reset("Emp")(2), e), [2]) - partial_trace(e, [2])) < 1e-12
end

@testset "The exponential of a generator on sites apart" begin
    # exp(t * L)(i, j) is the exponential of L(i, j): its superoperator, split into
    # superoperators of one site, takes the strings of its kets and of its bras
    fs = System(3, Fermion())
    t = 0.4
    for L in [-im * (dag(C) ⊗ C - C ⊗ dag(C)) + Dissipator(C ⊗ Id) + Dissipator(Id ⊗ C),
              Evolver(0.3 * (C ⊗ Id + dag(C) ⊗ N)) + Dissipator(Id ⊗ C)]
        for index in ((1, 3), (2, 1))
            @test superop(fs, exp(t * L)(index...)) ≈ exp(t * superop(fs, L(index...)))
        end
    end
    @test superop(fs, (Gate(C)^2)(2)) ≈ superop(fs, Gate(C)(2))^2
end

@testset "The dissipator of a channel on sites apart" begin
    # its channel minus half the anticommutator of the sum of the K†K, read off the trace
    sided = TMS.sided
    for (sys, Φ) in [(System(3, Qubit()), (Gate(X) ⊗ Reset("Up"))(1, 3)),
                     (System(3, Fermion()), (Gate(C) ⊗ Reset("Emp"))(1, 3))]
        S, ss = superop(sys, Φ), sys.sites
        m = TMS.kraus_sum(S, ss)
        i = one(m)
        @test superop(sys, Dissipator(Φ)) ≈ S - (sided(m, i, ss) + sided(i, m, ss)) / 2
    end
end

@testset "The trace on a system" begin
    # tr ρ is the row of vec(I) in the order of the sites, which a channel keeping the trace
    # leaves and a generator sends to zero
    fs = System(3, Fermion())
    t = transpose(foldl(kron, fill(vec(eye(2)), 3)))
    L = -im * (dag(C)(1) * C(3) + dag(C)(3) * C(1)) + Dissipator(C(1) + C(3))
    @test norm(t * superop(fs, L)) < 1e-12
    Φ = reset_gate(0.3, "Emp")(2) * Gate((C + dag(C))(1)) * depolarizing_gate(0.2, 2)(1, 3)
    @test t * superop(fs, Φ) ≈ t
end

@testset "A channel or a generator refused on a system" begin
    q2, q3 = System(2, Qubit(conserve = N)), System(3, Qubit())
    @test_throws "carries no definite charge of N" TMS.lower(q2, Dephase(X ⊗ X)(1, 2))
    @test_throws "acts on site 5, which the system does not have" TMS.lower(q3, Gate(X)(5))
    @test_throws "acts on site 4, which the system does not have" TMS.lower(
        q3, Dissipator(X(1) + X(4)))
    # an operator applied to a density matrix is checked as its gate, as it was written, not as
    # the factors its strings are laid into, which named Gate(F)(5) for C(9) in 1.x
    ρ = mix(State{Pure}(System(4, Fermion()), "Emp"))
    @test_throws "Gate(C)(9) acts on site 9, which the system does not have" apply(C(9), ρ)
end

@testset "The MPO of a channel or a generator on a density matrix" begin
    fs = System(3, Fermion())
    ρ = State{Mixed}(fs, "FullyMixed")
    h = dag(C)(1) * C(3) + dag(C)(3) * C(1)
    l = -im * (dag(C) ⊗ C - C ⊗ dag(C)) + Dissipator(C ⊗ Id) + Dissipator(Id ⊗ C)
    for a in [-im * h + Dissipator(C)(1) + 0.5 * Dissipator(dag(C)(2)), Gate(C(1) + C(3)),
              Reset("Occ")(2), Dephase(im * (C + dag(C)) ⊗ (C + dag(C)))(1, 3), exp(0.4 * l)(1, 3),
              reset_gate(0.3, "Emp")(2) * Gate(C)(1)]
        @test laid_mixed(ρ, a) ≈ superop(fs, a)
    end
    # an operator given for a density matrix is its evolver, and a vector takes time functions
    @test laid_mixed(ρ, h) ≈ superop(fs, Evolver(h))
    @test laid_mixed(ρ, [-im * h, Dissipator(C)(2)], [0.5, 2.]) ≈
        superop(fs, -0.5im * h + 2 * Dissipator(C)(2))
    # as few channels as a triangular MPO of the superoperator can have, on sites of dimension d²
    qs = System(5, Qubit())
    hq = sum(X(i) * X(i + 1) + 0.3 * Z(i) for i in 1:4)
    L = -im * hq + sum(Dissipator(sqrt(0.1) * Sm)(i) for i in 1:5) + Dissipator(X ⊗ X)(1, 3)
    ρq = State{Mixed}(qs, "FullyMixed")
    S = superop(qs, L)
    @test laid_mixed(ρq, L) ≈ S
    @test PreMPO(ρq, L).linkdims .- 1 == least_channels(S, fill(4, 5))
    # the order of the mixed index, read once for each kind of site
    sys = System([Fermion(conserve = N), Qubit(), Fermion(conserve = N), Boson(3, conserve = N),
                  Qubit()])
    @test TMS.block_orders(Mixed, sys) == [ TMS.mixed_elements(sys, i) for i in 1:5 ]
end

@testset "Superoperators made of the sides of their operators" begin
    # a gate or a dissipator of an operator of many terms, an evolver of a product that would
    # develop and a product of such gates are made of the sides Kρ and ρK† of their operators
    # multiplied: the superoperator itself, with the least channels; on qubits, on fermions with
    # their strings and on fermions conserving N weakly and strongly
    qs, fs = System(5, Qubit()), System(4, Fermion())
    cs, ss = System(4, Fermion(conserve = N)), System(4, Fermion(conserve = strong(N)))
    hq = sum(X(i) * X(i + 1) for i in 1:4) + sum(0.3 * Z(i) for i in 1:5)
    hop = sum(dag(C)(i) * C(i + 1) + dag(C)(i + 1) * C(i) for i in 1:3)
    loss = Dissipator(sum(C(i) for i in 1:4))
    for (sys, a) in [(qs, -im * hq + Dissipator(sum(Sm(i) for i in 1:5))), (qs, Dissipator(hq)),
                     (qs, Evolver(-im * hq^2)), (qs, Gate(hq)), (qs, Gate(hq^2)),
                     (fs, -im * hop + loss), (fs, Gate(hop)), (fs, Dissipator(hop + 0.4 * C(2))),
                     (cs, -im * hop + loss), (ss, -im * hop + Dissipator(sum(N(i) for i in 1:4)))]
        start = sys === qs || sys === fs ? "FullyMixed" : ["Occ", "Emp", "Occ", "Emp"]
        ρ = State{Mixed}(sys, start)
        @test any(TMS.super_by_blocks, TMS.addends(a))
        S = superop(sys, a)
        @test laid_mixed(ρ, a) ≈ S
        @test PreMPO(ρ, a).linkdims .- 1 == least_channels(S, fill(4, length(sys)))
    end
    # a gate of two terms keeps its terms, and a loss under a strong symmetry is refused as by
    # its terms
    @test !TMS.super_by_blocks(Gate(X(1) * X(2) + Z(3)))
    @test_throws "changes N by different amounts" PreMPO(
        State{Mixed}(ss, ["Occ", "Emp", "Occ", "Emp"]), loss)
    # K ρ K† of many terms applied to a density matrix
    k = hq + 0.5 * X(1) * Z(3)
    ρ = random_state(Mixed, qs, 2; rng = Xoshiro(20261016))
    @test ket_bra_matrix(apply(k, ρ)) ≈ lowered(qs, k) * ket_bra_matrix(ρ) * lowered(qs, k)'
end

@testset "Strings carried by the channels on a density matrix" begin
    # the strings of the kets and of the bras, laid on the blocks of a term only, carried by the
    # channels between them, the identity on the qubit
    sys = System([Fermion(), Qubit(), Fermion(), Fermion()])
    ρ = State{Mixed}(sys, "FullyMixed")
    h = dag(C)(1) * C(4) + dag(C)(4) * C(1)
    for a in [-im * h + sum(Dissipator(C)(j) for j in (1, 3, 4)) + Dissipator(dag(C)(1) * C(4)) +
                  Dissipator(X)(2),
              Reset("Occ")(4) * Gate(C(1) + C(3))]
        S = superop(sys, a)
        @test laid_mixed(ρ, a) ≈ S
        @test PreMPO(ρ, a).linkdims .- 1 == least_channels(S, fill(4, 4))
    end
end

@testset "The MPO applied to a density matrix, its sites charged or not" begin
    # the density matrix read back is the map applied to the one it was made from, written on
    # dense matrices
    h = dag(C)(1) * C(3) + dag(C)(3) * C(1)
    for site in (Fermion(), Fermion(conserve = N))
        sys = System(3, site)
        ρ = random_state(Mixed, sys, ["Occ", "Emp", "Occ"], 3)
        r = ket_bra_matrix(ρ)
        H, k = lowered(sys, h), lowered(sys, C(1) + C(3))
        d(c) = c * r * c' - (c' * c * r + r * c' * c) / 2
        after(a) = ket_bra_matrix(apply(make_mpo(ρ, a), ρ))
        @test after(Gate(C(1) + C(3))) ≈ k * r * k'
        @test after(-im * h + Dissipator(C)(1) + Dissipator(dag(C)(2))) ≈
            -im * (H * r - r * H) + d(lowered(sys, C(1))) + d(lowered(sys, dag(C)(2)))
        @test after(h) ≈ H * r + r * H
        ks = [ lowered(sys, Operator{1}("K", [1., 0.] * [ j == i ? 1. : 0. for j in 1:2 ]')(2))
               for i in 1:2 ]
        @test after(Reset("Emp")(2)) ≈ sum(k * r * k' for k in ks)
    end
end

@testset "A channel or a generator under a strong symmetry" begin
    ss = System(2, Fermion(conserve = strong(N)))
    ρ = mix(random_state(State{Pure}(ss, ["Occ", "Emp"]), 2))
    r = ket_bra_matrix(ρ)
    n1 = lowered(ss, N(1))
    # a jump commuting with the charge is what a strong symmetry asks for
    @test ket_bra_matrix(apply(make_mpo(ρ, Dissipator(N)(1)), ρ)) ≈
        n1 * r * n1 - (n1 * n1 * r + r * n1 * n1) / 2
    # one changing it in some parts and not in others would mix the sectors it keeps apart
    @test_throws("changes N by different amounts from one part to another, which conserving it " *
                 "strongly forbids", make_mpo(ρ, Reset("Emp")(1)))
    @test_throws "drop `strong` or weaken the state" make_mpo(ρ,
        -im * (dag(C)(1) * C(2) + dag(C)(2) * C(1)) + Dissipator(C)(1))
    # a jump of no definite charge is refused for what it lacks, weak or strong, its parts of
    # each charge changing it by amounts that differ
    for site in (Boson(3, conserve = N), Boson(3, conserve = strong(N)))
        ρb = mix(State{Pure}(System(2, site), "1"))
        @test_throws "so it has no definite flux" make_mpo(ρb, Gate(Q)(1))
        @test_throws "so it has no definite flux" make_mpo(ρb, Dissipator(Q)(1))
    end
end

@testset "Channels and generators applied to a density matrix" begin
    # the density matrix read back is the map applied, written on dense matrices: a gate as its
    # operator on both sides, any other factor as a dense block or through its MPO, on sites
    # charged or not, the states entangled so that no result vanishes
    h = dag(C) ⊗ C + dag(dag(C) ⊗ C)
    hp = dag(C)(1) * C(3) + dag(C)(3) * C(1)
    for site in (Fermion(), Fermion(conserve = N))
        sys = System(4, site)
        ρ = mix(random_state(State{Pure}(sys, ["Occ", "Emp", "Occ", "Emp"]), 3)) +
            mix(random_state(State{Pure}(sys, ["Emp", "Occ", "Occ", "Emp"]), 3))
        r = ket_bra_matrix(ρ)
        after(a) = ket_bra_matrix(apply(a, ρ))
        k = lowered(sys, C(1) * dag(C)(3) + 0.5 * N(2))
        @test after(Gate(C(1) * dag(C)(3) + 0.5 * N(2))) ≈ k * r * k'
        k = lowered(sys, C(2) * dag(C)(4))
        @test after(Gate(C)(2) * Gate(dag(C))(4)) ≈ k * r * k'
        # an operator applied to a density matrix is its gate
        u = lowered(sys, exp(-0.4im * h)(1, 3))
        @test after(exp(-0.4im * h)(1, 3)) ≈ u * r * u'
        ks = [ lowered(sys, Operator{1}("K", [1., 0.] * [ j == i ? 1. : 0. for j in 1:2 ]')(2))
               for i in 1:2 ]
        reset(x) = sum(k * x * k' for k in ks)
        @test after(Reset("Emp")(2)) ≈ reset(r)
        k = lowered(sys, C(1) + C(3))
        @test after(reset_gate(0.3, "Emp")(2) * Gate(C(1) + C(3))) ≈
            0.7 * k * r * k' + 0.3 * reset(k * r * k')
        ps = [ lowered(sys, Proj(N ⊗ N => λ)(1, 3)) for λ in (0., 1.) ]
        @test after(Dephase(N ⊗ N)(1, 3)) ≈ sum(p * r * p for p in ps)
        H, c1 = lowered(sys, hp), lowered(sys, C(1))
        @test after(-im * hp + Dissipator(C)(1)) ≈
            -im * (H * r - r * H) + c1 * r * c1' - (c1' * c1 * r + r * c1' * c1) / 2
        @test after(0.5 * Gate(Id)(1)) ≈ 0.5 * r
        # through its MPO on the sites it acts on alone, the first one left out
        kh = lowered(sys, dag(C)(2) * C(4) + dag(C)(4) * C(2))
        @test after(Gate(dag(C)(2) * C(4) + dag(C)(4) * C(2))) ≈ kh * r * kh'
        # inner(a, Φ, b) is tr(a† Φ(b))
        ρ2 = mix(random_state(State{Pure}(sys, ["Emp", "Occ", "Occ", "Emp"]), 2))
        @test inner(ρ2, Gate(C(1) + C(3)), ρ) ≈ LA.tr(ket_bra_matrix(ρ2)' * k * r * k')
    end
end

@testset "Noisy gates" begin
    # a mixture of gates leaves one qubit fully mixed, and two qubits mixed with the correlation
    # ZZ they keep, the purity tr ρ² being inner(ρ, ρ)
    ρ = apply((0.5 * Gate(Id) + 0.5 * Gate(X))(1), State{Mixed}(System(1, Qubit()), "Up"))
    @test norm([expect(ρ, X(1)), expect(ρ, Y(1)), expect(ρ, Z(1))]) < 1e-12
    @test inner(ρ, ρ) ≈ 0.5
    ρ = apply((0.5 * Gate(Id ⊗ Id) + 0.5 * Gate(X ⊗ X))(1, 2),
              State{Mixed}(System(2, Qubit()), "Up"))
    @test norm([expect1(ρ, X); expect1(ρ, Y); expect1(ρ, Z)]) < 1e-12
    @test inner(ρ, ρ) ≈ 0.5
    @test expect(ρ, Z(1) * Z(2)) ≈ 1
    # a power of a noisy gate goes into an MPO as well, and the two must agree
    ρ = mix(State{Pure}(System(2, Qubit()), ["Up", "+"]))
    for a in [(Gate(X)^2)(1), ((0.9 * Gate(Id) + 0.1 * Gate(X))^3)(1)]
        @test norm(apply(a, ρ) - apply(make_mpo(ρ, a), ρ)) < 1e-12
    end
    # a power that is not an integer is no channel, which 1.x applied
    @test_throws "Gate(H) has no power 0.5" Gate(H)^0.5
end

@testset "A product of gates is the operator it denotes" begin
    # its rightmost factor acts first, as in make_mpo and expect, and on a density matrix a gate
    # g acts as g ρ g†, what mixing the pure result gives, on qubits and beside a fermion
    q = State{Pure}(System(2, Qubit()), ["Up", "X+"])
    qf = State{Pure}(System([Qubit(), Fermion()]), ["X+", "Emp"])
    for (st, g) in [(q, Swap(1, 2) * X(1)), (q, X(1) * Swap(1, 2)), (q, Sp(2) * Sm(2)),
                    (qf, Sp(1) * Sm(1) * dag(C)(2)), (qf, dag(C)(2) * Sm(1) * Sp(1))]
        @test norm(apply(g, st) - apply(make_mpo(st, g), st)) < 1e-12
        @test norm(apply(g, mix(st)) - mix(apply(g, st))) < 1e-12
    end
end

@testset "Fermionic gates" begin
    # the two branches differ in the parity of sites 1 and 2, which makes the string acting
    # further right observable instead of a global phase
    fe = Fermion()
    sys = System(4, fe)
    st = normalize(State{Pure}(sys, ["1", "0", "1", "0"]) + State{Pure}(sys, ["0", "0", "1", "0"]))
    ρ = mix(st)
    # mixing the result of the pure application gives what applying it to ρ gives
    for a in [C(3), dag(C)(4) * C(3)]
        @test norm(mix(apply(a, st)) - apply(a, ρ)) < 1e-12
    end
    # a fermionic operator inside a gate takes its string, which 1.x built from the bare matrices
    for a in [Gate(C)(3), Gate(C ⊗ dag(C))(3, 4), Gate(C + dag(C))(3)]
        @test norm(apply(a, ρ) - apply(make_mpo(ρ, a), ρ)) < 1e-12
    end
    # the descent lays the strings, so that a dissipator of a fermionic operator and a function
    # of an odd operator, of several sites or on site 1, are applied, here against the dense
    # matrices of c₁ = c ⊗ 1 ⊗ 1 ⊗ 1 and c₃ = F ⊗ F ⊗ c ⊗ 1
    c, f, id = matrix(C, fe), matrix(F, fe), matrix(Id, fe)
    c1, c3 = kron(c, id, id, id), kron(f, f, c, id)
    r, v = ket_bra_matrix(ρ), raw_vector(st)
    @test ket_bra_matrix(apply(Dissipator(C)(3), ρ)) ≈
        c3 * r * c3' - (c3' * c3 * r + r * c3' * c3) / 2
    @test raw_vector(apply(exp(0.3 * (C ⊗ Id + Id ⊗ C))(1, 3), st)) ≈ exp(0.3 * (c1 + c3)) * v
    @test raw_vector(apply(exp(0.3 * (C + dag(C)))(1), st)) ≈ exp(0.3 * (c1 + c1')) * v
    # a gate of several sites defined by a matrix is kept whole, and a string covering only some
    # of its sites does not commute with it: moved past it in 1.x, it changed the sign of the
    # branch where the gate moves a fermion onto site 2
    sw = Operator{2}("Sw", [1. 0. 0. 0. ; 0. 0. 1. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 1.])
    a = dag(C)(1) * dag(C)(4) * sw(1, 2)
    @test norm(apply(a, st) - apply(dag(C)(1), apply(dag(C)(4), apply(sw(1, 2), st)))) < 1e-12
    @test norm(mix(apply(a, st)) - apply(a, ρ)) < 1e-12
    # a term of coefficient zero leaves the parity of a gate alone, where it made one a sum of
    # fermionic and non fermionic operators, and a gate whose terms all vanish makes the state
    # null, as a gate that annihilates it does
    g = 0.0
    @test norm(apply((g * C + dag(C))(4), st) - apply(dag(C)(4), st)) < 1e-12
    @test norm(apply((g * C + g * dag(C))(4), st)) == 0
    @test trace(apply(0 * X(1), mix(State{Pure}(System(2, Qubit()), "Up")))) == 0
    # a factor contributing no tensor, and Swap, defined by an expression, which developed would
    # be a sum, alone or beside a fermionic factor
    @test_ok apply(Id(1) * X(2), State{Pure}(System(2, Qubit()), "Up"))
    @test_ok apply(Swap(1, 2) * Swap(3, 4), State{Pure}(System(4, Qubit()), "Up"))
    sq = State{Pure}(System([Qubit(), Qubit(), Fermion()]), ["Up", "Dn", "0"])
    @test norm(apply(Swap(1, 2) * dag(C)(3), sq) - apply(Swap(1, 2), apply(dag(C)(3), sq))) < 1e-12
end

@testset "Tensor products of superoperators" begin
    # (A ⊗ B)(i, j) is A(i) * B(j) for channels too, an operator being its gate, so that a
    # channel of several sites is written once and placed anywhere: two qubits depolarized
    # together against the Pauli form, on sites apart and reversed, and through the MPO
    p = 0.4
    FM = Reset("FullyMixed")
    paulis = (Id, X, Y, Z)
    depolarizing2 = (1 - p) * Gate(Id ⊗ Id) + p * (FM ⊗ FM)
    pauli_form = (1 - p) * Gate(Id ⊗ Id) + p / 16 * sum(Gate(a ⊗ b) for a in paulis for b in paulis)
    ρ = mix(random_state(Pure, System(3, Qubit()), 4))
    for s in [(1, 2), (1, 3), (3, 1)]
        @test norm(apply(depolarizing2(s...), ρ) - apply(pauli_form(s...), ρ)) < 1e-12
    end
    @test norm(apply(make_mpo(ρ, depolarizing2(1, 3)), ρ) - apply(depolarizing2(1, 3), ρ)) < 1e-12
    @test norm(apply((X ⊗ Reset("Up"))(1, 3), ρ) - apply(Gate(X)(1) * Reset("Up")(3), ρ)) < 1e-12
    # on a charged system, where the difference of two equal states is exactly zero
    w = mix(State{Pure}(System(3, Qubit(conserve = N)), ["Up", "Dn", "Up"]))
    @test norm(apply((Gate(Sp) ⊗ Gate(Sm))(2, 3), w) - apply(Gate(Sp)(2) * Gate(Sm)(3), w)) < 1e-12
    # on fermions, the strings are laid factor by factor, in either order of the sites; C ⊗
    # Reset("Emp"), whose Reset has no definite parity after a fermionic site, is refused, see
    # "The tensor product of channels on fermions"
    f = mix(random_state(Pure, System(4, Fermion()), 4))
    for (g, s) in [(Gate(C) ⊗ Gate(dag(C)), (1, 3)), (Gate(C) ⊗ Gate(dag(C)), (3, 1)),
                   (FM ⊗ FM, (1, 3))]
        ref = apply(prod(o(i) for (o, i) in zip((g.subs...,), s)), f)
        @test norm(apply(g(s...), f) - ref) < 1e-12
    end
    @test norm(apply(make_mpo(f, depolarizing2(2, 4)), f) - apply(depolarizing2(2, 4), f)) < 1e-12
end

@testset "Depolarization" begin
    # the channel on any site, compared with the Pauli form on qubits, and the generator, whose
    # evolution for a time t is the channel of probability 1 - exp(-γt)
    p = 0.4
    paulis = (Id, X, Y, Z)
    q = mix(random_state(Pure, System(3, Qubit()), 4))
    one_site = (1 - 3p / 4) * Gate(Id) + p / 4 * (Gate(X) + Gate(Y) + Gate(Z))
    @test norm(apply(depolarizing_gate(p)(2), q) - apply(one_site(2), q)) < 1e-12
    two_sites = (1 - p) * Gate(Id ⊗ Id) + p / 16 * sum(Gate(a ⊗ b) for a in paulis for b in paulis)
    @test norm(apply(depolarizing_gate(p, 2)(1, 3), q) - apply(two_sites(1, 3), q)) < 1e-12
    s = State{Mixed}(System(2, Spin(1)), "1")
    @test expect(apply(depolarizing_gate(p)(1), s), Proj("1")(1)) ≈ 1 - p + p / 3
    γ, t = 0.7, 0.5
    for (ρ, n, sites) in [(mix(random_state(Pure, System(3, Spin(1)), 4)), 1, (2,)), (q, 2, (1, 3))]
        evolved = evolve(depolarizing_dissipator(γ, n)(sites...), t, ρ; nsteps = 20,
                         limits = Limits(maxdim = 100))
        @test norm(evolved - apply(depolarizing_gate(1 - exp(-γ * t), n)(sites...), ρ)) < 1e-10
    end
    strong_sys = System(2, Qubit(conserve = strong(N)))
    @test_throws "conserving it strongly forbids" apply(depolarizing_gate(p)(1),
                                                        State{Mixed}(strong_sys, "Up"))
end

@testset "Reset" begin
    # a reset overwrites its site whatever it held, mixed or not, to a named state or to a
    # density matrix, leaves the other sites as they were and keeps the trace
    st = State{Mixed}(System(3, Qubit()), ["Dn", "FullyMixed", "+"])
    @test expect1(apply(Reset("Up")(1), st), Z)[1] ≈ 1
    st3 = apply(Reset("Dn")(2), st)
    @test expect1(st3, Z)[2] ≈ -1
    @test expect1(st3, Z)[[1, 3]] ≈ expect1(st, Z)[[1, 3]]
    st4 = apply(Reset("FullyMixed")(1), st)
    @test expect1(st4, Z)[1] ≈ 0 atol = 1e-12
    @test trace(st4) ≈ 1
    p0 = 0.3
    st5 = apply(Reset([p0 0. ; 0. 1 - p0])(1), st)
    @test expect1(st5, Z)[1] ≈ 2p0 - 1
    @test trace(st5) ≈ 1
end

@testset "Reset channels" begin
    # reset with probability p, to a state, n times the same state or a state per site, and its
    # generator, whose evolution for a time t is the channel of probability 1 - exp(-γt)
    p = 0.4
    q = mix(random_state(Pure, System(3, Qubit()), 4))
    applied(g) = apply(g, q)
    @test norm(applied(reset_gate(p, "Up")(2)) -
               applied(((1 - p) * Gate(Id) + p * Reset("Up"))(2))) < 1e-12
    each = (1 - p) * Gate(Id ⊗ Id) + p * (Reset("Up") ⊗ Reset("Dn"))
    @test norm(applied(reset_gate(p, ["Up", "Dn"])(1, 3)) - applied(each(1, 3))) < 1e-12
    # a level is given as Level => 0, a number denoting no state
    by_level = applied(reset_gate(p, 2, Level => 0)(1, 3))
    @test norm(by_level - applied(reset_gate(p, ["Up", "Up"])(1, 3))) < 1e-12
    # a vector of numbers is the amplitudes of a single state, as for State
    by_vector = applied(reset_gate(p, [0., 1.])(2))
    @test norm(by_vector - applied(reset_gate(p, "Dn")(2))) < 1e-12
    γ, t = 0.7, 0.5
    evolved = evolve(reset_dissipator(γ, "Up")(2), t, q; nsteps = 20, limits = Limits(maxdim = 64))
    @test norm(evolved - applied(reset_gate(1 - exp(-γ * t), "Up")(2))) < 1e-10
    @test depolarizing_gate(p, 2) == reset_gate(p, 2, "FullyMixed")
end

@testset "Reset on fermions" begin
    # the reset of a fermionic site keeps the correlations of the other sites across it, which
    # its matrix, laid with no Jordan-Wigner string, changed in 1.x: compared with its Kraus
    # operators, each gate taking its strings
    kraus = Dict("Emp" => [(1, Gate(Id - N)), (1, Gate(C))],
                 "Occ" => [(1, Gate(N)), (1, Gate(dag(C)))],
                 "FullyMixed" => [(0.5, Gate(Id - N)), (0.5, Gate(C)), (0.5, Gate(N)),
                                  (0.5, Gate(dag(C)))])
    exact(s, i, ρ) = sum(c * apply(g(i), ρ) for (c, g) in kraus[s])
    q, γ, t = 0.3, 0.7, 0.5
    sys = System(4, Fermion())
    configs = [ collect(v) for v in Iterators.product(fill(["Emp", "Occ"], 4)...) ]
    function sector(k)
        s = sum(randn(ComplexF64) * State{Pure}(sys, v)
                for v in configs if count(==("Occ"), v) == k)
        return s / norm(s)
    end
    ρ = 0.7 * mix(sector(2)) + 0.3 * mix(sector(1))
    for s in ["Emp", "Occ", "FullyMixed"], i in 1:4
        @test norm(apply(Reset(s)(i), ρ) - exact(s, i, ρ)) < 1e-12
    end
    @test norm(apply(reset_gate(q, "Occ")(3), ρ) - ((1 - q) * ρ + q * exact("Occ", 3, ρ))) < 1e-12
    @test norm(apply(depolarizing_gate(q, 2)(2, 4), ρ) -
               ((1 - q) * ρ + q * exact("FullyMixed", 2, exact("FullyMixed", 4, ρ)))) < 1e-12
    @test norm(apply((Reset("Occ") ⊗ Reset("Emp"))(3, 2), ρ) -
               exact("Emp", 2, exact("Occ", 3, ρ))) < 1e-12
    @test norm(apply(Gate(C)(1) * Reset("Occ")(3), ρ) - apply(Gate(C)(1), exact("Occ", 3, ρ))) <
          1e-12
    @test norm(apply(exp(t * reset_dissipator(γ, "Emp"))(3), ρ) -
               apply(reset_gate(1 - exp(-γ * t), "Emp")(3), ρ)) < 1e-12
    # the generator, through its MPO, against the channel
    for (generator, channel) in
            [(reset_dissipator(γ, "Occ")(2), reset_gate(1 - exp(-γ * t), "Occ")(2)),
             (depolarizing_dissipator(γ, 2)(2, 4), depolarizing_gate(1 - exp(-γ * t), 2)(2, 4))]
        evolved = evolve(generator, t, ρ; nsteps = 20, limits = Limits(maxdim = 100))
        @test norm(evolved - apply(channel, ρ)) < 1e-10
    end
    # a qubit between the fermions takes no sign
    mixed_sys = System([Fermion(), Qubit(), Fermion(), Fermion()])
    m = mix((State{Pure}(mixed_sys, ["Occ", "Up", "Occ", "Emp"]) +
             State{Pure}(mixed_sys, ["Emp", "Dn", "Occ", "Occ"])) / sqrt(2))
    @test norm(apply(Reset("Emp")(3), m) - exact("Emp", 3, m)) < 1e-12
    # on an Electron, whose parity is that of both spins, the reset generator leaves the other
    # sites in their state, see "Resets and dephasings on a density matrix"
    esys = System(3, Electron())
    e = mix((State{Pure}(esys, ["Up", "Up", "Emp"]) + State{Pure}(esys, ["Emp", "Up", "Up"])) /
            sqrt(2))
    evolved = evolve(reset_dissipator(γ, "Dn")(2), t, e; nsteps = 20, limits = Limits(maxdim = 100))
    @test norm(partial_trace(evolved, [2]) - partial_trace(e, [2])) < 1e-10
    # the descent lays the strings of each parity, so that a reset to a state of no definite
    # parity, against its Kraus operators |v⟩⟨i|, and a mixture of a reset and a fermionic gate
    # are applied
    v = [1., 1.] / sqrt(2)
    ks = [ Operator{1}("K", v * [ j == i ? 1. : 0. for j in 1:2 ]') for i in 1:2 ]
    @test norm(apply(Reset(v)(2), ρ) - sum(apply(Gate(k)(2), ρ) for k in ks)) < 1e-12
    @test norm(apply((0.5 * Gate(C) + 0.5 * Reset("Emp"))(2), ρ) -
               (0.5 * apply(Gate(C)(2), ρ) + 0.5 * exact("Emp", 2, ρ))) < 1e-12
    # and a function of a reset of several sites, which was refused, its second factor having
    # no definite parity after a fermionic site: a reset being idempotent, its evolution is the
    # reset with probability 1 - exp(-γt), on sites apart as well
    @test norm(apply(exp(t * reset_dissipator(γ, 2, "Emp"))(2, 4), ρ) -
               apply(reset_gate(1 - exp(-γ * t), 2, "Emp")(2, 4), ρ)) < 1e-12
    @test norm(apply(exp(t * depolarizing_dissipator(γ, 2))(1, 3), ρ) -
               apply(depolarizing_gate(1 - exp(-γ * t), 2)(1, 3), ρ)) < 1e-12
end

@testset "Qubit decay and relaxation" begin
    # the decay goes from "Dn", the excitation N counts, to "Up", the coherences keeping
    # sqrt(1 - p); the relaxation of times T1 and T2 is that decay and a dephasing, which
    # commute, so that the generator for a time t gives the product of their channels
    p = 0.3
    ψ = State{Pure}(System(2, Qubit()), [[0.6, 0.8], "Up"])
    decayed = apply(amplitude_damping_gate(p)(1), mix(ψ))
    @test trace(decayed) ≈ 1
    @test expect(decayed, N(1)) ≈ 0.64 * (1 - p)
    @test expect(decayed, X(1)) ≈ 2 * 0.48 * sqrt(1 - p)
    γ, t, T1, T2 = 0.7, 0.5, 2., 1.5
    by_rate = evolve(amplitude_damping_dissipator(γ)(1), t, mix(ψ); nsteps = 20,
                     limits = Limits(maxdim = 16))
    @test norm(by_rate - apply(amplitude_damping_gate(1 - exp(-γ * t))(1), mix(ψ))) < 1e-10
    by_times = evolve(thermal_relaxation_dissipator(T1, T2)(1), t, mix(ψ); nsteps = 20,
                      limits = Limits(maxdim = 16))
    @test norm(by_times - apply(thermal_relaxation_gate(T1, T2, t)(1), mix(ψ))) < 1e-10
    # read on the channel
    relaxed = apply(thermal_relaxation_gate(T1, T2, t)(1), mix(ψ))
    @test expect(relaxed, N(1)) ≈ 0.64 * exp(-t / T1)
    @test expect(relaxed, X(1)) ≈ 2 * 0.48 * exp(-t / T2)
    @test_throws "no relaxation reaches" thermal_relaxation_gate(1., 3., t)
    @test_throws "must be positive" thermal_relaxation_dissipator(0., 1.)
    @test_throws "between 0 and 1" amplitude_damping_gate(1.5)
end

@testset "Dephasing" begin
    # Dephase(A) erases the coherences between the eigenspaces of A and keeps those within one;
    # its generator damps them all at one rate, where Dissipator(A) depends on the eigenvalues.
    # Level, the basis of the site, is the default of dephasing_gate and dephasing_dissipator
    b = mix(State{Pure}(System(2, Boson(4)), [[1., 1., 1., 0.] / sqrt(3), "0"]))
    @test norm(apply(Dephase(Level)(1), b) - apply(Dephase(N)(1), b)) < 1e-12
    @test abs(expect(apply(Dephase(N)(1), b), A(1))) < 1e-12
    @test expect(apply(Dephase(N)(1), b), N(1)) ≈ expect(b, N(1))
    @test norm(apply(make_mpo(b, dephasing_gate(0.3)(1)), b) - apply(dephasing_gate(0.3)(1), b)) <
          1e-12
    γ, t = 0.7, 0.5
    evolved = evolve(dephasing_dissipator(γ)(1), t, b; nsteps = 20, limits = Limits(maxdim = 100))
    @test norm(evolved - apply(dephasing_gate(1 - exp(-γ * t))(1), b)) < 1e-10
    q = mix(random_state(Pure, System(3, Qubit()), 4))
    by_rate = evolve(dephasing_dissipator(2, Z)(2), t, q; nsteps = 20, limits = Limits(maxdim = 64))
    by_dissipator = evolve(Dissipator(Z)(2), t, q; nsteps = 20, limits = Limits(maxdim = 64))
    @test norm(by_rate - by_dissipator) < 1e-10
    @test dephasing_gate(0.3) == dephasing_gate(0.3, Level)
    # on an electron, Ntot does not tell Up from Dn
    e = mix(State{Pure}(System(1, Electron()), [[0., 1., 1., 0.] / sqrt(2)]))
    coherence = named([ i == 2 && j == 3 ? 1. : 0. for i in 1:4, j in 1:4 ], "UpDn")
    @test expect(apply(Dephase(Ntot)(1), e), coherence(1)) ≈ 0.5
    @test abs(expect(apply(Dephase(Level)(1), e), coherence(1))) < 1e-12
    # of several sites, as a gate: the parity of two qubits of |++⟩ not read leaves mixed the two
    # Bell pairs of X1X2 = 1; in a generator, the dissipators of its projectors
    pp = State{Mixed}(System(2, Qubit()), "+")
    @test Dephase(Z ⊗ Z) isa GenericOp{IsChannel, 2}
    @test repr(Dephase(Z ⊗ Z)) == "Dephase(Z⊗Z)"
    dz = apply(Dephase(Z ⊗ Z)(1, 2), pp)
    @test real(expect(dz, X(1) * X(2))) ≈ 1
    @test abs(expect(dz, Z(1) * Z(2))) < 1e-12
    @test abs(expect(dz, X(1))) < 1e-12
    @test real(expect(apply(dephasing_gate(0.3, Z ⊗ Z)(1, 2), pp), X(1))) ≈ 0.7
    @test matrix(dephasing_dissipator(1.0, Z ⊗ Z), Qubit()) ≈
          matrix(sum(Dissipator(Proj(Z ⊗ Z => λ)) for λ in (-1, 1)), Qubit())
    # its matrix is split into factors of one site, so that its MPO gives the channel
    @test norm(apply(make_mpo(pp, Dephase(Z ⊗ Z)(1, 2)), pp) - dz) < 1e-12
    # diagonal in the charges, it passes a strong conservation
    strong_q = State{Mixed}(System(2, Qubit(conserve = strong(N))), "Up")
    @test trace(apply(dephasing_gate(0.3)(1), strong_q)) ≈ 1
    @test_throws "Dephase needs a Hermitian operator" apply(Dephase(Sp)(1), pp)
end

@testset "Exponentials of superoperators" begin
    # exp(t * L) is the channel of the evolution under L for a time t; its closed forms are in
    # "The exponential of a generator is the channel of its evolution"
    t = 0.5
    q = mix(random_state(Pure, System(3, Qubit()), 4))
    evolved = evolve(Dissipator(Sm)(2), t, q; nsteps = 20)
    @test norm(apply(exp(t * Dissipator(Sm))(2), q) - evolved) < 1e-10
    # of two sites, against tdvp on two sites, whose block holds the whole evolution, the
    # Hamiltonian written with ⊗
    h = X ⊗ X + 0.5 * (Z ⊗ Id + Id ⊗ Z)
    l2 = -im * h + Dissipator(Sm ⊗ Id) + Dissipator(Id ⊗ Sm)
    lp = -im * (X(1) * X(2) + 0.5 * (Z(1) + Z(2))) + Dissipator(Sm)(1) + Dissipator(Sm)(2)
    ρ = mix(random_state(Pure, System(2, Qubit()), 2))
    lim = Limits(cutoff = 0)
    evolved = evolve(lp, t, ρ; nsteps = 4, limits = lim)
    @test norm(apply(exp(t * l2)(1, 2), ρ; limits = lim) - evolved) < 1e-10
    # the identity through cosh and sinh, Evolver(Id) being twice the identity
    @test norm(apply(exp(0.25 * Evolver(Id))(1), q) - exp(0.5) * q) < 1e-12
    # a fermionic jump, its strings on both sides of ρ, has one: the amplitude damping of the
    # fermion, against its Kraus operators
    @test exp(t * Dissipator(C)) isa TMS.ExpOp{IsChannel}
    pd = 1 - exp(-t)
    fs = mix(random_state(Pure, System(3, Fermion()), 4))
    @test norm(apply(exp(t * Dissipator(C))(2), fs) -
               apply((Gate(Id - (1 - sqrt(1 - pd)) * N) + pd * Gate(C))(2), fs)) < 1e-12
    # the number of fermions has one, on a site conserving it: a dephasing at rate 1/2
    f = State{Mixed}(System(3, Fermion(conserve = N)), ["Occ", "Emp", "Occ"])
    f = apply(exp(-0.4im * (dag(C) ⊗ C + dag(dag(C) ⊗ C)))(1, 2), f)
    @test norm(apply(exp(t * Dissipator(N))(2), f) - evolve(Dissipator(N)(2), t, f; nsteps = 20)) <
          1e-10
    @test norm(apply(exp(t * Dissipator(N))(2), f) - apply(dephasing_gate(1 - exp(-t / 2))(2), f)) <
          1e-12
end

@testset "Larger channels on a density matrix" begin
    # a dephasing on 6 qubits, beyond dense_limit, through its MPO
    qs = System(6, Qubit())
    ρ = random_state(Mixed, qs, 3)
    r = ket_bra_matrix(ρ)
    z = lowered(qs, prod(Z(i) for i in 1:6))
    ps = [ (eye(64) + s * z) / 2 for s in (1, -1) ]
    @test ket_bra_matrix(apply(Dephase(reduce(⊗, fill(Z, 6)))(1, 2, 3, 4, 5, 6), ρ)) ≈
        sum(p * r * p for p in ps)
    # and a collective dissipator on 40 qubits, its result traceless, whose block of dimension
    # 4^40 an Int wrapped around to 0, taking it for dense, as 2^64 for a sum on 64 qubits
    @test !TMS.dense_enough(Mixed, System(32, Qubit()), 1:32)
    @test !TMS.dense_enough(Pure, System(64, Qubit()), 1:64)
    q40 = State{Mixed}(System(40, Qubit()), "Up")
    @test abs(trace(apply(Dissipator(sum(Sm(i) for i in 1:40)), q40))) < 1e-10
    @test norm(apply(sum(X(i) for i in 1:64), State{Pure}(System(64, Qubit()), "Up")))^2 ≈ 64
    # a hundred gates, applied one by one, never developed
    q5 = System(5, Qubit())
    ρ5 = mix(State{Pure}(q5, ["Dn", "Up", "Up", "Up", "Up"]))
    swaps = [ Swap(mod1(i, 4), mod1(i, 4) + 1) for i in 1:100 ]
    one_by_one = foldl((s, g) -> apply(g, s), reverse(swaps); init = ρ5)
    @test norm(apply(prod(swaps), ρ5) - one_by_one) < 1e-10
    # a dissipator of a sum or of a product, and a sum of gates, on more than two sites, go through
    # their MPO, 20 to 250 times faster than their dense blocks on 4 and 5 sites, as the pieces of
    # an operator; a channel given by its matrix stays dense
    r5 = random_state(Mixed, q5, 3)
    for d in (Dissipator(sum(Sm(i) for i in 1:5)), Dissipator(prod(X(i) for i in 1:5)),
              0.5 * Gate(Id(1)) + 0.5 * Gate(prod(Z(i) for i in 1:5)))
        @test all(s -> s isa TMS.WindowMPO, last(TMS.mixed_steps(q5, d)))
        @test norm(apply(d, r5) - apply(make_mpo(r5, d), r5)) < 1e-12
    end
    z4 = Dephase(reduce(⊗, fill(Z, 4)))(1, 2, 3, 4)
    @test only(last(TMS.mixed_steps(q5, z4))) isa TMS.GateBatch
    # and gates far apart leave the sites in their places
    ρf = random_state(Mixed, System(6, Fermion()), 3)
    @test TMS.ITensorMPS.siteinds(apply(Gate((C ⊗ dag(C))(1, 5)) * Gate(Swap(2, 6)), ρf).state) ==
        TMS.ITensorMPS.siteinds(ρf.state)
    # a product of gates is the gates applied one by one
    ρg = random_state(Mixed, System(8, Fermion()), 3)
    g1, g2 = Gate((C ⊗ dag(C))(2, 6)), Gate(Swap(3, 5))
    g3 = Gate(exp(-0.3im * (dag(C) ⊗ C + dag(dag(C) ⊗ C)))(4, 7))
    @test norm(apply(g1 * g2 * g3, ρg) - apply(g1, apply(g2, apply(g3, ρg)))) < 1e-12
    # bosons of dimension 8, whose gate on two sites is too large to be contracted: its MPO,
    # on the kets then on the bras
    bs = System(2, Boson(8))
    ρb = random_state(Mixed, bs, 2)
    rb = ket_bra_matrix(ρb)
    g = exp(-0.2im * (A ⊗ dag(A) + dag(A) ⊗ A))(1, 2)
    u = lowered(bs, g)
    @test ket_bra_matrix(apply(g, ρb)) ≈ u * rb * u'
end

@testset "Gates on a density matrix, on the kets then on the bras" begin
    # a gate of several sites acts on the kets, then on the bras, through its MPO made to act on
    # the mixed indices: one gate on consecutive sites, its MPO on the sites between otherwise; a
    # gate of one site as its superoperator; another channel on sites apart through its MPO
    function same(ρ, a, k)
        r = ket_bra_matrix(ρ)
        return norm(ket_bra_matrix(apply(a, ρ)) - k * r * k') ≤ 1e-12 * max(1, norm(k * r * k'))
    end
    kinds(ρ, a) = typeof.(last(TMS.mixed_steps(ρ.system, a)))
    qs = System(6, Qubit())
    ρ = random_state(Mixed, qs, 4)
    g4 = exp(-0.3im * (X ⊗ X + Y ⊗ Y + Z ⊗ Z))
    for op in [controlled(X)(2, 5), controlled(X)(3, 4), g4(1, 6), g4(2, 3), H(3),
               controlled(controlled(X))(1, 3, 6), controlled(X)(1, 5) * H(2) * g4(3, 4)]
        @test same(ρ, Gate(op), lowered(qs, op))
    end
    @test kinds(ρ, Gate(controlled(X)(2, 5))) == [TMS.WindowMPO, TMS.WindowMPO]
    @test kinds(ρ, Gate(controlled(X)(3, 4))) == [TMS.GateBatch]
    @test kinds(ρ, Gate(H(3))) == [TMS.GateBatch]
    ps = [ lowered(qs, Proj(Z ⊗ Z => λ)(1, 4)) for λ in (1., -1.) ]
    r = ket_bra_matrix(ρ)
    @test ket_bra_matrix(apply(Dephase(Z ⊗ Z)(1, 4), ρ)) ≈ sum(p * r * p for p in ps)
    @test kinds(ρ, Dephase(Z ⊗ Z)(1, 4)) == [TMS.WindowMPO]
    # fermions, charged or not, and bosons
    hop = exp(-0.4im * (dag(C) ⊗ C + dag(dag(C) ⊗ C)))
    for site in (Fermion(), Fermion(conserve = N))
        fs = System(5, site)
        ρf = mix(random_state(State{Pure}(fs, ["Occ", "Emp", "Occ", "Emp", "Occ"]), 3)) +
             mix(random_state(State{Pure}(fs, ["Emp", "Occ", "Occ", "Emp", "Occ"]), 3))
        for op in [hop(1, 2), hop(2, 4), (C ⊗ N)(1, 4), C(1) * dag(C)(3) + 0.5 * N(2)]
            @test same(ρf, Gate(op), lowered(fs, op))
        end
    end
    bs = System(4, Boson(4))
    ρb = random_state(Mixed, bs, 3)
    bh = exp(-0.2im * (dag(A) ⊗ A + dag(dag(A) ⊗ A)))
    for op in [bh(1, 2), bh(1, 4)]
        @test same(ρb, Gate(op), lowered(bs, op))
    end
end

@testset "Channels applied under a strong symmetry" begin
    ss = System(3, Fermion(conserve = strong(N)))
    ρ = mix(random_state(State{Pure}(ss, ["Occ", "Emp", "Occ"]), 2))
    r = ket_bra_matrix(ρ)
    k, n2 = lowered(ss, dag(C)(1) * C(3)), lowered(ss, N(2))
    @test ket_bra_matrix(apply(Gate(dag(C)(1) * C(3)), ρ)) ≈ k * r * k'
    @test ket_bra_matrix(apply(Dissipator(N)(2), ρ)) ≈ n2 * r * n2 - (n2 * n2 * r + r * n2 * n2) / 2
    @test_throws "conserving it strongly forbids" apply(Reset("Emp")(1), ρ)
    # a gate of several sites moving a charge between them, whose links on the side of the bras
    # take their starred names, where they kept those of the kets
    g = exp(-0.37im * (dag(C) ⊗ C + dag(dag(C) ⊗ C)))
    for op in (g(1, 2), g(2, 1), Swap(1, 2))
        kg = lowered(ss, op)
        @test ket_bra_matrix(apply(op, ρ)) ≈ kg * r * kg'
    end
    ρq = mix(State{Pure}(System(3, Qubit(conserve = strong(N))), ["Up", "Dn", "Up"]))
    @test expect(apply(Swap(1, 2), ρq), [N(1), N(2), N(3)]) ≈ [1, 0, 0]
end

@testset "A channel taking a density matrix to zero on charged sites" begin
    # the null state, as for a pure state, see observables.jl
    ρ = mix(State{Pure}(System(5, Fermion(conserve = N)), ["Emp", "Occ", "Emp", "Occ", "Emp"]))
    @test norm(apply(Gate((C ⊗ C)(1, 2)), ρ)) == 0
    @test norm(apply(Gate(C(1) + C(3)), ρ)) == 0
    @test norm(apply(Gate(dag(C)(1) * C(2)), ρ)) ≈ 1
    # and an operator given by a zero matrix, which ITensors failed to factorise into an MPO,
    # on a density matrix as on a pure state
    z = Operator{2}("Z0", zeros(4, 4))
    ψ = State{Pure}(System(4, Fermion(conserve = N)), ["Occ", "Emp", "Occ", "Emp"])
    @test norm(apply(z(1, 3), ψ)) == 0
    @test norm(apply(z(1, 2), mix(ψ))) == 0
end

@testset "Strings kept apart on a density matrix" begin
    # channels that are not gates between fermionic gates far apart, developed for their
    # superoperator, their strings on the kets and on the bras laid at the end, give the product
    # of the superoperators of their factors
    fs = System(4, Fermion())
    ss = fs.sites
    l = Dissipator(C ⊗ Id) + Evolver(-im * (dag(C) ⊗ C + dag(dag(C) ⊗ C)))
    for a in [Gate(C)(1) * Reset("Emp")(2) * Gate(dag(C))(4),
              Gate(C(4)) * exp(0.3 * l)(2, 3) * Gate(dag(C)(1)),
              Gate(C)(4) * reset_gate(0.4, "Occ")(3) * Gate(C)(1)]
        @test superop(fs, a) ≈ TMS.scalarcoef(a) * prod(superop(fs, f) for f in TMS.prodsubs(a))
    end
    # and the dissipator of such a product, its M read off the strings of both sides
    Φ = Gate(C)(1) * Reset("Emp")(2) * Gate(dag(C))(4)
    S = superop(fs, Φ)
    m = TMS.kraus_sum(S, ss)
    i = one(m)
    @test superop(fs, Dissipator(Φ)) ≈ S - (TMS.sided(m, i, ss) + TMS.sided(i, m, ss)) / 2
end

@testset "A factor holding the identity on a density matrix" begin
    # (Id + Z)(1) * Z(2) stands for Z(2) + Z(1) * Z(2) in the channel and the generator it makes,
    # its part of one site laid on that site, as for a pure state (compact.jl)
    ρ = random_state(Mixed, System(3, Qubit()), 4; rng = Xoshiro(20261013))
    a = ((Id + Z)^0.5)(1) * ((Id + Z)^0.5)(1) * Z(2)
    b = Z(2) + Z(1) * Z(2)
    for op in (Gate, Evolver)
        @test norm(apply(make_mpo(ρ, op(a)), ρ) - apply(make_mpo(ρ, op(b)), ρ)) < 1e-12
    end
end
