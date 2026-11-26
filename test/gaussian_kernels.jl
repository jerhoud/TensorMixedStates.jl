# The numerical kernels of the Gaussian backend (module GaussianKernels) against the exact
# evaluator of gaussian_dense.jl, up to 10 modes, and by consistency at 200 modes: numbers held
# by their logarithm, pfaffians, Fock states, supports of Gaussian states.

include("gaussian_dense.jl")
using TensorMixedStates.GaussianKernels
using .LA: qr, det, dot

const gk_rng = Random.Xoshiro(20261126)

# a random pure covariance of n modes and the parity of x
function gk_random_covariance(n, x)
    O = Matrix(qr(randn(gk_rng, 2n, 2n)).Q)
    if det(O) < 0
        O[:, 1] .*= -1
    end
    M = O * fock_covariance(x) * transpose(O)
    (M - transpose(M)) / 2
end

@testset "LogComplex" begin
    z1, z2 = 0.3 - 1.2im, -2.0 + 0.5im
    a, b = LogComplex(z1), LogComplex(z2)
    @test ComplexF64(a * b) ≈ z1 * z2
    @test ComplexF64(a / b) ≈ z1 / z2
    @test ComplexF64(a + b) ≈ z1 + z2
    @test ComplexF64(a - b) ≈ z1 - z2
    @test ComplexF64(conj(a)) ≈ conj(z1)
    @test ComplexF64(2im * a / 3) ≈ 2im * z1 / 3
    @test abs(a) ≈ abs(z1)
    @test LogComplex(-2.0).phase == -1
    @test iszero(LogComplex(0)) && iszero(LogComplex(-Inf, im)) && iszero(LogComplex(1.0, 0))
    @test iszero(a * zero(LogComplex)) && a + zero(LogComplex) == a && zero(LogComplex) + a == a
    @test iszero(a - a)
    @test ComplexF64(zero(LogComplex)) == 0
    @test_throws DivideError inv(zero(LogComplex))
    @test_throws DomainError LogComplex(Inf, 1)
    @test_throws DomainError LogComplex(NaN, 1)
    # beyond the range of Float64
    t = LogComplex(-1000.0, im)
    @test (t * t).logabs == -2000 && (t * t).phase == -1
    @test (t + t).logabs ≈ -1000 + log(2) && (t + t).phase == im
    @test (t / LogComplex(-3000.0, 1)).logabs == 2000
end

@testset "pfaffians" begin
    for m in 0:2:10, T in (Float64, ComplexF64)
        B = randn(gk_rng, T, m, m)
        A = B - transpose(B)
        @test ComplexF64(logpfaffian(A)) ≈ dense_pfaffian(A) rtol = 1e-12
    end
    @test iszero(logpfaffian(zeros(3, 3)))
    @test iszero(logpfaffian(zeros(4, 4)))
    B = randn(gk_rng, ComplexF64, 10, 10)
    A = B - transpose(B)
    @test ComplexF64(logpfaffian(A, (4, 2))) ≈ dense_pfaffian(A) rtol = 1e-12
    @test ComplexF64(logpfaffian(A, (10,))) ≈ dense_pfaffian(A) rtol = 1e-12
    A[1, 2] = A[2, 1] = 0
    @test_throws ArgumentError logpfaffian(A, (2,))
    # 2n = 400, |Pf| = 1e-600 beyond Float64: blocks t [0 1; -1 0] turned by O, Pf = det(O) tⁿ
    O = Matrix(qr(randn(gk_rng, 400, 400)).Q)
    A = O * (1e-3 * fock_covariance(zeros(Bool, 200))) * transpose(O)
    p = logpfaffian(A)
    @test p.logabs ≈ 200 * log(1e-3) rtol = 1e-12
    @test p.phase ≈ sign(det(O))
end

@testset "Fock states" begin
    for n in 1:4, k in 0:2^n-1
        x = dense_bits(k, n)
        @test fock_covariance(x) ≈ dense_covariance(dense_fock(x)) atol = 1e-15
        for b in 1:2n
            coef, y = majorana_on_fock(b, x)
            @test dense_majorana(b, dense_fock(x)) ≈ coef .* dense_fock(y)
        end
    end
end

@testset "supports and descriptions" begin
    for n in 1:10, _ in 1:(n <= 6 ? 20 : 5)
        x0 = rand(gk_rng, Bool, n)
        ψ = dense_gaussian(x0, gk_rng)
        g = describe(dense_covariance(ψ))
        p = abs2(ψ[dense_index(g.x)])
        @test nmodes(g) == n
        @test p >= 2.0^-n * (1 - 1e-12)
        @test exp(2 * g.amp.logabs) ≈ p rtol = 1e-12
        @test g.amp.phase == 1
        @test iseven(count(g.x)) == iseven(count(x0))
    end
    # mode 1 empty with probability 1/2, so that either choice is valid
    ψ = ComplexF64[1, 0, 0, 1] / sqrt(2)
    g = describe(dense_covariance(ψ))
    @test g.x in ([false, false], [true, true])
    @test exp(2 * g.amp.logabs) ≈ 1 / 2
    # 200 modes: the probability of the support by its chain of conditional probabilities
    # against its pfaffian, |⟨x|ψ⟩|² = 2⁻ⁿ |Pf(M + M_x)|
    n = 200
    M = gk_random_covariance(n, rand(gk_rng, Bool, n))
    g = describe(M)
    @test 2 * g.amp.logabs >= -n * log(2)
    @test 2 * g.amp.logabs ≈ logpfaffian(M + fock_covariance(g.x)).logabs - n * log(2) rtol = 1e-12
end
