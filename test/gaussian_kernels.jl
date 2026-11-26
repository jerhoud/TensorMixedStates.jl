# The numerical kernels of the Gaussian backend (module GaussianKernels) against the exact
# evaluator of gaussian_dense.jl, up to 10 modes, and by consistency at 200 modes: numbers held
# by their logarithm, pfaffians, Fock states, supports of Gaussian states.

include("gaussian_dense.jl")
using TensorMixedStates.GaussianKernels
using .LA: qr, det, dot, I

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

# the description of a dense Gaussian state ψ, with its phase
function gk_from_dense(ψ)
    g = describe(dense_covariance(ψ))
    PureGaussian(g.cov, g.x, LogComplex(ψ[dense_index(g.x)]))
end

# up to 4 distinct Majoranas in random order
gk_monomial(n) = randperm(gk_rng, 2n)[1:rand(gk_rng, 0:min(4, 2n))]

# the relative distance of two nonzero LogComplex
gk_dist(a, b) = abs(a.logabs - b.logabs) + abs(a.phase - b.phase)

@testset "overlaps, monomials and Wick" begin
    ewick = 0.0
    for n in 1:10, _ in 1:(n <= 6 ? 10 : 3)
        x1 = rand(gk_rng, Bool, n)
        x2 = rand(gk_rng, Bool, n)
        x2[1] = isodd(count(x1)) ⊻ isodd(count(x2[2:end]))
        x3 = copy(x2)
        x3[end] = !x3[end]
        ψ1, ψ2, ψ3 = dense_gaussian(x1, gk_rng), dense_gaussian(x2, gk_rng), dense_gaussian(x3, gk_rng)
        g1, g2, g3 = gk_from_dense(ψ1), gk_from_dense(ψ2), gk_from_dense(ψ3)
        @test ComplexF64(overlap(g1, g2)) ≈ dot(ψ1, ψ2) atol = 1e-12
        @test ComplexF64(overlap(g1, g1)) ≈ 1 atol = 1e-12
        @test iszero(overlap(g1, g3))
        for _ in 1:4
            α = gk_monomial(n)
            h = apply_monomial(α, g2)
            φ = dense_monomial(α, ψ2)
            @test h.cov ≈ dense_covariance(φ) atol = 1e-12
            @test ComplexF64(h.amp) ≈ φ[dense_index(h.x)] atol = 1e-12
            @test ComplexF64(cross_term(g1, α, g2)) ≈ dot(ψ1, φ) atol = 1e-12
            @test ComplexF64(cross_term(g1, α, g3)) ≈ dot(ψ1, dense_monomial(α, ψ3)) atol = 1e-12
            @test ComplexF64(expect_monomial(g1.cov, α)) ≈ dot(ψ1, dense_monomial(α, ψ1)) atol = 1e-12
        end
        W = transition_matrix(g1.cov, g2.cov)
        o = dot(ψ1, ψ2)
        for a in 1:2n, b in 1:2n
            ewick = max(ewick, abs(o * W[a, b] - dot(ψ1, dense_majorana(a, dense_majorana(b, ψ2)))))
        end
    end
    @test ewick < 1e-11
    # orthogonal states: Fock states, then a Gaussian unitary V on two of them
    for (x, y) in (([false, false], [true, true]), ([false, true, false, true], [true, false, true, false]))
        @test abs(overlap(describe(fock_covariance(x)), describe(fock_covariance(y)))) < 1e-15
    end
    n = 6
    x, y = rand(gk_rng, Bool, n), rand(gk_rng, Bool, n)
    y[1:2] .= .!x[1:2]
    y[3:end] .= x[3:end]
    V = dense_rotations(n, gk_rng)
    ψx, ψy = dense_rotate(dense_fock(x), V), dense_rotate(dense_fock(y), V)
    @test abs(ComplexF64(overlap(gk_from_dense(ψx), gk_from_dense(ψy)))) < 1e-12
    # nearly orthogonal: |⟨ψx|ψε⟩| = sin(θ/2) ≈ 5e-7
    ψε = dense_rotate(dense_rotate(dense_fock(y), [(1, 3, 1e-6)]), V)
    o = dot(ψx, ψε)
    @test abs(o) ≈ sin(5e-7)
    @test abs(ComplexF64(overlap(gk_from_dense(ψx), gk_from_dense(ψε))) - o) < 1e-8 * abs(o)
    # edge cases
    g = gk_from_dense(dense_gaussian(rand(gk_rng, Bool, 3), gk_rng))
    @test ComplexF64(expect_monomial(g.cov, Int[])) == 1
    @test cross_term(g, Int[], g) == overlap(g, g)
    @test_throws ArgumentError expect_monomial(g.cov, [1, 2, 1, 3])
    @test_throws DimensionMismatch overlap(g, describe(fock_covariance([false, false])))
end

@testset "overlaps at 200 modes, by consistency" begin
    n = 200
    gs = [describe(gk_random_covariance(n, zeros(Bool, n))) for _ in 1:3]
    gs = [PureGaussian(g.cov, g.x, LogComplex(g.amp.logabs, cis(2π * rand(gk_rng)))) for g in gs]
    o01, o12, o20 = overlap(gs[1], gs[2]), overlap(gs[2], gs[3]), overlap(gs[3], gs[1])
    @test gk_dist(o01, conj(overlap(gs[2], gs[1]))) < 1e-9
    # |⟨ψ₁|ψ₂⟩|² = 2⁻ⁿ |Pf(M₁ + M₂)|
    @test abs(2 * o12.logabs - (logpfaffian(gs[2].cov + gs[3].cov).logabs - n * log(2))) < 1e-9
    # the product of the three overlaps, Bravyi and Gosset, eq. (24), even states
    Id = Matrix{Float64}(I, 2n, 2n)
    T = [im*gs[1].cov -Id Id; Id im*gs[2].cov -Id; -Id Id im*gs[3].cov]
    @test gk_dist(o01 * o12 * o20, LogComplex(-n * log(4), im^n) * logpfaffian(T)) < 1e-9
    α = [7, 2, 150, 33]
    @test gk_dist(cross_term(gs[1], α, gs[2]), conj(cross_term(gs[2], reverse(α), gs[1]))) < 1e-9
end
