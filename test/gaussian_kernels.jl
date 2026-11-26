# The numerical kernels of the Gaussian backend (module GaussianKernels) against the exact
# evaluator of gaussian_dense.jl, up to 10 modes, and by consistency at 200 modes: numbers held
# by their logarithm, pfaffians, Fock states, supports of Gaussian states.

include("gaussian_dense.jl")
using TensorMixedStates.GaussianKernels
using .LA: qr, det, dot, I, eigen, Hermitian

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

gk_antisymmetric(m) = (B = randn(gk_rng, m, m); B - transpose(B))

@testset "evolution and the transport of phases" begin
    # propagator against exp(-iHt)
    for n in (1, 3, 5)
        ψ = dense_gaussian(rand(gk_rng, Bool, n), gk_rng)
        h = gk_antisymmetric(2n)
        φ = exp(-im * 0.7 * dense_hamiltonian(h)) * ψ
        R = propagator(h, 0.7)
        @test R * dense_covariance(ψ) * transpose(R) ≈ dense_covariance(φ) atol = 1e-12
        @test rotate(gk_from_dense(ψ), R).cov ≈ dense_covariance(φ) atol = 1e-12
        @test dense_vector(gk_from_dense(ψ)) ≈ ψ atol = 1e-12
    end
    # components born of one another by monomials, orthogonal ones, both parities, a h at each step
    for n in (2, 4, 6)
        ψ1 = dense_gaussian(zeros(Bool, n), gk_rng)
        xodd = zeros(Bool, n)
        xodd[end] = true
        ψs = [ψ1, dense_monomial(randperm(gk_rng, 2n)[1:4], ψ1), dense_monomial([2], ψ1),
              dense_fock(zeros(Bool, n)), dense_fock(vcat(trues(2), falses(n - 2))),
              dense_gaussian(xodd, gk_rng)]
        s = GaussianSum(gk_from_dense.(ψs), randn(gk_rng, ComplexF64, length(ψs)), gk_rng)
        for step in 1:4
            h, t = gk_antisymmetric(2n), rand(gk_rng)
            evolve!(s, propagator(h, t))
            ψs = [exp(-im * t * dense_hamiltonian(h)) * ψ for ψ in ψs]
            if step == 2
                fix!(s)
            end
        end
        fix!(s)
        # every component equal to the exact one up to one phase common to all
        z = [dot(ψ, dense_vector(g)) for (ψ, g) in zip(ψs, s.comps)]
        @test all(abs(zk - z[1]) < 1e-10 for zk in z)
        @test abs(z[1]) ≈ 1 atol = 1e-10
    end
    @test_throws DimensionMismatch GaussianSum([describe(fock_covariance([false]))], [1, 2], gk_rng)
    # a reflection changes the parity, it is not the rotation of a Gaussian unitary of H
    s = GaussianSum([describe(fock_covariance([false, false]))], [1], gk_rng)
    evolve!(s, [-1.0 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1])
    @test_throws ArgumentError fix!(s)
    @test_throws DimensionMismatch GaussianSum(describe.(fock_covariance.([[false], [false, false]])), [1, 1], gk_rng)
end

@testset "transport of phases at 200 modes, by consistency" begin
    n = 200
    G1 = describe(gk_random_covariance(n, zeros(Bool, n)))
    comps = [G1, apply_monomial([2, 7, 11, 30], G1), describe(fock_covariance(zeros(Bool, n))),
             describe(fock_covariance(vcat(trues(4), falses(n - 4)))), apply_monomial([1], G1),
             apply_monomial([3, 5, 9], G1)]
    G0 = [overlap(a, b) for a in comps, b in comps]
    # a transverse field Ising chain in Majoranas
    h = zeros(2n, 2n)
    for j in 1:n
        h[2j-1, 2j] = 0.7
        if j < n
            h[2j, 2j+1] = 1.0
        end
    end
    h = h - transpose(h)
    R = propagator(h, 0.01)
    s = GaussianSum(comps, ones(length(comps)), gk_rng)
    for step in 1:200
        evolve!(s, R)
        if step == 100
            fix!(s)
        end
    end
    fix!(s)
    for k in eachindex(comps), l in eachindex(comps)
        o = overlap(s.comps[k], s.comps[l])
        if abs(G0[k, l]) < 1e-12
            @test abs(o) < 1e-12
        else
            @test gk_dist(o, G0[k, l]) < 1e-9
        end
    end
    # between parities: ⟨UG₁|γ_b|U γ₁G₁⟩ = Σ_c R_bc ⟨G₁|γ_c γ₁|G₁⟩
    u = im .* G1.cov[:, 1]
    u[1] = 1
    w = R^200 * u
    b = argmax(abs.(w))
    @test gk_dist(cross_term(s.comps[1], [b], s.comps[5]), LogComplex(w[b])) < 1e-9
end

# the projector on mode j occupied or empty
gk_projector(j, occupied, ψ) = dense_site(occupied ? ComplexF64[0 0; 0 1] : ComplexF64[1 0; 0 0], j, ψ)

@testset "measurements of occupations" begin
    kept = changed = 0
    for n in 1:6, _ in 1:(n <= 4 ? 6 : 2), j in 1:n, occupied in (false, true)
        ψ = dense_gaussian(rand(gk_rng, Bool, n), gk_rng)
        g = gk_from_dense(ψ)
        φ = gk_projector(j, occupied, ψ)
        if n == 1
            # a single mode is in a Fock state
            @test_throws ArgumentError project(g, j, !g.x[1])
            continue
        end
        p, h = project(g, j, occupied)
        @test p ≈ dot(φ, φ) atol = 1e-12
        @test sqrt(p) .* dense_vector(h) ≈ φ atol = 1e-11
        @test exp(2 * h.amp.logabs) >= 2.0^-n * (1 - 1e-12)
        if h.x == g.x
            kept += 1
        else
            changed += 1
        end
    end
    @test kept > 0 && changed > 0
    # certainty, and an outcome of probability 1e-8
    x = [false, true, false]
    p, h = project(describe(fock_covariance(x)), 2, true)
    @test p == 1 && h.x == x
    @test_throws ArgumentError project(describe(fock_covariance(x)), 2, false)
    ψ = dense_rotate(dense_fock(zeros(Bool, 4)), [(1, 3, 2e-4)])
    φ = gk_projector(1, true, ψ)
    p, h = project(gk_from_dense(ψ), 1, true)
    @test p ≈ 1e-8 rtol = 1e-6
    @test sqrt(p) .* dense_vector(h) ≈ φ rtol = 1e-8
    # superpositions, after an evolution
    for n in (3, 5)
        ψ1 = dense_gaussian(zeros(Bool, n), gk_rng)
        ψs = [ψ1, dense_monomial(randperm(gk_rng, 2n)[1:2], ψ1), dense_monomial([3], ψ1),
              dense_fock(zeros(Bool, n)), dense_gaussian(vcat(true, falses(n - 1)), gk_rng)]
        a = randn(gk_rng, ComplexF64, length(ψs))
        s = GaussianSum(gk_from_dense.(ψs), a, gk_rng)
        h, t = gk_antisymmetric(2n), rand(gk_rng)
        evolve!(s, propagator(h, t))
        Ψ = exp(-im * t * dense_hamiltonian(h)) * sum(a .* ψs)
        @test ComplexF64(norm2(s)) ≈ dot(Ψ, Ψ) rtol = 1e-11
        project!(s, 2, true)
        Φ = gk_projector(2, true, Ψ)
        @test ComplexF64(norm2(s)) ≈ dot(Φ, Φ) rtol = 1e-11
        Φs = sum(ComplexF64(c) .* dense_vector(g) for (c, g) in zip(s.coefs, s.comps))
        z = dot(Φ, Φs) / dot(Φ, Φ)
        @test abs(z) ≈ 1 atol = 1e-10
        @test Φs ≈ z .* Φ atol = 1e-10
    end
    # a component whose mode 1 is empty, but for the rounding of its covariance, adds no spurious
    # weight of order √eps: 1 - M₁₂ ~ 1e-16 > 0 in some of the draws, which (1 - M₁₂)/2 would take as
    # a probability
    n = 4
    rounded = 0
    for _ in 1:6
        rots = [r for r in dense_rotations(n, gk_rng) if r[1] > 2 && r[2] > 2]
        ψs = [dense_gaussian(zeros(Bool, n), gk_rng), dense_rotate(dense_fock(zeros(Bool, n)), rots)]
        s = GaussianSum(gk_from_dense.(ψs), [1, 1], gk_rng)
        rounded += s.comps[2].cov[1, 2] < 1
        project!(s, 1, true)
        Φ = gk_projector(1, true, sum(ψs))
        @test ComplexF64(norm2(s)) ≈ dot(Φ, Φ) rtol = 1e-12
        Φs = sum(ComplexF64(c) .* dense_vector(g) for (c, g) in zip(s.coefs, s.comps))
        @test Φs ≈ (dot(Φ, Φs) / dot(Φ, Φ)) .* Φ atol = 1e-12
    end
    @test rounded > 0
    # a component that cannot be found occupied is dropped, all of them is an error
    s = GaussianSum([describe(fock_covariance([false, false])), describe(fock_covariance([true, true]))],
                    [1, 1], gk_rng)
    project!(s, 1, true)
    @test length(s.comps) == 1 && s.comps[1].x == [true, true]
    @test_throws ArgumentError project!(s, 2, false)
end

@testset "measurements at 200 modes, by consistency" begin
    n = 200
    g = describe(gk_random_covariance(n, zeros(Bool, n)))
    for j in (1, 77, 200), occupied in (false, true)
        # ⟨ψ|ψ′⟩ = ⟨ψ|Πψ⟩/√p = √p
        p, h = project(g, j, occupied)
        @test gk_dist(overlap(g, h), LogComplex(log(p) / 2, 1)) < 1e-10
    end
end
