# The evolutions in time and the thermal states: evolve by Tdvp and ApproxW, on cases whose exact
# solution is known, dissipative ones included, and against the exact dense evolution of small
# chains, pure and mixed, with the order at which each method converges; the evolvers that
# depend on time, the number of the first step and the limits of each step, the MPOs a PreMPO
# keeps, the parameters of the methods and what they refuse; the approximations WI and WII
# against their closed forms; and thermalize against e^{-βH}/Z. The searches are in searches.jl,
# the algorithms of one's own in extension.jl.
#
# The references are dense matrices written here with kron, sharing no code with the package,
# or exact values computed by the scripts of test/reference. A tolerance says how
# accurate a method is on its problem: as tight as the method allows, but not at the floating
# point floor, where rounding differs with the BLAS build.
#
# The phases Evolve and Thermalize are tested here as far as their computation goes, their
# resume in runner.jl; the gates and the channels applied to states are in channels.jl.

"""
the dense generator of ``\\rho \\mapsto -i[h, \\rho] + \\sum_l l \\rho l^\\dagger - \\{l^\\dagger
l, \\rho\\}/2`` on `vec(ρ)`, stacked by columns, so that ``vec(a \\rho b) = (b^T \\otimes a)
vec(\\rho)``
"""
function dense_lindbladian(h, ls = [])
    d = size(h, 1)
    id = Matrix{ComplexF64}(LA.I, d, d)
    g = -im * (kron(id, h) - kron(transpose(h), id))
    for l in ls
        g += kron(conj(l), l) - 0.5 * (kron(id, l' * l) + kron(transpose(l' * l), id))
    end
    return g
end

"""
the solution at `t1` of ``dy/dt = f(t, y)`` from `y0` at `t0`, by the classical Runge-Kutta
method of order 4 in `nsteps` steps: a reference independent of the integrators of the package
"""
function runge_kutta(f, y0, t0, t1, nsteps)
    dt = (t1 - t0) / nsteps
    y = y0
    for k in 0:nsteps-1
        s = t0 + k * dt
        k1 = f(s, y)
        k2 = f(s + dt / 2, y + dt / 2 * k1)
        k3 = f(s + dt / 2, y + dt / 2 * k2)
        k4 = f(s + dt, y + dt * k3)
        y = y + dt / 6 * (k1 + 2k2 + 2k3 + k4)
    end
    return y
end

@testset "Dissipative evolve" begin
    lim = Limits(maxdim = 10, cutoff = 1e-15)
    # a qubit damped towards Up while it precesses, from X+: X(t) = e^{-γt/2} cos 2t,
    # Y(t) = e^{-γt/2} sin 2t and Z(t) = 1 - e^{-γt}
    γ = 0.5
    ev = -im * Z(1) + Dissipator(sqrt(γ) * Sp)(1)
    ρ = State{Mixed}(System(2, Qubit()), "X+")
    exact = [exp(-γ / 2) * cos(2), exp(-γ / 2) * sin(2), 1 - exp(-γ)]
    for (algo, nsteps, tol) in [(Tdvp(), 10, 1e-12), (ApproxW(order = 4, w = 2), 10, 1e-12),
                                (ApproxW(order = 4, w = 1), 10, 2e-6),
                                (ApproxW(order = 1, w = 1), 100, 0.03)]
        st = evolve(ev, 1., ρ; algo, nsteps, limits = lim)
        @test real(expect(st, [X(1), Y(1), Z(1)]; complex = true)) ≈ exact atol = tol
        @test trace(st) ≈ 1 atol = 1e-12
    end

    # two qubits coupled, one dephased and the other damped, against the dense Lindbladian
    x, y, z, o = [0. 1. ; 1. 0.], [0. -im ; im 0.], [1. 0. ; 0. -1.], [1. 0. ; 0. 1.]
    h = kron(x, x) + 0.7 * kron(z, o)
    g = dense_lindbladian(h, [sqrt(0.3) * kron(z, o), sqrt(0.2) * kron(o, [0. 1. ; 0. 0.])])
    ψ = kron([1., 1.] / sqrt(2), [1., 0.])
    r = reshape(exp(g) * vec(ψ * ψ'), 4, 4)
    tr = LA.tr
    ref = real([tr(r * kron(x, o)), tr(r * kron(o, z)), tr(r * kron(y, x))])
    ev2 = -im * (X(1) * X(2) + 0.7 * Z(1)) + Dissipator(sqrt(0.3) * Z)(1) +
          Dissipator(sqrt(0.2) * Sp)(2)
    ρ2 = State{Mixed}(System(2, Qubit()), ["X+", "Up"])
    for (algo, tol) in [(Tdvp(), 1e-11), (ApproxW(order = 4, w = 2), 1e-7),
                        (ApproxW(order = 4, w = 1), 1e-7)]
        st = evolve(ev2, 1., ρ2; algo, nsteps = 20, limits = lim)
        @test real(expect(st, [X(1), Z(2), Y(1) * X(2)]; complex = true)) ≈ ref atol = tol
    end
end

@testset "A step of tdvp whatever the norm of the state" begin
    # KrylovKit takes its precision as absolute, and a vector below it for a fixed point: a
    # state of norm 1e-13 did not move, and a density matrix of trace one, whose norm is the
    # square root of its purity, lost precision as it neared the fully mixed state on many sites
    ρ = random_state(Mixed, System(4, Qubit()), 2; rng = Xoshiro(20261019))
    ev = -im * sum(X(i) * X(i + 1) + Z(i) * Z(i + 1) for i in 1:3) + Dissipator(sqrt(0.1) * Sm)(2)
    r1 = evolve(ev, 0.5, ρ; nsteps = 5)
    @test norm(1e13 * evolve(ev, 0.5, 1e-13 * ρ; nsteps = 5) - r1) < 1e-10 * norm(r1)
end

@testset "Each method against the exact evolution" begin
    # a chain of three qubits, coupled between neighbours and between the two ends across the
    # middle one, in fields along Z and X, evolved for a time 1 from a product state, pure and
    # under a decay and a dephasing, against the exponential of its dense Hamiltonian or
    # Lindbladian, with no truncation, so that what is left is the error of the method. Tdvp,
    # at full bond dimension, is exact; ApproxW converges at its order, its error divided by
    # about 2^order when the step is halved, WII being the more accurate
    n, t = 3, 1.
    sys = System(n, Qubit())
    x, y, z = [0. 1. ; 1. 0.], [0. -im ; im 0.], [1. 0. ; 0. -1.]
    sm = [0. 0. ; 1. 0.]
    h = sum(X(i) * X(i + 1) for i in 1:n-1) + 0.5 * Z(1) * Z(3) + 0.7 * sum(Z(i) for i in 1:n) +
        0.3 * X(1)
    hd = sum(on_qubit(x, i, n) * on_qubit(x, i + 1, n) for i in 1:n-1) +
         0.5 * on_qubit(z, 1, n) * on_qubit(z, 3, n) + 0.7 * sum(on_qubit(z, i, n) for i in 1:n) +
         0.3 * on_qubit(x, 1, n)
    jumps = [sqrt(0.4) * on_qubit(sm, 1, n), sqrt(0.3) * on_qubit(z, 3, n)]
    dissipators = Dissipator(sqrt(0.4) * Sm)(1) + Dissipator(sqrt(0.3) * Z)(3)
    ψ = State{Pure}(sys, ["Up", "+", "Dn"])
    v = dense_vector(ψ)
    exact_pure(hd) = dense_state(sys, exp(-im * t * hd) * v)
    exact_mixed(hd, ls) = dense_state(sys, reshape(exp(t * dense_lindbladian(hd, ls)) * vec(v * v'),
                                                   2^n, 2^n))
    lim = Limits(cutoff = 0)
    for (st, ev, exact, tol) in [(ψ, -im * h, exact_pure(hd), 1e-12),
                                 (mix(ψ), -im * h + dissipators, exact_mixed(hd, jumps), 1e-11)]
        @test norm(evolve(ev, t, st; nsteps = 4, limits = lim) - exact) < tol
        for order in 1:4
            err = Dict((w, nsteps) => norm(evolve(ev, t, st; algo = ApproxW(; order, w), nsteps,
                                                  limits = lim) - exact)
                       for w in 1:2, nsteps in (8, 16))
            for w in 1:2
                @test abs(log2(err[w, 8] / err[w, 16]) - order) < 0.4
            end
            @test err[2, 16] < err[1, 16]
        end
    end
    # terms of one site, whose exponential WII takes exactly, as it takes the dissipators of one
    # site, where WI is of order one
    h1 = 0.7 * sum(Z(i) for i in 1:n) + 0.3 * X(1) + 0.5 * Y(2)
    hd1 = 0.7 * sum(on_qubit(z, i, n) for i in 1:n) + 0.3 * on_qubit(x, 1, n) +
          0.5 * on_qubit(y, 2, n)
    for (st, ev, exact) in [(ψ, -im * h1, exact_pure(hd1)),
                            (mix(ψ), -im * h1 + dissipators, exact_mixed(hd1, jumps))]
        w2(w) = evolve(ev, t, st; nsteps = 1, algo = ApproxW(order = 1, w = w), limits = lim)
        @test norm(w2(2) - exact) < 1e-12
        @test norm(w2(1) - exact) > 1e-2
    end
end

@testset "Time dependent evolve" begin
    # h(t) = Z(1) + t * Z(2): two commuting terms of one site carrying different time functions,
    # so that each spin precesses at a rate known exactly and a function matched with the wrong
    # term shows at once. The middle of the step, and the two Gauss points of the Magnus
    # integrator of ApproxW of order 4, are exact on functions of degree one, and so are tdvp
    # and WII on these terms: what is left is rounding, and on a mixed state the cutoff of the
    # hundred MPOs ApproxW builds, about 1e-13
    terms = [-im * Z(1) => t -> 1., -im * Z(2) => t -> t]
    lim = Limits(maxdim = 10, cutoff = 1e-15)
    precessed = [cos(2), sin(2), cos(1), sin(1)]
    for (algo, tol) in [(Tdvp(), 1e-13), (ApproxW(order = 4, w = 2), 1e-12)], R in (Pure, Mixed)
        st = evolve(terms, 1., State{R}(System(2, Qubit()), "X+"); algo, nsteps = 100,
                    limits = lim)
        @test real(expect(st, [X(1), Y(1), X(2), Y(2)])) ≈ precessed atol = tol
    end
    # read from `time`: from 1, Z(2) turns by 3 over a time 1
    st = evolve(terms, 1., State{Pure}(System(2, Qubit()), "X+"); nsteps = 10, time = 1.)
    @test real(expect(st, X(2))) ≈ cos(3) atol = 1e-12
    # the vector given as its PreMPO, on a mixed state, so that each operator of the vector is
    # lifted to the mixed representation
    ρ = mix(State{Pure}(System(2, Qubit()), "X+"))
    st = evolve(PreMPO(ρ, terms), 1., ρ; nsteps = 100)
    @test real(expect(st, X(1))) ≈ cos(2.) atol = 1e-13
    @test real(expect(st, X(2))) ≈ cos(1.) atol = 1e-13
end

@testset "Time dependent evolvers against the exact evolution" begin
    # two qubits coupled, in two fields whose time functions do not commute with the coupling,
    # from a time 0.3, against a Runge-Kutta integration of the dense equation in many steps,
    # less than 1e-12 away from the exact solution. The functions are read at the middle of
    # each step by Tdvp, of order two then, exact otherwise on two sites, and by ApproxW of
    # orders 1 and 2, and from order 3 at the two Gauss points of the commutator-free Magnus
    # integrator of order 4: each method keeps its order, where the middle of the step left
    # order 4 at order two in 1.x
    x, z, o = [0. 1. ; 1. 0.], [1. 0. ; 0. -1.], [1. 0. ; 0. 1.]
    sm = [0. 0. ; 1. 0.]
    ops = [X(1) * X(2), Z(1) + Z(2), X(1) + X(2)]
    fs = [t -> 1., t -> cos(3t), t -> sin(2t)^2]
    dense_ops = [kron(x, x), kron(z, o) + kron(o, z), kron(x, o) + kron(o, x)]
    hd(t) = sum(f(t) * m for (f, m) in zip(fs, dense_ops))
    jumps = [sqrt(0.2) * kron(sm, o), sqrt(0.2) * kron(o, sm)]
    sys = System(2, Qubit())
    ψ = State{Pure}(sys, ["Up", "+"])
    v = dense_vector(ψ)
    t0, t = 0.3, 1.
    exact_pure = dense_state(sys, runge_kutta((s, y) -> -im * hd(s) * y, v, t0, t0 + t, 10000))
    lind = [-im * ops[1] + 0.2 * (Dissipator(Sm)(1) + Dissipator(Sm)(2)), -im * ops[2],
            -im * ops[3]]
    r = runge_kutta((s, y) -> dense_lindbladian(hd(s), jumps) * y, vec(v * v'), t0, t0 + t, 10000)
    exact_mixed = dense_state(sys, reshape(r, 4, 4))
    lim = Limits(cutoff = 0)
    for (st, ev, exact) in [(ψ, -im .* ops, exact_pure), (mix(ψ), lind, exact_mixed)]
        # the first function is constant, and so is its term
        terms = [ev[1], ev[2] => fs[2], ev[3] => fs[3]]
        err(algo, nsteps) =
            norm(evolve(terms, t, st; algo, nsteps, time = t0, limits = lim) - exact)
        @test abs(log2(err(Tdvp(), 16) / err(Tdvp(), 32)) - 2) < 0.4
        for order in 1:4
            @test abs(log2(err(ApproxW(; order), 16) / err(ApproxW(; order), 32)) - order) < 0.4
        end
    end
end

@testset "The terms of an evolver depending on time" begin
    # each term carries its function, which leaves no count to match: too few functions raised a
    # BoundsError on an internal vector in 1.x, and too many were ignored. Its vectors paired are
    # refused, alone, with a PreMPO, or for a single term
    hs = [-im * Z(1), -im * Z(2)]
    st = State{Pure}(System(2, Qubit()), "X+")
    msg = "where 1.x paired vectors, [A, B] => [f, g]"
    @test_throws msg evolve(hs => [t -> 1., t -> 1.], 0.1, st; nsteps = 1)
    @test_throws msg evolve(PreMPO(st, hs) => [t -> 1., t -> 1.], 0.1, st; nsteps = 1)
    @test_throws msg evolve(-im * Z(1) => [t -> 1.], 0.1, st; nsteps = 1)
    @test_throws "3 is not a term of an evolver" evolve([-im * Z(1), 3], 0.1, st; nsteps = 1)
    # terms with no function are their sum, given as they are or as their PreMPO
    ref = evolve(sum(hs), 0.4, st; nsteps = 2)
    @test norm(evolve(hs, 0.4, st; nsteps = 2) - ref) < 1e-12
    @test norm(evolve(PreMPO(st, hs), 0.4, st; nsteps = 2) - ref) < 1e-12
    # a single term, in a vector or not, and the PreMPO keeping the function
    for ev in (-im * Z(1) => t -> 2., [-im * Z(1) => t -> 2.], PreMPO(st, -im * Z(1) => t -> 2.))
        @test real(expect(evolve(ev, 0.4, st; nsteps = 1), X(1))) ≈ cos(1.6)
    end
    # which an Evolve phase takes as well
    run(evolver) = runTMS(SimData(phases = [CreateState{Pure}(2, Qubit(), "X+"),
                                            Evolve(; evolver, duration = 0.4, time_step = 0.1)]);
                          output = devnull)
    @test real(expect(run(-im * Z(1) => t -> 2.).state, X(1))) ≈ cos(1.6)
    @test real(expect(run([-im * Z(2), -im * Z(1) => t -> 2.]).state, X(1))) ≈ cos(1.6)
    @test_throws msg run([-im * Z(1)] => [t -> 1.])
end

@testset "Time functions take real values" begin
    # a complex value multiplied ρ A† by itself rather than by its conjugate in 1.x, which gave a
    # mixed state of complex trace; it is refused on pure states as well
    hs = -im * [X(1), -Y(1), Z(1) * Z(2) + X(2)]
    ψ = State{Pure}(System(2, Qubit()), ["Up", "X+"])
    for st in (ψ, mix(ψ)), algo in (Tdvp(), ApproxW(order = 2))
        @test_throws "time functions take real values" evolve(
            [hs[1] => t -> cis(3t), hs[2] => t -> 0., hs[3]], 1., st; nsteps = 1, algo)
    end
    # the same drive written with its real and imaginary parts: evolving the mixed state is
    # mixing the evolved pure state, the Hamiltonian not commuting at different times. Tdvp is
    # exact at full bond dimension, which a cutoff of zero keeps, and the tolerance of ApproxW
    # is its error measured here
    terms = [hs[1] => t -> cos(3t), hs[2] => t -> sin(3t), hs[3]]
    for (algo, lim, tol) in [(Tdvp(), Limits(cutoff = 0), 1e-11),
                             (ApproxW(order = 4), Limits(), 1e-7)]
        ev(st) = evolve(terms, 1., st; algo, nsteps = 50, limits = lim)
        @test norm(ev(mix(ψ)) - mix(ev(ψ))) < tol
    end
end

@testset "Time dependent terms of several sites" begin
    # the coefficient of a term goes once into the MPO, whatever the number of sites it spans:
    # laid on each of its pieces, a term of k sites took it to the power k in 1.x, which the
    # terms of one site above cannot show. Folding the coefficients into the terms gives the MPO
    # the time functions must give
    for R in (Pure, Mixed)
        ψ = random_state(Pure, System(4, Qubit()), 4)
        st = R == Pure ? ψ : mix(ψ)
        hs = [-im * X(1) * Z(2), -im * Z(1) * Y(2) * X(4), -im * Z(3)]
        c = [0.5, 0.3, 0.7]
        pre = PreMPO(st, hs)
        folded = PreMPO(st, sum(c .* hs))
        @test norm(prod(make_mpo(pre, c)) - prod(make_mpo(folded))) < 1e-12
        for w in (TMS.make_approx_W1, TMS.make_approx_W2)
            @test norm(prod(w(pre, 0.1, c)) - prod(w(folded, 0.1))) < 1e-12
        end
    end
    # and through an evolution: -i X(1)X(2)/2 over a time 1 turns Z(1) by an angle 1
    st = evolve(-im * X(1) * X(2) => t -> 0.5, 1., State{Pure}(System(2, Qubit()), "Up");
                nsteps = 10, limits = Limits(maxdim = 4, cutoff = 1e-14))
    @test real(expect(st, Z(1))) ≈ cos(1.) atol = 1e-12
end

@testset "A constant term in the approximation WII" begin
    # a constant term is laid on the first site as a delta, whose diagonal storage ITensors
    # could not add a tensor of another element type to, nor, on a charged system, add anything
    # to or exponentiate: WII failed on each of these operators in 1.x. A constant commutes with
    # everything, so all it does is multiply WII by its exponential
    for sys in (System(2, Qubit()), System(2, Qubit(conserve = N)))
        st = State{Pure}(sys, "Dn")
        w(op, coefs...) = prod(TMS.make_approx_W2(PreMPO(st, op), 0.1, coefs...))
        @test norm(w(-im * N(1) + 0.5 * Id(1)) - exp(0.05) * w(-im * N(1))) < 1e-12
        @test norm(w(-im * (N(2) + 0.5 * Id(1))) - exp(-0.05im) * w(-im * N(2))) < 1e-12
        @test norm(w([Id(1), -im * N(1)], [0.7, 0.3]) - exp(0.07) * w([-im * N(1)], [0.3])) < 1e-12
    end
end

@testset "WII keeps the products of terms that cross no same link" begin
    # WII is the approximation of Zaletel et al., whose blocks put the exponential of the terms
    # of one site around every piece, transports included, and take a closing and an opening
    # on the same site in both orders. On these operators it has a closed form, which the blocks
    # of 1.x before 1.5 missed by a term of order τ²
    q = Qubit()
    x, y, z, o = (matrix(a, q) for a in (X, Y, Z, Id))
    id8 = kron(o, o, o)
    τ, h, J, K = 0.1, 0.7, 0.4, 0.3
    w(st, op) = mpo_matrix(st, TMS.make_approx_W2(PreMPO(st, op), τ))
    st2, st3 = State{Pure}(System(2, q), "Up"), State{Pure}(System(3, q), "Up")
    # a term of one site that a coupling goes through, commuting with it
    @test norm(w(st3, h * Z(2) + J * X(1) * X(3)) -
               kron(o, exp(τ * h * z), o) * (id8 + τ * J * kron(x, o, x))) < 1e-12
    # one at the end of a coupling, anticommuting with its factor there
    @test norm(w(st2, h * Z(1) + J * X(1) * X(2)) -
               (kron(exp(τ * h * z), o) + J * sinh(τ * h) / h * kron(x, x))) < 1e-12
    # two couplings meeting on a site by anticommuting factors, whose two orders cancel
    @test norm(w(st3, J * X(1) * X(2) + K * Y(2) * Y(3)) -
               (id8 + τ * (J * kron(x, x, o) + K * kron(o, y, y)))) < 1e-12
end

@testset "The error of WII comes from the terms that cross a same link" begin
    # WII keeps, up to order τ³, every product of terms of which no two cross the same link:
    # the coefficient of τ² in exp(τH) - WII is half the sum of H_x H_y over the ordered pairs
    # of terms crossing a link in common, a term with itself included. It is read off by a
    # Richardson extrapolation, which removes the order τ³. The atoms X, Y and Z leave the
    # terms on their sites when PreMPO compacts them
    rng = Xoshiro(20261001)
    n = 5
    st = State{Pure}(System(n, Qubit()), "Up")
    pauli() = rand(rng, (X, Y, Z))
    terms = [ [ (randn(rng) * pauli()(i), i:i) for i in 1:n ];
              [ (randn(rng) * pauli()(i) * pauli()(j), i:j) for i in 1:n for j in i+1:n ];
              [ (randn(rng) * pauli()(i) * pauli()(j) * pauli()(k), i:k)
                for (i, j, k) in ((1, 2, 4), (2, 3, 5)) ] ]
    h = sum(first, terms)
    ds = [ mpo_matrix(st, make_mpo(st, t)) for (t, _) in terms ]
    crossing(a, b) = max(first(a), first(b)) < min(last(a), last(b))
    c2 = sum(ds[x] * ds[y] for x in eachindex(terms), y in eachindex(terms)
             if crossing(last(terms[x]), last(terms[y]))) / 2
    hd = mpo_matrix(st, make_mpo(st, h))
    err(τ) = exp(τ * hd) - mpo_matrix(st, TMS.make_approx_W2(PreMPO(st, h), τ))
    τ = 1e-3
    @test norm((8 * err(τ / 2) - err(τ)) / τ^2 - c2) < 1e-3 * norm(c2)
end

@testset "WI and WII of a factor holding the identity" begin
    # N is (1 - Z) / 2 where Z also opens terms, and the identity taken out moves part of the
    # operator to terms of one site: the approximations remain of the first order, their error
    # over a step shrinking four times when the step is halved
    n = 4
    q = Qubit()
    st = State{Pure}(System(n, q), "Up")
    h = sum(0.7^(j - i) * (N(i) * N(j) + Z(i) * Z(j)) for i in 1:n for j in i+1:n) +
        0.4 * sum(X(i) for i in 1:n)
    on(a, i) = on_qubit(matrix(a, q), i, n)
    hd = sum(0.7^(j - i) * (on(N, i) * on(N, j) + on(Z, i) * on(Z, j))
             for i in 1:n for j in i+1:n) + 0.4 * sum(on(X, i) for i in 1:n)
    for w in (TMS.make_approx_W1, TMS.make_approx_W2)
        err(τ) = norm(mpo_matrix(st, w(PreMPO(st, -im * h), τ)) - exp(-im * τ * hd))
        @test 3.5 < err(1e-3) / err(5e-4) < 4.5
    end
end

@testset "The first step and the limits of each step" begin
    # `cutoff` and `maxdim` may be given one value per step: a schedule does step by step
    # exactly what the same limits do one step at a time, a constant schedule is the plain
    # value, and one shorter than the steps keeps its last value
    n = 6
    sys = System(n, Qubit())
    h = -sum(Z(i) * Z(i + 1) for i in 1:n - 1) - sum(X(i) for i in 1:n)
    st0 = State{Pure}(sys, "X+")
    sched = [2, 4, 4, 8]
    lim(m) = Limits(cutoff = 1e-14, maxdim = m)
    for algo in (Tdvp(), ApproxW(order = 1))
        chained = foldl((st, m) -> evolve(-im * h, 0.1, st; nsteps = 1, algo, limits = lim(m)),
                        sched;
                        init = st0)
        scheduled = evolve(-im * h, 0.4, st0; algo, nsteps = 4, limits = lim(sched))
        @test maxlinkdim(scheduled) == maxlinkdim(chained)
        @test expect1(scheduled, Z) ≈ expect1(chained, Z)
        flat = evolve(-im * h, 0.4, st0; algo, nsteps = 4, limits = lim(4))
        for m in [[4, 4, 4, 4], [4]]
            @test expect1(evolve(-im * h, 0.4, st0; algo, nsteps = 4, limits = lim(m)), Z) ≈
                  expect1(flat, Z)
        end
        # step k takes the k-th value: the evolution resumed at step 3 is the whole one
        half = evolve(-im * h, 0.2, st0; algo, nsteps = 2, limits = lim(sched))
        @test norm(evolve(-im * h, 0.2, half; algo, nsteps = 2, step = 3, limits = lim(sched)) -
                   scheduled) < 1e-12
    end
    # on the ring of "An expansion before the first step", whose first step from a product
    # state needs the expansion before it, tdvp expands a state of bond dimension below
    # min(5, maxdim) before any step, whatever its number, the periods coming on top of that
    ring = -im * (sum(Z(i) * Z(i + 1) for i in 1:5) + Z(6) * Z(1) - sum(X(i) for i in 1:6))
    x1(step; algo = Tdvp()) =
        real(expect(evolve(ring, 0.1, st0; nsteps = 1, algo, step, limits = lim(8)), X(1)))
    @test x1(2) ≈ x1(1)
    @test x1(4; algo = Tdvp(expand_period = 3)) ≈ x1(1)
    # a mixed state made hermitian after step 3, which an evolution resumed at step 3 keeps,
    # and one started again at step 1 does not
    ρ0 = mix(State{Pure}(System(4, Qubit()), "X+"))
    ev = -im * (-sum(Z(i) * Z(i + 1) for i in 1:3) - sum(X(i) for i in 1:4)) +
         Dissipator(sqrt(0.2) * Sm)(2)
    for algo in (Tdvp(expand_period = 3, hermitianize_period = 3),
                 ApproxW(order = 2, hermitianize_period = 3))
        whole = evolve(ev, 0.4, ρ0; algo, nsteps = 4, limits = lim(16))
        half = evolve(ev, 0.2, ρ0; algo, nsteps = 2, limits = lim(16))
        rest = evolve(ev, 0.2, half; algo, nsteps = 2, step = 3, limits = lim(16))
        @test norm(rest - whole) < 1e-12
    end
    algo = ApproxW(order = 2, hermitianize_period = 3)
    whole = evolve(ev, 0.4, ρ0; algo, nsteps = 4, limits = lim(16))
    half = evolve(ev, 0.2, ρ0; algo, nsteps = 2, limits = lim(16))
    @test norm(evolve(ev, 0.2, half; algo, nsteps = 2, limits = lim(16)) - whole) > 1e-8
    # and time: an evolver that depends on time resumed at the time it had reached
    terms = [-im * h, -im * sum(Z(i) for i in 1:n) => t -> cos(4t)]
    for algo in (Tdvp(), ApproxW(order = 4))
        whole = evolve(terms, 0.4, st0; algo, nsteps = 4, limits = lim(16))
        half = evolve(terms, 0.2, st0; algo, nsteps = 2, limits = lim(16))
        @test norm(evolve(terms, 0.2, half; algo, nsteps = 2, time = 0.2, step = 3,
                          limits = lim(16)) - whole) < 1e-12
    end
    # and thermalize, the log traces of the two parts adding up to that of the whole
    ρ = State{Mixed}(System(4, Qubit()), "FullyMixed")
    hq = -sum(Z(i) * Z(i + 1) for i in 1:3) - sum(X(i) for i in 1:4)
    for algo in (Tdvp(), ApproxW(order = 2))
        @test maxlinkdim(last(thermalize(hq, 0.1, ρ; nsteps = 1, algo, limits = lim([2, 8])))) == 2
        @test maxlinkdim(last(thermalize(hq, 0.1, ρ; nsteps = 1, algo, step = 2,
                                         limits = lim([2, 8])))) > 2
        l, whole = thermalize(hq, 0.4, ρ; algo, nsteps = 4, limits = lim(sched))
        l1, half = thermalize(hq, 0.2, ρ; algo, nsteps = 2, limits = lim(sched))
        l2, rest = thermalize(hq, 0.2, half; algo, nsteps = 2, step = 3, limits = lim(sched))
        @test l1 + l2 ≈ l
        @test norm(rest - whole) < 1e-12
    end
    # an Evolve phase takes the step dividing its duration into the nearest whole number of
    # steps, ends where it was asked to, and skips a duration shorter than half a step
    function steps(duration, time_step)
        io = IOBuffer()
        sim = runTMS(SimData(time_start = 1., phases = [
            CreateState{Pure}(2, Qubit(), "X+"),
            Evolve(evolver = -im * Z(1), duration = duration, time_step = time_step,
                   measurements = Data("t") => :step)]); output = io)
        return sim, [ last(first(r)) for r in get(sim.data, "t", []) ], String(take!(io))
    end
    sim, times, log = steps(1., 0.3)
    @test sim.time == 2.
    @test times ≈ [4/3, 5/3, 2.]
    @test occursin("Taking a step of 0.3333333333333333 rather than 0.3 to cover 1.0 in 3 steps",
                   log)
    sim, times, _ = steps(-1., 0.5)
    @test sim.time == 0. && times ≈ [0.5, 0.]
    @test real(expect(sim.state, X(1))) ≈ cos(2.)
    sim, times, log = steps(0.1, 0.3)
    @test sim.time == 1. && isempty(times)
    @test occursin("Skipping an evolution of 0.1, shorter than half a step", log)
    @test_throws "the time_step of an Evolve cannot be zero" Evolve(
        evolver = -im * Z(1), duration = 1., time_step = 0)
    # its number of steps in place of its time step, one of the two
    @test_throws "an Evolve needs its number of steps" Evolve(evolver = -im * Z(1), duration = 1.)
    @test_throws "an Evolve takes nsteps or time_step, not both" Evolve(
        evolver = -im * Z(1), duration = 1., time_step = 0.1, nsteps = 10)
    evolved(; kw...) = runTMS(SimData(phases = [CreateState{Pure}(2, Qubit(), "X+"),
        Evolve(; evolver = -im * Z(1), duration = 1., kw...)]); output = devnull).state
    @test norm(evolved(nsteps = 3) - evolved(time_step = 0.3)) < 1e-14
end

@testset "Periods below one mean never" begin
    # `mod(step, 0)` raised a division by zero in 1.x and `mod(step, -2)` is zero on every
    # second step, so anything below one is read as never rather than as one of those
    due = TMS.due
    for period in (0, -1, -2, -3)
        @test !any(k -> due(period, k), 1:12)
    end
    @test filter(k -> due(1, k), 1:4) == [1, 2, 3, 4]
    @test filter(k -> due(3, k), 1:10) == [3, 6, 9]
    # end to end: no hermitianization with a period of zero, the state as by default, up to
    # rounding, the links of each MPO taking indices of their own; a negative period is refused
    # when the method is written, see "The methods and their parameters"
    ρ = mix(State{Pure}(System(3, Qubit()), "X+"))
    ev = -im * (X(1) * X(2) + Z(2) * Z(3)) + Dissipator(Sm)(1)
    plain = evolve(ev, 0.3, ρ; algo = ApproxW(order = 2), nsteps = 3)
    @test 1 - hermiticity(plain) > 1e-8
    @test norm(evolve(ev, 0.3, ρ; algo = ApproxW(order = 2, hermitianize_period = 0),
                      nsteps = 3) - plain) < 1e-14
    # the measurements of an Evolve phase every measurements_period steps, from its start
    sim = runTMS(SimData(phases = [CreateState{Pure}(2, Qubit(), "X+"),
                                   Evolve(evolver = -im * Z(1), duration = 1., time_step = 0.1,
                                          measurements = Data("p") => :step,
                                          measurements_period = 3)]);
                 output = devnull)
    @test [ last(r[2]) for r in sim.data["p"] ] == [3, 6, 9]
end

@testset "Hermitianizing a mixed state" begin
    # from order 2, the complex fractions of the step of ApproxW keep a density matrix
    # hermitian only up to the order of the approximation, which hermitianize_period corrects
    ρ = mix(State{Pure}(System(3, Qubit()), "X+"))
    ev = -im * (X(1) * X(2) + Z(2) * Z(3) + 0.5 * Z(1)) + Dissipator(sqrt(0.5) * Sm)(1)
    deviation(algo) = 1 - hermiticity(evolve(ev, 1., ρ; algo, nsteps = 4))
    @test deviation(ApproxW(order = 2)) > 1e-8
    @test deviation(ApproxW(order = 2, hermitianize_period = 1)) < 1e-14
    @test deviation(Tdvp(hermitianize_period = 1)) < 1e-14
end

@testset "An expansion before the first step" begin
    # a ring of six qubits from a product state, with tdvp: the first step, of bond dimension
    # one, left the tangent space through Z(6)Z(1), and the expansion, coming after it in 1.x,
    # left an error of order the time step, 0.016 on X(1), whatever expand_period. The values
    # are exact, computed by test/reference/ising_ed.jl
    ring = -im * (sum(Z(i) * Z(i + 1) for i in 1:5) + Z(6) * Z(1) - sum(X(i) for i in 1:6))
    for expand_period in (1, 2)
        st = evolve(ring, 1., State{Pure}(System(6, Qubit()), "X+"); algo = Tdvp(; expand_period),
                    nsteps = 10, limits = Limits(maxdim = 8, cutoff = 1e-14))
        @test real(expect(st, X(1))) ≈ 0.48881258418 atol = 1e-8
        @test real(expect(st, Z(1) * Z(6))) ≈ -0.51118741582 atol = 1e-8
    end
end

@testset "The expansion before a step of tdvp" begin
    # tdvp expands a state whose bond dimension is below min(5, maxdim) before any step, the
    # default Tdvp() included: on the ring above, the same values as with expand_period. The
    # warnings tested are those of TMS, the deprecations Julia 1.14 brings to the QR of NDTensors,
    # which the expansion goes through, left aside
    ring = -im * (sum(Z(i) * Z(i + 1) for i in 1:5) + Z(6) * Z(1) - sum(X(i) for i in 1:6))
    st = @test_logs without_deprecations(
        () -> evolve(ring, 1., State{Pure}(System(6, Qubit()), "X+"); nsteps = 10,
                     limits = Limits(maxdim = 8, cutoff = 1e-14)))
    @test real(expect(st, X(1))) ≈ 0.48881258418 atol = 1e-8
    # from FullyMixed on a system conserving N, where the Krylov vectors of ITensorMPS, truncated
    # to one more than the bond dimension, kept nothing and expand failed
    h = sum(X(i) * X(i + 1) + Y(i) * Y(i + 1) + Z(i) * Z(i + 1) for i in 1:3)
    ρ = State{Mixed}(System(4, Qubit(conserve = N)), "FullyMixed")
    @test_logs without_deprecations(
        () -> thermalize(h, 0.1, ρ; nsteps = 1, limits = Limits(maxdim = 16)))
    # on a mixed product state under a term of small weight, expand changed the state by 2.4e-5,
    # the vectors it adds not quite orthogonal to its basis: the state written again in the
    # enlarged bonds is exact, and so is the step
    ρx = mix(State{Pure}(System(2, Qubit()), "X+"))
    m = make_mpo(ρx, -im * (Z(1) + 1e-3 * Z(2)))
    @test norm(TMS.expanded(ρx, m) - ρx) < 1e-14
    st = @test_logs without_deprecations(
        () -> evolve(-im * (Z(1) + 1e-3 * Z(2)), 0.01, ρx; nsteps = 1,
                     limits = Limits(maxdim = 10)))
    @test real(expect(st, X(2))) ≈ cos(2e-5) atol = 1e-14
    # and on a mixed state of fermions, where the enlarged bonds miss about 1e-13 of the state a
    # bond, within the precision of the expansion, kept: from FullyMixed in 10 steps, against the
    # same thermal state in 80
    tv = sum(dag(C)(i) * C(i + 1) + dag(C)(i + 1) * C(i) + 0.5 * N(i) * N(i + 1) for i in 1:5)
    f = State{Mixed}(System(6, Fermion(conserve = N)), "FullyMixed")
    _, ρ10 = @test_logs without_deprecations(
        () -> thermalize(tv, 2., f; nsteps = 10, limits = Limits(maxdim = 64)))
    _, ρ80 = thermalize(tv, 2., f; nsteps = 80, limits = Limits(maxdim = 64))
    @test real(expect(ρ10, tv)) ≈ real(expect(ρ80, tv)) atol = 1e-6
    # the expansion is dropped, with a warning, when expand fails on charged indices, as on a
    # state the evolver maps onto itself: the step is then exact all the same
    up = mix(State{Pure}(System(4, Qubit(conserve = N)), "Up"))
    st = @test_logs (:warn, r"failed on charged indices \(BoundsError\)") without_deprecations(
        () -> evolve(-im * h, 0.1, up; nsteps = 1, limits = Limits(maxdim = 16)))
    @test real(expect(st, Z(1))) ≈ 1
end

@testset "Ising chain" begin
    # The reference values are exact, computed by diagonalizing the 64 dimensional Hilbert space
    # in test/reference/ising_ed.jl, which shares no code with what is tested here. The
    # tolerance is therefore the error of the evolution below, about 2e-8, and not the
    # precision of the references. The last one is analytic: the product of all X commutes with
    # the Hamiltonian and starts at 1, so it stays at 1
    ring = -im * (sum(Z(i) * Z(i + 1) for i in 1:5) + Z(6) * Z(1) - sum(X(i) for i in 1:6))
    st = evolve(ring, 1., State{Pure}(System(6, Qubit()), "X+"); algo = ApproxW(order = 4, w = 2),
                nsteps = 50, limits = Limits(maxdim = 8))
    re(ops) = real(expect(st, ops))
    @test re([X(i) for i in 1:6]) ≈ fill(0.48881258418, 6) atol = 1e-7
    @test re([Y(i) for i in 1:6]) ≈ zeros(6) atol = 1e-7
    @test re([Z(i) for i in 1:6]) ≈ zeros(6) atol = 1e-7
    @test re([Z(1)Z(2), Z(2)Z(3), Z(1)Z(6)]) ≈ fill(-0.51118741582, 3) atol = 1e-7
    @test re([Y(1)Y(2), Y(2)Y(3), Y(1)Y(6)]) ≈ fill(-0.2518341076, 3) atol = 1e-7
    @test re([X(1)X(2), X(2)X(3), X(1)X(6)]) ≈ fill(0.1123126174, 3) atol = 1e-7
    fours = [X(1)X(2)X(3)X(4), X(2)X(3)X(4)X(5), X(4)X(5)X(6)X(1)]
    @test re(fours) ≈ fill(0.1123126174, 3) atol = 1e-7
    @test entanglement_entropy(st, 3) ≈ 1.15220908566 atol = 1e-7
    @test re(X(1)X(2)X(3)X(4)X(5)X(6)) ≈ 1 atol = 1e-8
end

@testset "Free fermions with source" begin
    # The reference values are exact, from the dense Lindblad evolution of the 32 dimensional
    # Fock space in test/reference/fermion_lindblad.jl, which shares no code with what is tested
    # here. The tolerance is the error of the evolution below, about 8e-8
    ρ = evolve(-im * sum(dag(C)(i) * C(i+1) + dag(C)(i+1) * C(i) for i in 1:4) +
               Dissipator(sqrt(2 * 0.2) * dag(C))(3), 1., State{Mixed}(System(5, Fermion()), "0");
               nsteps = 20, limits = Limits(maxdim = 16))
    @test real(expect1(ρ, N)) ≈ [0.0125582080327, 0.063008590052, 0.1950187333854, 0.063008590052,
                                 0.0125582080327] atol = 1e-6
    @test expect(ρ, [dag(C)(3) * C(i) for i in 1:5]) ≈
          [-0.023529279887, -0.0783682258151im, 0.1950187333854, -0.0783682258151im,
           -0.023529279887] atol = 1e-6
    @test purity(ρ) ≈ 0.5256648672193 atol = 1e-6
end

@testset "Free bosons with source" begin
    # The reference values are exact, from the Gaussian moments of this quadratic Lindbladian
    # in test/reference/boson_gaussian.jl, which shares no code with what is tested here. The
    # evolution below is 6e-7 to 8e-7 away from them. At a bond dimension of 10 rather than 16,
    # which singular values survive depends on the rounding, and the answer on the order the
    # sums are taken in, 1e-5 being missed on one BLAS thread count in 1.x
    ρ = evolve(-im * sum(A(i) * dag(A)(i+1) + dag(A)(i) * A(i+1) for i in 1:3) +
               Dissipator(2 * sqrt(0.1) * dag(A))(2), 0.3, State{Mixed}(System(4, Boson(7)), "0");
               algo = ApproxW(order = 4, w = 2), nsteps = 3, limits = Limits(maxdim = 16))
    @test real(expect1(ρ, N; complex = true)) ≈
          [0.0036328913502, 0.1200236303285, 0.0035680174576, 4.86613106e-5] atol = 1e-5
    @test expect(ρ, [dag(A)(2) * A(i) for i in 1:4]; complex = true) ≈
          [-0.0180013634952im, 0.1200236303285, -0.0178675798942im, -0.0017865872242] atol = 1e-5
    @test purity(ρ) ≈ 0.7965328475916 atol = 1e-5
end

@testset "The Krylov parameters reach tdvp" begin
    # a single Krylov vector, never rebuilt, cannot hold the exponential of a step: the
    # precession is then off, where the default parameters give it exactly
    st0 = State{Pure}(System(2, Qubit()), "X+")
    lim = Limits(maxdim = 10, cutoff = 1e-15)
    evolved(krylov) = evolve(-im * Z(1), 1., st0; algo = Tdvp(; krylov), nsteps = 10, limits = lim)
    precession_error(st) = abs(expect(st, X(1)) - cos(2.)) + abs(expect(st, Y(1)) - sin(2.))
    @test precession_error(evolved(Krylov())) < 1e-13
    # KrylovKit warns at each step that it has not converged
    _, poor = Test.collect_test_logs(() -> evolved(Krylov(dim = 1, maxiter = 1)))
    @test precession_error(poor) > 1e-2
    # and so does an Evolve phase, KrylovKit warning outside of the log of the run
    _, sim = Test.collect_test_logs(() -> runTMS(SimData(phases = [
        CreateState{Pure}(2, Qubit(), "X+"),
        Evolve(evolver = -im * Z(1), duration = 1., time_step = 0.1, limits = lim,
               algo = Tdvp(krylov = Krylov(dim = 1, maxiter = 1)))]); output = devnull))
    @test precession_error(sim.state) > 1e-2
end

@testset "Evolving a state that carries charges" begin
    # the links of the MPO of the generator carry the charge each of its channels has
    # accumulated, and the answer must not notice. The mixed case uses dephasing, whose jump
    # commutes with the charge, so that it is a strong symmetry as well as a weak one
    h = -sum(dag(C)(i) * C(i+1) + dag(C)(i+1) * C(i) for i in 1:3)
    lind = -im * h + sum(Dissipator(sqrt(0.3) * N)(i) for i in 1:4)
    lim = Limits(cutoff = 1e-14, maxdim = 32)
    function start(site, mixed)
        sys = System(4, site)
        p(v) = State{Pure}(sys, v)
        s = (p(["Occ", "Emp", "Occ", "Emp"]) + 0.5 * p(["Emp", "Occ", "Occ", "Emp"])) / sqrt(1.25)
        return mixed ? mix(s) : s
    end
    for (op, mixed) in ((-im * h, false), (lind, true)), algo in (Tdvp(), ApproxW(order = 2))
        numbers(site) = real.(expect1(evolve(op, 0.2, start(site, mixed); nsteps = 1, algo,
                                             limits = lim), N))
        ref = numbers(Fermion())
        for site in (Fermion(conserve = N), Fermion(conserve = strong(N)))
            @test numbers(site) ≈ ref
        end
    end
    # WI and WII share one channel between the term not begun and the term finished, so they
    # need a generator of zero flux. One that moved the charge would not keep the state in its
    # sector, which is why this is a refusal and not a gap
    charged = start(Fermion(conserve = N), false)
    @test_throws "need an operator of zero flux" evolve(dag(C)(1), 0.1, charged; nsteps = 1,
                                                        algo = ApproxW(order = 2), limits = lim)
end

@testset "The MPOs kept in the PreMPO" begin
    # each method keeps the last MPOs it built, for the parameters they were built for, and a
    # second evolution under the same PreMPO and step takes them again
    sys = System(3, Qubit())
    st = State{Pure}(sys, "+")
    h = sum(Z(i) * Z(i + 1) for i in 1:2) + 0.5 * sum(X(i) for i in 1:3)
    pre = PreMPO(st, -im * h)
    @test isempty(pre.cache)
    a = evolve(pre, 0.2, st; nsteps = 2)
    mpo = last(pre.cache[:mpo])
    b = evolve(pre, 0.2, a; nsteps = 2)
    @test last(pre.cache[:mpo]) === mpo
    @test norm(b - evolve(-im * h, 0.4, st; nsteps = 4)) < 1e-12
    # those of ApproxW for its order, its approximation and its step
    w2 = ApproxW(order = 2)
    a = evolve(pre, 0.2, st; algo = w2, nsteps = 2)
    key, ws = pre.cache[:approx_W]
    @test key == (2, 2, 0.1)
    @test length(ws) == 2
    b = evolve(pre, 0.4, a; algo = w2, nsteps = 4)
    @test last(pre.cache[:approx_W]) === ws
    @test norm(b - evolve(-im * h, 0.6, st; algo = w2, nsteps = 6)) < 1e-12
    for (algo, dt) in [(ApproxW(order = 2), 0.2), (ApproxW(order = 2, w = 1), 0.1),
                       (ApproxW(order = 3), 0.1)]
        evolve(pre, dt, st; nsteps = 1, algo)
        @test first(pre.cache[:approx_W]) == (algo.order, algo.w, dt)
        @test last(pre.cache[:approx_W]) !== ws
    end
    # the methods keep theirs apart
    @test last(pre.cache[:mpo]) === mpo
    # an evolver that depends on time keeps nothing, its MPOs changing at each step
    pv = PreMPO(st, [-im * h, -im * X(1) => cos])
    for algo in (Tdvp(), ApproxW(order = 4))
        evolve(pv, 0.2, st; algo, nsteps = 2)
    end
    @test isempty(pv.cache)
    # thermalize keeps them too, for a step of -dβ/2
    ρ = State{Mixed}(sys, "FullyMixed")
    ph = PreMPO(ρ, h)
    l1, ρ1 = thermalize(ph, 0.4, ρ; nsteps = 2)
    mpo = last(ph.cache[:mpo])
    l2, ρ2 = thermalize(ph, 0.4, ρ; nsteps = 2)
    @test last(ph.cache[:mpo]) === mpo
    @test l1 ≈ l2
    @test norm(ρ1 - ρ2) < 1e-14
    thermalize(ph, 0.4, ρ; algo = w2, nsteps = 2)
    @test first(ph.cache[:approx_W]) == (2, 2, -0.1)
end

@testset "The methods and their parameters" begin
    @test Tdvp() == Tdvp(Krylov(), 0, 0)
    @test repr(Tdvp()) == "Tdvp(krylov = Krylov(), expand_period = 0, hermitianize_period = 0)"
    @test repr(Tdvp(krylov = Krylov(tol = 1e-10), expand_period = 2)) ==
          "Tdvp(krylov = Krylov(tol = 1.0e-10), expand_period = 2, hermitianize_period = 0)"
    @test repr(ApproxW(order = 4)) ==
          "ApproxW(order = 4, w = 2, apply_algo = \"naive\", hermitianize_period = 0)"
    @test ApproxW(order = 3, w = 1, apply_algo = "densitymatrix").apply_algo == "densitymatrix"
    @test Tdvp() isa TMS.EvolutionAlgorithm && ApproxW(order = 1) isa TMS.EvolutionAlgorithm
    # a period is a number of steps, 0 for never, a negative one refused as for every period
    @test_throws "expand_period is a number of steps, 0 for never, not -1" Tdvp(expand_period = -1)
    @test_throws "hermitianize_period is a number of steps" Tdvp(hermitianize_period = -2)
    @test_throws "hermitianize_period is a number of steps" ApproxW(order = 2,
                                                                   hermitianize_period = -3)
    @test_throws "measurements_period is a number of steps" Evolve(
        evolver = -im * Z(1), duration = 1., nsteps = 1, measurements_period = -1)
    @test_throws "measurements_period is a number of steps" Update(identity;
                                                                  measurements_period = -1)
    # the order and the approximation are checked when the method is written, not at its first
    # step, and the order has no default
    @test_throws UndefKeywordError ApproxW()
    for order in (0, 5)
        @test_throws "ApproxW has orders 1 to 4, not $order" ApproxW(; order)
    end
    @test_throws "ApproxW has the approximations WI and WII, w = 1 or 2, not 3" ApproxW(
        order = 2, w = 3)
    for alg in ("fit", "zipup", "exact")
        @test_throws "apply_algo is" ApproxW(order = 2, apply_algo = alg)
    end
    # the fractions of a step of each order add up to the whole step
    for c in TMS.w_approx_coefs
        @test sum(c) ≈ 1
    end
end

@testset "A single site" begin
    # which tdvp of ITensorMPS does not take: the step exponentiates the site whole, exact for an
    # evolver that does not depend on time or commutes with itself at all times
    one = System(1, Qubit())
    ψ = State{Pure}(one, "X+")
    @test real(expect(evolve(-im * Z(1), 0.5, ψ; nsteps = 5), X(1))) ≈ cos(1)
    @test real(expect(evolve(-im * Z(1) => t -> 2t, 0.5, ψ; nsteps = 5), X(1))) ≈ cos(0.5)
    # a decaying qubit, from Up
    ρ = evolve(Dissipator(sqrt(0.5) * Sm)(1), 1.0, mix(State{Pure}(one, "Up")); nsteps = 4)
    @test real(expect(ρ, Z(1))) ≈ 2exp(-0.5) - 1
    # and thermalize, from the fully mixed state
    l, ρβ = thermalize(Z(1), 0.7, State{Mixed}(one, "FullyMixed"); nsteps = 1)
    @test l ≈ log(cosh(0.7))
    @test real(expect(ρβ, Z(1))) ≈ -tanh(0.7)
end

@testset "What evolve refuses" begin
    st = State{Pure}(System(2, Qubit()), "Up")
    ρ = mix(st)
    # an MPO, which has lost the terms the methods build their MPOs from
    msg = "the solvers take an operator or its PreMPO, which keeps the MPOs they build, not an MPO"
    for s in (st, ρ)
        @test_throws msg evolve(make_mpo(s, -im * X(1)), 0.1, s; nsteps = 1)
    end
    # a PreMPO of another system, or of the other representation, alone or with its functions;
    # a system of the same sites is the same system
    other = State{Pure}(System(3, Qubit()), "Up")
    @test_throws "the PreMPO was prepared on another system than that of the state" evolve(
        PreMPO(other, -im * X(1)), 0.1, st; nsteps = 1)
    @test_throws "the PreMPO was prepared on another system than that of the state" evolve(
        PreMPO(other, -im * X(1) => t -> 1.), 0.1, st; nsteps = 1)
    @test_throws "the PreMPO was prepared on a state of the other representation" evolve(
        PreMPO(ρ, -im * X(1)), 0.1, st; nsteps = 1)
    @test_throws "the PreMPO was prepared on a state of the other representation" evolve(
        PreMPO(st, -im * X(1)), 0.1, ρ; nsteps = 1)
    # with its functions, where it ran into ITensors, or gave a broken state by ApproxW
    for algo in (Tdvp(), ApproxW(order = 1))
        @test_throws "the PreMPO was prepared on a state of the other representation" evolve(
            PreMPO(st, -im * X(1) => t -> 1.), 0.1, ρ; nsteps = 1, algo)
        @test_throws "the PreMPO was prepared on a state of the other representation" evolve(
            PreMPO(ρ, -im * X(1) => t -> 1.), 0.1, st; nsteps = 1, algo)
    end
    same = PreMPO(State{Pure}(System(2, Qubit()), "Dn"), -im * X(1))
    @test real(expect(evolve(same, 0.4, st; nsteps = 4), Z(1))) ≈ cos(0.8)
    @test_throws "3 is not an evolver" evolve(3, 0.1, st; nsteps = 1)
    for nsteps in (0, -3), algo in (Tdvp(), ApproxW(order = 2))
        @test_throws "an evolution takes at least one step, and nsteps is $nsteps" evolve(
            -im * X(1), 0.5, st; algo, nsteps)
    end
    # the options of ITensorMPS are not passed on: the truncation goes through `limits`
    @test_throws MethodError evolve(-im * Z(1), 0.1, st; nsteps = 1, maxdim = 4)
    @test_throws MethodError evolve(-im * Z(1), 0.1, st; nsteps = 1, algo = ApproxW(order = 1),
                                    maxdim = 4)
    # the number of steps or the time step, one of the two: the step adjusted as that of an
    # Evolve, and no step for a duration shorter than half of it
    @test_throws "an evolution needs its number of steps, nsteps = 100 for instance" evolve(
        -im * X(1), 0.5, st)
    @test_throws "an evolution takes nsteps or time_step, not both" evolve(
        -im * X(1), 0.5, st; nsteps = 5, time_step = 0.1)
    @test_throws "the time_step of an evolution cannot be zero" evolve(
        -im * X(1), 0.5, st; time_step = 0)
    ev(; kw...) = evolve(-im * X(1), 0.5, st; kw...)
    @test norm(ev(time_step = 0.1) - ev(nsteps = 5)) < 1e-14
    @test norm(ev(time_step = 0.3) - ev(nsteps = 2)) < 1e-14
    @test evolve(-im * X(1), 0.04, st; time_step = 0.1) === st
end

@testset "Thermal states" begin
    lim = Limits(cutoff = 1e-14, maxdim = 64)
    tr = LA.tr
    x, y, z = [0. 1. ; 1. 0.], [0. -im ; im 0.], [1. 0. ; 0. -1.]

    # from the fully mixed state, the thermal state of a transverse field Ising chain against
    # the exponential of its dense Hamiltonian, the trace included
    n, β = 4, 1.3
    sys = System(n, Qubit())
    h = -sum(Z(i) * Z(i+1) for i in 1:n-1) - 0.7 * sum(X(i) for i in 1:n)
    ρ0 = State{Mixed}(sys, "FullyMixed")
    l, ρ = thermalize(h, β, ρ0; nsteps = 13, limits = lim)
    @test keys(thermalize(h, β, ρ0; nsteps = 13, limits = lim)) == (:log_trace, :state)
    dense = -sum(on_qubit(z, i, n) * on_qubit(z, i+1, n) for i in 1:n-1) -
            0.7 * sum(on_qubit(x, i, n) for i in 1:n)
    g = exp(-β * dense)
    @test real(expect(ρ, h)) ≈ tr(g * dense) / tr(g) atol = 1e-6
    xx = on_qubit(x, 1, n) * on_qubit(x, 3, n)
    @test real(expect(ρ, X(1) * X(3))) ≈ tr(g * xx) / tr(g) atol = 1e-7
    @test purity(ρ) ≈ tr(g * g) / tr(g)^2 atol = 1e-7
    @test l ≈ log(tr(g) / 2^n) atol = 1e-7
    @test trace(ρ) ≈ 1
    @test norm(ρ - dense_state(sys, g / tr(g))) < 1e-6
    # the Hamiltonian prepared on the mixed state
    lp, ρp = thermalize(PreMPO(ρ0, h), β, ρ0; nsteps = 13, limits = lim)
    @test lp ≈ l
    @test norm(ρp - ρ) < 1e-12
    # the log traces of successive calls add up, those of steps of the same size giving the
    # same state
    l1, ρ1 = thermalize(h, 0.5, ρ0; nsteps = 5, limits = lim)
    l2, ρ2 = thermalize(h, 0.8, ρ1; nsteps = 8, limits = lim)
    l13, ρ13 = thermalize(h, β, ρ0; nsteps = 13, limits = lim)
    @test l1 + l2 ≈ l13
    @test norm(ρ2 - ρ13) < 1e-10
    # a state of any trace is normalized first
    l3, ρ3 = thermalize(h, β, 3 * ρ0; nsteps = 13, limits = lim)
    @test l3 ≈ l
    @test norm(ρ3 - ρ) < 1e-12
    # on sites of different dimensions, log Z - Σ log d_i
    l, ρ = thermalize(Z(1) + Sz(2) + 0.5 * Z(3), 0.7,
                      State{Mixed}(System([Qubit(), Spin(1), Qubit()]), "FullyMixed"); nsteps = 4)
    @test l ≈ log(2cosh(0.7)) + log(1 + 2cosh(0.7)) + log(2cosh(0.35)) - 2log(2) - log(3)
    @test real(expect(ρ, Sz(2))) ≈ -2sinh(0.7) / (1 + 2cosh(0.7))

    # ApproxW, its orders converging in imaginary time as in real time; the order 4, two of
    # whose sub steps go backwards, converges as well when dβ times the width of the spectrum
    # is small, about 1 here
    for order in 1:4
        err = [ norm(last(thermalize(h, β, ρ0; nsteps, algo = ApproxW(; order),
                                     limits = Limits(cutoff = 0))) - dense_state(sys, g / tr(g)))
                for nsteps in (8, 16) ]
        @test abs(log2(err[1] / err[2]) - order) < 0.4
    end
    # the trace one in its real part: from order 2, the complex fractions of the step keep the
    # state hermitian only up to the order of the approximation, see "Hermitianizing a mixed
    # state"
    for (order, nsteps, tol, etol) in [(1, 64, 0.05, 0.01), (2, 32, 2e-3, 3e-4),
                                       (3, 16, 1e-5, 1e-4), (4, 16, 1e-6, 2e-6)]
        lw, ρw = thermalize(h, β, ρ0; nsteps, algo = ApproxW(; order), limits = lim)
        @test lw ≈ log(tr(g) / 2^n) atol = tol
        @test real(expect(ρw, h; complex = true)) ≈ tr(g * dense) / tr(g) atol = etol
        @test real(trace(ρw)) ≈ 1
        @test abs(imag(trace(ρw))) < 1e-7
    end

    # free fermions against Fermi-Dirac, <c†_i c_j> = [(e^{βh} + 1)^{-1}]_{ji}
    n, β = 8, 1.5
    hop = zeros(n, n)
    for i in 1:n-1
        hop[i, i+1] = hop[i+1, i] = -1.
    end
    l, ρ = thermalize(-sum(dag(C)(i) * C(i+1) + dag(C)(i+1) * C(i) for i in 1:n-1), β,
                      State{Mixed}(System(n, Fermion()), "FullyMixed"); nsteps = 15, limits = lim)
    ε = LA.eigvals(hop)
    @test expect2(ρ, (dag(C), C)) ≈ transpose(inv(exp(β * hop) + LA.I)) atol = 5e-5
    @test l ≈ sum(log.(1 .+ exp.(-β * ε))) - n * log(2) atol = 1e-5

    # a state that commutes with the Hamiltonian gives the thermal state restricted to it: one
    # electron of mixed spin per site, under a strong conservation of their number, makes the
    # t-J model the Heisenberg model, hopping blocked
    n, β, J = 4, 1.1, 0.8
    hopping = -sum(dag(Cup)(i) * Cup(i+1) + dag(Cup)(i+1) * Cup(i) +
                   dag(Cdn)(i) * Cdn(i+1) + dag(Cdn)(i+1) * Cdn(i) for i in 1:n-1)
    heisenberg = J * sum(Sz(i) * Sz(i+1) + (Sp(i) * Sm(i+1) + Sm(i) * Sp(i+1)) / 2 for i in 1:n-1)
    sys = System(n, Tj(conserve = strong(Ntot)))
    l, ρ = thermalize(hopping + heisenberg, β, State{Mixed}(sys, "MixedSpin"); nsteps = 11,
                      limits = lim)
    s = [x / 2, y / 2, z / 2]
    dense = J * sum(real(sum(on_qubit(m, i, n) * on_qubit(m, i+1, n) for m in s)) for i in 1:n-1)
    g = exp(-β * dense)
    @test real(expect(ρ, heisenberg)) ≈ tr(g * dense) / tr(g) atol = 1e-6
    @test real(expect(ρ, Sz(1) * Sz(2))) ≈
          tr(g * on_qubit(s[3], 1, n) * on_qubit(s[3], 2, n)) / tr(g) atol = 1e-7
    @test l ≈ log(tr(g) / 2^n) atol = 1e-7

    # the energy is the derivative of the logarithm of the trace, the two computed apart, by
    # successive calls of one step each
    n = 4
    h = -sum(Z(i) * Z(i+1) for i in 1:n-1) - 0.7 * sum(X(i) for i in 1:n)
    ρ = State{Mixed}(System(n, Qubit()), "FullyMixed")
    ls, es = [0.], Float64[]
    for k in 1:20
        lk, ρ = thermalize(h, 0.05, ρ; nsteps = 1, limits = lim)
        push!(ls, ls[end] + lk)
        push!(es, real(expect(ρ, h)))
    end
    @test maximum(abs((ls[k+2] - ls[k]) / 0.1 + es[k]) for k in 1:19) < 5e-3
    # which a Thermalize phase gives as :beta and :log_trace, the time left alone; it normalizes
    # at its first step only, where thermalize does at each call
    sim = runTMS(SimData(time_start = 1.5, phases = [
        CreateState{Mixed}(n, Qubit(), "FullyMixed"),
        Thermalize(hamiltonian = h, beta = 1., beta_step = 0.05, limits = lim,
                   measurements = Data("th") => [:beta, :log_trace])]); output = devnull)
    rows = sim.data["th"]
    @test [ last(r[2]) for r in rows ] ≈ 0.05:0.05:1.
    @test [ last(r[3]) for r in rows ] ≈ ls[2:end]
    @test all(last(r[1]) == 1.5 for r in rows) && sim.time == 1.5
    @test sim.log_trace ≈ ls[end] && sim.beta == 1.
    # shorter than half a step, it evolves nothing and still normalizes and gives the values, 0,
    # where it gave neither
    sim = runTMS(SimData(phases = [CreateState{Mixed}(n, Qubit(), "FullyMixed"),
                                   Update(c -> (c.state = 2 * c.state)),
                                   Thermalize(hamiltonian = h, beta = 0.02, beta_step = 0.05)]);
                 output = devnull)
    @test sim.beta == 0 && sim.log_trace == 0 && trace(sim.state) ≈ 1
end

@testset "What thermalize refuses" begin
    n = 4
    h = -sum(Z(i) * Z(i+1) for i in 1:n-1) - 0.7 * sum(X(i) for i in 1:n)
    # a pure state, which e^{-βH/2} takes to its ground state
    st = State{Pure}(System(n, Qubit()), "Up")
    for hp in (h, PreMPO(st, h))
        @test_throws "a thermal state needs a mixed representation: mix the state first" thermalize(
            hp, 1., st; nsteps = 1)
    end
    ρ = State{Mixed}(System(n, Qubit()), "FullyMixed")
    @test_throws "not an MPO" thermalize(make_mpo(ρ, h), 1., ρ; nsteps = 1)
    @test_throws "the PreMPO was prepared on a state of the other representation" thermalize(
        PreMPO(st, h), 1., ρ; nsteps = 1)
    @test_throws "the PreMPO was prepared on another system than that of the state" thermalize(
        PreMPO(State{Mixed}(System(3, Qubit()), "FullyMixed"), Z(1)), 1., ρ; nsteps = 1)
    @test_throws "thermalize takes at least one step, and nsteps is 0" thermalize(
        h, 1., ρ; nsteps = 0)
    # the number of steps or the step of beta, as for evolve
    @test_throws "thermalize needs its number of steps" thermalize(h, 1., ρ)
    @test_throws "thermalize takes nsteps or beta_step, not both" thermalize(
        h, 1., ρ; nsteps = 2, beta_step = 0.5)
    lb, ρb = thermalize(h, 1., ρ; beta_step = 0.5)
    ln, ρn = thermalize(h, 1., ρ; nsteps = 2)
    @test lb ≈ ln && norm(ρb - ρn) < 1e-14
    @test first(thermalize(h, 0.1, ρ; beta_step = 0.5)) == 0
    # a Hamiltonian depending on time, refused by a message that spoke of a search, and its PreMPO
    @test_throws "is not an operator or its PreMPO" thermalize([Z(1) => t -> 1.], 1., ρ;
                                                                nsteps = 1)
    @test_throws "the PreMPO depends on time, which an evolution alone follows" thermalize(
        PreMPO(ρ, Z(1) => t -> 1.), 1., ρ; nsteps = 1)
    @test_throws MethodError thermalize(h, 1., ρ; nsteps = 1, maxdim = 4)
    # and a Thermalize phase, a step of zero or none when it is written, a pure state when it
    # runs, its number of steps given in place of its step
    @test_throws "the beta_step of a Thermalize cannot be zero" Thermalize(
        hamiltonian = h, beta = 1., beta_step = 0)
    @test_throws "a Thermalize needs its number of steps" Thermalize(hamiltonian = h, beta = 1.)
    thermal(; kw...) = runTMS(SimData(phases = [CreateState{Mixed}(n, Qubit(), "FullyMixed"),
        Thermalize(; hamiltonian = h, beta = 1., kw...)]); output = devnull).state
    @test norm(thermal(nsteps = 2) - thermal(beta_step = 0.5)) < 1e-14
    @test_throws "a thermal state needs a mixed representation" runTMS(SimData(phases = [
        CreateState{Pure}(n, Qubit(), "Up"),
        Thermalize(hamiltonian = h, beta = 1., beta_step = 0.5)]); output = devnull)
end

@testset "An evolution under its MPO truncated" begin
    # mpo_limits truncates the MPO of the evolver, see PreMPO, in its triangular form, which
    # ApproxW reads: given an operator, a PreMPO truncated already, or a PreMPO to truncate, the
    # same evolution, near the exact one; and so for thermalize
    n = 8
    sys = System(n, Qubit())
    h = sum(Z(i) * Z(j) / (j - i)^2 for i in 1:n for j in i+1:n) + sum(0.5 * X(i) for i in 1:n)
    ψ = random_state(Pure, sys, 4; rng = Xoshiro(20261016))
    lim = Limits(cutoff = 1e-6)
    pre = PreMPO(ψ, -im * h; limits = lim)
    @test pre.linkdims .- 1 == [1, 2, 3, 3, 3, 2, 1]
    for algo in (ApproxW(order = 2), Tdvp())
        a = raw_vector(evolve(pre, 0.2, ψ; nsteps = 4, algo))
        @test norm(a - raw_vector(evolve(-im * h, 0.2, ψ; nsteps = 4, algo))) < 1e-4
        @test raw_vector(evolve(-im * h, 0.2, ψ; nsteps = 4, algo, mpo_limits = lim)) ≈ a
        @test raw_vector(evolve(PreMPO(ψ, -im * h), 0.2, ψ; nsteps = 4, algo, mpo_limits = lim)) ≈ a
    end
    ρ = State{Mixed}(sys, "FullyMixed")
    l0, ρ0 = thermalize(h, 0.5, ρ; nsteps = 5)
    l1, ρ1 = thermalize(h, 0.5, ρ; nsteps = 5, mpo_limits = lim)
    @test abs(l1 - l0) < 1e-4
end
