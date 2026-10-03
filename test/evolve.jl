# Time evolution.
#
# Goes here: the evolution algorithms themselves (Tdvp, ApproxW) on cases whose exact
# solution is known, including dissipative ones. A physical study that happens to use
# an evolution belongs to algorithms.jl instead.

@testset "Simple evolve" begin
    # the tolerances say how accurate each algorithm is on this problem, so they are as
    # tight as the method allows. They must not be set at the floating point floor though:
    # the two that were at 1e-14 failed on Julia 1.10 by 1.3e-14, pure rounding noise that
    # differs with the BLAS build. An error of the algorithm would be orders of magnitude
    # larger than these bounds, so 1e-13 loses nothing
    for (algo, time_step, tol) in [
        (Tdvp(), 0.1, 1e-13),
        (ApproxW(order=1, w=1), 0.01, 0.03),
        (ApproxW(order=4, w=1), 0.01, 1e-10),
        (ApproxW(order=1, w=2), 0.01, 1e-13),
        (ApproxW(order=4, w=2), 0.01, 1e-13),
        ]
        @test_pm test_phases([
            CreateState{type}(2, Qubit(), "X+"),
            Evolve(; time_step, algo,
            limits = Limits(maxdim = 10, cutoff = 1e-15),
            evolver = -im * Z(1),
            duration = 1,
            final_measures = check([X(1), Y(1)], t->[cos(2t), sin(2t)], tol)
            )
            ])
        end
end

@testset "Multi evolve" begin
    for (algo, time_step, tol) in [
        (Tdvp(), 0.1, 1e-14),
        (ApproxW(order=1, w=1), 0.01, 0.03),
        (ApproxW(order=4, w=1), 0.01, 1e-10),
        (ApproxW(order=1, w=2), 0.01, 0.03),
        (ApproxW(order=4, w=2), 0.01, 1e-10),       
        ]
        @test_pm test_phases([
            CreateState{type}(2, Qubit(), ["X+", "Z-"]),
            Evolve(; time_step, algo,
            limits = Limits(maxdim = 10, cutoff = 1e-15),
            evolver = -im * Z(1) * Z(2),
            duration = 1,
            final_measures = check([X(1), Y(1), Z(2)], t->[cos(2t), -sin(2t), -1], tol)
            )
        ])
    end
end

@testset "Complex evolve" begin
    for (algo, time_step, tol) in [
        (Tdvp(), 0.1, 1e-14),
        (ApproxW(order=1, w=1), 0.01, 0.04),
        (ApproxW(order=4, w=1), 0.01, 3e-9),
        (ApproxW(order=1, w=2), 0.01, 0.04),
        (ApproxW(order=4, w=2), 0.01, 2e-9),       
    ]
        @test_pm test_phases([
            CreateState{type}(5, Qubit(), ["X+", "Z+", "Z-", "X+", "X+"]),
            Evolve(; algo, time_step,
            limits = Limits(maxdim = 10, cutoff = 1e-15),
            evolver = -im * (Z(2)Z(4) + Z(3)Z(1) + Z(3)Z(5)),
            duration = 1,
            final_measures = check([Y(1), Y(4), Y(5)], t->sin(2t) * [-1., 1., -1.], tol)
            )
            ])
        end
    end

@testset "Time dependent evolve" begin
    # h(t) = Z(1) + t * Z(2): two commuting one site terms carrying different time
    # functions, so each spin precesses at a rate that is known exactly and a coefficient
    # matched with the wrong term shows up at once. The midpoint rule the solvers use to
    # sample the time functions is exact on both of them, so what is left is the error of
    # the evolution algorithm itself. The tolerances below come from the error measured
    # in this exact configuration: at machine precision everywhere except for ApproxW on a
    # mixed state, where rebuilding the MPO at every step and truncating at the cutoff a
    # hundred times leaves about 1e-13.
    hs = [-im * Z(1), -im * Z(2)]
    coefs = [t -> 1.0, t -> t]
    for (algo, time_step, tol) in [
        (Tdvp(), 0.01, 1e-13),
        (ApproxW(order=4, w=2), 0.01, 1e-12),
        ]
        @test_pm test_phases([
            CreateState{type}(2, Qubit(), "X+"),
            Evolve(; time_step, algo,
            limits = Limits(maxdim = 10, cutoff = 1e-15),
            evolver = hs => coefs,
            duration = 1,
            final_measures = check([X(1), Y(1), X(2), Y(2)],
                t->[cos(2t), sin(2t), cos(t^2), sin(t^2)], tol)
            )
            ])
    end
    # the low level entry point documented in others.md, on a mixed state so that lifting
    # a pure evolver to the mixed representation is exercised on every term of the vector
    st = tdvp(hs, 1.0, mix(State{Pure}(System(2, Qubit()), "X+")); coefs, nsweeps = 100)
    @test expect(st, X(1)) ≈ cos(2.0) atol = 1e-13
    @test expect(st, X(2)) ≈ cos(1.0) atol = 1e-13
end

@testset "One time function per term" begin
    # too few raised a BoundsError on an internal vector, and too many were ignored
    hs = [-im * Z(1), -im * Z(2)]
    st = State{Pure}(System(2, Qubit()), "X+")
    @test_throws "takes as many time functions, got 1" tdvp(hs, 0.1, st; coefs = [t -> 1.0])
    @test_throws "takes as many time functions, got 3" tdvp(hs, 0.1, st;
        coefs = [t -> 1.0, t -> 1.0, t -> 1.0])
    # in an Evolve, refused when it is written rather than at its first step, and a single
    # term with its function, written without the vectors, which failed on a MethodError
    evolve(evolver) = Evolve(; duration = 0.2, time_step = 0.1, algo = Tdvp(), evolver)
    @test_throws "one function of time per term" evolve(hs => [t -> 1.0])
    @test_throws "one function of time per term" evolve(-im * Z(1) => [t -> 1.0])
    @test evolve(-im * Z(1) => (t -> 1.0)).evolver isa Pair{<:Vector, <:Vector}
    sim = runTMS(SimData(phases = [CreateState{Pure}(2, Qubit(), "X+"),
                                   evolve(-im * Z(1) => (t -> 2.0))]); output = devnull)
    @test real(expect(sim.state, X(1))) ≈ cos(0.8)
end

@testset "Time functions take real values" begin
    # a complex value multiplied rho A† by itself rather than by its conjugate, which gave a
    # mixed state of complex trace; it is refused on pure states as well
    hs = -im * [X(1), -Y(1), Z(1) * Z(2) + X(2)]
    ψ = State{Pure}(System(2, Qubit()), ["Up", "X+"])
    for st in (ψ, mix(ψ))
        @test_throws "time functions take real values" tdvp(hs, 1.0, st;
            coefs = [t -> cis(3t), t -> 0.0, t -> 1.0])
    end
    # the same drive written with its real and imaginary parts: evolving the mixed state is
    # mixing the evolved pure state, the hamiltonian not commuting at different times. Tdvp
    # is exact at full bond dimension, the tolerance of ApproxW is its error measured here
    coefs = [t -> cos(3t), t -> sin(3t), t -> 1.0]
    for (ev, tol) in [(st -> tdvp(hs, 1.0, st; coefs, nsweeps = 50), 1e-11),
                      (st -> approx_W(hs, 1.0, st; coefs, nsweeps = 50, order = 4), 1e-7)]
        @test norm(ev(mix(ψ)) - mix(ev(ψ))) < tol
    end
end

@testset "An evolution ends where it was asked to" begin
    # the step is adjusted to divide the duration: a duration of 1 in steps of 0.3 stopped at
    # 0.9, and one shorter than half a step ran no step at all
    function evolve(duration, time_step; algo = Tdvp())
        return runTMS(SimData(phases = [
                CreateState{Pure}(2, Qubit(), "Up"),
                Evolve(; duration, time_step, algo, evolver = -im * X(1),
                       measures = Data("m") => Z(1))]);
            output = devnull)
    end
    sim = evolve(1.0, 0.3)
    @test sim.time ≈ 1.0
    @test sim.data["m"]["Z(1)"]["times"] ≈ [1/3, 2/3, 1]
    @test expect(sim.state, Z(1)) ≈ cos(2.0) atol = 1e-12
    sim = evolve(0.1, 0.3)
    @test sim.time == 0
    @test !haskey(sim.data, "m")
    # the duration gives the direction: a step of the other sign made no step, while the time
    # went on to the end of the duration
    for algo in (Tdvp(), ApproxW(order = 4)), duration in (0.5, -0.5), time_step in (0.1, -0.1)
        sim = evolve(duration, time_step; algo)
        @test sim.time ≈ duration
        @test length(sim.data["m"]["Z(1)"]["times"]) == 5
        @test expect(sim.state, Z(1)) ≈ cos(2duration) atol = 1e-10
    end
    # and the solvers take at least one step, where none advanced the time of a simulation
    st = State{Pure}(System(2, Qubit()), "Up")
    for nsweeps in (0, -3)
        @test_throws "takes at least one step" tdvp(-im * X(1), 0.5, st; nsweeps)
        @test_throws "takes at least one step" approx_W(-im * X(1), 0.5, st; nsweeps, order = 2)
    end
end

@testset "Positions checked before a gate on a mixed state" begin
    # the operator is checked as it was written, not as the factors its string is prepared
    # into, which named Gate(F)(5) for C(9)
    st = mix(State{Pure}(System(4, Fermion()), "Emp"))
    @test_throws "C(9) acts on site 9, which the system does not have" apply(C(9), st)
end

@testset "Time dependent terms of several sites" begin
    # the coefficient of a term goes once into the MPO, whatever the number of sites it spans:
    # laid on each of its pieces, a term of k sites took it to the power k, which the one site
    # terms above cannot show. Folding the coefficients into the terms gives the MPO the time
    # functions must give
    for R in (Pure, Mixed)
        ψ = RandomState{Pure}(System(4, Qubit()), 4)
        st = R == Pure ? ψ : mix(ψ)
        hs = [-im * X(1) * Z(2), -im * Z(1) * Y(2) * X(4), -im * Z(3)]
        c = [0.5, 0.3, 0.7]
        pre = PreMPO(st, hs)
        folded = PreMPO(st, sum(c .* hs))
        @test norm(prod(make_mpo(pre, c)) - prod(make_mpo(folded))) < 1e-12
        @test norm(prod(make_approx_W1(pre, 0.1, c)) - prod(make_approx_W1(folded, 0.1))) < 1e-12
        @test norm(prod(make_approx_W2(pre, 0.1, c)) - prod(make_approx_W2(folded, 0.1))) < 1e-12
    end
    # and through an evolution: -i X(1)X(2)/2 over a time 1 turns Z(1) by an angle 1
    st = tdvp([-im * X(1) * X(2)], 1.0, State{Pure}(System(2, Qubit()), "Up");
              coefs = [t -> 0.5], nsweeps = 10, limits = Limits(maxdim = 4, cutoff = 1e-14))
    @test expect(st, Z(1)) ≈ cos(1.0) atol = 1e-12
end

@testset "A constant term in the approximation WII" begin
    # a constant term is laid on the first site as a delta, whose diagonal storage ITensors
    # could not add a tensor of another element type to, nor, on a charged system, add
    # anything to or exponentiate: WII failed on each of these operators. A constant commutes
    # with everything, so all it does is multiply WII by its exponential
    for sys in (System(2, Qubit()), System(2, Qubit(conserve = N)))
        st = State{Pure}(sys, "Dn")
        w(op, coefs...) = prod(make_approx_W2(PreMPO(st, op), 0.1, coefs...))
        @test norm(w(-im * N(1) + 0.5 * Id(1)) - exp(0.05) * w(-im * N(1))) < 1e-12
        @test norm(w(-im * (N(2) + 0.5 * Id(1))) - exp(-0.05im) * w(-im * N(2))) < 1e-12
        @test norm(w([Id(1), -im * N(1)], [0.7, 0.3]) - exp(0.07) * w([-im * N(1)], [0.3])) < 1e-12
    end
end

@testset "WII keeps the products of terms that cross no same link" begin
    # WII is the approximation of Zaletel et al., whose blocks put the exponential of the terms
    # of one site around every piece, transports included, and take a closing and an opening
    # on the same site in both orders. On these operators it has a closed form, which the
    # former blocks missed by a term of order τ²
    q = Qubit()
    x, y, z, o = (matrix(a, q) for a in (X, Y, Z, Id))
    id8 = kron(o, o, o)
    τ, h, J, K = 0.1, 0.7, 0.4, 0.3
    w(st, op) = dense(st, make_approx_W2(st, op, τ))
    st2, st3 = State{Pure}(System(2, q), "Up"), State{Pure}(System(3, q), "Up")
    # a term of one site that a coupling goes through, commuting with it
    @test norm(w(st3, h * Z(2) + J * X(1) * X(3)) - kron(o, exp(τ * h * z), o) * (id8 + τ * J * kron(x, o, x))) < 1e-12
    # one at the end of a coupling, anticommuting with its factor there
    @test norm(w(st2, h * Z(1) + J * X(1) * X(2)) - (kron(exp(τ * h * z), o) + J * sinh(τ * h) / h * kron(x, x))) < 1e-12
    # two couplings meeting on a site by anticommuting factors, whose two orders cancel
    @test norm(w(st3, J * X(1) * X(2) + K * Y(2) * Y(3)) - (id8 + τ * (J * kron(x, x, o) + K * kron(o, y, y)))) < 1e-12
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
              [ (randn(rng) * pauli()(i) * pauli()(j) * pauli()(k), i:k) for (i, j, k) in ((1, 2, 4), (2, 3, 5)) ] ]
    h = sum(first, terms)
    ds = [ dense(st, make_mpo(st, t)) for (t, _) in terms ]
    crossing(a, b) = max(first(a), first(b)) < min(last(a), last(b))
    c2 = sum(ds[x] * ds[y] for x in eachindex(terms), y in eachindex(terms)
             if crossing(last(terms[x]), last(terms[y]))) / 2
    hd = dense(st, make_mpo(st, h))
    err(τ) = exp(τ * hd) - dense(st, make_approx_W2(st, h, τ))
    τ = 1e-3
    @test norm((8 * err(τ / 2) - err(τ)) / τ^2 - c2) < 1e-3 * norm(c2)
end

@testset "Per sweep limits in an evolution" begin
    # `cutoff` and `maxdim` may be given one value per sweep. dmrg is handed the whole
    # schedule, but the evolution solvers drive their sweeps themselves and have to pick
    # the value out for each one, so what is checked is that a schedule does sweep by
    # sweep exactly what the same limits do one sweep at a time
    n = 6
    sys = System(n, Qubit())
    h = sum(-Z(i) * Z(i + 1) for i in 1:n - 1) - sum(1. * X(i) for i in 1:n)
    st0 = State{Pure}(sys, "X+")
    sched = [2, 4, 4, 8]
    # `approx_W` takes the order of the approximation, with no default, the same way the
    # `ApproxW` phase does; `tdvp` has nothing of the sort, hence the per solver arguments
    for (solver, opts) in [(tdvp, (;)), (approx_W, (; order = 1))]
        chained = foldl(sched; init = st0) do st, m
            solver(-im * h, 0.1, st; nsweeps = 1, limits = Limits(cutoff = 1e-14, maxdim = m),
                   opts...)
        end
        scheduled = solver(-im * h, 0.4, st0;
                           nsweeps = 4, limits = Limits(cutoff = 1e-14, maxdim = sched), opts...)
        @test maxlinkdim(scheduled) == maxlinkdim(chained)
        @test expect1(scheduled, Z) ≈ expect1(chained, Z)
        # a constant schedule is the plain value, and one shorter than the sweeps keeps
        # its last value for the rest of them
        flat = solver(-im * h, 0.4, st0; nsweeps = 4, limits = Limits(cutoff = 1e-14, maxdim = 4),
                      opts...)
        for m in [[4, 4, 4, 4], [4]]
            st = solver(-im * h, 0.4, st0;
                        nsweeps = 4, limits = Limits(cutoff = 1e-14, maxdim = m), opts...)
            @test expect1(st, Z) ≈ expect1(flat, Z)
        end
    end
end

@testset "The Krylov parameters reach tdvp" begin
    # a single Krylov vector, never rebuilt, cannot hold the exponential of a step: the
    # precession is then off, where the default parameters give it exactly
    st0 = State{Pure}(System(2, Qubit()), "X+")
    lim = Limits(maxdim = 10, cutoff = 1e-15)
    poor = Krylov(dim = 1, maxiter = 1)
    evolved(krylov) = tdvp(-im * Z(1), 1., st0; nsweeps = 10, limits = lim, krylov)
    precession_error(st) = abs(expect(st, X(1)) - cos(2.)) + abs(expect(st, Y(1)) - sin(2.))
    @test precession_error(evolved(Krylov())) < 1e-13
    @test precession_error(evolved(poor)) > 1e-2
    # and a Tdvp phase hands them on
    sim = runTMS(SimData(phases = [
            CreateState(type = Pure(), state = st0),
            Evolve(duration = 1., time_step = 0.1, algo = Tdvp(krylov = poor), limits = lim,
                   evolver = -im * Z(1))]);
        output = devnull)
    @test precession_error(sim.state) ≈ precession_error(evolved(poor))
end

@testset "The algorithms of the product of an MPO by a state" begin
    # without truncation, the two give the exact product
    n = 4
    h = sum(X(i) * X(i + 1) for i in 1:n - 1)
    up = State{Pure}(System(n, Qubit()), "Up")
    st = tdvp(-im * h, 0.5, up; nsweeps = 5, limits = Limits(maxdim = 16, cutoff = 1e-14))
    m = make_mpo(st, h)
    w = approx_W(-im * h, 0.5, up; order = 2, nsweeps = 5)
    @test norm(apply(m, st; apply_algo = "naive") - apply(m, st)) < 1e-12
    @test norm(approx_W(-im * h, 0.5, up; order = 2, nsweeps = 5, apply_algo = "naive") - w) < 1e-12
    # "fit" needs a number of sweeps of its own, "zipup" is missing for a state from the lowest
    # ITensorMPS accepted, and a name that is none of them is refused too
    for alg in ("fit", "zipup", "exact")
        @test_throws "apply_algo is" apply(m, st; apply_algo = alg)
        @test_throws "apply_algo is" approx_W(-im * h, 0.5, up; order = 2, apply_algo = alg)
    end
    # and an ApproxW phase hands it on
    @test_throws "apply_algo is" runTMS(SimData(phases = [
            CreateState(type = Pure(), state = up),
            Evolve(duration = 0.2, time_step = 0.1, algo = ApproxW(order = 2, apply_algo = "fit"),
                   evolver = -im * h)]);
        output = devnull)
end

@testset "tdvp and approx_W refuse the options of ITensorMPS" begin
    # they no longer pass on what they do not know: the truncation goes through `limits`
    st = State{Pure}(System(2, Qubit()), "Up")
    @test_throws MethodError tdvp(-im * Z(1), 0.1, st; maxdim = 4)
    @test_throws MethodError approx_W(-im * Z(1), 0.1, st; order = 1, maxdim = 4)
end

@testset "Noisy gates" begin
    @test_ok test_phases([
        CreateState{Mixed}(1, Qubit(),"Up"),
        Gates(
            gates = (0.5*Gate(Id)+0.5*Gate(X))(1),
            final_measures = [
                check([X(1), Y(1), Z(1)], [0, 0, 0], 1e-12),
                check(Trace2, 0.5, 1e-12)
            ]
        )
    ])
    @test_ok test_phases([
        CreateState{Mixed}(2, Qubit(),"Up"),
        Gates(
            gates = (0.5*Gate(Id⊗Id)+0.5*Gate(X⊗X))(1, 2),
            final_measures = [
                check([X, Y, Z], [[0, 0], [0, 0], [0, 0]], 1e-12),
                check(Trace2, 0.5, 1e-12),
                check(Z(1)Z(2), 1, 1e-12)
            ]
        )
    ])
    # a power of a noisy gate goes into an MPO as well, which simplifies it where applying it
    # does not, and the two must agree
    ρ = mix(State{Pure}(System(2, Qubit()), ["Up", "+"]))
    for a in [(Gate(X)^2)(1), (Left(X)^2)(2), ((0.9 * Gate(Id) + 0.1 * Gate(X))^3)(1),
              (Gate(H)^0.5)(2)]
        @test norm(apply(a, ρ) - apply(make_mpo(ρ, a), ρ)) < 1e-12
    end
end

@testset "A product of gates is the operator it denotes" begin
    # its rightmost factor acts first, as in make_mpo and expect: ITensorMPS applies a list of
    # tensors first to last, which is the product read backwards. On a density matrix a gate
    # g acts as g ρ g†, what mixing the pure result gives
    q = State{Pure}(System(2, Qubit()), ["Up", "X+"])
    qf = State{Pure}(System([Qubit(), Fermion()]), ["X+", "Emp"])
    for (st, g) in [(q, Swap(1, 2) * X(1)), (q, X(1) * Swap(1, 2)), (q, Sp(2) * Sm(2)),
                    (qf, Sp(1) * Sm(1) * dag(C)(2)), (qf, dag(C)(2) * Sm(1) * Sp(1))]
        @test norm(apply(g, st) - apply(make_mpo(st, g), st)) < 1e-12
        @test norm(apply(g, mix(st)) - mix(apply(g, st))) < 1e-12
    end
end

@testset "Fermionic gates" begin
    # apply places one local tensor per factor and has no way to build a Jordan-Wigner
    # string, so an operator that still needs one is simplified first, which is what
    # inserts it. Only a fermionic operator is simplified, so a gate defined by an
    # expression, such as Swap, is never expanded into a sum apply could not place.
    # The two branches below differ in the parity of sites 1 and 2, which is what makes
    # the string acting further right observable instead of a global phase.
    sys = System(4, Fermion())
    st = normalize(State{Pure}(sys, ["1", "0", "1", "0"]) +
                   State{Pure}(sys, ["0", "0", "1", "0"]))
    for a in [C(3), dag(C)(4), dag(C)(4) * C(3)]
        @test norm(apply(a, st) - apply(make_mpo(st, a), st)) < 1e-12
    end
    # the mixed representation takes the same path: the one site factors removeMulti
    # leaves behind are built by the Multi_F constructor, which is where the knowledge of
    # which side of the density matrix the string acts on already lives. Mixing the result
    # of the pure application must give what applying it to the mixed state gives.
    for a in [C(3), dag(C)(4) * C(3)]
        @test norm(mix(apply(a, st)) - apply(a, mix(st))) < 1e-12
    end
    # a fermionic operator inside a superoperator or a tensor product takes its string as
    # well, (C ⊗ dag(C))(3, 4) being C(3) * dag(C)(4) by definition. Both used to be built
    # from the bare matrices
    ρ = mix(st)
    for a in [Gate(C)(3), Left(C)(3), Right(C)(3), Gate(C ⊗ dag(C))(3, 4)]
        @test norm(apply(a, ρ) - apply(make_mpo(ρ, a), ρ)) < 1e-12
    end
    for a in [(dag(C) ⊗ Id ⊗ C)(1, 2, 3), (C ⊗ dag(C))(3, 4)]
        @test norm(apply(a, st) - apply(make_mpo(st, a), st)) < 1e-12
    end
    # the terms of an odd sum of one site share their string, so that it is a single product
    # and can be applied, however it is written and wherever it sits
    for a in [(C + dag(C))(3), C(3) + dag(C)(3), (2C + im * dag(C))(4), C(3) * (C + dag(C))(2)]
        @test norm(apply(a, st) - apply(make_mpo(st, a), st)) < 1e-12
    end
    @test norm(apply(Gate(C + dag(C))(3), ρ) - apply(make_mpo(ρ, Gate(C + dag(C))(3)), ρ)) < 1e-12
    # what has no gate to become is refused: a dissipator of a fermionic operator turns into
    # a sum, and a function of an odd operator of several sites mixes the two parities. One of
    # an even operator is a gate, checked against dense matrices in observables.jl
    @test_throws "makes it a sum" apply(Dissipator(C)(3), ρ)
    @test_throws "mixes the two parities" apply(exp(0.3 * (C ⊗ Id + Id ⊗ C))(1, 3), st)
    # on the first site as well, where it has no string: it is refused for its parity
    @test_throws "a function of an odd fermionic operator, which mixes the two parities" apply(
        exp(0.3 * (C + dag(C)))(1), st)
    @test_ok apply(exp(-0.3im * (dag(C) ⊗ C + C ⊗ dag(C)))(3, 4), st)
    # a gate of several sites defined by a matrix is kept whole, and a string covering only
    # some of its sites does not commute with it. It used to be moved past it, which changed
    # the sign of the branch where the gate moves a fermion onto site 2
    sw = Operator{2}("Sw", [1. 0. 0. 0. ; 0. 0. 1. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 1.], involution_op)
    a = dag(C)(1) * dag(C)(4) * sw(1, 2)
    @test norm(apply(a, st) - apply(dag(C)(1), apply(dag(C)(4), apply(sw(1, 2), st)))) < 1e-12
    @test norm(mix(apply(a, st)) - apply(a, ρ)) < 1e-12
    # a term of coefficient zero leaves the parity of a gate alone, where it made one a sum of
    # fermionic and non fermionic operators, and a gate whose terms all vanish makes the state
    # null, as a gate that annihilates it does, where it was refused
    g = 0.0
    @test norm(apply((g * C + dag(C))(4), st) - apply(dag(C)(4), st)) < 1e-12
    @test norm(apply((g * C + g * dag(C))(4), st)) == 0
    @test trace(apply(0 * X(1), mix(State{Pure}(System(2, Qubit()), "Up")))) == 0
    # a factor contributing no tensor, an identity or a Jordan-Wigner string, must not
    # leave the gate list untyped: ITensorMPS.product has no method for a Vector{Any}
    @test_ok apply(Id(1) * X(2), State{Pure}(System(2, Qubit()), "Up"))
    # a non fermionic gate must not be simplified: Swap is defined by an expression and
    # simplifying a product of them would make a sum, which apply cannot place
    @test_ok apply(Swap(1, 2) * Swap(3, 4), State{Pure}(System(4, Qubit()), "Up"))
    # nor is one beside a fermionic factor, which used to make the whole gate a sum, refused
    sq = State{Pure}(System([Qubit(), Qubit(), Fermion()]), ["Up", "Dn", "0"])
    @test norm(apply(Swap(1, 2) * dag(C)(3), sq) - apply(Swap(1, 2), apply(dag(C)(3), sq))) < 1e-12
end

@testset "Gates of sums and of mixed parities" begin
    # Gate takes an operator not placed on sites. The gate of a placed sum K, which a product
    # of K by an operator on mixed states takes, is Left(K) Right(K), what K does to the pure
    # state, strings included, and on one site it is the gate of the operator of that site
    build_gate = TensorMixedStates.build_gate
    @test_throws "which is placed on sites" Gate(X(1) + Z(1))
    @test_throws "which is placed on sites" Gate(X(1))
    @test build_gate(X(1) + Z(1)) == Gate(X + Z)(1)
    for (sys, K, obs) in [(System(3, Qubit()), (X(1) + Z(2)) / sqrt(2), [X(1), Z(2), X(1) * Y(3)]),
                          (System(3, Fermion()), (C(1) + dag(C)(3)) / sqrt(2),
                           [N(1), N(3), dag(C)(1) * C(3)])]
        ψ = RandomState{Pure}(sys, 2)
        kψ = apply(make_mpo(ψ, K), ψ)
        ρk = apply(make_mpo(mix(ψ), build_gate(K)), mix(ψ))
        @test [expect(ρk, o) for o in obs] ≈ [expect(kψ, o) for o in obs]
    end
    # a factor of no definite parity holds an odd part: on site 1 it needs no string and is
    # applied, further on it becomes a sum, which is refused for that reason
    sys = System(3, Fermion())
    ψ = normalize(State{Pure}(sys, ["Occ", "Emp", "Occ"]) + State{Pure}(sys, ["Emp", "Emp", "Occ"]))
    @test norm(apply((C + N)(1), ψ) - apply(make_mpo(ψ, (C + N)(1)), ψ)) < 1e-12
    @test_throws "makes it a sum" apply((C + N)(3), ψ)
    # an Evolver is Left + Right of its argument, a sum, which a gate cannot be
    @test_throws "cannot apply sums as gates" apply(Evolver(Z(1)), mix(ψ))
end

@testset "Periods below one mean never" begin
    # one rule for every period of the library: `measures_period`, `n_expand`,
    # `n_hermitianize`, and the `checkpoint_interval` covered in checkpoint.jl.
    # `mod(sweep, 0)` raised a division by zero and `mod(sweep, -2)` is zero on every
    # second sweep, so anything below one is read as never rather than as one of those
    due = TensorMixedStates.sweep_due
    for period in (0, -1, -2, -3)
        @test !any(sweep -> due(period, sweep), 1:12)
    end
    @test filter(sweep -> due(1, sweep), 1:4) == [1, 2, 3, 4]
    @test filter(sweep -> due(3, sweep), 1:10) == [3, 6, 9]

    # end to end, counted in a Data destination, which collects whatever the output of the
    # run is redirected to
    function count_measures(period)
        sim = runTMS(SimData(phases = [
                CreateState{Pure}(2, Qubit(), "Up"),
                Evolve(duration = 0.4, time_step = 0.1, algo = Tdvp(),
                       evolver = -im * Z(1), measures = Data("m") => Z,
                       measures_period = period)]);
            output = devnull)
        return haskey(sim.data, "m") ? length(sim.data["m"]["Z"]["times"]) : 0
    end
    @test count_measures(1) == 4
    @test count_measures(2) == 2
    @test count_measures(0) == 0
    @test count_measures(-2) == 0
end

@testset "Evolving a state that carries charges" begin
    strong = TensorMixedStates.strong
    h = -sum(dag(C)(i) * C(i+1) + dag(C)(i+1) * C(i) for i in 1:3)
    lind = -im * h + sum(Dissipator(sqrt(0.3) * N)(i) for i in 1:4)
    lim = Limits(cutoff = 1e-14, maxdim = 32)
    start(site, mixed) = begin
        sys = System(4, site)
        p(v) = State{Pure}(sys, v)
        s = (p(["Occ", "Emp", "Occ", "Emp"]) + 0.5 * p(["Emp", "Occ", "Occ", "Emp"])) / sqrt(1.25)
        return mixed ? mix(s) : s
    end
    charged = [Fermion(conserve = N), Fermion(conserve = strong(N))]

    # the links of the MPO of the generator carry the charge each of its channels has
    # accumulated, and the answer must not notice. The mixed case uses dephasing, whose jump
    # commutes with the charge, so that it is a strong symmetry as well as a weak one
    for (op, mixed) in ((-im * h, false), (lind, true))
        ref_t = real.(expect1(tdvp(op, 0.2, start(Fermion(), mixed); nsweeps = 1, limits = lim), N))
        ref_w = real.(expect1(approx_W(op, 0.2, start(Fermion(), mixed); limits = lim, order = 2), N))
        for site in charged
            @test real.(expect1(tdvp(op, 0.2, start(site, mixed); nsweeps = 1, limits = lim), N)) ≈ ref_t
            @test real.(expect1(approx_W(op, 0.2, start(site, mixed); limits = lim, order = 2), N)) ≈ ref_w
        end
    end

    # WI and WII share one channel between the term not begun and the term finished, so they
    # need a generator of zero flux. One that moved the charge would not keep the state in
    # its sector, which is why this is a refusal and not a gap
    @test_throws "need an operator of zero flux" approx_W(
        dag(C)(1), 0.1, start(Fermion(conserve = N), false); limits = lim, order = 2)
end
