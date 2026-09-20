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
    for solver in [tdvp, approx_W]
        chained = foldl(sched; init = st0) do st, m
            solver(-im * h, 0.1, st; nsweeps = 1, limits = Limits(cutoff = 1e-14, maxdim = m))
        end
        scheduled = solver(-im * h, 0.4, st0;
                           nsweeps = 4, limits = Limits(cutoff = 1e-14, maxdim = sched))
        @test maxlinkdim(scheduled) == maxlinkdim(chained)
        @test expect1(scheduled, Z) ≈ expect1(chained, Z)
        # a constant schedule is the plain value, and one shorter than the sweeps keeps
        # its last value for the rest of them
        flat = solver(-im * h, 0.4, st0; nsweeps = 4, limits = Limits(cutoff = 1e-14, maxdim = 4))
        for m in [[4, 4, 4, 4], [4]]
            st = solver(-im * h, 0.4, st0;
                        nsweeps = 4, limits = Limits(cutoff = 1e-14, maxdim = m))
            @test expect1(st, Z) ≈ expect1(flat, Z)
        end
    end
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
    # a factor contributing no tensor, an identity or a Jordan-Wigner string, must not
    # leave the gate list untyped: ITensorMPS.product has no method for a Vector{Any}
    @test_ok apply(Id(1) * X(2), State{Pure}(System(2, Qubit()), "Up"))
    # a gate built from a sum does not distribute and must say so rather than guess
    @test_throws ErrorException Gate(X(1) + Y(2))
    # a non fermionic gate must not be simplified: Swap is defined by an expression and
    # simplifying a product of them would make a sum, which apply cannot place
    @test_ok apply(Swap(1, 2) * Swap(3, 4), State{Pure}(System(4, Qubit()), "Up"))
end
