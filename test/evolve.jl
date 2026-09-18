# Time evolution.
#
# Goes here: the evolution algorithms themselves (Tdvp, ApproxW) on cases whose exact
# solution is known, including dissipative ones. A physical study that happens to use
# an evolution belongs to algorithms.jl instead.

@testset "Simple evolve" begin
    for (algo, time_step, tol) in [
        (Tdvp(), 0.1, 1e-14),
        (ApproxW(order=1, w=1), 0.01, 0.03),
        (ApproxW(order=4, w=1), 0.01, 1e-10),
        (ApproxW(order=1, w=2), 0.01, 1e-13),
        (ApproxW(order=4, w=2), 0.01, 1e-14),       
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
