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
