using TensorMixedStates
using Test

const state_vals = [
    (1, "Qubit", "Z+", [X(1), Y(1), Z(1)], [0, 0, 1]),
    (6, "Qubit", ["X+", "Y+", "Z+", "X-", "Y-", "Z-"], [X, Y, Z],
        [[1, 0, 0, -1, 0, 0], [0, 1, 0, 0, -1, 0], [0, 0, 1, 0, 0, -1]]),
    (3, "Qubit", ["X+", "Y-", "Z-"], [(X, Y), (Z, Z)],
        [[0 -1 0; 0 0 0; 0 0 -im], [1 0 0; 0 1 0; 0 0 1]]),
    (4, "Qubit", ["Z+", "X-", "Z-", "Y-"], [Z(1)Y(4), X(2)Z(3), Z(1)Z(3), X(2)Y(4)], [-1, 1, -1, 1])
]

function check_state_vals()
    for sv in state_vals
        n, s, st, m, v = sv
        sys = System(n, s)
        meas = Measure(m)
        statep = State(Pure, sys, st)
        statem = State(Mixed, sys, st)
        mp = last.(measure(statep, meas))
        mm = last.(measure(statem, meas))
        if mp ≈ v && mm ≈ v
            continue
        end
        error("state test $sv failed with values $mp and $mm instead of $v")
    end
    return true
end

const single_evolve_vals = [
    (2, "Qubit", "X+", -im * Z(1), [X(1), Y(1)], t->[cos(2t), sin(2t)])
]

const multi_evolve_vals = [
    (2, "Qubit", ["X+", "Z-"], -im*Z(1)Z(2), [X(1), Y(1), Z(2)], t->[cos(2t), -sin(2t), -1])
]

const complex_evolve_vals = [
    (5, "Qubit", ["X+", "Z+", "Z-", "X+", "X+"],
    -im*(Z(2)Z(4)+Z(3)Z(1)+Z(3)Z(5)), [Y(1), Y(4), Y(5)], t->sin(2t) * [-1., 1., -1.])
]

function check_evolve_vals(algo, vals; step, atol, kwargs...)
    for ev in vals
        n, s, st, h, m, f = ev
        sys = System(n, s)
        meas = Measure(m)
        statep = State(Pure, sys, st)
        prep = PreMPO(statep, h)
        statem = State(Mixed, sys, st)
        prem = PreMPO(statem, h)
        for t in step:step:1.0
            statep = algo(prep, step, statep; kwargs...)
            statem = algo(prem, step, statem; kwargs...)
            mp = last.(measure(statep, meas))
            mm = last.(measure(statem, meas))
            vf = f(t)
            if ≈(mp, vf; atol) && ≈(mm, vf; atol)
                continue
            end
            error("evolve test $ev failed with values $mp and $mm instead of $vf at t = $t")
        end
    end
    return true
end

@testset verbose=true "TensorMixedStates.jl" begin
    @testset "System building" begin
        @test (System(1, "Qubit"); true)
        @test (System(3, "Qubit"); true)
        @test (System([]); true)
        @test (System(["Qubit"]); true)
        @test (System(["Qubit", "Qubit", "Qubit"]); true)
    end
    @testset "State measuring" begin
        @test check_state_vals()
    end
    @testset "Tdvp single evolution" begin
        @test check_evolve_vals(tdvp, single_evolve_vals; step = 0.1, atol = 1e-8)
    end
    @testset "Tdvp multi evolution" begin
        @test check_evolve_vals(tdvp, multi_evolve_vals; step = 0.1, atol = 1e-8)
    end
    @testset "approx_W1 single evolution" begin
        @test check_evolve_vals(approx_W, single_evolve_vals; step = 0.01, atol = 0.05, order=1, w=1)
        @test check_evolve_vals(approx_W, single_evolve_vals; step = 0.01, atol = 5e-3, order=2, w=1)
        @test check_evolve_vals(approx_W, single_evolve_vals; step = 0.01, atol = 5e-5, order=3, w=1)
        @test check_evolve_vals(approx_W, single_evolve_vals; step = 0.01, atol = 5e-7, order=3, w=1)
    end
    @testset "approx_W2 single evolution" begin
        @test check_evolve_vals(approx_W, single_evolve_vals; step = 0.1, atol = 1e-8, order=1, w=2)
        @test check_evolve_vals(approx_W, single_evolve_vals; step = 0.1, atol = 1e-8, order=2, w=2)
        @test check_evolve_vals(approx_W, single_evolve_vals; step = 0.1, atol = 1e-8, order=3, w=2)
        @test check_evolve_vals(approx_W, single_evolve_vals; step = 0.1, atol = 1e-8, order=4, w=2)
    end
    @testset "approx_W1 multi evolution" begin
        @test check_evolve_vals(approx_W, multi_evolve_vals; step = 0.01, atol = 0.05, order=1, w=2)
        @test check_evolve_vals(approx_W, multi_evolve_vals; step = 0.01, atol = 5e-3, order=2, w=2)
        @test check_evolve_vals(approx_W, multi_evolve_vals; step = 0.01, atol = 5e-5, order=3, w=2)
        @test check_evolve_vals(approx_W, multi_evolve_vals; step = 0.01, atol = 5e-7, order=4, w=2)
    end
    @testset "approx_W2 multi evolution" begin
        @test check_evolve_vals(approx_W, multi_evolve_vals; step = 0.01, atol = 0.05, order=1, w=2)
        @test check_evolve_vals(approx_W, multi_evolve_vals; step = 0.01, atol = 5e-3, order=2, w=2)
        @test check_evolve_vals(approx_W, multi_evolve_vals; step = 0.01, atol = 5e-5, order=3, w=2)
        @test check_evolve_vals(approx_W, multi_evolve_vals; step = 0.01, atol = 5e-7, order=4, w=2)
    end
    @testset "approx_W1 complex evolution" begin
        @test check_evolve_vals(approx_W, complex_evolve_vals; step = 0.01, atol = 0.05, order=1, w=2)
        @test check_evolve_vals(approx_W, complex_evolve_vals; step = 0.01, atol = 5e-3, order=2, w=2)
        @test check_evolve_vals(approx_W, complex_evolve_vals; step = 0.01, atol = 5e-5, order=3, w=2)
        @test check_evolve_vals(approx_W, complex_evolve_vals; step = 0.01, atol = 5e-7, order=4, w=2)
    end
    @testset "approx_W2 complex evolution" begin
        @test check_evolve_vals(approx_W, complex_evolve_vals; step = 0.01, atol = 0.05, order=1, w=2)
        @test check_evolve_vals(approx_W, complex_evolve_vals; step = 0.01, atol = 5e-3, order=2, w=2)
        @test check_evolve_vals(approx_W, complex_evolve_vals; step = 0.01, atol = 5e-5, order=3, w=2)
        @test check_evolve_vals(approx_W, complex_evolve_vals; step = 0.01, atol = 5e-7, order=4, w=2)
    end
end
