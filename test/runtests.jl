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
        error("test $sv failed with values $mp and $mm")
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
end
