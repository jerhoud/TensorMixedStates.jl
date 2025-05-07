using TensorMixedStates, .Qubits, .Fermions, .Bosons, .Spins, .Electrons, .Tjs
using Test

"""
executes the given phases
"""
function test_phases(phases)
    try
        runTMS(SimData(;phases); output = devnull)
    catch
        println("test_phases failed for $phases")
        rethrow()
    end
end

check(a, b, tol=1e-8) = "" => Check("", a, b, tol)

"""
test whether a statement execute without throwing an exception
"""
macro test_ok(a)
    :(@test ($(esc(a)); true))
end

"""
test a statement (with `@test_ok`) with `type = Pure()` and `type = Mixed()` 
"""
macro test_pm(a)
    tp = :type
    quote
        @test_ok (($(esc(tp)) -> $(esc(a)))(Pure()))
        @test_ok (($(esc(tp)) -> $(esc(a)))(Mixed()))
    end
end

"""
checks whether measurements are equal between a random pure state and its computed mixed represenetation
"""
function check_mix(dims, sites, measures)
    for d in dims, (s, m) in zip(sites, measures)
        sys = System(10, s)
        stp = random_state(Pure(), sys, d)
        stm = mix(stp)
        meas = Measure(m...)
        mp = last.(measure(stp, meas))
        mm = last.(measure(stm, meas))
        if mp ≈ mm
            continue
        end
        error("mix check $d, $s, $m fails with $mp and $mm")
    end
end

@testset verbose=true "TensorMixedStates.jl" begin
    @testset "System building" begin
        @test_ok System(1, Qubit())
        @test_ok System(3, Qubit())
        @test_ok System([Qubit()])
        @test_ok System([Qubit(), Qubit(), Qubit()])
    end
    @testset "State building" begin
        @test_pm State(type, System(1, Qubit()), "Up")
        @test_pm State(type, 3, Qubit(), "Up")
        @test_pm State(type, [Qubit(), Fermion(), Boson(4)], ["Up", "Occ", "2"])
        @test_pm State(type, System(3, Qubit()), [1., 0.])
        @test_ok State(Pure(), System(3, Qubit()), [[1., 0.], [0, 1], [0, im]])
        @test_ok State(Mixed(), System(3, Qubit()), [[1 0 ; 0 0 ], [1 1; 1 1], [0 0 ; 0 1]])
        @test_ok State(Mixed(), System(3, Qubit()), "FullyMixed")
    end
    @testset "Simulation building" begin
       @test_ok Simulation(nothing)
       @test_pm Simulation(State(type, System(3, Qubit()), "Up"))
    end
    @testset "Qubit measuring" begin
        @test_pm test_phases(CreateState(type, 1, Qubit(), "Z+"; 
            final_measures = check([X(1), Y(1), Z(1)], [0, 0, 1])))
        @test_pm test_phases(CreateState(type, 6, Qubit(), ["X+", "Y+", "Z+", "X-", "Y-", "Z-"];
            final_measures = check([X, Y, Z], [[1, 0, 0, -1, 0, 0], [0, 1, 0, 0, -1, 0], [0, 0, 1, 0, 0, -1]])))
        @test_pm test_phases(CreateState(type, 3, Qubit(), ["X+", "Y-", "Z-"];
            final_measures = check([(X, Y), (Z, Z)], [[0 -1 0; 0 0 0; 0 0 -im], [1 0 0; 0 1 0; 0 0 1]])))
        @test_pm test_phases(CreateState(type, 4, Qubit(), ["Z+", "X-", "Z-", "Y-"];
            final_measures = check([Z(1)Y(4), X(2)Z(3), Z(1)Z(3), X(2)Y(4)], [-1, 1, -1, 1])))
    end
    @testset "Fermion measuring" begin
        @test_pm test_phases(CreateState(type, 1, Fermion(), "1";
            final_measures = check(N(1), 1)))
        @test_pm test_phases(CreateState(type, 2, Fermion(), ["0", "1"];
            final_measures = check(N, [0, 1])))
    end
    @testset "Boson measuring" begin
        @test_pm test_phases(CreateState(type, 4, Boson(4), ["0", "1", "2", "3"];
            final_measures = check(N, [0, 1, 2, 3])))
    end
    @testset "Spin measuring" begin
        @test_pm test_phases(CreateState(type, 4, Spin(3/2), ["-3/2", "-1/2", "1/2", "3/2"];
            final_measures = check([Sx, Sy, Sz], [[0, 0, 0, 0], [0, 0, 0, 0], [-3/2, -1/2, 1/2, 3/2]])))        
    end
    @testset "Electron measuring" begin
        @test_pm test_phases(CreateState(type, 4, Electron(), ["Emp", "Up", "Dn", "UpDn"];
            final_measures = check([Nup, Ndn, Nupdn, Ntot], [[0, 1, 0, 1], [0, 0, 1, 1], [0, 0, 0, 1], [0, 1, 1, 2]])))        
    end
    @testset "Tj measuring" begin
        @test_pm test_phases(CreateState(type, 3, Tj(), ["Emp", "Up", "Dn"];
            final_measures = check([Nup, Ndn, Ntot], [[0, 1, 0], [0, 0, 1], [0, 1, 1]])))        
    end
    @testset "State mixing" begin
        @test_ok check_mix(
            [1, 5, 15],
            [Qubit(), Fermion(), Boson(3), Spin(3/2), Electron(), Tj()],
            [[X, Y, Z], [N], [N], [Sx, Sy, Sz], [Nup, Ndn, Nupdn, Ntot, Sx, Sy, Sz], [Nup, Ndn, Ntot, Sx, Sy, Sz]])
    end
    @testset "Complete graphs" begin
        @test_pm test_phases(create_graph_state(type, complete_graph(4);
            final_measures = check([X, Y, Z, X(1)Z(2)Z(3)Z(4), (Y, Y)],
            [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], 1, [1 1 1 1; 1 1 1 1; 1 1 1 1; 1 1 1 1]])))
    end
    @testset "Simple evolve" begin
        for (algo, time_step, tol) in [
            (Tdvp(), 0.1, 1e-8),
            (ApproxW(order=1, w=1), 0.01, 0.05),
            (ApproxW(order=4, w=1), 0.01, 5e-7),
            (ApproxW(order=1, w=2), 0.01, 1e-8),
            (ApproxW(order=4, w=2), 0.01, 1e-8),       
        ]
            @test_pm test_phases([
                CreateState(type, 2, Qubit(), "X+"),
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
            (Tdvp(), 0.1, 1e-8),
            (ApproxW(order=1, w=1), 0.01, 0.05),
            (ApproxW(order=4, w=1), 0.01, 5e-7),
            (ApproxW(order=1, w=2), 0.01, 0.05),
            (ApproxW(order=4, w=2), 0.01, 5e-7),       
        ]
            @test_pm test_phases([
                CreateState(type, 2, Qubit(), ["X+", "Z-"]),
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
            (Tdvp(), 0.1, 1e-8),
            (ApproxW(order=1, w=1), 0.01, 0.05),
            (ApproxW(order=4, w=1), 0.01, 5e-7),
            (ApproxW(order=1, w=2), 0.01, 0.05),
            (ApproxW(order=4, w=2), 0.01, 5e-7),       
        ]
            @test_pm test_phases([
                CreateState(type, 5, Qubit(), ["X+", "Z+", "Z-", "X+", "X+"]),
                Evolve(; algo, time_step,
                    limits = Limits(maxdim = 10, cutoff = 1e-15),
                    evolver = -im * (Z(2)Z(4) + Z(3)Z(1) + Z(3)Z(5)),
                    duration = 1,
                    final_measures = check([Y(1), Y(4), Y(5)], t->sin(2t) * [-1., 1., -1.], tol)
                )
            ])
        end
    end
    @testset "Dmrg" begin
        @test_ok test_phases([
            CreateState(
                type = Pure(),
                system = System(5, Qubit()),
                randomize = 10,
            ),
            Dmrg(
                hamiltonian = sum(-Z(i) for i in 1:5),
                limits = Limits(maxdim = 10),
                nsweeps = 2,
                final_measures = check([X, Y, Z, Norm], [[0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [1, 1, 1, 1, 1], 1], 1e-7)
            )
        ])
    end
    @testset "GHZ" begin
        @test_ok begin
            sys = System(6, Qubit())
            ghz = (State(Mixed(), sys, "Up") + State(Mixed(), sys, "Dn")) / 2
            test_phases([
                CreateState(state = ghz),
                Partial_Trace(
                    keep_positions = [2, 3, 5],
                    final_measures = check([X, Y, Z, (Z, Z)], [[0, 0, 0], [0, 0, 0], [0, 0, 0], [1 1 1 ; 1 1 1 ; 1 1 1]])
                )
            ])
        end
    end
end
