using TensorMixedStates, .Qubits, .Fermions, .Bosons, .Spins, .Electrons, .Tjs, .Qbosons
using Test, Aqua

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
test a statement (with `@test_ok`) with `type = Pure` and `type = Mixed` 
"""
macro test_pm(a)
    tp = :type
    quote
        @test_ok (($(esc(tp)) -> $(esc(a)))(Pure))
        @test_ok (($(esc(tp)) -> $(esc(a)))(Mixed))
    end
end

"""
checks whether measurements are equal between a random pure state and its computed mixed represenetation
"""
function check_mix(dims, sites, measures)
    for d in dims, (s, m) in zip(sites, measures)
        sys = System(10, s)
        stp = RandomState{Pure}(sys, d)
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
    @testset "Aqua" begin
        Aqua.test_all(TensorMixedStates)
    end
    @testset "System building" begin
        @test_ok System(1, Qubit())
        @test_ok System(3, Qubit())
        @test_ok System([Qubit()])
        @test_ok System([Qubit(), Qubit(), Qubit()])
    end
    @testset "State building" begin
        @test_pm State{type}(System(1, Qubit()), "Up")
        @test_pm State{type}(3, Qubit(), "Up")
        @test_pm State{type}([Qubit(), Fermion(), Boson(4)], ["Up", "Occ", "2"])
        @test_pm State{type}(System(3, Qubit()), [1., 0.])
        @test_ok State{Pure}(System(3, Qubit()), [[1., 0.], [0, 1], [0, im]])
        @test_ok State{Mixed}(System(3, Qubit()), [[1 0 ; 0 0 ], [1 1; 1 1], [0 0 ; 0 1]])
        @test_ok State{Mixed}(System(3, Qubit()), "FullyMixed")
    end
    @testset "Simulation building" begin
       @test_ok Simulation(nothing)
       @test_pm Simulation(State{type}(System(3, Qubit()), "Up"))
    end
    @testset "Qubit measuring" begin
        @test_pm test_phases(CreateState{type}(1, Qubit(), "Z+"; 
            final_measures = check([X(1), Y(1), Z(1)], [0, 0, 1])))
        @test_pm test_phases(CreateState{type}(6, Qubit(), ["X+", "Y+", "Z+", "X-", "Y-", "Z-"];
            final_measures = check([X, Y, Z], [[1, 0, 0, -1, 0, 0], [0, 1, 0, 0, -1, 0], [0, 0, 1, 0, 0, -1]])))
        @test_pm test_phases(CreateState{type}(3, Qubit(), ["X+", "Y-", "Z-"];
            final_measures = check([(X, Y), (Z, Z)], [[0 -1 0; 0 0 0; 0 0 -im], [1 0 0; 0 1 0; 0 0 1]])))
        @test_pm test_phases(CreateState{type}(4, Qubit(), ["Z+", "X-", "Z-", "Y-"];
            final_measures = check([Z(1)Y(4), X(2)Z(3), Z(1)Z(3), X(2)Y(4)], [-1, 1, -1, 1])))
        @test norm(matrix(Sx^2+Sy^2+Sz^2-S2, Qubit()))≈0 atol=1e-12
    end
    @testset "Fermion measuring" begin
        @test_pm test_phases(CreateState{type}(1, Fermion(), "1";
            final_measures = check(N(1), 1)))
        @test_pm test_phases(CreateState{type}(2, Fermion(), ["0", "1"];
            final_measures = check(N, [0, 1])))
    end
    @testset "Boson measuring" begin
        @test_pm test_phases(CreateState{type}(4, Boson(4), ["0", "1", "2", "3"];
            final_measures = check(N, [0, 1, 2, 3])))
    end
    @testset "Spin measuring" begin
        @test_pm test_phases(CreateState{type}(4, Spin(3/2), ["-3/2", "-1/2", "1/2", "3/2"];
            final_measures = check([Sx, Sy, Sz], [[0, 0, 0, 0], [0, 0, 0, 0], [-3/2, -1/2, 1/2, 3/2]])))
        @test norm(matrix(Sx^2+Sy^2+Sz^2-S2,Spin(5/2)))≈0 atol=1e-12
        @test norm(matrix(Sx^2+Sy^2+Sz^2-S2,Spin(4)))≈0 atol=1e-12
        @test_pm test_phases(CreateState{type}(9, Spin(1),
                ["X-1", "X0", "X1", "Y-1", "Y0", "Y1", "Z-1", "Z0", "Z1"];
            final_measures = check([Sx, Sy, Sz],
                [[-1, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, -1, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, -1, 0, 1]])))
        @test_pm test_phases(CreateState{type}(12, Spin(3/2),
                ["X-3/2", "X-1/2", "X1/2", "X3/2", "Y-3/2", "Y-1/2", "Y1/2", "Y3/2", "Z-3/2", "Z-1/2", "Z1/2", "Z3/2"];
            final_measures = check([Sx, Sy, Sz],
                [[-3/2, -1/2, 1/2, 3/2, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, -3/2, -1/2, 1/2, 3/2, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0, -3/2, -1/2, 1/2, 3/2]])))
    end
    @testset "Electron measuring" begin
        @test_pm test_phases(CreateState{type}(4, Electron(), ["Emp", "Up", "Dn", "UpDn"];
            final_measures = check([Nup, Ndn, Nupdn, Ntot], [[0, 1, 0, 1], [0, 0, 1, 1], [0, 0, 0, 1], [0, 1, 1, 2]])))        
    end
    @testset "Tj measuring" begin
        @test_pm test_phases(CreateState{type}(3, Tj(), ["Emp", "Up", "Dn"];
            final_measures = check([Nup, Ndn, Ntot], [[0, 1, 0], [0, 0, 1], [0, 1, 1]])))        
    end
    @testset "Qboson measuring" begin
        @test_pm test_phases(CreateState{type}(4, Qboson(0.1, 4), ["0", "1", "2", "3"];
            final_measures = check(N, [0, 1, 2, 3])))
    end
    @testset "State mixing" begin
        @test_ok check_mix(
            [1, 5, 15],
            [Qubit(), Fermion(), Boson(3), Spin(3/2), Electron(), Tj()],
            [[X, Y, Z], [N], [N], [Sx, Sy, Sz], [Nup, Ndn, Nupdn, Ntot, Sx, Sy, Sz], [Nup, Ndn, Ntot, Sx, Sy, Sz]])
    end
    @testset "Complete graphs" begin
        @test_ok test_phases(create_graph_state(complete_graph(4);
            final_measures = check([X, Y, Z, X(1)Z(2)Z(3)Z(4), (Y, Y)],
            [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], 1, [1 1 1 1; 1 1 1 1; 1 1 1 1; 1 1 1 1]])))
        @test_ok test_phases([ create_graph_state(complete_graph(4)), ToMixed(;
            final_measures = check([X, Y, Z, X(1)Z(2)Z(3)Z(4), (Y, Y)],
            [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], 1, [1 1 1 1; 1 1 1 1; 1 1 1 1; 1 1 1 1]]))])

    end
    
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
    @testset "Dmrg" begin
        @test_ok test_phases([
            CreateState(
                type = Pure(),
                system = System(5, Qubit()),
                randomize = 10,
                ),
                GroundState(
                    hamiltonian = sum(-Z(i) for i in 1:5),
                    limits = Limits(maxdim = 10),
                    nsweeps = 2,
                    final_measures = check([X, Y, Z, Norm], [[0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [1, 1, 1, 1, 1], 1], 5e-8)
                    )
                    ])
                end
                @testset "GHZ" begin
                    @test_ok begin
                        sys = System(6, Qubit())
                        ghz = (State{Mixed}(sys, "Up") + State{Mixed}(sys, "Dn")) / 2
                        test_phases([
            CreateState(type = Mixed(), state = ghz),
            PartialTrace(
                keep_positions = [2, 3, 5],
                final_measures = check([X, Y, Z, (Z, Z)], [[0, 0, 0], [0, 0, 0], [0, 0, 0], [1 1 1 ; 1 1 1 ; 1 1 1]])
                )
                ])
            end
        end
    @testset "Ising chain" begin
        @test_ok test_phases([
            CreateState{Pure}(6, Qubit(), "X+"),
            Evolve(
        algo =  ApproxW(order = 4, w = 2),
        limits = Limits(maxdim = 8),
        duration = 1.0,
        time_step = 0.02,
        evolver =
            -im*(sum(Z(i)*Z(i+1) for i in 1:5)+Z(6)*Z(1)-sum(X(i) for i in 1:6)),
        final_measures = [
            check([X,Y,Z],[[0.48881258678,0.48881258678,0.48881258678,0.48881258678,0.48881258678,0.48881258678],[0.0,0,0,0,0,0],[0.0,0,0,0,0,0]],1e-7),
            check([Z(1)Z(2),Z(2)Z(3),Z(1)Z(6)],[-0.51118739903,-0.51118739903,-0.51118739903],1e-7),
            check([Y(1)Y(2),Y(2)Y(3),Y(1)Y(6)],[-0.25183410946,-0.25183410946,-0.25183410946],1e-7),
            check([X(1)X(2),X(2)X(3),X(1)X(6)],[0.11231262241,0.11231262241,0.11231262241],1e-7),
            check([X(1)X(2)X(3)X(4),X(2)X(3)X(4)X(5),X(4)X(5)X(6)X(1)],[0.11231262241,0.11231262241,0.11231262241],1e-7),
            check(EE(3), 1.15220908795, 1e-7),
            check(X(1)X(2)X(3)X(4)X(5)X(6),1.0,1e-8)
            ])
            ])
        end
    @testset "Free fermions with source" begin
        @test_ok test_phases([
            CreateState{Mixed}(5, Fermion(), "0"),
            Evolve(
                algo=Tdvp(),
                limits = Limits(maxdim = 16),
                duration = 1.0,
                time_step = 0.05,
                evolver =
                    -im * sum(dag(C)(i)*C(i+1)+dag(C)(i+1)*C(i) for i in 1:4) + Dissipator(sqrt(2*0.2)*dag(C))(3),
                final_measures = [
                    check(N, [0.125581972006622e-1,0.630086192053610e-1,.195018685802510,0.630086192053610e-1,0.125581972006622e-1], 1e-6),
                    check([dag(C)(3)*C(i) for i in 1:5],
                    [-0.235293024730342e-1, -0.783683050686146e-1*im,.195018685802510,-0.783683050686146e-1*im,-0.235293024730342e-1],1e-6),
                    check(Purity,.525664902939503,1e-6)
                ])
        ])
    end
    @testset "Free bosons with source" begin
        @test_ok test_phases([
            CreateState{Mixed}(4, Boson(7), "0"),
            Evolve(
                algo = ApproxW(order = 4, w = 2),
                limits = Limits(maxdim = 10),
                duration = 0.3,
                time_step = 0.1,
                evolver =
                    -im*sum(A(i)*dag(A)(i+1)+dag(A)(i)*A(i+1) for i in 1:3) + Dissipator(2*sqrt(0.1)*dag(A))(2),
                final_measures = [
                        check(N,[0.363288753916464e-2,.120023637053104,0.356800815401577e-2,0.486649287358378e-4],1e-6),
                        check([dag(A)(2)*A(i) for i in 1:4],[-0.180013492215079e-1*im,.120023637053104,-0.178675615179132e-1*im,-0.178658178615343e-2],1e-5),
                        check(Purity,.7965328508313,1e-6)
                    ])
        ])
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
    @testset "Steady state" begin
        @test_ok test_phases([
            CreateState{Mixed}(2, Qubit(), "+"),
            SteadyState(
                lindbladian = Dissipator(Sp)(1) + Dissipator(Sm)(2),
                nsweeps = 20,
                limits = Limits(cutoff = 1e-10, maxdim = 10),
                final_measures = check([X, Y, Z], [[0, 0], [0, 0], [1, -1]], 1e-2)
            )
        ])
    end
    @test_ok test_phases([
        CreateState(
            type = Mixed(),
            system = System(5, Qubit()),
            randomize = 10,
        ),
        SteadyState(
            lindbladian = -im * (-sum(Z(i)Z(i+1) for i in 1:4)) + sum(Dissipator(Sp)(i) for i in 1:5),
            limits = Limits(maxdim = 10, cutoff = 1e-10),
            nsweeps = 20,
            final_measures = check([X, Y, Z], [[0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [1, 1, 1, 1, 1]], 1e-6)
        )
    ])
end
