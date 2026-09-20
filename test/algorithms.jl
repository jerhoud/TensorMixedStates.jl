# Complete physical scenarios.
#
# Goes here: multi phase simulations reproducing a known physical result, ground state
# search, steady state, and the graph helpers used to build them. These are the slowest
# tests of the suite; keep the systems small.

@testset "Graph utilities" begin
    @test line_graph(4) == [(1, 2), (2, 3), (3, 4)]
    @test circle_graph(4) == [(1, 2), (2, 3), (3, 4), (4, 1)]
    @test complete_graph(4) == [(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]
    @test graph_base_size(circle_graph(7)) == 7
    # the snake runs down the columns: 1 4 5 on the first row, 2 3 6 on the second
    @test square_lattice(3, 2) == [(1, 2), (3, 4), (5, 6), (1, 4), (2, 3), (4, 5), (3, 6)]
    @test square_lattice(3) == square_lattice(3, 3)
    for (nx, ny) in [(2, 2), (3, 2), (4, 3), (10, 3), (1, 5), (5, 1)]
        g = square_lattice(nx, ny)
        @test graph_base_size(g) == nx * ny
        @test length(g) == ny * (nx - 1) + nx * (ny - 1)
        @test allunique(g)
        # this is what the snake buys: no bond spans more than 2ny - 1 sites, whatever
        # the length of the lattice along x
        @test maximum(abs(b - a) for (a, b) in g) ≤ max(1, 2ny - 1)
    end
end

@testset "Complete graphs" begin
    @test_ok test_phases(create_graph_state(complete_graph(4);
        final_measures = check([X, Y, Z, X(1)Z(2)Z(3)Z(4), (Y, Y)],
        [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], 1, [1 1 1 1; 1 1 1 1; 1 1 1 1; 1 1 1 1]])))
    @test_ok test_phases([ create_graph_state(complete_graph(4)), ToMixed(;
        final_measures = check([X, Y, Z, X(1)Z(2)Z(3)Z(4), (Y, Y)],
        [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], 1, [1 1 1 1; 1 1 1 1; 1 1 1 1; 1 1 1 1]]))])

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
                final_measures = check([X, Y, Z, Norm], [[0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [1, 1, 1, 1, 1], 1], 1e-7)
                )
                ])
            end
            @testset "GHZ" begin
                @test_ok begin
                    sys = System(6, Qubit())
                    # adding mixed states adds density matrices, so summing the two mixed
                    # product states gives their classical mixture and not a GHZ state.
                    # The superposition has to be made in pure representation, where the
                    # addition is on amplitudes.
                    ghz = mix((State{Pure}(sys, "Up") + State{Pure}(sys, "Dn")) / sqrt(2))
                    test_phases([
        CreateState(type = Mixed(), state = ghz,
            final_measures = check([Purity, prod(X(i) for i in 1:6)], [1, 1], 1e-10)),
        PartialTrace(
            keep_positions = [2, 3, 5],
            final_measures = check([X, Y, Z, (Z, Z)], [[0, 0, 0], [0, 0, 0], [0, 0, 0], [1 1 1 ; 1 1 1 ; 1 1 1]])
            )
            ])
        end
    end

@testset "Ising chain" begin
    # The reference values are exact, computed by diagonalizing the 64 dimensional Hilbert
    # space in test/reference/ising_ed.jl, which shares no code with what is tested here.
    # The tolerance is therefore the error of the evolution below, about 2e-8, and not the
    # precision of the references. The last one is analytic rather than computed: the
    # product of all X commutes with the hamiltonian and starts at 1, so it stays at 1.
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
        check([X,Y,Z],[[0.48881258418,0.48881258418,0.48881258418,0.48881258418,0.48881258418,0.48881258418],[0.0,0,0,0,0,0],[0.0,0,0,0,0,0]],1e-7),
        check([Z(1)Z(2),Z(2)Z(3),Z(1)Z(6)],[-0.51118741582,-0.51118741582,-0.51118741582],1e-7),
        check([Y(1)Y(2),Y(2)Y(3),Y(1)Y(6)],[-0.2518341076,-0.2518341076,-0.2518341076],1e-7),
        check([X(1)X(2),X(2)X(3),X(1)X(6)],[0.1123126174,0.1123126174,0.1123126174],1e-7),
        check([X(1)X(2)X(3)X(4),X(2)X(3)X(4)X(5),X(4)X(5)X(6)X(1)],[0.1123126174,0.1123126174,0.1123126174],1e-7),
        check(EE(3), 1.15220908566, 1e-7),
        check(X(1)X(2)X(3)X(4)X(5)X(6),1.0,1e-8)
        ])
        ])
    end

@testset "Free fermions with source" begin
    # The reference values are exact, from the dense Lindblad evolution of the 32
    # dimensional Fock space in test/reference/fermion_lindblad.jl, which shares no code
    # with what is tested here. The tolerance is the error of the evolution below, about
    # 8e-8.
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
                check(N, [0.0125582080327,0.063008590052,0.1950187333854,0.063008590052,0.0125582080327], 1e-6),
                check([dag(C)(3)*C(i) for i in 1:5],
                [-0.023529279887, -0.0783682258151im,0.1950187333854,-0.0783682258151im,-0.023529279887],1e-6),
                check(Purity,0.5256648672193,1e-6)
            ])
    ])
end

@testset "Free bosons with source" begin
    # Unlike the two testsets above, these reference values are recorded from a run of the
    # library rather than computed independently: with four sites of dimension 7 the Fock
    # space has 2401 states and the vectorized Liouvillian 2401^2, which puts a dense
    # reference out of reach. A reference built from the Gaussian moments of this quadratic
    # Lindbladian would close the gap and has not been written. Until then, read this
    # testset as a regression check on the behaviour of the day it was recorded.
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
            final_measures = check([X, Y, Z], [[0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [1, 1, 1, 1, 1]], 1e-2)
        )
    ])
    @test_ok test_phases([
        CreateState{Mixed}(4, Qubit(), "FullyMixed"),
        SteadyState(
            lindbladian =
                -im * (sum(X(i)X(i+1)+Y(i)Y(i+1) for i in 1:3))
                + Dissipator(Sp)(1) + Dissipator(Sm)(4),
            nsweeps = 200,
            limits = Limits(cutoff = 1e-10, maxdim = 10),
            final_measures = [
                check(Z, [0.05882352941176472, 0.0, 0.0,-0.05882352941176472], 5e-6),
                check([2(X(i)Y(i+1)-Y(i)X(i+1)) for i in 1:3], fill(0.9411764705882353, 3), 5e-6),
            ]
        )
    ])
end
