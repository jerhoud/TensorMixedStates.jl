# Complete physical scenarios.
#
# Goes here: multi phase simulations reproducing a known physical result, ground state
# search, steady state, and the graph helpers used to build them. These are the slowest
# tests of the suite; keep the systems small.

@testset "Graph utilities" begin
    @test line_graph(4) == [(1, 2), (2, 3), (3, 4)]
    @test circle_graph(4) == [(1, 2), (2, 3), (3, 4), (4, 1)]
    @test circle_graph(2) == [(1, 2), (2, 1)]
    @test_throws "a ring has two vertices at least" circle_graph(1)
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

@testset "A failed Check fails test_phases" begin
    # what every test_phases of a check relies on: a Check that fails stops the run with its
    # error, rather than being logged. The tolerance bounds the norm of the difference
    @test_throws "failed with values" redirect_stdout(devnull) do
        test_phases([CreateState{Pure}(2, Qubit(), "Up"),
                     Evolve(duration = 0.2, time_step = 0.1, algo = Tdvp(), evolver = -im * Z(1),
                            final_measures = check(Z(1), -1))])
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
        check(EntanglementEntropy(3), 1.15220908566, 1e-7),
        check(X(1)X(2)X(3)X(4)X(5)X(6),1.0,1e-8)
        ])
        ])
    end

@testset "An expansion before the first step" begin
    # the ring above from a product state, with tdvp: the first step, of bond dimension one,
    # left the tangent space through Z(6)Z(1), and the expansion, coming after it, left an
    # error of order the time step, 0.016 on X(1), whatever n_expand
    for n_expand in (1, 2)
        @test_ok test_phases([
            CreateState{Pure}(6, Qubit(), "X+"),
            Evolve(algo = Tdvp(; n_expand), limits = Limits(maxdim = 8, cutoff = 1e-14),
                   duration = 1.0, time_step = 0.1,
                   evolver = -im * (sum(Z(i) * Z(i + 1) for i in 1:5) + Z(6) * Z(1) -
                                    sum(X(i) for i in 1:6)),
                   final_measures = [check(X(1), 0.48881258418, 1e-8),
                                     check(Z(1)Z(6), -0.51118741582, 1e-8)])])
    end
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
    # The reference values are exact, from the Gaussian moments of this quadratic Lindbladian
    # in test/reference/boson_gaussian.jl, which shares no code with what is tested here: with
    # four sites of dimension 7 a dense reference, on a vectorized Liouvillian of 2401^2, is
    # out of reach. The evolution below is 6e-7 to 8e-7 away from them.
    #
    # The bond dimension is 16 rather than the 10 first used, and that is not a detail. At
    # 10 the truncation itself is unstable: which singular values survive depends on the
    # rounding, so the answer follows the order the sums happen to be taken in. Changing
    # nothing but the BLAS thread count on one machine moved this correlation by 4.4e-4
    # relative, an overall deviation swinging between 8.8e-7 and 8.0e-6, and the Windows
    # job, landing on a different thread count from the Linux and macOS ones, missed the
    # 1e-5 tolerance by one percent. At 16 the dependence is gone, the deviation being
    # 8.6e-7 whatever the thread count, so the tolerance measures the accuracy of the
    # method again rather than the arithmetic of the runner. It is left at 1e-5 and not
    # tightened, since 8.6e-7 leaves too little margin under 1e-6.
    @test_ok test_phases([
        CreateState{Mixed}(4, Boson(7), "0"),
        Evolve(
            algo = ApproxW(order = 4, w = 2),
            limits = Limits(maxdim = 16),
            duration = 0.3,
            time_step = 0.1,
            evolver =
                -im*sum(A(i)*dag(A)(i+1)+dag(A)(i)*A(i+1) for i in 1:3) + Dissipator(2*sqrt(0.1)*dag(A))(2),
            final_measures = [
                    check(N, [0.0036328913502, 0.1200236303285, 0.0035680174576, 4.86613106e-5], 1e-5),
                    check([dag(A)(2)*A(i) for i in 1:4],
                          [-0.0180013634952im, 0.1200236303285, -0.0178675798942im, -0.0017865872242], 1e-5),
                    check(Purity, 0.7965328475916, 1e-5)
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
            final_measures = check([X, Y, Z], [[0, 0], [0, 0], [1, -1]], 1e-10)
        )
    ])
    @test_ok test_phases([
        CreateState(
            type = Mixed(),
            system = System(5, Qubit()),
            randomize = 10,
        ),
        # the tolerances of these three searches hid an error of the default Krylov search of
        # ITensorMPS, 1e-5 on the first, which the default of steady_state now resolves
        SteadyState(
            lindbladian = -im * (-sum(Z(i)Z(i+1) for i in 1:4)) + sum(Dissipator(Sp)(i) for i in 1:5),
            limits = Limits(maxdim = 10, cutoff = 1e-10),
            nsweeps = 40,
            final_measures = check([X, Y, Z], [[0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [1, 1, 1, 1, 1]], 1e-8)
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
                check(Z, [0.05882352941176472, 0.0, 0.0,-0.05882352941176472], 1e-7),
                check([2(X(i)Y(i+1)-Y(i)X(i+1)) for i in 1:3], fill(0.9411764705882353, 3), 1e-7),
            ]
        )
    ])
    # the state found is a density matrix, of trace one, where the eigenvector of (L+)L has
    # norm one and a sign of its own
    ρ0 = mix(State{Pure}(System(2, Qubit()), "+"))
    _, ρ = steady_state(Dissipator(Sp)(1) + Dissipator(Sm)(2), ρ0;
                        nsweeps = 10, limits = Limits(cutoff = 1e-10, maxdim = 10))
    @test trace(ρ) ≈ 1
end

@testset "The options of steady_state" begin
    # it resumes at `first_sweep` as dmrg does, which is how a SteadyState phase continues
    # after a checkpoint, and takes a noise as GroundState does
    L = Dissipator(Sp)(1) + Dissipator(Sm)(2)
    ρ0 = State{Mixed}(System(2, Qubit()), "FullyMixed")
    lim = Limits(maxdim = 16)
    for first_sweep in 1:3
        obs = TensorMixedStates.ITensorMPS.DMRGObserver()
        steady_state(L, ρ0; nsweeps = 3, first_sweep, limits = lim, observer! = obs)
        @test length(obs.energies) == 4 - first_sweep
    end
    _, ρ = steady_state(L, ρ0; nsweeps = 2, limits = lim, noise = 1e-6)
    @test expect1(ρ, Z) ≈ [1, -1]
end

@testset "dmrg with no sweep left" begin
    # a search resumed after its last sweep is asked for none, on which ITensorMPS gave an
    # energy of 0: it is that of the state given
    h = -sum(Z(i) for i in 1:2) - 0.5 * X(1)
    st = State{Pure}(System(2, Qubit()), "Up")
    e, s = dmrg(h, st; nsweeps = 2, first_sweep = 3)
    @test e ≈ real(expect(st, h))
    @test s === st
end

@testset "The Krylov parameters reach dmrg" begin
    # a Krylov space of a single vector holds nothing but the state it starts from, so that
    # the search does not leave it: dmrg ends where it started
    n = 4
    h = sum(X(i) * X(i + 1) for i in 1:n - 1)
    up = State{Pure}(System(n, Qubit()), "Up")
    lim = Limits(maxdim = 16)
    @test first(dmrg(h, up; nsweeps = 3, limits = lim)) ≈ -3
    e, st = dmrg(h, up; nsweeps = 3, limits = lim, krylov = Krylov(dim = 1))
    @test e ≈ 0 atol = 1e-12
    @test expect1(st, Z) ≈ ones(n)
    # and steady_state, from the fully mixed state
    L = Dissipator(Sp)(1) + Dissipator(Sm)(2)
    ρ0 = State{Mixed}(System(2, Qubit()), "FullyMixed")
    _, ρ = steady_state(L, ρ0; nsweeps = 2, limits = lim, krylov = Krylov(dim = 1))
    @test expect1(ρ, Z) ≈ [0, 0] atol = 1e-12
    # the phases hand them on: the state stays up through both searches, where the ground
    # state of `h` has no magnetization and the steady state of `L` a spin down on site 2
    sim = runTMS(SimData(phases = [
            CreateState(type = Pure(), state = up),
            GroundState(hamiltonian = h, limits = lim, nsweeps = 3, krylov = Krylov(dim = 1)),
            ToMixed(),
            SteadyState(lindbladian = L, limits = lim, nsweeps = 2, krylov = Krylov(dim = 1))]);
        output = devnull)
    @test expect1(sim.state, Z) ≈ ones(n)
end

@testset "dmrg and steady_state refuse the options of ITensorMPS" begin
    # they no longer pass on what they do not know, `outputlevel` included
    up = State{Pure}(System(2, Qubit()), "Up")
    ρ0 = State{Mixed}(System(2, Qubit()), "FullyMixed")
    @test_throws MethodError dmrg(X(1) * X(2), up; outputlevel = 1)
    @test_throws MethodError steady_state(Dissipator(Sm)(1), ρ0; outputlevel = 1)
end

@testset "Dmrg of a hamiltonian on a mixed state" begin
    # it would minimise ρ ↦ Hρ + ρH, whose lowest eigenvector is neither the ground state nor a
    # density matrix, and is refused. A superoperator given as such is left to the caller
    h = -Z(1) * Z(2) - 0.5 * (X(1) + X(2))
    ρ = mix(RandomState{Pure}(System(2, Qubit()), 2))
    lim = Limits(maxdim = 4)
    @test_throws "ground state of a pure state" dmrg(h, ρ; nsweeps = 2, limits = lim)
    @test_throws "ground state of a pure state" runTMS(SimData(phases = [
        CreateState{Mixed}(2, Qubit(), "+"),
        GroundState(hamiltonian = h, nsweeps = 2, limits = lim)]); output = devnull)
    @test_ok dmrg(sum(Left(Z)(i) + Right(Z)(i) for i in 1:2), ρ; nsweeps = 2, limits = lim)
end

@testset "Ground and steady states on a charged system" begin
    strong = TensorMixedStates.strong
    h = -sum(dag(C)(i) * C(i+1) + dag(C)(i+1) * C(i) for i in 1:3)
    lind = -im * h + sum(Dissipator(sqrt(0.3) * N)(i) for i in 1:4)
    lim = Limits(cutoff = 1e-12, maxdim = 32)
    start(site) = begin
        sys = System(4, site)
        p(v) = State{Pure}(sys, v)
        return (p(["Occ", "Emp", "Occ", "Emp"]) + 0.5 * p(["Emp", "Occ", "Occ", "Emp"])) / sqrt(1.25)
    end

    # dmrg searches inside the sector its starting state lives in, and steady_state builds
    # the MPO of `(L+)L`, whose flux is zero whenever that of `L` is. Both must land where
    # the dense computation does
    e = first(dmrg(h, start(Fermion()); nsweeps = 3, limits = lim))
    # the two modes of negative energy filled, -2cos(π/5) - 2cos(2π/5)
    @test e ≈ -sqrt(5)
    z = real(first(steady_state(lind, mix(start(Fermion())); nsweeps = 2, limits = lim)))
    for site in (Fermion(conserve = N), Fermion(conserve = strong(N)))
        @test first(dmrg(h, start(site); nsweeps = 3, limits = lim)) ≈ e
        @test real(first(steady_state(lind, mix(start(site)); nsweeps = 2, limits = lim))) ≈ z atol = 1e-10
    end
end

@testset "Weakening between two phases" begin
    strong = TensorMixedStates.strong
    # dephasing commutes with the number of particles and allows a strong symmetry, loss does
    # not. Weakening in between must give what a weak symmetry gives from the start
    h = -sum(dag(C)(i) * C(i+1) + dag(C)(i+1) * C(i) for i in 1:3)
    dephasing = -im * h + sum(Dissipator(sqrt(0.3) * N)(i) for i in 1:4)
    loss = -im * h + Dissipator(sqrt(0.2) * C)(2)
    evolve(ev) = Evolve(algo = Tdvp(), limits = Limits(cutoff = 1e-12, maxdim = 32),
                        duration = 0.5, time_step = 0.05, evolver = ev)
    run(site, middle) = runTMS(SimData(phases = [
            CreateState{Mixed}(4, site, ["Occ", "Emp", "Occ", "Emp"]),
            evolve(dephasing), middle..., evolve(loss) ]); output = devnull).state
    numbers(s) = [ real(trace(s)); real.(expect1(s, N)); expect(s, dag(C)(1) * C(3)) ]
    reference = numbers(run(Fermion(conserve = N), []))

    s = run(Fermion(conserve = strong(N)), [Weaken()])
    @test repr(symmetries(s.system)) == "N"
    @test numbers(s) ≈ reference
    s = run(Fermion(conserve = strong(N)), [Weaken(target = ())])
    @test !TensorMixedStates.is_charged(s.system)
    @test numbers(s) ≈ reference

    @test_throws "drop `strong`" run(Fermion(conserve = strong(N)), [])
end
