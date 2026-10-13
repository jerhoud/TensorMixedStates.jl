# The precompilation workload, small computations run when the package is precompiled, so that a
# first script does not wait for them to compile; set TMS_SKIP_PRECOMPILE_WORKLOAD=true to skip it
# during development, never for a release.

using PrecompileTools: @setup_workload, @compile_workload

if get(ENV, "TMS_SKIP_PRECOMPILE_WORKLOAD", "false") ≠ "true"
    @setup_workload begin
        # four sites, for a bond away from both ends, as most are: its tensors have more indices,
        # and their contractions compile apart
        system = System(4, Qubit())
        h = -sum(Z(i) * Z(i + 1) for i in 1:3) - 0.5 * sum(X(i) for i in 1:4)
        lim = Limits(maxdim = 4)
        @compile_workload begin
            # silent: the warnings of a computation, and the line runTMS ends with
            Logging.with_logger(Logging.NullLogger()) do
                # pure and mixed states, of real tensors and of complex ones, product states and
                # entangled ones: each takes paths of its own
                ψ = random_state(Pure, system, 4; rng = Xoshiro(1))
                up = State{Pure}(system, "Up")
                ρ = mix(ψ)
                ρup = State{Mixed}(system, "Up")
                m = Measure(X, (X, Y), Z(1) * Z(2), h, Purity, Trace, EntanglementEntropy(1))
                for st in (ψ, up, ρ, ρup)
                    measure(st, m; time = 0.)
                    expect(st, h)
                end
                expect1(ψ, Z)
                expect2(ψ, (Z, Z))
                for st in (ψ, up)
                    evolve(-im * h, 0.1, st; limits = lim)
                    evolve(-im * h, 0.1, st; algo = ApproxW(order = 2), limits = lim)
                end
                dissipative = -im * h + Dissipator(Sm)(1)
                for st in (ρ, ρup)
                    evolve(dissipative, 0.1, st; limits = lim)
                    evolve(dissipative, 0.1, st; algo = ApproxW(order = 2), limits = lim)
                end
                ground_state(h, ψ; algo = Dmrg(nsweeps = 2), limits = lim)
                steady_state(dissipative, ρup; algo = Dmrg(nsweeps = 1), limits = lim)
                thermalize(h, 0.1, State{Mixed}(system, "FullyMixed"); limits = lim)
                apply(X(1) * Z(2), ψ)
                apply(Gate(X)(1), ρ)
                # a conserved quantity makes the tensors block sparse, on a density matrix
                # strongly or weakly, and fermions bring their Jordan-Wigner strings
                alternating = ["Occ", "Emp", "Occ", "Emp"]
                hopping = -sum(dag(C)(i) * C(i + 1) + dag(C)(i + 1) * C(i) for i in 1:3)
                ψf = State{Pure}(System(4, Fermion(conserve = N)), alternating)
                _, ψg = ground_state(hopping, ψf; algo = Dmrg(nsweeps = 2), limits = lim)
                for st in (ψf, ψg)
                    evolve(-im * hopping, 0.1, st; limits = lim)
                end
                expect2(ψg, (dag(C), C))
                lindblad = -im * hopping + Dissipator(N)(1)
                ρf = State{Mixed}(System(4, Fermion(conserve = strong(N))), alternating)
                ρf = evolve(lindblad, 0.1, ρf; limits = lim)
                ρf = evolve(lindblad + Dissipator(C)(1), 0.1, weaken(ρf); limits = lim)
                measure(ρf, Measure(N, Purity, Trace, EntanglementEntropy(1)); time = 0.)
                # a run with its checkpoints, and a run whose output is redirected
                mktempdir() do dir
                    cd(dir) do
                        phases = [CreateState{Mixed}(4, Qubit(), "Up"),
                                  Evolve(evolver = dissipative, duration = 0.2, time_step = 0.1,
                                         limits = lim, measurements = "data" => [Z, Purity])]
                        runTMS(SimData(; name = "workload", phases); checkpoint_interval = 1e-9)
                        runTMS(SimData(; name = "workload", phases); output = devnull)
                    end
                end
            end
        end
    end
end
