# The precompilation workload, small computations run when the package is precompiled, so that a
# first script does not wait for them to compile; set TMS_SKIP_PRECOMPILE_WORKLOAD=true to skip it
# during development, never for a release.

using PrecompileTools: @setup_workload, @compile_workload

if get(ENV, "TMS_SKIP_PRECOMPILE_WORKLOAD", "false") ≠ "true"
    @setup_workload begin
        system = System(3, Qubit())
        h = -Z(1) * Z(2) - Z(2) * Z(3) - 0.5 * (X(1) + X(2) + X(3))
        lim = Limits(maxdim = 4)
        @compile_workload begin
            # silent: the warnings of a computation, and the line runTMS ends with
            Logging.with_logger(Logging.NullLogger()) do
                # pure and mixed states, of real tensors and of complex ones, which take paths of
                # their own
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
                evolve(dissipative, 0.1, ρup; limits = lim)
                evolve(dissipative, 0.1, ρ; algo = ApproxW(order = 2), limits = lim)
                ground_state(h, ψ; algo = Dmrg(nsweeps = 2), limits = lim)
                steady_state(dissipative, ρup; algo = Dmrg(nsweeps = 1), limits = lim)
                thermalize(h, 0.1, State{Mixed}(system, "FullyMixed"); limits = lim)
                apply(X(1) * Z(2), ψ)
                apply(Gate(X)(1), ρ)
                # a run with its checkpoints, and a run whose output is redirected
                mktempdir() do dir
                    cd(dir) do
                        phases = [CreateState{Mixed}(3, Qubit(), "Up"),
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
