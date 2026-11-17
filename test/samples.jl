# The samples of a simulation: Realizations, the loop running each sample from the branch point,
# their generators, and their checkpoints, resumption and prolongation.
#
# Goes here: src/runner/Samples.jl and what it changes in the runner. The tests run in temporary
# directories, which they make the working directory while they run.

using TensorMixedStates: Realization, disorder_rng, sample_rng

# a chain of four qubits, its couplings drawn by the generator of the sample, evolved three steps
function quench(c; measurements = "steps" => [:sample, :step, X(1), Z])
    J = [1 + 0.5 * randn(c.rng) for _ in 1:3]
    H = -sum(J[i] * Z(i) * Z(i + 1) for i in 1:3) - sum(X(i) for i in 1:4)
    return Evolve(evolver = -im * H, duration = 0.3, time_step = 0.1,
                  limits = Limits(maxdim = 16), measurements = measurements)
end

sampled(name; samples = 3, phases = [CreateState{Pure}(4, Qubit(), "Up"), Realizations(quench)],
        final_measurements = "end" => [:sample, X(1), Z(2)]) =
    SimData(; name, seed = 7, samples, phases, final_measurements)

@testset "A sample runs the phases of its realization from the branch point" begin
    in_tmp() do
        draw = "draw" => c -> rand(c.rng)
        sim = runTMS(sampled("s"; final_measurements = Data("end") => [:sample, X(1), Z(2), draw]))
        rows = sim.data["end"]
        @test [ r["sample"] for r in rows ] == [1, 2, 3]
        for k in 1:3
            hand = runTMS(SimData(name = "h$k",
                                  phases = [CreateState{Pure}(4, Qubit(), "Up"),
                                            quench(Realization(k, disorder_rng(7, k));
                                                   measurements = [])],
                                  final_measurements = Data("end") => [X(1), Z(2)]))
            @test rows[k]["X(1)"] == hand.data["end"][1]["X(1)"]
            @test rows[k]["Z(2)"] == hand.data["end"][1]["Z(2)"]
            # the phases draw with the generator of the sample, the disorder with its own
            @test rows[k]["draw"] == rand(sample_rng(7, k))
        end
        @test rows[1]["X(1)"] ≠ rows[2]["X(1)"]
        @test length(readlines("s/steps")) == 1 + 3 * 3
        log = read("s/log", String)
        # the branch point once, the phases and steps of the first sample only
        @test count("Starting phase", log) == 2
        @test count("step 3, time", log) == 1
        @test occursin("Sample 3 of 3 done", log)
        @test isfile("s/branch.h5")
    end
end

@testset "Samples need no phase before the branch point" begin
    in_tmp() do
        sim = runTMS(sampled("s"; samples = 2,
                             phases = [Realizations(c -> [CreateState{Pure}(4, Qubit(), "Up"),
                                                          quench(c)])],
                             final_measurements = Data("end") => [:sample, X(1)]))
        @test [ r["sample"] for r in sim.data["end"] ] == [1, 2]
    end
end

@testset "The ground state of a disordered chain, sample by sample" begin
    in_tmp() do
        function search(c)
            J = [1 + 0.5 * randn(c.rng) for _ in 1:3]
            H = -sum(J[i] * Z(i) * Z(i + 1) for i in 1:3) - sum(X(i) for i in 1:4)
            return GroundState(hamiltonian = H, algo = Dmrg(nsweeps = 20, tol = 1e-10),
                               limits = Limits(maxdim = 16), measurements = "sweeps" => [:energy])
        end
        sim = runTMS(sampled("g"; phases = [CreateState{Pure}(4, Qubit(), "Up"),
                                            Realizations(search)],
                             final_measurements = Data("e") => [:sample, :energy]))
        energies = [ r["energy"] for r in sim.data["e"] ]
        @test length(unique(energies)) == 3
    end
end

@testset "Samples stopped anywhere resume as the uninterrupted ones" begin
    in_tmp() do
        ref = runTMS(sampled("ref"))
        sim = run_through(sampled("chk"); max_time = -1)
        for f in ("steps", "end")
            @test read("chk/$f", String) == read("ref/$f", String)
        end
        @test norm(sim.state - ref.state) < 1e-12
        log = read("chk/log", String)
        @test occursin("Resuming from the checkpoint: sample 3", log)
        # the branch point once
        @test count("Starting phase \"Creating state\"", log) == 1
    end
end

@testset "A finished run is prolonged by more samples" begin
    in_tmp() do
        runTMS(sampled("ref"; samples = 4))
        withenv(() -> runTMS(sampled("p")), "TMS_SAMPLES" => "2")
        @test length(readlines("p/end")) == 1 + 2
        withenv(() -> runTMS(sampled("p")), "TMS_SAMPLES" => "4")
        for f in ("steps", "end")
            @test read("p/$f", String) == read("ref/$f", String)
        end
        # fewer samples than those done only end the one the checkpoint is in
        withenv(() -> runTMS(sampled("p")), "TMS_SAMPLES" => "3")
        for f in ("steps", "end")
            @test read("p/$f", String) == read("ref/$f", String)
        end
        log = read("p/log", String)
        @test occursin("beyond the 3 asked for", log)
        @test count("Starting phase", log) == 2
        # the same number again does nothing more
        runTMS(sampled("p"; samples = 4))
        @test read("p/end", String) == read("ref/end", String)
        @test_throws "runTMS cannot read" withenv(() -> runTMS(sampled("p")), "TMS_SAMPLES" => "x")
        @test_throws "1 or more" withenv(() -> runTMS(sampled("p")), "TMS_SAMPLES" => "0")
    end
end

@testset "Samples refused" begin
    up = CreateState{Pure}(4, Qubit(), "Up")
    @test_throws "last phase" SimData(samples = 2, phases = [Realizations(quench), up])
    @test_throws "needs samples" SimData(phases = [up, Realizations(quench)])
    @test_throws "ends with the phase Realizations" SimData(samples = 2, phases = [up])
    @test_throws "1 or more" SimData(samples = 0, phases = [up, Realizations(quench)])
    @test_throws "at most typemax" SimData(samples = 2, seed = big(2)^70,
                                           phases = [up, Realizations(quench)])
    @test_throws "returns a Realizations" SimData(samples = 2, phases = [up,
                                                  Realizations(c -> Realizations(quench))])
    @test_throws "overwrite" SimData(samples = 2,
                                     phases = [up, Realizations(c -> SaveState(file = "f.h5"))])
    @test_throws "changes nothing" SimData(samples = 2, phases = [up, Realizations(c -> (c.x = 1))])
    @test_throws "not c.state" SimData(samples = 2, phases = [up, Realizations(c -> c.state)])
    @test_throws "only a simulation with samples" SimData(phases = [up],
                                                          final_measurements = "f" => [:sample])
    in_tmp() do
        other(c) = c.sample == 1 ? quench(c) : Gates(gates = X(1))
        sd = sampled("o"; phases = [up, Realizations(other)])
        @test_throws "phase 1 of sample 2 is Gates" runTMS(sd)
    end
end
