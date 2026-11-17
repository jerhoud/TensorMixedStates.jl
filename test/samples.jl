# The samples of a simulation: Realizations, the loop running each sample from the branch point,
# their generators, and their checkpoints, resumption and prolongation.
#
# Goes here: src/runner/Samples.jl and what it changes in the runner. The tests run in temporary
# directories, which they make the working directory while they run.

using TensorMixedStates: Realization

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
                                            quench(Realization(k, Xoshiro(UInt64[7, k]));
                                                   measurements = [])],
                                  final_measurements = Data("end") => [X(1), Z(2)]))
            @test rows[k]["X(1)"] == hand.data["end"][1]["X(1)"]
            @test rows[k]["Z(2)"] == hand.data["end"][1]["Z(2)"]
            # the phases draw with the generator of the sample, the disorder with its own, both
            # seeded by the seed of the run and the sample
            @test rows[k]["draw"] == rand(Xoshiro(UInt64[7, k, 1]))
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
        @test_throws "runTMS cannot read" withenv(() -> runTMS(sampled("p")), "TMS_SAMPLES" => "0")
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
    @test SimData(samples = Int32(2), phases = [up, Realizations(quench)]).samples === 2
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

@testset "Mixed samples, a state named before the branch point, Data, stops and prolongation" begin
    in_tmp() do
        function noisy(c)
            J = [1 + 0.5 * randn(c.rng) for _ in 1:2]
            H = -sum(J[i] * Z(i) * Z(i + 1) for i in 1:2) - sum(X(i) for i in 1:3)
            return [ToMixed(),
                    Evolve(evolver = -im * H + Dissipator(sqrt(0.3) * Sm)(1), duration = 0.2,
                           time_step = 0.1, algo = ApproxW(order = 2), limits = Limits(maxdim = 16),
                           measurements = Data("steps") => [:sample, :step, Z(1), Purity])]
        end
        phases = [CreateState{Pure}(3, Qubit(), "Up"), CopyState(:state => :init),
                  Realizations(noisy)]
        sd(name) = SimData(; name, seed = 3, samples = 3, phases,
                           final_measurements = Data("end") => [:sample, Z(2), Fidelity(:init)])
        ref = runTMS(sd("ref"))
        withenv(() -> run_through(sd("chk"); max_time = -1), "TMS_SAMPLES" => "2")
        chk = run_through(sd("chk"); max_time = -1)
        for name in ("steps", "end")
            @test chk.data[name] == ref.data[name]
        end
        @test length(ref.data["end"]) == 3
        @test norm(chk.state - ref.state) < 1e-12
    end
end

@testset "Errors, stops and checks of a sample" begin
    in_tmp() do
        up = CreateState{Pure}(4, Qubit(), "Up")
        ref = runTMS(sampled("ref"))
        # an error while Realizations draws names the sample
        failing(c) = c.sample == 2 ? error("no draw") : quench(c)
        @test_throws "no draw" runTMS(sampled("e"; phases = [up, Realizations(failing)]))
        @test occursin("in sample 2, outside its phases", read("e/error", String))
        # a stop there names where the checkpoint is, and the run resumes from it
        interrupted = Ref(true)
        function stopping(c)
            if c.sample == 2 && interrupted[]
                interrupted[] = false
                throw(InterruptException())
            end
            return quench(c)
        end
        sd = sampled("i"; phases = [up, Realizations(stopping)])
        @test stopped(runTMS(sd))
        @test occursin("Stopping before phase 2 of sample 2", read("i/log", String))
        @test !stopped(runTMS(sd))
        @test read("i/end", String) == read("ref/end", String)
        # the phases of each sample are checked before it runs
        reading(c) = Gates(gates = X(1),
                           final_measurements = Data("g") => [c.sample == 1 ? X(1) : :energy])
        @test_throws "reads energy" runTMS(sampled("r"; phases = [up, Realizations(reading)]))
        # a state saved by each sample under its own name
        saving(c) = [quench(c), SaveState(file = "states.h5", statename = "s$(c.sample)")]
        runTMS(sampled("v"; samples = 2, phases = [up, Realizations(saving)]))
        @test load_state("v/states.h5", "s1") isa State{Pure}
        @test load_state("v/states.h5", "s2") isa State{Pure}
        # the number of samples of the environment, for a simulation without any
        withenv(() -> runTMS(SimData(name = "n", phases = [up])), "TMS_SAMPLES" => "5")
        @test occursin("TMS_SAMPLES is not used", read("n/log", String))
    end
end
