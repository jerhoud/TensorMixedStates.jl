# The samples of a simulation: Realizations, the loop running each sample from the branch point,
# their generators, their averaged outputs and EachSample, and their checkpoints, resumption and
# prolongation.
#
# Goes here: src/runner/Samples.jl, the averages of measurements/Outputs.jl and what they change
# in the runner. The tests run in temporary
# directories, which they make the working directory while they run.

using TensorMixedStates: Realization

# the seeds of a sample written again from SplitMix64, apart from the code under test
function mix(x::UInt64)
    x += 0x9e3779b97f4a7c15
    x = (x ⊻ (x >> 30)) * 0xbf58476d1ce4e5b9
    x = (x ⊻ (x >> 27)) * 0x94d049bb133111eb
    return x ⊻ (x >> 31)
end
seeded(seed, sample, stream) =
    Xoshiro(mix(mix(mix(UInt64(seed)) ⊻ UInt64(sample)) ⊻ UInt64(stream)))

# a phase of a single call writing by itself between its steps
struct Kicks <: TensorMixedStates.AbstractPhase
    name::String
    time_start
    final_measurements
end

function TensorMixedStates.run_phase(p::Kicks, sim, done)
    for k in done + 1:3
        output(sim, "u" => "u" => c -> 1.0 * c.sample + k)
        sim = TensorMixedStates.step_done(sim, p, k)
    end
    return sim
end

# a chain of four qubits, its couplings drawn by the generator of the sample, evolved three steps
function quench(c; measurements = "steps" => [:step, X(1), Z])
    J = [1 + 0.5 * randn(c.rng) for _ in 1:3]
    H = -sum(J[i] * Z(i) * Z(i + 1) for i in 1:3) - sum(X(i) for i in 1:4)
    return Evolve(evolver = -im * H, duration = 0.3, time_step = 0.1,
                  limits = Limits(maxdim = 16), measurements = measurements)
end

sampled(name; samples = 3, phases = [CreateState{Pure}(4, Qubit(), "Up"), Realizations(quench)],
        final_measurements = "end" => [X(1), Z(2)]) =
    SimData(; name, seed = 7, samples, phases, final_measurements)

@testset "A sample runs the phases of its realization from the branch point" begin
    in_tmp() do
        draw = "draw" => c -> rand(c.rng)
        lists = [Data("end") => [X(1), Z(2)], EachSample(Data("each")) => [X(1), Z(2), draw]]
        sim = runTMS(sampled("s"; final_measurements = lists))
        rows = sim.data["each"]
        @test [ r["sample"] for r in rows ] == [1, 2, 3]
        for k in 1:3
            hand = runTMS(SimData(name = "h$k",
                                  phases = [CreateState{Pure}(4, Qubit(), "Up"),
                                            quench(Realization(k, seeded(7, k, 0));
                                                   measurements = [])],
                                  final_measurements = Data("end") => [X(1), Z(2)]))
            @test rows[k]["X(1)"] == hand.data["end"][1]["X(1)"]
            @test rows[k]["Z(2)"] == hand.data["end"][1]["Z(2)"]
            # the phases draw with the generator of the sample, the disorder with its own, both
            # seeded by the seed of the run and the sample
            @test rows[k]["draw"] == rand(seeded(7, k, 1))
        end
        @test rows[1]["X(1)"] ≠ rows[2]["X(1)"]
        # the averages of the samples
        x = [ r["X(1)"] for r in rows ]
        m = sum(x) / 3
        sd = sqrt(sum((x .- m) .^ 2) / 2)
        avg = only(sim.data["end"])
        @test avg["N"] == 3 && avg["N_eff"] ≈ 3
        @test avg["X(1)"] ≈ m rtol = 1e-14
        @test avg["std(X(1))"] ≈ sd rtol = 1e-10
        @test avg["err(X(1))"] ≈ sd / sqrt(3) rtol = 1e-10
        # a row for each step, the header first
        @test length(readlines("s/steps")) == 1 + 3
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
                             final_measurements = EachSample(Data("end")) => X(1)))
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
                             final_measurements = [Data("e") => :energy,
                                                   EachSample(Data("each")) => :energy]))
        energies = [ r["energy"] for r in sim.data["each"] ]
        @test length(unique(energies)) == 3
        # the final energies averaged whatever the sweeps each search made
        @test only(sim.data["e"])["energy"] ≈ sum(energies) / 3 rtol = 1e-14
        @test split(readlines("g/sweeps")[2], '\t')[2] == "3"
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
        @test length(readlines("p/end")) == 2
        @test split(readlines("p/end")[2], '\t')[2] == "2"
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
    @test_throws "has no samples" SimData(phases = [up], final_measurements = EachSample("f") => X)
    in_tmp() do
        @test_throws "both averaged and for each sample" runTMS(sampled("b";
            final_measurements = ["f" => X(1), EachSample("f") => X(1)]))
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
                           measurements = Data("steps") => [:step, Z(1), Purity])]
        end
        phases = [CreateState{Pure}(3, Qubit(), "Up"), CopyState(:state => :init),
                  Realizations(noisy)]
        sd(name) = SimData(; name, seed = 3, samples = 3, phases,
                           final_measurements = Data("end") => [Z(2), Fidelity(:init)])
        ref = runTMS(sd("ref"))
        withenv(() -> run_through(sd("chk"); max_time = -1), "TMS_SAMPLES" => "2")
        chk = run_through(sd("chk"); max_time = -1)
        for name in ("steps", "end")
            @test chk.data[name] == ref.data[name]
        end
        @test only(ref.data["end"])["N"] == 3
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

@testset "Averages: their precision, the rows before the samples, the times of the samples" begin
    in_tmp() do
        up = CreateState{Pure}(4, Qubit(), "Up")
        before = Gates(gates = X(1),
                       final_measurements = ["pre" => Z(1), EachSample("each") => Z(1)])
        # a deviation small against the mean, which the sums of squares would lose
        drift = "x" => c -> 10 + 1e-7 * c.sample
        sim = runTMS(SimData(name = "a", samples = 3,
                             phases = [up, before, Realizations(c -> Gates(gates = X(2)))],
                             final_measurements = Data("x") => drift))
        x = only(sim.data["x"])
        @test x["N"] == 3
        @test x["x"] ≈ 10 + 2e-7 rtol = 1e-15
        @test x["std(x)"] ≈ 1e-7 rtol = 1e-6
        @test x["err(x)"] ≈ 1e-7 / sqrt(3) rtol = 1e-6
        # a row before the samples, shared by all, with no error
        pre = readlines("a/pre")
        @test split(pre[1], '\t') == ["time", "N", "N_eff", "Z(1)", "err(Z(1))", "std(Z(1))"]
        @test split(pre[2], '\t')[[2, 4, 5, 6]] == ["3", "-1", "0", "0"]
        @test split(readlines("a/each")[1], '\t') == ["time", "sample", "Z(1)"]
        @test split(readlines("a/each")[2], '\t')[2] == "0"
        # samples whose rows have other times cannot be averaged
        steps(c) = Evolve(evolver = -im * X(1), duration = 0.2,
                          time_step = c.sample == 1 ? 0.1 : 0.05, measurements = "t" => Z(1))
        @test_throws "the same times" runTMS(SimData(name = "t", samples = 2,
                                                     phases = [up, Realizations(steps)]))
    end
end

@testset "Two simulations of the same seed draw the same disorder" begin
    in_tmp() do
        disorder(c) = [1 + 0.5 * randn(c.rng) for _ in 1:3]
        keep(J) = Update(c -> c.J = J; provides = [:J])
        a = runTMS(SimData(name = "a", seed = 5, samples = 3,
                           phases = [CreateState{Pure}(4, Qubit(), "Up"),
                                     Realizations(c -> [keep(disorder(c)), Gates(gates = X(1))])],
                           final_measurements = EachSample(Data("J")) => :J))
        # other phases, drawing with the generator of the sample, and fewer samples
        random = CreateState{Pure}(4, Qubit(), "Up"; randomize = 2)
        b = runTMS(SimData(name = "b", seed = 5, samples = 2,
                           phases = [Realizations(c -> [random, keep(disorder(c))])],
                           final_measurements = EachSample(Data("J")) => :J))
        Ja = [ r["J"] for r in a.data["J"] ]
        @test [ r["J"] for r in b.data["J"] ] == Ja[1:2]
        @test Ja[1] ≠ Ja[2]
    end
end

@testset "Averages of JSON, complex and vector values, rows some samples lack, one sample" begin
    in_tmp() do
        up = CreateState{Pure}(3, Qubit(), "Up")
        z = "z" => c -> complex(c.sample, -c.sample)
        v = "v" => c -> [1.0, c.sample]
        gate = Realizations(c -> Gates(gates = X(1)))
        sim = runTMS(SimData(name = "j", samples = 3, phases = [up, gate],
                             final_measurements = ["f.jsonl" => [z, v], Data("f") => [z, v]]))
        row = TensorMixedStates.JSON.parse(only(readlines("j/f.jsonl")))
        @test row["N"] == 3 && row["Re(z)"] ≈ 2 && row["Im(z)"] ≈ -2 && row["std(Re(z))"] ≈ 1
        @test row["v[1]"] == 1 && row["std(v[1])"] == 0 && row["err(v[2])"] ≈ 1 / sqrt(3)
        @test only(sim.data["f"])["v[2]"] ≈ 2
        # the first sample has fewer steps: the later rows average the samples that have them
        short(c) = Evolve(evolver = -im * (0.5 * c.sample * X(1)),
                          duration = c.sample == 1 ? 0.1 : 0.3, time_step = 0.1,
                          measurements = [Data("t") => Z(1), EachSample(Data("e")) => Z(1)])
        sim = runTMS(SimData(name = "s", samples = 3, phases = [up, Realizations(short)]))
        t, e = sim.data["t"], sim.data["e"]
        @test [ r["N"] for r in t ] == [3, 2, 2]
        @test t[3]["Z(1)"] ≈ sum(r["Z(1)"] for r in e if r["time"] == t[3]["time"]) / 2
        # one sample has no deviation
        one = runTMS(SimData(name = "o", samples = 1, phases = [up, gate],
                             final_measurements = Data("x") => X(1)))
        @test isnan(only(one.data["x"])["std(X(1))"])
    end
end

@testset "Averages refused: values, columns, lists, destinations" begin
    in_tmp() do
        up = CreateState{Pure}(3, Qubit(), "Up")
        gate = Realizations(c -> Gates(gates = X(1)))
        ending(name, final) = runTMS(SimData(; name, samples = 2, phases = [up, gate],
                                             final_measurements = final))
        sometimes = "m" => c -> c.sample == 1 ? missing : 1.0
        @test_throws "missing in some samples only" ending("m", "m" => sometimes)
        @test_throws "cannot be averaged" ending("s", "s" => "s" => c -> "text")
        @test_throws "other columns" ending("v", "v" => "v" => c -> ones(c.sample))
        @test_throws "two columns" ending("d", "d" => [Z, "Z(1)" => Z(1)])
        @test_throws "one the averages write" ending("n", "n" => "N" => X(1))
        # a list written by an Update, to an averaged text file holding another
        lists = Realizations(c -> Update(d -> (output(d, "g" => X(1)); output(d, "g" => Z(1)))))
        @test_throws "holds the list" runTMS(SimData(name = "l", samples = 2, phases = [up, lists]))
        # a destination both averaged and written for each sample, without a directory or by an
        # Update
        both = SimData(samples = 2, phases = [up, gate],
                       final_measurements = [Data("x") => X(1), EachSample(Data("x")) => X(1)])
        @test_throws "both averaged" runTMS(both; output = devnull)
        each = Realizations(c -> Update(d -> output(d, EachSample("g") => X(1))))
        @test_throws "both averaged" runTMS(SimData(name = "b", samples = 2, phases = [up, each],
                                                    final_measurements = "g" => X(1)))
        # EachSample written by an Update in a simulation without samples
        alone = Update(d -> output(d, EachSample("g") => X(1)))
        @test_throws "has no samples" runTMS(SimData(name = "a", phases = [up, alone]))
    end
end

@testset "A phase of one's own writing between its steps resumes as the uninterrupted one" begin
    in_tmp() do
        kicks = Realizations(c -> Kicks("kicks", nothing, []))
        sd(name) = SimData(; name, samples = 2,
                           phases = [CreateState{Pure}(2, Qubit(), "Up"), kicks])
        runTMS(sd("ref"))
        run_through(sd("chk"); max_time = -1)
        @test read("chk/u", String) == read("ref/u", String)
        @test length(readlines("ref/u")) == 4
    end
end

@testset "Averages keep the digits of a deviation small against the mean, and seeds stay" begin
    in_tmp() do
        x(k) = 1e8 + 1e-3 * k + 1e-4 * k^2
        sim = runTMS(SimData(name = "p", samples = 6,
                             phases = [CreateState{Pure}(2, Qubit(), "Up"),
                                       Realizations(c -> Gates(gates = X(1)))],
                             final_measurements = Data("x") => "x" => c -> x(c.sample)))
        xs = big.(x.(1:6))
        sd = sqrt(sum((xs .- sum(xs) / 6) .^ 2) / 5)
        @test only(sim.data["x"])["std(x)"] ≈ sd rtol = 1e-12
        # the first draw of the disorder of sample 1 of seed 7, the same in Julia 1.10 and 1.13
        @test rand(TensorMixedStates.disorder_rng(7, 1)) == 0.284828998012419
    end
end
