# The samples of a simulation: Realizations, the loop running each sample from the branch point,
# their generators, their averaged outputs and EachSample, and their checkpoints, resumption and
# prolongation.
#
# Goes here: src/runner/Samples.jl, the averages of measurements/Outputs.jl and what they change
# in the runner. The tests run in temporary
# directories, which they make the working directory while they run.

using TensorMixedStates: Realization

# the seeds of a sample written again from SplitMix64, apart from the code under test
function mix64(x::UInt64)
    x += 0x9e3779b97f4a7c15
    x = (x ⊻ (x >> 30)) * 0xbf58476d1ce4e5b9
    x = (x ⊻ (x >> 27)) * 0x94d049bb133111eb
    return x ⊻ (x >> 31)
end
seeded(seed, sample, stream) =
    Xoshiro(mix64(mix64(mix64(UInt64(seed)) ⊻ UInt64(sample)) ⊻ UInt64(stream)))

# a phase of a single call writing by itself between its steps
struct WritingKicks <: TensorMixedStates.AbstractPhase
    name::String
    time_start
    final_measurements
end

function TensorMixedStates.run_phase(p::WritingKicks, sim, done)
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
        # the branch point once, in the directory of the simulation, then the phases and steps of
        # the first sample only, in the directory of its part
        @test count("Starting phase", read("s/log", String)) == 1
        log = read("s/parts/1/log", String)
        @test count("Starting phase", log) == 1
        @test count("step 3, time", log) == 1
        @test occursin("Sample 3 done", log)
        @test isfile("s/branch.h5") && isempty(readdir("s/pending"))
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
        lists = ["end" => [X(1), Z(2)], EachSample("each") => X(1)]
        ref = runTMS(sampled("ref"; final_measurements = lists))
        sim = run_through(sampled("chk"; final_measurements = lists); max_time = -1)
        for f in ("steps", "end", "each")
            @test read("chk/$f", String) == read("ref/$f", String)
        end
        @test norm(sim.state - ref.state) < 1e-12
        @test occursin("Resuming from the checkpoint: sample 3", read("chk/parts/1/log", String))
        # the branch point once
        @test count("Starting phase \"Creating state\"", read("chk/log", String)) == 1
        # the tolerance reached at the same sample, ±1 values reaching 0.6 at the fourth
        v = "v" => "v" => c -> isodd(c.sample) ? 1. : -1.
        tol(name) = SimData(; name, seed = 7, samples = 6, tol = "v" => 0.6, min_samples = 2,
                            phases = [CreateState{Pure}(4, Qubit(), "Up"), Realizations(quench)],
                            final_measurements = ["end" => [X(1), Z(2)], v])
        runTMS(tol("tref"))
        run_through(tol("tchk"); max_time = -1)
        @test TensorMixedStates.read_merged(abspath("tref")).first == 4
        @test TensorMixedStates.read_merged(abspath("tchk")).first == 4
        for f in ("steps", "end", "v")
            @test read("tchk/$f", String) == read("tref/$f", String)
        end
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
        # fewer samples than those done change nothing
        withenv(() -> runTMS(sampled("p")), "TMS_SAMPLES" => "3")
        for f in ("steps", "end")
            @test read("p/$f", String) == read("ref/$f", String)
        end
        @test count("Starting phase", read("p/log", String)) == 1
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
        @test occursin("in sample 2, outside its phases", read("e/parts/1/error", String))
        # the sample is noted failed: the other parts stop, and a new run tries it again
        failed = TensorMixedStates.read_queue(abspath("e")).failed
        @test collect(keys(failed)) == [2] && occursin("no draw", failed[2])
        # and told in the error file and the log of the simulation
        @test startswith(read("e/error", String), "ERROR in sample 2: ") &&
              occursin("no draw", read("e/error", String))
        @test occursin("ERROR in sample 2: ", read("e/log", String))
        @test_throws "no draw" runTMS(sampled("e"; phases = [up, Realizations(failing)]))
        # a run that ends removes the error file
        runTMS(sampled("e"))
        @test !isfile("e/error")
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
        @test occursin("Stopping before phase 2 of sample 2", read("i/parts/1/log", String))
        @test !stopped(runTMS(sd))
        @test read("i/end", String) == read("ref/end", String)
        # the phases of each sample are checked before it runs
        reading(c) = Gates(gates = X(1),
                           final_measurements = Data("g") => [c.sample == 1 ? X(1) : :energy])
        @test_throws "reads energy" runTMS(sampled("r"; phases = [up, Realizations(reading)]))
        # a state saved by each sample under its own name
        saving(c) = [quench(c), SaveState(file = "state$(c.sample).h5")]
        runTMS(sampled("v"; samples = 2, phases = [up, Realizations(saving)]))
        @test load_state("v/state1.h5", "state") isa State{Pure}
        @test load_state("v/state2.h5", "state") isa State{Pure}
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
        times = SimData(name = "t", samples = 3, phases = [up, Realizations(steps)])
        @test_throws "sample 2 cannot be merged" runTMS(times)
        # the next run merges first, and fails again before running another sample
        @test_throws "the same times" runTMS(times)
        @test readdir("t/pending") == ["sample-2"]
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
        @test_throws "one the averages or the samples write" ending("n", "n" => "N" => X(1))
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
        kicks = Realizations(c -> WritingKicks("kicks", nothing, []))
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

@testset "A part a process left is taken over, and the file stop of the simulation stops all" begin
    in_tmp() do
        ref = runTMS(sampled("ref"; samples = 4))
        # a part stopped in a sample, then its file running naming a process that has ended
        function stop_in_two(d)
            if d.sample == 2
                touch("s/stop")
            end
        end
        stopping = Realizations(c -> [Update(stop_in_two), quench(c)])
        once(name) = sampled(name; samples = 4,
                             phases = [CreateState{Pure}(4, Qubit(), "Up"), stopping])
        @test stopped(runTMS(once("s")))
        @test TensorMixedStates.read_queue(abspath("s")).current == Dict(1 => 2)
        q = run(`$(Base.julia_cmd()) --startup-file=no -e "exit()"`; wait = false)
        ended = getpid(q)
        wait(q)
        write("s/parts/1/running", "$(gethostname()) $ended 2026-11-18T12:00:00\n")
        sim = runTMS(once("s"))
        @test !stopped(sim)
        log = read("s/parts/1/log", String)
        @test occursin("process $ended on $(gethostname()), ended without removing running", log)
        @test occursin("Resuming from the checkpoint: sample 2", log)
        @test readdir("s/parts") == ["1"]
        for f in ("steps", "end")
            @test read("s/$f", String) == read("ref/$f", String)
        end
        # a live part is not taken: another one is made
        p = run(`$(Base.julia_cmd()) --startup-file=no -e "sleep(300)"`; wait = false)
        try
            write("s/parts/1/running", "$(gethostname()) $(getpid(p)) 2026-11-18T12:00:00\n")
            withenv(() -> runTMS(once("s")), "TMS_SAMPLES" => "5")
            @test readdir("s/parts") == ["1", "2"]
            @test TensorMixedStates.read_merged(abspath("s")).first == 5
            @test_throws "may still be going on" cleanTMS("s")
        finally
            kill(p, Base.SIGKILL)
        end
    end
end

@testset "Two processes run the samples of a simulation as one does" begin
    in_tmp() do
        program = """
            using TensorMixedStates, .Qubits
            # each sample waits for the second process, so that both run samples
            two = endswith(ARGS[1], "two")
            second = joinpath(ARGS[1], "parts", "2")
            waiting = Update(c -> (t = time(); while two && !isdir(second) &&
                                                         time() - t < 600; sleep(0.1); end);
                             provides = Symbol[])
            function quench(c)
                J = [1 + 0.5 * randn(c.rng) for _ in 1:3]
                H = -sum(J[i] * Z(i) * Z(i + 1) for i in 1:3) - sum(X(i) for i in 1:4)
                return [waiting, Evolve(evolver = -im * H, duration = 0.3, time_step = 0.1,
                                        limits = Limits(maxdim = 16),
                                        measurements = ["steps" => [X(1), Z],
                                                        EachSample("each") => Z(2)])]
            end
            # with a tolerance, ±1 values reaching 0.6 at the fourth sample
            tol = ARGS[1] in ("tone", "ttwo") ? "v" => 0.6 : []
            runTMS(SimData(name = ARGS[1], seed = 7, samples = 8, tol = tol, min_samples = 2,
                           phases = [CreateState{Pure}(4, Qubit(), "Up"), Realizations(quench)],
                           final_measurements = ["end.json" => [X(1), Z(2)],
                                                 EachSample("each.jsonl") => [X(1)],
                                                 "v" => "v" => c -> isodd(c.sample) ? 1. : -1.]))
            """
        write("prog.jl", program)
        julia = `$(Base.julia_cmd()) --project=$(Base.active_project()) --startup-file=no prog.jl`
        run(`$julia one`)
        ps = [ run(pipeline(`$julia two`; stdout = devnull, stderr = devnull); wait = false)
               for _ in 1:2 ]
        foreach(wait, ps)
        @test all(success, ps)
        @test readdir("two/parts") == ["1", "2"]
        @test all(occursin("Sample", read("two/parts/$p/log", String)) for p in 1:2)
        for f in ("steps", "each", "end.json", "each.jsonl")
            @test read("two/$f", String) == read("one/$f", String)
        end
        @test TensorMixedStates.read_merged(abspath("two")).first == 8
        @test isempty(readdir("two/pending"))
        # with a tolerance, the same samples merged whatever the processes
        run(`$julia tone`)
        ps = [ run(pipeline(`$julia ttwo`; stdout = devnull, stderr = devnull); wait = false)
               for _ in 1:2 ]
        foreach(wait, ps)
        @test all(success, ps)
        for f in ("steps", "each", "end.json", "each.jsonl", "v")
            @test read("ttwo/$f", String) == read("tone/$f", String)
        end
        @test TensorMixedStates.read_merged(abspath("ttwo")).first == 4
    end
end

@testset "A run without a directory gives the same Data, and a sample reads those merged" begin
    in_tmp() do
        look(d) = (d.seen = haskey(d.data, "x") ? Float64(only(d.data["x"])["N"]) : 0.)
        seeing = Realizations(c -> [Update(look; provides = [:seen]), quench(c)])
        sd(name) = sampled(name; phases = [CreateState{Pure}(4, Qubit(), "Up"), seeing],
                           final_measurements = [Data("x") => [X(1), Z(2)],
                                                 EachSample(Data("seen")) => :seen])
        a = runTMS(sd("a"))
        b = runTMS(sd("b"); output = devnull)
        @test [ r["seen"] for r in a.data["seen"] ] == [0, 1, 2]
        for name in ("x", "seen")
            @test a.data[name] == b.data[name]
        end
    end
end

@testset "Samples given back, a launch joining while waiting, an edited program" begin
    in_tmp() do
        up = CreateState{Pure}(4, Qubit(), "Up")
        ref = runTMS(sampled("ref"; phases = [up, Realizations(c -> [Update(d -> nothing),
                                                                     quench(c)])]))
        # a stop in sample 3, whose part is then given fewer samples
        once = Ref(true)
        function stop_in_three(d)
            if d.sample == 3 && once[]
                once[] = false
                touch("l/stop")
            end
        end
        sd = sampled("l"; phases = [up, Realizations(c -> [Update(stop_in_three), quench(c)])])
        @test stopped(runTMS(sd))
        withenv(() -> runTMS(sd), "TMS_SAMPLES" => "2")
        q = TensorMixedStates.read_queue(abspath("l"))
        @test q.free == [3] && isempty(q.current)
        @test TensorMixedStates.read_merged(abspath("l")).first == 2
        runTMS(sd)
        for f in ("steps", "end")
            @test read("l/$f", String) == read("ref/$f", String)
        end
        # the checkpoint of sample 3 resumed once, then removed when the sample was given back
        @test count("Resuming", read("l/parts/1/log", String)) == 1
        # an edited program, refused, then accepted and kept
        write("l/prog.jl", "# another program")
        @test_throws "written by another program" withenv(() -> runTMS(sd), "TMS_SAMPLES" => "4")
        withenv(() -> runTMS(sd), "TMS_SAMPLES" => "4", "TMS_ACCEPT_PROGRAM" => "1")
        @test occursin("TMS_ACCEPT_PROGRAM=1 allows", read("l/log", String)) &&
              occursin("the 3 samples merged before", read("l/log", String))
        withenv(() -> runTMS(sd), "TMS_SAMPLES" => "5")
        @test TensorMixedStates.read_merged(abspath("l")).first == 5
        # a launch waiting for the phases before the branch point, which another process runs
        p = run(`$(Base.julia_cmd()) --startup-file=no -e "sleep(300)"`; wait = false)
        try
            dir = abspath("w")
            mkpath(dir)
            write("w/running", "$(gethostname()) $(getpid(p)) 2026-11-18T12:00:00\n")
            launch(first) = TensorMixedStates.with_lock(dir) do
                TensorMixedStates.launch_samples!(dir, 2, nothing; first)
            end
            @test launch(true).first == :wait
            touch("w/stop")
            @test launch(true).first == :wait
            @test launch(false).first == :stop
            @test isfile("w/stop")
            touch("w/restart")
            @test_throws "holds the file restart" launch(true)
        finally
            kill(p, Base.SIGKILL)
        end
    end
end

@testset "The samples stop at the first whose average reaches the tolerance" begin
    in_tmp() do
        up = CreateState{Pure}(3, Qubit(), "Up")
        gate = Realizations(c -> Gates(gates = X(1)))
        value = "v" => c -> isodd(c.sample) ? 1. : -1.
        # the first number of samples, from `least` on, whose error is within `atol`
        function reaching(atol, least)
            for n in max(least, 2):10^4
                x = big.([ isodd(k) ? 1 : -1 for k in 1:n ])
                m = sum(x) / n
                if sqrt(sum((x .- m) .^ 2) / (n - 1) / n) ≤ atol
                    return n
                end
            end
        end
        sd(name, tol; min_samples = 10) =
            SimData(; name, samples = 400, tol, min_samples, phases = [up, gate],
                    final_measurements = [Data("v") => value, "f" => value])
        n = reaching(0.25, 10)
        @test n == 17
        sim = runTMS(sd("t", Data("v") => 0.25))
        @test only(sim.data["v"])["N"] == n
        @test TensorMixedStates.read_queue(abspath("t")).reached == n
        @test split(readlines("t/f")[2], '\t')[2] == string(n)
        @test only(runTMS(sd("d", Data("v") => 0.25); output = devnull).data["v"])["N"] == n
        @test startswith(TensorMixedStates.samples_line(abspath("t"), n),
                         "$n samples merged, within the tolerance")
        runTMS(SimData(; name = "n", samples = 12, tol = Data("v") => 1e-6,
                       phases = [up, gate], final_measurements = [Data("v") => value]))
        @test TensorMixedStates.samples_line(abspath("n"), 12; tolerance = true) ==
              "12 of 12 samples merged, the tolerance not reached"
        # every row holds min_samples: one written by two samples never reaches it
        w = Realizations(c -> Gates(gates = X(1), final_measurements =
                                        c.sample ≤ 2 ? [Data("w") => "w" => c -> 1.] : []))
        sim = runTMS(SimData(; name = "w", samples = 30, tol = Data("w") => 0.1,
                             phases = [up, w], final_measurements = [Data("v") => value]))
        @test only(sim.data["v"])["N"] == 30
        # more samples change nothing, a tighter tolerance goes on
        withenv(() -> runTMS(sd("t", Data("v") => 0.25)), "TMS_SAMPLES" => "500")
        @test TensorMixedStates.read_merged(abspath("t")).first == n
        sim = runTMS(sd("t", "f" => (atol = 0.1,)))
        @test only(sim.data["v"])["N"] == reaching(0.1, 10)
        # at least min_samples, and a file
        sim = runTMS(sd("m", "f" => (atol = 0.25, rtol = 0.); min_samples = 30))
        @test only(sim.data["v"])["N"] == reaching(0.25, 30) == 30
        # refused
        @test_throws "has no samples" SimData(phases = [up], tol = "f" => 1e-3,
                                              final_measurements = "f" => X(1))
        @test_throws "which no list of measurements writes" sd("x", "g" => 0.1)
        @test_throws "of zero or above" sd("x", "f" => -1)
        @test_throws "of zero or above" sd("x", "f" => NaN)
        @test_throws "which is not averaged" sd("x", EachSample("f") => 0.1)
        @test_throws "which is not averaged" sd("x", "" => 0.1)
        @test_throws "a tolerance is" sd("x", "f" => (atol = "a",))
        @test_throws "min_samples is 1 or more" sd("x", "f" => 0.1; min_samples = 0)
        @test_throws "names a Data(name) alone" runTMS(sd("x", "f" => 0.1); output = devnull)
        # a sample failed beyond the one reaching the tolerance stops no part
        T = TensorMixedStates
        dir = mkpath(abspath("q"))
        T.write_queue(dir, T.Queue(40, 17, 19, Int[], Dict(2 => 18), Dict(18 => "boom")))
        @test isnothing(T.take_sample!(dir, 1))
        T.write_queue(dir, T.Queue(40, nothing, 19, Int[], Dict{Int, Int}(), Dict(18 => "boom")))
        @test_throws "sample 18 failed" T.take_sample!(dir, 1)
    end
end

@testset "The files of the samples name what wrote them, refused where this one cannot read" begin
    in_tmp() do
        T = TensorMixedStates
        sd = SimData(; name = "v", samples = 3,
                     phases = [CreateState{Pure}(3, Qubit(), "Up"),
                               Realizations(c -> Gates(gates = X(1)))],
                     final_measurements = "f" => "x" => c -> 1. * c.sample)
        runTMS(sd)
        file = abspath("v/merged.jls")
        kept = read(file)
        header = T.written_by()
        @test String(kept[1:ncodeunits(header)]) == header
        body = kept[ncodeunits(header) + 1:end]
        write(file, [codeunits(replace(header, "Julia $VERSION" => "Julia 1.9.0")); body])
        @test_throws "go on with Julia 1.9" runTMS(sd)
        write(file, [codeunits(replace(header, "TensorMixedStates $(T.samples_file_version)" =>
                                               "TensorMixedStates 0")); body])
        @test_throws "the layout of another version" runTMS(sd)
        write(file, kept)
        # another version of TensorMixedStates goes on, noted in the log
        T.HDF5.h5open(abspath("v/branch.h5"), "r+") do f
            T.HDF5.delete_object(f, "written_by")
            f["written_by"] = replace(header, r"version .*" => "version 2.1.0")
        end
        withenv(() -> runTMS(sd), "TMS_SAMPLES" => "4")
        @test occursin("those before by 2.1.0", read("v/log", String))
        @test T.read_merged(abspath("v")).first == 4
    end
end

@testset "The names the samples reserve are free in a simulation without samples" begin
    in_tmp() do
        up = CreateState{Pure}(2, Qubit(), "Up")
        give = Update(c -> (c.sample = 7); provides = [:sample])
        runTMS(SimData(; name = "a", phases = [up, give],
                       final_measurements = ["queue.json" => :sample, "f" => "weight" => X(1)]))
        @test isfile("a/queue.json") && isfile("a/f")
        r = Realizations(c -> Gates(gates = X(1)))
        @test_throws "a name the run reserves" SimData(; name = "b", samples = 2,
                                                         phases = [up, give, r])
        sd(name, final) = SimData(; name, samples = 2, phases = [up, r], final_measurements = final)
        @test_throws "a file of the simulation directory" runTMS(sd("b", "queue.json" => X(1)))
        for (k, list) in enumerate((EachSample("e") => "weight" => X(1),
                                    EachSample("e") => "sample" => X(1),
                                    EachSample("e") => [X(1), :sample],
                                    "f" => "trajectory" => X(1)))
            @test_throws "the samples write" runTMS(sd("c$k", list))
        end
        runTMS(sd("d", EachSample("e") => [:sample, X(1)]))
        @test split(readlines("d/e")[1], '\t') == ["time", "sample", "X(1)"]
    end
end

@testset "The lock of a simulation is released on an error and taken as soon as it is free" begin
    in_tmp() do
        T = TensorMixedStates
        dir = mkpath(abspath("l"))
        @test_throws "inside" T.with_lock(() -> error("inside"), dir)
        t = time()
        @test T.with_lock(() -> 1, dir) == 1
        @test time() - t < 1
        released = Ref(0.)
        holder = @async T.with_lock(dir) do
            sleep(0.5)
            released[] = time()
        end
        sleep(0.1)
        got = T.with_lock(time, dir)
        wait(holder)
        @test released[] ≤ got < released[] + 1
        # the file running read once, absent when its process removed it
        @test isnothing(T.running_marker(dir)) && !T.live_run(dir)
        write(joinpath(dir, "running"), "host 1 date\n")
        @test T.running_marker(dir) == "host 1 date\n"
    end
end

@testset "The lists to one EachSample are gathered behind one column sample" begin
    in_tmp() do
        up = CreateState{Pure}(2, Qubit(), "Up")
        r = Realizations(c -> Gates(gates = X(1)))
        sd(name, final) = SimData(; name, samples = 2, phases = [up, r], final_measurements = final)
        runTMS(sd("a", [EachSample("e") => X(1), EachSample("e") => Z(1)]))
        @test split(readlines("a/e")[1], '\t') == ["time", "sample", "X(1)", "Z(1)"]
        @test_throws "put :sample first" runTMS(sd("b", EachSample("e") => [X(1), :sample]))
        # a value sample given by an Update does not make samples
        give = Update(c -> (c.sample = 7); provides = [:sample])
        each = Update(c -> output(c, EachSample("e") => X(1)))
        @test_throws "has no samples" runTMS(SimData(; name = "c", phases = [up, give, each]))
    end
end
