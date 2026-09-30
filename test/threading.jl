# The threading of the tensor contractions.
#
# Goes here: set_threading, ThreadingState and threading_settings, and the threading field of
# SimData. These are settings of the process: every test puts back the threading it found, so
# that the groups after it run as they would have.

"""
    quietly(f)

`f()`, with what it prints thrown away: ITensors prints a warning when block sparse
multithreading is switched on while Julia has a single thread.
"""
quietly(f) = redirect_stdout(f, devnull)

@testset "ThreadingState" begin
    now = threading_settings()
    s = ThreadingState()
    @test (s.blas, s.strided, s.blocksparse) == (now.blas, now.strided, now.blocksparse)
    # a field left out takes the value in force
    @test ThreadingState(blas = 3) == ThreadingState(3, s.strided, s.blocksparse)
    # printed as the call that builds it
    t = ThreadingState(blas = 2, strided = 1, blocksparse = true)
    @test repr(t) == "ThreadingState(blas = 2, strided = 1, blocksparse = true)"
    @test eval(Meta.parse(repr(t))) == t
end

@testset "set_threading" begin
    start = ThreadingState()
    try
        # every change returns the threading it replaces, which set_threading takes back
        @test set_threading(:dense) == start
        dense = ThreadingState()
        @test dense == ThreadingState(blas = start.blas, strided = 1, blocksparse = false)
        @test quietly(() -> set_threading(:blocks)) == dense
        blocks = ThreadingState(blas = 1, strided = 1, blocksparse = true)
        @test ThreadingState() == blocks
        @test set_threading(dense) == blocks
        @test ThreadingState() == dense

        # :dense leaves BLAS as it is, but gives it the threads Julia starts it with when it
        # comes after :blocks, which put it on one
        set_threading(ThreadingState(blas = 2))
        set_threading(:dense)
        @test threading_settings().blas == 2
        quietly(() -> set_threading(:blocks))
        set_threading(:dense)
        @test threading_settings().blas == TensorMixedStates.default_blas_threads()

        # the mode suited to a system, a state or a simulation: blocks only when something is
        # conserved and Julia has threads to run them on
        set_threading(System(2, Qubit()))
        @test !ThreadingState().blocksparse
        set_threading(State{Pure}(System(2, Fermion(conserve = N)), "Occ"))
        @test ThreadingState().blocksparse == (Threads.nthreads() > 1)
        set_threading(Simulation(State{Pure}(System(2, Qubit()), "Up"); output = devnull))
        @test !ThreadingState().blocksparse

        @test_throws ErrorException set_threading(:blocksparse)
    finally
        quietly(() -> set_threading(start))
    end
end

@testset "The threading of a simulation" begin
    start = ThreadingState()
    phases = [CreateState{Pure}(2, Qubit(), "Up")]
    @test SimData(; phases).threading == :dense
    @test_throws ErrorException SimData(; phases, threading = :blocksparse)
    try
        mktempdir() do dir
            cd(dir) do
                # a mode that does not depend on the state is set from the start, which the
                # stamp records, and the threading found is put back at the end
                runTMS(SimData(name = "dense", phases = phases))
                stamp = readlines(joinpath("dense", "stamp"))
                @test "Threading :dense" in stamp
                @test "Strided threads 1" in stamp
                @test ThreadingState() == start

                runTMS(SimData(name = "left", threading = nothing, phases = phases))
                @test "Threading nothing" in readlines(joinpath("left", "stamp"))
                @test ThreadingState() == start

                # :auto chooses before each phase from the system of the state, and logs each
                # change: blocks once the state conserving something exists, dense once the
                # conserved quantity is dropped
                hopping = sum(dag(C)(i)C(i + 1) + dag(C)(i + 1)C(i) for i in 1:3)
                evolve = Evolve(duration = 0.1, time_step = 0.1, algo = Tdvp(), evolver = -im * hopping)
                auto = [CreateState{Mixed}(4, Fermion(conserve = N), ["Occ", "Emp", "Occ", "Emp"]),
                        evolve, Weaken(target = ()), evolve]
                runTMS(SimData(name = "auto", threading = :auto, phases = auto))
                changes = filter(l -> startswith(l, "Threading set to"),
                                 readlines(joinpath("auto", "log")))
                if Threads.nthreads() > 1
                    @test length(changes) == 2
                    @test startswith(changes[1], "Threading set to :blocks")
                    @test startswith(changes[2], "Threading set to :dense")
                else
                    @test isempty(changes)
                end
                @test ThreadingState() == start
            end
        end
    finally
        quietly(() -> set_threading(start))
    end
end
