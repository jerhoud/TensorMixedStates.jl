# How the tensor contractions of ITensors are threaded: threading_settings reads it,
# ThreadingState holds what a program can change of it, set_threading applies it.

export set_threading, threading_settings, ThreadingState

"""
    threading_settings()

how the tensor contractions are threaded, as a named tuple:

- `julia`: the threads of Julia, on which block sparse multithreading runs
- `gc`: the threads of the garbage collector
- `blas`, `blas_library`: the threads of BLAS, which runs the products of dense matrices, and
  its library, OpenBLAS unless MKL was loaded
- `strided`: the threads of Strided, which ITensors uses for the dense permutations
- `blocksparse`: whether ITensors multithreads the contractions of block sparse tensors

`julia` and `gc` are fixed when Julia starts, by `julia --threads=N --gcthreads=M`; the others
can be changed at any time, see `ThreadingState` and `set_threading`.
"""
threading_settings() =
    (julia = Threads.nthreads(), gc = Threads.ngcthreads(), blas = BLAS.get_num_threads(),
     blas_library = join((basename(l.libname) for l in BLAS.get_config().loaded_libs), ", "),
     strided = ITensors.NDTensors.Strided.get_num_threads(),
     blocksparse = ITensors.using_threaded_blocksparse())

"""
    threading_mode(system)

the mode of `set_threading` suited to `system`: `:blocks` when it conserves something and
Julia has several threads, `:dense` otherwise
"""
threading_mode(system::System) =
    isempty(symmetries(system).names) || Threads.nthreads() == 1 ? :dense : :blocks

"""
    ThreadingState(; blas, strided, blocksparse)

how the tensor contractions are threaded, in what a program can change, see
`threading_settings`:

- `blas`: the threads of BLAS, which runs the products of dense matrices
- `strided`: the threads of Strided, which ITensors uses for the dense permutations
- `blocksparse`: whether ITensors multithreads the contractions of block sparse tensors

A field left out takes the value in force, so that `ThreadingState()` is the threading in
force. `set_threading` applies one, returning the one it replaces.

# Examples

    ThreadingState()                          # the threading in force
    set_threading(ThreadingState(blas = 2))   # BLAS on two threads, the rest as it is
"""
struct ThreadingState
    blas::Int
    strided::Int
    blocksparse::Bool
end

ThreadingState(; blas = BLAS.get_num_threads(),
               strided = ITensors.NDTensors.Strided.get_num_threads(),
               blocksparse = ITensors.using_threaded_blocksparse()) =
    ThreadingState(blas, strided, blocksparse)

show(io::IO, s::ThreadingState) =
    print(io, "ThreadingState(blas = ", s.blas, ", strided = ", s.strided,
          ", blocksparse = ", s.blocksparse, ")")

"""
    default_blas_threads()

the threads Julia gives BLAS when it starts: those an environment variable of OpenBLAS asks
for, or else one for every two threads of the processor, one for each on an Apple processor
"""
function default_blas_threads()
    for name in ("OPENBLAS_NUM_THREADS", "GOTO_NUM_THREADS", "OMP_NUM_THREADS")
        n = tryparse(Int, get(ENV, name, ""))
        if !isnothing(n)
            return max(1, n)
        end
    end
    # the threads given to the process, which Julia 1.10 does not count
    cpus = isdefined(Sys, :EFFECTIVE_CPU_THREADS) ? Sys.EFFECTIVE_CPU_THREADS : Sys.CPU_THREADS
    return Sys.isapple() && Sys.ARCH === :aarch64 ? max(1, cpus) : max(1, cpus ÷ 2)
end

"""
    set_threading(threading::ThreadingState)
    set_threading(mode)
    set_threading(system)
    set_threading(state)
    set_threading(sim)

set how ITensors threads the tensor contractions, and return the `ThreadingState` it
replaces. The threading is given as a `ThreadingState`, or as one of two modes:

- `:dense`: no block sparse multithreading, BLAS on its threads, Strided on a single one, as
  ITensors recommends: started with Julia, Strided has as many threads as Julia, which
  compete with those of BLAS and can slow down dense contractions considerably. BLAS keeps its
  threads, except after `:blocks`, where it gets back those Julia starts it with. `SimData`
  applies this mode by default, and a program calling the functions of TMS directly should
  start with `set_threading(:dense)`;
- `:blocks`: block sparse multithreading on the threads of Julia, BLAS and Strided on a single
  thread, for a system that conserves something. Julia has to be started with several
  threads, `julia -t N`: on a single one, everything runs on one core, and ITensors warns.

Given a system, a state or a simulation, the mode is `:blocks` when the system conserves
something and Julia has several threads, `:dense` otherwise; which is faster depends on the
calculation, see [Threads and performance](@ref). These settings of the process last until
they are changed; the threads of Julia and of its garbage collector, fixed when Julia starts,
are not among them.

# Examples

    set_threading(:dense)
    old = set_threading(System(20, Fermion(conserve = N)))   # :blocks, given several threads
    set_threading(old)                                       # back to what it was
    set_threading(ThreadingState(blas = 2))
"""
function set_threading(threading::ThreadingState)
    old = ThreadingState()
    if !threading.blocksparse
        ITensors.enable_threaded_blocksparse(false)
    end
    ITensors.NDTensors.Strided.set_num_threads(threading.strided)
    BLAS.set_num_threads(threading.blas)
    # enabled once BLAS and Strided are set: switching it on warns while they have several threads
    if threading.blocksparse
        ITensors.enable_threaded_blocksparse(true)
    end
    return old
end

function set_threading(mode::Symbol)
    if mode ∉ (:dense, :blocks)
        error("unknown threading mode :$mode, expected :dense or :blocks")
    end
    blas = mode == :blocks ? 1 :
           ITensors.using_threaded_blocksparse() ? default_blas_threads() : BLAS.get_num_threads()
    return set_threading(ThreadingState(; blas, strided = 1, blocksparse = mode == :blocks))
end

set_threading(system::System) = set_threading(threading_mode(system))

set_threading(state::AbstractState) = set_threading(state.system)

set_threading(sim::Simulation) = set_threading(sim.state)
