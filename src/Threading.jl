# How the tensor contractions of ITensors are threaded: threading_settings, which reads the
# settings of the process they depend on, and set_threading, which applies the dense or the
# block sparse mode, chosen by hand or from what a system conserves.

export set_threading, threading_settings

"""
    blas_before_blocks

the threads BLAS had before the `:blocks` mode of `set_threading` put it on a single one,
which `:dense` gives back, and `nothing` outside that mode.
"""
const blas_before_blocks = Ref{Union{Nothing, Int}}(nothing)

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
can be changed at any time, see `set_threading`.
"""
threading_settings() =
    (julia = Threads.nthreads(), gc = Threads.ngcthreads(), blas = BLAS.get_num_threads(),
     blas_library = join((basename(l.libname) for l in BLAS.get_config().loaded_libs), ", "),
     strided = ITensors.NDTensors.Strided.get_num_threads(),
     blocksparse = ITensors.using_threaded_blocksparse())

"""
    threading_mode(system)

the mode of `set_threading` suited to `system`: `:blocks` when it conserves something, its
tensors being then block sparse, and Julia has several threads to run the blocks on, `:dense`
otherwise.
"""
threading_mode(system::System) =
    isempty(symmetries(system).names) || Threads.nthreads() == 1 ? :dense : :blocks

"""
    ThreadingState

the threading of the process as `set_threading` found it, which `set_threading` takes back to
restore it: the settings `threading_settings` returns, and the threads `:dense` would give
back to BLAS.
"""
struct ThreadingState
    settings::NamedTuple
    blas_before_blocks::Union{Nothing, Int}
end

show(io::IO, s::ThreadingState) = print(io, "ThreadingState", s.settings)

"""
    save_threading()

the threading of the process now, as a `ThreadingState`.
"""
save_threading() = ThreadingState(threading_settings(), blas_before_blocks[])

"""
    set_threading(mode)
    set_threading(system)
    set_threading(state)
    set_threading(sim)
    set_threading(old)

set how ITensors threads the tensor contractions, in one of two modes, and return the threading
in force before, which `set_threading(old)` puts back:

- `:dense`: no block sparse multithreading, BLAS running each product of dense matrices on its
  threads, and Strided, which ITensors uses for the dense permutations, on a single one, as
  ITensors recommends. Started with Julia, Strided has as many threads as Julia has, which
  compete with those of BLAS and can slow down a calculation on dense tensors considerably.
  `SimData` applies this mode by default, and a program calling the functions of TMS directly
  should start with `set_threading(:dense)`;
- `:blocks`: block sparse multithreading on the threads of Julia, BLAS and Strided on a single
  thread, for the many products of blocks of a system that conserves something. Julia has to
  be started with several threads, `julia --threads=N` or `julia -t N`: on a single one,
  everything runs on one core, and ITensors warns. `:dense` gives BLAS its threads back.

Given a system, a state or a simulation, the mode is `:blocks` when the system conserves
something and Julia has several threads, `:dense` otherwise: which of the two is faster
depends on the calculation and on the BLAS library, see [Threads and performance](@ref).
These are settings of the process, which last until they are changed; the threads of Julia and
of its garbage collector are not among them, being fixed when Julia starts. The settings now
in force are read by `threading_settings`.

# Examples

    set_threading(:dense)
    old = set_threading(System(20, Fermion(conserve = N)))   # :blocks, given several threads
    set_threading(old)                                       # back to what it was
"""
function set_threading(mode::Symbol)
    if mode ∉ (:dense, :blocks)
        error("unknown threading mode :$mode, expected :dense or :blocks")
    end
    old = save_threading()
    ITensors.NDTensors.Strided.set_num_threads(1)
    if mode == :dense
        ITensors.enable_threaded_blocksparse(false)
        if !isnothing(blas_before_blocks[])
            BLAS.set_num_threads(blas_before_blocks[])
            blas_before_blocks[] = nothing
        end
    else
        if isnothing(blas_before_blocks[])
            blas_before_blocks[] = BLAS.get_num_threads()
        end
        # before block sparse multithreading is enabled, which warns about it otherwise
        BLAS.set_num_threads(1)
        ITensors.enable_threaded_blocksparse(true)
    end
    return old
end

function set_threading(old::ThreadingState)
    current = save_threading()
    s = old.settings
    # disabled first and enabled last: enabling warns while BLAS or Strided have several threads
    ITensors.enable_threaded_blocksparse(false)
    ITensors.NDTensors.Strided.set_num_threads(s.strided)
    BLAS.set_num_threads(s.blas)
    blas_before_blocks[] = old.blas_before_blocks
    if s.blocksparse
        ITensors.enable_threaded_blocksparse(true)
    end
    return current
end

set_threading(system::System) = set_threading(threading_mode(system))

set_threading(state::State) = set_threading(state.system)

set_threading(sim::Simulation) = set_threading(sim.state)
