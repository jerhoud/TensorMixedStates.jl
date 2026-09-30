# How the tensor contractions of ITensors are threaded: threading_settings, which reads the
# settings of the process they depend on, and set_threading, which applies the dense or the
# block sparse mode, chosen by hand or from what a system conserves.

export set_threading, threading_settings

"""
    initial_threads

the threads of BLAS and of Strided when the package was loaded, the choices of Julia, of
Strided and of the user, which the `:dense` mode of `set_threading` puts back after `:blocks`
has set them to one.
"""
const initial_threads = Ref((blas = 1, strided = 1))

# Strided, a dependency, has set its threads by then
function __init__()
    initial_threads[] = (blas = BLAS.get_num_threads(),
                         strided = ITensors.NDTensors.Strided.get_num_threads())
end

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
    set_threading(mode)
    set_threading(system)
    set_threading(state)
    set_threading(sim)

set how ITensors threads the tensor contractions, in one of two modes, and return the settings
applied, see `threading_settings`:

- `:dense`: no block sparse multithreading, BLAS and Strided, which ITensors uses for the
  dense permutations, on the threads they had when the package was loaded, for the large
  products of dense matrices;
- `:blocks`: block sparse multithreading on the threads of Julia, BLAS and Strided on a single
  thread as ITensors recommends, for the many products of blocks of a system that conserves
  something. Julia must have been started with several threads, `julia --threads=N`.

Given a system, a state or a simulation, the mode is `:blocks` when the system conserves
something and Julia has several threads, `:dense` otherwise: which of the two is faster
depends on the calculation and on the BLAS library, see [Threads and performance](@ref).

These are settings of the process, which last until they are changed. The threads of Julia
and of its garbage collector are not among them, being fixed when Julia starts: see
[Threads and performance](@ref).

# Examples

    set_threading(:dense)
    set_threading(System(20, Fermion(conserve = N)))   # :blocks, given several Julia threads
"""
function set_threading(mode::Symbol)
    strided = ITensors.NDTensors.Strided
    if mode == :dense
        ITensors.enable_threaded_blocksparse(false)
        strided.set_num_threads(initial_threads[].strided)
        BLAS.set_num_threads(initial_threads[].blas)
    elseif mode == :blocks
        if Threads.nthreads() == 1
            error("the :blocks mode runs the blocks on the threads of Julia, which has a " *
                  "single one: start it with --threads=N")
        end
        # set before block sparse multithreading is enabled, which warns about both otherwise
        strided.set_num_threads(1)
        BLAS.set_num_threads(1)
        ITensors.enable_threaded_blocksparse(true)
    else
        error("unknown threading mode :$mode, expected :dense or :blocks")
    end
    return threading_settings()
end

set_threading(system::System) = set_threading(threading_mode(system))

set_threading(state::State) = set_threading(state.system)

set_threading(sim::Simulation) = set_threading(sim.state)

"""
    restore_threading(settings)

put back the settings `threading_settings` returned, but for the threads of Julia and of its
garbage collector, which cannot change.
"""
function restore_threading(s)
    # disabled first and enabled last: enabling warns while BLAS or Strided have several threads
    ITensors.enable_threaded_blocksparse(false)
    ITensors.NDTensors.Strided.set_num_threads(s.strided)
    BLAS.set_num_threads(s.blas)
    if s.blocksparse
        ITensors.enable_threaded_blocksparse(true)
    end
end
