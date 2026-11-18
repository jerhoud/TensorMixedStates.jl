# The parts of a simulation with samples, which processes run at once, each in a directory of
# its own: the lock of the directory of the simulation, the queue of its samples, and the merge of
# the samples ended into its outputs, in the order of the samples.

"""
    with_lock(f, dir)

`f()` holding the lock of the directory `dir`, its file `lock`, which a process holds a short
while: one left by a process that died is passed over after a minute, five for one of another
machine, see `mkpidlock`
"""
with_lock(f::Function, dir::String) =
    mkpidlock(f, joinpath(dir, "lock"); stale_age = 60, poll_interval = 1)

"""
    part_dir(dir, part)

the directory of the part numbered `part` of the simulation of directory `dir`
"""
part_dir(dir::String, part::Int) = joinpath(dir, "parts", string(part))

"""
    Queue(samples, next, current, failed)

the queue of the samples of a simulation, its file `queue.json`, read and written under its
lock: the samples asked for, the next one no part has taken, the sample each part runs, by the
number of the part, and the samples failed, with the part whose file `error` says why
"""
mutable struct Queue
    samples::Int
    next::Int
    current::Dict{Int, Int}
    failed::Dict{Int, Int}
end

queue_file(dir::String) = joinpath(dir, "queue.json")

function read_queue(dir::String)
    q = JSON.parsefile(queue_file(dir))
    return Queue(q["samples"], q["next"],
                 Dict{Int, Int}(parse(Int, p) => k for (p, k) in q["current"]),
                 Dict{Int, Int}(parse(Int, k) => p for (k, p) in q["failed"]))
end

function write_queue(dir::String, q::Queue)
    tmp = queue_file(dir) * ".tmp"
    open(tmp, "w") do f
        JSON.print(f, Dict("samples" => q.samples, "next" => q.next,
                           "current" => Dict(string(p) => k for (p, k) in q.current),
                           "failed" => Dict(string(k) => p for (k, p) in q.failed)))
    end
    sync_file(tmp)
    Base.Filesystem.rename(tmp, queue_file(dir))
    return nothing
end

"""
    take_sample!(dir, part)

the sample the part `part` of the simulation of directory `dir` runs next: the one it holds, or
else the next one no part has taken, `nothing` when none is left; refused once a sample failed
"""
function take_sample!(dir::String, part::Int)
    return with_lock(dir) do
        q = read_queue(dir)
        if haskey(q.current, part)
            return q.current[part]
        elseif !isempty(q.failed)
            k = minimum(keys(q.failed))
            error("sample $k failed, as the file error of $(part_dir(dir, q.failed[k])) says: " *
                  "the samples stop")
        elseif q.next > q.samples
            return nothing
        end
        k = q.next
        q.next += 1
        q.current[part] = k
        write_queue(dir, q)
        return k
    end
end

"""
    sample_failed!(dir, sample, part)

note in the queue of the simulation of directory `dir` that `sample` failed, the file `error` of
the part `part` saying why, so that the other parts stop, see `take_sample!`
"""
function sample_failed!(dir::String, sample::Int, part::Int)
    with_lock(dir) do
        q = read_queue(dir)
        q.failed[sample] = part
        write_queue(dir, q)
    end
    return nothing
end

"""
    sample_ended!(dir, part, format)

note in the queue of the simulation of directory `dir` that the part `part` ended its sample,
whose rows are written, see `write_pending`, and merge the samples ended, see `merge_samples!`;
the first sample merged and the last, `nothing` when none was
"""
function sample_ended!(dir::String, part::Int, format::String)
    return with_lock(dir) do
        q = read_queue(dir)
        delete!(q.current, part)
        write_queue(dir, q)
        return merge_samples!(dir, part, format)
    end
end

"""
    merge_samples!(dir, part, format)

under the lock of the directory `dir` of a simulation, its samples ended added to its outputs in
their order, see `replay!`, from the first not merged up to the first not ended: its files
written, then `merged.jls`, then the rows of the samples merged removed, so that a crash at any
point leaves a merge to go on with. It returns the first sample merged and the last, `nothing`
when none was. A sample that cannot be merged is noted failed, the file `error` of the part
`part`, which merges, saying why, see `sample_failed!`.
"""
function merge_samples!(dir::String, part::Int, format::String)
    m, marks = read_merged(dir)
    k = m + 1
    if !isfile(pending_file(dir, k))
        return nothing
    end
    o = reopen!(Outputs(; dir, format, sampled = true), marks)
    try
        while isfile(pending_file(dir, k))
            replay!(o, read_pending(dir, k))
            k += 1
        end
    catch
        q = read_queue(dir)
        q.failed[k] = part
        write_queue(dir, q)
        # the averaged files left as the merge before wrote them, the rows written to the others
        # cut by the next merge, which reopens them at their marks
        foreach(close, (h for h in values(o.opened) if !(h isa AveragedOpened)))
        rethrow()
    end
    marks = resume_marks(o)
    close!(o)
    write_merged(dir, k - 1, marks)
    # those a crash left after the merge before as well
    for f in readdir(joinpath(dir, "pending"))
        n = tryparse(Int, chopprefix(f, "sample-"))
        if startswith(f, "sample-") && !isnothing(n) && n < k
            rm(joinpath(dir, "pending", f); force = true)
        end
    end
    return m + 1 => k - 1
end
