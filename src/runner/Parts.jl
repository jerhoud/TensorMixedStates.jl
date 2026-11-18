# The parts of a simulation with samples, which processes run at once, each in a directory of
# its own: the lock of the directory of the simulation, the queue of its samples, and the merge of
# the samples ended into its outputs, in the order of the samples.

"""
    with_lock(f, dir)

`f()` holding the lock of the directory `dir`, its file `lock`: one left by a process that died
is passed over after a minute, five for another machine, see `mkpidlock`
"""
function with_lock(f::Function, dir::String)
    path = joinpath(dir, "lock")
    # tried every tenth of a second, a stale lock looked for every ten seconds: the wait of
    # mkpidlock misses the removal of the file by another machine on a network file system, and
    # the look for a stale lock warns of a clock skew each time
    tries, held = 0, false
    while held === false
        if tries > 0
            sleep(0.1)
        end
        held = trymkpidlock(path; stale_age = tries % 100 == 0 ? 60 : 0, refresh = 30)
        tries += 1
    end
    try
        return f()
    finally
        close(held)
    end
end

# the directory of the part numbered `part` of the simulation of directory `dir`
part_dir(dir::String, part::Int) = joinpath(dir, "parts", string(part))

"""
    Queue(samples, reached, next, free, current, failed)

the queue of the samples of a simulation, its file `queue.json`: the samples asked for, the
sample whose merge reached the tolerance, `nothing` before, the next one never taken, those
given back, the sample each part holds, by part, and the samples failed, with what went wrong
"""
mutable struct Queue
    samples::Int
    reached::Union{Nothing, Int}
    next::Int
    free::Vector{Int}
    current::Dict{Int, Int}
    failed::Dict{Int, String}
end

queue_file(dir::String) = joinpath(dir, "queue.json")

function read_queue(dir::String)
    q = JSON.parsefile(queue_file(dir))
    return Queue(q["samples"], q["reached"], q["next"], Int[ k for k in q["free"] ],
                 Dict{Int, Int}(parse(Int, p) => k for (p, k) in q["current"]),
                 Dict{Int, String}(parse(Int, k) => m for (k, m) in q["failed"]))
end

function write_queue(dir::String, q::Queue)
    tmp = queue_file(dir) * ".tmp"
    open(tmp, "w") do f
        JSON.print(f, Dict("samples" => q.samples, "reached" => q.reached, "next" => q.next,
                           "free" => sort(q.free),
                           "current" => Dict(string(p) => k for (p, k) in q.current),
                           "failed" => Dict(string(k) => m for (k, m) in q.failed)))
    end
    sync_file(tmp)
    Base.Filesystem.rename(tmp, queue_file(dir))
    return nothing
end

# the last sample the parts run: those asked for, or that reaching the tolerance
last_sample(q::Queue) = min(q.samples, something(q.reached, q.samples))

# the samples failed among those the parts run, see last_sample
failed_samples(q::Queue) = filter(≤(last_sample(q)), collect(keys(q.failed)))

# the error of a queue with a sample failed, naming the first
function failed_error(q::Queue)
    k = minimum(failed_samples(q))
    return ErrorException("sample $k failed, the samples stop: $(q.failed[k])")
end

# the sample a part takes next, given back or never taken, up to the last one
function next_sample(q::Queue)
    free = filter(≤(last_sample(q)), q.free)
    return !isempty(free) ? minimum(free) : q.next ≤ last_sample(q) ? q.next : nothing
end

"""
    take_sample!(dir, part; stop = false)

the sample the part `part` of the simulation of directory `dir` runs next: the one it holds, or
else the next one, `nothing` when none is left, `:stop` with `stop` when one is left, which it
does not take; refused once a sample failed. A sample held beyond the last one is given back,
see `last_sample`.
"""
function take_sample!(dir::String, part::Int; stop::Bool = false)
    return with_lock(dir) do
        q = read_queue(dir)
        k = get(q.current, part, nothing)
        if !isnothing(k) && k > last_sample(q)
            delete!(q.current, part)
            push!(q.free, k)
            write_queue(dir, q)
            k = nothing
        end
        if !isnothing(k)
            return stop ? :stop : k
        elseif !isempty(failed_samples(q))
            throw(failed_error(q))
        end
        k = next_sample(q)
        if isnothing(k) || stop
            return isnothing(k) ? nothing : :stop
        end
        filter!(≠(k), q.free)
        q.next = max(q.next, k + 1)
        q.current[part] = k
        write_queue(dir, q)
        return k
    end
end

# `sample` noted failed with `message` in the queue `q` of the simulation of directory `dir`, and
# when the parts run it, in its log and in its file error, which names the first sample failed
function note_failed!(dir::String, q::Queue, sample::Int, message::AbstractString)
    q.failed[sample] = message
    write_queue(dir, q)
    if sample ≤ last_sample(q)
        open(io -> println(io, "\nERROR in sample $sample: $message"), joinpath(dir, "log"), "a")
        k = minimum(failed_samples(q))
        write(joinpath(dir, "error"), "ERROR in sample $k: $(q.failed[k])\n")
    end
    return nothing
end

"""
    sample_failed!(dir, sample, message)

note that `sample` failed with `message`, so that the other parts stop, see `take_sample!`
"""
sample_failed!(dir::String, sample::Int, message::AbstractString) =
    with_lock(() -> note_failed!(dir, read_queue(dir), sample, message), dir)

"""
    sample_ended!(dir, part)

note that the part `part` ended its sample, whose rows are written, see `write_pending`
"""
function sample_ended!(dir::String, part::Int)
    with_lock(dir) do
        q = read_queue(dir)
        delete!(q.current, part)
        write_queue(dir, q)
    end
    return nothing
end

"""
    merge_samples!(dir, merging; limit = 100)

under the lock of the directory `dir` of a simulation, at most `limit` of its samples ended added
to its outputs in their order, see `replay!`, from the first not merged up to the first not
ended, or the first reaching the tolerance, `merging = (; format, tol, min_samples)`, see
`within`, which the queue then notes: its files written, then `merged.jls`, then the rows merged
removed, so that a crash leaves a merge to go on with. It returns `(first, last, reached)`, the
samples merged and whether the tolerance was reached, `nothing` for none. A sample that cannot
be merged is noted failed, see `sample_failed!`.
"""
function merge_samples!(dir::String, merging; limit::Int = 100)
    (; format, tol, min_samples) = merging
    m, marks = read_merged(dir)
    k = m + 1
    q = read_queue(dir)
    if !isnothing(q.reached) || !isfile(pending_file(dir, k))
        return nothing
    elseif !isempty(tol) && within(marks, dir, format, tol, min_samples)
        # reached by a merge that stopped before the queue noted it
        q.reached = m
        write_queue(dir, q)
        return nothing
    end
    o = reopen!(Outputs(; dir, format, sampled = true), marks)
    reached = false
    try
        while !reached && k ≤ m + limit && isfile(pending_file(dir, k))
            replay!(o, read_pending(dir, k))
            reached = !isempty(tol) && within(o.opened, dir, format, tol, min_samples)
            k += 1
        end
    catch e
        # the averaged files left as the merge before wrote them, the rows written to the others
        # cut by the next merge, which reopens them at their marks
        foreach(close, (h for h in values(o.opened) if !(h isa AveragedOpened)))
        if e isa InterruptException
            rethrow()
        end
        message = "sample $k cannot be merged: $(sprint(showerror, e))"
        note_failed!(dir, read_queue(dir), k, message)
        error(message)
    end
    marks = resume_marks(o)
    close!(o)
    write_merged(dir, k - 1, marks)
    if reached
        q.reached = k - 1
        write_queue(dir, q)
    end
    # those a crash left after the merge before as well
    for f in readdir(joinpath(dir, "pending"))
        n = tryparse(Int, chopprefix(f, "sample-"))
        if startswith(f, "sample-") && !isnothing(n) && n < k
            rm(joinpath(dir, "pending", f); force = true)
        end
    end
    return (first = m + 1, last = k - 1, reached)
end

"""
    merge_all!(dir, merging)

the samples ended of the simulation of directory `dir` merged, see `merge_samples!`, the lock
taken anew for every `limit` of them, so that it is held a short while; `(first, last, reached)`
as `merge_samples!` returns it, `nothing` for none
"""
function merge_all!(dir::String, merging)
    merged = nothing
    while true
        step = with_lock(() -> merge_samples!(dir, merging), dir)
        if isnothing(step)
            return merged
        end
        merged = isnothing(merged) ? step : (; merged.first, step.last, step.reached)
    end
end
