# Helpers shared by every test group. Nothing here runs a test by itself: this file is
# included before the groups, outside of the global testset, so that the macros it
# defines exist by the time the group files are parsed.

"""
executes the given phases
"""
function test_phases(phases)
    try
        runTMS(SimData(;phases); output = devnull)
    catch
        println("test_phases failed for $phases")
        rethrow()
    end
end

check(a, b, tol=1e-8) = "" => Check("", a, b, tol)

"""
ask the simulation being run to stop, as a user does by creating the file `stop` in its
directory: the one directory of the working directory that holds the marker `running`,
since a measurement is not told the directory of its simulation, and `runTMS` leaves the
working directory as it is
"""
function request_stop()
    for d in readdir()
        if isfile(joinpath(d, "running"))
            touch(joinpath(d, "stop"))
        end
    end
end

"""
test whether a statement execute without throwing an exception
"""
macro test_ok(a)
    :(@test ($(esc(a)); true))
end

"""
test a statement (with `@test_ok`) with `type = Pure` and `type = Mixed` 
"""
macro test_pm(a)
    tp = :type
    quote
        @test_ok (($(esc(tp)) -> $(esc(a)))(Pure))
        @test_ok (($(esc(tp)) -> $(esc(a)))(Mixed))
    end
end

"""
the vector of a pure state in the full basis, site 1 most significant, normalised. Only usable
on small systems.
"""
function dense_vector(state::State{Pure})
    n = length(state)
    idx = [ SysIndex{Pure}(state.system, k) for k in 1:n ]
    psi = vec(Array(reduce(*, [state.state[k] for k in 1:n]), reverse(idx)...))
    return psi / norm(psi)
end

"""
the matrix, on `n` sites `site`, of the product of operators of one site `factors`, pairs
`op => i` in the order of the product, written with explicit Jordan-Wigner strings: a
fermionic operator on site `i` carries `F` on every site before it. Independent of the
fermionic machinery of the package, so it can be used as a reference for it.
"""
function jw_matrix(site, n, factors)
    mf = matrix(F, site)
    id = matrix(Id, site)
    placed(a, j) = foldl(kron, [ k == j ? matrix(a, site) : (k < j && isfermionic(a) ? mf : id)
                                 for k in 1:n ])
    return prod(placed(a, j) for (a, j) in factors)
end

"""
exact one particle correlation matrix <c^dag_i c_j> of a pure state, see `jw_matrix`
"""
function exact_fermionic_correlations(state::State{Pure}, site)
    n = length(state)
    psi = dense_vector(state)
    return [ psi' * jw_matrix(site, n, [dag(C) => i, C => j]) * psi for i in 1:n, j in 1:n ]
end

"""
checks whether measurements are equal between a random pure state and its computed mixed represenetation
"""
function check_mix(dims, sites, measurements)
    for d in dims, (s, m) in zip(sites, measurements)
        sys = System(10, s)
        stp = RandomState{Pure}(sys, d)
        stm = mix(stp)
        meas = Measure(m...)
        mp = last.(measure(stp, meas))
        mm = last.(measure(stm, meas))
        if mp ≈ mm
            continue
        end
        error("mix check $d, $s, $m fails with $mp and $mm")
    end
end

"""
the dense matrix of the MPO `mpo` on the sites of `state`, the first site varying slowest
"""
function dense(state::State{R}, mpo) where R
    t = prod(mpo)
    s = [ SysIndex{R}(state.system, k) for k in 1:length(state) ]
    d = prod(dim, s)
    return reshape(Array(t, reverse([ x' for x in s ])..., reverse(s)...), d, d)
end
