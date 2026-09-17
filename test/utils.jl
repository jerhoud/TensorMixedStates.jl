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
exact one particle correlation matrix <c^dag_i c_j> of a pure state, computed in the full
Fock basis with an explicit Jordan-Wigner string. Independent of the fermionic machinery
of the package, so it can be used as a reference for it. Only usable on small systems.
"""
function exact_fermionic_correlations(state::State{Pure}, site)
    n = length(state)
    mc = matrix(C, site)
    mf = matrix(F, site)
    id = matrix(Id, site)
    c(j) = foldl(kron, [ k < j ? mf : (k == j ? mc : id) for k in 1:n ])
    idx = [ SysIndex{Pure}(state.system, k) for k in 1:n ]
    psi = reshape(Array(reduce(*, [state.state[k] for k in 1:n]), reverse(idx)...), dim(site)^n)
    psi /= norm(psi)
    return [ psi' * (c(i)' * c(j)) * psi for i in 1:n, j in 1:n ]
end

"""
checks whether measurements are equal between a random pure state and its computed mixed represenetation
"""
function check_mix(dims, sites, measures)
    for d in dims, (s, m) in zip(sites, measures)
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
