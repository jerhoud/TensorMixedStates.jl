# Getting data in and out: saving a state to a file and reading it back, with a round trip
# check whenever possible, and printing what has to read back.
#
# Resetting a site is tested in channels.jl, the partial trace in observables.jl and the output
# of a simulation in outputs.jl and runner.jl.

# Sites whose fields are not all numbers, which version 1 of the state file format could not
# carry, writing every field as a Float64. The fields cover every kind a state file accepts but
# the Conservation of the built in sites, and the reference file of version 2 holds a Kindly.
struct Kindly <: AbstractSite
    conserve::Symbol
    on::Bool
    label::String
    n::Int
    x::Float64
    void::Nothing
end

TensorMixedStates.dim(::Kindly) = 2

# a field of a kind no state file can carry, to check that the refusal names it
struct Unkindly <: AbstractSite
    range::UnitRange{Int}
end

TensorMixedStates.dim(::Unkindly) = 2

# a site type of a module inside another one, which a state file named by the last name of its
# module alone and could not find again
module Enclosing
    using TensorMixedStates
    struct Enclosed <: AbstractSite end
    TensorMixedStates.dim(::Enclosed) = 2
end

# fields a state file wrote as their printed form and could not read back: a float of another
# width prints as "0.1f0", and an integer beyond `Int` overflowed
struct Widely <: AbstractSite
    y::Float32
    big::UInt64
end

TensorMixedStates.dim(::Widely) = 2

# a site whose inner constructor replaces the one taking its fields, by which a state file and
# weaken rebuild a site
struct Guarded <: AbstractSite
    conserve::Conservation
    Guarded(; conserve = ()) = new(Conservation(new(Conservation()), conserve))
end

TensorMixedStates.dim(::Guarded) = 2

@def_operators Guarded [ N = [0. 0. ; 0. 1.] ]

@testset "Weakening beside a site whose conserve field is not a Conservation" begin
    # such a site conserves nothing, and weakening, which the trace of a state
    # conserving something strongly goes through, leaves it as it is
    sys = System([Kindly(:sz, true, "a label", 3, 1.5, nothing), Fermion(conserve = strong(N))])
    ρ = mix(State{Pure}(sys, Any[[1., 0.], "Occ"]))
    @test trace(ρ) ≈ 1
    @test expect(ρ, N(2)) ≈ 1
    @test_ok weaken(ρ, ())
end

@testset "A site rebuilt from its fields must keep the constructor taking them" begin
    sys = System([Guarded(conserve = strong(N)), Guarded(conserve = strong(N))])
    st = State{Pure}(sys, Any[[0., 1.], [1., 0.]])
    @test expect(st, N(1)) ≈ 1
    kept = "Guarded must keep its constructor Guarded(conserve)"
    @test_throws "$kept for its conserved" weaken(st, ())
    mktempdir() do dir
        file = joinpath(dir, "guarded.h5")
        save_state(file, "g", st)
        @test_throws "$kept for a state file" load_state(file, "g")
    end
end

@testset "Saving and loading states" begin
    dir = mktempdir()
    file = joinpath(dir, "states.h5")

    # pure state round trip, read back on a system of the same sites, so it compares with the
    # state saved
    sys = System(4, Qubit())
    stp = State{Pure}(sys, ["Up", "Dn", "+", "i"])
    save_state(file, "pure", stp)
    lp = load_state(file, "pure")
    @test lp isa State{Pure}
    @test length(lp) == 4
    @test expect1(lp, Z) ≈ expect1(stp, Z)
    @test expect1(lp, X) ≈ expect1(stp, X)
    @test lp.system == sys
    @test inner(stp, lp) ≈ 1

    # mixed state round trip, in the same file
    stm = mix(stp)
    save_state(file, "mixed", stm)
    lm = load_state(file, "mixed")
    @test lm isa State{Mixed}
    @test expect1(lm, Z) ≈ expect1(stm, Z)
    @test trace(lm) ≈ 1
    @test inner(stm, lm) ≈ 1

    # both states are still readable
    @test expect1(load_state(file, "pure"), Z) ≈ expect1(stp, Z)
    @test inner(stp, load_state(file, "pure")) ≈ 1

    # saving under a name already used replaces the state
    save_state(file, "pure", State{Pure}(System(2, Qubit()), "Dn"))
    lp2 = load_state(file, "pure")
    @test length(lp2) == 2
    @test expect1(lp2, Z) ≈ [-1, -1]
    @test inner(State{Pure}(System(2, Qubit()), "Dn"), lp2) ≈ 1

    # sites with parameters are rebuilt identically
    sysp = System([Qubit(), Boson(4), Spin(3/2), Qboson(0.1, 3)])
    stq = State{Pure}(sysp, ["Up", "2", "1/2", "1"])
    save_state(file, "params", stq)
    lq = load_state(file, "params")
    @test lq.system.sites == sysp.sites

    # a loaded state can still be used for further computations
    @test trace(mix(lq)) ≈ 1

    # a site type of a module inside another one is rebuilt as well, its module being found
    # from the root one
    nested = State{Pure}(System(2, Enclosing.Enclosed()), ["0", "1"])
    save_state(file, "nested", nested)
    @test load_state(file, "nested").system.sites == nested.system.sites

    # SaveState and LoadState phases, the check fails if the state is not restored
    file2 = joinpath(dir, "phases.h5")
    @test_ok runTMS(SimData(phases = [
        CreateState{Pure}(3, Qubit(), ["Up", "Dn", "Up"]),
        SaveState(file = file2, statename = "st"),
        CreateState{Pure}(1, Qubit(), "Dn"),
        LoadState(file = file2, statename = "st", final_measurements = "" => [
            "z" => c -> Check(real(expect1(c.state, Z)), [1, -1, 1]; tol = 1e-12)])]);
        output = devnull)
end

@testset "A state file with its conserved quantities unsorted" begin
    # version 2 recorded them as a string, in the order of their declaration before they were
    # sorted: the state read back sat on sites no longer equal to those built now. The file is
    # the version 2 reference, its quantities put out of order
    file = joinpath(mktempdir(), "old.h5")
    TMS.HDF5.h5open(joinpath(@__DIR__, "reference", "state_v2.h5"), "r") do src
        TMS.HDF5.h5open(file, "w") do dst
            TMS.HDF5.copy_object(src, "electrons", dst, "st")
            TMS.HDF5.delete_object(dst["st"], "params")
            dst["st"]["params"] = fill("Ntot:0,1,1,2;2Sz:0,1,-1,0", 2)
        end
    end
    st = load_state(file, "st")
    @test st.system == System(2, Electron(conserve = (Ntot, 2Sz)))
end

@testset "A state read back compares with the one saved" begin
    # read back, the state lands on a system equal to the one it was saved from,
    # whose indices it shares, so the two compare without being put on a system
    dir = mktempdir()
    file = joinpath(dir, "onto.h5")
    sys = System(3, Qubit())
    a = random_state(Pure, sys, 4)
    save_state(file, "a", a)
    b = load_state(file, "a")
    @test b.system == sys
    @test inner(a, b) ≈ norm(a)^2
    @test fidelity(a, b) ≈ 1
    @test expect1(b, Z) ≈ expect1(a, Z)

    # and the same for a mixed state
    m = mix(a)
    save_state(file, "m", m)
    lm = load_state(file, "m")
    @test lm.system == sys
    @test inner(m, lm) ≈ inner(m, m)
    @test hs_fidelity(m, lm) ≈ 1

    # with the orthogonality centre the file records, which renaming the indices forgot: a
    # resumed run took other steps, moving the centre from elsewhere
    c = State(a, TMS.ITensorMPS.orthogonalize(a.state, 2))
    save_state(file, "c", c)
    @test TMS.ITensorMPS.ortho_lims(load_state(file, "c").state) == 2:2
end

@testset "Site fields of every kind" begin
    dir = mktempdir()
    file = joinpath(dir, "fields.h5")

    site = Kindly(:sz, true, "a label", 3, 1.5, nothing)
    st = State{Pure}(System(2, site), [1., 0.])
    save_state(file, "kinds", st)
    l = load_state(file, "kinds")

    # each field comes back as what it was, not as the Float64 version 1 would have made of it
    back = l.system.sites[1]
    @test back.conserve === :sz
    @test back.on === true
    @test back.label == "a label"
    @test back.n === 3
    @test back.x === 1.5
    @test back.void === nothing
    @test l.system.sites == st.system.sites

    # the first field is named conserve without being one, being a Symbol and not a
    # Conservation. A site carrying such a field conserves nothing
    @test TMS.conserved(site) == Conservation()

    # a field of an unsupported kind is refused by a message naming the site and the field,
    # rather than by a MethodError raised by convert somewhere inside HDF5
    bad = State{Pure}(System(2, Unkindly(1:3)), [1., 0.])
    @test_throws "its field range is a" save_state(file, "bad", bad)

    # and the refusal comes before the file is touched, so saving a state that cannot be
    # written over a name already in use leaves what was there intact
    @test_throws "its field range is a" save_state(file, "kinds", bad)
    @test load_state(file, "kinds").system.sites == st.system.sites

    wide = Widely(0.1f0, typemax(UInt64))
    save_state(file, "wide", State{Pure}(System(2, wide), [1., 0.]))
    @test load_state(file, "wide").system.sites[1] === wide
end

@testset "Reading a version 1 state file" begin
    # written by the released 1.3.0, see reference/make_state_v1.jl. Version 1 wrote every
    # site field as a Float64, and `load_state` has to go on reading it
    file = joinpath(@__DIR__, "reference", "state_v1.h5")

    p = load_state(file, "pure")
    @test p isa State{Pure}
    @test length(p) == 4
    @test p.system.sites == [Qubit(), Boson(4), Spin(3/2), Qboson(0.1, 3)]
    @test expect(p, Z(1)) ≈ 1
    @test expect(p, N(2)) ≈ 2
    @test expect(p, Sz(3)) ≈ 0.5
    @test inner(State{Pure}(p.system, ["Up", "2", "1/2", "1"]), p) ≈ 1

    m = load_state(file, "mixed")
    @test m isa State{Mixed}
    @test m.system.sites == p.system.sites
    @test trace(m) ≈ 1
    @test expect(m, Z(1)) ≈ 1
    @test inner(mix(State{Pure}(m.system, ["Up", "2", "1/2", "1"])), m) ≈ 1

    # a file of an unknown version is still refused, by a message naming what is accepted.
    # The group needs nothing but its version attribute, the check coming first
    future = joinpath(mktempdir(), "future.h5")
    TMS.HDF5.h5open(future, "cw") do h
        g = TMS.HDF5.create_group(h, "s")
        TMS.HDF5.attributes(g)["version"] = 99
    end
    @test_throws "expected one of 1, 2, 3" load_state(future, "s")
    # and so it is when the sites alone are read
    @test_throws "expected one of 1, 2, 3" TMS.saved_sites(future, "s")
end

@testset "Reading a version 2 state file" begin
    # written by the released 1.6.0, see reference/make_state_v2.jl: a file save_state writes
    # stays readable by every later version
    file = joinpath(@__DIR__, "reference", "state_v2.h5")

    e = load_state(file, "electrons")
    @test e isa State{Pure}
    @test e.system.sites == fill(Electron(conserve = (Ntot, 2Sz)), 2)
    @test real(expect1(e, Nup)) ≈ [1, 0]
    @test real(expect1(e, Ndn)) ≈ [0, 1]
    @test inner(State{Pure}(e.system, ["Up", "Dn"]), e) ≈ 1

    ρ = load_state(file, "fermions_strong")
    @test ρ isa State{Mixed}
    @test ρ.system.sites == fill(Fermion(conserve = strong(N)), 2)
    @test real(expect1(ρ, N)) ≈ [0.8, 0.2]
    @test expect(ρ, dag(C)(1) * C(2)) ≈ 0.4
    @test trace(ρ) ≈ 1
    p(v) = State{Pure}(ρ.system, v)
    ψ = (p(["Occ", "Emp"]) + 0.5 * p(["Emp", "Occ"])) / sqrt(1.25)
    @test inner(mix(ψ), ρ) ≈ 1

    k = load_state(file, "kindly")
    @test k.system.sites[1] == Kindly(:none, true, "a", 3, 0.5, nothing)
    @test expect(k, Z(2)) ≈ 1
    @test inner(State{Pure}(k.system, ["1", "Up"]), k) ≈ 1
end

@testset "Reading a version 3 state file" begin
    # written by version 2 before its release, see reference/make_state_v3.jl: what the sites
    # conserve is written typed, apart from their other fields
    file = joinpath(@__DIR__, "reference", "state_v3.h5")

    c = load_state(file, "charged")
    @test c isa State{Pure}
    @test c.system.sites == [Boson(4, conserve = N), Qubit(), Electron(conserve = (Ntot, 2Sz))]
    @test inner(State{Pure}(c.system, ["2", "Up", "Dn"]), c) ≈ 1

    ρ = load_state(file, "parity_strong")
    @test ρ isa State{Mixed}
    @test ρ.system.sites == fill(Fermion(conserve = strong(parity(N))), 2)
    @test trace(ρ) ≈ 1
    p(v) = State{Pure}(ρ.system, v)
    ψ = (p(["Emp", "Emp"]) + 0.5 * p(["Occ", "Occ"])) / sqrt(1.25)
    @test inner(mix(ψ), ρ) ≈ 1
end

@testset "Saving a state that carries charges" begin
    dir = mktempdir()
    file = joinpath(dir, "charged.h5")

    # what a site conserves is written with its strength, so a state read back draws the same
    # indices and is the same state
    for (name, site) in (("weak", Fermion(conserve = N)), ("strong", Fermion(conserve = strong(N))))
        sys = System(3, site)
        stp = State{Pure}(sys, ["Occ", "Emp", "Occ"])
        save_state(file, name, stp)
        lp = load_state(file, name)
        @test lp isa State{Pure}
        @test lp.system.sites == sys.sites
        @test flux(lp.state) == flux(stp.state)
        @test expect1(lp, N) ≈ expect1(stp, N)
        @test inner(stp, lp) ≈ 1

        stm = mix(stp)
        save_state(file, name * "-mixed", stm)
        lm = load_state(file, name * "-mixed")
        @test lm isa State{Mixed}
        @test flux(lm.state) == flux(stm.state)
        @test trace(lm) ≈ 1
        @test expect1(lm, N) ≈ expect1(stm, N)
        @test inner(stm, lm) ≈ 1
    end
end

@testset "Printing Limits" begin
    # printed as the call that builds it, the fields left nothing omitted, so that
    # the text evaluates back to the same fields
    @test repr(Limits()) == "Limits()"
    @test repr(Limits(cutoff = 1e-10, maxdim = 50)) == "Limits(cutoff = 1.0e-10, maxdim = 50)"
    # a value given is printed, the default of a function included
    @test repr(Limits(cutoff = eps())) == "Limits(cutoff = $(eps()))"
    for l in [Limits(), Limits(cutoff = 1e-10, maxdim = 50), Limits(maxdim = [10, 20, 50]),
              Limits(cutoff = 1e-14, maxdim = 100, mindim = 10)]
        @test eval(Meta.parse(repr(l))) == l
    end
end

@testset "Loading a partial trace of a charged system" begin
    # the qubit kept conserves nothing: its partial trace is densified, and a file holds it as
    # any state conserving nothing
    file = joinpath(mktempdir(), "traced.h5")
    sq = System([Qubit(), Fermion(conserve = N)])
    ρ = partial_trace(mix(State{Pure}(sq, ["+", "Occ"])); keep = [1])
    save_state(file, "traced", ρ)
    lt = load_state(file, "traced")
    @test lt.system == System(1, Qubit())
    @test trace(lt) ≈ 1
    @test expect(lt, X(1)) ≈ expect(ρ, X(1)) ≈ 1
end
