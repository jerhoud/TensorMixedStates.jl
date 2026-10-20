# The hints Julia prints after an error, for a program written for 1.x: the names 2.0 removed or
# renamed, the keywords it moved, a function given a Simulation, and the forms of 1.x of a few
# calls.
#
# Goes here: src/Hints.jl. The hints are part of the message of the error, as Julia shows it.

"""
    message(f)

the message of the error `f()` throws, hints included, as Julia shows it
"""
function message(f)
    try
        f()
    catch e
        return sprint(showerror, e)
    end
    error("no error was thrown")
end

@testset "The names of 1.x" begin
    # every name of the table, looked up in a module that defines none of them
    m = Module()
    for (name, h) in TMS.renamed_names
        msg = message(() -> Core.eval(m, name))
        @test occursin("TensorMixedStates 2.0: $h", msg)
    end
    @test occursin("tdvp became evolve(evolver, duration, state; algo = Tdvp(), nsteps)",
                   message(() -> Core.eval(m, :(tdvp(1, 2, 3)))))
    @test occursin("Upgrading from 1.x", message(() -> Core.eval(m, :HermiticityError)))
    # a name that is no name of 1.x has no hint
    @test !occursin("TensorMixedStates", message(() -> Core.eval(m, :undefined_name)))
    # nor a static parameter named as a gate of 1.x, which Julia tells apart from 1.11 on
    if VERSION ≥ v"1.11"
        parameter(x::Union{Int, Vector{T}}) where T = T
        @test !occursin("TensorMixedStates", message(() -> parameter(1)))
    end
    # a solver of ITensorMPS called on the objects of this package
    st = State{Pure}(System(2, Qubit()), "Up")
    @test occursin("dmrg became ground_state", message(() -> TMS.ITensorMPS.dmrg(Z(1), st)))
end

@testset "The keywords of 1.x" begin
    phases = [CreateState{Pure}(2, Qubit(), "Up")]
    h = Z(1) * Z(2)
    st = State{Pure}(System(2, Qubit()), "Up")
    ρ = mix(st)
    algo = Dmrg(nsweeps = 2)
    L = Dissipator(Sm)(1)
    cases = [
        (() -> GroundState(hamiltonian = h, nsweeps = 3), "algo = Dmrg(nsweeps"),
        (() -> SteadyState(lindbladian = L, tol = 1e-5), "algo = Dmrg(nsweeps"),
        (() -> ground_state(h, st; nsweeps = 3), "algo = Dmrg(nsweeps"),
        (() -> steady_state(L, ρ; algo, first_sweep = 2), "phase resumes its search"),
        (() -> evolve(-im * h, 0.1, st; order = 2), "ApproxW(order"),
        (() -> evolve(-im * h, 0.1, st; time_start = 1.), "time_start became time"),
        (() -> evolve(-im * h, 0.1, st; coefs = [t -> 1]), "[A, B] => [f, g]"),
        (() -> thermalize(h, 0.1, ρ; first_step = 2), "first_step became step"),
        (() -> SimData(; phases, checkpoint_interval = 10), "runTMS(sim_data; checkpoint_interval"),
        (() -> SimData(; phases, max_time = 10), "TMS_MAX_TIME"),
        (() -> SimData(; phases, data_format = "%g"), "data_format became format"),
        (() -> runTMS(SimData(; phases); restart = true), "cleanTMS(sim_data)"),
        (() -> CreateState{Pure}(2, Qubit(), "Up"; seed = 1), "SimData(seed"),
        (() -> expect(st, Z(1); normalize = false), "multiply the value by trace(state)"),
        # those 1.6 had, which the development of 2.0 had renamed before
        (() -> SimData(; phases, final_measures = []), "final_measures became final_measurements"),
        (() -> GroundState(hamiltonian = h, algo, measures = "d" => [Z]),
         "measures became measurements"),
        (() -> GroundState(hamiltonian = h, tolerance = 1e-8), "algo = Dmrg(nsweeps"),
        (() -> Dmrg(nsweeps = 2, hamiltonian = h), "the phase is GroundState(hamiltonian"),
        (() -> Tdvp(n_expand = 2), "n_expand became expand_period"),
        (() -> ApproxW(order = 2, n_hermitianize = 1), "hermitianize_period"),
        (() -> evolve(-im * h, 0.1, st; nsweeps = 2), "nsweeps became nsteps"),
        (() -> evolve(-im * h, 0.1, st; first_sweep = 2), "first_sweep became step"),
        (() -> steady_state(L, ρ; algo, alg = "naive"), "alg was removed: the MPO of L†L"),
        (() -> steady_state(L, ρ; algo, mpo_algo = "naive"), "mpo_algo was removed"),
        (() -> SteadyState(lindbladian = L, algo, mpo_algo = "naive"), "mpo_algo was removed"),
        (() -> partial_trace(ρ, [1]; keepers = true), "partial_trace(state; keep = positions)"),
        (() -> PartialTrace(trace_positions = [1]), "trace_positions became trace_out"),
        (() -> PartialTrace(keep_positions = [1]), "keep_positions became keep"),
        (() -> load_state("state.h5", "s"; system = st.system), "systems with the same sites"),
        (() -> controlled(X; type = 1), "type was removed"),
    ]
    for (f, h) in cases
        @test occursin("TensorMixedStates 2.0", message(f)) && occursin(h, message(f))
    end
    # a keyword that was never one has no hint
    @test !occursin("TensorMixedStates 2.0", message(() -> evolve(-im * h, 0.1, st; tolerance = 1)))
end

@testset "Calls of 1.x" begin
    st = State{Pure}(System(2, Qubit()), "Up")
    sim = runTMS(SimData(phases = [CreateState{Pure}(2, Qubit(), "Up")]); output = devnull)
    @test occursin("the functions take a state, sim.state", message(() -> apply(X(1), sim)))
    # not for a function taking the simulation, whose other argument is wrong
    @test !occursin("sim.state", message(() -> output(sim, Z(1))))
    @test occursin("a Simulation is made by runTMS", message(() -> Simulation(st)))
    @test occursin("measure(state, list; time = t)", message(() -> measure(st, [Z(1)], 0.3)))
    @test occursin("returns the Check", message(() -> Check("name", Z(1), 1.)))
    @test occursin("takes no type and no sites",
                   message(() -> Operator{1}("A", [1. 0.; 0. 1.], 1, Qubit())))
    @test occursin("output(Simulation(sim, sim.state; energy = e), list)",
                   message(() -> output(sim, [Z(1)]; energy = 1.)))
    @test occursin("variance(state, op)", message(() -> variance(Z(1), st)))
    @test occursin("SubRenyi2(1:k)", message(() -> SubRenyi2(1)))
    @test occursin("renyi2(state, 1:k)", message(() -> renyi2(st, 1)))
    @test occursin("matrix(X, site)", message(() -> matrix("X", Qubit())))
    @test occursin("systems with the same sites are equal", message(() -> State(st.system, st)))
    @test occursin("partial_trace(state; trace_out = positions)",
                   message(() -> partial_trace(mix(st), [1])))
    # the sites of PartialTrace are asked for by an error of their own, which leaves the
    # keywords of 1.x to their hints
    @test occursin("PartialTrace takes the sites to keep", message(() -> PartialTrace()))
    # the steps a list of measurements reads
    @test occursin("the :sweep of 1.x is :step", message(() -> SimData(phases = [
        CreateState{Pure}(2, Qubit(), "Up"),
        Gates(gates = X(1), final_measurements = "d" => [:sweep])])))
    # a local state given as a number, as 1.x took it, is the basis state of that number
    @test occursin("is Level => 1", message(() -> State{Pure}(st.system, 1)))
    @test occursin("is Proj(Level => 1)", message(() -> Proj(1)))
    @test occursin("is Reset(Level => 1)", message(() -> Reset(1)))
end

@testset "data_to_frame without DataFrames" begin
    # in a process that has not loaded DataFrames, the error says to load it
    script = """
        using TensorMixedStates
        try
            data_to_frame(Dict())
        catch e
            print(sprint(showerror, e))
        end
        """
    cmd = `$(Base.julia_cmd()) --project=$(Base.active_project()) --startup-file=no -e $script`
    out = read(cmd, String)
    @test occursin("data_to_frame needs the DataFrames package: run using DataFrames", out)
end
