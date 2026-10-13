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
    cases = [
        (() -> GroundState(hamiltonian = h, nsweeps = 3), "algo = Dmrg(nsweeps"),
        (() -> SteadyState(lindbladian = Dissipator(Sm)(1), tol = 1e-5), "algo = Dmrg(nsweeps"),
        (() -> ground_state(h, st; nsweeps = 3), "algo = Dmrg(nsweeps"),
        (() -> steady_state(Dissipator(Sm)(1), ρ; algo, first_sweep = 2), "phase resumes its search"),
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
    @test occursin("a Simulation is made by runTMS", message(() -> Simulation(st)))
    @test occursin("measure(state, list; time = t)", message(() -> measure(st, [Z(1)], 0.3)))
    @test occursin("returns the Check", message(() -> Check("name", Z(1), 1.)))
    @test occursin("takes no type and no sites",
                   message(() -> Operator{1}("A", [1. 0.; 0. 1.], 1, Qubit())))
end
