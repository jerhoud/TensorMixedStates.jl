# The hints Julia prints after an error: for a program written for 1.x, the names that 2.0 removed
# or renamed and the keywords it moved, a function given a Simulation where it takes a state, and
# data_to_frame without DataFrames. They are registered when the package is loaded.

"""
    renamed_names

the names of 1.x that 2.0 removed or renamed, with what to write instead, which Julia shows after
the `UndefVarError` of a program written for 1.x
"""
const renamed_names = Dict{Symbol, String}(
    :tdvp => "tdvp became evolve(evolver, duration, state; algo = Tdvp(), nsteps)",
    :approx_W => "approx_W became evolve(evolver, duration, state; " *
                 "algo = ApproxW(order = 4), nsteps)",
    :dmrg => "dmrg became ground_state(hamiltonian, state; algo = Dmrg(nsweeps = 10))",
    :thermal_state => "thermal_state became thermalize(hamiltonian, beta, state; nsteps), which " *
                      "returns (log_trace, state)",
    :RandomState => "RandomState became random_state(Pure, system, linkdims), " *
                    "random_state(Mixed, system, linkdims) or random_state(state, linkdims)",
    :StateFunc => "StateFunc was removed: a list of measurements takes \"name\" => c -> f(c.state)",
    :TimeFunc => "TimeFunc was removed: a list of measurements takes \"name\" => c -> f(c.time)",
    :RealValue => "RealValue was removed: write \"name\" => c -> real(expect(c.state, op))",
    :ImaginaryValue => "ImaginaryValue was removed: write " *
                       "\"name\" => c -> imag(expect(c.state, op))",
    :ComplexValue => "ComplexValue was removed: the value of an operator that is not " *
                     "Hermitian is complex",
    :trace2 => "trace2 became purity",
    :EE => "EE became EntanglementEntropy",
    :Linkdim => "Linkdim became MaxLinkdim",
    :Mutual_Info_Renyi2 => "Mutual_Info_Renyi2 became MutualInfoRenyi2, of a block of sites: " *
                           "MutualInfoRenyi2(1:k) for the first k",
    :DataToFrame => "DataToFrame became data_to_frame",
    :Trace2 => "Trace2 became Purity",
    :TraceError => "TraceError became TraceDeviation",
    :HermiticityError => "HermiticityError became HermiticityDeviation",
    :Basis => "Basis became Level",
    :SetState => "SetState became Reset, a channel: γ * Dissipator(Reset(state)) resets at the " *
                 "rate γ in an evolver",
    :Left => "Left was removed: write the channel Gate(A), ρ → AρA†, or the generator " *
             "Evolver(K), ρ → Kρ + ρK†",
    :Right => "Right was removed: write the channel Gate(A), ρ → AρA†, or the generator " *
              "Evolver(K), ρ → Kρ + ρK†",
    :TdvpObserver => "the observers were removed: measure in an Evolve phase, or call evolve " *
                     "step by step and measure between the calls",
    :ApproxWObserver => "the observers were removed: measure in an Evolve phase, or call evolve " *
                        "step by step and measure between the calls",
    :ThermalObserver => "the observers were removed: measure in a Thermalize phase, or call " *
                        "thermalize step by step and measure between the calls",
    :DmrgObserver => "the observers were removed: measure in a GroundState phase, and stop a " *
                     "search on its energy with Dmrg(tol = …)",
    :get_sim_file => "get_sim_file was removed: write output(sim, destination => measurements), " *
                     "or a Destination of one's own",
    :close_sim_files => "close_sim_files was removed: runTMS closes the files of a simulation",
    :log_msg => "log_msg became log_message",
    :compact => "compact was removed: the MPOs are compacted when they are built; keep a PreMPO " *
                "or a Measure to build them once",
    :simplify => "simplify was removed: an operator is lowered on the system it acts on",
    :OpType => "the types of operators were removed: write @def_operators T [name = matrix, …] " *
               "and Operator{N}(name, definition)",
    :plain_op => "the types of operators were removed: write @def_operators T [name = matrix, …]",
    :fermionic_op => "the types of operators were removed: write @def_operators T " *
                     "[name = matrix, …], the parity of an operator is found from " *
                     "fermion_parities of the site",
    :selfadjoint_op => "the types of operators were removed: write @def_operators T " *
                       "[name = matrix, …]",
    :involution_op => "the types of operators were removed: write @def_operators T " *
                      "[name = matrix, …]",
    :isfermionic => "isfermionic was removed: the parity of an operator is found from " *
                    "fermion_parities of its sites",
    :hasfermionic => "hasfermionic was removed: the parity of an operator is found from " *
                     "fermion_parities of its sites",
    :has_fermionic => "has_fermionic was removed: the parity of an operator is found from " *
                      "fermion_parities of its sites",
    :conserve_string => "conserve_string became conserving(T; conserve), in the keyword " *
                        "constructor of a site type",
    :state => "state(site, name) is no longer exported: TensorMixedStates.local_state(site, name)",
    :identity_operator => "identity_operator became Id",
    :Index => "Index(site) was removed: the index of site i is siteind(Pure, system, i)",
    :SysIndex => "SysIndex was removed: the index of site i is siteind(Pure, system, i)",
    :make_approx_W1 => "make_approx_W1 is no longer exported: TensorMixedStates.make_approx_W1",
    :make_approx_W2 => "make_approx_W2 is no longer exported: TensorMixedStates.make_approx_W2",
    :run_steps => "run_steps was removed: a phase of one's own made of steps, see " *
                  "TensorMixedStates.prepare_phase, or Update(f; nsteps)",
    :resume_step => "resume_step was removed: a phase of one's own made of steps resumes by " *
                    "itself, see TensorMixedStates.prepare_phase",
    :resume_time => "resume_time was removed: a phase of one's own made of steps resumes by " *
                    "itself, see TensorMixedStates.prepare_phase",
    :committed_time => "committed_time was removed: a phase of one's own made of steps resumes " *
                       "by itself, see TensorMixedStates.prepare_phase",
    :phase_system => "phase_system was removed: systems compare by their content",
    :Phases => "Phases was removed: a phase is a subtype of TensorMixedStates.AbstractPhase",
    :Algo => "Algo became EvolutionAlgorithm, and SearchAlgorithm for the searches",
    :create_graph_state => "create_graph_state was removed: write " *
                           "CreateState(type = Pure(), state = graph_state(graph))",
    :Fup => "Fup was removed: write parity(Nup)",
    :Fdn => "Fdn was removed: write parity(Ndn)",
    :Swap => "Swap is no longer in Qubits: TensorMixedStates exports it, for two sites of any " *
             "one type",
    :controlled => "controlled is no longer in Qubits: TensorMixedStates exports it, for a " *
                   "first site of any type",
    :Zd => "Zd became Z of Qudits",
    :Xd => "Xd became X of Qudits",
    :Hd => "Hd became H of Qudits",
    :Sumd => "Sumd became Sum of Qudits",
)

"""
    moved_keywords

the keywords of 1.x that 2.0 moved, by the function or the type that took them, with what to
write instead, which Julia shows after the `MethodError` of a call written for 1.x
"""
const moved_keywords = let
    search_algo = "the search method is algo = Dmrg(nsweeps = …, tol = …, noise = …, krylov = …)"
    evolve_algo = "the method is algo = Tdvp(krylov = …, expand_period = …, " *
                  "hermitianize_period = …) or ApproxW(order = …, w = …, apply_algo = …, " *
                  "hermitianize_period = …)"
    observers = "the observers were removed: call the function step by step and measure between " *
                "the calls, or measure in a phase"
    d = Dict{Tuple{Any, Symbol}, String}()
    for f in (GroundState, SteadyState, ground_state, steady_state),
        k in (:nsweeps, :tol, :noise, :krylov)
        d[(f, k)] = search_algo
    end
    for f in (ground_state, steady_state)
        d[(f, :first_sweep)] = "first_sweep was removed: a GroundState or SteadyState phase " *
                               "resumes its search itself"
        d[(f, Symbol("observer!"))] = observers
    end
    for f in (evolve, thermalize)
        for k in (:order, :w, :apply_algo, :expand_period, :hermitianize_period, :krylov)
            d[(f, k)] = evolve_algo
        end
        d[(f, :first_step)] = "first_step became step"
        d[(f, Symbol("observer!"))] = observers
    end
    for k in (:n_expand, :n_hermitianize)
        d[(evolve, k)] = evolve_algo
    end
    d[(evolve, :nsweeps)] = "nsweeps became nsteps, the number of steps"
    d[(evolve, :first_sweep)] = "first_sweep became step: time is then the time the evolution " *
                                "starts at, duration and nsteps those of the steps left"
    d[(evolve, :coefs)] = "an evolver depending on time is written [A, B] => [f, g]"
    d[(evolve, :time_start)] = "time_start became time"
    d[(thermalize, :log_trace)] = "the log_trace of successive calls of thermalize add up"
    for k in (:checkpoint_interval, :max_time, :threading)
        d[(SimData, k)] = "$k is a keyword of runTMS, runTMS(sim_data; $k = …), or the " *
                          "environment variable TMS_$(uppercase(string(k)))"
    end
    for k in (:time_format, :data_format)
        d[(SimData, k)] = "$k became format, of every number written"
    end
    d[(runTMS, :restart)] = "restart was removed: put a file restart in the directory of the " *
                            "simulation, or call cleanTMS(sim_data) before runTMS"
    d[(runTMS, :clean)] = "clean was removed: call cleanTMS(sim_data)"
    d[(CreateState, :seed)] = "the seed is that of the simulation, SimData(seed = …), and " *
                              "the draws of the phases use its generator"
    for T in (CreateState, LoadState, SaveState, ToMixed, Weaken, PartialTrace, Gates, Evolve,
              GroundState, SteadyState, SimData)
        d[(T, :final_measures)] = "final_measures became final_measurements"
    end
    for T in (Evolve, GroundState, SteadyState)
        d[(T, :measures)] = "measures became measurements"
        d[(T, :measures_period)] = "measures_period became measurements_period"
    end
    for T in (GroundState, SteadyState)
        d[(T, :tolerance)] = search_algo
    end
    d[(Dmrg, :hamiltonian)] = "Dmrg is now the method of a search: the phase is " *
                              "GroundState(hamiltonian = …, algo = Dmrg(nsweeps = …))"
    d[(Tdvp, :n_expand)] = "n_expand became expand_period"
    for T in (Tdvp, ApproxW)
        d[(T, :n_hermitianize)] = "n_hermitianize became hermitianize_period"
    end
    d[(PartialTrace, :trace_positions)] = "trace_positions became positions"
    d[(PartialTrace, :keep_positions)] = "keep_positions = p became positions = p, keep = true"
    for (f, k) in ((steady_state, :alg), (steady_state, :mpo_algo), (SteadyState, :mpo_algo))
        d[(f, k)] = "$k was removed: the MPO of L†L is built from that of L with as few " *
                    "channels as its rank, mpo_limits truncating it if given"
    end
    d[(partial_trace, :keepers)] = "keepers became keep"
    d[(load_state, :system)] = "system was removed: systems with the same sites are equal"
    for f in (named, controlled)
        d[(f, :type)] = "type was removed with the types of operators: the parity of an " *
                        "operator is found from fermion_parities of its sites"
    end
    d[(expect, :normalize)] = "normalize was removed: multiply the value by trace(state)"
    d
end

"""
    keyword_names(kwargs)

the names of the keywords of a call that failed, as Julia hands them to a hint, pairs of a name
and a value or of a name and a type depending on its version
"""
keyword_names(kwargs) = Symbol[ first(k) for k in kwargs ]

"""
    takes_simulation(f)

whether a method of `f` takes a `Simulation`, the error of a call to it being then no
`Simulation` given for a state
"""
takes_simulation(f) =
    any(m -> any(t -> t isa Type && t <: Simulation, Base.unwrap_unionall(m.sig).parameters[2:end]),
        methods(f))

"""
    method_hint(io, e, argtypes, kwargs)

print what to write instead after the `MethodError` `e` of a call written for 1.x: a keyword it
moved, a `Simulation` where a state is taken, a solver of ITensorMPS called on objects of this
package, and the forms of 1.x of `measure`, `Check`, `Operator`, `Simulation`, `output`,
`variance`, the Rényi entropies of a block, `matrix` and `State`
"""
function method_hint(io::IO, e::MethodError, argtypes, kwargs)
    f = e.f
    for k in keyword_names(kwargs)
        h = get(moved_keywords, (f, k), nothing)
        if !isnothing(h)
            print(io, "\nTensorMixedStates 2.0: ", h)
            return nothing
        end
    end
    types = collect(argtypes)
    if f === data_to_frame && isnothing(Base.get_extension(@__MODULE__, :DataFramesExt))
        print(io, "\ndata_to_frame needs the DataFrames package: run using DataFrames")
    elseif f === data_to_frame && length(types) == 1 && types[1] <: Simulation
        print(io, "\ndata_to_frame takes the rows of a Data: data_to_frame(sim.data[\"name\"])")
    elseif f === Simulation && length(types) == 1 && types[1] <: AbstractState
        print(io, "\nTensorMixedStates 2.0: a Simulation is made by runTMS; outside a run, call " *
                  "the functions on the state")
    elseif f === output && !isempty(kwargs)
        print(io, "\nTensorMixedStates 2.0: output measures with the values of the simulation, " *
                  "output(Simulation(sim, sim.state; energy = e), list)")
    elseif any(t -> t isa Type && t <: Simulation, types) && !takes_simulation(f)
        print(io, "\nTensorMixedStates 2.0: the functions take a state, sim.state, and " *
                  "Simulation(sim, state) gives the simulation another one")
    elseif f === ITensorMPS.tdvp || f === ITensorMPS.dmrg
        print(io, "\nTensorMixedStates 2.0: ", renamed_names[nameof(f)])
    elseif f === measure && length(types) == 3
        print(io, "\nTensorMixedStates 2.0: the time is a keyword, measure(state, list; time = t)")
    elseif f === Check && !isempty(types) && types[1] <: AbstractString
        print(io, "\nTensorMixedStates 2.0: a function of the context returns the Check, " *
                  "\"name\" => c -> Check(expect(c.state, op), reference; tol = 1e-8)")
    elseif f isa Type && f <: Operator && length(types) > 2
        print(io, "\nTensorMixedStates 2.0: Operator{N}(name, definition) takes no type and no " *
                  "sites")
    elseif f === variance && length(types) == 2 && types[2] <: AbstractState
        print(io, "\nTensorMixedStates 2.0: the state comes first, variance(state, op)")
    elseif (f === SubRenyi2 || f === MutualInfoRenyi2) && length(types) == 1 &&
           types[1] <: Integer
        print(io, "\nTensorMixedStates 2.0: a block of sites is given by its positions, " *
                  "$(nameof(f))(1:k) for the first k")
    elseif (f === renyi2 || f === mutual_info_renyi2) && length(types) == 2 &&
           types[2] <: Integer
        print(io, "\nTensorMixedStates 2.0: a block of sites is given by its positions, " *
                  "$(nameof(f))(state, 1:k) for the first k")
    elseif f === matrix && !isempty(types) && types[1] <: AbstractString
        print(io, "\nTensorMixedStates 2.0: matrix takes the operator, matrix(X, site)")
    elseif f === State && length(types) == 2 && types[1] <: System && types[2] <: AbstractState
        print(io, "\nTensorMixedStates 2.0: State(system, state) was removed: systems with the " *
                  "same sites are equal")
    end
    return nothing
end

function __init__()
    Base.Experimental.register_error_hint(UndefVarError) do io, e
        h = get(renamed_names, e.var, nothing)
        if !isnothing(h)
            print(io, "\nTensorMixedStates 2.0: ", h,
                  " (see Upgrading from 1.x in its documentation)")
        end
    end
    Base.Experimental.register_error_hint(method_hint, MethodError)
end
