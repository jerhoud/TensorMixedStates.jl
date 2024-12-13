export Limits, no_limits
export Dmrg, Tdvp, Gates, Evolve, ToMixed, SaveState, LoadState, CreateState

struct Trotter
  order::Int
  correct_trace::Int
  correct_hermitian::Int
end

@kwdef struct Output
  file::String = ""
  state_info::Bool = false
  time_format::String = "%g"
  data_format::String = "%12.8f"
  observables::Vector{ProdLit} = []
  expect::Vector{Function} = []
  correl::Vector{Tuple{Function, Function}} = []
  entanglement_entropy::Vector{Int} = []
  show_ee_spectrum::Int = 0
  checks::Vector{Pair{ProdLit, Function}} = []
end

const no_output = Output()

@kwdef struct Limits
  cutoff::Float64
  maxdim::Int
end

const no_limits = Limits(cutoff = 0, maxdim = typemax(Int))

@kwdef struct Gates
  name::String
  output::Output
  gates::ProdLit
  limits::Limits
end

@kwdef struct Tdvp
  name::String = "tdvp evolution"
  output::Output = no_output
  limits::Limits
  duration::Float64
  tau::Float64
  hamiltonian::Lits
  output_periodicity::Int = 1
  data_output::Output = no_output
  expand::Int = 0
end

@kwdef struct Dmrg
  name::String = "Dmrg optimization"
  output::Output = no_output
  hamiltonian::Lits
  cutoff::Float64
  maxdim::Vector{Int}
  nsweep::Int
  output_periodicity::Int = 1
  data_output::Output = no_output
end



@kwdef struct LoadState
  name::String = "Loading state"
  output::Output = no_output
  file::String
  object_name::String = "state"
  limits::Limits = no_limits
end

@kwdef struct CreateState
  name::String = "Creating state"
  output::Output = no_output
  type::Union{Type{Pure}, Type{Mixed}}
  system
  state::Union{Nothing, String, Vector{String}} = nothing
  randomize::Int = 0
end

@kwdef struct SimData
  name::String=""
  description::String=""
  restart::Bool=false
  clean::Bool=false
  debug::Bool=false
  phases
end