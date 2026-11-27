# The states of the Gaussian backend (GaussianState) against dense states: product states made by
# TMS, carried to the basis of the modes and lifted with the ancilla, Ψ = |0⟩ψ_even + i|1⟩ψ_odd.

if !isdefined(Main, :dense_majorana)
    include("gaussian_dense.jl")
end
if !isdefined(Main, :gc_to_modes)
    include("gaussian_tms.jl")
end

const gs_rng = Random.Xoshiro(20261128)

# the dense vector of a product state of TMS, in the basis of the modes, lifted with the ancilla
function gs_lifted(sys, ψ; axes = Dict{Int, Symbol}())
    W = reduce(kron, [gc_kind(s) == :qubit ? gc_axis_unitary(TensorMixedStates.axis_vector(get(axes, k, :Z))) :
                      Matrix{ComplexF64}(I, dim(s), dim(s)) for (k, s) in enumerate(sys.sites)])
    x = (W' * ψ)[gc_to_modes(sys)]
    M = trailing_zeros(length(x))
    even = [iseven(count(dense_bits(k - 1, M))) ? x[k] : 0 for k in eachindex(x)]
    vcat(even, im .* (x - even))
end

using .LA: dot, eigvals, eigen, Hermitian, I

gs_random_vector(d) = (v = randn(gs_rng, ComplexF64, d); v / LA.norm(v))

@testset "creation from product states" begin
    for (sites, axes) in (([Fermion(), Electron(), Qubit(), Fermion()], Dict(3 => :X)),
                          ([Qubit(), Qubit(), Qubit(), Qubit()], Dict(1 => :X, 2 => :Y)),
                          ([Qubit(), Fermion(), Qubit()], Dict{Int, Symbol}()),
                          ([Electron()], Dict{Int, Symbol}()))
        sys = System(sites)
        for _ in 1:5
            vs = [gs_random_vector(dim(s)) for s in sites]
            g = TensorMixedStates.GaussianState(State{Pure}(sys, vs); axes)
            Ψ = gs_lifted(sys, reduce(kron, vs); axes)
            @test g.cov ≈ dense_covariance(Ψ) atol = 1e-12
            @test g.cov * transpose(g.cov) ≈ I atol = 1e-12
            @test !TensorMixedStates.definite(g)
        end
    end
    # states of definite parity, Fock states and eigenstates of the axes
    sys = System([Fermion(), Electron(), Qubit()])
    g = TensorMixedStates.GaussianState(State{Pure}(sys, ["Occ", "Up", "Dn"]))
    @test TensorMixedStates.definite(g)
    @test g.cov ≈ dense_covariance(gs_lifted(sys, reduce(kron, [[0, 1.0], [0, 1.0, 0, 0], [0, 1.0]]))) atol = 1e-12
    # refused: a state that is not a product, a mixed state with coherences between occupations
    sys = System([Qubit() for _ in 1:4])
    @test_throws ArgumentError TensorMixedStates.GaussianState(random_state(Pure, sys, 2))
    @test_throws "coherences" TensorMixedStates.GaussianState(State{Mixed}(sys, "X+"))
end

# a random Gaussian state of sys (axes Z): with a definite parity, an even Gaussian of the modes; else
# a Gaussian Ψ of the modes and the ancilla, whose physical state is ψ_even - iψ_odd; the GaussianState
# holding it and its physical vector in the basis of TMS
function gs_random(sys; definite::Bool)
    M = sum(gc_kind(s) == :fermion ? length(fermion_species(s)) : 1 for s in sys.sites)
    if definite
        x = dense_gaussian(zeros(Bool, M), gs_rng)
        Ψ = vcat(x, zeros(ComplexF64, 2^M))
    else
        Ψ = dense_gaussian(zeros(Bool, M + 1), gs_rng)
        x = Ψ[1:2^M] - im * Ψ[2^M+1:end]
    end
    ψ = similar(x)
    ψ[gc_to_modes(sys)] = x
    layout = TensorMixedStates.MajoranaLayout(sys)
    g = TensorMixedStates.GaussianState{Pure}(sys, Dict{Int, Symbol}(), layout, dense_covariance(Ψ), nothing,
                                              definite ? [1] : nothing, 1.0, 64)
    g, ψ
end

# the reduced density matrix of the sites kept of the vector ψ (in the basis of TMS), from ExactState,
# which lays the fermionic signs of a State
gs_partial(sys, ψ, kept) = reduced_density_matrix(TensorMixedStates.ExactState{Pure}(sys, ψ / LA.norm(ψ), Val(:raw)), kept)

gs_entropy(ρ) = (λ = max.(eigvals(Hermitian(ρ)), 0); -sum(x * log(x) for x in λ if x > 0; init = 0.0))

@testset "values and entropies" begin
    for (sites, ops) in (([Fermion() for _ in 1:5], (C, dag(C), N)),
                         ([Qubit() for _ in 1:5], (X, Y, Z)),
                         ([Fermion(), Qubit(), Fermion(), Qubit()], nothing),
                         ([Electron(), Fermion(), Electron()], nothing)),
        definite in (true, false)
        sys = System(sites)
        g, ψ = gs_random(sys; definite)
        n = length(sites)
        dims = [dim(s) for s in sites]
        opsof(k) = gc_kind(sites[k]) == :qubit ? (X, Y, Z) : sites[k] isa Electron ? (Cup, Cdn, Nup) : (C, dag(C), N)
        for _ in 1:5
            ks = rand(gs_rng, 1:n, rand(gs_rng, 1:3))
            Os = [rand(gs_rng, opsof(k)) for k in ks]
            op = prod(o(k) for (o, k) in zip(Os, ks))
            O = prod(gc_op(sys, o, k) for (o, k) in zip(Os, ks))
            @test TensorMixedStates.expect_value(g, op) ≈ dot(ψ, O * ψ) atol = 1e-12
        end
        for kept in ([1], [1, 2], [2, 3], [1, 3], [n - 1, n], [n])
            ρ = gs_partial(sys, ψ, kept)
            @test reduced_density_matrix(g, kept) ≈ ρ atol = 1e-12
            @test von_neumann_entropy(g, kept) ≈ gs_entropy(ρ) atol = 1e-10
        end
        for cut in 1:n
            @test entanglement_entropy(g, cut) ≈ gs_entropy(gs_partial(sys, ψ, 1:cut)) atol = 1e-10
        end
    end
end

# a random product state of sys, its GaussianState and its vector in the basis of TMS
function gs_product(sys; axes = Dict{Int, Symbol}())
    vs = [gs_random_vector(dim(s)) for s in sys.sites]
    TensorMixedStates.GaussianState(State{Pure}(sys, vs); axes), reduce(kron, vs)
end

# the values of random operators on the sites and the entropies of the cuts, on g and on ψ
function gs_same(g, ψ, sys; nops = 6)
    n = length(sys.sites)
    opsof(k) = gc_kind(sys.sites[k]) == :qubit ? (X, Y, Z) : sys.sites[k] isa Electron ? (Cup, Cdn, Nup) : (C, dag(C), N)
    ok = true
    for _ in 1:nops
        ks = rand(gs_rng, 1:n, rand(gs_rng, 1:2))
        Os = [rand(gs_rng, opsof(k)) for k in ks]
        v = TensorMixedStates.expect_value(g, prod(o(k) for (o, k) in zip(Os, ks)))
        ok &= isapprox(v, dot(ψ, prod(gc_op(sys, o, k) for (o, k) in zip(Os, ks)) * ψ) / dot(ψ, ψ); atol = 1e-10)
    end
    for cut in 1:n-1
        ok &= isapprox(entanglement_entropy(g, cut), gs_entropy(gs_partial(sys, ψ, 1:cut)); atol = 1e-9)
    end
    for kept in ([2], [1, 3], [2, n])
        if !allunique(kept) || maximum(kept) > n
            continue
        end
        ok &= isapprox(von_neumann_entropy(g, kept), gs_entropy(gs_partial(sys, ψ, kept)); atol = 1e-9)
    end
    ok
end

@testset "evolution" begin
    # fermions with hopping and pairing; qubits XX + YY with fields; Ising along X with fields on both ends
    sys1 = System([Fermion() for _ in 1:5])
    H1 = sum(-(dag(C)(j) * C(j + 1) + dag(C)(j + 1) * C(j)) + 0.4 * (C(j) * C(j + 1) + dag(C)(j + 1) * dag(C)(j)) for j in 1:4) +
         0.3 * sum(N(j) for j in 1:5)
    H1d(sys) = sum(-(gc_op(sys, dag(C), j) * gc_op(sys, C, j + 1) + gc_op(sys, dag(C), j + 1) * gc_op(sys, C, j)) +
                   0.4 * (gc_op(sys, C, j) * gc_op(sys, C, j + 1) + gc_op(sys, dag(C), j + 1) * gc_op(sys, dag(C), j)) for j in 1:4) +
                   0.3 * sum(gc_op(sys, N, j) for j in 1:5)
    sys2 = System([Qubit() for _ in 1:5])
    q(sys, O, k) = gc_op(sys, O, k)
    H2 = sum(X(j) * X(j + 1) + Y(j) * Y(j + 1) for j in 1:4) + 0.5 * Z(2) - 0.3 * Z(4)
    H2d(sys) = sum(q(sys, X, j) * q(sys, X, j + 1) + q(sys, Y, j) * q(sys, Y, j + 1) for j in 1:4) + 0.5 * q(sys, Z, 2) - 0.3 * q(sys, Z, 4)
    H3 = -sum(Z(j) * Z(j + 1) for j in 1:4) - 0.7 * sum(X(j) for j in 1:5) + 0.3 * Z(1) - 0.2 * Z(5)
    H3d(sys) = -sum(q(sys, Z, j) * q(sys, Z, j + 1) for j in 1:4) - 0.7 * sum(q(sys, X, j) for j in 1:5) + 0.3 * q(sys, Z, 1) - 0.2 * q(sys, Z, 5)
    for (sys, H, Hd, axes) in ((sys1, H1, H1d, Dict{Int, Symbol}()), (sys2, H2, H2d, Dict{Int, Symbol}()),
                               (sys2, H3, H3d, Dict(k => :X for k in 1:5)))
        g, ψ = gs_product(sys; axes)
        pre = TensorMixedStates.prepared_evolver(g, -im * H)
        dt = 0.13
        for step in 1:4
            g = TensorMixedStates.evolve_step(Tdvp(), pre, g, (step - 1) * dt, dt, step, Limits())
        end
        ψ = exp(-im * 4dt * Hd(sys)) * ψ
        @test gs_same(g, ψ, sys)
    end
    # refused: an evolver that is not unitary
    g, _ = gs_product(sys1)
    @test_throws ArgumentError TensorMixedStates.prepared_evolver(g, -im * H1 - 0.1 * N(2))
end

@testset "ground states" begin
    sys = System([Fermion() for _ in 1:5])
    H = sum(-(dag(C)(j) * C(j + 1) + dag(C)(j + 1) * C(j)) + 0.4 * (C(j) * C(j + 1) + dag(C)(j + 1) * dag(C)(j)) for j in 1:4) +
        0.3 * sum(N(j) for j in 1:5)
    Hd = sum(-(gc_op(sys, dag(C), j) * gc_op(sys, C, j + 1) + gc_op(sys, dag(C), j + 1) * gc_op(sys, C, j)) +
             0.4 * (gc_op(sys, C, j) * gc_op(sys, C, j + 1) + gc_op(sys, dag(C), j + 1) * gc_op(sys, dag(C), j)) for j in 1:4) +
         0.3 * sum(gc_op(sys, N, j) for j in 1:5)
    E = eigen(Hermitian(Hd))
    g, _ = gs_product(sys)
    r = ground_state(H, g)
    @test r.energy ≈ E.values[1] atol = 1e-12
    @test gs_same(r.state, E.vectors[:, 1], sys)
    sysq = System([Qubit() for _ in 1:5])
    Hq = -sum(Z(j) * Z(j + 1) for j in 1:4) - 0.7 * sum(X(j) for j in 1:5) + 0.3 * Z(1)
    Hqd = -sum(gc_op(sysq, Z, j) * gc_op(sysq, Z, j + 1) for j in 1:4) - 0.7 * sum(gc_op(sysq, X, j) for j in 1:5) + 0.3 * gc_op(sysq, Z, 1)
    Eq = eigen(Hermitian(Hqd))
    gq, _ = gs_product(sysq; axes = Dict(k => :X for k in 1:5))
    rq = ground_state(Hq, gq)
    @test rq.energy ≈ Eq.values[1] atol = 1e-12
    @test gs_same(rq.state, Eq.vectors[:, 1], sysq)
    # refused: a charged system, states to avoid
    gc_ = TensorMixedStates.GaussianState(State{Pure}(System([Fermion(conserve = N) for _ in 1:3]), ["Occ", "Emp", "Emp"]))
    @test_throws ArgumentError ground_state(sum(dag(C)(j) * C(j + 1) + dag(C)(j + 1) * C(j) for j in 1:2), gc_)
    @test_throws ArgumentError ground_state(H, g; orthogonal_to = [r.state])
end

@testset "gates" begin
    sys = System([Qubit() for _ in 1:4])
    m(O, sites...) = TensorMixedStates.matrix(O, (sys.sites[k] for k in sites)...)
    embed2(U, k) = reduce(kron, [j == k ? U : Matrix{ComplexF64}(I, 2, 2) for j in 1:4 if j != k + 1])
    g, ψ = gs_product(sys)
    gates = [(Rx(0.3)(1), gc_embed(sys, 1, m(Rx(0.3), 1))), (Rz(0.7)(2), gc_embed(sys, 2, m(Rz(0.7), 2))),
             (Rxx(0.4)(2, 3) * Ryy(0.4)(2, 3), embed2(m(Rxx(0.4), 2, 3) * m(Ryy(0.4), 2, 3), 2)),
             (Ry(0.5)(4), gc_embed(sys, 4, m(Ry(0.5), 4))), (X(2), gc_op(sys, X, 2)),
             (Z(2) * X(3), gc_op(sys, Z, 2) * gc_op(sys, X, 3))]
    for (G_, Gd) in gates
        g = apply(G_, g)
        ψ = Gd * ψ
        @test gs_same(g, ψ, sys)
    end
    g2 = apply(2 * X(1), g)
    @test TensorMixedStates.trace(g2) ≈ 4 * TensorMixedStates.trace(g)
    # a fermionic gate written as an exponential, on a mixed system
    sys = System([Fermion(), Qubit(), Fermion()])
    g, ψ = gs_product(sys)
    A = -0.3im * (dag(C) ⊗ C + dag(dag(C) ⊗ C))
    g = apply(exp(A)(1, 3), g)
    Ad = -0.3im * (gc_op(sys, dag(C), 1) * gc_op(sys, C, 3) + (gc_op(sys, dag(C), 1) * gc_op(sys, C, 3))')
    @test gs_same(g, exp(Ad) * ψ, sys)
end

@testset "files and conversion to State" begin
    for (sites, axes) in (([Fermion(), Electron(), Qubit(), Fermion()], Dict(3 => :X)),
                          ([Qubit(), Qubit(), Qubit()], Dict(1 => :Y)),
                          ([Electron()], Dict{Int, Symbol}()))
        sys = System(sites)
        g, ψ = gs_product(sys; axes)
        # the vector of the conversion against the state: the same values and entropies
        st = State(g)
        v = TensorMixedStates.ExactState(st).state
        @test abs(dot(v, ψ)) ≈ 1 atol = 1e-10
        @test gs_same(g, v, sys)
        # files keep the covariance, the axes, the parity and the trace
        mktempdir() do dir
            file = joinpath(dir, "g.h5")
            save_state(file, "g", g)
            h = load_state(file, "g")
            @test h isa TensorMixedStates.GaussianState{Pure}
            @test h.cov == g.cov && h.axes == g.axes && h.parities == g.parities && h.trace == g.trace
        end
    end
    # an entangled state: a ground state converted
    sys = System([Fermion() for _ in 1:4])
    H = sum(-(dag(C)(j) * C(j + 1) + dag(C)(j + 1) * C(j)) + 0.4 * (C(j) * C(j + 1) + dag(C)(j + 1) * dag(C)(j)) for j in 1:3)
    g, _ = gs_product(sys)
    r = ground_state(H, g)
    v = TensorMixedStates.ExactState(State(r.state)).state
    @test gs_same(r.state, v, sys)
    @test entanglement_entropy(r.state, 2) > 0.1
end

# the numbers of the files of measurements of two runs, which must have the same files and lines
function gs_same_outputs(a, b; atol = 1e-10)
    files = sort(filter(f -> endswith(f, ".dat"), readdir(a)))
    ok = !isempty(files) && files == sort(filter(f -> endswith(f, ".dat"), readdir(b)))
    for f in files
        la, lb = readlines(joinpath(a, f)), readlines(joinpath(b, f))
        ok &= length(la) == length(lb)
        for (x, y) in zip(la, lb)
            ta, tb = split(x), split(y)
            ok &= length(ta) == length(tb)
            for (u, v) in zip(ta, tb)
                p, q = tryparse(Float64, u), tryparse(Float64, v)
                ok &= isnothing(p) ? u == v : !isnothing(q) && isapprox(p, q; atol)
            end
        end
    end
    ok
end

# a program run on GaussianState and on ExactState: the same measurements, the same final state
function gs_programs_agree(phases; kwargs...)
    ok = true
    in_tmp() do
        g = runTMS(SimData(; name = "g", backend = TensorMixedStates.GaussianState, phases, kwargs...))
        e = runTMS(SimData(; name = "e", backend = TensorMixedStates.ExactState, phases, kwargs...))
        ok &= g.state isa TensorMixedStates.GaussianState
        ok &= gs_same_outputs("g", "e")
        v, w = TensorMixedStates.ExactState(State(g.state)).state, e.state.state
        ok &= isapprox(abs(dot(v, w)), LA.norm(v) * LA.norm(w); atol = 1e-10) && isapprox(LA.norm(v), LA.norm(w); rtol = 1e-10)
    end
    ok
end

@testset "programs against ExactState" begin
    # fermions: product state, evolution with pairing, gates, ground state
    for _ in 1:2
        n = 6
        sys = System(n, Fermion())
        t, Δ, μ = randn(gs_rng, 3)
        H = sum(t * (dag(C)(j) * C(j + 1) + dag(C)(j + 1) * C(j)) + Δ * (C(j) * C(j + 1) + dag(C)(j + 1) * dag(C)(j)) for j in 1:n-1) +
            μ * sum(N(j) for j in 1:n)
        meas = ["f.dat" => [N, dag(C)(1) * C(4), EntanglementEntropy(3)]]
        phases = [CreateState{Pure}(sys, [gs_random_vector(2) for _ in 1:n]; final_measurements = meas),
                  Evolve(evolver = -im * H, duration = 0.6, time_step = 0.2, measurements = meas),
                  Gates(gates = exp(-0.4im * (dag(C) ⊗ C + dag(dag(C) ⊗ C)))(2, 3) * (C + dag(C))(5), final_measurements = meas),
                  GroundState(hamiltonian = H + 0.3 * N(1), algo = Dmrg(nsweeps = 10), final_measurements = meas)]
        @test gs_programs_agree(phases)
    end
    # qubits along Z: Ising J XX - h Z from X+, gates, XX + YY ground state
    for _ in 1:2
        n = 6
        sys = System(n, Qubit())
        J, h = randn(gs_rng, 2)
        H = J * sum(X(j) * X(j + 1) for j in 1:n-1) - h * sum(Z(j) for j in 1:n)
        meas = ["q.dat" => [X, Z, X(1) * X(2), Y(2) * Y(3), EntanglementEntropy(3)]]
        phases = [CreateState{Pure}(sys, "X+"; final_measurements = meas),
                  Evolve(evolver = -im * H, duration = 0.6, time_step = 0.2, measurements = meas),
                  Gates(gates = Rz(0.3)(2) * Rxx(0.2)(3, 4) * Ryy(0.2)(3, 4) * X(1) * Rx(0.5)(1), final_measurements = meas),
                  GroundState(hamiltonian = sum(X(j) * X(j + 1) + Y(j) * Y(j + 1) for j in 1:n-1) + 0.3 * Z(2) - 0.1 * Z(5),
                              algo = Dmrg(nsweeps = 10), final_measurements = meas)]
        @test gs_programs_agree(phases)
    end
    # qubits then fermions, the last qubit coupled to the fermions: quadratic in this order only (with
    # the fermions first, σ⁺c is the parity of the fermions times a hopping)
    for _ in 1:2
        sys = System([Qubit(), Qubit(), Fermion(), Fermion()])
        a, b, c = randn(gs_rng, 3)
        H = a * (dag(C)(3) * C(4) + dag(C)(4) * C(3)) + b * (Sp(2) * C(4) + dag(C)(4) * Sm(2)) +
            (X(1) * X(2) + Y(1) * Y(2)) + c * (N(3) + Z(1))
        meas = ["m.dat" => [N(3), N(4), Z(1), Z(2), dag(C)(4) * Sm(2), EntanglementEntropy(2)]]
        phases = [CreateState{Pure}(sys, [gs_random_vector(2) for _ in 1:4]; final_measurements = meas),
                  Evolve(evolver = -im * H, duration = 0.6, time_step = 0.2, measurements = meas)]
        @test gs_programs_agree(phases)
    end
    @test_throws ArgumentError TensorMixedStates.gaussian_generator(System([Fermion(), Fermion(), Qubit()]),
                                                                    Sp(3) * C(2) + dag(C)(2) * Sm(3))
end

@testset "checks from the reviews of lot 3" begin
    # a monomial gate multiplies the norm by its coefficient
    sys = System([Spin(1/2), Spin(1/2)])
    g, _ = gs_product(sys)
    @test TensorMixedStates.trace(apply(Sz(1), g)) ≈ TensorMixedStates.trace(g) / 4
    @test TensorMixedStates.trace(apply(Sx(1) * Sz(2), g)) ≈ TensorMixedStates.trace(g) / 16
    # eigenstates of the axes have a definite parity, whose entropies take the Gaussian formula
    g = TensorMixedStates.GaussianState(State{Pure}(System(30, Qubit()), "X+"); axes = :X)
    @test TensorMixedStates.definite(g) && abs(entanglement_entropy(g, 15)) < 1e-12
    # the forms of the axes, and the sites without modes refused
    st = State{Pure}(System([Qubit(), Fermion(), Qubit()]), [gs_random_vector(2) for _ in 1:3])
    g1 = TensorMixedStates.GaussianState(st; axes = :X)
    @test g1.cov == TensorMixedStates.GaussianState(st; axes = [:X, nothing, :X]).cov ==
          TensorMixedStates.GaussianState(st; axes = Dict(1 => :X, 3 => :X)).cov
    @test_throws ArgumentError TensorMixedStates.GaussianState(st; axes = [:X, :X])
    @test_throws ArgumentError TensorMixedStates.GaussianState(State{Pure}(System([Qubit(), Boson(3)]), ["Up", "0"]))
    # values of no definite flux on a charged system, as on a State
    gq = TensorMixedStates.GaussianState(State{Pure}(System([Fermion(conserve = N) for _ in 1:3]), ["Occ", "Emp", "Emp"]))
    @test TensorMixedStates.expect_value(gq, C(1) + dag(C)(1)) ≈ 0 atol = 1e-14
    # truncate, the cuts
    @test truncate(g1; limits = Limits(maxdim = 2)) === g1
    @test_throws ErrorException entanglement_entropy(g1, 0)
    @test_throws ErrorException entanglement_entropy(g1, 4)
    # a Hermitian operator given by a matrix with complex entries, whose expansion leaves rounding
    sys = System([Qubit() for _ in 1:3])
    A = Operator{1}("A", ComplexF64[0.25 0.5+0.25im; 0.5-0.25im -0.75])
    H = A(1) + X(1) * X(2) + Y(1) * Y(2) + 0.3 * Z(3) + X(2) * X(3)
    Hd = gc_op(sys, A, 1) + gc_op(sys, X, 1) * gc_op(sys, X, 2) + gc_op(sys, Y, 1) * gc_op(sys, Y, 2) +
         0.3 * gc_op(sys, Z, 3) + gc_op(sys, X, 2) * gc_op(sys, X, 3)
    g, ψ = gs_product(sys)
    E = eigen(Hermitian(Hd))
    r = ground_state(H, g)
    @test r.energy ≈ E.values[1] atol = 1e-12
    pre = TensorMixedStates.prepared_evolver(g, -im * H)
    g2 = TensorMixedStates.evolve_step(Tdvp(), pre, g, 0.0, 0.3, 1, Limits())
    @test gs_same(g2, exp(-0.3im * Hd) * ψ, sys)
    # the conversion to State is refused beyond 16 modes
    @test_throws ArgumentError State(TensorMixedStates.GaussianState(State{Pure}(System(17, Qubit()), "Up")))
end

@testset "superpositions" begin
    # gates that are not Gaussian make several components: against the dense state, values, entropies
    sys = System([Qubit() for _ in 1:4])
    g, ψ = gs_product(sys)
    m(O, sites...) = TensorMixedStates.matrix(O, (sys.sites[k] for k in sites)...)
    embed2(U, k) = reduce(kron, [j == k ? U : Matrix{ComplexF64}(I, 2, 2) for j in 1:4 if j != k + 1])
    # T, H and S·Rz multiply the components by 2 each, CNOT and CZ by 4 (D4: no merging): 8 at most here
    for (G_, Gd) in ((CNOT(1, 2), embed2(m(CNOT, 1, 2), 1)), (0.5 * CZ(3, 4), 0.5 * embed2(m(CZ, 3, 4), 3)))
        g1, ψ1 = gs_product(sys)
        g1 = apply(G_, g1)
        @test gs_same(g1, Gd * ψ1, sys)
        @test TensorMixedStates.trace(g1) ≈ real(dot(Gd * ψ1, Gd * ψ1)) rtol = 1e-10
    end
    steps = [(T(2), gc_embed(sys, 2, m(T, 2))), (H(3), gc_embed(sys, 3, m(H, 3))),
             (Rxx(0.4)(2, 3), embed2(m(Rxx(0.4), 2, 3), 2)),
             (S(4) * Rz(0.3)(1), gc_embed(sys, 4, m(S, 4)) * gc_embed(sys, 1, m(Rz(0.3), 1)))]
    for (G_, Gd) in steps
        g = apply(G_, g)
        ψ = Gd * ψ
        @test gs_same(g, ψ, sys)
        @test TensorMixedStates.trace(g) ≈ real(dot(ψ, ψ)) rtol = 1e-10
    end
    @test TensorMixedStates.ncomponents(g) == 8
    # an evolution of a superposition, its files, its conversion to State
    Hxy = sum(X(j) * X(j + 1) + Y(j) * Y(j + 1) for j in 1:3) + 0.4 * Z(2)
    Hd = sum(gc_op(sys, X, j) * gc_op(sys, X, j + 1) + gc_op(sys, Y, j) * gc_op(sys, Y, j + 1) for j in 1:3) + 0.4 * gc_op(sys, Z, 2)
    pre = TensorMixedStates.prepared_evolver(g, -im * Hxy)
    for step in 1:3
        g = TensorMixedStates.evolve_step(Tdvp(), pre, g, 0.0, 0.2, step, Limits())
    end
    ψ = exp(-0.6im * Hd) * ψ
    @test gs_same(g, ψ, sys)
    v = TensorMixedStates.ExactState(State(g)).state
    @test abs(dot(v, ψ)) ≈ LA.norm(v) * LA.norm(ψ) rtol = 1e-10
    mktempdir() do dir
        file = joinpath(dir, "g.h5")
        save_state(file, "g", g)
        h = load_state(file, "g")
        @test gs_same(h, ψ, sys) && TensorMixedStates.ncomponents(h) == TensorMixedStates.ncomponents(g)
    end
    # fermions: a quartic operator applied, a product of odd operators
    sys = System([Fermion(), Electron(), Fermion()])
    g, ψ = gs_product(sys)
    g = apply(Nupdn(2) + 0.3 * dag(C)(1) * C(3), g)
    ψ = (gc_op(sys, Nupdn, 2) + 0.3 * gc_op(sys, dag(C), 1) * gc_op(sys, C, 3)) * ψ
    @test gs_same(g, ψ, sys)
    # max_terms: beyond, a clear refusal
    sys = System([Qubit() for _ in 1:4])
    g = TensorMixedStates.GaussianState(State{Pure}(sys, "Up"); max_terms = 4)
    g = apply(H(2) * H(3), g)
    @test TensorMixedStates.ncomponents(g) == 4
    @test_throws ArgumentError apply(H(1), g)
end

# the projected vector, normalised
gs_projected(P, ψ) = (φ = P * ψ; φ / LA.norm(φ))

@testset "measurements" begin
    # projections on occupations: exact, no new component, one or several components
    sys = System([Fermion(), Electron(), Qubit()])
    g, ψ = gs_product(sys)
    h = TensorMixedStates.project_on(g, N(1))
    @test TensorMixedStates.ncomponents(h) == 1
    @test gs_same(h, gs_projected(gc_op(sys, N, 1), ψ), sys)
    # the levels of an Electron, products of the occupations of its two modes
    up = TensorMixedStates.matrix(Nup, sys.sites[2]) * (I - TensorMixedStates.matrix(Ndn, sys.sites[2]))
    r = collapse(g, Level(2) => 1)
    P = gc_embed(sys, 2, up)
    @test r.probability ≈ real(dot(ψ, P * ψ)) / real(dot(ψ, ψ)) atol = 1e-12
    @test TensorMixedStates.ncomponents(r.state) == 1
    @test gs_same(r.state, gs_projected(P, ψ), sys)
    # a superposition projected on an occupation
    g2 = apply(T(3) * H(3), g)
    ψ2 = gc_embed(sys, 3, TensorMixedStates.matrix(T, sys.sites[3]) * TensorMixedStates.matrix(H, sys.sites[3])) * ψ
    h2 = TensorMixedStates.project_on(g2, (Id(1) - Z(3)) / 2)
    @test TensorMixedStates.ncomponents(h2) == TensorMixedStates.ncomponents(g2)
    @test gs_same(h2, gs_projected((I - gc_op(sys, Z, 3)) / 2, ψ2), sys)
    # another projector, as a polynomial
    h3 = TensorMixedStates.project_on(g, (Id(1) + X(3)) / 2)
    @test gs_same(h3, gs_projected((I + gc_op(sys, X, 3)) / 2, ψ), sys)
    # probabilities and a collapse on a result
    pr = probabilities(g, N(1) + N(3))
    for (λ, p) in pr
        Pλ = sum(P for (v, P) in ((0, (I - gc_op(sys, N, 1)) * (I - gc_op(sys, N, 3))),
                                  (1, gc_op(sys, N, 1) * (I - gc_op(sys, N, 3)) + (I - gc_op(sys, N, 1)) * gc_op(sys, N, 3)),
                                  (2, gc_op(sys, N, 1) * gc_op(sys, N, 3))) if v == round(Int, λ))
        @test p ≈ real(dot(ψ, Pλ * ψ)) / real(dot(ψ, ψ)) atol = 1e-12
    end
    @test_throws ArgumentError TensorMixedStates.project_count(g, [1, 3], nothing, 2; limits = Limits())
    # samples: a Fock state gives its occupations, any state values in range
    sysf = System([Fermion() for _ in 1:4])
    gf = TensorMixedStates.GaussianState(State{Pure}(sysf, ["Occ", "Emp", "Occ", "Occ"]))
    @test sample(gf) == [1, 0, 1, 1]
    gr, _ = gs_product(System([Fermion(), Qubit(), Electron()]))
    s = sample(apply(H(2), gr))
    @test length(s) == 3 && 0 <= s[1] <= 1 && 0 <= s[2] <= 1 && 0 <= s[3] <= 3
end

@testset "evolution depending on time, sectors of the parity" begin
    # H(t) = XY + cos(2t) Σ Z, against the evolution of an ExactState (its own Magnus integrator)
    sys = System([Qubit() for _ in 1:4])
    g, _ = gs_product(sys)
    e = TensorMixedStates.ExactState(State(g))
    terms = [-im * sum(X(j) * X(j + 1) + Y(j) * Y(j + 1) for j in 1:3), -im * sum(Z(j) for j in 1:4) => t -> cos(2t)]
    g = evolve(terms, 1.2, g; nsteps = 3)
    e = evolve(terms, 1.2, e; nsteps = 3)
    @test gs_same(g, e.state, sys)
    @test_throws ArgumentError evolve([-im * Z(1), -im * Z(2) => t -> exp(im * t)], 0.3, g; nsteps = 1)
    # the quench of ising_quench, pure: a ring along X from X+, against ExactState
    n = 6
    sys = System(n, Qubit())
    ring = sum(Z(j) * Z(j + 1) for j in 1:n-1) + Z(n) * Z(1)
    H = ring - 0.7 * sum(X(j) for j in 1:n)
    g = TensorMixedStates.GaussianState(State{Pure}(sys, "X+"); axes = :X)
    @test TensorMixedStates.definite(g)
    e = TensorMixedStates.ExactState(State{Pure}(sys, "X+"))
    g = evolve(-im * H, 1.0, g; nsteps = 5)
    e = evolve(-im * H, 1.0, e; nsteps = 5)
    @test gs_same(g, e.state, sys)
    @test TensorMixedStates.expect_value(g, Z(1) * Z(4)) ≈ TensorMixedStates.expect_value(e, Z(1) * Z(4)) atol = 1e-10
    # its ground state, from the two sectors, against the lowest level
    Hd = sum(gc_op(sys, Z, j) * gc_op(sys, Z, mod1(j + 1, n)) for j in 1:n) - 0.7 * sum(gc_op(sys, X, j) for j in 1:n)
    E = eigen(Hermitian(Hd))
    r = ground_state(H, g)
    @test r.energy ≈ E.values[1] atol = 1e-10
    @test gs_same(r.state, E.vectors[:, 1], sys)
    # a ring on a state of no definite parity: its two sectors under their own Gaussian unitaries
    vs = [gs_random_vector(2) for _ in 1:n]
    gm = TensorMixedStates.GaussianState(State{Pure}(sys, vs); axes = :X)
    em = TensorMixedStates.ExactState(State{Pure}(sys, vs))
    gm = evolve(-im * H, 1.0, gm; nsteps = 4)
    em = evolve(-im * H, 1.0, em; nsteps = 4)
    @test TensorMixedStates.ncomponents(gm) == 2 && gm.parities == [1, -1]
    @test gs_same(gm, em.state, sys)
    @test TensorMixedStates.expect_value(gm, Z(2)) ≈ TensorMixedStates.expect_value(em, Z(2)) atol = 1e-10
    # depending on time, refused for now
    @test_throws ArgumentError evolve([-im * H, -im * X(1) => t -> cos(t)],
                                      0.1, TensorMixedStates.GaussianState(State{Pure}(sys, vs); axes = :X); nsteps = 1)
end

@testset "programs with superpositions against ExactState" begin
    # qubits: gates that are not Gaussian, an evolution depending on time, a projection, from X+
    for _ in 1:2
        n = 5
        sys = System(n, Qubit())
        a, b = randn(gs_rng, 2)
        meas = ["s.dat" => [X, Z, X(1) * X(2), Y(2) * Y(4), EntanglementEntropy(2)]]
        evolver = [-im * sum(X(j) * X(j + 1) + Y(j) * Y(j + 1) for j in 1:n-1), -im * a * sum(Z(j) for j in 1:n) => t -> cos(b * t)]
        phases = [CreateState{Pure}(sys, "X+"; final_measurements = meas),
                  Gates(gates = T(2) * H(3), final_measurements = meas),
                  Evolve(evolver = evolver, duration = 0.6, time_step = 0.2, measurements = meas),
                  Gates(gates = (Id(1) + Z(4)) / 2 * S(5), final_measurements = meas)]
        @test gs_programs_agree(phases)
    end
    # fermions: a polynomial applied, pairing depending on time, a projection on an occupation
    for _ in 1:2
        n = 5
        sys = System(n, Fermion())
        a, b = randn(gs_rng, 2)
        meas = ["f.dat" => [N, dag(C)(1) * C(3), EntanglementEntropy(2)]]
        evolver = [-im * sum(dag(C)(j) * C(j + 1) + dag(C)(j + 1) * C(j) for j in 1:n-1),
                   -im * a * sum(C(j) * C(j + 1) + dag(C)(j + 1) * dag(C)(j) for j in 1:n-1) => t -> sin(b * t)]
        phases = [CreateState{Pure}(sys, [gs_random_vector(2) for _ in 1:n]; final_measurements = meas),
                  Gates(gates = dag(C)(1) * C(3) + 0.5 * N(2), final_measurements = meas),
                  Evolve(evolver = evolver, duration = 0.6, time_step = 0.2, measurements = meas),
                  Gates(gates = N(4), final_measurements = meas)]
        @test gs_programs_agree(phases)
    end
end

@testset "trajectories with projective and linear jumps" begin
    # the no-jump evolution not unitary (Σ L†L not constant): monitored occupations, losses and gains of
    # fermions, losses on qubits behind their strings, against ExactState{Trajectory} with the same draws
    function both(sys, st, ev, T; seed = 1, nsteps = 8)
        g = TensorMixedStates.GaussianState(State{Trajectory}(sys, st))
        e = TensorMixedStates.ExactState(State{Trajectory}(sys, st))
        jg, je = [], []
        for i in 1:3
            g = evolve(ev, T / 3, g; nsteps, time = (i - 1) * T / 3, rng = Random.Xoshiro(seed + i), jumps = jg)
            e = evolve(ev, T / 3, e; nsteps, time = (i - 1) * T / 3, rng = Random.Xoshiro(seed + i), jumps = je)
        end
        g, e, jg, je
    end
    fs = System(4, Fermion())
    hop = sum(dag(C)(j) * C(j + 1) + dag(C)(j + 1) * C(j) for j in 1:3) + 0.3 * N(2)
    qs = System(5, Qubit())
    xy = sum(X(j) * X(j + 1) + Y(j) * Y(j + 1) for j in 1:4) + 0.4 * sum(Z(j) for j in 1:5)
    cases = [(fs, ["Occ", "Emp", "Occ", "Emp"], -im * hop + sum(Dissipator(N)(j) for j in 1:4), [N(1), N(3), C(1) * dag(C)(2)]),
             (fs, ["Occ", "Occ", "Occ", "Emp"], -im * hop + sum(0.5 * Dissipator(C)(j) for j in 1:4), [N(2), N(4), C(2) * dag(C)(3)]),
             (fs, "Emp", -im * hop + 0.4 * Dissipator(dag(C))(2) + 0.6 * Dissipator(C)(3) + Dissipator(N)(1), [N(1), N(2), N(4)]),
             (qs, ["Dn", "Up", "Dn", "Dn", "Up"], -im * xy + sum(0.5 * Dissipator(Sp)(j) for j in 1:5) + 0.3 * Dissipator(N)(4),
              [Z(1), Z(3), X(2) * X(3), Z(5)])]
    for (sys, st, ev, obs) in cases, seed in (1, 2)
        g, e, jg, je = both(sys, st, ev, 1.2; seed)
        @test length(jg) == length(je)
        @test all(a.time ≈ b.time && a.kraus == b.kraus for (a, b) in zip(jg, je))
        @test TensorMixedStates.trace(g) ≈ TensorMixedStates.trace(e) rtol = 1e-9
        for o in obs
            @test TensorMixedStates.expect_value(g, o) ≈ TensorMixedStates.expect_value(e, o) atol = 1e-9
        end
    end
    # the jumps of the review of 5c: a projector with its rate, a gate exp(A) not unitary, a null jump, rates
    # depending on time, sites of both kinds, one site; unitary jumps that are not Gaussian when the evolution
    # without jumps is unitary (the polynomial way, the components multiplied)
    mixed = System([Fermion(), Qubit(), Fermion()])
    more = [(fs, ["Occ", "Emp", "Occ", "Emp"], -im * hop + 0.5 * Dissipator(Id - N)(2) + Dissipator(N)(3), [N(1), N(2), N(4)]),
            (fs, ["Occ", "Emp", "Occ", "Emp"], -im * hop + 3.0 * Dissipator(exp(0.5 * N))(2) + Dissipator(C)(1), [N(1), N(3)]),
            (fs, ["Occ", "Emp", "Occ", "Emp"], -im * hop + Dissipator(C * C)(1) + Dissipator(C)(3), [N(2), N(3)]),
            (fs, ["Occ", "Occ", "Emp", "Emp"], [-im * hop, sum(Dissipator(N)(j) for j in 1:4) => t -> 1 + 0.5 * sin(3t)],
             [N(1), N(4), C(2) * dag(C)(3)]),
            (mixed, ["Occ", "Dn", "Occ"], -im * (N(1) + 0.7 * Z(2) + N(3)) + Dissipator(C)(1) + 0.5 * Dissipator(Sp)(2) +
             Dissipator(C)(3), [N(1), Z(2), N(3)]),
            (System(1, Fermion()), "Occ", -im * N(1) + 0.8 * Dissipator(C)(1), [N(1)]),
            (qs, "X+", -im * xy + 0.4 * Dissipator(T)(2), [X(1), Z(2), X(3) * X(4)])]
    for (sys, st, ev, obs) in more
        g, e, jg, je = both(sys, st, ev, 1.2)
        @test length(jg) == length(je)
        @test all(a.time ≈ b.time for (a, b) in zip(jg, je))
        @test TensorMixedStates.trace(g) ≈ TensorMixedStates.trace(e) rtol = 1e-9
        for o in obs
            @test TensorMixedStates.expect_value(g, o) ≈ TensorMixedStates.expect_value(e, o) atol = 1e-9
        end
    end
    # long steps under strong losses, one step: the vacuum stays empty, a fermion lost as in ExactState
    g = evolve(-im * hop + 10 * Dissipator(C)(1), 10.0, TensorMixedStates.GaussianState(State{Trajectory}(fs, "Emp"));
               nsteps = 1, rng = Random.Xoshiro(1))
    @test TensorMixedStates.trace(g) ≈ 1 atol = 1e-12
    g, e, jg, je = both(fs, ["Emp", "Emp", "Emp", "Occ"], -im * hop + 10 * Dissipator(C)(1), 30.0; nsteps = 1)
    @test length(jg) == length(je)
    @test TensorMixedStates.trace(g) ≈ TensorMixedStates.trace(e) rtol = 1e-8
    @test TensorMixedStates.expect_value(g, N(2)) ≈ TensorMixedStates.expect_value(e, N(2)) atol = 1e-9
    # a linear operator on a state it nearly annihilates, exact by the projection: Sp after Ry(1e-10)
    q1 = System(2, Qubit())
    g = apply(Sp(1), apply(Ry(1e-10)(1), TensorMixedStates.GaussianState(State{Pure}(q1, "Up"))))
    @test abs(TensorMixedStates.expect_value(g, X(1))) < 1e-14
    # refused: a jump that does not keep the state Gaussian under an evolution that is not unitary, a
    # superposition under it, a state of both sectors of the parity on a ring under it
    gf = TensorMixedStates.GaussianState(State{Trajectory}(fs, ["Occ", "Emp", "Occ", "Emp"]))
    @test_throws "keep a state Gaussian" TensorMixedStates.prepared_evolver(gf, [-im * hop + Dissipator(Id + 0.3 * C * dag(C))(1) +
                                                                                Dissipator(N)(2)])
    gq = apply(T(2), TensorMixedStates.GaussianState(State{Trajectory}(qs, "X+")))
    @test TensorMixedStates.ncomponents(gq) == 2
    @test_throws "superposition" evolve(-im * xy + Dissipator(N)(1), 0.2, gq; nsteps = 1)
    rs = System(4, Qubit())
    ring = sum(Z(j) * Z(j + 1) for j in 1:3) + Z(4) * Z(1) - sum(X(j) for j in 1:4)
    gr = TensorMixedStates.GaussianState(State{Trajectory}(rs, "Dn"); axes = :X)
    @test_throws "each sector" evolve(-im * ring + Dissipator((Id - X) / 2)(1), 0.2, gr; nsteps = 1)
end

@testset "thermalize of Gaussian trajectories" begin
    # e^{-βH/2}ψ normalised and log ‖·‖², against ExactState{Trajectory}: fermions, qubits, a ring of qubits from a
    # state of definite parity (compiled in its sector), and a thermalization in two calls
    fs = System(4, Fermion())
    hf = sum(dag(C)(j) * C(j + 1) + dag(C)(j + 1) * C(j) for j in 1:3) + 0.4 * N(1) - 0.3 * N(4)
    qs = System(4, Qubit())
    hq = sum(X(j) * X(j + 1) + 0.6 * Y(j) * Y(j + 1) for j in 1:3) + 0.5 * sum(Z(j) for j in 1:4)
    hr = sum(X(j) * X(mod1(j + 1, 4)) + Y(j) * Y(mod1(j + 1, 4)) for j in 1:4) + 0.3 * sum(Z(j) for j in 1:4)
    for (sys, st, h, obs) in ((fs, ["Occ", "Emp", "Occ", "Emp"], hf, [N(1), N(3), C(1) * dag(C)(2)]),
                              (qs, "X+", hq, [X(1), Z(2), X(2) * X(3)]),
                              (qs, ["Up", "Dn", "Up", "Dn"], hr, [Z(1), X(1) * X(2), Y(2) * Y(3)]))
        g = thermalize(h, 1.5, TensorMixedStates.GaussianState(State{Trajectory}(sys, st)); nsteps = 3)
        e = thermalize(h, 1.5, TensorMixedStates.ExactState(State{Trajectory}(sys, st)); nsteps = 3)
        @test g.log_trace ≈ e.log_trace atol = 1e-9
        @test TensorMixedStates.trace(g.state) ≈ 1 atol = 1e-12
        for o in obs
            @test TensorMixedStates.expect_value(g.state, o) ≈ TensorMixedStates.expect_value(e.state, o) atol = 1e-9
        end
        g1 = thermalize(h, 0.75, TensorMixedStates.GaussianState(State{Trajectory}(sys, st)); nsteps = 1)
        g2 = thermalize(h, 0.75, g1.state; nsteps = 1, step = 2)
        @test g1.log_trace + g2.log_trace ≈ g.log_trace atol = 1e-9
    end
    # refused: a pure state, a superposition
    @test_throws "trajectories" thermalize(hf, 1.0, TensorMixedStates.GaussianState(State{Pure}(fs, "Emp")); nsteps = 1)
    gq = apply(T(2), TensorMixedStates.GaussianState(State{Trajectory}(qs, "X+")))
    @test_throws "superposition" thermalize(hq, 1.0, gq; nsteps = 1)
end

@testset "cross terms of nearly orthogonal components" begin
    # X+ along Z has ⟨Z₁⟩ = 0: T(1) makes G and Z₁G, orthogonal by structure (C of `pair_transition`); Ry(θ)
    # on qubit 1 first makes them nearly orthogonal (θ = 1e-7, the precise way) or not (θ = 0.3, the
    # transition matrix)
    sys = System([Qubit() for _ in 1:4])
    m(O, k) = gc_embed(sys, k, TensorMixedStates.matrix(O, sys.sites[k]))
    for θ in (0.0, 1e-7, 0.3)
        g = TensorMixedStates.GaussianState(State{Pure}(sys, "X+"))
        ψ = reduce(kron, fill(ComplexF64[1, 1] / sqrt(2), 4))
        g = apply(T(1) * H(2) * Ry(θ)(1), g)
        ψ = m(T, 1) * m(H, 2) * m(Ry(θ), 1) * ψ
        @test gs_same(g, ψ, sys)
        P = TensorMixedStates.pair_transitions(g.sum)
        pairs = [P[k, l] for k in axes(P, 1), l in axes(P, 2) if k < l]
        if θ == 0.3
            @test any(p -> !isnothing(p) && isempty(p[3]), pairs)
        elseif θ == 0
            @test all(!isnothing, pairs) && any(p -> !isempty(p[3]), pairs)
        else
            @test any(isnothing, pairs)
        end
    end
end

@testset "fixes from the reviews of lot 4" begin
    n = 4
    sys = System(n, Qubit())
    q(O, k) = gc_op(sys, O, k)
    ringxy = sum(X(j) * X(mod1(j + 1, n)) + Y(j) * Y(mod1(j + 1, n)) for j in 1:n)
    ringxyd = sum(q(X, j) * q(X, mod1(j + 1, n)) + q(Y, j) * q(Y, mod1(j + 1, n)) for j in 1:n)
    st = State{Pure}(sys, ["Up", "Dn", "Up", "Dn"])
    g = TensorMixedStates.GaussianState(st)
    e = TensorMixedStates.ExactState(st)
    # a ring with a term changing the parity: refused, not reduced to a sector
    @test_throws ArgumentError evolve([-im * ringxy, -im * X(1) => cos], 0.5, g; nsteps = 2)
    # a ring depending on time, the parity kept
    terms = [-im * ringxy, -im * Z(2) => cos]
    @test gs_same(evolve(terms, 0.5, g; nsteps = 2), evolve(terms, 0.5, e; nsteps = 2).state, sys)
    # one evolver for states of the two parities, each in its sector
    pre = TensorMixedStates.prepared_evolver(g, -im * ringxy)
    go = apply(X(1), g)
    go = TensorMixedStates.evolve_step(Tdvp(), pre, go, 0.0, 0.4, 1, Limits())
    @test gs_same(go, exp(-0.4im * ringxyd) * (q(X, 1) * e.state), sys)
    # an evolution depending on time unitary as a whole, not term by term
    sysf = System(3, Fermion())
    gf, ψf = gs_product(sysf)
    ef = TensorMixedStates.ExactState{Pure}(sysf, ψf, Val(:raw))
    tf = [-im * dag(C)(1) * C(2) => t -> cis(t), -im * dag(C)(2) * C(1) => t -> cis(-t)]
    @test gs_same(evolve(tf, 0.6, gf; nsteps = 2), evolve(tf, 0.6, ef; nsteps = 2).state, sysf)
    # null states: refused when measured, as for an ExactState
    g0 = apply(Id(1) - Z(1), TensorMixedStates.GaussianState(State{Pure}(sys, "Up")))
    @test TensorMixedStates.trace(g0) == 0 && TensorMixedStates.norm(g0) == 0
    @test_throws ErrorException TensorMixedStates.expect_value(g0, Z(2))
    @test TensorMixedStates.trace(apply(N(1) - N(1), gf)) == 0
    @test TensorMixedStates.trace(apply(C(2), TensorMixedStates.GaussianState(State{Pure}(sysf, "Emp")))) == 0
    # a projection keeps a definite parity
    @test TensorMixedStates.definite(TensorMixedStates.project_on(g, (Id(1) + Z(1)) / 2))
    # a collapse on a sum whose components conserve the number
    gu = TensorMixedStates.GaussianState(State{Pure}(sys, "Up"))
    eu = TensorMixedStates.ExactState(State{Pure}(sys, "Up"))
    # an open chain: a ring needs a definite parity, which a sum has not
    Hn = sum(X(j) * X(j + 1) + Y(j) * Y(j + 1) for j in 1:n-1) + 0.4 * Z(2)
    gu = evolve(-im * Hn, 0.7, apply(H(1), gu); nsteps = 2)
    eu = evolve(-im * Hn, 0.7, apply(H(1), eu); nsteps = 2)
    rg, re = collapse(gu, 3 => 1), collapse(eu, 3 => 1)
    @test rg.probability ≈ re.probability atol = 1e-12
    @test gs_same(rg.state, re.state.state, sys)
    # gates exp(A), A not quadratic, as polynomials
    gp, ψp = gs_product(sys)
    gp = apply(Rx(0.3)(2) * exp(-0.3im * (Z ⊗ Z))(1, 2), gp)
    ψp = q(Rx(0.3), 2) * exp(-0.3im * q(Z, 1) * q(Z, 2)) * ψp
    @test gs_same(gp, ψp, sys)
    # the log of a degenerate level of the other sector, not returned, is not written
    sys5 = System(5, Qubit())
    H5 = sum(Z(j) * Z(mod1(j + 1, 5)) for j in 1:5) - 0.7 * sum(X(j) for j in 1:5)
    g5 = TensorMixedStates.GaussianState(State{Pure}(sys5, "X+"); axes = :X)
    @test_logs ground_state(H5, g5)
    # sample of one site from the probabilities, without projecting
    gx = apply(T(2), TensorMixedStates.GaussianState(State{Pure}(System(3, Qubit()), "X+"); axes = :X, max_terms = 2))
    @test sample(gx, 1) in (0, 1)
    # the sites of a term refused in a sector
    m = try
        TensorMixedStates.prepared_evolver(TensorMixedStates.GaussianState(State{Pure}(System(6, Qubit()), "Up")),
                                           -im * Z(1) * Z(2) * Z(3) * Z(4))
        ""
    catch err
        err.msg
    end
    @test occursin("sites 1, 2, 3, 4", m)
end

@testset "trajectories with jumps of constant Σ L†L" begin
    # ising_quench: a ring along X from X+, Z jumps, against ExactState with the same draws
    n = 6
    sys = System(n, Qubit())
    H = sum(Z(j) * Z(j + 1) for j in 1:n-1) + Z(n) * Z(1) - 0.7 * sum(X(j) for j in 1:n)
    L = -im * H + sum(0.4 * Dissipator(Z)(j) for j in 1:n)
    for seed in 1:3
        g = TensorMixedStates.GaussianState(State{Trajectory}(sys, "X+"); axes = :X)
        e = TensorMixedStates.ExactState(State{Trajectory}(sys, "X+"))
        @test g isa TensorMixedStates.GaussianState{Trajectory}
        jg, je = [], []
        g = evolve(L, 2.0, g; nsteps = 20, rng = Random.Xoshiro(seed), jumps = jg)
        e = evolve(L, 2.0, e; nsteps = 20, rng = Random.Xoshiro(seed), jumps = je)
        @test [j.dissipator for j in jg] == [j.dissipator for j in je]
        @test all(isapprox(a.time, b.time; atol = 1e-12) for (a, b) in zip(jg, je))
        @test TensorMixedStates.trace(g) ≈ real(TensorMixedStates.trace(e)) atol = 1e-9
        @test gs_same(g, e.state, sys)
    end
    # fermions, Majorana jumps C + C† and σ⁺, σ⁻ pairs of the same rate (each L†L is not constant, their
    # sum is), one depending on time, against ExactState
    sys = System(5, Fermion())
    Hf = sum(dag(C)(j) * C(j + 1) + dag(C)(j + 1) * C(j) + 0.3 * (C(j) * C(j + 1) + dag(C)(j + 1) * dag(C)(j)) for j in 1:4)
    terms = [-im * Hf + 0.3 * Dissipator(C + dag(C))(2), 0.5 * (Dissipator(C)(4) + Dissipator(dag(C))(4)) => t -> 1 + t]
    g, ψ = gs_product(sys)
    e = TensorMixedStates.ExactState(State(g))
    e = TensorMixedStates.ExactState{Trajectory}(sys, e.state, Val(:raw))
    g = TensorMixedStates.as_representation(Trajectory, g)
    jg, je = [], []
    g = evolve(terms, 1.5, g; nsteps = 10, rng = Random.Xoshiro(4), jumps = jg)
    e = evolve(terms, 1.5, e; nsteps = 10, rng = Random.Xoshiro(4), jumps = je)
    @test !isempty(jg) && [j.dissipator for j in jg] == [j.dissipator for j in je]
    @test TensorMixedStates.trace(g) ≈ real(TensorMixedStates.trace(e)) rtol = 1e-9
    @test gs_same(g, e.state, sys)
    # a program: evolution with jumps, a Pauli channel drawn, measurements, against ExactState with the
    # same draws
    sys = System(5, Qubit())
    meas = EachSample(Data("t")) => [X, Z, X(1) * X(2), EntanglementEntropy(2)]
    phases = [CreateState{Trajectory}(sys, "X+"; final_measurements = meas),
              Evolve(evolver = -im * sum(X(j) * X(j + 1) + Y(j) * Y(j + 1) for j in 1:4) + sum(0.3 * Dissipator(Z)(j) for j in 1:5),
                     duration = 1.0, time_step = 0.1, measurements = meas),
              Gates(gates = depolarizing_channel(0.6)(3), final_measurements = meas)]
    for seed in 1:3
        run(b) = runTMS(SimData(; name = "t", backend = b, phases, seed, samples = 1); output = devnull)
        rg, re = in_tmp(() -> run(TensorMixedStates.GaussianState)), in_tmp(() -> run(TensorMixedStates.ExactState))
        @test rg.state isa TensorMixedStates.GaussianState{Trajectory}
        tg, te = rg.data["t"], re.data["t"]
        @test length(tg) == length(te) == 12
        same(p, q) = first(p) == first(q) && isapprox(last(p), last(q); atol = 1e-10)
        @test all(length(a) == length(b) && all(same(p, q) for (p, q) in zip(a, b)) for (a, b) in zip(tg, te))
        v, w = TensorMixedStates.ExactState(State(rg.state)).state, re.state.state
        @test abs(dot(v, w)) ≈ LA.norm(v) * LA.norm(w) atol = 1e-10
    end
    gt = TensorMixedStates.GaussianState(State{Trajectory}(sys, "Up"))
    @test_throws ErrorException TensorMixedStates.as_representation(Pure, gt)
    # files keep the representation
    mktempdir() do dir
        file = joinpath(dir, "t.h5")
        save_state(file, "t", gt)
        @test load_state(file, "t") isa TensorMixedStates.GaussianState{Trajectory}
    end
end

@testset "sectors of the parity of ising_quench" begin
    # from |↓…↓⟩ (no definite parity along X), the ring is refused with its remedy; each sector, projected,
    # has a definite parity, and the mean of the two trajectories drawn with the same generator is the
    # trajectory of ExactState from |↓…↓⟩ for the even observables
    n = 6
    sys = System(n, Qubit())
    J = 1 .+ 0.3 * randn(gs_rng, n)
    H = sum(J[j] * Z(j) * Z(j + 1) for j in 1:n-1) + J[n] * Z(n) * Z(1) - sum(X(j) for j in 1:n)
    ev = -im * H + sum(0.1 * Dissipator(Z)(j) for j in 1:n)
    g = TensorMixedStates.GaussianState(State{Trajectory}(sys, "Dn"); axes = :X)
    @test !TensorMixedStates.definite(g)
    P = prod(X(j) for j in 1:n)
    sectors = [TensorMixedStates.project_on(g, (Id(1) + s * P) / 2) for s in (1, -1)]
    @test all(TensorMixedStates.definite(x) for x in sectors)
    @test [TensorMixedStates.state_parity(x) for x in sectors] == [1, -1]
    for seed in 1:2
        gs = [evolve(ev, 2.0, x; nsteps = 10, rng = Random.Xoshiro(seed)) for x in sectors]
        e = evolve(ev, 2.0, TensorMixedStates.ExactState(State{Trajectory}(sys, "Dn")); nsteps = 10, rng = Random.Xoshiro(seed))
        for o in (X(2), X(1) * X(4), Y(2) * Y(3), Z(1) * Z(4), Z(3) * Z(6))
            @test sum(TensorMixedStates.expect_value(x, o) for x in gs) / 2 ≈ TensorMixedStates.expect_value(e, o) atol = 1e-10
        end
    end
    # the whole state from |↓…↓⟩, its two branches together: the same jumps as ExactState, odd values (⟨Z⟩)
    # included, the parities of the branches exchanged by each jump
    for seed in 1:3
        jg, je = [], []
        x = evolve(ev, 2.0, g; nsteps = 10, rng = Random.Xoshiro(seed), jumps = jg)
        e = evolve(ev, 2.0, TensorMixedStates.ExactState(State{Trajectory}(sys, "Dn")); nsteps = 10,
                   rng = Random.Xoshiro(seed), jumps = je)
        @test [j.dissipator for j in jg] == [j.dissipator for j in je]
        @test x.parities == (iseven(length(jg)) ? [1, -1] : [-1, 1])
        @test gs_same(x, e.state, sys)
        for k in (1, 4)
            @test TensorMixedStates.expect_value(x, Z(k)) ≈ TensorMixedStates.expect_value(e, Z(k)) atol = 1e-9
        end
        @test TensorMixedStates.trace(x) ≈ real(TensorMixedStates.trace(e)) atol = 1e-9
    end
    # files keep the parities of a sum
    x = evolve(-im * H, 0.3, TensorMixedStates.GaussianState(State{Pure}(sys, "Dn"); axes = :X); nsteps = 1)
    mktempdir() do dir
        file = joinpath(dir, "b.h5")
        save_state(file, "b", x)
        @test load_state(file, "b").parities == x.parities == [1, -1]
    end
    # a ground state with zero modes in the system and the ancilla decoupled: its parity is that of the
    # ancilla, the entropies by the covariance
    sysf = System(4, Fermion())
    Hf = sum(-(dag(C)(j) * C(j + 1) + dag(C)(j + 1) * C(j)) for j in 2:3)
    r = ground_state(Hf, TensorMixedStates.GaussianState(State{Pure}(sysf, "Emp")))
    σ = only(r.state.parities)
    @test r.state.cov[1, 2] == σ
    v = TensorMixedStates.ExactState(State(r.state)).state
    @test gs_same(r.state, v, sysf)
end

@testset "fixes from the reviews of lot 5" begin
    # ground states with a term in the parity, which lands on the ancilla: energy and state against the dense
    for (sys, axes, H) in ((System(4, Qubit()), :X, -sum(X(j) for j in 1:4) + 0.5 * prod(X(j) for j in 1:4)),
                           (System(4, Qubit()), :X, -sum(X(j) for j in 1:4) + 2 * prod(X(j) for j in 1:4)),
                           (System(2, Qubit()), :Z, Z(1) * Z(2)),
                           (System(2, Qubit()), :Z, Z(1) * Z(2) + 0.3 * X(1) * X(2)))
        g = TensorMixedStates.GaussianState(State{Pure}(sys, "Up"); axes)
        e = TensorMixedStates.ExactState(State{Pure}(sys, "Up"))
        r, re = ground_state(H, g), ground_state(H, e)
        @test r.energy ≈ re.energy atol = 1e-10
        @test real(TensorMixedStates.expect_value(r.state, H)) ≈ re.energy atol = 1e-10
    end
    # a ring on a sum of unknown parity, its components all in one sector (the other sector empty), and on
    # states of one parity whose parity was lost by gates coupling the ancilla (the other sector noise)
    n = 5
    sys = System(n, Qubit())
    ring = sum(Z(j) * Z(j + 1) for j in 1:n-1) + Z(n) * Z(1) - 0.6 * sum(X(j) for j in 1:n)
    start = State{Pure}(sys, ["Up"; fill("X+", n - 1)])
    prep(st) = apply(Z(2) * Z(3) + 0.5 * Z(4) * Z(5), TensorMixedStates.project_on(st, 0.5 * (Id(1) + X(1)); limits = Limits()))
    g = prep(TensorMixedStates.GaussianState(start; axes = :X))
    e = prep(TensorMixedStates.ExactState(start))
    @test isnothing(g.parities)
    g, e = evolve(-im * ring, 0.5, g; nsteps = 2), evolve(-im * ring, 0.5, e; nsteps = 2)
    @test gs_same(g, e.state, sys)
    xp = State{Pure}(sys, "X+")
    for gates in ((exp(-im * π / 2 * Z)(1),), (exp(-0.37im * Z)(1), exp(0.37im * Z)(1)))
        g = foldl((st, x) -> apply(x, st), gates; init = TensorMixedStates.GaussianState(xp; axes = :X))
        e = foldl((st, x) -> apply(x, st), gates; init = TensorMixedStates.ExactState(xp))
        g, e = evolve(-im * ring, 0.5, g; nsteps = 2), evolve(-im * ring, 0.5, e; nsteps = 2)
        @test TensorMixedStates.ncomponents(g) == 1
        @test gs_same(g, e.state, sys)
    end
    # the rotations kept by a trajectory: one per parity, whatever its jumps
    gt = TensorMixedStates.GaussianState(State{Trajectory}(sys, "X+"); axes = :X)
    pre = TensorMixedStates.with_rng(TensorMixedStates.prepared_evolver(gt, -im * ring + sum(0.5 * Dissipator(Z)(j) for j in 1:n)),
                                     Random.Xoshiro(1), [])
    for k in 0:39
        gt = TensorMixedStates.evolve_step(Tdvp(), pre, gt, 0.1k, 0.1, k + 1, Limits())
    end
    @test !isempty(pre.jumps) && length(pre.pre.rotations) <= 2
    # refusals with their remedy
    gp = TensorMixedStates.GaussianState(xp; axes = :X)
    @test_throws "is a generator" apply(Dissipator(Z)(1), gp)
    @test_throws "is a generator" apply(Dissipator(Z)(1), TensorMixedStates.as_representation(Trajectory, gp))
    @test_throws "draw it on a trajectory" apply(depolarizing_channel(0.2)(1), gp)
    @test_throws "not searched" ground_state(ring, TensorMixedStates.as_representation(Trajectory, gp))
    @test_throws "not a pure or a mixed state" TensorMixedStates.as_representation(Mixed, TensorMixedStates.as_representation(Trajectory, gp))
    gm = TensorMixedStates.GaussianState(State{Pure}(sys, "Up"); axes = :X)
    @test_throws "project it on each sector" TensorMixedStates.prepared_evolver(gm, [-im * ring, -im * X(1) => t -> cos(t)])
end

@testset "cross terms of orthogonal components" begin
    # a sum of 7 components on qubits and fermions, some exactly orthogonal after the evolution (their
    # overlap matrix singular: NaN by the pfaffian of SkewLinearAlgebra 1.1), against ExactState
    sys = System([Qubit(), Qubit(), Fermion(), Fermion()])
    st = State{Pure}(sys, ["Up", "Dn", "Occ", "Emp"])
    H = dag(C)(3) * C(4) + dag(C)(4) * C(3) + 0.3 * Z(1) + 0.2 * X(1) * X(2)
    gate = X(1) + 0.5 * C(3) + 0.2 * dag(C)(4) * C(3)
    g = evolve(-im * H, 0.8, apply(gate, TensorMixedStates.GaussianState(st)); nsteps = 3)
    e = evolve(-im * H, 0.8, apply(gate, TensorMixedStates.ExactState(st)); nsteps = 3)
    @test TensorMixedStates.ncomponents(g) > 1
    for o in (X(2), X(2) * N(4), Z(1), N(3), C(3))
        @test TensorMixedStates.expect_value(g, o) ≈ TensorMixedStates.expect_value(e, o) atol = 1e-10
    end
    @test gs_same(g, e.state, sys)
end

@testset "cross terms of the branches" begin
    # the two branches of a ring, in the two sectors of the ancilla: their pair orthogonal by structure
    # (C of `pair_transition`, one pair of modes), the terms leaving c(α)G_l in the other sector zero, every
    # value as the double sum of the cross terms by overlaps
    GK = TensorMixedStates.GaussianKernels
    n = 6
    sys = System(n, Qubit())
    H = sum(Z(j) * Z(j + 1) for j in 1:n-1) + Z(n) * Z(1) - sum(X(j) for j in 1:n)
    x = evolve(-im * H, 0.7, TensorMixedStates.GaussianState(State{Pure}(sys, "Dn"); axes = :X); nsteps = 2)
    s = x.sum
    @test x.parities == [1, -1]
    @test size(TensorMixedStates.pair_transitions(s)[1, 2][3], 2) == 2
    function by_overlaps(α)
        v = sum(conj(ak) * al * (k == l ? GK.expect_monomial(gk.cov, α) : GK.cross_term(gk, α, gl))
                for (k, (gk, ak)) in enumerate(zip(s.comps, s.coefs)), (l, (gl, al)) in enumerate(zip(s.comps, s.coefs)))
        ComplexF64(v) / x.trace
    end
    for o in (Z(1), Z(3), X(2), Z(1) * Z(4), X(1) * Z(2), Y(2) * Z(5), Id(1))
        p = TensorMixedStates.polynomial(sys, x.layout, TensorMixedStates.lower(TensorMixedStates.weaken(sys, ()), o))
        @test TensorMixedStates.poly_value(x, p) ≈ sum((c * by_overlaps(α) for (α, c) in p.terms); init = 0im) atol = 1e-12
    end
end

@testset "cross terms of components orthogonal by structure" begin
    # T gates on three qubits of X+ along Z, then an XY chain: 8 components, pairs orthogonal by one to three
    # pairs of modes, spread over the chain; every value as the double sum of the cross terms by overlaps
    GK = TensorMixedStates.GaussianKernels
    n = 6
    sys = System(n, Qubit())
    H = sum(X(j) * X(j + 1) + Y(j) * Y(j + 1) for j in 1:n-1) + 0.3 * sum(Z(j) for j in 1:n)
    x = TensorMixedStates.GaussianState(State{Pure}(sys, "X+"))
    for j in (1, 3, 5)
        x = apply(T(j), x)
    end
    x = evolve(-im * H, 2.0, x; nsteps = 1)
    s = x.sum
    GK.fix!(s)
    P = TensorMixedStates.pair_transitions(s)
    widths = [size(P[k, l][3], 2) for k in axes(P, 1), l in axes(P, 2) if k < l && !isnothing(P[k, l])]
    @test length(widths) == 28 && Set(widths) == Set([2, 4, 6])
    function by_overlaps(α)
        v = sum(conj(ak) * al * (k == l ? GK.expect_monomial(gk.cov, α) : GK.cross_term(gk, α, gl))
                for (k, (gk, ak)) in enumerate(zip(s.comps, s.coefs)), (l, (gl, al)) in enumerate(zip(s.comps, s.coefs)))
        ComplexF64(v) / x.trace
    end
    for o in (Z(1), Z(4), X(1) * X(2), X(1) * X(6), Y(2) * Y(5), Z(1) * Z(3) * Z(5), Id(1))
        p = TensorMixedStates.polynomial(sys, x.layout, TensorMixedStates.lower(TensorMixedStates.weaken(sys, ()), o))
        @test TensorMixedStates.poly_value(x, p) ≈ sum((c * by_overlaps(α) for (α, c) in p.terms); init = 0im) atol = 1e-12
    end
end

@testset "parity of a product state judged on the amplitude" begin
    # an amplitude of 6e-7 in the other sector (a probability of 3.6e-13) makes no definite parity: kept by
    # the values and by the branches, against ExactState; an eigenstate of the parity is definite
    sys = System(3, Qubit())
    st = State{Pure}(sys, [ComplexF64[1, 6e-7] / sqrt(1 + 3.6e-13), ComplexF64[1, 0], ComplexF64[1, 0]])
    g = TensorMixedStates.GaussianState(st)
    @test isnothing(g.parities)
    @test TensorMixedStates.definite(TensorMixedStates.GaussianState(State{Pure}(sys, "Up")))
    @test TensorMixedStates.definite(TensorMixedStates.GaussianState(State{Pure}(sys, "X+"); axes = :X))
    gate = 0.5 * Id(1) + X(2)
    x, e = apply(gate, g), apply(gate, TensorMixedStates.ExactState(st))
    for o in (X(1), X(1) * X(2), Z(1))
        @test TensorMixedStates.expect_value(x, o) ≈ TensorMixedStates.expect_value(e, o) atol = 1e-13
    end
    n = 4
    ring = System(n, Qubit())
    H = sum(Z(j) * Z(j + 1) for j in 1:n-1) + Z(n) * Z(1) - sum(X(j) for j in 1:n)
    v = ComplexF64[1 + 7e-7, 1 - 7e-7] / sqrt(2 * (1 + 4.9e-13))
    st = State{Pure}(ring, [[v]; fill(ComplexF64[1, 1] / sqrt(2), n - 1)])
    x = evolve(-im * H, 0.5, TensorMixedStates.GaussianState(st; axes = :X); nsteps = 2)
    e = evolve(-im * H, 0.5, TensorMixedStates.ExactState(st); nsteps = 2)
    for o in (Z(1), Z(1) * Z(3), X(2), Y(1) * Z(2))
        @test TensorMixedStates.expect_value(x, o) ≈ TensorMixedStates.expect_value(e, o) atol = 1e-9
    end
end

@testset "mixed states" begin
    # GaussianState{Mixed} against ExactState{Mixed}: products commuting with the parity and pure states of a
    # definite parity made mixed, then gates, products of Majoranas, projections; values, entropies,
    # measurements, files; the refusals
    TS = TensorMixedStates
    fs = System(4, Fermion())
    qs = System(4, Qubit())
    hop(a, b) = exp(-0.4im * (dag(C) ⊗ C + dag(dag(C) ⊗ C)))(a, b)
    xy(a, b) = exp(-0.3im * (X ⊗ X + 0.6 * Y ⊗ Y))(a, b)
    fgates = [hop(1, 2), hop(2, 3), (C + dag(C))(3), hop(3, 4), exp(-0.7im * N)(2), hop(1, 3)]
    qgates = [xy(1, 2), xy(2, 3), X(1), xy(3, 4), exp(-0.5im * Z)(2), X(2) * Y(3)]
    # with the total parity, which a mixed state need not have
    fobs = [N(1), N(2), N(4), C(1) * dag(C)(2), dag(C)(1) * C(3), C(2) * C(3), N(1) * N(3), C(1) + dag(C)(1),
            prod(F(j) for j in 1:4)]
    qobs = [Z(1), Z(3), X(1) * X(2), Y(2) * X(3), Z(1) * Z(4), X(1), Y(3), prod(Z(j) for j in 1:4), prod(X(j) for j in 1:4)]
    # occupied along X: Y and Z are the Majoranas, X the parity
    yz(a, b) = exp(-0.3im * (Z ⊗ Z + 0.6 * Y ⊗ Y))(a, b)
    xgates = [yz(1, 2), yz(2, 3), X(1), yz(3, 4), exp(-0.5im * X)(2), X(2) * Y(3)]
    applied(gates, st) = foldl((x, a) -> apply(a, x), gates; init = st)
    near(g, e, obs; atol = 1e-12) = TS.trace(g) ≈ real(TS.trace(e)) && all(abs(TS.expect_value(g, o) - TS.expect_value(e, o)) <= atol for o in obs)
    # the occupation of site 2 as a projector and as a measured operator
    cases = [(fs, ["FullyMixed", "Occ", "FullyMixed", "Emp"], fgates, fobs, Dict{Int, Symbol}(), N(2), N(2) => 1),
             (qs, ["FullyMixed", "Up", "Dn", "FullyMixed"], qgates, qobs, Dict{Int, Symbol}(), (Id(2) - Z(2)) / 2, Z(2) => -1),
             (qs, ["X+", "FullyMixed", "X-", "X+"], xgates, qobs, Dict(k => :X for k in 1:4), (Id(2) - X(2)) / 2, X(2) => -1)]
    for (sys, st, gates, obs, axes, proj, meas) in cases
        g = TS.GaussianState(State{Mixed}(sys, st); axes)
        e = TS.ExactState(State{Mixed}(sys, st))
        @test g isa TS.GaussianState{Mixed} && near(g, e, obs)
        g, e = applied(gates, g), applied(gates, e)
        @test near(g, e, obs; atol = 1e-10)
        for kept in ([1, 2], [2, 3], [1, 2, 3, 4])
            @test TS.reduced_density_matrix(g, kept) ≈ TS.reduced_density_matrix(e, kept) atol = 1e-10
            @test von_neumann_entropy(g, kept) ≈ von_neumann_entropy(e, kept) atol = 1e-10
        end
        # projections, measurements given and drawn
        @test near(apply(2 * proj, g), apply(2 * proj, e), obs; atol = 1e-10)
        r = collapse(g, meas)
        x = collapse(e, meas)
        @test r.probability ≈ x.probability atol = 1e-10
        @test near(r.state, x.state, obs; atol = 1e-10)
        d = collapse(g, first(meas); rng = Random.Xoshiro(3))
        @test d.probability ≈ collapse(e, first(meas) => d.result).probability atol = 1e-10
        if isempty(axes)
            @test length(sample(g; rng = Random.Xoshiro(4))) == 4
        else
            # the levels along Z are not occupations along X: projected, the mixed state would not commute with the parity
            @test_throws "not Gaussian on a mixed state" sample(g; rng = Random.Xoshiro(4))
        end
        # files
        mktempdir() do dir
            file = joinpath(dir, "m.h5")
            save_state(file, "m", g)
            h = load_state(file, "m")
            @test h isa TS.GaussianState{Mixed} && h.cov == g.cov && TS.trace(h) == TS.trace(g)
        end
    end
    # the modes of a fermion taken out and put back, C ρ C† and C† ρ C
    g = TS.GaussianState(State{Mixed}(fs, ["FullyMixed", "Occ", "FullyMixed", "Emp"]))
    e = TS.ExactState(State{Mixed}(fs, ["FullyMixed", "Occ", "FullyMixed", "Emp"]))
    g, e = applied(fgates[1:2], g), applied(fgates[1:2], e)
    for op in (C(2), dag(C)(3), 0.5 * dag(C)(1))
        @test near(apply(op, g), apply(op, e), fobs; atol = 1e-10)
    end
    @test near(TS.project_on(g, N(3)), TS.project_on(e, N(3); limits = Limits()), fobs; atol = 1e-10)
    @test TS.trace(2 * g) ≈ 2 * TS.trace(g) && TS.trace(g / 4) ≈ TS.trace(g) / 4
    # sites of several modes, and a fermionic gate between a fermion and a qubit
    sys = System([Fermion(), Electron(), Qubit(), Fermion()])
    st = ["FullyMixed", "FullyMixed", "Dn", "Occ"]
    g = TS.GaussianState(State{Mixed}(sys, st))
    e = TS.ExactState(State{Mixed}(sys, st))
    A = -0.3im * (dag(C) ⊗ C + dag(dag(C) ⊗ C))
    gates = [exp(A)(1, 4), exp(-0.2im * (dag(Cup) ⊗ C + dag(dag(Cup) ⊗ C)))(2, 4)]
    obs = [N(1), Nup(2), Ndn(2), Z(3), N(4), C(1) * dag(C)(4), Cup(2) * dag(C)(4)]
    @test near(applied(gates, g), applied(gates, e), obs; atol = 1e-10)
    # pure states of a definite parity, even and odd, made mixed
    for (sys, st, gates, obs) in ((fs, ["Occ", "Emp", "Occ", "Emp"], fgates, fobs), (qs, ["Up", "Dn", "Dn", "Dn"], qgates, qobs))
        p = applied(gates, TS.GaussianState(State{Pure}(sys, st)))
        g = TS.as_representation(Mixed, p)
        e = applied(gates, TS.ExactState(State{Mixed}(sys, st)))
        @test g isa TS.GaussianState{Mixed} && near(g, e, obs; atol = 1e-10)
        @test von_neumann_entropy(g, 1:4) ≈ 0 atol = 1e-10
        @test near(apply(gates[1], g), apply(gates[1], e), obs; atol = 1e-10)
    end
    # refused
    gq = TS.GaussianState(State{Mixed}(qs, "FullyMixed"))
    gf = TS.GaussianState(State{Mixed}(fs, "FullyMixed"))
    @test_throws "coherences" TS.GaussianState(State{Mixed}(qs, ["Up", "X+", "Up", "Up"]))
    @test_throws "odd in the Majoranas" apply(Ry(0.3)(1), gq)
    @test_throws "not this operator" apply(T(2), gq)
    @test_throws "not anti-Hermitian" apply(exp(0.5 * N)(2), gf)
    @test_throws "not Gaussian on a mixed state" TS.project_on(gq, (Id(1) + X(1)) / 2)
    @test_throws "several Kraus" apply(depolarizing_channel(0.2)(1), gq)
    @test_throws "one Gaussian component" TS.as_representation(Mixed, apply(T(2), TS.GaussianState(State{Pure}(qs, "X+"))))
    @test_throws "no definite parity" TS.as_representation(Mixed, TS.GaussianState(State{Pure}(qs, "X+")))
    @test_throws "does not evolve" evolve(-im * (X(1) * X(2)), 0.1, gq; nsteps = 1)
    @test_throws "operator space" entanglement_entropy(gq, 2)
    @test_throws "does not convert" State(gq)
    @test_throws "real number" 2im * gq
    @test_throws "does not thermalize" thermalize(Z(1), 1.0, gq; nsteps = 1)
    @test_throws "searches a pure state" ground_state(Z(1), gq)
end
