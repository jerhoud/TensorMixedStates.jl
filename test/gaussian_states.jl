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
            @test !g.definite
        end
    end
    # states of definite parity, Fock states and eigenstates of the axes
    sys = System([Fermion(), Electron(), Qubit()])
    g = TensorMixedStates.GaussianState(State{Pure}(sys, ["Occ", "Up", "Dn"]))
    @test g.definite
    @test g.cov ≈ dense_covariance(gs_lifted(sys, reduce(kron, [[0, 1.0], [0, 1.0, 0, 0], [0, 1.0]]))) atol = 1e-12
    # refused: a state that is not a product, a mixed state
    sys = System([Qubit() for _ in 1:4])
    @test_throws ArgumentError TensorMixedStates.GaussianState(random_state(Pure, sys, 2))
    @test_throws ArgumentError TensorMixedStates.GaussianState(State{Mixed}(sys, "Up"))
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
    g = TensorMixedStates.GaussianState{Pure}(sys, Dict{Int, Symbol}(), layout, dense_covariance(Ψ), definite, 1.0)
    g, ψ
end

# the reduced density matrix of the sites kept of the vector ψ (the first site the slowest)
function gs_partial(ψ, dims, kept)
    n = length(dims)
    T = reshape(ψ, reverse(dims)...)
    rows = reverse([n + 1 - k for k in kept])
    A = reshape(permutedims(T, vcat(rows, setdiff(1:n, rows))), prod(dims[kept]), :)
    A * A'
end

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
            ρ = gs_partial(ψ, dims, kept)
            @test reduced_density_matrix(g, kept) ≈ ρ atol = 1e-12
            @test von_neumann_entropy(g, kept) ≈ gs_entropy(ρ) atol = 1e-10
        end
        for cut in 0:n
            @test entanglement_entropy(g, cut) ≈ gs_entropy(gs_partial(ψ, dims, 1:cut)) atol = 1e-10
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
    dims = [dim(s) for s in sys.sites]
    for cut in 1:n-1
        ok &= isapprox(entanglement_entropy(g, cut), gs_entropy(gs_partial(ψ / LA.norm(ψ), dims, 1:cut)); atol = 1e-9)
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
    # refused: an evolver depending on time, one that is not unitary
    g, _ = gs_product(sys1)
    @test_throws ArgumentError TensorMixedStates.prepared_evolver(g, [-im * N(1), -im * N(2) => t -> t])
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
    @test_throws ArgumentError apply(H(2), g)
    @test_throws ArgumentError apply(CNOT(2, 3), g)
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
            @test h.cov == g.cov && h.axes == g.axes && h.definite == g.definite && h.trace == g.trace
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
function gs_programs_agree(phases)
    ok = true
    in_tmp() do
        g = runTMS(SimData(; name = "g", backend = TensorMixedStates.GaussianState, phases))
        e = runTMS(SimData(; name = "e", backend = TensorMixedStates.ExactState, phases))
        ok &= g.state isa TensorMixedStates.GaussianState
        ok &= gs_same_outputs("g", "e")
        v = TensorMixedStates.ExactState(State(g.state)).state
        ok &= isapprox(abs(dot(v, e.state.state)), LA.norm(e.state.state); atol = 1e-10)
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
