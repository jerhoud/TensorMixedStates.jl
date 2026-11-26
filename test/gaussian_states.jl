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
    g = TensorMixedStates.GaussianState{Pure}(sys, Dict{Int, Symbol}(), layout, dense_covariance(Ψ), nothing, definite, 1.0, 64)
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

@testset "checks from the reviews of lot 3" begin
    # a monomial gate multiplies the norm by its coefficient
    sys = System([Spin(1/2), Spin(1/2)])
    g, _ = gs_product(sys)
    @test TensorMixedStates.trace(apply(Sz(1), g)) ≈ TensorMixedStates.trace(g) / 4
    @test TensorMixedStates.trace(apply(Sx(1) * Sz(2), g)) ≈ TensorMixedStates.trace(g) / 16
    # eigenstates of the axes have a definite parity, whose entropies take the Gaussian formula
    g = TensorMixedStates.GaussianState(State{Pure}(System(30, Qubit()), "X+"); axes = :X)
    @test g.definite && abs(entanglement_entropy(g, 15)) < 1e-12
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
