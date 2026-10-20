# The Jordan-Wigner strings: the operations where they matter, on the sites of every kind of
# fermion and with a site holding none between them, on pure and mixed states of no definite
# parity, which make every sign count, checked against dense matrices written from the matrices
# of one site alone. The descent and the compaction that lay the strings are tested in
# operators.jl and compact.jl, and each operation in its own group.

"""
the matrix on the whole of `system` of the matrix `m` of the sites `on`, in their order, written
in the Jordan-Wigner basis of those sites alone: each element of `m`, a product of elements of
one site, taken to the product of the same elements placed on the whole system, see `jw_placed`
"""
function jw_embedded(system, on, m)
    sites = [ system[k] for k in on ]
    ds = dim.(sites)
    parity(j, x) = fermion_parities(sites[j])[x]
    levels = CartesianIndices(Tuple(reverse(ds)))
    D = prod(dim(system[k]) for k in 1:length(system))
    out = zeros(ComplexF64, D, D)
    for X in axes(m, 1), Y in axes(m, 2)
        if m[X, Y] == 0
            continue
        end
        x, y = reverse(Tuple(levels[X])), reverse(Tuple(levels[Y]))
        p = [ xor(parity(j, x[j]), parity(j, y[j])) for j in eachindex(sites) ]
        # the placed product differs from the element by the sign each element of one site takes
        # going past the strings of those after it
        sign = prod((-1)^(sum(p[j+1:end]; init = 0) * parity(j, y[j])) for j in eachindex(sites))
        units = map(eachindex(sites)) do j
            u = zeros(ds[j], ds[j])
            u[x[j], y[j]] = 1
            return u
        end
        out += m[X, Y] * sign * prod(jw_placed(system, on[j], units[j]) for j in eachindex(sites))
    end
    return out
end

"""
the density matrix of the sites `kept`, in increasing order, of the dense density matrix `ρ` of
`system`: its elements are the values on `ρ` of those of the sites kept, see `jw_embedded`
"""
function jw_reduced(system, ρ, kept)
    d = prod(dim(system[k]) for k in kept)
    function unit(x, y)
        u = zeros(d, d)
        u[x, y] = 1
        return u
    end
    return [ LA.tr(ρ * jw_embedded(system, kept, unit(y, x))) for x in 1:d, y in 1:d ]
end

"""
the operators of one site the net draws on `site`, at position `k`: odd, even, odd again, and
one of no definite parity given by a random matrix of `rng`
"""
function jw_singles(site, k, rng)
    d = dim(site)
    generic = Operator{1}("M$k", randn(rng, ComplexF64, d, d))
    if site isa Fermion
        return [C, dag(C), N, C + 2dag(C), generic]
    elseif site isa Electron || site isa Tj
        return [Cup, dag(Cdn), Ntot, Cup + dag(Cdn), generic]
    elseif site isa Qubit
        return [X, Z, generic]
    end
    return [generic]
end

"""
a random product of `rng` on `system`, of one to three factors on sites in any order and apart:
operators of one site from `ops`, the operators of each site, tensor products of two, and at
most one operator of two sites given by a random matrix, on two sites of 9 states at most, as
one develops into many terms. The product, its adjoint and its dense matrix, see `jw_placed`
and `jw_embedded`
"""
function jw_product(system, ops, rng)
    n = length(system)
    placed(a, k) = jw_placed(system, k, matrix(a, system[k]))
    generic = false
    factors = map(1:rand(rng, 1:3)) do _
        i = rand(rng, 1:n)
        a = rand(rng, ops[i])
        r = rand(rng, generic ? (1:2) : (1:3))
        if r == 1
            return a(i), dag(a)(i), placed(a, i)
        end
        j = rand(rng, setdiff(1:n, i))
        if r == 2 || dim(system[i]) * dim(system[j]) > 9
            b = rand(rng, ops[j])
            return (a ⊗ b)(i, j), dag(a ⊗ b)(i, j), placed(a, i) * placed(b, j)
        end
        generic = true
        d = dim(system[i]) * dim(system[j])
        g = randn(rng, ComplexF64, d, d)
        G = Operator{2}("G$i$j", g)
        return G(i, j), dag(G)(i, j), jw_embedded(system, [i, j], g)
    end
    return prod(f[1] for f in factors), prod(f[2] for f in reverse(factors)),
           prod(f[3] for f in factors)
end

@testset "Jordan-Wigner strings on $name" for (name, sys) in [
        ("fermions", System(4, Fermion())), ("electrons", System(3, Electron())),
        ("t-J sites", System(3, Tj())),
        ("a qubit among fermions", System([Fermion(), Qubit(), Electron(), Fermion()]))]
    rng = Xoshiro(20261015)
    n = length(sys)
    ψ = random_state(Pure, sys, 4; rng)
    ρ = random_state(Mixed, sys, 4; rng)
    w = raw_vector(ψ)
    v = w / norm(w)
    R0 = ket_bra_matrix(ρ)
    R = R0 / LA.tr(R0)
    ops = [ jw_singles(sys[k], k, rng) for k in 1:n ]
    placed(a, k) = jw_placed(sys, k, matrix(a, sys[k]))
    fermionic = [ k for k in 1:n if any(isodd, fermion_parities(sys[k])) ]
    c(k) = sys[k] isa Fermion ? C : Cup

    @testset "Products of operators" begin
        for _ in 1:8
            o, o_dag, d = jw_product(sys, ops, rng)
            @test expect(ψ, o) ≈ v' * d * v atol = 1e-10 rtol = 1e-10
            @test expect(ρ, o) ≈ LA.tr(R * d) atol = 1e-10 rtol = 1e-10
            @test mpo_matrix(ψ, o) ≈ d atol = 1e-10 rtol = 1e-10
            @test raw_vector(apply(o, ψ)) ≈ d * w atol = 1e-10 rtol = 1e-10
            @test ket_bra_matrix(apply(o, ρ)) ≈ d * R0 * d' atol = 1e-10 rtol = 1e-10
            @test ket_bra_matrix(apply(Gate(o), ρ)) ≈ d * R0 * d' atol = 1e-10 rtol = 1e-10
            @test inner(ψ, o, ψ) ≈ w' * d * w atol = 1e-10 rtol = 1e-10
            h, dh = o + o_dag, d + d'
            var = real(v' * dh * dh * v - (v' * dh * v)^2)
            @test variance(ψ, h) ≈ var atol = 1e-10 rtol = 1e-10
            @test last(only(measure(ψ, [h]))) ≈ real(v' * dh * v) atol = 1e-10 rtol = 1e-10
            # the bound of the width of the spectrum the search of excited states reads
            λ = LA.eigvals(LA.Hermitian(dh))
            @test sum(TMS.term_width(sys, t) for t in TMS.lower(sys, h); init = 0.0) ≥
                  last(λ) - first(λ) - 1e-10
        end
        # functions of operators of two sites, applied as gates
        for _ in 1:4
            i, j = randperm(rng, n)[1:2]
            a, b = rand(rng, ops[i]), rand(rng, ops[j])
            g = exp(0.3 * (a ⊗ b) - 0.3 * dag(a ⊗ b))(i, j)
            d = placed(a, i) * placed(b, j)
            u = exp(0.3 * d - 0.3 * d')
            @test raw_vector(apply(g, ψ)) ≈ u * w atol = 1e-10 rtol = 1e-10
            @test ket_bra_matrix(apply(g, ρ)) ≈ u * R0 * u' atol = 1e-10 rtol = 1e-10
        end
    end

    @testset "Operators of one site, and pairs" begin
        defined(a, k) =
            try
                matrix(a, sys[k])
                true
            catch
                false
            end
        for k in 1:n, a in ops[k]
            for (st, value) in ((ψ, d -> v' * d * v), (ρ, d -> LA.tr(R * d)))
                e = expect1(st, a)
                @test [ ismissing(x) for x in e ] == [ !defined(a, j) for j in 1:n ]
                @test all(e[j] ≈ value(placed(a, j)) for j in 1:n if defined(a, j))
            end
            b = rand(rng, ops[rand(rng, 1:n)])
            e = expect2(ψ, (a, b))
            pair(i, j) = i == j ? jw_placed(sys, i, matrix(a, sys[i]) * matrix(b, sys[i])) :
                                  placed(a, i) * placed(b, j)
            @test all(e[i, j] ≈ v' * pair(i, j) * v for i in 1:n, j in 1:n
                      if defined(a, i) && defined(b, j))
        end
    end

    @testset "Parts of the state" begin
        S2(r) = -log(real(LA.tr(r * r)))
        SvN(r) = -sum(x * log(x) for x in LA.eigvals(LA.Hermitian((r + r') / 2)) if x > 0;
                      init = 0.0)
        full = v * v'
        for m in 1:2^n - 2
            kept = [ k for k in 1:n if isodd(m >> (k - 1)) ]
            rest = setdiff(1:n, kept)
            r, rm = jw_reduced(sys, full, kept), jw_reduced(sys, R, kept)
            @test reduced_density_matrix(ψ, kept) ≈ r atol = 1e-10
            @test reduced_density_matrix(ρ, kept) ≈ rm atol = 1e-10
            pt = ket_bra_matrix(partial_trace(ρ; keep = kept))
            @test pt / LA.tr(pt) ≈ rm atol = 1e-10
            @test renyi2(ψ, kept) ≈ S2(r) atol = 1e-10
            @test renyi2(ρ, kept) ≈ S2(rm) atol = 1e-10
            @test mutual_info_renyi2(ψ, kept) ≈ S2(r) + S2(jw_reduced(sys, full, rest)) atol = 1e-10
            @test mutual_info_renyi2(ρ, kept) ≈ S2(rm) + S2(jw_reduced(sys, R, rest)) - S2(R) atol =
                  1e-10
            @test vonneumann_entropy(ψ, kept) ≈ SvN(r) atol = 1e-10
            if kept == 1:last(kept)
                @test first(entanglement_entropy(ψ, last(kept))) ≈ SvN(r) atol = 1e-10
            end
        end
    end

    @testset "Measuring odd operators" begin
        # a Majorana operator on each fermionic site, a Hermitian product of two, and a Hermitian
        # operator of two sites given by its matrix, placed in reverse order, on 8 states at most
        i, j = first(fermionic), last(fermionic)
        γ(k) = c(k) + dag(c(k))
        η(k) = c(k) - dag(c(k))
        measured = Any[ (γ(k), (k,), placed(γ(k), k)) for k in fermionic ]
        push!(measured, (γ(i) ⊗ η(j), (i, j), placed(γ(i), i) * placed(η(j), j)))
        d = dim(sys[j]) * dim(sys[i])
        if d ≤ 8
            g = randn(rng, ComplexF64, d, d)
            g += g'
            push!(measured, (Operator{2}("H$j$i", g), (j, i), jw_embedded(sys, [j, i], g)))
        end
        for (a, on, dh) in measured
            h = a(on...)
            λs, U = LA.eigen(LA.Hermitian(dh))
            projector(λ) = (u = U[:, abs.(λs .- λ) .< 1e-8]; u * u')
            for (λ, p) in probabilities(ψ, h)
                @test p ≈ real(v' * projector(λ) * v) atol = 1e-10
            end
            for (λ, p) in probabilities(ρ, h)
                @test p ≈ real(LA.tr(R * projector(λ))) atol = 1e-10
            end
            λ = first(argmax(last, probabilities(ψ, h)))
            P = projector(λ)
            u = raw_vector(last(collapse(ψ, h => λ)))
            @test u / norm(u) ≈ P * v / norm(P * v) atol = 1e-10
            r = ket_bra_matrix(last(collapse(ρ, h => λ)))
            @test r / LA.tr(r) ≈ P * R * P / LA.tr(P * R * P) atol = 1e-10
            ps = [ projector(x) for x in unique(round.(λs; digits = 8)) ]
            @test ket_bra_matrix(apply(Dephase(a)(on...), ρ)) ≈ sum(p * R0 * p for p in ps) atol =
                  1e-10
        end
    end

    @testset "Channels of one site" begin
        for k in fermionic
            K = placed(c(k), k)
            @test ket_bra_matrix(apply(Dissipator(c(k))(k), ρ)) ≈
                  K * R0 * K' - (K' * K * R0 + R0 * K' * K) / 2 atol = 1e-10
            for s in (sys[k] isa Fermion ? ("Emp", "Occ") : ("Emp", "Up"))
                l = findfirst(!iszero, TMS.local_state(sys[k], s))
                kraus = map(1:dim(sys[k])) do x
                    e = zeros(dim(sys[k]), dim(sys[k]))
                    e[l, x] = 1
                    return jw_placed(sys, k, e)
                end
                after = apply(Reset(s)(k), ρ)
                @test ket_bra_matrix(after) ≈ sum(e * R0 * e' for e in kraus) atol = 1e-10
                # its docstring: an odd operator of another site takes the value of ⟨F⟩_σ F_k O
                for j in setdiff(fermionic, k)
                    taken = (-1)^fermion_parities(sys[k])[l] *
                            LA.tr(R * placed(F, k) * placed(c(j), j))
                    @test expect(after, c(j)(j)) ≈ taken atol = 1e-10
                end
            end
        end
    end

    @testset "Evolution" begin
        i, j = first(fermionic), last(fermionic)
        ci, cj = c(i), c(j)
        H = dag(ci)(i) * cj(j) + dag(cj)(j) * ci(i) +
            0.7 * (ci(i) * cj(j) + dag(cj)(j) * dag(ci)(i)) + 0.4 * (ci + dag(ci))(i)
        Ci, Cj = placed(ci, i), placed(cj, j)
        dH = Ci' * Cj + Cj' * Ci + 0.7 * (Ci * Cj + Cj' * Ci') + 0.4 * (Ci + Ci')
        @test raw_vector(evolve(-im * H, 0.3, ψ; nsteps = 3)) ≈ exp(-0.3im * dH) * w atol = 1e-9
        # a generator, by its MPO applied and by an evolution
        L = -im * H + Dissipator(cj)(j)
        loss = Cj' * Cj
        @test ket_bra_matrix(apply(make_mpo(ρ, L), ρ)) ≈
              -im * (dH * R0 - R0 * dH) + Cj * R0 * Cj' - (loss * R0 + R0 * loss) / 2 atol = 1e-10
        D = size(R0, 1)
        if D ≤ 32
            id = Matrix{Float64}(LA.I, D, D)
            sup = -im * (kron(id, dH) - kron(transpose(dH), id)) + kron(conj(Cj), Cj) -
                  (kron(id, loss) + kron(transpose(loss), id)) / 2
            evolved = ket_bra_matrix(evolve(L, 0.3, ρ; nsteps = 3))
            @test vec(evolved) ≈ exp(0.3 * sup) * vec(R0) atol = 1e-9
        end
    end
end

@testset "Slater determinants" begin
    # the product of the creation operators of the orbitals, written with their strings, on the
    # vacuum, up to a global phase
    rng = Xoshiro(20261015)
    orbitals(n, k) = Matrix(LA.qr(randn(rng, ComplexF64, n, k)).Q)[:, 1:k]
    vacuum(sys) = foldl(kron, [ (e = zeros(dim(sys[i])); e[1] = 1; e) for i in 1:length(sys) ])
    created(sys, c, Φ, sites) =
        prod(sum(Φ[p, k] * jw_placed(sys, s, matrix(dag(c), sys[s])) for (p, s) in enumerate(sites))
             for k in axes(Φ, 2))
    fe = System(4, Fermion())
    Φ = orbitals(4, 2)
    ref = created(fe, C, Φ, 1:4) * vacuum(fe)
    ψ = raw_vector(slater_state(fe, Φ))
    @test abs(ref' * ψ) ≈ norm(ref) * norm(ψ) atol = 1e-10
    el = System(3, Electron())
    up, dn = orbitals(3, 1), orbitals(3, 2)
    ref = created(el, Cup, up, 1:3) * created(el, Cdn, dn, 1:3) * vacuum(el)
    ψ = raw_vector(slater_state(el, up, dn))
    @test abs(ref' * ψ) ≈ norm(ref) * norm(ψ) atol = 1e-10
end
