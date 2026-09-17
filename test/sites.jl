# Per site type behaviour.
#
# Goes here: for each site type (Qubit, Fermion, Boson, Spin, Electron, Tj, Qboson),
# that its named local states and local operators have the expected values, and that
# the mixed representation of a state agrees with the pure one. One testset per site
# type, so a new site type gets a new testset here.

@testset "Qubit measuring" begin
    @test_pm test_phases(CreateState{type}(1, Qubit(), "Z+"; 
        final_measures = check([X(1), Y(1), Z(1)], [0, 0, 1])))
    @test_pm test_phases(CreateState{type}(6, Qubit(), ["X+", "Y+", "Z+", "X-", "Y-", "Z-"];
        final_measures = check([X, Y, Z], [[1, 0, 0, -1, 0, 0], [0, 1, 0, 0, -1, 0], [0, 0, 1, 0, 0, -1]])))
    @test_pm test_phases(CreateState{type}(3, Qubit(), ["X+", "Y-", "Z-"];
        final_measures = check([(X, Y), (Z, Z)], [[0 -1 0; 0 0 0; 0 0 -im], [1 0 0; 0 1 0; 0 0 1]])))
    @test_pm test_phases(CreateState{type}(4, Qubit(), ["Z+", "X-", "Z-", "Y-"];
        final_measures = check([Z(1)Y(4), X(2)Z(3), Z(1)Z(3), X(2)Y(4)], [-1, 1, -1, 1])))
    @test norm(matrix(Sx^2+Sy^2+Sz^2-S2, Qubit()))≈0 atol=1e-12
end

@testset "Fermion measuring" begin
    @test_pm test_phases(CreateState{type}(1, Fermion(), "1";
        final_measures = check(N(1), 1)))
    @test_pm test_phases(CreateState{type}(2, Fermion(), ["0", "1"];
        final_measures = check(N, [0, 1])))
end

@testset "Boson measuring" begin
    @test_pm test_phases(CreateState{type}(4, Boson(4), ["0", "1", "2", "3"];
        final_measures = check(N, [0, 1, 2, 3])))
end

@testset "Spin measuring" begin
    @test_pm test_phases(CreateState{type}(4, Spin(3/2), ["-3/2", "-1/2", "1/2", "3/2"];
        final_measures = check([Sx, Sy, Sz], [[0, 0, 0, 0], [0, 0, 0, 0], [-3/2, -1/2, 1/2, 3/2]])))
    @test norm(matrix(Sx^2+Sy^2+Sz^2-S2,Spin(5/2)))≈0 atol=1e-12
    @test norm(matrix(Sx^2+Sy^2+Sz^2-S2,Spin(4)))≈0 atol=1e-12
    @test_pm test_phases(CreateState{type}(9, Spin(1),
            ["X-1", "X0", "X1", "Y-1", "Y0", "Y1", "Z-1", "Z0", "Z1"];
        final_measures = check([Sx, Sy, Sz],
            [[-1, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, -1, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, -1, 0, 1]])))
    @test_pm test_phases(CreateState{type}(12, Spin(3/2),
            ["X-3/2", "X-1/2", "X1/2", "X3/2", "Y-3/2", "Y-1/2", "Y1/2", "Y3/2", "Z-3/2", "Z-1/2", "Z1/2", "Z3/2"];
        final_measures = check([Sx, Sy, Sz],
            [[-3/2, -1/2, 1/2, 3/2, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, -3/2, -1/2, 1/2, 3/2, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, -3/2, -1/2, 1/2, 3/2]])))
end

@testset "Electron measuring" begin
    @test_pm test_phases(CreateState{type}(4, Electron(), ["Emp", "Up", "Dn", "UpDn"];
        final_measures = check([Nup, Ndn, Nupdn, Ntot], [[0, 1, 0, 1], [0, 0, 1, 1], [0, 0, 0, 1], [0, 1, 1, 2]])))        
end

@testset "Tj measuring" begin
    @test_pm test_phases(CreateState{type}(3, Tj(), ["Emp", "Up", "Dn"];
        final_measures = check([Nup, Ndn, Ntot], [[0, 1, 0], [0, 0, 1], [0, 1, 1]])))        
end

@testset "Qboson measuring" begin
    @test_pm test_phases(CreateState{type}(4, Qboson(0.1, 4), ["0", "1", "2", "3"];
        final_measures = check(N, [0, 1, 2, 3])))
end

@testset "State mixing" begin
    @test_ok check_mix(
        [1, 5, 15],
        [Qubit(), Fermion(), Boson(3), Spin(3/2), Electron(), Tj()],
        [[X, Y, Z], [N], [N], [Sx, Sy, Sz], [Nup, Ndn, Nupdn, Ntot, Sx, Sy, Sz], [Nup, Ndn, Ntot, Sx, Sy, Sz]])
end

@testset "Mix ordering" begin
    sys = System(2, Qubit())
    ux2 = Id ⊗ X
    psi = State{Pure}(sys, "0")
    rho_from_pure = mix(apply(ux2(1, 2), psi))
    rho_direct = apply(ux2(1, 2), mix(psi))

    @test expect(rho_from_pure, Z(1)) ≈ 1
    @test expect(rho_direct, Z(1)) ≈ 1
    @test expect(rho_from_pure, Z(2)) ≈ -1
    @test expect(rho_direct, Z(2)) ≈ -1
end
