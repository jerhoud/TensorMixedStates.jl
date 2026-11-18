# The representation Trajectory, a state held as one quantum trajectory: its types, the algebra it
# shares with Pure, its files, and what it refuses

const Trajectory = TensorMixedStates.Trajectory
const Ket = TensorMixedStates.Ket

@testset "Trajectory is a ket, as Pure, and holds its states as Pure does" begin
    @test Trajectory <: Ket && Pure <: Ket && !(Mixed <: Ket) && Ket <: TensorMixedStates.PM
    system = System(4, Fermion(conserve = strong(N)))
    @test all(siteind(Trajectory, system, i) == siteind(Pure, system, i) for i in 1:4)
    st = State{Trajectory}(system, [ isodd(i) ? "Occ" : "Emp" for i in 1:4 ])
    pure = State{Pure}(system, [ isodd(i) ? "Occ" : "Emp" for i in 1:4 ])
    @test st isa State{Trajectory} && TensorMixedStates.representation(st) === Trajectory
    overlap(a, b) = TensorMixedStates.ITensorMPS.inner(a.state, b.state)
    @test overlap(st, pure) ≈ 1
    mktempdir() do dir
        file = joinpath(dir, "t.h5")
        save_state(file, "t", st)
        back = load_state(file, "t")
        @test back isa State{Trajectory} && overlap(back, st) ≈ 1
    end
    @test_throws "ExactState holds Pure or Mixed states" ExactState{Trajectory}(System(2, Qubit()),
                                                                                "Up")
end
