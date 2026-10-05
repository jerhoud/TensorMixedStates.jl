# The precompilation workload, a small simulation run when the package is precompiled; set
# TMS_SKIP_PRECOMPILE_WORKLOAD=true to skip it during development.

using PrecompileTools: @compile_workload

# skip this workload during development (set TMS_SKIP_PRECOMPILE_WORKLOAD=true)
# to speed up rebuilds; keep it enabled for releases so end users get the
# usual precompiled-code benefit on first use
if get(ENV, "TMS_SKIP_PRECOMPILE_WORKLOAD", "false") != "true"
    @compile_workload begin
        s = System(3, Qubit())
        stp = RandomState{Pure}(s, 4)
        stpm = mix(stp)
        stm = State{Mixed}(s, ["FullyMixed", "Y+", "Z+"])
        sstp = string(stp)
        sstm = string(stm)
        m = Measure(X, Y(1), Z(2)Y(1) + Z(3)Y(2), (X, Y), Y, (Z, Z))
        sm = string(m)
        mp = measure(stp, m)
        mm = measure(stm, m)
        tdvp(-im*(Y(2)+2Z(1)X(3)), 0.1, stp; limits = Limits(maxdim = 3))
        approx_W(-im*(Y(2)+2Z(1)X(3)), 0.1, stp; order = 1, limits = Limits(maxdim = 3))
        # measured: with the workload as it stood, the first `tdvp` went from 30 s to
        # 3.6 s while the first `dmrg` stayed at 7.8 and the first `apply` at 5.0, for
        # want of being here. The mixed `apply` is a path of its own, the operator being
        # wrapped in a `Gate` on the way in
        dmrg(Z(1)Z(2) + Z(2)Z(3), stp; nsweeps = 2, limits = Limits(maxdim = 3))
        apply(X(1) * H(2), stp; limits = Limits(maxdim = 3))
        apply(X(1), stm; limits = Limits(maxdim = 3))
        # a product state has real tensors, which take paths of their own through the
        # solvers, and runTMS has its own as well. Measured: the first dmrg on a product state
        # went from 22 s to 5 s, the first tdvp of a mixed one from 14 s to 4 s and the first
        # runTMS from 7 s to 0.2 s
        stpr = State{Pure}(s, "Up")
        stmr = State{Mixed}(s, "Up")
        measure(stpr, m)
        measure(stmr, m)
        dmrg(Z(1)Z(2) + X(2), stpr; nsweeps = 2, limits = Limits(maxdim = 3))
        tdvp(-im*(Y(2)+2Z(1)X(3)), 0.1, stpr; limits = Limits(maxdim = 3))
        tdvp(-im*X(1) + Dissipator(Sm)(2), 0.1, stmr; limits = Limits(maxdim = 3))
        mktempdir() do dir
            cd(dir) do
                runTMS(SimData(name = "workload", phases = [CreateState{Mixed}(3, Qubit(), "Up"),
                    Evolve(duration = 0.2, time_step = 0.1, algo = Tdvp(), limits = Limits(maxdim = 3),
                           evolver = -im * X(1) + Dissipator(Sm)(2),
                           measurements = "data" => [Z, Purity])]); output = devnull)
            end
        end
    end
end