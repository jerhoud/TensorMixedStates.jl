# The precompilation workload, a small simulation run when the package is precompiled; set
# TMS_SKIP_PRECOMPILE_WORKLOAD=true to skip it during development, never for a release.

using PrecompileTools: @compile_workload

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
        # dmrg and apply need calls of their own, those of tdvp not compiling them; the mixed
        # `apply` is a path of its own, the operator being wrapped in a `Gate`
        dmrg(Z(1)Z(2) + Z(2)Z(3), stp; nsweeps = 2, limits = Limits(maxdim = 3))
        apply(X(1) * H(2), stp; limits = Limits(maxdim = 3))
        apply(X(1), stm; limits = Limits(maxdim = 3))
        # a product state has real tensors, which take paths of their own through the
        # solvers, and runTMS has its own as well
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