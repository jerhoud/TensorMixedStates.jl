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
        tdvp(-im*(Y(2)+2Z(1)X(3)), 0.1, stp; maxdim = 3)
        approx_W(-im*(Y(2)+2Z(1)X(3)), 0.1, stp; order = 1, maxdim = 3)
    end
end