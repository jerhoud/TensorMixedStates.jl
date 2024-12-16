using PrecompileTools: @compile_workload

@compile_workload begin
    s = System(3, "Qubit")
    stp = State(Pure, s, "+")
    stm = State(Mixed, s, "Y+")
    sstp = string(stp)
    sstm = string(stm)
    m = Measure(X, Y(1), Z(2)Y(1) + Z(3)Y(2), (X, Y))
    mp = measure(stp, m)
    mm = measure(stm, m)
end