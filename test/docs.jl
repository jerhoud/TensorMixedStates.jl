# The documentation shipped with the package against the code.
#
# Goes here: what a page states that the code can check, as the options of the quick reference,
# each field of Limits, Dmrg, Tdvp, ApproxW and Krylov written at its default.

@testset "the options of the quick reference" begin
    page = read(joinpath(pkgdir(TensorMixedStates), "docs", "src", "quickref.md"), String)
    lines = [ m.match for m in eachmatch(r"^(Limits|Dmrg|Tdvp|ApproxW|Krylov)\(.*$"m, page) ]
    calls = Dict(ex.args[1] => ex for ex in Meta.parse.(lines))
    @test length(lines) == 5
    @test sort!(collect(keys(calls))) == [:ApproxW, :Dmrg, :Krylov, :Limits, :Tdvp]
    required = Dict(:Dmrg => :nsweeps, :ApproxW => :order)
    for (name, call) in calls
        T = getfield(TensorMixedStates, name)
        given = call.args[2:end]
        @test [ kw.args[1] for kw in given ] == collect(fieldnames(T))
        written = eval(call)
        kept = filter(kw -> kw.args[1] == get(required, name, nothing), given)
        least = eval(Expr(:call, name, kept...))
        @test [ f => getfield(written, f) for f in fieldnames(T) ] ==
              [ f => getfield(least, f) for f in fieldnames(T) ]
        if haskey(required, name)
            @test_throws UndefKeywordError(required[name]) eval(Expr(:call, name))
        end
    end
end
