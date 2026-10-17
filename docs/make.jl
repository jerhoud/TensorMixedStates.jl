using Documenter
using TensorMixedStates, .Qubits, .Qudits, .Fermions, .Electrons, .Spins, .Bosons, .Tjs, .Qbosons,
    .Circuits

# context in which the doctests of the docstrings are run
DocMeta.setdocmeta!(
    TensorMixedStates,
    :DocTestSetup,
    quote
        using TensorMixedStates
        using .Qubits, .Qudits, .Fermions, .Electrons, .Spins, .Bosons, .Tjs, .Qbosons, .Circuits
    end;
    recursive = true,
)

makedocs(
    sitename = "TensorMixedStates",
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true", repolink = "https://github.com/jerhoud/TensorMixedStates.jl"),
    modules = [TensorMixedStates],
    pages = [
        "index.md",
        "quickref.md",
        "correspondence.md",
        "manual.md",
        "accuracy.md",
        "simulations.md",
        "conservation.md",
        "performance.md",
        "extending.md",
        "upgrading.md",
        "Reference" => [
            "sites.md",
            "operators.md",
            "states.md",
            "algorithms.md",
            "measurements.md",
            "highlevel.md",
            "others.md"
        ]
    ],
    checkdocs=:exports,
)

# llms.txt, the map of the documentation that language models read, and llms-full.txt, the
# guides in full, written into the build so that they are published at the root of the site
const site = "https://jerhoud.github.io/TensorMixedStates.jl/stable"
const guides = [
    ("quickref.md", "Quick reference", "the whole library on one page, with its conventions"),
    ("correspondence.md", "Coming from QuTiP or ITensors", "their functions next to those of TMS"),
    ("manual.md", "Manual", "sites, states, operators, algorithms and measurements, with examples"),
    ("accuracy.md", "Choosing and checking the parameters", "the truncation, the time step, the sweeps"),
    ("simulations.md", "High level interface", "simulations as phases, their files and checkpoints"),
    ("conservation.md", "Conserving a quantity", "weak and strong symmetries, and what they save"),
    ("performance.md", "Threads and performance", "threads, BLAS, long sums"),
    ("extending.md", "Extending TMS", "site types, phases, methods, destinations and states of one's own"),
    ("upgrading.md", "Upgrading from 1.x", "the names, keywords and defaults that changed in 2.0"),
]
let build = joinpath(@__DIR__, "build"), src = joinpath(@__DIR__, "src")
    open(joinpath(build, "llms.txt"), "w") do io
        println(io, "# TensorMixedStates\n")
        println(io, "> A Julia library simulating closed and open quantum systems, pure states and density ",
                "matrices, with matrix product states: ground and excited states, Hamiltonian and ",
                "Lindbladian evolution, steady and thermal states, gates and noisy channels, conserved ",
                "quantities, and simulations written as phases with checkpoints.\n")
        println(io, "## Guides\n")
        for (file, title, what) in guides
            println(io, "- [$title]($site/$(splitext(file)[1])/): $what")
        end
        println(io, "\n## Optional\n")
        println(io, "- [The guides in full]($site/llms-full.txt): the pages above in one file")
    end
    open(joinpath(build, "llms-full.txt"), "w") do io
        for (file, _, _) in guides
            println(io, read(joinpath(src, file), String))
        end
    end
end

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
# /dev/ follows the branch dev, where the next version is made, rather than the default branch
# of the repository, main, which only receives the releases; /stable/ is built from the tags
deploydocs(
    repo = "github.com/jerhoud/TensorMixedStates.jl.git",
    devbranch = "dev",
)
