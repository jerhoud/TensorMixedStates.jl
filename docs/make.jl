using Documenter
using TensorMixedStates, .Qubits, .Fermions, .Electrons, .Spins, .Bosons, .Tjs 

makedocs(
    sitename = "TensorMixedStates",
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true", repolink = "https://github.com/jerhoud/TensorMixedStates.jl"),
    modules = [TensorMixedStates],
    pages = [
        "index.md",
        "manual.md",
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
    checkdocs=:none,
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/jerhoud/TensorMixedStates.jl.git"
)
