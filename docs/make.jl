using Documenter
using TensorMixedStates, .Qubits, .Qudits, .Fermions, .Electrons, .Spins, .Bosons, .Tjs, .Qbosons

# context in which the doctests of the docstrings are run
DocMeta.setdocmeta!(
    TensorMixedStates,
    :DocTestSetup,
    quote
        using TensorMixedStates
        using .Qubits, .Qudits, .Fermions, .Electrons, .Spins, .Bosons, .Tjs, .Qbosons
    end;
    recursive = true,
)

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
    checkdocs=:exports,
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/jerhoud/TensorMixedStates.jl.git"
)
