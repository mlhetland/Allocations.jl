using Documenter
using Allocations

DocMeta.setdocmeta!(Allocations, :DocTestSetup, :(using Allocations); recursive=true)

makedocs(
    sitename = "Allocations",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        mathengine = Documenter.KaTeX(
            Dict(
                :macros => Dict(
                    "\\le"  => "\\leqslant",
                    "\\leq" => "\\leqslant",
                    "\\ge"  => "\\geqslant",
                    "\\geq" => "\\geqslant",
                    ["$c" => "{\\char $(codepoint(c))}" for c='A':'Z']...
                ),
            ),
        ),
    ),
    modules = [Allocations]
)

if haskey(ENV, "GITHUB_REPOSITORY")
    deploydocs(
        repo = "github.com/" * ENV["GITHUB_REPOSITORY"] * ".git"
    )
end
