using Documenter
using Allocations

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

repo = get(ENV, "GITHUB_REPOSITORY", nothing) ≠ nothing &&
deploydocs(
    repo = "github.com/$repo.git"
)

; # Suppress `false` in REPL
