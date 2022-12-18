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

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
