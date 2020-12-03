"""
    Allocations

Algorithms for fair allocation of indivisible items.

# Basic use

To specify an allocation problem instance, create a `Valuation`:

    julia> V = Valuation([1 2 3; 2 3 1])
    Additive{Array{Int64,2}} with 2 agents and 3 items:
     1  2  3
     2  3  1

`Valuation` is an abstract class, and `Valuation(X::Matrix)` is an alias for
`Additive(X)`. Once you have a valuation, you can use an allocation function
(ones `alloc_...`), e.g., for finding a maximin (MM) allocation:

    julia> res = alloc_mm(V);

These functions take a `Valuation` as input and return a named tuple with the
field `alloc` referring to an `Allocation`:

    julia> A = res.alloc
    Allocation with 2 agents and 3 items:
      1 => {1, 3}
      2 => {2}

The bundles of each agent is available through the `bundle` function:

    julia> bundle(A, 1)
    Set{Int64} with 2 elements:
      3
      1

!!! note

    Bundles should not be modified directly, as the `Allocation` also
    maintains an inverse mapping, from items to agents. Rather, use the
    `give!` and `deny!` functions.

Some allocation functions may produce other results as well, such as
properties of the allocation that are naturally computed as part of the
allocation process. For the maximin case, the objective value (the minimum
utility, which is being maximized) is available as `mm`:

    julia> res.mm
    3.0

The allocation functions also permit a matrix allocation as a shortcut,
implicitly creating an `Additive`. For example, you can find a maximum Nash
welfare (MNW) allocation as follows:

    julia> alloc_mnw([1 1 2 3; 2 1 2 3]).alloc
    Allocation with 2 agents and 4 items:
      1 => {2, 4}
      2 => {1, 3}

"""
module Allocations

import JuMP
using JuMP: Model, optimizer_with_attributes, objective_value, @variable,
            @objective, @constraint, fix, optimize!, termination_status, MOI
using Cbc

"""
    MIP_SOLVER

The (factory for) the JuMP optimizer to be used (by default) for mixed-integer
programming. Initially set to `Cbc.Optimizer`, with a `logLevel` of `0`. This
can be overridden either by setting `MIP_SOLVER` to another value (e.g., using
the JuMP function `optimizer_with_attributes`) or by passing it directly to
the appropriate allocation functions.
"""
MIP_SOLVER = nothing
# Not using the const Ref trick to avoid global variables, as we'll need to
# use type Any anyway.

function __init__()
    global MIP_SOLVER = optimizer_with_attributes(
        Cbc.Optimizer, "logLevel" => 0
    )
end

include("exports.jl")
include("types.jl")
include("util.jl")
include("checks.jl")
include("measures.jl")
include("mip.jl")

end # module
