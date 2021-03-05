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

The allocation functions also permit a matrix argument as a shortcut,
implicitly creating an `Additive`. For example, you can find a maximum Nash
welfare (MNW) allocation as follows:

    julia> alloc_mnw([1 1 2 3; 2 1 2 3]).alloc
    Allocation with 2 agents and 4 items:
      1 => {2, 4}
      2 => {1, 3}

Several allocation functions use mixed-integer linear programming via
[JuMP](https://jump.dev). Depending on the choice of MIP solver, solving even
moderately-sized instances may take a significant amount of time. Choosing a
different solver (from the default `Cbc.Optimizer`) may speed things up
considerably. For example, with the appropriate license, one could use use
[Gurobi](https://www.gurobi.com) as follows:

    Allocations.conf.MIP_SOLVER = Gurobi.Optimizer

It is also possible to supply the `Optimizer` (or other optimizer factories,
e.g., constructed using `optimizer_with_attributes`) as the `solver` keyword
argument to the relevant allocation functions.
"""
module Allocations

import JuMP
using JuMP: Model, optimizer_with_attributes, objective_value, @variable,
            @objective, @constraint, fix, optimize!, termination_status, MOI
using LightGraphs
using Cbc
using Random

include("exports.jl")
include("conf.jl")
include("types.jl")
include("util.jl")
include("checks.jl")
include("measures.jl")
include("mip.jl")
include("reductions.jl")
include("algorithms.jl")

end # module
