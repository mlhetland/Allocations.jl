module Allocations

import JuMP
using JuMP: Model, optimizer_with_attributes, objective_value, @variable,
            @objective, @constraint, fix, optimize!, termination_status, MOI
using Graphs
using HiGHS
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
include("precompile.jl")

end # module
