module Allocations

import Base: reduce
# Temporary (cf. https://github.com/mlhetland/Allocations.jl-private/issues/66)

import JuMP
using JuMP: Model, optimizer_with_attributes, objective_value, objective_sense,
            @variable, @objective, @constraint, @build_constraint, delete, fix,
            optimize!, termination_status, set_attribute, callback_node_status,
            callback_value, MOI
using Distributions
using Graphs
using HiGHS
using Random
using Random: default_rng

include("exports.jl")
include("conf.jl")
include("smallbitset.jl")
include("matroids.jl")
include("types.jl")
include("util.jl")
include("checks.jl")
include("matroid_algorithms.jl")
include("matroid_checks.jl")
include("matroid_util.jl")
include("measures.jl")
include("mip.jl")
include("reductions.jl")
include("algorithms.jl")
include("data.jl")
include("deprecated.jl")
include("precompile.jl")

end # module
