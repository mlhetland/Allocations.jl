mutable struct Conf
    MIP_SOLVER
    MIP_SUCCESS
end

"""
    conf

A struct with fields for global configuration of the `Allocations` module.


# Fields

    MIP_SOLVER :: Any

The (factory for) the JuMP optimizer to be used (by default) for mixed-integer
programming. Initially set to `Cbc.Optimizer`, with a `logLevel` of `0`. This
can be overridden either by setting `MIP_SOLVER` to another value (e.g., using
the JuMP function `optimizer_with_attributes`) or by passing the solver
directly to the appropriate allocation functions.

    MIP_SUCCESS :: Any

Container of acceptable MIP statuses. By default, has the value
`[MOI.OPTIMAL]`.
"""
const conf = Conf(nothing, nothing)

function __init__()
    conf.MIP_SOLVER = optimizer_with_attributes(
        Cbc.Optimizer, "logLevel" => 0
    )
    conf.MIP_SUCCESS = [MOI.OPTIMAL]
end


