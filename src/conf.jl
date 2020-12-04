mutable struct Conf
    MIP_SOLVER
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
"""
const conf = Conf(nothing)

function __init__()
    conf.MIP_SOLVER = optimizer_with_attributes(
        Cbc.Optimizer, "logLevel" => 0
    )
end


