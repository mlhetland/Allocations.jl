"""
    nash_welfare(V, A; nonzero=true)

Compute the Nash welfare of the allocation `A`, given the valuation `V`, i.e.,
the product of the individual agent utilities resulting from `A`. The
`nonzero` keyword indicates that agents with a utility of zero are left out.
If no agents with nonzero utility exist, the result os zero.
"""
function nash_welfare(V, A; nonzero=true)
    x = [value(V, i, bundle(A, i)) for i in agents(V)]
    if nonzero
        x = x[x .> 0]
    end
    return isempty(x) ? 0.0 : prod(x)
end
