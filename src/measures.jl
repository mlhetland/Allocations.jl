"""
    nash_welfare(V, A; nonzero=true)

Compute the Nash welfare of the allocation `A`, given the valuation `V`, i.e.,
the product of the individual agent utilities resulting from `A`. The
`nonzero` keyword indicates that agents with a utility of zero are left out.
If no agents with nonzero utility exist, the result is zero. To avoid overflow
with large utilities, the product is performed using floating-point arithmetic,
even if the utilities are integral.
"""
function nash_welfare(V, A; nonzero=true)
    x = Float64[value(V, i, A) for i in agents(V)]
    if nonzero
        x = x[x .> sqrt(eps(Float64))] # More than approximately zero
    end
    return isempty(x) ? 0.0 : prod(x)
end
