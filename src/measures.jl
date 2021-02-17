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


"""
    prop_alpha(V, A)

Compute the fraction of proportionality guaranteed to every agent, that is, what
fraction each agent is guaranteed to get of `1/n` of their value for the grand
bundle `M`.
"""
function prop_alpha(V, A)
    na(V) * minimum(value(V, i, A) / value(V, i, items(V)) for i in agents(V))
end
