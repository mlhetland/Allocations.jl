# Social welfare functions


"""
    utility(V, A)

Compute the utilitarian welfare of the allocation `A`, given the profile `V`,
i.e., the sum of the individual agent utilities (i.e., the bundle values)
resulting from `A`.
"""
utility(V, A) = sum(value(V, i, A) for i in agents(V))


"""
    nash_welfare(V, A; nonzero=true)

Compute the Nash welfare of the allocation `A`, given the profile `V`, i.e., the
product of the individual agent utilities resulting from `A`. The `nonzero`
keyword indicates that agents with a utility of zero are left out. If no agents
with nonzero utility exist, the result is zero. To avoid overflow with large
utilities, the product is performed using floating-point arithmetic, even if the
utilities are integral.
"""
function nash_welfare(V, A; nonzero=true)
    x = Float64[value(V, i, A) for i in agents(V)]
    if nonzero
        x = x[x .> sqrt(eps(Float64))] # More than approximately zero
    end
    return isempty(x) ? 0.0 : prod(x)
end


# Approximation factors


"""
    prop_alpha(V, A)

Compute the fraction of proportionality guaranteed to every agent, that is, what
fraction each agent is guaranteed to get of `1/n` of their value for the grand
bundle `M`.
"""
function prop_alpha(V, A)
    na(V) * minimum(value(V, i, A) / value(V, i, items(V)) for i in agents(V))
end


# General approximation factor for EF-like properties (EF_ = EF, EF1, EFX, ...),
# where the value of the other agent's bundle is somehow modified, e.g., by
# removing some item. This modified value is provided by the function value_:
function _ef_alpha(V, A, value_)
    N = agents(V)
    frac = [value(V, i, A) / value_(V, i, bundle(A, j)) for i in N,
            j in N if i !== j]
    frac[isnan.(frac)] .= Inf
    return minimum(frac)
end


"""
    ef_alpha(V, A)

Find the approximation factor for *envy-freeness* in the allocation `A` with
the valuation profile `V`, i.e., how close to envy-freeness each agent is
guaranteed to get.

If every agent values all other agents' bundles as `0`, the `alpha` returned by
this function will be `Inf`. This is the case even if the agents also value
their own bundles as `0`.
"""
ef_alpha(V, A) = _ef_alpha(V, A, value)


"""
    ef1_alpha(V, A)

Find the approximation factor for *envy-freeness up to one item* in the
allocation `A` with the valuation profile `V`, i.e., how close to EF1 each agent
is guaranteed to get.

If every agent values all other agents' bundles as `0`, the `alpha` returned by
this function will be `Inf`. This is the case even if the agents also value
their own bundles as `0`.
"""
ef1_alpha(V, A) = _ef_alpha(V, A, value_1)


"""
    efx_alpha(V, A)

Find the approximation factor for *envy-freeness up to any item* in the
allocation `A` with the valuation profile `V`, i.e., how close to EFX each agent
is guaranteed to get.

If every agent values all other agents' bundles as `0`, the `alpha` returned by
this function will be `Inf`. This is the case even if the agents also value
their own bundles as `0`.
"""
efx_alpha(V, A) = _ef_alpha(V, A, value_x)
