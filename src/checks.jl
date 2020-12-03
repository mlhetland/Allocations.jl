# Checks for various fairness criteria.


# General check for EF-like properties (EF_ = EF, EF1, EFX, ...), where the
# value of the other agent's bundle is somehow modified, e.g., by removing
# some item. This modified value is provided by the function value_:
function check_ef_(V, A, value_)

    N = agents(V)

    for i in N, j in N

        i !== j || continue

        if value(V, i, bundle(A, i)) < value_(V, i, bundle(A, j))
            return false
        end

    end

    return true

end


"""
    check_ef(V, A)

Check whether the allocation `A` is *envy-free* for the valuation `V`, i.e.,
if no agent strictly prefers another agent's bundle.
"""
check_ef(V, A) = check_ef_(V, A, value)


"""
    check_ef1(V, A)

Check whether the allocation `A` is *envy-free up to one item* for the
valuation `V`, i.e., if no agent strictly prefers another agent's bundle,
given that an appropriate (e.g., the most valuable) item is removed.
"""
check_ef1(V, A) = check_ef_(V, A, value_1)


"""
    check_efx(V, A)

Check whether the allocation `A` is *envy-free up to any item* for the
valuation `V`, i.e., if no agent strictly prefers another agent's bundle,
given that an appropriate (e.g., the least valuable) item is removed.
"""
check_efx(V, A) = check_ef_(V, A, value_x)
