# Checks for various fairness criteria.


"""
    check_complete(A)

Check that the allocation is complete, or effective, in the sense that each item
has been allocated to at least one agent.
"""
function check_complete(A)
    for g in items(A)
        length(owners(A, g)) >= 1 || return false
    end
    return true
end


"""
    check_partition(A)

Check that the allocation is a partition, i.e., that each item has been
allocated to exactly one agent.
"""
function check_partition(A)
    for g in items(A)
        length(owners(A, g)) == 1 || return false
    end
    return true
end


# General check for EF-like properties (EF_ = EF, EF1, EFX, ...), where the
# value of the other agent's bundle is somehow modified, e.g., by removing
# some item. This modified value is provided by the function value_:
function check_ef_(V, A, value_)

    N = agents(V)

    for i in N, j in N

        i !== j || continue

        if value(V, i, A) < value_(V, i, bundle(A, j))
            return false
        end

    end

    return true

end


"""
    check_ef(V, A)

Check whether the allocation `A` is *envy-free* for the profile `V`, i.e., if no
agent strictly prefers another agent's bundle.
"""
check_ef(V, A) = check_ef_(V, A, value)


"""
    check_ef1(V, A)

Check whether the allocation `A` is *envy-free up to one item* for the profile
`V`, i.e., if no agent strictly prefers another agent's bundle, given that an
appropriate (e.g., the most valuable) item is removed.
"""
check_ef1(V, A) = check_ef_(V, A, value_1)


"""
    check_efx(V, A)

Check whether the allocation `A` is *envy-free up to any item* for the profile
`V`, i.e., if no agent strictly prefers another agent's bundle, given that an
appropriate (e.g., the least valuable) item is removed.
"""
check_efx(V, A) = check_ef_(V, A, value_x)


"""
    check(V, A, C)

Check that the allocation `A` obeys the `Constraint` `C`, given the profile `V`.
"""
function check end


"""
    check(V, A, C::Nothing)

Trivial check that `A` satisfies a null-constraint. Always returns true.
"""
check(V, A, ::Nothing) = true


"""
    check(V, A, C::Counts)

Check whether the allocation `A` respects the cardinality constraints `C`.
"""
function check(V, A, C::Counts)

    for c in C

        counts = zeros(Int, na(A))

        for g in c, i in owners(A, g)
            counts[i] += 1
        end

        maximum(counts) <= c.threshold || return false

    end

    return true

end


"""
    check(V, A, C::Conflicts)

Check whether the allocation `A` respects the item conflicts `C`.
"""
function check(V, A, C::Conflicts)

    G = C.graph

    for e in edges(G)
        isempty(owners(A, src(e)) ∩ owners(A, dst(e))) || return false
    end

    return true

end


"""
    check(V, A, C::MatroidConstraints)

Check whether the allocation `A` respects the matroid constraints `C`.
"""
function check(_, A, C::MatroidConstraints)
    Ms = C.matroids

    @assert length(Ms) == na(A) "Invalid number of matroid constraints"

    for (i, B) in A
        if !is_indep(Ms[i], B)
            return false
        end
    end

    return true
end


"""
    check(V, A, C::MatroidConstraint)

Check whether the allocation `A` respects the matroid constraint `C`.
"""
function check(_, A, C::MatroidConstraint)
    M = C.matroid

    for (_, B) in A
        if !is_indep(M, B)
            return false
        end
    end

    return true
end
