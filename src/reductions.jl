"""
    reduce(V::Additive, C::Counts{OrderedCategory}, i, B)

Reduce the instance given by the pair (V, C) to a new instance by giving the
supplied agent, `i`, the supplied bundle, `B`. Returns a reduction, where the
transformation, in addition to converting the allocation to one for the original
instance, allocates `B` to `i`.
"""
function reduce(V::Additive, C::Counts{OrderedCategory}, i, B)
    N, M = agents(V), items(V)
    n′, m′ = na(V) - 1, ni(V) - length(B)

    Vs = zeros(n′, m′)
    λg = Int[]

    for g in M
        g ∈ B && continue

        push!(λg, g)
        g′ = length(λg)

        Vs[1:n′,g′] = [value(V, j, g) for j in N if i != j]
    end

    # Create new ordered categories
    Cs = OrderedCategory[]
    index = 1
    for c in C
        newlength = length(c ∩ λg)
        push!(Cs, OrderedCategory(index, newlength, c.threshold))
        index += newlength
    end

    λi = [j for j in agents(V) if i != j]
    V′, C′ = Additive(Vs), Counts(Cs)

    return Reduction(V′, C′, λi, λg, A -> revert(λg, i, B, A))
end


"""
    reduce(V::Submodular, i, B)

Reduce the instance given by `V` to a new instance by giving the specified
bundle, `B`, to agent `i`. Returns a reduction, where the transformation, in
addition to converting the allocation to one for the original instance,
allocates `B` to `i`.
"""
function reduce(V::Submodular, i, B)

    N, M = agents(V), items(V)

    # Translation table
    λg = [g for g in M if g ∉ B]
    λi = [j for j in N if j != i]

    Vf′ = [B -> value(V, j, Set(λg[g] for g in B)) for j in λi]
    V′ = Submodular(Vf′, length(λg))

    return Reduction(V′, λi, λg, (A) -> revert(λg, i, B, A))

end


"""
    revert(translate::Vector{Int}, i, B, A::Allocation)

Convert an allocation for a reduced instance to one for the original instance,
including giving the removed bundle, `B`, to the removed agent, `i`.
"""
function revert(λg::Vector{Int}, i, B, A::Allocation)
    A′ = Allocation(na(A) + 1, ni(A) + length(B))

    for j in 1:na(A)
        j′ = j + (j >= i)
        give!(A′, j′, [λg[g] for g in bundle(A, j)])
    end

    give!(A′, i, B)

    return A′
end


"""
    reduce(V::Additive, C::Vector{OrderedCategory}, α)

Reduce an ordered instance by normalizing the values and giving any agent that
value an individual item greater than or equal to α the item and any low value
items required to reduce to a valid instance. This reduction is performed
recursively until no more such items exist. The reduction does not decrease the
MMS guarantee of any remaining agents and all agents that are allocated a
bundle in the reduction is guaranteed to value their bundle at least α of their
MMS guarantee.
"""
function reduce(V::Additive, C::Counts{OrderedCategory}, α::Real)
    N, n, M = agents(V), na(V), items(V)

    if n == 1 return reduce(V, C, 1, M) end

    V = normalize(V)

    # Find any agent who values a single item at least α and give that item to
    # the agent along with the least valuable items in the remaining bundles
    # that must be given to the agent to guarantee a valid instance.
    for i in N, g in M
        if value(V, i, g) >= α
            bundle = union(
                Set{Int}(g),
                [c[end - required(c, n) + (g in c) + 1:end] for c in C]...)

            R = reduce(V, C, i, bundle)

            V, C = profile(R), constraint(R)
            # Recursive application, as the removed items may cause items worth
            # less than α to be worth α or more after a new normalization.
            R′ = reduce(V, C, α)

            # Combine the reductions
            return chain(R, R′)
        end
    end

    # Return an empty reduction
    return Reduction(V, C)
end


"""
    reduce(V::Profile, α::Real)

Produce a reduced instance by giving an item to any agent that values it at `α`
or more. This reduction is performed repeatedly, until no such item exists.
"""
function reduce(V::Profile, α::Real)
    N, M = agents(V), items(V)

    for i in N, g in M
        if value(V, i, g) ≥ α
            R = reduce(V, i, Set(g))
            V = profile(R)

            # Recursive application on the new instance
            R′ = reduce(V, α)

            # Combine the reductions
            return chain(R, R′)
        end
    end

    return Reduction(V)
end


"""
    order(V::Additive, C::Counts)

Create an ordered instance for the given weights and categories. The items are
reorded such that each category has a continous range of indices for its items.
Returns a reduction, with a transformation that converts an allocation to one in
the original instance where each agent gets at least the same value as in the
ordered instance.
"""
function order(V::Additive, C::Counts)
    N = agents(V)
    Vo = zeros(na(V), ni(V))
    Co = OrderedCategory[]

    m = 1
    for c in C
        for i in N
            Vo[i,m:m+length(c)-1] = sort([value(V, i, g) for g in c], rev=true)
        end

        push!(Co, OrderedCategory(m, length(c), threshold(c)))
        m += length(c)
    end

    V′, C′ = Additive(Vo), Counts(Co)

    return Reduction(V′, C′, agents(V), items(V), A -> revert(V, C, C′, A))
end


"""
    revert(V::Additive, C::Counts, C′::Counts, A)

Convert an allocation for the ordered instance (`C′`) to one for the original
instance `(V, C)`.
"""
function revert(V::Additive, C::Counts, C′::Counts, A)
    A′ = Allocation(na(A), ni(A))

    for (unordered_category, ordered_category) in zip(C, C′)
        items = [g for g in unordered_category]
        for g in ordered_category
            i = owner(A, g)

            # Select the good in the category with the highest value
            g′ = maximum((g) -> (value(V, i, g), g), items)[2]

            give!(A′, i, g′)
            items = setdiff(items, g′)
        end
    end

    return A′
end
