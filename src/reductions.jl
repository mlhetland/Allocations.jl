import Base: reduce
# Temporary (cf. https://github.com/mlhetland/Allocations.jl-private/issues/66)


"""
    reduceutil(V::Profile, assignments::Pair...)

Utility function that given valuations and a collection of assignments of
bundles to agents (`i => B`), creates a reduced instance, translation tables
from the reduced instance and a function to convert an allocation in the reduced
instance to one in the original instance -- including the given assignements.
The function returns a `Reduction` object without any constraints.
"""
function reduceutil(V::Profile, assignments::Pair...)
    # Create translation tables
    λi = collect(agents(V))
    λg = collect(items(V))

    # Replace all agents and items assigned with 0
    for (i, B) in assignments
        λi[i] = 0
        for g in B
            λg[g] = 0
        end
    end

    # Remove assigned agents and items
    filter!(!iszero, λi)
    filter!(!iszero, λg)

    # Gather relevant valuations
    V′ = reducevaluation(V, λi, λg)

    return Reduction(V′, λi, λg, (A) -> revert(λi, λg, assignments, A))
end


"""
    reducevaluation(V::Additive, λi, λg)

Utility function that given additive valuations prior to a reduction and
translation tables for the reduction, returns new additive valuations for the
reduced instance. The new valuations are as prior to the reduction, except for
missing items/agents and changes in item/agent numbers.
"""
reducevaluation(V::Additive, λi, λg) = Additive([value(V, i, g) for i=λi, g=λg])


"""
    reducevaluation(V::Submodular, λi, λg)

Utility function that given submodular valuations prior to a reduction and
translation tables for the reduction, returns new submodular valuations for the
reduced instance. The new valuations are as prior to the reduction, except for
missing items/agents and changes in item/agent numbers. That is, the new
valuation functions work by translating the item numbers to what they would be
prior to the reduction and calling the valuation function of the agent prior to
the reduction.
"""
function reducevaluation(V::Submodular, λi, λg)
    Vf′ = [B -> value(V, i, [λg[g] for g in B]) for i in λi]
    return Submodular(Vf′, length(λg))
end


"""
    reduce(V::Valuation, assignment::Pair...)

Reduce the instance given to a new instance where the involved agents and
bundles in the assignments are removed. Returns new valuations and a function
that turns an allocation in the reduced instance into one for the original
instance, including giving the supplied agent the supplied bundle.
"""
reduce(V::Profile, assignments::Pair...) = reduceutil(V, assignments...)


"""
    reduce(V::Additive, C::Counts{OrderedCategory}, i, B)

Reduce the instance given by the pair (V, C) to a new instance by giving the
supplied agent, `i`, the supplied bundle, `B`. Returns a reduction, where the
transformation, in addition to converting the allocation to one for the original
instance, allocates `B` to `i`.
"""
function reduce(V::Additive, C::Counts{OrderedCategory}, i, B)
    reduction = reduceutil(V, i => B)

    # Create new ordered categories
    C′ = OrderedCategory[]
    index = 1
    for c in C
        newlength = length(c ∩ reduction.λg)
        push!(C′, OrderedCategory(index, newlength, c.threshold))
        index += newlength
    end

    return Reduction(reduction, Counts(C′))
end


"""
    reduce(V::Submodular, i, B)

Reduce the instance given by `V` to a new instance by giving the specified
bundle, `B`, to agent `i`. Returns a reduction, where the transformation, in
addition to converting the allocation to one for the original instance,
allocates `B` to `i`.
"""
reduce(V::Submodular, i, B) = reduceutil(V, i => B)


"""
    revert(λi, λg, assignments, A)

Convert an allocation for a reduced instance to one for the original instance,
including giving the removed bundles to the removed agents.
"""
function revert(λi, λg, assignments, A)
    A′ = Allocation(
        na(A) + length(assignments),
        ni(A) + sum(el -> length(el.second), assignments)
    )

    for i in 1:na(A)
        give!(A′, λi[i], [λg[g] for g in bundle(A, i)])
    end

    for (i, B) in assignments give!(A′, i, B) end

    return A′
end


"""
    reduce(V::Additive, F::Function...)

Reduce an instance V by repeatedly applying the functions f ∈ F to find bundles
to be allocated. The functions in F are expected to return either a pair, `(i,
B)`, consisting of an agent `i` and the bundle `B` to be assigned to agent `i`,
or the value `nothing` if the function couldn't find a valid bundle-agent-pair.
The functions are called in prioritized order and the instance is reduced and
normalized between each invocation. The functions are invoked with the valuation
matrix.
"""
function reduce(V::Additive, F::Function...)

    if na(V) == 0 return Reduction(V) end
    if na(V) == 1 return reduce(V, 1 => items(V)) end

    V = normalize(V)

    for f in F
        res = f(V)
        if res !== nothing
            i, B = res
            R₁ = reduce(V, i => B)
            R₂ = reduce(R₁.profile, F...)
            return chain(R₁, R₂)
        end
    end

    # No reduction to apply
    return Reduction(V)
end


"""
    reduce(V::Additive, α::Real; greedy::Bool=true)

Reduce an ordered instance by normalizing the values and giving any agent that
value an individual item greater than or equal to α the item.  This reduction is
performed recursively until no more such items exist. The reduction does not
decrease the MMS guarantee of any remaining agents and all agents that are
allocated a bundle in the reduction is guaranteed to value their bundle at least
α of their MMS guarantee. The agent-item pairs are either selected greedily or
by finding a maximum matching between agents and such items.
"""
function reduce(V::Additive, α::Real; greedy::Bool=true)
    if greedy
        function find_combo(V)
            i = findfirst(i -> value(V, i, 1) ≥ α, agents(V))
            return i === nothing ? nothing : i => [1]
        end

        return reduce(V, find_combo)
    end

    if na(V) == 0 return Reduction(V) end
    if na(V) == 1 return reduce(V, 1 => items(V)) end

    V = normalize(V)
    X = BitArray(value(V, i, g) ≥ α for i = agents(V), g = items(V))
    x = bipartite_matching(X)

    allocations = [i => [g] for (i, g) in x]

    if isempty(allocations) return Reduction(V) end

    # If no agent remains after allocating the matching, give one of the agents
    # all the remaining goods.
    if length(allocations) == na(V)
        append!(allocations[end].second,
                filter(g -> !any(g == g′ for (i, g′) in x), items(V)))
    end

    R₁ = reduce(V, allocations...)
    R₂ = reduce(R₁.profile, α, greedy=false)
    return chain(R₁, R₂)
end


"""
    reduce(V::Additive, C::Counts{OrderedCategory}, α)

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

            V, C = R.profile, R.constraint
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
            R = reduce(V, i, [g])
            V = R.profile

            # Recursive application on the new instance
            R′ = reduce(V, α)

            # Combine the reductions
            return chain(R, R′)
        end
    end

    return Reduction(V)
end


"""
    order(V::Additive)

Create an ordered instance for the given weights. The weights are reordered for
each agent such that item 1 is worth the most and item m is worth the least.
Returns new additive valuations and a function to convert an allocation in the
ordered instance into one for the original instance.
"""
order(V::Additive) = Reduction(
        Additive(sort(V.values, dims=2, rev=true)),
        agents(V), items(V), A -> revert(V, A)
    )


"""
    revert(V::Additive, A)

Convert an allocation for the ordered instance to one for the original instance.
"""
function revert(V::Additive, A)
    A′ = Allocation(na(A), ni(A))

    M = Set(items(V))
    for g in items(V)
        i = owner(A, g)

        # Select the remaining good with the highest value
        g′ = maximum(g -> (value(V, i, g), g), M)[2]

        give!(A′, i, g′)
        M = setdiff(M, g′)
    end

    return A′
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

        push!(Co, OrderedCategory(m, length(c), c.threshold))
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
