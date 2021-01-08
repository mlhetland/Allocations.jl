"""
    reduce(V::Additive, C::Vector{OrderedCategory}, agent, removedbundle)

Reduce the instance given by the pair (V, C) to a new instance by giving the
supplied agent the supplied bundle, `B`. Returns an additive valuation, a set
of ordered categories for the new reduced instance and a function that turns an
allocation in the reduced instance into one for the original instance,
including giving the supplied agent the supplied bundle.
"""
function reduce(V::Additive, C::Vector{OrderedCategory}, agent, B)
    N, M = agents(V), items(V)
    n′, m′ = na(V) - 1, ni(V) - length(B)

    Vs = zeros(n′, m′)
    translate = zeros(Int, m′)

    Δg = 0

    for g in M
        if g in B 
            Δg += 1
            continue
        end

        g′ = g - Δg
        translate[g′] = g

        Vs[1:n′,g′] = [value(V, i, g) for i in N if i != agent]
    end

    # Create new ordered categories
    Cs = OrderedCategory[]
    index = 1
    for c in C
        newlength = length(c ∩ translate)
        push!(Cs, OrderedCategory(index, newlength, c.threshold))
        index += newlength
    end

    return Additive(Vs), Cs, (A) -> revert(translate, agent, B, A)
end


"""
    revert(translate::Vector{Int}, agent, removedbundle, A::Allocation)

Convert an allocation for a reduced instance to one for the original instance,
including giving the removed bundle, `B`, to the removed agent.
"""
function revert(translate::Vector{Int}, agent, B, A::Allocation)
    A′ = Allocation(na(A) + 1, ni(A) + length(B))
    for i in 1:na(A)
        i′ = i + (i >= agent)
        give!(A′, i′, [translate[g] for g in bundle(A, i)])
    end

    give!(A′, agent, B)

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
function reduce(V::Additive, C::Vector{OrderedCategory}, α::Float64)
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

            V, C, convert = reduce(V, C, i, bundle)

            # Recursive application, as the removed items may cause items worth
            # less than α to be worth α or more after a new normalization.
            V, C, prev_convert = reduce(V, C, α)

            # The converters have to be applied in reverse.
            return V, C, (A) -> convert(prev_convert(A))
        end
    end

    return V, C, (A) -> A
end


"""
    order(V::Additive, C::Counts)

Create an ordered instance for the given weights and categories. The items are
reorded such that each category has a continous range of indices for its items.
Returns new additive valutations, an array of `OrderedCategory` objects and
a function that converts an allocation in the ordered instance to one for the
original instance.
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

    return Additive(Vo), Co, (A) -> revert(V, C, Co, A)
end


"""
    revert(V::Additive, C::Counts, Co::Vector{OrderedCategory}, A)

Convert an allocation for the ordered instance to one for the original instance.
"""
function revert(V::Additive, C::Counts, Co::Vector{OrderedCategory}, A)
    A′ = Allocation(na(A), ni(A))
 
    for (unordered_category, ordered_category) in zip(C, Co)
        items = copy(unordered_category.members)
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
