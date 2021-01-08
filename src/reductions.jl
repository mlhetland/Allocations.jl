"""
    reduce(V::Additive, C::Vector{OrderedCategory}, agent, removedbundle)

Reduce the instance given by the pair (V, C) to a new instance by giving the
supplied agent the supplied bundle. Returns an additive valuation, a set of
ordered categories for the new reduced instance and a function that turns an
allocation in the reduced instance into one for the original instance,
including giving the supplied agent the supplied bundle.
"""
function reduce(V::Additive, C::Vector{OrderedCategory}, agent, removedbundle)
    N, M = agents(V), items(V)
    n, m = na(V), ni(V)

    Vs = zeros(n - 1, m - length(removedbundle))
    translate = zeros(Int, m - length(removedbundle))

    itemcounterchange = 0

    for j in M
        if j in removedbundle
            itemcounterchange += 1
            continue
        end

        newj = j - itemcounterchange
        translate[newj] = j

        Vs[1:n-1,newj] = [value(V, i, j) for i in N if i != agent]
    end

    # Create new ordered categories
    Cs = OrderedCategory[]
    index = 1
    for c in C
        newlength = length(c ∩ translate)
        push!(Cs, OrderedCategory(index, newlength, c.threshold))
        index += newlength
    end

    return Additive(Vs), Cs, (A) -> revert(translate, agent, removedbundle, A)
end


"""
    revert(translate::Vector{Int}, agent, removedbundle, A::Allocation)

Convert an allocation for a reduced instance to one for the original instance,
including giving the removed bundle to the removed agent.
"""
function revert(translate::Vector{Int}, agent, removedbundle, A::Allocation)
    A′ = Allocation(na(A) + 1, ni(A) + length(removedbundle))
    for i in 1:na(A)
        new_i = i + (i >= agent)
        for j in bundle(A, i)
            give!(A′, new_i, translate[j])
        end
    end

    for j in removedbundle
        give!(A′, agent, j)
    end

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

    if n == 1 return reduce(V, C, i, M) end

    V = normalize(V)

    # Find any agent who values a single item at least α and give that item to
    # the agent along with the least valuable items in the remaining bundles
    # that must be given to the agent to guarantee a valid instance.
    for i in N, j in M
        if value(V, i, j) >= α
            bundle = union(
                Set{Int}(j),
                [c[end - required(c, n) + (j in c) + 1:end] for c in C]...)
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

    itemcounter = 1
    for c in C
        for i in N
            Vo[i,itemcounter:itemcounter+length(c)-1] = sort([value(V, i, j) for j in c], rev=true)
        end

        push!(Co, OrderedCategory(itemcounter, length(c), threshold(c)))
        itemcounter += length(c)
    end

    return Additive(Vo), Co, (A) -> revert(V, C, Co, A)
end


"""
    revert(V::Additive, C::Counts, Co::Vector{OrderedCategory}, A:Allocation)

Convert an allocation for the ordered instance to one for the original instance.
"""
function revert(V::Additive, C::Counts, Co::Vector{OrderedCategory}, A::Allocation)
    A′ = Allocation(na(A), ni(A))
 
    for (orig, new) in zip(C, Co)
        items = copy(orig.members)
        for j in new 
            i = owner(A, j)
            item = maximum((el) -> (value(V, i, el), el), items)
            give!(A′, i, item[2])
            items = setdiff(items, item[2])
        end
    end

    return A′
end
