"""
    alloc_half_mms(V::Additive, C::Counts)

Find a 1/2-approximate MMS allocation that obeys the constraints imposed by C.
"""
function alloc_half_mms(V::Additive, C::Counts)
    V, C, convert = create_ordered_instance(V, C)
    return convert(alloc_half_mms_big_item_reduction(V, C))
end


"""
    alloc_half_mms_big_item_reduction(V::Additive, C::Vector{OrderedCategory})

First part of `half_mms`. Takes an ordered instance, normalizes the values and
gives any agent that value an individual item greater than or equal to 1/2 the
item and any low value items required to reduce to a valid instance. After
performing this reduction, the function calls itself recursively. If there is
no such agent-item pairs, `alloc_half_mms_bag_filling` is called.
"""
function alloc_half_mms_big_item_reduction(V::Additive, C::Vector{OrderedCategory})
    N, n, M = agents(V), na(V), items(V)

    if n == 1
        A = Allocation(1, ni(V))
        for j in M
            give!(A, 1, j)
        end
        return A
    end

    V = normalize(V)

    # Find any agent who values a single item at least 0.5 and give that item
    # to the agent along with the least valuable items in the remaining bundles
    # that must be given to the agent to guarantee a valid instance.
    for i in N, j in M
		if value(V, i, j) >= 0.5
			bundle = union(Set{Int}(j), [c[end - required(c, n) + (j in c) + 1:end] for c in C]...)
			V, C, convert = reduce_instance(V, C, i, bundle)
			return convert(alloc_half_mms_big_item_reduction(V, C))
		end
    end

    return alloc_half_mms_bag_filling(V, C)
end


"""
	alloc_half_mms_bag_filling(V::Additive, C::Vector{OrderedCategory})

The second part of `half_mms`. Takes a normalized ordered instance without any
items worth 1/2 or more to any agent and produces a valid
1/2-MMS allocation. The algorithm creates a bundle of the ⌊length(category)/n⌋
lowest-valued items in each category. Repeatedly, it converts each of these to
the highest-valued remaining item in the category until it either runs out of
items to convert or an agent values the bundle at least 1/2. If the procedure
runs out of items to convert, it adds the highest-valued remaining
item in each category, in order, to get ⌈length(category)/n⌉ items from each category.
After each such item is added, the value is again checked for each agent. As
long as the initial ordered instance was normalized, the procedure will never
run out of items to add and have no agent value the bundle at least 1/2.
"""
function alloc_half_mms_bag_filling(V::Additive, C::Vector{OrderedCategory})
    N, n = agents(V), na(V)

    if n == 1
        A = Allocation(1, ni(V))
        for j in items(V)
            give!(A, 1, j)
        end
        return A 
    end

    # Fill the bundle with the floor(k_h/n) least valuable items
    bundle = union(Set{Int}(), [c[end - floor_n(c, n) + 1:end] for c in C]...)

    # Convert lower value items to higher value items
    categories = Iterators.Stateful(C)
    category, converted = popfirst!(categories), 0
    while all(i -> value(V, i, bundle) < 0.5, N)
        if converted == floor_n(category, n) 
            isempty(categories) && break

            category, converted = popfirst!(categories), 0
            continue
        end

        # Convert by removing the lowest valued remaining item and adding the
        # highest valued item not in the bundle.
        converted += 1
        setdiff!(bundle, Set(category[end-converted+1]))
        bundle = bundle ∪ category[converted]
    end

    # Add higher value items
    categories = Iterators.Stateful(C)
    while all((i) -> value(V, i, bundle) < 0.5, N)
        category = popfirst!(categories)

        if length(category) % n > 0
            bundle = bundle ∪ category[ceil_n(category, n)]
        end
    end

    V, C, convert = reduce_instance(V, C, findfirst(i -> value(V, i, bundle) >= 0.5, N), bundle)
    return convert(alloc_half_mms_bag_filling(V, C))
end
