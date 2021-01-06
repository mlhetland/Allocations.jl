"""
    half_mms(V::Additive, C::Counts)

Find a 1/2-approximate MMS allocation that obeys the constraints imposed by C.
"""
function half_mms(V::Additive, C::Counts)
    V, C, convert = create_ordered_instance(V, C)
    return convert(subprocedure1(V, C))
end


"""
    subprocedure1(V::Additive, C::Array{OrderedCategory, 1})

First part of `half_mms`. Takes an ordered instance, normalizes the values and
gives any agent that value an individual item greater than or equal to 1/2 the
item and any low value items required to reduce to a valid instance. After
performing this reduction, the function calls itself recursively. If there is
no such agent-item pairs, `subprocedure2` is called.
"""
function subprocedure1(V::Additive, C::Array{OrderedCategory, 1})
    N, n, M = agents(V), na(V), items(V)

    if n == 1
        allocation = Allocation(1, ni(V))
        for j in M
            give!(allocation, 1, j)
        end
        return allocation
    end

    V = normalize(V)

    # Find any agent who values a single item at least 0.5 and give that item
    # to the agent along with the least valuable items in the remaining bundles
    # that must be given to the agent to guarantee a valid instance.
    for i in N
        for j in M
            if value(V, i, j) >= 0.5
                bundle = union(Set{Int}(j), [category[end - required(category, n) + (j in category) + 1:end] for category in C]...)
                V, C, converter = reduce_instance(V, C, i, bundle)
                return converter(subprocedure1(V, C))
            end
        end
    end

    return subprocedure2(V, C)
end


"""
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
function subprocedure2(V::Additive, C::Array{OrderedCategory, 1})
    N, n = agents(V), na(V)

    if n == 1
        allocation = Allocation(1, ni(V))
        for j in items(V)
            give!(allocation, 1, j)
        end
        return allocation
    end

    # Fill the bundle with the floor(k_h/n) least valuable items
    bundle = union(Set{Int}(), [category[end - floor_n(category, n) + 1:end] for category in C]...)

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

    V, C, converter = reduce_instance(V, C, findfirst(i -> value(V, i, bundle) >= 0.5, N), bundle)
    return converter(subprocedure2(V, C))
end
