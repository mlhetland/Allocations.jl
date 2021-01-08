"""
    alloc_half_mms(V::Additive, C::Counts)

Find a 1/2-approximate MMS allocation that obeys the constraints imposed by C.
"""
function alloc_half_mms(V::Additive, C::Counts)
    V, C, convert = order(V, C)

    # Normalize and allocate items worth 0.5 or more.
    V, C, convert2 = reduce(V, C, 0.5)

    return convert(convert2(alloc_half_mms_bag_filling(V, C)))
end




"""
    alloc_half_mms_bag_filling(V::Additive, C::Vector{OrderedCategory})

Create a valid 1/2-MMS allocation from a normalized ordered instance that does
not contain any items worth 1/2 or more to any agent. The algorithm creates a
bundle of the ``⌊length(category)/n⌋`` lowest-valued items in each category.
Repeatedly, it converts each of these to the highest-valued remaining item in
the category until it either runs out of items to convert or an agent values
the bundle at least 1/2. If the procedure runs out of items to convert, it adds
the highest-valued remaining item in each category, in order, to get
``⌈length(category)/n⌉`` items from each category. After each such item is
added, the value is again checked for each agent. As long as the initial
ordered instance was normalized and without items worth 1/2 or more, the bundle
created will always be worth more than 1/2 to one of the remaining agents
before the procedure runs out of items to add to it or convert from low- to 
high-valued.
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

    V, C, convert = reduce(V, C, findfirst(i -> value(V, i, bundle) >= 0.5, N), bundle)
    return convert(alloc_half_mms_bag_filling(V, C))
end
