"""
    alloc_half_mms(V::Additive, C::Counts)

Find a 1/2-approximate MMS allocation that obeys the constraints imposed by C.
First the instance is reduced to an ordered normalized instance where each good
is worth less than 1/2. While there are more than one agent remaining, the
algorithm creates a bundle with the ``⌊length(category)/n⌋`` lowest-valued
items in each category. Repeatedly, it converts each of these to the
highest-valued remaining item in the category until it either runs out of items
to convert or an agent values the bundle at least 1/2. If the procedure runs
out of items to convert, it adds the highest-valued remaining item in each
category, in order, to get ``⌈length(category)/n⌉`` items from each category.
After each such item is added, the value is again checked for each agent. Since
the instance was ordered normalized and without items worth 1/2 or more, the
bundle created will always be worth more than 1/2 to one of the remaining
agents before the procedure runs out of items to add to it or convert from low-
to high-valued.
"""
function alloc_half_mms(V::Additive, C::Counts)
    V, C, convert = order(V, C)

    # Normalize and allocate items worth 0.5 or more.
    V, C, convert2 = reduce(V, C, 0.5)

    N, n = agents(V), na(V)
    converts = [convert, convert2]

    while n > 1
        # Fill the bundle with the floor(k_h/n) least valuable items
        B = union(Set{Int}(), [c[end - floor_n(c, n) + 1:end] for c in C]...)

        # Convert lower value items to higher value items
        categories = Iterators.Stateful(C)
        c, converted = popfirst!(categories), 0
        while all(i -> value(V, i, B) < 0.5, N)
            if converted == floor_n(c, n)
                isempty(categories) && break

                c, converted = popfirst!(categories), 0
                continue
            end

            # Convert by removing the lowest valued remaining item and adding
            # the  highest valued item not in the bundle.
            converted += 1
            setdiff!(B, Set(c[end-converted+1]))
            B = B ∪ c[converted]
        end

        # Add higher value items
        categories = Iterators.Stateful(C)
        while all((i) -> value(V, i, B) < 0.5, N)
            c = popfirst!(categories)

            if length(c) % n > 0
                B = B ∪ c[ceil_n(c, n)]
            end
        end

        V, C, convert = reduce(V, C, findfirst(i -> value(V, i, B) >= 0.5, N), B)
        push!(converts, convert)
        N, n = agents(V), na(V)
    end

    # Use an empty allocation if there are no more agents remaining
    if n == 0
        A = Allocation(0, 0)
    else
        A = Allocation(1, ni(V))
        give!(A, 1, items(V))
    end

    # Apply all the collected converters in order
    for convert in Iterators.reverse(converts)
        A = convert(A)
    end

    return A
end

