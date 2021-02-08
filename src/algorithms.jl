"""
    alloc_rand(V)

A straithgtforward lottery that allocates the items randomly to the agents. For
each item, its agent is selected uniformly at random.
"""
function alloc_rand(V)

    A = Allocation(na(V), ni(V))

    N, M = agents(A), items(A)

    for g in items(A)
        give!(A, rand(N), g)
    end

    return (alloc=A,)

end

# [The Randomized Coloring Procedure with
# Symmetry-Breaking](https://doi.org/10.1007/978-3-540-70575-8_26)
"""
    alloc_rand(V, C::Conflicts)

"""
function alloc_rand(V, C::Conflicts)

    # ...

    # return (alloc=alloc,)

end
