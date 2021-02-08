"""
    alloc_rand(V)

A straightforward lottery that allocates the items randomly to the agents. For
each item, its agent is selected uniformly at random. The valuation `V` is not
used, other than to determine the number of agents and items.
"""
alloc_rand(V) = alloc_rand(na(V), ni(V))


"""
    alloc_rand(n::Int, m::Int)

Same as `alloc_rand(V)`, for `n` agents and `m` items.
"""
function alloc_rand(n::Int, m::Int)

    A = Allocation(n, m)

    N, M = agents(A), items(A)

    for g in items(A)
        give!(A, rand(N), g)
    end

    return (alloc=A,)

end


"""
    alloc_rand(V, C::Conflicts)

Allocate items to agents randomly, respecting the item conflicts. Uses the
[randomized coloring procedure with
symmetry-breaking](https://doi.org/10.1007/978-3-540-70575-8_26) of Pemmaraju
and Srinivasan, which works as follows:

1. Give the items random priorities, corresponding to a permutation selected
   uniformly at ramdom.
2. Tentatively allocate each item randomly to an agent, without concern for the
   item conflicts.
3. If an agent has received conflicting items, it keeps the highest-priority
   item (i.e., earliest in the permutation), and the others are reallocated
   arbitrarily.

This final arbitrary reallocation is also performed randomly in this
implementation, by going through the items in random order, allocating each to a
randomly selected agent among those able to receive it.

The valuation `V` is not used, other than to determine the number of agents and
items.

For this algorithm to function properly, the maximum degree of the conflict
graph should be strictly less than the number of agents.
"""
alloc_rand(V, C::Conflicts) = alloc_rand(na(V), ni(V), C)


"""
    alloc_rand(n::Int, m::Int, C::Conflicts)

Same as `alloc_rand(V, C)`, for `n` agents and `m` items.
"""
function alloc_rand(n::Int, m::Int, C::Conflicts)

    G = graph(C)

    @assert Δ(G) < n

    # Random permutation of the items, for priority:
    π = randperm(m)

    # Tentative allocation:
    A = alloc_rand(n, m).alloc

    removed = Set{Int}()

    # In every conflict, remove lower-priority items:
    for e in edges(G)

        g, h = src(e), dst(e)

        (owned(A, g) && owned(A, h)) || continue

        i, j = owner(A, g), owner(A, h)

        i == j || continue

        g = π[g] < π[h] ? h : g

        deny!(A, i, g)
        push!(removed, g)

    end

    # Allocate each item to some agent able to receive it:
    for g in shuffle(collect(removed))

        conflicted = neighbors(G, g)

        for i in randperm(n)

            isempty(conflicted ∩ bundle(A, i)) || continue

            give!(A, i, g)

            break

        end

    end

    return (alloc=A,)

end
