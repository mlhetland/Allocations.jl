using DataStructures

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


"""
    alloc_rand(V)

A straightforward lottery that allocates the items randomly to the agents. For
each item, its agent is selected uniformly at random. The valuation `V` is not
used, other than to determine the number of agents and items. The return value
is a named tuple with the field `alloc` (the `Allocation`).
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


"""
    alloc_bkv18_1(V)

The first algorithm (**Alg-Identical**) described by Barman, Krishnamurty and
Vaish in their 2018 paper [Greedy Algorithms for Maximizing Nash Social
Welfare](https://doi.org/10.1145/3355902). The algorithm finds a
1.061-approximate MNW allocation when agents have identical valuations, i.e.,
for any agents `i`, `j` and item `g`, `value(V, i, g) == value(V, j, g)`. (This
approximation ratio applies to the geometric mean of agent utilities, not the
raw product.) The result will also be envy-free up to any item (EFX).

The algorithm follows a straightforward greedy allocation procedure, where in
each iteration, the most valuable item is allocated to the agent with the lowest
utility.
"""
function alloc_bkv18_1(V)

    n, m = na(V), ni(V)
    v = [value(V, 1, g) for g = 1:m]

    pq = PriorityQueue{Int, Int}(i => 0 for i = 1:n)

    A = Allocation(n, m)

    for g in sortperm(v, rev=true)
        i, u = peek(pq)
        give!(A, i, g)
        pq[i] = u + v[g]
    end

    return (alloc = A,)

end


# The original descriptions of the following algorithm leaves a great latitude
# in how one actually implements it. The gist of the procedure is a local
# search, where an allocation is incrementally improved along paths in a graph,
# which has edges added or removed as items are reallocated. Here one might,
# e.g., use one of the structures developed for dynamically updating transitive
# closures. Because the graph in question is one with the agents as its
# vertices, which means it will probably be quite small, any asymptotic
# improvements could quite possibly be swamped by overhead. For now, things are
# kept simple, with a fresh traversal for each agent pair.

"""
    alloc_bkv18_2(V; randpri=true, complete=false)
    alloc_hpps20_1(V; randpri=true, complete=false) # alias

The second algorithm (**Alg-Binary**) described by Barman, Krishnamurty and
Vaish in their 2018 paper [Greedy Algorithms for Maximizing Nash Social
Welfare](https://doi.org/10.1145/3355902). The algorithm finds MNW allocations
in polynomial time for binary additive valuations, i.e., where each agent values
any given object at 0 or 1 (e.g., an `Additive{BitMatrix}`). It also works in a
more general setting, where `value(V, i, S)`, for any given `i`, is a concave
function of the number of items `g` in `S` for which `value(V, i, g) == 1`.

The original algorithm builds on an initial allocation, but does not specify
what this allocation should be. It also does not deal with the case where one or
more agents ends up with zero utility; in fact the procedure will not work even
if we start with two or more agents with zero utility in the intial allocation.
The strategy followed here is the same as that of Caragiannis et al.
(https://doi.org/10.1145/3355902), where a maximum cardinality set of agents
achieving positive utility is found using bipartite matching (with no fairness
considerations). The remaining items are randomly allocated to agents among
these that value them, if any. Remaining agents and items are ignored by the
procedure.

Following the algorithm of Barman et al., the tie-breaking procedure (Algorithm
1) of Halpern et al. (http://arxiv.org/abs/2007.06073) is used, where the MNW
allocation is transformed into the lexically greatest MNW, according to some
ordering of the agents, providing group-strategyproofness (GSP) in addition to
the EF1 and PO guarantees that follow from MNW. By default, the agent
ordering/priority is random; if this randomization is turned off, the default
ordering is used, with agent `1` receiving the highest priority, etc.

!!! note

    Despite the use of randomization here, by default, this is the
    *deterministic* procedure of Halpern et al. They also describe a randomized
    procedure, which functions in an entirely different manner.

Finally, if the `complete` argument is set to `true`, the allocation is
completed with `fill_even!` (which means that some agents that must necessarily
get a utility of zero can still receive items valued zero, if that evens out the
bundle cardinalities). Note that this undermines the GSP guarantee, which
requires that these items be discarded. The return value is a named tuple with
the fields `alloc` (the `Allocation`) and `mnw` (the Nash welfare, ignoring
agents with zero utility).
"""
function alloc_bkv18_2(V; randpri=true, complete=false)

    n, m = na(V), ni(V)

    # Tentative allocation, where O[g] is the owner of g. Will assign positive
    # utility to as many agents as possible (the set N, found through maximum
    # bipartite matching betewen agents and items they value), and will allocate
    # items only to them.
    O = zeros(Int, m)

    N = Int[]
    for (i, g) in bipartite_matching(V)
        O[g] = i
        push!(N, i)
    end

    M = Int[]
    for g in 1:m
        # To simplify the logic below, and ideally to reduce the number of
        # iterations, we assign items only to agents who value them, if any. If
        # no agents value an item, it will not participate in the procedure
        # anyway, and is left unallocated until the end of the function.
        fans = [i for i in N if value(V, i, g) > 0]
        isempty(fans) && continue # Leave g unallocated for now
        push!(M, g)
        O[g] == 0 && (O[g] = rand(fans))
    end

    # U[i]: The utility of agent i
    U = zeros(Int, n)
    for g = M
        i = O[g]
        U[i] += value(V, i, g)
    end

    # Pre-allocation
    G = zeros(Int, n, n)
    Π = zeros(Int, n, n)

    # Update the graph G based on current ownership.
    function update()

        # A form of adjacency matrix, where any G[i, j] > 0 is an arbitrary item
        # we may transfer from j to i.
        G .= 0
        for i in N, g in M
            (value(V, i, g) == 0 || O[g] == i) && continue
            G[i, O[g]] = g
        end

    end

    # We should never need more than this (as shown by Barman, Krishnamurty and
    # Vaish). We'll add one iteration, which should never be reached, as a
    # failsafe (with an assertion, below the loop):
    maxiter = ceil(Int, 2m*(n+1)*log(n*m))

    it = 0
    for outer it = 1:maxiter + 1

        update()

        # Predecessor matrix, where any Π[i, j] > 0 is the predecessor of j
        # along the path (sequence of possible transfers) from i.
        Π .= 0
        for i in N, j in N
            G[i, j] > 0 && (Π[i, j] = i)
        end

        # Straightforward version of Warshall's algorithm, using only Π and
        # restricting nodes of interest to N, i.e., the agents with nonzero
        # utility.
        for k in N, i in N, j in N
            (i == j || Π[i, j] > 0) && continue
            if Π[i, k] > 0 && Π[k, j] > 0
                Π[i, j] = Π[k, j]
            end
        end

        # The augmentation procedure of Barman et al., except that we only care
        # about agents with nonzero utility. Only items that are valued by at
        # least one agent will be used in Π. Such items are initially allocated
        # to agents that value them, and will remain so after each update, which
        # means that transfering the item will necessarily increment U[i] and
        # decrement U[j].
        s = t = 0
        best = 0
        for i in N, j in N
            @assert U[i] ≥ 1 && U[j] ≥ 1
            (Π[i, j] > 0 && U[j] > 1) || continue
            change = (U[i] + 1)*(U[j] - 1)/(U[i]*U[j])
            if change > best
                best = change
                s, t = i, j
            end
        end

        # No improvement, so we're done
        best > 1 || break

        # There was an improvement, and the transfer of one unit of utility.
        U[t] -= 1; U[s] += 1

        # Pass items back along path from s to t
        while t ≠ s
            p = Π[s, t]
            @assert O[G[p, t]] == t
            O[G[p, t]] = p
            t = p
        end

    end

    @assert it ≤ maxiter # Make sure we found optimum

    # We now have an MNW allocation, but want to applie the tiebreaker from
    # HPPS20(1), ending up with the maximum MNW allocation in lexicographic
    # order, given some priority ordering of the agents.

    # Normally, the priorities are randomized (and function primarily to prevent
    # group strategyproofness), but that may be turned off, in which case agents
    # are prioritized by their index.
    if randpri
        shuffle!(N) # Random priorities
    else
        sort!(N)    # If i < j, i is prioritized
    end

    # Find a path (in the same graph as before) from i (with priority p) to some
    # node j with lower priority and utility U[i] + 1. If one is found, pass
    # items back to increment U[i] and decrement U[j]. (For simplicity, we
    # increment/decrement along the path as well, though that cancels out.)
    #
    # We could have use Warshall here as well, but we'll potentially be
    # modifying the graph between every start node, so a straight DFS is
    # probably less wasteful.
    function search(p, i, u, seen=falses(n))
        seen[i] = true
        for (q, j) in enumerate(N)
            (G[i, j] > 0 && !seen[j]) || continue
            target = (q > p && U[j] == u + 1)
            if target || search(q, j, u, seen)
                @assert O[G[i, j]] == j
                O[G[i, j]] = i
                U[j] -= 1; U[i] += 1
                return true
            end
        end
        return false
    end

    # In order of priority, p:
    for (p, i) in enumerate(N)
        update()
        search(p, i, U[i])
    end

    # We now have the appropriate MNW allocation, represented by the ownership
    # vector O, and create the resulting Allocation.

    A = Allocation(n, m)
    for g in M
        give!(A, O[g], g)
    end

    if complete
        # Distribute the unvalued items as well, somewhat evenly. Note that this
        # ruins the group strategyproofness guarantees of HPPS20.
        fill_even!(A)
    end

    return (alloc=A, mnw=nash_welfare(V, A))

end

# Rather than BKV18(2) with HPPS20(1) for tie-breaking, one can view the
# function as HPPS20(1) with BKV18(2) for finding the initial allocation.
const alloc_hpps20_1 = alloc_bkv18_2
