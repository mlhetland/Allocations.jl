using DataStructures

"""
    alloc_half_mms(V::Additive, C::Counts; α=0.5)

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

Another approximation ratio, `α`, can be supplied. If `α ≤ 0.5` the algorithm is
guaranteed to succeed. Otherwise, the method will try to find an allocation with
an approximation ratio of `α`, but may fail. In the latter case, the results
will indicate a failure by setting `res.fail` to `true`.
"""
function alloc_half_mms(V::Additive, C::Counts; α=0.5)
    R = order(V, C)

    # Normalize and allocate items worth α or more.
    R′ = reduce(valuations(R), constraints(R), α)

    R = chain(R, R′)

    V, C = valuations(R), constraints(R)
    N, n = agents(V), na(V)

    while n > 1
        # Fill the bundle with the floor(k_h/n) least valuable items
        B = union(Set{Int}(), [c[end - floor_n(c, n) + 1:end] for c in C]...)

        # Convert lower value items to higher value items
        categories = Iterators.Stateful(C)
        c, converted = popfirst!(categories), 0
        while all(i -> value(V, i, B) < α, N)
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
        while all((i) -> value(V, i, B) < α, N)
            # If there are no more categories, the supplied α does not work
            if isempty(categories)
                return (alloc=nothing, fail=true)
            end

            c = popfirst!(categories)

            if length(c) % n > 0
                B = B ∪ c[ceil_n(c, n)]
            end
        end

        R′ = reduce(V, C, findfirst(i -> value(V, i, B) ≥ α, N), B)
        R = chain(R, R′)
        V, C = valuations(R), constraints(R)
        N, n = agents(V), na(V)
    end

    # Use an empty allocation if there are no more agents remaining
    if n == 0
        A = Allocation()
    else
        # If there is not enough value remaining for the last agent, the
        # supplied α does not work
        if value(V, 1, items(V)) < α
            return (alloc=nothing, fail=true)
        end

        A = Allocation(V)
        give!(A, 1, items(V))
    end

    # Create an allocation in the original instance
    return (alloc=transform(R, A), fail=false)

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
    alloc_bkv18_1(V; randpri=true)

The first algorithm (**Alg-Identical**) described by Barman, Krishnamurty and
Vaish in their 2018 paper [Greedy Algorithms for Maximizing Nash Social
Welfare](https://doi.org/10.1145/3355902). The algorithm finds a
1.061-approximate MNW allocation when agents have identical valuations, i.e.,
for any agents `i`, `j` and item `g`, `value(V, i, g) == value(V, j, g)`. (This
approximation ratio applies to the geometric mean of agent utilities, not the
raw product.) The result will also be envy-free up to any item (EFX).

The algorithm follows a straightforward greedy allocation procedure, where in
each iteration, the most valuable item is allocated to the agent with the lowest
utility. By default, ties are broken by giving the agents random priorities; if
`randpri` is set to false, they are instead broken lexicographically (as
specified by Barman et al.), so that the agent with the lower index is
preferred.
"""
function alloc_bkv18_1(V; randpri=true)

    N, M = agents(V), items(V)
    v = [value(V, 1, g) for g in M]

    π = randpri ? shuffle(N) : N

    pq = PriorityQueue(i => (0, p) for (i, p) in zip(N, π))

    A = Allocation(V)

    for g in sortperm(v, rev=true)
        i, (u, p) = peek(pq)
        give!(A, i, g)
        pq[i] = (u + v[g], p)
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
allocation is transformed into the lexicographically greatest MNW, according to
some ordering of the agents, providing group-strategyproofness (GSP) in addition
to the EF1 and PO guarantees that follow from MNW. By default, the agent
ordering/priority is random; if this randomization is turned off (with `randpri`
set to false), the default ordering is used, with agent `1` receiving the
highest priority, etc.

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


"""
    alloc_ghss18_4(V::Submodular, MMSs)

The fourth algorithm (**Algorithm 4**) described by Ghodsi et al. in the 2018
paper [Fair allocation of Indivisible Goods: Improvements and
Generalizations](https://arxiv.org/abs/1704.00222). The algorithm finds a
1/3-approximate MMS allocation for a given submodular instance and corresponding
maximin shares for the agents (`MMSs[i]` should be the MMS of agent `i`). If the
supplied maximin shares, are higher than the actual maximin shares, the method
may fail. In that case, this will be indicated in the result, where `res.fail`
will be set to true and `res.agent` will be set to the agent last considered
when the method failed to improve. If the maximin shares are unknown, use
`alloc_ghss18_4b`.
"""
function alloc_ghss18_4(V::Submodular, MMSs)

    @assert na(V) == length(MMSs) "There must be one MMS per agent"

    # Scale valuations so that each agent's MMS is 1
    Vo = V
    V = Submodular([B -> (value(Vo, i, B) / MMSs[i]) for i in agents(Vo)], ni(Vo))

    # Allocate all items worth at least 1/3
    R = reduce(V, 1/3)
    V = valuations(R)

    A = Allocation(V)

    # In case all agents have received an item worth at least 1/3
    if na(V) == 0
        A = transform(R, A)

        # The algorithm does not specify what to do in this situation. We
        # therefore simply allocate the goods randomly to the agents.
        fill_random!(A)

        return (alloc=A, fail=false, agent=0)
    end

    # Create an arbitrary allocation with the remaining items by for each item
    # allocating it to a random agent
    fill_random!(A)

    N, M = agents(V), items(V)

    # All valuations from now on will be wrapped by the ceiling function for 2/3
    V′ = Submodular([B -> min(2/3, value(V, i, B)) for i in N], ni(V))

    # Continue until all agents have a bundle valued at least 1/3
    while any(i -> value(V′, i, bundle(A, i)) < 1/3, N)
        # Want to improve the bundle of the worst of agent
        i = argmin([value(V′, i, bundle(A, i)) for i in N])

        found = false
        for g in M
            j = owner(A, g)
            i == j && continue

            Bᵢ, Bⱼ = bundle(A, i), bundle(A, j)

            # Only the contribution of agents `i` and `j` to `ex(A)` (the
            # function defined by Ghodsi et al.) will change. Thus, we do not
            # need to compute the entire value and can simply check the
            # difference in the contributions of agents `i` and `j` before and
            # after moving `g`.
            v = value(V′, i, Bᵢ) + value(V′, j, Bⱼ)
            v′ = value(V′, i, Bᵢ ∪ g) + value(V′, j, symdiff(Bⱼ, g))

            # If the change in total value is more than 1/3m when moving g from
            # i to j, perform the move.
            if v′ ≥ v + 1/(3*ni(V))
                found = true
                deny!(A, j, g)
                give!(A, i, g)
                break
            end
        end

        # If no sufficently large improvement can be made, then the MMS values
        # are incorrect (too large)
        if !found
            return (alloc=nothing, fail=true, agent=item(R, i))
        end
    end

    # Convert the allocation to one for the original instance
    return (alloc=transform(R, A), fail=false, agent=0)

end


"""
    alloc_ghss18_4b(V::Submodular; a=3, x_warn=true)

A variation on the fourth algorithm (**Algorithm 4**) described by Ghodsi et al.
in the 2018 paper [Fair allocation of Indivisible Goods: Improvements and
Generalizations](https://arxiv.org/abs/1704.00222). The algorithm finds a
1/3-approximate MMS allocation for a given submodular instance. The method
starts by overestimating the MMS of each agent and slowly decreasing the MMS of
specific agents until `alloc_ghss18_4` returns an allocation.

The amount that the MMS of an agent should be reduced by in each iteration is
not specified by Ghodsi et al. One can show that if the factor is `1/(1 + 1/x)`,
where `x ≥ 3n - 1`, then the algorithm will successfully find a 1/3-approximate
MMS allocation. One way to show this, is to modify Lemma 4.6 in their paper to
assume that each of the bundles `Sᵢ` is valued at least `1/(1 + 1/x)`. Using
this modified version of Lemma 4.6, one can modify the proof of Theorem 4.7 to
show that as long as `x ≥ 3n - 1`, the change in expectance from moving an item
is at least `1/(3m)`. The value of `x` used in this implementation is `x = an`,
where the keyword argument `a` is set to `3` by default (i.e., `x = 3n`). If `a`
is set so that `x < 3n - 1` a warning will be given. The warning can be turned
off by setting `x_warn` to `false`.
"""
function alloc_ghss18_4b(V::Submodular; a=3, x_warn=true)

    if x_warn && a * na(V) < 3 * na(V) - 1
        @warn("The value of `a` may be too small")
    end

    MMSs = Vector{Float64}([value(V, i, items(V)) for i in agents(V)])

    res = alloc_ghss18_4(V, MMSs)
    while res.fail
        MMSs[res.agent] /= (1 + 1/(a * na(V)))
        res = alloc_ghss18_4(V, MMSs)
    end

    return (alloc=res.alloc,)

end


"""
    alloc_bb18_3(V::Additive, C::Counts; a=3, ghss18_4b_warn=true)

The 1/3-approximate MMS-allocation under cardinality constraints algorithm
(Section 5) described by Biswas and Barman in their 2018 paper [Fair Division
Under Cardinality Constraints] (https://doi.org/10.24963/ijcai.2018/13). Finds a
1/3-approximate MMS allocation for an instance of the fair allocation problem
under cardinality constraints by converting the additive instance under
cardinality constraints to a submodular instance without cardinality
constraints. The allocation is then found by using the method of Ghodsi et al.
(`alloc_ghhs18_4b`), with possible reallocation of items to satisfy the
constraints. Both keyword arguments, `a` and `ghss18_4b_warn`, are passed
directly to `alloc_ghhs18_4b` as respectively the keyword arguments `a` and
`x_warn`. See [`alloc_ghhs18_4b`](@ref) for documentation on how to use them.
"""
function alloc_bb18_3(V::Additive, C::Counts; a=3, ghss18_4b_warn=true)

    # The submodular valuation function construction for V and C
    function submodular_valuation(i, B)
        total = 0

        for c in C
            # Extract the overlap between B and the category. We do not care
            # about which items overlap, only their value.
            overlap = [value(V, i, g) for g in c.members ∩ B]

            # An agent cannot get value from more items than the threshold of
            # the given category.
            n_items = min(threshold(c), length(overlap))

            # Give the agent the `n_items` most valuable items in the overlap.
            total += sum(partialsort!(overlap, 1:n_items, rev=true))
        end

        return total
    end

    N = agents(V)

    # Create a submodular instance and use the method of Ghodsi et al.
    V′ = Submodular([B -> submodular_valuation(i, B) for i in N], ni(V))
    res = alloc_ghss18_4b(V′, a=a, x_warn=ghss18_4b_warn)
    A = res.alloc

    # Easy way to find overlap between agent's bundle and a category
    overlap(i, c) = bundle(A, i) ∩ c.members

    # It is possible that the 1/3-approximate MMS allocation produced by
    # `alloc_ghss18_4b` does not adhere to the cardinality constraints. If this
    # is the case, then we for each non-adhering bundle, give away the
    # corresponding agent's least-preferred items in the bundle from any
    # category for which the bundle breaks the category's threshold. This does
    # not result in a decrease in the agent's perceived valuation of their
    # bundle based on the submodular valuations.
    for i in N, c in C
        B = overlap(i, c)
        # If there are more items than allowed
        while length(B) > threshold(c)
            # Give away the item in c that `i` prefers the least
            g = minimum(g -> (value(V, i, g), g), B)[2]

            # Give the item away to any agent that is allowed more items from
            # the category
            i′ = findfirst(i′ -> length(overlap(i′, c)) < threshold(c), N)

            deny!(A, i, g)
            give!(A, i′, g)
            delete!(B, g)
        end
    end

    return (alloc=A,)

end
